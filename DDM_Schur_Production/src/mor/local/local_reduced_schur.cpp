#include "local_reduced_schur.hpp"

#include "ddm_schur/common.hpp"
#include "sipg_core.hpp"
#include "linear_solvers.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

namespace mor::local {

// Runtime state owned by the production interface solver. The sparse matrix
// is fixed for a constant Backward-Euler time step, so the expensive PARDISO
// factorization is constructed once and reused by all subsequent solves.
struct LocalReducedSchurSolver::AugmentedDirectData {
    std::unique_ptr<SubdomainDirectSolver> factor;
    std::vector<int> rankOffsets;
    int dimension = 0;
    int solveThreads = 1;
    int workThreads = 1;
    int solveCalls = 0;
    double couplingSymmetryRelativeError = 0.0;
    double validatedRelativeResidual =
        std::numeric_limits<double>::infinity();
};

namespace {

using Clock = std::chrono::steady_clock;

// Reused HBM subdomains store their large numerical payload only in the
// representative template. Instances retain their own global/interface maps.
const SubdomainModel& templatePayload(const Model& model,
                                      const SubdomainModel& instance)
{
    if (!instance.templateReused) {
        return instance;
    }
    for (const SubdomainModel& candidate : model.subdomains) {
        if (candidate.templateId == instance.templateId
            && !candidate.templateReused) {
            return candidate;
        }
    }
    throw std::runtime_error(
        "[Production Schur] Reused template payload is missing.");
}

struct UpperCsr {
    std::vector<int> rowPtr;
    std::vector<int> colInd;
    std::vector<double> values;
};

struct ColumnValue {
    int column = 0;
    double value = 0.0;
};

// Assemble the symmetric upper triangle of
//
//   [ A_GammaGamma       A_Gamma,r ]
//   [ A_r,Gamma          A_rr      ].
//
// The local interior coordinates remain separate by subdomain. Consequently
// two subdomains communicate only through shared physical interface rows, and
// their local block assembly can proceed concurrently. PARDISO later performs
// the global sparse elimination needed to enforce interface continuity.
UpperCsr buildAugmentedUpperCsr(const Model& model,
                                std::vector<int>& rankOffsets,
                                double& couplingSymmetryRelativeError,
                                int assemblyThreads)
{
    rankOffsets.assign(model.subdomains.size() + 1, 0);
    for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
        rankOffsets[slot + 1] = rankOffsets[slot]
            + model.subdomains[slot].rank;
    }
    if (rankOffsets.back() != model.totalLocalRank) {
        throw std::runtime_error(
            "[Augmented direct] Local rank offsets are inconsistent.");
    }

    const int dimension = model.interfaceDofs + model.totalLocalRank;
    if (dimension <= 0) {
        throw std::runtime_error("[Augmented direct] Empty reduced system.");
    }

    // Pass 1 counts raw entries. Duplicate contributions are allowed here and
    // are compacted after each row is sorted.
    std::vector<std::size_t> rowCounts(
        static_cast<std::size_t>(dimension), 0);
    for (const InterfaceEntry& entry : model.interfaceEntries) {
        if (entry.row < 0 || entry.row >= model.interfaceDofs
            || entry.column < 0 || entry.column >= model.interfaceDofs) {
            throw std::runtime_error(
                "[Augmented direct] Interface entry is out of range.");
        }
        if (entry.row <= entry.column && entry.value != 0.0) {
            ++rowCounts[static_cast<std::size_t>(entry.row)];
        }
    }

    double maximumCoupling = 0.0;
    double maximumCouplingDifference = 0.0;
    for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
        const SubdomainModel& local = model.subdomains[slot];
        const SubdomainModel& data = templatePayload(model, local);
        const int rankOffset = model.interfaceDofs + rankOffsets[slot];

        for (int row = 0; row < local.rank; ++row) {
            for (int column = row; column < local.rank; ++column) {
                const double value = 0.5 * (
                    data.reducedInterior[static_cast<std::size_t>(
                        row * local.rank + column)]
                    + data.reducedInterior[static_cast<std::size_t>(
                        column * local.rank + row)]);
                if (value != 0.0) {
                    ++rowCounts[static_cast<std::size_t>(rankOffset + row)];
                }
            }
        }

        for (std::size_t localGamma = 0;
             localGamma < local.interfaceIndices.size(); ++localGamma) {
            const int gamma = local.interfaceIndices[localGamma];
            if (gamma < 0 || gamma >= model.interfaceDofs) {
                throw std::runtime_error(
                    "[Augmented direct] Local interface index is out of range.");
            }
            for (int mode = 0; mode < local.rank; ++mode) {
                const double interfaceInterior =
                    data.reducedInterfaceInterior[
                        localGamma * static_cast<std::size_t>(local.rank)
                        + static_cast<std::size_t>(mode)];
                const double interiorInterface =
                    data.reducedInteriorInterface[
                        static_cast<std::size_t>(mode)
                            * local.localInterfaceDofs
                        + localGamma];
                maximumCoupling = std::max(maximumCoupling,
                    std::max(std::abs(interfaceInterior),
                             std::abs(interiorInterface)));
                maximumCouplingDifference = std::max(
                    maximumCouplingDifference,
                    std::abs(interfaceInterior - interiorInterface));
                if (0.5 * (interfaceInterior + interiorInterface) != 0.0) {
                    ++rowCounts[static_cast<std::size_t>(gamma)];
                }
            }
        }
    }
    couplingSymmetryRelativeError = maximumCoupling > 0.0
        ? maximumCouplingDifference / maximumCoupling : 0.0;

    std::vector<std::size_t> rawRowPtr(
        static_cast<std::size_t>(dimension + 1), 0);
    for (int row = 0; row < dimension; ++row) {
        rawRowPtr[static_cast<std::size_t>(row + 1)] =
            rawRowPtr[static_cast<std::size_t>(row)]
            + rowCounts[static_cast<std::size_t>(row)];
    }
    if (rawRowPtr.back()
        > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error(
            "[Augmented direct] CSR exceeds 32-bit PARDISO indices.");
    }

    // Pass 2 writes every raw contribution into its preallocated row.
    std::vector<ColumnValue> rawEntries(rawRowPtr.back());
    std::vector<std::size_t> cursor = rawRowPtr;
    const auto append = [&](int row, int column, double value) {
        if (value == 0.0) {
            return;
        }
        const std::size_t target = cursor[static_cast<std::size_t>(row)]++;
        rawEntries[target] = {column, value};
    };

    for (const InterfaceEntry& entry : model.interfaceEntries) {
        if (entry.row <= entry.column) {
            append(entry.row, entry.column, entry.value);
        }
    }
    for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
        const SubdomainModel& local = model.subdomains[slot];
        const SubdomainModel& data = templatePayload(model, local);
        const int rankOffset = model.interfaceDofs + rankOffsets[slot];

        for (int row = 0; row < local.rank; ++row) {
            for (int column = row; column < local.rank; ++column) {
                append(rankOffset + row, rankOffset + column, 0.5 * (
                    data.reducedInterior[static_cast<std::size_t>(
                        row * local.rank + column)]
                    + data.reducedInterior[static_cast<std::size_t>(
                        column * local.rank + row)]));
            }
        }
        for (std::size_t localGamma = 0;
             localGamma < local.interfaceIndices.size(); ++localGamma) {
            const int gamma = local.interfaceIndices[localGamma];
            for (int mode = 0; mode < local.rank; ++mode) {
                append(gamma, rankOffset + mode, 0.5 * (
                    data.reducedInterfaceInterior[
                        localGamma * static_cast<std::size_t>(local.rank)
                        + static_cast<std::size_t>(mode)]
                    + data.reducedInteriorInterface[
                        static_cast<std::size_t>(mode)
                            * local.localInterfaceDofs
                        + localGamma]));
            }
        }
    }
    for (int row = 0; row < dimension; ++row) {
        if (cursor[static_cast<std::size_t>(row)]
            != rawRowPtr[static_cast<std::size_t>(row + 1)]) {
            throw std::runtime_error(
                "[Augmented direct] CSR row count mismatch.");
        }
    }

    // Rows are independent: sorting is parallel, while deterministic row
    // order preserves reproducible cache and timing diagnostics.
#pragma omp parallel for num_threads(assemblyThreads) if(assemblyThreads > 1) schedule(static)
    for (int row = 0; row < dimension; ++row) {
        auto begin = rawEntries.begin() + static_cast<std::ptrdiff_t>(
            rawRowPtr[static_cast<std::size_t>(row)]);
        auto end = rawEntries.begin() + static_cast<std::ptrdiff_t>(
            rawRowPtr[static_cast<std::size_t>(row + 1)]);
        std::sort(begin, end,
            [](const ColumnValue& left, const ColumnValue& right) {
                return left.column < right.column;
            });
    }

    std::vector<std::size_t> compactCounts(
        static_cast<std::size_t>(dimension), 0);
    for (int row = 0; row < dimension; ++row) {
        const std::size_t begin = rawRowPtr[static_cast<std::size_t>(row)];
        const std::size_t end = rawRowPtr[static_cast<std::size_t>(row + 1)];
        std::size_t read = begin;
        std::size_t write = begin;
        while (read < end) {
            const int column = rawEntries[read].column;
            double value = 0.0;
            while (read < end && rawEntries[read].column == column) {
                value += rawEntries[read].value;
                ++read;
            }
            if (value != 0.0) {
                if (!std::isfinite(value) || column < row
                    || column >= dimension) {
                    throw std::runtime_error(
                        "[Augmented direct] Invalid CSR entry.");
                }
                rawEntries[write++] = {column, value};
            }
        }
        compactCounts[static_cast<std::size_t>(row)] = write - begin;
    }

    UpperCsr result;
    result.rowPtr.assign(static_cast<std::size_t>(dimension + 1), 0);
    for (int row = 0; row < dimension; ++row) {
        const std::size_t next = static_cast<std::size_t>(
            result.rowPtr[static_cast<std::size_t>(row)])
            + compactCounts[static_cast<std::size_t>(row)];
        if (next > static_cast<std::size_t>(
                std::numeric_limits<int>::max())) {
            throw std::runtime_error(
                "[Augmented direct] Compacted CSR exceeds 32-bit indices.");
        }
        result.rowPtr[static_cast<std::size_t>(row + 1)] =
            static_cast<int>(next);
    }
    result.colInd.resize(static_cast<std::size_t>(result.rowPtr.back()));
    result.values.resize(result.colInd.size());

#pragma omp parallel for num_threads(assemblyThreads) if(assemblyThreads > 1) schedule(static)
    for (int row = 0; row < dimension; ++row) {
        const std::size_t source = rawRowPtr[static_cast<std::size_t>(row)];
        const std::size_t target = static_cast<std::size_t>(
            result.rowPtr[static_cast<std::size_t>(row)]);
        const std::size_t count = compactCounts[static_cast<std::size_t>(row)];
        for (std::size_t entry = 0; entry < count; ++entry) {
            result.colInd[target + entry] = rawEntries[source + entry].column;
            result.values[target + entry] = rawEntries[source + entry].value;
        }
    }
    return result;
}

double vectorDot(const std::vector<double>& left,
                 const std::vector<double>& right)
{
    if (left.size() != right.size()) {
        throw std::runtime_error("[Production Schur] Dot-product size mismatch.");
    }
    double value = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        value += left[i] * right[i];
    }
    return value;
}

} // namespace

DenseSymmetricFactor factorDenseSymmetric(
    const std::vector<double>& matrix, int size)
{
    if (size <= 0
        || matrix.size() != static_cast<std::size_t>(size * size)) {
        throw std::runtime_error(
            "[Local ROM] Invalid dense symmetric factor dimensions.");
    }

    DenseSymmetricFactor factor;
    factor.size = size;
    factor.lower.assign(matrix.size(), 0.0);
    double diagonalScale = 0.0;
    for (int row = 0; row < size; ++row) {
        diagonalScale = std::max(diagonalScale,
            std::abs(matrix[static_cast<std::size_t>(row * size + row)]));
    }
    const double threshold = 128.0 * std::numeric_limits<double>::epsilon()
        * std::max(std::numeric_limits<double>::min(), diagonalScale);

    // Prefer Cholesky because the projected thermal operator should be SPD.
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column <= row; ++column) {
            double value = 0.5 * (
                matrix[static_cast<std::size_t>(row * size + column)]
                + matrix[static_cast<std::size_t>(column * size + row)]);
            for (int k = 0; k < column; ++k) {
                value -= factor.lower[static_cast<std::size_t>(
                    row * size + k)]
                    * factor.lower[static_cast<std::size_t>(
                        column * size + k)];
            }
            if (row == column) {
                if (!(value > threshold) || !std::isfinite(value)) {
                    factor.cholesky = false;
                    break;
                }
                factor.lower[static_cast<std::size_t>(
                    row * size + column)] = std::sqrt(value);
            } else {
                factor.lower[static_cast<std::size_t>(
                    row * size + column)] = value
                    / factor.lower[static_cast<std::size_t>(
                        column * size + column)];
            }
        }
        if (!factor.cholesky) {
            break;
        }
    }
    if (factor.cholesky) {
        return factor;
    }

    // LDL^T is retained as a defensive path for small projected matrices
    // whose roundoff-level asymmetry defeats a strict Cholesky pivot test.
    factor.lower.assign(matrix.size(), 0.0);
    factor.diagonal.assign(static_cast<std::size_t>(size), 0.0);
    for (int row = 0; row < size; ++row) {
        factor.lower[static_cast<std::size_t>(row * size + row)] = 1.0;
        double pivot = matrix[static_cast<std::size_t>(row * size + row)];
        for (int k = 0; k < row; ++k) {
            const double entry = factor.lower[static_cast<std::size_t>(
                row * size + k)];
            pivot -= entry * entry
                * factor.diagonal[static_cast<std::size_t>(k)];
        }
        if (!(std::abs(pivot) > threshold) || !std::isfinite(pivot)) {
            throw std::runtime_error(
                "[Local ROM] Dense factor encountered a tiny pivot.");
        }
        factor.diagonal[static_cast<std::size_t>(row)] = pivot;
        for (int next = row + 1; next < size; ++next) {
            double value = 0.5 * (
                matrix[static_cast<std::size_t>(next * size + row)]
                + matrix[static_cast<std::size_t>(row * size + next)]);
            for (int k = 0; k < row; ++k) {
                value -= factor.lower[static_cast<std::size_t>(
                    next * size + k)]
                    * factor.lower[static_cast<std::size_t>(
                        row * size + k)]
                    * factor.diagonal[static_cast<std::size_t>(k)];
            }
            factor.lower[static_cast<std::size_t>(next * size + row)] =
                value / pivot;
        }
    }
    return factor;
}

void solveDenseSymmetric(const DenseSymmetricFactor& factor,
                         std::vector<double>& rightHandSide)
{
    const int size = factor.size;
    if (rightHandSide.size() != static_cast<std::size_t>(size)) {
        throw std::runtime_error("[Local ROM] Dense RHS has the wrong size.");
    }
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column < row; ++column) {
            rightHandSide[static_cast<std::size_t>(row)] -=
                factor.lower[static_cast<std::size_t>(row * size + column)]
                * rightHandSide[static_cast<std::size_t>(column)];
        }
        if (factor.cholesky) {
            rightHandSide[static_cast<std::size_t>(row)] /=
                factor.lower[static_cast<std::size_t>(row * size + row)];
        }
    }
    if (!factor.cholesky) {
        for (int row = 0; row < size; ++row) {
            rightHandSide[static_cast<std::size_t>(row)] /=
                factor.diagonal[static_cast<std::size_t>(row)];
        }
    }
    for (int row = size - 1; row >= 0; --row) {
        for (int column = row + 1; column < size; ++column) {
            rightHandSide[static_cast<std::size_t>(row)] -=
                factor.lower[static_cast<std::size_t>(column * size + row)]
                * rightHandSide[static_cast<std::size_t>(column)];
        }
        if (factor.cholesky) {
            rightHandSide[static_cast<std::size_t>(row)] /=
                factor.lower[static_cast<std::size_t>(row * size + row)];
        }
    }
}

void solveDenseSymmetricMultiple(
    const DenseSymmetricFactor& factor,
    std::vector<double>& rightHandSides,
    int rightHandSideCount)
{
    const int size = factor.size;
    if (rightHandSideCount <= 0
        || rightHandSides.size() != static_cast<std::size_t>(
            size * rightHandSideCount)) {
        throw std::runtime_error(
            "[Local ROM] Dense multi-RHS dimensions are invalid.");
    }
#ifdef USE_MKL_PARDISO
    // The block is row-major size-by-rightHandSideCount. Two BLAS-3
    // triangular solves are substantially faster than scalar RHS loops.
    const CBLAS_DIAG diagonal = factor.cholesky
        ? CblasNonUnit : CblasUnit;
    cblas_dtrsm(
        CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans,
        diagonal, size, rightHandSideCount, 1.0,
        factor.lower.data(), size,
        rightHandSides.data(), rightHandSideCount);
    if (!factor.cholesky) {
        for (int row = 0; row < size; ++row) {
            const double inverse =
                1.0 / factor.diagonal[static_cast<std::size_t>(row)];
            double* values = rightHandSides.data()
                + static_cast<std::size_t>(row * rightHandSideCount);
            for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
                values[rhs] *= inverse;
            }
        }
    }
    cblas_dtrsm(
        CblasRowMajor, CblasLeft, CblasLower, CblasTrans,
        diagonal, size, rightHandSideCount, 1.0,
        factor.lower.data(), size,
        rightHandSides.data(), rightHandSideCount);
#else
    for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
        std::vector<double> column(static_cast<std::size_t>(size), 0.0);
        for (int row = 0; row < size; ++row) {
            column[static_cast<std::size_t>(row)] = rightHandSides[
                static_cast<std::size_t>(row * rightHandSideCount + rhs)];
        }
        solveDenseSymmetric(factor, column);
        for (int row = 0; row < size; ++row) {
            rightHandSides[static_cast<std::size_t>(
                row * rightHandSideCount + rhs)] =
                column[static_cast<std::size_t>(row)];
        }
    }
#endif
}

LocalReducedSchurSolver::LocalReducedSchurSolver(
    const Model& model,
    const Mesh& mesh,
    const CaseConfig& physics,
    const ddm_schur::InterfacePartition& partition,
    const ddm_schur::Options& options,
    int matrixFreeInterfaceThreshold,
    const std::filesystem::path& outputDirectory)
    : model_(model)
{
    (void)mesh;
    (void)physics;
    (void)partition;
    (void)matrixFreeInterfaceThreshold;
    if (options.interfaceKrylov != "augmented-direct") {
        throw std::runtime_error(
            "[Production Schur] Only augmented-direct is available.");
    }
    initializeAugmentedDirect(options, outputDirectory);
}

LocalReducedSchurSolver::~LocalReducedSchurSolver() = default;

void LocalReducedSchurSolver::initializeAugmentedDirect(
    const ddm_schur::Options& options,
    const std::filesystem::path& outputDirectory)
{
    augmentedDirect_ = std::make_unique<AugmentedDirectData>();
    augmentedDirect_->dimension = model_.interfaceDofs
        + model_.totalLocalRank;
    augmentedDirect_->solveThreads = std::max(
        1, options.localPardisoThreads);
    augmentedDirect_->workThreads = std::max(
        1, options.localSolveThreads);

    const auto assemblyStart = Clock::now();
    UpperCsr matrix = buildAugmentedUpperCsr(
        model_, augmentedDirect_->rankOffsets,
        augmentedDirect_->couplingSymmetryRelativeError,
        std::max(1, options.localSolveThreads));
    interfaceNonzeros_ = matrix.values.size();
    assemblySeconds_ = std::chrono::duration<double>(
        Clock::now() - assemblyStart).count();

    const auto factorStart = Clock::now();
    {
        ScopedDirectSolverMklThreads mklThreads(
            augmentedDirect_->solveThreads);
        augmentedDirect_->factor = std::make_unique<SubdomainDirectSolver>(
            augmentedDirect_->dimension,
            matrix.rowPtr, matrix.colInd, matrix.values);
    }
    factorizationSeconds_ = std::chrono::duration<double>(
        Clock::now() - factorStart).count();

    // Keep a compact machine-readable factor report. It is small enough for
    // production runs and contains no temperature field data.
    std::filesystem::create_directories(outputDirectory);
    std::ofstream summary(outputDirectory / "augmented_direct_summary.csv");
    summary
        << "dimension,interface_dofs,reduced_dofs,nonzeros,"
           "coupling_symmetry_relative_error,assembly_seconds,"
           "factor_seconds,symbolic_seconds,numerical_seconds,"
           "factor_memory_bytes,factor_threads\n"
        << std::setprecision(17)
        << augmentedDirect_->dimension << ',' << model_.interfaceDofs << ','
        << model_.totalLocalRank << ',' << interfaceNonzeros_ << ','
        << augmentedDirect_->couplingSymmetryRelativeError << ','
        << assemblySeconds_ << ',' << factorizationSeconds_ << ','
        << augmentedDirect_->factor->symbolicAnalysisSeconds() << ','
        << augmentedDirect_->factor->numericalFactorizationSeconds() << ','
        << augmentedDirect_->factor->memoryBytes() << ','
        << augmentedDirect_->solveThreads << '\n';

    std::cout << "[Augmented direct] dimension="
              << augmentedDirect_->dimension
              << ", nnz=" << interfaceNonzeros_
              << ", assembly=" << assemblySeconds_
              << " s, factor=" << factorizationSeconds_ << " s\n";
}

double LocalReducedSchurSolver::symbolicAnalysisSeconds() const
{
    return augmentedDirect_ && augmentedDirect_->factor
        ? augmentedDirect_->factor->symbolicAnalysisSeconds() : 0.0;
}

double LocalReducedSchurSolver::numericalFactorizationSeconds() const
{
    return augmentedDirect_ && augmentedDirect_->factor
        ? augmentedDirect_->factor->numericalFactorizationSeconds() : 0.0;
}

int LocalReducedSchurSolver::symbolicAnalysisCalls() const
{
    return augmentedDirect_ && augmentedDirect_->factor
        ? augmentedDirect_->factor->symbolicAnalysisCalls() : 0;
}

int LocalReducedSchurSolver::numericalFactorizationCalls() const
{
    return augmentedDirect_ && augmentedDirect_->factor
        ? augmentedDirect_->factor->numericalFactorizationCalls() : 0;
}

std::size_t LocalReducedSchurSolver::factorMemoryBytes() const
{
    return augmentedDirect_ && augmentedDirect_->factor
        ? augmentedDirect_->factor->memoryBytes() : 0;
}

SolveResult LocalReducedSchurSolver::solve(
    const std::vector<double>& globalRhs,
    const std::vector<double>* interfaceInitialGuess,
    double interfaceToleranceOverride)
{
    (void)interfaceInitialGuess;
    (void)interfaceToleranceOverride;
    if (globalRhs.size() != static_cast<std::size_t>(model_.globalDofs)) {
        throw std::runtime_error("[Local ROM] Global RHS has the wrong size.");
    }

    std::vector<double> interfaceRhs(
        static_cast<std::size_t>(model_.interfaceDofs), 0.0);
    for (int gamma = 0; gamma < model_.interfaceDofs; ++gamma) {
        interfaceRhs[static_cast<std::size_t>(gamma)] = globalRhs[
            static_cast<std::size_t>(model_.interfaceGlobalDofs[
                static_cast<std::size_t>(gamma)])];
    }

    const auto projectionStart = Clock::now();
    std::vector<std::vector<double>> projectedInteriorRhs(
        model_.subdomains.size());
    for (std::size_t slot = 0; slot < model_.subdomains.size(); ++slot) {
        const SubdomainModel& local = model_.subdomains[slot];
        const SubdomainModel& data = templatePayload(model_, local);
        std::vector<double>& reduced = projectedInteriorRhs[slot];
        reduced.assign(static_cast<std::size_t>(local.rank), 0.0);
        for (int mode = 0; mode < local.rank; ++mode) {
            for (int row = 0; row < local.interiorDofs; ++row) {
                reduced[static_cast<std::size_t>(mode)] +=
                    data.basis[static_cast<std::size_t>(
                        mode * local.interiorDofs + row)]
                    * globalRhs[static_cast<std::size_t>(
                        local.interiorGlobalDofs[
                            static_cast<std::size_t>(row)])];
            }
        }
    }
    const double projectionSeconds = std::chrono::duration<double>(
        Clock::now() - projectionStart).count();

    SolveResult result = solveWithReducedRhs(
        std::move(projectedInteriorRhs), std::move(interfaceRhs));
    result.timing.localReducedAssemblySeconds += projectionSeconds;
    return result;
}

SolveResult LocalReducedSchurSolver::solveReducedRhs(
    const std::vector<std::vector<double>>& projectedInteriorRhs,
    const std::vector<double>& interfaceRhs,
    const std::vector<double>* interfaceInitialGuess,
    double interfaceToleranceOverride)
{
    (void)interfaceInitialGuess;
    (void)interfaceToleranceOverride;
    return solveWithReducedRhs(projectedInteriorRhs, interfaceRhs);
}

SolveResult LocalReducedSchurSolver::solveWithReducedRhs(
    std::vector<std::vector<double>> projectedInteriorRhs,
    std::vector<double> interfaceRhs)
{
    if (!augmentedDirect_ || !augmentedDirect_->factor) {
        throw std::runtime_error(
            "[Production Schur] Augmented factor is not initialized.");
    }
    if (interfaceRhs.size() != static_cast<std::size_t>(model_.interfaceDofs)
        || projectedInteriorRhs.size() != model_.subdomains.size()) {
        throw std::runtime_error("[Local ROM] Reduced RHS dimensions are invalid.");
    }
    for (std::size_t slot = 0; slot < model_.subdomains.size(); ++slot) {
        if (projectedInteriorRhs[slot].size()
            != static_cast<std::size_t>(model_.subdomains[slot].rank)) {
            throw std::runtime_error(
                "[Local ROM] Reduced interior RHS rank is invalid.");
        }
    }

    const auto totalStart = Clock::now();
    const auto assemblyStart = Clock::now();

    // Reference-state terms are stored separately from the affine reduced
    // coordinates. Subtract A*T_ref before placing each block in the coupled
    // RHS so the solved q_i remains a correction around that reference.
    std::vector<double> condensed = std::move(interfaceRhs);
    for (std::size_t slot = 0; slot < model_.subdomains.size(); ++slot) {
        const SubdomainModel& local = model_.subdomains[slot];
        const SubdomainModel& data = templatePayload(model_, local);
        std::vector<double>& reduced = projectedInteriorRhs[slot];
        const bool hasInteriorReferenceImage = std::any_of(
            data.interiorReferenceImage.begin(),
            data.interiorReferenceImage.end(),
            [](double value) { return value != 0.0; });
        if (hasInteriorReferenceImage) {
            for (int mode = 0; mode < local.rank; ++mode) {
                for (int row = 0; row < local.interiorDofs; ++row) {
                    reduced[static_cast<std::size_t>(mode)] -=
                        data.basis[static_cast<std::size_t>(
                            mode * local.interiorDofs + row)]
                        * data.interiorReferenceImage[
                            static_cast<std::size_t>(row)];
                }
            }
        }
        for (std::size_t localGamma = 0;
             localGamma < local.interfaceIndices.size(); ++localGamma) {
            condensed[static_cast<std::size_t>(
                local.interfaceIndices[localGamma])] -=
                data.interfaceReferenceImage[localGamma];
        }
    }

    std::vector<double> augmentedRhs(
        static_cast<std::size_t>(augmentedDirect_->dimension), 0.0);
    std::copy(condensed.begin(), condensed.end(), augmentedRhs.begin());
    for (std::size_t slot = 0; slot < model_.subdomains.size(); ++slot) {
        const int offset = model_.interfaceDofs
            + augmentedDirect_->rankOffsets[slot];
        std::copy(projectedInteriorRhs[slot].begin(),
                  projectedInteriorRhs[slot].end(),
                  augmentedRhs.begin() + offset);
    }

    SolveResult result;
    result.timing.localReducedAssemblySeconds =
        std::chrono::duration<double>(Clock::now() - assemblyStart).count();

    const auto solveStart = Clock::now();
    std::vector<double> augmentedSolution;
    {
        ScopedDirectSolverMklThreads mklThreads(
            augmentedDirect_->solveThreads);
        augmentedDirect_->factor->solve(
            augmentedRhs, augmentedSolution);
    }
    result.timing.interfaceSolveSeconds =
        std::chrono::duration<double>(Clock::now() - solveStart).count();

    result.interfaceTemperature.assign(
        augmentedSolution.begin(),
        augmentedSolution.begin() + model_.interfaceDofs);
    result.localReducedCoordinates.resize(model_.subdomains.size());
    for (std::size_t slot = 0; slot < model_.subdomains.size(); ++slot) {
        const int begin = model_.interfaceDofs
            + augmentedDirect_->rankOffsets[slot];
        const int end = model_.interfaceDofs
            + augmentedDirect_->rankOffsets[slot + 1];
        result.localReducedCoordinates[slot].assign(
            augmentedSolution.begin() + begin,
            augmentedSolution.begin() + end);
    }

    // A residual audit every 25 solves catches corrupted factors or assembly
    // drift while avoiding an additional sparse product on every time step.
    constexpr int residualAuditInterval = 25;
    if (augmentedDirect_->solveCalls % residualAuditInterval == 0) {
        std::vector<double> interfaceResidual = condensed;
        std::vector<std::vector<double>> interiorResidual =
            projectedInteriorRhs;
        for (const InterfaceEntry& entry : model_.interfaceEntries) {
            interfaceResidual[static_cast<std::size_t>(entry.row)] -=
                entry.value * result.interfaceTemperature[
                    static_cast<std::size_t>(entry.column)];
        }
        for (std::size_t slot = 0; slot < model_.subdomains.size(); ++slot) {
            const SubdomainModel& local = model_.subdomains[slot];
            const SubdomainModel& data = templatePayload(model_, local);
            const std::vector<double>& coordinates =
                result.localReducedCoordinates[slot];
            for (std::size_t localGamma = 0;
                 localGamma < local.interfaceIndices.size(); ++localGamma) {
                const int gamma = local.interfaceIndices[localGamma];
                for (int mode = 0; mode < local.rank; ++mode) {
                    interfaceResidual[static_cast<std::size_t>(gamma)] -=
                        data.reducedInterfaceInterior[
                            localGamma * static_cast<std::size_t>(local.rank)
                            + static_cast<std::size_t>(mode)]
                        * coordinates[static_cast<std::size_t>(mode)];
                    interiorResidual[slot][static_cast<std::size_t>(mode)] -=
                        data.reducedInteriorInterface[
                            static_cast<std::size_t>(mode)
                                * local.localInterfaceDofs
                            + localGamma]
                        * result.interfaceTemperature[
                            static_cast<std::size_t>(gamma)];
                }
            }
            for (int row = 0; row < local.rank; ++row) {
                for (int column = 0; column < local.rank; ++column) {
                    interiorResidual[slot][static_cast<std::size_t>(row)] -=
                        data.reducedInterior[static_cast<std::size_t>(
                            row * local.rank + column)]
                        * coordinates[static_cast<std::size_t>(column)];
                }
            }
        }

        double residualSquared = vectorDot(
            interfaceResidual, interfaceResidual);
        double rhsSquared = vectorDot(condensed, condensed);
        for (std::size_t slot = 0;
             slot < interiorResidual.size(); ++slot) {
            residualSquared += vectorDot(
                interiorResidual[slot], interiorResidual[slot]);
            rhsSquared += vectorDot(
                projectedInteriorRhs[slot], projectedInteriorRhs[slot]);
        }
        augmentedDirect_->validatedRelativeResidual =
            std::sqrt(std::max(0.0, residualSquared))
            / std::max(std::sqrt(std::max(0.0, rhsSquared)),
                       std::numeric_limits<double>::min());
    }
    for (const double value : augmentedSolution) {
        if (!std::isfinite(value)) {
            augmentedDirect_->validatedRelativeResidual =
                std::numeric_limits<double>::infinity();
            break;
        }
    }
    ++augmentedDirect_->solveCalls;

    result.timing.interfaceIterations = 1;
    result.timing.interfaceMatvecs = 0;
    result.timing.interfaceInitialRelativeResidual = 1.0;
    result.timing.interfaceRelativeResidual =
        augmentedDirect_->validatedRelativeResidual;
    result.timing.interfaceKrylovActual = "augmented-direct";
    result.status = std::isfinite(
            augmentedDirect_->validatedRelativeResidual)
            && augmentedDirect_->validatedRelativeResidual <= 1.0e-8
        ? "success" : "augmented_direct_residual_failed";

    // Reconstructing the full vector is currently needed by residual and
    // summary calculations, but the production entry point never writes it.
    const auto reconstructionStart = Clock::now();
    result.temperature.assign(static_cast<std::size_t>(model_.globalDofs), 0.0);
    for (int gamma = 0; gamma < model_.interfaceDofs; ++gamma) {
        result.temperature[static_cast<std::size_t>(
            model_.interfaceGlobalDofs[static_cast<std::size_t>(gamma)])] =
            result.interfaceTemperature[static_cast<std::size_t>(gamma)];
    }
    const int domainCount = static_cast<int>(model_.subdomains.size());
#pragma omp parallel for num_threads(augmentedDirect_->workThreads) \
    if(augmentedDirect_->workThreads > 1) schedule(static)
    for (int slot = 0; slot < domainCount; ++slot) {
        const SubdomainModel& local =
            model_.subdomains[static_cast<std::size_t>(slot)];
        const SubdomainModel& data = templatePayload(model_, local);
        const std::vector<double>& coordinates =
            result.localReducedCoordinates[static_cast<std::size_t>(slot)];
        for (int row = 0; row < local.interiorDofs; ++row) {
            double value = data.referenceInterior[static_cast<std::size_t>(row)];
            for (int mode = 0; mode < local.rank; ++mode) {
                value += data.basis[static_cast<std::size_t>(
                    mode * local.interiorDofs + row)]
                    * coordinates[static_cast<std::size_t>(mode)];
            }
            result.temperature[static_cast<std::size_t>(
                local.interiorGlobalDofs[static_cast<std::size_t>(row)])] =
                value;
        }
    }
    result.timing.fullFieldReconstructionSeconds =
        std::chrono::duration<double>(
            Clock::now() - reconstructionStart).count();
    result.timing.totalSeconds =
        std::chrono::duration<double>(Clock::now() - totalStart).count();
    return result;
}

} // namespace mor::local
