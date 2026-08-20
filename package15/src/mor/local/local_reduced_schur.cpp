#include "local_reduced_schur.hpp"

#include "ddm_schur/interface_operator.hpp"
#include "ddm_schur/schur_fgmres.hpp"
#include "ddm_schur/schur_proxy.hpp"
#include "port_core_solver.hpp"
#include "sipg_core.hpp"
#include "linear_solvers.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>

namespace mor::local {

struct LocalReducedSchurSolver::SparseFactor {
    explicit SparseFactor(int size, const std::vector<MatrixEntry>& entries)
        : solver(size, entries)
    {
    }

    SubdomainDirectSolver solver;
};

struct LocalReducedSchurSolver::AugmentedDirectData {
    std::unique_ptr<SubdomainDirectSolver> factor;
    std::unique_ptr<GeneralSparseDirectSolver> generalFactor;
    std::vector<int> rankOffsets;
    int dimension = 0;
    int solveThreads = 1;
    int solveCalls = 0;
    double couplingSymmetryRelativeError = 0.0;
    double validatedRelativeResidual =
        std::numeric_limits<double>::infinity();
};

std::vector<MatrixEntry> buildAugmentedGeneralEntries(
    const Model& model, std::vector<int>& rankOffsets,
    double& couplingSymmetryRelativeError)
{
    rankOffsets.assign(model.subdomains.size() + 1, 0);
    for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
        rankOffsets[slot + 1] = rankOffsets[slot] + model.subdomains[slot].rank;
    }
    if (rankOffsets.back() != model.totalLocalRank) {
        throw std::runtime_error("[Augmented direct] Local rank offsets are inconsistent.");
    }
    std::vector<MatrixEntry> entries;
    entries.reserve(model.interfaceEntries.size() +
        static_cast<std::size_t>(model.totalLocalRank) * 32);
    for (const InterfaceEntry& entry : model.interfaceEntries) {
        entries.push_back({entry.row, entry.column, entry.value});
    }
    double maximumCoupling = 0.0;
    double maximumDifference = 0.0;
    for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
        const SubdomainModel& local = model.subdomains[slot];
        // Petrov projection is rebuilt per local model, so it must not borrow
        // a Galerkin template payload.
        const SubdomainModel& data = local;
        const int offset = model.interfaceDofs + rankOffsets[slot];
        for (int row = 0; row < local.rank; ++row) {
            for (int column = 0; column < local.rank; ++column) {
                const double value = data.reducedInterior[static_cast<std::size_t>(
                    row * local.rank + column)];
                if (value != 0.0) entries.push_back({offset + row, offset + column, value});
            }
        }
        for (std::size_t gammaLocal = 0; gammaLocal < local.interfaceIndices.size(); ++gammaLocal) {
            const int gamma = local.interfaceIndices[gammaLocal];
            for (int mode = 0; mode < local.rank; ++mode) {
                const double gammaInterior = data.reducedInterfaceInterior[
                    gammaLocal * static_cast<std::size_t>(local.rank) + mode];
                const double interiorGamma = data.reducedInteriorInterface[
                    static_cast<std::size_t>(mode) * local.localInterfaceDofs + gammaLocal];
                if (gammaInterior != 0.0) entries.push_back({gamma, offset + mode, gammaInterior});
                if (interiorGamma != 0.0) entries.push_back({offset + mode, gamma, interiorGamma});
                maximumCoupling = std::max(maximumCoupling,
                    std::max(std::abs(gammaInterior), std::abs(interiorGamma)));
                maximumDifference = std::max(maximumDifference,
                    std::abs(gammaInterior - interiorGamma));
            }
        }
    }
    couplingSymmetryRelativeError = maximumCoupling > 0.0
        ? maximumDifference / maximumCoupling : 0.0;
    return entries;
}

struct LocalReducedSchurSolver::MatrixFreeData {
    struct CoarseVector {
        bool dense = false;
        std::vector<int> indices;
        std::vector<double> values;

        double dot(const std::vector<double>& vector) const
        {
            double result = 0.0;
            if (dense) {
                if (values.size() != vector.size()) {
                    throw std::runtime_error(
                        "[Local ROM] Dense coarse vector size mismatch.");
                }
                for (std::size_t i = 0; i < values.size(); ++i) {
                    result += values[i] * vector[i];
                }
            } else {
                for (std::size_t i = 0; i < indices.size(); ++i) {
                    result += values[i]
                        * vector[static_cast<std::size_t>(indices[i])];
                }
            }
            return result;
        }

        void addScaled(double coefficient, std::vector<double>& vector) const
        {
            if (dense) {
                if (values.size() != vector.size()) {
                    throw std::runtime_error(
                        "[Local ROM] Dense coarse vector size mismatch.");
                }
                for (std::size_t i = 0; i < values.size(); ++i) {
                    vector[i] += coefficient * values[i];
                }
            } else {
                for (std::size_t i = 0; i < indices.size(); ++i) {
                    vector[static_cast<std::size_t>(indices[i])] +=
                        coefficient * values[i];
                }
            }
        }

        std::vector<double> expand(std::size_t size) const
        {
            if (dense) {
                if (values.size() != size) {
                    throw std::runtime_error(
                        "[Local ROM] Dense coarse vector size mismatch.");
                }
                return values;
            }
            std::vector<double> vector(size, 0.0);
            for (std::size_t i = 0; i < indices.size(); ++i) {
                vector[static_cast<std::size_t>(indices[i])] = values[i];
            }
            return vector;
        }
    };
    struct LocalApplyWorkspace {
        std::vector<double> interfaceInput;
        std::vector<double> reducedProduct;
        std::vector<double> correction;
        std::vector<double> recoveryCoordinates;
    };

    ddm_schur::Options options;
    std::vector<int> rowPtr;
    std::vector<int> colInd;
    std::vector<double> values;
    std::unique_ptr<ddm_schur::SchurProxyPreconditioner> proxy;
    std::vector<CoarseVector> coarseBasis;
    std::vector<CoarseVector> coarseImages;
    DenseSymmetricFactor coarseFactor;
    DenseSymmetricFactor predictorFactor;
    int geometricCoarseDimension = 0;
    int operatorCoarseDimension = 0;
    bool operatorCoarseCacheHit = false;
    double operatorCoarseSetupSeconds = 0.0;
    double operatorCoarseCacheLoadSeconds = 0.0;
    double operatorCoarseCacheSaveSeconds = 0.0;
    std::vector<LocalApplyWorkspace> localApplyWorkspaces;
    std::vector<double> preconditionerCoarseRhs;
    std::vector<double> preconditionerLeftRhs;
    std::vector<double> preconditionerCoarseCorrection;
    std::vector<double> preconditionerCoarseImage;
    std::vector<double> preconditionerProjectedResidual;
    std::vector<double> preconditionerLocalCorrection;
    ddm_schur::MatrixFreeFgmresWorkspace fgmresWorkspace;
    std::vector<double> previousInterfaceSolution;
    std::vector<double> predictorSolution;
    std::vector<double> predictorProduct;
    std::vector<double> predictorResidual;
    std::vector<double> predictorCoarseRhs;
    int matvecs = 0;
    double coarseSolveSeconds = 0.0;
};

namespace {

constexpr std::uint64_t operatorCoarseCacheMagic =
    UINT64_C(0x4f50434f41525345); // "OPCOARSE"
constexpr int operatorCoarseCacheVersion = 3;

template <typename T>
void appendHash(std::uint64_t& hash, const T& value)
{
    static_assert(std::is_trivially_copyable_v<T>);
    const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
    for (std::size_t i = 0; i < sizeof(T); ++i) {
        hash ^= bytes[i];
        hash *= UINT64_C(1099511628211);
    }
}

template <typename T>
void appendHash(std::uint64_t& hash, const std::vector<T>& values)
{
    static_assert(std::is_trivially_copyable_v<T>);
    appendHash(hash, values.size());
    const auto* bytes = reinterpret_cast<const unsigned char*>(values.data());
    const std::size_t byteCount = values.size() * sizeof(T);
    for (std::size_t i = 0; i < byteCount; ++i) {
        hash ^= bytes[i];
        hash *= UINT64_C(1099511628211);
    }
}

template <typename T>
void writeBinary(std::ofstream& output, const T& value)
{
    static_assert(std::is_trivially_copyable_v<T>);
    output.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template <typename T>
void writeBinaryVector(std::ofstream& output, const std::vector<T>& values)
{
    static_assert(std::is_trivially_copyable_v<T>);
    writeBinary(output, static_cast<std::uint64_t>(values.size()));
    if (!values.empty()) {
        output.write(reinterpret_cast<const char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(T)));
    }
}

template <typename T>
bool readBinary(std::ifstream& input, T& value)
{
    static_assert(std::is_trivially_copyable_v<T>);
    input.read(reinterpret_cast<char*>(&value), sizeof(T));
    return static_cast<bool>(input);
}

template <typename T>
bool readBinaryVector(std::ifstream& input,
                      std::vector<T>& values,
                      std::uint64_t maximumCount)
{
    static_assert(std::is_trivially_copyable_v<T>);
    std::uint64_t count = 0;
    if (!readBinary(input, count) || count > maximumCount) {
        return false;
    }
    values.resize(static_cast<std::size_t>(count));
    if (!values.empty()) {
        input.read(reinterpret_cast<char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(T)));
    }
    return static_cast<bool>(input);
}

const SubdomainModel& templatePayload(const Model& model,
                                      const SubdomainModel& instance)
{
    if (!instance.templateReused) {
        return instance;
    }
    for (const SubdomainModel& candidate : model.subdomains) {
        if (candidate.templateId == instance.templateId && !candidate.templateReused) {
            return candidate;
        }
    }
    throw std::runtime_error("[Local ROM] Reused template payload is missing.");
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

UpperCsr buildAugmentedUpperCsr(const Model& model,
                                std::vector<int>& rankOffsets,
                                double& couplingSymmetryRelativeError,
                                int assemblyThreads)
{
    // Unknown ordering is [full physical interface, local ROM coordinates].
    // Keeping the interface unreduced avoids a data-trained port basis while
    // PARDISO performs the equivalent sparse elimination internally.
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
        throw std::runtime_error(
            "[Augmented direct] Empty reduced system.");
    }

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
            "[Augmented direct] CSR exceeds the 32-bit PARDISO index range.");
    }
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
                        "[Augmented direct] Invalid assembled CSR entry.");
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
    result.colInd.resize(
        static_cast<std::size_t>(result.rowPtr.back()));
    result.values.resize(result.colInd.size());
#pragma omp parallel for num_threads(assemblyThreads) if(assemblyThreads > 1) schedule(static)
    for (int row = 0; row < dimension; ++row) {
        const std::size_t source = rawRowPtr[static_cast<std::size_t>(row)];
        const std::size_t target = static_cast<std::size_t>(
            result.rowPtr[static_cast<std::size_t>(row)]);
        const std::size_t count = compactCounts[static_cast<std::size_t>(row)];
        for (std::size_t entry = 0; entry < count; ++entry) {
            result.colInd[target + entry] =
                rawEntries[source + entry].column;
            result.values[target + entry] =
                rawEntries[source + entry].value;
        }
    }
    return result;
}

double vectorDot(const std::vector<double>& left,
                 const std::vector<double>& right)
{
    double value = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        value += left[i] * right[i];
    }
    return value;
}

} // namespace

DenseSymmetricFactor factorDenseSymmetric(const std::vector<double>& matrix, int size)
{
    if (size <= 0 || matrix.size() != static_cast<std::size_t>(size * size)) {
        throw std::runtime_error("[Local ROM] Invalid dense symmetric factor dimensions.");
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
    for (int row = 0; row < size; ++row) {
        for (int column = 0; column <= row; ++column) {
            double value = 0.5 * (matrix[static_cast<std::size_t>(row * size + column)]
                + matrix[static_cast<std::size_t>(column * size + row)]);
            for (int k = 0; k < column; ++k) {
                value -= factor.lower[static_cast<std::size_t>(row * size + k)]
                    * factor.lower[static_cast<std::size_t>(column * size + k)];
            }
            if (row == column) {
                if (!(value > threshold) || !std::isfinite(value)) {
                    factor.cholesky = false;
                    break;
                }
                factor.lower[static_cast<std::size_t>(row * size + column)] = std::sqrt(value);
            } else {
                factor.lower[static_cast<std::size_t>(row * size + column)] = value
                    / factor.lower[static_cast<std::size_t>(column * size + column)];
            }
        }
        if (!factor.cholesky) {
            break;
        }
    }
    if (factor.cholesky) {
        return factor;
    }

    factor.lower.assign(matrix.size(), 0.0);
    factor.diagonal.assign(static_cast<std::size_t>(size), 0.0);
    for (int row = 0; row < size; ++row) {
        factor.lower[static_cast<std::size_t>(row * size + row)] = 1.0;
        double pivot = matrix[static_cast<std::size_t>(row * size + row)];
        for (int k = 0; k < row; ++k) {
            const double entry = factor.lower[static_cast<std::size_t>(row * size + k)];
            pivot -= entry * entry * factor.diagonal[static_cast<std::size_t>(k)];
        }
        if (!(std::abs(pivot) > threshold) || !std::isfinite(pivot)) {
            throw std::runtime_error("[Local ROM] Reduced matrix factorization encountered a tiny pivot.");
        }
        factor.diagonal[static_cast<std::size_t>(row)] = pivot;
        for (int next = row + 1; next < size; ++next) {
            double value = 0.5 * (matrix[static_cast<std::size_t>(next * size + row)]
                + matrix[static_cast<std::size_t>(row * size + next)]);
            for (int k = 0; k < row; ++k) {
                value -= factor.lower[static_cast<std::size_t>(next * size + k)]
                    * factor.lower[static_cast<std::size_t>(row * size + k)]
                    * factor.diagonal[static_cast<std::size_t>(k)];
            }
            factor.lower[static_cast<std::size_t>(next * size + row)] = value / pivot;
        }
    }
    return factor;
}

void solveDenseSymmetric(const DenseSymmetricFactor& factor,
                         std::vector<double>& rightHandSide)
{
    const int size = factor.size;
    if (rightHandSide.size() != static_cast<std::size_t>(size)) {
        throw std::runtime_error("[Local ROM] Dense solve RHS has the wrong size.");
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
            "[Local ROM] Dense multi-RHS solve dimensions are invalid.");
    }
#ifdef USE_MKL_PARDISO
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
    for (int row = 0; row < size; ++row) {
        double* values = rightHandSides.data()
            + static_cast<std::size_t>(row * rightHandSideCount);
        for (int column = 0; column < row; ++column) {
            const double coefficient = factor.lower[
                static_cast<std::size_t>(row * size + column)];
            const double* prior = rightHandSides.data()
                + static_cast<std::size_t>(
                    column * rightHandSideCount);
            for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
                values[rhs] -= coefficient * prior[rhs];
            }
        }
        if (factor.cholesky) {
            const double inverse = 1.0 / factor.lower[
                static_cast<std::size_t>(row * size + row)];
            for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
                values[rhs] *= inverse;
            }
        }
    }
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
    for (int row = size - 1; row >= 0; --row) {
        double* values = rightHandSides.data()
            + static_cast<std::size_t>(row * rightHandSideCount);
        for (int column = row + 1; column < size; ++column) {
            const double coefficient = factor.lower[
                static_cast<std::size_t>(column * size + row)];
            const double* next = rightHandSides.data()
                + static_cast<std::size_t>(
                    column * rightHandSideCount);
            for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
                values[rhs] -= coefficient * next[rhs];
            }
        }
        if (factor.cholesky) {
            const double inverse = 1.0 / factor.lower[
                static_cast<std::size_t>(row * size + row)];
            for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
                values[rhs] *= inverse;
            }
        }
    }
#endif
}

LocalReducedSchurSolver::LocalReducedSchurSolver(const Model& model)
    : model_(model)
{
    localFactors_.reserve(model.subdomains.size());
    for (const SubdomainModel& local : model.subdomains) {
        const SubdomainModel& data = templatePayload(model, local);
        localFactors_.push_back(factorDenseSymmetric(data.reducedInterior, local.rank));
    }

    initializeExplicit();
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
    const auto localFactorStart = std::chrono::steady_clock::now();
    const bool coupledDirect = options.interfaceKrylov == "port-core"
        || options.interfaceKrylov == "augmented-direct";
    if (!coupledDirect) {
        localFactors_.reserve(model.subdomains.size());
        for (const SubdomainModel& local : model.subdomains) {
            const SubdomainModel& data = templatePayload(model, local);
            localFactors_.push_back(
                factorDenseSymmetric(data.reducedInterior, local.rank));
        }
    }
    factorizationSeconds_ = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - localFactorStart).count();

    if (matrixFreeInterfaceThreshold < 0) {
        throw std::runtime_error(
            "[Local ROM] Matrix-free interface threshold must be nonnegative.");
    }
    if (options.interfaceKrylov == "port-core") {
        initializePortCore(partition, options, outputDirectory);
    } else if (options.interfaceKrylov == "augmented-direct") {
        initializeAugmentedDirect(options, outputDirectory);
    } else if (model.interfaceDofs > matrixFreeInterfaceThreshold) {
        initializeMatrixFree(mesh, physics, partition, options, outputDirectory);
    } else {
        initializeExplicit();
    }
}

void LocalReducedSchurSolver::initializeAugmentedDirect(
    const ddm_schur::Options& options,
    const std::filesystem::path& outputDirectory)
{
    augmentedDirect_ = std::make_unique<AugmentedDirectData>();
    augmentedDirect_->dimension = model_.interfaceDofs
        + model_.totalLocalRank;
    augmentedDirect_->solveThreads = std::max(
        1, options.localPardisoThreads);

    const auto assemblyStart = std::chrono::steady_clock::now();
    UpperCsr matrix;
    std::vector<MatrixEntry> generalEntries;
    if (model_.usesPetrovTestSpace) {
        generalEntries = buildAugmentedGeneralEntries(
            model_, augmentedDirect_->rankOffsets,
            augmentedDirect_->couplingSymmetryRelativeError);
        interfaceNonzeros_ = generalEntries.size();
    } else {
        matrix = buildAugmentedUpperCsr(
            model_, augmentedDirect_->rankOffsets,
            augmentedDirect_->couplingSymmetryRelativeError,
            std::max(1, options.localSolveThreads));
        interfaceNonzeros_ = matrix.values.size();
    }
    assemblySeconds_ = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - assemblyStart).count();

    const auto factorStart = std::chrono::steady_clock::now();
    {
        ScopedDirectSolverMklThreads mklThreads(
            augmentedDirect_->solveThreads);
        if (model_.usesPetrovTestSpace) {
            augmentedDirect_->generalFactor =
                std::make_unique<GeneralSparseDirectSolver>(
                    augmentedDirect_->dimension, generalEntries);
        } else {
            augmentedDirect_->factor =
                std::make_unique<SubdomainDirectSolver>(
                    augmentedDirect_->dimension,
                    matrix.rowPtr, matrix.colInd, matrix.values);
        }
    }
    const double factorSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - factorStart).count();
    factorizationSeconds_ += factorSeconds;
    interfaceSolver_ = "augmented-pardiso";

    std::filesystem::create_directories(outputDirectory);
    std::ofstream summary(
        outputDirectory / "augmented_direct_summary.csv");
    summary
        << "dimension,interface_dofs,reduced_dofs,nonzeros,"
           "coupling_symmetry_relative_error,assembly_seconds,"
           "factor_seconds,symbolic_seconds,numerical_seconds,"
           "factor_memory_bytes,factor_threads\n"
        << std::setprecision(17)
        << augmentedDirect_->dimension << ',' << model_.interfaceDofs << ','
        << model_.totalLocalRank << ',' << interfaceNonzeros_ << ','
        << augmentedDirect_->couplingSymmetryRelativeError << ','
        << assemblySeconds_ << ',' << factorSeconds << ','
        << (model_.usesPetrovTestSpace
                ? augmentedDirect_->generalFactor->symbolicAnalysisSeconds()
                : augmentedDirect_->factor->symbolicAnalysisSeconds()) << ','
        << (model_.usesPetrovTestSpace
                ? augmentedDirect_->generalFactor->numericalFactorizationSeconds()
                : augmentedDirect_->factor->numericalFactorizationSeconds()) << ','
        << (model_.usesPetrovTestSpace
                ? augmentedDirect_->generalFactor->memoryBytes()
                : augmentedDirect_->factor->memoryBytes()) << ','
        << augmentedDirect_->solveThreads << '\n';
    std::cout << "[Augmented direct] dimension="
              << augmentedDirect_->dimension
              << ", nnz=" << interfaceNonzeros_
              << ", assembly=" << assemblySeconds_
              << " s, factor=" << factorSeconds << " s\n";
}

void LocalReducedSchurSolver::initializePortCore(
    const ddm_schur::InterfacePartition& partition,
    const ddm_schur::Options& options,
    const std::filesystem::path& outputDirectory)
{
    const auto setupStart = std::chrono::steady_clock::now();
    portCore_ = std::make_unique<PortCoreSolver>(
        model_, partition, options.localSolveThreads,
        options.localPardisoThreads, outputDirectory,
        options.portCoreCacheEnabled
            ? std::filesystem::path(options.portCoreCachePath)
            : std::filesystem::path{});
    interfaceSolver_ = "port-core-direct";
    interfaceNonzeros_ = portCore_->nonzeros();
    factorizationSeconds_ += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - setupStart).count();
}

void LocalReducedSchurSolver::initializeExplicit()
{
    constexpr int denseInterfaceLimit = 2000;
    const bool useSparse = model_.interfaceDofs > denseInterfaceLimit;
    const auto assemblyStart = std::chrono::steady_clock::now();
    std::vector<double> dense;
    std::vector<MatrixEntry> sparse;
    if (useSparse) {
        interfaceSolver_ = "sparse-pardiso";
        std::size_t reserve = model_.interfaceEntries.size();
        for (const SubdomainModel& local : model_.subdomains) {
            reserve += local.interfaceIndices.size()
                * (local.interfaceIndices.size() + 1) / 2;
        }
        sparse.reserve(reserve);
        for (const InterfaceEntry& entry : model_.interfaceEntries) {
            sparse.push_back({entry.row, entry.column, entry.value});
        }
    } else {
        interfaceSolver_ = "dense-llt";
        dense.assign(static_cast<std::size_t>(model_.interfaceDofs)
            * static_cast<std::size_t>(model_.interfaceDofs), 0.0);
        for (const InterfaceEntry& entry : model_.interfaceEntries) {
            dense[static_cast<std::size_t>(entry.row)
                * static_cast<std::size_t>(model_.interfaceDofs)
                + static_cast<std::size_t>(entry.column)] += entry.value;
        }
    }

    for (std::size_t slot = 0; slot < model_.subdomains.size(); ++slot) {
        const SubdomainModel& local = model_.subdomains[slot];
        const SubdomainModel& data = templatePayload(model_, local);
        for (std::size_t localColumn = 0;
             localColumn < local.interfaceIndices.size(); ++localColumn) {
            const int gammaColumn = local.interfaceIndices[localColumn];
            std::vector<double> eliminated(static_cast<std::size_t>(local.rank), 0.0);
            for (int mode = 0; mode < local.rank; ++mode) {
                eliminated[static_cast<std::size_t>(mode)] =
                    data.reducedInteriorInterface[static_cast<std::size_t>(
                        mode * local.localInterfaceDofs) + localColumn];
            }
            solveDenseSymmetric(localFactors_[slot], eliminated);
            for (std::size_t localRow = 0;
                 localRow < local.interfaceIndices.size(); ++localRow) {
                const int gammaRow = local.interfaceIndices[localRow];
                if (useSparse && gammaRow > gammaColumn) {
                    continue;
                }
                double correction = 0.0;
                for (int mode = 0; mode < local.rank; ++mode) {
                    correction += data.reducedInterfaceInterior[static_cast<std::size_t>(
                        localRow * static_cast<std::size_t>(local.rank) + mode)]
                        * eliminated[static_cast<std::size_t>(mode)];
                }
                if (useSparse) {
                    sparse.push_back({gammaRow, gammaColumn, -correction});
                } else {
                    dense[static_cast<std::size_t>(gammaRow)
                        * static_cast<std::size_t>(model_.interfaceDofs)
                        + static_cast<std::size_t>(gammaColumn)] -= correction;
                }
            }
        }
    }
    assemblySeconds_ = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - assemblyStart).count();

    const auto factorStart = std::chrono::steady_clock::now();
    if (useSparse) {
        interfaceNonzeros_ = sparse.size();
        sparseFactor_ = std::make_unique<SparseFactor>(model_.interfaceDofs, sparse);
    } else {
        interfaceNonzeros_ = dense.size();
        schurFactor_ = factorDenseSymmetric(dense, model_.interfaceDofs);
    }
    factorizationSeconds_ += std::chrono::duration<double>(
        std::chrono::steady_clock::now() - factorStart).count();
}

void LocalReducedSchurSolver::initializeMatrixFree(
    const Mesh& mesh,
    const CaseConfig& physics,
    const ddm_schur::InterfacePartition& partition,
    const ddm_schur::Options& options,
    const std::filesystem::path& outputDirectory)
{
    if (options.interfaceKrylov != "fgmres"
        && options.interfaceKrylov != "pcg") {
        throw std::runtime_error(
            "[Local ROM] Interface Krylov solver must be fgmres or pcg.");
    }
    if (model_.interfaceDofs
            != static_cast<int>(partition.interfaceGlobalDofs.size())
        || model_.interfaceGlobalDofs != partition.interfaceGlobalDofs) {
        throw std::runtime_error(
            "[Local ROM] Matrix-free Schur partition/model ordering mismatch.");
    }
    const auto setupStart = std::chrono::steady_clock::now();
    interfaceSolver_ = options.interfaceKrylov == "pcg"
        ? "matrix-free-pcg-protected" : "matrix-free-fgmres";
    matrixFree_ = std::make_unique<MatrixFreeData>();
    matrixFree_->options = options;
    matrixFree_->localApplyWorkspaces.resize(model_.subdomains.size());

    std::vector<InterfaceEntry> entries = model_.interfaceEntries;
    std::sort(entries.begin(), entries.end(), [](const InterfaceEntry& left,
                                                 const InterfaceEntry& right) {
        return left.row < right.row
            || (left.row == right.row && left.column < right.column);
    });
    matrixFree_->rowPtr.assign(
        static_cast<std::size_t>(model_.interfaceDofs + 1), 0);
    for (std::size_t begin = 0; begin < entries.size();) {
        std::size_t end = begin + 1;
        double value = entries[begin].value;
        while (end < entries.size()
               && entries[end].row == entries[begin].row
               && entries[end].column == entries[begin].column) {
            value += entries[end].value;
            ++end;
        }
        if (value != 0.0) {
            if (entries[begin].row < 0 || entries[begin].row >= model_.interfaceDofs
                || entries[begin].column < 0
                || entries[begin].column >= model_.interfaceDofs) {
                throw std::runtime_error(
                    "[Local ROM] Matrix-free interface entry is out of range.");
            }
            matrixFree_->colInd.push_back(entries[begin].column);
            matrixFree_->values.push_back(value);
            ++matrixFree_->rowPtr[static_cast<std::size_t>(entries[begin].row + 1)];
        }
        begin = end;
    }
    std::partial_sum(matrixFree_->rowPtr.begin(), matrixFree_->rowPtr.end(),
                     matrixFree_->rowPtr.begin());
    interfaceNonzeros_ = matrixFree_->values.size();

    if (options.proxyEnabled) {
        matrixFree_->proxy = std::make_unique<ddm_schur::SchurProxyPreconditioner>(
            mesh, physics, partition,
            options.proxyHighConductivityThreshold,
            options.proxyUseMaterialConnectivity,
            options.proxyRing,
            options.proxyBlockSize,
            options.localPardisoThreads,
            options.proxyValidateBlockEquivalence,
            options.proxyCacheEnabled,
            options.proxyCachePath,
            outputDirectory,
            [&](const std::vector<double>& input, std::vector<double>& output) {
                applyMatrixFree(input, output, false);
            },
            [&](const std::vector<int>& probeColors,
                int firstColor,
                int rightHandSides,
                const ddm_schur::ExactSchurResponseConsumer& consumeResponse) {
                std::vector<double> probe(
                    static_cast<std::size_t>(model_.interfaceDofs), 0.0);
                for (int rhs = 0; rhs < rightHandSides; ++rhs) {
                    const int color = firstColor + rhs;
                    for (int gamma = 0; gamma < model_.interfaceDofs; ++gamma) {
                        probe[static_cast<std::size_t>(gamma)] =
                            probeColors[static_cast<std::size_t>(gamma)] == color
                                ? 1.0 : 0.0;
                    }
                    std::vector<double> response;
                    applyMatrixFree(probe, response, false);
                    consumeResponse(rhs, response);
                }
            });
    }

    if (!options.proxyDisableCoarse) {
        const auto appendDomainModes = [&](const SubdomainModel& local) {
            const std::size_t count = local.interfaceIndices.size();
            if (count == 0) {
                return;
            }
            MatrixFreeData::CoarseVector constant;
            constant.indices = local.interfaceIndices;
            constant.values.assign(count, 1.0 / std::sqrt(static_cast<double>(count)));
            matrixFree_->coarseBasis.push_back(std::move(constant));

            std::vector<std::vector<double>> linearModes;
            const auto appendAxis = [&](int axis) {
                double mean = 0.0;
                for (int global : local.interfaceGlobalDofs) {
                    const Vec3& point = mesh.nodes[static_cast<std::size_t>(global)].p;
                    mean += axis == 0 ? point.x : (axis == 1 ? point.y : point.z);
                }
                mean /= static_cast<double>(count);
                std::vector<double> values;
                values.reserve(count);
                for (int global : local.interfaceGlobalDofs) {
                    const Vec3& point = mesh.nodes[static_cast<std::size_t>(global)].p;
                    values.push_back(
                        (axis == 0 ? point.x : (axis == 1 ? point.y : point.z)) - mean);
                }
                for (const std::vector<double>& previous : linearModes) {
                    const double projection = vectorDot(values, previous);
                    for (std::size_t i = 0; i < count; ++i) {
                        values[i] -= projection * previous[i];
                    }
                }
                const double norm = std::sqrt(std::max(0.0, vectorDot(values, values)));
                if (norm <= 1.0e-12 * std::sqrt(static_cast<double>(count))) {
                    return;
                }
                for (double& value : values) {
                    value /= norm;
                }
                MatrixFreeData::CoarseVector linear;
                linear.indices = local.interfaceIndices;
                linear.values = values;
                matrixFree_->coarseBasis.push_back(std::move(linear));
                linearModes.push_back(std::move(values));
            };
            if (options.coarseLinearXY) {
                appendAxis(0);
                appendAxis(1);
            }
            if (options.coarseLinearZ) {
                appendAxis(2);
            }
        };
        for (const SubdomainModel& local : model_.subdomains) {
            appendDomainModes(local);
        }
    }

    if (options.interfaceOperatorCoarseRank < 0
        || options.interfaceOperatorCoarseSweeps <= 0) {
        throw std::runtime_error(
            "[Local ROM] Interface operator coarse rank/sweeps are invalid.");
    }
    const int rawGeometricDimension =
        static_cast<int>(matrixFree_->coarseBasis.size());
    matrixFree_->geometricCoarseDimension = rawGeometricDimension;
    const bool buildOperatorCoarse = options.interfaceOperatorCoarseRank > 0;
    if (buildOperatorCoarse
        && (!matrixFree_->proxy || rawGeometricDimension == 0)) {
        throw std::runtime_error(
            "[Local ROM] Interface operator coarse modes require the Schur "
            "proxy and geometric coarse modes.");
    }

    std::uint64_t operatorSignature = UINT64_C(1469598103934665603);
    if (buildOperatorCoarse) {
        constexpr int operatorCoarsePolicyVersion = 2;
        appendHash(operatorSignature, operatorCoarsePolicyVersion);
        appendHash(operatorSignature, model_.globalDofs);
        appendHash(operatorSignature, model_.interfaceDofs);
        appendHash(operatorSignature, model_.totalLocalRank);
        appendHash(operatorSignature, model_.fingerprints.mesh);
        appendHash(operatorSignature, model_.fingerprints.system);
        appendHash(operatorSignature, model_.fingerprints.interfaceOrdering);
        appendHash(operatorSignature, options.coarseLinearXY);
        appendHash(operatorSignature, options.coarseLinearZ);
        appendHash(operatorSignature, options.interfaceOperatorCoarseRank);
        appendHash(operatorSignature, options.interfaceOperatorCoarseSweeps);
        appendHash(operatorSignature, matrixFree_->proxy->valueHash());
        appendHash(operatorSignature, matrixFree_->rowPtr);
        appendHash(operatorSignature, matrixFree_->colInd);
        appendHash(operatorSignature, matrixFree_->values);
        const auto appendSampled = [&](const std::vector<double>& values) {
            appendHash(operatorSignature, values.size());
            if (values.empty()) {
                return;
            }
            const std::size_t stride = std::max<std::size_t>(
                1, values.size() / 4096);
            for (std::size_t i = 0; i < values.size(); i += stride) {
                appendHash(operatorSignature, values[i]);
            }
            appendHash(operatorSignature, values.back());
        };
        for (const SubdomainModel& local : model_.subdomains) {
            const SubdomainModel& data = templatePayload(model_, local);
            appendHash(operatorSignature, local.subdomain);
            appendHash(operatorSignature, local.rank);
            appendHash(operatorSignature, local.boundaryInterfaceFingerprint);
            appendHash(operatorSignature, local.interfaceIndices);
            appendHash(operatorSignature, data.reducedInterior);
            appendSampled(data.reducedInteriorInterface);
            appendSampled(data.reducedInterfaceInterior);
        }
    }

    const std::filesystem::path operatorCachePath =
        options.interfaceOperatorCoarseCachePath;
    const auto loadOperatorCache = [&]() {
        if (!buildOperatorCoarse || operatorCachePath.empty()
            || !std::filesystem::exists(operatorCachePath)) {
            return false;
        }
        const auto start = std::chrono::steady_clock::now();
        bool valid = false;
        try {
            std::ifstream input(operatorCachePath, std::ios::binary);
            std::uint64_t magic = 0;
            int version = 0;
            std::uint64_t signature = 0;
            int interfaceDofs = 0;
            int geometricDimension = 0;
            int requestedRank = 0;
            int requestedSweeps = 0;
            int basisCount = 0;
            valid = input
                && readBinary(input, magic)
                && readBinary(input, version)
                && readBinary(input, signature)
                && readBinary(input, interfaceDofs)
                && readBinary(input, geometricDimension)
                && readBinary(input, requestedRank)
                && readBinary(input, requestedSweeps)
                && readBinary(input, basisCount)
                && magic == operatorCoarseCacheMagic
                && version == operatorCoarseCacheVersion
                && signature == operatorSignature
                && interfaceDofs == model_.interfaceDofs
                && geometricDimension == rawGeometricDimension
                && requestedRank == options.interfaceOperatorCoarseRank
                && requestedSweeps == options.interfaceOperatorCoarseSweeps
                && basisCount >= geometricDimension
                && basisCount <= geometricDimension + requestedRank;
            std::vector<MatrixFreeData::CoarseVector> basis;
            std::vector<MatrixFreeData::CoarseVector> images;
            if (valid) {
                basis.resize(static_cast<std::size_t>(basisCount));
                images.resize(static_cast<std::size_t>(basisCount));
                for (MatrixFreeData::CoarseVector& value : basis) {
                    valid = readBinary(input, value.dense)
                        && readBinaryVector(
                            input, value.indices,
                            static_cast<std::uint64_t>(model_.interfaceDofs))
                        && readBinaryVector(
                            input, value.values,
                            static_cast<std::uint64_t>(model_.interfaceDofs))
                        && (value.dense
                            ? value.indices.empty()
                                && value.values.size() == static_cast<std::size_t>(
                                    model_.interfaceDofs)
                            : value.indices.size() == value.values.size());
                    if (!valid) break;
                }
                for (MatrixFreeData::CoarseVector& value : images) {
                    valid = readBinary(input, value.dense)
                        && readBinaryVector(
                            input, value.indices,
                            static_cast<std::uint64_t>(model_.interfaceDofs))
                        && readBinaryVector(
                            input, value.values,
                            static_cast<std::uint64_t>(model_.interfaceDofs))
                        && (value.dense
                            ? value.indices.empty()
                                && value.values.size() == static_cast<std::size_t>(
                                    model_.interfaceDofs)
                            : value.indices.size() == value.values.size());
                    if (!valid) break;
                }
            }
            const auto readFactor = [&](DenseSymmetricFactor& factor) {
                return readBinary(input, factor.size)
                    && readBinary(input, factor.cholesky)
                    && readBinaryVector(
                        input, factor.lower,
                        static_cast<std::uint64_t>(basisCount)
                            * static_cast<std::uint64_t>(basisCount))
                    && readBinaryVector(
                        input, factor.diagonal,
                        static_cast<std::uint64_t>(basisCount))
                    && factor.size == basisCount
                    && factor.lower.size() == static_cast<std::size_t>(
                        basisCount * basisCount)
                    && (factor.cholesky || factor.diagonal.size()
                        == static_cast<std::size_t>(basisCount));
            };
            DenseSymmetricFactor factor;
            DenseSymmetricFactor predictorFactor;
            if (valid) {
                valid = readFactor(factor) && readFactor(predictorFactor);
            }
            if (valid) {
                for (const MatrixFreeData::CoarseVector* group :
                     {basis.data(), images.data()}) {
                    for (int column = 0; column < basisCount && valid; ++column) {
                        const MatrixFreeData::CoarseVector& value = group[column];
                        valid = value.dense || std::all_of(
                                value.indices.begin(), value.indices.end(),
                                [&](int index) {
                                    return index >= 0
                                        && index < model_.interfaceDofs;
                                });
                    }
                }
            }
            if (valid) {
                matrixFree_->coarseBasis = std::move(basis);
                matrixFree_->coarseImages = std::move(images);
                matrixFree_->coarseFactor = std::move(factor);
                matrixFree_->predictorFactor = std::move(predictorFactor);
                matrixFree_->operatorCoarseDimension =
                    basisCount - geometricDimension;
            }
        } catch (const std::exception&) {
            valid = false;
        }
        matrixFree_->operatorCoarseCacheLoadSeconds =
            std::chrono::duration<double>(
                std::chrono::steady_clock::now() - start).count();
        return valid;
    };

    matrixFree_->operatorCoarseCacheHit = loadOperatorCache();
    if (!matrixFree_->operatorCoarseCacheHit) {
        const auto operatorSetupStart = std::chrono::steady_clock::now();
        const auto sparseToDense = [&](const MatrixFreeData::CoarseVector& sparse) {
            return sparse.expand(static_cast<std::size_t>(model_.interfaceDofs));
        };
        const auto appendOrthonormal = [&](std::vector<double> candidate) {
            const double originalNorm = std::sqrt(std::max(
                0.0, vectorDot(candidate, candidate)));
            if (!(originalNorm > 0.0) || !std::isfinite(originalNorm)) {
                return false;
            }
            for (int pass = 0; pass < 2; ++pass) {
                for (const MatrixFreeData::CoarseVector& existing :
                     matrixFree_->coarseBasis) {
                    const double projection = existing.dot(candidate);
                    existing.addScaled(-projection, candidate);
                }
            }
            const double norm = std::sqrt(std::max(
                0.0, vectorDot(candidate, candidate)));
            if (!(norm > 1.0e-10 * originalNorm) || !std::isfinite(norm)) {
                return false;
            }
            const double inverseNorm = 1.0 / norm;
            double maximum = 0.0;
            for (double& value : candidate) {
                value *= inverseNorm;
                maximum = std::max(maximum, std::abs(value));
            }
            MatrixFreeData::CoarseVector sparse;
            const double dropTolerance = 1.0e-14 * maximum;
            std::size_t retained = 0;
            for (double value : candidate) {
                retained += std::abs(value) > dropTolerance ? 1 : 0;
            }
            if (retained * 2 > candidate.size()) {
                sparse.dense = true;
                sparse.values = std::move(candidate);
                matrixFree_->coarseBasis.push_back(std::move(sparse));
                return true;
            }
            sparse.indices.reserve(retained);
            sparse.values.reserve(retained);
            for (int gamma = 0; gamma < model_.interfaceDofs; ++gamma) {
                const double value = candidate[static_cast<std::size_t>(gamma)];
                if (std::abs(value) > dropTolerance) {
                    sparse.indices.push_back(gamma);
                    sparse.values.push_back(value);
                }
            }
            matrixFree_->coarseBasis.push_back(std::move(sparse));
            return true;
        };

        if (buildOperatorCoarse) {
            std::vector<MatrixFreeData::CoarseVector> geometric =
                std::move(matrixFree_->coarseBasis);
            matrixFree_->coarseBasis.clear();
            matrixFree_->coarseBasis.reserve(static_cast<std::size_t>(
                rawGeometricDimension + options.interfaceOperatorCoarseRank));
            for (const MatrixFreeData::CoarseVector& value : geometric) {
                if (!appendOrthonormal(sparseToDense(value))) {
                    throw std::runtime_error(
                        "[Local ROM] Geometric coarse basis lost rank during "
                        "operator-space orthogonalization.");
                }
            }
            matrixFree_->geometricCoarseDimension =
                static_cast<int>(matrixFree_->coarseBasis.size());
            int sourceBegin = 0;
            int sourceEnd = matrixFree_->geometricCoarseDimension;
            for (int sweep = 0;
                 sweep < options.interfaceOperatorCoarseSweeps
                    && matrixFree_->operatorCoarseDimension
                        < options.interfaceOperatorCoarseRank;
                 ++sweep) {
                std::vector<std::pair<double, std::vector<double>>> candidates;
                candidates.reserve(static_cast<std::size_t>(
                    std::max(0, sourceEnd - sourceBegin)
                    + (sweep == 0 ? std::min(16,
                        std::max(4, options.interfaceOperatorCoarseRank)) : 0)));
                const auto appendFilteredCandidate =
                    [&](std::vector<double> input) {
                        std::vector<double> exactImage;
                        std::vector<double> proxyCorrection;
                        applyMatrixFree(input, exactImage, false);
                        matrixFree_->proxy->solve(exactImage, proxyCorrection);
                        for (std::size_t i = 0; i < input.size(); ++i) {
                            input[i] -= proxyCorrection[i];
                        }
                        const double score = std::sqrt(std::max(
                            0.0, vectorDot(input, input)));
                        if (score > 0.0 && std::isfinite(score)) {
                            candidates.emplace_back(score, std::move(input));
                        }
                    };
                for (int source = sourceBegin; source < sourceEnd; ++source) {
                    appendFilteredCandidate(sparseToDense(
                        matrixFree_->coarseBasis[static_cast<std::size_t>(source)]));
                }
                if (sweep == 0) {
                    const int probeCount = std::min(
                        16, std::max(4, options.interfaceOperatorCoarseRank));
                    const double scale = 1.0 / std::sqrt(
                        static_cast<double>(model_.interfaceDofs));
                    for (int probe = 0; probe < probeCount; ++probe) {
                        std::vector<double> input(
                            static_cast<std::size_t>(model_.interfaceDofs), 0.0);
                        for (int gamma = 0; gamma < model_.interfaceDofs; ++gamma) {
                            std::uint64_t state =
                                (static_cast<std::uint64_t>(gamma) + 1)
                                    * UINT64_C(0x9e3779b97f4a7c15)
                                + (static_cast<std::uint64_t>(probe) + 1)
                                    * UINT64_C(0xbf58476d1ce4e5b9);
                            state ^= state >> 30;
                            state *= UINT64_C(0xbf58476d1ce4e5b9);
                            state ^= state >> 27;
                            state *= UINT64_C(0x94d049bb133111eb);
                            state ^= state >> 31;
                            input[static_cast<std::size_t>(gamma)] =
                                (state & 1) == 0 ? scale : -scale;
                        }
                        appendFilteredCandidate(std::move(input));
                    }
                }
                std::stable_sort(
                    candidates.begin(), candidates.end(),
                    [](const auto& left, const auto& right) {
                        return left.first > right.first;
                    });
                const int acceptedBegin =
                    static_cast<int>(matrixFree_->coarseBasis.size());
                const int targetAfterSweep = std::min(
                    options.interfaceOperatorCoarseRank,
                    (options.interfaceOperatorCoarseRank * (sweep + 1)
                        + options.interfaceOperatorCoarseSweeps - 1)
                        / options.interfaceOperatorCoarseSweeps);
                for (auto& candidate : candidates) {
                    if (matrixFree_->operatorCoarseDimension
                        >= targetAfterSweep) {
                        break;
                    }
                    if (appendOrthonormal(std::move(candidate.second))) {
                        ++matrixFree_->operatorCoarseDimension;
                    }
                }
                sourceBegin = acceptedBegin;
                sourceEnd = static_cast<int>(matrixFree_->coarseBasis.size());
                if (sourceBegin == sourceEnd) {
                    break;
                }
            }
            matrixFree_->proxy->resetRuntimeCounters();
        }

        const int dimension = static_cast<int>(matrixFree_->coarseBasis.size());
        if (dimension > 0) {
            std::vector<double> coarseMatrix(
                static_cast<std::size_t>(dimension * dimension), 0.0);
            matrixFree_->coarseImages.clear();
            matrixFree_->coarseImages.reserve(static_cast<std::size_t>(dimension));
            for (int column = 0; column < dimension; ++column) {
                std::vector<double> basis = sparseToDense(
                    matrixFree_->coarseBasis[static_cast<std::size_t>(column)]);
                std::vector<double> image;
                applyMatrixFree(basis, image, false);
                for (int row = 0; row < dimension; ++row) {
                    const MatrixFreeData::CoarseVector& test =
                        matrixFree_->coarseBasis[static_cast<std::size_t>(row)];
                    const double value = test.dot(image);
                    coarseMatrix[static_cast<std::size_t>(
                        row * dimension + column)] = value;
                }
                MatrixFreeData::CoarseVector sparseImage;
                sparseImage.dense = true;
                sparseImage.values = std::move(image);
                matrixFree_->coarseImages.push_back(std::move(sparseImage));
            }
            matrixFree_->coarseFactor = factorDenseSymmetric(
                coarseMatrix, dimension);
            if (matrixFree_->operatorCoarseDimension > 0) {
                std::vector<double> predictorMatrix(
                    static_cast<std::size_t>(dimension * dimension), 0.0);
                for (int row = 0; row < dimension; ++row) {
                    for (int column = row; column < dimension; ++column) {
                        const double value = matrixFree_->coarseImages[
                            static_cast<std::size_t>(row)].dot(
                                matrixFree_->coarseImages[
                                    static_cast<std::size_t>(column)].values);
                        predictorMatrix[static_cast<std::size_t>(
                            row * dimension + column)] = value;
                        predictorMatrix[static_cast<std::size_t>(
                            column * dimension + row)] = value;
                    }
                }
                try {
                    matrixFree_->predictorFactor = factorDenseSymmetric(
                        predictorMatrix, dimension);
                } catch (const std::runtime_error&) {
                    double scale = 0.0;
                    for (int row = 0; row < dimension; ++row) {
                        scale = std::max(scale, predictorMatrix[
                            static_cast<std::size_t>(row * dimension + row)]);
                    }
                    const double shift = 1.0e-12 * scale;
                    for (int row = 0; row < dimension; ++row) {
                        predictorMatrix[static_cast<std::size_t>(
                            row * dimension + row)] += shift;
                    }
                    matrixFree_->predictorFactor = factorDenseSymmetric(
                        predictorMatrix, dimension);
                }
            }
        }
        matrixFree_->operatorCoarseSetupSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - operatorSetupStart).count();

        if (buildOperatorCoarse && !operatorCachePath.empty()) {
            const auto cacheSaveStart = std::chrono::steady_clock::now();
            std::filesystem::path temporary = operatorCachePath;
            temporary += ".tmp";
            try {
                if (!operatorCachePath.parent_path().empty()) {
                    std::filesystem::create_directories(
                        operatorCachePath.parent_path());
                }
                std::ofstream output(
                    temporary, std::ios::binary | std::ios::trunc);
                const int basisCount =
                    static_cast<int>(matrixFree_->coarseBasis.size());
                writeBinary(output, operatorCoarseCacheMagic);
                writeBinary(output, operatorCoarseCacheVersion);
                writeBinary(output, operatorSignature);
                writeBinary(output, model_.interfaceDofs);
                writeBinary(output, matrixFree_->geometricCoarseDimension);
                writeBinary(output, options.interfaceOperatorCoarseRank);
                writeBinary(output, options.interfaceOperatorCoarseSweeps);
                writeBinary(output, basisCount);
                for (const MatrixFreeData::CoarseVector& value :
                     matrixFree_->coarseBasis) {
                    writeBinary(output, value.dense);
                    writeBinaryVector(output, value.indices);
                    writeBinaryVector(output, value.values);
                }
                for (const MatrixFreeData::CoarseVector& value :
                     matrixFree_->coarseImages) {
                    writeBinary(output, value.dense);
                    writeBinaryVector(output, value.indices);
                    writeBinaryVector(output, value.values);
                }
                writeBinary(output, matrixFree_->coarseFactor.size);
                writeBinary(output, matrixFree_->coarseFactor.cholesky);
                writeBinaryVector(output, matrixFree_->coarseFactor.lower);
                writeBinaryVector(output, matrixFree_->coarseFactor.diagonal);
                writeBinary(output, matrixFree_->predictorFactor.size);
                writeBinary(output, matrixFree_->predictorFactor.cholesky);
                writeBinaryVector(output, matrixFree_->predictorFactor.lower);
                writeBinaryVector(output, matrixFree_->predictorFactor.diagonal);
                output.close();
                if (!output) {
                    throw std::runtime_error("cache write failed");
                }
                std::error_code error;
                std::filesystem::rename(temporary, operatorCachePath, error);
                if (error) {
                    std::filesystem::remove(operatorCachePath, error);
                    error.clear();
                    std::filesystem::rename(temporary, operatorCachePath, error);
                }
                if (error) {
                    throw std::runtime_error(error.message());
                }
            } catch (const std::exception& error) {
                std::error_code ignored;
                std::filesystem::remove(temporary, ignored);
                std::cerr << "[Local ROM] Warning: interface operator coarse "
                          << "cache was not saved: " << error.what() << '\n';
            }
            matrixFree_->operatorCoarseCacheSaveSeconds =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - cacheSaveStart).count();
        }
    }

    const int coarseDimension = static_cast<int>(matrixFree_->coarseBasis.size());

    assemblySeconds_ = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - setupStart).count();
    factorizationSeconds_ += matrixFree_->proxy
        ? matrixFree_->proxy->symbolicSeconds()
            + matrixFree_->proxy->numericalSeconds()
        : 0.0;
    std::cout << "[Local ROM] matrix-free full-interface Schur: dofs="
              << model_.interfaceDofs << ", A_GammaGamma nnz="
              << interfaceNonzeros_ << ", coarse=" << coarseDimension
              << " (geometry=" << matrixFree_->geometricCoarseDimension
              << ", operator=" << matrixFree_->operatorCoarseDimension
              << ", operator cache="
              << (matrixFree_->operatorCoarseCacheHit ? "hit" : "miss") << ')'
              << ", proxy colors=" << proxyColors() << '\n';
}

LocalReducedSchurSolver::~LocalReducedSchurSolver() = default;

void LocalReducedSchurSolver::applyMatrixFree(
    const std::vector<double>& input,
    std::vector<double>& output,
    bool countMatvec)
{
    if (!matrixFree_
        || input.size() != static_cast<std::size_t>(model_.interfaceDofs)) {
        throw std::runtime_error(
            "[Local ROM] Invalid matrix-free Schur apply dimensions.");
    }
    output.resize(input.size());
    const bool parallelInterfaceMultiply =
        matrixFree_->options.localSolveThreads > 1
        && model_.interfaceDofs >= 4096;
#pragma omp parallel for num_threads(matrixFree_->options.localSolveThreads) \
    if(parallelInterfaceMultiply) schedule(static)
    for (int row = 0; row < model_.interfaceDofs; ++row) {
        double value = 0.0;
        for (int offset = matrixFree_->rowPtr[static_cast<std::size_t>(row)];
             offset < matrixFree_->rowPtr[static_cast<std::size_t>(row + 1)];
             ++offset) {
            value += matrixFree_->values[static_cast<std::size_t>(offset)]
                * input[static_cast<std::size_t>(
                    matrixFree_->colInd[static_cast<std::size_t>(offset)])];
        }
        output[static_cast<std::size_t>(row)] = value;
    }

    const int domainCount = static_cast<int>(model_.subdomains.size());
#pragma omp parallel for num_threads(matrixFree_->options.localSolveThreads) \
    if(matrixFree_->options.localSolveThreads > 1) schedule(static)
    for (int slot = 0; slot < domainCount; ++slot) {
        const SubdomainModel& local = model_.subdomains[static_cast<std::size_t>(slot)];
        const SubdomainModel& data = templatePayload(model_, local);
        MatrixFreeData::LocalApplyWorkspace& workspace =
            matrixFree_->localApplyWorkspaces[static_cast<std::size_t>(slot)];
        std::vector<double>& interfaceInput = workspace.interfaceInput;
        interfaceInput.resize(local.interfaceIndices.size());
        for (std::size_t localGamma = 0;
             localGamma < local.interfaceIndices.size(); ++localGamma) {
            interfaceInput[localGamma] = input[static_cast<std::size_t>(
                local.interfaceIndices[localGamma])];
        }
        std::vector<double>& reducedProduct = workspace.reducedProduct;
        reducedProduct.resize(static_cast<std::size_t>(local.rank));
        for (int mode = 0; mode < local.rank; ++mode) {
            double value = 0.0;
            for (std::size_t localGamma = 0;
                 localGamma < local.interfaceIndices.size(); ++localGamma) {
                value += data.reducedInteriorInterface[static_cast<std::size_t>(
                    mode * local.localInterfaceDofs) + localGamma]
                    * interfaceInput[localGamma];
            }
            reducedProduct[static_cast<std::size_t>(mode)] = value;
        }
        solveDenseSymmetric(
            localFactors_[static_cast<std::size_t>(slot)], reducedProduct);
        std::vector<double>& correction = workspace.correction;
        correction.resize(local.interfaceIndices.size());
        for (std::size_t localGamma = 0;
             localGamma < local.interfaceIndices.size(); ++localGamma) {
            double value = 0.0;
            for (int mode = 0; mode < local.rank; ++mode) {
                value += data.reducedInterfaceInterior[static_cast<std::size_t>(
                    localGamma * static_cast<std::size_t>(local.rank) + mode)]
                    * reducedProduct[static_cast<std::size_t>(mode)];
            }
            correction[localGamma] = value;
        }
    }
    // Junction DOFs can be touched by more than one subdomain. Keep the
    // accumulation deterministic and free of concurrent global writes.
    for (int slot = 0; slot < domainCount; ++slot) {
        const SubdomainModel& local =
            model_.subdomains[static_cast<std::size_t>(slot)];
        const std::vector<double>& correction =
            matrixFree_->localApplyWorkspaces[
                static_cast<std::size_t>(slot)].correction;
        for (std::size_t localGamma = 0;
             localGamma < local.interfaceIndices.size(); ++localGamma) {
            output[static_cast<std::size_t>(
                local.interfaceIndices[localGamma])] -= correction[localGamma];
        }
    }
    if (countMatvec) {
        ++matrixFree_->matvecs;
    }
}

void LocalReducedSchurSolver::applyMatrixFreePreconditioner(
    const std::vector<double>& residual,
    std::vector<double>& result)
{
    if (!matrixFree_
        || residual.size() != static_cast<std::size_t>(model_.interfaceDofs)) {
        throw std::runtime_error(
            "[Local ROM] Invalid matrix-free preconditioner dimensions.");
    }
    const auto solveCoarse = [&](std::vector<double>& rhs) {
        if (rhs.empty()) {
            return;
        }
        const auto start = std::chrono::steady_clock::now();
        solveDenseSymmetric(matrixFree_->coarseFactor, rhs);
        matrixFree_->coarseSolveSeconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
    };
    const auto applyLevelOne = [&](const std::vector<double>& rhs,
                                   std::vector<double>& correction) {
        if (matrixFree_->proxy) {
            matrixFree_->proxy->solve(rhs, correction);
        } else {
            correction = rhs;
        }
    };

    const std::size_t coarseDimension = matrixFree_->coarseBasis.size();
    if (coarseDimension == 0) {
        applyLevelOne(residual, result);
        return;
    }
    if (matrixFree_->operatorCoarseDimension == 0) {
        std::vector<double>& coarseRhs = matrixFree_->preconditionerCoarseRhs;
        coarseRhs.assign(coarseDimension, 0.0);
        for (std::size_t coarse = 0; coarse < coarseDimension; ++coarse) {
            coarseRhs[coarse] = matrixFree_->coarseBasis[coarse].dot(residual);
        }
        solveCoarse(coarseRhs);
        std::vector<double>& coarseCorrection =
            matrixFree_->preconditionerCoarseCorrection;
        std::vector<double>& coarseImage =
            matrixFree_->preconditionerCoarseImage;
        coarseCorrection.assign(residual.size(), 0.0);
        coarseImage.assign(residual.size(), 0.0);
        for (std::size_t coarse = 0; coarse < coarseDimension; ++coarse) {
            matrixFree_->coarseBasis[coarse].addScaled(
                coarseRhs[coarse], coarseCorrection);
            matrixFree_->coarseImages[coarse].addScaled(
                coarseRhs[coarse], coarseImage);
        }
        std::vector<double>& projectedResidual =
            matrixFree_->preconditionerProjectedResidual;
        projectedResidual.resize(residual.size());
        for (std::size_t i = 0; i < residual.size(); ++i) {
            projectedResidual[i] = residual[i] - coarseImage[i];
        }
        std::vector<double>& localCorrection =
            matrixFree_->preconditionerLocalCorrection;
        applyLevelOne(projectedResidual, localCorrection);
        std::vector<double>& leftRhs = matrixFree_->preconditionerLeftRhs;
        leftRhs.assign(coarseDimension, 0.0);
        for (std::size_t coarse = 0; coarse < coarseDimension; ++coarse) {
            leftRhs[coarse] = matrixFree_->coarseImages[coarse].dot(
                localCorrection);
        }
        solveCoarse(leftRhs);
        result.resize(residual.size());
        for (std::size_t i = 0; i < result.size(); ++i) {
            result[i] = coarseCorrection[i] + localCorrection[i];
        }
        for (std::size_t coarse = 0; coarse < coarseDimension; ++coarse) {
            matrixFree_->coarseBasis[coarse].addScaled(
                -leftRhs[coarse], result);
        }
        return;
    }

    std::vector<double>& coarseRhs = matrixFree_->preconditionerCoarseRhs;
    coarseRhs.assign(coarseDimension, 0.0);
    const bool parallelCoarse = matrixFree_->options.localSolveThreads > 1
        && residual.size() >= 4096 && coarseDimension >= 8;
#pragma omp parallel for num_threads(matrixFree_->options.localSolveThreads) \
    if(parallelCoarse) schedule(static)
    for (int coarse = 0; coarse < static_cast<int>(coarseDimension); ++coarse) {
        const MatrixFreeData::CoarseVector& basis =
            matrixFree_->coarseBasis[static_cast<std::size_t>(coarse)];
        coarseRhs[static_cast<std::size_t>(coarse)] = basis.dot(residual);
    }
    solveCoarse(coarseRhs);
    std::vector<double>& coarseCorrection =
        matrixFree_->preconditionerCoarseCorrection;
    std::vector<double>& coarseImage =
        matrixFree_->preconditionerCoarseImage;
    coarseCorrection.assign(residual.size(), 0.0);
    coarseImage.assign(residual.size(), 0.0);
    for (std::size_t coarse = 0; coarse < coarseDimension; ++coarse) {
        const double coefficient = coarseRhs[coarse];
        const MatrixFreeData::CoarseVector& basis =
            matrixFree_->coarseBasis[coarse];
        if (!basis.dense) {
            basis.addScaled(coefficient, coarseCorrection);
        }
        const MatrixFreeData::CoarseVector& image =
            matrixFree_->coarseImages[coarse];
        if (!image.dense) {
            image.addScaled(coefficient, coarseImage);
        }
    }
#pragma omp parallel for num_threads(matrixFree_->options.localSolveThreads) \
    if(parallelCoarse) schedule(static)
    for (int row = 0; row < static_cast<int>(residual.size()); ++row) {
        double basisValue = 0.0;
        double imageValue = 0.0;
        for (std::size_t coarse = 0; coarse < coarseDimension; ++coarse) {
            const double coefficient = coarseRhs[coarse];
            const MatrixFreeData::CoarseVector& basis =
                matrixFree_->coarseBasis[coarse];
            if (basis.dense) {
                basisValue += coefficient
                    * basis.values[static_cast<std::size_t>(row)];
            }
            const MatrixFreeData::CoarseVector& image =
                matrixFree_->coarseImages[coarse];
            if (image.dense) {
                imageValue += coefficient
                    * image.values[static_cast<std::size_t>(row)];
            }
        }
        coarseCorrection[static_cast<std::size_t>(row)] += basisValue;
        coarseImage[static_cast<std::size_t>(row)] += imageValue;
    }
    std::vector<double>& projectedResidual =
        matrixFree_->preconditionerProjectedResidual;
    projectedResidual.resize(residual.size());
    for (std::size_t i = 0; i < residual.size(); ++i) {
        projectedResidual[i] = residual[i] - coarseImage[i];
    }
    std::vector<double>& localCorrection =
        matrixFree_->preconditionerLocalCorrection;
    applyLevelOne(projectedResidual, localCorrection);

    std::vector<double>& leftRhs = matrixFree_->preconditionerLeftRhs;
    leftRhs.assign(coarseDimension, 0.0);
#pragma omp parallel for num_threads(matrixFree_->options.localSolveThreads) \
    if(parallelCoarse) schedule(static)
    for (int coarse = 0; coarse < static_cast<int>(coarseDimension); ++coarse) {
        const MatrixFreeData::CoarseVector& image =
            matrixFree_->coarseImages[static_cast<std::size_t>(coarse)];
        leftRhs[static_cast<std::size_t>(coarse)] = image.dot(localCorrection);
    }
    solveCoarse(leftRhs);
    result.resize(residual.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = coarseCorrection[i] + localCorrection[i];
    }
    for (std::size_t coarse = 0; coarse < coarseDimension; ++coarse) {
        const MatrixFreeData::CoarseVector& basis =
            matrixFree_->coarseBasis[coarse];
        if (!basis.dense) {
            basis.addScaled(-leftRhs[coarse], result);
        }
    }
#pragma omp parallel for num_threads(matrixFree_->options.localSolveThreads) \
    if(parallelCoarse) schedule(static)
    for (int row = 0; row < static_cast<int>(result.size()); ++row) {
        double value = 0.0;
        for (std::size_t coarse = 0; coarse < coarseDimension; ++coarse) {
            const MatrixFreeData::CoarseVector& basis =
                matrixFree_->coarseBasis[coarse];
            if (basis.dense) {
                value += leftRhs[coarse]
                    * basis.values[static_cast<std::size_t>(row)];
            }
        }
        result[static_cast<std::size_t>(row)] -= value;
    }
}

double LocalReducedSchurSolver::symbolicAnalysisSeconds() const
{
    if (portCore_) {
        return portCore_->symbolicSeconds();
    }
    if (augmentedDirect_) {
        return augmentedDirect_->generalFactor
            ? augmentedDirect_->generalFactor->symbolicAnalysisSeconds()
            : augmentedDirect_->factor->symbolicAnalysisSeconds();
    }
    if (sparseFactor_) {
        return sparseFactor_->solver.symbolicAnalysisSeconds();
    }
    return matrixFree_ && matrixFree_->proxy
        ? matrixFree_->proxy->symbolicSeconds() : 0.0;
}

double LocalReducedSchurSolver::numericalFactorizationSeconds() const
{
    if (portCore_) {
        return portCore_->numericalSeconds();
    }
    if (augmentedDirect_) {
        return augmentedDirect_->generalFactor
            ? augmentedDirect_->generalFactor->numericalFactorizationSeconds()
            : augmentedDirect_->factor->numericalFactorizationSeconds();
    }
    if (sparseFactor_) {
        return sparseFactor_->solver.numericalFactorizationSeconds();
    }
    return matrixFree_ && matrixFree_->proxy
        ? matrixFree_->proxy->numericalSeconds() : factorizationSeconds_;
}

int LocalReducedSchurSolver::symbolicAnalysisCalls() const
{
    if (portCore_) {
        return portCore_->symbolicCalls();
    }
    if (augmentedDirect_) {
        return augmentedDirect_->generalFactor
            ? augmentedDirect_->generalFactor->symbolicAnalysisCalls()
            : augmentedDirect_->factor->symbolicAnalysisCalls();
    }
    if (sparseFactor_) {
        return sparseFactor_->solver.symbolicAnalysisCalls();
    }
    return matrixFree_ && matrixFree_->proxy
        && matrixFree_->proxy->symbolicSeconds() > 0.0 ? 1 : 0;
}

int LocalReducedSchurSolver::numericalFactorizationCalls() const
{
    if (portCore_) {
        return portCore_->numericalCalls();
    }
    if (augmentedDirect_) {
        return augmentedDirect_->generalFactor
            ? augmentedDirect_->generalFactor->numericalFactorizationCalls()
            : augmentedDirect_->factor->numericalFactorizationCalls();
    }
    if (sparseFactor_) {
        return sparseFactor_->solver.numericalFactorizationCalls();
    }
    return matrixFree_ && matrixFree_->proxy
        && matrixFree_->proxy->numericalSeconds() > 0.0 ? 1 : 0;
}

std::size_t LocalReducedSchurSolver::factorMemoryBytes() const
{
    if (portCore_) {
        return portCore_->memoryBytes();
    }
    if (augmentedDirect_) {
        return augmentedDirect_->generalFactor
            ? augmentedDirect_->generalFactor->memoryBytes()
            : augmentedDirect_->factor->memoryBytes();
    }
    if (sparseFactor_) {
        return sparseFactor_->solver.memoryBytes();
    }
    if (matrixFree_) {
        std::size_t bytes = matrixFree_->rowPtr.capacity() * sizeof(int)
            + matrixFree_->colInd.capacity() * sizeof(int)
            + matrixFree_->values.capacity() * sizeof(double)
            + proxyMemoryBytes()
            + matrixFree_->coarseFactor.lower.capacity() * sizeof(double)
            + matrixFree_->coarseFactor.diagonal.capacity() * sizeof(double)
            + matrixFree_->predictorFactor.lower.capacity() * sizeof(double)
            + matrixFree_->predictorFactor.diagonal.capacity() * sizeof(double);
        for (const MatrixFreeData::CoarseVector& value : matrixFree_->coarseBasis) {
            bytes += value.indices.capacity() * sizeof(int)
                + value.values.capacity() * sizeof(double);
        }
        for (const MatrixFreeData::CoarseVector& value : matrixFree_->coarseImages) {
            bytes += value.indices.capacity() * sizeof(int)
                + value.values.capacity() * sizeof(double);
        }
        return bytes;
    }
    return schurFactor_.lower.capacity() * sizeof(double)
        + schurFactor_.diagonal.capacity() * sizeof(double);
}

int LocalReducedSchurSolver::coarseDimension() const
{
    return matrixFree_
        ? static_cast<int>(matrixFree_->coarseBasis.size()) : 0;
}

int LocalReducedSchurSolver::geometricCoarseDimension() const
{
    return matrixFree_ ? matrixFree_->geometricCoarseDimension : 0;
}

int LocalReducedSchurSolver::operatorCoarseDimension() const
{
    return matrixFree_ ? matrixFree_->operatorCoarseDimension : 0;
}

bool LocalReducedSchurSolver::operatorCoarseCacheHit() const
{
    return matrixFree_ && matrixFree_->operatorCoarseCacheHit;
}

double LocalReducedSchurSolver::operatorCoarseSetupSeconds() const
{
    return matrixFree_ ? matrixFree_->operatorCoarseSetupSeconds : 0.0;
}

double LocalReducedSchurSolver::operatorCoarseCacheLoadSeconds() const
{
    return matrixFree_ ? matrixFree_->operatorCoarseCacheLoadSeconds : 0.0;
}

double LocalReducedSchurSolver::operatorCoarseCacheSaveSeconds() const
{
    return matrixFree_ ? matrixFree_->operatorCoarseCacheSaveSeconds : 0.0;
}

int LocalReducedSchurSolver::proxyColors() const
{
    return matrixFree_ && matrixFree_->proxy ? matrixFree_->proxy->colors() : 0;
}

int LocalReducedSchurSolver::proxyProbingApplies() const
{
    return matrixFree_ && matrixFree_->proxy
        ? matrixFree_->proxy->probingSchurApplies() : 0;
}

std::size_t LocalReducedSchurSolver::proxyNonzeros() const
{
    return matrixFree_ && matrixFree_->proxy ? matrixFree_->proxy->nnz() : 0;
}

double LocalReducedSchurSolver::proxySetupSeconds() const
{
    return matrixFree_ && matrixFree_->proxy
        ? matrixFree_->proxy->setupSeconds() : 0.0;
}

double LocalReducedSchurSolver::proxySymbolicSeconds() const
{
    return matrixFree_ && matrixFree_->proxy
        ? matrixFree_->proxy->symbolicSeconds() : 0.0;
}

double LocalReducedSchurSolver::proxyNumericalSeconds() const
{
    return matrixFree_ && matrixFree_->proxy
        ? matrixFree_->proxy->numericalSeconds() : 0.0;
}

std::size_t LocalReducedSchurSolver::proxyMemoryBytes() const
{
    return matrixFree_ && matrixFree_->proxy
        ? matrixFree_->proxy->memoryBytes() : 0;
}

bool LocalReducedSchurSolver::proxyMatrixCacheHit() const
{
    return matrixFree_ && matrixFree_->proxy
        && matrixFree_->proxy->matrixCacheHit();
}

bool LocalReducedSchurSolver::proxyFactorCacheHit() const
{
    return matrixFree_ && matrixFree_->proxy
        && matrixFree_->proxy->factorCacheHit();
}

bool LocalReducedSchurSolver::portCoreCacheHit() const
{
    return portCore_ && portCore_->cacheHit();
}

double LocalReducedSchurSolver::portCoreCacheLoadSeconds() const
{
    return portCore_ ? portCore_->cacheLoadSeconds() : 0.0;
}

double LocalReducedSchurSolver::portCoreCacheSaveSeconds() const
{
    return portCore_ ? portCore_->cacheSaveSeconds() : 0.0;
}

std::size_t LocalReducedSchurSolver::portCoreCacheBytes() const
{
    return portCore_ ? portCore_->cacheBytes() : 0;
}

double LocalReducedSchurSolver::portCorePartitionSeconds() const
{
    return portCore_ ? portCore_->partitionSeconds() : 0.0;
}

double LocalReducedSchurSolver::portCoreCouplingAssemblySeconds() const
{
    return portCore_ ? portCore_->couplingAssemblySeconds() : 0.0;
}

double LocalReducedSchurSolver::portCoreLeafCsrSeconds() const
{
    return portCore_ ? portCore_->leafCsrSeconds() : 0.0;
}

double LocalReducedSchurSolver::portCoreLeafFactorSeconds() const
{
    return portCore_ ? portCore_->leafFactorSeconds() : 0.0;
}

double LocalReducedSchurSolver::portCoreEliminationSeconds() const
{
    return portCore_ ? portCore_->eliminationSeconds() : 0.0;
}

double LocalReducedSchurSolver::portCoreMultiRhsPortSeconds() const
{
    return portCore_ ? portCore_->multiRhsPortSeconds() : 0.0;
}

double LocalReducedSchurSolver::portCoreSchurProductPortSeconds() const
{
    return portCore_ ? portCore_->schurProductPortSeconds() : 0.0;
}

double LocalReducedSchurSolver::portCoreCoreAccumulationSeconds() const
{
    return portCore_ ? portCore_->coreAccumulationSeconds() : 0.0;
}

double LocalReducedSchurSolver::portCoreCoreCsrSeconds() const
{
    return portCore_ ? portCore_->coreCsrSeconds() : 0.0;
}

double LocalReducedSchurSolver::portCoreCoreFactorSeconds() const
{
    return portCore_ ? portCore_->coreFactorSeconds() : 0.0;
}

SolveResult LocalReducedSchurSolver::solve(
    const std::vector<double>& globalRhs,
    const std::vector<double>* interfaceInitialGuess,
    double interfaceToleranceOverride)
{
    if (globalRhs.size() != static_cast<std::size_t>(model_.globalDofs)) {
        throw std::runtime_error("[Local ROM] Global RHS has the wrong size.");
    }
    const int gammaSize = model_.interfaceDofs;
    std::vector<double> interfaceRhs(static_cast<std::size_t>(gammaSize), 0.0);
    for (int gamma = 0; gamma < gammaSize; ++gamma) {
        interfaceRhs[static_cast<std::size_t>(gamma)] = globalRhs[static_cast<std::size_t>(
            model_.interfaceGlobalDofs[static_cast<std::size_t>(gamma)])];
    }

    std::vector<std::vector<double>> projectedInteriorRhs(model_.subdomains.size());
    const auto assemblyStart = std::chrono::steady_clock::now();
    for (std::size_t slot = 0; slot < model_.subdomains.size(); ++slot) {
        const SubdomainModel& local = model_.subdomains[slot];
        const SubdomainModel& data = templatePayload(model_, local);
        std::vector<double>& reducedRhs = projectedInteriorRhs[slot];
        reducedRhs.assign(static_cast<std::size_t>(local.rank), 0.0);
        for (int mode = 0; mode < local.rank; ++mode) {
            for (int row = 0; row < local.interiorDofs; ++row) {
                const std::vector<double>& test = data.testBasis.empty()
                    ? data.basis : data.testBasis;
                reducedRhs[static_cast<std::size_t>(mode)] +=
                    test[static_cast<std::size_t>(mode * local.interiorDofs + row)]
                    * globalRhs[static_cast<std::size_t>(
                        local.interiorGlobalDofs[static_cast<std::size_t>(row)])];
            }
        }
    }
    const double projectionSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - assemblyStart).count();
    SolveResult result = solveWithReducedRhs(
        std::move(projectedInteriorRhs), std::move(interfaceRhs),
        interfaceInitialGuess, interfaceToleranceOverride);
    result.timing.localReducedAssemblySeconds += projectionSeconds;
    return result;
}

SolveResult LocalReducedSchurSolver::solveReducedRhs(
    const std::vector<std::vector<double>>& projectedInteriorRhs,
    const std::vector<double>& interfaceRhs,
    const std::vector<double>* interfaceInitialGuess,
    double interfaceToleranceOverride)
{
    return solveWithReducedRhs(
        projectedInteriorRhs, interfaceRhs, interfaceInitialGuess,
        interfaceToleranceOverride);
}

SolveResult LocalReducedSchurSolver::solveWithReducedRhs(
    std::vector<std::vector<double>> projectedInteriorRhs,
    std::vector<double> interfaceRhs,
    const std::vector<double>* interfaceInitialGuess,
    double interfaceToleranceOverride)
{
    if (interfaceRhs.size() != static_cast<std::size_t>(model_.interfaceDofs)
        || projectedInteriorRhs.size() != model_.subdomains.size()) {
        throw std::runtime_error("[Local ROM] Reduced RHS has the wrong dimensions.");
    }
    for (std::size_t slot = 0; slot < model_.subdomains.size(); ++slot) {
        if (projectedInteriorRhs[slot].size()
            != static_cast<std::size_t>(model_.subdomains[slot].rank)) {
            throw std::runtime_error(
                "[Local ROM] Reduced interior RHS has the wrong rank.");
        }
    }
    const auto totalStart = std::chrono::steady_clock::now();
    const int gammaSize = model_.interfaceDofs;
    std::vector<double> condensed = std::move(interfaceRhs);
    const auto assemblyStart = std::chrono::steady_clock::now();
    for (std::size_t slot = 0; slot < model_.subdomains.size(); ++slot) {
        const SubdomainModel& local = model_.subdomains[slot];
        const SubdomainModel& data = templatePayload(model_, local);
        std::vector<double>& reducedRhs = projectedInteriorRhs[slot];
        const bool hasInteriorReferenceImage = std::any_of(
            data.interiorReferenceImage.begin(),
            data.interiorReferenceImage.end(),
            [](double value) { return value != 0.0; });
        if (hasInteriorReferenceImage) {
            for (int mode = 0; mode < local.rank; ++mode) {
                for (int row = 0; row < local.interiorDofs; ++row) {
                    const std::vector<double>& test = data.testBasis.empty()
                        ? data.basis : data.testBasis;
                    reducedRhs[static_cast<std::size_t>(mode)] -=
                        test[static_cast<std::size_t>(
                            mode * local.interiorDofs + row)]
                        * data.interiorReferenceImage[static_cast<std::size_t>(row)];
                }
            }
        }
        for (std::size_t localGamma = 0;
             localGamma < local.interfaceIndices.size(); ++localGamma) {
            const int gamma = local.interfaceIndices[localGamma];
            condensed[static_cast<std::size_t>(gamma)] -=
                data.interfaceReferenceImage[localGamma];
        }
        if (!portCore_ && !augmentedDirect_) {
            std::vector<double> eliminated = reducedRhs;
            solveDenseSymmetric(localFactors_[slot], eliminated);
            for (std::size_t localGamma = 0;
                 localGamma < local.interfaceIndices.size(); ++localGamma) {
                const int gamma = local.interfaceIndices[localGamma];
                for (int mode = 0; mode < local.rank; ++mode) {
                    condensed[static_cast<std::size_t>(gamma)] -=
                        data.reducedInterfaceInterior[static_cast<std::size_t>(
                            localGamma * static_cast<std::size_t>(local.rank)
                            + mode)]
                        * eliminated[static_cast<std::size_t>(mode)];
                }
            }
        }
    }
    SolveResult result;
    result.timing.localReducedAssemblySeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - assemblyStart).count();

    const auto interfaceStart = std::chrono::steady_clock::now();
    std::vector<std::vector<double>> coupledDirectCoordinates;
    if (portCore_) {
        PortCoreSolveResult direct = portCore_->solve(
            projectedInteriorRhs, condensed);
        result.interfaceTemperature = std::move(direct.interfaceTemperature);
        coupledDirectCoordinates = std::move(direct.localReducedCoordinates);
        result.timing.interfaceIterations = 1;
        result.timing.interfaceMatvecs = 0;
        result.timing.interfaceInitialRelativeResidual = 1.0;
        result.timing.interfaceRelativeResidual =
            direct.reducedRelativeResidual;
        result.timing.interfaceKrylovActual = "port-core-direct";
        result.timing.portForwardSolveSeconds = direct.portForwardSeconds;
        result.timing.portCoreSolveSeconds = direct.coreSolveSeconds;
        result.timing.portBackSubstitutionSeconds =
            direct.portBackSubstitutionSeconds;
        result.status = std::isfinite(direct.reducedRelativeResidual)
                && direct.reducedRelativeResidual <= 1.0e-8
            ? "success" : "port_core_residual_failed";
    } else if (augmentedDirect_) {
        // The matrix is fixed for constant dt, so online work is only RHS
        // assembly plus a solve with the factor built in initializeAugmentedDirect().
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
        std::vector<double> augmentedSolution;
        {
            ScopedDirectSolverMklThreads mklThreads(
                augmentedDirect_->solveThreads);
            if (augmentedDirect_->generalFactor) {
                augmentedDirect_->generalFactor->solve(
                    augmentedRhs, augmentedSolution);
            } else {
                augmentedDirect_->factor->solve(
                    augmentedRhs, augmentedSolution);
            }
        }
        result.interfaceTemperature.assign(
            augmentedSolution.begin(),
            augmentedSolution.begin() + model_.interfaceDofs);
        coupledDirectCoordinates.resize(model_.subdomains.size());
        for (std::size_t slot = 0; slot < model_.subdomains.size(); ++slot) {
            const int begin = model_.interfaceDofs
                + augmentedDirect_->rankOffsets[slot];
            const int end = model_.interfaceDofs
                + augmentedDirect_->rankOffsets[slot + 1];
            coupledDirectCoordinates[slot].assign(
                augmentedSolution.begin() + begin,
                augmentedSolution.begin() + end);
        }

        // Audit the original unsymmetrized reduced equations periodically; this
        // catches assembly drift without adding a full matvec to every step.
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
                    coupledDirectCoordinates[slot];
                for (std::size_t localGamma = 0;
                     localGamma < local.interfaceIndices.size(); ++localGamma) {
                    const int gamma = local.interfaceIndices[localGamma];
                    for (int mode = 0; mode < local.rank; ++mode) {
                        interfaceResidual[static_cast<std::size_t>(gamma)] -=
                            data.reducedInterfaceInterior[
                                localGamma
                                    * static_cast<std::size_t>(local.rank)
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
        for (double value : augmentedSolution) {
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
    } else if (matrixFree_) {
        ddm_schur::Options krylovOptions = matrixFree_->options;
        if (interfaceToleranceOverride > 0.0) {
            krylovOptions.relativeTolerance = interfaceToleranceOverride;
        }
        const int matvecsBefore = matrixFree_->matvecs;
        const double proxySecondsBefore = matrixFree_->proxy
            ? matrixFree_->proxy->solveSeconds() : 0.0;
        const double coarseSecondsBefore = matrixFree_->coarseSolveSeconds;
        const auto applyOperator =
            [&](const std::vector<double>& input, std::vector<double>& output) {
                applyMatrixFree(input, output, true);
            };
        const auto applyPreconditioner =
            [&](const std::vector<double>& input, std::vector<double>& output) {
                applyMatrixFreePreconditioner(input, output);
            };
        const std::vector<double>* cachedInitialProduct = nullptr;
        if (interfaceInitialGuess
            && interfaceInitialGuess->size()
                == matrixFree_->previousInterfaceSolution.size()
            && matrixFree_->fgmresWorkspace.product.size()
                == interfaceInitialGuess->size()
            && std::equal(
                interfaceInitialGuess->begin(), interfaceInitialGuess->end(),
                matrixFree_->previousInterfaceSolution.begin())) {
            cachedInitialProduct = &matrixFree_->fgmresWorkspace.product;
        }
        const std::vector<double>* krylovInitialSolution =
            interfaceInitialGuess;
        const std::vector<double>* krylovInitialProduct =
            cachedInitialProduct;
        if (matrixFree_->options.interfaceKrylov == "fgmres"
            && matrixFree_->options.interfaceOperatorCoarsePredictor
            && matrixFree_->operatorCoarseDimension > 0
            && !matrixFree_->coarseBasis.empty()) {
            const auto predictorStart = std::chrono::steady_clock::now();
            std::vector<double>& predictorSolution =
                matrixFree_->predictorSolution;
            std::vector<double>& predictorProduct =
                matrixFree_->predictorProduct;
            std::vector<double>& predictorResidual =
                matrixFree_->predictorResidual;
            predictorSolution.assign(condensed.size(), 0.0);
            if (interfaceInitialGuess) {
                predictorSolution = *interfaceInitialGuess;
            }
            if (cachedInitialProduct) {
                predictorProduct = *cachedInitialProduct;
            } else if (interfaceInitialGuess) {
                applyMatrixFree(predictorSolution, predictorProduct, true);
            } else {
                predictorProduct.assign(condensed.size(), 0.0);
            }
            predictorResidual.resize(condensed.size());
            for (std::size_t i = 0; i < condensed.size(); ++i) {
                predictorResidual[i] = condensed[i] - predictorProduct[i];
            }
            const double rhsNorm = std::sqrt(std::max(
                0.0, vectorDot(condensed, condensed)));
            const double initialResidualNorm = std::sqrt(std::max(
                0.0, vectorDot(predictorResidual, predictorResidual)));
            const double initialRelativeResidual = rhsNorm > 0.0
                ? initialResidualNorm / rhsNorm : 0.0;
            std::vector<double>& coarseRhs =
                matrixFree_->predictorCoarseRhs;
            coarseRhs.assign(matrixFree_->coarseBasis.size(), 0.0);
            for (std::size_t coarse = 0;
                 coarse < matrixFree_->coarseBasis.size(); ++coarse) {
                const MatrixFreeData::CoarseVector& image =
                    matrixFree_->coarseImages[coarse];
                coarseRhs[coarse] = image.dot(predictorResidual);
            }
            const auto coarseStart = std::chrono::steady_clock::now();
            solveDenseSymmetric(matrixFree_->predictorFactor, coarseRhs);
            matrixFree_->coarseSolveSeconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - coarseStart).count();
            for (std::size_t coarse = 0;
                 coarse < matrixFree_->coarseBasis.size(); ++coarse) {
                const double coefficient = coarseRhs[coarse];
                const MatrixFreeData::CoarseVector& basis =
                    matrixFree_->coarseBasis[coarse];
                basis.addScaled(coefficient, predictorSolution);
                const MatrixFreeData::CoarseVector& image =
                    matrixFree_->coarseImages[coarse];
                image.addScaled(coefficient, predictorProduct);
            }
            for (std::size_t i = 0; i < condensed.size(); ++i) {
                predictorResidual[i] = condensed[i] - predictorProduct[i];
            }
            const double predictedResidualNorm = std::sqrt(std::max(
                0.0, vectorDot(predictorResidual, predictorResidual)));
            const double predictedRelativeResidual = rhsNorm > 0.0
                ? predictedResidualNorm / rhsNorm : 0.0;
            const bool acceptPredictor = std::isfinite(predictedRelativeResidual)
                && predictedRelativeResidual < initialRelativeResidual;
            result.timing.interfacePredictorApplied = true;
            result.timing.interfacePredictorAccepted = acceptPredictor;
            result.timing.interfacePredictorInitialRelativeResidual =
                initialRelativeResidual;
            result.timing.interfacePredictorRelativeResidual =
                predictedRelativeResidual;
            if (acceptPredictor) {
                krylovInitialSolution = &predictorSolution;
                krylovInitialProduct = &predictorProduct;
            } else if (!cachedInitialProduct && interfaceInitialGuess) {
                for (std::size_t coarse = 0;
                     coarse < matrixFree_->coarseImages.size(); ++coarse) {
                    const double coefficient = coarseRhs[coarse];
                    const MatrixFreeData::CoarseVector& image =
                        matrixFree_->coarseImages[coarse];
                    image.addScaled(-coefficient, predictorProduct);
                }
                krylovInitialProduct = &predictorProduct;
            }
            result.timing.interfacePredictorSeconds =
                std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - predictorStart).count();
        }
        ddm_schur::MatrixFreeFgmresResult krylovResult;
        if (matrixFree_->options.interfaceKrylov == "pcg") {
            krylovResult = ddm_schur::solveMatrixFreePcg(
                condensed,
                krylovOptions,
                applyOperator,
                applyPreconditioner,
                interfaceInitialGuess);
            if (!krylovResult.converged) {
                const int pcgIterations = krylovResult.iterations;
                const double initialResidual =
                    krylovResult.initialRelativeResidual;
                const std::string fallbackReason =
                    krylovResult.fallbackReason;
                ddm_schur::MatrixFreeFgmresResult fallback =
                    ddm_schur::solveMatrixFreeFgmres(
                        condensed,
                        krylovOptions,
                        applyOperator,
                        applyPreconditioner,
                        &krylovResult.solution,
                        &matrixFree_->fgmresWorkspace);
                fallback.iterations += pcgIterations;
                fallback.initialRelativeResidual = initialResidual;
                fallback.actualSolver = "pcg->fgmres";
                fallback.fallbackTriggered = true;
                fallback.fallbackReason = fallbackReason;
                fallback.operatorSeconds += krylovResult.operatorSeconds;
                fallback.preconditionerSeconds +=
                    krylovResult.preconditionerSeconds;
                fallback.orthogonalizationSeconds +=
                    krylovResult.orthogonalizationSeconds;
                fallback.vectorUpdateSeconds +=
                    krylovResult.vectorUpdateSeconds;
                krylovResult = std::move(fallback);
            }
        } else {
            krylovResult = ddm_schur::solveMatrixFreeFgmres(
                condensed,
                krylovOptions,
                applyOperator,
                applyPreconditioner,
                krylovInitialSolution,
                &matrixFree_->fgmresWorkspace,
                krylovInitialProduct);
        }
        const bool cacheInterfaceProduct = krylovResult.converged
            && krylovResult.solutionProductAvailable
            && (krylovResult.actualSolver == "fgmres"
                || krylovResult.actualSolver == "pcg->fgmres")
            && matrixFree_->fgmresWorkspace.product.size()
                == krylovResult.solution.size();
        result.interfaceTemperature = std::move(krylovResult.solution);
        if (cacheInterfaceProduct) {
            matrixFree_->previousInterfaceSolution =
                result.interfaceTemperature;
        }
        result.timing.interfaceIterations = krylovResult.iterations;
        result.timing.interfaceInitialRelativeResidual =
            krylovResult.initialRelativeResidual;
        result.timing.interfaceRelativeResidual = krylovResult.relativeResidual;
        result.timing.interfaceKrylovActual = krylovResult.actualSolver;
        result.timing.interfaceKrylovFallback =
            krylovResult.fallbackTriggered;
        result.timing.interfaceKrylovFallbackReason =
            krylovResult.fallbackReason;
        result.timing.interfaceMatvecs = matrixFree_->matvecs - matvecsBefore;
        result.timing.proxySolveSeconds = matrixFree_->proxy
            ? matrixFree_->proxy->solveSeconds() - proxySecondsBefore : 0.0;
        result.timing.coarseSolveSeconds =
            matrixFree_->coarseSolveSeconds - coarseSecondsBefore;
        result.timing.interfaceOperatorSeconds =
            krylovResult.operatorSeconds;
        result.timing.interfacePreconditionerSeconds =
            krylovResult.preconditionerSeconds;
        result.timing.interfaceOrthogonalizationSeconds =
            krylovResult.orthogonalizationSeconds;
        result.timing.interfaceVectorUpdateSeconds =
            krylovResult.vectorUpdateSeconds;
        result.status = krylovResult.converged
            ? "success" : "interface_krylov_not_converged";
    } else if (sparseFactor_) {
        sparseFactor_->solver.solve(condensed, result.interfaceTemperature);
    } else {
        result.interfaceTemperature = condensed;
        solveDenseSymmetric(schurFactor_, result.interfaceTemperature);
    }
    result.timing.interfaceSolveSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - interfaceStart).count();

    const auto recoveryStart = std::chrono::steady_clock::now();
    const int domainCount = static_cast<int>(model_.subdomains.size());
    const bool parallelRecovery = matrixFree_
        && matrixFree_->options.localSolveThreads > 1;
    std::vector<std::vector<double>> explicitRecoveryCoordinates;
    if (portCore_ || augmentedDirect_) {
        explicitRecoveryCoordinates = std::move(coupledDirectCoordinates);
    } else if (!matrixFree_) {
        explicitRecoveryCoordinates.resize(model_.subdomains.size());
    }
    if (!portCore_ && !augmentedDirect_) {
#pragma omp parallel for num_threads(matrixFree_ ? matrixFree_->options.localSolveThreads : 1) \
    if(parallelRecovery) schedule(static)
        for (int slot = 0; slot < domainCount; ++slot) {
            const SubdomainModel& local =
                model_.subdomains[static_cast<std::size_t>(slot)];
            const SubdomainModel& data = templatePayload(model_, local);
            std::vector<double>& coordinates = matrixFree_
                ? matrixFree_->localApplyWorkspaces[static_cast<std::size_t>(slot)]
                      .recoveryCoordinates
                : explicitRecoveryCoordinates[static_cast<std::size_t>(slot)];
            coordinates = projectedInteriorRhs[static_cast<std::size_t>(slot)];
            for (int mode = 0; mode < local.rank; ++mode) {
                for (std::size_t localGamma = 0;
                     localGamma < local.interfaceIndices.size(); ++localGamma) {
                    const int gamma = local.interfaceIndices[localGamma];
                    coordinates[static_cast<std::size_t>(mode)] -=
                        data.reducedInteriorInterface[static_cast<std::size_t>(
                            mode * local.localInterfaceDofs)
                            + localGamma]
                        * result.interfaceTemperature[
                            static_cast<std::size_t>(gamma)];
                }
            }
            solveDenseSymmetric(
                localFactors_[static_cast<std::size_t>(slot)], coordinates);
        }
    }
    result.timing.localRecoverySeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - recoveryStart).count();
    result.localReducedCoordinates.resize(static_cast<std::size_t>(domainCount));
    for (int slot = 0; slot < domainCount; ++slot) {
        result.localReducedCoordinates[static_cast<std::size_t>(slot)] =
            matrixFree_
                ? matrixFree_->localApplyWorkspaces[static_cast<std::size_t>(slot)]
                      .recoveryCoordinates
                : explicitRecoveryCoordinates[static_cast<std::size_t>(slot)];
    }

    const auto reconstructionStart = std::chrono::steady_clock::now();
    result.temperature.assign(static_cast<std::size_t>(model_.globalDofs), 0.0);
    for (int gamma = 0; gamma < gammaSize; ++gamma) {
        result.temperature[static_cast<std::size_t>(
            model_.interfaceGlobalDofs[static_cast<std::size_t>(gamma)])] =
            result.interfaceTemperature[static_cast<std::size_t>(gamma)];
    }
#pragma omp parallel for num_threads(matrixFree_ ? matrixFree_->options.localSolveThreads : 1) \
    if(parallelRecovery) schedule(static)
    for (int slot = 0; slot < domainCount; ++slot) {
        const SubdomainModel& local = model_.subdomains[static_cast<std::size_t>(slot)];
        const SubdomainModel& data = templatePayload(model_, local);
        const std::vector<double>& coordinates = matrixFree_
            ? matrixFree_->localApplyWorkspaces[static_cast<std::size_t>(slot)]
                  .recoveryCoordinates
            : explicitRecoveryCoordinates[static_cast<std::size_t>(slot)];
        // Preserve the per-row accumulation order used by exact steady-state reuse.
        for (int row = 0; row < local.interiorDofs; ++row) {
            double value = data.referenceInterior[static_cast<std::size_t>(row)];
            for (int mode = 0; mode < local.rank; ++mode) {
                value += data.basis[static_cast<std::size_t>(
                    mode * local.interiorDofs + row)]
                    * coordinates[static_cast<std::size_t>(mode)];
            }
            result.temperature[static_cast<std::size_t>(
                local.interiorGlobalDofs[static_cast<std::size_t>(row)])]
                = value;
        }
    }
    result.timing.fullFieldReconstructionSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - reconstructionStart).count();
    result.timing.totalSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - totalStart).count();
    if (result.status == "not_run") {
        result.status = "success";
    }
    return result;
}

} // namespace mor::local
