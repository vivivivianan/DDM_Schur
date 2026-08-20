// Operator-derived local Block-Arnoldi construction. The implementation solves
// with retained local interior factors, performs deterministic block
// orthogonalization/deflation, and truncates by numerical rank. Production uses
// one moment block (M1); no snapshot matrix or POD training data is consumed.

#include "block_arnoldi.hpp"

#include "linear_solvers.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

namespace mor::transient {
namespace {

using Clock = std::chrono::steady_clock;

double secondsSince(const Clock::time_point& start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

double dotColumns(const double* left, const double* right, int rows)
{
#ifdef USE_MKL_PARDISO
    return cblas_ddot(rows, left, 1, right, 1);
#else
    double result = 0.0;
    for (int row = 0; row < rows; ++row) result += left[row] * right[row];
    return result;
#endif
}

double normColumn(const double* column, int rows)
{
    return std::sqrt(std::max(0.0, dotColumns(column, column, rows)));
}

void axpyColumn(double coefficient,
                const double* source,
                double* target,
                int rows)
{
#ifdef USE_MKL_PARDISO
    cblas_daxpy(rows, coefficient, source, 1, target, 1);
#else
    for (int row = 0; row < rows; ++row) {
        target[row] += coefficient * source[row];
    }
#endif
}

void scaleColumn(double coefficient, double* column, int rows)
{
#ifdef USE_MKL_PARDISO
    cblas_dscal(rows, coefficient, column, 1);
#else
    for (int row = 0; row < rows; ++row) column[row] *= coefficient;
#endif
}

void swapColumns(std::vector<double>& matrix, int rows, int left, int right)
{
    if (left == right) return;
    double* a = matrix.data() + static_cast<std::size_t>(left) * rows;
    double* b = matrix.data() + static_cast<std::size_t>(right) * rows;
    for (int row = 0; row < rows; ++row) std::swap(a[row], b[row]);
}

void orthogonalizeAgainstBasis(std::vector<double>& block,
                               int rows,
                               int columns,
                               const std::vector<double>& basis,
                               int basisColumns)
{
    if (basisColumns == 0 || columns == 0) return;
#ifdef USE_MKL_PARDISO
    std::vector<double> coefficients(
        static_cast<std::size_t>(basisColumns * columns), 0.0);
    for (int pass = 0; pass < 2; ++pass) {
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
            basisColumns, columns, rows, 1.0,
            basis.data(), rows, block.data(), rows,
            0.0, coefficients.data(), basisColumns);
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
            rows, columns, basisColumns, -1.0,
            basis.data(), rows, coefficients.data(), basisColumns,
            1.0, block.data(), rows);
    }
#else
    for (int pass = 0; pass < 2; ++pass) {
        for (int column = 0; column < columns; ++column) {
            double* candidate = block.data() + static_cast<std::size_t>(column) * rows;
            for (int mode = 0; mode < basisColumns; ++mode) {
                const double* q = basis.data() + static_cast<std::size_t>(mode) * rows;
                axpyColumn(-dotColumns(q, candidate, rows), q, candidate, rows);
            }
        }
    }
#endif
}

struct OrthonormalBlock {
    std::vector<double> values;
    int columns = 0;
    int deflated = 0;
    double residual = 0.0;
    double orthogonalityError = 0.0;
};

OrthonormalBlock rankRevealingMgs(std::vector<double> candidates,
                                  int rows,
                                  int columns,
                                  const std::vector<double>& existing,
                                  int existingColumns,
                                  double tolerance)
{
    OrthonormalBlock result;
    double largestCandidateNorm = 0.0;
    for (int column = 0; column < columns; ++column) {
        largestCandidateNorm = std::max(largestCandidateNorm, normColumn(
            candidates.data() + static_cast<std::size_t>(column) * rows, rows));
    }
    orthogonalizeAgainstBasis(
        candidates, rows, columns, existing, existingColumns);
    double originalSquared = 0.0;
    double largestNorm = 0.0;
    for (int column = 0; column < columns; ++column) {
        const double norm = normColumn(
            candidates.data() + static_cast<std::size_t>(column) * rows, rows);
        originalSquared += norm * norm;
        largestNorm = std::max(largestNorm, norm);
    }
    // If the existing Krylov basis already represents this block, projection
    // can leave only roundoff. Scale by the unprojected block so that noise is
    // deflated instead of promoted to additional Arnoldi modes.
    const double thresholdScale = existingColumns > 0
        ? largestCandidateNorm : largestNorm;
    const double threshold = tolerance * std::max(1.0e-300, thresholdScale);
#ifdef USE_MKL_PARDISO
    std::vector<lapack_int> pivots(static_cast<std::size_t>(columns), 0);
    std::vector<double> tau(static_cast<std::size_t>(
        std::min(rows, columns)), 0.0);
    const lapack_int factorInfo = LAPACKE_dgeqp3(
        LAPACK_COL_MAJOR, static_cast<lapack_int>(rows),
        static_cast<lapack_int>(columns), candidates.data(),
        static_cast<lapack_int>(rows), pivots.data(), tau.data());
    if (factorInfo != 0) {
        throw std::runtime_error(
            "Block Arnoldi rank-revealing QR factorization failed.");
    }
    const int maximumRank = std::min(rows, columns);
    while (result.columns < maximumRank) {
        const double diagonal = std::abs(candidates[static_cast<std::size_t>(
            result.columns * rows + result.columns)]);
        if (!(diagonal > threshold) || !std::isfinite(diagonal)) break;
        ++result.columns;
    }
    result.deflated = columns - result.columns;
    double remainingSquared = 0.0;
    for (int column = result.columns; column < columns; ++column) {
        const int lastRow = std::min(column, rows - 1);
        for (int row = result.columns; row <= lastRow; ++row) {
            const double value = candidates[static_cast<std::size_t>(
                column * rows + row)];
            remainingSquared += value * value;
        }
    }
    result.residual = std::sqrt(remainingSquared)
        / std::max(1.0e-300, std::sqrt(originalSquared));
    if (result.columns > 0) {
        const lapack_int basisInfo = LAPACKE_dorgqr(
            LAPACK_COL_MAJOR, static_cast<lapack_int>(rows),
            static_cast<lapack_int>(result.columns),
            static_cast<lapack_int>(result.columns), candidates.data(),
            static_cast<lapack_int>(rows), tau.data());
        if (basisInfo != 0) {
            throw std::runtime_error(
                "Block Arnoldi orthonormal basis generation failed.");
        }
    }
    candidates.resize(static_cast<std::size_t>(rows)
        * static_cast<std::size_t>(result.columns));
    result.values = std::move(candidates);

    double maximumDot = 0.0;
    if (result.columns > 0) {
        std::vector<double> gram(static_cast<std::size_t>(
            result.columns * result.columns), 0.0);
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
            result.columns, result.columns, rows, 1.0,
            result.values.data(), rows, result.values.data(), rows,
            0.0, gram.data(), result.columns);
        for (int column = 0; column < result.columns; ++column) {
            for (int row = 0; row < result.columns; ++row) {
                maximumDot = std::max(maximumDot, std::abs(
                    gram[static_cast<std::size_t>(column * result.columns + row)]
                    - (row == column ? 1.0 : 0.0)));
            }
        }
        if (existingColumns > 0) {
            std::vector<double> cross(static_cast<std::size_t>(
                existingColumns * result.columns), 0.0);
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                existingColumns, result.columns, rows, 1.0,
                existing.data(), rows, result.values.data(), rows,
                0.0, cross.data(), existingColumns);
            for (double value : cross) {
                maximumDot = std::max(maximumDot, std::abs(value));
            }
        }
    }
    result.orthogonalityError = maximumDot;
    return result;
#else
    for (int selected = 0; selected < columns; ++selected) {
        int pivot = selected;
        double pivotNorm = -1.0;
        for (int column = selected; column < columns; ++column) {
            const double norm = normColumn(candidates.data()
                + static_cast<std::size_t>(column) * rows, rows);
            if (norm > pivotNorm) {
                pivotNorm = norm;
                pivot = column;
            }
        }
        if (!(pivotNorm > threshold) || !std::isfinite(pivotNorm)) break;
        swapColumns(candidates, rows, selected, pivot);
        double* q = candidates.data() + static_cast<std::size_t>(selected) * rows;
        scaleColumn(1.0 / pivotNorm, q, rows);
        ++result.columns;
        for (int column = selected + 1; column < columns; ++column) {
            double* candidate = candidates.data()
                + static_cast<std::size_t>(column) * rows;
            for (int pass = 0; pass < 2; ++pass) {
                axpyColumn(-dotColumns(q, candidate, rows), q, candidate, rows);
            }
        }
    }
    result.deflated = columns - result.columns;
    double remainingSquared = 0.0;
    for (int column = result.columns; column < columns; ++column) {
        const double norm = normColumn(candidates.data()
            + static_cast<std::size_t>(column) * rows, rows);
        remainingSquared += norm * norm;
    }
    result.residual = std::sqrt(remainingSquared)
        / std::max(1.0e-300, std::sqrt(originalSquared));

    // Accepted columns have already been pivoted and normalized into the
    // leading part of candidates. Reuse that multi-RHS allocation instead of
    // holding a second full n-by-block copy, which is material for RRAM26.
    candidates.resize(static_cast<std::size_t>(rows)
        * static_cast<std::size_t>(result.columns));
    result.values = std::move(candidates);

    double maximumDot = 0.0;
    for (int column = 0; column < result.columns; ++column) {
        const double* q = result.values.data()
            + static_cast<std::size_t>(column) * rows;
        maximumDot = std::max(maximumDot,
            std::abs(dotColumns(q, q, rows) - 1.0));
        for (int previous = 0; previous < column; ++previous) {
            maximumDot = std::max(maximumDot, std::abs(dotColumns(
                result.values.data() + static_cast<std::size_t>(previous) * rows,
                q, rows)));
        }
        for (int previous = 0; previous < existingColumns; ++previous) {
            maximumDot = std::max(maximumDot, std::abs(dotColumns(
                existing.data() + static_cast<std::size_t>(previous) * rows,
                q, rows)));
        }
    }
    result.orthogonalityError = maximumDot;
    return result;
#endif
}

void capacityBlockProduct(const SparseMatrix& capacity,
                          const double* basisBlock,
                          int blockColumns,
                          std::vector<double>& product)
{
    const int rows = capacity.size();
    product.assign(static_cast<std::size_t>(rows * blockColumns), 0.0);
    parallelFor(static_cast<std::size_t>(rows), [&](std::size_t row) {
        for (int offset = capacity.rowPtr[row];
             offset < capacity.rowPtr[row + 1]; ++offset) {
            const int column = capacity.colInd[static_cast<std::size_t>(offset)];
            const double value = capacity.values[static_cast<std::size_t>(offset)];
            for (int block = 0; block < blockColumns; ++block) {
                product[static_cast<std::size_t>(block) * rows + row] -= value
                    * basisBlock[static_cast<std::size_t>(block) * rows + column];
            }
        }
    });
}

struct SecondMomentSelection {
    int candidateColumns = 0;
    int selectedColumns = 0;
    double retainedEnergy = 1.0;
};

SecondMomentSelection selectSecondMomentColumns(
    std::vector<double>& rhs,
    int rows,
    int columns,
    double energyTarget,
    int maximumColumns)
{
    SecondMomentSelection selection;
    selection.candidateColumns = columns;
    selection.selectedColumns = columns;
    if (columns == 0 || (energyTarget == 1.0 && maximumColumns == 0)) {
        return selection;
    }

    struct ColumnEnergy {
        double squared = 0.0;
        int column = 0;
    };
    std::vector<ColumnEnergy> energies(static_cast<std::size_t>(columns));
    double totalSquared = 0.0;
    for (int column = 0; column < columns; ++column) {
        const double squared = dotColumns(rhs.data()
                + static_cast<std::size_t>(column) * rows,
            rhs.data() + static_cast<std::size_t>(column) * rows, rows);
        energies[static_cast<std::size_t>(column)] = {
            std::max(0.0, squared), column};
        totalSquared += std::max(0.0, squared);
    }
    std::stable_sort(energies.begin(), energies.end(),
        [](const ColumnEnergy& left, const ColumnEnergy& right) {
            return left.squared > right.squared;
        });

    const int limit = maximumColumns > 0
        ? std::min(maximumColumns, columns) : columns;
    const double requiredSquared = energyTarget * totalSquared;
    double retainedSquared = 0.0;
    int selected = 0;
    while (selected < limit && (selected == 0 || retainedSquared < requiredSquared)) {
        retainedSquared += energies[static_cast<std::size_t>(selected)].squared;
        ++selected;
    }
    selection.selectedColumns = selected;
    selection.retainedEnergy = totalSquared > 0.0
        ? std::min(1.0, retainedSquared / totalSquared) : 1.0;

    std::vector<double> selectedRhs(static_cast<std::size_t>(rows) * selected);
    for (int selectedColumn = 0; selectedColumn < selected; ++selectedColumn) {
        const int sourceColumn = energies[static_cast<std::size_t>(selectedColumn)].column;
        std::copy_n(rhs.data() + static_cast<std::size_t>(sourceColumn) * rows,
            rows, selectedRhs.data() + static_cast<std::size_t>(selectedColumn) * rows);
    }
    rhs = std::move(selectedRhs);
    return selection;
}

} // namespace

BlockArnoldiResult buildBlockArnoldiBasis(
    const ThermalDescriptorSystem& system,
    int moments,
    double expansionPoint,
    double rankTolerance,
    double secondMomentEnergy,
    int secondMomentMaximumColumns,
    int directSolverThreads)
{
    const auto totalStart = Clock::now();
    if (moments <= 0 || system.sourceChannels <= 0) {
        throw std::runtime_error(
            "Block Arnoldi requires positive moments and source channels.");
    }
    if (expansionPoint != 0.0) {
        throw std::runtime_error(
            "Stage 2C.1 implements the validated s0=0 expansion only.");
    }
    if (!(rankTolerance > 0.0) || !std::isfinite(rankTolerance)) {
        throw std::runtime_error("Block Arnoldi rank tolerance must be positive.");
    }
    if (!(secondMomentEnergy > 0.0) || secondMomentEnergy > 1.0
        || !std::isfinite(secondMomentEnergy)) {
        throw std::runtime_error(
            "Block Arnoldi second-moment energy must lie in (0, 1].");
    }
    if (secondMomentMaximumColumns < 0) {
        throw std::runtime_error(
            "Block Arnoldi second-moment maximum columns must be nonnegative.");
    }
    if (directSolverThreads <= 0) {
        throw std::runtime_error(
            "Block Arnoldi direct-solver threads must be positive.");
    }
    BlockArnoldiResult result;
    result.rows = system.dofs;
    result.blockSize = system.sourceChannels;
    result.moments = moments;
    result.expansionPoint = expansionPoint;
    result.rankTolerance = rankTolerance;

    // Factor K_II once and reuse it for every block moment. The caller controls
    // inner MKL threads because several subdomain Arnoldi builds may run at once.
    ScopedDirectSolverMklThreads directThreads(directSolverThreads);
    const std::vector<MatrixEntry> entries =
        sparseMatrixEntries(system.conductivity);
    SubdomainDirectSolver factor(system.dofs, entries);
    result.timing.symbolicAnalysisSeconds = factor.symbolicAnalysisSeconds();
    result.timing.numericalFactorizationSeconds =
        factor.numericalFactorizationSeconds();
    result.timing.symbolicAnalysisCalls = factor.symbolicAnalysisCalls();
    result.timing.numericalFactorizationCalls = factor.numericalFactorizationCalls();

    auto orthStart = Clock::now();
    const std::vector<double> noExistingBasis;
    OrthonormalBlock compressedInput = rankRevealingMgs(
        system.input, system.dofs, system.sourceChannels,
        noExistingBasis, 0, rankTolerance * 1.0e-2);
    const double inputCompressionSeconds = secondsSince(orthStart);
    result.timing.orthogonalizationSeconds += inputCompressionSeconds;
    if (compressedInput.columns <= 0) {
        throw std::runtime_error(
            "Block Arnoldi input block is numerically rank deficient.");
    }

    std::vector<double> initialResponse;
    auto solveStart = Clock::now();
    factor.solveMultiple(
        compressedInput.values, compressedInput.columns, initialResponse,
        directSolverThreads);
    double solveSeconds = secondsSince(solveStart);
    result.timing.multiRhsSolveSeconds += solveSeconds;
    std::vector<double> reference;
    solveStart = Clock::now();
    factor.solve(system.boundaryRhs, reference);
    result.timing.multiRhsSolveSeconds += secondsSince(solveStart);
    result.referenceTemperature = std::move(reference);

    orthStart = Clock::now();
    OrthonormalBlock block = rankRevealingMgs(
        std::move(initialResponse), system.dofs, compressedInput.columns,
        result.basis, 0, rankTolerance);
    double orthSeconds = secondsSince(orthStart);
    result.timing.orthogonalizationSeconds += orthSeconds;
    result.basis = std::move(block.values);
    result.rank = block.columns;
    const int initialDeflated = system.sourceChannels
        - compressedInput.columns + block.deflated;
    result.history.push_back({
        1, system.sourceChannels, block.columns, result.rank, initialDeflated,
        block.orthogonalityError, block.residual, solveSeconds,
        inputCompressionSeconds + orthSeconds,
        result.basis.size() * sizeof(double)});
    int previousBlockStart = 0;
    int previousBlockColumns = block.columns;

    for (int moment = 2; moment <= moments && previousBlockColumns > 0; ++moment) {
        std::vector<double> rhs;
        capacityBlockProduct(system.capacity,
            result.basis.data() + static_cast<std::size_t>(previousBlockStart)
                * system.dofs,
            previousBlockColumns, rhs);
        SecondMomentSelection selection;
        if (moment == 2) {
            selection = selectSecondMomentColumns(rhs, system.dofs,
                previousBlockColumns, secondMomentEnergy,
                secondMomentMaximumColumns);
        } else {
            selection.candidateColumns = previousBlockColumns;
            selection.selectedColumns = previousBlockColumns;
        }
        std::vector<double> response;
        solveStart = Clock::now();
        factor.solveMultiple(
            rhs, selection.selectedColumns, response, directSolverThreads);
        solveSeconds = secondsSince(solveStart);
        result.timing.multiRhsSolveSeconds += solveSeconds;
        orthStart = Clock::now();
        block = rankRevealingMgs(
            std::move(response), system.dofs, selection.selectedColumns,
            result.basis, result.rank, rankTolerance);
        orthSeconds = secondsSince(orthStart);
        result.timing.orthogonalizationSeconds += orthSeconds;
        previousBlockStart = result.rank;
        previousBlockColumns = block.columns;
        result.basis.reserve(result.basis.size() + block.values.size());
        result.basis.insert(result.basis.end(),
            std::make_move_iterator(block.values.begin()),
            std::make_move_iterator(block.values.end()));
        result.rank += block.columns;
        result.history.push_back({
            moment, selection.selectedColumns, block.columns, result.rank,
            selection.candidateColumns - selection.selectedColumns + block.deflated,
            block.orthogonalityError, block.residual, solveSeconds, orthSeconds,
            result.basis.size() * sizeof(double)});
    }
    result.timing.totalSeconds = secondsSince(totalStart);
    return result;
}

} // namespace mor::transient
