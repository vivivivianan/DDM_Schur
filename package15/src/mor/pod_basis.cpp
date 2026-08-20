#include "pod_basis.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

namespace mor {
namespace {

struct EigenSystem {
    std::vector<double> values;
    std::vector<double> vectors;
};

EigenSystem jacobiSymmetric(std::vector<double> a, int n)
{
    EigenSystem result;
    result.vectors.assign(static_cast<std::size_t>(n * n), 0.0);
    for (int i = 0; i < n; ++i) {
        result.vectors[static_cast<std::size_t>(i * n + i)] = 1.0;
    }
    if (n == 0) {
        return result;
    }
    const int maximumSweeps = std::max(30, 8 * n * n);
    for (int sweep = 0; sweep < maximumSweeps; ++sweep) {
        int p = 0;
        int q = 0;
        double largest = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const double value = std::abs(a[static_cast<std::size_t>(i * n + j)]);
                if (value > largest) {
                    largest = value;
                    p = i;
                    q = j;
                }
            }
        }
        double diagonalScale = 0.0;
        for (int i = 0; i < n; ++i) {
            diagonalScale = std::max(diagonalScale,
                std::abs(a[static_cast<std::size_t>(i * n + i)]));
        }
        if (largest <= 32.0 * std::numeric_limits<double>::epsilon()
            * std::max(1.0, diagonalScale)) {
            break;
        }
        const double app = a[static_cast<std::size_t>(p * n + p)];
        const double aqq = a[static_cast<std::size_t>(q * n + q)];
        const double apq = a[static_cast<std::size_t>(p * n + q)];
        const double angle = 0.5 * std::atan2(2.0 * apq, aqq - app);
        const double c = std::cos(angle);
        const double s = std::sin(angle);
        for (int k = 0; k < n; ++k) {
            if (k == p || k == q) {
                continue;
            }
            const double akp = a[static_cast<std::size_t>(k * n + p)];
            const double akq = a[static_cast<std::size_t>(k * n + q)];
            const double newKp = c * akp - s * akq;
            const double newKq = s * akp + c * akq;
            a[static_cast<std::size_t>(k * n + p)] = newKp;
            a[static_cast<std::size_t>(p * n + k)] = newKp;
            a[static_cast<std::size_t>(k * n + q)] = newKq;
            a[static_cast<std::size_t>(q * n + k)] = newKq;
        }
        a[static_cast<std::size_t>(p * n + p)] =
            c * c * app - 2.0 * s * c * apq + s * s * aqq;
        a[static_cast<std::size_t>(q * n + q)] =
            s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[static_cast<std::size_t>(p * n + q)] = 0.0;
        a[static_cast<std::size_t>(q * n + p)] = 0.0;
        for (int k = 0; k < n; ++k) {
            const double vkp = result.vectors[static_cast<std::size_t>(k * n + p)];
            const double vkq = result.vectors[static_cast<std::size_t>(k * n + q)];
            result.vectors[static_cast<std::size_t>(k * n + p)] = c * vkp - s * vkq;
            result.vectors[static_cast<std::size_t>(k * n + q)] = s * vkp + c * vkq;
        }
    }
    result.values.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        result.values[static_cast<std::size_t>(i)] =
            a[static_cast<std::size_t>(i * n + i)];
    }
    std::vector<int> order(static_cast<std::size_t>(n));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int left, int right) {
        return result.values[static_cast<std::size_t>(left)]
            > result.values[static_cast<std::size_t>(right)];
    });
    EigenSystem sorted;
    sorted.values.resize(static_cast<std::size_t>(n));
    sorted.vectors.resize(static_cast<std::size_t>(n * n));
    for (int column = 0; column < n; ++column) {
        const int original = order[static_cast<std::size_t>(column)];
        sorted.values[static_cast<std::size_t>(column)] =
            result.values[static_cast<std::size_t>(original)];
        for (int row = 0; row < n; ++row) {
            sorted.vectors[static_cast<std::size_t>(row * n + column)] =
                result.vectors[static_cast<std::size_t>(row * n + original)];
        }
    }
    return sorted;
}

EigenSystem solveSymmetric(std::vector<double> matrix, int n)
{
#ifdef USE_MKL_PARDISO
    // A symmetric matrix has identical row-major and column-major storage.
    // LAPACKE returns eigenvectors as columns; copy them into the row-major
    // (row,mode) layout used by the POD assembly below and reverse the
    // ascending LAPACK eigenvalue order.
    std::vector<double> eigenvalues(static_cast<std::size_t>(n), 0.0);
    const lapack_int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U',
        static_cast<lapack_int>(n), matrix.data(), static_cast<lapack_int>(n),
        eigenvalues.data());
    if (info != 0) {
        throw std::runtime_error("[MOR] MKL symmetric Gram eigensolve failed with info="
            + std::to_string(info));
    }
    EigenSystem result;
    result.values.resize(static_cast<std::size_t>(n));
    result.vectors.resize(static_cast<std::size_t>(n * n));
    for (int mode = 0; mode < n; ++mode) {
        const int original = n - 1 - mode;
        result.values[static_cast<std::size_t>(mode)] =
            eigenvalues[static_cast<std::size_t>(original)];
        for (int row = 0; row < n; ++row) {
            result.vectors[static_cast<std::size_t>(row * n + mode)] =
                matrix[static_cast<std::size_t>(row + original * n)];
        }
    }
    return result;
#else
    return jacobiSymmetric(std::move(matrix), n);
#endif
}

} // namespace

std::vector<double> symmetricEigenvalues(std::vector<double> matrix, int n)
{
    return solveSymmetric(std::move(matrix), n).values;
}

PodResult buildGramPod(const SnapshotDatabase& snapshots,
                       int requestedRank,
                       double discardedEnergyTolerance,
                       double relativeSingularTolerance)
{
    if (snapshots.rows <= 0 || snapshots.cases.empty()) {
        throw std::runtime_error("[MOR] Cannot build POD from an empty snapshot database.");
    }
    const int rows = snapshots.rows;
    const int columns = static_cast<int>(snapshots.cases.size());
    if (static_cast<int>(snapshots.values.size()) != rows * columns) {
        throw std::runtime_error("[MOR] Snapshot database dimensions are inconsistent.");
    }
    PodResult result;
    result.rows = rows;
    result.columns = columns;

    const auto gramStart = std::chrono::steady_clock::now();
    std::vector<double> gram(static_cast<std::size_t>(columns * columns), 0.0);
#ifdef USE_MKL_PARDISO
    // X is stored column-major (rows x columns).  Form the snapshot Gram
    // matrix with one BLAS-3 operation instead of scalar triple loops.
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
        columns, columns, rows, 1.0,
        snapshots.values.data(), rows,
        snapshots.values.data(), rows,
        0.0, gram.data(), columns);
#else
    for (int left = 0; left < columns; ++left) {
        for (int right = left; right < columns; ++right) {
            double value = 0.0;
            const std::size_t leftOffset = static_cast<std::size_t>(left * rows);
            const std::size_t rightOffset = static_cast<std::size_t>(right * rows);
            for (int row = 0; row < rows; ++row) {
                value += snapshots.values[leftOffset + static_cast<std::size_t>(row)]
                    * snapshots.values[rightOffset + static_cast<std::size_t>(row)];
            }
            gram[static_cast<std::size_t>(left * columns + right)] = value;
            gram[static_cast<std::size_t>(right * columns + left)] = value;
        }
    }
#endif
    result.gramSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - gramStart).count();

    const auto eigenStart = std::chrono::steady_clock::now();
    EigenSystem eig = solveSymmetric(std::move(gram), columns);
    result.eigenSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - eigenStart).count();
    result.singularValues.resize(static_cast<std::size_t>(columns), 0.0);
    double totalEnergy = 0.0;
    for (int i = 0; i < columns; ++i) {
        const double eigenvalue = std::max(0.0, eig.values[static_cast<std::size_t>(i)]);
        result.singularValues[static_cast<std::size_t>(i)] = std::sqrt(eigenvalue);
        totalEnergy += eigenvalue;
    }
    const double sigma0 = result.singularValues.empty() ? 0.0 : result.singularValues.front();
    result.numericalRank = 0;
    for (double sigma : result.singularValues) {
        if (sigma0 > 0.0 && sigma > relativeSingularTolerance * sigma0) {
            ++result.numericalRank;
        }
    }
    int energyRank = result.numericalRank;
    bool energyRankFound = false;
    double cumulative = 0.0;
    result.retainedEnergy.resize(static_cast<std::size_t>(columns), 0.0);
    for (int i = 0; i < columns; ++i) {
        cumulative += result.singularValues[static_cast<std::size_t>(i)]
            * result.singularValues[static_cast<std::size_t>(i)];
        result.retainedEnergy[static_cast<std::size_t>(i)] = totalEnergy > 0.0
            ? cumulative / totalEnergy : 1.0;
        if (!energyRankFound && i < result.numericalRank
            && 1.0 - result.retainedEnergy[static_cast<std::size_t>(i)]
                <= discardedEnergyTolerance) {
            energyRank = i + 1;
            energyRankFound = true;
        }
    }
    result.selectedRank = requestedRank > 0
        ? std::min(requestedRank, result.numericalRank)
        : energyRank;
    result.selectedRank = std::min(result.selectedRank, rows);
    if (result.selectedRank <= 0) {
        throw std::runtime_error("[MOR] All interface snapshots are numerically rank deficient.");
    }

    const auto basisStart = std::chrono::steady_clock::now();
    result.basis.assign(static_cast<std::size_t>(rows * result.selectedRank), 0.0);
#ifdef USE_MKL_PARDISO
    std::vector<double> coefficients(static_cast<std::size_t>(
        columns * result.selectedRank), 0.0);
    for (int mode = 0; mode < result.selectedRank; ++mode) {
        const double sigma = result.singularValues[static_cast<std::size_t>(mode)];
        for (int snapshot = 0; snapshot < columns; ++snapshot) {
            coefficients[static_cast<std::size_t>(mode * columns + snapshot)] =
                eig.vectors[static_cast<std::size_t>(snapshot * columns + mode)] / sigma;
        }
    }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        rows, result.selectedRank, columns, 1.0,
        snapshots.values.data(), rows,
        coefficients.data(), columns,
        0.0, result.basis.data(), rows);
#else
    for (int mode = 0; mode < result.selectedRank; ++mode) {
        const double sigma = result.singularValues[static_cast<std::size_t>(mode)];
        for (int snapshot = 0; snapshot < columns; ++snapshot) {
            const double coefficient =
                eig.vectors[static_cast<std::size_t>(snapshot * columns + mode)] / sigma;
            const std::size_t snapshotOffset = static_cast<std::size_t>(snapshot * rows);
            const std::size_t modeOffset = static_cast<std::size_t>(mode * rows);
            for (int row = 0; row < rows; ++row) {
                result.basis[modeOffset + static_cast<std::size_t>(row)] += coefficient
                    * snapshots.values[snapshotOffset + static_cast<std::size_t>(row)];
            }
        }
    }
#endif
    for (int mode = 0; mode < result.selectedRank; ++mode) {
        // Two-pass modified Gram-Schmidt removes the small loss of
        // orthogonality introduced by the Gram formulation.
        const std::size_t modeOffset = static_cast<std::size_t>(mode * rows);
        for (int pass = 0; pass < 2; ++pass) {
            for (int previous = 0; previous < mode; ++previous) {
                const std::size_t previousOffset = static_cast<std::size_t>(previous * rows);
#ifdef USE_MKL_PARDISO
                const double projection = cblas_ddot(rows,
                    result.basis.data() + previousOffset, 1,
                    result.basis.data() + modeOffset, 1);
                cblas_daxpy(rows, -projection,
                    result.basis.data() + previousOffset, 1,
                    result.basis.data() + modeOffset, 1);
#else
                double projection = 0.0;
                for (int row = 0; row < rows; ++row) {
                    projection += result.basis[previousOffset + static_cast<std::size_t>(row)]
                        * result.basis[modeOffset + static_cast<std::size_t>(row)];
                }
                for (int row = 0; row < rows; ++row) {
                    result.basis[modeOffset + static_cast<std::size_t>(row)] -= projection
                        * result.basis[previousOffset + static_cast<std::size_t>(row)];
                }
#endif
            }
        }
#ifdef USE_MKL_PARDISO
        const double norm = cblas_dnrm2(
            rows, result.basis.data() + modeOffset, 1);
#else
        double norm = 0.0;
        for (int row = 0; row < rows; ++row) {
            const double value = result.basis[modeOffset + static_cast<std::size_t>(row)];
            norm += value * value;
        }
        norm = std::sqrt(norm);
#endif
        if (!(norm > 0.0)) {
            throw std::runtime_error("[MOR] POD basis orthogonalization failed.");
        }
#ifdef USE_MKL_PARDISO
        cblas_dscal(rows, 1.0 / norm,
            result.basis.data() + modeOffset, 1);
#else
        for (int row = 0; row < rows; ++row) {
            result.basis[modeOffset + static_cast<std::size_t>(row)] /= norm;
        }
#endif
    }
    double maximumOrthogonalityError = 0.0;
    for (int left = 0; left < result.selectedRank; ++left) {
        for (int right = 0; right <= left; ++right) {
#ifdef USE_MKL_PARDISO
            const double product = cblas_ddot(
                rows,
                result.basis.data() + static_cast<std::size_t>(left * rows), 1,
                result.basis.data() + static_cast<std::size_t>(right * rows), 1);
#else
            double product = 0.0;
            for (int row = 0; row < rows; ++row) {
                product += result.basis[static_cast<std::size_t>(row + left * rows)]
                    * result.basis[static_cast<std::size_t>(row + right * rows)];
            }
#endif
            maximumOrthogonalityError = std::max(
                maximumOrthogonalityError,
                std::abs(product - (left == right ? 1.0 : 0.0)));
        }
    }
    if (maximumOrthogonalityError > 1.0e-8) {
        throw std::runtime_error(
            "[MOR] POD basis failed the orthogonality gate: max|V^T V-I|="
            + std::to_string(maximumOrthogonalityError));
    }
    result.basisSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - basisStart).count();
    return result;
}

} // namespace mor
