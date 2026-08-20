#include "port_eigensolver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace mor::transient {
namespace {

double dot(const double* left, const double* right, int rows)
{
    long double value = 0.0L;
    for (int row = 0; row < rows; ++row) {
        value += static_cast<long double>(left[row]) * right[row];
    }
    return static_cast<double>(value);
}

double norm(const double* values, int rows)
{
    return std::sqrt(std::max(0.0, dot(values, values, rows)));
}

int orthonormalize(std::vector<double>& columns, int rows, int count,
                   double tolerance)
{
    std::vector<double> accepted;
    accepted.reserve(columns.size());
    int rank = 0;
    for (int column = 0; column < count; ++column) {
        std::vector<double> candidate(
            columns.begin() + static_cast<std::ptrdiff_t>(column * rows),
            columns.begin() + static_cast<std::ptrdiff_t>((column + 1) * rows));
        const double original = norm(candidate.data(), rows);
        if (!(original > 0.0)) continue;
        for (int pass = 0; pass < 2; ++pass) {
            for (int prior = 0; prior < rank; ++prior) {
                const double* basis = accepted.data()
                    + static_cast<std::size_t>(prior * rows);
                const double coefficient =
                    dot(basis, candidate.data(), rows);
                for (int row = 0; row < rows; ++row) {
                    candidate[static_cast<std::size_t>(row)] -=
                        coefficient * basis[row];
                }
            }
        }
        const double candidateNorm = norm(candidate.data(), rows);
        if (!(candidateNorm > tolerance * original)) continue;
        for (double& value : candidate) value /= candidateNorm;
        accepted.insert(
            accepted.end(), candidate.begin(), candidate.end());
        ++rank;
    }
    columns = std::move(accepted);
    return rank;
}

void symmetricJacobi(std::vector<double>& matrix, int size,
                     std::vector<double>& eigenvalues,
                     std::vector<double>& eigenvectors)
{
    eigenvectors.assign(static_cast<std::size_t>(size * size), 0.0);
    for (int row = 0; row < size; ++row) {
        eigenvectors[static_cast<std::size_t>(row * size + row)] = 1.0;
    }
    const int maximumSweeps = std::max(16, 10 * size * size);
    for (int sweep = 0; sweep < maximumSweeps; ++sweep) {
        int p = 0;
        int q = 0;
        double maximum = 0.0;
        for (int row = 0; row < size; ++row) {
            for (int column = row + 1; column < size; ++column) {
                const double value = std::abs(matrix[
                    static_cast<std::size_t>(row * size + column)]);
                if (value > maximum) {
                    maximum = value;
                    p = row;
                    q = column;
                }
            }
        }
        double diagonalScale = 0.0;
        for (int row = 0; row < size; ++row) {
            diagonalScale = std::max(diagonalScale,
                std::abs(matrix[static_cast<std::size_t>(
                    row * size + row)]));
        }
        if (maximum <= 64.0 * std::numeric_limits<double>::epsilon()
                * std::max(1.0, diagonalScale)) {
            break;
        }
        const double app =
            matrix[static_cast<std::size_t>(p * size + p)];
        const double aqq =
            matrix[static_cast<std::size_t>(q * size + q)];
        const double apq =
            matrix[static_cast<std::size_t>(p * size + q)];
        const double angle =
            0.5 * std::atan2(2.0 * apq, aqq - app);
        const double cosine = std::cos(angle);
        const double sine = std::sin(angle);
        for (int k = 0; k < size; ++k) {
            const double mkp =
                matrix[static_cast<std::size_t>(k * size + p)];
            const double mkq =
                matrix[static_cast<std::size_t>(k * size + q)];
            matrix[static_cast<std::size_t>(k * size + p)] =
                cosine * mkp - sine * mkq;
            matrix[static_cast<std::size_t>(k * size + q)] =
                sine * mkp + cosine * mkq;
        }
        for (int k = 0; k < size; ++k) {
            const double mpk =
                matrix[static_cast<std::size_t>(p * size + k)];
            const double mqk =
                matrix[static_cast<std::size_t>(q * size + k)];
            matrix[static_cast<std::size_t>(p * size + k)] =
                cosine * mpk - sine * mqk;
            matrix[static_cast<std::size_t>(q * size + k)] =
                sine * mpk + cosine * mqk;
        }
        for (int k = 0; k < size; ++k) {
            const double vkp =
                eigenvectors[static_cast<std::size_t>(k * size + p)];
            const double vkq =
                eigenvectors[static_cast<std::size_t>(k * size + q)];
            eigenvectors[static_cast<std::size_t>(k * size + p)] =
                cosine * vkp - sine * vkq;
            eigenvectors[static_cast<std::size_t>(k * size + q)] =
                sine * vkp + cosine * vkq;
        }
    }
    eigenvalues.resize(static_cast<std::size_t>(size));
    for (int row = 0; row < size; ++row) {
        eigenvalues[static_cast<std::size_t>(row)] =
            matrix[static_cast<std::size_t>(row * size + row)];
    }
}

} // namespace

MatrixFreeEigenResult solveLargestSymmetricEigenpairs(
    int dimension,
    const SymmetricOperatorApply& apply,
    const MatrixFreeEigenOptions& options)
{
    if (dimension <= 0 || options.requestedEigenpairs <= 0
        || options.maximumIterations <= 0
        || !(options.relativeTolerance > 0.0)
        || !(options.deflationTolerance > 0.0)) {
        throw std::runtime_error(
            "[Optimal port] Invalid matrix-free eigen options.");
    }
    MatrixFreeEigenResult result;
    const auto wallStart = std::chrono::steady_clock::now();
    const auto wallLimitReached = [&]() {
        return options.maximumWallSeconds > 0.0
            && std::chrono::duration<double>(
                std::chrono::steady_clock::now() - wallStart).count()
                >= options.maximumWallSeconds;
    };
    result.dimension = dimension;
    const int requested =
        std::min(dimension, options.requestedEigenpairs);
    int block = std::min(dimension,
        requested + std::max(0, options.oversamplingVectors));
    std::vector<double> basis(
        static_cast<std::size_t>(dimension * block), 0.0);
    for (int column = 0; column < block; ++column) {
        for (int row = 0; row < dimension; ++row) {
            basis[static_cast<std::size_t>(
                column * dimension + row)] =
                std::sin(0.7548776662466927
                    * static_cast<double>((row + 1) * (column + 2)))
                + std::cos(0.5698402909980532
                    * static_cast<double>((row + 3) * (column + 1)));
        }
    }
    block = orthonormalize(
        basis, dimension, block, options.deflationTolerance);
    if (block == 0) {
        result.status = "initial_block_deflated";
        return result;
    }

    for (int iteration = 0;
         iteration < options.maximumIterations; ++iteration) {
        std::vector<double> image(
            static_cast<std::size_t>(dimension * block), 0.0);
        for (int column = 0; column < block; ++column) {
            std::vector<double> input(
                basis.begin()
                    + static_cast<std::ptrdiff_t>(column * dimension),
                basis.begin()
                    + static_cast<std::ptrdiff_t>((column + 1) * dimension));
            std::vector<double> output;
            apply(input, output);
            ++result.operatorApplies;
            if (wallLimitReached()) {
                result.status = "time_limit";
                result.iterations = iteration;
                return result;
            }
            if (output.size() != static_cast<std::size_t>(dimension)) {
                throw std::runtime_error(
                    "[Optimal port] Eigen operator output size mismatch.");
            }
            std::copy(output.begin(), output.end(),
                image.begin()
                    + static_cast<std::ptrdiff_t>(column * dimension));
        }
        block = orthonormalize(
            image, dimension, block, options.deflationTolerance);
        if (block == 0) {
            result.status = "zero_operator";
            result.iterations = iteration + 1;
            return result;
        }
        basis = std::move(image);

        std::vector<double> applied(
            static_cast<std::size_t>(dimension * block), 0.0);
        for (int column = 0; column < block; ++column) {
            std::vector<double> input(
                basis.begin()
                    + static_cast<std::ptrdiff_t>(column * dimension),
                basis.begin()
                    + static_cast<std::ptrdiff_t>((column + 1) * dimension));
            std::vector<double> output;
            apply(input, output);
            ++result.operatorApplies;
            if (wallLimitReached()) {
                result.status = "time_limit";
                result.iterations = iteration;
                return result;
            }
            std::copy(output.begin(), output.end(),
                applied.begin()
                    + static_cast<std::ptrdiff_t>(column * dimension));
        }
        std::vector<double> rayleigh(
            static_cast<std::size_t>(block * block), 0.0);
        for (int row = 0; row < block; ++row) {
            for (int column = 0; column < block; ++column) {
                rayleigh[static_cast<std::size_t>(
                    row * block + column)] =
                    dot(basis.data()
                            + static_cast<std::size_t>(row * dimension),
                        applied.data()
                            + static_cast<std::size_t>(column * dimension),
                        dimension);
            }
        }
        for (int row = 0; row < block; ++row) {
            for (int column = 0; column < row; ++column) {
                const double average = 0.5
                    * (rayleigh[static_cast<std::size_t>(
                        row * block + column)]
                       + rayleigh[static_cast<std::size_t>(
                        column * block + row)]);
                rayleigh[static_cast<std::size_t>(
                    row * block + column)] = average;
                rayleigh[static_cast<std::size_t>(
                    column * block + row)] = average;
            }
        }
        std::vector<double> values;
        std::vector<double> vectors;
        symmetricJacobi(rayleigh, block, values, vectors);
        std::vector<int> order(static_cast<std::size_t>(block));
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(),
            [&](int left, int right) {
                return values[static_cast<std::size_t>(left)]
                    > values[static_cast<std::size_t>(right)];
            });

        std::vector<double> ritz(
            static_cast<std::size_t>(dimension * block), 0.0);
        std::vector<double> ritzApplied(
            static_cast<std::size_t>(dimension * block), 0.0);
        for (int selected = 0; selected < block; ++selected) {
            const int eigenColumn = order[static_cast<std::size_t>(selected)];
            for (int source = 0; source < block; ++source) {
                const double coefficient = vectors[
                    static_cast<std::size_t>(
                        source * block + eigenColumn)];
                for (int row = 0; row < dimension; ++row) {
                    ritz[static_cast<std::size_t>(
                        selected * dimension + row)] += coefficient
                        * basis[static_cast<std::size_t>(
                            source * dimension + row)];
                    ritzApplied[static_cast<std::size_t>(
                        selected * dimension + row)] += coefficient
                        * applied[static_cast<std::size_t>(
                            source * dimension + row)];
                }
            }
        }
        std::vector<double> residuals(
            static_cast<std::size_t>(block), 0.0);
        for (int selected = 0; selected < block; ++selected) {
            const double value = values[static_cast<std::size_t>(
                order[static_cast<std::size_t>(selected)])];
            double* residual = ritzApplied.data()
                + static_cast<std::size_t>(selected * dimension);
            const double* vector = ritz.data()
                + static_cast<std::size_t>(selected * dimension);
            const double imageNorm = norm(residual, dimension);
            for (int row = 0; row < dimension; ++row) {
                residual[row] -= value * vector[row];
            }
            residuals[static_cast<std::size_t>(selected)] =
                norm(residual, dimension)
                / std::max(std::numeric_limits<double>::epsilon(),
                    imageNorm);
        }
        basis = std::move(ritz);
        result.iterations = iteration + 1;
        const int retained = std::min(requested, block);
        result.rank = retained;
        result.eigenvalues.resize(static_cast<std::size_t>(retained));
        result.residuals.resize(static_cast<std::size_t>(retained));
        result.eigenvectors.resize(
            static_cast<std::size_t>(dimension * retained));
        bool converged = retained == requested;
        for (int selected = 0; selected < retained; ++selected) {
            result.eigenvalues[static_cast<std::size_t>(selected)] =
                std::max(0.0, values[static_cast<std::size_t>(
                    order[static_cast<std::size_t>(selected)])]);
            result.residuals[static_cast<std::size_t>(selected)] =
                residuals[static_cast<std::size_t>(selected)];
            std::copy_n(
                basis.begin()
                    + static_cast<std::ptrdiff_t>(selected * dimension),
                dimension,
                result.eigenvectors.begin()
                    + static_cast<std::ptrdiff_t>(selected * dimension));
            converged = converged
                && result.residuals[static_cast<std::size_t>(selected)]
                    <= options.relativeTolerance;
        }
        result.converged = converged;
        result.status = converged ? "success" : "iterating";
        if (converged) return result;
    }
    result.status = "maximum_iterations";
    return result;
}

} // namespace mor::transient
