#include "reduced_schur_model.hpp"

#include "ddm_schur/schur_fgmres.hpp"
#include "pod_basis.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

namespace mor {
namespace {

double dotColumn(const std::vector<double>& columnMajor,
                 int rows,
                 int column,
                 const std::vector<double>& vector)
{
#ifdef USE_MKL_PARDISO
    return cblas_ddot(rows,
        columnMajor.data() + static_cast<std::size_t>(column * rows), 1,
        vector.data(), 1);
#else
    double result = 0.0;
    const std::size_t offset = static_cast<std::size_t>(column * rows);
    for (int row = 0; row < rows; ++row) {
        result += columnMajor[offset + static_cast<std::size_t>(row)]
            * vector[static_cast<std::size_t>(row)];
    }
    return result;
#endif
}

bool factorCholesky(const std::vector<double>& matrix,
                    int n,
                    std::vector<double>& lower)
{
    lower.assign(static_cast<std::size_t>(n * n), 0.0);
    double diagonalScale = 0.0;
    for (int i = 0; i < n; ++i) {
        diagonalScale = std::max(diagonalScale,
            std::abs(matrix[static_cast<std::size_t>(i * n + i)]));
    }
    const double threshold = 128.0 * std::numeric_limits<double>::epsilon()
        * std::max(std::numeric_limits<double>::min(), diagonalScale);
    for (int row = 0; row < n; ++row) {
        for (int column = 0; column <= row; ++column) {
            double value = 0.5 * (matrix[static_cast<std::size_t>(row * n + column)]
                + matrix[static_cast<std::size_t>(column * n + row)]);
            for (int k = 0; k < column; ++k) {
                value -= lower[static_cast<std::size_t>(row * n + k)]
                    * lower[static_cast<std::size_t>(column * n + k)];
            }
            if (row == column) {
                if (!(value > threshold) || !std::isfinite(value)) {
                    lower.clear();
                    return false;
                }
                lower[static_cast<std::size_t>(row * n + column)] = std::sqrt(value);
            } else {
                lower[static_cast<std::size_t>(row * n + column)] = value
                    / lower[static_cast<std::size_t>(column * n + column)];
            }
        }
    }
    return true;
}

void factorLdlt(const std::vector<double>& matrix,
                int n,
                std::vector<double>& lower,
                std::vector<double>& diagonal)
{
    lower.assign(static_cast<std::size_t>(n * n), 0.0);
    diagonal.assign(static_cast<std::size_t>(n), 0.0);
    double scale = 0.0;
    for (double value : matrix) {
        scale = std::max(scale, std::abs(value));
    }
    const double threshold = 128.0 * std::numeric_limits<double>::epsilon()
        * std::max(std::numeric_limits<double>::min(), scale);
    for (int row = 0; row < n; ++row) {
        lower[static_cast<std::size_t>(row * n + row)] = 1.0;
        double pivot = matrix[static_cast<std::size_t>(row * n + row)];
        for (int k = 0; k < row; ++k) {
            const double value = lower[static_cast<std::size_t>(row * n + k)];
            pivot -= value * value * diagonal[static_cast<std::size_t>(k)];
        }
        if (!(std::abs(pivot) > threshold) || !std::isfinite(pivot)) {
            throw std::runtime_error(
                "[MOR] Dense Cholesky failed and dense LDLT encountered a zero/tiny pivot at row "
                + std::to_string(row));
        }
        diagonal[static_cast<std::size_t>(row)] = pivot;
        for (int next = row + 1; next < n; ++next) {
            double value = 0.5 * (matrix[static_cast<std::size_t>(next * n + row)]
                + matrix[static_cast<std::size_t>(row * n + next)]);
            for (int k = 0; k < row; ++k) {
                value -= lower[static_cast<std::size_t>(next * n + k)]
                    * lower[static_cast<std::size_t>(row * n + k)]
                    * diagonal[static_cast<std::size_t>(k)];
            }
            lower[static_cast<std::size_t>(next * n + row)] = value / pivot;
        }
    }
}

std::vector<double> solveSymmetricFactor(const std::vector<double>& lower,
                                         const std::vector<double>& diagonal,
                                         bool cholesky,
                                         std::vector<double> rhs)
{
    const int n = static_cast<int>(rhs.size());
    // L y = b.
    for (int row = 0; row < n; ++row) {
        for (int column = 0; column < row; ++column) {
            rhs[static_cast<std::size_t>(row)] -=
                lower[static_cast<std::size_t>(row * n + column)]
                * rhs[static_cast<std::size_t>(column)];
        }
        if (cholesky) {
            rhs[static_cast<std::size_t>(row)] /=
                lower[static_cast<std::size_t>(row * n + row)];
        }
    }
    if (!cholesky) {
        for (int row = 0; row < n; ++row) {
            rhs[static_cast<std::size_t>(row)] /= diagonal[static_cast<std::size_t>(row)];
        }
    }
    // L^T x = y.
    for (int row = n - 1; row >= 0; --row) {
        for (int column = row + 1; column < n; ++column) {
            rhs[static_cast<std::size_t>(row)] -=
                lower[static_cast<std::size_t>(column * n + row)]
                * rhs[static_cast<std::size_t>(column)];
        }
        if (cholesky) {
            rhs[static_cast<std::size_t>(row)] /=
                lower[static_cast<std::size_t>(row * n + row)];
        }
    }
    return rhs;
}

} // namespace

ReducedSchurModel buildReducedSchurModel(ddm_schur::DdmSchurSolver& solver,
                                         const PodResult& pod,
                                         int rank,
                                         const std::vector<double>& referenceInterface,
                                         const Fingerprints& fingerprints)
{
    rank = std::min(rank, pod.selectedRank);
    if (rank <= 0 || pod.rows != solver.interfaceDofs()
        || static_cast<int>(referenceInterface.size()) != solver.interfaceDofs()) {
        throw std::runtime_error("[MOR] Invalid dimensions while constructing reduced Schur model.");
    }
    ReducedSchurModel model;
    model.interfaceDofs = pod.rows;
    model.rank = rank;
    model.referenceInterface = referenceInterface;
    model.fingerprints = fingerprints;
    model.singularValues.assign(pod.singularValues.begin(),
        pod.singularValues.begin() + rank);
    model.basis.assign(pod.basis.begin(),
        pod.basis.begin() + static_cast<std::ptrdiff_t>(rank * pod.rows));
    model.schurBasis.assign(static_cast<std::size_t>(rank * pod.rows), 0.0);
    model.reducedOperator.assign(static_cast<std::size_t>(rank * rank), 0.0);
    std::vector<double> column(static_cast<std::size_t>(pod.rows), 0.0);
    for (int mode = 0; mode < rank; ++mode) {
        std::copy_n(model.basis.begin() + static_cast<std::ptrdiff_t>(mode * pod.rows),
                    pod.rows, column.begin());
        std::vector<double> image;
        const auto applyStart = std::chrono::steady_clock::now();
        solver.applyExactSchur(column, image);
        model.exactSchurApplySeconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - applyStart).count();
        std::copy(image.begin(), image.end(),
                  model.schurBasis.begin() + static_cast<std::ptrdiff_t>(mode * pod.rows));
        const auto assemblyStart = std::chrono::steady_clock::now();
        for (int test = 0; test < rank; ++test) {
            model.reducedOperator[static_cast<std::size_t>(test * rank + mode)] =
                dotColumn(model.basis, pod.rows, test, image);
        }
        model.reducedAssemblySeconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - assemblyStart).count();
    }

    double skewSquared = 0.0;
    double normSquared = 0.0;
    std::vector<double> symmetricPart(static_cast<std::size_t>(rank * rank), 0.0);
    for (int row = 0; row < rank; ++row) {
        for (int columnIndex = 0; columnIndex < rank; ++columnIndex) {
            const double value = model.reducedOperator[
                static_cast<std::size_t>(row * rank + columnIndex)];
            const double transpose = model.reducedOperator[
                static_cast<std::size_t>(columnIndex * rank + row)];
            const double skew = value - transpose;
            skewSquared += skew * skew;
            normSquared += value * value;
            symmetricPart[static_cast<std::size_t>(row * rank + columnIndex)] =
                0.5 * (value + transpose);
        }
    }
    model.symmetryError = std::sqrt(skewSquared) /
        std::max(1.0e-300, std::sqrt(normSquared));
    const std::vector<double> eigenvalues = symmetricEigenvalues(std::move(symmetricPart), rank);
    if (!eigenvalues.empty()) {
        model.maximumSymmetricEigenvalue = eigenvalues.front();
        model.minimumSymmetricEigenvalue = eigenvalues.back();
        model.conditionEstimate = std::abs(model.maximumSymmetricEigenvalue)
            / std::max(1.0e-300, std::abs(model.minimumSymmetricEigenvalue));
    }
    return model;
}

ReducedSchurModel truncateModel(const ReducedSchurModel& source, int rank)
{
    rank = std::min(rank, source.rank);
    if (rank <= 0) {
        throw std::runtime_error("[MOR] Requested truncation rank is invalid.");
    }
    ReducedSchurModel target = source;
    target.rank = rank;
    target.basis.resize(static_cast<std::size_t>(rank * source.interfaceDofs));
    target.schurBasis.resize(static_cast<std::size_t>(rank * source.interfaceDofs));
    target.singularValues.resize(static_cast<std::size_t>(rank));
    target.reducedOperator.assign(static_cast<std::size_t>(rank * rank), 0.0);
    for (int row = 0; row < rank; ++row) {
        for (int column = 0; column < rank; ++column) {
            target.reducedOperator[static_cast<std::size_t>(row * rank + column)] =
                source.reducedOperator[static_cast<std::size_t>(row * source.rank + column)];
        }
    }
    std::vector<double> symmetricPart(static_cast<std::size_t>(rank * rank), 0.0);
    double skewSquared = 0.0;
    double normSquared = 0.0;
    for (int row = 0; row < rank; ++row) {
        for (int column = 0; column < rank; ++column) {
            const double value = target.reducedOperator[static_cast<std::size_t>(row * rank + column)];
            const double transpose = target.reducedOperator[static_cast<std::size_t>(column * rank + row)];
            symmetricPart[static_cast<std::size_t>(row * rank + column)] = 0.5 * (value + transpose);
            skewSquared += (value - transpose) * (value - transpose);
            normSquared += value * value;
        }
    }
    target.symmetryError = std::sqrt(skewSquared) / std::max(1.0e-300, std::sqrt(normSquared));
    const std::vector<double> eigenvalues = symmetricEigenvalues(std::move(symmetricPart), rank);
    target.maximumSymmetricEigenvalue = eigenvalues.front();
    target.minimumSymmetricEigenvalue = eigenvalues.back();
    target.conditionEstimate = std::abs(eigenvalues.front())
        / std::max(1.0e-300, std::abs(eigenvalues.back()));
    return target;
}

ReducedSchurOnlineSolver::ReducedSchurOnlineSolver(
    ddm_schur::DdmSchurSolver& solver,
    const ReducedSchurModel& model)
    : solver_(solver), model_(model)
{
    if (model_.interfaceDofs != solver_.interfaceDofs()) {
        throw std::runtime_error("[MOR] Loaded basis does not match the Schur interface dimension.");
    }
    solver_.applyExactSchur(model_.referenceInterface, referenceSchurImage_);
    choleskyFactor_ = factorCholesky(
        model_.reducedOperator, model_.rank, denseFactor_);
    if (!choleskyFactor_) {
        std::cerr << "[MOR] Warning: reduced Schur Cholesky factorization failed; "
                     "using the explicit symmetric LDLT fallback.\n";
        factorLdlt(model_.reducedOperator, model_.rank,
                   denseFactor_, diagonalFactor_);
    }
}

const char* ReducedSchurOnlineSolver::factorizationType() const
{
    return choleskyFactor_ ? "cholesky" : "ldlt";
}

SolveResult ReducedSchurOnlineSolver::solve(const std::vector<double>& globalRhs,
                                            bool corrected)
{
    SolveResult result;
    const auto totalStart = std::chrono::steady_clock::now();
    const auto condensationStart = std::chrono::steady_clock::now();
    const std::vector<double> condensed = solver_.condensedRhs(globalRhs);
    result.timing.condensationSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - condensationStart).count();

    const auto projectionStart = std::chrono::steady_clock::now();
    std::vector<double> shifted(condensed.size(), 0.0);
    for (std::size_t i = 0; i < condensed.size(); ++i) {
        shifted[i] = condensed[i] - referenceSchurImage_[i];
    }
    std::vector<double> reducedRhs(static_cast<std::size_t>(model_.rank), 0.0);
    for (int mode = 0; mode < model_.rank; ++mode) {
        reducedRhs[static_cast<std::size_t>(mode)] =
            dotColumn(model_.basis, model_.interfaceDofs, mode, shifted);
    }
    result.timing.projectionSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - projectionStart).count();

    const auto solveStart = std::chrono::steady_clock::now();
    const std::vector<double> coefficients = solveSymmetricFactor(
        denseFactor_, diagonalFactor_, choleskyFactor_, std::move(reducedRhs));
    result.timing.reducedSolveSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - solveStart).count();
    const auto reconstructionStart = std::chrono::steady_clock::now();
    result.interfaceTemperature = model_.referenceInterface;
    for (int mode = 0; mode < model_.rank; ++mode) {
        const double coefficient = coefficients[static_cast<std::size_t>(mode)];
        const std::size_t offset = static_cast<std::size_t>(mode * model_.interfaceDofs);
        for (int row = 0; row < model_.interfaceDofs; ++row) {
            result.interfaceTemperature[static_cast<std::size_t>(row)] += coefficient
                * model_.basis[offset + static_cast<std::size_t>(row)];
        }
    }
    result.timing.interfaceReconstructionSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - reconstructionStart).count();

    if (corrected) {
        const auto correctionStart = std::chrono::steady_clock::now();
        ddm_schur::SolveResult correctedResult =
            solver_.solveCorrection(globalRhs, result.interfaceTemperature);
        result.timing.correctionSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - correctionStart).count();
        result.interfaceTemperature = std::move(correctedResult.interfaceSolution);
        result.temperature = std::move(correctedResult.temperature);
        result.correctionIterations = correctedResult.report.iterations;
        result.error.interfaceRelativeResidual =
            correctedResult.report.interfaceRelativeResidual;
        result.status = correctedResult.report.status;
    } else {
        std::vector<double> image;
        solver_.applyExactSchur(result.interfaceTemperature, image);
        double residualNormSquared = 0.0;
        double rhsNormSquared = 0.0;
        for (std::size_t i = 0; i < condensed.size(); ++i) {
            const double residual = condensed[i] - image[i];
            residualNormSquared += residual * residual;
            rhsNormSquared += condensed[i] * condensed[i];
        }
        result.error.interfaceRelativeResidual = std::sqrt(residualNormSquared)
            / std::max(1.0e-300, std::sqrt(rhsNormSquared));
        const auto recoveryStart = std::chrono::steady_clock::now();
        result.temperature = solver_.recover(globalRhs, result.interfaceTemperature);
        result.timing.recoverySeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - recoveryStart).count();
        result.status = "success";
    }
    result.timing.totalSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - totalStart).count();
    return result;
}

} // namespace mor
