#include "schur_fgmres.hpp"
#include "schur_operator.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace ddm_schur {
namespace {

double dotVector(const std::vector<double>& a, const std::vector<double>& b)
{
    double value = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        value += a[i] * b[i];
    }
    return value;
}

double normVector(const std::vector<double>& value)
{
    return std::sqrt(std::max(0.0, dotVector(value, value)));
}

void axpy(double alpha, const std::vector<double>& x, std::vector<double>& y)
{
    for (std::size_t i = 0; i < x.size(); ++i) {
        y[i] += alpha * x[i];
    }
}

std::vector<double> solveLeastSquaresCoefficients(
    const std::vector<std::vector<double>>& hessenberg,
    const std::vector<double>& rotatedRhs,
    int used)
{
    std::vector<double> coefficients(static_cast<std::size_t>(used), 0.0);
    for (int row = used - 1; row >= 0; --row) {
        double value = rotatedRhs[static_cast<std::size_t>(row)];
        for (int col = row + 1; col < used; ++col) {
            value -= hessenberg[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]
                * coefficients[static_cast<std::size_t>(col)];
        }
        const double diagonal = hessenberg[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
        if (std::abs(diagonal) <= 1.0e-30) {
            throw std::runtime_error("[Schur] FGMRES triangular solve broke down.");
        }
        coefficients[static_cast<std::size_t>(row)] = value / diagonal;
    }
    return coefficients;
}

struct FgmresResult {
    std::vector<double> solution;
    int iterations = 0;
    double relativeResidual = std::numeric_limits<double>::infinity();
    bool converged = false;
};

FgmresResult solveFgmres(SchurOperator& schur,
                         const std::vector<double>& rhs,
                         const Options& options)
{
    FgmresResult result;
    const int n = static_cast<int>(rhs.size());
    result.solution.assign(rhs.size(), 0.0);
    if (n == 0) {
        result.relativeResidual = 0.0;
        result.converged = true;
        return result;
    }

    const double rhsNorm = normVector(rhs);
    if (rhsNorm == 0.0) {
        result.relativeResidual = 0.0;
        result.converged = true;
        return result;
    }

    const int restart = std::max(1, std::min(options.restart, n));
    std::vector<double> product;
    schur.apply(result.solution, product);
    std::vector<double> residual(rhs.size(), 0.0);
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        residual[i] = rhs[i] - product[i];
    }

    while (result.iterations < options.maxIterations) {
        const double beta = normVector(residual);
        result.relativeResidual = beta / rhsNorm;
        if (result.relativeResidual <= options.relativeTolerance) {
            result.converged = true;
            break;
        }

        std::vector<std::vector<double>> basis(static_cast<std::size_t>(restart + 1),
                                                std::vector<double>(rhs.size(), 0.0));
        std::vector<std::vector<double>> preconditioned(static_cast<std::size_t>(restart),
                                                         std::vector<double>(rhs.size(), 0.0));
        std::vector<std::vector<double>> hessenberg(static_cast<std::size_t>(restart + 1),
                                                     std::vector<double>(static_cast<std::size_t>(restart), 0.0));
        std::vector<double> cosine(static_cast<std::size_t>(restart), 0.0);
        std::vector<double> sine(static_cast<std::size_t>(restart), 0.0);
        std::vector<double> rotatedRhs(static_cast<std::size_t>(restart + 1), 0.0);
        rotatedRhs[0] = beta;
        for (std::size_t i = 0; i < residual.size(); ++i) {
            basis[0][i] = residual[i] / beta;
        }

        int used = 0;
        bool estimatedConverged = false;
        for (int column = 0; column < restart && result.iterations < options.maxIterations; ++column) {
            schur.applyBlockPreconditioner(basis[static_cast<std::size_t>(column)],
                                           preconditioned[static_cast<std::size_t>(column)]);
            std::vector<double> work;
            schur.apply(preconditioned[static_cast<std::size_t>(column)], work);

            for (int row = 0; row <= column; ++row) {
                hessenberg[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] =
                    dotVector(basis[static_cast<std::size_t>(row)], work);
                axpy(-hessenberg[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)],
                     basis[static_cast<std::size_t>(row)], work);
            }
            // A second MGS pass is important for long, unrestarted interface solves.
            // It keeps the Arnoldi basis orthogonal when the Schur operator is ill-conditioned.
            for (int row = 0; row <= column; ++row) {
                const double correction = dotVector(basis[static_cast<std::size_t>(row)], work);
                hessenberg[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] += correction;
                axpy(-correction, basis[static_cast<std::size_t>(row)], work);
            }
            const double nextNorm = normVector(work);
            hessenberg[static_cast<std::size_t>(column + 1)][static_cast<std::size_t>(column)] = nextNorm;
            if (nextNorm > 10.0 * std::numeric_limits<double>::epsilon()) {
                for (std::size_t i = 0; i < work.size(); ++i) {
                    basis[static_cast<std::size_t>(column + 1)][i] = work[i] / nextNorm;
                }
            }

            for (int row = 0; row < column; ++row) {
                const double upper = hessenberg[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)];
                const double lower = hessenberg[static_cast<std::size_t>(row + 1)][static_cast<std::size_t>(column)];
                hessenberg[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] =
                    cosine[static_cast<std::size_t>(row)] * upper
                    + sine[static_cast<std::size_t>(row)] * lower;
                hessenberg[static_cast<std::size_t>(row + 1)][static_cast<std::size_t>(column)] =
                    -sine[static_cast<std::size_t>(row)] * upper
                    + cosine[static_cast<std::size_t>(row)] * lower;
            }

            const double diagonal = hessenberg[static_cast<std::size_t>(column)][static_cast<std::size_t>(column)];
            const double subdiagonal = hessenberg[static_cast<std::size_t>(column + 1)][static_cast<std::size_t>(column)];
            const double magnitude = std::hypot(diagonal, subdiagonal);
            cosine[static_cast<std::size_t>(column)] = magnitude == 0.0 ? 1.0 : diagonal / magnitude;
            sine[static_cast<std::size_t>(column)] = magnitude == 0.0 ? 0.0 : subdiagonal / magnitude;
            hessenberg[static_cast<std::size_t>(column)][static_cast<std::size_t>(column)] = magnitude;
            hessenberg[static_cast<std::size_t>(column + 1)][static_cast<std::size_t>(column)] = 0.0;

            const double currentRhs = rotatedRhs[static_cast<std::size_t>(column)];
            rotatedRhs[static_cast<std::size_t>(column)] =
                cosine[static_cast<std::size_t>(column)] * currentRhs;
            rotatedRhs[static_cast<std::size_t>(column + 1)] =
                -sine[static_cast<std::size_t>(column)] * currentRhs;

            ++result.iterations;
            used = column + 1;
            result.relativeResidual =
                std::abs(rotatedRhs[static_cast<std::size_t>(column + 1)]) / rhsNorm;
            if (result.relativeResidual <= options.relativeTolerance) {
                const std::vector<double> trialCoefficients =
                    solveLeastSquaresCoefficients(hessenberg, rotatedRhs, used);
                std::vector<double> trialSolution = result.solution;
                for (int i = 0; i < used; ++i) {
                    axpy(trialCoefficients[static_cast<std::size_t>(i)],
                         preconditioned[static_cast<std::size_t>(i)],
                         trialSolution);
                }
                std::vector<double> trialProduct;
                schur.apply(trialSolution, trialProduct);
                std::vector<double> trialResidual(rhs.size(), 0.0);
                for (std::size_t i = 0; i < rhs.size(); ++i) {
                    trialResidual[i] = rhs[i] - trialProduct[i];
                }
                if (normVector(trialResidual) / rhsNorm <= options.relativeTolerance) {
                    estimatedConverged = true;
                    break;
                }
            }
            if (nextNorm == 0.0) {
                break;
            }
        }

        const std::vector<double> coefficients =
            solveLeastSquaresCoefficients(hessenberg, rotatedRhs, used);
        for (int i = 0; i < used; ++i) {
            axpy(coefficients[static_cast<std::size_t>(i)],
                 preconditioned[static_cast<std::size_t>(i)], result.solution);
        }

        schur.apply(result.solution, product);
        for (std::size_t i = 0; i < rhs.size(); ++i) {
            residual[i] = rhs[i] - product[i];
        }
        result.relativeResidual = normVector(residual) / rhsNorm;
        std::cout << "[Schur] FGMRES cycle: used=" << used
                  << ", total=" << result.iterations
                  << ", estimated="
                  << std::abs(rotatedRhs[static_cast<std::size_t>(used)]) / rhsNorm
                  << ", true=" << result.relativeResidual << '\n';
        if (result.relativeResidual <= options.relativeTolerance) {
            result.converged = true;
            break;
        }
        if (estimatedConverged && used == 0) {
            break;
        }
    }
    return result;
}

} // namespace

MatrixFreeFgmresResult solveMatrixFreeFgmres(
    const std::vector<double>& rhs,
    const Options& options,
    const MatrixFreeApply& applyOperator,
    const MatrixFreeApply& applyPreconditioner,
    const std::vector<double>* initialSolution,
    MatrixFreeFgmresWorkspace* workspace,
    const std::vector<double>* initialProduct)
{
    MatrixFreeFgmresResult result;
    const int n = static_cast<int>(rhs.size());
    result.solution.assign(rhs.size(), 0.0);
    if (initialSolution) {
        if (initialSolution->size() != rhs.size()) {
            throw std::runtime_error(
                "[Local ROM] Matrix-free FGMRES initial solution has the wrong size.");
        }
        result.solution = *initialSolution;
    }
    if (n == 0) {
        result.relativeResidual = 0.0;
        result.converged = true;
        return result;
    }
    if (!applyOperator || !applyPreconditioner || options.maxIterations <= 0
        || options.restart <= 0 || !(options.relativeTolerance > 0.0)) {
        throw std::runtime_error("[Local ROM] Invalid matrix-free FGMRES options.");
    }

    const double rhsNorm = normVector(rhs);
    if (rhsNorm == 0.0) {
        result.relativeResidual = 0.0;
        result.converged = true;
        return result;
    }

    MatrixFreeFgmresWorkspace ownedWorkspace;
    MatrixFreeFgmresWorkspace& storage = workspace
        ? *workspace : ownedWorkspace;
    storage.product.resize(rhs.size());
    storage.residual.resize(rhs.size());
    auto& product = storage.product;
    auto& residual = storage.residual;
    auto timedApply = [&](const MatrixFreeApply& apply,
                          const std::vector<double>& input,
                          std::vector<double>& output,
                          double& seconds) {
        const auto start = std::chrono::steady_clock::now();
        apply(input, output);
        seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
    };
    const std::vector<double>* initialImage = &product;
    if (initialProduct) {
        if (!initialSolution || initialProduct->size() != rhs.size()) {
            throw std::runtime_error(
                "[Local ROM] Matrix-free FGMRES initial product has the wrong size.");
        }
        initialImage = initialProduct;
    } else {
        timedApply(
            applyOperator, result.solution, product, result.operatorSeconds);
    }
    result.solutionProductAvailable = true;
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        residual[i] = rhs[i] - (*initialImage)[i];
    }
    result.initialRelativeResidual = normVector(residual) / rhsNorm;
    result.relativeResidual = result.initialRelativeResidual;
    if (result.relativeResidual <= options.relativeTolerance) {
        if (initialImage != &product) {
            product = *initialImage;
        }
        result.converged = true;
        return result;
    }

    const int restart = std::max(1, std::min(options.restart, n));
    storage.basis.resize(static_cast<std::size_t>(restart + 1));
    for (std::vector<double>& value : storage.basis) {
        value.resize(rhs.size());
    }
    storage.preconditioned.resize(static_cast<std::size_t>(restart));
    for (std::vector<double>& value : storage.preconditioned) {
        value.resize(rhs.size());
    }
    storage.hessenberg.resize(static_cast<std::size_t>(restart + 1));
    for (std::vector<double>& value : storage.hessenberg) {
        value.resize(static_cast<std::size_t>(restart));
    }
    storage.cosine.resize(static_cast<std::size_t>(restart));
    storage.sine.resize(static_cast<std::size_t>(restart));
    storage.rotatedRhs.resize(static_cast<std::size_t>(restart + 1));
    storage.work.resize(rhs.size());
    storage.trialSolution.resize(rhs.size());
    storage.trialProduct.resize(rhs.size());
    storage.trialResidual.resize(rhs.size());
    auto& basis = storage.basis;
    auto& preconditioned = storage.preconditioned;
    auto& hessenberg = storage.hessenberg;
    auto& cosine = storage.cosine;
    auto& sine = storage.sine;
    auto& rotatedRhs = storage.rotatedRhs;
    auto& work = storage.work;
    auto& trialSolution = storage.trialSolution;
    auto& trialProduct = storage.trialProduct;
    auto& trialResidual = storage.trialResidual;
    while (result.iterations < options.maxIterations) {
        const double beta = normVector(residual);
        result.relativeResidual = beta / rhsNorm;
        if (result.relativeResidual <= options.relativeTolerance) {
            result.converged = true;
            break;
        }

        for (std::vector<double>& row : hessenberg) {
            std::fill(row.begin(), row.end(), 0.0);
        }
        std::fill(cosine.begin(), cosine.end(), 0.0);
        std::fill(sine.begin(), sine.end(), 0.0);
        std::fill(rotatedRhs.begin(), rotatedRhs.end(), 0.0);
        rotatedRhs[0] = beta;
        for (std::size_t i = 0; i < residual.size(); ++i) {
            basis[0][i] = residual[i] / beta;
        }

        int used = 0;
        bool estimatedConverged = false;
        double trialRelativeResidual =
            std::numeric_limits<double>::infinity();
        for (int column = 0;
             column < restart && result.iterations < options.maxIterations;
             ++column) {
            timedApply(
                applyPreconditioner,
                basis[static_cast<std::size_t>(column)],
                preconditioned[static_cast<std::size_t>(column)],
                result.preconditionerSeconds);
            timedApply(
                applyOperator,
                preconditioned[static_cast<std::size_t>(column)], work,
                result.operatorSeconds);
            if (work.size() != rhs.size()) {
                throw std::runtime_error(
                    "[Local ROM] Matrix-free Schur apply returned the wrong size.");
            }

            const auto orthogonalizationStart =
                std::chrono::steady_clock::now();
            for (int row = 0; row <= column; ++row) {
                hessenberg[static_cast<std::size_t>(row)]
                          [static_cast<std::size_t>(column)] =
                    dotVector(basis[static_cast<std::size_t>(row)], work);
                axpy(-hessenberg[static_cast<std::size_t>(row)]
                                  [static_cast<std::size_t>(column)],
                     basis[static_cast<std::size_t>(row)], work);
            }
            for (int row = 0; row <= column; ++row) {
                const double correction =
                    dotVector(basis[static_cast<std::size_t>(row)], work);
                hessenberg[static_cast<std::size_t>(row)]
                          [static_cast<std::size_t>(column)] += correction;
                axpy(-correction, basis[static_cast<std::size_t>(row)], work);
            }
            const double nextNorm = normVector(work);
            hessenberg[static_cast<std::size_t>(column + 1)]
                      [static_cast<std::size_t>(column)] = nextNorm;
            if (nextNorm > 10.0 * std::numeric_limits<double>::epsilon()) {
                for (std::size_t i = 0; i < work.size(); ++i) {
                    basis[static_cast<std::size_t>(column + 1)][i] = work[i] / nextNorm;
                }
            } else {
                std::fill(
                    basis[static_cast<std::size_t>(column + 1)].begin(),
                    basis[static_cast<std::size_t>(column + 1)].end(), 0.0);
            }
            result.orthogonalizationSeconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now()
                - orthogonalizationStart).count();

            for (int row = 0; row < column; ++row) {
                const double upper = hessenberg[static_cast<std::size_t>(row)]
                                              [static_cast<std::size_t>(column)];
                const double lower = hessenberg[static_cast<std::size_t>(row + 1)]
                                              [static_cast<std::size_t>(column)];
                hessenberg[static_cast<std::size_t>(row)]
                          [static_cast<std::size_t>(column)] =
                    cosine[static_cast<std::size_t>(row)] * upper
                    + sine[static_cast<std::size_t>(row)] * lower;
                hessenberg[static_cast<std::size_t>(row + 1)]
                          [static_cast<std::size_t>(column)] =
                    -sine[static_cast<std::size_t>(row)] * upper
                    + cosine[static_cast<std::size_t>(row)] * lower;
            }

            const double diagonal = hessenberg[static_cast<std::size_t>(column)]
                                               [static_cast<std::size_t>(column)];
            const double subdiagonal = hessenberg[static_cast<std::size_t>(column + 1)]
                                                  [static_cast<std::size_t>(column)];
            const double magnitude = std::hypot(diagonal, subdiagonal);
            cosine[static_cast<std::size_t>(column)] =
                magnitude == 0.0 ? 1.0 : diagonal / magnitude;
            sine[static_cast<std::size_t>(column)] =
                magnitude == 0.0 ? 0.0 : subdiagonal / magnitude;
            hessenberg[static_cast<std::size_t>(column)]
                      [static_cast<std::size_t>(column)] = magnitude;
            hessenberg[static_cast<std::size_t>(column + 1)]
                      [static_cast<std::size_t>(column)] = 0.0;

            const double currentRhs = rotatedRhs[static_cast<std::size_t>(column)];
            rotatedRhs[static_cast<std::size_t>(column)] =
                cosine[static_cast<std::size_t>(column)] * currentRhs;
            rotatedRhs[static_cast<std::size_t>(column + 1)] =
                -sine[static_cast<std::size_t>(column)] * currentRhs;

            ++result.iterations;
            used = column + 1;
            result.relativeResidual =
                std::abs(rotatedRhs[static_cast<std::size_t>(column + 1)]) / rhsNorm;
            if (result.relativeResidual <= options.relativeTolerance) {
                const auto vectorUpdateStart =
                    std::chrono::steady_clock::now();
                const std::vector<double> trialCoefficients =
                    solveLeastSquaresCoefficients(hessenberg, rotatedRhs, used);
                trialSolution = result.solution;
                for (int i = 0; i < used; ++i) {
                    axpy(trialCoefficients[static_cast<std::size_t>(i)],
                         preconditioned[static_cast<std::size_t>(i)], trialSolution);
                }
                result.vectorUpdateSeconds += std::chrono::duration<double>(
                    std::chrono::steady_clock::now()
                    - vectorUpdateStart).count();
                timedApply(
                    applyOperator, trialSolution, trialProduct,
                    result.operatorSeconds);
                for (std::size_t i = 0; i < rhs.size(); ++i) {
                    trialResidual[i] = rhs[i] - trialProduct[i];
                }
                trialRelativeResidual = normVector(trialResidual) / rhsNorm;
                if (trialRelativeResidual <= options.relativeTolerance) {
                    estimatedConverged = true;
                    break;
                }
            }
            if (nextNorm == 0.0) {
                break;
            }
        }

        if (used == 0) {
            break;
        }
        if (estimatedConverged) {
            result.solution.swap(trialSolution);
            product.swap(trialProduct);
            residual.swap(trialResidual);
            result.relativeResidual = trialRelativeResidual;
        } else {
            const auto vectorUpdateStart = std::chrono::steady_clock::now();
            const std::vector<double> coefficients =
                solveLeastSquaresCoefficients(hessenberg, rotatedRhs, used);
            for (int i = 0; i < used; ++i) {
                axpy(coefficients[static_cast<std::size_t>(i)],
                     preconditioned[static_cast<std::size_t>(i)], result.solution);
            }
            result.vectorUpdateSeconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - vectorUpdateStart).count();

            timedApply(
                applyOperator, result.solution, product,
                result.operatorSeconds);
            for (std::size_t i = 0; i < rhs.size(); ++i) {
                residual[i] = rhs[i] - product[i];
            }
            result.relativeResidual = normVector(residual) / rhsNorm;
        }
        std::cout << "[Local ROM] FGMRES cycle: used=" << used
                  << ", total=" << result.iterations
                  << ", estimated="
                  << std::abs(rotatedRhs[static_cast<std::size_t>(used)]) / rhsNorm
                  << ", true=" << result.relativeResidual << '\n';
        if (result.relativeResidual <= options.relativeTolerance) {
            result.converged = true;
            break;
        }
    }
    return result;
}

MatrixFreeFgmresResult solveMatrixFreePcg(
    const std::vector<double>& rhs,
    const Options& options,
    const MatrixFreeApply& applyOperator,
    const MatrixFreeApply& applyPreconditioner,
    const std::vector<double>* initialSolution)
{
    MatrixFreeFgmresResult result;
    result.actualSolver = "pcg";
    result.solution.assign(rhs.size(), 0.0);
    if (initialSolution) {
        if (initialSolution->size() != rhs.size()) {
            throw std::runtime_error(
                "[Local ROM] Matrix-free PCG initial solution has the wrong size.");
        }
        result.solution = *initialSolution;
    }
    if (!applyOperator || !applyPreconditioner || options.maxIterations <= 0
        || !(options.relativeTolerance > 0.0)) {
        throw std::runtime_error("[Local ROM] Invalid matrix-free PCG options.");
    }
    if (rhs.empty()) {
        result.relativeResidual = 0.0;
        result.converged = true;
        return result;
    }
    const double rhsNorm = normVector(rhs);
    if (rhsNorm == 0.0) {
        result.solution.assign(rhs.size(), 0.0);
        result.relativeResidual = 0.0;
        result.converged = true;
        return result;
    }

    auto timedApply = [&](const MatrixFreeApply& apply,
                          const std::vector<double>& input,
                          std::vector<double>& output,
                          double& seconds) {
        const auto start = std::chrono::steady_clock::now();
        apply(input, output);
        seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
    };
    std::vector<double> product;
    timedApply(
        applyOperator, result.solution, product, result.operatorSeconds);
    if (product.size() != rhs.size()) {
        throw std::runtime_error(
            "[Local ROM] Matrix-free Schur apply returned the wrong size.");
    }
    std::vector<double> residual(rhs.size(), 0.0);
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        residual[i] = rhs[i] - product[i];
    }
    result.initialRelativeResidual = normVector(residual) / rhsNorm;
    result.relativeResidual = result.initialRelativeResidual;
    if (result.relativeResidual <= options.relativeTolerance) {
        result.converged = true;
        return result;
    }

    std::vector<double> preconditioned;
    timedApply(
        applyPreconditioner, residual, preconditioned,
        result.preconditionerSeconds);
    if (preconditioned.size() != rhs.size()) {
        throw std::runtime_error(
            "[Local ROM] Matrix-free preconditioner returned the wrong size.");
    }
    double residualPreconditioned = dotVector(residual, preconditioned);
    if (!(residualPreconditioned > 0.0)
        || !std::isfinite(residualPreconditioned)) {
        result.fallbackReason = "nonpositive_preconditioner_curvature";
        return result;
    }
    std::vector<double> direction = preconditioned;

    while (result.iterations < options.maxIterations) {
        std::vector<double> image;
        timedApply(
            applyOperator, direction, image, result.operatorSeconds);
        if (image.size() != rhs.size()) {
            throw std::runtime_error(
                "[Local ROM] Matrix-free Schur apply returned the wrong size.");
        }
        const double curvature = dotVector(direction, image);
        if (!(curvature > 0.0) || !std::isfinite(curvature)) {
            result.fallbackReason = "nonpositive_operator_curvature";
            break;
        }
        const double alpha = residualPreconditioned / curvature;
        if (!std::isfinite(alpha)) {
            result.fallbackReason = "nonfinite_pcg_step";
            break;
        }
        const auto vectorUpdateStart = std::chrono::steady_clock::now();
        axpy(alpha, direction, result.solution);
        axpy(-alpha, image, residual);
        result.vectorUpdateSeconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - vectorUpdateStart).count();
        ++result.iterations;

        const bool exactCheck = result.iterations % 50 == 0
            || normVector(residual) / rhsNorm <= options.relativeTolerance;
        if (exactCheck) {
            timedApply(
                applyOperator, result.solution, product,
                result.operatorSeconds);
            for (std::size_t i = 0; i < rhs.size(); ++i) {
                residual[i] = rhs[i] - product[i];
            }
            result.relativeResidual = normVector(residual) / rhsNorm;
            if (result.relativeResidual <= options.relativeTolerance) {
                result.converged = true;
                break;
            }
            // A true-residual replacement is a PCG restart.  This prevents
            // accumulated recursive-residual drift from corrupting the next
            // conjugacy coefficient.
            timedApply(
                applyPreconditioner, residual, preconditioned,
                result.preconditionerSeconds);
            residualPreconditioned = dotVector(residual, preconditioned);
            if (!(residualPreconditioned > 0.0)
                || !std::isfinite(residualPreconditioned)) {
                result.fallbackReason = "nonpositive_preconditioner_curvature";
                break;
            }
            direction = preconditioned;
            continue;
        }

        timedApply(
            applyPreconditioner, residual, preconditioned,
            result.preconditionerSeconds);
        const double nextResidualPreconditioned =
            dotVector(residual, preconditioned);
        if (!(nextResidualPreconditioned > 0.0)
            || !std::isfinite(nextResidualPreconditioned)) {
            result.fallbackReason = "nonpositive_preconditioner_curvature";
            break;
        }
        const double beta = nextResidualPreconditioned
            / residualPreconditioned;
        if (!std::isfinite(beta)) {
            result.fallbackReason = "nonfinite_pcg_step";
            break;
        }
        const auto directionUpdateStart = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < direction.size(); ++i) {
            direction[i] = preconditioned[i] + beta * direction[i];
        }
        result.vectorUpdateSeconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - directionUpdateStart).count();
        residualPreconditioned = nextResidualPreconditioned;
    }

    if (!result.converged) {
        applyOperator(result.solution, product);
        for (std::size_t i = 0; i < rhs.size(); ++i) {
            residual[i] = rhs[i] - product[i];
        }
        result.relativeResidual = normVector(residual) / rhsNorm;
        if (result.relativeResidual <= options.relativeTolerance) {
            result.converged = true;
            result.fallbackReason.clear();
        } else if (result.fallbackReason.empty()) {
            result.fallbackReason = "maximum_iterations";
        }
    }
    return result;
}

struct DdmSchurSolver::Impl {
    Options options;
    SchurOperator schur;
    Report setup;

    Impl(const Mesh& mesh,
         const SparseMatrix& system,
         const CaseConfig& physics,
         Options requested)
        : options(requested),
          schur(mesh,
                system,
                physics,
                requested.coarseLinearXY,
                requested.coarseLinearZ,
                requested.coarseGlobalQuadraticZ,
                requested.coarseInterfacePatches,
                requested.coarseInterfaceLinearXY,
                requested.coarseEnergyAdaptive,
                requested.energyMaxModesPerDomain,
                requested.energySubspaceIterations,
                requested.energyEigenvalueThreshold,
                requested.coarseGlobalSlow,
                requested.globalSlowModes,
                requested.globalSlowSubspaceDimension,
                requested.proxyEnabled,
                requested.proxyDisableCoarse,
                requested.proxyDiagnostics,
                requested.proxyHighConductivityThreshold,
                requested.proxyUseMaterialConnectivity,
                requested.proxyRing,
                requested.proxyProbeColumns,
                requested.proxyBlockSize,
                requested.proxyValidateBlockEquivalence,
                requested.localSolveThreads,
                requested.localPardisoThreads,
                requested.proxyCacheEnabled,
                requested.proxyCachePath,
                requested.proxyOutputDirectory)
    {
        setup.domains = schur.domains();
        setup.totalDofs = schur.totalDofs();
        setup.interfaceDofs = schur.interfaceDofs();
        setup.interiorDofs = schur.interiorDofs();
        setup.coarseDimension = schur.coarseDimension();
        setup.interfacePatchCount = schur.interfacePatchCount();
        setup.energyCandidateModes = schur.energyCandidateModes();
        setup.energyEigenSetupSeconds = schur.energyEigenSetupSeconds();
        setup.energySelectedEigenvalueMin = schur.energySelectedEigenvalueMin();
        setup.energySelectedEigenvalueMax = schur.energySelectedEigenvalueMax();
        setup.globalSlowCandidateDimension = schur.globalSlowCandidateDimension();
        setup.globalSlowSetupSeconds = schur.globalSlowSetupSeconds();
        setup.globalSlowEstimatedLambdaMax = schur.globalSlowEstimatedLambdaMax();
        setup.globalSlowSelectedRitzMin = schur.globalSlowSelectedRitzMin();
        setup.globalSlowSelectedRitzMax = schur.globalSlowSelectedRitzMax();
        setup.proxyGraphEdges = schur.proxyGraphEdges();
        setup.proxyProbeColumns = schur.proxyProbeColumns();
        setup.proxySchurApplies = schur.proxySchurApplies();
        setup.proxyDiagnosticsSeconds = schur.proxyDiagnosticsSeconds();
        setup.proxyRing3Coverage = schur.proxyRing3Coverage();
        setup.proxyRing3OperatorError = schur.proxyRing3OperatorError();
        setup.proxyRing3EstimatedNnz = schur.proxyRing3EstimatedNnz();
        setup.proxyRing3MemoryEstimateBytes = schur.proxyRing3MemoryEstimateBytes();
        setup.proxyRecommended = schur.proxyRecommended();
        setup.proxyNnz = schur.proxyNnz();
        setup.proxyDensity = schur.proxyDensity();
        setup.proxyColors = schur.proxyColors();
        setup.proxyProbingSchurApplies = schur.proxyProbingSchurApplies();
        setup.proxyProbingBlockSize = schur.proxyProbingBlockSize();
        setup.proxyProbingBlockCalls = schur.proxyProbingBlockCalls();
        setup.proxyValidationSchurApplies = schur.proxyValidationSchurApplies();
        setup.proxySymbolicCalls = schur.proxySymbolicCalls();
        setup.proxyNumericalCalls = schur.proxyNumericalCalls();
        setup.proxySetupSeconds = schur.proxySetupSeconds();
        setup.proxySymbolicSeconds = schur.proxySymbolicSeconds();
        setup.proxyNumericalSeconds = schur.proxyNumericalSeconds();
        setup.proxySymmetryError = schur.proxySymmetryError();
        setup.proxyMinimumTestRayleigh = schur.proxyMinimumTestRayleigh();
        setup.proxyDiagonalShift = schur.proxyDiagonalShift();
        setup.proxyDiagonalCompensation = schur.proxyDiagonalCompensation();
        setup.proxyValueHash = schur.proxyValueHash();
        setup.proxyBlockMaximumDifference = schur.proxyBlockMaximumDifference();
        setup.proxyBlockRelativeDifference = schur.proxyBlockRelativeDifference();
        setup.proxyMemoryBytes = schur.proxyMemoryBytes();
        setup.proxyMatrixCacheHit = schur.proxyMatrixCacheHit();
        setup.proxyFactorCacheHit = schur.proxyFactorCacheHit();
        setup.localFactorizationSeconds = schur.localFactorizationSeconds();
        setup.localSymbolicAnalysisSeconds = schur.localSymbolicAnalysisSeconds();
        setup.localNumericalFactorizationSeconds = schur.localNumericalFactorizationSeconds();
        setup.localSymbolicAnalysisCalls = schur.localSymbolicAnalysisCalls();
        setup.localNumericalFactorizationCalls = schur.localNumericalFactorizationCalls();
        setup.memoryBytes = schur.memoryBytes();
        setup.subdomainPerformance = schur.subdomainPerformance();
        setup.status = "ready";
    }
};

DdmSchurSolver::DdmSchurSolver(const Mesh& mesh,
                               const SparseMatrix& system,
                               const CaseConfig& physics,
                               Options options)
{
    if (options.maxIterations <= 0 || options.restart <= 0 || options.relativeTolerance <= 0.0) {
        throw std::runtime_error("[Schur] Invalid FGMRES options.");
    }
    if (options.localSolveThreads <= 0 || options.localPardisoThreads <= 0) {
        throw std::runtime_error("[Schur] Local OpenMP/PARDISO thread counts must be positive.");
    }
    if (options.proxyCacheEnabled && options.proxyCachePath.empty()) {
        throw std::runtime_error("[Schur proxy] Cache is enabled but no cache path was provided.");
    }
    if (options.coarseGlobalQuadraticZ && !options.coarseLinearZ) {
        throw std::runtime_error("[Schur] Global quadratic-z enrichment requires linear-z coarse modes.");
    }
    if (options.coarseInterfacePatches && options.coarseGlobalQuadraticZ) {
        throw std::runtime_error(
            "[Schur] Interface-patch coarse space cannot use global quadratic-z enrichment.");
    }
    if (options.coarseEnergyAdaptive
        && (options.coarseInterfacePatches || options.coarseGlobalQuadraticZ)) {
        throw std::runtime_error(
            "[Schur] Energy-adaptive coarse space cannot be combined with patch or quadratic-z modes.");
    }
    if (options.coarseGlobalSlow
        && (options.coarseEnergyAdaptive
            || options.coarseInterfacePatches
            || options.coarseGlobalQuadraticZ)) {
        throw std::runtime_error(
            "[Schur] Global slow coarse space cannot be combined with local spectral, patch, or quadratic-z modes.");
    }
    if (options.coarseGlobalSlow
        && (options.globalSlowModes <= 0
            || options.globalSlowSubspaceDimension < options.globalSlowModes)) {
        throw std::runtime_error("[Schur] Invalid global slow-mode options.");
    }
    if (options.coarseEnergyAdaptive
        && (options.energyMaxModesPerDomain <= 0
            || options.energySubspaceIterations <= 0
            || options.energyEigenvalueThreshold <= 0.0)) {
        throw std::runtime_error("[Schur] Invalid energy-adaptive coarse-space options.");
    }
    if (options.proxyDiagnostics
        && (!(options.proxyHighConductivityThreshold > 0.0)
            || options.proxyProbeColumns <= 0
            || options.proxyRing < 1
            || options.proxyRing > 3)) {
        throw std::runtime_error("[Schur proxy] Invalid diagnostic options.");
    }
    if (options.proxyEnabled && (options.proxyRing != 1
        || !(options.proxyHighConductivityThreshold > 0.0)
        || options.proxyBlockSize <= 0)) {
        throw std::runtime_error("[Schur proxy] Invalid or locality-unapproved proxy options.");
    }
    if (options.proxyDisableCoarse && !options.proxyEnabled) {
        throw std::runtime_error("[Schur proxy] Coarse disabling is only valid for proxy-only runs.");
    }
    const auto setupStart = std::chrono::steady_clock::now();
    impl_ = std::make_unique<Impl>(mesh, system, physics, options);
    impl_->setup.setupSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - setupStart).count();

    std::cout << std::setprecision(6)
              << "[Schur] domains: " << impl_->setup.domains << '\n'
              << "[Schur] total dofs: " << impl_->setup.totalDofs << '\n'
              << "[Schur] interface dofs: " << impl_->setup.interfaceDofs << '\n'
              << "[Schur] interior dofs: " << impl_->setup.interiorDofs << '\n'
              << "[Schur] coarse dimension: " << impl_->setup.coarseDimension << '\n'
              << "[Schur] interface boundary-side patches: "
              << impl_->setup.interfacePatchCount << '\n'
              << "[Schur] energy candidate modes: "
              << impl_->setup.energyCandidateModes << '\n'
              << "[Schur] energy selected eigenvalue range: ["
              << impl_->setup.energySelectedEigenvalueMin << ", "
              << impl_->setup.energySelectedEigenvalueMax << "]\n"
              << "[Schur] energy eigen setup time: "
              << impl_->setup.energyEigenSetupSeconds << " s\n"
              << "[Schur] global slow candidate dimension: "
              << impl_->setup.globalSlowCandidateDimension << '\n'
              << "[Schur] global slow estimated lambda max: "
              << impl_->setup.globalSlowEstimatedLambdaMax << '\n'
              << "[Schur] global slow selected Ritz range: ["
              << impl_->setup.globalSlowSelectedRitzMin << ", "
              << impl_->setup.globalSlowSelectedRitzMax << "]\n"
              << "[Schur] global slow setup time: "
              << impl_->setup.globalSlowSetupSeconds << " s\n"
              << "[Schur proxy] diagnostic probe columns: "
              << impl_->setup.proxyProbeColumns << '\n'
              << "[Schur proxy] exact diagnostic applies: "
              << impl_->setup.proxySchurApplies << '\n'
              << "[Schur proxy] ring-3 coverage/error: "
              << impl_->setup.proxyRing3Coverage << " / "
              << impl_->setup.proxyRing3OperatorError << '\n'
              << "[Schur proxy] diagnostic time: "
              << impl_->setup.proxyDiagnosticsSeconds << " s\n"
              << "[Schur proxy] nnz/colors/probes: "
              << impl_->setup.proxyNnz << " / "
              << impl_->setup.proxyColors << " / "
              << impl_->setup.proxyProbingSchurApplies << '\n'
              << "[Schur proxy] probing block size/calls: "
              << impl_->setup.proxyProbingBlockSize << " / "
              << impl_->setup.proxyProbingBlockCalls << '\n'
              << "[Schur proxy] scalar validation applies/max/relative difference: "
              << impl_->setup.proxyValidationSchurApplies << " / "
              << impl_->setup.proxyBlockMaximumDifference << " / "
              << impl_->setup.proxyBlockRelativeDifference << '\n'
              << "[Schur proxy] phase 11/22: "
              << impl_->setup.proxySymbolicCalls << " / "
              << impl_->setup.proxyNumericalCalls << '\n'
              << "[Schur proxy] symbolic/numerical time: "
              << impl_->setup.proxySymbolicSeconds << " / "
              << impl_->setup.proxyNumericalSeconds << " s\n"
              << "[Schur] local symbolic analyses: "
              << impl_->setup.localSymbolicAnalysisCalls << '\n'
              << "[Schur] local numerical factorizations: "
              << impl_->setup.localNumericalFactorizationCalls << '\n'
              << "[Schur] symbolic analysis time: "
              << impl_->setup.localSymbolicAnalysisSeconds << " s\n"
              << "[Schur] numerical factorization time: "
              << impl_->setup.localNumericalFactorizationSeconds << " s\n";
}

DdmSchurSolver::~DdmSchurSolver() = default;
DdmSchurSolver::DdmSchurSolver(DdmSchurSolver&&) noexcept = default;
DdmSchurSolver& DdmSchurSolver::operator=(DdmSchurSolver&&) noexcept = default;

const Report& DdmSchurSolver::setupReport() const { return impl_->setup; }

int DdmSchurSolver::interfaceDofs() const { return impl_->schur.interfaceDofs(); }

const std::vector<int>& DdmSchurSolver::interfaceGlobalDofs() const
{
    return impl_->schur.interfaceGlobalDofs();
}

std::vector<double> DdmSchurSolver::condensedRhs(const std::vector<double>& globalRhs)
{
    return impl_->schur.condensedRhs(globalRhs);
}

void DdmSchurSolver::applyExactSchur(const std::vector<double>& interfaceVector,
                                     std::vector<double>& result)
{
    impl_->schur.apply(interfaceVector, result);
}

std::vector<double> DdmSchurSolver::recover(
    const std::vector<double>& globalRhs,
    const std::vector<double>& interfaceSolution)
{
    return impl_->schur.recover(globalRhs, interfaceSolution);
}

SolveResult DdmSchurSolver::solve(const std::vector<double>& globalRhs)
{
    SolveResult result;
    const auto totalStart = std::chrono::steady_clock::now();
    const std::vector<SubdomainPerformance> domainPerformanceBefore =
        impl_->schur.subdomainPerformance();
    const int proxySolvesBefore = impl_->schur.proxySolveCalls();
    const double proxySecondsBefore = impl_->schur.proxySolveSeconds();

    const int callsBeforeCondensation = impl_->schur.localSolveCalls();
    const double localSecondsBeforeCondensation = impl_->schur.localSolveSeconds();
    const auto rhsStart = std::chrono::steady_clock::now();
    const std::vector<double> condensed = impl_->schur.condensedRhs(globalRhs);
    const double condensedRhsSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - rhsStart).count();
    const int condensedSolveCalls =
        impl_->schur.localSolveCalls() - callsBeforeCondensation;
    const double condensedLocalSeconds =
        impl_->schur.localSolveSeconds() - localSecondsBeforeCondensation;

    const auto globalSetupStart = std::chrono::steady_clock::now();
    const bool preparedGlobalSlow = impl_->schur.prepareGlobalSlow(condensed);
    const double deferredGlobalSetupSeconds = preparedGlobalSlow
        ? std::chrono::duration<double>(
            std::chrono::steady_clock::now() - globalSetupStart).count()
        : 0.0;
    if (preparedGlobalSlow) {
        impl_->setup.setupSeconds += deferredGlobalSetupSeconds;
        impl_->setup.coarseDimension = impl_->schur.coarseDimension();
        impl_->setup.globalSlowCandidateDimension = impl_->schur.globalSlowCandidateDimension();
        impl_->setup.globalSlowSetupSeconds = impl_->schur.globalSlowSetupSeconds();
        impl_->setup.globalSlowEstimatedLambdaMax = impl_->schur.globalSlowEstimatedLambdaMax();
        impl_->setup.globalSlowSelectedRitzMin = impl_->schur.globalSlowSelectedRitzMin();
        impl_->setup.globalSlowSelectedRitzMax = impl_->schur.globalSlowSelectedRitzMax();
        impl_->setup.memoryBytes = impl_->schur.memoryBytes();
        impl_->schur.resetRuntimeCounters();
        std::cout << "[Schur] RHS-seeded global harmonic coarse dimension: "
                  << impl_->setup.coarseDimension << '\n'
                  << "[Schur] RHS-seeded global harmonic Ritz range: ["
                  << impl_->setup.globalSlowSelectedRitzMin << ", "
                  << impl_->setup.globalSlowSelectedRitzMax << "]\n"
                  << "[Schur] RHS-seeded global harmonic setup time: "
                  << impl_->setup.globalSlowSetupSeconds << " s\n";
    }
    result.report = impl_->setup;
    result.report.condensedRhsSeconds = condensedRhsSeconds;

    const int solvesBeforeInterface = impl_->schur.localSolveCalls();
    const int matvecsBeforeInterface = impl_->schur.matvecCalls();
    const double localSecondsBeforeInterface = impl_->schur.localSolveSeconds();
    const double coarseSecondsBeforeInterface = impl_->schur.coarseSolveSeconds();
    const double schurApplySecondsBeforeInterface = impl_->schur.schurApplySeconds();

    const auto interfaceStart = std::chrono::steady_clock::now();
    FgmresResult interfaceResult = solveFgmres(impl_->schur, condensed, impl_->options);
    result.report.interfaceSolveSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - interfaceStart).count();
    result.report.fgmresSeconds = result.report.interfaceSolveSeconds;
    result.report.iterations = interfaceResult.iterations;
    result.report.interfaceRelativeResidual = interfaceResult.relativeResidual;

    const auto recoveryStart = std::chrono::steady_clock::now();
    result.interfaceSolution = interfaceResult.solution;
    result.temperature = impl_->schur.recover(globalRhs, result.interfaceSolution);
    result.report.recoverySeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - recoveryStart).count();
    result.report.totalSolveSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - totalStart).count()
        - deferredGlobalSetupSeconds;
    result.report.totalSeconds = result.report.setupSeconds + result.report.totalSolveSeconds;
    result.report.localSolveCalls = condensedSolveCalls
        + impl_->schur.localSolveCalls() - solvesBeforeInterface;
    result.report.schurMatvecs = impl_->schur.matvecCalls() - matvecsBeforeInterface;
    result.report.localSolveSeconds = condensedLocalSeconds
        + impl_->schur.localSolveSeconds() - localSecondsBeforeInterface;
    result.report.schurApplySeconds = impl_->schur.schurApplySeconds()
        - schurApplySecondsBeforeInterface;
    result.report.subdomainPerformance = impl_->schur.subdomainPerformance();
    for (std::size_t slot = 0; slot < result.report.subdomainPerformance.size(); ++slot) {
        result.report.subdomainPerformance[slot].phase33Calls -=
            domainPerformanceBefore[slot].phase33Calls;
        result.report.subdomainPerformance[slot].phase33Seconds -=
            domainPerformanceBefore[slot].phase33Seconds;
    }
    result.report.coarseSolveSeconds =
        impl_->schur.coarseSolveSeconds() - coarseSecondsBeforeInterface;
    result.report.proxySolveCalls = impl_->schur.proxySolveCalls() - proxySolvesBefore;
    result.report.proxySolveSeconds = impl_->schur.proxySolveSeconds() - proxySecondsBefore;
    result.report.status = interfaceResult.converged ? "success" : "max_iterations";

    std::cout << std::setprecision(12)
              << "[Schur] FGMRES iterations: " << result.report.iterations << '\n'
              << "[Schur] residual: " << result.report.interfaceRelativeResidual << '\n'
              << "[Schur] local solve time: " << result.report.localSolveSeconds << " s\n"
              << "[Schur] coarse solve time: " << result.report.coarseSolveSeconds << " s\n"
              << "[Schur proxy] phase-33 calls/time: "
              << result.report.proxySolveCalls << " / "
              << result.report.proxySolveSeconds << " s\n"
              << "[Schur] total FGMRES time: " << result.report.fgmresSeconds << " s\n"
              << "[Schur] total solve time: " << result.report.totalSolveSeconds << " s\n"
              << "[Schur] total time (setup + solve): " << result.report.totalSeconds << " s\n";
    return result;
}

SolveResult DdmSchurSolver::solveCorrection(
    const std::vector<double>& globalRhs,
    const std::vector<double>& interfaceInitialGuess)
{
    if (static_cast<int>(interfaceInitialGuess.size()) != impl_->schur.interfaceDofs()) {
        throw std::runtime_error("[Schur MOR] Interface initial guess has the wrong size.");
    }
    SolveResult result;
    const auto totalStart = std::chrono::steady_clock::now();
    const int proxySolvesBefore = impl_->schur.proxySolveCalls();
    const double proxySecondsBefore = impl_->schur.proxySolveSeconds();
    const int solvesBefore = impl_->schur.localSolveCalls();
    const int matvecsBefore = impl_->schur.matvecCalls();
    const double localBefore = impl_->schur.localSolveSeconds();
    const double coarseBefore = impl_->schur.coarseSolveSeconds();

    const auto rhsStart = std::chrono::steady_clock::now();
    const std::vector<double> condensed = impl_->schur.condensedRhs(globalRhs);
    std::vector<double> initialImage;
    impl_->schur.apply(interfaceInitialGuess, initialImage);
    std::vector<double> correctionRhs(condensed.size(), 0.0);
    for (std::size_t i = 0; i < condensed.size(); ++i) {
        correctionRhs[i] = condensed[i] - initialImage[i];
    }
    const double condensedSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - rhsStart).count();

    const auto interfaceStart = std::chrono::steady_clock::now();
    // The validated Stage-1 FGMRES implementation and its true-residual gate
    // are called verbatim, now on the exact correction equation.
    const double condensedNorm = std::max(1.0e-300, normVector(condensed));
    const double correctionNorm = normVector(correctionRhs);
    FgmresResult correction;
    correction.solution.assign(correctionRhs.size(), 0.0);
    if (correctionNorm / condensedNorm <= impl_->options.relativeTolerance) {
        correction.converged = true;
        correction.relativeResidual = correctionNorm / condensedNorm;
    } else {
        Options correctionOptions = impl_->options;
        // solveFgmres keeps its validated gate ||r_k||/||r_0||.  Scaling the
        // requested tolerance makes that unchanged gate exactly equivalent to
        // the conventional initial-guess criterion ||r_k||/||g|| <= tol.
        correctionOptions.relativeTolerance = std::min(0.99,
            impl_->options.relativeTolerance * condensedNorm / correctionNorm);
        correction = solveFgmres(impl_->schur, correctionRhs, correctionOptions);
    }
    result.interfaceSolution = interfaceInitialGuess;
    for (std::size_t i = 0; i < result.interfaceSolution.size(); ++i) {
        result.interfaceSolution[i] += correction.solution[i];
    }
    const double interfaceSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - interfaceStart).count();

    std::vector<double> finalImage;
    impl_->schur.apply(result.interfaceSolution, finalImage);
    std::vector<double> finalResidual(condensed.size(), 0.0);
    for (std::size_t i = 0; i < condensed.size(); ++i) {
        finalResidual[i] = condensed[i] - finalImage[i];
    }
    const double rhsNorm = condensedNorm;

    const auto recoveryStart = std::chrono::steady_clock::now();
    result.temperature = impl_->schur.recover(globalRhs, result.interfaceSolution);
    const double recoverySeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - recoveryStart).count();

    result.report = impl_->setup;
    result.report.condensedRhsSeconds = condensedSeconds;
    result.report.interfaceSolveSeconds = interfaceSeconds;
    result.report.fgmresSeconds = interfaceSeconds;
    result.report.recoverySeconds = recoverySeconds;
    result.report.totalSolveSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - totalStart).count();
    result.report.totalSeconds = result.report.setupSeconds + result.report.totalSolveSeconds;
    result.report.iterations = correction.iterations;
    result.report.interfaceRelativeResidual = normVector(finalResidual) / rhsNorm;
    result.report.localSolveCalls = impl_->schur.localSolveCalls() - solvesBefore;
    result.report.schurMatvecs = impl_->schur.matvecCalls() - matvecsBefore;
    result.report.localSolveSeconds = impl_->schur.localSolveSeconds() - localBefore;
    result.report.coarseSolveSeconds = impl_->schur.coarseSolveSeconds() - coarseBefore;
    result.report.proxySolveCalls = impl_->schur.proxySolveCalls() - proxySolvesBefore;
    result.report.proxySolveSeconds = impl_->schur.proxySolveSeconds() - proxySecondsBefore;
    result.report.status = correction.converged ? "success" : "max_iterations";
    return result;
}

} // namespace ddm_schur
