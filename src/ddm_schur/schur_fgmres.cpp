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

struct DdmSchurSolver::Impl {
    Options options;
    SchurOperator schur;
    Report setup;

    Impl(const Mesh& mesh, const SparseMatrix& system, Options requested)
        : options(requested),
          schur(mesh,
                system,
                requested.coarseLinearXY,
                requested.coarseLinearZ,
                requested.coarseGlobalQuadraticZ,
                requested.coarseInterfacePatches,
                requested.coarseInterfaceLinearXY)
    {
        setup.domains = schur.domains();
        setup.totalDofs = schur.totalDofs();
        setup.interfaceDofs = schur.interfaceDofs();
        setup.interiorDofs = schur.interiorDofs();
        setup.coarseDimension = schur.coarseDimension();
        setup.interfacePatchCount = schur.interfacePatchCount();
        setup.localFactorizationSeconds = schur.localFactorizationSeconds();
        setup.localSymbolicAnalysisSeconds = schur.localSymbolicAnalysisSeconds();
        setup.localNumericalFactorizationSeconds = schur.localNumericalFactorizationSeconds();
        setup.localSymbolicAnalysisCalls = schur.localSymbolicAnalysisCalls();
        setup.localNumericalFactorizationCalls = schur.localNumericalFactorizationCalls();
        setup.memoryBytes = schur.memoryBytes();
        setup.status = "ready";
    }
};

DdmSchurSolver::DdmSchurSolver(const Mesh& mesh, const SparseMatrix& system, Options options)
{
    if (options.maxIterations <= 0 || options.restart <= 0 || options.relativeTolerance <= 0.0) {
        throw std::runtime_error("[Schur] Invalid FGMRES options.");
    }
    if (options.coarseGlobalQuadraticZ && !options.coarseLinearZ) {
        throw std::runtime_error("[Schur] Global quadratic-z enrichment requires linear-z coarse modes.");
    }
    if (options.coarseInterfacePatches && options.coarseGlobalQuadraticZ) {
        throw std::runtime_error(
            "[Schur] Interface-patch coarse space cannot use global quadratic-z enrichment.");
    }
    const auto setupStart = std::chrono::steady_clock::now();
    impl_ = std::make_unique<Impl>(mesh, system, options);
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

SolveResult DdmSchurSolver::solve(const std::vector<double>& globalRhs)
{
    SolveResult result;
    result.report = impl_->setup;
    const int solvesBefore = impl_->schur.localSolveCalls();
    const int matvecsBefore = impl_->schur.matvecCalls();
    const double localSecondsBefore = impl_->schur.localSolveSeconds();
    const double coarseSecondsBefore = impl_->schur.coarseSolveSeconds();
    const auto totalStart = std::chrono::steady_clock::now();

    const auto rhsStart = std::chrono::steady_clock::now();
    const std::vector<double> condensed = impl_->schur.condensedRhs(globalRhs);
    result.report.condensedRhsSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - rhsStart).count();

    const auto interfaceStart = std::chrono::steady_clock::now();
    FgmresResult interfaceResult = solveFgmres(impl_->schur, condensed, impl_->options);
    result.report.interfaceSolveSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - interfaceStart).count();
    result.report.fgmresSeconds = result.report.interfaceSolveSeconds;
    result.report.iterations = interfaceResult.iterations;
    result.report.interfaceRelativeResidual = interfaceResult.relativeResidual;

    const auto recoveryStart = std::chrono::steady_clock::now();
    result.temperature = impl_->schur.recover(globalRhs, interfaceResult.solution);
    result.report.recoverySeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - recoveryStart).count();
    result.report.totalSolveSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - totalStart).count();
    result.report.totalSeconds = result.report.setupSeconds + result.report.totalSolveSeconds;
    result.report.localSolveCalls = impl_->schur.localSolveCalls() - solvesBefore;
    result.report.schurMatvecs = impl_->schur.matvecCalls() - matvecsBefore;
    result.report.localSolveSeconds = impl_->schur.localSolveSeconds() - localSecondsBefore;
    result.report.coarseSolveSeconds = impl_->schur.coarseSolveSeconds() - coarseSecondsBefore;
    result.report.status = interfaceResult.converged ? "success" : "max_iterations";

    std::cout << std::setprecision(12)
              << "[Schur] FGMRES iterations: " << result.report.iterations << '\n'
              << "[Schur] residual: " << result.report.interfaceRelativeResidual << '\n'
              << "[Schur] local solve time: " << result.report.localSolveSeconds << " s\n"
              << "[Schur] coarse solve time: " << result.report.coarseSolveSeconds << " s\n"
              << "[Schur] total FGMRES time: " << result.report.fgmresSeconds << " s\n"
              << "[Schur] total solve time: " << result.report.totalSolveSeconds << " s\n"
              << "[Schur] total time (setup + solve): " << result.report.totalSeconds << " s\n";
    return result;
}

} // namespace ddm_schur
