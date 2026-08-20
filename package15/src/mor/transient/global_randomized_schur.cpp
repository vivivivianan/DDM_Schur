#include "global_randomized_schur.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

namespace mor::transient {
namespace {

using Clock = std::chrono::steady_clock;

double secondsSince(const Clock::time_point& start)
{
    return std::chrono::duration<double>(
        Clock::now() - start).count();
}

double dot(const std::vector<double>& left,
           const std::vector<double>& right)
{
    if (left.size() != right.size()) {
        throw std::runtime_error(
            "[Global randomized] Dot-product size mismatch.");
    }
    long double value = 0.0L;
    for (std::size_t row = 0; row < left.size(); ++row) {
        value += static_cast<long double>(left[row]) * right[row];
    }
    return static_cast<double>(value);
}

double norm(const std::vector<double>& values)
{
    return std::sqrt(std::max(0.0, dot(values, values)));
}

double weightedDot(const std::vector<double>& left,
                   const std::vector<double>& right,
                   const std::vector<double>& metric)
{
    if (left.size() != right.size()
        || left.size() != metric.size()) {
        throw std::runtime_error(
            "[Global randomized] Weighted dot-product size mismatch.");
    }
    long double value = 0.0L;
    for (std::size_t row = 0; row < left.size(); ++row) {
        value += static_cast<long double>(left[row])
            * metric[row] * right[row];
    }
    return static_cast<double>(value);
}

double weightedDot(const double* left,
                   const double* right,
                   const std::vector<double>& metric,
                   int rows)
{
    long double value = 0.0L;
    for (int row = 0; row < rows; ++row) {
        value += static_cast<long double>(left[row])
            * metric[static_cast<std::size_t>(row)] * right[row];
    }
    return static_cast<double>(value);
}

struct SolveResult {
    std::vector<double> solution;
    int iterations = 0;
    double residual = std::numeric_limits<double>::infinity();
    bool converged = false;
    bool timedOut = false;
};

SolveResult solveGlobalSchurPcg(
    const std::vector<double>& rightHandSide,
    LocalPortReducedSchurSolver& localPreconditioner,
    const ReducedDynamicSchurOperator& schur,
    int maximumIterations,
    double tolerance,
    const Clock::time_point& deadline,
    GlobalRandomizedSchurDiagnostics& diagnostics)
{
    const int rows = schur.size();
    SolveResult result;
    result.solution.assign(static_cast<std::size_t>(rows), 0.0);
    const double rightHandSideNorm = norm(rightHandSide);
    if (!(rightHandSideNorm > 0.0)) {
        result.residual = 0.0;
        result.converged = true;
        return result;
    }

    const std::vector<double>& diagonal = schur.diagonal();
    double diagonalScale = 0.0;
    for (double value : diagonal) {
        diagonalScale = std::max(diagonalScale, std::abs(value));
    }
    const double diagonalFloor =
        1024.0 * std::numeric_limits<double>::epsilon()
        * std::max(std::numeric_limits<double>::min(),
                   diagonalScale);

    const auto applyPreconditioner =
        [&](const std::vector<double>& residual,
            std::vector<double>& output) {
            // Symmetric two-level preconditioner.  The frozen M8.9 local
            // space is used only to accelerate the global solve; it is not
            // inserted into the global-only randomized trial space.
            std::vector<double> local;
            localPreconditioner.localGalerkinInterfaceResponse(
                residual, local);
            std::vector<double> localImage;
            schur.apply(local, localImage);
            ++diagnostics.globalSchurApplyCount;
            std::vector<double> dualComplement(residual.size(), 0.0);
            for (std::size_t row = 0; row < residual.size(); ++row) {
                dualComplement[row] = residual[row] - localImage[row];
            }
            std::vector<double> jacobi(residual.size(), 0.0);
            for (std::size_t row = 0; row < residual.size(); ++row) {
                if (!(std::abs(diagonal[row]) > diagonalFloor)
                    || !std::isfinite(diagonal[row])) {
                    throw std::runtime_error(
                        "[Global randomized] Exact Schur Jacobi "
                        "preconditioner has a tiny/nonfinite diagonal.");
                }
                jacobi[row] = dualComplement[row] / diagonal[row];
            }
            std::vector<double> primalComplement;
            localPreconditioner.projectSchurEnergyComplement(
                jacobi, primalComplement);
            ++diagnostics.globalSchurApplyCount;
            output.resize(residual.size());
            for (std::size_t row = 0; row < residual.size(); ++row) {
                output[row] = local[row] + primalComplement[row];
            }
            ++diagnostics.preconditionerApplyCount;
        };

    std::vector<double> residual = rightHandSide;
    std::vector<double> preconditioned;
    applyPreconditioner(residual, preconditioned);
    std::vector<double> direction = preconditioned;
    double residualPreconditioned =
        dot(residual, preconditioned);
    if (!(residualPreconditioned > 0.0)
        || !std::isfinite(residualPreconditioned)) {
        return result;
    }
    for (int iteration = 0;
         iteration < maximumIterations; ++iteration) {
        if (Clock::now() >= deadline) {
            result.timedOut = true;
            break;
        }
        std::vector<double> image;
        schur.apply(direction, image);
        ++diagnostics.globalSchurApplyCount;
        const double denominator = dot(direction, image);
        if (!(denominator > 0.0) || !std::isfinite(denominator)) {
            break;
        }
        const double alpha = residualPreconditioned / denominator;
        for (int row = 0; row < rows; ++row) {
            result.solution[static_cast<std::size_t>(row)]
                += alpha * direction[static_cast<std::size_t>(row)];
            residual[static_cast<std::size_t>(row)]
                -= alpha * image[static_cast<std::size_t>(row)];
        }
        result.iterations = iteration + 1;
        result.residual = norm(residual) / rightHandSideNorm;
        if (result.residual <= tolerance) {
            result.converged = true;
            break;
        }
        if (Clock::now() >= deadline) {
            result.timedOut = true;
            break;
        }
        std::vector<double> nextPreconditioned;
        applyPreconditioner(residual, nextPreconditioned);
        const double nextResidualPreconditioned =
            dot(residual, nextPreconditioned);
        if (!(nextResidualPreconditioned > 0.0)
            || !std::isfinite(nextResidualPreconditioned)) {
            break;
        }
        const double beta =
            nextResidualPreconditioned / residualPreconditioned;
        for (int row = 0; row < rows; ++row) {
            direction[static_cast<std::size_t>(row)] =
                nextPreconditioned[static_cast<std::size_t>(row)]
                + beta * direction[static_cast<std::size_t>(row)];
        }
        residualPreconditioned = nextResidualPreconditioned;
    }

    // The gate always uses a newly evaluated, unpreconditioned Schur
    // residual rather than the PCG recurrence residual.
    std::vector<double> finalImage;
    schur.apply(result.solution, finalImage);
    ++diagnostics.globalSchurApplyCount;
    for (std::size_t row = 0; row < finalImage.size(); ++row) {
        finalImage[row] -= rightHandSide[row];
    }
    result.residual = norm(finalImage) / rightHandSideNorm;
    result.converged = result.residual <= tolerance;
    return result;
}

bool appendWeightedColumn(
    std::vector<double>& basis,
    std::vector<double> candidate,
    const std::vector<double>& metric,
    int rows,
    double relativeTolerance)
{
    const double original =
        std::sqrt(std::max(0.0,
            weightedDot(candidate, candidate, metric)));
    if (!(original > 0.0) || !std::isfinite(original)) {
        return false;
    }
    const int rank = static_cast<int>(
        basis.size() / static_cast<std::size_t>(rows));
#ifdef USE_MKL_PARDISO
    std::vector<double> weightedCandidate(
        static_cast<std::size_t>(rows), 0.0);
    std::vector<double> coefficients(
        static_cast<std::size_t>(rank), 0.0);
    for (int pass = 0; pass < 2; ++pass) {
        for (int row = 0; row < rows; ++row) {
            weightedCandidate[static_cast<std::size_t>(row)] =
                metric[static_cast<std::size_t>(row)]
                * candidate[static_cast<std::size_t>(row)];
        }
        if (rank > 0) {
            cblas_dgemv(
                CblasColMajor, CblasTrans, rows, rank,
                1.0, basis.data(), rows, weightedCandidate.data(), 1,
                0.0, coefficients.data(), 1);
            cblas_dgemv(
                CblasColMajor, CblasNoTrans, rows, rank,
                -1.0, basis.data(), rows, coefficients.data(), 1,
                1.0, candidate.data(), 1);
        }
    }
#else
    for (int pass = 0; pass < 2; ++pass) {
        for (int mode = 0; mode < rank; ++mode) {
            const double* basisColumn = basis.data()
                + static_cast<std::size_t>(mode * rows);
            const double coefficient =
                weightedDot(
                    basisColumn, candidate.data(), metric, rows);
            for (int row = 0; row < rows; ++row) {
                candidate[static_cast<std::size_t>(row)]
                    -= coefficient
                    * basisColumn[row];
            }
        }
    }
#endif
    const double remaining =
        std::sqrt(std::max(0.0,
            weightedDot(candidate, candidate, metric)));
    if (!(remaining > relativeTolerance * original)
        || !std::isfinite(remaining)) {
        return false;
    }
    for (double& value : candidate) value /= remaining;
    basis.insert(basis.end(), candidate.begin(), candidate.end());
    return true;
}

double weightedOrthogonalityError(
    const std::vector<double>& basis,
    const std::vector<double>& metric,
    int rows,
    int rank)
{
#ifdef USE_MKL_PARDISO
    if (rank == 0) return 0.0;
    std::vector<double> weightedBasis(basis.size(), 0.0);
    for (int mode = 0; mode < rank; ++mode) {
        const double* source = basis.data()
            + static_cast<std::size_t>(mode * rows);
        double* target = weightedBasis.data()
            + static_cast<std::size_t>(mode * rows);
        for (int row = 0; row < rows; ++row) {
            target[row] =
                metric[static_cast<std::size_t>(row)] * source[row];
        }
    }
    std::vector<double> gram(
        static_cast<std::size_t>(rank * rank), 0.0);
    cblas_dgemm(
        CblasColMajor, CblasTrans, CblasNoTrans,
        rank, rank, rows, 1.0, basis.data(), rows,
        weightedBasis.data(), rows, 0.0, gram.data(), rank);
    double maximum = 0.0;
    for (int left = 0; left < rank; ++left) {
        for (int right = 0; right < rank; ++right) {
            maximum = std::max(maximum, std::abs(
                gram[static_cast<std::size_t>(left + right * rank)]
                - (left == right ? 1.0 : 0.0)));
        }
    }
    return maximum;
#else
    double maximum = 0.0;
    for (int left = 0; left < rank; ++left) {
        const double* leftColumn = basis.data()
            + static_cast<std::size_t>(left * rows);
        for (int right = 0; right < rank; ++right) {
            const double* rightColumn = basis.data()
                + static_cast<std::size_t>(right * rows);
            const double gram =
                weightedDot(leftColumn, rightColumn, metric, rows);
            maximum = std::max(maximum, std::abs(
                gram - (left == right ? 1.0 : 0.0)));
        }
    }
    return maximum;
#endif
}

} // namespace

GlobalRandomizedSchurResult buildGlobalRandomizedSchurPortSpace(
    const LocalPortModel& localModel,
    LocalPortReducedSchurSolver& localPreconditioner,
    const ReducedDynamicSchurOperator& exactSchur,
    const std::vector<double>& interfaceMetricDiagonal,
    const GlobalRandomizedSchurOptions& options)
{
    const auto start = Clock::now();
    GlobalRandomizedSchurResult result;
    auto& model = result.model;
    auto& diagnostics = result.diagnostics;
    diagnostics.globalInterfaceDofs = exactSchur.size();
    diagnostics.requestedRank = options.requestedRank;
    diagnostics.seed = options.seed;
    diagnostics.composition = options.composition;
    diagnostics.baselineWorkingSetBytes =
        globalCoarseCurrentWorkingSetBytes();
    diagnostics.peakWorkingSetBytes =
        diagnostics.baselineWorkingSetBytes;

    if (options.requestedRank <= 0
        || options.requestedRank > exactSchur.size()
        || options.innerMaximumIterations <= 0
        || !(options.innerTolerance > 0.0)
        || !(options.orthogonalityTolerance > 0.0)
        || !(options.deflationTolerance > 0.0)
        || !(options.maximumBasisTimeSeconds > 0.0)
        || (options.composition != "global-only"
            && options.composition != "augment-local")
        || localModel.fullInterfaceDofs != exactSchur.size()
        || interfaceMetricDiagonal.size()
            != static_cast<std::size_t>(exactSchur.size())) {
        throw std::runtime_error(
            "[Global randomized] Invalid options or dimensions.");
    }

    std::vector<double> metric = interfaceMetricDiagonal;
    double maximumMetric = 0.0;
    for (double value : metric) {
        if (value < 0.0 || !std::isfinite(value)) {
            throw std::runtime_error(
                "[Global randomized] Interface metric must be "
                "nonnegative finite.");
        }
        maximumMetric = std::max(maximumMetric, value);
    }
    if (!(maximumMetric > 0.0)) {
        throw std::runtime_error(
            "[Global randomized] Interface metric is identically zero.");
    }
    const double metricFloor =
        128.0 * std::numeric_limits<double>::epsilon()
        * maximumMetric;
    for (double& value : metric) {
        value = std::max(value, metricFloor);
    }

    model.method = "global-randomized-schur";
    model.fullInterfaceDofs = exactSchur.size();
    model.interfaceGlobalDofs = localModel.interfaceGlobalDofs;
    model.rank = 0;
    model.snapshotUsed = false;
    model.fomUsedForBasis = false;
    model.podUsed = false;
    model.svdUsed = false;
    model.basis.reserve(static_cast<std::size_t>(
        exactSchur.size() * options.requestedRank));
    model.schurImages.reserve(static_cast<std::size_t>(
        exactSchur.size() * options.requestedRank));

    std::mt19937_64 generator(options.seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    bool allSolvesConverged = true;
    const auto deadline = start
        + std::chrono::duration_cast<Clock::duration>(
            std::chrono::duration<double>(
                options.maximumBasisTimeSeconds));
    for (int column = 0; column < options.requestedRank; ++column) {
        if (Clock::now() >= deadline) {
            diagnostics.status = "global_randomized_port_failed";
            diagnostics.failureReason =
                "basis_build_time_gate_failed";
            break;
        }
        std::vector<double> omega(
            static_cast<std::size_t>(exactSchur.size()), 0.0);
        for (double& value : omega) value = normal(generator);
        const auto solveStart = Clock::now();
        SolveResult solve = solveGlobalSchurPcg(
            omega, localPreconditioner, exactSchur,
            options.innerMaximumIterations,
            options.innerTolerance, deadline, diagnostics);
        const double solveSeconds = secondsSince(solveStart);
        ++diagnostics.globalRhsCount;
        diagnostics.globalInnerIterations += solve.iterations;
        diagnostics.maximumInnerIterations = std::max(
            diagnostics.maximumInnerIterations, solve.iterations);
        diagnostics.totalSolveTimeSeconds += solveSeconds;
        diagnostics.maximumSolveTimeSeconds = std::max(
            diagnostics.maximumSolveTimeSeconds, solveSeconds);
        diagnostics.maximumTargetSolveResidual = std::max(
            diagnostics.maximumTargetSolveResidual, solve.residual);
        diagnostics.schurResidual = std::max(
            diagnostics.schurResidual, solve.residual);
        allSolvesConverged = allSolvesConverged && solve.converged;
        if (solve.timedOut) {
            diagnostics.status = "global_randomized_port_failed";
            diagnostics.failureReason =
                "basis_build_time_gate_failed";
            break;
        }
        if (!solve.converged) {
            diagnostics.status =
                "global_randomized_target_solve_gate_failed";
            break;
        }
        if (!appendWeightedColumn(
                model.basis, std::move(solve.solution), metric,
                exactSchur.size(), options.deflationTolerance)) {
            ++diagnostics.deflatedColumns;
            continue;
        }
        ++model.rank;
        const double* accepted = model.basis.data()
            + static_cast<std::size_t>(
                (model.rank - 1) * exactSchur.size());
        std::vector<double> acceptedColumn(
            accepted, accepted + exactSchur.size());
        std::vector<double> image;
        exactSchur.apply(acceptedColumn, image);
        ++diagnostics.globalSchurApplyCount;
        model.schurImages.insert(
            model.schurImages.end(), image.begin(), image.end());
        diagnostics.peakWorkingSetBytes = std::max(
            diagnostics.peakWorkingSetBytes,
            globalCoarseCurrentWorkingSetBytes());
        std::cout
            << "[Global randomized] RHS " << (column + 1)
            << '/' << options.requestedRank
            << ": iterations=" << solve.iterations
            << ", true residual=" << solve.residual
            << ", accepted rank=" << model.rank
            << ", solve=" << solveSeconds << " s"
            << std::endl;
        if (secondsSince(start) >= options.maximumBasisTimeSeconds) {
            diagnostics.status = "global_randomized_port_failed";
            diagnostics.failureReason =
                "basis_build_time_gate_failed";
            break;
        }
    }

    diagnostics.acceptedRank = model.rank;
    diagnostics.meanSolveTimeSeconds =
        diagnostics.globalRhsCount > 0
        ? diagnostics.totalSolveTimeSeconds
            / diagnostics.globalRhsCount
        : 0.0;
    diagnostics.orthogonalityError =
        weightedOrthogonalityError(
            model.basis, metric, exactSchur.size(), model.rank);
    diagnostics.compressionRatio =
        model.rank > 0
        ? static_cast<double>(exactSchur.size()) / model.rank
        : 0.0;
    diagnostics.basisBuildTimeSeconds = secondsSince(start);
    diagnostics.peakWorkingSetBytes = std::max(
        diagnostics.peakWorkingSetBytes,
        globalCoarseCurrentWorkingSetBytes());
    diagnostics.peakIncrementalMemoryBytes =
        diagnostics.peakWorkingSetBytes
            > diagnostics.baselineWorkingSetBytes
        ? diagnostics.peakWorkingSetBytes
            - diagnostics.baselineWorkingSetBytes
        : 0;

    if (!diagnostics.failureReason.empty()) {
        diagnostics.status = "global_randomized_port_failed";
    } else if (!allSolvesConverged) {
        if (diagnostics.status == "not_run") {
            diagnostics.status =
                "global_randomized_target_solve_gate_failed";
            diagnostics.failureReason =
                "target_solve_residual_gate_failed";
        }
    } else if (model.rank != options.requestedRank) {
        diagnostics.status =
            "global_randomized_rank_deflation_gate_failed";
        diagnostics.failureReason =
            "requested_rank_not_accepted";
    } else if (!(diagnostics.maximumTargetSolveResidual
                     <= options.innerTolerance)) {
        diagnostics.status =
            "global_randomized_target_solve_gate_failed";
        diagnostics.failureReason =
            "target_solve_residual_gate_failed";
    } else if (!(diagnostics.schurResidual <= 1.0e-8)) {
        diagnostics.status =
            "global_randomized_schur_residual_gate_failed";
        diagnostics.failureReason =
            "schur_residual_gate_failed";
    } else if (!(diagnostics.orthogonalityError
                     <= options.orthogonalityTolerance)) {
        diagnostics.status =
            "global_randomized_orthogonality_gate_failed";
        diagnostics.failureReason =
            "weighted_orthogonality_gate_failed";
    } else {
        diagnostics.status = "passed";
    }
    return result;
}

void writeGlobalRandomizedSchurDiagnostics(
    const GlobalRandomizedSchurResult& result,
    const std::string& caseName,
    const std::filesystem::path& outputDirectory)
{
    std::filesystem::create_directories(outputDirectory);
    const auto& row = result.diagnostics;
    std::ofstream output(
        outputDirectory
        / "milestone8_global_randomized_diagnostics.csv");
    output
        << "case,composition,global_interface_dofs,"
        "requested_global_rank,global_port_rank,compression_ratio,"
        "seed,deflated_columns,orthogonality_error,schur_residual,"
        "target_solve_residual,basis_build_time_s,"
        "global_schur_apply_count,pardiso_phase33_calls,"
        "global_rhs_count,global_inner_iterations,"
        "maximum_inner_iterations,mean_solve_time_s,"
        "maximum_solve_time_s,total_solve_time_s,"
        "preconditioner_apply_count,baseline_working_set_bytes,"
        "peak_working_set_bytes,peak_incremental_memory_bytes,"
        "inner_solver_requested,inner_solver_actual,"
        "snapshot_used,fom_used_for_basis,pod_used,svd_used,"
        "training_waveform_used,status,failure_reason\n"
        << std::setprecision(17)
        << caseName << ',' << row.composition << ','
        << row.globalInterfaceDofs << ','
        << row.requestedRank << ',' << row.acceptedRank << ','
        << row.compressionRatio << ',' << row.seed << ','
        << row.deflatedColumns << ',' << row.orthogonalityError << ','
        << row.schurResidual << ','
        << row.maximumTargetSolveResidual << ','
        << row.basisBuildTimeSeconds << ','
        << row.globalSchurApplyCount << ','
        << row.pardisoPhase33Calls << ','
        << row.globalRhsCount << ','
        << row.globalInnerIterations << ','
        << row.maximumInnerIterations << ','
        << row.meanSolveTimeSeconds << ','
        << row.maximumSolveTimeSeconds << ','
        << row.totalSolveTimeSeconds << ','
        << row.preconditionerApplyCount << ','
        << row.baselineWorkingSetBytes << ','
        << row.peakWorkingSetBytes << ','
        << row.peakIncrementalMemoryBytes << ','
        << row.innerSolverRequested << ','
        << row.innerSolverActual << ','
        << (row.snapshotUsed ? 1 : 0) << ','
        << (row.fomUsedForBasis ? 1 : 0) << ','
        << (row.podUsed ? 1 : 0) << ','
        << (row.svdUsed ? 1 : 0) << ','
        << (row.trainingWaveformUsed ? 1 : 0) << ','
        << row.status << ',' << row.failureReason << '\n';
}

} // namespace mor::transient
