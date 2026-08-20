#include "reduced_schur_workflow.hpp"

#include "ddm_schur/schur_fgmres.hpp"
#include "model_io.hpp"
#include "mor_diagnostics.hpp"
#include "pod_basis.hpp"
#include "reduced_schur_model.hpp"
#include "source_parameterization.hpp"
#include "deployment_response_model.hpp"
#include "interior_model_serialization.hpp"
#include "local_interior_basis.hpp"
#include "sipg_core.hpp"
#include "config_io.hpp"
#include "mesh_loader.hpp"
#include "fem_assembly.hpp"
#include "linear_solvers.hpp"
#include "diagnostics_io.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>

namespace mor {
namespace {

double vectorNorm(const std::vector<double>& values)
{
    double squared = 0.0;
    for (double value : values) {
        squared += value * value;
    }
    return std::sqrt(squared);
}

ErrorMetrics calculateErrors(const SparseMatrix& system,
                             const std::vector<double>& rhs,
                             const SolveResult& reduced,
                             const std::vector<double>* truth,
                             const std::vector<double>* truthInterface)
{
    ErrorMetrics metrics = reduced.error;
    const std::vector<double> product = system.multiply(reduced.temperature);
    std::vector<double> residual(rhs.size(), 0.0);
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        residual[i] = rhs[i] - product[i];
    }
    metrics.globalRelativeResidual = vectorNorm(residual)
        / std::max(1.0e-300, vectorNorm(rhs));
    if (truthInterface != nullptr
        && truthInterface->size() == reduced.interfaceTemperature.size()) {
        double errorSquared = 0.0;
        double truthSquared = 0.0;
        for (std::size_t i = 0; i < truthInterface->size(); ++i) {
            const double error = reduced.interfaceTemperature[i] - (*truthInterface)[i];
            errorSquared += error * error;
            truthSquared += (*truthInterface)[i] * (*truthInterface)[i];
        }
        metrics.interfaceRelativeL2 = std::sqrt(errorSquared)
            / std::max(1.0e-300, std::sqrt(truthSquared));
    }
    if (truth != nullptr && truth->size() == reduced.temperature.size()) {
        double errorSquared = 0.0;
        double truthSquared = 0.0;
        double absoluteSum = 0.0;
        double maximum = 0.0;
        for (std::size_t i = 0; i < truth->size(); ++i) {
            const double error = reduced.temperature[i] - (*truth)[i];
            errorSquared += error * error;
            truthSquared += (*truth)[i] * (*truth)[i];
            absoluteSum += std::abs(error);
            maximum = std::max(maximum, std::abs(error));
        }
        metrics.relativeL2 = std::sqrt(errorSquared)
            / std::max(1.0e-300, std::sqrt(truthSquared));
        metrics.maximumAbsolute = maximum;
        metrics.meanAbsolute = absoluteSum / static_cast<double>(truth->size());
        const double reducedMaximum = *std::max_element(
            reduced.temperature.begin(), reduced.temperature.end());
        const double truthMaximum = *std::max_element(truth->begin(), truth->end());
        metrics.maximumTemperatureError = std::abs(reducedMaximum - truthMaximum);
    }
    return metrics;
}

std::vector<double> nominalPowers(const SourceParameterization& sources)
{
    std::vector<double> powers;
    powers.reserve(sources.channels.size());
    for (const SourceChannel& channel : sources.channels) {
        powers.push_back(channel.nominalPowerW);
    }
    return powers;
}

int energySelectedRank(const PodResult& pod, double tolerance)
{
    for (int i = 0; i < pod.numericalRank; ++i) {
        if (1.0 - pod.retainedEnergy[static_cast<std::size_t>(i)] <= tolerance) {
            return i + 1;
        }
    }
    return pod.numericalRank;
}

struct CaseRow {
    int rank = 0;
    int caseIndex = 0;
    std::string split;
    std::string family;
    int correctionIterations = 0;
    ErrorMetrics error;
    OnlineTiming timing;
    double fomSeconds = 0.0;
    std::string status;
};

struct FullOrderSample {
    std::vector<double> temperature;
    std::vector<double> interfaceTemperature;
    std::string status;
};

void writeCaseRows(const std::vector<CaseRow>& rows,
                   const std::filesystem::path& path)
{
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("[MOR] Cannot write case result table.");
    }
    out << "rank,case,split,family,correction_iterations,interface_relative_residual,"
        << "global_relative_residual,interface_relative_l2,full_field_relative_l2,maximum_absolute_error,mean_absolute_error,"
        << "maximum_temperature_error,interface_jump_rms,relative_heat_flux_imbalance,"
        << "rhs_seconds,condensation_seconds,projection_seconds,"
        << "reduced_solve_seconds,interface_reconstruction_seconds,recovery_seconds,correction_seconds,total_online_seconds,"
        << "fom_seconds,status\n";
    out << std::setprecision(17);
    for (const CaseRow& row : rows) {
        out << row.rank << ',' << row.caseIndex << ',' << row.split << ',' << row.family << ','
            << row.correctionIterations << ',' << row.error.interfaceRelativeResidual << ','
            << row.error.globalRelativeResidual << ',' << row.error.interfaceRelativeL2 << ','
            << row.error.relativeL2 << ','
            << row.error.maximumAbsolute << ',' << row.error.meanAbsolute << ','
            << row.error.maximumTemperatureError << ','
            << row.error.interfaceJump << ',' << row.error.heatFluxImbalance << ','
            << row.timing.rhsSeconds << ','
            << row.timing.condensationSeconds << ',' << row.timing.projectionSeconds << ','
            << row.timing.reducedSolveSeconds << ','
            << row.timing.interfaceReconstructionSeconds << ','
            << row.timing.recoverySeconds << ','
            << row.timing.correctionSeconds << ',' << row.timing.totalSeconds << ','
            << row.fomSeconds << ',' << row.status << '\n';
    }
}

void writeSingularValues(const PodResult& pod,
                         const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "index,singular_value,relative_singular_value,retained_energy\n";
    out << std::setprecision(17);
    const double first = pod.singularValues.empty() ? 1.0 : pod.singularValues.front();
    for (std::size_t i = 0; i < pod.singularValues.size(); ++i) {
        out << i << ',' << pod.singularValues[i] << ','
            << pod.singularValues[i] / std::max(1.0e-300, first) << ','
            << pod.retainedEnergy[i] << '\n';
    }
}

} // namespace

WorkflowResult runReducedSchurWorkflow(
    const Mesh& mesh,
    const SparseMatrix& system,
    const CaseConfig& physics,
    const std::vector<double>& assembledSource,
    const std::vector<double>& heatOnlySource,
    const std::vector<double>& fixedAdjust,
    const ddm_schur::Options& schurOptions,
    const Options& options,
    const std::filesystem::path& outputDirectory)
{
    if (options.mode != "pure" && options.mode != "corrected") {
        throw std::runtime_error("--mor-mode must be pure or corrected.");
    }
    if (!options.exactInteriorRecovery) {
        throw std::runtime_error(
            "Stage 2A.1 offline construction and corrected mode require the existing exact local recovery.");
    }
    if (options.interiorMode != "pardiso"
        && options.interiorMode != "exact-response"
        && options.interiorMode != "compressed-rb") {
        throw std::runtime_error(
            "--mor-interior-mode must be pardiso, exact-response, or compressed-rb.");
    }
    if (options.storagePrecision != "float64" && options.storagePrecision != "float32") {
        throw std::runtime_error("--mor-storage-precision must be float64 or float32.");
    }
    if (options.snapshotSolver != "auto"
        && options.snapshotSolver != "direct"
        && options.snapshotSolver != "schur") {
        throw std::runtime_error(
            "--mor-snapshot-solver must be auto, direct, or schur.");
    }
    if (!options.generateModel && options.loadModelPath.empty()) {
        throw std::runtime_error(
            "reduced-schur requires --mor-generate-model or --mor-load-model <directory>.");
    }
    std::filesystem::create_directories(outputDirectory);
    const auto workflowStart = std::chrono::steady_clock::now();

    const auto sourceStart = std::chrono::steady_clock::now();
    const SourceParameterization sources = buildSourceParameterization(
        mesh, physics, assembledSource, heatOnlySource, fixedAdjust);
    const double sourceSetupSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - sourceStart).count();
    if (sources.channels.empty()) {
        throw std::runtime_error("[MOR] No independent heat-source channels were found.");
    }
    writeSourceChannels(sources, outputDirectory / "mor_source_channels.csv");

    std::string snapshotSolver = options.snapshotSolver;
    if (snapshotSolver == "auto") {
#ifdef USE_MKL_PARDISO
        snapshotSolver = "direct";
#else
        snapshotSolver = "schur";
#endif
    }
#ifndef USE_MKL_PARDISO
    if (snapshotSolver == "direct") {
        throw std::runtime_error(
            "[MOR] Direct snapshot generation requires an MKL PARDISO build.");
    }
#endif
    const bool needsFomBackend = options.generateModel || options.compareFom;
    const bool needsIterativeStage1 = options.mode == "corrected"
        || (snapshotSolver == "schur" && needsFomBackend);
    ddm_schur::Options effectiveSchurOptions = schurOptions;
    if (!needsIterativeStage1) {
        // Pure Reduced Schur needs local factors for exact S applies,
        // condensation, and recovery, but it never applies the Stage 1
        // iterative preconditioner.  Avoid constructing an unused proxy.
        effectiveSchurOptions.proxyEnabled = false;
        effectiveSchurOptions.proxyDiagnostics = false;
    }

    const auto schurSetupStart = std::chrono::steady_clock::now();
    ddm_schur::DdmSchurSolver solver(mesh, system, physics, effectiveSchurOptions);
    const double schurSetupWallSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - schurSetupStart).count();
    const Fingerprints fingerprints = computeFingerprints(
        mesh, system, physics, solver.interfaceGlobalDofs());

    const std::vector<ParameterCase> trainingCases = makeParameterCases(
        sources, "training", options.trainingCases, options.seed, true);
    const std::vector<ParameterCase> validationCases = makeParameterCases(
        sources, "validation", options.validationCases, options.seed + 104729ULL, false);
    const std::vector<ParameterCase> testCases = makeParameterCases(
        sources, "test", options.testCases, options.seed + 224737ULL, false);
    writeParameterCases(trainingCases, outputDirectory / "mor_training_cases.csv");
    writeParameterCases(validationCases, outputDirectory / "mor_validation_cases.csv");
    writeParameterCases(testCases, outputDirectory / "mor_test_cases.csv");

    std::unique_ptr<SubdomainDirectSolver> directFom;
    double snapshotFomSetupSeconds = 0.0;
    if (snapshotSolver == "direct" && needsFomBackend) {
        const auto directSetupStart = std::chrono::steady_clock::now();
        directFom = std::make_unique<SubdomainDirectSolver>(
            system.size(), sparseMatrixEntries(system));
        snapshotFomSetupSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - directSetupStart).count();
    }
    const auto solveFom = [&](const std::vector<double>& rhs) {
        FullOrderSample sample;
        if (directFom) {
            directFom->solve(rhs, sample.temperature);
            sample.interfaceTemperature.resize(
                static_cast<std::size_t>(solver.interfaceDofs()));
            const std::vector<int>& interfaceDofs = solver.interfaceGlobalDofs();
            for (int row = 0; row < solver.interfaceDofs(); ++row) {
                sample.interfaceTemperature[static_cast<std::size_t>(row)] =
                    sample.temperature[static_cast<std::size_t>(
                        interfaceDofs[static_cast<std::size_t>(row)])];
            }
            sample.status = "success";
        } else {
            ddm_schur::SolveResult result = solver.solve(rhs);
            sample.temperature = std::move(result.temperature);
            sample.interfaceTemperature = std::move(result.interfaceSolution);
            sample.status = result.report.status;
        }
        return sample;
    };

    ReducedSchurModel selectedModel;
    ReducedSchurModel sweepBaseModel;
    SnapshotDatabase trainingSnapshots;
    PodResult pod;
    double snapshotSeconds = 0.0;
    double reducedOperatorSeconds = 0.0;
    double serializationSeconds = 0.0;
    int snapshotCount = 0;
    if (options.generateModel) {
        const auto snapshotStart = std::chrono::steady_clock::now();
        FullOrderSample reference = solveFom(sources.referenceRhs);
        if (reference.status != "success") {
            throw std::runtime_error("[MOR] Reference FOM solve did not converge.");
        }
        std::vector<std::vector<double>> unitResponses;
        unitResponses.reserve(sources.channels.size());
        for (std::size_t channel = 0; channel < sources.channels.size(); ++channel) {
            std::vector<double> powers(sources.channels.size(), 0.0);
            powers[channel] = 1.0;
            FullOrderSample unit = solveFom(composeRhs(sources, powers));
            if (unit.status != "success") {
                throw std::runtime_error("[MOR] Unit-channel FOM solve did not converge.");
            }
            std::vector<double> response(unit.interfaceTemperature.size(), 0.0);
            for (std::size_t row = 0; row < response.size(); ++row) {
                response[row] = unit.interfaceTemperature[row]
                    - reference.interfaceTemperature[row];
            }
            unitResponses.push_back(std::move(response));
        }
        trainingSnapshots.rows = solver.interfaceDofs();
        trainingSnapshots.cases = trainingCases;
        trainingSnapshots.values.assign(static_cast<std::size_t>(trainingSnapshots.rows)
            * trainingCases.size(), 0.0);
        for (std::size_t column = 0; column < trainingCases.size(); ++column) {
            for (std::size_t channel = 0; channel < sources.channels.size(); ++channel) {
                const double power = trainingCases[column].powersW[channel];
                if (power == 0.0) {
                    continue;
                }
                const auto& response = unitResponses[channel];
                for (int row = 0; row < trainingSnapshots.rows; ++row) {
                    trainingSnapshots.values[column * static_cast<std::size_t>(trainingSnapshots.rows)
                        + static_cast<std::size_t>(row)] += power
                        * response[static_cast<std::size_t>(row)];
                }
            }
        }
        snapshotCount = static_cast<int>(trainingCases.size());
        snapshotSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - snapshotStart).count();

        pod = buildGramPod(trainingSnapshots, static_cast<int>(trainingCases.size()),
                           options.energyTolerance, options.singularValueTolerance);
        // With fixed A and b=b0+Bp, no response space can have a physical
        // dimension larger than rank(B).  Gram eigenvalues square the
        // condition number, so roundoff can otherwise manufacture tiny modes.
        pod.numericalRank = std::min(pod.numericalRank,
            static_cast<int>(sources.channels.size()));
        pod.selectedRank = std::min(pod.selectedRank, pod.numericalRank);
        writeSingularValues(pod, outputDirectory / "mor_singular_values.csv");
        const int selectedRank = options.rank > 0
            ? std::min(options.rank, pod.numericalRank)
            : energySelectedRank(pod, options.energyTolerance);
        const auto reducedStart = std::chrono::steady_clock::now();
        sweepBaseModel = buildReducedSchurModel(
            solver, pod, pod.selectedRank, reference.interfaceTemperature, fingerprints);
        sweepBaseModel.sourceChannels = static_cast<int>(sources.channels.size());
        selectedModel = truncateModel(sweepBaseModel, selectedRank);
        selectedModel.sourceChannels = static_cast<int>(sources.channels.size());
        reducedOperatorSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - reducedStart).count();
        if (!options.saveModelPath.empty()) {
            const auto serializationStart = std::chrono::steady_clock::now();
            saveModel(selectedModel, options.saveModelPath);
            writeSourceChannels(sources, options.saveModelPath / "source_channels.csv");
            writeParameterCases(trainingCases,
                                options.saveModelPath / "training_cases.csv");
            writeParameterCases(validationCases,
                                options.saveModelPath / "validation_cases.csv");
            writeParameterCases(testCases,
                                options.saveModelPath / "test_cases.csv");
            writeSingularValues(pod,
                                options.saveModelPath / "singular_values.csv");
            serializationSeconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - serializationStart).count();
        }
    } else {
        selectedModel = loadModel(options.loadModelPath, fingerprints,
                                  solver.interfaceDofs(),
                                  static_cast<int>(sources.channels.size()));
        snapshotCount = 0;
        sweepBaseModel = selectedModel;
    }

    std::cout << std::setprecision(12)
              << "[MOR] snapshot/truth FOM backend: " << snapshotSolver << '\n'
              << "[MOR] Stage 1 iterative proxy constructed: "
              << (needsIterativeStage1 ? "yes" : "no") << '\n'
              << "[MOR] independent source channels: " << sources.channels.size() << '\n'
              << "[MOR] snapshots/numerical rank/selected rank: "
              << snapshotCount << '/' << (pod.numericalRank > 0 ? pod.numericalRank : selectedModel.rank)
              << '/' << selectedModel.rank << '\n'
              << "[MOR] reduced Schur symmetry error: " << selectedModel.symmetryError << '\n'
              << "[MOR] reduced symmetric eigenvalue range: ["
              << selectedModel.minimumSymmetricEigenvalue << ", "
              << selectedModel.maximumSymmetricEigenvalue << "]\n";

    const auto onlineSolverSetupStart = std::chrono::steady_clock::now();
    ReducedSchurOnlineSolver selectedOnline(solver, selectedModel);
    const double onlineSolverSetupSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - onlineSolverSetupStart).count();

    std::vector<int> ranks;
    if (options.generateModel) {
        for (int requested : options.rankSweep) {
            ranks.push_back(std::min(requested, pod.numericalRank));
        }
    }
    ranks.push_back(selectedModel.rank);
    ranks.erase(std::remove_if(ranks.begin(), ranks.end(), [](int rank) { return rank <= 0; }),
                ranks.end());
    std::sort(ranks.begin(), ranks.end());
    ranks.erase(std::unique(ranks.begin(), ranks.end()), ranks.end());

    std::vector<ParameterCase> evaluationCases;
    if (options.generateModel) {
        evaluationCases = trainingCases;
    }
    const std::size_t trainingEvaluationCount = evaluationCases.size();
    evaluationCases.insert(evaluationCases.end(), validationCases.begin(), validationCases.end());
    evaluationCases.insert(evaluationCases.end(), testCases.begin(), testCases.end());
    std::vector<std::vector<double>> truthTemperatures(evaluationCases.size());
    std::vector<std::vector<double>> truthInterfaces(evaluationCases.size());
    std::vector<char> truthAvailable(evaluationCases.size(), 0);
    std::vector<double> truthSeconds(evaluationCases.size(), 0.0);
    for (std::size_t index = 0; index < trainingEvaluationCount; ++index) {
        truthInterfaces[index] = selectedModel.referenceInterface;
        const std::size_t offset = index
            * static_cast<std::size_t>(trainingSnapshots.rows);
        for (int row = 0; row < trainingSnapshots.rows; ++row) {
            truthInterfaces[index][static_cast<std::size_t>(row)] +=
                trainingSnapshots.values[offset + static_cast<std::size_t>(row)];
        }
        const std::vector<double> rhs = composeRhs(sources, evaluationCases[index].powersW);
        truthTemperatures[index] = solver.recover(rhs, truthInterfaces[index]);
        truthAvailable[index] = 1;
    }
    if (options.compareFom) {
        for (std::size_t index = trainingEvaluationCount;
             index < evaluationCases.size(); ++index) {
            const auto truthStart = std::chrono::steady_clock::now();
            const std::vector<double> rhs = composeRhs(
                sources, evaluationCases[index].powersW);
            FullOrderSample truth = solveFom(rhs);
            truthSeconds[index] = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - truthStart).count();
            if (truth.status != "success") {
                throw std::runtime_error("[MOR] Validation/test FOM solve did not converge.");
            }
            truthTemperatures[index] = std::move(truth.temperature);
            truthInterfaces[index] = std::move(truth.interfaceTemperature);
            truthAvailable[index] = 1;
        }
    }
    // Snapshot/truth construction is complete.  Release a global direct
    // factor before rank-sweep recovery so the reported online path retains
    // only the Stage 1 local/proxy factors used in production.
    directFom.reset();

    const bool deploymentRequested = options.interiorMode != "pardiso"
        || options.precomputePowerResponse || options.compareInteriorModes;
    DeploymentResponseModel deploymentModel;
    bool deploymentAvailable = false;
    double deploymentSerializationSeconds = 0.0;
    if (deploymentRequested) {
        if (options.generateModel) {
            deploymentModel = buildDeploymentResponseModel(
                mesh, sources, solver.interfaceGlobalDofs(), selectedModel, selectedOnline);
            if (options.interiorMode == "compressed-rb" || options.compareInteriorModes) {
                buildCompressedInteriorBases(
                    deploymentModel, options.interiorRank,
                    options.interiorEnergyTolerance,
                    options.interiorSingularValueTolerance);
            }
            if (!options.saveModelPath.empty()) {
                const std::string storedMode = options.interiorMode == "compressed-rb"
                    ? "compressed-rb" : "exact-response";
                saveDeploymentResponseModel(deploymentModel, options.saveModelPath,
                    storedMode, options.storagePrecision);
                deploymentSerializationSeconds = deploymentModel.serializationSeconds;
            }
            deploymentAvailable = true;
        } else {
            deploymentModel = loadDeploymentResponseModel(
                options.loadModelPath, options.interiorMode,
                &fingerprints, system.size(), solver.interfaceDofs(),
                selectedModel.rank, static_cast<int>(sources.channels.size()));
            deploymentAvailable = true;
        }
    }

    std::vector<CaseRow> caseRows;
    std::ofstream rankOut(outputDirectory / "mor_rank_sweep.csv");
    rankOut << "requested_or_clamped_rank,cases,average_online_seconds,average_fom_seconds,"
        << "maximum_interface_residual,maximum_interface_relative_l2,maximum_global_residual,maximum_relative_l2,"
        << "maximum_absolute_error,break_even_queries\n";
    rankOut << std::setprecision(17);
    double selectedAverageOnline = 0.0;
    double selectedAverageFom = 0.0;
    for (int rank : ranks) {
        ReducedSchurModel rankModel = rank == selectedModel.rank
            ? selectedModel : truncateModel(sweepBaseModel, rank);
        rankModel.sourceChannels = static_cast<int>(sources.channels.size());
        ReducedSchurOnlineSolver online(solver, rankModel);
        double onlineSum = 0.0;
        double fomSum = 0.0;
        int fomCount = 0;
        double maxInterface = 0.0;
        double maxInterfaceL2 = 0.0;
        double maxGlobal = 0.0;
        double maxL2 = 0.0;
        double maxAbsolute = 0.0;
        for (std::size_t index = 0; index < evaluationCases.size(); ++index) {
            const auto rhsStart = std::chrono::steady_clock::now();
            const std::vector<double> rhs = composeRhs(sources, evaluationCases[index].powersW);
            const double rhsSeconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - rhsStart).count();
            const bool useDeployment = deploymentAvailable
                && rank == selectedModel.rank
                && options.interiorMode != "pardiso"
                && options.mode == "pure";
            SolveResult reduced = useDeployment
                ? evaluateDeploymentResponse(deploymentModel,
                    evaluationCases[index].powersW, options.interiorMode)
                : online.solve(rhs, false);
            if (!useDeployment) {
                reduced.timing.rhsSeconds = rhsSeconds;
                reduced.timing.totalSeconds += rhsSeconds;
            }
            if (useDeployment) {
                const std::vector<double> condensed = solver.condensedRhs(rhs);
                std::vector<double> image;
                solver.applyExactSchur(reduced.interfaceTemperature, image);
                std::vector<double> interfaceResidual(condensed.size(), 0.0);
                for (std::size_t row = 0; row < condensed.size(); ++row) {
                    interfaceResidual[row] = condensed[row] - image[row];
                }
                reduced.error.interfaceRelativeResidual = vectorNorm(interfaceResidual)
                    / std::max(1.0e-300, vectorNorm(condensed));
            }
            reduced.error = calculateErrors(system, rhs, reduced,
                truthAvailable[index] ? &truthTemperatures[index] : nullptr,
                truthAvailable[index] ? &truthInterfaces[index] : nullptr);
            const InterfacePhysicsMetrics interfacePhysics =
                calculateInterfacePhysicsMetrics(mesh, physics, reduced.temperature);
            reduced.error.interfaceJump = interfacePhysics.temperatureJumpRms;
            reduced.error.heatFluxImbalance = interfacePhysics.relativeFluxImbalance;
            CaseRow row;
            row.rank = rank;
            row.caseIndex = evaluationCases[index].index;
            row.split = evaluationCases[index].split;
            row.family = evaluationCases[index].family;
            row.error = reduced.error;
            row.timing = reduced.timing;
            row.fomSeconds = truthSeconds[index];
            row.status = reduced.status;
            caseRows.push_back(row);
            onlineSum += reduced.timing.totalSeconds;
            fomSum += truthSeconds[index];
            if (truthSeconds[index] > 0.0) {
                ++fomCount;
            }
            maxInterface = std::max(maxInterface, reduced.error.interfaceRelativeResidual);
            maxInterfaceL2 = std::max(maxInterfaceL2, reduced.error.interfaceRelativeL2);
            maxGlobal = std::max(maxGlobal, reduced.error.globalRelativeResidual);
            maxL2 = std::max(maxL2, reduced.error.relativeL2);
            maxAbsolute = std::max(maxAbsolute, reduced.error.maximumAbsolute);
        }
        const double caseCount = static_cast<double>(std::max<std::size_t>(1, evaluationCases.size()));
        const double averageOnline = onlineSum / caseCount;
        const double averageFom = fomCount > 0
            ? fomSum / static_cast<double>(fomCount) : 0.0;
        if (rank == selectedModel.rank) {
            selectedAverageOnline = averageOnline;
            selectedAverageFom = averageFom;
        }
        const double offline = sourceSetupSeconds + schurSetupWallSeconds + snapshotSeconds
            + pod.gramSeconds + pod.eigenSeconds + pod.basisSeconds
            + reducedOperatorSeconds + serializationSeconds + onlineSolverSetupSeconds
            + snapshotFomSetupSeconds
            + deploymentModel.responseConstructionSeconds
            + deploymentModel.compressionSeconds + deploymentSerializationSeconds;
        const double breakEven = averageFom > averageOnline
            ? offline / (averageFom - averageOnline) : -1.0;
        rankOut << rank << ',' << evaluationCases.size() << ',' << averageOnline << ','
                << averageFom << ',' << maxInterface << ',' << maxInterfaceL2 << ','
                << maxGlobal << ',' << maxL2 << ','
                << maxAbsolute << ',' << breakEven << '\n';
    }
    writeCaseRows(caseRows, outputDirectory / "mor_case_results.csv");
    if (deploymentAvailable) {
        writeCaseRows(caseRows, outputDirectory / "local_interior_accuracy_by_case.csv");
        std::ofstream subdomainRanks(outputDirectory / "local_interior_subdomain_ranks.csv");
        subdomainRanks << "subdomain,interior_dofs,direct_power_channels,rank,numerical_rank,retained_energy\n";
        subdomainRanks << std::setprecision(17);
        for (const InteriorResponseBlock& block : deploymentModel.interiors) {
            const double energy = block.rank > 0
                && static_cast<std::size_t>(block.rank) <= block.retainedEnergy.size()
                ? block.retainedEnergy[static_cast<std::size_t>(block.rank - 1)] : 1.0;
            subdomainRanks << block.subdomain << ',' << block.globalDofs.size() << ','
                << block.directPowerChannels.size() << ',' << block.rank << ','
                << block.singularValues.size() << ',' << energy << '\n';
        }
        std::ofstream memory(outputDirectory / "local_interior_memory_breakdown.csv");
        const std::uint64_t interfaceBytes = static_cast<std::uint64_t>(
            deploymentModel.interfacePowerResponse.size() * sizeof(double));
        std::uint64_t exactInteriorBytes = 0;
        std::uint64_t compressedInteriorBytes = 0;
        for (const InteriorResponseBlock& block : deploymentModel.interiors) {
            exactInteriorBytes += static_cast<std::uint64_t>(
                block.exactResponse.size() * sizeof(double));
            compressedInteriorBytes += static_cast<std::uint64_t>(
                (block.localBasis.size() + block.localCoordinateMap.size()) * sizeof(double));
        }
        memory << "component,bytes\n"
            << "interface_power_response," << interfaceBytes << '\n'
            << "exact_interior_responses," << exactInteriorBytes << '\n'
            << "compressed_interior_data," << compressedInteriorBytes << '\n'
            << "deployment_files," << deploymentModel.deploymentFileBytes << '\n'
            << "complete_model_files," << deploymentModel.completeModelFileBytes << '\n'
            << "peak_working_set," << peakWorkingSetBytes() << '\n';
        std::ofstream sizes(outputDirectory / "local_interior_model_sizes.csv");
        sizes << "interior_mode,storage_precision,deployment_file_bytes,complete_model_file_bytes\n"
              << deploymentModel.interiorMode << ',' << deploymentModel.storagePrecision << ','
              << deploymentModel.deploymentFileBytes << ','
              << deploymentModel.completeModelFileBytes << '\n';
        std::ofstream onlineTiming(outputDirectory / "local_interior_online_timing.csv");
        onlineTiming << "rank,case,split,interface_seconds,interior_seconds,total_compute_seconds\n"
                     << std::setprecision(17);
        for (const CaseRow& row : caseRows) {
            if (row.rank == selectedModel.rank) {
                onlineTiming << row.rank << ',' << row.caseIndex << ',' << row.split << ','
                    << row.timing.interfaceReconstructionSeconds << ','
                    << row.timing.recoverySeconds << ',' << row.timing.totalSeconds << '\n';
            }
        }
        double stage2aPardisoAverage = -1.0;
        if (options.compareInteriorModes && !evaluationCases.empty()) {
            double sum = 0.0;
            for (const ParameterCase& parameterCase : evaluationCases) {
                const std::vector<double> rhs = composeRhs(sources, parameterCase.powersW);
                sum += selectedOnline.solve(rhs, false).timing.totalSeconds;
            }
            stage2aPardisoAverage = sum / static_cast<double>(evaluationCases.size());
        }
        const double additionalOffline = deploymentModel.responseConstructionSeconds
            + deploymentModel.compressionSeconds + deploymentSerializationSeconds;
        std::ofstream breakEven(outputDirectory / "local_interior_break_even.csv");
        breakEven << "reference,reference_online_seconds,new_online_seconds,"
                     "additional_offline_seconds,break_even_rhs\n"
                  << std::setprecision(17);
        const double directBreakEven = selectedAverageFom > selectedAverageOnline
            ? additionalOffline / (selectedAverageFom - selectedAverageOnline) : -1.0;
        breakEven << "monolithic_pardiso_reuse," << selectedAverageFom << ','
                  << selectedAverageOnline << ',' << additionalOffline << ','
                  << directBreakEven << '\n';
        const double stage2aBreakEven = stage2aPardisoAverage > selectedAverageOnline
            ? additionalOffline / (stage2aPardisoAverage - selectedAverageOnline) : -1.0;
        breakEven << "stage2a_pardiso_interior," << stage2aPardisoAverage << ','
                  << selectedAverageOnline << ',' << additionalOffline << ','
                  << stage2aBreakEven << '\n';
    }
    struct SplitAggregate {
        int count = 0;
        double online = 0.0;
        double maximumInterfaceL2 = 0.0;
        double maximumFullL2 = 0.0;
        double maximumAbsolute = 0.0;
        double maximumGlobalResidual = 0.0;
    };
    std::map<std::pair<int, std::string>, SplitAggregate> splitAggregates;
    for (const CaseRow& row : caseRows) {
        SplitAggregate& aggregate = splitAggregates[{row.rank, row.split}];
        ++aggregate.count;
        aggregate.online += row.timing.totalSeconds;
        aggregate.maximumInterfaceL2 = std::max(
            aggregate.maximumInterfaceL2, row.error.interfaceRelativeL2);
        aggregate.maximumFullL2 = std::max(
            aggregate.maximumFullL2, row.error.relativeL2);
        aggregate.maximumAbsolute = std::max(
            aggregate.maximumAbsolute, row.error.maximumAbsolute);
        aggregate.maximumGlobalResidual = std::max(
            aggregate.maximumGlobalResidual, row.error.globalRelativeResidual);
    }
    std::ofstream splitOut(outputDirectory / "mor_rank_split_summary.csv");
    splitOut << "rank,split,cases,average_online_seconds,maximum_interface_relative_l2,"
        << "maximum_full_field_relative_l2,maximum_absolute_error,maximum_global_residual\n";
    splitOut << std::setprecision(17);
    for (const auto& entry : splitAggregates) {
        const SplitAggregate& aggregate = entry.second;
        splitOut << entry.first.first << ',' << entry.first.second << ',' << aggregate.count << ','
                 << aggregate.online / static_cast<double>(std::max(1, aggregate.count)) << ','
                 << aggregate.maximumInterfaceL2 << ',' << aggregate.maximumFullL2 << ','
                 << aggregate.maximumAbsolute << ',' << aggregate.maximumGlobalResidual << '\n';
    }

    const std::vector<double> nominal = nominalPowers(sources);
    const auto nominalRhsStart = std::chrono::steady_clock::now();
    const std::vector<double> nominalRhs = composeRhs(sources, nominal);
    const double nominalRhsSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - nominalRhsStart).count();
    const bool nominalUsesDeployment = deploymentAvailable
        && options.interiorMode != "pardiso" && options.mode == "pure";
    SolveResult nominalResult = nominalUsesDeployment
        ? evaluateDeploymentResponse(deploymentModel, nominal, options.interiorMode)
        : selectedOnline.solve(nominalRhs, options.mode == "corrected");
    if (!nominalUsesDeployment) {
        nominalResult.timing.rhsSeconds = nominalRhsSeconds;
        nominalResult.timing.totalSeconds += nominalRhsSeconds;
    } else {
        const std::vector<double> condensed = solver.condensedRhs(nominalRhs);
        std::vector<double> image;
        solver.applyExactSchur(nominalResult.interfaceTemperature, image);
        std::vector<double> interfaceResidual(condensed.size(), 0.0);
        for (std::size_t row = 0; row < condensed.size(); ++row) {
            interfaceResidual[row] = condensed[row] - image[row];
        }
        nominalResult.error.interfaceRelativeResidual = vectorNorm(interfaceResidual)
            / std::max(1.0e-300, vectorNorm(condensed));
    }
    nominalResult.error = calculateErrors(system, nominalRhs, nominalResult, nullptr, nullptr);
    const InterfacePhysicsMetrics nominalInterfacePhysics =
        calculateInterfacePhysicsMetrics(mesh, physics, nominalResult.temperature);
    nominalResult.error.interfaceJump = nominalInterfacePhysics.temperatureJumpRms;
    nominalResult.error.heatFluxImbalance = nominalInterfacePhysics.relativeFluxImbalance;

    const double offlineSeconds = sourceSetupSeconds + schurSetupWallSeconds + snapshotSeconds
        + pod.gramSeconds + pod.eigenSeconds + pod.basisSeconds
        + reducedOperatorSeconds + serializationSeconds + onlineSolverSetupSeconds
        + snapshotFomSetupSeconds + deploymentModel.responseConstructionSeconds
        + deploymentModel.compressionSeconds + deploymentSerializationSeconds;
    std::ofstream summary(outputDirectory / "mor_summary.json");
    summary << std::setprecision(17)
        << "{\n"
        << "  \"source_channels\": " << sources.channels.size() << ",\n"
        << "  \"snapshot_count\": " << snapshotCount << ",\n"
        << "  \"numerical_rank\": "
        << (pod.numericalRank > 0 ? pod.numericalRank : selectedModel.rank) << ",\n"
        << "  \"selected_rank\": " << selectedModel.rank << ",\n"
        << "  \"mode\": \"" << options.mode << "\",\n"
        << "  \"interior_mode\": \"" << options.interiorMode << "\",\n"
        << "  \"snapshot_fom_backend\": \"" << snapshotSolver << "\",\n"
        << "  \"stage1_iterative_proxy_constructed\": "
        << (needsIterativeStage1 ? "true" : "false") << ",\n"
        << "  \"dense_factorization\": \"" << selectedOnline.factorizationType() << "\",\n"
        << "  \"reduced_operator_symmetry_error\": "
        << selectedModel.symmetryError << ",\n"
        << "  \"reduced_operator_minimum_symmetric_eigenvalue\": "
        << selectedModel.minimumSymmetricEigenvalue << ",\n"
        << "  \"reduced_operator_maximum_symmetric_eigenvalue\": "
        << selectedModel.maximumSymmetricEigenvalue << ",\n"
        << "  \"reduced_operator_condition_estimate\": "
        << selectedModel.conditionEstimate << ",\n"
        << "  \"source_setup_seconds\": " << sourceSetupSeconds << ",\n"
        << "  \"schur_setup_seconds\": " << schurSetupWallSeconds << ",\n"
        << "  \"snapshot_fom_factorization_seconds\": "
        << snapshotFomSetupSeconds << ",\n"
        << "  \"snapshot_seconds\": " << snapshotSeconds << ",\n"
        << "  \"snapshot_storage_bytes\": "
        << trainingSnapshots.values.size() * sizeof(double) << ",\n"
        << "  \"gram_seconds\": " << pod.gramSeconds << ",\n"
        << "  \"pod_eigen_seconds\": " << pod.eigenSeconds << ",\n"
        << "  \"pod_basis_seconds\": " << pod.basisSeconds << ",\n"
        << "  \"reduced_operator_seconds\": " << reducedOperatorSeconds << ",\n"
        << "  \"exact_schur_basis_apply_seconds\": "
        << selectedModel.exactSchurApplySeconds << ",\n"
        << "  \"reduced_matrix_assembly_seconds\": "
        << selectedModel.reducedAssemblySeconds << ",\n"
        << "  \"model_serialization_seconds\": " << serializationSeconds << ",\n"
        << "  \"interior_response_construction_seconds\": "
        << deploymentModel.responseConstructionSeconds << ",\n"
        << "  \"interior_compression_seconds\": "
        << deploymentModel.compressionSeconds << ",\n"
        << "  \"deployment_serialization_seconds\": "
        << deploymentSerializationSeconds << ",\n"
        << "  \"deployment_model_file_bytes\": "
        << deploymentModel.deploymentFileBytes << ",\n"
        << "  \"online_solver_setup_seconds\": " << onlineSolverSetupSeconds << ",\n"
        << "  \"peak_working_set_bytes\": " << peakWorkingSetBytes() << ",\n"
        << "  \"offline_total_seconds\": " << offlineSeconds << ",\n"
        << "  \"nominal_online_seconds\": " << nominalResult.timing.totalSeconds << ",\n"
        << "  \"nominal_rhs_seconds\": " << nominalResult.timing.rhsSeconds << ",\n"
        << "  \"nominal_condensation_seconds\": "
        << nominalResult.timing.condensationSeconds << ",\n"
        << "  \"nominal_projection_seconds\": "
        << nominalResult.timing.projectionSeconds << ",\n"
        << "  \"nominal_reduced_solve_seconds\": "
        << nominalResult.timing.reducedSolveSeconds << ",\n"
        << "  \"nominal_interface_reconstruction_seconds\": "
        << nominalResult.timing.interfaceReconstructionSeconds << ",\n"
        << "  \"nominal_exact_interior_recovery_seconds\": "
        << nominalResult.timing.recoverySeconds << ",\n"
        << "  \"nominal_exact_correction_seconds\": "
        << nominalResult.timing.correctionSeconds << ",\n"
        << "  \"nominal_interface_residual\": "
        << nominalResult.error.interfaceRelativeResidual << ",\n"
        << "  \"nominal_global_residual\": " << nominalResult.error.globalRelativeResidual << ",\n"
        << "  \"nominal_interface_jump_rms\": " << nominalResult.error.interfaceJump << ",\n"
        << "  \"nominal_relative_heat_flux_imbalance\": "
        << nominalResult.error.heatFluxImbalance << ",\n"
        << "  \"nominal_correction_iterations\": " << nominalResult.correctionIterations << "\n"
        << "}\n";

    std::cout << "[MOR] offline setup time: " << offlineSeconds << " s\n"
              << "[MOR] nominal online time: " << nominalResult.timing.totalSeconds << " s\n"
              << "[MOR] nominal interface/global residual: "
              << nominalResult.error.interfaceRelativeResidual << " / "
              << nominalResult.error.globalRelativeResidual << '\n'
              << "[MOR] correction iterations: " << nominalResult.correctionIterations << '\n';

    WorkflowResult workflow;
    workflow.nominalTemperature = std::move(nominalResult.temperature);
    workflow.selectedRank = selectedModel.rank;
    workflow.snapshotCount = snapshotCount;
    workflow.sourceChannels = static_cast<int>(sources.channels.size());
    workflow.correctionIterations = nominalResult.correctionIterations;
    workflow.setupSeconds = offlineSeconds;
    workflow.nominalOnlineSeconds = nominalResult.timing.totalSeconds;
    workflow.nominalGlobalResidual = nominalResult.error.globalRelativeResidual;
    workflow.nominalInterfaceResidual = nominalResult.error.interfaceRelativeResidual;
    workflow.status = nominalResult.status;
    (void)workflowStart;
    return workflow;
}

} // namespace mor
