#include "transient_workflow.hpp"

#include "block_arnoldi.hpp"
#include "thermal_descriptor_system.hpp"
#include "transient_input_waveform.hpp"
#include "transient_reduced_model.hpp"
#include "transient_rom_diagnostics.hpp"
#include "transient_rom_serialization.hpp"
#include "transient_time_integrator.hpp"

#include "sipg_core.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace mor::transient {
namespace {

using Clock = std::chrono::steady_clock;

PowerWaveform waveformFor(const Options& options,
                          const TransientReducedModel& model,
                          double timeStep,
                          int timeSteps,
                          std::uint64_t seedOffset = 0)
{
    return options.inputPath.empty()
        ? makeBuiltinWaveform(options.waveform, model.nominalPowersW,
            timeStep, timeSteps, options.seed + seedOffset)
        : loadPowerWaveformCsv(options.inputPath, model.sourceChannels);
}

void writePureMaximumCurve(const PowerWaveform& waveform,
                           const ReducedTrajectory& trajectory,
                           const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "waveform,step,time_s,rom_maximum_k\n" << std::setprecision(17);
    for (int step = 0; step <= trajectory.steps; ++step) {
        out << waveform.name << ',' << step << ','
            << trajectory.times[static_cast<std::size_t>(step)] << ','
            << trajectory.maximumTemperature[static_cast<std::size_t>(step)] << '\n';
    }
}

void writeFinalField(const TransientReducedModel& model,
                     const ReducedTrajectory& trajectory,
                     const std::string& outputMode,
                     const std::filesystem::path& outputDirectory)
{
    if (outputMode != "full-field" && outputMode != "selected-dofs"
        && outputMode != "selected-regions") return;
    std::vector<double> coordinates(static_cast<std::size_t>(model.rank), 0.0);
    std::copy_n(trajectory.states.end() - model.rank,
                static_cast<std::size_t>(model.rank), coordinates.begin());
    if (outputMode == "selected-regions") {
        const std::vector<double> temperature = reconstructTemperature(model, coordinates);
        int maximumSubdomain = -1;
        for (const DeploymentDof& dof : model.deploymentDofs) {
            maximumSubdomain = std::max(maximumSubdomain, dof.subdomain);
        }
        std::vector<double> sum(static_cast<std::size_t>(maximumSubdomain + 1), 0.0);
        std::vector<double> maximum(static_cast<std::size_t>(maximumSubdomain + 1),
                                    -std::numeric_limits<double>::infinity());
        std::vector<std::size_t> count(static_cast<std::size_t>(maximumSubdomain + 1), 0);
        for (std::size_t row = 0; row < temperature.size(); ++row) {
            const std::size_t region = static_cast<std::size_t>(
                model.deploymentDofs[row].subdomain);
            sum[region] += temperature[row];
            maximum[region] = std::max(maximum[region], temperature[row]);
            ++count[region];
        }
        std::ofstream out(outputDirectory / "transient_selected_regions.csv");
        out << "subdomain,dof_count,average_temperature_k,maximum_temperature_k\n"
            << std::setprecision(17);
        for (int region = 0; region <= maximumSubdomain; ++region) {
            const std::size_t entries = count[static_cast<std::size_t>(region)];
            if (entries == 0) continue;
            out << region << ',' << entries << ','
                << sum[static_cast<std::size_t>(region)] / static_cast<double>(entries)
                << ',' << maximum[static_cast<std::size_t>(region)] << '\n';
        }
        return;
    }
    std::ofstream out(outputDirectory / (outputMode == "full-field"
        ? "transient_final_temperature.csv" : "transient_selected_dofs.csv"));
    out << "dof,x_m,y_m,z_m,subdomain,temperature_k\n" << std::setprecision(17);
    const std::size_t count = outputMode == "selected-dofs"
        ? std::min<std::size_t>(16, model.deploymentDofs.size())
        : model.deploymentDofs.size();
    for (std::size_t row = 0; row < count; ++row) {
        const DeploymentDof& dof = model.deploymentDofs[row];
        double temperature = model.referenceTemperature[row];
        for (int mode = 0; mode < model.rank; ++mode) {
            temperature += model.basis[row + static_cast<std::size_t>(mode)
                * model.globalDofs] * coordinates[static_cast<std::size_t>(mode)];
        }
        out << row << ',' << dof.x << ',' << dof.y << ',' << dof.z << ','
            << dof.subdomain << ',' << temperature << '\n';
    }
}

void writeOnlineTiming(const TransientOnlineTiming& timing,
                       int cases,
                       double averageSeconds,
                       const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "cases,load_seconds,reduced_factorization_seconds,time_stepping_seconds,"
        << "reconstruction_seconds,output_seconds,total_seconds,average_seconds,"
        << "peak_working_set_bytes,ldlt_fallback\n" << std::setprecision(17)
        << cases << ',' << timing.loadSeconds << ','
        << timing.reducedFactorizationSeconds << ',' << timing.timeSteppingSeconds << ','
        << timing.reconstructionSeconds << ',' << timing.outputSeconds << ','
        << timing.totalSeconds << ',' << averageSeconds << ','
        << timing.peakWorkingSetBytes << ',' << (timing.usedLdltFallback ? 1 : 0) << '\n';
}

void writeModelSummary(const TransientReducedModel& model,
                       const TransientAccuracySummary* accuracy,
                       const std::vector<DcConsistencyRow>* dc,
                       double offlineSeconds,
                       std::uint64_t modelBytes,
                       const std::filesystem::path& path)
{
    const double dcWorst = dc == nullptr ? 0.0
        : std::max_element(dc->begin(), dc->end(), [](const auto& a, const auto& b) {
            return a.relativeL2 < b.relativeL2;
        })->relativeL2;
    std::ofstream out(path);
    out << "global_dofs,source_channels,moments,rank,block_size,expansion_point,"
        << "basis_orthogonality_error,capacity_min_eigenvalue,conductivity_min_eigenvalue,"
        << "reference_residual,dc_worst_relative_l2,space_time_relative_l2,"
        << "maximum_absolute_k,speedup,offline_seconds,model_bytes,status\n"
        << std::setprecision(17)
        << model.globalDofs << ',' << model.sourceChannels << ',' << model.moments << ','
        << model.rank << ',' << model.blockSize << ',' << model.expansionPoint << ','
        << model.basisOrthogonalityError << ','
        << model.capacityDiagnostic.minimumEigenvalue << ','
        << model.conductivityDiagnostic.minimumEigenvalue << ','
        << model.referenceResidual << ',' << dcWorst << ','
        << (accuracy == nullptr ? 0.0 : accuracy->spaceTimeRelativeL2) << ','
        << (accuracy == nullptr ? 0.0 : accuracy->maximumAbsolute) << ','
        << (accuracy == nullptr ? 0.0 : accuracy->speedup) << ','
        << offlineSeconds << ',' << modelBytes << ','
        << (accuracy == nullptr ? "not_compared"
            : ((accuracy->spaceTimeRelativeL2 < 1.0e-4
                && accuracy->maximumAbsolute < 0.1) ? "success" : "accuracy_failed"))
        << '\n';
}

void writeOfflineTiming(const ThermalDescriptorSystem& descriptor,
                        const TransientReducedModel& model,
                        double serialization,
                        double total,
                        const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "descriptor_assembly_seconds,symbolic_analysis_seconds,"
        << "numerical_factorization_seconds,multi_rhs_solve_seconds,"
        << "orthogonalization_seconds,reduced_projection_seconds,"
        << "serialization_seconds,total_offline_seconds,peak_working_set_bytes,"
        << "symbolic_calls,numerical_calls\n" << std::setprecision(17)
        << descriptor.assemblySeconds << ','
        << model.arnoldiTiming.symbolicAnalysisSeconds << ','
        << model.arnoldiTiming.numericalFactorizationSeconds << ','
        << model.arnoldiTiming.multiRhsSolveSeconds << ','
        << model.arnoldiTiming.orthogonalizationSeconds << ','
        << model.projectionSeconds << ',' << serialization << ',' << total << ','
        << peakWorkingSetBytes() << ','
        << model.arnoldiTiming.symbolicAnalysisCalls << ','
        << model.arnoldiTiming.numericalFactorizationCalls << '\n';
}

std::uint64_t residentModelBytes(const TransientReducedModel& model)
{
    return static_cast<std::uint64_t>(
        (model.basis.size() + model.reducedCapacity.size()
         + model.reducedConductivity.size() + model.reducedInput.size()
         + model.reducedBoundary.size() + model.referenceTemperature.size()
         + model.nominalPowersW.size() + model.minimumPowersW.size()
         + model.maximumPowersW.size()) * sizeof(double)
        + (model.sourceSubdomains.size() + model.sourceDomainEntities.size())
            * sizeof(int)
        + model.deploymentDofs.size() * sizeof(DeploymentDof));
}

} // namespace

void runTransientBlockArnoldiWorkflow(
    const Mesh& mesh,
    const CaseConfig& physics,
    const Options& options,
    const std::filesystem::path& outputDirectory)
{
    std::filesystem::create_directories(outputDirectory);
    if (options.outputMode != "full-field" && options.outputMode != "max-temperature"
        && options.outputMode != "selected-dofs"
        && options.outputMode != "selected-regions") {
        throw std::runtime_error("Unknown transient output mode: " + options.outputMode);
    }
    const auto offlineStart = Clock::now();
    ThermalDescriptorSystem descriptor = assembleThermalDescriptorSystem(
        mesh, physics, options.massType);
    writeTransientMatrixDiagnostics(
        descriptor, outputDirectory / "transient_matrix_diagnostics.csv");

    TransientReducedModel model;
    double serializationSeconds = 0.0;
    std::filesystem::path modelPath = options.savePath;
    if (options.generate) {
        BlockArnoldiResult arnoldi = buildBlockArnoldiBasis(
            descriptor, options.moments, options.expansionPoint,
            options.rankTolerance, options.secondMomentEnergy,
            options.secondMomentMaximumColumns);
        model = buildTransientReducedModel(descriptor, std::move(arnoldi));
        if (!modelPath.empty()) {
            const auto serializationStart = Clock::now();
            saveTransientReducedModel(model, modelPath);
            serializationSeconds = std::chrono::duration<double>(
                Clock::now() - serializationStart).count();
        }
    } else {
        if (options.loadPath.empty()) {
            throw std::runtime_error(
                "Transient workflow requires --mor-transient-generate or --mor-transient-load.");
        }
        modelPath = options.loadPath;
        model = loadTransientReducedModel(modelPath, &descriptor.fingerprints);
    }
    const double offlineSeconds = std::chrono::duration<double>(
        Clock::now() - offlineStart).count();
    writeArnoldiHistory(model, outputDirectory / "transient_arnoldi_rank_history.csv");
    writeOfflineTiming(descriptor, model, serializationSeconds, offlineSeconds,
                       outputDirectory / "transient_offline_timing.csv");

    const double dt = std::isfinite(options.timeStep)
        ? options.timeStep : physics.timeStep;
    const double endTime = std::isfinite(options.endTime)
        ? options.endTime : dt * std::max(1, physics.timeSteps);
    if (!(dt > 0.0) || !(endTime > 0.0)) {
        throw std::runtime_error(
            "Transient time step and end time must be positive.");
    }
    const int timeSteps = std::max(1, static_cast<int>(std::llround(endTime / dt)));
    const PowerWaveform waveform = waveformFor(options, model, dt, timeSteps);
    std::vector<double> initial = model.referenceTemperature;
    if (std::isfinite(options.initialTemperature)) {
        std::fill(initial.begin(), initial.end(), options.initialTemperature);
    }
    const std::vector<double> initialCoordinates =
        projectInitialConditionCWeighted(descriptor, model, initial);
    const double projectionError = initialProjectionOrthogonalityError(
        descriptor, model, initial, initialCoordinates);
    {
        std::ofstream out(outputDirectory / "transient_initial_projection.csv");
        out << "projection,c_weighted_orthogonality_error,initial_temperature_k\n"
            << std::setprecision(17) << "C_weighted," << projectionError << ','
            << (std::isfinite(options.initialTemperature)
                ? options.initialTemperature : model.referenceTemperature.front()) << '\n';
    }
    const bool computeMaximumTemperature = options.compareFom
        || options.outputMode == "max-temperature"
        || options.outputMode == "full-field";
    ReducedTrajectory trajectory = integrateReducedModel(
        model, waveform, dt, timeSteps, options.integrator,
        initialCoordinates, computeMaximumTemperature);
    const auto outputStart = Clock::now();
    if (computeMaximumTemperature) {
        writePureMaximumCurve(waveform, trajectory,
                              outputDirectory / "transient_max_temperature_curves.csv");
    }
    writeFinalField(model, trajectory, options.outputMode, outputDirectory);
    trajectory.timing.outputSeconds = std::chrono::duration<double>(
        Clock::now() - outputStart).count();
    trajectory.timing.totalSeconds += trajectory.timing.outputSeconds;
    writeOnlineTiming(trajectory.timing, 1, trajectory.timing.totalSeconds,
                      outputDirectory / "transient_online_timing.csv");

    const std::vector<DcConsistencyRow> dc =
        evaluateDcConsistency(descriptor, model);
    writeDcConsistency(dc, outputDirectory / "transient_dc_consistency.csv");
    std::vector<TransientAccuracySummary> accuracyRows;
    if (options.compareFom) {
        accuracyRows.push_back(compareWithFullOrder(
            descriptor, model, waveform, trajectory, dt, timeSteps,
            options.integrator, initial));
        writeAccuracyByTime(accuracyRows.front(),
                            outputDirectory / "transient_accuracy_by_time.csv");
        writeAccuracyByWaveform(accuracyRows,
                            outputDirectory / "transient_accuracy_by_waveform.csv");
        writeMaximumTemperatureCurves(accuracyRows.front(),
                            outputDirectory / "transient_max_temperature_curves.csv");
    }
    const std::uint64_t modelBytes = modelPath.empty()
        ? residentModelBytes(model) : transientModelFileBytes(modelPath);
    writeModelSummary(model,
        accuracyRows.empty() ? nullptr : &accuracyRows.front(), &dc,
        offlineSeconds, modelBytes,
        outputDirectory / "transient_block_arnoldi_summary.csv");
    {
        std::ofstream memory(outputDirectory / "transient_memory.csv");
        memory << "model_bytes,basis_bytes,peak_working_set_bytes,fom_factor_memory_bytes\n"
            << modelBytes << ',' << model.basis.size() * sizeof(double) << ','
            << peakWorkingSetBytes() << ','
            << (accuracyRows.empty() ? 0 : accuracyRows.front().fomFactorMemoryBytes) << '\n';
    }
    {
        const double difference = accuracyRows.empty() ? 0.0
            : accuracyRows.front().fomTimeSteppingSeconds
                - accuracyRows.front().romTimeSteppingSeconds;
        std::ofstream out(outputDirectory / "transient_break_even.csv");
        out << "offline_seconds,fom_seconds_per_waveform,rom_seconds_per_waveform,"
            << "break_even_waveforms\n" << std::setprecision(17)
            << offlineSeconds << ','
            << (accuracyRows.empty() ? 0.0 : accuracyRows.front().fomTimeSteppingSeconds)
            << ',' << trajectory.timing.timeSteppingSeconds << ','
            << (difference > 0.0 ? offlineSeconds / difference : -1.0) << '\n';
    }
    std::cout << "[Stage 2C.1] dofs=" << model.globalDofs
              << ", channels=" << model.sourceChannels
              << ", moments=" << model.moments << ", rank=" << model.rank
              << ", orthogonality=" << model.basisOrthogonalityError
              << ", online=" << trajectory.timing.totalSeconds << " s\n";
}

void runTransientDeploymentOnly(
    const std::filesystem::path& modelDirectory,
    const Options& options,
    const std::filesystem::path& outputDirectory)
{
    std::filesystem::create_directories(outputDirectory);
    if (options.outputMode != "full-field" && options.outputMode != "max-temperature"
        && options.outputMode != "selected-dofs"
        && options.outputMode != "selected-regions") {
        throw std::runtime_error("Unknown transient output mode: " + options.outputMode);
    }
    const auto totalStart = Clock::now();
    const auto loadStart = Clock::now();
    const TransientReducedModel model =
        loadTransientReducedModel(modelDirectory, nullptr);
    const double loadSeconds = std::chrono::duration<double>(
        Clock::now() - loadStart).count();
    if (std::isfinite(options.initialTemperature)) {
        throw std::runtime_error(
            "Pure transient deployment supports the serialized reference initial state; "
            "C-weighted nonzero initial projection requires the descriptor workflow.");
    }
    const double dt = std::isfinite(options.timeStep) ? options.timeStep : 1.0;
    const double endTime = std::isfinite(options.endTime) ? options.endTime : 10.0 * dt;
    if (!(dt > 0.0) || !(endTime > 0.0)
        || options.deploymentRhsCount <= 0) {
        throw std::runtime_error(
            "Pure transient deployment requires positive dt, end time, and case count.");
    }
    const int timeSteps = std::max(1, static_cast<int>(std::llround(endTime / dt)));
    std::vector<double> initial(static_cast<std::size_t>(model.rank), 0.0);
    TransientOnlineTiming aggregate;
    ReducedTrajectory representative;
    const bool computeMaximumTemperature = options.outputMode == "max-temperature"
        || options.outputMode == "full-field";
    for (int caseIndex = 0; caseIndex < options.deploymentRhsCount; ++caseIndex) {
        const PowerWaveform waveform = waveformFor(
            options, model, dt, timeSteps, static_cast<std::uint64_t>(caseIndex) * 7919ULL);
        ReducedTrajectory trajectory = integrateReducedModel(
            model, waveform, dt, timeSteps, options.integrator, initial,
            computeMaximumTemperature);
        aggregate.reducedFactorizationSeconds += trajectory.timing.reducedFactorizationSeconds;
        aggregate.timeSteppingSeconds += trajectory.timing.timeSteppingSeconds;
        aggregate.reconstructionSeconds += trajectory.timing.reconstructionSeconds;
        aggregate.totalSeconds += trajectory.timing.totalSeconds;
        aggregate.usedLdltFallback = aggregate.usedLdltFallback
            || trajectory.timing.usedLdltFallback;
        if (caseIndex == 0) representative = std::move(trajectory);
    }
    const auto outputStart = Clock::now();
    const PowerWaveform representativeWaveform = waveformFor(
        options, model, dt, timeSteps, 0);
    if (computeMaximumTemperature) {
        writePureMaximumCurve(representativeWaveform, representative,
                              outputDirectory / "transient_max_temperature_curves.csv");
    }
    writeFinalField(model, representative, options.outputMode, outputDirectory);
    aggregate.outputSeconds = std::chrono::duration<double>(
        Clock::now() - outputStart).count();
    aggregate.loadSeconds = loadSeconds;
    aggregate.totalSeconds = std::chrono::duration<double>(
        Clock::now() - totalStart).count();
    aggregate.peakWorkingSetBytes = peakWorkingSetBytes();
    writeOnlineTiming(aggregate, options.deploymentRhsCount,
        aggregate.totalSeconds / static_cast<double>(options.deploymentRhsCount),
        outputDirectory / "transient_online_timing.csv");
    std::cout << "[Stage 2C.1 pure deployment] rank=" << model.rank
              << ", cases=" << options.deploymentRhsCount
              << ", average=" << aggregate.totalSeconds
                    / static_cast<double>(options.deploymentRhsCount)
              << " s, peak=" << aggregate.peakWorkingSetBytes << " bytes\n";
}

} // namespace mor::transient
