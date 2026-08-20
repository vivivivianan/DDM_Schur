#include "parametric_workflow.hpp"

#include "affine_fem_components.hpp"
#include "parametric_reduced_schur.hpp"
#include "ddm_schur/schur_fgmres.hpp"
#include "ddm_schur/interface_operator.hpp"
#include "linear_solvers.hpp"
#include "config_io.hpp"
#include "mesh_loader.hpp"
#include "fem_assembly.hpp"
#include "diagnostics_io.hpp"
#include "mor/model_io.hpp"
#include "mor/mor_diagnostics.hpp"
#include "sipg_core.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>

namespace mor::parametric {
namespace {

using Clock = std::chrono::steady_clock;

double elapsed(const Clock::time_point& start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

struct ProcessMemorySnapshot {
    std::uint64_t workingSetBytes = 0;
    std::uint64_t privateBytes = 0;
    std::uint64_t committedBytes = 0;
};

ProcessMemorySnapshot processMemorySnapshot()
{
    ProcessMemorySnapshot snapshot;
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX counters{};
    counters.cb = sizeof(counters);
    if (GetProcessMemoryInfo(GetCurrentProcess(),
            reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&counters),
            sizeof(counters)) != 0) {
        snapshot.workingSetBytes = counters.WorkingSetSize;
        snapshot.privateBytes = counters.PrivateUsage;
        snapshot.committedBytes = counters.PagefileUsage;
    }
#endif
    return snapshot;
}

double mib(std::uint64_t bytes)
{
    return static_cast<double>(bytes) / (1024.0 * 1024.0);
}

void writeMemoryTimelineHeader(std::ostream& out)
{
    out << "stage,case_index,working_set_mib,private_bytes_mib,committed_bytes_mib,"
        << "model_resident_mib,workspace_mib,reconstruction_mib,result_cache_mib,"
        << "io_buffer_mib,allocator_retained_mib\n" << std::setprecision(17);
}

void writeMemoryTimelineRow(
    std::ostream& out,
    const std::string& stage,
    int caseIndex,
    const ParametricModelMemoryBreakdown* modelMemory,
    const ParametricRomWorkspace* workspace)
{
    const ProcessMemorySnapshot process = processMemorySnapshot();
    const std::uint64_t modelBytes = modelMemory == nullptr
        ? 0U : static_cast<std::uint64_t>(modelMemory->totalBytes());
    const std::uint64_t workspaceBytes = workspace == nullptr
        ? 0U : static_cast<std::uint64_t>(workspace->workspaceBytes());
    const std::uint64_t reconstructionBytes = workspace == nullptr
        ? 0U : static_cast<std::uint64_t>(workspace->reconstructionBytes());
    const std::uint64_t ioBytes = workspace == nullptr
        ? 0U : static_cast<std::uint64_t>(
            workspace->outputScratch.capacity() * sizeof(char));
    const std::uint64_t categorized = modelBytes + workspaceBytes
        + reconstructionBytes + ioBytes;
    const std::uint64_t allocatorRetained = process.privateBytes > categorized
        ? process.privateBytes - categorized : 0U;
    out << stage << ',' << caseIndex << ',' << mib(process.workingSetBytes) << ','
        << mib(process.privateBytes) << ',' << mib(process.committedBytes) << ','
        << mib(modelBytes) << ',' << mib(workspaceBytes) << ','
        << mib(reconstructionBytes) << ",0," << mib(ioBytes) << ','
        << mib(allocatorRetained) << '\n';
    out.flush();
}

double norm(const std::vector<double>& vector)
{
    double squared = 0.0;
    for (double value : vector) {
        squared += value * value;
    }
    return std::sqrt(squared);
}

struct ResidualRatios {
    double global = 0.0;
    double interfaceEquation = 0.0;
};

ResidualRatios calculateAffineResidual(
    const AffineFemComponents& affine,
    double parameterValue,
    const std::vector<double>& powers,
    const std::vector<double>& temperature,
    const std::vector<int>& interfaceGlobalDofs)
{
    std::vector<double> image = affine.matrixConstant.multiply(temperature);
    const std::vector<double> linearImage = affine.matrixLinear.multiply(temperature);
    for (std::size_t row = 0; row < image.size(); ++row) {
        image[row] += parameterValue * linearImage[row];
    }
    for (std::size_t group = 0; group < affine.matrixHarmonic.size(); ++group) {
        const std::vector<double> harmonicImage =
            affine.matrixHarmonic[group].multiply(temperature);
        const double coefficient = harmonicTheta(affine.parameter, parameterValue, group);
        for (std::size_t row = 0; row < image.size(); ++row) {
            image[row] += coefficient * harmonicImage[row];
        }
    }
    const std::vector<double> rhs = composeRhs(affine, parameterValue, powers);
    std::vector<double> residual(rhs.size(), 0.0);
    for (std::size_t row = 0; row < rhs.size(); ++row) {
        residual[row] = rhs[row] - image[row];
    }
    double interfaceResidualSquared = 0.0;
    double interfaceRhsSquared = 0.0;
    for (int global : interfaceGlobalDofs) {
        const double residualValue = residual[static_cast<std::size_t>(global)];
        const double rhsValue = rhs[static_cast<std::size_t>(global)];
        interfaceResidualSquared += residualValue * residualValue;
        interfaceRhsSquared += rhsValue * rhsValue;
    }
    ResidualRatios ratios;
    ratios.global = norm(residual) / std::max(1.0e-300, norm(rhs));
    ratios.interfaceEquation = std::sqrt(interfaceResidualSquared)
        / std::max(1.0e-300, std::sqrt(interfaceRhsSquared));
    return ratios;
}

std::vector<double> nominalPowers(const SourceParameterization& sources)
{
    std::vector<double> result;
    result.reserve(sources.channels.size());
    for (const SourceChannel& channel : sources.channels) {
        result.push_back(channel.nominalPowerW);
    }
    return result;
}

void appendUnique(std::vector<double>& values, double value)
{
    for (double existing : values) {
        if (std::abs(existing - value)
            <= 1.0e-13 * std::max({1.0, std::abs(existing), std::abs(value)})) {
            return;
        }
    }
    values.push_back(value);
}

std::vector<double> trainingParameters(const AffineParameter& parameter, int count)
{
    if (count < 3) {
        throw std::runtime_error("Stage 2B.1 requires at least three training parameter points.");
    }
    std::vector<double> result;
    appendUnique(result, parameter.reference);
    appendUnique(result, parameter.minimum);
    appendUnique(result, parameter.maximum);
    while (static_cast<int>(result.size()) < count) {
        std::vector<double> sorted = result;
        std::sort(sorted.begin(), sorted.end());
        double widestGap = -1.0;
        double midpoint = parameter.reference;
        for (std::size_t index = 1; index < sorted.size(); ++index) {
            const double gap = sorted[index] - sorted[index - 1];
            if (gap > widestGap) {
                widestGap = gap;
                midpoint = 0.5 * (sorted[index - 1] + sorted[index]);
            }
        }
        appendUnique(result, midpoint);
    }
    return result;
}

std::vector<double> validationParameters(const AffineParameter& parameter,
                                         const std::vector<double>& training,
                                         int count)
{
    std::vector<double> sorted = training;
    std::sort(sorted.begin(), sorted.end());
    std::vector<double> result;
    for (std::size_t index = 1; index < sorted.size()
         && static_cast<int>(result.size()) < count; ++index) {
        appendUnique(result, 0.5 * (sorted[index - 1] + sorted[index]));
    }
    for (int index = 0; static_cast<int>(result.size()) < count; ++index) {
        const double fraction = (static_cast<double>(index) + 0.5)
            / static_cast<double>(count);
        appendUnique(result, parameter.minimum
            + fraction * (parameter.maximum - parameter.minimum));
    }
    return result;
}

std::vector<double> testParameters(const AffineParameter& parameter,
                                   int count,
                                   std::uint64_t seed)
{
    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<double> uniform(
        parameter.minimum, parameter.maximum);
    std::vector<double> result;
    while (static_cast<int>(result.size()) < count) {
        appendUnique(result, uniform(generator));
    }
    return result;
}

void writeParameterPoints(const std::vector<double>& values,
                          const std::string& split,
                          const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "index,split,parameter_value\n" << std::setprecision(17);
    for (std::size_t index = 0; index < values.size(); ++index) {
        out << index << ',' << split << ',' << values[index] << '\n';
    }
}

struct OfflineTiming {
    double affineAssembly = 0.0;
    double fomAssembly = 0.0;
    double fomSymbolic = 0.0;
    double fomNumerical = 0.0;
    double multiRhsSolve = 0.0;
    double snapshotStorage = 0.0;
    double podQr = 0.0;
    double projectionPreparation = 0.0;
    double interfaceProjection = 0.0;
    double localProjection = 0.0;
    double serialization = 0.0;
    double total = 0.0;
    std::size_t peakMemory = 0;
};

void appendSnapshotCase(SnapshotDatabase& database,
                        const std::string& family,
                        const double* values)
{
    ParameterCase parameterCase;
    parameterCase.index = static_cast<int>(database.cases.size());
    parameterCase.split = "training";
    parameterCase.family = family;
    database.cases.push_back(std::move(parameterCase));
    database.values.insert(database.values.end(), values,
                           values + database.rows);
}

ParametricTrainingSnapshots generateSnapshots(
    const Mesh& mesh,
    const CaseConfig& physics,
    const AffineFemComponents& affine,
    const ddm_schur::InterfacePartition& partition,
    const std::vector<double>& parameters,
    OfflineTiming& timing)
{
    if (parameters.empty()
        || std::abs(parameters.front() - affine.parameter.reference) > 1.0e-12) {
        throw std::runtime_error("The reference parameter must be the first training point.");
    }
    ParametricTrainingSnapshots snapshots;
    snapshots.interfaceSnapshots.rows = static_cast<int>(partition.interfaceGlobalDofs.size());
    snapshots.interiors.resize(partition.domains.size());
    for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
        snapshots.interiors[slot].rows = static_cast<int>(
            partition.domains[slot].interiorGlobalDofs.size());
    }
    const int size = static_cast<int>(mesh.nodes.size());
    const int channels = static_cast<int>(affine.sources.channels.size());
    const std::vector<double> zeroPowers(static_cast<std::size_t>(channels), 0.0);

    for (std::size_t parameterIndex = 0;
         parameterIndex < parameters.size(); ++parameterIndex) {
        const double parameterValue = parameters[parameterIndex];
        const auto assemblyStart = Clock::now();
        DirectParametricSystem direct = assembleDirectParametricSystem(
            mesh, physics, affine.parameter, parameterValue);
        timing.fomAssembly += elapsed(assemblyStart);

        const auto factorStart = Clock::now();
        SubdomainDirectSolver solver(
            direct.matrix.size(), sparseMatrixEntries(direct.matrix));
        (void)factorStart;
        timing.fomSymbolic += solver.symbolicAnalysisSeconds();
        timing.fomNumerical += solver.numericalFactorizationSeconds();

        const std::vector<double> baseRhs = composeRhs(
            affine, parameterValue, zeroPowers);
        std::vector<double> rhs(static_cast<std::size_t>(size * (channels + 1)), 0.0);
        std::copy(baseRhs.begin(), baseRhs.end(), rhs.begin());
        for (int channel = 0; channel < channels; ++channel) {
            const std::size_t offset = static_cast<std::size_t>(channel + 1)
                * static_cast<std::size_t>(size);
            std::copy(baseRhs.begin(), baseRhs.end(), rhs.begin() + offset);
            const std::vector<double>& load = affine.sources.channels[
                static_cast<std::size_t>(channel)].rhsPerWatt;
            for (int row = 0; row < size; ++row) {
                rhs[offset + static_cast<std::size_t>(row)] +=
                    load[static_cast<std::size_t>(row)];
            }
        }
        const auto solveStart = Clock::now();
        std::vector<double> solution;
        solver.solveMultiple(rhs, channels + 1, solution);
        timing.multiRhsSolve += elapsed(solveStart);

        if (parameterIndex == 0) {
            snapshots.referenceTemperature.assign(
                solution.begin(), solution.begin() + size);
        }
        const auto storeStart = Clock::now();
        std::vector<double> globalColumn(static_cast<std::size_t>(size), 0.0);
        for (int column = 0; column <= channels; ++column) {
            const std::size_t offset = static_cast<std::size_t>(column)
                * static_cast<std::size_t>(size);
            for (int row = 0; row < size; ++row) {
                globalColumn[static_cast<std::size_t>(row)] = column == 0
                    ? solution[offset + static_cast<std::size_t>(row)]
                        - snapshots.referenceTemperature[static_cast<std::size_t>(row)]
                    : solution[offset + static_cast<std::size_t>(row)]
                        - solution[static_cast<std::size_t>(row)];
            }
            std::vector<double> interfaceColumn(
                partition.interfaceGlobalDofs.size(), 0.0);
            for (std::size_t row = 0;
                 row < partition.interfaceGlobalDofs.size(); ++row) {
                interfaceColumn[row] = globalColumn[static_cast<std::size_t>(
                    partition.interfaceGlobalDofs[row])];
            }
            appendSnapshotCase(snapshots.interfaceSnapshots,
                column == 0 ? "parameter_reference" : "unit_power_channel",
                interfaceColumn.data());
            for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
                const std::vector<int>& dofs =
                    partition.domains[slot].interiorGlobalDofs;
                std::vector<double> localColumn(dofs.size(), 0.0);
                for (std::size_t row = 0; row < dofs.size(); ++row) {
                    localColumn[row] = globalColumn[static_cast<std::size_t>(dofs[row])];
                }
                appendSnapshotCase(snapshots.interiors[slot],
                    column == 0 ? "parameter_reference" : "unit_power_channel",
                    localColumn.data());
            }
        }
        timing.snapshotStorage += elapsed(storeStart);
        std::cout << "[Stage 2B.1] snapshots mu=" << parameterValue
                  << ", multi-RHS=" << channels + 1
                  << ", symbolic=" << solver.symbolicAnalysisSeconds()
                  << " s, numerical=" << solver.numericalFactorizationSeconds()
                  << " s\n";
    }
    return snapshots;
}

struct AccuracyRow {
    std::string split;
    int parameterIndex = 0;
    int powerCase = 0;
    std::string family;
    int powerChannel = -1;
    double parameterValue = 0.0;
    double interfaceRelativeL2 = 0.0;
    double relativeL2 = 0.0;
    double meanAbsolute = 0.0;
    double maximumAbsolute = 0.0;
    double maximumTemperatureError = 0.0;
    double interfaceJump = 0.0;
    double fluxImbalance = 0.0;
    double globalResidual = 0.0;
    double fomSeconds = 0.0;
    ParametricOnlineTiming timing;
};

AccuracyRow calculateAccuracy(const Mesh& mesh,
                              const CaseConfig& physics,
                              const DirectParametricSystem& direct,
                              const std::vector<double>& rhs,
                              const std::vector<double>& truth,
                              const ParametricReducedModel& model,
                              const ParametricOnlineResult& reduced)
{
    AccuracyRow row;
    double errorSquared = 0.0;
    double truthSquared = 0.0;
    double absoluteSum = 0.0;
    for (std::size_t index = 0; index < truth.size(); ++index) {
        const double error = reduced.temperature[index] - truth[index];
        errorSquared += error * error;
        truthSquared += truth[index] * truth[index];
        absoluteSum += std::abs(error);
        row.maximumAbsolute = std::max(row.maximumAbsolute, std::abs(error));
    }
    row.relativeL2 = std::sqrt(errorSquared)
        / std::max(1.0e-300, std::sqrt(truthSquared));
    row.meanAbsolute = absoluteSum / static_cast<double>(truth.size());
    row.maximumTemperatureError = std::abs(
        *std::max_element(reduced.temperature.begin(), reduced.temperature.end())
        - *std::max_element(truth.begin(), truth.end()));
    double interfaceError = 0.0;
    double interfaceTruth = 0.0;
    for (int global : model.interfaceGlobalDofs) {
        const double error = reduced.temperature[static_cast<std::size_t>(global)]
            - truth[static_cast<std::size_t>(global)];
        interfaceError += error * error;
        interfaceTruth += truth[static_cast<std::size_t>(global)]
            * truth[static_cast<std::size_t>(global)];
    }
    row.interfaceRelativeL2 = std::sqrt(interfaceError)
        / std::max(1.0e-300, std::sqrt(interfaceTruth));
    const std::vector<double> image = direct.matrix.multiply(reduced.temperature);
    std::vector<double> residual(rhs.size(), 0.0);
    for (std::size_t index = 0; index < rhs.size(); ++index) {
        residual[index] = rhs[index] - image[index];
    }
    row.globalResidual = norm(residual) / std::max(1.0e-300, norm(rhs));
    const InterfacePhysicsMetrics interfacePhysics =
        calculateInterfacePhysicsMetrics(mesh, physics, reduced.temperature);
    row.interfaceJump = interfacePhysics.temperatureJumpRms;
    row.fluxImbalance = interfacePhysics.relativeFluxImbalance;
    row.timing = reduced.timing;
    return row;
}

std::vector<std::vector<double>> validationPowerCases(
    const SourceParameterization& sources,
    bool allUnitChannels,
    std::uint64_t seed,
    std::vector<std::string>& families,
    std::vector<int>& powerChannels)
{
    const int channels = static_cast<int>(sources.channels.size());
    std::vector<std::vector<double>> result;
    result.push_back(nominalPowers(sources));
    families.push_back("nominal");
    powerChannels.push_back(-1);
    if (allUnitChannels) {
        for (int channel = 0; channel < channels; ++channel) {
            std::vector<double> powers(static_cast<std::size_t>(channels), 0.0);
            powers[static_cast<std::size_t>(channel)] =
                sources.channels[static_cast<std::size_t>(channel)].nominalPowerW;
            result.push_back(std::move(powers));
            families.push_back("single_power_channel");
            powerChannels.push_back(channel);
        }
    }
    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.5);
    for (int sample = 0; sample < 2; ++sample) {
        std::vector<double> powers(static_cast<std::size_t>(channels), 0.0);
        for (int channel = 0; channel < channels; ++channel) {
            powers[static_cast<std::size_t>(channel)] = uniform(generator)
                * sources.channels[static_cast<std::size_t>(channel)].nominalPowerW;
        }
        result.push_back(std::move(powers));
        families.push_back("seeded_multi_hotspot");
        powerChannels.push_back(-1);
    }
    return result;
}

void evaluateParameterSet(const Mesh& mesh,
                          const CaseConfig& physics,
                          const AffineFemComponents& affine,
                          const ParametricReducedModel& model,
                          const std::vector<double>& parameters,
                          const std::string& split,
                          std::uint64_t seed,
                          std::vector<AccuracyRow>& rows)
{
    const int size = static_cast<int>(mesh.nodes.size());
    for (std::size_t parameterIndex = 0;
         parameterIndex < parameters.size(); ++parameterIndex) {
        const double value = parameters[parameterIndex];
        std::vector<std::string> families;
        std::vector<int> channels;
        const std::vector<std::vector<double>> powers = validationPowerCases(
            affine.sources, parameterIndex == 0, seed + parameterIndex,
            families, channels);
        DirectParametricSystem direct = assembleDirectParametricSystem(
            mesh, physics, affine.parameter, value);
        const auto factorStart = Clock::now();
        SubdomainDirectSolver solver(
            direct.matrix.size(), sparseMatrixEntries(direct.matrix));
        std::vector<double> rhs(
            static_cast<std::size_t>(size) * powers.size(), 0.0);
        for (std::size_t column = 0; column < powers.size(); ++column) {
            const std::vector<double> columnRhs = composeRhs(affine, value, powers[column]);
            std::copy(columnRhs.begin(), columnRhs.end(),
                rhs.begin() + static_cast<std::ptrdiff_t>(column * size));
        }
        std::vector<double> truth;
        solver.solveMultiple(rhs, static_cast<int>(powers.size()), truth);
        const double fomSeconds = elapsed(factorStart);
        for (std::size_t column = 0; column < powers.size(); ++column) {
            ParametricOnlineResult reduced = solveParametricRom(
                model, value, powers[column], false);
            std::vector<double> columnTruth(
                truth.begin() + static_cast<std::ptrdiff_t>(column * size),
                truth.begin() + static_cast<std::ptrdiff_t>((column + 1) * size));
            std::vector<double> columnRhs(
                rhs.begin() + static_cast<std::ptrdiff_t>(column * size),
                rhs.begin() + static_cast<std::ptrdiff_t>((column + 1) * size));
            CaseConfig currentPhysics = physicsAtParameter(physics, affine.parameter, value);
            AccuracyRow row = calculateAccuracy(
                mesh, currentPhysics, direct, columnRhs, columnTruth, model, reduced);
            row.split = split;
            row.parameterIndex = static_cast<int>(parameterIndex);
            row.powerCase = static_cast<int>(column);
            row.family = families[column];
            row.powerChannel = channels[column];
            row.parameterValue = value;
            row.fomSeconds = fomSeconds / static_cast<double>(powers.size());
            rows.push_back(std::move(row));
        }
    }
}

void writeAccuracy(const std::vector<AccuracyRow>& rows,
                   const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "split,parameter_index,power_case,family,power_channel,parameter_value,"
        << "interface_relative_l2,full_field_relative_l2,mean_absolute_error,"
        << "maximum_absolute_error,maximum_temperature_error,interface_jump_rms,"
        << "relative_flux_imbalance,full_global_residual,fom_seconds,rom_seconds\n"
        << std::setprecision(17);
    for (const AccuracyRow& row : rows) {
        out << row.split << ',' << row.parameterIndex << ',' << row.powerCase << ','
            << row.family << ',' << row.powerChannel << ',' << row.parameterValue << ','
            << row.interfaceRelativeL2 << ',' << row.relativeL2 << ','
            << row.meanAbsolute << ',' << row.maximumAbsolute << ','
            << row.maximumTemperatureError << ',' << row.interfaceJump << ','
            << row.fluxImbalance << ',' << row.globalResidual << ','
            << row.fomSeconds << ',' << row.timing.totalSeconds << '\n';
    }
}

void writePowerAccuracy(const std::vector<AccuracyRow>& rows,
                        const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "power_channel,parameter_value,full_field_relative_l2,maximum_absolute_error,"
        << "maximum_temperature_error,global_residual\n" << std::setprecision(17);
    for (const AccuracyRow& row : rows) {
        if (row.powerChannel >= 0) {
            out << row.powerChannel << ',' << row.parameterValue << ','
                << row.relativeL2 << ',' << row.maximumAbsolute << ','
                << row.maximumTemperatureError << ',' << row.globalResidual << '\n';
        }
    }
}

void writeParameterAccuracy(const std::vector<AccuracyRow>& rows,
                            const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "split,parameter_value,cases,maximum_relative_l2,maximum_absolute_error,"
        << "maximum_temperature_error,maximum_global_residual\n" << std::setprecision(17);
    for (const std::string split : {std::string("validation"), std::string("test")}) {
        std::vector<double> parameters;
        for (const AccuracyRow& row : rows) {
            if (row.split == split) appendUnique(parameters, row.parameterValue);
        }
        for (double parameter : parameters) {
            int cases = 0;
            double relative = 0.0, absolute = 0.0, maximum = 0.0, residual = 0.0;
            for (const AccuracyRow& row : rows) {
                if (row.split != split || row.parameterValue != parameter) continue;
                ++cases;
                relative = std::max(relative, row.relativeL2);
                absolute = std::max(absolute, row.maximumAbsolute);
                maximum = std::max(maximum, row.maximumTemperatureError);
                residual = std::max(residual, row.globalResidual);
            }
            out << split << ',' << parameter << ',' << cases << ',' << relative << ','
                << absolute << ',' << maximum << ',' << residual << '\n';
        }
    }
}

void writeOfflineTiming(const OfflineTiming& timing,
                        const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "affine_component_assembly_seconds,fom_matrix_assembly_seconds,"
        << "fom_symbolic_seconds,fom_numerical_seconds,multi_rhs_snapshot_solve_seconds,"
        << "snapshot_storage_seconds,pod_qr_seconds,projection_preparation_seconds,"
        << "interface_reduced_projection_seconds,local_reduced_projection_seconds,"
        << "total_reduced_projection_seconds,"
        << "serialization_seconds,total_offline_seconds,peak_memory_bytes\n"
        << std::setprecision(17)
        << timing.affineAssembly << ',' << timing.fomAssembly << ','
        << timing.fomSymbolic << ',' << timing.fomNumerical << ','
        << timing.multiRhsSolve << ',' << timing.snapshotStorage << ','
        << timing.podQr << ',' << timing.projectionPreparation << ','
        << timing.interfaceProjection << ',' << timing.localProjection << ','
        << timing.projectionPreparation + timing.interfaceProjection
            + timing.localProjection << ','
        << timing.serialization << ',' << timing.total << ',' << timing.peakMemory << '\n';
}

void writeOnlineTimingHeader(std::ostream& out);
void writeOnlineTimingRow(std::ostream& out,
                          std::size_t index,
                          const ParametricOnlineResult& row);

void writeOnlineTiming(const std::vector<ParametricOnlineResult>& rows,
                       const std::filesystem::path& path)
{
    std::ofstream out(path);
    writeOnlineTimingHeader(out);
    for (std::size_t index = 0; index < rows.size(); ++index) {
        writeOnlineTimingRow(out, index, rows[index]);
    }
}

void writeOnlineTimingHeader(std::ostream& out)
{
    out << "case,parameter_value,coefficient_seconds,local_assembly_seconds,"
        << "local_factorization_seconds,reduced_schur_assembly_seconds,reduced_solve_seconds,"
        << "reconstruction_seconds,output_io_seconds,total_seconds,extrapolated\n"
        << std::setprecision(17);
}

void writeOnlineTimingRow(std::ostream& out,
                          std::size_t index,
                          const ParametricOnlineResult& row)
{
    out << index << ',' << row.parameterValue << ','
        << row.timing.coefficientSeconds << ',' << row.timing.localAssemblySeconds << ','
        << row.timing.localFactorSeconds << ','
        << row.timing.reducedSchurAssemblySeconds << ','
        << row.timing.reducedSolveSeconds << ','
        << row.timing.reconstructionSeconds << ',' << row.timing.outputSeconds << ','
        << row.timing.totalSeconds << ',' << (row.extrapolated ? 1 : 0) << '\n';
}

void writeTemperature(const ParametricReducedModel& model,
                      const std::vector<double>& temperature,
                      const std::filesystem::path& path)
{
    const auto start = Clock::now();
    std::ofstream out(path);
    out << "dof,x_m,y_m,z_m,subdomain,temperature_k\n" << std::setprecision(17);
    for (std::size_t row = 0; row < temperature.size(); ++row) {
        const DeploymentDof& dof = model.dofs[row];
        out << row << ',' << dof.x << ',' << dof.y << ',' << dof.z << ','
            << dof.subdomain << ',' << temperature[row] << '\n';
    }
    (void)start;
}

void writeRankSweep(const ParametricReducedModel& model,
                    const std::vector<AccuracyRow>& rows,
                    const std::filesystem::path& path)
{
    double worstRelative = 0.0;
    double worstAbsolute = 0.0;
    double worstMaximum = 0.0;
    for (const AccuracyRow& row : rows) {
        worstRelative = std::max(worstRelative, row.relativeL2);
        worstAbsolute = std::max(worstAbsolute, row.maximumAbsolute);
        worstMaximum = std::max(worstMaximum, row.maximumTemperatureError);
    }
    int localMaximum = 0;
    int localTotal = 0;
    for (const ParametricLocalBlock& block : model.locals) {
        localMaximum = std::max(localMaximum, block.rank);
        localTotal += block.rank;
    }
    std::ofstream out(path);
    out << "interface_rank,maximum_local_rank,total_local_rank,maximum_relative_l2,"
        << "maximum_absolute_error,maximum_temperature_error\n"
        << std::setprecision(17) << model.interfaceRank << ',' << localMaximum << ','
        << localTotal << ',' << worstRelative << ',' << worstAbsolute << ','
        << worstMaximum << '\n';
}

} // namespace

WorkflowResult runParametricReducedSchurWorkflow(
    const Mesh& mesh,
    const SparseMatrix& referenceSystem,
    const CaseConfig& physics,
    const std::vector<double>& assembledSource,
    const std::vector<double>& heatOnlySource,
    const std::vector<double>& fixedAdjust,
    const ddm_schur::Options& schurOptions,
    const mor::Options& options,
    const std::filesystem::path& outputDirectory)
{
    if (!options.parametricGenerate && options.parametricLoadPath.empty()) {
        throw std::runtime_error(
            "Stage 2B.1 requires --mor-parametric-generate or --mor-parametric-load <directory>.");
    }
    if (options.parametricMode != "pure" && options.parametricMode != "corrected") {
        throw std::runtime_error("--mor-parametric-mode must be pure or corrected.");
    }
    std::filesystem::create_directories(outputDirectory);
    const auto totalStart = Clock::now();
    OfflineTiming offline;
    const AffineFemComponents affine = buildAffineFemComponents(
        mesh, physics, referenceSystem, assembledSource, heatOnlySource,
        fixedAdjust, options);
    offline.affineAssembly = affine.assemblySeconds;
    const ddm_schur::InterfacePartition partition =
        ddm_schur::buildInterfacePartition(mesh, referenceSystem);
    const Fingerprints fingerprints = computeFingerprints(
        mesh, referenceSystem, physics, partition.interfaceGlobalDofs);

    const std::vector<double> training = trainingParameters(
        affine.parameter, options.parameterTrainingCount);
    const std::vector<double> validation = validationParameters(
        affine.parameter, training, options.parameterValidationCount);
    const std::vector<double> test = testParameters(
        affine.parameter, options.parameterTestCount, options.seed + 424243ULL);
    writeParameterPoints(training, "training", outputDirectory / "training_parameters.csv");
    writeParameterPoints(validation, "validation", outputDirectory / "validation_parameters.csv");
    writeParameterPoints(test, "test", outputDirectory / "test_parameters.csv");

    std::vector<AffineValidationRow> affineValidation;
    const std::vector<double> nominal = nominalPowers(affine.sources);
    for (double value : {affine.parameter.minimum,
                         affine.parameter.reference,
                         affine.parameter.maximum}) {
        affineValidation.push_back(validateAffineOperator(
            mesh, physics, affine, value, nominal));
    }
    writeAffineValidation(affineValidation,
        outputDirectory / "parametric_affine_validation.csv");
    if (options.parametricAffineValidationOnly) {
        WorkflowResult validationOnly;
        validationOnly.nominalTemperature.assign(
            mesh.nodes.size(), physics.initialTemperature);
        validationOnly.sourceChannels = static_cast<int>(affine.sources.channels.size());
        validationOnly.setupSeconds = elapsed(totalStart);
        validationOnly.status = "success";
        return validationOnly;
    }

    ParametricReducedModel model;
    if (options.parametricGenerate) {
        ParametricTrainingSnapshots snapshots = generateSnapshots(
            mesh, physics, affine, partition, training, offline);
        model = buildParametricReducedModel(
            mesh, affine, partition, snapshots, options);
        model.fingerprints = fingerprints;
        model.sourceDefinitionHash = fingerprints.sources;
        offline.podQr = model.basisSeconds;
        offline.projectionPreparation = model.projectionPreparationSeconds;
        offline.interfaceProjection = model.interfaceProjectionSeconds;
        offline.localProjection = model.localProjectionSeconds;
        const std::filesystem::path modelPath = options.parametricSavePath.empty()
            ? outputDirectory / "mor_models" / "stage2b1"
            : options.parametricSavePath;
        const auto serializationStart = Clock::now();
        saveParametricModel(model, modelPath);
        offline.serialization = elapsed(serializationStart);
        model.fileBytes = std::filesystem::file_size(modelPath / "model.bin");
        std::cout << "[Stage 2B.1] model saved to " << modelPath.string() << '\n';
    } else {
        model = loadParametricModel(
            options.parametricLoadPath, &fingerprints, &affine);
        truncateParametricModel(model, options.interfaceRank, options.localRank);
    }
    offline.total = elapsed(totalStart);
    offline.peakMemory = peakWorkingSetBytes();
    writeOfflineTiming(offline, outputDirectory / "parametric_offline_timing.csv");
    {
        std::ofstream memory(outputDirectory / "parametric_memory.csv");
        memory << "model_file_bytes,offline_peak_working_set_bytes\n"
               << model.fileBytes << ',' << offline.peakMemory << '\n';
    }

    std::vector<AccuracyRow> accuracy;
    evaluateParameterSet(mesh, physics, affine, model, validation,
                         "validation", options.seed + 101ULL, accuracy);
    evaluateParameterSet(mesh, physics, affine, model, test,
                         "test", options.seed + 211ULL, accuracy);
    writeAccuracy(accuracy, outputDirectory / "parametric_accuracy_by_case.csv");
    writeParameterAccuracy(accuracy,
        outputDirectory / "parametric_accuracy_by_parameter.csv");
    writePowerAccuracy(accuracy,
        outputDirectory / "parametric_accuracy_by_power_channel.csv");
    writeRankSweep(model, accuracy, outputDirectory / "parametric_rank_sweep.csv");

    double onlineParameter = std::isfinite(options.parameterValue)
        ? options.parameterValue : affine.parameter.reference;
    std::vector<double> onlinePowers = options.onlinePowersW.empty()
        ? nominal : options.onlinePowersW;
    ParametricOnlineResult nominalResult = solveParametricRom(
        model, onlineParameter, onlinePowers, options.allowExtrapolation);

    std::ofstream onlineTiming(outputDirectory / "parametric_online_timing.csv");
    writeOnlineTimingHeader(onlineTiming);
    ParametricRomWorkspace timingWorkspace;
    timingWorkspace.initialize(model);
    std::mt19937_64 generator(options.seed + 991ULL);
    std::uniform_real_distribution<double> parameterDistribution(
        affine.parameter.minimum, affine.parameter.maximum);
    std::uniform_real_distribution<double> powerDistribution(0.0, 1.5);
    std::vector<double> powers(static_cast<std::size_t>(model.sourceChannels), 0.0);
    double averageOnline = 0.0;
    for (int sample = 0; sample < 100; ++sample) {
        for (int channel = 0; channel < model.sourceChannels; ++channel) {
            powers[static_cast<std::size_t>(channel)] = powerDistribution(generator)
                * model.nominalPowersW[static_cast<std::size_t>(channel)];
        }
        const ParametricOnlineResult timing = solveParametricRomInWorkspace(
            model, parameterDistribution(generator), powers, false,
            timingWorkspace);
        averageOnline += timing.timing.totalSeconds;
        writeOnlineTimingRow(onlineTiming, static_cast<std::size_t>(sample), timing);
    }
    averageOnline /= 100.0;
    onlineTiming.close();
    const auto nominalOutputStart = Clock::now();
    writeTemperature(model, nominalResult.temperature,
                     outputDirectory / "parametric_temperature_nodes.csv");
    nominalResult.timing.outputSeconds = elapsed(nominalOutputStart);

    int zeroIterations = 0;
    int romIterations = 0;
    double correctionSeconds = 0.0;
    double correctedResidual = std::numeric_limits<double>::quiet_NaN();
    if (options.parametricMode == "corrected") {
        const auto correctionStart = Clock::now();
        DirectParametricSystem direct = assembleDirectParametricSystem(
            mesh, physics, affine.parameter, onlineParameter);
        CaseConfig currentPhysics = physicsAtParameter(
            physics, affine.parameter, onlineParameter);
        ddm_schur::DdmSchurSolver exact(
            mesh, direct.matrix, currentPhysics, schurOptions);
        const std::vector<double> rhs = composeRhs(
            affine, onlineParameter, onlinePowers);
        const ddm_schur::SolveResult zero = exact.solve(rhs);
        const ddm_schur::SolveResult corrected = exact.solveCorrection(
            rhs, nominalResult.interfaceTemperature);
        zeroIterations = zero.report.iterations;
        romIterations = corrected.report.iterations;
        correctedResidual = corrected.report.interfaceRelativeResidual;
        correctionSeconds = elapsed(correctionStart);
        nominalResult.temperature = corrected.temperature;
        nominalResult.interfaceTemperature = corrected.interfaceSolution;
        nominalResult.status = corrected.report.status;
    }

    const ResidualRatios nominalResidual = calculateAffineResidual(
        affine, onlineParameter, onlinePowers, nominalResult.temperature,
        model.interfaceGlobalDofs);
    double worstRelative = 0.0, worstAbsolute = 0.0, worstMaximum = 0.0;
    int worstParameter = -1, worstChannel = -1;
    for (const AccuracyRow& row : accuracy) {
        worstRelative = std::max(worstRelative, row.relativeL2);
        worstMaximum = std::max(worstMaximum, row.maximumTemperatureError);
        if (row.maximumAbsolute >= worstAbsolute) {
            worstAbsolute = row.maximumAbsolute;
            worstParameter = row.parameterIndex;
            worstChannel = row.powerChannel;
        }
    }
    {
        std::ofstream summary(outputDirectory / "parametric_stage2b1_summary.csv");
        summary << "parameter,units,region_id,subdomain,training_points,interface_rank,"
            << "maximum_local_rank,source_channels,worst_relative_l2,worst_maximum_absolute_error,"
            << "worst_maximum_temperature_error,worst_parameter_index,worst_power_channel,"
            << "average_online_seconds,offline_seconds,model_bytes,nominal_global_residual,"
            << "nominal_interface_equation_residual,status\n"
            << std::setprecision(17);
        int localRank = 0;
        for (const auto& block : model.locals) localRank = std::max(localRank, block.rank);
        summary << model.parameter.name << ',' << model.parameter.units << ','
            << model.parameter.regionId << ',' << model.parameter.subdomain << ','
            << training.size() << ',' << model.interfaceRank << ',' << localRank << ','
            << model.sourceChannels << ',' << worstRelative << ',' << worstAbsolute << ','
            << worstMaximum << ',' << worstParameter << ',' << worstChannel << ','
            << averageOnline << ',' << offline.total << ',' << model.fileBytes << ','
            << nominalResidual.global << ',' << nominalResidual.interfaceEquation << ','
            << nominalResult.status << '\n';
    }
    {
        std::ofstream corrected(outputDirectory / "parametric_corrected_comparison.csv");
        corrected << "mode,zero_initial_guess_iterations,rom_initial_guess_iterations,"
            << "correction_seconds,true_interface_residual\n"
            << options.parametricMode << ',' << zeroIterations << ',' << romIterations << ','
            << correctionSeconds << ',' << correctedResidual << '\n';
    }
    {
        std::ofstream breakEven(outputDirectory / "parametric_break_even.csv");
        breakEven << "parameter_count,rhs_per_parameter,fom_total_seconds,rom_total_seconds,speedup\n";
        double averageFom = 0.0;
        for (const AccuracyRow& row : accuracy) averageFom += row.fomSeconds;
        averageFom /= std::max<std::size_t>(1, accuracy.size());
        for (int parameters : {1, 5, 10, 50}) {
            for (int rhs : {1, 10, 50, 100}) {
                const double fom = parameters * rhs * averageFom;
                const double rom = offline.total + parameters * rhs * averageOnline;
                breakEven << parameters << ',' << rhs << ',' << fom << ',' << rom << ','
                    << fom / std::max(1.0e-300, rom) << '\n';
            }
        }
    }
    {
        std::ofstream report(outputDirectory / "parametric_stage2b1_report.md");
        report << "# Stage 2B.1 parametric local RB / Reduced Schur\n\n"
            << "Implemented one active scalar matrix parameter (`" << model.parameter.name
            << "`, " << model.parameter.units << ") selected by subdomain/region ID. "
            << "The exact separated operator uses `A0+mu*A1` and, for material parameters "
            << "touching SIPG interfaces, analytic harmonic-penalty coefficient groups. "
            << "All components are assembled "
            << "from the existing FEM boundary/volume quadrature, with no finite differences.\n\n"
            << "The pure online path loads only dense reduced blocks and bases. It does not "
            << "construct a sparse FEM matrix, exact Schur operator, proxy, FGMRES, or PARDISO.\n\n"
            << "- training parameter points: " << training.size() << "\n"
            << "- interface rank: " << model.interfaceRank << "\n"
            << "- source channels: " << model.sourceChannels << "\n"
            << "- worst full-field relative L2: " << worstRelative << "\n"
            << "- worst maximum absolute error: " << worstAbsolute << " K\n"
            << "- worst maximum-temperature error: " << worstMaximum << " K\n"
            << "- 100-case average pure online time: " << averageOnline << " s\n"
            << "- serialized model: " << model.fileBytes << " bytes\n"
            << "- offline peak working set: " << offline.peakMemory << " bytes\n"
            << "- nominal global/interface-equation residual: "
            << nominalResidual.global << "/" << nominalResidual.interfaceEquation << "\n"
            << "- corrected zero/ROM initial iterations: " << zeroIterations << "/"
            << romIterations << "\n";
    }

    WorkflowResult workflow;
    workflow.nominalTemperature = nominalResult.temperature;
    workflow.selectedRank = model.interfaceRank;
    workflow.snapshotCount = static_cast<int>(
        training.size() * static_cast<std::size_t>(model.sourceChannels + 1));
    workflow.sourceChannels = model.sourceChannels;
    workflow.correctionIterations = romIterations;
    workflow.setupSeconds = offline.total;
    workflow.nominalOnlineSeconds = nominalResult.timing.totalSeconds;
    workflow.nominalGlobalResidual = nominalResidual.global;
    workflow.nominalInterfaceResidual = std::isfinite(correctedResidual)
        ? correctedResidual : nominalResidual.interfaceEquation;
    workflow.status = nominalResult.status;
    return workflow;
}

void runParametricDeploymentOnly(
    const std::filesystem::path& modelDirectory,
    const std::filesystem::path& outputDirectory,
    const mor::Options& options)
{
    std::filesystem::create_directories(outputDirectory);
    std::ofstream memoryTimeline(
        outputDirectory / "stage2b1_1_memory_timeline.csv");
    writeMemoryTimelineHeader(memoryTimeline);
    writeMemoryTimelineRow(memoryTimeline, "startup", -2, nullptr, nullptr);

    const auto loadStart = Clock::now();
    ParametricReducedModel model = loadParametricModel(modelDirectory);
    truncateParametricModel(model, options.interfaceRank, options.localRank);
    if (!options.saveModelPath.empty()) {
        saveParametricModel(model, options.saveModelPath);
        model.fileBytes = std::filesystem::file_size(
            options.saveModelPath / "model.bin");
    }
    const double loadSeconds = elapsed(loadStart);
    const ParametricModelMemoryBreakdown modelMemory =
        parametricModelMemoryBreakdown(model);
    writeMemoryTimelineRow(
        memoryTimeline, "model_loaded", -1, &modelMemory, nullptr);
    if (options.deploymentRhsCount <= 0) {
        throw std::runtime_error("--mor-deployment-rhs-count must be positive.");
    }

    ParametricRomWorkspace workspace;
    workspace.initialize(model);
    writeMemoryTimelineRow(
        memoryTimeline, "before_case_1", 0, &modelMemory, &workspace);

    std::mt19937_64 generator(options.seed);
    std::uniform_real_distribution<double> parameterDistribution(
        model.parameter.minimum, model.parameter.maximum);
    std::vector<double> powers(static_cast<std::size_t>(model.sourceChannels), 0.0);
    std::ofstream onlineTiming(outputDirectory / "parametric_online_timing.csv");
    onlineTiming.rdbuf()->pubsetbuf(
        workspace.outputScratch.data(),
        static_cast<std::streamsize>(workspace.outputScratch.size()));
    writeOnlineTimingHeader(onlineTiming);
    double representativeParameter = 0.0;
    bool representativeExtrapolated = false;
    std::string representativeStatus = "not_run";
    double accumulatedComputeSeconds = 0.0;
    double outputSeconds = 0.0;
    const auto batchStart = Clock::now();
    for (int index = 0; index < options.deploymentRhsCount; ++index) {
        const double value = std::isfinite(options.parameterValue)
            ? options.parameterValue
            : (options.deploymentRhsCount == 1
                ? model.parameter.reference : parameterDistribution(generator));
        if (!options.onlinePowersW.empty()) {
            std::copy(options.onlinePowersW.begin(), options.onlinePowersW.end(),
                      powers.begin());
        } else {
            std::copy(model.nominalPowersW.begin(), model.nominalPowersW.end(),
                      powers.begin());
            if (options.deploymentRhsCount > 1) {
                for (int channel = 0; channel < model.sourceChannels; ++channel) {
                    std::uniform_real_distribution<double> powerDistribution(
                        model.minimumPowersW[static_cast<std::size_t>(channel)],
                        model.maximumPowersW[static_cast<std::size_t>(channel)]);
                    powers[static_cast<std::size_t>(channel)] =
                        powerDistribution(generator);
                }
            }
        }
        ParametricOnlineResult result = solveParametricRomInWorkspace(
            model, value, powers, options.allowExtrapolation, workspace);
        accumulatedComputeSeconds += result.timing.totalSeconds;
        if (index == 0) {
            representativeParameter = result.parameterValue;
            representativeExtrapolated = result.extrapolated;
            representativeStatus = result.status;
            const auto outputStart = Clock::now();
            writeTemperature(model, workspace.fullTemperature,
                             outputDirectory / "parametric_temperature_nodes.csv");
            result.timing.outputSeconds = elapsed(outputStart);
            outputSeconds += result.timing.outputSeconds;
        }
        writeOnlineTimingRow(
            onlineTiming, static_cast<std::size_t>(index), result);
        const int completed = index + 1;
        if (completed == 1 || completed == 10
            || completed == 50 || completed == 100) {
            writeMemoryTimelineRow(
                memoryTimeline, "after_case_" + std::to_string(completed),
                completed, &modelMemory, &workspace);
        }
    }
    const double batchWallSeconds = elapsed(batchStart);
    onlineTiming.flush();
    onlineTiming.close();
    const double averageComputeSeconds = accumulatedComputeSeconds
        / static_cast<double>(options.deploymentRhsCount);
    const std::size_t deploymentPeakBytes = peakWorkingSetBytes();
    std::ofstream summary(outputDirectory / "parametric_deployment_summary.csv");
    summary << "deployment_cases,representative_parameter_value,source_channels,interface_rank,"
        << "load_seconds,average_compute_seconds,batch_compute_seconds,batch_wall_seconds,"
        << "output_seconds,end_to_end_seconds,model_bytes,peak_working_set_bytes,"
        << "extrapolated,status\n" << std::setprecision(17)
        << options.deploymentRhsCount << ',' << representativeParameter << ','
        << model.sourceChannels << ',' << model.interfaceRank << ',' << loadSeconds << ','
        << averageComputeSeconds << ',' << accumulatedComputeSeconds << ','
        << batchWallSeconds << ',' << outputSeconds << ','
        << loadSeconds + batchWallSeconds << ','
        << model.fileBytes << ',' << deploymentPeakBytes << ','
        << (representativeExtrapolated ? 1 : 0) << ',' << representativeStatus << '\n';
    summary.flush();

    const ProcessMemorySnapshot finalMemory = processMemorySnapshot();
    std::ofstream batchSummary(
        outputDirectory / "stage2b1_1_batch_summary.csv");
    batchSummary << "deployment_cases,baseline_peak_working_set_mib,"
        << "optimized_peak_working_set_mib,target_peak_working_set_mib,"
        << "preferred_peak_working_set_mib,baseline_average_online_seconds,"
        << "optimized_average_online_seconds,target_average_online_seconds,"
        << "load_seconds,batch_wall_seconds,reconstruction_output_seconds,"
        << "model_file_mib,model_resident_mib,interface_basis_mib,local_bases_mib,"
        << "reduced_blocks_mib,sources_mib,metadata_mib,workspace_mib,"
        << "reconstruction_mib,io_buffer_mib,final_private_mib,"
        << "result_cache_mib,mathematics_changed,status\n" << std::setprecision(17)
        << options.deploymentRhsCount << ",1206.272," << mib(deploymentPeakBytes)
        << ",1024,800,0.08336," << averageComputeSeconds << ",0.09,"
        << loadSeconds << ',' << batchWallSeconds << ',' << outputSeconds << ','
        << mib(model.fileBytes) << ',' << mib(modelMemory.totalBytes()) << ','
        << mib(modelMemory.interfaceBasisBytes) << ','
        << mib(modelMemory.localBasisBytes) << ','
        << mib(modelMemory.reducedBlockBytes) << ','
        << mib(modelMemory.sourceBytes) << ',' << mib(modelMemory.metadataBytes) << ','
        << mib(workspace.workspaceBytes()) << ','
        << mib(workspace.reconstructionBytes()) << ','
        << mib(workspace.outputScratch.capacity()) << ','
        << mib(finalMemory.privateBytes) << ",0,0,"
        << ((mib(deploymentPeakBytes) < 1024.0
                && averageComputeSeconds <= 0.09) ? "pass" : "review") << '\n';
    batchSummary.flush();
    writeMemoryTimelineRow(
        memoryTimeline, "after_output", options.deploymentRhsCount + 1,
        &modelMemory, &workspace);
    memoryTimeline.close();

    std::cout << "[Stage 2B.1.1 deployment] load=" << loadSeconds
              << " s, cases=" << options.deploymentRhsCount
              << ", average pure online=" << averageComputeSeconds
              << " s, model=" << model.fileBytes
              << " bytes, peak working set=" << deploymentPeakBytes << " bytes\n";
}

} // namespace mor::parametric
