#include "deployment_response_model.hpp"

#include "reduced_schur_model.hpp"
#include "source_parameterization.hpp"
#include "interior_model_serialization.hpp"
#include "local_interior_basis.hpp"
#include "sipg_core.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <fstream>
#include <random>
#include <stdexcept>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

namespace mor {
namespace {

class Fnv1a {
public:
    template <typename T>
    void add(const T& value)
    {
        const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
        for (std::size_t i = 0; i < sizeof(T); ++i) {
            value_ ^= static_cast<std::uint64_t>(bytes[i]);
            value_ *= 1099511628211ULL;
        }
    }
    std::uint64_t value() const { return value_; }
private:
    std::uint64_t value_ = 1469598103934665603ULL;
};

void matrixVector(const std::vector<double>& matrix,
                  int rows,
                  int columns,
                  const std::vector<double>& x,
                  std::vector<double>& y)
{
    if (static_cast<int>(matrix.size()) < rows * columns
        || static_cast<int>(x.size()) != columns) {
        throw std::runtime_error("[ROM Deployment] Dense response dimensions are inconsistent.");
    }
    y.assign(static_cast<std::size_t>(rows), 0.0);
    if (rows == 0 || columns == 0) {
        return;
    }
#ifdef USE_MKL_PARDISO
    cblas_dgemv(CblasColMajor, CblasNoTrans, rows, columns, 1.0,
        matrix.data(), rows, x.data(), 1, 0.0, y.data(), 1);
#else
    for (int column = 0; column < columns; ++column) {
        const double coefficient = x[static_cast<std::size_t>(column)];
        if (coefficient == 0.0) {
            continue;
        }
        const std::size_t offset = static_cast<std::size_t>(column * rows);
        for (int row = 0; row < rows; ++row) {
            y[static_cast<std::size_t>(row)] += coefficient
                * matrix[offset + static_cast<std::size_t>(row)];
        }
    }
#endif
}

std::uint64_t hashDofs(const std::vector<int>& dofs)
{
    Fnv1a hash;
    const std::uint64_t count = static_cast<std::uint64_t>(dofs.size());
    hash.add(count);
    for (int dof : dofs) {
        hash.add(dof);
    }
    return hash.value();
}

} // namespace

DeploymentResponseModel buildDeploymentResponseModel(
    const Mesh& mesh,
    const SourceParameterization& sources,
    const std::vector<int>& interfaceGlobalDofs,
    const ReducedSchurModel& interfaceModel,
    ReducedSchurOnlineSolver& onlineSolver)
{
    const auto constructionStart = std::chrono::steady_clock::now();
    DeploymentResponseModel model;
    model.globalDofs = static_cast<int>(mesh.nodes.size());
    model.interfaceDofs = static_cast<int>(interfaceGlobalDofs.size());
    model.interfaceRank = interfaceModel.rank;
    model.sourceChannels = static_cast<int>(sources.channels.size());
    model.fingerprints = interfaceModel.fingerprints;
    model.interfaceGlobalDofs = interfaceGlobalDofs;
    model.dofs.resize(mesh.nodes.size());
    model.nominalPowers.resize(sources.channels.size());
    model.minimumPowers.resize(sources.channels.size());
    model.maximumPowers.resize(sources.channels.size());
    model.sourceSubdomains.resize(sources.channels.size());
    model.sourceDomainEntities.resize(sources.channels.size());
    for (std::size_t channel = 0; channel < sources.channels.size(); ++channel) {
        model.nominalPowers[channel] = sources.channels[channel].nominalPowerW;
        model.minimumPowers[channel] = sources.channels[channel].minimumPowerW;
        model.maximumPowers[channel] = sources.channels[channel].maximumPowerW;
        model.sourceSubdomains[channel] = sources.channels[channel].subdomain;
        model.sourceDomainEntities[channel] = sources.channels[channel].domainEntity;
    }

    Fnv1a globalOrdering;
    const std::uint64_t globalCount = static_cast<std::uint64_t>(mesh.nodes.size());
    globalOrdering.add(globalCount);
    int subdomains = 0;
    for (std::size_t dof = 0; dof < mesh.nodes.size(); ++dof) {
        const Node& node = mesh.nodes[dof];
        DeploymentDof& descriptor = model.dofs[dof];
        descriptor.x = node.p.x;
        descriptor.y = node.p.y;
        descriptor.z = node.p.z;
        descriptor.subdomain = node.subdomain;
        descriptor.sourceVertex = node.sourceVertex;
        globalOrdering.add(dof);
        globalOrdering.add(node.subdomain);
        globalOrdering.add(node.sourceVertex);
        globalOrdering.add(node.p.x);
        globalOrdering.add(node.p.y);
        globalOrdering.add(node.p.z);
        subdomains = std::max(subdomains, node.subdomain + 1);
    }
    model.globalDofOrderingHash = globalOrdering.value();

    std::vector<char> isInterface(mesh.nodes.size(), 0);
    for (int dof : interfaceGlobalDofs) {
        if (dof < 0 || dof >= model.globalDofs
            || isInterface[static_cast<std::size_t>(dof)]) {
            throw std::runtime_error("[MOR Interior] Invalid or duplicate interface DOF ordering.");
        }
        isInterface[static_cast<std::size_t>(dof)] = 1;
    }
    model.interiors.resize(static_cast<std::size_t>(subdomains));
    for (int domain = 0; domain < subdomains; ++domain) {
        model.interiors[static_cast<std::size_t>(domain)].subdomain = domain;
    }
    for (int dof = 0; dof < model.globalDofs; ++dof) {
        if (!isInterface[static_cast<std::size_t>(dof)]) {
            const int domain = mesh.nodes[static_cast<std::size_t>(dof)].subdomain;
            model.interiors[static_cast<std::size_t>(domain)].globalDofs.push_back(dof);
        }
    }
    for (InteriorResponseBlock& block : model.interiors) {
        block.dofOrderingHash = hashDofs(block.globalDofs);
        for (int channel = 0; channel < model.sourceChannels; ++channel) {
            bool hasDirectInteriorLoad = false;
            const std::vector<double>& load =
                sources.channels[static_cast<std::size_t>(channel)].rhsPerWatt;
            for (int dof : block.globalDofs) {
                if (load[static_cast<std::size_t>(dof)] != 0.0) {
                    hasDirectInteriorLoad = true;
                    break;
                }
            }
            if (hasDirectInteriorLoad) {
                block.directPowerChannels.push_back(channel);
            }
        }
    }

    const SolveResult reference = onlineSolver.solve(sources.referenceRhs, false);
    if (reference.status != "success"
        || static_cast<int>(reference.temperature.size()) != model.globalDofs
        || static_cast<int>(reference.interfaceTemperature.size()) != model.interfaceDofs) {
        throw std::runtime_error("[MOR Interior] Reference Reduced Schur recovery failed.");
    }
    model.interfaceReference = reference.interfaceTemperature;
    model.interfacePowerResponse.assign(
        static_cast<std::size_t>(model.interfaceDofs * model.sourceChannels), 0.0);
    model.reducedInputOperator.assign(
        static_cast<std::size_t>(model.interfaceRank * model.sourceChannels), 0.0);
    for (InteriorResponseBlock& block : model.interiors) {
        const int rows = static_cast<int>(block.globalDofs.size());
        block.referenceTemperature.resize(static_cast<std::size_t>(rows));
        block.exactResponse.assign(
            static_cast<std::size_t>(rows * model.sourceChannels), 0.0);
        for (int row = 0; row < rows; ++row) {
            block.referenceTemperature[static_cast<std::size_t>(row)] =
                reference.temperature[static_cast<std::size_t>(
                    block.globalDofs[static_cast<std::size_t>(row)])];
        }
    }

    std::cout << "[MOR Interior] constructing exact response columns: "
              << model.sourceChannels << '\n';
    for (int channel = 0; channel < model.sourceChannels; ++channel) {
        std::vector<double> powers(static_cast<std::size_t>(model.sourceChannels), 0.0);
        powers[static_cast<std::size_t>(channel)] = 1.0;
        const SolveResult unit = onlineSolver.solve(composeRhs(sources, powers), false);
        if (unit.status != "success") {
            throw std::runtime_error("[MOR Interior] Unit response recovery failed for channel "
                + std::to_string(channel));
        }
        const std::size_t interfaceOffset = static_cast<std::size_t>(
            channel * model.interfaceDofs);
        for (int row = 0; row < model.interfaceDofs; ++row) {
            model.interfacePowerResponse[interfaceOffset + static_cast<std::size_t>(row)] =
                unit.interfaceTemperature[static_cast<std::size_t>(row)]
                - reference.interfaceTemperature[static_cast<std::size_t>(row)];
        }
        for (int mode = 0; mode < model.interfaceRank; ++mode) {
            double coefficient = 0.0;
#ifdef USE_MKL_PARDISO
            coefficient = cblas_ddot(model.interfaceDofs,
                interfaceModel.basis.data()
                    + static_cast<std::size_t>(mode * model.interfaceDofs), 1,
                model.interfacePowerResponse.data() + interfaceOffset, 1);
#else
            for (int row = 0; row < model.interfaceDofs; ++row) {
                coefficient += interfaceModel.basis[static_cast<std::size_t>(
                    mode * model.interfaceDofs + row)]
                    * model.interfacePowerResponse[interfaceOffset
                        + static_cast<std::size_t>(row)];
            }
#endif
            model.reducedInputOperator[static_cast<std::size_t>(
                channel * model.interfaceRank + mode)] = coefficient;
        }
        for (InteriorResponseBlock& block : model.interiors) {
            const int rows = static_cast<int>(block.globalDofs.size());
            const std::size_t offset = static_cast<std::size_t>(channel * rows);
            for (int row = 0; row < rows; ++row) {
                const int dof = block.globalDofs[static_cast<std::size_t>(row)];
                block.exactResponse[offset + static_cast<std::size_t>(row)] =
                    unit.temperature[static_cast<std::size_t>(dof)]
                    - reference.temperature[static_cast<std::size_t>(dof)];
            }
        }
        if ((channel + 1) % 10 == 0 || channel + 1 == model.sourceChannels) {
            std::cout << "[MOR Interior] exact response " << (channel + 1)
                      << '/' << model.sourceChannels << '\n';
        }
    }
    model.responseConstructionSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - constructionStart).count();
    return model;
}

SolveResult evaluateDeploymentResponse(
    const DeploymentResponseModel& model,
    const std::vector<double>& powersW,
    const std::string& interiorMode,
    int rankOverride)
{
    if (static_cast<int>(powersW.size()) != model.sourceChannels) {
        throw std::runtime_error("[ROM Deployment] Power vector size mismatch.");
    }
    if (interiorMode != "exact-response" && interiorMode != "compressed-rb") {
        throw std::runtime_error("[ROM Deployment] Interior mode must be exact-response or compressed-rb.");
    }
    SolveResult result;
    const auto totalStart = std::chrono::steady_clock::now();
    const auto interfaceStart = std::chrono::steady_clock::now();
    std::vector<double> interfaceRise;
    matrixVector(model.interfacePowerResponse, model.interfaceDofs,
                 model.sourceChannels, powersW, interfaceRise);
    result.interfaceTemperature = model.interfaceReference;
    for (int row = 0; row < model.interfaceDofs; ++row) {
        result.interfaceTemperature[static_cast<std::size_t>(row)] +=
            interfaceRise[static_cast<std::size_t>(row)];
    }
    result.timing.interfaceReconstructionSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - interfaceStart).count();

    const auto interiorStart = std::chrono::steady_clock::now();
    result.temperature.assign(static_cast<std::size_t>(model.globalDofs), 0.0);
    for (int row = 0; row < model.interfaceDofs; ++row) {
        result.temperature[static_cast<std::size_t>(
            model.interfaceGlobalDofs[static_cast<std::size_t>(row)])] =
            result.interfaceTemperature[static_cast<std::size_t>(row)];
    }
    for (const InteriorResponseBlock& block : model.interiors) {
        const int rows = static_cast<int>(block.globalDofs.size());
        std::vector<double> rise;
        if (interiorMode == "exact-response") {
            matrixVector(block.exactResponse, rows, model.sourceChannels, powersW, rise);
        } else {
            int rank = block.rank;
            if (rankOverride > 0) {
                rank = std::min(rank, rankOverride);
            }
            if (rank <= 0
                || static_cast<int>(block.localBasis.size()) < rows * rank
                || static_cast<int>(block.localCoordinateMap.size())
                    < block.rank * model.sourceChannels) {
                throw std::runtime_error("[ROM Deployment] Compressed local response is missing.");
            }
            std::vector<double> coordinates(static_cast<std::size_t>(rank), 0.0);
            for (int channel = 0; channel < model.sourceChannels; ++channel) {
                const double power = powersW[static_cast<std::size_t>(channel)];
                const std::size_t offset = static_cast<std::size_t>(channel * block.rank);
                for (int mode = 0; mode < rank; ++mode) {
                    coordinates[static_cast<std::size_t>(mode)] += power
                        * block.localCoordinateMap[offset + static_cast<std::size_t>(mode)];
                }
            }
            // POD/SVD modes are stored contiguously in descending importance,
            // so a rank audit can use the leading columns in place.
            matrixVector(block.localBasis, rows, rank, coordinates, rise);
        }
        for (int row = 0; row < rows; ++row) {
            const int dof = block.globalDofs[static_cast<std::size_t>(row)];
            result.temperature[static_cast<std::size_t>(dof)] =
                block.referenceTemperature[static_cast<std::size_t>(row)]
                + rise[static_cast<std::size_t>(row)];
        }
    }
    result.timing.recoverySeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - interiorStart).count();
    result.timing.totalSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - totalStart).count();
    result.status = "success";
    return result;
}

DeploymentRunResult runDeploymentOnly(
    const std::filesystem::path& modelDirectory,
    const std::filesystem::path& outputDirectory,
    const Options& options)
{
    const int rhsCount = options.deploymentRhsCount;
    if (rhsCount <= 0) {
        throw std::runtime_error("--mor-deployment-rhs-count must be positive.");
    }
    const auto endToEndStart = std::chrono::steady_clock::now();
    std::filesystem::create_directories(outputDirectory);
    DeploymentResponseModel model = loadDeploymentResponseModel(
        modelDirectory, "");
    std::string mode = options.interiorMode == "pardiso"
        ? model.interiorMode : options.interiorMode;
    const bool exactRepresentationAvailable = model.interiorMode == "exact-response";
    if (mode == "compressed-rb" && model.interiorMode == "exact-response") {
        int buildRank = options.interiorRank;
        if (buildRank > 0 && !options.interiorRankSweep.empty()) {
            buildRank = std::max(buildRank, *std::max_element(
                options.interiorRankSweep.begin(), options.interiorRankSweep.end()));
        }
        buildCompressedInteriorBases(model, buildRank,
            options.interiorEnergyTolerance,
            options.interiorSingularValueTolerance);
        if (!options.saveModelPath.empty()) {
            saveDeploymentResponseModel(model,
                std::filesystem::absolute(options.saveModelPath).lexically_normal(),
                "compressed-rb", options.storagePrecision);
        }
    } else if (mode == "exact-response" && model.interiorMode == "exact-response"
               && !options.saveModelPath.empty()) {
        saveDeploymentResponseModel(model,
            std::filesystem::absolute(options.saveModelPath).lexically_normal(),
            "exact-response", options.storagePrecision);
    } else if (mode != model.interiorMode) {
        throw std::runtime_error(
            "[ROM Deployment] Requested representation is not stored and cannot be derived.");
    }

    std::mt19937_64 generator(options.seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::ofstream timingOut(outputDirectory / "local_interior_online_timing.csv");
    timingOut << "case,interface_seconds,interior_seconds,compute_seconds\n"
              << std::setprecision(17);
    DeploymentRunResult result;
    result.rhsCount = rhsCount;
    result.setupSeconds = model.loadSeconds;
    result.deploymentFileBytes = model.deploymentFileBytes;
    std::vector<std::vector<double>> benchmarkPowers;
    benchmarkPowers.reserve(static_cast<std::size_t>(rhsCount));
    for (int sample = 0; sample < rhsCount; ++sample) {
        std::vector<double> powers = model.nominalPowers;
        if (sample > 0) {
            for (int channel = 0; channel < model.sourceChannels; ++channel) {
                const double lo = model.minimumPowers[static_cast<std::size_t>(channel)];
                const double hi = model.maximumPowers[static_cast<std::size_t>(channel)];
                if (sample % 7 == 1) {
                    powers[static_cast<std::size_t>(channel)] =
                        channel == sample % model.sourceChannels
                        ? model.nominalPowers[static_cast<std::size_t>(channel)] : 0.0;
                } else if (sample % 7 == 2) {
                    powers[static_cast<std::size_t>(channel)] =
                        (channel == sample % model.sourceChannels
                         || channel == (sample * 5 + 1) % model.sourceChannels)
                        ? model.nominalPowers[static_cast<std::size_t>(channel)] : 0.0;
                } else if (sample % 7 == 3) {
                    powers[static_cast<std::size_t>(channel)] =
                        ((channel + sample) % 2 == 0 ? 1.4 : 0.2)
                        * model.nominalPowers[static_cast<std::size_t>(channel)];
                } else {
                    powers[static_cast<std::size_t>(channel)] = lo + uniform(generator) * (hi - lo);
                }
            }
        }
        benchmarkPowers.push_back(powers);
        SolveResult evaluated = evaluateDeploymentResponse(model, powers, mode);
        if (sample == 0) {
            result.nominalTemperature = evaluated.temperature;
        }
        result.averageComputeSeconds += evaluated.timing.totalSeconds;
        result.averageInterfaceSeconds += evaluated.timing.interfaceReconstructionSeconds;
        result.averageInteriorSeconds += evaluated.timing.recoverySeconds;
        timingOut << sample << ',' << evaluated.timing.interfaceReconstructionSeconds
                  << ',' << evaluated.timing.recoverySeconds << ','
                  << evaluated.timing.totalSeconds << '\n';
    }
    result.averageComputeSeconds /= static_cast<double>(rhsCount);
    result.averageInterfaceSeconds /= static_cast<double>(rhsCount);
    result.averageInteriorSeconds /= static_cast<double>(rhsCount);

    if (exactRepresentationAvailable && !model.interiors.empty()
        && !model.interiors.front().localBasis.empty()) {
        std::ofstream accuracy(outputDirectory / "local_interior_accuracy_by_case.csv");
        accuracy << "case,rank,relative_l2,maximum_absolute_error,hotspot_maximum_error\n"
                 << std::setprecision(17);
        std::ofstream rankSweep(outputDirectory / "local_interior_rank_sweep.csv");
        rankSweep << "rank,cases,average_compute_seconds,maximum_relative_l2,"
                     "maximum_absolute_error,maximum_hotspot_error,worst_unit_channel_relative_l2\n"
                  << std::setprecision(17);
        std::vector<int> ranks = options.interiorRankSweep;
        if (options.interiorRank > 0) {
            ranks.push_back(options.interiorRank);
        } else {
            // Rank zero means: use the independently energy-selected rank stored
            // in each subdomain block.
            ranks.push_back(0);
        }
        std::sort(ranks.begin(), ranks.end());
        ranks.erase(std::unique(ranks.begin(), ranks.end()), ranks.end());
        for (int rank : ranks) {
            if (rank < 0) {
                continue;
            }
            double maximumRelative = 0.0;
            double maximumAbsolute = 0.0;
            double maximumHotspot = 0.0;
            double worstUnitRelative = 0.0;
            double computeSeconds = 0.0;
            int caseIndex = 0;
            auto compare = [&](const std::vector<double>& powers, bool unitChannel) {
                const SolveResult truth = evaluateDeploymentResponse(
                    model, powers, "exact-response");
                const SolveResult reduced = evaluateDeploymentResponse(
                    model, powers, "compressed-rb", rank);
                long double errorSquared = 0.0L;
                long double truthSquared = 0.0L;
                double maximum = 0.0;
                for (std::size_t dof = 0; dof < truth.temperature.size(); ++dof) {
                    const double error = reduced.temperature[dof] - truth.temperature[dof];
                    errorSquared += static_cast<long double>(error) * error;
                    truthSquared += static_cast<long double>(truth.temperature[dof])
                        * truth.temperature[dof];
                    maximum = std::max(maximum, std::abs(error));
                }
                const double relative = std::sqrt(static_cast<double>(errorSquared
                    / std::max(std::numeric_limits<long double>::min(), truthSquared)));
                const double reducedMaximum = *std::max_element(
                    reduced.temperature.begin(), reduced.temperature.end());
                const double truthMaximum = *std::max_element(
                    truth.temperature.begin(), truth.temperature.end());
                const double hotspot = std::abs(reducedMaximum - truthMaximum);
                maximumRelative = std::max(maximumRelative, relative);
                maximumAbsolute = std::max(maximumAbsolute, maximum);
                maximumHotspot = std::max(maximumHotspot, hotspot);
                if (unitChannel) {
                    worstUnitRelative = std::max(worstUnitRelative, relative);
                }
                computeSeconds += reduced.timing.totalSeconds;
                if (rank == options.interiorRank || options.interiorRank <= 0) {
                    accuracy << caseIndex << ',' << rank << ',' << relative << ','
                             << maximum << ',' << hotspot << '\n';
                }
                ++caseIndex;
            };
            for (const std::vector<double>& powers : benchmarkPowers) {
                compare(powers, false);
            }
            for (int channel = 0; channel < model.sourceChannels; ++channel) {
                std::vector<double> unit(static_cast<std::size_t>(model.sourceChannels), 0.0);
                unit[static_cast<std::size_t>(channel)] = 1.0;
                compare(unit, true);
            }
            rankSweep << rank << ',' << caseIndex << ','
                      << computeSeconds / static_cast<double>(std::max(1, caseIndex)) << ','
                      << maximumRelative << ',' << maximumAbsolute << ','
                      << maximumHotspot << ',' << worstUnitRelative << '\n';
        }
    }

    const auto ioStart = std::chrono::steady_clock::now();
    std::ofstream temperatureOut(outputDirectory / "temperature_rom_deployment_nodes.csv");
    temperatureOut << "global_dof,subdomain,source_vertex,x,y,z,temperature\n"
                   << std::setprecision(17);
    for (int dof = 0; dof < model.globalDofs; ++dof) {
        const DeploymentDof& descriptor = model.dofs[static_cast<std::size_t>(dof)];
        temperatureOut << dof << ',' << descriptor.subdomain << ','
                       << descriptor.sourceVertex << ',' << descriptor.x << ','
                       << descriptor.y << ',' << descriptor.z << ','
                       << result.nominalTemperature[static_cast<std::size_t>(dof)] << '\n';
    }
    temperatureOut.close();
    result.outputIoSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - ioStart).count();
    result.peakWorkingSetBytes = peakWorkingSetBytes();
    result.endToEndSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - endToEndStart).count();

    int minimumRank = std::numeric_limits<int>::max();
    int maximumRank = 0;
    for (const InteriorResponseBlock& block : model.interiors) {
        minimumRank = std::min(minimumRank, block.rank);
        maximumRank = std::max(maximumRank, block.rank);
    }
    if (mode == "exact-response") {
        minimumRank = 0;
        maximumRank = 0;
    }
    std::ofstream summary(outputDirectory / "deployment_summary.json");
    summary << std::setprecision(17)
        << "{\n"
        << "  \"model\": \"" << modelDirectory.generic_string() << "\",\n"
        << "  \"interior_mode\": \"" << mode << "\",\n"
        << "  \"global_dofs\": " << model.globalDofs << ",\n"
        << "  \"interface_rank\": " << model.interfaceRank << ",\n"
        << "  \"minimum_local_rank\": " << minimumRank << ",\n"
        << "  \"maximum_local_rank\": " << maximumRank << ",\n"
        << "  \"rhs_count\": " << rhsCount << ",\n"
        << "  \"model_file_bytes\": " << model.deploymentFileBytes << ",\n"
        << "  \"setup_seconds\": " << result.setupSeconds << ",\n"
        << "  \"compression_setup_seconds\": " << model.compressionSeconds << ",\n"
        << "  \"average_online_compute_seconds\": " << result.averageComputeSeconds << ",\n"
        << "  \"average_interface_seconds\": " << result.averageInterfaceSeconds << ",\n"
        << "  \"average_interior_reconstruction_seconds\": "
        << result.averageInteriorSeconds << ",\n"
        << "  \"output_io_seconds\": " << result.outputIoSeconds << ",\n"
        << "  \"end_to_end_seconds\": " << result.endToEndSeconds << ",\n"
        << "  \"peak_working_set_bytes\": " << result.peakWorkingSetBytes << ",\n"
        << "  \"fem_assembly_initialized\": false,\n"
        << "  \"schur_initialized\": false,\n"
        << "  \"fgmres_initialized\": false,\n"
        << "  \"pardiso_initialized\": false\n"
        << "}\n";

    std::ofstream memory(outputDirectory / "local_interior_memory_breakdown.csv");
    const std::uint64_t interfaceBytes = static_cast<std::uint64_t>(
        (model.interfaceReference.size() + model.interfacePowerResponse.size())
        * sizeof(double));
    std::uint64_t exactInteriorBytes = 0;
    std::uint64_t compressedInteriorBytes = 0;
    for (const InteriorResponseBlock& block : model.interiors) {
        exactInteriorBytes += static_cast<std::uint64_t>(
            block.exactResponse.size() * sizeof(double));
        compressedInteriorBytes += static_cast<std::uint64_t>(
            (block.localBasis.size() + block.localCoordinateMap.size()) * sizeof(double));
    }
    memory << "component,bytes\n"
           << "interface_response," << interfaceBytes << '\n'
           << "exact_interior_response_in_memory," << exactInteriorBytes << '\n'
           << "compressed_interior_response_in_memory," << compressedInteriorBytes << '\n'
           << "deployment_files," << model.deploymentFileBytes << '\n'
           << "peak_working_set," << result.peakWorkingSetBytes << '\n';
    std::ofstream sizes(outputDirectory / "local_interior_model_sizes.csv");
    sizes << "interior_mode,deployment_file_bytes,complete_model_file_bytes,storage_precision\n"
          << mode << ',' << model.deploymentFileBytes << ','
          << model.completeModelFileBytes << ',' << model.storagePrecision << '\n';
    std::ofstream ranksOut(outputDirectory / "local_interior_subdomain_ranks.csv");
    ranksOut << "subdomain,interior_dofs,direct_power_channels,selected_rank,"
                "available_singular_values,retained_energy\n"
             << std::setprecision(17);
    for (const InteriorResponseBlock& block : model.interiors) {
        const double retained = block.rank > 0
            && static_cast<std::size_t>(block.rank) <= block.retainedEnergy.size()
            ? block.retainedEnergy[static_cast<std::size_t>(block.rank - 1)] : 1.0;
        ranksOut << block.subdomain << ',' << block.globalDofs.size() << ','
                 << block.directPowerChannels.size() << ',' << block.rank << ','
                 << block.singularValues.size() << ',' << retained << '\n';
    }

    std::cout << std::setprecision(12)
              << "[ROM Deployment]\n"
              << "model: " << modelDirectory.string() << '\n'
              << "interior mode: " << mode << '\n'
              << "global DOFs: " << model.globalDofs << '\n'
              << "interface rank: " << model.interfaceRank << '\n'
              << "local ranks: " << minimumRank << ".." << maximumRank << '\n'
              << "model file size: " << model.deploymentFileBytes << " bytes\n"
              << "setup time: " << result.setupSeconds << " s\n"
              << "peak memory: " << result.peakWorkingSetBytes << " bytes\n"
              << "online evaluation time: " << result.averageComputeSeconds << " s/RHS\n"
              << "reconstruction time: " << result.averageInteriorSeconds << " s/RHS\n"
              << "output I/O time: " << result.outputIoSeconds << " s\n"
              << "FEM/Schur/FGMRES/PARDISO initialized: no/no/no/no\n";
    if (!options.reportIoTime) {
        std::cout << "[ROM Deployment] I/O was measured separately; use --mor-report-io-time "
                     "to request it explicitly in benchmark commands.\n";
    }
    return result;
}

} // namespace mor
