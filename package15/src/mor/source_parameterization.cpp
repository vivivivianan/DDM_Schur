#include "source_parameterization.hpp"

#include "sipg_core.hpp"
#include "config_io.hpp"
#include "mesh_loader.hpp"
#include "fem_assembly.hpp"
#include "linear_solvers.hpp"
#include "diagnostics_io.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>
#include <stdexcept>

namespace mor {
namespace {

std::vector<std::vector<double>> assembleUnitHeatChannels(
    const Mesh& mesh,
    const CaseConfig& physics)
{
    const std::size_t channelCount = physics.heatSources.size();
    std::vector<double> selectedVolumes(channelCount, 0.0);
    for (const Tet& tet : mesh.tets) {
        const double volume = tetVolume(mesh, tet);
        for (std::size_t channel = 0; channel < channelCount; ++channel) {
            if (tetMatchesHeatSource(tet, physics.heatSources[channel])) {
                selectedVolumes[channel] += volume;
            }
        }
    }
    std::vector<std::vector<double>> channels(channelCount,
        std::vector<double>(mesh.nodes.size(), 0.0));
    const std::vector<TetQuadraturePoint> quadrature = makeTetQuadrature();
    for (const Tet& tet : mesh.tets) {
        std::vector<std::size_t> matching;
        for (std::size_t channel = 0; channel < channelCount; ++channel) {
            if (tetMatchesHeatSource(tet, physics.heatSources[channel])) {
                matching.push_back(channel);
            }
        }
        if (matching.empty()) {
            continue;
        }
        const ElementGeometry geometry = elementGeometry(mesh, tet);
        std::array<double, 10> unitDensityLoad{};
        for (const TetQuadraturePoint& point : quadrature) {
            const double weight = point.weight * geometry.detJ;
            for (int local = 0; local < 10; ++local) {
                unitDensityLoad[static_cast<std::size_t>(local)] += weight
                    * point.shape[static_cast<std::size_t>(local)];
            }
        }
        for (std::size_t channel : matching) {
            if (!(selectedVolumes[channel] > 0.0)) {
                continue;
            }
            const double inverseVolume = 1.0 / selectedVolumes[channel];
            for (int local = 0; local < 10; ++local) {
                channels[channel][static_cast<std::size_t>(
                    tet.dof[static_cast<std::size_t>(local)])] += inverseVolume
                    * unitDensityLoad[static_cast<std::size_t>(local)];
            }
        }
    }
    return channels;
}

} // namespace

SourceParameterization buildSourceParameterization(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::vector<double>& assembledSource,
    const std::vector<double>& heatOnlySource,
    const std::vector<double>& fixedAdjust)
{
    if (assembledSource.size() != heatOnlySource.size()
        || assembledSource.size() != mesh.nodes.size()) {
        throw std::runtime_error("[MOR] Source vector dimensions are inconsistent.");
    }
    SourceParameterization result;
    std::vector<double> referenceBeforeDirichlet(assembledSource.size(), 0.0);
    for (std::size_t i = 0; i < assembledSource.size(); ++i) {
        referenceBeforeDirichlet[i] = assembledSource[i] - heatOnlySource[i];
    }
    result.referenceRhs = referenceBeforeDirichlet;
    applyDirichletRhs(mesh, fixedAdjust, result.referenceRhs);

    std::vector<std::vector<double>> unitHeatChannels =
        assembleUnitHeatChannels(mesh, physics);
    result.channels.reserve(physics.heatSources.size());
    for (std::size_t channelIndex = 0;
         channelIndex < physics.heatSources.size(); ++channelIndex) {
        const HeatSource& heatSource = physics.heatSources[channelIndex];
        if (!(heatSource.heatRateW > 0.0)) {
            throw std::runtime_error("[MOR] Every parameterized heat-source channel must have positive nominal power.");
        }
        std::vector<double> unitCase = std::move(unitHeatChannels[channelIndex]);
        for (std::size_t i = 0; i < unitCase.size(); ++i) {
            unitCase[i] += referenceBeforeDirichlet[i];
        }
        applyDirichletRhs(mesh, fixedAdjust, unitCase);

        SourceChannel channel;
        channel.index = static_cast<int>(channelIndex);
        channel.subdomain = heatSource.subdomain;
        channel.domainEntity = heatSource.domainEntity;
        channel.nominalPowerW = heatSource.heatRateW * physics.thermalSourceScale;
        channel.minimumPowerW = 0.0;
        channel.maximumPowerW = 1.5 * channel.nominalPowerW;
        channel.rhsPerWatt.resize(unitCase.size());
        for (std::size_t i = 0; i < unitCase.size(); ++i) {
            channel.rhsPerWatt[i] = unitCase[i] - result.referenceRhs[i];
        }
        result.channels.push_back(std::move(channel));
    }
    return result;
}

std::vector<double> composeRhs(const SourceParameterization& sources,
                               const std::vector<double>& powersW)
{
    if (powersW.size() != sources.channels.size()) {
        throw std::runtime_error("[MOR] Parameter vector has the wrong number of source channels.");
    }
    std::vector<double> rhs = sources.referenceRhs;
    for (std::size_t channel = 0; channel < sources.channels.size(); ++channel) {
        if (powersW[channel] == 0.0) {
            continue;
        }
        const auto& load = sources.channels[channel].rhsPerWatt;
        for (std::size_t row = 0; row < rhs.size(); ++row) {
            rhs[row] += powersW[channel] * load[row];
        }
    }
    return rhs;
}

std::vector<ParameterCase> makeParameterCases(const SourceParameterization& sources,
                                              const std::string& split,
                                              int count,
                                              std::uint64_t seed,
                                              bool includeUnitChannels)
{
    if (count < 0) {
        throw std::runtime_error("[MOR] Case count must be nonnegative.");
    }
    const int channels = static_cast<int>(sources.channels.size());
    const int target = includeUnitChannels ? std::max(count, channels) : count;
    std::vector<ParameterCase> result;
    result.reserve(static_cast<std::size_t>(target));
    if (includeUnitChannels) {
        for (int channel = 0; channel < channels; ++channel) {
            ParameterCase parameterCase;
            parameterCase.index = static_cast<int>(result.size());
            parameterCase.split = split;
            parameterCase.family = "unit_channel";
            parameterCase.powersW.assign(static_cast<std::size_t>(channels), 0.0);
            parameterCase.powersW[static_cast<std::size_t>(channel)] = 1.0;
            result.push_back(std::move(parameterCase));
        }
    }
    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    const int splitOffset = channels > 0
        ? static_cast<int>(seed % static_cast<std::uint64_t>(channels)) : 0;
    const double splitScale = 0.75
        + 0.05 * static_cast<double>((seed / 11ULL) % 9ULL);
    while (static_cast<int>(result.size()) < target) {
        ParameterCase parameterCase;
        parameterCase.index = static_cast<int>(result.size());
        parameterCase.split = split;
        parameterCase.powersW.assign(static_cast<std::size_t>(channels), 0.0);
        const int family = channels == 0 ? 0 : parameterCase.index % 7;
        if (family == 0) {
            parameterCase.family = "single_source";
            const int selected = (parameterCase.index + splitOffset)
                % std::max(1, channels);
            if (channels > 0) {
                parameterCase.powersW[static_cast<std::size_t>(selected)] =
                    splitScale
                    * sources.channels[static_cast<std::size_t>(selected)].nominalPowerW;
            }
        } else if (family == 1) {
            parameterCase.family = "double_source";
            for (int channel = 0; channel < channels; ++channel) {
                if (channel == (parameterCase.index + splitOffset) % channels
                    || channel == (parameterCase.index * 5 + 1 + splitOffset) % channels) {
                    parameterCase.powersW[static_cast<std::size_t>(channel)] =
                        splitScale
                        * sources.channels[static_cast<std::size_t>(channel)].nominalPowerW;
                }
            }
        } else if (family == 2) {
            parameterCase.family = "same_total_power";
            double nominalTotal = 0.0;
            for (const SourceChannel& channel : sources.channels) {
                nominalTotal += channel.nominalPowerW;
            }
            std::vector<double> weights(static_cast<std::size_t>(channels), 0.0);
            double weightSum = 0.0;
            for (double& weight : weights) {
                weight = 0.05 + uniform(generator);
                weightSum += weight;
            }
            for (int channel = 0; channel < channels; ++channel) {
                parameterCase.powersW[static_cast<std::size_t>(channel)] =
                    nominalTotal * weights[static_cast<std::size_t>(channel)]
                    / std::max(1.0e-300, weightSum);
            }
        } else if (family == 3) {
            parameterCase.family = "layer_pattern";
            for (int channel = 0; channel < channels; ++channel) {
                if (sources.channels[static_cast<std::size_t>(channel)].subdomain % 2
                    == (parameterCase.index + splitOffset) % 2) {
                    parameterCase.powersW[static_cast<std::size_t>(channel)] =
                        sources.channels[static_cast<std::size_t>(channel)].nominalPowerW;
                }
            }
        } else if (family == 4) {
            parameterCase.family = "high_low";
            for (int channel = 0; channel < channels; ++channel) {
                const double factor = ((channel + parameterCase.index + splitOffset) % 2 == 0)
                    ? 1.4 : 0.2;
                parameterCase.powersW[static_cast<std::size_t>(channel)] = factor
                    * sources.channels[static_cast<std::size_t>(channel)].nominalPowerW;
            }
        } else if (family == 5) {
            parameterCase.family = "uniform_power";
            const double uniformFactor = 0.65
                + 0.05 * static_cast<double>((seed / 17ULL) % 7ULL);
            for (int channel = 0; channel < channels; ++channel) {
                parameterCase.powersW[static_cast<std::size_t>(channel)] =
                    uniformFactor
                    * sources.channels[static_cast<std::size_t>(channel)].nominalPowerW;
            }
        } else {
            parameterCase.family = "seeded_random";
            for (int channel = 0; channel < channels; ++channel) {
                const SourceChannel& source = sources.channels[static_cast<std::size_t>(channel)];
                parameterCase.powersW[static_cast<std::size_t>(channel)] =
                    source.minimumPowerW
                    + uniform(generator) * (source.maximumPowerW - source.minimumPowerW);
            }
        }
        result.push_back(std::move(parameterCase));
    }
    return result;
}

void writeSourceChannels(const SourceParameterization& sources,
                         const std::filesystem::path& path)
{
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("[MOR] Cannot write source-channel metadata: " + path.string());
    }
    out << "channel,subdomain,domain_entity,nominal_power_w,minimum_power_w,maximum_power_w\n";
    out << std::setprecision(17);
    for (const SourceChannel& channel : sources.channels) {
        out << channel.index << ',' << channel.subdomain << ',' << channel.domainEntity << ','
            << channel.nominalPowerW << ',' << channel.minimumPowerW << ','
            << channel.maximumPowerW << '\n';
    }
}

void writeParameterCases(const std::vector<ParameterCase>& cases,
                         const std::filesystem::path& path)
{
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("[MOR] Cannot write parameter cases: " + path.string());
    }
    std::size_t channels = cases.empty() ? 0 : cases.front().powersW.size();
    out << "case,split,family";
    for (std::size_t channel = 0; channel < channels; ++channel) {
        out << ",power_w_" << channel;
    }
    out << '\n' << std::setprecision(17);
    for (const ParameterCase& parameterCase : cases) {
        out << parameterCase.index << ',' << parameterCase.split << ',' << parameterCase.family;
        for (double power : parameterCase.powersW) {
            out << ',' << power;
        }
        out << '\n';
    }
}

} // namespace mor
