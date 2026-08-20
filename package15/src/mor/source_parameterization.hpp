#pragma once

#include "types.hpp"

#include <filesystem>
#include <vector>

struct CaseConfig;
struct Mesh;

namespace mor {

struct SourceParameterization {
    std::vector<double> referenceRhs;
    std::vector<SourceChannel> channels;
};

SourceParameterization buildSourceParameterization(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::vector<double>& assembledSource,
    const std::vector<double>& heatOnlySource,
    const std::vector<double>& fixedAdjust);

std::vector<double> composeRhs(const SourceParameterization& sources,
                               const std::vector<double>& powersW);

std::vector<ParameterCase> makeParameterCases(const SourceParameterization& sources,
                                              const std::string& split,
                                              int count,
                                              std::uint64_t seed,
                                              bool includeUnitChannels);

void writeSourceChannels(const SourceParameterization& sources,
                         const std::filesystem::path& path);
void writeParameterCases(const std::vector<ParameterCase>& cases,
                         const std::filesystem::path& path);

} // namespace mor
