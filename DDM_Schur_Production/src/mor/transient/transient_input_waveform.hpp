#pragma once

// Time-dependent source-channel interface. Waveform generation is independent
// of ROM construction, preventing online deployment inputs from becoming data.

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace mor::transient {

struct PowerWaveform {
    std::string name;
    int channels = 0;
    std::vector<double> times;
    // Step-major: values[step * channels + channel].
    std::vector<double> values;

    std::vector<double> sample(double time) const;
};

PowerWaveform makeBuiltinWaveform(
    const std::string& name,
    const std::vector<double>& nominalPowersW,
    double timeStep,
    int timeSteps,
    std::uint64_t seed = 20260721ULL);

PowerWaveform loadPowerWaveformCsv(
    const std::filesystem::path& path,
    int expectedChannels);

std::vector<std::string> builtinWaveformNames();

} // namespace mor::transient
