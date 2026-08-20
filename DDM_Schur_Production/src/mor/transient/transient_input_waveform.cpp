// Deterministic deployment waveforms. Production uses asynchronous hotspots;
// this module maps normalized channel amplitudes to physical source powers and
// writes only compact waveform metadata needed to reproduce a run.

#include "transient_input_waveform.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>

namespace mor::transient {
namespace {

std::vector<std::string> split(const std::string& line)
{
    std::vector<std::string> result;
    std::istringstream stream(line);
    std::string item;
    while (std::getline(stream, item, ',')) result.push_back(item);
    return result;
}

} // namespace

std::vector<double> PowerWaveform::sample(double time) const
{
    if (times.empty() || values.size() != times.size()
        * static_cast<std::size_t>(channels)) {
        throw std::runtime_error("Transient waveform is empty or inconsistent.");
    }
    if (time <= times.front()) {
        return std::vector<double>(values.begin(), values.begin() + channels);
    }
    if (time >= times.back()) {
        return std::vector<double>(values.end() - channels, values.end());
    }
    const auto upper = std::upper_bound(times.begin(), times.end(), time);
    const std::size_t right = static_cast<std::size_t>(
        std::distance(times.begin(), upper));
    const std::size_t left = right - 1;
    const double fraction = (time - times[left]) / (times[right] - times[left]);
    std::vector<double> result(static_cast<std::size_t>(channels), 0.0);
    for (int channel = 0; channel < channels; ++channel) {
        const double a = values[left * channels + static_cast<std::size_t>(channel)];
        const double b = values[right * channels + static_cast<std::size_t>(channel)];
        result[static_cast<std::size_t>(channel)] = a + fraction * (b - a);
    }
    return result;
}

std::vector<std::string> builtinWaveformNames()
{
    return {
        "single_step", "multi_step", "rectangular_pulse", "piecewise_multilevel",
        "variable_duty_cycle", "mixed_frequency", "asynchronous_hotspots",
        "remote_heating", "unseen_waveform", "unseen_channel_combination",
        "same_average_power"
    };
}

PowerWaveform makeBuiltinWaveform(
    const std::string& name,
    const std::vector<double>& nominalPowersW,
    double timeStep,
    int timeSteps,
    std::uint64_t seed)
{
    if (!(timeStep > 0.0) || timeSteps <= 0 || nominalPowersW.empty()) {
        throw std::runtime_error("Builtin waveform dimensions are invalid.");
    }
    const std::vector<std::string> names = builtinWaveformNames();
    if (std::find(names.begin(), names.end(), name) == names.end()) {
        throw std::runtime_error("Unknown builtin transient waveform: " + name);
    }
    PowerWaveform waveform;
    waveform.name = name;
    waveform.channels = static_cast<int>(nominalPowersW.size());
    waveform.times.resize(static_cast<std::size_t>(timeSteps + 1));
    waveform.values.assign(static_cast<std::size_t>(timeSteps + 1)
        * nominalPowersW.size(), 0.0);
    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<double> random(0.15, 1.25);
    for (int step = 0; step <= timeSteps; ++step) {
        const double t = static_cast<double>(step) * timeStep;
        const double fraction = static_cast<double>(step)
            / static_cast<double>(timeSteps);
        waveform.times[static_cast<std::size_t>(step)] = t;
        for (int channel = 0; channel < waveform.channels; ++channel) {
            const double nominal = nominalPowersW[static_cast<std::size_t>(channel)];
            double factor = 0.0;
            if (name == "single_step") {
                factor = channel == 0 && step > 0 ? 1.0 : 0.0;
            } else if (name == "multi_step") {
                factor = step > 0 ? (channel < std::min(4, waveform.channels) ? 1.0 : 0.25) : 0.0;
            } else if (name == "rectangular_pulse") {
                factor = fraction >= 0.2 && fraction <= 0.65
                    ? (channel % 3 == 0 ? 1.2 : 0.15) : 0.0;
            } else if (name == "piecewise_multilevel") {
                factor = fraction < 0.25 ? 0.15
                    : (fraction < 0.55 ? 1.0 : (fraction < 0.8 ? 0.45 : 1.25));
            } else if (name == "variable_duty_cycle") {
                const int period = std::max(2, 3 + channel % 7);
                const int duty = 1 + (step / std::max(1, period)) % period;
                factor = step % period < duty ? 1.0 : 0.05;
            } else if (name == "mixed_frequency") {
                constexpr double pi = 3.14159265358979323846;
                factor = 0.6 + 0.25 * std::sin(2.0 * pi * fraction * (1 + channel % 3))
                    + 0.12 * std::sin(2.0 * pi * fraction * (5 + channel % 5));
            } else if (name == "asynchronous_hotspots") {
                const double start = 0.08 * static_cast<double>(channel % 8);
                const double stop = std::min(1.0, start + 0.28 + 0.03 * (channel % 4));
                factor = fraction >= start && fraction <= stop ? 1.1 : 0.0;
            } else if (name == "remote_heating") {
                factor = channel == waveform.channels - 1 && step > 0 ? 1.0 : 0.0;
            } else if (name == "unseen_waveform") {
                factor = fraction < 0.18 ? 0.0
                    : (fraction < 0.42 ? 0.85 : (fraction < 0.73 ? 0.22 : 1.12));
                factor *= 0.8 + 0.2 * random(generator);
            } else if (name == "unseen_channel_combination") {
                factor = ((channel * 7 + 3) % 11 < 4)
                    ? (step > 0 ? 0.7 + 0.4 * random(generator) : 0.0) : 0.0;
            } else if (name == "same_average_power") {
                // Even/odd channel groups receive phase-shifted pulses with
                // the same time-average power but different distributions.
                factor = (step + (channel % 2) * 2) % 4 < 2 ? 1.2 : 0.2;
            }
            waveform.values[static_cast<std::size_t>(step) * waveform.channels
                + static_cast<std::size_t>(channel)] = std::max(0.0, factor) * nominal;
        }
    }
    return waveform;
}

PowerWaveform loadPowerWaveformCsv(
    const std::filesystem::path& path,
    int expectedChannels)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read transient input CSV: " + path.string());
    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Transient input CSV is empty.");
    const std::vector<std::string> header = split(line);
    if (header.size() != static_cast<std::size_t>(expectedChannels + 1)) {
        throw std::runtime_error("Transient input CSV channel count mismatch.");
    }
    PowerWaveform waveform;
    waveform.name = path.stem().string();
    waveform.channels = expectedChannels;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = split(line);
        if (fields.size() != header.size()) {
            throw std::runtime_error("Transient input CSV row width mismatch.");
        }
        waveform.times.push_back(std::stod(fields[0]));
        for (int channel = 0; channel < expectedChannels; ++channel) {
            waveform.values.push_back(std::stod(fields[static_cast<std::size_t>(channel + 1)]));
        }
    }
    if (waveform.times.empty()) throw std::runtime_error("Transient input CSV has no data rows.");
    return waveform;
}

} // namespace mor::transient
