#pragma once

#include "transient_time_integrator.hpp"

#include <filesystem>

namespace mor::transient {

struct DcConsistencyRow {
    std::string caseName;
    double relativeL2 = 0.0;
    double maximumAbsolute = 0.0;
    double fomMaximumTemperature = 0.0;
    double romMaximumTemperature = 0.0;
};

std::vector<DcConsistencyRow> evaluateDcConsistency(
    const ThermalDescriptorSystem& descriptor,
    const TransientReducedModel& model);

double initialProjectionOrthogonalityError(
    const ThermalDescriptorSystem& descriptor,
    const TransientReducedModel& model,
    const std::vector<double>& initialTemperature,
    const std::vector<double>& coordinates);

void writeArnoldiHistory(const TransientReducedModel& model,
                         const std::filesystem::path& path);
void writeAccuracyByTime(const TransientAccuracySummary& accuracy,
                         const std::filesystem::path& path,
                         bool append = false);
void writeAccuracyByWaveform(const std::vector<TransientAccuracySummary>& rows,
                             const std::filesystem::path& path);
void writeMaximumTemperatureCurves(const TransientAccuracySummary& accuracy,
                                   const std::filesystem::path& path,
                                   bool append = false);
void writeDcConsistency(const std::vector<DcConsistencyRow>& rows,
                        const std::filesystem::path& path);

} // namespace mor::transient
