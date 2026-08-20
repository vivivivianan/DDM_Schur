#pragma once

#include "transient_input_waveform.hpp"
#include "transient_reduced_model.hpp"

#include <string>
#include <vector>

namespace mor::transient {

struct TransientOnlineTiming {
    double loadSeconds = 0.0;
    double reducedFactorizationSeconds = 0.0;
    double timeSteppingSeconds = 0.0;
    double reconstructionSeconds = 0.0;
    double outputSeconds = 0.0;
    double totalSeconds = 0.0;
    std::size_t peakWorkingSetBytes = 0;
    bool usedLdltFallback = false;
};

struct ReducedTrajectory {
    int rank = 0;
    int steps = 0;
    double timeStep = 0.0;
    std::vector<double> times;
    // Step-major reduced coordinates, including step zero.
    std::vector<double> states;
    std::vector<double> maximumTemperature;
    TransientOnlineTiming timing;
};

struct TransientAccuracyRow {
    std::string waveform;
    int step = 0;
    double time = 0.0;
    double relativeL2 = 0.0;
    double maximumAbsolute = 0.0;
    double fomMaximumTemperature = 0.0;
    double romMaximumTemperature = 0.0;
    double fullResidual = 0.0;
    double reducedResidual = 0.0;
    double energyBalanceError = 0.0;
};

struct TransientAccuracySummary {
    std::string waveform;
    double spaceTimeRelativeL2 = 0.0;
    double maximumAbsolute = 0.0;
    double maximumTemperatureRmse = 0.0;
    double peakTemperatureError = 0.0;
    double peakTimeError = 0.0;
    double finalRelativeL2 = 0.0;
    double maximumFullResidual = 0.0;
    double maximumReducedResidual = 0.0;
    double maximumEnergyBalanceError = 0.0;
    double fomFactorizationSeconds = 0.0;
    double fomTimeSteppingSeconds = 0.0;
    double romTimeSteppingSeconds = 0.0;
    double romReconstructionSeconds = 0.0;
    double speedup = 0.0;
    std::size_t fomFactorMemoryBytes = 0;
    std::vector<TransientAccuracyRow> rows;
};

ReducedTrajectory integrateReducedModel(
    const TransientReducedModel& model,
    const PowerWaveform& waveform,
    double timeStep,
    int timeSteps,
    const std::string& integrator,
    const std::vector<double>& initialCoordinates,
    bool computeMaximumTemperature = true);

TransientAccuracySummary compareWithFullOrder(
    const ThermalDescriptorSystem& descriptor,
    const TransientReducedModel& model,
    const PowerWaveform& waveform,
    const ReducedTrajectory& reduced,
    double timeStep,
    int timeSteps,
    const std::string& integrator,
    const std::vector<double>& initialTemperature);

} // namespace mor::transient
