#pragma once

#include "block_arnoldi.hpp"

#include <string>
#include <vector>

namespace mor::transient {

struct ReducedMatrixDiagnostic {
    double symmetryError = 0.0;
    double minimumEigenvalue = 0.0;
    double maximumEigenvalue = 0.0;
    bool choleskySucceeded = false;
};

struct TransientReducedModel {
    int formatVersion = 1;
    int globalDofs = 0;
    int rank = 0;
    int blockSize = 0;
    int moments = 0;
    int sourceChannels = 0;
    double expansionPoint = 0.0;
    double rankTolerance = 1.0e-10;
    std::string massType = "consistent";
    std::vector<double> basis;
    std::vector<double> reducedCapacity;
    std::vector<double> reducedConductivity;
    std::vector<double> reducedInput;
    std::vector<double> reducedBoundary;
    std::vector<double> referenceTemperature;
    std::vector<double> nominalPowersW;
    std::vector<double> minimumPowersW;
    std::vector<double> maximumPowersW;
    std::vector<int> sourceSubdomains;
    std::vector<int> sourceDomainEntities;
    std::vector<DeploymentDof> deploymentDofs;
    std::vector<ArnoldiHistoryRow> arnoldiHistory;
    TransientFingerprints fingerprints;
    ReducedMatrixDiagnostic capacityDiagnostic;
    ReducedMatrixDiagnostic conductivityDiagnostic;
    double basisOrthogonalityError = 0.0;
    double referenceResidual = 0.0;
    double projectionSeconds = 0.0;
    BlockArnoldiTiming arnoldiTiming;
};

TransientReducedModel buildTransientReducedModel(
    const ThermalDescriptorSystem& descriptor,
    BlockArnoldiResult arnoldi);

std::vector<double> projectInitialConditionCWeighted(
    const ThermalDescriptorSystem& descriptor,
    const TransientReducedModel& model,
    const std::vector<double>& initialTemperature);

std::vector<double> reconstructTemperature(
    const TransientReducedModel& model,
    const std::vector<double>& coordinates);

} // namespace mor::transient
