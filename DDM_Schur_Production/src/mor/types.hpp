#pragma once

// Common MOR fingerprints, affine channel descriptions, timing records, and
// validation metrics. Structures are deliberately plain for stable binary/CSV
// serialization across cold and warm production runs.

#include <cstdint>
#include <limits>
#include <filesystem>
#include <string>
#include <vector>

namespace mor {

struct Options {
    bool generateModel = false;
    std::filesystem::path loadModelPath;
    std::filesystem::path saveModelPath;
    int rank = 0;
    double energyTolerance = 1.0e-10;
    double singularValueTolerance = 1.0e-12;
    int trainingCases = 60;
    int validationCases = 20;
    int testCases = 20;
    std::uint64_t seed = 20260721ULL;
    std::string mode = "pure";
    std::string snapshotSolver = "auto";
    bool exactInteriorRecovery = true;
    bool compareFom = false;
    std::vector<int> rankSweep{5, 10, 20, 40, 60, 80, 100};
    std::string interiorMode = "pardiso";
    int interiorRank = 0;
    double interiorEnergyTolerance = 1.0e-8;
    double interiorSingularValueTolerance = 1.0e-12;
    std::vector<int> interiorRankSweep{5, 10, 20, 40, 60, 80, 100, 125};
    bool deploymentOnly = false;
    bool precomputePowerResponse = false;
    std::string storagePrecision = "float64";
    bool reportIoTime = false;
    bool compareInteriorModes = false;
    int deploymentRhsCount = 1;

    // Stage 2B.1: exactly one matrix parameter plus arbitrary source powers.
    bool parametricGenerate = false;
    bool parametricDeploymentOnly = false;
    std::filesystem::path parametricLoadPath;
    std::filesystem::path parametricSavePath;
    std::string matrixParameter;
    int parameterSubdomain = -1;
    int parameterRegionId = -1;
    double parameterMinimum = std::numeric_limits<double>::quiet_NaN();
    double parameterMaximum = std::numeric_limits<double>::quiet_NaN();
    double parameterReference = std::numeric_limits<double>::quiet_NaN();
    double parameterValue = std::numeric_limits<double>::quiet_NaN();
    int parameterTrainingCount = 5;
    int parameterValidationCount = 4;
    int parameterTestCount = 4;
    int interfaceRank = 0;
    int localRank = 0;
    std::string parametricMode = "pure";
    bool allowExtrapolation = false;
    bool parametricAffineValidationOnly = false;
    std::vector<double> onlinePowersW;
};

struct SourceChannel {
    int index = -1;
    int subdomain = -1;
    int domainEntity = -1;
    double nominalPowerW = 0.0;
    double minimumPowerW = 0.0;
    double maximumPowerW = 0.0;
    std::vector<double> rhsPerWatt;
};

struct ParameterCase {
    int index = -1;
    std::string split;
    std::string family;
    std::vector<double> powersW;
};

struct Fingerprints {
    std::uint64_t mesh = 0;
    std::uint64_t system = 0;
    std::uint64_t interfaceOrdering = 0;
    std::uint64_t sources = 0;
};

struct SnapshotDatabase {
    int rows = 0;
    std::vector<ParameterCase> cases;
    // Column-major: value(row, column) = values[column * rows + row].
    std::vector<double> values;
};

struct PodResult {
    int rows = 0;
    int columns = 0;
    int numericalRank = 0;
    int selectedRank = 0;
    std::vector<double> singularValues;
    std::vector<double> retainedEnergy;
    std::vector<double> basis;
    double gramSeconds = 0.0;
    double eigenSeconds = 0.0;
    double basisSeconds = 0.0;
};

struct ReducedSchurModel {
    int interfaceDofs = 0;
    int rank = 0;
    int sourceChannels = 0;
    std::vector<double> referenceInterface;
    std::vector<double> basis;
    std::vector<double> schurBasis;
    std::vector<double> reducedOperator;
    std::vector<double> singularValues;
    Fingerprints fingerprints;
    double symmetryError = 0.0;
    double minimumSymmetricEigenvalue = 0.0;
    double maximumSymmetricEigenvalue = 0.0;
    double conditionEstimate = 0.0;
    double exactSchurApplySeconds = 0.0;
    double reducedAssemblySeconds = 0.0;
};

struct OnlineTiming {
    double rhsSeconds = 0.0;
    double condensationSeconds = 0.0;
    double projectionSeconds = 0.0;
    double reducedSolveSeconds = 0.0;
    double interfaceReconstructionSeconds = 0.0;
    double recoverySeconds = 0.0;
    double correctionSeconds = 0.0;
    double totalSeconds = 0.0;
};

struct ErrorMetrics {
    double interfaceRelativeResidual = 0.0;
    double interfaceRelativeL2 = 0.0;
    double globalRelativeResidual = 0.0;
    double relativeL2 = 0.0;
    double maximumAbsolute = 0.0;
    double meanAbsolute = 0.0;
    double maximumTemperatureError = 0.0;
    double heatFluxImbalance = 0.0;
    double interfaceJump = 0.0;
};

struct SolveResult {
    std::vector<double> interfaceTemperature;
    std::vector<double> temperature;
    OnlineTiming timing;
    ErrorMetrics error;
    int correctionIterations = 0;
    std::string status = "not_run";
};

struct WorkflowResult {
    std::vector<double> nominalTemperature;
    int selectedRank = 0;
    int snapshotCount = 0;
    int sourceChannels = 0;
    int correctionIterations = 0;
    double setupSeconds = 0.0;
    double nominalOnlineSeconds = 0.0;
    double nominalGlobalResidual = 0.0;
    double nominalInterfaceResidual = 0.0;
    std::string status = "not_run";
};

struct DeploymentDof {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    int subdomain = -1;
    int sourceVertex = -1;
};

struct InteriorResponseBlock {
    int subdomain = -1;
    int rank = 0;
    std::uint64_t dofOrderingHash = 0;
    std::vector<int> globalDofs;
    std::vector<int> directPowerChannels;
    std::vector<double> referenceTemperature;
    // Column-major: global power channel columns, local interior DOF rows.
    std::vector<double> exactResponse;
    // Column-major U (local DOFs x rank) and R (rank x power channels).
    std::vector<double> localBasis;
    std::vector<double> localCoordinateMap;
    std::vector<double> singularValues;
    std::vector<double> retainedEnergy;
};

struct DeploymentResponseModel {
    int formatVersion = 2;
    int globalDofs = 0;
    int interfaceDofs = 0;
    int interfaceRank = 0;
    int sourceChannels = 0;
    std::string interiorMode = "exact-response";
    std::string storagePrecision = "float64";
    Fingerprints fingerprints;
    std::uint64_t globalDofOrderingHash = 0;
    std::vector<DeploymentDof> dofs;
    std::vector<int> interfaceGlobalDofs;
    std::vector<double> interfaceReference;
    std::vector<double> interfacePowerResponse;
    std::vector<double> reducedInputOperator;
    std::vector<double> nominalPowers;
    std::vector<double> minimumPowers;
    std::vector<double> maximumPowers;
    std::vector<int> sourceSubdomains;
    std::vector<int> sourceDomainEntities;
    std::vector<InteriorResponseBlock> interiors;
    double responseConstructionSeconds = 0.0;
    double compressionSeconds = 0.0;
    double serializationSeconds = 0.0;
    double loadSeconds = 0.0;
    std::uint64_t deploymentFileBytes = 0;
    std::uint64_t completeModelFileBytes = 0;
};

struct DeploymentRunResult {
    int rhsCount = 0;
    double setupSeconds = 0.0;
    double averageComputeSeconds = 0.0;
    double averageInterfaceSeconds = 0.0;
    double averageInteriorSeconds = 0.0;
    double outputIoSeconds = 0.0;
    double endToEndSeconds = 0.0;
    std::size_t peakWorkingSetBytes = 0;
    std::uint64_t deploymentFileBytes = 0;
    std::vector<double> nominalTemperature;
};

} // namespace mor
