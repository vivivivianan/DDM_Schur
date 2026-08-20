#pragma once

// Compatibility types for the removed global randomized Schur experiment.

#include "global_interface_coarse.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace mor::transient {

struct GlobalRandomizedSchurOptions {
    int requestedRank = 50;
    std::uint64_t seed = 12345;
    int innerMaximumIterations = 1000;
    double innerTolerance = 1.0e-10;
    double orthogonalityTolerance = 1.0e-10;
    double deflationTolerance = 1.0e-12;
    double maximumBasisTimeSeconds = 600.0;
    std::string composition = "global-only";
};

struct GlobalRandomizedSchurDiagnostics {
    int globalInterfaceDofs = 0;
    int requestedRank = 0;
    int acceptedRank = 0;
    int deflatedColumns = 0;
    std::uint64_t seed = 0;
    std::string composition;
    std::string innerSolverRequested = "matrix-free-pcg";
    std::string innerSolverActual = "matrix-free-pcg";
    int globalSchurApplyCount = 0;
    int preconditionerApplyCount = 0;
    int pardisoPhase33Calls = 0;
    int globalRhsCount = 0;
    int globalInnerIterations = 0;
    int maximumInnerIterations = 0;
    double meanSolveTimeSeconds = 0.0;
    double maximumSolveTimeSeconds = 0.0;
    double totalSolveTimeSeconds = 0.0;
    double maximumTargetSolveResidual = 0.0;
    double schurResidual = 0.0;
    double orthogonalityError = 0.0;
    double basisBuildTimeSeconds = 0.0;
    double compressionRatio = 0.0;
    std::size_t baselineWorkingSetBytes = 0;
    std::size_t peakWorkingSetBytes = 0;
    std::size_t peakIncrementalMemoryBytes = 0;
    bool snapshotUsed = false;
    bool fomUsedForBasis = false;
    bool podUsed = false;
    bool svdUsed = false;
    bool trainingWaveformUsed = false;
    std::string status = "not_run";
    std::string failureReason;
};

struct GlobalRandomizedSchurResult {
    GlobalInterfaceCoarseModel model;
    GlobalRandomizedSchurDiagnostics diagnostics;
};

GlobalRandomizedSchurResult buildGlobalRandomizedSchurPortSpace(
    const LocalPortModel& localModel,
    LocalPortReducedSchurSolver& localPreconditioner,
    const ReducedDynamicSchurOperator& exactSchur,
    const std::vector<double>& interfaceMetricDiagonal,
    const GlobalRandomizedSchurOptions& options);

void writeGlobalRandomizedSchurDiagnostics(
    const GlobalRandomizedSchurResult& result,
    const std::string& caseName,
    const std::filesystem::path& outputDirectory);

} // namespace mor::transient
