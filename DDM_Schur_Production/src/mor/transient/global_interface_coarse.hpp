#pragma once

// Compatibility types for a removed global spectral-coarse experiment.
// No production implementation is linked; see removed_research_methods.cpp.

#include "local_port_reduced_schur.hpp"
#include "transfer_operator.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace mor::transient {

struct GlobalInterfaceCoarseOptions {
    std::string inverseMode = "schur-jacobi";
    bool explicitReference = false;
    int requestedRank = 4;
    int candidateDimension = 0;
    int maximumIterations = 12;
    int innerMaximumIterations = 1000;
    int historyRank = 64;
    int krylovSweeps = 2;
    double innerTolerance = 1.0e-10;
    double ritzTolerance = 1.0e-8;
    double orthogonalityTolerance = 1.0e-10;
    double deflationTolerance = 1.0e-12;
    std::vector<int> selectedInterfaceIds;
};

struct GlobalInterfaceCoarseDiagnostics {
    int requestedRank = 0;
    int acceptedRank = 0;
    int candidateRank = 0;
    int iterations = 0;
    int rawInputColumns = 0;
    int rawBoundaryColumns = 0;
    int rawHistoryColumns = 0;
    int compressedHistoryColumns = 0;
    int rawGraphColumns = 0;
    int independentSeedColumns = 0;
    int deflatedSeedColumns = 0;
    int schurApplyCount = 0;
    int preconditionerApplyCount = 0;
    int innerSolveCount = 0;
    int innerSolveIterations = 0;
    int metricFlooredRows = 0;
    double metricRawMinimum = 0.0;
    double metricRawMaximum = 0.0;
    double metricFloor = 0.0;
    double maximumInnerSolveResidual = 0.0;
    double innerSolveSeconds = 0.0;
    double symmetryError = 0.0;
    double maximumRitzResidual = 0.0;
    double maximumProjectedRitzResidual = 0.0;
    double maximumInverseMetricRitzResidual = 0.0;
    double maximumLocalCoarseOrthogonality = 0.0;
    double localCoarseOrthogonalityAbsolute = 0.0;
    double localCoarseOrthogonalityRelative = 0.0;
    double coarseSchurGramError = 0.0;
    double coarseMassGramError = 0.0;
    double projectorIdempotenceError = 0.0;
    double projectorSSymmetryError = 0.0;
    double initialProbeResidual = 0.0;
    double finalProbeResidual = 0.0;
    double basisSeconds = 0.0;
    std::size_t baselineWorkingSetBytes = 0;
    std::size_t peakWorkingSetBytes = 0;
    std::size_t peakIncrementalMemoryBytes = 0;
    std::string requestedSpectralOperator;
    std::string actualSpectralOperator;
    std::string jacobiRole;
    std::string eigenvalueMapping;
    std::string innerSolverRequested;
    std::string innerSolverActual;
    bool innerSolveConverged = true;
    bool explicitReferenceRequested = false;
    bool explicitReferenceRan = false;
    int explicitReferenceDimension = 0;
    double explicitReferenceFullResidual = 0.0;
    double explicitReferenceProjectedResidual = 0.0;
    double explicitEigenvalueRelativeError = 0.0;
    double explicitMaximumPrincipalAngleRadians = 0.0;
    double explicitReferenceSeconds = 0.0;
    std::string explicitReferenceStatus = "not_requested";
    std::string status = "not_run";
};

struct GlobalInterfaceCoarseModel {
    int formatVersion = 1;
    std::string method = "schur-spectral-prototype";
    int fullInterfaceDofs = 0;
    int rank = 0;
    std::vector<int> interfaceGlobalDofs;
    std::vector<int> selectedInterfaceIds;
    // Column-major fullInterfaceDofs x rank.
    std::vector<double> basis;
    std::vector<double> schurImages;
    std::vector<double> eigenvalues;
    std::vector<double> ritzResiduals;
    std::vector<double> projectedRitzResiduals;
    std::vector<double> inverseMetricRitzResiduals;
    std::vector<double> participationRatios;
    std::uint64_t localBasisFingerprint = 0;
    bool snapshotUsed = false;
    bool fomUsedForBasis = false;
    bool podUsed = false;
    bool svdUsed = false;
    GlobalInterfaceCoarseDiagnostics diagnostics;
};

GlobalInterfaceCoarseModel buildGlobalInterfaceCoarsePrototype(
    const LocalPortModel& localModel,
    LocalPortReducedSchurSolver& localSolver,
    const ReducedDynamicSchurOperator& exactSchur,
    const std::vector<PortPatch>& physicalPatches,
    const std::vector<double>& interfaceMetricDiagonal,
    const GeneralizedTransferSourceBlocks& sources,
    const GlobalInterfaceCoarseOptions& options);

void writeGlobalInterfaceCoarseDiagnostics(
    const GlobalInterfaceCoarseModel& model,
    const std::string& caseName,
    const std::filesystem::path& outputDirectory);

std::size_t globalCoarseCurrentWorkingSetBytes();

} // namespace mor::transient
