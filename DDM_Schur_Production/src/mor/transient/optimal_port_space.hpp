#pragma once

// Compatibility declarations for removed Steklov, optimal-transfer, and
// residual-Krylov port methods. Their implementations remain in the research
// repository only.

#include "local_port_reduced_schur.hpp"
#include "transfer_operator.hpp"

#include <filesystem>
#include <string>
#include <vector>

struct Mesh;
namespace ddm_schur { struct InterfacePartition; }
namespace mor { namespace local { struct Model; } }

namespace mor::transient {

// M8.1 direct-patch Steklov baseline.  M8.2 extends this surface with the
// matrix-free transfer operator while keeping the projected Schur solver and
// Local Block Arnoldi path unchanged.
struct SteklovPortOptions {
    int requestedRank = 16;
    int directRowLimit = 2048;
    int inverseIterations = 24;
    double relativeTolerance = 1.0e-10;
    std::string innerProduct = "penalty-weighted-mass";
};

LocalPortModel buildSteklovSchurPortModel(
    const Mesh& mesh,
    const ddm_schur::InterfacePartition& partition,
    const local::Model& dynamicModel,
    const std::vector<double>& traceMassDiagonal,
    const std::vector<double>& penaltyMassDiagonal,
    const SteklovPortOptions& options);

struct OptimalTransferPortOptions {
    int requestedRank = 16;
    // Negative means that requestedRank is the total port budget. A
    // nonnegative value fixes the enrichment count after mandatory
    // deflation; used by the interface-16 rank-4 pilot.
    int requestedTransferRank = -1;
    std::string rankMode = "fixed";
    std::filesystem::path rankFile;
    double eigenvalueTolerance = 1.0e-2;
    int minimumRank = 1;
    int maximumRank = 64;
    std::string innerProduct = "penalty-weighted-mass";
    int oversamplingLayers = 0;
    int eigensolverMaximumIterations = 40;
    double eigensolverTolerance = 1.0e-8;
    double relativeDeflationTolerance = 1.0e-12;
    std::string ablationMode = "constant-geometry-trace";
    std::string sourceMode = "trace-only";
    std::vector<int> selectedInterfaceIds;
    bool targetSolverPilotPreflight = false;
    bool columnConsistencyCheck = false;
    // A positive value is a hard wall-time budget for the pilot build.
    // Normal basis builds leave this disabled.
    double maximumPilotSeconds = 0.0;
    double maximumTargetSetupSeconds = 300.0;
    double maximumMeanTargetSolveSeconds = 5.0;
    double maximumOperatorCheckSeconds = 60.0;
    double maximumTargetResidual = 1.0e-9;
    double maximumWeightedAdjointError = 1.0e-8;
    std::size_t maximumIncrementalWorkspaceBytes =
        UINT64_C(2147483648);
    PatchTransferOptions innerSolver;
};

struct OptimalPortInterfaceDiagnostics {
    int interfaceId = -1;
    int targetRows = 0;
    int sourceRows = 0;
    int mandatoryRank = 0;
    int requestedTransferRank = 0;
    int convergedTransferRank = 0;
    int totalPortRank = 0;
    int traceSourceRows = 0;
    int inputSourceRows = 0;
    int boundarySourceRows = 0;
    int historySourceRows = 0;
    int eigenIterations = 0;
    int eigenOperatorApplies = 0;
    bool eigenConverged = false;
    std::string eigenStatus;
    double transferIndicator = 0.0;
    double gramMinimum = 0.0;
    double gramMaximum = 0.0;
    double gramRegularization = 0.0;
    double adjointRelativeError = 0.0;
    double explicitColumnReferenceError = 0.0;
    double maximumEigenpairResidual = 0.0;
    double targetSolvePreflightSeconds = 0.0;
    double operatorCheckSeconds = 0.0;
    double transferApplySeconds = 0.0;
    double transposeApplySeconds = 0.0;
    std::string pilotStatus = "not_requested";
    std::size_t peakWorkspaceBytes = 0;
    std::size_t eigensolverWorkspaceBytes = 0;
    std::size_t finalBasisStorageBytes = 0;
    std::vector<double> eigenvalues;
    PatchInnerSolverStatistics innerSolver;
};

struct OptimalPortBuildResult {
    LocalPortModel model;
    std::vector<OptimalPortInterfaceDiagnostics> interfaces;
    double patchSetupSeconds = 0.0;
    double gramAssemblySeconds = 0.0;
    double mandatoryModeSeconds = 0.0;
    double eigenSolveSeconds = 0.0;
    double orthogonalizationSeconds = 0.0;
    double totalSeconds = 0.0;
};

OptimalPortBuildResult buildOptimalTransferPortModel(
    const Mesh& mesh,
    const ddm_schur::InterfacePartition& partition,
    const local::Model& dynamicModel,
    const std::vector<double>& traceMassDiagonal,
    const std::vector<double>& penaltyMassDiagonal,
    const std::vector<double>& input,
    int sourceChannels,
    const std::vector<double>& boundaryLoad,
    const std::vector<double>& condensedHistory,
    int historyChannels,
    const OptimalTransferPortOptions& options,
    // The formal multi-interface/rank pilot passes one immutable operator
    // here so its local factors and selected Jacobi diagonal are constructed
    // once. Normal builds may omit it.
    const ReducedDynamicSchurOperator* sharedSchur = nullptr);

struct ResidualKrylovPortOptions {
    std::string basisMethod = "residual-krylov";
    std::string innerProduct = "penalty-weighted-mass";
    int oversamplingLayers = 0;
    int maximumEnrichmentRank = 4;
    int maximumSweeps = 2;
    int blockSize = 4;
    double residualTolerance = 1.0e-4;
    double relativeDeflationTolerance = 1.0e-12;
    std::string probeMode = "operator-geometry";
    // M8.6 compresses the operator-side Local Block Arnoldi history load
    // columns before any S_tt^{-1} particular response is requested.
    std::string historyCompressionMethod = "none";
    int historyCompressionRank = 0;
    double historyCompressionTolerance = 1.0e-12;
    std::vector<int> selectedInterfaceIds;
    PatchTransferOptions innerSolver;
};

struct ResidualKrylovInterfaceDiagnostics {
    int interfaceId = -1;
    int targetRows = 0;
    int sourceRows = 0;
    int constantRank = 0;
    int geometryRank = 0;
    int inputRank = 0;
    int boundaryRank = 0;
    int historyRank = 0;
    int rawHistoryChannels = 0;
    int activeHistoryChannels = 0;
    int requestedHistoryRank = 0;
    int compressedHistoryRank = 0;
    int deflatedHistoryChannels = 0;
    int historyTargetRightHandSides = 0;
    int mandatoryRankTotal = 0;
    int requestedRandomizedRank = 0;
    int acceptedRandomizedRank = 0;
    int rawProbeColumns = 0;
    int independentProbeColumns = 0;
    int deflatedProbeColumns = 0;
    int probeBlockRank = 0;
    int requestedEnrichmentRank = 0;
    int acceptedEnrichmentRank = 0;
    int enrichmentSweeps = 0;
    int targetSolveCount = 0;
    int schurApplyCount = 0;
    bool deflationLimited = false;
    double initialMaximumProbeResidual = 0.0;
    double finalMaximumProbeResidual = 0.0;
    double residualReductionFactor = 1.0;
    double historyCompressionRelativeError = 0.0;
    double historyCompressionSeconds = 0.0;
    double weightedAdjointError = 0.0;
    double totalSeconds = 0.0;
    std::size_t historyCompressionWorkspaceBytes = 0;
    std::size_t peakIncrementalMemoryBytes = 0;
    std::uint64_t historyCompressionFingerprint = 0;
    std::string historyCompressionMethod = "none";
    std::string status = "not_run";
    PatchInnerSolverStatistics innerSolver;
};

struct ResidualKrylovBuildResult {
    LocalPortModel model;
    std::vector<ResidualKrylovInterfaceDiagnostics> interfaces;
    double mandatoryModeSeconds = 0.0;
    double probeSetupSeconds = 0.0;
    double enrichmentSeconds = 0.0;
    double orthogonalizationSeconds = 0.0;
    double totalSeconds = 0.0;
};

ResidualKrylovBuildResult buildResidualKrylovPortModel(
    const Mesh& mesh,
    const ddm_schur::InterfacePartition& partition,
    const local::Model& dynamicModel,
    const std::vector<double>& traceMassDiagonal,
    const std::vector<double>& penaltyMassDiagonal,
    const std::vector<double>& input,
    int sourceChannels,
    const std::vector<double>& boundaryLoad,
    const std::vector<double>& condensedHistory,
    int historyChannels,
    const ResidualKrylovPortOptions& options,
    const ReducedDynamicSchurOperator* sharedSchur = nullptr,
    // Optional operator-only transfer basis.  Its columns are inserted after
    // mandatory modes and before residual enrichment, and are reorthogonalized
    // in the target metric so duplicate directions are deterministically
    // deflated.
    const LocalPortModel* initialTransferModel = nullptr);

} // namespace mor::transient
