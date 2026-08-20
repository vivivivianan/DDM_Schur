#pragma once

#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>

struct CaseConfig;
struct Mesh;

namespace mor::transient {

struct Options {
    bool generate = false;
    std::filesystem::path loadPath;
    std::filesystem::path savePath;
    int moments = 2;
    double expansionPoint = 0.0;
    double rankTolerance = 1.0e-10;
    // Optional second-moment preselection.  The default retains every first
    // moment direction, preserving the original Block Arnoldi construction.
    double secondMomentEnergy = 1.0;
    int secondMomentMaximumColumns = 0;
    // Non-data-driven local interior bases.  "standard" preserves the
    // validated M3 path; "shift-invert" changes only the Krylov resolvent.
    std::string basisType = "standard";
    // "galerkin" keeps W=V.  "petrov" uses the fixed-step CN local test
    // space W=(C_II/dt+0.5K_II)V without changing the Arnoldi trial basis.
    std::string projection = "galerkin";
    double basisShift = 0.0;
    // Adds C_II*1 to every local initial block. This is deliberately an
    // enrichment only; it does not alter the physical boundary conditions.
    bool boundaryAwareBasis = false;
    // ROM-only local residual enrichment.  A pilot uses the current reduced
    // Dynamic Schur solution only; it never invokes a FOM transient solve or
    // consumes trajectory snapshots.
    bool residualEnrichment = false;
    int residualTopN = 4;
    int residualEnrichmentRounds = 1;
    int residualPilotSteps = 200;
    std::string massType = "consistent";
    double timeStep = std::numeric_limits<double>::quiet_NaN();
    double endTime = std::numeric_limits<double>::quiet_NaN();
    std::string integrator = "backward-euler";
    std::filesystem::path inputPath;
    std::string waveform = "multi_step";
    std::string outputMode = "max-temperature";
    // Validation remains the compatibility default. Production mode disables
    // the eager FOM reference while retaining the exact residual gate.
    bool compareFom = true;
    // Run the factor-reused FOM comparator without emitting its final field.
    bool compareFomSummaryOnly = false;
    // Execute only the monolithic SIPG transient system and write Tmin/Tmax
    // at every stored time. No Arnoldi basis or reduced Schur system is built.
    bool fomOnly = false;
    // Form local and interface history RHS blocks from cached reduced state
    // instead of projecting the expanded global RHS back onto each basis.
    bool nativeReducedHistory = true;
    std::string interfaceInitialGuess = "previous";
    double fullResidualTolerance = 1.0e-5;
    // When false, retain and report an out-of-tolerance ROM state instead of
    // replacing the time step with a monolithic FOM solve.
    bool fullResidualFallback = false;
    double initialTemperature = std::numeric_limits<double>::quiet_NaN();
    std::string initialMode = "ambient";
    bool reuseIdenticalSubdomains = false;
    // global-fom preserves the historical full-order steady trace solves.
    // operator-coarse uses parallel local elimination, a geometric Schur
    // coarse correction, and source-particular PCG traces, so cold
    // construction needs no global factor.
    std::string constructionTraceMode = "global-fom";
    // Validation-only entry point: reuse the Local Block Arnoldi interior
    // basis and run a full-interface reduced Dynamic Schur solve.
    bool localInteriorArnoldiReducedSchurValidation = false;
    bool skipReducedSchurValidation = false;
    bool sourceAlignedInterfaceValidation = false;
    int interfaceExcitationRank = 0;
    int matrixFreeInterfaceThreshold = 20000;
    int interfaceMaxIterations = 1000;
    int interfaceRestart = 30;
    double interfaceTolerance = 1.0e-10;
    // Zero disables the production-only loose-first/strict-retry policy.
    double adaptiveInterfaceTolerance = 0.0;
    bool coarseLinearXY = true;
    bool coarseLinearZ = true;
    bool proxyEnabled = true;
    bool proxyDisableCoarse = false;
    double proxyHighConductivityThreshold = 20.0;
    bool proxyUseMaterialConnectivity = true;
    int proxyRing = 1;
    int proxyProbeColumns = 16;
    int proxyBlockSize = 64;
    bool proxyValidateBlockEquivalence = false;
    int localSolveThreads = 1;
    int localPardisoThreads = 1;
    std::string interfaceKrylov = "fgmres";
    bool portCoreCacheEnabled = false;
    std::filesystem::path portCoreCachePath;
    bool proxyCacheEnabled = false;
    std::filesystem::path proxyCachePath;
    int interfaceOperatorCoarseRank = 0;
    int interfaceOperatorCoarseSweeps = 2;
    bool interfaceOperatorCoarsePredictor = true;
    std::filesystem::path interfaceOperatorCoarseCachePath;
    bool portReduction = false;
    std::string portBasisMethod = "port-pod";
    int optimalPortRank = 16;
    std::string optimalPortRankMode = "fixed";
    std::filesystem::path optimalPortRankFile;
    double optimalPortEigenvalueTolerance = 1.0e-2;
    int optimalPortMinimumRank = 1;
    int optimalPortMaximumRank = 64;
    std::string optimalPortInnerProduct = "penalty-weighted-mass";
    int optimalPortOversamplingLayers = 0;
    std::string optimalPortInnerSolver = "auto";
    double optimalPortInnerTolerance = 1.0e-10;
    int optimalPortInnerMaximumIterations = 1000;
    int optimalPortInnerRefinementMaximumIterations = 3;
    double optimalPortInnerRefinementTolerance = 1.0e-10;
    double optimalPortEigenTolerance = 1.0e-8;
    int optimalPortEigenMaximumIterations = 40;
    std::string optimalPortAblation = "constant-geometry-trace";
    std::string optimalPortSourceMode = "trace-only";
    bool optimalPortTopologyAudit = false;
    bool optimalPortBasisPilot = false;
    bool optimalPortTargetSolverComparison = false;
    bool optimalPortWoodburyPilot = false;
    bool optimalPortRefinementValidation = false;
    bool optimalPortRepresentativeInterfacePilot = false;
    bool optimalPortMaximumInterfaceRefinementPilot = false;
    bool optimalPortAllInterfaceBasis = false;
    std::filesystem::path optimalPortTopologyAuditCsv;
    int residualKrylovMaximumRank = 4;
    int residualKrylovMaximumSweeps = 2;
    double residualKrylovTolerance = 1.0e-4;
    int residualKrylovBlockSize = 4;
    std::string residualKrylovProbeMode = "operator-geometry";
    std::string residualKrylovInnerSolver = "woodbury-exact";
    bool residualKrylovRepresentativePilot = false;
    bool residualKrylovAllInterfaceBasis = false;
    int randomizedPortRank = 8;
    int randomizedPortOversampling = 5;
    int randomizedPortPowerIterations = 1;
    std::uint64_t randomizedPortSeed = 12345;
    bool randomizedPortCompareOptimal = false;
    bool randomizedPortRepresentativePilot = false;
    std::string historyCompressionMethod = "none";
    int historyCompressionRank = 0;
    double historyCompressionTolerance = 1.0e-12;
    bool historyCompressionMaximumInterfacePilot = false;
    bool milestone8ProductionBasisOnly = false;
    bool milestone8AdaptiveProduction = false;
    bool adaptivePortLocalPilot = false;
    std::vector<int> adaptivePortInterfaceIds;
    bool globalInterfaceCoarsePrototype = false;
    std::string globalInterfaceCoarseInverseMode =
        "schur-jacobi";
    bool globalInterfaceCoarseExplicitReference = false;
    int globalInterfaceCoarseRank = 4;
    int globalInterfaceCoarseCandidateDimension = 0;
    int globalInterfaceCoarseMaximumIterations = 12;
    int globalInterfaceCoarseInnerMaximumIterations = 1000;
    int globalInterfaceCoarseKrylovSweeps = 2;
    double globalInterfaceCoarseTolerance = 1.0e-8;
    double globalInterfaceCoarseInnerTolerance = 1.0e-10;
    std::vector<int> globalInterfaceCoarseInterfaceIds;
    bool globalRandomizedSchur = false;
    int globalRandomizedRank = 50;
    std::uint64_t globalRandomizedSeed = 12345;
    int globalRandomizedInnerMaximumIterations = 1000;
    double globalRandomizedInnerTolerance = 1.0e-10;
    std::string globalRandomizedComposition = "global-only";
    bool projectionDiagnosis = false;
    bool fluxOperatorAudit = false;
    bool fluxAwarePort = false;
    std::string fluxAwareFluxType = "both";
    std::vector<int> projectionInterfaceIds;
    int localPortRank = 16;
    double localPortEnergyTolerance = 1.0e-8;
    std::filesystem::path localPortRankFile;
    double localPortTemperatureWeight = 1.0;
    double localPortFluxWeight = 1.0;
    double localPortResidualWeight = 1.0;
    int localPortEnrichmentRounds = 0;
    bool localPortCorrected = false;
    int deploymentRhsCount = 1;
    std::string basisOrthogonalization = "euclidean";
    std::uint64_t seed = 20260721ULL;
};

void runTransientBlockArnoldiWorkflow(
    const Mesh& mesh,
    const CaseConfig& physics,
    const Options& options,
    const std::filesystem::path& outputDirectory);

void runTransientDeploymentOnly(
    const std::filesystem::path& modelDirectory,
    const Options& options,
    const std::filesystem::path& outputDirectory);

} // namespace mor::transient
