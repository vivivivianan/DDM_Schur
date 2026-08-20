#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

struct Mesh;
struct SparseMatrix;
struct CaseConfig;

namespace ddm_schur {

struct Options {
    int maxIterations = 500;
    int restart = 50;
    double relativeTolerance = 1.0e-12;
    bool coarseLinearXY = true;
    bool coarseLinearZ = true;
    bool coarseGlobalQuadraticZ = false;
    bool coarseInterfacePatches = false;
    bool coarseInterfaceLinearXY = false;
    bool coarseEnergyAdaptive = false;
    int energyMaxModesPerDomain = 8;
    int energySubspaceIterations = 12;
    double energyEigenvalueThreshold = 0.1;
    bool coarseGlobalSlow = false;
    int globalSlowModes = 6;
    int globalSlowSubspaceDimension = 100;
    bool proxyDiagnostics = false;
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
    std::string portCoreCachePath;
    bool proxyCacheEnabled = false;
    std::string proxyCachePath;
    int interfaceOperatorCoarseRank = 0;
    int interfaceOperatorCoarseSweeps = 2;
    bool interfaceOperatorCoarsePredictor = true;
    std::string interfaceOperatorCoarseCachePath;
    std::string proxyOutputDirectory;
};

struct SubdomainPerformance {
    int subdomainId = -1;
    int interiorDofs = 0;
    int interfaceDofs = 0;
    double phase11Seconds = 0.0;
    double phase22Seconds = 0.0;
    int phase33Calls = 0;
    double phase33Seconds = 0.0;
    std::size_t factorMemoryBytes = 0;
};

struct Report {
    int domains = 0;
    int totalDofs = 0;
    int interfaceDofs = 0;
    int interiorDofs = 0;
    int iterations = 0;
    int coarseDimension = 0;
    int interfacePatchCount = 0;
    int energyCandidateModes = 0;
    int globalSlowCandidateDimension = 0;
    int schurMatvecs = 0;
    int localSolveCalls = 0;
    int localSymbolicAnalysisCalls = 0;
    int localNumericalFactorizationCalls = 0;
    int proxyGraphEdges = 0;
    int proxyProbeColumns = 0;
    int proxySchurApplies = 0;
    double setupSeconds = 0.0;
    double energyEigenSetupSeconds = 0.0;
    double energySelectedEigenvalueMin = 0.0;
    double energySelectedEigenvalueMax = 0.0;
    double globalSlowSetupSeconds = 0.0;
    double globalSlowEstimatedLambdaMax = 0.0;
    double globalSlowSelectedRitzMin = 0.0;
    double globalSlowSelectedRitzMax = 0.0;
    double proxyDiagnosticsSeconds = 0.0;
    double proxyRing3Coverage = 0.0;
    double proxyRing3OperatorError = 0.0;
    std::size_t proxyRing3EstimatedNnz = 0;
    std::size_t proxyRing3MemoryEstimateBytes = 0;
    bool proxyRecommended = false;
    std::size_t proxyNnz = 0;
    double proxyDensity = 0.0;
    int proxyColors = 0;
    int proxyProbingSchurApplies = 0;
    int proxyProbingBlockSize = 0;
    int proxyProbingBlockCalls = 0;
    int proxyValidationSchurApplies = 0;
    int proxySymbolicCalls = 0;
    int proxyNumericalCalls = 0;
    int proxySolveCalls = 0;
    double proxySetupSeconds = 0.0;
    double proxySymbolicSeconds = 0.0;
    double proxyNumericalSeconds = 0.0;
    double proxySolveSeconds = 0.0;
    bool proxyMatrixCacheHit = false;
    bool proxyFactorCacheHit = false;
    double proxySymmetryError = 0.0;
    double proxyMinimumTestRayleigh = 0.0;
    double proxyDiagonalShift = 0.0;
    double proxyDiagonalCompensation = 0.0;
    std::uint64_t proxyValueHash = 0;
    double proxyBlockMaximumDifference = 0.0;
    double proxyBlockRelativeDifference = 0.0;
    std::size_t proxyMemoryBytes = 0;
    double localFactorizationSeconds = 0.0;
    double localSymbolicAnalysisSeconds = 0.0;
    double localNumericalFactorizationSeconds = 0.0;
    double localSolveSeconds = 0.0;
    double schurApplySeconds = 0.0;
    double coarseSolveSeconds = 0.0;
    double condensedRhsSeconds = 0.0;
    double interfaceSolveSeconds = 0.0;
    double fgmresSeconds = 0.0;
    double recoverySeconds = 0.0;
    double totalSolveSeconds = 0.0;
    double totalSeconds = 0.0;
    double interfaceRelativeResidual = 0.0;
    std::size_t memoryBytes = 0;
    std::vector<SubdomainPerformance> subdomainPerformance;
    std::string status = "not_run";
};

struct SolveResult {
    std::vector<double> temperature;
    // Exposed for the optional Stage-2 interface ROM.  The Stage-1 solve and
    // recovery paths are unchanged; this is the converged Schur unknown that
    // was already computed internally.
    std::vector<double> interfaceSolution;
    Report report;
};

} // namespace ddm_schur
