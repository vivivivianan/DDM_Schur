#pragma once

#include "common.hpp"

#include <memory>
#include <vector>

namespace ddm_schur {

class SchurOperator {
public:
    SchurOperator(const Mesh& mesh,
                  const SparseMatrix& system,
                  const CaseConfig& physics,
                  bool coarseLinearXY,
                  bool coarseLinearZ,
                  bool coarseGlobalQuadraticZ,
                  bool coarseInterfacePatches,
                  bool coarseInterfaceLinearXY,
                  bool coarseEnergyAdaptive,
                  int energyMaxModesPerDomain,
                  int energySubspaceIterations,
                  double energyEigenvalueThreshold,
                  bool coarseGlobalSlow,
                  int globalSlowModes,
                  int globalSlowSubspaceDimension,
                  bool proxyEnabled,
                  bool proxyDisableCoarse,
                  bool proxyDiagnostics,
                  double proxyHighConductivityThreshold,
                  bool proxyUseMaterialConnectivity,
                  int proxyRing,
                  int proxyProbeColumns,
                  int proxyBlockSize,
                  bool proxyValidateBlockEquivalence,
                  int localSolveThreads,
                  int localPardisoThreads,
                  bool proxyCacheEnabled,
                  const std::string& proxyCachePath,
                  const std::string& proxyOutputDirectory);
    ~SchurOperator();

    SchurOperator(SchurOperator&&) noexcept;
    SchurOperator& operator=(SchurOperator&&) noexcept;
    SchurOperator(const SchurOperator&) = delete;
    SchurOperator& operator=(const SchurOperator&) = delete;

    int domains() const;
    int totalDofs() const;
    int interfaceDofs() const;
    int interiorDofs() const;
    int coarseDimension() const;
    int interfacePatchCount() const;
    int energyCandidateModes() const;
    double energyEigenSetupSeconds() const;
    double energySelectedEigenvalueMin() const;
    double energySelectedEigenvalueMax() const;
    int globalSlowCandidateDimension() const;
    double globalSlowSetupSeconds() const;
    double globalSlowEstimatedLambdaMax() const;
    double globalSlowSelectedRitzMin() const;
    double globalSlowSelectedRitzMax() const;
    int proxyGraphEdges() const;
    int proxyProbeColumns() const;
    int proxySchurApplies() const;
    double proxyDiagnosticsSeconds() const;
    double proxyRing3Coverage() const;
    double proxyRing3OperatorError() const;
    std::size_t proxyRing3EstimatedNnz() const;
    std::size_t proxyRing3MemoryEstimateBytes() const;
    bool proxyRecommended() const;
    std::size_t proxyNnz() const;
    double proxyDensity() const;
    int proxyColors() const;
    int proxyProbingSchurApplies() const;
    int proxyProbingBlockSize() const;
    int proxyProbingBlockCalls() const;
    int proxyValidationSchurApplies() const;
    int proxySymbolicCalls() const;
    int proxyNumericalCalls() const;
    int proxySolveCalls() const;
    double proxySetupSeconds() const;
    double proxySymbolicSeconds() const;
    double proxyNumericalSeconds() const;
    double proxySolveSeconds() const;
    double proxySymmetryError() const;
    double proxyMinimumTestRayleigh() const;
    double proxyDiagonalShift() const;
    double proxyDiagonalCompensation() const;
    std::uint64_t proxyValueHash() const;
    double proxyBlockMaximumDifference() const;
    double proxyBlockRelativeDifference() const;
    std::size_t proxyMemoryBytes() const;
    bool proxyMatrixCacheHit() const;
    bool proxyFactorCacheHit() const;
    double localFactorizationSeconds() const;
    double localSymbolicAnalysisSeconds() const;
    double localNumericalFactorizationSeconds() const;
    double localSolveSeconds() const;
    double schurApplySeconds() const;
    double coarseSolveSeconds() const;
    std::size_t memoryBytes() const;
    int localSolveCalls() const;
    int localSymbolicAnalysisCalls() const;
    int localNumericalFactorizationCalls() const;
    int matvecCalls() const;
    std::vector<SubdomainPerformance> subdomainPerformance() const;
    const std::vector<int>& interfaceGlobalDofs() const;

    std::vector<double> condensedRhs(const std::vector<double>& globalRhs);
    bool prepareGlobalSlow(const std::vector<double>& arnoldiSeed);
    void resetRuntimeCounters();
    void apply(const std::vector<double>& interfaceVector, std::vector<double>& result);
    void applyBlockPreconditioner(const std::vector<double>& residual, std::vector<double>& result);
    std::vector<double> recover(const std::vector<double>& globalRhs,
                                const std::vector<double>& interfaceSolution);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace ddm_schur
