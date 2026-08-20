#pragma once

#include "interface_operator.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <vector>

struct CaseConfig;
struct Mesh;

namespace ddm_schur {

struct ProxyDistanceMetric {
    int distance = 0;
    std::size_t entryCount = 0;
    std::size_t nonzeroCount = 0;
    double meanAbsolute = 0.0;
    double meanNonzeroAbsolute = 0.0;
    double maximumAbsolute = 0.0;
    double energyFraction = 0.0;
};

struct ProxyRingMetric {
    int ring = 0;
    double meanColumnNnz = 0.0;
    std::size_t estimatedNnz = 0;
    double estimatedDensity = 0.0;
    double energyCoverage = 0.0;
    double randomOperatorErrorMean = 0.0;
    double randomOperatorErrorMaximum = 0.0;
    std::size_t csrMemoryEstimateBytes = 0;
};

struct ProxyDiagnosticsResult {
    int interfaceDofs = 0;
    std::size_t graphEdges = 0;
    std::size_t materialStrongEdges = 0;
    int crossDomainHighKComponents = 0;
    int probeColumns = 0;
    int exactSchurApplies = 0;
    double setupSeconds = 0.0;
    bool proxyRecommended = false;
    std::vector<ProxyDistanceMetric> distanceMetrics;
    std::vector<ProxyRingMetric> ringMetrics;
};

using ExactSchurApply = std::function<void(
    const std::vector<double>&, std::vector<double>&)>;
using ExactSchurResponseConsumer = std::function<void(
    int, const std::vector<double>&)>;
using ExactSchurApplyBlock = std::function<void(
    const std::vector<int>&, int, int, const ExactSchurResponseConsumer&)>;

ProxyDiagnosticsResult runSchurProxyDiagnostics(
    const Mesh& mesh,
    const CaseConfig& physics,
    const InterfacePartition& partition,
    double highConductivityThreshold,
    bool useMaterialConnectivity,
    int requestedProbeColumns,
    const std::filesystem::path& outputDirectory,
    const ExactSchurApply& applyExactSchur);

class SchurProxyPreconditioner {
public:
    SchurProxyPreconditioner(const Mesh& mesh,
                             const CaseConfig& physics,
                             const InterfacePartition& partition,
                             double highConductivityThreshold,
                             bool useMaterialConnectivity,
                             int ring,
                             int probingBlockSize,
                             int pardisoThreads,
                             bool validateBlockEquivalence,
                             bool cacheEnabled,
                             const std::filesystem::path& cachePath,
                             const std::filesystem::path& outputDirectory,
                             const ExactSchurApply& applyExactSchur,
                             const ExactSchurApplyBlock& applyExactSchurBlock);
    ~SchurProxyPreconditioner();

    SchurProxyPreconditioner(SchurProxyPreconditioner&&) noexcept;
    SchurProxyPreconditioner& operator=(SchurProxyPreconditioner&&) noexcept;
    SchurProxyPreconditioner(const SchurProxyPreconditioner&) = delete;
    SchurProxyPreconditioner& operator=(const SchurProxyPreconditioner&) = delete;

    void solve(const std::vector<double>& rhs, std::vector<double>& solution);
    void resetRuntimeCounters();
    int ring() const;
    std::size_t nnz() const;
    double density() const;
    int colors() const;
    int probingSchurApplies() const;
    int probingBlockSize() const;
    int probingBlockCalls() const;
    int validationSchurApplies() const;
    double setupSeconds() const;
    double symbolicSeconds() const;
    double numericalSeconds() const;
    double solveSeconds() const;
    int solveCalls() const;
    int symbolicCalls() const;
    int numericalCalls() const;
    double symmetryError() const;
    double minimumTestRayleigh() const;
    double diagonalShift() const;
    double diagonalCompensation() const;
    std::uint64_t valueHash() const;
    double blockMaximumDifference() const;
    double blockRelativeDifference() const;
    std::size_t memoryBytes() const;
    bool matrixCacheHit() const;
    bool factorCacheHit() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace ddm_schur
