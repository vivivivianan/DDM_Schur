#pragma once

#include "local_subdomain_model.hpp"

#include <cstddef>
#include <filesystem>
#include <memory>
#include <vector>

namespace ddm_schur { struct InterfacePartition; }

namespace mor::local {

struct PortCoreSolveResult {
    std::vector<double> interfaceTemperature;
    std::vector<std::vector<double>> localReducedCoordinates;
    double portForwardSeconds = 0.0;
    double coreSolveSeconds = 0.0;
    double portBackSubstitutionSeconds = 0.0;
    double reducedRelativeResidual = 0.0;
};

// Exact block elimination for a fixed full-interface Local ROM system.  Trace
// unknowns that belong to one physical port are eliminated independently;
// junction traces and all local ROM coordinates remain in a small global core.
class PortCoreSolver {
public:
    PortCoreSolver(const Model& model,
                   const ddm_schur::InterfacePartition& partition,
                   int localSolveThreads,
                   int corePardisoThreads,
                   const std::filesystem::path& outputDirectory,
                   const std::filesystem::path& cachePath = {});
    ~PortCoreSolver();

    PortCoreSolveResult solve(
        const std::vector<std::vector<double>>& projectedInteriorRhs,
        const std::vector<double>& interfaceRhs);

    int portCount() const;
    int separatorDofs() const;
    int coreDimension() const;
    std::size_t nonzeros() const;
    std::size_t memoryBytes() const;
    double setupSeconds() const;
    double symbolicSeconds() const;
    double numericalSeconds() const;
    int symbolicCalls() const;
    int numericalCalls() const;
    bool cacheHit() const;
    double cacheLoadSeconds() const;
    double cacheSaveSeconds() const;
    std::size_t cacheBytes() const;
    double partitionSeconds() const;
    double leafCsrSeconds() const;
    double leafFactorSeconds() const;
    double couplingAssemblySeconds() const;
    double eliminationSeconds() const;
    double multiRhsPortSeconds() const;
    double schurProductPortSeconds() const;
    double coreAccumulationSeconds() const;
    double coreCsrSeconds() const;
    double coreFactorSeconds() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace mor::local
