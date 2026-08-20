#include "../sipg_core.hpp"
#include "interface_operator.hpp"
#include "local_solver.hpp"
#include "schur_operator.hpp"
#include "schur_proxy.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <exception>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <tuple>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_MKL_PARDISO
#include <mkl_lapacke.h>
#endif

namespace ddm_schur {

namespace {

struct CoarseVector {
    std::vector<int> indices;
    std::vector<double> values;
};

void appendInterfaceStencilDofs(const Tet& tet, std::vector<int>& dofs)
{
    dofs.insert(dofs.end(), tet.dof.begin(), tet.dof.end());
}

} // namespace

struct SchurOperator::Impl {
    InterfacePartition partition;
    SparseMatrix interfaceMatrix;
    std::vector<LocalSolver> interiorSolvers;
    std::vector<LocalSolver> blockSolvers;
    LocalSolver coarseSolver;
    std::vector<CoarseVector> coarseBasis;
    std::vector<CoarseVector> coarseImages;
    double factorSeconds = 0.0;
    double symbolicSeconds = 0.0;
    double numericalSeconds = 0.0;
    double localSeconds = 0.0;
    double schurApplySeconds = 0.0;
    double coarseSeconds = 0.0;
    int solveCalls = 0;
    int symbolicCalls = 0;
    int numericalCalls = 0;
    int matvecs = 0;
    int patchCount = 0;
    int energyCandidates = 0;
    double energySeconds = 0.0;
    double energyEigenMin = 0.0;
    double energyEigenMax = 0.0;
    int globalSlowCandidates = 0;
    double globalSlowSeconds = 0.0;
    double globalSlowLambdaMax = 0.0;
    double globalSlowRitzMin = 0.0;
    double globalSlowRitzMax = 0.0;
    bool globalSlowPending = false;
    int pendingGlobalSlowModes = 0;
    int pendingGlobalSlowSubspaceDimension = 0;
    ProxyDiagnosticsResult proxyDiagnostics;
    std::unique_ptr<SchurProxyPreconditioner> proxySolver;
    std::vector<SubdomainPerformance> domainPerformance;
    int localSolveThreads = 1;
    int localPardisoThreads = 1;
    bool blockSolversBuilt = false;
    std::size_t bytes = 0;

    Impl(const Mesh& mesh,
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
         bool buildProxy,
         bool proxyDisableCoarse,
         bool runProxyDiagnostics,
         double proxyHighConductivityThreshold,
         bool proxyUseMaterialConnectivity,
         int proxyRing,
         int proxyProbeColumns,
         int proxyBlockSize,
         bool proxyValidateBlockEquivalence,
         int requestedLocalSolveThreads,
         int requestedLocalPardisoThreads,
         bool proxyCacheEnabled,
         const std::string& proxyCachePath,
         const std::string& proxyOutputDirectory)
        : partition(buildInterfacePartition(mesh, system)),
          interfaceMatrix(static_cast<int>(partition.interfaceGlobalDofs.size())),
          localSolveThreads(std::max(1, requestedLocalSolveThreads)),
          localPardisoThreads(std::max(1, requestedLocalPardisoThreads))
    {
        for (const Entry& entry : partition.interfaceEntries) {
            interfaceMatrix.add(entry.row, entry.col, entry.value);
        }
        interfaceMatrix.finalizeCsr();

        interiorSolvers.reserve(partition.domains.size());
        blockSolversBuilt = !buildProxy || coarseEnergyAdaptive;
        if (blockSolversBuilt) {
            blockSolvers.reserve(partition.domains.size());
        }
        domainPerformance.resize(partition.domains.size());
        for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
            const DomainBlocks& domain = partition.domains[slot];
            interiorSolvers.emplace_back(
                static_cast<int>(domain.interiorGlobalDofs.size()), domain.interiorEntries,
                localPardisoThreads);
            if (blockSolversBuilt) {
                blockSolvers.emplace_back(
                    static_cast<int>(domain.interiorGlobalDofs.size()
                        + domain.interfaceGlobalDofs.size()),
                    domain.fullBlockEntries,
                    localPardisoThreads);
            }
            SubdomainPerformance& performance = domainPerformance[slot];
            performance.subdomainId = domain.domainId;
            performance.interiorDofs = static_cast<int>(domain.interiorGlobalDofs.size());
            performance.interfaceDofs = static_cast<int>(domain.interfaceGlobalDofs.size());
        }
        globalSlowPending = coarseGlobalSlow;
        pendingGlobalSlowModes = globalSlowModes;
        pendingGlobalSlowSubspaceDimension = globalSlowSubspaceDimension;
        if (runProxyDiagnostics) {
            proxyDiagnostics = runSchurProxyDiagnostics(
                mesh,
                physics,
                partition,
                proxyHighConductivityThreshold,
                proxyUseMaterialConnectivity,
                proxyProbeColumns,
                proxyOutputDirectory,
                [&](const std::vector<double>& input, std::vector<double>& output) {
                    applyRaw(input, output, false);
                });
        }
        if (buildProxy) {
            proxySolver = std::make_unique<SchurProxyPreconditioner>(
                mesh,
                physics,
                partition,
                proxyHighConductivityThreshold,
                proxyUseMaterialConnectivity,
                proxyRing,
                proxyBlockSize,
                localPardisoThreads,
                proxyValidateBlockEquivalence,
                proxyCacheEnabled,
                proxyCachePath,
                proxyOutputDirectory,
                [&](const std::vector<double>& input, std::vector<double>& output) {
                    applyRaw(input, output, false);
                },
                [&](const std::vector<int>& probeColors,
                    int firstColor,
                    int rightHandSides,
                    const ExactSchurResponseConsumer& consumeResponse) {
                    applyRawColorBlock(
                        probeColors, firstColor, rightHandSides, consumeResponse, false);
                });
        }
        if (!proxyDisableCoarse) {
            buildCoarseCorrection(mesh,
                                  coarseLinearXY,
                                  coarseLinearZ,
                                  coarseGlobalQuadraticZ,
                                  coarseInterfacePatches,
                                  coarseInterfaceLinearXY,
                                  coarseEnergyAdaptive,
                                  energyMaxModesPerDomain,
                                  energySubspaceIterations,
                                  energyEigenvalueThreshold);
        }
        for (const LocalSolver& solver : interiorSolvers) {
            symbolicSeconds += solver.symbolicAnalysisSeconds();
            numericalSeconds += solver.numericalFactorizationSeconds();
            symbolicCalls += solver.symbolicAnalysisCalls();
            numericalCalls += solver.numericalFactorizationCalls();
        }
        for (const LocalSolver& solver : blockSolvers) {
            symbolicSeconds += solver.symbolicAnalysisSeconds();
            numericalSeconds += solver.numericalFactorizationSeconds();
            symbolicCalls += solver.symbolicAnalysisCalls();
            numericalCalls += solver.numericalFactorizationCalls();
        }
        factorSeconds = symbolicSeconds + numericalSeconds;
        for (std::size_t slot = 0; slot < domainPerformance.size(); ++slot) {
            SubdomainPerformance& performance = domainPerformance[slot];
            performance.phase11Seconds = interiorSolvers[slot].symbolicAnalysisSeconds();
            performance.phase22Seconds = interiorSolvers[slot].numericalFactorizationSeconds();
            performance.factorMemoryBytes = interiorSolvers[slot].memoryBytes();
            if (blockSolversBuilt) {
                performance.phase11Seconds += blockSolvers[slot].symbolicAnalysisSeconds();
                performance.phase22Seconds += blockSolvers[slot].numericalFactorizationSeconds();
                performance.factorMemoryBytes += blockSolvers[slot].memoryBytes();
            }
        }

        refreshMemoryBytes();
        solveCalls = 0;
        matvecs = 0;
        localSeconds = 0.0;
        schurApplySeconds = 0.0;
        coarseSeconds = 0.0;
        for (SubdomainPerformance& performance : domainPerformance) {
            performance.phase33Calls = 0;
            performance.phase33Seconds = 0.0;
        }
    }

    double solveDomain(std::size_t slot,
                       LocalSolver& solver,
                       const std::vector<double>& rhs,
                       std::vector<double>& solution)
    {
        const auto start = std::chrono::steady_clock::now();
        solver.solve(rhs, solution);
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        domainPerformance[slot].phase33Seconds += seconds;
        ++domainPerformance[slot].phase33Calls;
        return seconds;
    }

    double solveDomainMultiple(std::size_t slot,
                               LocalSolver& solver,
                               const std::vector<double>& rhs,
                               int rightHandSides,
                               std::vector<double>& solution)
    {
        const auto start = std::chrono::steady_clock::now();
        solver.solveMultiple(rhs, rightHandSides, solution);
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        domainPerformance[slot].phase33Seconds += seconds;
        ++domainPerformance[slot].phase33Calls;
        return seconds;
    }

    void solveLocal(LocalSolver& solver,
                    const std::vector<double>& rhs,
                    std::vector<double>& solution)
    {
        const auto start = std::chrono::steady_clock::now();
        solver.solve(rhs, solution);
        localSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        ++solveCalls;
    }

    void solveLocalMultiple(LocalSolver& solver,
                            const std::vector<double>& rhs,
                            int rightHandSides,
                            std::vector<double>& solution)
    {
        const auto start = std::chrono::steady_clock::now();
        solver.solveMultiple(rhs, rightHandSides, solution);
        localSeconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        ++solveCalls;
    }

    void solveCoarse(const std::vector<double>& rhs, std::vector<double>& solution)
    {
        const auto start = std::chrono::steady_clock::now();
        coarseSolver.solve(rhs, solution);
        coarseSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }

    void applyRaw(const std::vector<double>& interfaceVector,
                  std::vector<double>& result,
                  bool countMatvec)
    {
        const auto applyStart = std::chrono::steady_clock::now();
        result = interfaceMatrix.multiply(interfaceVector);
        const int domainCount = static_cast<int>(partition.domains.size());
        std::vector<std::vector<double>> corrections(static_cast<std::size_t>(domainCount));
        std::vector<double> solveTimes(static_cast<std::size_t>(domainCount), 0.0);
        std::vector<std::exception_ptr> errors(static_cast<std::size_t>(domainCount));
#pragma omp parallel for num_threads(localSolveThreads) if(localSolveThreads > 1) schedule(static)
        for (int slotIndex = 0; slotIndex < domainCount; ++slotIndex) {
            try {
                const std::size_t slot = static_cast<std::size_t>(slotIndex);
                const DomainBlocks& domain = partition.domains[slot];
                std::vector<double> interiorProduct(domain.interiorGlobalDofs.size(), 0.0);
                for (std::size_t row = 0; row < domain.interiorInterfaceRows.size(); ++row) {
                    for (const auto& entry : domain.interiorInterfaceRows[row]) {
                        interiorProduct[row] += entry.second
                            * interfaceVector[static_cast<std::size_t>(entry.first)];
                    }
                }
                std::vector<double> eliminated;
                solveTimes[slot] = solveDomain(
                    slot, interiorSolvers[slot], interiorProduct, eliminated);
                std::vector<double>& localCorrection = corrections[slot];
                localCorrection.assign(domain.interfaceGlobalDofs.size(), 0.0);
                for (std::size_t localGamma = 0;
                     localGamma < domain.interfaceGlobalDofs.size(); ++localGamma) {
                    for (const auto& entry : domain.interfaceInteriorRows[localGamma]) {
                        localCorrection[localGamma] += entry.second
                            * eliminated[static_cast<std::size_t>(entry.first)];
                    }
                }
            } catch (...) {
                errors[static_cast<std::size_t>(slotIndex)] = std::current_exception();
            }
        }
        for (const std::exception_ptr& error : errors) {
            if (error) {
                std::rethrow_exception(error);
            }
        }
        localSeconds += std::accumulate(solveTimes.begin(), solveTimes.end(), 0.0);
        solveCalls += domainCount;
        for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
            const DomainBlocks& domain = partition.domains[slot];
            for (std::size_t localGamma = 0;
                 localGamma < domain.interfaceGlobalDofs.size(); ++localGamma) {
                const int gamma = partition.globalToInterface[
                    static_cast<std::size_t>(domain.interfaceGlobalDofs[localGamma])];
                result[static_cast<std::size_t>(gamma)] -= corrections[slot][localGamma];
            }
        }
        if (countMatvec) {
            ++matvecs;
            schurApplySeconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - applyStart).count();
        }
    }

    void applyRawMultiple(const std::vector<double>& interfaceVectors,
                          int rightHandSides,
                          std::vector<double>& result,
                          bool countMatvec)
    {
        const int interfaceCount = static_cast<int>(partition.interfaceGlobalDofs.size());
        if (rightHandSides <= 0
            || interfaceVectors.size() != static_cast<std::size_t>(interfaceCount)
                * static_cast<std::size_t>(rightHandSides)) {
            throw std::runtime_error("[Schur] Invalid matrix-free multi-RHS dimensions.");
        }

        result.assign(interfaceVectors.size(), 0.0);
        for (int row = 0; row < interfaceCount; ++row) {
            for (int offset = interfaceMatrix.rowPtr[static_cast<std::size_t>(row)];
                 offset < interfaceMatrix.rowPtr[static_cast<std::size_t>(row + 1)];
                 ++offset) {
                const int column = interfaceMatrix.colInd[static_cast<std::size_t>(offset)];
                const double value = interfaceMatrix.values[static_cast<std::size_t>(offset)];
                for (int rhs = 0; rhs < rightHandSides; ++rhs) {
                    const std::size_t vectorOffset = static_cast<std::size_t>(rhs)
                        * static_cast<std::size_t>(interfaceCount);
                    result[vectorOffset + static_cast<std::size_t>(row)] += value
                        * interfaceVectors[vectorOffset + static_cast<std::size_t>(column)];
                }
            }
        }

        for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
            const DomainBlocks& domain = partition.domains[slot];
            const int interiorCount = static_cast<int>(domain.interiorGlobalDofs.size());
            std::vector<double> interiorProducts(
                static_cast<std::size_t>(interiorCount)
                    * static_cast<std::size_t>(rightHandSides),
                0.0);
            for (int row = 0; row < interiorCount; ++row) {
                for (const auto& entry :
                     domain.interiorInterfaceRows[static_cast<std::size_t>(row)]) {
                    for (int rhs = 0; rhs < rightHandSides; ++rhs) {
                        const std::size_t interfaceOffset = static_cast<std::size_t>(rhs)
                            * static_cast<std::size_t>(interfaceCount);
                        const std::size_t interiorOffset = static_cast<std::size_t>(rhs)
                            * static_cast<std::size_t>(interiorCount);
                        interiorProducts[interiorOffset + static_cast<std::size_t>(row)] +=
                            entry.second * interfaceVectors[
                                interfaceOffset + static_cast<std::size_t>(entry.first)];
                    }
                }
            }
            std::vector<double> eliminated;
            solveLocalMultiple(
                interiorSolvers[slot], interiorProducts, rightHandSides, eliminated);
            for (std::size_t localGamma = 0;
                 localGamma < domain.interfaceGlobalDofs.size(); ++localGamma) {
                const int gamma = partition.globalToInterface[static_cast<std::size_t>(
                    domain.interfaceGlobalDofs[localGamma])];
                for (int rhs = 0; rhs < rightHandSides; ++rhs) {
                    const std::size_t interfaceOffset = static_cast<std::size_t>(rhs)
                        * static_cast<std::size_t>(interfaceCount);
                    const std::size_t interiorOffset = static_cast<std::size_t>(rhs)
                        * static_cast<std::size_t>(interiorCount);
                    double correction = 0.0;
                    for (const auto& entry : domain.interfaceInteriorRows[localGamma]) {
                        correction += entry.second * eliminated[
                            interiorOffset + static_cast<std::size_t>(entry.first)];
                    }
                    result[interfaceOffset + static_cast<std::size_t>(gamma)] -= correction;
                }
            }
        }
        if (countMatvec) {
            matvecs += rightHandSides;
        }
    }

    // Exact matrix-free Schur action for distance-2 coloring probes.  The
    // input columns contain only zero/one values determined by probeColors,
    // so materializing an interfaceCount*nrhs dense input would duplicate
    // information and can add hundreds of MiB at block=64.  The local phase-33
    // calls still receive all nrhs columns together; only the read-only probe
    // representation is compressed.
    void applyRawColorBlock(const std::vector<int>& probeColors,
                            int firstColor,
                            int rightHandSides,
                            const ExactSchurResponseConsumer& consumeResponse,
                            bool countMatvec)
    {
        const int interfaceCount = static_cast<int>(partition.interfaceGlobalDofs.size());
        if (rightHandSides <= 0
            || probeColors.size() != static_cast<std::size_t>(interfaceCount)) {
            throw std::runtime_error("[Schur] Invalid colored multi-RHS probe dimensions.");
        }

        std::vector<std::vector<std::pair<int, double>>> interfaceProducts(
            static_cast<std::size_t>(rightHandSides));
        for (int row = 0; row < interfaceCount; ++row) {
            for (int offset = interfaceMatrix.rowPtr[static_cast<std::size_t>(row)];
                 offset < interfaceMatrix.rowPtr[static_cast<std::size_t>(row + 1)];
                 ++offset) {
                const int column = interfaceMatrix.colInd[static_cast<std::size_t>(offset)];
                const int rhs = probeColors[static_cast<std::size_t>(column)] - firstColor;
                if (rhs >= 0 && rhs < rightHandSides) {
                    interfaceProducts[static_cast<std::size_t>(rhs)].push_back({
                        row, interfaceMatrix.values[static_cast<std::size_t>(offset)]});
                }
            }
        }

        std::vector<std::vector<double>> eliminatedByDomain(partition.domains.size());
        const int domainCount = static_cast<int>(partition.domains.size());
        std::vector<double> solveTimes(static_cast<std::size_t>(domainCount), 0.0);
        std::vector<std::exception_ptr> errors(static_cast<std::size_t>(domainCount));
#pragma omp parallel for num_threads(localSolveThreads) if(localSolveThreads > 1) schedule(static)
        for (int slotIndex = 0; slotIndex < domainCount; ++slotIndex) {
            try {
                const std::size_t slot = static_cast<std::size_t>(slotIndex);
                const DomainBlocks& domain = partition.domains[slot];
                const int interiorCount = static_cast<int>(domain.interiorGlobalDofs.size());
                std::vector<double> interiorProducts(
                    static_cast<std::size_t>(interiorCount)
                        * static_cast<std::size_t>(rightHandSides),
                    0.0);
                for (int row = 0; row < interiorCount; ++row) {
                    for (const auto& entry :
                         domain.interiorInterfaceRows[static_cast<std::size_t>(row)]) {
                        const int rhs = probeColors[static_cast<std::size_t>(entry.first)]
                            - firstColor;
                        if (rhs >= 0 && rhs < rightHandSides) {
                            interiorProducts[static_cast<std::size_t>(rhs)
                                * static_cast<std::size_t>(interiorCount)
                                + static_cast<std::size_t>(row)] += entry.second;
                        }
                    }
                }
                solveTimes[slot] = solveDomainMultiple(
                    slot, interiorSolvers[slot], interiorProducts, rightHandSides,
                    eliminatedByDomain[slot]);
            } catch (...) {
                errors[static_cast<std::size_t>(slotIndex)] = std::current_exception();
            }
        }
        for (const std::exception_ptr& error : errors) {
            if (error) {
                std::rethrow_exception(error);
            }
        }
        localSeconds += std::accumulate(solveTimes.begin(), solveTimes.end(), 0.0);
        solveCalls += domainCount;

        std::vector<double> result(static_cast<std::size_t>(interfaceCount), 0.0);
        for (int rhs = 0; rhs < rightHandSides; ++rhs) {
            std::fill(result.begin(), result.end(), 0.0);
            for (const auto& entry : interfaceProducts[static_cast<std::size_t>(rhs)]) {
                result[static_cast<std::size_t>(entry.first)] += entry.second;
            }
            for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
                const DomainBlocks& domain = partition.domains[slot];
                const int interiorCount = static_cast<int>(domain.interiorGlobalDofs.size());
                const std::vector<double>& eliminated = eliminatedByDomain[slot];
            for (std::size_t localGamma = 0;
                 localGamma < domain.interfaceGlobalDofs.size(); ++localGamma) {
                const int gamma = partition.globalToInterface[static_cast<std::size_t>(
                    domain.interfaceGlobalDofs[localGamma])];
                    const std::size_t interiorOffset = static_cast<std::size_t>(rhs)
                        * static_cast<std::size_t>(interiorCount);
                    double correction = 0.0;
                    for (const auto& entry : domain.interfaceInteriorRows[localGamma]) {
                        correction += entry.second * eliminated[
                            interiorOffset + static_cast<std::size_t>(entry.first)];
                    }
                    result[static_cast<std::size_t>(gamma)] -= correction;
                }
            }
            consumeResponse(rhs, result);
        }
        if (countMatvec) {
            matvecs += rightHandSides;
        }
    }

    static double weightedDot(const std::vector<double>& left,
                              const std::vector<double>& right,
                              const std::vector<double>& weight)
    {
        double value = 0.0;
        for (std::size_t i = 0; i < left.size(); ++i) {
            value += weight[i] * left[i] * right[i];
        }
        return value;
    }

    static void weightedOrthonormalize(std::vector<std::vector<double>>& vectors,
                                       const std::vector<double>& weight)
    {
        for (std::size_t col = 0; col < vectors.size(); ++col) {
            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t previous = 0; previous < col; ++previous) {
                    const double projection = weightedDot(
                        vectors[col], vectors[previous], weight);
                    for (std::size_t row = 0; row < vectors[col].size(); ++row) {
                        vectors[col][row] -= projection * vectors[previous][row];
                    }
                }
            }
            const double normSquared = weightedDot(vectors[col], vectors[col], weight);
            if (!(normSquared > 1.0e-28) || !std::isfinite(normSquared)) {
                throw std::runtime_error(
                    "[Schur] Energy eigenspace lost rank during weighted orthogonalization.");
            }
            const double inverseNorm = 1.0 / std::sqrt(normSquared);
            for (double& value : vectors[col]) {
                value *= inverseNorm;
            }
        }
    }

    static double deterministicInitialValue(int globalDof, int mode)
    {
        std::uint64_t state = static_cast<std::uint64_t>(globalDof + 1)
            * UINT64_C(0x9E3779B97F4A7C15);
        state ^= static_cast<std::uint64_t>(mode + 1) * UINT64_C(0xD1B54A32D192ED03);
        state ^= state >> 30;
        state *= UINT64_C(0xBF58476D1CE4E5B9);
        state ^= state >> 27;
        state *= UINT64_C(0x94D049BB133111EB);
        state ^= state >> 31;
        const double unit = static_cast<double>(state >> 11)
            * (1.0 / 9007199254740992.0);
        return 2.0 * unit - 1.0;
    }

    void applyLocalSchur(std::size_t slot,
                         const std::vector<int>& globalGammaToLocal,
                         const std::vector<double>& localInterface,
                         std::vector<double>& result)
    {
        const DomainBlocks& domain = partition.domains[slot];
        const int interiorCount = static_cast<int>(domain.interiorGlobalDofs.size());
        const int interfaceCount = static_cast<int>(domain.interfaceGlobalDofs.size());
        std::vector<double> interiorProduct(static_cast<std::size_t>(interiorCount), 0.0);
        for (int row = 0; row < interiorCount; ++row) {
            for (const auto& entry : domain.interiorInterfaceRows[static_cast<std::size_t>(row)]) {
                const int local = globalGammaToLocal[static_cast<std::size_t>(entry.first)];
                if (local >= 0) {
                    interiorProduct[static_cast<std::size_t>(row)] +=
                        entry.second * localInterface[static_cast<std::size_t>(local)];
                }
            }
        }
        std::vector<double> eliminated;
        solveLocal(interiorSolvers[slot], interiorProduct, eliminated);

        result.assign(static_cast<std::size_t>(interfaceCount), 0.0);
        for (const Entry& entry : domain.fullBlockEntries) {
            if (entry.row >= interiorCount && entry.col >= interiorCount) {
                result[static_cast<std::size_t>(entry.row - interiorCount)] +=
                    entry.value * localInterface[static_cast<std::size_t>(entry.col - interiorCount)];
            }
        }
        for (int row = 0; row < interfaceCount; ++row) {
            for (const auto& entry : domain.interfaceInteriorRows[static_cast<std::size_t>(row)]) {
                result[static_cast<std::size_t>(row)] -=
                    entry.second * eliminated[static_cast<std::size_t>(entry.first)];
            }
        }
    }

    void buildEnergyCoarseBasis(int maxModesPerDomain,
                                int subspaceIterations,
                                double eigenvalueThreshold)
    {
        const auto start = std::chrono::steady_clock::now();
        double selectedMin = std::numeric_limits<double>::infinity();
        double selectedMax = 0.0;

        for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
            const DomainBlocks& domain = partition.domains[slot];
            const int interiorCount = static_cast<int>(domain.interiorGlobalDofs.size());
            const int interfaceCount = static_cast<int>(domain.interfaceGlobalDofs.size());
            if (interfaceCount == 0) {
                continue;
            }
            const int selectedLimit = std::min(maxModesPerDomain, interfaceCount);
            const int candidateCount = std::min(interfaceCount, selectedLimit + 4);
            energyCandidates += candidateCount;

            std::vector<int> globalGammaToLocal(partition.interfaceGlobalDofs.size(), -1);
            std::vector<int> interfaceIndices(static_cast<std::size_t>(interfaceCount), -1);
            for (int local = 0; local < interfaceCount; ++local) {
                const int gamma = partition.globalToInterface[static_cast<std::size_t>(
                    domain.interfaceGlobalDofs[static_cast<std::size_t>(local)])];
                interfaceIndices[static_cast<std::size_t>(local)] = gamma;
                globalGammaToLocal[static_cast<std::size_t>(gamma)] = local;
            }

            std::vector<double> massDiagonal(static_cast<std::size_t>(interfaceCount), 0.0);
            for (const Entry& entry : domain.fullBlockEntries) {
                if (entry.row == entry.col && entry.row >= interiorCount) {
                    massDiagonal[static_cast<std::size_t>(entry.row - interiorCount)] +=
                        std::abs(entry.value);
                }
            }
            double positiveMean = 0.0;
            int positiveCount = 0;
            for (double value : massDiagonal) {
                if (value > 0.0 && std::isfinite(value)) {
                    positiveMean += value;
                    ++positiveCount;
                }
            }
            positiveMean = positiveCount > 0
                ? positiveMean / static_cast<double>(positiveCount)
                : 1.0;
            // Use a scalar-lumped algebraic interface mass.  Per-DOF Jacobi
            // scaling would normalize away the coefficient contrast that the
            // adaptive eigenproblem is meant to detect.
            for (double& value : massDiagonal) {
                value = positiveMean;
            }

            std::vector<std::vector<double>> subspace(
                static_cast<std::size_t>(candidateCount),
                std::vector<double>(static_cast<std::size_t>(interfaceCount), 0.0));
            for (int row = 0; row < interfaceCount; ++row) {
                subspace[0][static_cast<std::size_t>(row)] = 1.0;
                const int globalDof = domain.interfaceGlobalDofs[static_cast<std::size_t>(row)];
                for (int col = 1; col < candidateCount; ++col) {
                    subspace[static_cast<std::size_t>(col)][static_cast<std::size_t>(row)] =
                        deterministicInitialValue(globalDof, col);
                }
            }
            weightedOrthonormalize(subspace, massDiagonal);

            for (int iteration = 0; iteration < subspaceIterations; ++iteration) {
                std::vector<std::vector<double>> next(
                    static_cast<std::size_t>(candidateCount),
                    std::vector<double>(static_cast<std::size_t>(interfaceCount), 0.0));
                for (int col = 0; col < candidateCount; ++col) {
                    std::vector<double> fullRhs(
                        static_cast<std::size_t>(interiorCount + interfaceCount), 0.0);
                    for (int row = 0; row < interfaceCount; ++row) {
                        fullRhs[static_cast<std::size_t>(interiorCount + row)] =
                            massDiagonal[static_cast<std::size_t>(row)]
                            * subspace[static_cast<std::size_t>(col)][static_cast<std::size_t>(row)];
                    }
                    std::vector<double> fullSolution;
                    solveLocal(blockSolvers[slot], fullRhs, fullSolution);
                    for (int row = 0; row < interfaceCount; ++row) {
                        next[static_cast<std::size_t>(col)][static_cast<std::size_t>(row)] =
                            fullSolution[static_cast<std::size_t>(interiorCount + row)];
                    }
                }
                weightedOrthonormalize(next, massDiagonal);
                subspace = std::move(next);
            }

            std::vector<std::vector<double>> images(
                static_cast<std::size_t>(candidateCount));
            for (int col = 0; col < candidateCount; ++col) {
                applyLocalSchur(slot,
                                globalGammaToLocal,
                                subspace[static_cast<std::size_t>(col)],
                                images[static_cast<std::size_t>(col)]);
            }
            std::vector<double> projected(
                static_cast<std::size_t>(candidateCount * candidateCount), 0.0);
            for (int col = 0; col < candidateCount; ++col) {
                for (int row = 0; row <= col; ++row) {
                    double value = 0.0;
                    for (int local = 0; local < interfaceCount; ++local) {
                        value += subspace[static_cast<std::size_t>(row)][static_cast<std::size_t>(local)]
                            * images[static_cast<std::size_t>(col)][static_cast<std::size_t>(local)];
                    }
                    projected[static_cast<std::size_t>(row + col * candidateCount)] = value;
                    projected[static_cast<std::size_t>(col + row * candidateCount)] = value;
                }
            }
            std::vector<double> eigenvalues(static_cast<std::size_t>(candidateCount), 0.0);
#ifdef USE_MKL_PARDISO
            const lapack_int info = LAPACKE_dsyev(
                LAPACK_COL_MAJOR,
                'V',
                'U',
                static_cast<lapack_int>(candidateCount),
                projected.data(),
                static_cast<lapack_int>(candidateCount),
                eigenvalues.data());
            if (info != 0) {
                throw std::runtime_error(
                    "[Schur] Local generalized energy eigensolve failed in LAPACKE_dsyev.");
            }
#else
            throw std::runtime_error(
                "[Schur] Energy-adaptive coarse space requires MKL LAPACKE.");
#endif

            std::vector<std::vector<double>> selectedLocalModes;
            for (int eigen = 0; eigen < candidateCount
                 && static_cast<int>(selectedLocalModes.size()) < selectedLimit;
                 ++eigen) {
                const double lambda = std::max(0.0, eigenvalues[static_cast<std::size_t>(eigen)]);
                if (!selectedLocalModes.empty() && lambda > eigenvalueThreshold) {
                    break;
                }
                std::vector<double> mode(static_cast<std::size_t>(interfaceCount), 0.0);
                for (int col = 0; col < candidateCount; ++col) {
                    const double coefficient = projected[static_cast<std::size_t>(
                        col + eigen * candidateCount)];
                    for (int row = 0; row < interfaceCount; ++row) {
                        mode[static_cast<std::size_t>(row)] += coefficient
                            * subspace[static_cast<std::size_t>(col)][static_cast<std::size_t>(row)];
                    }
                }
                for (int pass = 0; pass < 2; ++pass) {
                    for (const CoarseVector& previous : coarseBasis) {
                        double projection = 0.0;
                        for (std::size_t i = 0; i < previous.indices.size(); ++i) {
                            const int local = globalGammaToLocal[static_cast<std::size_t>(
                                previous.indices[i])];
                            if (local >= 0) {
                                projection += mode[static_cast<std::size_t>(local)]
                                    * previous.values[i];
                            }
                        }
                        if (projection != 0.0) {
                            for (std::size_t i = 0; i < previous.indices.size(); ++i) {
                                const int local = globalGammaToLocal[static_cast<std::size_t>(
                                    previous.indices[i])];
                                if (local >= 0) {
                                    mode[static_cast<std::size_t>(local)] -= projection
                                        * previous.values[i];
                                }
                            }
                        }
                    }
                }
                double normSquared = 0.0;
                for (double value : mode) {
                    normSquared += value * value;
                }
                if (!(normSquared > 1.0e-24)) {
                    continue;
                }
                const double inverseNorm = 1.0 / std::sqrt(normSquared);
                CoarseVector coarse;
                coarse.indices = interfaceIndices;
                coarse.values.reserve(static_cast<std::size_t>(interfaceCount));
                for (double& value : mode) {
                    value *= inverseNorm;
                    coarse.values.push_back(value);
                }
                selectedLocalModes.push_back(std::move(mode));
                coarseBasis.push_back(std::move(coarse));
                selectedMin = std::min(selectedMin, lambda);
                selectedMax = std::max(selectedMax, lambda);
            }
        }

        if (std::isfinite(selectedMin)) {
            energyEigenMin = selectedMin;
            energyEigenMax = selectedMax;
        }
        energySeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
    }

    void applyOneLevelBlock(const std::vector<double>& residual,
                            std::vector<double>& result)
    {
        if (!blockSolversBuilt) {
            throw std::runtime_error(
                "[Schur] Full block solvers were not built for the selected proxy path.");
        }
        result.assign(residual.size(), 0.0);
        const int domainCount = static_cast<int>(partition.domains.size());
        std::vector<std::vector<double>> localSolutions(static_cast<std::size_t>(domainCount));
        std::vector<double> solveTimes(static_cast<std::size_t>(domainCount), 0.0);
        std::vector<std::exception_ptr> errors(static_cast<std::size_t>(domainCount));
#pragma omp parallel for num_threads(localSolveThreads) if(localSolveThreads > 1) schedule(static)
        for (int slotIndex = 0; slotIndex < domainCount; ++slotIndex) {
            try {
                const std::size_t slot = static_cast<std::size_t>(slotIndex);
                const DomainBlocks& domain = partition.domains[slot];
                const std::size_t interiorCount = domain.interiorGlobalDofs.size();
                std::vector<double> localRhs(
                    interiorCount + domain.interfaceGlobalDofs.size(), 0.0);
                for (std::size_t localGamma = 0;
                     localGamma < domain.interfaceGlobalDofs.size(); ++localGamma) {
                    const int gamma = partition.globalToInterface[static_cast<std::size_t>(
                        domain.interfaceGlobalDofs[localGamma])];
                    localRhs[interiorCount + localGamma] =
                        residual[static_cast<std::size_t>(gamma)];
                }
                solveTimes[slot] = solveDomain(
                    slot, blockSolvers[slot], localRhs, localSolutions[slot]);
            } catch (...) {
                errors[static_cast<std::size_t>(slotIndex)] = std::current_exception();
            }
        }
        for (const std::exception_ptr& error : errors) {
            if (error) {
                std::rethrow_exception(error);
            }
        }
        localSeconds += std::accumulate(solveTimes.begin(), solveTimes.end(), 0.0);
        solveCalls += domainCount;
        for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
            const DomainBlocks& domain = partition.domains[slot];
            const std::size_t interiorCount = domain.interiorGlobalDofs.size();
            for (std::size_t localGamma = 0;
                 localGamma < domain.interfaceGlobalDofs.size();
                 ++localGamma) {
                const int gamma = partition.globalToInterface[static_cast<std::size_t>(
                    domain.interfaceGlobalDofs[localGamma])];
                result[static_cast<std::size_t>(gamma)] =
                    localSolutions[slot][interiorCount + localGamma];
            }
        }
    }

    void applyLevelOne(const std::vector<double>& residual,
                       std::vector<double>& result)
    {
        if (proxySolver) {
            proxySolver->solve(residual, result);
        } else {
            applyOneLevelBlock(residual, result);
        }
    }

    static double euclideanDot(const std::vector<double>& left,
                               const std::vector<double>& right)
    {
        double value = 0.0;
        for (std::size_t i = 0; i < left.size(); ++i) {
            value += left[i] * right[i];
        }
        return value;
    }

    static double euclideanNormalize(std::vector<double>& vector)
    {
        const double norm = std::sqrt(std::max(0.0, euclideanDot(vector, vector)));
        if (norm > 0.0) {
            const double inverseNorm = 1.0 / norm;
            for (double& value : vector) {
                value *= inverseNorm;
            }
        }
        return norm;
    }

    void buildGlobalSlowCoarseBasis(int requestedModes,
                                    int requestedSubspaceDimension,
                                    const std::vector<double>& arnoldiSeed)
    {
        const auto start = std::chrono::steady_clock::now();
        const int interfaceCount = static_cast<int>(partition.interfaceGlobalDofs.size());
        const int arnoldiLimit = std::min(interfaceCount, requestedSubspaceDimension);
        if (interfaceCount == 0 || arnoldiLimit == 0) {
            return;
        }

        // Build a temporary copy of the already selected fixed coarse
        // correction.  The harmonic analysis must target the actual volume-xyz
        // baseline preconditioner, not the raw one-level operator.  This dense
        // coarse factor is setup-only and does not alter the final balanced
        // correction or the reusable PARDISO factors.
        const int fixedDimension = static_cast<int>(coarseBasis.size());
        std::vector<CoarseVector> fixedImages;
        std::vector<double> fixedFactor;
        std::vector<double> fixedEnergyDiagonal;
        if (fixedDimension > 0) {
            fixedImages.reserve(static_cast<std::size_t>(fixedDimension));
            fixedFactor.assign(
                static_cast<std::size_t>(fixedDimension * fixedDimension), 0.0);
            for (int col = 0; col < fixedDimension; ++col) {
                std::vector<double> basis(static_cast<std::size_t>(interfaceCount), 0.0);
                const CoarseVector& sparseBasis = coarseBasis[static_cast<std::size_t>(col)];
                for (std::size_t i = 0; i < sparseBasis.indices.size(); ++i) {
                    basis[static_cast<std::size_t>(sparseBasis.indices[i])] =
                        sparseBasis.values[i];
                }
                std::vector<double> image;
                applyRaw(basis, image, false);
                CoarseVector sparseImage;
                for (int gamma = 0; gamma < interfaceCount; ++gamma) {
                    const double value = image[static_cast<std::size_t>(gamma)];
                    if (value != 0.0) {
                        sparseImage.indices.push_back(gamma);
                        sparseImage.values.push_back(value);
                    }
                }
                for (int row = 0; row < fixedDimension; ++row) {
                    const CoarseVector& rowBasis = coarseBasis[static_cast<std::size_t>(row)];
                    double value = 0.0;
                    for (std::size_t i = 0; i < rowBasis.indices.size(); ++i) {
                        value += rowBasis.values[i]
                            * image[static_cast<std::size_t>(rowBasis.indices[i])];
                    }
                    fixedFactor[static_cast<std::size_t>(row + col * fixedDimension)] = value;
                }
                fixedImages.push_back(std::move(sparseImage));
            }
            for (int col = 0; col < fixedDimension; ++col) {
                for (int row = 0; row < col; ++row) {
                    const double average = 0.5 * (
                        fixedFactor[static_cast<std::size_t>(row + col * fixedDimension)]
                        + fixedFactor[static_cast<std::size_t>(col + row * fixedDimension)]);
                    fixedFactor[static_cast<std::size_t>(row + col * fixedDimension)] = average;
                    fixedFactor[static_cast<std::size_t>(col + row * fixedDimension)] = average;
                }
            }
            fixedEnergyDiagonal.resize(static_cast<std::size_t>(fixedDimension));
            for (int coarse = 0; coarse < fixedDimension; ++coarse) {
                fixedEnergyDiagonal[static_cast<std::size_t>(coarse)] =
                    fixedFactor[static_cast<std::size_t>(coarse + coarse * fixedDimension)];
            }
#ifdef USE_MKL_PARDISO
            const lapack_int factorInfo = LAPACKE_dpotrf(
                LAPACK_COL_MAJOR,
                'U',
                static_cast<lapack_int>(fixedDimension),
                fixedFactor.data(),
                static_cast<lapack_int>(fixedDimension));
            if (factorInfo != 0) {
                throw std::runtime_error(
                    "[Schur] Fixed coarse factorization failed during global harmonic analysis.");
            }
#else
            throw std::runtime_error("[Schur] Global slow coarse space requires MKL LAPACKE.");
#endif
        }

        auto solveFixedCoarse = [&](std::vector<double>& rhs) {
            if (fixedDimension == 0) {
                return;
            }
#ifdef USE_MKL_PARDISO
            const lapack_int solveInfo = LAPACKE_dpotrs(
                LAPACK_COL_MAJOR,
                'U',
                static_cast<lapack_int>(fixedDimension),
                1,
                fixedFactor.data(),
                static_cast<lapack_int>(fixedDimension),
                rhs.data(),
                static_cast<lapack_int>(fixedDimension));
            if (solveInfo != 0) {
                throw std::runtime_error(
                    "[Schur] Fixed coarse solve failed during global harmonic analysis.");
            }
#endif
        };

        auto applyFixedBalanced = [&](const std::vector<double>& residual,
                                      std::vector<double>& result) {
            if (fixedDimension == 0) {
                applyLevelOne(residual, result);
                return;
            }
            std::vector<double> coarseRhs(static_cast<std::size_t>(fixedDimension), 0.0);
            for (int coarse = 0; coarse < fixedDimension; ++coarse) {
                const CoarseVector& basis = coarseBasis[static_cast<std::size_t>(coarse)];
                for (std::size_t i = 0; i < basis.indices.size(); ++i) {
                    coarseRhs[static_cast<std::size_t>(coarse)] += basis.values[i]
                        * residual[static_cast<std::size_t>(basis.indices[i])];
                }
            }
            solveFixedCoarse(coarseRhs);
            std::vector<double> coarseCorrection(
                static_cast<std::size_t>(interfaceCount), 0.0);
            std::vector<double> coarseImage(
                static_cast<std::size_t>(interfaceCount), 0.0);
            for (int coarse = 0; coarse < fixedDimension; ++coarse) {
                const double coefficient = coarseRhs[static_cast<std::size_t>(coarse)];
                const CoarseVector& basis = coarseBasis[static_cast<std::size_t>(coarse)];
                for (std::size_t i = 0; i < basis.indices.size(); ++i) {
                    coarseCorrection[static_cast<std::size_t>(basis.indices[i])] +=
                        coefficient * basis.values[i];
                }
                const CoarseVector& image = fixedImages[static_cast<std::size_t>(coarse)];
                for (std::size_t i = 0; i < image.indices.size(); ++i) {
                    coarseImage[static_cast<std::size_t>(image.indices[i])] +=
                        coefficient * image.values[i];
                }
            }
            std::vector<double> projectedResidual(
                static_cast<std::size_t>(interfaceCount), 0.0);
            for (int gamma = 0; gamma < interfaceCount; ++gamma) {
                projectedResidual[static_cast<std::size_t>(gamma)] =
                    residual[static_cast<std::size_t>(gamma)]
                    - coarseImage[static_cast<std::size_t>(gamma)];
            }
            std::vector<double> localCorrection;
            applyLevelOne(projectedResidual, localCorrection);

            std::vector<double> leftRhs(static_cast<std::size_t>(fixedDimension), 0.0);
            for (int coarse = 0; coarse < fixedDimension; ++coarse) {
                const CoarseVector& image = fixedImages[static_cast<std::size_t>(coarse)];
                for (std::size_t i = 0; i < image.indices.size(); ++i) {
                    leftRhs[static_cast<std::size_t>(coarse)] += image.values[i]
                        * localCorrection[static_cast<std::size_t>(image.indices[i])];
                }
            }
            solveFixedCoarse(leftRhs);
            std::vector<double> leftCorrection(
                static_cast<std::size_t>(interfaceCount), 0.0);
            for (int coarse = 0; coarse < fixedDimension; ++coarse) {
                const double coefficient = leftRhs[static_cast<std::size_t>(coarse)];
                const CoarseVector& basis = coarseBasis[static_cast<std::size_t>(coarse)];
                for (std::size_t i = 0; i < basis.indices.size(); ++i) {
                    leftCorrection[static_cast<std::size_t>(basis.indices[i])] +=
                        coefficient * basis.values[i];
                }
            }
            result.resize(static_cast<std::size_t>(interfaceCount));
            for (int gamma = 0; gamma < interfaceCount; ++gamma) {
                result[static_cast<std::size_t>(gamma)] =
                    coarseCorrection[static_cast<std::size_t>(gamma)]
                    + localCorrection[static_cast<std::size_t>(gamma)]
                    - leftCorrection[static_cast<std::size_t>(gamma)];
            }
        };

        auto applyRightPreconditionedSchur = [&](const std::vector<double>& vector,
                                                 std::vector<double>& result) {
            std::vector<double> correction;
            applyFixedBalanced(vector, correction);
            applyRaw(correction, result, false);
        };

        if (static_cast<int>(arnoldiSeed.size()) != interfaceCount) {
            throw std::runtime_error(
                "[Schur] Global harmonic Arnoldi seed has the wrong size.");
        }
        std::vector<double> initial = arnoldiSeed;
        if (!(euclideanNormalize(initial) > 1.0e-24)) {
            throw std::runtime_error("[Schur] Global harmonic Arnoldi seed lost rank.");
        }

        std::vector<std::vector<double>> arnoldiBasis;
        arnoldiBasis.reserve(static_cast<std::size_t>(arnoldiLimit));
        arnoldiBasis.push_back(std::move(initial));
        std::vector<double> hessenberg(
            static_cast<std::size_t>(arnoldiLimit * arnoldiLimit), 0.0);
        int used = 0;
        double finalSubdiagonal = 0.0;
        for (int col = 0; col < arnoldiLimit; ++col) {
            std::vector<double> work;
            applyRightPreconditionedSchur(
                arnoldiBasis[static_cast<std::size_t>(col)], work);
            for (int pass = 0; pass < 2; ++pass) {
                for (int row = 0; row <= col; ++row) {
                    const double projection = euclideanDot(
                        arnoldiBasis[static_cast<std::size_t>(row)], work);
                    hessenberg[static_cast<std::size_t>(row + col * arnoldiLimit)] +=
                        projection;
                    for (int gamma = 0; gamma < interfaceCount; ++gamma) {
                        work[static_cast<std::size_t>(gamma)] -= projection
                            * arnoldiBasis[static_cast<std::size_t>(row)][static_cast<std::size_t>(gamma)];
                    }
                }
            }
            finalSubdiagonal = euclideanNormalize(work);
            used = col + 1;
            if (col + 1 >= arnoldiLimit || !(finalSubdiagonal > 1.0e-13)) {
                break;
            }
            hessenberg[static_cast<std::size_t>(col + 1 + col * arnoldiLimit)] =
                finalSubdiagonal;
            arnoldiBasis.push_back(std::move(work));
        }
        if (used < requestedModes) {
            throw std::runtime_error("[Schur] Global harmonic Arnoldi space is too small.");
        }
        globalSlowCandidates = used;

        std::vector<double> harmonicLeft(
            static_cast<std::size_t>(used * used), 0.0);
        std::vector<double> harmonicRight(
            static_cast<std::size_t>(used * used), 0.0);
        for (int col = 0; col < used; ++col) {
            for (int row = 0; row < used; ++row) {
                double value = 0.0;
                for (int inner = 0; inner < used; ++inner) {
                    value += hessenberg[static_cast<std::size_t>(inner + row * arnoldiLimit)]
                        * hessenberg[static_cast<std::size_t>(inner + col * arnoldiLimit)];
                }
                if (row == used - 1 && col == used - 1) {
                    value += finalSubdiagonal * finalSubdiagonal;
                }
                harmonicLeft[static_cast<std::size_t>(row + col * used)] = value;
                harmonicRight[static_cast<std::size_t>(row + col * used)] =
                    hessenberg[static_cast<std::size_t>(col + row * arnoldiLimit)];
            }
        }

        std::vector<double> alphaReal(static_cast<std::size_t>(used), 0.0);
        std::vector<double> alphaImaginary(static_cast<std::size_t>(used), 0.0);
        std::vector<double> beta(static_cast<std::size_t>(used), 0.0);
        std::vector<double> rightEigenvectors(
            static_cast<std::size_t>(used * used), 0.0);
#ifdef USE_MKL_PARDISO
        std::vector<double> unusedLeft(1, 0.0);
        const lapack_int eigenInfo = LAPACKE_dggev(
            LAPACK_COL_MAJOR,
            'N',
            'V',
            static_cast<lapack_int>(used),
            harmonicLeft.data(),
            static_cast<lapack_int>(used),
            harmonicRight.data(),
            static_cast<lapack_int>(used),
            alphaReal.data(),
            alphaImaginary.data(),
            beta.data(),
            unusedLeft.data(),
            1,
            rightEigenvectors.data(),
            static_cast<lapack_int>(used));
        if (eigenInfo != 0) {
            throw std::runtime_error(
                "[Schur] Global harmonic Ritz generalized eigensolve failed.");
        }
#else
        throw std::runtime_error("[Schur] Global slow coarse space requires MKL LAPACKE.");
#endif

        struct HarmonicGroup {
            int index = -1;
            int width = 1;
            double magnitude = 0.0;
        };
        std::vector<HarmonicGroup> groups;
        double candidateMaximum = 0.0;
        for (int eigen = 0; eigen < used; ++eigen) {
            const double denominator = beta[static_cast<std::size_t>(eigen)];
            if (!(std::abs(denominator) > 1.0e-14)) {
                continue;
            }
            const double real = alphaReal[static_cast<std::size_t>(eigen)] / denominator;
            const double imaginary = alphaImaginary[static_cast<std::size_t>(eigen)] / denominator;
            const double magnitude = std::hypot(real, imaginary);
            if (!std::isfinite(magnitude)) {
                continue;
            }
            candidateMaximum = std::max(candidateMaximum, magnitude);
            const double complexTolerance = 1.0e-10 * (1.0 + std::abs(real));
            if (imaginary > complexTolerance && eigen + 1 < used) {
                groups.push_back({eigen, 2, magnitude});
                ++eigen;
            } else if (imaginary >= -complexTolerance) {
                groups.push_back({eigen, 1, magnitude});
            }
        }
        std::sort(groups.begin(), groups.end(),
                  [](const HarmonicGroup& left, const HarmonicGroup& right) {
                      return left.magnitude < right.magnitude;
                  });

        double selectedMin = std::numeric_limits<double>::infinity();
        double selectedMax = 0.0;
        double maximumPhysicsSOverlap = 0.0;
        double maximumHarmonicSOverlap = 0.0;
        std::vector<std::vector<double>> selectedHarmonicModes;
        std::vector<std::vector<double>> selectedHarmonicImages;
        selectedHarmonicModes.reserve(static_cast<std::size_t>(requestedModes));
        selectedHarmonicImages.reserve(static_cast<std::size_t>(requestedModes));
        int selected = 0;
        for (const HarmonicGroup& group : groups) {
            for (int part = 0; part < group.width && selected < requestedModes; ++part) {
                const int vectorColumn = group.index + part;
                std::vector<double> residualMode(
                    static_cast<std::size_t>(interfaceCount), 0.0);
                for (int col = 0; col < used; ++col) {
                    const double coefficient = rightEigenvectors[static_cast<std::size_t>(
                        col + vectorColumn * used)];
                    for (int gamma = 0; gamma < interfaceCount; ++gamma) {
                        residualMode[static_cast<std::size_t>(gamma)] += coefficient
                            * arnoldiBasis[static_cast<std::size_t>(col)][static_cast<std::size_t>(gamma)];
                    }
                }
                std::vector<double> mode;
                applyFixedBalanced(residualMode, mode);
                std::vector<double> modeImage;
                applyRaw(mode, modeImage, false);

                // Energy orthogonalization is performed against the fixed
                // physics basis first and then against previously accepted
                // harmonic modes.  Updating mode and S*mode together avoids
                // an extra Schur apply on each reorthogonalization pass and
                // preserves the matrix-free operator exactly.
                for (int pass = 0; pass < 2; ++pass) {
                    if (fixedDimension > 0) {
                        std::vector<double> coefficients(
                            static_cast<std::size_t>(fixedDimension), 0.0);
                        for (int coarse = 0; coarse < fixedDimension; ++coarse) {
                            const CoarseVector& image =
                                fixedImages[static_cast<std::size_t>(coarse)];
                            for (std::size_t i = 0; i < image.indices.size(); ++i) {
                                coefficients[static_cast<std::size_t>(coarse)] +=
                                    image.values[i]
                                    * mode[static_cast<std::size_t>(image.indices[i])];
                            }
                        }
                        solveFixedCoarse(coefficients);
                        for (int coarse = 0; coarse < fixedDimension; ++coarse) {
                            const double coefficient =
                                coefficients[static_cast<std::size_t>(coarse)];
                            const CoarseVector& basis =
                                coarseBasis[static_cast<std::size_t>(coarse)];
                            for (std::size_t i = 0; i < basis.indices.size(); ++i) {
                                mode[static_cast<std::size_t>(basis.indices[i])] -=
                                    coefficient * basis.values[i];
                            }
                            const CoarseVector& image =
                                fixedImages[static_cast<std::size_t>(coarse)];
                            for (std::size_t i = 0; i < image.indices.size(); ++i) {
                                modeImage[static_cast<std::size_t>(image.indices[i])] -=
                                    coefficient * image.values[i];
                            }
                        }
                    }
                    for (std::size_t previous = 0;
                         previous < selectedHarmonicModes.size(); ++previous) {
                        const double coefficient = euclideanDot(
                            selectedHarmonicModes[previous], modeImage);
                        for (int gamma = 0; gamma < interfaceCount; ++gamma) {
                            mode[static_cast<std::size_t>(gamma)] -= coefficient
                                * selectedHarmonicModes[previous][static_cast<std::size_t>(gamma)];
                            modeImage[static_cast<std::size_t>(gamma)] -= coefficient
                                * selectedHarmonicImages[previous][static_cast<std::size_t>(gamma)];
                        }
                    }
                }

                const double energy = euclideanDot(mode, modeImage);
                if (!(energy > 1.0e-24) || !std::isfinite(energy)) {
                    continue;
                }
                const double inverseEnergyNorm = 1.0 / std::sqrt(energy);
                for (int gamma = 0; gamma < interfaceCount; ++gamma) {
                    mode[static_cast<std::size_t>(gamma)] *= inverseEnergyNorm;
                    modeImage[static_cast<std::size_t>(gamma)] *= inverseEnergyNorm;
                }

                for (int coarse = 0; coarse < fixedDimension; ++coarse) {
                    const CoarseVector& image = fixedImages[static_cast<std::size_t>(coarse)];
                    double overlap = 0.0;
                    for (std::size_t i = 0; i < image.indices.size(); ++i) {
                        overlap += image.values[i]
                            * mode[static_cast<std::size_t>(image.indices[i])];
                    }
                    const double fixedEnergy =
                        fixedEnergyDiagonal[static_cast<std::size_t>(coarse)];
                    if (fixedEnergy > 0.0) {
                        maximumPhysicsSOverlap = std::max(
                            maximumPhysicsSOverlap,
                            std::abs(overlap) / std::sqrt(fixedEnergy));
                    }
                }
                for (const std::vector<double>& previousImage : selectedHarmonicImages) {
                    maximumHarmonicSOverlap = std::max(
                        maximumHarmonicSOverlap,
                        std::abs(euclideanDot(mode, previousImage)));
                }

                CoarseVector coarse;
                coarse.indices.reserve(static_cast<std::size_t>(interfaceCount));
                coarse.values.reserve(static_cast<std::size_t>(interfaceCount));
                for (int gamma = 0; gamma < interfaceCount; ++gamma) {
                    const double value = mode[static_cast<std::size_t>(gamma)];
                    if (value != 0.0) {
                        coarse.indices.push_back(gamma);
                        coarse.values.push_back(value);
                    }
                }
                coarseBasis.push_back(std::move(coarse));
                selectedHarmonicModes.push_back(mode);
                selectedHarmonicImages.push_back(std::move(modeImage));
                selectedMin = std::min(selectedMin, group.magnitude);
                selectedMax = std::max(selectedMax, group.magnitude);
                ++selected;
            }
            if (selected >= requestedModes) {
                break;
            }
        }
        if (selected != requestedModes) {
            throw std::runtime_error(
                "[Schur] Global harmonic analysis did not produce enough independent modes.");
        }
        const double sOrthogonalityTolerance = 1.0e-8;
        if (maximumPhysicsSOverlap > sOrthogonalityTolerance
            || maximumHarmonicSOverlap > sOrthogonalityTolerance) {
            throw std::runtime_error(
                "[Schur] Global harmonic S-orthogonalization verification failed.");
        }
        std::cout << "[Schur] global harmonic S-orthogonality: physics="
                  << maximumPhysicsSOverlap
                  << ", harmonic=" << maximumHarmonicSOverlap << '\n';
        globalSlowLambdaMax = candidateMaximum;
        globalSlowRitzMin = selectedMin;
        globalSlowRitzMax = selectedMax;
        globalSlowSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
    }

    void buildCoarseCorrection(const Mesh& mesh,
                               bool coarseLinearXY,
                               bool coarseLinearZ,
                               bool coarseGlobalQuadraticZ,
                               bool coarseInterfacePatches,
                               bool coarseInterfaceLinearXY,
                               bool coarseEnergyAdaptive,
                               int energyMaxModesPerDomain,
                               int energySubspaceIterations,
                               double energyEigenvalueThreshold)
    {
        const bool addLinearXY = coarseInterfacePatches
            ? coarseInterfaceLinearXY
            : coarseLinearXY;
        const bool addLinearZ = !coarseInterfacePatches && coarseLinearZ;
        const bool addGlobalQuadraticZ = !coarseInterfacePatches
            && coarseGlobalQuadraticZ;
        std::vector<double> globalQuadraticZ(partition.interfaceGlobalDofs.size(), 0.0);

        if (coarseInterfacePatches) {
            using PatchKey = std::tuple<int, int, int, int>;
            std::map<PatchKey, std::vector<int>> patchDofs;
            for (const InterfaceFace& face : mesh.interfaceFaces) {
                if (face.leftTet < 0 || face.rightTet < 0) {
                    continue;
                }
                const Tet& left = mesh.tets[static_cast<std::size_t>(face.leftTet)];
                const Tet& right = mesh.tets[static_cast<std::size_t>(face.rightTet)];
                // A physical interface contributes one patch on each
                // subdomain boundary.  Keeping the two trace sides separate
                // lets the coarse space represent both mean and jump errors.
                const PatchKey leftKey{left.subdomain,
                                       right.subdomain,
                                       face.leftBoundaryEntity,
                                       face.rightBoundaryEntity};
                const PatchKey rightKey{right.subdomain,
                                        left.subdomain,
                                        face.rightBoundaryEntity,
                                        face.leftBoundaryEntity};
                // SIPG normal-flux couplings make every P2 DOF in the
                // boundary-tet stencil part of the algebraic Schur interface,
                // not only the six DOFs with nonzero geometric face trace.
                appendInterfaceStencilDofs(left, patchDofs[leftKey]);
                appendInterfaceStencilDofs(right, patchDofs[rightKey]);
            }

            // Make physical patches a disjoint partition of interface DOFs.
            // Shared edge/corner DOFs are assigned to the first deterministic
            // (domain pair, entity pair) key, avoiding dependent patch modes.
            std::vector<int> patchOwner(partition.interfaceGlobalDofs.size(), -1);
            for (auto& patch : patchDofs) {
                std::vector<int>& dofs = patch.second;
                std::sort(dofs.begin(), dofs.end());
                dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
                std::vector<int> ownedDofs;
                ownedDofs.reserve(dofs.size());
                for (int globalDof : dofs) {
                    const int gamma = partition.globalToInterface[static_cast<std::size_t>(globalDof)];
                    if (gamma >= 0 && patchOwner[static_cast<std::size_t>(gamma)] < 0) {
                        patchOwner[static_cast<std::size_t>(gamma)] = patchCount;
                        ownedDofs.push_back(globalDof);
                    }
                }
                if (ownedDofs.empty()) {
                    continue;
                }

                std::vector<int> interfaceIndices;
                interfaceIndices.reserve(ownedDofs.size());
                for (int globalDof : ownedDofs) {
                    interfaceIndices.push_back(
                        partition.globalToInterface[static_cast<std::size_t>(globalDof)]);
                }
                const std::size_t count = ownedDofs.size();
                CoarseVector constant;
                constant.indices = interfaceIndices;
                constant.values.assign(count, 1.0 / std::sqrt(static_cast<double>(count)));
                coarseBasis.push_back(std::move(constant));

                if (coarseInterfaceLinearXY) {
                    std::vector<std::vector<double>> localLinearModes;
                    for (int axis = 0; axis < 2; ++axis) {
                        double mean = 0.0;
                        for (int globalDof : ownedDofs) {
                            const Vec3& point = mesh.nodes[static_cast<std::size_t>(globalDof)].p;
                            mean += axis == 0 ? point.x : point.y;
                        }
                        mean /= static_cast<double>(count);
                        std::vector<double> values;
                        values.reserve(count);
                        for (int globalDof : ownedDofs) {
                            const Vec3& point = mesh.nodes[static_cast<std::size_t>(globalDof)].p;
                            values.push_back((axis == 0 ? point.x : point.y) - mean);
                        }
                        for (const std::vector<double>& previous : localLinearModes) {
                            double projection = 0.0;
                            for (std::size_t i = 0; i < count; ++i) {
                                projection += values[i] * previous[i];
                            }
                            for (std::size_t i = 0; i < count; ++i) {
                                values[i] -= projection * previous[i];
                            }
                        }
                        double normSquared = 0.0;
                        for (double value : values) {
                            normSquared += value * value;
                        }
                        if (normSquared > 1.0e-24) {
                            const double inverseNorm = 1.0 / std::sqrt(normSquared);
                            for (double& value : values) {
                                value *= inverseNorm;
                            }
                            CoarseVector linear;
                            linear.indices = interfaceIndices;
                            linear.values = values;
                            coarseBasis.push_back(std::move(linear));
                            localLinearModes.push_back(std::move(values));
                        }
                    }
                }
                ++patchCount;
            }
        }

        for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
            const DomainBlocks& domain = partition.domains[slot];
            if (domain.interfaceGlobalDofs.empty()) {
                continue;
            }
            if (coarseInterfacePatches) {
                continue;
            }

            const auto appendSubdomainModes = [&](const std::vector<int>& subdomainDofs) {
                const std::size_t count = subdomainDofs.size();
                if (count == 0) {
                    return;
                }
                std::vector<int> interfaceIndices;
                interfaceIndices.reserve(count);
                for (int globalDof : subdomainDofs) {
                    interfaceIndices.push_back(
                        partition.globalToInterface[static_cast<std::size_t>(globalDof)]);
                }

                CoarseVector constant;
                constant.indices = interfaceIndices;
                constant.values.assign(count, 1.0 / std::sqrt(static_cast<double>(count)));
                coarseBasis.push_back(std::move(constant));

                // Build a numerically stable local basis for span{1,x,y,z}.
                // Centering removes the constant component.  Modified
                // Gram-Schmidt removes correlations among coordinate modes on
                // thin or tilted interfaces without changing their span.
                std::vector<std::vector<double>> localLinearModes;
                const auto coordinate = [&](int globalDof, int axis) {
                    const Vec3& point = mesh.nodes[static_cast<std::size_t>(globalDof)].p;
                    return axis == 0 ? point.x : (axis == 1 ? point.y : point.z);
                };
                const auto appendLinearMode = [&](int axis) {
                    double mean = 0.0;
                    for (int globalDof : subdomainDofs) {
                        mean += coordinate(globalDof, axis);
                    }
                    mean /= static_cast<double>(count);

                    std::vector<double> values;
                    values.reserve(count);
                    for (int globalDof : subdomainDofs) {
                        values.push_back(coordinate(globalDof, axis) - mean);
                    }

                    for (const std::vector<double>& previous : localLinearModes) {
                        double projection = 0.0;
                        for (std::size_t i = 0; i < count; ++i) {
                            projection += values[i] * previous[i];
                        }
                        for (std::size_t i = 0; i < count; ++i) {
                            values[i] -= projection * previous[i];
                        }
                    }

                    double normSquared = 0.0;
                    for (double value : values) {
                        normSquared += value * value;
                    }
                    if (normSquared > 1.0e-24) {
                        const double inverseNorm = 1.0 / std::sqrt(normSquared);
                        for (double& value : values) {
                            value *= inverseNorm;
                        }
                        CoarseVector linear;
                        linear.indices = interfaceIndices;
                        linear.values = values;
                        coarseBasis.push_back(std::move(linear));
                        localLinearModes.push_back(std::move(values));
                        return true;
                    }
                    return false;
                };

                if (addLinearXY) {
                    appendLinearMode(0);
                    appendLinearMode(1);
                }
                const bool zModeAdded = addLinearZ && appendLinearMode(2);

                if (addGlobalQuadraticZ && zModeAdded) {
                    double zMean = 0.0;
                    for (int globalDof : subdomainDofs) {
                        zMean += coordinate(globalDof, 2);
                    }
                    zMean /= static_cast<double>(count);
                    std::vector<double> quadraticValues;
                    quadraticValues.reserve(count);
                    double quadraticMean = 0.0;
                    for (int globalDof : subdomainDofs) {
                        const double delta = coordinate(globalDof, 2) - zMean;
                        quadraticValues.push_back(delta * delta);
                        quadraticMean += delta * delta;
                    }
                    quadraticMean /= static_cast<double>(count);
                    for (double& value : quadraticValues) {
                        value -= quadraticMean;
                    }
                    for (const std::vector<double>& linearMode : localLinearModes) {
                        double projection = 0.0;
                        for (std::size_t i = 0; i < count; ++i) {
                            projection += quadraticValues[i] * linearMode[i];
                        }
                        for (std::size_t i = 0; i < count; ++i) {
                            quadraticValues[i] -= projection * linearMode[i];
                        }
                    }

                    double quadraticNormSquared = 0.0;
                    for (double value : quadraticValues) {
                        quadraticNormSquared += value * value;
                    }
                    if (quadraticNormSquared > 1.0e-24) {
                        for (std::size_t i = 0; i < count; ++i) {
                            globalQuadraticZ[static_cast<std::size_t>(interfaceIndices[i])] +=
                                quadraticValues[i];
                        }
                    }
                }
            };

            appendSubdomainModes(domain.interfaceGlobalDofs);
        }

        double globalQuadraticNormSquared = 0.0;
        for (double value : globalQuadraticZ) {
            globalQuadraticNormSquared += value * value;
        }
        if (addGlobalQuadraticZ && globalQuadraticNormSquared > 1.0e-24) {
            const double inverseNorm = 1.0 / std::sqrt(globalQuadraticNormSquared);
            CoarseVector quadratic;
            for (std::size_t i = 0; i < globalQuadraticZ.size(); ++i) {
                if (globalQuadraticZ[i] != 0.0) {
                    quadratic.indices.push_back(static_cast<int>(i));
                    quadratic.values.push_back(globalQuadraticZ[i] * inverseNorm);
                }
            }
            coarseBasis.push_back(std::move(quadratic));
        }

        if (coarseEnergyAdaptive) {
            // Adaptive eigenvectors enrich, rather than replace, the fixed
            // near-kernel space.  buildEnergyCoarseBasis orthogonalizes every
            // selected local mode against these existing volume modes.
            buildEnergyCoarseBasis(energyMaxModesPerDomain,
                                   energySubspaceIterations,
                                   energyEigenvalueThreshold);
        }
        const int coarseDim = static_cast<int>(coarseBasis.size());
        if (coarseDim == 0) {
            return;
        }

        std::vector<Entry> coarseEntries;
        std::vector<std::vector<std::pair<int, double>>> basisByDof(
            partition.interfaceGlobalDofs.size());
        for (int coarse = 0; coarse < coarseDim; ++coarse) {
            const CoarseVector& basis = coarseBasis[static_cast<std::size_t>(coarse)];
            for (std::size_t i = 0; i < basis.indices.size(); ++i) {
                basisByDof[static_cast<std::size_t>(basis.indices[i])].push_back(
                    {coarse, basis.values[i]});
            }
        }
        std::vector<double> coarseColumn(static_cast<std::size_t>(coarseDim), 0.0);
        std::vector<int> touchedRows;
        touchedRows.reserve(static_cast<std::size_t>(coarseDim));
        std::vector<unsigned char> rowTouched(static_cast<std::size_t>(coarseDim), 0);
        coarseImages.reserve(static_cast<std::size_t>(coarseDim));
        for (int col = 0; col < coarseDim; ++col) {
            std::vector<double> basis(partition.interfaceGlobalDofs.size(), 0.0);
            const CoarseVector& sparseBasis = coarseBasis[static_cast<std::size_t>(col)];
            for (std::size_t i = 0; i < sparseBasis.indices.size(); ++i) {
                basis[static_cast<std::size_t>(sparseBasis.indices[i])] = sparseBasis.values[i];
            }
            std::vector<double> image;
            applyRaw(basis, image, false);
            CoarseVector sparseImage;
            for (std::size_t i = 0; i < image.size(); ++i) {
                if (image[i] != 0.0) {
                    sparseImage.indices.push_back(static_cast<int>(i));
                    sparseImage.values.push_back(image[i]);
                }
            }
            for (std::size_t i = 0; i < sparseImage.indices.size(); ++i) {
                const int dof = sparseImage.indices[i];
                const double imageValue = sparseImage.values[i];
                for (const auto& incidence : basisByDof[static_cast<std::size_t>(dof)]) {
                    const int row = incidence.first;
                    if (rowTouched[static_cast<std::size_t>(row)] == 0) {
                        rowTouched[static_cast<std::size_t>(row)] = 1;
                        touchedRows.push_back(row);
                    }
                    coarseColumn[static_cast<std::size_t>(row)] +=
                        incidence.second * imageValue;
                }
            }
            std::sort(touchedRows.begin(), touchedRows.end());
            for (int row : touchedRows) {
                const double value = coarseColumn[static_cast<std::size_t>(row)];
                if (value != 0.0) {
                    coarseEntries.push_back({row, col, value});
                }
                coarseColumn[static_cast<std::size_t>(row)] = 0.0;
                rowTouched[static_cast<std::size_t>(row)] = 0;
            }
            touchedRows.clear();
            coarseImages.push_back(std::move(sparseImage));
        }
        coarseSolver = LocalSolver(coarseDim, coarseEntries);
    }

    void rebuildCoarseFactor()
    {
        const int coarseDim = static_cast<int>(coarseBasis.size());
        if (coarseDim == 0) {
            return;
        }

        std::vector<Entry> coarseEntries;
        std::vector<std::vector<std::pair<int, double>>> basisByDof(
            partition.interfaceGlobalDofs.size());
        for (int coarse = 0; coarse < coarseDim; ++coarse) {
            const CoarseVector& basis = coarseBasis[static_cast<std::size_t>(coarse)];
            for (std::size_t i = 0; i < basis.indices.size(); ++i) {
                basisByDof[static_cast<std::size_t>(basis.indices[i])].push_back(
                    {coarse, basis.values[i]});
            }
        }
        std::vector<double> coarseColumn(static_cast<std::size_t>(coarseDim), 0.0);
        std::vector<int> touchedRows;
        touchedRows.reserve(static_cast<std::size_t>(coarseDim));
        std::vector<unsigned char> rowTouched(static_cast<std::size_t>(coarseDim), 0);
        coarseImages.clear();
        coarseImages.reserve(static_cast<std::size_t>(coarseDim));
        for (int col = 0; col < coarseDim; ++col) {
            std::vector<double> basis(partition.interfaceGlobalDofs.size(), 0.0);
            const CoarseVector& sparseBasis = coarseBasis[static_cast<std::size_t>(col)];
            for (std::size_t i = 0; i < sparseBasis.indices.size(); ++i) {
                basis[static_cast<std::size_t>(sparseBasis.indices[i])] = sparseBasis.values[i];
            }
            std::vector<double> image;
            applyRaw(basis, image, false);
            CoarseVector sparseImage;
            for (std::size_t i = 0; i < image.size(); ++i) {
                if (image[i] != 0.0) {
                    sparseImage.indices.push_back(static_cast<int>(i));
                    sparseImage.values.push_back(image[i]);
                }
            }
            for (std::size_t i = 0; i < sparseImage.indices.size(); ++i) {
                const int dof = sparseImage.indices[i];
                const double imageValue = sparseImage.values[i];
                for (const auto& incidence : basisByDof[static_cast<std::size_t>(dof)]) {
                    const int row = incidence.first;
                    if (rowTouched[static_cast<std::size_t>(row)] == 0) {
                        rowTouched[static_cast<std::size_t>(row)] = 1;
                        touchedRows.push_back(row);
                    }
                    coarseColumn[static_cast<std::size_t>(row)] +=
                        incidence.second * imageValue;
                }
            }
            std::sort(touchedRows.begin(), touchedRows.end());
            for (int row : touchedRows) {
                const double value = coarseColumn[static_cast<std::size_t>(row)];
                if (value != 0.0) {
                    coarseEntries.push_back({row, col, value});
                }
                coarseColumn[static_cast<std::size_t>(row)] = 0.0;
                rowTouched[static_cast<std::size_t>(row)] = 0;
            }
            touchedRows.clear();
            coarseImages.push_back(std::move(sparseImage));
        }
        coarseSolver = LocalSolver(coarseDim, coarseEntries);
    }

    bool prepareGlobalSlow(const std::vector<double>& arnoldiSeed)
    {
        if (!globalSlowPending) {
            return false;
        }
        buildGlobalSlowCoarseBasis(pendingGlobalSlowModes,
                                   pendingGlobalSlowSubspaceDimension,
                                   arnoldiSeed);
        rebuildCoarseFactor();
        refreshMemoryBytes();
        globalSlowPending = false;
        return true;
    }

    void refreshMemoryBytes()
    {
        bytes = partitionMemoryBytes(partition)
            + interfaceMatrix.rowPtr.size() * sizeof(int)
            + interfaceMatrix.colInd.size() * sizeof(int)
            + interfaceMatrix.values.size() * sizeof(double);
        for (const LocalSolver& solver : interiorSolvers) {
            bytes += solver.memoryBytes();
        }
        for (const LocalSolver& solver : blockSolvers) {
            bytes += solver.memoryBytes();
        }
        bytes += coarseSolver.memoryBytes()
            + coarseBasis.size() * sizeof(CoarseVector);
        if (proxySolver) {
            bytes += proxySolver->memoryBytes();
        }
        for (const CoarseVector& basis : coarseBasis) {
            bytes += basis.indices.size() * sizeof(int) + basis.values.size() * sizeof(double);
        }
        for (const CoarseVector& image : coarseImages) {
            bytes += image.indices.size() * sizeof(int) + image.values.size() * sizeof(double);
        }
    }

    void resetRuntimeCounters()
    {
        solveCalls = 0;
        matvecs = 0;
        localSeconds = 0.0;
        schurApplySeconds = 0.0;
        coarseSeconds = 0.0;
        for (SubdomainPerformance& performance : domainPerformance) {
            performance.phase33Calls = 0;
            performance.phase33Seconds = 0.0;
        }
        if (proxySolver) {
            proxySolver->resetRuntimeCounters();
        }
    }

    void applyCoarse(const std::vector<double>& residual,
                     std::vector<double>& correction,
                     std::vector<double>* image)
    {
        correction.assign(residual.size(), 0.0);
        if (image != nullptr) {
            image->assign(residual.size(), 0.0);
        }
        if (coarseBasis.empty()) {
            return;
        }
        std::vector<double> coarseRhs(coarseBasis.size(), 0.0);
        for (std::size_t coarse = 0; coarse < coarseBasis.size(); ++coarse) {
            for (std::size_t i = 0; i < coarseBasis[coarse].indices.size(); ++i) {
                coarseRhs[coarse] += coarseBasis[coarse].values[i]
                    * residual[static_cast<std::size_t>(coarseBasis[coarse].indices[i])];
            }
        }
        std::vector<double> coefficients;
        solveCoarse(coarseRhs, coefficients);
        for (std::size_t coarse = 0; coarse < coarseBasis.size(); ++coarse) {
            for (std::size_t i = 0; i < coarseBasis[coarse].indices.size(); ++i) {
                correction[static_cast<std::size_t>(coarseBasis[coarse].indices[i])] +=
                    coarseBasis[coarse].values[i] * coefficients[coarse];
            }
            if (image != nullptr) {
                const CoarseVector& coarseImage = coarseImages[coarse];
                for (std::size_t i = 0; i < coarseImage.indices.size(); ++i) {
                    (*image)[static_cast<std::size_t>(coarseImage.indices[i])] +=
                        coarseImage.values[i] * coefficients[coarse];
                }
            }
        }
    }

    void applyLeftCoarseProjection(const std::vector<double>& vector,
                                   std::vector<double>& correction)
    {
        correction.assign(vector.size(), 0.0);
        if (coarseBasis.empty()) {
            return;
        }

        // coarseImages[j] = S * Z_j.  The assembled Schur complement is
        // symmetric, so (S * Z_j)^T * vector = Z_j^T * S * vector.  Reusing
        // these setup-time images applies Q*S without another Schur matvec.
        std::vector<double> coarseRhs(coarseBasis.size(), 0.0);
        for (std::size_t coarse = 0; coarse < coarseBasis.size(); ++coarse) {
            const CoarseVector& image = coarseImages[coarse];
            for (std::size_t i = 0; i < image.indices.size(); ++i) {
                coarseRhs[coarse] += image.values[i]
                    * vector[static_cast<std::size_t>(image.indices[i])];
            }
        }

        std::vector<double> coefficients;
        solveCoarse(coarseRhs, coefficients);
        for (std::size_t coarse = 0; coarse < coarseBasis.size(); ++coarse) {
            for (std::size_t i = 0; i < coarseBasis[coarse].indices.size(); ++i) {
                correction[static_cast<std::size_t>(coarseBasis[coarse].indices[i])] +=
                    coarseBasis[coarse].values[i] * coefficients[coarse];
            }
        }
    }
};

SchurOperator::SchurOperator(const Mesh& mesh,
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
                             const std::string& proxyOutputDirectory)
    : impl_(std::make_unique<Impl>(mesh,
                                   system,
                                   physics,
                                   coarseLinearXY,
                                   coarseLinearZ,
                                   coarseGlobalQuadraticZ,
                                   coarseInterfacePatches,
                                   coarseInterfaceLinearXY,
                                   coarseEnergyAdaptive,
                                   energyMaxModesPerDomain,
                                   energySubspaceIterations,
                                   energyEigenvalueThreshold,
                                   coarseGlobalSlow,
                                   globalSlowModes,
                                   globalSlowSubspaceDimension,
                                   proxyEnabled,
                                   proxyDisableCoarse,
                                   proxyDiagnostics,
                                   proxyHighConductivityThreshold,
                                   proxyUseMaterialConnectivity,
                                   proxyRing,
                                   proxyProbeColumns,
                                   proxyBlockSize,
                                   proxyValidateBlockEquivalence,
                                   localSolveThreads,
                                   localPardisoThreads,
                                   proxyCacheEnabled,
                                   proxyCachePath,
                                   proxyOutputDirectory)) {}

SchurOperator::~SchurOperator() = default;
SchurOperator::SchurOperator(SchurOperator&&) noexcept = default;
SchurOperator& SchurOperator::operator=(SchurOperator&&) noexcept = default;

int SchurOperator::domains() const { return static_cast<int>(impl_->partition.domains.size()); }
int SchurOperator::totalDofs() const { return impl_->partition.totalDofs; }
int SchurOperator::interfaceDofs() const { return static_cast<int>(impl_->partition.interfaceGlobalDofs.size()); }
int SchurOperator::interiorDofs() const { return totalDofs() - interfaceDofs(); }
int SchurOperator::coarseDimension() const { return static_cast<int>(impl_->coarseBasis.size()); }
int SchurOperator::interfacePatchCount() const { return impl_->patchCount; }
int SchurOperator::energyCandidateModes() const { return impl_->energyCandidates; }
double SchurOperator::energyEigenSetupSeconds() const { return impl_->energySeconds; }
double SchurOperator::energySelectedEigenvalueMin() const { return impl_->energyEigenMin; }
double SchurOperator::energySelectedEigenvalueMax() const { return impl_->energyEigenMax; }
int SchurOperator::globalSlowCandidateDimension() const { return impl_->globalSlowCandidates; }
double SchurOperator::globalSlowSetupSeconds() const { return impl_->globalSlowSeconds; }
double SchurOperator::globalSlowEstimatedLambdaMax() const { return impl_->globalSlowLambdaMax; }
double SchurOperator::globalSlowSelectedRitzMin() const { return impl_->globalSlowRitzMin; }
double SchurOperator::globalSlowSelectedRitzMax() const { return impl_->globalSlowRitzMax; }
int SchurOperator::proxyGraphEdges() const
{
    return static_cast<int>(std::min<std::size_t>(
        impl_->proxyDiagnostics.graphEdges,
        static_cast<std::size_t>(std::numeric_limits<int>::max())));
}
int SchurOperator::proxyProbeColumns() const { return impl_->proxyDiagnostics.probeColumns; }
int SchurOperator::proxySchurApplies() const { return impl_->proxyDiagnostics.exactSchurApplies; }
double SchurOperator::proxyDiagnosticsSeconds() const { return impl_->proxyDiagnostics.setupSeconds; }
double SchurOperator::proxyRing3Coverage() const
{
    return impl_->proxyDiagnostics.ringMetrics.empty()
        ? 0.0 : impl_->proxyDiagnostics.ringMetrics.back().energyCoverage;
}
double SchurOperator::proxyRing3OperatorError() const
{
    return impl_->proxyDiagnostics.ringMetrics.empty()
        ? 0.0 : impl_->proxyDiagnostics.ringMetrics.back().randomOperatorErrorMaximum;
}
std::size_t SchurOperator::proxyRing3EstimatedNnz() const
{
    return impl_->proxyDiagnostics.ringMetrics.empty()
        ? 0 : impl_->proxyDiagnostics.ringMetrics.back().estimatedNnz;
}
std::size_t SchurOperator::proxyRing3MemoryEstimateBytes() const
{
    return impl_->proxyDiagnostics.ringMetrics.empty()
        ? 0 : impl_->proxyDiagnostics.ringMetrics.back().csrMemoryEstimateBytes;
}
bool SchurOperator::proxyRecommended() const { return impl_->proxyDiagnostics.proxyRecommended; }
std::size_t SchurOperator::proxyNnz() const
{
    return impl_->proxySolver ? impl_->proxySolver->nnz() : 0;
}
double SchurOperator::proxyDensity() const
{
    return impl_->proxySolver ? impl_->proxySolver->density() : 0.0;
}
int SchurOperator::proxyColors() const
{
    return impl_->proxySolver ? impl_->proxySolver->colors() : 0;
}
int SchurOperator::proxyProbingSchurApplies() const
{
    return impl_->proxySolver ? impl_->proxySolver->probingSchurApplies() : 0;
}
int SchurOperator::proxyProbingBlockSize() const
{
    return impl_->proxySolver ? impl_->proxySolver->probingBlockSize() : 0;
}
int SchurOperator::proxyProbingBlockCalls() const
{
    return impl_->proxySolver ? impl_->proxySolver->probingBlockCalls() : 0;
}
int SchurOperator::proxyValidationSchurApplies() const
{
    return impl_->proxySolver ? impl_->proxySolver->validationSchurApplies() : 0;
}
int SchurOperator::proxySymbolicCalls() const
{
    return impl_->proxySolver ? impl_->proxySolver->symbolicCalls() : 0;
}
int SchurOperator::proxyNumericalCalls() const
{
    return impl_->proxySolver ? impl_->proxySolver->numericalCalls() : 0;
}
int SchurOperator::proxySolveCalls() const
{
    return impl_->proxySolver ? impl_->proxySolver->solveCalls() : 0;
}
double SchurOperator::proxySetupSeconds() const
{
    return impl_->proxySolver ? impl_->proxySolver->setupSeconds() : 0.0;
}
double SchurOperator::proxySymbolicSeconds() const
{
    return impl_->proxySolver ? impl_->proxySolver->symbolicSeconds() : 0.0;
}
double SchurOperator::proxyNumericalSeconds() const
{
    return impl_->proxySolver ? impl_->proxySolver->numericalSeconds() : 0.0;
}
double SchurOperator::proxySolveSeconds() const
{
    return impl_->proxySolver ? impl_->proxySolver->solveSeconds() : 0.0;
}
double SchurOperator::proxySymmetryError() const
{
    return impl_->proxySolver ? impl_->proxySolver->symmetryError() : 0.0;
}
double SchurOperator::proxyMinimumTestRayleigh() const
{
    return impl_->proxySolver ? impl_->proxySolver->minimumTestRayleigh() : 0.0;
}
double SchurOperator::proxyDiagonalShift() const
{
    return impl_->proxySolver ? impl_->proxySolver->diagonalShift() : 0.0;
}
double SchurOperator::proxyDiagonalCompensation() const
{
    return impl_->proxySolver ? impl_->proxySolver->diagonalCompensation() : 0.0;
}
std::uint64_t SchurOperator::proxyValueHash() const
{
    return impl_->proxySolver ? impl_->proxySolver->valueHash() : 0;
}
double SchurOperator::proxyBlockMaximumDifference() const
{
    return impl_->proxySolver ? impl_->proxySolver->blockMaximumDifference() : 0.0;
}
double SchurOperator::proxyBlockRelativeDifference() const
{
    return impl_->proxySolver ? impl_->proxySolver->blockRelativeDifference() : 0.0;
}
std::size_t SchurOperator::proxyMemoryBytes() const
{
    return impl_->proxySolver ? impl_->proxySolver->memoryBytes() : 0;
}
bool SchurOperator::proxyMatrixCacheHit() const
{
    return impl_->proxySolver && impl_->proxySolver->matrixCacheHit();
}
bool SchurOperator::proxyFactorCacheHit() const
{
    return impl_->proxySolver && impl_->proxySolver->factorCacheHit();
}
double SchurOperator::localFactorizationSeconds() const { return impl_->factorSeconds; }
double SchurOperator::localSymbolicAnalysisSeconds() const { return impl_->symbolicSeconds; }
double SchurOperator::localNumericalFactorizationSeconds() const { return impl_->numericalSeconds; }
double SchurOperator::localSolveSeconds() const { return impl_->localSeconds; }
double SchurOperator::schurApplySeconds() const { return impl_->schurApplySeconds; }
double SchurOperator::coarseSolveSeconds() const { return impl_->coarseSeconds; }
std::size_t SchurOperator::memoryBytes() const { return impl_->bytes; }
int SchurOperator::localSolveCalls() const { return impl_->solveCalls; }
int SchurOperator::localSymbolicAnalysisCalls() const { return impl_->symbolicCalls; }
int SchurOperator::localNumericalFactorizationCalls() const { return impl_->numericalCalls; }
int SchurOperator::matvecCalls() const { return impl_->matvecs; }
std::vector<SubdomainPerformance> SchurOperator::subdomainPerformance() const
{
    return impl_->domainPerformance;
}
const std::vector<int>& SchurOperator::interfaceGlobalDofs() const
{
    return impl_->partition.interfaceGlobalDofs;
}

std::vector<double> SchurOperator::condensedRhs(const std::vector<double>& globalRhs)
{
    if (static_cast<int>(globalRhs.size()) != totalDofs()) {
        throw std::runtime_error("[Schur] Global right-hand side has the wrong size.");
    }
    std::vector<double> result(static_cast<std::size_t>(interfaceDofs()), 0.0);
    for (int gamma = 0; gamma < interfaceDofs(); ++gamma) {
        result[static_cast<std::size_t>(gamma)] =
            globalRhs[static_cast<std::size_t>(impl_->partition.interfaceGlobalDofs[static_cast<std::size_t>(gamma)])];
    }
    const int domainCount = static_cast<int>(impl_->partition.domains.size());
    std::vector<std::vector<double>> corrections(static_cast<std::size_t>(domainCount));
    std::vector<double> solveTimes(static_cast<std::size_t>(domainCount), 0.0);
    std::vector<std::exception_ptr> errors(static_cast<std::size_t>(domainCount));
#pragma omp parallel for num_threads(impl_->localSolveThreads) if(impl_->localSolveThreads > 1) schedule(static)
    for (int slotIndex = 0; slotIndex < domainCount; ++slotIndex) {
        try {
            const std::size_t slot = static_cast<std::size_t>(slotIndex);
            const DomainBlocks& domain = impl_->partition.domains[slot];
            std::vector<double> localRhs(domain.interiorGlobalDofs.size(), 0.0);
            for (std::size_t i = 0; i < domain.interiorGlobalDofs.size(); ++i) {
                localRhs[i] = globalRhs[static_cast<std::size_t>(domain.interiorGlobalDofs[i])];
            }
            std::vector<double> eliminated;
            solveTimes[slot] = impl_->solveDomain(
                slot, impl_->interiorSolvers[slot], localRhs, eliminated);
            corrections[slot].assign(domain.interfaceGlobalDofs.size(), 0.0);
            for (std::size_t localGamma = 0;
                 localGamma < domain.interfaceGlobalDofs.size(); ++localGamma) {
                for (const auto& entry : domain.interfaceInteriorRows[localGamma]) {
                    corrections[slot][localGamma] += entry.second
                        * eliminated[static_cast<std::size_t>(entry.first)];
                }
            }
        } catch (...) {
            errors[static_cast<std::size_t>(slotIndex)] = std::current_exception();
        }
    }
    for (const std::exception_ptr& error : errors) {
        if (error) {
            std::rethrow_exception(error);
        }
    }
    impl_->localSeconds += std::accumulate(solveTimes.begin(), solveTimes.end(), 0.0);
    impl_->solveCalls += domainCount;
    for (std::size_t slot = 0; slot < impl_->partition.domains.size(); ++slot) {
        const DomainBlocks& domain = impl_->partition.domains[slot];
        for (std::size_t localGamma = 0;
             localGamma < domain.interfaceGlobalDofs.size(); ++localGamma) {
            const int gamma = impl_->partition.globalToInterface[static_cast<std::size_t>(
                domain.interfaceGlobalDofs[localGamma])];
            result[static_cast<std::size_t>(gamma)] -= corrections[slot][localGamma];
        }
    }
    return result;
}

bool SchurOperator::prepareGlobalSlow(const std::vector<double>& arnoldiSeed)
{
    return impl_->prepareGlobalSlow(arnoldiSeed);
}

void SchurOperator::resetRuntimeCounters()
{
    impl_->resetRuntimeCounters();
}

void SchurOperator::apply(const std::vector<double>& interfaceVector, std::vector<double>& result)
{
    if (static_cast<int>(interfaceVector.size()) != interfaceDofs()) {
        throw std::runtime_error("[Schur] Interface vector has the wrong size.");
    }
    impl_->applyRaw(interfaceVector, result, true);
}

void SchurOperator::applyBlockPreconditioner(const std::vector<double>& residual,
                                             std::vector<double>& result)
{
    if (static_cast<int>(residual.size()) != interfaceDofs()) {
        throw std::runtime_error("[Schur] Preconditioner vector has the wrong size.");
    }
    std::vector<double> coarseCorrection;
    std::vector<double> coarseImage;
    impl_->applyCoarse(residual, coarseCorrection, &coarseImage);
    std::vector<double> projectedResidual(residual.size(), 0.0);
    for (std::size_t i = 0; i < residual.size(); ++i) {
        projectedResidual[i] = residual[i] - coarseImage[i];
    }

    std::vector<double> localCorrection;
    impl_->applyLevelOne(projectedResidual, localCorrection);

    std::vector<double> leftCoarseProjection;
    impl_->applyLeftCoarseProjection(localCorrection, leftCoarseProjection);
    result.resize(residual.size());
    for (std::size_t i = 0; i < residual.size(); ++i) {
        result[i] = coarseCorrection[i] + localCorrection[i] - leftCoarseProjection[i];
    }
}

std::vector<double> SchurOperator::recover(const std::vector<double>& globalRhs,
                                           const std::vector<double>& interfaceSolution)
{
    if (static_cast<int>(globalRhs.size()) != totalDofs()
        || static_cast<int>(interfaceSolution.size()) != interfaceDofs()) {
        throw std::runtime_error("[Schur] Recovery vector has the wrong size.");
    }
    std::vector<double> solution(static_cast<std::size_t>(totalDofs()), 0.0);
    for (int gamma = 0; gamma < interfaceDofs(); ++gamma) {
        solution[static_cast<std::size_t>(impl_->partition.interfaceGlobalDofs[static_cast<std::size_t>(gamma)])] =
            interfaceSolution[static_cast<std::size_t>(gamma)];
    }
    const int domainCount = static_cast<int>(impl_->partition.domains.size());
    std::vector<std::vector<double>> localSolutions(static_cast<std::size_t>(domainCount));
    std::vector<double> solveTimes(static_cast<std::size_t>(domainCount), 0.0);
    std::vector<std::exception_ptr> errors(static_cast<std::size_t>(domainCount));
#pragma omp parallel for num_threads(impl_->localSolveThreads) if(impl_->localSolveThreads > 1) schedule(static)
    for (int slotIndex = 0; slotIndex < domainCount; ++slotIndex) {
        try {
            const std::size_t slot = static_cast<std::size_t>(slotIndex);
            const DomainBlocks& domain = impl_->partition.domains[slot];
            std::vector<double> localRhs(domain.interiorGlobalDofs.size(), 0.0);
            for (std::size_t row = 0; row < domain.interiorGlobalDofs.size(); ++row) {
                localRhs[row] = globalRhs[static_cast<std::size_t>(
                    domain.interiorGlobalDofs[row])];
                for (const auto& entry : domain.interiorInterfaceRows[row]) {
                    localRhs[row] -= entry.second
                        * interfaceSolution[static_cast<std::size_t>(entry.first)];
                }
            }
            solveTimes[slot] = impl_->solveDomain(
                slot, impl_->interiorSolvers[slot], localRhs, localSolutions[slot]);
        } catch (...) {
            errors[static_cast<std::size_t>(slotIndex)] = std::current_exception();
        }
    }
    for (const std::exception_ptr& error : errors) {
        if (error) {
            std::rethrow_exception(error);
        }
    }
    impl_->localSeconds += std::accumulate(solveTimes.begin(), solveTimes.end(), 0.0);
    impl_->solveCalls += domainCount;
    for (std::size_t slot = 0; slot < impl_->partition.domains.size(); ++slot) {
        const DomainBlocks& domain = impl_->partition.domains[slot];
        for (std::size_t i = 0; i < domain.interiorGlobalDofs.size(); ++i) {
            solution[static_cast<std::size_t>(domain.interiorGlobalDofs[i])] =
                localSolutions[slot][i];
        }
    }
    return solution;
}

} // namespace ddm_schur
