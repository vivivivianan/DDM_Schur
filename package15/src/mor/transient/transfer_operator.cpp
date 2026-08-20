#include "transfer_operator.hpp"

#include "ddm_schur/interface_operator.hpp"
#include "sipg_core.hpp"
#include "linear_solvers.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace mor::transient {
namespace {

using Clock = std::chrono::steady_clock;

double elapsed(const Clock::time_point& start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

double dot(const std::vector<double>& left,
           const std::vector<double>& right)
{
    long double value = 0.0L;
    for (std::size_t row = 0; row < left.size(); ++row) {
        value += static_cast<long double>(left[row]) * right[row];
    }
    return static_cast<double>(value);
}

double dot(const double* left, const double* right, int rows)
{
    long double value = 0.0L;
    for (int row = 0; row < rows; ++row) {
        value += static_cast<long double>(left[row]) * right[row];
    }
    return static_cast<double>(value);
}

double norm(const std::vector<double>& values)
{
    return std::sqrt(std::max(0.0, dot(values, values)));
}

double stableNorm(const std::vector<double>& values)
{
#ifdef USE_MKL_PARDISO
    // MKL's xNRM2 uses a scaled sum of squares and avoids both overflow and
    // the serial cost of a hand-written norm in the iterative baseline.
    return cblas_dnrm2(
        static_cast<MKL_INT>(values.size()), values.data(), 1);
#else
    // LAPACK-style scaled sum of squares. Interface penalty units can make
    // a perfectly finite vector overflow in a naive dot(v,v).
    double scale = 0.0;
    double sum = 1.0;
    for (double value : values) {
        const double magnitude = std::abs(value);
        if (magnitude == 0.0) continue;
        if (!std::isfinite(magnitude)) {
            return std::numeric_limits<double>::infinity();
        }
        if (scale < magnitude) {
            const double ratio = scale / magnitude;
            sum = 1.0 + sum * ratio * ratio;
            scale = magnitude;
        } else {
            const double ratio = magnitude / scale;
            sum += ratio * ratio;
        }
    }
    return scale == 0.0 ? 0.0 : scale * std::sqrt(sum);
#endif
}

double relativeDifference(const std::vector<double>& candidate,
                          const std::vector<double>& reference)
{
    if (candidate.size() != reference.size()) {
        return std::numeric_limits<double>::infinity();
    }
    std::vector<double> difference = candidate;
    for (std::size_t row = 0; row < difference.size(); ++row) {
        difference[row] -= reference[row];
    }
    return norm(difference)
        / std::max(1.0e-300, norm(reference));
}

template <typename T>
void hashValue(std::uint64_t& hash, const T& value)
{
    const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
    for (std::size_t byte = 0; byte < sizeof(T); ++byte) {
        hash ^= bytes[byte];
        hash *= UINT64_C(1099511628211);
    }
}

std::uint64_t fingerprint(const std::vector<int>& values)
{
    std::uint64_t result = UINT64_C(1469598103934665603);
    hashValue(result, values.size());
    for (int value : values) hashValue(result, value);
    return result;
}

std::uint64_t fingerprint(const std::vector<double>& values)
{
    std::uint64_t result = UINT64_C(1469598103934665603);
    hashValue(result, values.size());
    for (double value : values) hashValue(result, value);
    return result;
}

std::vector<double> condensedUnitInputs(
    const local::Model& model,
    const std::vector<double>& input,
    int sourceChannels)
{
    if (sourceChannels < 0
        || input.size() != static_cast<std::size_t>(
            model.globalDofs * sourceChannels)) {
        throw std::runtime_error(
            "[Optimal port] Unit-input matrix dimensions are invalid.");
    }
    std::vector<double> result(static_cast<std::size_t>(
        model.interfaceDofs * sourceChannels), 0.0);
    for (int channel = 0; channel < sourceChannels; ++channel) {
        for (int gamma = 0; gamma < model.interfaceDofs; ++gamma) {
            const int global =
                model.interfaceGlobalDofs[static_cast<std::size_t>(gamma)];
            result[static_cast<std::size_t>(
                channel * model.interfaceDofs + gamma)] =
                input[static_cast<std::size_t>(
                    channel * model.globalDofs + global)];
        }
    }
    for (const auto& subdomain : model.subdomains) {
        const local::DenseSymmetricFactor factor =
            local::factorDenseSymmetric(
                subdomain.reducedInterior, subdomain.rank);
        for (int channel = 0; channel < sourceChannels; ++channel) {
            std::vector<double> reduced(
                static_cast<std::size_t>(subdomain.rank), 0.0);
            for (int mode = 0; mode < subdomain.rank; ++mode) {
                for (int row = 0; row < subdomain.interiorDofs; ++row) {
                    const int global = subdomain.interiorGlobalDofs[
                        static_cast<std::size_t>(row)];
                    reduced[static_cast<std::size_t>(mode)] +=
                        subdomain.basis[static_cast<std::size_t>(
                            mode * subdomain.interiorDofs + row)]
                        * input[static_cast<std::size_t>(
                            channel * model.globalDofs + global)];
                }
            }
            local::solveDenseSymmetric(factor, reduced);
            for (int gamma = 0;
                 gamma < subdomain.localInterfaceDofs; ++gamma) {
                double correction = 0.0;
                for (int mode = 0; mode < subdomain.rank; ++mode) {
                    correction += subdomain.reducedInterfaceInterior[
                        static_cast<std::size_t>(
                            gamma * subdomain.rank + mode)]
                        * reduced[static_cast<std::size_t>(mode)];
                }
                const int globalGamma = subdomain.interfaceIndices[
                    static_cast<std::size_t>(gamma)];
                result[static_cast<std::size_t>(
                    channel * model.interfaceDofs + globalGamma)] -=
                    correction;
            }
        }
    }
    return result;
}

} // namespace

GeneralizedTransferSourceBlocks buildGeneralizedTransferSourceBlocks(
    const local::Model& model,
    const std::vector<double>& input,
    int sourceChannels,
    const std::vector<double>& boundaryLoad,
    const std::vector<double>& condensedHistory,
    int historyChannels)
{
    if (historyChannels < 0
        || condensedHistory.size() != static_cast<std::size_t>(
            model.interfaceDofs * historyChannels)) {
        throw std::runtime_error(
            "[Optimal port] Condensed history dimensions are invalid.");
    }
    const int boundaryChannels = boundaryLoad.empty() ? 0 : 1;
    GeneralizedTransferSourceBlocks blocks;
    blocks.interfaceDofs = model.interfaceDofs;
    blocks.inputChannels = sourceChannels;
    blocks.boundaryChannels = boundaryChannels;
    blocks.historyChannels = historyChannels;
    blocks.input = condensedUnitInputs(
        model, input, sourceChannels);
    blocks.boundary = condensedUnitInputs(
        model, boundaryLoad, boundaryChannels);
    blocks.history = condensedHistory;
    blocks.inputFingerprint = fingerprint(input);
    blocks.boundaryFingerprint = fingerprint(boundaryLoad);
    blocks.historyFingerprint = fingerprint(condensedHistory);
    return blocks;
}

struct PatchTransferOperator::PatchAlgebra {
    struct Entry {
        int row = 0;
        int column = 0;
        double value = 0.0;
    };
    struct LocalCorrection {
        std::size_t subdomainSlot = 0;
        int rank = 0;
        int correctionOffset = 0;
        std::vector<int> targetLocalDofs;
        std::vector<int> targetRows;
        std::vector<int> sourceLocalDofs;
        std::vector<int> sourceRows;
    };

    std::vector<Entry> targetEntries;
    std::vector<Entry> crossEntries;
    std::vector<MatrixEntry> sparseTargetEntries;
    std::vector<LocalCorrection> localCorrections;
};

struct PatchTransferOperator::SparseTargetFactor {
    SparseTargetFactor(int dimension,
                       const std::vector<MatrixEntry>& entries,
                       bool symmetric,
                       std::string& fallbackReason)
    {
        if (symmetric) {
            try {
                spd = std::make_unique<SubdomainDirectSolver>(
                    dimension, entries);
                return;
            } catch (const std::exception& error) {
                fallbackReason = std::string("pardiso_spd_failed:")
                    + error.what();
            }
        }
        general = std::make_unique<GeneralSparseDirectSolver>(
            dimension, entries);
    }

    bool symmetric() const
    {
        return static_cast<bool>(spd);
    }

    void solve(const std::vector<double>& rightHandSide,
               bool transpose,
               std::vector<double>& solution) const
    {
        if (spd) {
            spd->solve(rightHandSide, solution);
        } else if (transpose) {
            general->solveTranspose(rightHandSide, solution);
        } else {
            general->solve(rightHandSide, solution);
        }
    }

    void solveMultiple(const std::vector<double>& rightHandSides,
                       int rightHandSideCount,
                       bool transpose,
                       std::vector<double>& solutions) const
    {
        if (spd) {
            spd->solveMultiple(
                rightHandSides, rightHandSideCount, solutions);
        } else {
            general->solveMultiple(
                rightHandSides, rightHandSideCount, transpose, solutions);
        }
    }

    double factorizationSeconds() const
    {
        if (spd) {
            return spd->symbolicAnalysisSeconds()
                + spd->numericalFactorizationSeconds();
        }
        return general->symbolicAnalysisSeconds()
            + general->numericalFactorizationSeconds();
    }

    std::size_t memoryBytes() const
    {
        return spd ? spd->memoryBytes() : general->memoryBytes();
    }

    std::unique_ptr<SubdomainDirectSolver> spd;
    std::unique_ptr<GeneralSparseDirectSolver> general;
};

std::vector<PortPatch> buildOptimalPortPatches(
    const Mesh& mesh,
    const ddm_schur::InterfacePartition& partition,
    int oversamplingLayers)
{
    if (oversamplingLayers < 0 || oversamplingLayers > 2) {
        throw std::runtime_error(
            "[Optimal port] Oversampling layers must be 0, 1, or 2.");
    }
    std::map<std::pair<int, int>, int> interfaceIds;
    for (std::size_t index = 0; index < mesh.interfaceSummaries.size(); ++index) {
        const auto& summary = mesh.interfaceSummaries[index];
        interfaceIds[std::minmax(
            summary.leftSubdomain, summary.rightSubdomain)] =
            static_cast<int>(index);
    }
    std::map<std::pair<int, int>, std::set<int>> candidateTargets;
    std::map<int, const ddm_schur::DomainBlocks*> domains;
    std::map<int, std::set<int>> neighbors;
    for (const auto& domain : partition.domains) {
        domains[domain.domainId] = &domain;
        for (const auto& neighbor : domain.interfaceGlobalDofsByNeighbor) {
            const auto pair = std::minmax(domain.domainId, neighbor.first);
            auto& target = candidateTargets[pair];
            neighbors[domain.domainId].insert(neighbor.first);
            neighbors[neighbor.first].insert(domain.domainId);
            for (int global : neighbor.second) {
                const int gamma =
                    partition.globalToInterface[static_cast<std::size_t>(global)];
                if (gamma >= 0) target.insert(gamma);
            }
        }
    }

    std::vector<PortPatch> patches;
    int fallbackId = static_cast<int>(mesh.interfaceSummaries.size());
    for (const auto& item : candidateTargets) {
        PortPatch patch;
        const auto found = interfaceIds.find(item.first);
        patch.interfaceId =
            found == interfaceIds.end() ? fallbackId++ : found->second;
        patch.leftSubdomain = item.first.first;
        patch.rightSubdomain = item.first.second;
        patch.target.assign(item.second.begin(), item.second.end());
        patches.push_back(std::move(patch));
    }
    std::sort(patches.begin(), patches.end(),
        [](const PortPatch& left, const PortPatch& right) {
            return left.interfaceId < right.interfaceId;
        });

    std::vector<int> owner(partition.interfaceGlobalDofs.size(), -1);
    for (PortPatch& patch : patches) {
        std::vector<int> target;
        target.reserve(patch.target.size());
        for (int gamma : patch.target) {
            if (owner[static_cast<std::size_t>(gamma)] < 0) {
                owner[static_cast<std::size_t>(gamma)] = patch.interfaceId;
                target.push_back(gamma);
            }
        }
        patch.target = std::move(target);
    }
    patches.erase(std::remove_if(patches.begin(), patches.end(),
        [](const PortPatch& patch) { return patch.target.empty(); }),
        patches.end());
    if (std::find(owner.begin(), owner.end(), -1) != owner.end()) {
        throw std::runtime_error(
            "[Optimal port] Target supports do not cover the SIPG trace.");
    }

    for (PortPatch& patch : patches) {
        std::set<int> patchDomains{
            patch.leftSubdomain, patch.rightSubdomain};
        std::set<int> frontier = patchDomains;
        for (int layer = 0; layer < oversamplingLayers; ++layer) {
            std::set<int> next;
            for (int domain : frontier) {
                for (int neighbor : neighbors[domain]) {
                    if (patchDomains.insert(neighbor).second) {
                        next.insert(neighbor);
                    }
                }
            }
            frontier = std::move(next);
        }
        patch.patchSubdomains.assign(
            patchDomains.begin(), patchDomains.end());

        std::set<int> source;
        for (int domainId : patch.patchSubdomains) {
            const auto found = domains.find(domainId);
            if (found == domains.end()) {
                throw std::runtime_error(
                    "[Optimal port] Patch subdomain is absent from the partition.");
            }
            for (const auto& neighbor :
                 found->second->interfaceGlobalDofsByNeighbor) {
                // Only the outer boundary of the oversampled patch is a
                // source. Interfaces internal to the patch are not prescribed
                // trace data and must not enter T_p's source space.
                if (patchDomains.count(neighbor.first) != 0) continue;
                for (int global : neighbor.second) {
                    const int gamma =
                        partition.globalToInterface[
                            static_cast<std::size_t>(global)];
                    if (gamma >= 0) source.insert(gamma);
                }
            }
        }
        for (int gamma : patch.target) source.erase(gamma);
        patch.source.assign(source.begin(), source.end());
        patch.targetFingerprint = fingerprint(patch.target);
        patch.sourceFingerprint = fingerprint(patch.source);
    }
    return patches;
}

std::vector<PortTopologyAudit> auditOptimalPortTopology(
    const Mesh& mesh,
    const CaseConfig& physics,
    const ddm_schur::InterfacePartition& partition,
    const std::vector<int>& sourceSubdomains,
    int oversamplingLayers,
    int requestedRank,
    const std::string& requestedInnerSolver,
    int directRowLimit)
{
    if (requestedRank <= 0 || directRowLimit <= 0) {
        throw std::runtime_error(
            "[Optimal port] Invalid topology-audit rank/solver limit.");
    }
    const std::vector<PortPatch> patches = buildOptimalPortPatches(
        mesh, partition, oversamplingLayers);
    std::set<std::pair<int, int>> internalBoundaryEntities;
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        const int left = mesh.tets[
            static_cast<std::size_t>(face.leftTet)].subdomain;
        const int right = mesh.tets[
            static_cast<std::size_t>(face.rightTet)].subdomain;
        internalBoundaryEntities.insert(
            {left, face.leftBoundaryEntity});
        internalBoundaryEntities.insert(
            {right, face.rightBoundaryEntity});
    }
    std::map<int, std::set<int>> externalEntitiesBySubdomain;
    for (const BoundaryFace& face : mesh.boundaryFaces) {
        if (internalBoundaryEntities.count(
                {face.subdomain, face.boundaryEntity}) == 0) {
            externalEntitiesBySubdomain[face.subdomain].insert(
                face.boundaryEntity);
        }
    }

    std::vector<PortTopologyAudit> result;
    result.reserve(patches.size());
    for (const PortPatch& patch : patches) {
        PortTopologyAudit row;
        row.interfaceId = patch.interfaceId;
        row.leftSubdomain = patch.leftSubdomain;
        row.rightSubdomain = patch.rightSubdomain;
        row.targetDofs = static_cast<int>(patch.target.size());
        row.sourceDofs = static_cast<int>(patch.source.size());
        row.sourceEmpty = patch.source.empty();
        const std::set<int> patchDomains(
            patch.patchSubdomains.begin(), patch.patchSubdomains.end());
        for (int subdomain : sourceSubdomains) {
            if (patchDomains.count(subdomain) != 0) {
                ++row.heatSourceChannelCount;
            }
        }
        std::set<std::tuple<std::string, int, int>> boundaryChannels;
        for (int subdomain : patch.patchSubdomains) {
            const auto& external = externalEntitiesBySubdomain[subdomain];
            for (const BoundaryCondition& condition :
                 physics.dirichletConditions) {
                if ((condition.subdomain < 0
                        || condition.subdomain == subdomain)
                    && external.count(condition.boundaryEntity) != 0) {
                    boundaryChannels.insert({
                        "dirichlet", subdomain,
                        condition.boundaryEntity});
                }
            }
            for (const ConvectionCondition& condition :
                 physics.convectionConditions) {
                if ((condition.subdomain < 0
                        || condition.subdomain == subdomain)
                    && external.count(condition.boundaryEntity) != 0) {
                    boundaryChannels.insert({
                        "convection", subdomain,
                        condition.boundaryEntity});
                }
            }
        }
        row.externalBoundaryChannelCount =
            static_cast<int>(boundaryChannels.size());

        int geometryModes = 1;
        if (!patch.target.empty()) {
            Vec3 minimum = mesh.nodes[static_cast<std::size_t>(
                partition.interfaceGlobalDofs[
                    static_cast<std::size_t>(patch.target.front())])].p;
            Vec3 maximum = minimum;
            for (int gamma : patch.target) {
                const Vec3& point = mesh.nodes[static_cast<std::size_t>(
                    partition.interfaceGlobalDofs[
                        static_cast<std::size_t>(gamma)])].p;
                minimum.x = std::min(minimum.x, point.x);
                minimum.y = std::min(minimum.y, point.y);
                minimum.z = std::min(minimum.z, point.z);
                maximum.x = std::max(maximum.x, point.x);
                maximum.y = std::max(maximum.y, point.y);
                maximum.z = std::max(maximum.z, point.z);
            }
            const double scale = std::max({
                1.0, std::abs(minimum.x), std::abs(minimum.y),
                std::abs(minimum.z), std::abs(maximum.x),
                std::abs(maximum.y), std::abs(maximum.z)});
            const double tolerance =
                256.0 * std::numeric_limits<double>::epsilon() * scale;
            geometryModes += maximum.x - minimum.x > tolerance ? 1 : 0;
            geometryModes += maximum.y - minimum.y > tolerance ? 1 : 0;
            geometryModes += maximum.z - minimum.z > tolerance ? 1 : 0;
        }
        row.mandatoryModeCount =
            geometryModes + row.heatSourceChannelCount;

        const std::size_t target =
            static_cast<std::size_t>(row.targetDofs);
        const std::size_t source =
            static_cast<std::size_t>(row.sourceDofs);
        const std::size_t rank = static_cast<std::size_t>(
            std::min(requestedRank, row.targetDofs));
        const std::size_t block = std::min(
            target, rank + static_cast<std::size_t>(4));
        row.eigensolverWorkspaceBytes =
            4 * target * block * sizeof(double);
        const bool direct = requestedInnerSolver == "direct"
            || (requestedInnerSolver == "auto"
                && row.targetDofs <= directRowLimit);
        if (direct) {
            row.innerSolverWorkspaceBytes =
                2 * target * target * sizeof(double);
        } else {
            constexpr std::size_t restart = 100;
            row.innerSolverWorkspaceBytes =
                (2 * restart + 8) * target * sizeof(double);
        }
        row.transferWorkspaceBytes =
            (5 * target + 3 * source) * sizeof(double);
        row.basisStorageEstimateBytes =
            target * rank * sizeof(double);
        row.estimatedWorkspaceBytes =
            row.eigensolverWorkspaceBytes
            + row.innerSolverWorkspaceBytes
            + row.transferWorkspaceBytes;
        result.push_back(row);
    }
    return result;
}

ReducedDynamicSchurOperator::ReducedDynamicSchurOperator(
    const local::Model& model,
    bool exactJacobiDiagonal)
    : model_(model)
{
    factors_.resize(model_.subdomains.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (std::ptrdiff_t slot = 0;
         slot < static_cast<std::ptrdiff_t>(model_.subdomains.size());
         ++slot) {
        const auto& subdomain =
            model_.subdomains[static_cast<std::size_t>(slot)];
        factors_[static_cast<std::size_t>(slot)] =
            local::factorDenseSymmetric(
                subdomain.reducedInterior, subdomain.rank);
    }

    std::vector<local::InterfaceEntry> entries = model_.interfaceEntries;
    std::sort(entries.begin(), entries.end(),
        [](const auto& left, const auto& right) {
            return left.row < right.row
                || (left.row == right.row
                    && left.column < right.column);
        });
    rowPtr_.assign(static_cast<std::size_t>(model_.interfaceDofs + 1), 0);
    for (std::size_t begin = 0; begin < entries.size();) {
        std::size_t end = begin + 1;
        double value = entries[begin].value;
        while (end < entries.size()
               && entries[end].row == entries[begin].row
               && entries[end].column == entries[begin].column) {
            value += entries[end].value;
            ++end;
        }
        if (value != 0.0) {
            columns_.push_back(entries[begin].column);
            values_.push_back(value);
            ++rowPtr_[static_cast<std::size_t>(entries[begin].row + 1)];
        }
        begin = end;
    }
    std::partial_sum(rowPtr_.begin(), rowPtr_.end(), rowPtr_.begin());

    diagonal_.assign(static_cast<std::size_t>(model_.interfaceDofs), 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int row = 0; row < model_.interfaceDofs; ++row) {
        for (int entry = rowPtr_[static_cast<std::size_t>(row)];
             entry < rowPtr_[static_cast<std::size_t>(row + 1)]; ++entry) {
            if (columns_[static_cast<std::size_t>(entry)] == row) {
                diagonal_[static_cast<std::size_t>(row)] +=
                    values_[static_cast<std::size_t>(entry)];
            }
        }
    }
    if (!exactJacobiDiagonal) {
        return;
    }
    constexpr int diagonalBlockSize = 64;
    for (std::size_t slot = 0; slot < model_.subdomains.size(); ++slot) {
        const auto& subdomain = model_.subdomains[slot];
        for (int firstGamma = 0;
             firstGamma < subdomain.localInterfaceDofs;
             firstGamma += diagonalBlockSize) {
            const int block = std::min(
                diagonalBlockSize,
                subdomain.localInterfaceDofs - firstGamma);
            std::vector<double> columns(static_cast<std::size_t>(
                subdomain.rank * block), 0.0);
            for (int mode = 0; mode < subdomain.rank; ++mode) {
                for (int local = 0; local < block; ++local) {
                    columns[static_cast<std::size_t>(
                        mode * block + local)] =
                        subdomain.reducedInteriorInterface[
                        static_cast<std::size_t>(
                            mode * subdomain.localInterfaceDofs
                            + firstGamma + local)];
                }
            }
            local::solveDenseSymmetricMultiple(
                factors_[slot], columns, block);
            if (firstGamma == 0) {
                std::vector<double> scalar(
                    static_cast<std::size_t>(subdomain.rank), 0.0);
                for (int mode = 0; mode < subdomain.rank; ++mode) {
                    scalar[static_cast<std::size_t>(mode)] =
                        subdomain.reducedInteriorInterface[
                            static_cast<std::size_t>(
                                mode * subdomain.localInterfaceDofs)];
                }
                local::solveDenseSymmetric(factors_[slot], scalar);
                double difference = 0.0;
                double magnitude = 0.0;
                for (int mode = 0; mode < subdomain.rank; ++mode) {
                    const double delta =
                        columns[static_cast<std::size_t>(mode * block)]
                        - scalar[static_cast<std::size_t>(mode)];
                    difference = std::max(difference, std::abs(delta));
                    magnitude = std::max(
                        magnitude,
                        std::abs(scalar[static_cast<std::size_t>(mode)]));
                }
                if (difference > 1.0e-10
                        * std::max(1.0, magnitude)) {
                    throw std::runtime_error(
                        "[Optimal port] Batched exact Schur diagonal "
                        "solve disagrees with scalar solve.");
                }
            }
            for (int local = 0; local < block; ++local) {
                const int gamma = firstGamma + local;
                double correction = 0.0;
                for (int mode = 0; mode < subdomain.rank; ++mode) {
                    correction += subdomain.reducedInterfaceInterior[
                        static_cast<std::size_t>(
                            gamma * subdomain.rank + mode)]
                        * columns[static_cast<std::size_t>(
                            mode * block + local)];
                }
                const int globalGamma =
                    subdomain.interfaceIndices[
                        static_cast<std::size_t>(gamma)];
                diagonal_[static_cast<std::size_t>(globalGamma)]
                    -= correction;
            }
        }
    }
}

void ReducedDynamicSchurOperator::solveReducedInterior(
    std::size_t subdomainSlot,
    std::vector<double>& rightHandSide) const
{
    if (subdomainSlot >= factors_.size()) {
        throw std::runtime_error(
            "[Optimal port] Reduced-interior factor slot is invalid.");
    }
    local::solveDenseSymmetric(
        factors_[subdomainSlot], rightHandSide);
}

void ReducedDynamicSchurOperator::solveReducedInteriorMultiple(
    std::size_t subdomainSlot,
    std::vector<double>& rightHandSides,
    int rightHandSideCount) const
{
    if (subdomainSlot >= factors_.size()) {
        throw std::runtime_error(
            "[Optimal port] Reduced-interior factor slot is invalid.");
    }
    local::solveDenseSymmetricMultiple(
        factors_[subdomainSlot], rightHandSides, rightHandSideCount);
}

void ReducedDynamicSchurOperator::apply(
    const std::vector<double>& input,
    std::vector<double>& output) const
{
    if (input.size() != static_cast<std::size_t>(model_.interfaceDofs)) {
        throw std::runtime_error(
            "[Optimal port] Schur input size mismatch.");
    }
    output.assign(input.size(), 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int row = 0; row < model_.interfaceDofs; ++row) {
        for (int entry = rowPtr_[static_cast<std::size_t>(row)];
             entry < rowPtr_[static_cast<std::size_t>(row + 1)]; ++entry) {
            output[static_cast<std::size_t>(row)] +=
                values_[static_cast<std::size_t>(entry)]
                * input[static_cast<std::size_t>(
                    columns_[static_cast<std::size_t>(entry)])];
        }
    }
    std::vector<std::vector<double>> corrections(
        model_.subdomains.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (std::ptrdiff_t slot = 0;
         slot < static_cast<std::ptrdiff_t>(model_.subdomains.size());
         ++slot) {
        const auto& subdomain =
            model_.subdomains[static_cast<std::size_t>(slot)];
        std::vector<double> eliminated(
            static_cast<std::size_t>(subdomain.rank), 0.0);
        for (int mode = 0; mode < subdomain.rank; ++mode) {
            for (int gamma = 0;
                 gamma < subdomain.localInterfaceDofs; ++gamma) {
                eliminated[static_cast<std::size_t>(mode)] +=
                    subdomain.reducedInteriorInterface[
                        static_cast<std::size_t>(
                            mode * subdomain.localInterfaceDofs + gamma)]
                    * input[static_cast<std::size_t>(
                        subdomain.interfaceIndices[
                            static_cast<std::size_t>(gamma)])];
            }
        }
        local::solveDenseSymmetric(
            factors_[static_cast<std::size_t>(slot)], eliminated);
        corrections[static_cast<std::size_t>(slot)].assign(
            static_cast<std::size_t>(
                subdomain.localInterfaceDofs), 0.0);
        for (int gamma = 0;
             gamma < subdomain.localInterfaceDofs; ++gamma) {
            double correction = 0.0;
            for (int mode = 0; mode < subdomain.rank; ++mode) {
                correction += subdomain.reducedInterfaceInterior[
                    static_cast<std::size_t>(
                        gamma * subdomain.rank + mode)]
                    * eliminated[static_cast<std::size_t>(mode)];
            }
            corrections[static_cast<std::size_t>(slot)][
                static_cast<std::size_t>(gamma)] = correction;
        }
    }
    for (std::size_t slot = 0;
         slot < model_.subdomains.size(); ++slot) {
        const auto& subdomain = model_.subdomains[slot];
        for (int gamma = 0;
             gamma < subdomain.localInterfaceDofs; ++gamma) {
            output[static_cast<std::size_t>(
                subdomain.interfaceIndices[
                    static_cast<std::size_t>(gamma)])] -=
                corrections[slot][static_cast<std::size_t>(gamma)];
        }
    }
}

void ReducedDynamicSchurOperator::applyTranspose(
    const std::vector<double>& input,
    std::vector<double>& output) const
{
    if (input.size() != static_cast<std::size_t>(model_.interfaceDofs)) {
        throw std::runtime_error(
            "[Optimal port] Schur transpose input size mismatch.");
    }
    output.assign(input.size(), 0.0);
    std::vector<std::vector<double>> sparseTranspose(
        static_cast<std::size_t>(
            std::max(1,
#ifdef _OPENMP
                omp_get_max_threads()
#else
                1
#endif
                )),
        std::vector<double>(input.size(), 0.0));
#ifdef _OPENMP
#pragma omp parallel
    {
        const int thread = omp_get_thread_num();
#pragma omp for schedule(static)
#endif
    for (int row = 0; row < model_.interfaceDofs; ++row) {
        for (int entry = rowPtr_[static_cast<std::size_t>(row)];
             entry < rowPtr_[static_cast<std::size_t>(row + 1)]; ++entry) {
            sparseTranspose[
#ifdef _OPENMP
                static_cast<std::size_t>(thread)
#else
                0
#endif
                ][static_cast<std::size_t>(
                columns_[static_cast<std::size_t>(entry)])] +=
                values_[static_cast<std::size_t>(entry)]
                * input[static_cast<std::size_t>(row)];
        }
    }
#ifdef _OPENMP
    }
#endif
    for (const std::vector<double>& localOutput : sparseTranspose) {
        for (std::size_t row = 0; row < output.size(); ++row) {
            output[row] += localOutput[row];
        }
    }
    std::vector<std::vector<double>> corrections(
        model_.subdomains.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (std::ptrdiff_t slot = 0;
         slot < static_cast<std::ptrdiff_t>(model_.subdomains.size());
         ++slot) {
        const auto& subdomain =
            model_.subdomains[static_cast<std::size_t>(slot)];
        std::vector<double> eliminated(
            static_cast<std::size_t>(subdomain.rank), 0.0);
        for (int mode = 0; mode < subdomain.rank; ++mode) {
            for (int gamma = 0;
                 gamma < subdomain.localInterfaceDofs; ++gamma) {
                eliminated[static_cast<std::size_t>(mode)] +=
                    subdomain.reducedInterfaceInterior[
                        static_cast<std::size_t>(
                            gamma * subdomain.rank + mode)]
                    * input[static_cast<std::size_t>(
                        subdomain.interfaceIndices[
                            static_cast<std::size_t>(gamma)])];
            }
        }
        local::solveDenseSymmetric(
            factors_[static_cast<std::size_t>(slot)], eliminated);
        corrections[static_cast<std::size_t>(slot)].assign(
            static_cast<std::size_t>(
                subdomain.localInterfaceDofs), 0.0);
        for (int gamma = 0;
             gamma < subdomain.localInterfaceDofs; ++gamma) {
            double correction = 0.0;
            for (int mode = 0; mode < subdomain.rank; ++mode) {
                correction += subdomain.reducedInteriorInterface[
                    static_cast<std::size_t>(
                        mode * subdomain.localInterfaceDofs + gamma)]
                    * eliminated[static_cast<std::size_t>(mode)];
            }
            corrections[static_cast<std::size_t>(slot)][
                static_cast<std::size_t>(gamma)] = correction;
        }
    }
    for (std::size_t slot = 0;
         slot < model_.subdomains.size(); ++slot) {
        const auto& subdomain = model_.subdomains[slot];
        for (int gamma = 0;
             gamma < subdomain.localInterfaceDofs; ++gamma) {
            output[static_cast<std::size_t>(
                subdomain.interfaceIndices[
                    static_cast<std::size_t>(gamma)])] -=
                corrections[slot][static_cast<std::size_t>(gamma)];
        }
    }
}

PatchTransferOperator::PatchTransferOperator(
    const ReducedDynamicSchurOperator& schur,
    PortPatch patch,
    const PatchTransferOptions& options)
    : schur_(schur), patch_(std::move(patch)), options_(options)
{
    statistics_.requestedSolver = options_.innerSolver;
    if (patch_.target.empty()) {
        throw std::runtime_error(
            "[Optimal port] Transfer patch has an empty target.");
    }
    if (options_.innerSolver != "auto"
        && options_.innerSolver != "direct"
        && options_.innerSolver != "pcg"
        && options_.innerSolver != "fgmres"
        && options_.innerSolver != "iterative-schur"
        && options_.innerSolver != "assembled-dense"
        && options_.innerSolver != "woodbury-exact") {
        throw std::runtime_error(
            "[Optimal port] Inner solver must be iterative-schur, "
            "assembled-dense, woodbury-exact, or a supported legacy alias.");
    }
    if (options_.directRowLimit <= 0
        || !(options_.relativeTolerance > 0.0)
        || options_.maximumIterations <= 0
        || options_.refinementMaximumIterations < 0
        || options_.refinementMaximumIterations > 3
        || !(options_.refinementTolerance > 0.0)) {
        throw std::runtime_error(
            "[Optimal port] Invalid inner-solver options.");
    }
    const bool legacyDirect = options_.innerSolver == "direct";
    const bool legacyIterative =
        options_.innerSolver == "pcg"
        || options_.innerSolver == "fgmres";
    assembledDense_ = options_.innerSolver == "assembled-dense"
        || legacyDirect
        || (options_.innerSolver == "auto"
            && targetRows() <= options_.directRowLimit);
    woodbury_ = options_.innerSolver == "woodbury-exact"
        || (options_.innerSolver == "auto"
            && targetRows() > options_.directRowLimit);
    forceFgmres_ = options_.innerSolver == "fgmres";
    if (legacyDirect && targetRows() > options_.directRowLimit) {
        throw std::runtime_error(
            "[Optimal port] Legacy direct S_tt row limit exceeded; "
            "use assembled-dense explicitly for a reference run.");
    }
    direct_ = assembledDense_ || woodbury_;

    const auto start = Clock::now();
    initializePatchAlgebra();
    if (assembledDense_) {
        initializeAssembledDense();
    } else if (woodbury_) {
        initializeWoodbury();
    } else {
        (void)legacyIterative;
        statistics_.solver = forceFgmres_
            ? "matrix-free-fgmres" : "matrix-free-pcg";
        statistics_.actualSolver = statistics_.solver;
        initializeIterativePreconditioner();
    }
    statistics_.setupSeconds = elapsed(start);
    statistics_.status = "ready";
}

PatchTransferOperator::~PatchTransferOperator() = default;

void PatchTransferOperator::initializePatchAlgebra()
{
    const auto start = Clock::now();
    algebra_ = std::make_unique<PatchAlgebra>();
    const local::Model& model = schur_.model();
    std::vector<int> targetPosition(
        static_cast<std::size_t>(model.interfaceDofs), -1);
    std::vector<int> sourcePosition(
        static_cast<std::size_t>(model.interfaceDofs), -1);
    for (int row = 0; row < targetRows(); ++row) {
        targetPosition[static_cast<std::size_t>(
            patch_.target[static_cast<std::size_t>(row)])] = row;
    }
    for (int row = 0; row < sourceRows(); ++row) {
        sourcePosition[static_cast<std::size_t>(
            patch_.source[static_cast<std::size_t>(row)])] = row;
    }

    std::map<std::pair<int, int>, double> targetValues;
    std::map<std::pair<int, int>, double> crossValues;
    for (const local::InterfaceEntry& entry : model.interfaceEntries) {
        if (entry.row < 0 || entry.row >= model.interfaceDofs
            || entry.column < 0 || entry.column >= model.interfaceDofs) {
            throw std::runtime_error(
                "[Optimal port] Interface entry index is invalid.");
        }
        const int targetRow =
            targetPosition[static_cast<std::size_t>(entry.row)];
        const int targetColumn =
            targetPosition[static_cast<std::size_t>(entry.column)];
        const int sourceColumn =
            sourcePosition[static_cast<std::size_t>(entry.column)];
        if (targetRow >= 0 && targetColumn >= 0) {
            targetValues[{targetRow, targetColumn}] += entry.value;
        }
        if (targetRow >= 0 && sourceColumn >= 0) {
            crossValues[{targetRow, sourceColumn}] += entry.value;
        }
    }
    long double baseAsymmetry = 0.0L;
    long double baseMagnitude = 0.0L;
    for (const auto& item : targetValues) {
        if (item.second == 0.0) continue;
        algebra_->targetEntries.push_back({
            item.first.first, item.first.second, item.second});
        algebra_->sparseTargetEntries.push_back({
            item.first.first, item.first.second, item.second});
        const auto transpose = targetValues.find({
            item.first.second, item.first.first});
        const double transposeValue =
            transpose == targetValues.end() ? 0.0 : transpose->second;
        const double difference = item.second - transposeValue;
        baseAsymmetry +=
            static_cast<long double>(difference) * difference;
        baseMagnitude +=
            static_cast<long double>(item.second) * item.second;
    }
    for (const auto& item : crossValues) {
        if (item.second != 0.0) {
            algebra_->crossEntries.push_back({
                item.first.first, item.first.second, item.second});
        }
    }

    long double couplingAsymmetry = 0.0L;
    long double couplingMagnitude = 0.0L;
    int correctionOffset = 0;
    for (std::size_t slot = 0;
         slot < model.subdomains.size(); ++slot) {
        const local::SubdomainModel& subdomain =
            model.subdomains[slot];
        PatchAlgebra::LocalCorrection correction;
        correction.subdomainSlot = slot;
        correction.rank = subdomain.rank;
        for (int gamma = 0;
             gamma < subdomain.localInterfaceDofs; ++gamma) {
            const int globalGamma = subdomain.interfaceIndices[
                static_cast<std::size_t>(gamma)];
            const int targetRow =
                targetPosition[static_cast<std::size_t>(globalGamma)];
            const int sourceRow =
                sourcePosition[static_cast<std::size_t>(globalGamma)];
            if (targetRow >= 0) {
                correction.targetLocalDofs.push_back(gamma);
                correction.targetRows.push_back(targetRow);
                for (int mode = 0; mode < subdomain.rank; ++mode) {
                    const double left =
                        subdomain.reducedInterfaceInterior[
                            static_cast<std::size_t>(
                                gamma * subdomain.rank + mode)];
                    const double right =
                        subdomain.reducedInteriorInterface[
                            static_cast<std::size_t>(
                                mode * subdomain.localInterfaceDofs
                                + gamma)];
                    const double difference = left - right;
                    couplingAsymmetry +=
                        static_cast<long double>(difference) * difference;
                    couplingMagnitude +=
                        static_cast<long double>(left) * left;
                }
            }
            if (sourceRow >= 0) {
                correction.sourceLocalDofs.push_back(gamma);
                correction.sourceRows.push_back(sourceRow);
            }
        }
        if (!correction.targetRows.empty()) {
            correction.correctionOffset = correctionOffset;
            correctionOffset += correction.rank;
            algebra_->localCorrections.push_back(
                std::move(correction));
        }
    }
    statistics_.aTtDimension = targetRows();
    statistics_.aTtNonzeros = algebra_->targetEntries.size();
    statistics_.reducedCorrectionRank = correctionOffset;
    const double baseRelative = std::sqrt(
        static_cast<double>(baseAsymmetry))
        / std::max(1.0e-300,
            std::sqrt(static_cast<double>(baseMagnitude)));
    const double couplingRelative = std::sqrt(
        static_cast<double>(couplingAsymmetry))
        / std::max(1.0e-300,
            std::sqrt(static_cast<double>(couplingMagnitude)));
    statistics_.relativeAsymmetry =
        std::max(baseRelative, couplingRelative);
    statistics_.aTtAssemblySeconds = elapsed(start);

    const auto referenceStart = Clock::now();
    std::vector<double> targetProbe(
        patch_.target.size(), 0.0);
    std::vector<double> full(
        static_cast<std::size_t>(schur_.size()), 0.0);
    for (std::size_t row = 0; row < targetProbe.size(); ++row) {
        targetProbe[row] = std::sin(
            0.319 * static_cast<double>(row + 1));
        full[static_cast<std::size_t>(patch_.target[row])] =
            targetProbe[row];
    }
    std::vector<double> localTarget;
    applyTargetBlock(targetProbe, localTarget, false);
    std::vector<double> fullImage;
    schur_.apply(full, fullImage);
    std::vector<double> referenceTarget(
        patch_.target.size(), 0.0);
    for (std::size_t row = 0;
         row < patch_.target.size(); ++row) {
        referenceTarget[row] = fullImage[static_cast<std::size_t>(
            patch_.target[row])];
    }
    statistics_.referenceTargetActionError =
        relativeDifference(localTarget, referenceTarget);
    if (!patch_.source.empty()) {
        std::vector<double> sourceProbe(
            patch_.source.size(), 0.0);
        full.assign(
            static_cast<std::size_t>(schur_.size()), 0.0);
        for (std::size_t row = 0; row < sourceProbe.size(); ++row) {
            sourceProbe[row] = std::cos(
                0.211 * static_cast<double>(row + 1));
            full[static_cast<std::size_t>(patch_.source[row])] =
                sourceProbe[row];
        }
        std::vector<double> localCross;
        applyCrossBlock(sourceProbe, localCross, false);
        schur_.apply(full, fullImage);
        std::vector<double> referenceCross(
            patch_.target.size(), 0.0);
        for (std::size_t row = 0;
             row < patch_.target.size(); ++row) {
            referenceCross[row] = fullImage[static_cast<std::size_t>(
                patch_.target[row])];
        }
        statistics_.referenceCrossActionError =
            relativeDifference(localCross, referenceCross);

        full.assign(
            static_cast<std::size_t>(schur_.size()), 0.0);
        for (std::size_t row = 0; row < targetProbe.size(); ++row) {
            full[static_cast<std::size_t>(patch_.target[row])] =
                targetProbe[row];
        }
        std::vector<double> localTranspose;
        applyCrossBlock(targetProbe, localTranspose, true);
        schur_.applyTranspose(full, fullImage);
        std::vector<double> referenceTranspose(
            patch_.source.size(), 0.0);
        for (std::size_t row = 0;
             row < patch_.source.size(); ++row) {
            referenceTranspose[row] =
                fullImage[static_cast<std::size_t>(
                    patch_.source[row])];
        }
        statistics_.referenceCrossTransposeError =
            relativeDifference(localTranspose, referenceTranspose);
    }
    statistics_.referenceActionCheckSeconds =
        elapsed(referenceStart);
    if (statistics_.referenceTargetActionError > 1.0e-10
        || statistics_.referenceCrossActionError > 1.0e-10
        || statistics_.referenceCrossTransposeError > 1.0e-10) {
        throw std::runtime_error(
            "[Optimal port] Patch-local reduced-block action "
            "does not match the exact global Schur reference.");
    }
}

void PatchTransferOperator::initializeAssembledDense()
{
    const int rows = targetRows();
    const auto assemblyStart = Clock::now();
    std::vector<double> matrix(
        static_cast<std::size_t>(rows) * rows, 0.0);
    for (const PatchAlgebra::Entry& entry :
         algebra_->targetEntries) {
        matrix[static_cast<std::size_t>(
            entry.row * rows + entry.column)] += entry.value;
    }
    std::size_t peakTemporary = matrix.capacity() * sizeof(double);
    const local::Model& model = schur_.model();
    for (const PatchAlgebra::LocalCorrection& correction :
         algebra_->localCorrections) {
        const local::SubdomainModel& subdomain =
            model.subdomains[correction.subdomainSlot];
        std::vector<double> left(
            static_cast<std::size_t>(rows * correction.rank), 0.0);
        std::vector<double> solved(
            static_cast<std::size_t>(correction.rank * rows), 0.0);
        for (std::size_t item = 0;
             item < correction.targetRows.size(); ++item) {
            const int gamma = correction.targetLocalDofs[item];
            const int targetRow = correction.targetRows[item];
            for (int mode = 0; mode < correction.rank; ++mode) {
                left[static_cast<std::size_t>(
                    targetRow * correction.rank + mode)] =
                    subdomain.reducedInterfaceInterior[
                        static_cast<std::size_t>(
                            gamma * correction.rank + mode)];
                solved[static_cast<std::size_t>(
                    mode * rows + targetRow)] =
                    subdomain.reducedInteriorInterface[
                        static_cast<std::size_t>(
                            mode * subdomain.localInterfaceDofs
                            + gamma)];
            }
        }
        schur_.solveReducedInteriorMultiple(
            correction.subdomainSlot, solved, rows);
#ifdef USE_MKL_PARDISO
        cblas_dgemm(
            CblasRowMajor, CblasNoTrans, CblasNoTrans,
            rows, rows, correction.rank,
            -1.0, left.data(), correction.rank,
            solved.data(), rows,
            1.0, matrix.data(), rows);
#else
        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < rows; ++column) {
                long double value = 0.0L;
                for (int mode = 0; mode < correction.rank; ++mode) {
                    value += static_cast<long double>(left[
                        static_cast<std::size_t>(
                            row * correction.rank + mode)])
                        * solved[static_cast<std::size_t>(
                            mode * rows + column)];
                }
                matrix[static_cast<std::size_t>(
                    row * rows + column)] -=
                    static_cast<double>(value);
            }
        }
#endif
        peakTemporary = std::max(
            peakTemporary,
            matrix.capacity() * sizeof(double)
            + left.capacity() * sizeof(double)
            + solved.capacity() * sizeof(double));
    }
    statistics_.aTtAssemblySeconds += elapsed(assemblyStart);

    long double asymmetry = 0.0L;
    long double magnitude = 0.0L;
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < rows; ++column) {
            const double value = matrix[static_cast<std::size_t>(
                row * rows + column)];
            const double transpose = matrix[static_cast<std::size_t>(
                column * rows + row)];
            const double difference = value - transpose;
            asymmetry +=
                static_cast<long double>(difference) * difference;
            magnitude += static_cast<long double>(value) * value;
        }
    }
    statistics_.relativeAsymmetry = std::max(
        statistics_.relativeAsymmetry,
        std::sqrt(static_cast<double>(asymmetry))
            / std::max(1.0e-300,
                std::sqrt(static_cast<double>(magnitude))));
    const auto factorStart = Clock::now();
    if (statistics_.relativeAsymmetry <= 1.0e-10) {
        factor_ = local::factorDenseSymmetric(matrix, rows);
        statistics_.actualSolver = factor_.cholesky
            ? "assembled-dense-llt" : "assembled-dense-ldlt";
        if (!factor_.cholesky) {
            statistics_.fallbackReason =
                "dense_llt_failed_using_ldlt";
        }
        statistics_.solver = statistics_.actualSolver;
        statistics_.aTtFactorBytes =
            factor_.lower.capacity() * sizeof(double)
            + factor_.diagonal.capacity() * sizeof(double);
    } else {
        std::vector<MatrixEntry> entries;
        entries.reserve(matrix.size());
        for (int row = 0; row < rows; ++row) {
            for (int column = 0; column < rows; ++column) {
                const double value = matrix[static_cast<std::size_t>(
                    row * rows + column)];
                if (value != 0.0) {
                    entries.push_back({row, column, value});
                }
            }
        }
        std::string fallback = "assembled_dense_asymmetric";
        sparseFactor_ = std::make_unique<SparseTargetFactor>(
            rows, entries, false, fallback);
        statistics_.fallbackReason = fallback;
        statistics_.actualSolver = "assembled-dense-pardiso-general";
        statistics_.solver = statistics_.actualSolver;
        statistics_.aTtFactorBytes = sparseFactor_->memoryBytes();
    }
    statistics_.aTtFactorizationSeconds = elapsed(factorStart);
    statistics_.factorBytes = statistics_.aTtFactorBytes;
    statistics_.peakIncrementalMemoryBytes =
        std::max(peakTemporary, statistics_.factorBytes);
}

void PatchTransferOperator::initializeWoodbury()
{
    const int rows = targetRows();
    const int rank = statistics_.reducedCorrectionRank;
    if (rank <= 0) {
        throw std::runtime_error(
            "[Optimal port] Woodbury target has no reduced correction.");
    }
    const bool numericallySymmetric =
        statistics_.relativeAsymmetry <= 1.0e-10;
    const auto factorStart = Clock::now();
    std::string aFallback;
    sparseFactor_ = std::make_unique<SparseTargetFactor>(
        rows, algebra_->sparseTargetEntries,
        numericallySymmetric, aFallback);
#ifdef USE_MKL_PARDISO
    // SubdomainDirectSolver and GeneralSparseDirectSolver both configure
    // PARDISO iparm[7]=2. This is recorded separately from the outer
    // exact-Schur iterative refinement implemented below.
    statistics_.pardisoInternalRefinementSteps = 2;
#endif
    statistics_.aTtFactorizationSeconds = elapsed(factorStart);
    statistics_.aTtFactorBytes = sparseFactor_->memoryBytes();
    if (!aFallback.empty()) {
        statistics_.fallbackReason = aFallback;
    }

    woodburyU_.assign(
        static_cast<std::size_t>(rows * rank), 0.0);
    woodburyV_.assign(
        static_cast<std::size_t>(rows * rank), 0.0);
    std::vector<double> h(
        static_cast<std::size_t>(rank * rank), 0.0);
    const local::Model& model = schur_.model();
    for (const PatchAlgebra::LocalCorrection& correction :
         algebra_->localCorrections) {
        const local::SubdomainModel& subdomain =
            model.subdomains[correction.subdomainSlot];
        const int offset = correction.correctionOffset;
        for (int row = 0; row < correction.rank; ++row) {
            for (int column = 0; column < correction.rank; ++column) {
                h[static_cast<std::size_t>(
                    (offset + row) * rank + offset + column)] =
                    subdomain.reducedInterior[static_cast<std::size_t>(
                        row * correction.rank + column)];
            }
        }
        for (std::size_t item = 0;
             item < correction.targetRows.size(); ++item) {
            const int gamma = correction.targetLocalDofs[item];
            const int targetRow = correction.targetRows[item];
            for (int mode = 0; mode < correction.rank; ++mode) {
                woodburyU_[static_cast<std::size_t>(
                    (offset + mode) * rows + targetRow)] =
                    subdomain.reducedInterfaceInterior[
                        static_cast<std::size_t>(
                            gamma * correction.rank + mode)];
                woodburyV_[static_cast<std::size_t>(
                    (offset + mode) * rows + targetRow)] =
                    subdomain.reducedInteriorInterface[
                        static_cast<std::size_t>(
                            mode * subdomain.localInterfaceDofs
                            + gamma)];
            }
        }
    }

    const auto wStart = Clock::now();
    sparseFactor_->solveMultiple(
        woodburyU_, rank, false, woodburyW_);
    if (numericallySymmetric && sparseFactor_->symmetric()) {
        woodburyWTranspose_.clear();
    } else {
        sparseFactor_->solveMultiple(
            woodburyV_, rank, true, woodburyWTranspose_);
    }
    statistics_.wSetupSeconds = elapsed(wStart);
    statistics_.wBytes =
        (woodburyW_.capacity()
         + woodburyWTranspose_.capacity()) * sizeof(double);

    const auto qAssemblyStart = Clock::now();
    std::vector<double> q = h;
#ifdef USE_MKL_PARDISO
    cblas_dgemm(
        CblasRowMajor, CblasNoTrans, CblasTrans,
        rank, rank, rows,
        -1.0, woodburyV_.data(), rows,
        woodburyW_.data(), rows,
        1.0, q.data(), rank);
#else
    for (int row = 0; row < rank; ++row) {
        for (int column = 0; column < rank; ++column) {
            long double value = 0.0L;
            for (int target = 0; target < rows; ++target) {
                value += static_cast<long double>(
                    woodburyV_[static_cast<std::size_t>(
                        row * rows + target)])
                    * woodburyW_[static_cast<std::size_t>(
                        column * rows + target)];
            }
            q[static_cast<std::size_t>(
                row * rank + column)] -= static_cast<double>(value);
        }
    }
#endif
    statistics_.qAssemblySeconds = elapsed(qAssemblyStart);
    statistics_.qDimension = rank;
    long double qAsymmetry = 0.0L;
    long double qMagnitude = 0.0L;
    for (int row = 0; row < rank; ++row) {
        for (int column = 0; column < rank; ++column) {
            const double value = q[static_cast<std::size_t>(
                row * rank + column)];
            const double transpose = q[static_cast<std::size_t>(
                column * rank + row)];
            const double difference = value - transpose;
            qAsymmetry +=
                static_cast<long double>(difference) * difference;
            qMagnitude += static_cast<long double>(value) * value;
        }
    }
    const double qRelativeAsymmetry =
        std::sqrt(static_cast<double>(qAsymmetry))
        / std::max(1.0e-300,
            std::sqrt(static_cast<double>(qMagnitude)));
    const auto qFactorStart = Clock::now();
    if (qRelativeAsymmetry <= 1.0e-10) {
        qFactor_ = local::factorDenseSymmetric(q, rank);
        double minimumFactorDiagonal =
            std::numeric_limits<double>::infinity();
        double maximumFactorDiagonal = 0.0;
        for (int row = 0; row < rank; ++row) {
            const double value = std::abs(
                qFactor_.cholesky
                    ? qFactor_.lower[static_cast<std::size_t>(
                        row * rank + row)]
                    : qFactor_.diagonal[
                        static_cast<std::size_t>(row)]);
            minimumFactorDiagonal =
                std::min(minimumFactorDiagonal, value);
            maximumFactorDiagonal =
                std::max(maximumFactorDiagonal, value);
        }
        statistics_.qMinimumAbsoluteFactorDiagonal =
            std::isfinite(minimumFactorDiagonal)
                ? minimumFactorDiagonal : 0.0;
        statistics_.qMaximumAbsoluteFactorDiagonal =
            maximumFactorDiagonal;
        statistics_.qFactorDiagonalRatio =
            maximumFactorDiagonal
            / std::max(std::numeric_limits<double>::min(),
                statistics_.qMinimumAbsoluteFactorDiagonal);
        statistics_.qBytes =
            qFactor_.lower.capacity() * sizeof(double)
            + qFactor_.diagonal.capacity() * sizeof(double);
        statistics_.actualSolver =
            std::string("woodbury-exact:pardiso-")
            + (sparseFactor_->symmetric() ? "spd" : "general")
            + (qFactor_.cholesky ? "+q-llt" : "+q-ldlt");
        if (!qFactor_.cholesky) {
            if (!statistics_.fallbackReason.empty()) {
                statistics_.fallbackReason += ';';
            }
            statistics_.fallbackReason +=
                "q_llt_failed_using_ldlt";
        }
    } else {
        std::vector<MatrixEntry> entries;
        entries.reserve(q.size());
        for (int row = 0; row < rank; ++row) {
            for (int column = 0; column < rank; ++column) {
                const double value = q[static_cast<std::size_t>(
                    row * rank + column)];
                if (value != 0.0) {
                    entries.push_back({row, column, value});
                }
            }
        }
        std::string qFallback = "woodbury_q_asymmetric";
        qSparseFactor_ = std::make_unique<SparseTargetFactor>(
            rank, entries, false, qFallback);
        if (!statistics_.fallbackReason.empty()) {
            statistics_.fallbackReason += ';';
        }
        statistics_.fallbackReason += qFallback;
        statistics_.qBytes = qSparseFactor_->memoryBytes();
        statistics_.actualSolver =
            "woodbury-exact:pardiso-general+q-general";
    }
    woodburyQ_ = std::move(q);
    statistics_.qBytes +=
        woodburyQ_.capacity() * sizeof(double);
    statistics_.qFactorizationSeconds = elapsed(qFactorStart);
    statistics_.solver = statistics_.actualSolver;
    statistics_.factorBytes =
        statistics_.aTtFactorBytes
        + statistics_.wBytes
        + statistics_.qBytes;
    statistics_.peakIncrementalMemoryBytes =
        std::max(
            statistics_.factorBytes
                + (woodburyU_.capacity()
                   + woodburyV_.capacity()
                   + h.capacity()) * sizeof(double),
            statistics_.factorBytes
                + (woodburyU_.capacity()
                   + woodburyV_.capacity()
                   + static_cast<std::size_t>(8 * rows)
                   + static_cast<std::size_t>(4 * rank))
                    * sizeof(double));
}

void PatchTransferOperator::initializeIterativePreconditioner()
{
    preconditionerDiagonal_.clear();
    preconditionerDiagonal_.reserve(patch_.target.size());
    double maximum = 0.0;
    bool nonPositive = false;
    for (int gamma : patch_.target) {
        const double value =
            schur_.diagonal()[static_cast<std::size_t>(gamma)];
        if (!std::isfinite(value)) {
            throw std::runtime_error(
                "[Optimal port] S_tt Jacobi diagonal is not finite.");
        }
        nonPositive = nonPositive || !(value > 0.0);
        maximum = std::max(maximum, std::abs(value));
        preconditionerDiagonal_.push_back(std::abs(value));
    }
    if (!(maximum > 0.0)) {
        throw std::runtime_error(
            "[Optimal port] S_tt Jacobi diagonal is identically zero.");
    }
    const double floor = 128.0 * std::numeric_limits<double>::epsilon()
        * maximum;
    for (double& value : preconditionerDiagonal_) {
        value = std::max(value, floor);
    }
    if (nonPositive && !forceFgmres_) {
        forceFgmres_ = true;
        statistics_.fallbackReason = "non_positive_jacobi_diagonal";
        statistics_.solver = "matrix-free-fgmres";
        statistics_.actualSolver = "matrix-free-fgmres";
    }
    statistics_.diagonalShift = floor;
    statistics_.factorBytes =
        preconditionerDiagonal_.capacity() * sizeof(double);
}

void PatchTransferOperator::applyTargetBlock(
    const std::vector<double>& input,
    std::vector<double>& output,
    bool transpose) const
{
    if (input.size() != patch_.target.size()) {
        throw std::runtime_error(
            "[Optimal port] S_tt input size mismatch.");
    }
    output.assign(input.size(), 0.0);
    for (const PatchAlgebra::Entry& entry :
         algebra_->targetEntries) {
        if (transpose) {
            output[static_cast<std::size_t>(entry.column)] +=
                entry.value * input[static_cast<std::size_t>(entry.row)];
        } else {
            output[static_cast<std::size_t>(entry.row)] +=
                entry.value * input[static_cast<std::size_t>(entry.column)];
        }
    }
    const local::Model& model = schur_.model();
    for (const PatchAlgebra::LocalCorrection& correction :
         algebra_->localCorrections) {
        const local::SubdomainModel& subdomain =
            model.subdomains[correction.subdomainSlot];
        std::vector<double> eliminated(
            static_cast<std::size_t>(correction.rank), 0.0);
        for (std::size_t item = 0;
             item < correction.targetRows.size(); ++item) {
            const int gamma = correction.targetLocalDofs[item];
            const int targetRow = correction.targetRows[item];
            for (int mode = 0; mode < correction.rank; ++mode) {
                eliminated[static_cast<std::size_t>(mode)] +=
                    (transpose
                        ? subdomain.reducedInterfaceInterior[
                            static_cast<std::size_t>(
                                gamma * correction.rank + mode)]
                        : subdomain.reducedInteriorInterface[
                            static_cast<std::size_t>(
                                mode * subdomain.localInterfaceDofs
                                + gamma)])
                    * input[static_cast<std::size_t>(targetRow)];
            }
        }
        schur_.solveReducedInterior(
            correction.subdomainSlot, eliminated);
        for (std::size_t item = 0;
             item < correction.targetRows.size(); ++item) {
            const int gamma = correction.targetLocalDofs[item];
            const int targetRow = correction.targetRows[item];
            long double value = 0.0L;
            for (int mode = 0; mode < correction.rank; ++mode) {
                value += static_cast<long double>(
                    transpose
                        ? subdomain.reducedInteriorInterface[
                            static_cast<std::size_t>(
                                mode * subdomain.localInterfaceDofs
                                + gamma)]
                        : subdomain.reducedInterfaceInterior[
                            static_cast<std::size_t>(
                                gamma * correction.rank + mode)])
                    * eliminated[static_cast<std::size_t>(mode)];
            }
            output[static_cast<std::size_t>(targetRow)] -=
                static_cast<double>(value);
        }
    }
}

void PatchTransferOperator::applyCrossBlock(
    const std::vector<double>& input,
    std::vector<double>& output,
    bool transpose) const
{
    const std::size_t expected = transpose
        ? patch_.target.size() : patch_.source.size();
    if (input.size() != expected) {
        throw std::runtime_error(
            "[Optimal port] S_to input size mismatch.");
    }
    output.assign(
        transpose ? patch_.source.size() : patch_.target.size(), 0.0);
    for (const PatchAlgebra::Entry& entry :
         algebra_->crossEntries) {
        if (transpose) {
            output[static_cast<std::size_t>(entry.column)] +=
                entry.value * input[static_cast<std::size_t>(entry.row)];
        } else {
            output[static_cast<std::size_t>(entry.row)] +=
                entry.value * input[static_cast<std::size_t>(entry.column)];
        }
    }
    const local::Model& model = schur_.model();
    for (const PatchAlgebra::LocalCorrection& correction :
         algebra_->localCorrections) {
        if (correction.sourceRows.empty()) continue;
        const local::SubdomainModel& subdomain =
            model.subdomains[correction.subdomainSlot];
        std::vector<double> eliminated(
            static_cast<std::size_t>(correction.rank), 0.0);
        if (transpose) {
            for (std::size_t item = 0;
                 item < correction.targetRows.size(); ++item) {
                const int gamma = correction.targetLocalDofs[item];
                const int targetRow = correction.targetRows[item];
                for (int mode = 0; mode < correction.rank; ++mode) {
                    eliminated[static_cast<std::size_t>(mode)] +=
                        subdomain.reducedInterfaceInterior[
                            static_cast<std::size_t>(
                                gamma * correction.rank + mode)]
                        * input[static_cast<std::size_t>(targetRow)];
                }
            }
        } else {
            for (std::size_t item = 0;
                 item < correction.sourceRows.size(); ++item) {
                const int gamma = correction.sourceLocalDofs[item];
                const int sourceRow = correction.sourceRows[item];
                for (int mode = 0; mode < correction.rank; ++mode) {
                    eliminated[static_cast<std::size_t>(mode)] +=
                        subdomain.reducedInteriorInterface[
                            static_cast<std::size_t>(
                                mode * subdomain.localInterfaceDofs
                                + gamma)]
                        * input[static_cast<std::size_t>(sourceRow)];
                }
            }
        }
        schur_.solveReducedInterior(
            correction.subdomainSlot, eliminated);
        if (transpose) {
            for (std::size_t item = 0;
                 item < correction.sourceRows.size(); ++item) {
                const int gamma = correction.sourceLocalDofs[item];
                const int sourceRow = correction.sourceRows[item];
                long double value = 0.0L;
                for (int mode = 0; mode < correction.rank; ++mode) {
                    value += static_cast<long double>(
                        subdomain.reducedInteriorInterface[
                            static_cast<std::size_t>(
                                mode * subdomain.localInterfaceDofs
                                + gamma)])
                        * eliminated[static_cast<std::size_t>(mode)];
                }
                output[static_cast<std::size_t>(sourceRow)] -=
                    static_cast<double>(value);
            }
        } else {
            for (std::size_t item = 0;
                 item < correction.targetRows.size(); ++item) {
                const int gamma = correction.targetLocalDofs[item];
                const int targetRow = correction.targetRows[item];
                long double value = 0.0L;
                for (int mode = 0; mode < correction.rank; ++mode) {
                    value += static_cast<long double>(
                        subdomain.reducedInterfaceInterior[
                            static_cast<std::size_t>(
                                gamma * correction.rank + mode)])
                        * eliminated[static_cast<std::size_t>(mode)];
                }
                output[static_cast<std::size_t>(targetRow)] -=
                    static_cast<double>(value);
            }
        }
    }
}

bool PatchTransferOperator::solveFgmres(
    const std::vector<double>& rightHandSide,
    bool transpose,
    std::vector<double>& solution,
    int& iterations,
    double& relativeResidual)
{
    const int rows = static_cast<int>(rightHandSide.size());
    // A moderately wide restart is important for the highly conditioned
    // target Schur blocks.  Workspace remains O(restart * target_rows) and is
    // released after every solve/port.
    const int restart = std::min(100, rows);
    const double rightNorm = norm(rightHandSide);
    solution.assign(rightHandSide.size(), 0.0);
    std::vector<double> residual = rightHandSide;
    relativeResidual = rightNorm > 0.0 ? 1.0 : 0.0;
    iterations = 0;
    while (iterations < options_.maximumIterations
           && relativeResidual > options_.relativeTolerance) {
        const int columns = std::min(
            restart, options_.maximumIterations - iterations);
        const double beta = norm(residual);
        if (!(beta > 0.0)) {
            relativeResidual = 0.0;
            break;
        }
        std::vector<double> basis(
            static_cast<std::size_t>((columns + 1) * rows), 0.0);
        std::vector<double> preconditionedBasis(
            static_cast<std::size_t>(columns * rows), 0.0);
        std::vector<double> hessenberg(
            static_cast<std::size_t>((columns + 1) * columns), 0.0);
        std::vector<double> cosine(
            static_cast<std::size_t>(columns), 0.0);
        std::vector<double> sine(
            static_cast<std::size_t>(columns), 0.0);
        std::vector<double> projected(
            static_cast<std::size_t>(columns + 1), 0.0);
        projected[0] = beta;
        for (int row = 0; row < rows; ++row) {
            basis[static_cast<std::size_t>(row)] =
                residual[static_cast<std::size_t>(row)] / beta;
        }

        int accepted = 0;
        for (int column = 0; column < columns; ++column) {
            const double* vector = basis.data()
                + static_cast<std::size_t>(column * rows);
            double* z = preconditionedBasis.data()
                + static_cast<std::size_t>(column * rows);
            for (int row = 0; row < rows; ++row) {
                z[row] = vector[row]
                    / preconditionerDiagonal_[static_cast<std::size_t>(row)];
            }
            std::vector<double> input(
                z, z + static_cast<std::ptrdiff_t>(rows));
            std::vector<double> image;
            applyTargetBlock(input, image, transpose);
            for (int pass = 0; pass < 2; ++pass) {
                for (int prior = 0; prior <= column; ++prior) {
                    const double* priorVector = basis.data()
                        + static_cast<std::size_t>(prior * rows);
                    const double coefficient =
                        dot(priorVector, image.data(), rows);
                    hessenberg[static_cast<std::size_t>(
                        column * (columns + 1) + prior)] += coefficient;
                    for (int row = 0; row < rows; ++row) {
                        image[static_cast<std::size_t>(row)] -=
                            coefficient * priorVector[row];
                    }
                }
            }
            const double nextNorm = norm(image);
            hessenberg[static_cast<std::size_t>(
                column * (columns + 1) + column + 1)] = nextNorm;
            if (nextNorm > 0.0) {
                double* nextVector = basis.data()
                    + static_cast<std::size_t>((column + 1) * rows);
                for (int row = 0; row < rows; ++row) {
                    nextVector[row] =
                        image[static_cast<std::size_t>(row)] / nextNorm;
                }
            }
            for (int prior = 0; prior < column; ++prior) {
                double& upper = hessenberg[static_cast<std::size_t>(
                    column * (columns + 1) + prior)];
                double& lower = hessenberg[static_cast<std::size_t>(
                    column * (columns + 1) + prior + 1)];
                const double rotated =
                    cosine[static_cast<std::size_t>(prior)] * upper
                    + sine[static_cast<std::size_t>(prior)] * lower;
                lower =
                    -sine[static_cast<std::size_t>(prior)] * upper
                    + cosine[static_cast<std::size_t>(prior)] * lower;
                upper = rotated;
            }
            double& diagonal = hessenberg[static_cast<std::size_t>(
                column * (columns + 1) + column)];
            double& subdiagonal = hessenberg[static_cast<std::size_t>(
                column * (columns + 1) + column + 1)];
            const double rotationNorm = std::hypot(diagonal, subdiagonal);
            if (!(rotationNorm > 0.0)) {
                accepted = column;
                break;
            }
            cosine[static_cast<std::size_t>(column)] =
                diagonal / rotationNorm;
            sine[static_cast<std::size_t>(column)] =
                subdiagonal / rotationNorm;
            diagonal = rotationNorm;
            subdiagonal = 0.0;
            projected[static_cast<std::size_t>(column + 1)] =
                -sine[static_cast<std::size_t>(column)]
                * projected[static_cast<std::size_t>(column)];
            projected[static_cast<std::size_t>(column)] *=
                cosine[static_cast<std::size_t>(column)];
            accepted = column + 1;
            ++iterations;
            relativeResidual =
                std::abs(projected[static_cast<std::size_t>(column + 1)])
                / rightNorm;
            if (relativeResidual <= options_.relativeTolerance
                || !(nextNorm > 0.0)) {
                break;
            }
        }
        if (accepted == 0) break;
        std::vector<double> coefficients(
            static_cast<std::size_t>(accepted), 0.0);
        for (int row = accepted - 1; row >= 0; --row) {
            double value = projected[static_cast<std::size_t>(row)];
            for (int column = row + 1; column < accepted; ++column) {
                value -= hessenberg[static_cast<std::size_t>(
                    column * (columns + 1) + row)]
                    * coefficients[static_cast<std::size_t>(column)];
            }
            const double diagonal = hessenberg[static_cast<std::size_t>(
                row * (columns + 1) + row)];
            if (!(std::abs(diagonal)
                    > std::numeric_limits<double>::min())) {
                break;
            }
            coefficients[static_cast<std::size_t>(row)] =
                value / diagonal;
        }
        for (int column = 0; column < accepted; ++column) {
            const double coefficient =
                coefficients[static_cast<std::size_t>(column)];
            const double* z = preconditionedBasis.data()
                + static_cast<std::size_t>(column * rows);
            for (int row = 0; row < rows; ++row) {
                solution[static_cast<std::size_t>(row)] +=
                    coefficient * z[row];
            }
        }
        std::vector<double> image;
        applyTargetBlock(solution, image, transpose);
        for (int row = 0; row < rows; ++row) {
            residual[static_cast<std::size_t>(row)] =
                rightHandSide[static_cast<std::size_t>(row)]
                - image[static_cast<std::size_t>(row)];
        }
        relativeResidual = norm(residual) / rightNorm;
    }
    return relativeResidual <= options_.relativeTolerance;
}

void PatchTransferOperator::solveSparseA(
    const std::vector<double>& rightHandSide,
    bool transpose,
    std::vector<double>& solution) const
{
    if (!sparseFactor_) {
        throw std::runtime_error(
            "[Optimal port] Sparse target factor is missing.");
    }
    sparseFactor_->solve(rightHandSide, transpose, solution);
}

void PatchTransferOperator::solveQ(
    std::vector<double>& rightHandSide,
    bool transpose) const
{
    if (qSparseFactor_) {
        std::vector<double> solution;
        qSparseFactor_->solve(
            rightHandSide, transpose, solution);
        rightHandSide = std::move(solution);
    } else {
        // The symmetric LLT/LDLT factor is identical for Q and Q^T.
        local::solveDenseSymmetric(qFactor_, rightHandSide);
    }
}

double PatchTransferOperator::sparseAResidual(
    const std::vector<double>& rightHandSide,
    const std::vector<double>& solution,
    bool transpose) const
{
    std::vector<double> image(rightHandSide.size(), 0.0);
    for (const PatchAlgebra::Entry& entry :
         algebra_->targetEntries) {
        if (transpose) {
            image[static_cast<std::size_t>(entry.column)] +=
                entry.value
                * solution[static_cast<std::size_t>(entry.row)];
        } else {
            image[static_cast<std::size_t>(entry.row)] +=
                entry.value
                * solution[static_cast<std::size_t>(entry.column)];
        }
    }
    for (std::size_t row = 0; row < image.size(); ++row) {
        image[row] = rightHandSide[row] - image[row];
    }
    return stableNorm(image)
        / std::max(std::numeric_limits<double>::epsilon(),
            stableNorm(rightHandSide));
}

double PatchTransferOperator::qResidual(
    const std::vector<double>& rightHandSide,
    const std::vector<double>& solution,
    bool transpose,
    std::vector<double>* outputResidual) const
{
    const int rank = statistics_.qDimension;
    if (rank <= 0
        || rightHandSide.size() != static_cast<std::size_t>(rank)
        || solution.size() != static_cast<std::size_t>(rank)
        || woodburyQ_.size()
            != static_cast<std::size_t>(rank * rank)) {
        return std::numeric_limits<double>::infinity();
    }
    std::vector<double> residual = rightHandSide;
    for (int row = 0; row < rank; ++row) {
        long double value = 0.0L;
        for (int column = 0; column < rank; ++column) {
            value += static_cast<long double>(
                woodburyQ_[static_cast<std::size_t>(
                    (transpose ? column * rank + row
                               : row * rank + column))])
                * solution[static_cast<std::size_t>(column)];
        }
        residual[static_cast<std::size_t>(row)] -=
            static_cast<double>(value);
    }
    const double relative = stableNorm(residual)
        / std::max(std::numeric_limits<double>::epsilon(),
            stableNorm(rightHandSide));
    if (outputResidual != nullptr) {
        *outputResidual = std::move(residual);
    }
    return relative;
}

double PatchTransferOperator::exactTargetResidual(
    const std::vector<double>& rightHandSide,
    const std::vector<double>& solution,
    bool transpose,
    std::vector<double>* residual) const
{
    std::vector<double> image;
    applyTargetBlock(solution, image, transpose);
    for (std::size_t row = 0; row < image.size(); ++row) {
        image[row] = rightHandSide[row] - image[row];
    }
    const double relative = stableNorm(image)
        / std::max(std::numeric_limits<double>::epsilon(),
            stableNorm(rightHandSide));
    if (residual != nullptr) {
        *residual = std::move(image);
    }
    return relative;
}

void PatchTransferOperator::solveWoodburyOnce(
    const std::vector<double>& rightHandSide,
    bool transpose,
    std::vector<double>& solution)
{
    std::vector<double> aInverse;
    solveSparseA(rightHandSide, transpose, aInverse);
    statistics_.aTtSolveRelativeResidual = std::max(
        statistics_.aTtSolveRelativeResidual,
        sparseAResidual(
            rightHandSide, aInverse, transpose));

    const int rows = targetRows();
    const int rank = statistics_.reducedCorrectionRank;
    std::vector<double> reduced(
        static_cast<std::size_t>(rank), 0.0);
    const std::vector<double>& projection =
        transpose ? woodburyU_ : woodburyV_;
    for (int mode = 0; mode < rank; ++mode) {
        reduced[static_cast<std::size_t>(mode)] = dot(
            projection.data()
                + static_cast<std::size_t>(mode * rows),
            aInverse.data(), rows);
    }
    const std::vector<double> qRightHandSide = reduced;
    solveQ(reduced, transpose);
    const double qRelativeResidual = qResidual(
        qRightHandSide, reduced, transpose);
    statistics_.qSolvePreRefinementResidual = std::max(
        statistics_.qSolvePreRefinementResidual,
        qRelativeResidual);
    statistics_.qSolveRelativeResidual = std::max(
        statistics_.qSolveRelativeResidual,
        qRelativeResidual);

    const std::vector<double>& lift =
        transpose && !woodburyWTranspose_.empty()
            ? woodburyWTranspose_ : woodburyW_;
    std::vector<double> correction(
        static_cast<std::size_t>(rows), 0.0);
    for (int mode = 0; mode < rank; ++mode) {
        const double coefficient =
            reduced[static_cast<std::size_t>(mode)];
        const double* column = lift.data()
            + static_cast<std::size_t>(mode * rows);
        for (int row = 0; row < rows; ++row) {
            correction[static_cast<std::size_t>(row)] +=
                coefficient * column[row];
        }
    }
    const double aInverseNorm = stableNorm(aInverse);
    const double correctionNorm = stableNorm(correction);
    solution = std::move(aInverse);
    for (int row = 0; row < rows; ++row) {
        solution[static_cast<std::size_t>(row)] +=
            correction[static_cast<std::size_t>(row)];
    }
    statistics_.woodburyCancellationFactor = std::max(
        statistics_.woodburyCancellationFactor,
        (aInverseNorm + correctionNorm)
        / std::max(std::numeric_limits<double>::epsilon(),
            stableNorm(solution)));
}

void PatchTransferOperator::solveWoodburyMultipleOnce(
    const std::vector<double>& rightHandSides,
    int rightHandSideCount, bool transpose,
    std::vector<double>& solutions)
{
    const int rows = targetRows();
    const int rank = statistics_.reducedCorrectionRank;
    if (!sparseFactor_ || rightHandSideCount <= 0
        || rightHandSides.size() != static_cast<std::size_t>(
            rows * rightHandSideCount)) {
        throw std::runtime_error(
            "[Residual Krylov] Invalid Woodbury multi-RHS block.");
    }
    std::vector<double> aInverse;
    sparseFactor_->solveMultiple(
        rightHandSides, rightHandSideCount, transpose, aInverse);
    const std::vector<double>& projection =
        transpose ? woodburyU_ : woodburyV_;
    std::vector<double> reduced(static_cast<std::size_t>(
        rank * rightHandSideCount), 0.0);
    for (int mode = 0; mode < rank; ++mode) {
        const double* projectionColumn = projection.data()
            + static_cast<std::size_t>(mode * rows);
        for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
            reduced[static_cast<std::size_t>(
                mode * rightHandSideCount + rhs)] =
                dot(projectionColumn,
                    aInverse.data()
                        + static_cast<std::size_t>(rhs * rows),
                    rows);
        }
    }
    const std::vector<double> qRightHandSides = reduced;
    if (qSparseFactor_) {
        std::vector<double> packed(static_cast<std::size_t>(
            rank * rightHandSideCount), 0.0);
        for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
            for (int mode = 0; mode < rank; ++mode) {
                packed[static_cast<std::size_t>(
                    rhs * rank + mode)] =
                    reduced[static_cast<std::size_t>(
                        mode * rightHandSideCount + rhs)];
            }
        }
        std::vector<double> packedSolution;
        qSparseFactor_->solveMultiple(
            packed, rightHandSideCount,
            transpose, packedSolution);
        for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
            for (int mode = 0; mode < rank; ++mode) {
                reduced[static_cast<std::size_t>(
                    mode * rightHandSideCount + rhs)] =
                    packedSolution[static_cast<std::size_t>(
                        rhs * rank + mode)];
            }
        }
    } else {
        local::solveDenseSymmetricMultiple(
            qFactor_, reduced, rightHandSideCount);
    }
    const std::vector<double>& lift =
        transpose && !woodburyWTranspose_.empty()
            ? woodburyWTranspose_ : woodburyW_;
    solutions = aInverse;
    for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
        for (int mode = 0; mode < rank; ++mode) {
            const double coefficient =
                reduced[static_cast<std::size_t>(
                    mode * rightHandSideCount + rhs)];
            const double* liftColumn = lift.data()
                + static_cast<std::size_t>(mode * rows);
            for (int row = 0; row < rows; ++row) {
                solutions[static_cast<std::size_t>(
                    rhs * rows + row)] +=
                    coefficient * liftColumn[row];
            }
        }
    }
    for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
        std::vector<double> original(
            rightHandSides.begin()
                + static_cast<std::ptrdiff_t>(rhs * rows),
            rightHandSides.begin()
                + static_cast<std::ptrdiff_t>((rhs + 1) * rows));
        std::vector<double> aColumn(
            aInverse.begin()
                + static_cast<std::ptrdiff_t>(rhs * rows),
            aInverse.begin()
                + static_cast<std::ptrdiff_t>((rhs + 1) * rows));
        statistics_.aTtSolveRelativeResidual = std::max(
            statistics_.aTtSolveRelativeResidual,
            sparseAResidual(
                original, aColumn, transpose));
        std::vector<double> qRhs(
            static_cast<std::size_t>(rank), 0.0);
        std::vector<double> qSolution(
            static_cast<std::size_t>(rank), 0.0);
        for (int mode = 0; mode < rank; ++mode) {
            qRhs[static_cast<std::size_t>(mode)] =
                qRightHandSides[static_cast<std::size_t>(
                    mode * rightHandSideCount + rhs)];
            qSolution[static_cast<std::size_t>(mode)] =
                reduced[static_cast<std::size_t>(
                    mode * rightHandSideCount + rhs)];
        }
        const double qRelative = qResidual(
            qRhs, qSolution, transpose);
        statistics_.qSolvePreRefinementResidual = std::max(
            statistics_.qSolvePreRefinementResidual,
            qRelative);
        statistics_.qSolveRelativeResidual = std::max(
            statistics_.qSolveRelativeResidual,
            qRelative);
    }
}

void PatchTransferOperator::solveTargetMultiple(
    const std::vector<double>& rightHandSides,
    int rightHandSideCount, bool transpose,
    std::vector<double>& solutions)
{
    if (!woodbury_) {
        throw std::runtime_error(
            "[Residual Krylov] Block target solve currently requires "
            "woodbury-exact.");
    }
    const auto start = Clock::now();
    const int rows = targetRows();
    solveWoodburyMultipleOnce(
        rightHandSides, rightHandSideCount,
        transpose, solutions);
    std::vector<double> residuals(rightHandSides.size(), 0.0);
    std::vector<double> relative(
        static_cast<std::size_t>(rightHandSideCount), 0.0);
    auto updateResiduals = [&]() {
        double maximum = 0.0;
        for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
            std::vector<double> original(
                rightHandSides.begin()
                    + static_cast<std::ptrdiff_t>(rhs * rows),
                rightHandSides.begin()
                    + static_cast<std::ptrdiff_t>((rhs + 1) * rows));
            std::vector<double> solution(
                solutions.begin()
                    + static_cast<std::ptrdiff_t>(rhs * rows),
                solutions.begin()
                    + static_cast<std::ptrdiff_t>((rhs + 1) * rows));
            std::vector<double> residual;
            relative[static_cast<std::size_t>(rhs)] =
                exactTargetResidual(
                    original, solution, transpose, &residual);
            std::copy(
                residual.begin(), residual.end(),
                residuals.begin()
                    + static_cast<std::ptrdiff_t>(rhs * rows));
            maximum = std::max(
                maximum,
                relative[static_cast<std::size_t>(rhs)]);
        }
        return maximum;
    };
    double maximumResidual = updateResiduals();
    statistics_.woodburyPreRefinementResidual = std::max(
        statistics_.woodburyPreRefinementResidual,
        maximumResidual);
    int refinements = 0;
    if (maximumResidual > options_.refinementTolerance
        && options_.refinementMaximumIterations > 0) {
        ++statistics_.refinementTriggeredSolveCalls;
    }
    const auto refinementStart = Clock::now();
    while (maximumResidual > options_.refinementTolerance
        && refinements
            < options_.refinementMaximumIterations) {
        std::vector<double> corrections;
        solveWoodburyMultipleOnce(
            residuals, rightHandSideCount,
            transpose, corrections);
        for (std::size_t row = 0;
             row < solutions.size(); ++row) {
            solutions[row] += corrections[row];
        }
        ++refinements;
        maximumResidual = updateResiduals();
    }
    statistics_.refinementSeconds += elapsed(refinementStart);
    statistics_.refinementIterations = std::max(
        statistics_.refinementIterations, refinements);
    statistics_.woodburyPostRefinementResidual = std::max(
        statistics_.woodburyPostRefinementResidual,
        maximumResidual);
    statistics_.refinementConverged =
        statistics_.refinementConverged
        && maximumResidual <= options_.refinementTolerance;
    const double seconds = elapsed(start);
    ++statistics_.solveCalls;
    statistics_.solveRightHandSides += rightHandSideCount;
    statistics_.totalIterations += refinements;
    statistics_.maximumIterations = std::max(
        statistics_.maximumIterations, refinements);
    statistics_.totalSolveSeconds += seconds;
    statistics_.maximumSolveSeconds = std::max(
        statistics_.maximumSolveSeconds, seconds);
    statistics_.finalRelativeResidual = maximumResidual;
    statistics_.maximumRelativeResidual = std::max(
        statistics_.maximumRelativeResidual,
        maximumResidual);
    statistics_.status =
        maximumResidual <= 1.0e-9
        ? "success" : "direct_target_residual_failed";
}

void PatchTransferOperator::solveTarget(
    std::vector<double>& rightHandSide,
    bool transpose)
{
    const auto start = Clock::now();
    const std::vector<double> originalRightHandSide = rightHandSide;
    const double rhsNorm = direct_
        ? stableNorm(rightHandSide) : norm(rightHandSide);
    int iterations = 0;
    double relativeResidual = 0.0;
    bool directResidualFailed = false;
    bool exactResidualAvailable = false;
    if (assembledDense_ && sparseFactor_) {
        const std::vector<double> rhs = rightHandSide;
        sparseFactor_->solve(rhs, transpose, rightHandSide);
    } else if (assembledDense_) {
        local::solveDenseSymmetric(factor_, rightHandSide);
    } else if (woodbury_) {
        solveWoodburyOnce(
            originalRightHandSide, transpose, rightHandSide);
        std::vector<double> residual;
        const double preRefinementResidual =
            exactTargetResidual(
                originalRightHandSide, rightHandSide,
                transpose, &residual);
        std::vector<double> residualHistory(4,
            std::numeric_limits<double>::quiet_NaN());
        residualHistory[0] = preRefinementResidual;
        double maximumCorrectionRelativeNorm = 0.0;
        int refinementIterations = 0;
        const auto refinementStart = Clock::now();
        relativeResidual = preRefinementResidual;
        if (relativeResidual > options_.refinementTolerance
            && options_.refinementMaximumIterations > 0) {
            ++statistics_.refinementTriggeredSolveCalls;
        }
        while (relativeResidual > options_.refinementTolerance
               && refinementIterations
                    < options_.refinementMaximumIterations) {
            std::vector<double> correction;
            solveWoodburyOnce(residual, transpose, correction);
            const double solutionNorm = stableNorm(rightHandSide);
            maximumCorrectionRelativeNorm = std::max(
                maximumCorrectionRelativeNorm,
                stableNorm(correction)
                / std::max(
                    std::numeric_limits<double>::epsilon(),
                    solutionNorm));
            for (std::size_t row = 0;
                 row < rightHandSide.size(); ++row) {
                rightHandSide[row] += correction[row];
            }
            ++refinementIterations;
            relativeResidual = exactTargetResidual(
                originalRightHandSide, rightHandSide,
                transpose, &residual);
            residualHistory[static_cast<std::size_t>(
                refinementIterations)] = relativeResidual;
        }
        statistics_.refinementSeconds +=
            elapsed(refinementStart);
        statistics_.woodburyPostRefinementResidual = std::max(
            statistics_.woodburyPostRefinementResidual,
            relativeResidual);
        statistics_.refinementIterations = std::max(
            statistics_.refinementIterations,
            refinementIterations);
        statistics_.refinementCorrectionRelativeNorm = std::max(
            statistics_.refinementCorrectionRelativeNorm,
            maximumCorrectionRelativeNorm);
        statistics_.refinementConverged =
            statistics_.refinementConverged
            && relativeResidual <= options_.refinementTolerance;
        if (preRefinementResidual
                >= statistics_.woodburyPreRefinementResidual) {
            statistics_.woodburyPreRefinementResidual =
                preRefinementResidual;
            statistics_.refinementResidual0 = residualHistory[0];
            statistics_.refinementResidual1 = residualHistory[1];
            statistics_.refinementResidual2 = residualHistory[2];
            statistics_.refinementResidual3 = residualHistory[3];
            statistics_.refinementReductionFactor =
                relativeResidual
                / std::max(std::numeric_limits<double>::epsilon(),
                    preRefinementResidual);
        }
        exactResidualAvailable = true;
    } else if (rhsNorm > 0.0 && forceFgmres_) {
        const std::vector<double> rhs = rightHandSide;
        std::vector<double> solution;
        if (!solveFgmres(
                rhs, transpose, solution, iterations,
                relativeResidual)) {
            statistics_.status = "fgmres_maximum_iterations";
            throw std::runtime_error(
                "[Optimal port] S_tt fallback FGMRES did not converge.");
        }
        rightHandSide = std::move(solution);
    } else if (rhsNorm > 0.0) {
        const std::vector<double> rhs = rightHandSide;
        std::vector<double> solution(rhs.size(), 0.0);
        std::vector<double> residual = rhs;
        std::vector<double> preconditioned(rhs.size(), 0.0);
        for (std::size_t row = 0; row < rhs.size(); ++row) {
            preconditioned[row] =
                residual[row] / preconditionerDiagonal_[row];
        }
        std::vector<double> direction = preconditioned;
        double residualPreconditioned = dot(residual, preconditioned);
        relativeResidual = norm(residual) / rhsNorm;
        bool pcgBreakdown = false;
        for (; iterations < options_.maximumIterations
               && relativeResidual > options_.relativeTolerance;
             ++iterations) {
            std::vector<double> image;
            applyTargetBlock(direction, image, transpose);
            const double curvature = dot(direction, image);
            const double scale =
                std::max(std::numeric_limits<double>::min(),
                    norm(direction) * norm(image));
            if (!(curvature > 64.0
                    * std::numeric_limits<double>::epsilon() * scale)
                || !std::isfinite(curvature)) {
                pcgBreakdown = true;
                break;
            }
            const double alpha = residualPreconditioned / curvature;
            for (std::size_t row = 0; row < rhs.size(); ++row) {
                solution[row] += alpha * direction[row];
                residual[row] -= alpha * image[row];
            }
            relativeResidual = norm(residual) / rhsNorm;
            if (relativeResidual <= options_.relativeTolerance) {
                ++iterations;
                break;
            }
            for (std::size_t row = 0; row < rhs.size(); ++row) {
                preconditioned[row] =
                    residual[row] / preconditionerDiagonal_[row];
            }
            const double next = dot(residual, preconditioned);
            const double beta = next / residualPreconditioned;
            for (std::size_t row = 0; row < rhs.size(); ++row) {
                direction[row] =
                    preconditioned[row] + beta * direction[row];
            }
            residualPreconditioned = next;
        }
        if (pcgBreakdown
            || relativeResidual > options_.relativeTolerance) {
            const int pcgIterations = iterations;
            int fgmresIterations = 0;
            statistics_.fallbackReason = pcgBreakdown
                ? "pcg_non_positive_curvature"
                : "pcg_maximum_iterations";
            statistics_.solver = "matrix-free-pcg->fgmres";
            statistics_.actualSolver = "matrix-free-fgmres";
            if (!solveFgmres(
                    rhs, transpose, solution, fgmresIterations,
                    relativeResidual)) {
                statistics_.status = "fgmres_maximum_iterations";
                throw std::runtime_error(
                    "[Optimal port] S_tt PCG fallback FGMRES did not converge.");
            }
            iterations = pcgIterations + fgmresIterations;
        }
        rightHandSide = std::move(solution);
    }
    if (rhsNorm > 0.0 && !exactResidualAvailable) {
        std::vector<double> image;
        applyTargetBlock(rightHandSide, image, transpose);
        for (std::size_t row = 0; row < image.size(); ++row) {
            image[row] -= originalRightHandSide[row];
        }
        relativeResidual = (direct_
            ? stableNorm(image) : norm(image)) / rhsNorm;
    }
    if (rhsNorm > 0.0
        && (assembledDense_ || woodbury_)
        && relativeResidual > std::max(
            1.0e-9, 10.0 * options_.relativeTolerance)) {
        statistics_.status = "direct_target_residual_failed";
        statistics_.finalRelativeResidual = relativeResidual;
        statistics_.maximumRelativeResidual = std::max(
            statistics_.maximumRelativeResidual,
            relativeResidual);
        directResidualFailed = true;
    }
    const double seconds = elapsed(start);
    ++statistics_.solveCalls;
    ++statistics_.solveRightHandSides;
    statistics_.totalIterations += iterations;
    statistics_.maximumIterations =
        std::max(statistics_.maximumIterations, iterations);
    statistics_.totalSolveSeconds += seconds;
    statistics_.maximumSolveSeconds =
        std::max(statistics_.maximumSolveSeconds, seconds);
    statistics_.finalRelativeResidual = relativeResidual;
    statistics_.maximumRelativeResidual = std::max(
        statistics_.maximumRelativeResidual, relativeResidual);
    statistics_.status = directResidualFailed
        ? "direct_target_residual_failed"
        : (statistics_.fallbackReason.empty()
            ? "success"
            : (assembledDense_ || woodbury_
                ? "success_reported_fallback"
                : "success_fallback_fgmres"));
}

void PatchTransferOperator::apply(
    const std::vector<double>& source,
    std::vector<double>& target)
{
    const auto start = Clock::now();
    if (source.size() != patch_.source.size()) {
        throw std::runtime_error(
            "[Optimal port] Transfer source size mismatch.");
    }
    target.assign(patch_.target.size(), 0.0);
    if (source.empty()) return;
    formTraceRightHandSide(source, target);
    solveTarget(target, false);
    statistics_.transferApplySeconds += elapsed(start);
}

void PatchTransferOperator::applyTranspose(
    const std::vector<double>& target,
    std::vector<double>& source)
{
    const auto start = Clock::now();
    if (target.size() != patch_.target.size()) {
        throw std::runtime_error(
            "[Optimal port] Transfer transpose target size mismatch.");
    }
    source.assign(patch_.source.size(), 0.0);
    if (source.empty()) return;
    std::vector<double> solved = target;
    solveTarget(solved, true);
    applyTraceRightHandSideTranspose(solved, source);
    statistics_.transposeApplySeconds += elapsed(start);
}

void PatchTransferOperator::applyMultiple(
    const std::vector<double>& sources,
    int rightHandSideCount,
    std::vector<double>& targets)
{
    const auto start = Clock::now();
    const int sourceRowsValue = sourceRows();
    const int targetRowsValue = targetRows();
    if (rightHandSideCount <= 0
        || sources.size() != static_cast<std::size_t>(
            sourceRowsValue * rightHandSideCount)) {
        throw std::runtime_error(
            "[Optimal port] Transfer multi-source dimensions are invalid.");
    }
    std::vector<double> rightHandSides(static_cast<std::size_t>(
        targetRowsValue * rightHandSideCount), 0.0);
    for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
        std::vector<double> source(
            sources.begin() + static_cast<std::ptrdiff_t>(
                rhs * sourceRowsValue),
            sources.begin() + static_cast<std::ptrdiff_t>(
                (rhs + 1) * sourceRowsValue));
        std::vector<double> target;
        formTraceRightHandSide(source, target);
        std::copy(
            target.begin(), target.end(),
            rightHandSides.begin() + static_cast<std::ptrdiff_t>(
                rhs * targetRowsValue));
    }
    solveTargetMultiple(
        rightHandSides, rightHandSideCount, false, targets);
    statistics_.transferApplySeconds += elapsed(start);
}

void PatchTransferOperator::applyTransposeMultiple(
    const std::vector<double>& targets,
    int rightHandSideCount,
    std::vector<double>& sources)
{
    const auto start = Clock::now();
    const int sourceRowsValue = sourceRows();
    const int targetRowsValue = targetRows();
    if (rightHandSideCount <= 0
        || targets.size() != static_cast<std::size_t>(
            targetRowsValue * rightHandSideCount)) {
        throw std::runtime_error(
            "[Optimal port] Transfer multi-target dimensions are invalid.");
    }
    std::vector<double> solved;
    solveTargetMultiple(
        targets, rightHandSideCount, true, solved);
    sources.assign(static_cast<std::size_t>(
        sourceRowsValue * rightHandSideCount), 0.0);
    for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
        std::vector<double> target(
            solved.begin() + static_cast<std::ptrdiff_t>(
                rhs * targetRowsValue),
            solved.begin() + static_cast<std::ptrdiff_t>(
                (rhs + 1) * targetRowsValue));
        std::vector<double> source;
        applyTraceRightHandSideTranspose(target, source);
        std::copy(
            source.begin(), source.end(),
            sources.begin() + static_cast<std::ptrdiff_t>(
                rhs * sourceRowsValue));
    }
    statistics_.transposeApplySeconds += elapsed(start);
}

void PatchTransferOperator::solveTargetResponse(
    const std::vector<double>& rightHandSide,
    std::vector<double>& target)
{
    if (rightHandSide.size() != patch_.target.size()) {
        throw std::runtime_error(
            "[Optimal port] Particular target RHS size mismatch.");
    }
    target = rightHandSide;
    solveTarget(target, false);
}

void PatchTransferOperator::solveTargetResponses(
    const std::vector<double>& rightHandSides,
    int rightHandSideCount,
    std::vector<double>& targets)
{
    if (rightHandSideCount <= 0
        || rightHandSides.size()
            != patch_.target.size()
                * static_cast<std::size_t>(
                    rightHandSideCount)) {
        throw std::runtime_error(
            "[Residual Krylov] Particular target block dimensions "
            "are invalid.");
    }
    solveTargetMultiple(
        rightHandSides, rightHandSideCount,
        false, targets);
}

void PatchTransferOperator::solveTargetResponsesTranspose(
    const std::vector<double>& rightHandSides,
    int rightHandSideCount,
    std::vector<double>& targets)
{
    if (rightHandSideCount <= 0
        || rightHandSides.size()
            != patch_.target.size()
                * static_cast<std::size_t>(
                    rightHandSideCount)) {
        throw std::runtime_error(
            "[Randomized port] Transpose target block dimensions "
            "are invalid.");
    }
    solveTargetMultiple(
        rightHandSides, rightHandSideCount,
        true, targets);
}

void PatchTransferOperator::solveTargetResponseTranspose(
    const std::vector<double>& rightHandSide,
    std::vector<double>& target)
{
    if (rightHandSide.size() != patch_.target.size()) {
        throw std::runtime_error(
            "[Optimal port] Particular transpose target RHS size mismatch.");
    }
    target = rightHandSide;
    solveTarget(target, true);
}

void PatchTransferOperator::applyTargetAction(
    const std::vector<double>& input, std::vector<double>& output) const
{
    applyTargetBlock(input, output, false);
}

void PatchTransferOperator::applyTargetActionTranspose(
    const std::vector<double>& input, std::vector<double>& output) const
{
    applyTargetBlock(input, output, true);
}

void PatchTransferOperator::formTraceRightHandSide(
    const std::vector<double>& source,
    std::vector<double>& target) const
{
    applyCrossBlock(source, target, false);
    for (double& value : target) value = -value;
}

void PatchTransferOperator::applyTraceRightHandSideTranspose(
    const std::vector<double>& target,
    std::vector<double>& source) const
{
    applyCrossBlock(target, source, true);
    for (double& value : source) value = -value;
}

std::vector<double> PatchTransferOperator::explicitMatrix()
{
    std::vector<double> matrix(
        static_cast<std::size_t>(targetRows() * sourceRows()), 0.0);
    for (int column = 0; column < sourceRows(); ++column) {
        std::vector<double> input(
            static_cast<std::size_t>(sourceRows()), 0.0);
        input[static_cast<std::size_t>(column)] = 1.0;
        std::vector<double> output;
        apply(input, output);
        for (int row = 0; row < targetRows(); ++row) {
            matrix[static_cast<std::size_t>(
                row * sourceRows() + column)] =
                output[static_cast<std::size_t>(row)];
        }
    }
    return matrix;
}

GeneralizedPatchTransferOperator::GeneralizedPatchTransferOperator(
    const ReducedDynamicSchurOperator& schur,
    PortPatch patch,
    GeneralizedTransferSourceBlocks blocks,
    bool includeTrace,
    bool includeInput,
    bool includeBoundary,
    bool includeHistory,
    const PatchTransferOptions& options)
    : schur_(schur),
      blocks_(std::move(blocks)),
      includeTrace_(includeTrace),
      includeInput_(includeInput),
      includeBoundary_(includeBoundary),
      includeHistory_(includeHistory),
      targetSolver_(schur, std::move(patch), options)
{
    if (blocks_.interfaceDofs != schur_.size()
        || blocks_.inputChannels < 0
        || blocks_.boundaryChannels < 0
        || blocks_.historyChannels < 0
        || blocks_.input.size() != static_cast<std::size_t>(
            blocks_.interfaceDofs * blocks_.inputChannels)
        || blocks_.boundary.size() != static_cast<std::size_t>(
            blocks_.interfaceDofs * blocks_.boundaryChannels)
        || blocks_.history.size() != static_cast<std::size_t>(
            blocks_.interfaceDofs * blocks_.historyChannels)) {
        throw std::runtime_error(
            "[Optimal port] Generalized source block dimensions are invalid.");
    }
}

int GeneralizedPatchTransferOperator::sourceRows() const
{
    return traceRows() + inputRows() + boundaryRows() + historyRows();
}

PatchInnerSolverStatistics
GeneralizedPatchTransferOperator::statistics() const
{
    PatchInnerSolverStatistics result = targetSolver_.statistics();
    result.transferApplySeconds += transferApplySeconds_;
    result.transposeApplySeconds += transposeApplySeconds_;
    return result;
}

void GeneralizedPatchTransferOperator::addTraceRightHandSide(
    const double* source,
    std::vector<double>& target) const
{
    std::vector<double> traceSource(
        source, source + targetSolver_.patch().source.size());
    std::vector<double> traceRightHandSide;
    targetSolver_.formTraceRightHandSide(
        traceSource, traceRightHandSide);
    for (std::size_t row = 0; row < target.size(); ++row) {
        target[row] += traceRightHandSide[row];
    }
}

void GeneralizedPatchTransferOperator::addColumns(
    const std::vector<double>& columns,
    int channels,
    const double* coefficients,
    std::vector<double>& target) const
{
    for (int channel = 0; channel < channels; ++channel) {
        const double coefficient =
            coefficients[static_cast<std::size_t>(channel)];
        if (coefficient == 0.0) continue;
        const double* column = columns.data()
            + static_cast<std::size_t>(channel * blocks_.interfaceDofs);
        for (std::size_t row = 0;
             row < targetSolver_.patch().target.size(); ++row) {
            target[row] += coefficient * column[
                static_cast<std::size_t>(
                    targetSolver_.patch().target[row])];
        }
    }
}

void GeneralizedPatchTransferOperator::transposeColumns(
    const std::vector<double>& columns,
    int channels,
    const std::vector<double>& target,
    double* coefficients) const
{
    for (int channel = 0; channel < channels; ++channel) {
        const double* column = columns.data()
            + static_cast<std::size_t>(channel * blocks_.interfaceDofs);
        long double value = 0.0L;
        for (std::size_t row = 0;
             row < targetSolver_.patch().target.size(); ++row) {
            value += static_cast<long double>(column[
                static_cast<std::size_t>(
                    targetSolver_.patch().target[row])]) * target[row];
        }
        coefficients[static_cast<std::size_t>(channel)] =
            static_cast<double>(value);
    }
}

void GeneralizedPatchTransferOperator::apply(
    const std::vector<double>& source,
    std::vector<double>& target)
{
    const auto start = Clock::now();
    if (source.size() != static_cast<std::size_t>(sourceRows())) {
        throw std::runtime_error(
            "[Optimal port] Generalized transfer source size mismatch.");
    }
    target.assign(static_cast<std::size_t>(targetRows()), 0.0);
    std::size_t offset = 0;
    if (includeTrace_) {
        addTraceRightHandSide(source.data() + offset, target);
        offset += static_cast<std::size_t>(traceRows());
    }
    if (includeInput_) {
        addColumns(blocks_.input, blocks_.inputChannels,
                   source.data() + offset, target);
        offset += static_cast<std::size_t>(inputRows());
    }
    if (includeBoundary_) {
        addColumns(blocks_.boundary, blocks_.boundaryChannels,
                   source.data() + offset, target);
        offset += static_cast<std::size_t>(boundaryRows());
    }
    if (includeHistory_) {
        addColumns(blocks_.history, blocks_.historyChannels,
                   source.data() + offset, target);
    }
    targetSolver_.solveTargetResponse(target, target);
    transferApplySeconds_ += elapsed(start);
}

void GeneralizedPatchTransferOperator::applyTranspose(
    const std::vector<double>& target,
    std::vector<double>& source)
{
    const auto start = Clock::now();
    if (target.size() != static_cast<std::size_t>(targetRows())) {
        throw std::runtime_error(
            "[Optimal port] Generalized transfer transpose size mismatch.");
    }
    std::vector<double> solved;
    targetSolver_.solveTargetResponseTranspose(target, solved);
    source.assign(static_cast<std::size_t>(sourceRows()), 0.0);
    std::size_t offset = 0;
    if (includeTrace_) {
        std::vector<double> traceTranspose;
        targetSolver_.applyTraceRightHandSideTranspose(
            solved, traceTranspose);
        for (std::size_t row = 0; row < traceTranspose.size(); ++row) {
            source[offset + row] = traceTranspose[row];
        }
        offset += static_cast<std::size_t>(traceRows());
    }
    if (includeInput_) {
        transposeColumns(blocks_.input, blocks_.inputChannels, solved,
                         source.data() + offset);
        offset += static_cast<std::size_t>(inputRows());
    }
    if (includeBoundary_) {
        transposeColumns(blocks_.boundary, blocks_.boundaryChannels, solved,
                         source.data() + offset);
        offset += static_cast<std::size_t>(boundaryRows());
    }
    if (includeHistory_) {
        transposeColumns(blocks_.history, blocks_.historyChannels, solved,
                         source.data() + offset);
    }
    transposeApplySeconds_ += elapsed(start);
}

void GeneralizedPatchTransferOperator::applyMultiple(
    const std::vector<double>& sources,
    int rightHandSideCount,
    std::vector<double>& targets)
{
    const auto start = Clock::now();
    if (rightHandSideCount <= 0
        || sources.size() != static_cast<std::size_t>(
            sourceRows() * rightHandSideCount)) {
        throw std::runtime_error(
            "[Randomized port] Generalized multi-source dimensions "
            "are invalid.");
    }
    std::vector<double> rightHandSides(static_cast<std::size_t>(
        targetRows() * rightHandSideCount), 0.0);
    for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
        const double* source = sources.data()
            + static_cast<std::size_t>(rhs * sourceRows());
        std::vector<double> target(
            static_cast<std::size_t>(targetRows()), 0.0);
        std::size_t offset = 0;
        if (includeTrace_) {
            addTraceRightHandSide(source + offset, target);
            offset += static_cast<std::size_t>(traceRows());
        }
        if (includeInput_) {
            addColumns(
                blocks_.input, blocks_.inputChannels,
                source + offset, target);
            offset += static_cast<std::size_t>(inputRows());
        }
        if (includeBoundary_) {
            addColumns(
                blocks_.boundary, blocks_.boundaryChannels,
                source + offset, target);
            offset += static_cast<std::size_t>(boundaryRows());
        }
        if (includeHistory_) {
            addColumns(
                blocks_.history, blocks_.historyChannels,
                source + offset, target);
        }
        std::copy(
            target.begin(), target.end(),
            rightHandSides.begin() + static_cast<std::ptrdiff_t>(
                rhs * targetRows()));
    }
    targetSolver_.solveTargetResponses(
        rightHandSides, rightHandSideCount, targets);
    transferApplySeconds_ += elapsed(start);
}

void GeneralizedPatchTransferOperator::applyTransposeMultiple(
    const std::vector<double>& targets,
    int rightHandSideCount,
    std::vector<double>& sources)
{
    const auto start = Clock::now();
    if (rightHandSideCount <= 0
        || targets.size() != static_cast<std::size_t>(
            targetRows() * rightHandSideCount)) {
        throw std::runtime_error(
            "[Randomized port] Generalized multi-target dimensions "
            "are invalid.");
    }
    std::vector<double> solved;
    targetSolver_.solveTargetResponsesTranspose(
        targets, rightHandSideCount, solved);
    sources.assign(static_cast<std::size_t>(
        sourceRows() * rightHandSideCount), 0.0);
    for (int rhs = 0; rhs < rightHandSideCount; ++rhs) {
        std::vector<double> target(
            solved.begin() + static_cast<std::ptrdiff_t>(
                rhs * targetRows()),
            solved.begin() + static_cast<std::ptrdiff_t>(
                (rhs + 1) * targetRows()));
        double* source = sources.data()
            + static_cast<std::size_t>(rhs * sourceRows());
        std::size_t offset = 0;
        if (includeTrace_) {
            std::vector<double> traceTranspose;
            targetSolver_.applyTraceRightHandSideTranspose(
                target, traceTranspose);
            std::copy(
                traceTranspose.begin(), traceTranspose.end(),
                source + offset);
            offset += static_cast<std::size_t>(traceRows());
        }
        if (includeInput_) {
            transposeColumns(
                blocks_.input, blocks_.inputChannels,
                target, source + offset);
            offset += static_cast<std::size_t>(inputRows());
        }
        if (includeBoundary_) {
            transposeColumns(
                blocks_.boundary, blocks_.boundaryChannels,
                target, source + offset);
            offset += static_cast<std::size_t>(boundaryRows());
        }
        if (includeHistory_) {
            transposeColumns(
                blocks_.history, blocks_.historyChannels,
                target, source + offset);
        }
    }
    transposeApplySeconds_ += elapsed(start);
}

void GeneralizedPatchTransferOperator::solveTargetResponse(
    const std::vector<double>& rightHandSide,
    std::vector<double>& target)
{
    targetSolver_.solveTargetResponse(rightHandSide, target);
}

} // namespace mor::transient
