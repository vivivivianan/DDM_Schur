#include "port_core_solver.hpp"

#include "ddm_schur/interface_operator.hpp"
#include "sipg_core.hpp"
#include "linear_solvers.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <system_error>
#include <type_traits>
#include <utility>

namespace mor::local {
namespace {

using Clock = std::chrono::steady_clock;
struct UpperCsr {
    std::vector<int> rowPtr;
    std::vector<int> colInd;
    std::vector<double> values;
};

constexpr std::uint64_t portCoreCacheMagic =
    UINT64_C(0x504f5254434f5245); // "PORTCORE"
constexpr int portCoreCacheVersion = 1;

double elapsed(const Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

template <typename T>
void appendHash(std::uint64_t& hash, const T& value);

template <typename T>
void appendHash(std::uint64_t& hash, const std::vector<T>& values);

const SubdomainModel& templatePayload(const Model& model,
                                      const SubdomainModel& instance)
{
    if (!instance.templateReused) {
        return instance;
    }
    for (const SubdomainModel& candidate : model.subdomains) {
        if (candidate.templateId == instance.templateId
            && !candidate.templateReused) {
            return candidate;
        }
    }
    throw std::runtime_error("[Port core] Reused template payload is missing.");
}

std::uint64_t portCoreOperatorSignature(
    const Model& model,
    const ddm_schur::InterfacePartition& partition)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    appendHash(hash, portCoreCacheVersion);
    appendHash(hash, model.formatVersion);
    appendHash(hash, model.globalDofs);
    appendHash(hash, model.interfaceDofs);
    appendHash(hash, model.totalLocalRank);
    appendHash(hash, model.fingerprints.mesh);
    appendHash(hash, model.fingerprints.system);
    appendHash(hash, model.fingerprints.interfaceOrdering);
    appendHash(hash, model.interfaceGlobalDofs);
    appendHash(hash, model.interfaceEntries.size());
    for (const InterfaceEntry& entry : model.interfaceEntries) {
        appendHash(hash, entry.row);
        appendHash(hash, entry.column);
        appendHash(hash, entry.value);
    }
    appendHash(hash, model.subdomains.size());
    for (const SubdomainModel& local : model.subdomains) {
        const SubdomainModel& data = templatePayload(model, local);
        appendHash(hash, local.subdomain);
        appendHash(hash, local.rank);
        appendHash(hash, local.interfaceIndices);
        appendHash(hash, data.reducedInterior);
        appendHash(hash, data.reducedInteriorInterface);
        appendHash(hash, data.reducedInterfaceInterior);
    }
    appendHash(hash, partition.domains.size());
    for (const ddm_schur::DomainBlocks& domain : partition.domains) {
        appendHash(hash, domain.domainId);
        appendHash(hash, domain.interfaceGlobalDofsByNeighbor.size());
        for (const auto& [neighbor, globalDofs] :
             domain.interfaceGlobalDofsByNeighbor) {
            appendHash(hash, neighbor);
            appendHash(hash, globalDofs);
        }
    }
    return hash;
}

std::size_t packedUpperSize(int size)
{
    if (size <= 0) {
        throw std::runtime_error("[Port core] Invalid packed matrix size.");
    }
    return static_cast<std::size_t>(size)
        * static_cast<std::size_t>(size + 1) / 2;
}

std::size_t packedUpperIndex(int size, int row, int column)
{
    if (row > column) {
        std::swap(row, column);
    }
    if (row < 0 || column >= size) {
        throw std::runtime_error("[Port core] Core matrix index is out of range.");
    }
    const std::size_t rowOffset = static_cast<std::size_t>(row)
        * static_cast<std::size_t>(2 * size - row + 1) / 2;
    return rowOffset + static_cast<std::size_t>(column - row);
}

void addPackedUpper(std::vector<double>& packed,
                    int size,
                    int row,
                    int column,
                    double value)
{
    if (value != 0.0) {
        packed[packedUpperIndex(size, row, column)] += value;
    }
}

UpperCsr buildUpperCsr(int size, std::vector<MatrixEntry> entries)
{
    for (MatrixEntry& entry : entries) {
        if (entry.row > entry.col) {
            std::swap(entry.row, entry.col);
        }
        if (entry.row < 0 || entry.col >= size) {
            throw std::runtime_error(
                "[Port core] Sparse block entry is out of range.");
        }
    }
    std::sort(entries.begin(), entries.end(),
        [](const MatrixEntry& left, const MatrixEntry& right) {
            return left.row < right.row
                || (left.row == right.row && left.col < right.col);
        });

    UpperCsr result;
    result.rowPtr.assign(static_cast<std::size_t>(size + 1), 0);
    result.colInd.reserve(entries.size());
    result.values.reserve(entries.size());
    std::size_t cursor = 0;
    for (int row = 0; row < size; ++row) {
        while (cursor < entries.size() && entries[cursor].row < row) {
            ++cursor;
        }
        while (cursor < entries.size() && entries[cursor].row == row) {
            const int column = entries[cursor].col;
            double value = 0.0;
            while (cursor < entries.size()
                   && entries[cursor].row == row
                   && entries[cursor].col == column) {
                value += entries[cursor].value;
                ++cursor;
            }
            if (value != 0.0) {
                result.colInd.push_back(column);
                result.values.push_back(value);
            }
        }
        result.rowPtr[static_cast<std::size_t>(row + 1)] =
            static_cast<int>(result.values.size());
    }
    return result;
}

UpperCsr buildUpperCsr(int size, const std::vector<double>& packed)
{
    if (packed.size() != packedUpperSize(size)) {
        throw std::runtime_error(
            "[Port core] Packed core matrix size mismatch.");
    }
    UpperCsr result;
    result.rowPtr.assign(static_cast<std::size_t>(size + 1), 0);
    for (int row = 0; row < size; ++row) {
        for (int column = row; column < size; ++column) {
            const double value = packed[packedUpperIndex(
                size, row, column)];
            if (value != 0.0) {
                result.colInd.push_back(column);
                result.values.push_back(value);
            }
        }
        result.rowPtr[static_cast<std::size_t>(row + 1)] =
            static_cast<int>(result.values.size());
    }
    return result;
}

bool isValidUpperCsr(int size, const UpperCsr& matrix)
{
    if (size < 0
        || matrix.rowPtr.size() != static_cast<std::size_t>(size + 1)
        || matrix.rowPtr.front() != 0
        || matrix.rowPtr.back() < 0
        || static_cast<std::size_t>(matrix.rowPtr.back())
            != matrix.colInd.size()
        || matrix.colInd.size() != matrix.values.size()) {
        return false;
    }
    for (int row = 0; row < size; ++row) {
        const int begin = matrix.rowPtr[static_cast<std::size_t>(row)];
        const int end = matrix.rowPtr[static_cast<std::size_t>(row + 1)];
        if (begin < 0 || begin > end
            || static_cast<std::size_t>(end) > matrix.colInd.size()) {
            return false;
        }
        int previous = row - 1;
        for (int entry = begin; entry < end; ++entry) {
            const int column = matrix.colInd[static_cast<std::size_t>(entry)];
            if (column < row || column >= size || column <= previous
                || !std::isfinite(
                    matrix.values[static_cast<std::size_t>(entry)])) {
                return false;
            }
            previous = column;
        }
    }
    return true;
}

bool isStrictlyIncreasingInRange(const std::vector<int>& values,
                                 int upperBound)
{
    int previous = -1;
    for (int value : values) {
        if (value < 0 || value >= upperBound || value <= previous) {
            return false;
        }
        previous = value;
    }
    return true;
}

template <typename T>
void appendHash(std::uint64_t& hash, const T& value)
{
    static_assert(std::is_trivially_copyable_v<T>);
    const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
    for (std::size_t i = 0; i < sizeof(T); ++i) {
        hash ^= bytes[i];
        hash *= UINT64_C(1099511628211);
    }
}

template <typename T>
void appendHash(std::uint64_t& hash, const std::vector<T>& values)
{
    static_assert(std::is_trivially_copyable_v<T>);
    appendHash(hash, values.size());
    const auto* bytes = reinterpret_cast<const unsigned char*>(values.data());
    const std::size_t bytesCount = values.size() * sizeof(T);
    for (std::size_t i = 0; i < bytesCount; ++i) {
        hash ^= bytes[i];
        hash *= UINT64_C(1099511628211);
    }
}

template <typename T>
void writeBinary(std::ofstream& output, const T& value)
{
    static_assert(std::is_trivially_copyable_v<T>);
    output.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template <typename T>
void writeBinaryVector(std::ofstream& output, const std::vector<T>& values)
{
    static_assert(std::is_trivially_copyable_v<T>);
    writeBinary(output, static_cast<std::uint64_t>(values.size()));
    if (!values.empty()) {
        output.write(reinterpret_cast<const char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(T)));
    }
}

template <typename T>
bool readBinary(std::ifstream& input, T& value)
{
    static_assert(std::is_trivially_copyable_v<T>);
    input.read(reinterpret_cast<char*>(&value), sizeof(T));
    return static_cast<bool>(input);
}

template <typename T>
bool readBinaryVector(std::ifstream& input,
                      std::vector<T>& values,
                      std::uint64_t maximumCount)
{
    static_assert(std::is_trivially_copyable_v<T>);
    std::uint64_t count = 0;
    if (!readBinary(input, count) || count > maximumCount) {
        return false;
    }
    values.resize(static_cast<std::size_t>(count));
    if (!values.empty()) {
        input.read(reinterpret_cast<char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(T)));
    }
    return static_cast<bool>(input);
}

double vectorNormSquared(const std::vector<double>& values)
{
    double result = 0.0;
    for (double value : values) {
        result += value * value;
    }
    return result;
}

class ScopedMklThreads {
public:
    explicit ScopedMklThreads(int threads)
    {
#ifdef USE_MKL_PARDISO
        previous_ = mkl_get_max_threads();
        mkl_set_num_threads_local(std::max(1, threads));
#else
        (void)threads;
#endif
    }

    ~ScopedMklThreads()
    {
#ifdef USE_MKL_PARDISO
        mkl_set_num_threads_local(previous_);
#endif
    }

private:
#ifdef USE_MKL_PARDISO
    int previous_ = 1;
#endif
};

} // namespace

struct PortCoreSolver::Impl {
    struct PortBlock {
        int leftDomain = -1;
        int rightDomain = -1;
        int membershipDofs = 0;
        std::vector<int> leafGamma;
        std::vector<int> activeCoreIndices;
        std::vector<int> matrixRowPtr;
        std::vector<int> matrixColInd;
        std::vector<double> matrixValues;
        // Column-major leafDofs x activeCoreIndices matrix D_p^{-1} E_p.
        std::vector<double> solvedCoupling;
        std::vector<double> coreUpdate;
        std::unique_ptr<SubdomainDirectSolver> factor;
        std::size_t matrixNonzeros = 0;
        std::vector<double> rhs;
        std::vector<double> solution;
        std::vector<double> coreContribution;
        std::vector<double> activeCoreSolution;
        double multiRhsSeconds = 0.0;
        double schurProductSeconds = 0.0;
    };

    const Model& model;
    int threads = 1;
    int coreThreads = 1;
    std::vector<PortBlock> ports;
    std::vector<int> rankOffsets;
    std::vector<int> separatorGamma;
    std::vector<int> gammaToSeparator;
    std::vector<int> gammaToPort;
    std::vector<int> gammaToLeaf;
    std::unique_ptr<SubdomainDirectSolver> coreFactor;
    std::size_t coreNonzeros = 0;
    double couplingSymmetryRelativeError = 0.0;
    double setupTime = 0.0;
    double symbolicTime = 0.0;
    double numericalTime = 0.0;
    int symbolicCalls = 0;
    int numericalCalls = 0;
    bool artifactCacheHit = false;
    double artifactCacheLoadTime = 0.0;
    double artifactCacheSaveTime = 0.0;
    std::size_t artifactCacheBytes = 0;
    double partitionTime = 0.0;
    double leafCsrTime = 0.0;
    double leafFactorTime = 0.0;
    double couplingAssemblyTime = 0.0;
    double eliminationTime = 0.0;
    int eliminationWorkers = 1;
    int eliminationMklThreads = 1;
    double multiRhsPortTime = 0.0;
    double schurProductPortTime = 0.0;
    double coreAccumulationTime = 0.0;
    double coreCsrTime = 0.0;
    double coreFactorTime = 0.0;
    double validatedRelativeResidual = 0.0;
    int solveCalls = 0;

    Impl(const Model& inputModel,
         const ddm_schur::InterfacePartition& partition,
         int localSolveThreads,
         int corePardisoThreads,
         const std::filesystem::path& outputDirectory,
         const std::filesystem::path& cachePath)
        : model(inputModel), threads(std::max(1, localSolveThreads)),
          coreThreads(std::max(1, corePardisoThreads))
    {
        const auto setupStart = Clock::now();
        if (model.interfaceDofs
                != static_cast<int>(partition.interfaceGlobalDofs.size())
            || model.interfaceGlobalDofs != partition.interfaceGlobalDofs) {
            throw std::runtime_error(
                "[Port core] Interface partition/model ordering mismatch.");
        }

        rankOffsets.resize(model.subdomains.size() + 1, 0);
        for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
            rankOffsets[slot + 1] = rankOffsets[slot]
                + model.subdomains[slot].rank;
        }
        if (rankOffsets.back() != model.totalLocalRank) {
            throw std::runtime_error("[Port core] Total local rank mismatch.");
        }

        std::map<std::pair<int, int>, int> portByDomains;
        for (const ddm_schur::DomainBlocks& domain : partition.domains) {
            for (const auto& [neighbor, unused] :
                 domain.interfaceGlobalDofsByNeighbor) {
                (void)unused;
                const auto key = std::minmax(domain.domainId, neighbor);
                if (portByDomains.find(key) == portByDomains.end()) {
                    const int id = static_cast<int>(portByDomains.size());
                    portByDomains.emplace(key, id);
                }
            }
        }
        ports.resize(portByDomains.size());
        for (const auto& [domains, id] : portByDomains) {
            ports[static_cast<std::size_t>(id)].leftDomain = domains.first;
            ports[static_cast<std::size_t>(id)].rightDomain = domains.second;
        }

        std::vector<std::vector<int>> memberships(
            static_cast<std::size_t>(model.interfaceDofs));
        for (const ddm_schur::DomainBlocks& domain : partition.domains) {
            for (const auto& [neighbor, globalDofs] :
                 domain.interfaceGlobalDofsByNeighbor) {
                const auto key = std::minmax(domain.domainId, neighbor);
                const int port = portByDomains.at(key);
                for (int global : globalDofs) {
                    if (global < 0
                        || global >= static_cast<int>(partition.globalToInterface.size())) {
                        throw std::runtime_error(
                            "[Port core] Physical port DOF is out of range.");
                    }
                    const int gamma = partition.globalToInterface[
                        static_cast<std::size_t>(global)];
                    if (gamma < 0 || gamma >= model.interfaceDofs) {
                        throw std::runtime_error(
                            "[Port core] Physical port contains a non-interface DOF.");
                    }
                    memberships[static_cast<std::size_t>(gamma)].push_back(port);
                }
            }
        }

        std::vector<int> initialOwner(
            static_cast<std::size_t>(model.interfaceDofs), -1);
        std::vector<char> separator(
            static_cast<std::size_t>(model.interfaceDofs), 0);
        for (int gamma = 0; gamma < model.interfaceDofs; ++gamma) {
            std::vector<int>& owned = memberships[static_cast<std::size_t>(gamma)];
            std::sort(owned.begin(), owned.end());
            owned.erase(std::unique(owned.begin(), owned.end()), owned.end());
            if (owned.size() == 1) {
                initialOwner[static_cast<std::size_t>(gamma)] = owned.front();
                ++ports[static_cast<std::size_t>(owned.front())].membershipDofs;
            } else {
                separator[static_cast<std::size_t>(gamma)] = 1;
            }
        }

        // Cover every direct A_GammaGamma edge between two physical ports with
        // a separator endpoint.  A deterministic maximum-degree vertex-cover
        // heuristic retains substantially more independent trace unknowns than
        // promoting both endpoints while preserving exact block separation.
        std::set<std::pair<int, int>> crossPortEdgeSet;
        for (const InterfaceEntry& entry : model.interfaceEntries) {
            if (entry.row == entry.column || entry.value == 0.0) {
                continue;
            }
            const int left = initialOwner[static_cast<std::size_t>(entry.row)];
            const int right = initialOwner[static_cast<std::size_t>(entry.column)];
            if (left >= 0 && right >= 0 && left != right) {
                crossPortEdgeSet.emplace(
                    std::min(entry.row, entry.column),
                    std::max(entry.row, entry.column));
            }
        }
        std::vector<std::pair<int, int>> crossPortEdges(
            crossPortEdgeSet.begin(), crossPortEdgeSet.end());
        std::vector<std::vector<int>> incident(
            static_cast<std::size_t>(model.interfaceDofs));
        std::vector<int> degree(
            static_cast<std::size_t>(model.interfaceDofs), 0);
        for (int edge = 0; edge < static_cast<int>(crossPortEdges.size()); ++edge) {
            const auto [left, right] = crossPortEdges[static_cast<std::size_t>(edge)];
            incident[static_cast<std::size_t>(left)].push_back(edge);
            incident[static_cast<std::size_t>(right)].push_back(edge);
            ++degree[static_cast<std::size_t>(left)];
            ++degree[static_cast<std::size_t>(right)];
        }
        std::priority_queue<std::pair<int, int>> candidates;
        for (int gamma = 0; gamma < model.interfaceDofs; ++gamma) {
            if (degree[static_cast<std::size_t>(gamma)] > 0) {
                candidates.emplace(degree[static_cast<std::size_t>(gamma)], gamma);
            }
        }
        std::vector<char> edgeCovered(crossPortEdges.size(), 0);
        std::size_t remainingEdges = crossPortEdges.size();
        while (remainingEdges > 0) {
            if (candidates.empty()) {
                throw std::runtime_error(
                    "[Port core] Cross-port vertex cover stalled.");
            }
            const auto [candidateDegree, gamma] = candidates.top();
            candidates.pop();
            if (candidateDegree != degree[static_cast<std::size_t>(gamma)]
                || candidateDegree == 0) {
                continue;
            }
            separator[static_cast<std::size_t>(gamma)] = 1;
            degree[static_cast<std::size_t>(gamma)] = 0;
            for (int edge : incident[static_cast<std::size_t>(gamma)]) {
                if (edgeCovered[static_cast<std::size_t>(edge)] != 0) {
                    continue;
                }
                edgeCovered[static_cast<std::size_t>(edge)] = 1;
                --remainingEdges;
                const auto [left, right] =
                    crossPortEdges[static_cast<std::size_t>(edge)];
                const int neighbor = left == gamma ? right : left;
                if (degree[static_cast<std::size_t>(neighbor)] > 0) {
                    --degree[static_cast<std::size_t>(neighbor)];
                    candidates.emplace(
                        degree[static_cast<std::size_t>(neighbor)], neighbor);
                }
            }
        }

        gammaToSeparator.assign(static_cast<std::size_t>(model.interfaceDofs), -1);
        gammaToPort.assign(static_cast<std::size_t>(model.interfaceDofs), -1);
        gammaToLeaf.assign(static_cast<std::size_t>(model.interfaceDofs), -1);
        for (int gamma = 0; gamma < model.interfaceDofs; ++gamma) {
            if (separator[static_cast<std::size_t>(gamma)] != 0
                || initialOwner[static_cast<std::size_t>(gamma)] < 0) {
                gammaToSeparator[static_cast<std::size_t>(gamma)] =
                    static_cast<int>(separatorGamma.size());
                separatorGamma.push_back(gamma);
                continue;
            }
            const int port = initialOwner[static_cast<std::size_t>(gamma)];
            PortBlock& block = ports[static_cast<std::size_t>(port)];
            gammaToPort[static_cast<std::size_t>(gamma)] = port;
            gammaToLeaf[static_cast<std::size_t>(gamma)] =
                static_cast<int>(block.leafGamma.size());
            block.leafGamma.push_back(gamma);
        }

        const int separatorSize = static_cast<int>(separatorGamma.size());
        const int coreSize = separatorSize + model.totalLocalRank;
        if (coreSize <= 0) {
            throw std::runtime_error("[Port core] Empty global core.");
        }
        const std::uint64_t operatorSignature =
            portCoreOperatorSignature(model, partition);
        partitionTime = elapsed(setupStart);

        UpperCsr coreCsr;
        const auto loadArtifacts = [&]() {
            if (cachePath.empty() || !std::filesystem::exists(cachePath)) {
                return false;
            }
            const auto loadStart = Clock::now();
            try {
                std::ifstream input(cachePath, std::ios::binary);
                std::uint64_t magic = 0;
                int version = 0;
                std::uint64_t signature = 0;
                int cachedInterfaceDofs = 0;
                int cachedSeparatorSize = 0;
                int cachedRank = 0;
                int cachedCoreSize = 0;
                int cachedPorts = 0;
                double cachedSymmetry = 0.0;
                bool valid = static_cast<bool>(input)
                    && readBinary(input, magic)
                    && readBinary(input, version)
                    && readBinary(input, signature)
                    && readBinary(input, cachedInterfaceDofs)
                    && readBinary(input, cachedSeparatorSize)
                    && readBinary(input, cachedRank)
                    && readBinary(input, cachedCoreSize)
                    && readBinary(input, cachedPorts)
                    && readBinary(input, cachedSymmetry)
                    && magic == portCoreCacheMagic
                    && version == portCoreCacheVersion
                    && signature == operatorSignature
                    && cachedInterfaceDofs == model.interfaceDofs
                    && cachedSeparatorSize == separatorSize
                    && cachedRank == model.totalLocalRank
                    && cachedCoreSize == coreSize
                    && cachedPorts == static_cast<int>(ports.size())
                    && std::isfinite(cachedSymmetry);

                struct CachedPort {
                    int leftDomain = -1;
                    int rightDomain = -1;
                    int membershipDofs = 0;
                    std::vector<int> leafGamma;
                    std::vector<int> activeCoreIndices;
                    UpperCsr matrix;
                    std::vector<double> solvedCoupling;
                };
                std::vector<CachedPort> cached(ports.size());
                for (std::size_t port = 0; valid && port < cached.size(); ++port) {
                    CachedPort& block = cached[port];
                    const std::uint64_t leafCount =
                        static_cast<std::uint64_t>(ports[port].leafGamma.size());
                    const std::uint64_t maximumLeafNonzeros =
                        leafCount * (leafCount + 1) / 2;
                    valid = readBinary(input, block.leftDomain)
                        && readBinary(input, block.rightDomain)
                        && readBinary(input, block.membershipDofs)
                        && readBinaryVector(input, block.leafGamma, leafCount)
                        && readBinaryVector(input, block.activeCoreIndices,
                                            static_cast<std::uint64_t>(coreSize))
                        && readBinaryVector(input, block.matrix.rowPtr,
                                            leafCount + 1)
                        && readBinaryVector(input, block.matrix.colInd,
                                            maximumLeafNonzeros)
                        && readBinaryVector(input, block.matrix.values,
                                            maximumLeafNonzeros)
                        && readBinaryVector(input, block.solvedCoupling,
                            leafCount * static_cast<std::uint64_t>(coreSize));
                    if (!valid) {
                        break;
                    }
                    const std::uint64_t expectedSolved = leafCount
                        * static_cast<std::uint64_t>(
                            block.activeCoreIndices.size());
                    valid = block.leftDomain == ports[port].leftDomain
                        && block.rightDomain == ports[port].rightDomain
                        && block.membershipDofs == ports[port].membershipDofs
                        && block.leafGamma == ports[port].leafGamma
                        && block.solvedCoupling.size() == expectedSolved
                        && isStrictlyIncreasingInRange(
                            block.activeCoreIndices, coreSize)
                        && isValidUpperCsr(
                            static_cast<int>(leafCount), block.matrix);
                }
                const std::uint64_t maximumCoreNonzeros =
                    static_cast<std::uint64_t>(coreSize)
                    * static_cast<std::uint64_t>(coreSize + 1) / 2;
                valid = valid
                    && readBinaryVector(input, coreCsr.rowPtr,
                        static_cast<std::uint64_t>(coreSize + 1))
                    && readBinaryVector(input, coreCsr.colInd,
                        maximumCoreNonzeros)
                    && readBinaryVector(input, coreCsr.values,
                        maximumCoreNonzeros)
                    && isValidUpperCsr(coreSize, coreCsr);
                if (!valid) {
                    artifactCacheLoadTime = elapsed(loadStart);
                    std::cout << "[Port core cache] Validation failed; rebuilding.\n";
                    return false;
                }
                for (std::size_t port = 0; port < ports.size(); ++port) {
                    ports[port].activeCoreIndices =
                        std::move(cached[port].activeCoreIndices);
                    ports[port].matrixRowPtr =
                        std::move(cached[port].matrix.rowPtr);
                    ports[port].matrixColInd =
                        std::move(cached[port].matrix.colInd);
                    ports[port].matrixValues =
                        std::move(cached[port].matrix.values);
                    ports[port].solvedCoupling =
                        std::move(cached[port].solvedCoupling);
                    ports[port].matrixNonzeros =
                        ports[port].matrixValues.size();
                }
                couplingSymmetryRelativeError = cachedSymmetry;
                coreNonzeros = coreCsr.values.size();
                artifactCacheLoadTime = elapsed(loadStart);
                std::error_code sizeError;
                artifactCacheBytes = static_cast<std::size_t>(
                    std::filesystem::file_size(cachePath, sizeError));
                if (sizeError) {
                    artifactCacheBytes = 0;
                }
                return true;
            } catch (const std::exception& error) {
                artifactCacheLoadTime = elapsed(loadStart);
                std::cout << "[Port core cache] Load failed: "
                          << error.what() << "; rebuilding.\n";
                return false;
            }
        };

        artifactCacheHit = loadArtifacts();
        std::vector<std::map<int, std::map<int, double>>> portCouplings;
        std::vector<double> corePacked;
        if (!artifactCacheHit) {
            const auto assemblyStart = Clock::now();
            std::vector<std::vector<MatrixEntry>> portEntries(ports.size());
            portCouplings.resize(ports.size());
            corePacked.assign(packedUpperSize(coreSize), 0.0);

            for (const InterfaceEntry& entry : model.interfaceEntries) {
                if (entry.row > entry.column || entry.value == 0.0) {
                    continue;
                }
                const int rowPort = gammaToPort[
                    static_cast<std::size_t>(entry.row)];
                const int colPort = gammaToPort[
                    static_cast<std::size_t>(entry.column)];
                const int rowSeparator = gammaToSeparator[
                    static_cast<std::size_t>(entry.row)];
                const int colSeparator = gammaToSeparator[
                    static_cast<std::size_t>(entry.column)];
                if (rowPort >= 0 && colPort >= 0) {
                    if (rowPort != colPort) {
                        throw std::runtime_error(
                            "[Port core] Cross-port leaf coupling escaped separator detection.");
                    }
                    portEntries[static_cast<std::size_t>(rowPort)].push_back({
                        gammaToLeaf[static_cast<std::size_t>(entry.row)],
                        gammaToLeaf[static_cast<std::size_t>(entry.column)],
                        entry.value});
                } else if (rowSeparator >= 0 && colSeparator >= 0) {
                    addPackedUpper(corePacked, coreSize,
                        rowSeparator, colSeparator, entry.value);
                } else {
                    const int port = rowPort >= 0 ? rowPort : colPort;
                    const int leaf = rowPort >= 0
                        ? gammaToLeaf[static_cast<std::size_t>(entry.row)]
                        : gammaToLeaf[static_cast<std::size_t>(entry.column)];
                    const int core = rowSeparator >= 0
                        ? rowSeparator : colSeparator;
                    portCouplings[static_cast<std::size_t>(port)][core][leaf]
                        += entry.value;
                }
            }

            double maximumCoupling = 0.0;
            double maximumCouplingDifference = 0.0;
            for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
                const SubdomainModel& local = model.subdomains[slot];
                const SubdomainModel& data = templatePayload(model, local);
                const int rankOffset = separatorSize + rankOffsets[slot];
                for (int row = 0; row < local.rank; ++row) {
                    for (int column = row; column < local.rank; ++column) {
                        const double left = data.reducedInterior[
                            static_cast<std::size_t>(
                                row * local.rank + column)];
                        const double right = data.reducedInterior[
                            static_cast<std::size_t>(
                                column * local.rank + row)];
                        addPackedUpper(corePacked, coreSize,
                            rankOffset + row, rankOffset + column,
                            0.5 * (left + right));
                    }
                }
                for (std::size_t localGamma = 0;
                     localGamma < local.interfaceIndices.size(); ++localGamma) {
                    const int gamma = local.interfaceIndices[localGamma];
                    const int port = gammaToPort[
                        static_cast<std::size_t>(gamma)];
                    const int separatorIndex = gammaToSeparator[
                        static_cast<std::size_t>(gamma)];
                    for (int mode = 0; mode < local.rank; ++mode) {
                        const double interfaceInterior =
                            data.reducedInterfaceInterior[
                                localGamma
                                    * static_cast<std::size_t>(local.rank)
                                + static_cast<std::size_t>(mode)];
                        const double interiorInterface =
                            data.reducedInteriorInterface[
                                static_cast<std::size_t>(mode)
                                    * local.localInterfaceDofs
                                + localGamma];
                        maximumCoupling = std::max(maximumCoupling,
                            std::max(std::abs(interfaceInterior),
                                     std::abs(interiorInterface)));
                        maximumCouplingDifference = std::max(
                            maximumCouplingDifference,
                            std::abs(interfaceInterior - interiorInterface));
                        const double value =
                            0.5 * (interfaceInterior + interiorInterface);
                        const int core = rankOffset + mode;
                        if (port >= 0) {
                            const int leaf = gammaToLeaf[
                                static_cast<std::size_t>(gamma)];
                            portCouplings[static_cast<std::size_t>(port)]
                                [core][leaf] += value;
                        } else {
                            addPackedUpper(corePacked, coreSize,
                                separatorIndex, core, value);
                        }
                    }
                }
            }
            couplingSymmetryRelativeError = maximumCoupling > 0.0
                ? maximumCouplingDifference / maximumCoupling : 0.0;
            couplingAssemblyTime = elapsed(assemblyStart);

            const auto csrStart = Clock::now();
            std::exception_ptr csrFailure;
#pragma omp parallel for num_threads(threads) if(threads > 1) schedule(dynamic)
            for (int port = 0; port < static_cast<int>(ports.size()); ++port) {
                try {
                    PortBlock& block = ports[static_cast<std::size_t>(port)];
                    if (block.leafGamma.empty()) {
                        continue;
                    }
                    UpperCsr csr = buildUpperCsr(
                        static_cast<int>(block.leafGamma.size()),
                        std::move(portEntries[static_cast<std::size_t>(port)]));
                    block.matrixRowPtr = std::move(csr.rowPtr);
                    block.matrixColInd = std::move(csr.colInd);
                    block.matrixValues = std::move(csr.values);
                    block.matrixNonzeros = block.matrixValues.size();
                } catch (...) {
#pragma omp critical(port_core_csr_failure)
                    {
                        if (!csrFailure) {
                            csrFailure = std::current_exception();
                        }
                    }
                }
            }
            if (csrFailure) {
                std::rethrow_exception(csrFailure);
            }
            leafCsrTime = elapsed(csrStart);
        }

        const auto factorLeavesStart = Clock::now();
        std::exception_ptr factorFailure;
#pragma omp parallel for num_threads(threads) if(threads > 1) schedule(dynamic)
        for (int port = 0; port < static_cast<int>(ports.size()); ++port) {
            try {
                PortBlock& block = ports[static_cast<std::size_t>(port)];
                if (block.leafGamma.empty()) {
                    continue;
                }
                ScopedMklSingleThread mklThreads;
                block.factor = std::make_unique<SubdomainDirectSolver>(
                    static_cast<int>(block.leafGamma.size()),
                    block.matrixRowPtr, block.matrixColInd,
                    block.matrixValues);
            } catch (...) {
#pragma omp critical(port_core_factor_failure)
                {
                    if (!factorFailure) {
                        factorFailure = std::current_exception();
                    }
                }
            }
        }
        if (factorFailure) {
            std::rethrow_exception(factorFailure);
        }
        leafFactorTime = elapsed(factorLeavesStart);

        if (!artifactCacheHit) {
            // Form each exact leaf Schur update E_p^T D_p^{-1} E_p. The
            // solved coupling is retained for online back substitution.
            std::uint64_t denseCouplingDoubles = 0;
            for (std::size_t port = 0; port < ports.size(); ++port) {
                PortBlock& block = ports[port];
                const auto& sparseCoupling = portCouplings[port];
                block.activeCoreIndices.reserve(sparseCoupling.size());
                for (const auto& [core, unused] : sparseCoupling) {
                    (void)unused;
                    block.activeCoreIndices.push_back(core);
                }
                denseCouplingDoubles +=
                    static_cast<std::uint64_t>(block.leafGamma.size())
                    * static_cast<std::uint64_t>(
                        block.activeCoreIndices.size());
            }
            constexpr std::uint64_t hybridThresholdDoubles =
                UINT64_C(100000000);
            if (denseCouplingDoubles >= hybridThresholdDoubles
                && threads >= 4) {
                eliminationWorkers = std::max(1, std::min({
                    static_cast<int>(ports.size()), threads / 2, 4}));
                eliminationMklThreads =
                    std::max(1, threads / eliminationWorkers);
            } else {
                eliminationWorkers = std::max(1, std::min(
                    static_cast<int>(ports.size()), threads));
                eliminationMklThreads = 1;
            }
            std::vector<int> eliminationOrder(ports.size(), 0);
            std::iota(eliminationOrder.begin(), eliminationOrder.end(), 0);
            std::sort(eliminationOrder.begin(), eliminationOrder.end(),
                [&](int left, int right) {
                    const PortBlock& a = ports[static_cast<std::size_t>(left)];
                    const PortBlock& b = ports[static_cast<std::size_t>(right)];
                    const long double aWork =
                        static_cast<long double>(a.leafGamma.size())
                        * static_cast<long double>(a.activeCoreIndices.size())
                        * static_cast<long double>(a.activeCoreIndices.size());
                    const long double bWork =
                        static_cast<long double>(b.leafGamma.size())
                        * static_cast<long double>(b.activeCoreIndices.size())
                        * static_cast<long double>(b.activeCoreIndices.size());
                    return aWork > bWork;
                });
            const auto eliminationStart = Clock::now();
            std::exception_ptr couplingFailure;
#pragma omp parallel for num_threads(eliminationWorkers) if(eliminationWorkers > 1) schedule(dynamic)
            for (int task = 0; task < static_cast<int>(ports.size()); ++task) {
                try {
                    const int port = eliminationOrder[
                        static_cast<std::size_t>(task)];
                    PortBlock& block = ports[static_cast<std::size_t>(port)];
                    if (block.leafGamma.empty()) {
                        continue;
                    }
                    const auto& sparseCoupling =
                        portCouplings[static_cast<std::size_t>(port)];
                    const int leafCount =
                        static_cast<int>(block.leafGamma.size());
                    const int activeCount =
                        static_cast<int>(block.activeCoreIndices.size());
                    if (activeCount == 0) {
                        continue;
                    }
                    std::vector<double> coupling(
                        static_cast<std::size_t>(leafCount)
                            * static_cast<std::size_t>(activeCount),
                        0.0);
                    for (int active = 0; active < activeCount; ++active) {
                        const int core = block.activeCoreIndices[
                            static_cast<std::size_t>(active)];
                        for (const auto& [leaf, value] :
                             sparseCoupling.at(core)) {
                            coupling[static_cast<std::size_t>(active)
                                * static_cast<std::size_t>(leafCount)
                                + static_cast<std::size_t>(leaf)] += value;
                        }
                    }
                    ScopedMklThreads mklThreads(eliminationMklThreads);
                    const auto multiRhsStart = Clock::now();
                    block.factor->solveMultiple(
                        coupling, activeCount, block.solvedCoupling,
                        eliminationMklThreads);
                    block.multiRhsSeconds = elapsed(multiRhsStart);
                    block.coreUpdate.assign(
                        static_cast<std::size_t>(activeCount)
                            * static_cast<std::size_t>(activeCount),
                        0.0);
                    const auto productStart = Clock::now();
#ifdef USE_MKL_PARDISO
                    cblas_dgemm(
                        CblasColMajor, CblasTrans, CblasNoTrans,
                        activeCount, activeCount, leafCount, 1.0,
                        coupling.data(), leafCount,
                        block.solvedCoupling.data(), leafCount,
                        0.0, block.coreUpdate.data(), activeCount);
#else
                    for (int column = 0; column < activeCount; ++column) {
                        for (int row = 0; row < activeCount; ++row) {
                            double value = 0.0;
                            for (int leaf = 0; leaf < leafCount; ++leaf) {
                                value += coupling[static_cast<std::size_t>(row)
                                             * leafCount + leaf]
                                    * block.solvedCoupling[
                                        static_cast<std::size_t>(column)
                                            * leafCount + leaf];
                            }
                            block.coreUpdate[
                                static_cast<std::size_t>(column)
                                    * activeCount + row] = value;
                        }
                    }
#endif
                    block.schurProductSeconds = elapsed(productStart);
                } catch (...) {
#pragma omp critical(port_core_coupling_failure)
                    {
                        if (!couplingFailure) {
                            couplingFailure = std::current_exception();
                        }
                    }
                }
            }
            if (couplingFailure) {
                std::rethrow_exception(couplingFailure);
            }
            eliminationTime = elapsed(eliminationStart);
            for (const PortBlock& block : ports) {
                multiRhsPortTime += block.multiRhsSeconds;
                schurProductPortTime += block.schurProductSeconds;
            }

            const auto accumulationStart = Clock::now();
            for (PortBlock& block : ports) {
                const int activeCount =
                    static_cast<int>(block.activeCoreIndices.size());
                for (int row = 0; row < activeCount; ++row) {
                    for (int column = row; column < activeCount; ++column) {
                        const double value = 0.5 * (
                            block.coreUpdate[static_cast<std::size_t>(column)
                                * activeCount + row]
                            + block.coreUpdate[static_cast<std::size_t>(row)
                                * activeCount + column]);
                        addPackedUpper(corePacked, coreSize,
                            block.activeCoreIndices[
                                static_cast<std::size_t>(row)],
                            block.activeCoreIndices[
                                static_cast<std::size_t>(column)],
                            -value);
                    }
                }
                std::vector<double>().swap(block.coreUpdate);
            }
            coreAccumulationTime = elapsed(accumulationStart);

            const auto coreCsrStart = Clock::now();
            coreCsr = buildUpperCsr(coreSize, corePacked);
            std::vector<double>().swap(corePacked);
            coreNonzeros = coreCsr.values.size();
            coreCsrTime = elapsed(coreCsrStart);
        }

        const auto coreFactorStart = Clock::now();
        {
            ScopedMklThreads mklThreads(coreThreads);
            coreFactor = std::make_unique<SubdomainDirectSolver>(
                coreSize, coreCsr.rowPtr, coreCsr.colInd, coreCsr.values);
        }
        coreFactorTime = elapsed(coreFactorStart);

        const auto saveArtifacts = [&]() {
            if (artifactCacheHit || cachePath.empty()) {
                return;
            }
            const auto saveStart = Clock::now();
            std::filesystem::path temporary = cachePath;
            temporary += ".tmp";
            try {
                if (!cachePath.parent_path().empty()) {
                    std::filesystem::create_directories(
                        cachePath.parent_path());
                }
                std::error_code cleanupError;
                std::filesystem::remove(temporary, cleanupError);
                std::ofstream output(temporary,
                    std::ios::binary | std::ios::trunc);
                if (!output) {
                    throw std::runtime_error(
                        "cannot open temporary cache for writing");
                }
                writeBinary(output, portCoreCacheMagic);
                writeBinary(output, portCoreCacheVersion);
                writeBinary(output, operatorSignature);
                writeBinary(output, model.interfaceDofs);
                writeBinary(output, separatorSize);
                writeBinary(output, model.totalLocalRank);
                writeBinary(output, coreSize);
                writeBinary(output, static_cast<int>(ports.size()));
                writeBinary(output, couplingSymmetryRelativeError);
                for (const PortBlock& block : ports) {
                    writeBinary(output, block.leftDomain);
                    writeBinary(output, block.rightDomain);
                    writeBinary(output, block.membershipDofs);
                    writeBinaryVector(output, block.leafGamma);
                    writeBinaryVector(output, block.activeCoreIndices);
                    writeBinaryVector(output, block.matrixRowPtr);
                    writeBinaryVector(output, block.matrixColInd);
                    writeBinaryVector(output, block.matrixValues);
                    writeBinaryVector(output, block.solvedCoupling);
                }
                writeBinaryVector(output, coreCsr.rowPtr);
                writeBinaryVector(output, coreCsr.colInd);
                writeBinaryVector(output, coreCsr.values);
                output.close();
                if (!output) {
                    throw std::runtime_error("cache write did not complete");
                }
                std::error_code replaceError;
                std::filesystem::remove(cachePath, replaceError);
                replaceError.clear();
                std::filesystem::rename(
                    temporary, cachePath, replaceError);
                if (replaceError) {
                    throw std::runtime_error(
                        "cannot publish cache: " + replaceError.message());
                }
                std::error_code sizeError;
                artifactCacheBytes = static_cast<std::size_t>(
                    std::filesystem::file_size(cachePath, sizeError));
                if (sizeError) {
                    artifactCacheBytes = 0;
                }
            } catch (const std::exception& error) {
                std::error_code cleanupError;
                std::filesystem::remove(temporary, cleanupError);
                std::cout << "[Port core cache] Save failed: "
                          << error.what() << '\n';
            }
            artifactCacheSaveTime = elapsed(saveStart);
        };
        saveArtifacts();

        for (PortBlock& block : ports) {
            std::vector<int>().swap(block.matrixRowPtr);
            std::vector<int>().swap(block.matrixColInd);
            std::vector<double>().swap(block.matrixValues);
        }
        std::vector<int>().swap(coreCsr.rowPtr);
        std::vector<int>().swap(coreCsr.colInd);
        std::vector<double>().swap(coreCsr.values);

        for (const PortBlock& block : ports) {
            if (block.factor) {
                symbolicTime += block.factor->symbolicAnalysisSeconds();
                numericalTime += block.factor->numericalFactorizationSeconds();
                symbolicCalls += block.factor->symbolicAnalysisCalls();
                numericalCalls += block.factor->numericalFactorizationCalls();
            }
        }
        symbolicTime += coreFactor->symbolicAnalysisSeconds();
        numericalTime += coreFactor->numericalFactorizationSeconds();
        symbolicCalls += coreFactor->symbolicAnalysisCalls();
        numericalCalls += coreFactor->numericalFactorizationCalls();
        setupTime = elapsed(setupStart);

        std::filesystem::create_directories(outputDirectory);
        std::ofstream partitionOutput(outputDirectory / "port_core_partition.csv");
        partitionOutput
            << "port,left_domain,right_domain,membership_dofs,leaf_dofs,"
               "separator_columns,active_core_columns,block_nonzeros,"
               "factor_memory_bytes,multi_rhs_seconds,"
               "schur_product_seconds\n";
        for (std::size_t port = 0; port < ports.size(); ++port) {
            const PortBlock& block = ports[port];
            const int separatorColumns = static_cast<int>(std::count_if(
                block.activeCoreIndices.begin(), block.activeCoreIndices.end(),
                [separatorSize](int value) { return value < separatorSize; }));
            partitionOutput << port << ',' << block.leftDomain << ','
                << block.rightDomain << ',' << block.membershipDofs << ','
                << block.leafGamma.size() << ',' << separatorColumns << ','
                << block.activeCoreIndices.size() << ','
                << block.matrixNonzeros << ','
                << (block.factor ? block.factor->memoryBytes() : 0) << ','
                << block.multiRhsSeconds << ','
                << block.schurProductSeconds << '\n';
        }

        std::ofstream summaryOutput(outputDirectory / "port_core_summary.csv");
        summaryOutput << std::setprecision(16)
            << "physical_ports,leaf_dofs,separator_dofs,reduced_core_dofs,"
               "core_dimension,core_nonzeros,coupling_symmetry_relative_error,"
               "cache_hit,cache_load_seconds,cache_save_seconds,cache_bytes,"
               "partition_seconds,coupling_assembly_seconds,leaf_csr_seconds,"
               "leaf_factor_seconds,elimination_seconds,elimination_workers,"
               "elimination_mkl_threads,multi_rhs_port_seconds,"
               "schur_product_port_seconds,core_accumulation_seconds,"
               "core_csr_seconds,core_factor_seconds,setup_seconds,"
               "symbolic_seconds,numerical_seconds,memory_bytes\n";
        const int leafDofs = std::accumulate(
            ports.begin(), ports.end(), 0,
            [](int sum, const PortBlock& block) {
                return sum + static_cast<int>(block.leafGamma.size());
            });
        summaryOutput << ports.size() << ',' << leafDofs << ','
            << separatorSize << ',' << model.totalLocalRank << ','
            << coreSize << ',' << coreNonzeros << ','
            << couplingSymmetryRelativeError << ','
            << (artifactCacheHit ? 1 : 0) << ','
            << artifactCacheLoadTime << ',' << artifactCacheSaveTime << ','
            << artifactCacheBytes << ',' << partitionTime << ','
            << couplingAssemblyTime << ',' << leafCsrTime << ','
            << leafFactorTime << ',' << eliminationTime << ','
            << eliminationWorkers << ',' << eliminationMklThreads << ','
            << multiRhsPortTime << ',' << schurProductPortTime << ','
            << coreAccumulationTime << ',' << coreCsrTime << ','
            << coreFactorTime << ',' << setupTime << ',' << symbolicTime << ','
            << numericalTime << ',' << memoryBytes() << '\n';

        std::cout << "[Port core] ports=" << ports.size()
                  << ", leaves=" << leafDofs
                  << ", separator=" << separatorSize
                  << ", reduced=" << model.totalLocalRank
                  << ", core=" << coreSize
                  << ", coupling symmetry="
                  << couplingSymmetryRelativeError
                  << ", cache=" << (artifactCacheHit ? "hit" : "miss")
                  << ", setup=" << setupTime << " s\n";
    }

    int coreDimension() const
    {
        return static_cast<int>(separatorGamma.size()) + model.totalLocalRank;
    }

    std::size_t memoryBytes() const
    {
        std::size_t bytes = separatorGamma.capacity() * sizeof(int)
            + gammaToSeparator.capacity() * sizeof(int)
            + gammaToPort.capacity() * sizeof(int)
            + gammaToLeaf.capacity() * sizeof(int)
            + rankOffsets.capacity() * sizeof(int)
            + (coreFactor ? coreFactor->memoryBytes() : 0);
        for (const PortBlock& block : ports) {
            bytes += block.leafGamma.capacity() * sizeof(int)
                + block.activeCoreIndices.capacity() * sizeof(int)
                + block.solvedCoupling.capacity() * sizeof(double)
                + block.coreUpdate.capacity() * sizeof(double)
                + block.rhs.capacity() * sizeof(double)
                + block.solution.capacity() * sizeof(double)
                + block.coreContribution.capacity() * sizeof(double)
                + block.activeCoreSolution.capacity() * sizeof(double)
                + (block.factor ? block.factor->memoryBytes() : 0);
        }
        return bytes;
    }

    PortCoreSolveResult solve(
        const std::vector<std::vector<double>>& projectedInteriorRhs,
        const std::vector<double>& interfaceRhs)
    {
        if (interfaceRhs.size() != static_cast<std::size_t>(model.interfaceDofs)
            || projectedInteriorRhs.size() != model.subdomains.size()) {
            throw std::runtime_error("[Port core] Reduced RHS dimensions are invalid.");
        }
        for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
            if (projectedInteriorRhs[slot].size()
                != static_cast<std::size_t>(model.subdomains[slot].rank)) {
                throw std::runtime_error(
                    "[Port core] Local reduced RHS rank is invalid.");
            }
        }

        PortCoreSolveResult result;
        result.interfaceTemperature.assign(
            static_cast<std::size_t>(model.interfaceDofs), 0.0);
        const auto forwardStart = Clock::now();
        std::exception_ptr solveFailure;
#pragma omp parallel for num_threads(threads) if(threads > 1) schedule(dynamic)
        for (int port = 0; port < static_cast<int>(ports.size()); ++port) {
            try {
                PortBlock& block = ports[static_cast<std::size_t>(port)];
                if (block.leafGamma.empty()) {
                    continue;
                }
                block.rhs.resize(block.leafGamma.size());
                for (std::size_t leaf = 0; leaf < block.leafGamma.size(); ++leaf) {
                    block.rhs[leaf] = interfaceRhs[static_cast<std::size_t>(
                        block.leafGamma[leaf])];
                }
                ScopedMklSingleThread mklThreads;
                block.factor->solve(block.rhs, block.solution);
            } catch (...) {
#pragma omp critical(port_core_solve_failure)
                {
                    if (!solveFailure) {
                        solveFailure = std::current_exception();
                    }
                }
            }
        }
        if (solveFailure) {
            std::rethrow_exception(solveFailure);
        }
        result.portForwardSeconds = elapsed(forwardStart);

        std::vector<double> coreRhs(
            static_cast<std::size_t>(coreDimension()), 0.0);
        for (std::size_t separator = 0;
             separator < separatorGamma.size(); ++separator) {
            coreRhs[separator] = interfaceRhs[static_cast<std::size_t>(
                separatorGamma[separator])];
        }
        for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
            const int offset = static_cast<int>(separatorGamma.size())
                + rankOffsets[slot];
            for (int mode = 0; mode < model.subdomains[slot].rank; ++mode) {
                coreRhs[static_cast<std::size_t>(offset + mode)] =
                    projectedInteriorRhs[slot][static_cast<std::size_t>(mode)];
            }
        }
        // Each port forms its dense contribution independently; only the final
        // scatter into the shared core RHS is serialized.
#pragma omp parallel for num_threads(threads) if(threads > 1) schedule(dynamic)
        for (int port = 0; port < static_cast<int>(ports.size()); ++port) {
            PortBlock& block = ports[static_cast<std::size_t>(port)];
            const int leafCount = static_cast<int>(block.leafGamma.size());
            const int activeCount =
                static_cast<int>(block.activeCoreIndices.size());
            block.coreContribution.assign(
                static_cast<std::size_t>(activeCount), 0.0);
            if (leafCount == 0 || activeCount == 0) {
                continue;
            }
#ifdef USE_MKL_PARDISO
            ScopedMklSingleThread mklThreads;
            cblas_dgemv(
                CblasColMajor, CblasTrans, leafCount, activeCount, 1.0,
                block.solvedCoupling.data(), leafCount,
                block.rhs.data(), 1, 0.0,
                block.coreContribution.data(), 1);
#else
            for (int active = 0; active < activeCount; ++active) {
                for (int leaf = 0; leaf < leafCount; ++leaf) {
                    block.coreContribution[static_cast<std::size_t>(active)]
                        += block.solvedCoupling[static_cast<std::size_t>(
                                active * leafCount + leaf)]
                        * block.rhs[static_cast<std::size_t>(leaf)];
                }
            }
#endif
        }
        for (const PortBlock& block : ports) {
            for (int active = 0;
                 active < static_cast<int>(block.activeCoreIndices.size());
                 ++active) {
                coreRhs[static_cast<std::size_t>(
                    block.activeCoreIndices[static_cast<std::size_t>(active)])]
                    -= block.coreContribution[static_cast<std::size_t>(active)];
            }
        }

        const auto coreStart = Clock::now();
        std::vector<double> coreSolution;
        {
            ScopedMklThreads mklThreads(coreThreads);
            coreFactor->solve(coreRhs, coreSolution);
        }
        result.coreSolveSeconds = elapsed(coreStart);

        for (std::size_t separator = 0;
             separator < separatorGamma.size(); ++separator) {
            result.interfaceTemperature[static_cast<std::size_t>(
                separatorGamma[separator])] = coreSolution[separator];
        }
        result.localReducedCoordinates.resize(model.subdomains.size());
        for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
            const int offset = static_cast<int>(separatorGamma.size())
                + rankOffsets[slot];
            result.localReducedCoordinates[slot].assign(
                coreSolution.begin() + offset,
                coreSolution.begin() + offset + model.subdomains[slot].rank);
        }

        const auto backStart = Clock::now();
#pragma omp parallel for num_threads(threads) if(threads > 1) schedule(dynamic)
        for (int port = 0; port < static_cast<int>(ports.size()); ++port) {
            PortBlock& block = ports[static_cast<std::size_t>(port)];
            const int leafCount = static_cast<int>(block.leafGamma.size());
            const int activeCount =
                static_cast<int>(block.activeCoreIndices.size());
            block.activeCoreSolution.resize(static_cast<std::size_t>(activeCount));
            for (int active = 0; active < activeCount; ++active) {
                block.activeCoreSolution[static_cast<std::size_t>(active)] =
                    coreSolution[static_cast<std::size_t>(
                        block.activeCoreIndices[static_cast<std::size_t>(active)])];
            }
#ifdef USE_MKL_PARDISO
            if (leafCount > 0 && activeCount > 0) {
                ScopedMklSingleThread mklThreads;
                cblas_dgemv(
                    CblasColMajor, CblasNoTrans, leafCount, activeCount, -1.0,
                    block.solvedCoupling.data(), leafCount,
                    block.activeCoreSolution.data(), 1, 1.0,
                    block.solution.data(), 1);
            }
#else
            for (int active = 0; active < activeCount; ++active) {
                const double coefficient =
                    block.activeCoreSolution[static_cast<std::size_t>(active)];
                for (int leaf = 0; leaf < leafCount; ++leaf) {
                    block.solution[static_cast<std::size_t>(leaf)] -=
                        block.solvedCoupling[static_cast<std::size_t>(
                            active * leafCount + leaf)] * coefficient;
                }
            }
#endif
            for (int leaf = 0; leaf < leafCount; ++leaf) {
                result.interfaceTemperature[static_cast<std::size_t>(
                    block.leafGamma[static_cast<std::size_t>(leaf)])] =
                        block.solution[static_cast<std::size_t>(leaf)];
            }
        }
        result.portBackSubstitutionSeconds = elapsed(backStart);

        // Measure against the original, unsymmetrized reduced equations.  A
        // nonzero value therefore also exposes any projection coupling drift.
        // The outer full-order residual gate still runs every step, so this
        // reduced audit is periodic rather than duplicating that cost 100 times.
        constexpr int residualAuditInterval = 25;
        if (solveCalls % residualAuditInterval == 0) {
            std::vector<double> interfaceResidual = interfaceRhs;
            std::vector<std::vector<double>> interiorResidual =
                projectedInteriorRhs;
            for (const InterfaceEntry& entry : model.interfaceEntries) {
                interfaceResidual[static_cast<std::size_t>(entry.row)] -=
                    entry.value * result.interfaceTemperature[
                        static_cast<std::size_t>(entry.column)];
            }
            for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
                const SubdomainModel& local = model.subdomains[slot];
                const SubdomainModel& data = templatePayload(model, local);
                const std::vector<double>& coordinates =
                    result.localReducedCoordinates[slot];
                for (std::size_t localGamma = 0;
                     localGamma < local.interfaceIndices.size(); ++localGamma) {
                    const int gamma = local.interfaceIndices[localGamma];
                    for (int mode = 0; mode < local.rank; ++mode) {
                        interfaceResidual[static_cast<std::size_t>(gamma)] -=
                            data.reducedInterfaceInterior[
                                localGamma * static_cast<std::size_t>(local.rank)
                                + static_cast<std::size_t>(mode)]
                            * coordinates[static_cast<std::size_t>(mode)];
                        interiorResidual[slot][static_cast<std::size_t>(mode)] -=
                            data.reducedInteriorInterface[
                                static_cast<std::size_t>(mode)
                                    * local.localInterfaceDofs
                                + localGamma]
                            * result.interfaceTemperature[
                                static_cast<std::size_t>(gamma)];
                    }
                }
                for (int row = 0; row < local.rank; ++row) {
                    for (int column = 0; column < local.rank; ++column) {
                        interiorResidual[slot][static_cast<std::size_t>(row)] -=
                            data.reducedInterior[static_cast<std::size_t>(
                                row * local.rank + column)]
                            * coordinates[static_cast<std::size_t>(column)];
                    }
                }
            }
            double residualSquared = vectorNormSquared(interfaceResidual);
            double rhsSquared = vectorNormSquared(interfaceRhs);
            for (std::size_t slot = 0; slot < interiorResidual.size(); ++slot) {
                residualSquared += vectorNormSquared(interiorResidual[slot]);
                rhsSquared += vectorNormSquared(projectedInteriorRhs[slot]);
            }
            validatedRelativeResidual = std::sqrt(residualSquared)
                / std::max(std::sqrt(rhsSquared),
                           std::numeric_limits<double>::min());
        }
        for (double value : result.interfaceTemperature) {
            if (!std::isfinite(value)) {
                validatedRelativeResidual =
                    std::numeric_limits<double>::infinity();
                break;
            }
        }
        result.reducedRelativeResidual = validatedRelativeResidual;
        ++solveCalls;
        return result;
    }
};

PortCoreSolver::PortCoreSolver(
    const Model& model,
    const ddm_schur::InterfacePartition& partition,
    int localSolveThreads,
    int corePardisoThreads,
    const std::filesystem::path& outputDirectory,
    const std::filesystem::path& cachePath)
    : impl_(std::make_unique<Impl>(
          model, partition, localSolveThreads, corePardisoThreads,
          outputDirectory, cachePath))
{
}

PortCoreSolver::~PortCoreSolver() = default;

PortCoreSolveResult PortCoreSolver::solve(
    const std::vector<std::vector<double>>& projectedInteriorRhs,
    const std::vector<double>& interfaceRhs)
{
    return impl_->solve(projectedInteriorRhs, interfaceRhs);
}

int PortCoreSolver::portCount() const
{
    return static_cast<int>(impl_->ports.size());
}

int PortCoreSolver::separatorDofs() const
{
    return static_cast<int>(impl_->separatorGamma.size());
}

int PortCoreSolver::coreDimension() const
{
    return impl_->coreDimension();
}

std::size_t PortCoreSolver::nonzeros() const
{
    std::size_t result = impl_->coreNonzeros;
    for (const Impl::PortBlock& block : impl_->ports) {
        result += block.matrixNonzeros;
    }
    return result;
}

std::size_t PortCoreSolver::memoryBytes() const
{
    return impl_->memoryBytes();
}

double PortCoreSolver::setupSeconds() const
{
    return impl_->setupTime;
}

double PortCoreSolver::symbolicSeconds() const
{
    return impl_->symbolicTime;
}

double PortCoreSolver::numericalSeconds() const
{
    return impl_->numericalTime;
}

int PortCoreSolver::symbolicCalls() const
{
    return impl_->symbolicCalls;
}

int PortCoreSolver::numericalCalls() const
{
    return impl_->numericalCalls;
}

bool PortCoreSolver::cacheHit() const
{
    return impl_->artifactCacheHit;
}

double PortCoreSolver::cacheLoadSeconds() const
{
    return impl_->artifactCacheLoadTime;
}

double PortCoreSolver::cacheSaveSeconds() const
{
    return impl_->artifactCacheSaveTime;
}

std::size_t PortCoreSolver::cacheBytes() const
{
    return impl_->artifactCacheBytes;
}

double PortCoreSolver::partitionSeconds() const
{
    return impl_->partitionTime;
}

double PortCoreSolver::leafCsrSeconds() const
{
    return impl_->leafCsrTime;
}

double PortCoreSolver::leafFactorSeconds() const
{
    return impl_->leafFactorTime;
}

double PortCoreSolver::couplingAssemblySeconds() const
{
    return impl_->couplingAssemblyTime;
}

double PortCoreSolver::eliminationSeconds() const
{
    return impl_->eliminationTime;
}

double PortCoreSolver::multiRhsPortSeconds() const
{
    return impl_->multiRhsPortTime;
}

double PortCoreSolver::schurProductPortSeconds() const
{
    return impl_->schurProductPortTime;
}

double PortCoreSolver::coreAccumulationSeconds() const
{
    return impl_->coreAccumulationTime;
}

double PortCoreSolver::coreCsrSeconds() const
{
    return impl_->coreCsrTime;
}

double PortCoreSolver::coreFactorSeconds() const
{
    return impl_->coreFactorTime;
}

} // namespace mor::local
