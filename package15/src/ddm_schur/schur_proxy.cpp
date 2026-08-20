#include "../sipg_core.hpp"
#include "schur_proxy.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

namespace ddm_schur {
namespace {

struct DisjointSet {
    std::vector<int> parent;
    std::vector<unsigned char> rank;

    explicit DisjointSet(int size)
        : parent(static_cast<std::size_t>(size)), rank(static_cast<std::size_t>(size), 0)
    {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int find(int value)
    {
        int root = value;
        while (parent[static_cast<std::size_t>(root)] != root) {
            root = parent[static_cast<std::size_t>(root)];
        }
        while (parent[static_cast<std::size_t>(value)] != value) {
            const int next = parent[static_cast<std::size_t>(value)];
            parent[static_cast<std::size_t>(value)] = root;
            value = next;
        }
        return root;
    }

    void unite(int left, int right)
    {
        left = find(left);
        right = find(right);
        if (left == right) {
            return;
        }
        if (rank[static_cast<std::size_t>(left)] < rank[static_cast<std::size_t>(right)]) {
            std::swap(left, right);
        }
        parent[static_cast<std::size_t>(right)] = left;
        if (rank[static_cast<std::size_t>(left)] == rank[static_cast<std::size_t>(right)]) {
            ++rank[static_cast<std::size_t>(left)];
        }
    }
};

struct FaceKey {
    std::array<int, 3> vertices{};

    bool operator==(const FaceKey& other) const { return vertices == other.vertices; }
};

struct FaceKeyHash {
    std::size_t operator()(const FaceKey& key) const
    {
        std::size_t seed = 0;
        for (int value : key.vertices) {
            seed ^= static_cast<std::size_t>(value + 0x9e3779b9)
                + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct GraphCsr {
    std::vector<int> rowPtr;
    std::vector<int> colInd;
    std::vector<int> edgeId;
    std::vector<std::pair<int, int>> edges;
    std::size_t undirectedEdges = 0;
};

double materialConductivity(const Material& material)
{
    return std::max({material.conductivityX,
                     material.conductivityY,
                     material.conductivityZ});
}

void addUndirectedEdge(std::vector<std::pair<int, int>>& edges, int left, int right)
{
    if (left < 0 || right < 0 || left == right) {
        return;
    }
    if (right < left) {
        std::swap(left, right);
    }
    edges.push_back({left, right});
}

void addClique(std::vector<std::pair<int, int>>& edges, std::vector<int> vertices)
{
    std::sort(vertices.begin(), vertices.end());
    vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        for (std::size_t j = i + 1; j < vertices.size(); ++j) {
            addUndirectedEdge(edges, vertices[i], vertices[j]);
        }
    }
}

std::vector<int> interfaceDofsForTet(const Tet& tet,
                                     const InterfacePartition& partition)
{
    std::vector<int> result;
    result.reserve(tet.dof.size());
    for (int globalDof : tet.dof) {
        if (globalDof >= 0
            && globalDof < static_cast<int>(partition.globalToInterface.size())) {
            const int gamma = partition.globalToInterface[static_cast<std::size_t>(globalDof)];
            if (gamma >= 0) {
                result.push_back(gamma);
            }
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

std::size_t appendMaterialConnectivity(
    const Mesh& mesh,
    const CaseConfig& physics,
    const InterfacePartition& partition,
    double threshold,
    std::vector<std::pair<int, int>>& edges,
    int& crossDomainComponents)
{
    crossDomainComponents = 0;
    if (!(threshold > 0.0) || !std::isfinite(threshold) || mesh.tets.empty()) {
        return 0;
    }

    std::vector<unsigned char> highK(mesh.tets.size(), 0);
    for (std::size_t tetIndex = 0; tetIndex < mesh.tets.size(); ++tetIndex) {
        highK[tetIndex] = materialConductivity(
            materialForTet(physics, mesh.tets[tetIndex])) >= threshold ? 1 : 0;
    }

    DisjointSet components(static_cast<int>(mesh.tets.size()));
    std::unordered_map<FaceKey, int, FaceKeyHash> faceOwner;
    faceOwner.reserve(mesh.tets.size());
    constexpr int localFaces[4][3] = {
        {1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}
    };
    for (int tetIndex = 0; tetIndex < static_cast<int>(mesh.tets.size()); ++tetIndex) {
        if (highK[static_cast<std::size_t>(tetIndex)] == 0) {
            continue;
        }
        const Tet& tet = mesh.tets[static_cast<std::size_t>(tetIndex)];
        for (const auto& localFace : localFaces) {
            FaceKey key{{tet.v[static_cast<std::size_t>(localFace[0])],
                         tet.v[static_cast<std::size_t>(localFace[1])],
                         tet.v[static_cast<std::size_t>(localFace[2])]}};
            std::sort(key.vertices.begin(), key.vertices.end());
            const auto insertion = faceOwner.emplace(key, tetIndex);
            if (!insertion.second) {
                components.unite(tetIndex, insertion.first->second);
            }
        }
    }
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        if (face.leftTet >= 0 && face.rightTet >= 0
            && highK[static_cast<std::size_t>(face.leftTet)] != 0
            && highK[static_cast<std::size_t>(face.rightTet)] != 0) {
            components.unite(face.leftTet, face.rightTet);
        }
    }

    std::unordered_map<int, std::vector<int>> componentDofs;
    for (int tetIndex = 0; tetIndex < static_cast<int>(mesh.tets.size()); ++tetIndex) {
        if (highK[static_cast<std::size_t>(tetIndex)] == 0) {
            continue;
        }
        std::vector<int> dofs = interfaceDofsForTet(
            mesh.tets[static_cast<std::size_t>(tetIndex)], partition);
        if (!dofs.empty()) {
            std::vector<int>& target = componentDofs[components.find(tetIndex)];
            target.insert(target.end(), dofs.begin(), dofs.end());
        }
    }

    const std::size_t before = edges.size();
    for (auto& item : componentDofs) {
        std::vector<int>& dofs = item.second;
        std::sort(dofs.begin(), dofs.end());
        dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
        std::set<int> domains;
        Vec3 minimum{std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()};
        Vec3 maximum{-std::numeric_limits<double>::max(),
                     -std::numeric_limits<double>::max(),
                     -std::numeric_limits<double>::max()};
        for (int gamma : dofs) {
            const Node& node = mesh.nodes[static_cast<std::size_t>(
                partition.interfaceGlobalDofs[static_cast<std::size_t>(gamma)])];
            domains.insert(node.subdomain);
            minimum.x = std::min(minimum.x, node.p.x);
            minimum.y = std::min(minimum.y, node.p.y);
            minimum.z = std::min(minimum.z, node.p.z);
            maximum.x = std::max(maximum.x, node.p.x);
            maximum.y = std::max(maximum.y, node.p.y);
            maximum.z = std::max(maximum.z, node.p.z);
        }
        if (domains.size() < 2 || dofs.size() < 2) {
            continue;
        }
        ++crossDomainComponents;
        const std::array<double, 3> extent{
            maximum.x - minimum.x,
            maximum.y - minimum.y,
            maximum.z - minimum.z};
        const int axis = static_cast<int>(std::distance(
            extent.begin(), std::max_element(extent.begin(), extent.end())));
        std::sort(dofs.begin(), dofs.end(), [&](int left, int right) {
            const Vec3& a = mesh.nodes[static_cast<std::size_t>(
                partition.interfaceGlobalDofs[static_cast<std::size_t>(left)])].p;
            const Vec3& b = mesh.nodes[static_cast<std::size_t>(
                partition.interfaceGlobalDofs[static_cast<std::size_t>(right)])].p;
            const double av = axis == 0 ? a.x : (axis == 1 ? a.y : a.z);
            const double bv = axis == 0 ? b.x : (axis == 1 ? b.y : b.z);
            return std::tie(av, a.x, a.y, a.z, left)
                 < std::tie(bv, b.x, b.y, b.z, right);
        });
        // A sparse physical chain augments graph distance without replacing a
        // long conductive component by a single coarse/constant unknown.
        for (std::size_t i = 1; i < dofs.size(); ++i) {
            addUndirectedEdge(edges, dofs[i - 1], dofs[i]);
        }
    }
    return edges.size() - before;
}

GraphCsr buildInterfaceGraph(const Mesh& mesh,
                             const CaseConfig& physics,
                             const InterfacePartition& partition,
                             double highKThreshold,
                             bool useMaterialConnectivity,
                             std::size_t& materialStrongEdges,
                             int& crossDomainComponents)
{
    const int interfaceCount = static_cast<int>(partition.interfaceGlobalDofs.size());
    std::vector<std::pair<int, int>> edges;
    edges.reserve(partition.interfaceEntries.size());
    for (const Entry& entry : partition.interfaceEntries) {
        if (entry.value != 0.0) {
            addUndirectedEdge(edges, entry.row, entry.col);
        }
    }
    for (const Tet& tet : mesh.tets) {
        const std::vector<int> dofs = interfaceDofsForTet(tet, partition);
        if (dofs.size() > 1) {
            addClique(edges, dofs);
        }
    }
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        if (face.leftTet < 0 || face.rightTet < 0) {
            continue;
        }
        std::vector<int> dofs = interfaceDofsForTet(
            mesh.tets[static_cast<std::size_t>(face.leftTet)], partition);
        const std::vector<int> right = interfaceDofsForTet(
            mesh.tets[static_cast<std::size_t>(face.rightTet)], partition);
        dofs.insert(dofs.end(), right.begin(), right.end());
        addClique(edges, std::move(dofs));
    }
    materialStrongEdges = 0;
    crossDomainComponents = 0;
    if (useMaterialConnectivity) {
        materialStrongEdges = appendMaterialConnectivity(
            mesh, physics, partition, highKThreshold, edges, crossDomainComponents);
    }

    std::sort(edges.begin(), edges.end());
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    GraphCsr graph;
    graph.undirectedEdges = edges.size();
    graph.edges = edges;
    graph.rowPtr.assign(static_cast<std::size_t>(interfaceCount + 1), 0);
    for (const auto& edge : edges) {
        ++graph.rowPtr[static_cast<std::size_t>(edge.first + 1)];
        ++graph.rowPtr[static_cast<std::size_t>(edge.second + 1)];
    }
    for (int row = 0; row < interfaceCount; ++row) {
        graph.rowPtr[static_cast<std::size_t>(row + 1)] +=
            graph.rowPtr[static_cast<std::size_t>(row)];
    }
    graph.colInd.assign(graph.rowPtr.back(), -1);
    graph.edgeId.assign(graph.rowPtr.back(), -1);
    std::vector<int> cursor = graph.rowPtr;
    for (int edgeIndex = 0; edgeIndex < static_cast<int>(edges.size()); ++edgeIndex) {
        const auto& edge = edges[static_cast<std::size_t>(edgeIndex)];
        int& leftCursor = cursor[static_cast<std::size_t>(edge.first)];
        graph.colInd[static_cast<std::size_t>(leftCursor)] = edge.second;
        graph.edgeId[static_cast<std::size_t>(leftCursor++)] = edgeIndex;
        int& rightCursor = cursor[static_cast<std::size_t>(edge.second)];
        graph.colInd[static_cast<std::size_t>(rightCursor)] = edge.first;
        graph.edgeId[static_cast<std::size_t>(rightCursor++)] = edgeIndex;
    }
    return graph;
}

std::vector<unsigned char> distancesToRing3(const GraphCsr& graph, int source)
{
    const int count = static_cast<int>(graph.rowPtr.size()) - 1;
    constexpr unsigned char farDistance = 4;
    std::vector<unsigned char> distance(static_cast<std::size_t>(count), farDistance);
    std::queue<int> frontier;
    distance[static_cast<std::size_t>(source)] = 0;
    frontier.push(source);
    while (!frontier.empty()) {
        const int row = frontier.front();
        frontier.pop();
        const unsigned char next = static_cast<unsigned char>(
            distance[static_cast<std::size_t>(row)] + 1);
        if (next > 3) {
            continue;
        }
        for (int offset = graph.rowPtr[static_cast<std::size_t>(row)];
             offset < graph.rowPtr[static_cast<std::size_t>(row + 1)]; ++offset) {
            const int neighbor = graph.colInd[static_cast<std::size_t>(offset)];
            if (distance[static_cast<std::size_t>(neighbor)] > next) {
                distance[static_cast<std::size_t>(neighbor)] = next;
                frontier.push(neighbor);
            }
        }
    }
    return distance;
}

double deterministicCoefficient(int item, int trial)
{
    std::uint64_t state = static_cast<std::uint64_t>(item + 1)
        * UINT64_C(0x9E3779B97F4A7C15);
    state ^= static_cast<std::uint64_t>(trial + 1) * UINT64_C(0xD1B54A32D192ED03);
    state ^= state >> 30;
    state *= UINT64_C(0xBF58476D1CE4E5B9);
    state ^= state >> 27;
    state *= UINT64_C(0x94D049BB133111EB);
    state ^= state >> 31;
    const double unit = static_cast<double>(state >> 11)
        * (1.0 / 9007199254740992.0);
    return 2.0 * unit - 1.0;
}

std::vector<int> colorDistanceTwo(const GraphCsr& graph)
{
    const int count = static_cast<int>(graph.rowPtr.size()) - 1;
    std::vector<int> colors(static_cast<std::size_t>(count), -1);
    std::vector<int> marker(static_cast<std::size_t>(count), -1);
    int maximumColor = -1;
    for (int vertex = 0; vertex < count; ++vertex) {
        const auto markVertex = [&](int neighbor) {
            if (neighbor < vertex) {
                const int color = colors[static_cast<std::size_t>(neighbor)];
                if (color >= 0) {
                    marker[static_cast<std::size_t>(color)] = vertex;
                }
            }
        };
        markVertex(vertex);
        for (int offset = graph.rowPtr[static_cast<std::size_t>(vertex)];
             offset < graph.rowPtr[static_cast<std::size_t>(vertex + 1)]; ++offset) {
            const int neighbor = graph.colInd[static_cast<std::size_t>(offset)];
            markVertex(neighbor);
            for (int second = graph.rowPtr[static_cast<std::size_t>(neighbor)];
                 second < graph.rowPtr[static_cast<std::size_t>(neighbor + 1)]; ++second) {
                markVertex(graph.colInd[static_cast<std::size_t>(second)]);
            }
        }
        int color = 0;
        while (marker[static_cast<std::size_t>(color)] == vertex) {
            ++color;
        }
        colors[static_cast<std::size_t>(vertex)] = color;
        maximumColor = std::max(maximumColor, color);
    }
    return colors;
}

constexpr std::uint64_t proxyCacheMagic = UINT64_C(0x5350584341434831); // SPXCACH1
constexpr std::uint32_t proxyCacheVersion = 1;

void appendHashBytes(std::uint64_t& hash, const void* data, std::size_t bytes)
{
    const auto* begin = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < bytes; ++i) {
        hash ^= static_cast<std::uint64_t>(begin[i]);
        hash *= UINT64_C(1099511628211);
    }
}

template <class T>
void appendHashValue(std::uint64_t& hash, const T& value)
{
    appendHashBytes(hash, &value, sizeof(value));
}

std::uint64_t proxyOperatorFingerprint(const Mesh& mesh,
                                       const InterfacePartition& partition,
                                       double highConductivityThreshold,
                                       bool useMaterialConnectivity,
                                       int ring)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    appendHashValue(hash, highConductivityThreshold);
    appendHashValue(hash, useMaterialConnectivity);
    appendHashValue(hash, ring);
    appendHashValue(hash, partition.totalDofs);
    const std::size_t interfaceCount = partition.interfaceGlobalDofs.size();
    appendHashValue(hash, interfaceCount);
    for (int value : partition.interfaceGlobalDofs) {
        appendHashValue(hash, value);
    }
    auto hashEntries = [&](const std::vector<Entry>& entries) {
        const std::size_t count = entries.size();
        appendHashValue(hash, count);
        for (const Entry& entry : entries) {
            appendHashValue(hash, entry.row);
            appendHashValue(hash, entry.col);
            appendHashValue(hash, entry.value);
        }
    };
    hashEntries(partition.interfaceEntries);
    for (const DomainBlocks& domain : partition.domains) {
        appendHashValue(hash, domain.domainId);
        hashEntries(domain.interiorEntries);
        for (const auto& row : domain.interiorInterfaceRows) {
            const std::size_t count = row.size();
            appendHashValue(hash, count);
            for (const auto& entry : row) {
                appendHashValue(hash, entry.first);
                appendHashValue(hash, entry.second);
            }
        }
        for (const auto& row : domain.interfaceInteriorRows) {
            const std::size_t count = row.size();
            appendHashValue(hash, count);
            for (const auto& entry : row) {
                appendHashValue(hash, entry.first);
                appendHashValue(hash, entry.second);
            }
        }
    }
    const std::size_t tetCount = mesh.tets.size();
    appendHashValue(hash, tetCount);
    for (const Tet& tet : mesh.tets) {
        appendHashValue(hash, tet.subdomain);
        appendHashValue(hash, tet.domainEntity);
        for (int dof : tet.dof) {
            appendHashValue(hash, dof);
        }
    }
    return hash;
}

struct ProxyCacheData {
    std::uint64_t fingerprint = 0;
    int interfaceDofs = 0;
    GraphCsr graph;
    std::vector<int> colors;
    std::vector<double> diagonal;
    std::vector<double> edgeValues;
    double minimumRayleigh = 0.0;
    double diagonalShift = 0.0;
    double diagonalCompensation = 0.0;
    std::uint64_t valueHash = 0;
};

template <class T>
bool readPod(std::ifstream& in, T& value)
{
    return static_cast<bool>(in.read(reinterpret_cast<char*>(&value), sizeof(value)));
}

template <class T>
void writePod(std::ofstream& out, const T& value)
{
    out.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

template <class T>
bool readVector(std::ifstream& in, std::vector<T>& values)
{
    std::uint64_t count = 0;
    if (!readPod(in, count) || count > UINT64_C(2000000000)) {
        return false;
    }
    values.resize(static_cast<std::size_t>(count));
    return count == 0 || static_cast<bool>(in.read(
        reinterpret_cast<char*>(values.data()),
        static_cast<std::streamsize>(count * sizeof(T))));
}

template <class T>
void writeVector(std::ofstream& out, const std::vector<T>& values)
{
    const std::uint64_t count = static_cast<std::uint64_t>(values.size());
    writePod(out, count);
    if (!values.empty()) {
        out.write(reinterpret_cast<const char*>(values.data()),
                  static_cast<std::streamsize>(values.size() * sizeof(T)));
    }
}

bool loadProxyCache(const std::filesystem::path& path,
                    std::uint64_t expectedFingerprint,
                    ProxyCacheData& data)
{
    std::ifstream in(path, std::ios::binary);
    std::uint64_t magic = 0;
    std::uint32_t version = 0;
    if (!in || !readPod(in, magic) || !readPod(in, version)
        || magic != proxyCacheMagic || version != proxyCacheVersion
        || !readPod(in, data.fingerprint)
        || data.fingerprint != expectedFingerprint
        || !readPod(in, data.interfaceDofs)
        || !readVector(in, data.graph.rowPtr)
        || !readVector(in, data.graph.colInd)
        || !readVector(in, data.graph.edgeId)) {
        return false;
    }
    std::vector<int> edgePairs;
    if (!readVector(in, edgePairs) || edgePairs.size() % 2 != 0) {
        return false;
    }
    data.graph.edges.resize(edgePairs.size() / 2);
    for (std::size_t i = 0; i < data.graph.edges.size(); ++i) {
        data.graph.edges[i] = {edgePairs[2 * i], edgePairs[2 * i + 1]};
    }
    data.graph.undirectedEdges = data.graph.edges.size();
    const bool valid = readVector(in, data.colors)
        && readVector(in, data.diagonal)
        && readVector(in, data.edgeValues)
        && readPod(in, data.minimumRayleigh)
        && readPod(in, data.diagonalShift)
        && readPod(in, data.diagonalCompensation)
        && readPod(in, data.valueHash)
        && data.interfaceDofs >= 0
        && data.diagonal.size() == static_cast<std::size_t>(data.interfaceDofs)
        && data.colors.size() == data.diagonal.size()
        && data.edgeValues.size() == data.graph.edges.size();
    if (!valid) {
        return false;
    }
    std::uint64_t valueHash = UINT64_C(1469598103934665603);
    for (double value : data.diagonal) {
        appendHashBytes(valueHash, &value, sizeof(value));
    }
    for (double value : data.edgeValues) {
        appendHashBytes(valueHash, &value, sizeof(value));
    }
    return valueHash == data.valueHash;
}

void saveProxyCache(const std::filesystem::path& path, const ProxyCacheData& data)
{
    if (path.empty()) {
        return;
    }
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("[Schur proxy] Cannot create proxy cache file.");
    }
    writePod(out, proxyCacheMagic);
    writePod(out, proxyCacheVersion);
    writePod(out, data.fingerprint);
    writePod(out, data.interfaceDofs);
    writeVector(out, data.graph.rowPtr);
    writeVector(out, data.graph.colInd);
    writeVector(out, data.graph.edgeId);
    std::vector<int> edgePairs;
    edgePairs.reserve(2 * data.graph.edges.size());
    for (const auto& edge : data.graph.edges) {
        edgePairs.push_back(edge.first);
        edgePairs.push_back(edge.second);
    }
    writeVector(out, edgePairs);
    writeVector(out, data.colors);
    writeVector(out, data.diagonal);
    writeVector(out, data.edgeValues);
    writePod(out, data.minimumRayleigh);
    writePod(out, data.diagonalShift);
    writePod(out, data.diagonalCompensation);
    writePod(out, data.valueHash);
    out.close();
    if (!out) {
        throw std::runtime_error("[Schur proxy] Failed while writing proxy cache file.");
    }
    std::error_code error;
    std::filesystem::remove(path, error);
    error.clear();
    std::filesystem::rename(temporary, path, error);
    if (error) {
        throw std::runtime_error("[Schur proxy] Cannot install proxy cache file: "
                                 + error.message());
    }
}

class ScopedProxyMklThreads {
public:
    explicit ScopedProxyMklThreads(int threads)
    {
#ifdef USE_MKL_PARDISO
        previous_ = mkl_get_max_threads();
        mkl_set_num_threads_local(std::max(1, threads));
#else
        (void)threads;
#endif
    }
    ~ScopedProxyMklThreads()
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

class ProxyPardiso {
public:
    ProxyPardiso(int size,
                 const std::vector<std::pair<int, int>>& edges,
                 const std::vector<double>& edgeValues,
                 const std::vector<double>& diagonal,
                 int pardisoThreads)
    {
#ifdef USE_MKL_PARDISO
        pardisoThreads_ = std::max(1, pardisoThreads);
        ScopedProxyMklThreads threads(pardisoThreads_);
        n_ = static_cast<MKL_INT>(size);
        rowPtr_.assign(static_cast<std::size_t>(size + 1), 1);
        for (const auto& edge : edges) {
            ++rowPtr_[static_cast<std::size_t>(edge.first + 1)];
        }
        for (int row = 0; row < size; ++row) {
            rowPtr_[static_cast<std::size_t>(row + 1)] +=
                rowPtr_[static_cast<std::size_t>(row)];
        }
        colInd_.reserve(static_cast<std::size_t>(size) + edges.size());
        values_.reserve(static_cast<std::size_t>(size) + edges.size());
        std::size_t edgeIndex = 0;
        for (int row = 0; row < size; ++row) {
            colInd_.push_back(static_cast<MKL_INT>(row + 1));
            values_.push_back(diagonal[static_cast<std::size_t>(row)]);
            while (edgeIndex < edges.size()
                   && edges[edgeIndex].first == row) {
                colInd_.push_back(static_cast<MKL_INT>(edges[edgeIndex].second + 1));
                values_.push_back(edgeValues[edgeIndex]);
                ++edgeIndex;
            }
        }
        perm_.assign(static_cast<std::size_t>(size), 0);
        iparm_[0] = 1;
        iparm_[1] = 2;
        // The proxy is an FGMRES preconditioner, not the accepted solution.
        // The outer true-residual gate provides accuracy; avoid phase-33
        // iterative refinement inside every preconditioner application.
        iparm_[7] = 0;
        iparm_[9] = 13;
        iparm_[10] = 1;
        iparm_[12] = 1;
        iparm_[17] = -1;
        iparm_[18] = -1;
        std::vector<double> rhs(static_cast<std::size_t>(size), 0.0);
        std::vector<double> solution(static_cast<std::size_t>(size), 0.0);
        MKL_INT nrhs = 1;
        MKL_INT error = 0;
        MKL_INT phase = 11;
        auto start = std::chrono::steady_clock::now();
        pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase, &n_,
                values_.data(), rowPtr_.data(), colInd_.data(), perm_.data(),
                &nrhs, iparm_, &msglvl_, rhs.data(), solution.data(), &error);
        symbolicSeconds_ = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        symbolicCalls_ = 1;
        if (error != 0) {
            release();
            throw std::runtime_error(
                "[Schur proxy] PARDISO phase 11 failed with error "
                + std::to_string(error));
        }
        pardisoMemoryBytes_ = static_cast<std::size_t>(std::max<MKL_INT>(
            0, iparm_[14] + iparm_[15] + iparm_[16])) * 1024;
        constexpr std::size_t factorMemorySafetyLimit =
            static_cast<std::size_t>(12) * 1024 * 1024 * 1024;
        if (pardisoMemoryBytes_ > factorMemorySafetyLimit) {
            release();
            throw std::runtime_error(
                "[Schur proxy] PARDISO phase-11 memory estimate exceeds 12 GiB safety limit.");
        }
        phase = 22;
        error = 0;
        start = std::chrono::steady_clock::now();
        pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase, &n_,
                values_.data(), rowPtr_.data(), colInd_.data(), perm_.data(),
                &nrhs, iparm_, &msglvl_, rhs.data(), solution.data(), &error);
        numericalSeconds_ = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        numericalCalls_ = 1;
        if (error != 0) {
            release();
            throw std::runtime_error(
                "[Schur proxy] PARDISO phase 22 failed with error "
                + std::to_string(error));
        }
        factorized_ = true;
        pardisoMemoryBytes_ = std::max(
            pardisoMemoryBytes_,
            static_cast<std::size_t>(std::max<MKL_INT>(
                0, iparm_[14] + iparm_[15] + iparm_[16])) * 1024);
#else
        (void)size;
        (void)edges;
        (void)edgeValues;
        (void)diagonal;
        (void)pardisoThreads;
        throw std::runtime_error("[Schur proxy] MKL PARDISO is required.");
#endif
    }

    ~ProxyPardiso()
    {
#ifdef USE_MKL_PARDISO
        release();
#endif
    }

    void solve(const std::vector<double>& rhs, std::vector<double>& solution)
    {
#ifdef USE_MKL_PARDISO
        std::lock_guard<std::mutex> guard(mutex_);
        ScopedProxyMklThreads threads(pardisoThreads_);
        if (!factorized_ || rhs.size() != static_cast<std::size_t>(n_)) {
            throw std::runtime_error("[Schur proxy] Invalid phase-33 solve.");
        }
        // PARDISO phase 33 overwrites every solution entry. Preserve capacity
        // across Krylov iterations without an unnecessary full-vector clear.
        solution.resize(rhs.size());
        MKL_INT phase = 33;
        MKL_INT nrhs = 1;
        MKL_INT error = 0;
        const auto start = std::chrono::steady_clock::now();
        pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase, &n_,
                values_.data(), rowPtr_.data(), colInd_.data(), perm_.data(),
                &nrhs, iparm_, &msglvl_, const_cast<double*>(rhs.data()),
                solution.data(), &error);
        solveSeconds_ += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        ++solveCalls_;
        if (error != 0) {
            throw std::runtime_error(
                "[Schur proxy] PARDISO phase 33 failed with error "
                + std::to_string(error));
        }
#else
        (void)rhs;
        (void)solution;
#endif
    }

    double symbolicSeconds() const { return symbolicSeconds_; }
    double numericalSeconds() const { return numericalSeconds_; }
    double solveSeconds() const { return solveSeconds_; }
    int symbolicCalls() const { return symbolicCalls_; }
    int numericalCalls() const { return numericalCalls_; }
    int solveCalls() const { return solveCalls_; }
    void resetRuntimeCounters()
    {
        std::lock_guard<std::mutex> guard(mutex_);
        solveSeconds_ = 0.0;
        solveCalls_ = 0;
    }
    std::size_t memoryBytes() const
    {
#ifdef USE_MKL_PARDISO
        return values_.size() * sizeof(double)
            + rowPtr_.size() * sizeof(MKL_INT)
            + colInd_.size() * sizeof(MKL_INT)
            + perm_.size() * sizeof(MKL_INT)
            + pardisoMemoryBytes_;
#else
        return 0;
#endif
    }

private:
#ifdef USE_MKL_PARDISO
    void release() noexcept
    {
        if (n_ <= 0) {
            return;
        }
        ScopedProxyMklThreads threads(pardisoThreads_);
        MKL_INT phase = -1;
        MKL_INT nrhs = 1;
        MKL_INT error = 0;
        double dummy = 0.0;
        pardiso(pt_, &maxfct_, &mnum_, &mtype_, &phase, &n_,
                values_.data(), rowPtr_.data(), colInd_.data(), perm_.data(),
                &nrhs, iparm_, &msglvl_, &dummy, &dummy, &error);
        n_ = 0;
        factorized_ = false;
    }

    MKL_INT n_ = 0;
    std::vector<MKL_INT> rowPtr_;
    std::vector<MKL_INT> colInd_;
    std::vector<MKL_INT> perm_;
    std::vector<double> values_;
    void* pt_[64] = {};
    MKL_INT iparm_[64] = {};
    MKL_INT maxfct_ = 1;
    MKL_INT mnum_ = 1;
    MKL_INT mtype_ = 2;
    MKL_INT msglvl_ = 0;
    bool factorized_ = false;
    int pardisoThreads_ = 1;
    std::size_t pardisoMemoryBytes_ = 0;
#endif
    double symbolicSeconds_ = 0.0;
    double numericalSeconds_ = 0.0;
    double solveSeconds_ = 0.0;
    int symbolicCalls_ = 0;
    int numericalCalls_ = 0;
    int solveCalls_ = 0;
    mutable std::mutex mutex_;
};

struct CachedProxyFactor {
    std::uint64_t fingerprint = 0;
    std::shared_ptr<ProxyPardiso> solver;
};

std::mutex proxyFactorCacheMutex;
std::shared_ptr<CachedProxyFactor> proxyFactorCache;

} // namespace

struct SchurProxyPreconditioner::Impl {
    int ring = 0;
    std::size_t nonzeros = 0;
    double matrixDensity = 0.0;
    int colorCount = 0;
    int probingApplies = 0;
    int probingBlockSize = 1;
    int probingBlockCalls = 0;
    int validationApplies = 0;
    double setupTime = 0.0;
    double matrixSymmetryError = 0.0;
    double minimumRayleigh = 0.0;
    double globalDiagonalShift = 0.0;
    double totalDiagonalCompensation = 0.0;
    std::uint64_t matrixValueHash = 0;
    double blockMaximumDifference = 0.0;
    double blockRelativeDifference = 0.0;
    std::size_t graphMemoryBytes = 0;
    std::shared_ptr<ProxyPardiso> solver;
    bool matrixCacheHit = false;
    bool factorCacheHit = false;
};

SchurProxyPreconditioner::SchurProxyPreconditioner(
    const Mesh& mesh,
    const CaseConfig& physics,
    const InterfacePartition& partition,
    double highConductivityThreshold,
    bool useMaterialConnectivity,
    int requestedRing,
    int requestedBlockSize,
    int requestedPardisoThreads,
    bool validateBlockEquivalence,
    bool cacheEnabled,
    const std::filesystem::path& cachePath,
    const std::filesystem::path& outputDirectory,
    const ExactSchurApply& applyExactSchur,
    const ExactSchurApplyBlock& applyExactSchurBlock)
    : impl_(std::make_unique<Impl>())
{
    const auto setupStart = std::chrono::steady_clock::now();
    if (requestedRing != 1) {
        throw std::runtime_error(
            "[Schur proxy] Only the locality-approved 1-ring proxy may be constructed; "
            "2/3-ring patterns require a new memory review.");
    }
    impl_->ring = requestedRing;
    if (requestedBlockSize <= 0) {
        throw std::runtime_error("[Schur proxy] Probing block size must be positive.");
    }
    impl_->probingBlockSize = requestedBlockSize;
    const int count = static_cast<int>(partition.interfaceGlobalDofs.size());
    const std::uint64_t operatorFingerprint = proxyOperatorFingerprint(
        mesh, partition, highConductivityThreshold,
        useMaterialConnectivity, requestedRing);
    ProxyCacheData cacheData;
    GraphCsr graph;
    std::vector<int> colors;
    std::vector<double> diagonal;
    std::vector<double> edgeValues;
    if (cacheEnabled && !cachePath.empty()
        && loadProxyCache(cachePath, operatorFingerprint, cacheData)) {
        if (cacheData.interfaceDofs != count) {
            throw std::runtime_error(
                "[Schur proxy] Cache fingerprint matched but dimensions are inconsistent.");
        }
        graph = std::move(cacheData.graph);
        colors = std::move(cacheData.colors);
        diagonal = std::move(cacheData.diagonal);
        edgeValues = std::move(cacheData.edgeValues);
        impl_->minimumRayleigh = cacheData.minimumRayleigh;
        impl_->globalDiagonalShift = cacheData.diagonalShift;
        impl_->totalDiagonalCompensation = cacheData.diagonalCompensation;
        impl_->matrixValueHash = cacheData.valueHash;
        impl_->matrixCacheHit = true;
    } else {
        if (cacheEnabled && !cachePath.empty() && std::filesystem::exists(cachePath)) {
            std::cout << "[Schur proxy] cache fingerprint mismatch; rebuilding proxy.\n";
        }
        std::size_t materialEdges = 0;
        int materialComponents = 0;
        graph = buildInterfaceGraph(
            mesh, physics, partition, highConductivityThreshold,
            useMaterialConnectivity, materialEdges, materialComponents);
        colors = colorDistanceTwo(graph);
        diagonal.assign(static_cast<std::size_t>(count), 0.0);
        edgeValues.assign(graph.edges.size(), 0.0);
    }
    impl_->graphMemoryBytes = graph.rowPtr.size() * sizeof(int)
        + graph.colInd.size() * sizeof(int)
        + graph.edgeId.size() * sizeof(int)
        + graph.edges.size() * sizeof(std::pair<int, int>);

    impl_->colorCount = colors.empty()
        ? 0 : 1 + *std::max_element(colors.begin(), colors.end());
    if (!impl_->matrixCacheHit) {
    std::vector<std::vector<int>> columnsByColor(
        static_cast<std::size_t>(impl_->colorCount));
    for (int column = 0; column < count; ++column) {
        columnsByColor[static_cast<std::size_t>(colors[static_cast<std::size_t>(column)])]
            .push_back(column);
    }

    std::vector<double> omittedAbsolute(static_cast<std::size_t>(count), 0.0);
    std::vector<int> supportMarker(static_cast<std::size_t>(count), -1);
    for (int colorBegin = 0; colorBegin < impl_->colorCount;
         colorBegin += impl_->probingBlockSize) {
        const int activeBlock = std::min(
            impl_->probingBlockSize, impl_->colorCount - colorBegin);
        auto accumulateResponse = [&](int blockColumn,
                                      const std::vector<double>& response) {
            if (response.size() != static_cast<std::size_t>(count)) {
                throw std::runtime_error(
                    "[Schur proxy] Matrix-free block response has the wrong size.");
            }
            const int color = colorBegin + blockColumn;
            const std::vector<int>& columns =
                columnsByColor[static_cast<std::size_t>(color)];
            for (int column : columns) {
                diagonal[static_cast<std::size_t>(column)] =
                    response[static_cast<std::size_t>(column)];
                supportMarker[static_cast<std::size_t>(column)] = color;
                for (int offset = graph.rowPtr[static_cast<std::size_t>(column)];
                     offset < graph.rowPtr[static_cast<std::size_t>(column + 1)]; ++offset) {
                    const int row = graph.colInd[static_cast<std::size_t>(offset)];
                    supportMarker[static_cast<std::size_t>(row)] = color;
                    const int edge = graph.edgeId[static_cast<std::size_t>(offset)];
                    edgeValues[static_cast<std::size_t>(edge)] +=
                        0.5 * response[static_cast<std::size_t>(row)];
                }
            }
            for (int row = 0; row < count; ++row) {
                if (supportMarker[static_cast<std::size_t>(row)] != color) {
                    omittedAbsolute[static_cast<std::size_t>(row)] +=
                        std::abs(response[static_cast<std::size_t>(row)]);
                }
            }
        };

        if (impl_->probingBlockSize == 1) {
            std::vector<double> probes(static_cast<std::size_t>(count), 0.0);
            for (int column : columnsByColor[static_cast<std::size_t>(colorBegin)]) {
                probes[static_cast<std::size_t>(column)] = 1.0;
            }
            std::vector<double> responses;
            applyExactSchur(probes, responses);
            accumulateResponse(0, responses);
        } else {
            int consumedResponses = 0;
            applyExactSchurBlock(
                colors, colorBegin, impl_->probingBlockSize,
                [&](int blockColumn, const std::vector<double>& response) {
                    // Keep phase 33 at the requested nrhs even for the final
                    // partial color block.  Extra all-zero columns are not
                    // Schur probes and are intentionally excluded from the
                    // exact-vector count and proxy accumulation.
                    if (blockColumn < activeBlock) {
                        accumulateResponse(blockColumn, response);
                        ++consumedResponses;
                    }
                });
            if (consumedResponses != activeBlock) {
                throw std::runtime_error(
                    "[Schur proxy] Matrix-free block apply returned the wrong column count.");
            }
        }
        ++impl_->probingBlockCalls;
        impl_->probingApplies += activeBlock;
    }

    if (validateBlockEquivalence && impl_->probingBlockSize > 1) {
        std::vector<double> scalarDiagonal(static_cast<std::size_t>(count), 0.0);
        std::vector<double> scalarEdges(graph.edges.size(), 0.0);
        std::vector<double> scalarOmitted(static_cast<std::size_t>(count), 0.0);
        std::vector<int> scalarSupport(static_cast<std::size_t>(count), -1);
        std::vector<double> scalarProbe(static_cast<std::size_t>(count), 0.0);
        for (int color = 0; color < impl_->colorCount; ++color) {
            const std::vector<int>& columns =
                columnsByColor[static_cast<std::size_t>(color)];
            for (int column : columns) {
                scalarProbe[static_cast<std::size_t>(column)] = 1.0;
            }
            std::vector<double> scalarResponse;
            applyExactSchur(scalarProbe, scalarResponse);
            ++impl_->validationApplies;
            for (int column : columns) {
                scalarDiagonal[static_cast<std::size_t>(column)] =
                    scalarResponse[static_cast<std::size_t>(column)];
                scalarSupport[static_cast<std::size_t>(column)] = color;
                for (int offset = graph.rowPtr[static_cast<std::size_t>(column)];
                     offset < graph.rowPtr[static_cast<std::size_t>(column + 1)]; ++offset) {
                    const int row = graph.colInd[static_cast<std::size_t>(offset)];
                    scalarSupport[static_cast<std::size_t>(row)] = color;
                    scalarEdges[static_cast<std::size_t>(
                        graph.edgeId[static_cast<std::size_t>(offset)])] +=
                        0.5 * scalarResponse[static_cast<std::size_t>(row)];
                }
            }
            for (int row = 0; row < count; ++row) {
                if (scalarSupport[static_cast<std::size_t>(row)] != color) {
                    scalarOmitted[static_cast<std::size_t>(row)] +=
                        std::abs(scalarResponse[static_cast<std::size_t>(row)]);
                }
            }
            for (int column : columns) {
                scalarProbe[static_cast<std::size_t>(column)] = 0.0;
            }
        }

        double referenceMaximum = 0.0;
        auto compareValues = [&](const std::vector<double>& blocked,
                                 const std::vector<double>& scalar) {
            for (std::size_t i = 0; i < blocked.size(); ++i) {
                impl_->blockMaximumDifference = std::max(
                    impl_->blockMaximumDifference,
                    std::abs(blocked[i] - scalar[i]));
                referenceMaximum = std::max(referenceMaximum, std::abs(scalar[i]));
            }
        };
        compareValues(diagonal, scalarDiagonal);
        compareValues(edgeValues, scalarEdges);
        compareValues(omittedAbsolute, scalarOmitted);
        impl_->blockRelativeDifference = impl_->blockMaximumDifference
            / std::max(1.0e-300, referenceMaximum);
        if (!(impl_->blockRelativeDifference <= 1.0e-10)) {
            throw std::runtime_error(
                "[Schur proxy] Multi-RHS proxy differs from scalar probing beyond tolerance: max="
                + std::to_string(impl_->blockMaximumDifference)
                + ", relative=" + std::to_string(impl_->blockRelativeDifference) + ".");
        }
    }

    double meanPositiveDiagonal = 0.0;
    int positiveDiagonalCount = 0;
    double minimumDiagonal = std::numeric_limits<double>::max();
    for (double value : diagonal) {
        minimumDiagonal = std::min(minimumDiagonal, value);
        if (value > 0.0 && std::isfinite(value)) {
            meanPositiveDiagonal += value;
            ++positiveDiagonalCount;
        }
    }
    meanPositiveDiagonal = positiveDiagonalCount > 0
        ? meanPositiveDiagonal / positiveDiagonalCount : 1.0;
    const double safety = std::max(1.0e-30, 1.0e-12 * meanPositiveDiagonal);
    if (!(minimumDiagonal > 0.0) || !std::isfinite(minimumDiagonal)) {
        impl_->globalDiagonalShift = -std::min(0.0, minimumDiagonal) + safety;
        for (double& value : diagonal) {
            value += impl_->globalDiagonalShift;
        }
    }

    for (int row = 0; row < count; ++row) {
        diagonal[static_cast<std::size_t>(row)] +=
            omittedAbsolute[static_cast<std::size_t>(row)];
        impl_->totalDiagonalCompensation +=
            omittedAbsolute[static_cast<std::size_t>(row)];
    }

    impl_->minimumRayleigh = std::numeric_limits<double>::max();
    for (int trial = 0; trial < 8; ++trial) {
        std::vector<double> vector(static_cast<std::size_t>(count), 0.0);
        for (int row = 0; row < count; ++row) {
            vector[static_cast<std::size_t>(row)] = deterministicCoefficient(row, trial);
        }
        long double rayleigh = 0.0;
        long double normSquared = 0.0;
        for (int row = 0; row < count; ++row) {
            const double value = vector[static_cast<std::size_t>(row)];
            rayleigh += diagonal[static_cast<std::size_t>(row)] * value * value;
            normSquared += value * value;
        }
        for (std::size_t edge = 0; edge < graph.edges.size(); ++edge) {
            rayleigh += 2.0 * edgeValues[edge]
                * vector[static_cast<std::size_t>(graph.edges[edge].first)]
                * vector[static_cast<std::size_t>(graph.edges[edge].second)];
        }
        impl_->minimumRayleigh = std::min(
            impl_->minimumRayleigh,
            static_cast<double>(rayleigh / normSquared));
    }
    if (!(impl_->minimumRayleigh > 0.0)) {
        const double additionalShift = -impl_->minimumRayleigh + safety;
        impl_->globalDiagonalShift += additionalShift;
        for (double& value : diagonal) {
            value += additionalShift;
        }
        impl_->minimumRayleigh += additionalShift;
        if (!(impl_->minimumRayleigh > 0.0)) {
            throw std::runtime_error(
                "[Schur proxy] Stabilized proxy failed the positive Rayleigh test.");
        }
    }

    // Bitwise fingerprint is useful for reproducibility diagnostics.  The
    // explicit scalar comparison above is the numerical equivalence gate,
    // because MKL may round multi-RHS solves differently in the last bits.
    std::uint64_t valueHash = UINT64_C(1469598103934665603);
    auto appendHash = [&](double value) {
        std::uint64_t bits = 0;
        static_assert(sizeof(bits) == sizeof(value), "Unexpected double size.");
        std::memcpy(&bits, &value, sizeof(bits));
        for (int byte = 0; byte < 8; ++byte) {
            valueHash ^= (bits >> (8 * byte)) & UINT64_C(0xff);
            valueHash *= UINT64_C(1099511628211);
        }
    };
    for (double value : diagonal) {
        appendHash(value);
    }
    for (double value : edgeValues) {
        appendHash(value);
    }
    impl_->matrixValueHash = valueHash;
    if (cacheEnabled && !cachePath.empty()) {
        ProxyCacheData outputCache;
        outputCache.fingerprint = operatorFingerprint;
        outputCache.interfaceDofs = count;
        outputCache.graph = graph;
        outputCache.colors = colors;
        outputCache.diagonal = diagonal;
        outputCache.edgeValues = edgeValues;
        outputCache.minimumRayleigh = impl_->minimumRayleigh;
        outputCache.diagonalShift = impl_->globalDiagonalShift;
        outputCache.diagonalCompensation = impl_->totalDiagonalCompensation;
        outputCache.valueHash = impl_->matrixValueHash;
        saveProxyCache(cachePath, outputCache);
    }
    }

    impl_->nonzeros = static_cast<std::size_t>(count) + 2 * graph.edges.size();
    const long double dense = static_cast<long double>(count)
        * static_cast<long double>(count);
    impl_->matrixDensity = dense > 0.0
        ? static_cast<double>(impl_->nonzeros / dense) : 0.0;
    std::uint64_t factorFingerprint = operatorFingerprint;
    appendHashValue(factorFingerprint, impl_->matrixValueHash);
    appendHashValue(factorFingerprint, std::max(1, requestedPardisoThreads));
    if (cacheEnabled) {
        std::lock_guard<std::mutex> guard(proxyFactorCacheMutex);
        if (proxyFactorCache
            && proxyFactorCache->fingerprint == factorFingerprint
            && proxyFactorCache->solver) {
            impl_->solver = proxyFactorCache->solver;
            impl_->factorCacheHit = true;
            impl_->solver->resetRuntimeCounters();
        } else {
            impl_->solver = std::make_shared<ProxyPardiso>(
                count, graph.edges, edgeValues, diagonal,
                requestedPardisoThreads);
            proxyFactorCache = std::make_shared<CachedProxyFactor>();
            proxyFactorCache->fingerprint = factorFingerprint;
            proxyFactorCache->solver = impl_->solver;
        }
    } else {
        impl_->solver = std::make_shared<ProxyPardiso>(
            count, graph.edges, edgeValues, diagonal,
            requestedPardisoThreads);
    }
    impl_->setupTime = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - setupStart).count();

    if (!outputDirectory.empty()) {
        std::filesystem::create_directories(outputDirectory);
        std::ofstream out(outputDirectory / "schur_proxy_matrix_diagnostics.csv");
        out << "ring,interface_dofs,proxy_nnz,proxy_density,number_of_colors,"
            << "probing_schur_applies,probing_block_size,probing_block_calls,"
            << "validation_schur_applies,block_maximum_difference,block_relative_difference,"
            << "proxy_symbolic_calls,proxy_numerical_calls,"
            << "proxy_symbolic_seconds,proxy_numerical_seconds,symmetry_error,"
            << "minimum_test_rayleigh,diagonal_shift,diagonal_compensation,"
            << "proxy_memory_bytes,proxy_value_hash,matrix_cache_hit,factor_cache_hit,setup_seconds\n";
        out << std::setprecision(16)
            << impl_->ring << ',' << count << ',' << impl_->nonzeros << ','
            << impl_->matrixDensity << ',' << impl_->colorCount << ','
            << impl_->probingApplies << ',' << impl_->probingBlockSize << ','
            << impl_->probingBlockCalls << ',' << impl_->validationApplies << ','
            << impl_->blockMaximumDifference << ','
            << impl_->blockRelativeDifference << ','
            << impl_->solver->symbolicCalls() << ','
            << impl_->solver->numericalCalls() << ','
            << impl_->solver->symbolicSeconds() << ','
            << impl_->solver->numericalSeconds() << ','
            << impl_->matrixSymmetryError << ',' << impl_->minimumRayleigh << ','
            << impl_->globalDiagonalShift << ',' << impl_->totalDiagonalCompensation << ','
            << memoryBytes() << ',' << impl_->matrixValueHash << ','
            << (impl_->matrixCacheHit ? 1 : 0) << ','
            << (impl_->factorCacheHit ? 1 : 0) << ','
            << impl_->setupTime << '\n';
    }
    std::cout << "[Schur proxy] factorized ring-1 proxy: nnz=" << impl_->nonzeros
              << ", colors=" << impl_->colorCount
              << ", block=" << impl_->probingBlockSize
              << ", block calls=" << impl_->probingBlockCalls
              << ", matrix cache=" << (impl_->matrixCacheHit ? "hit" : "miss")
              << ", factor cache=" << (impl_->factorCacheHit ? "hit" : "miss")
              << ", symbolic=" << impl_->solver->symbolicSeconds() << " s"
              << ", numerical=" << impl_->solver->numericalSeconds() << " s"
              << ", memory=" << static_cast<double>(memoryBytes()) / (1024.0 * 1024.0)
              << " MiB\n";
}

SchurProxyPreconditioner::~SchurProxyPreconditioner() = default;
SchurProxyPreconditioner::SchurProxyPreconditioner(
    SchurProxyPreconditioner&&) noexcept = default;
SchurProxyPreconditioner& SchurProxyPreconditioner::operator=(
    SchurProxyPreconditioner&&) noexcept = default;

void SchurProxyPreconditioner::solve(
    const std::vector<double>& rhs, std::vector<double>& solution)
{
    impl_->solver->solve(rhs, solution);
}
void SchurProxyPreconditioner::resetRuntimeCounters()
{
    impl_->solver->resetRuntimeCounters();
}
int SchurProxyPreconditioner::ring() const { return impl_->ring; }
std::size_t SchurProxyPreconditioner::nnz() const { return impl_->nonzeros; }
double SchurProxyPreconditioner::density() const { return impl_->matrixDensity; }
int SchurProxyPreconditioner::colors() const { return impl_->colorCount; }
int SchurProxyPreconditioner::probingSchurApplies() const { return impl_->probingApplies; }
int SchurProxyPreconditioner::probingBlockSize() const { return impl_->probingBlockSize; }
int SchurProxyPreconditioner::probingBlockCalls() const { return impl_->probingBlockCalls; }
int SchurProxyPreconditioner::validationSchurApplies() const { return impl_->validationApplies; }
double SchurProxyPreconditioner::setupSeconds() const { return impl_->setupTime; }
double SchurProxyPreconditioner::symbolicSeconds() const
{
    return impl_->factorCacheHit ? 0.0 : impl_->solver->symbolicSeconds();
}
double SchurProxyPreconditioner::numericalSeconds() const
{
    return impl_->factorCacheHit ? 0.0 : impl_->solver->numericalSeconds();
}
double SchurProxyPreconditioner::solveSeconds() const { return impl_->solver->solveSeconds(); }
int SchurProxyPreconditioner::solveCalls() const { return impl_->solver->solveCalls(); }
int SchurProxyPreconditioner::symbolicCalls() const
{
    return impl_->factorCacheHit ? 0 : impl_->solver->symbolicCalls();
}
int SchurProxyPreconditioner::numericalCalls() const
{
    return impl_->factorCacheHit ? 0 : impl_->solver->numericalCalls();
}
double SchurProxyPreconditioner::symmetryError() const { return impl_->matrixSymmetryError; }
double SchurProxyPreconditioner::minimumTestRayleigh() const { return impl_->minimumRayleigh; }
double SchurProxyPreconditioner::diagonalShift() const { return impl_->globalDiagonalShift; }
double SchurProxyPreconditioner::diagonalCompensation() const { return impl_->totalDiagonalCompensation; }
std::uint64_t SchurProxyPreconditioner::valueHash() const { return impl_->matrixValueHash; }
double SchurProxyPreconditioner::blockMaximumDifference() const { return impl_->blockMaximumDifference; }
double SchurProxyPreconditioner::blockRelativeDifference() const { return impl_->blockRelativeDifference; }
std::size_t SchurProxyPreconditioner::memoryBytes() const
{
    return impl_->graphMemoryBytes + (impl_->solver ? impl_->solver->memoryBytes() : 0);
}
bool SchurProxyPreconditioner::matrixCacheHit() const { return impl_->matrixCacheHit; }
bool SchurProxyPreconditioner::factorCacheHit() const { return impl_->factorCacheHit; }

ProxyDiagnosticsResult runSchurProxyDiagnostics(
    const Mesh& mesh,
    const CaseConfig& physics,
    const InterfacePartition& partition,
    double highConductivityThreshold,
    bool useMaterialConnectivity,
    int requestedProbeColumns,
    const std::filesystem::path& outputDirectory,
    const ExactSchurApply& applyExactSchur)
{
    const auto start = std::chrono::steady_clock::now();
    ProxyDiagnosticsResult result;
    result.interfaceDofs = static_cast<int>(partition.interfaceGlobalDofs.size());
    if (result.interfaceDofs == 0) {
        return result;
    }
    if (requestedProbeColumns <= 0) {
        throw std::runtime_error("[Schur proxy] Probe-column count must be positive.");
    }

    GraphCsr graph = buildInterfaceGraph(
        mesh, physics, partition, highConductivityThreshold,
        useMaterialConnectivity, result.materialStrongEdges,
        result.crossDomainHighKComponents);
    result.graphEdges = graph.undirectedEdges;

    std::vector<double> interfaceConductivity(
        static_cast<std::size_t>(result.interfaceDofs), 0.0);
    for (const Tet& tet : mesh.tets) {
        const double conductivity = materialConductivity(materialForTet(physics, tet));
        for (int globalDof : tet.dof) {
            if (globalDof < 0
                || globalDof >= static_cast<int>(partition.globalToInterface.size())) {
                continue;
            }
            const int gamma = partition.globalToInterface[static_cast<std::size_t>(globalDof)];
            if (gamma >= 0) {
                interfaceConductivity[static_cast<std::size_t>(gamma)] = std::max(
                    interfaceConductivity[static_cast<std::size_t>(gamma)], conductivity);
            }
        }
    }

    std::vector<int> ordered(static_cast<std::size_t>(result.interfaceDofs));
    std::iota(ordered.begin(), ordered.end(), 0);
    std::sort(ordered.begin(), ordered.end(), [&](int left, int right) {
        const int leftGlobal = partition.interfaceGlobalDofs[static_cast<std::size_t>(left)];
        const int rightGlobal = partition.interfaceGlobalDofs[static_cast<std::size_t>(right)];
        const Node& a = mesh.nodes[static_cast<std::size_t>(leftGlobal)];
        const Node& b = mesh.nodes[static_cast<std::size_t>(rightGlobal)];
        return std::tie(interfaceConductivity[static_cast<std::size_t>(left)],
                        a.subdomain, a.p.z, a.p.y, a.p.x, left)
             < std::tie(interfaceConductivity[static_cast<std::size_t>(right)],
                        b.subdomain, b.p.z, b.p.y, b.p.x, right);
    });
    const int probeCount = std::min(requestedProbeColumns, result.interfaceDofs);
    std::vector<int> probes;
    probes.reserve(static_cast<std::size_t>(probeCount));
    for (int probe = 0; probe < probeCount; ++probe) {
        const std::size_t position = probeCount == 1 ? 0
            : static_cast<std::size_t>(
                (static_cast<long double>(probe) * (ordered.size() - 1))
                / static_cast<long double>(probeCount - 1));
        probes.push_back(ordered[position]);
    }
    std::sort(probes.begin(), probes.end());
    probes.erase(std::unique(probes.begin(), probes.end()), probes.end());
    result.probeColumns = static_cast<int>(probes.size());

    std::vector<std::vector<unsigned char>> distances;
    std::vector<std::vector<double>> columns;
    distances.reserve(probes.size());
    columns.reserve(probes.size());
    for (int probe : probes) {
        distances.push_back(distancesToRing3(graph, probe));
        std::vector<double> unit(static_cast<std::size_t>(result.interfaceDofs), 0.0);
        unit[static_cast<std::size_t>(probe)] = 1.0;
        std::vector<double> image;
        applyExactSchur(unit, image);
        columns.push_back(std::move(image));
        ++result.exactSchurApplies;
    }

    std::array<long double, 5> absoluteSum{};
    std::array<long double, 5> nonzeroAbsoluteSum{};
    std::array<long double, 5> energy{};
    std::array<double, 5> maximum{};
    std::array<std::size_t, 5> count{};
    std::array<std::size_t, 5> nonzeroCount{};
    std::array<long double, 3> neighborhoodSum{};
    constexpr double nonzeroTolerance = 1.0e-30;
    for (std::size_t probe = 0; probe < probes.size(); ++probe) {
        for (int row = 0; row < result.interfaceDofs; ++row) {
            const unsigned char distance = distances[probe][static_cast<std::size_t>(row)];
            const int bin = std::min(4, static_cast<int>(distance));
            const double absolute = std::abs(columns[probe][static_cast<std::size_t>(row)]);
            absoluteSum[static_cast<std::size_t>(bin)] += absolute;
            energy[static_cast<std::size_t>(bin)] += absolute * absolute;
            maximum[static_cast<std::size_t>(bin)] = std::max(
                maximum[static_cast<std::size_t>(bin)], absolute);
            ++count[static_cast<std::size_t>(bin)];
            if (absolute > nonzeroTolerance) {
                nonzeroAbsoluteSum[static_cast<std::size_t>(bin)] += absolute;
                ++nonzeroCount[static_cast<std::size_t>(bin)];
            }
            for (int ring = 1; ring <= 3; ++ring) {
                if (distance <= ring) {
                    neighborhoodSum[static_cast<std::size_t>(ring - 1)] += 1.0;
                }
            }
        }
    }
    const long double totalEnergy = std::accumulate(energy.begin(), energy.end(), 0.0L);
    for (int bin = 0; bin < 5; ++bin) {
        ProxyDistanceMetric metric;
        metric.distance = bin;
        metric.entryCount = count[static_cast<std::size_t>(bin)];
        metric.nonzeroCount = nonzeroCount[static_cast<std::size_t>(bin)];
        metric.meanAbsolute = metric.entryCount == 0 ? 0.0
            : static_cast<double>(absoluteSum[static_cast<std::size_t>(bin)]
                / static_cast<long double>(metric.entryCount));
        metric.meanNonzeroAbsolute = metric.nonzeroCount == 0 ? 0.0
            : static_cast<double>(nonzeroAbsoluteSum[static_cast<std::size_t>(bin)]
                / static_cast<long double>(metric.nonzeroCount));
        metric.maximumAbsolute = maximum[static_cast<std::size_t>(bin)];
        metric.energyFraction = totalEnergy > 0.0
            ? static_cast<double>(energy[static_cast<std::size_t>(bin)] / totalEnergy)
            : 0.0;
        result.distanceMetrics.push_back(metric);
    }

    constexpr int randomTrials = 3;
    for (int ring = 1; ring <= 3; ++ring) {
        ProxyRingMetric metric;
        metric.ring = ring;
        metric.meanColumnNnz = probes.empty() ? 0.0
            : static_cast<double>(neighborhoodSum[static_cast<std::size_t>(ring - 1)]
                / static_cast<long double>(probes.size()));
        const long double estimated = metric.meanColumnNnz
            * static_cast<long double>(result.interfaceDofs);
        const long double denseLimit = static_cast<long double>(result.interfaceDofs)
            * static_cast<long double>(result.interfaceDofs);
        metric.estimatedNnz = static_cast<std::size_t>(std::min(estimated, denseLimit));
        metric.estimatedDensity = denseLimit > 0.0
            ? static_cast<double>(estimated / denseLimit) : 0.0;
        long double coveredEnergy = 0.0;
        for (int distance = 0; distance <= ring; ++distance) {
            coveredEnergy += energy[static_cast<std::size_t>(distance)];
        }
        metric.energyCoverage = totalEnergy > 0.0
            ? static_cast<double>(coveredEnergy / totalEnergy) : 0.0;
        metric.csrMemoryEstimateBytes =
            static_cast<std::size_t>(result.interfaceDofs + 1) * sizeof(int)
            + metric.estimatedNnz * (sizeof(int) + sizeof(double));

        double errorSum = 0.0;
        for (int trial = 0; trial < randomTrials; ++trial) {
            long double exactNormSquared = 0.0;
            long double errorNormSquared = 0.0;
            for (int row = 0; row < result.interfaceDofs; ++row) {
                double exact = 0.0;
                double truncated = 0.0;
                for (std::size_t probe = 0; probe < probes.size(); ++probe) {
                    const double contribution = deterministicCoefficient(
                        static_cast<int>(probe), trial)
                        * columns[probe][static_cast<std::size_t>(row)];
                    exact += contribution;
                    if (distances[probe][static_cast<std::size_t>(row)] <= ring) {
                        truncated += contribution;
                    }
                }
                exactNormSquared += exact * exact;
                const double difference = exact - truncated;
                errorNormSquared += difference * difference;
            }
            const double relativeError = exactNormSquared > 0.0
                ? std::sqrt(static_cast<double>(errorNormSquared / exactNormSquared))
                : 0.0;
            errorSum += relativeError;
            metric.randomOperatorErrorMaximum = std::max(
                metric.randomOperatorErrorMaximum, relativeError);
        }
        metric.randomOperatorErrorMean = errorSum / randomTrials;
        result.ringMetrics.push_back(metric);
    }

    // "Recommended" means that at least the cheapest candidate is viable.
    // Larger rings remain subject to their own fill/memory stop conditions.
    const ProxyRingMetric& ring1 = result.ringMetrics.front();
    constexpr std::size_t maximumRecommendedCsrBytes =
        static_cast<std::size_t>(8) * 1024 * 1024 * 1024;
    result.proxyRecommended = ring1.energyCoverage >= 0.90
        && ring1.randomOperatorErrorMaximum <= 0.25
        && ring1.csrMemoryEstimateBytes <= maximumRecommendedCsrBytes;
    result.setupSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();

    if (!outputDirectory.empty()) {
        std::filesystem::create_directories(outputDirectory);
        std::ofstream distanceOut(outputDirectory / "schur_proxy_distance_decay.csv");
        distanceOut << "distance_bin,entry_count,nonzero_count,mean_abs_Sij,"
                    << "mean_nonzero_abs_Sij,max_abs_Sij,energy_fraction\n";
        distanceOut << std::setprecision(16);
        for (const ProxyDistanceMetric& metric : result.distanceMetrics) {
            distanceOut << (metric.distance < 4 ? std::to_string(metric.distance) : "gt3") << ','
                        << metric.entryCount << ',' << metric.nonzeroCount << ','
                        << metric.meanAbsolute << ',' << metric.meanNonzeroAbsolute << ','
                        << metric.maximumAbsolute << ',' << metric.energyFraction << '\n';
        }

        std::ofstream ringOut(outputDirectory / "schur_proxy_ring_diagnostics.csv");
        ringOut << "ring,interface_dofs,graph_edges,material_strong_edges,"
                << "cross_domain_high_k_components,probe_columns,exact_schur_applies,"
                << "mean_column_nnz,estimated_nnz,estimated_density,energy_coverage,"
                << "sampled_random_operator_error_mean,sampled_random_operator_error_max,"
                << "csr_memory_estimate_bytes,proxy_recommended,diagnostic_seconds\n";
        ringOut << std::setprecision(16);
        for (const ProxyRingMetric& metric : result.ringMetrics) {
            ringOut << metric.ring << ',' << result.interfaceDofs << ','
                    << result.graphEdges << ',' << result.materialStrongEdges << ','
                    << result.crossDomainHighKComponents << ',' << result.probeColumns << ','
                    << result.exactSchurApplies << ',' << metric.meanColumnNnz << ','
                    << metric.estimatedNnz << ',' << metric.estimatedDensity << ','
                    << metric.energyCoverage << ',' << metric.randomOperatorErrorMean << ','
                    << metric.randomOperatorErrorMaximum << ','
                    << metric.csrMemoryEstimateBytes << ','
                    << (result.proxyRecommended ? 1 : 0) << ',' << result.setupSeconds << '\n';
        }

        std::ofstream probeOut(outputDirectory / "schur_proxy_probe_columns.csv");
        probeOut << "probe,gamma,global_dof,subdomain,x,y,z,conductivity,"
                 << "ring1_nnz,ring2_nnz,ring3_nnz,column_l2\n";
        probeOut << std::setprecision(16);
        for (std::size_t probe = 0; probe < probes.size(); ++probe) {
            const int gamma = probes[probe];
            const int globalDof = partition.interfaceGlobalDofs[static_cast<std::size_t>(gamma)];
            const Node& node = mesh.nodes[static_cast<std::size_t>(globalDof)];
            std::array<int, 3> ringNnz{};
            long double normSquared = 0.0;
            for (int row = 0; row < result.interfaceDofs; ++row) {
                normSquared += columns[probe][static_cast<std::size_t>(row)]
                    * columns[probe][static_cast<std::size_t>(row)];
                for (int ring = 1; ring <= 3; ++ring) {
                    if (distances[probe][static_cast<std::size_t>(row)] <= ring) {
                        ++ringNnz[static_cast<std::size_t>(ring - 1)];
                    }
                }
            }
            probeOut << probe << ',' << gamma << ',' << globalDof << ','
                     << node.subdomain << ',' << node.p.x << ',' << node.p.y << ','
                     << node.p.z << ',' << interfaceConductivity[static_cast<std::size_t>(gamma)]
                     << ',' << ringNnz[0] << ',' << ringNnz[1] << ',' << ringNnz[2]
                     << ',' << std::sqrt(static_cast<double>(normSquared)) << '\n';
        }
    }

    std::cout << "[Schur proxy] interface graph: dofs=" << result.interfaceDofs
              << ", edges=" << result.graphEdges
              << ", material strong edges=" << result.materialStrongEdges
              << ", high-k cross-domain components="
              << result.crossDomainHighKComponents << '\n';
    for (const ProxyRingMetric& metric : result.ringMetrics) {
        std::cout << "[Schur proxy] ring=" << metric.ring
                  << ", estimated nnz=" << metric.estimatedNnz
                  << ", coverage=" << metric.energyCoverage
                  << ", sampled random error=" << metric.randomOperatorErrorMaximum
                  << ", CSR memory="
                  << static_cast<double>(metric.csrMemoryEstimateBytes) / (1024.0 * 1024.0)
                  << " MiB\n";
    }
    std::cout << "[Schur proxy] construction recommendation: "
              << (result.proxyRecommended ? "continue" : "stop after locality diagnostics")
              << '\n';
    return result;
}

} // namespace ddm_schur
