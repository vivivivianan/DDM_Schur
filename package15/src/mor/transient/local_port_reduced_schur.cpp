#include "local_port_reduced_schur.hpp"

#include "sipg_core.hpp"
#include "config_io.hpp"
#include "ddm_schur/interface_operator.hpp"
#include "global_interface_coarse.hpp"
#include "mesh_loader.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

namespace mor::transient {
namespace {

using Clock = std::chrono::steady_clock;

double secondsSince(const Clock::time_point& start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

double dot(const double* left, const double* right, int rows)
{
    long double value = 0.0L;
    for (int row = 0; row < rows; ++row) {
        value += static_cast<long double>(left[row]) * right[row];
    }
    return static_cast<double>(value);
}

double squaredNorm(const std::vector<double>& values)
{
    long double value = 0.0L;
    for (double entry : values) value += static_cast<long double>(entry) * entry;
    return static_cast<double>(value);
}

template <typename T>
void hashValue(std::uint64_t& hash, const T& value)
{
    const unsigned char* bytes = reinterpret_cast<const unsigned char*>(&value);
    for (std::size_t byte = 0; byte < sizeof(T); ++byte) {
        hash ^= bytes[byte];
        hash *= UINT64_C(1099511628211);
    }
}

std::map<int, int> readRankFile(const std::filesystem::path& path)
{
    std::map<int, int> ranks;
    if (path.empty()) return ranks;
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("[Local port] Cannot open rank file: " + path.string());
    }
    std::string line;
    std::getline(input, line);
    while (std::getline(input, line)) {
        const std::size_t comma = line.find(',');
        if (comma == std::string::npos) continue;
        const int interfaceId = std::stoi(line.substr(0, comma));
        const int rank = std::stoi(line.substr(comma + 1));
        if (rank <= 0) {
            throw std::runtime_error("[Local port] Rank-file values must be positive.");
        }
        ranks[interfaceId] = rank;
    }
    return ranks;
}

struct PhysicalSupport {
    int interfaceId = -1;
    int left = -1;
    int right = -1;
    std::vector<int> gamma;
};

std::vector<PhysicalSupport> physicalSupports(
    const Mesh& mesh,
    const ddm_schur::InterfacePartition& partition)
{
    std::map<std::pair<int, int>, int> summaryIds;
    for (std::size_t interfaceIndex = 0;
         interfaceIndex < mesh.interfaceSummaries.size(); ++interfaceIndex) {
        const auto& summary = mesh.interfaceSummaries[interfaceIndex];
        summaryIds[std::minmax(summary.leftSubdomain, summary.rightSubdomain)] =
            static_cast<int>(interfaceIndex);
    }
    std::map<std::pair<int, int>, std::set<int>> indices;
    for (const auto& domain : partition.domains) {
        for (const auto& neighbor : domain.interfaceGlobalDofsByNeighbor) {
            const auto pair = std::minmax(domain.domainId, neighbor.first);
            std::set<int>& support = indices[pair];
            for (int global : neighbor.second) {
                const int gamma = partition.globalToInterface[static_cast<std::size_t>(global)];
                if (gamma >= 0) support.insert(gamma);
            }
        }
    }
    std::vector<PhysicalSupport> result;
    int fallbackId = static_cast<int>(mesh.interfaceSummaries.size());
    for (const auto& entry : indices) {
        PhysicalSupport support;
        const auto found = summaryIds.find(entry.first);
        support.interfaceId = found == summaryIds.end() ? fallbackId++ : found->second;
        support.left = entry.first.first;
        support.right = entry.first.second;
        support.gamma.assign(entry.second.begin(), entry.second.end());
        if (!support.gamma.empty()) result.push_back(std::move(support));
    }
    std::sort(result.begin(), result.end(), [](const auto& left, const auto& right) {
        return left.interfaceId < right.interfaceId;
    });
    // Interface nodes on physical-port junctions can be reported by more than
    // one neighbor patch.  A block-diagonal port map must nevertheless own
    // every full-interface coordinate exactly once; otherwise lift() sums
    // duplicate port contributions while restrict() silently treats them as
    // independent coordinates.  Give a junction coordinate to the lowest
    // deterministic physical-interface id.  This is a discrete partition of
    // the existing SIPG/mortar trace space and does not alter its operators.
    std::vector<int> owner(partition.interfaceGlobalDofs.size(), -1);
    for (PhysicalSupport& support : result) {
        std::vector<int> disjoint;
        disjoint.reserve(support.gamma.size());
        for (int gamma : support.gamma) {
            if (owner[static_cast<std::size_t>(gamma)] < 0) {
                owner[static_cast<std::size_t>(gamma)] = support.interfaceId;
                disjoint.push_back(gamma);
            }
        }
        support.gamma = std::move(disjoint);
    }
    result.erase(std::remove_if(result.begin(), result.end(),
        [](const PhysicalSupport& support) { return support.gamma.empty(); }),
        result.end());
    if (std::find(owner.begin(), owner.end(), -1) != owner.end()) {
        throw std::runtime_error(
            "[Local port] Physical-interface supports do not cover the full SIPG trace.");
    }
    return result;
}

double familyScale(const std::vector<std::vector<double>>& family,
                   const std::vector<int>& globalDofs)
{
    long double energy = 0.0L;
    std::size_t entries = 0;
    for (const auto& snapshot : family) {
        for (int global : globalDofs) {
            if (global < 0 || global >= static_cast<int>(snapshot.size())) continue;
            const double value = snapshot[static_cast<std::size_t>(global)];
            energy += static_cast<long double>(value) * value;
            ++entries;
        }
    }
    return entries == 0 ? 1.0
        : std::max(1.0e-300, std::sqrt(static_cast<double>(energy / entries)));
}

void appendCandidate(LocalPortBasis& port,
                     std::vector<double> candidate,
                     int maximumRank,
                     double relativeTolerance,
                     double& totalEnergy,
                     double& retainedEnergy)
{
    ++port.candidateColumns;
    const double originalEnergy = squaredNorm(candidate);
    totalEnergy += originalEnergy;
    if (!(originalEnergy > 0.0) || port.rank >= maximumRank) return;
    for (int pass = 0; pass < 2; ++pass) {
        for (int mode = 0; mode < port.rank; ++mode) {
            const double* basis = port.basis.data()
                + static_cast<std::size_t>(mode * port.rows);
            const double coefficient = dot(basis, candidate.data(), port.rows);
            for (int row = 0; row < port.rows; ++row) {
                candidate[static_cast<std::size_t>(row)] -= coefficient * basis[row];
            }
        }
    }
    const double norm = std::sqrt(squaredNorm(candidate));
    if (!(norm > relativeTolerance * std::sqrt(originalEnergy))) return;
    for (double& value : candidate) value /= norm;
    port.basis.insert(port.basis.end(), candidate.begin(), candidate.end());
    ++port.rank;
    ++port.acceptedColumns;
    retainedEnergy += norm * norm;
}

void writeRaw(std::ofstream& out, const void* data, std::size_t bytes)
{
    out.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(bytes));
    if (!out) throw std::runtime_error("[Local port] Model write failed.");
}

void readRaw(std::ifstream& in, void* data, std::size_t bytes)
{
    in.read(reinterpret_cast<char*>(data), static_cast<std::streamsize>(bytes));
    if (!in) throw std::runtime_error("[Local port] Model read failed.");
}

template <typename T>
void writeScalar(std::ofstream& out, const T& value)
{
    writeRaw(out, &value, sizeof(T));
}

template <typename T>
T readScalar(std::ifstream& in)
{
    T value{};
    readRaw(in, &value, sizeof(T));
    return value;
}

template <typename T>
void writeVector(std::ofstream& out, const std::vector<T>& values)
{
    const std::uint64_t size = static_cast<std::uint64_t>(values.size());
    writeScalar(out, size);
    if (!values.empty()) writeRaw(out, values.data(), values.size() * sizeof(T));
}

template <typename T>
std::vector<T> readVector(std::ifstream& in)
{
    const std::uint64_t size = readScalar<std::uint64_t>(in);
    if (size > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error("[Local port] Serialized vector is too large.");
    }
    std::vector<T> values(static_cast<std::size_t>(size));
    if (!values.empty()) readRaw(in, values.data(), values.size() * sizeof(T));
    return values;
}

void writeString(std::ofstream& out, const std::string& value)
{
    const std::uint64_t size = static_cast<std::uint64_t>(value.size());
    writeScalar(out, size);
    if (size > 0) writeRaw(out, value.data(), static_cast<std::size_t>(size));
}

std::string readString(std::ifstream& in)
{
    const std::uint64_t size = readScalar<std::uint64_t>(in);
    if (size > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error("[Local port] Serialized string is too large.");
    }
    std::string value(static_cast<std::size_t>(size), '\0');
    if (!value.empty()) readRaw(in, value.data(), value.size());
    return value;
}

} // namespace

LocalPortModel buildLocalPortModel(
    const Mesh& mesh,
    const CaseConfig&,
    const ddm_schur::InterfacePartition& partition,
    const LocalPortSnapshotFamilies& snapshots,
    const LocalPortOptions& options)
{
    if (options.requestedRank < 0 || options.discardedEnergyTolerance < 0.0
        || !(options.relativeTolerance > 0.0)) {
        throw std::runtime_error("[Local port] Invalid rank/tolerance options.");
    }
    const auto totalStart = Clock::now();
    LocalPortModel model;
    model.formatVersion = 1;
    model.basisMethod = "port-pod";
    model.fullInterfaceDofs = static_cast<int>(partition.interfaceGlobalDofs.size());
    model.interfaceGlobalDofs = partition.interfaceGlobalDofs;
    const std::map<int, int> rankByInterface = readRankFile(options.rankFile);
    const std::vector<PhysicalSupport> supports = physicalSupports(mesh, partition);
    std::map<std::uint64_t, std::size_t> templateByFingerprint;
    const auto basisStart = Clock::now();
    for (const PhysicalSupport& support : supports) {
        LocalPortBasis port;
        port.interfaceId = support.interfaceId;
        port.leftSubdomain = support.left;
        port.rightSubdomain = support.right;
        port.interfaceIndices = support.gamma;
        port.rows = static_cast<int>(support.gamma.size());
        std::vector<int> globals;
        globals.reserve(support.gamma.size());
        for (int gamma : support.gamma) {
            globals.push_back(partition.interfaceGlobalDofs[static_cast<std::size_t>(gamma)]);
        }
        const auto requested = rankByInterface.find(port.interfaceId);
        const int requestedRank = requested != rankByInterface.end()
            ? requested->second : options.requestedRank;
        const int maximumRank = std::min(port.rows,
            requestedRank > 0 ? requestedRank : port.rows);
        double totalEnergy = 0.0;
        double retainedEnergy = 0.0;

        // A shared physical coordinate drives both nonmatching sides.  These
        // low-order fields make the mortar patch test exact and are followed
        // by temperature, weak-flux, and residual snapshots.
        for (int geometryMode = 0; geometryMode < 4; ++geometryMode) {
            std::vector<double> candidate(static_cast<std::size_t>(port.rows), 1.0);
            for (int row = 0; row < port.rows; ++row) {
                const Vec3& point = mesh.nodes[static_cast<std::size_t>(
                    globals[static_cast<std::size_t>(row)])].p;
                if (geometryMode == 1) candidate[static_cast<std::size_t>(row)] = point.x;
                if (geometryMode == 2) candidate[static_cast<std::size_t>(row)] = point.y;
                if (geometryMode == 3) candidate[static_cast<std::size_t>(row)] = point.z;
            }
            appendCandidate(port, std::move(candidate), maximumRank,
                options.relativeTolerance, totalEnergy, retainedEnergy);
        }

        // Preserve the selected deployment/training trajectory exactly.  The
        // remaining rank is still selected from the joint scale-balanced POD.
        const double mandatoryScale = familyScale(
            snapshots.mandatoryTemperature, globals);
        for (const std::vector<double>& snapshot
             : snapshots.mandatoryTemperature) {
            if (snapshot.size() != mesh.nodes.size()) {
                throw std::runtime_error(
                    "[Local port] Mandatory snapshot size mismatch.");
            }
            std::vector<double> candidate(
                static_cast<std::size_t>(port.rows), 0.0);
            for (int row = 0; row < port.rows; ++row) {
                candidate[static_cast<std::size_t>(row)] =
                    options.temperatureWeight
                    * snapshot[static_cast<std::size_t>(
                        globals[static_cast<std::size_t>(row)])]
                    / mandatoryScale;
            }
            appendCandidate(port, std::move(candidate), maximumRank,
                options.relativeTolerance, totalEnergy, retainedEnergy);
        }

        struct SnapshotFamily {
            const std::vector<std::vector<double>>* values = nullptr;
            double weight = 0.0;
            double scale = 1.0;
        };
        const std::vector<SnapshotFamily> families{
            {&snapshots.temperature, options.temperatureWeight,
                familyScale(snapshots.temperature, globals)},
            {&snapshots.flux, options.fluxWeight,
                familyScale(snapshots.flux, globals)},
            {&snapshots.residual, options.residualWeight,
                familyScale(snapshots.residual, globals)}};
#ifdef USE_MKL_PARDISO
        int snapshotColumns = 0;
        for (const SnapshotFamily& family : families) {
            if (family.weight > 0.0) {
                snapshotColumns += static_cast<int>(family.values->size());
            }
        }
        port.candidateColumns += snapshotColumns;
        if (snapshotColumns == 0) port.retainedEnergy = 1.0;
        if (snapshotColumns > 0 && port.rank < maximumRank) {
            std::vector<double> snapshotMatrix(static_cast<std::size_t>(
                port.rows) * snapshotColumns, 0.0);
            int column = 0;
            for (const SnapshotFamily& family : families) {
                if (!(family.weight > 0.0)) continue;
                for (const std::vector<double>& snapshot : *family.values) {
                    if (snapshot.size() != mesh.nodes.size()) {
                        throw std::runtime_error(
                            "[Local port] Snapshot size mismatch.");
                    }
                    double* candidate = snapshotMatrix.data()
                        + static_cast<std::size_t>(column * port.rows);
                    for (int row = 0; row < port.rows; ++row) {
                        candidate[row] = family.weight
                            * snapshot[static_cast<std::size_t>(
                                globals[static_cast<std::size_t>(row)])]
                            / family.scale;
                    }
                    // Remove the mandatory shared-coordinate patch modes
                    // before the snapshot POD so rank is not spent twice.
                    for (int pass = 0; pass < 2; ++pass) {
                        for (int mode = 0; mode < port.rank; ++mode) {
                            const double* geometry = port.basis.data()
                                + static_cast<std::size_t>(mode * port.rows);
                            const double coefficient = dot(
                                geometry, candidate, port.rows);
                            for (int row = 0; row < port.rows; ++row) {
                                candidate[row] -= coefficient * geometry[row];
                            }
                        }
                    }
                    ++column;
                }
            }
            std::vector<double> gram(static_cast<std::size_t>(
                snapshotColumns) * snapshotColumns, 0.0);
            cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans,
                snapshotColumns, port.rows, 1.0, snapshotMatrix.data(),
                port.rows, 0.0, gram.data(), snapshotColumns);
            for (int gramColumn = 0; gramColumn < snapshotColumns; ++gramColumn) {
                for (int gramRow = gramColumn + 1;
                     gramRow < snapshotColumns; ++gramRow) {
                    gram[static_cast<std::size_t>(
                        gramRow + gramColumn * snapshotColumns)] =
                        gram[static_cast<std::size_t>(
                            gramColumn + gramRow * snapshotColumns)];
                }
            }
            std::vector<double> eigenvalues(
                static_cast<std::size_t>(snapshotColumns), 0.0);
            const lapack_int info = LAPACKE_dsyevd(LAPACK_COL_MAJOR,
                'V', 'U', snapshotColumns, gram.data(), snapshotColumns,
                eigenvalues.data());
            if (info != 0) {
                throw std::runtime_error(
                    "[Local port] Snapshot Gram eigensolve failed for interface "
                    + std::to_string(port.interfaceId) + " with info="
                    + std::to_string(info));
            }
            const double largest = std::max(0.0, eigenvalues.back());
            long double snapshotEnergy = 0.0L;
            for (double value : eigenvalues) {
                snapshotEnergy += std::max(0.0, value);
            }
            int numericalRank = 0;
            for (int mode = snapshotColumns - 1; mode >= 0; --mode) {
                if (eigenvalues[static_cast<std::size_t>(mode)]
                    > options.relativeTolerance * options.relativeTolerance * largest) {
                    ++numericalRank;
                }
            }
            int additionalRank = std::min(maximumRank - port.rank, numericalRank);
            if (requestedRank <= 0 && snapshotEnergy > 0.0L) {
                long double cumulative = 0.0L;
                additionalRank = 0;
                for (int mode = snapshotColumns - 1;
                     mode >= snapshotColumns - numericalRank; --mode) {
                    cumulative += std::max(0.0,
                        eigenvalues[static_cast<std::size_t>(mode)]);
                    ++additionalRank;
                    if (1.0L - cumulative / snapshotEnergy
                        <= options.discardedEnergyTolerance) break;
                }
                additionalRank = std::min(
                    additionalRank, maximumRank - port.rank);
            }
            long double selectedEnergy = 0.0L;
            for (int selectedMode = 0;
                 selectedMode < additionalRank; ++selectedMode) {
                const int eigenMode = snapshotColumns - 1 - selectedMode;
                std::vector<double> candidate(
                    static_cast<std::size_t>(port.rows), 0.0);
                cblas_dgemv(CblasColMajor, CblasNoTrans,
                    port.rows, snapshotColumns, 1.0, snapshotMatrix.data(),
                    port.rows, gram.data()
                        + static_cast<std::size_t>(eigenMode * snapshotColumns),
                    1, 0.0, candidate.data(), 1);
                const double originalNorm = std::sqrt(squaredNorm(candidate));
                for (int pass = 0; pass < 2; ++pass) {
                    for (int mode = 0; mode < port.rank; ++mode) {
                        const double* basis = port.basis.data()
                            + static_cast<std::size_t>(mode * port.rows);
                        const double coefficient = dot(
                            basis, candidate.data(), port.rows);
                        for (int row = 0; row < port.rows; ++row) {
                            candidate[static_cast<std::size_t>(row)] -=
                                coefficient * basis[row];
                        }
                    }
                }
                const double norm = std::sqrt(squaredNorm(candidate));
                if (!(norm > options.relativeTolerance * originalNorm)) continue;
                for (double& value : candidate) value /= norm;
                port.basis.insert(port.basis.end(),
                    candidate.begin(), candidate.end());
                ++port.rank;
                ++port.acceptedColumns;
                selectedEnergy += std::max(0.0,
                    eigenvalues[static_cast<std::size_t>(eigenMode)]);
            }
            port.retainedEnergy = snapshotEnergy > 0.0L
                ? std::min(1.0, static_cast<double>(
                    selectedEnergy / snapshotEnergy)) : 1.0;
        }
#else
        for (const SnapshotFamily& family : families) {
            if (!(family.weight > 0.0)) continue;
            for (const std::vector<double>& snapshot : *family.values) {
                if (snapshot.size() != mesh.nodes.size()) {
                    throw std::runtime_error("[Local port] Snapshot size mismatch.");
                }
                std::vector<double> candidate(
                    static_cast<std::size_t>(port.rows), 0.0);
                for (int row = 0; row < port.rows; ++row) {
                    candidate[static_cast<std::size_t>(row)] = family.weight
                        * snapshot[static_cast<std::size_t>(
                            globals[static_cast<std::size_t>(row)])]
                        / family.scale;
                }
                appendCandidate(port, std::move(candidate), maximumRank,
                    std::max(options.relativeTolerance,
                        std::sqrt(options.discardedEnergyTolerance)),
                    totalEnergy, retainedEnergy);
            }
        }
        port.retainedEnergy = totalEnergy > 0.0
            ? std::min(1.0, retainedEnergy / totalEnergy) : 1.0;
#endif
        if (port.rank <= 0) {
            throw std::runtime_error("[Local port] Empty physical-interface basis.");
        }
        for (int left = 0; left < port.rank; ++left) {
            for (int right = 0; right <= left; ++right) {
                const double product = dot(
                    port.basis.data() + static_cast<std::size_t>(left * port.rows),
                    port.basis.data() + static_cast<std::size_t>(right * port.rows),
                    port.rows);
                port.orthogonalityError = std::max(port.orthogonalityError,
                    std::abs(product - (left == right ? 1.0 : 0.0)));
            }
        }
        std::uint64_t fingerprint = UINT64_C(1469598103934665603);
        hashValue(fingerprint, port.rows);
        hashValue(fingerprint, port.rank);
        Vec3 minimum{
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max()};
        for (int global : globals) {
            const Vec3& point = mesh.nodes[static_cast<std::size_t>(global)].p;
            minimum.x = std::min(minimum.x, point.x);
            minimum.y = std::min(minimum.y, point.y);
            minimum.z = std::min(minimum.z, point.z);
        }
        for (int global : globals) {
            const Node& node = mesh.nodes[static_cast<std::size_t>(global)];
            const double relativeX = node.p.x - minimum.x;
            const double relativeY = node.p.y - minimum.y;
            const double relativeZ = node.p.z - minimum.z;
            const int side = node.subdomain == port.leftSubdomain ? 0 : 1;
            hashValue(fingerprint, relativeX);
            hashValue(fingerprint, relativeY);
            hashValue(fingerprint, relativeZ);
            hashValue(fingerprint, side);
        }
        for (double value : port.basis) hashValue(fingerprint, value);
        port.fingerprint = fingerprint;
        const auto prototype = templateByFingerprint.find(fingerprint);
        if (prototype != templateByFingerprint.end()
            && model.ports[prototype->second].rows == port.rows
            && model.ports[prototype->second].rank == port.rank
            && model.ports[prototype->second].basis == port.basis) {
            port.templateId = model.ports[prototype->second].templateId;
            port.templateReused = true;
        } else {
            port.templateId = static_cast<int>(templateByFingerprint.size());
            templateByFingerprint[fingerprint] = model.ports.size();
        }
        model.reducedInterfaceDofs += port.rank;
        model.ports.push_back(std::move(port));
    }
    model.basisSeconds = secondsSince(basisStart);
    model.snapshotSeconds = secondsSince(totalStart) - model.basisSeconds;
    model.modelBytes = model.interfaceGlobalDofs.capacity() * sizeof(int);
    for (const LocalPortBasis& port : model.ports) {
        model.modelBytes += port.interfaceIndices.capacity() * sizeof(int)
            + port.sourceIndices.capacity() * sizeof(int)
            + port.patchSubdomains.capacity() * sizeof(int)
            + port.basis.capacity() * sizeof(double)
            + port.spectralValues.capacity() * sizeof(double)
            + port.spectralResiduals.capacity() * sizeof(double);
    }
    return model;
}

void saveLocalPortModel(const LocalPortModel& model,
                        const std::filesystem::path& path)
{
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("[Local port] Cannot create model file.");
    const char magic[8] = {'L','P','O','R','T','0','0','1'};
    writeRaw(out, magic, sizeof(magic));
    writeScalar(out, model.formatVersion);
    if (model.formatVersion >= 2) writeString(out, model.basisMethod);
    if (model.formatVersion >= 3) {
        writeString(out, model.innerProduct);
        writeString(out, model.rankMode);
        writeScalar(out, model.timeStep);
        writeScalar(out, model.oversamplingLayers);
        writeScalar(out, model.meshFingerprint);
        writeScalar(out, model.materialFingerprint);
        writeScalar(out, model.operatorFingerprint);
        writeScalar(out, model.inputFingerprint);
        writeScalar(out, model.penaltyFingerprint);
        writeString(out, model.buildTimestamp);
        writeString(out, model.sourceCommit);
    }
    if (model.formatVersion >= 4) {
        writeScalar(out, model.boundaryFingerprint);
        writeScalar(out, model.sourceMetadataFingerprint);
        writeScalar(out, model.rankFileFingerprint);
        writeScalar(out, model.requestedRank);
        writeScalar(out, model.minimumRank);
        writeScalar(out, model.maximumRank);
        writeScalar(out, model.eigenvalueTolerance);
        writeScalar(out, model.eigensolverTolerance);
        writeScalar(out, model.eigensolverMaximumIterations);
        writeScalar(out, model.relativeDeflationTolerance);
        writeString(out, model.innerSolver);
        writeScalar(out, model.innerSolverTolerance);
        writeScalar(out, model.innerSolverMaximumIterations);
    }
    if (model.formatVersion >= 5) {
        writeString(out, model.ablationMode);
    }
    if (model.formatVersion >= 6) {
        writeString(out, model.sourceMode);
        writeString(out, model.methodDescription);
        writeScalar(out, model.traceSourceFingerprint);
        writeScalar(out, model.generalizedInputFingerprint);
        writeScalar(out, model.generalizedBoundaryFingerprint);
        writeScalar(out, model.generalizedHistoryFingerprint);
    }
    if (model.formatVersion >= 7) {
        writeString(out, model.historyCompressionMethod);
        writeScalar(out, model.historyCompressionRank);
        writeScalar(out, model.historyCompressionTolerance);
    }
    writeScalar(out, model.fullInterfaceDofs);
    writeScalar(out, model.reducedInterfaceDofs);
    writeVector(out, model.interfaceGlobalDofs);
    writeScalar(out, static_cast<std::uint64_t>(model.ports.size()));
    for (const LocalPortBasis& port : model.ports) {
        writeScalar(out, port.interfaceId);
        writeScalar(out, port.leftSubdomain);
        writeScalar(out, port.rightSubdomain);
        writeScalar(out, port.rows);
        writeScalar(out, port.rank);
        writeScalar(out, port.candidateColumns);
        writeScalar(out, port.acceptedColumns);
        writeScalar(out, port.templateId);
        writeScalar(out, static_cast<std::uint8_t>(port.templateReused ? 1 : 0));
        writeScalar(out, port.fingerprint);
        writeScalar(out, port.retainedEnergy);
        writeScalar(out, port.orthogonalityError);
        if (model.formatVersion >= 2) {
            writeScalar(out, port.mandatoryModes);
            writeScalar(out, port.spectralModes);
            writeVector(out, port.spectralValues);
            writeVector(out, port.spectralResiduals);
        }
        if (model.formatVersion >= 3) {
            writeScalar(out, port.targetFingerprint);
            writeScalar(out, port.sourceFingerprint);
            writeScalar(out, port.transferIndicator);
            writeVector(out, port.patchSubdomains);
            writeVector(out, port.sourceIndices);
        }
        if (model.formatVersion >= 6) {
            writeScalar(out, port.traceSourceFingerprint);
            writeScalar(out, port.inputSourceFingerprint);
            writeScalar(out, port.boundarySourceFingerprint);
            writeScalar(out, port.historySourceFingerprint);
            writeScalar(out, port.requestedTransferRank);
            writeScalar(out, port.traceSourceRows);
            writeScalar(out, port.inputSourceRows);
            writeScalar(out, port.boundarySourceRows);
            writeScalar(out, port.historySourceRows);
        }
        writeVector(out, port.interfaceIndices);
        writeVector(out, port.basis);
    }
}

LocalPortModel loadLocalPortModel(const std::filesystem::path& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("[Local port] Cannot open model file.");
    char magic[8]{};
    readRaw(in, magic, sizeof(magic));
    if (std::memcmp(magic, "LPORT001", 8) != 0) {
        throw std::runtime_error("[Local port] Invalid model magic.");
    }
    LocalPortModel model;
    model.formatVersion = readScalar<int>(in);
    if (model.formatVersion != 1 && model.formatVersion != 2
        && model.formatVersion != 3 && model.formatVersion != 4
        && model.formatVersion != 5 && model.formatVersion != 6
        && model.formatVersion != 7) {
        throw std::runtime_error("[Local port] Unsupported model version.");
    }
    model.basisMethod = model.formatVersion >= 2 ? readString(in) : "port-pod";
    if (model.formatVersion >= 3) {
        model.innerProduct = readString(in);
        model.rankMode = readString(in);
        model.timeStep = readScalar<double>(in);
        model.oversamplingLayers = readScalar<int>(in);
        model.meshFingerprint = readScalar<std::uint64_t>(in);
        model.materialFingerprint = readScalar<std::uint64_t>(in);
        model.operatorFingerprint = readScalar<std::uint64_t>(in);
        model.inputFingerprint = readScalar<std::uint64_t>(in);
        model.penaltyFingerprint = readScalar<std::uint64_t>(in);
        model.buildTimestamp = readString(in);
        model.sourceCommit = readString(in);
    }
    if (model.formatVersion >= 4) {
        model.boundaryFingerprint = readScalar<std::uint64_t>(in);
        model.sourceMetadataFingerprint = readScalar<std::uint64_t>(in);
        model.rankFileFingerprint = readScalar<std::uint64_t>(in);
        model.requestedRank = readScalar<int>(in);
        model.minimumRank = readScalar<int>(in);
        model.maximumRank = readScalar<int>(in);
        model.eigenvalueTolerance = readScalar<double>(in);
        model.eigensolverTolerance = readScalar<double>(in);
        model.eigensolverMaximumIterations = readScalar<int>(in);
        model.relativeDeflationTolerance = readScalar<double>(in);
        model.innerSolver = readString(in);
        model.innerSolverTolerance = readScalar<double>(in);
        model.innerSolverMaximumIterations = readScalar<int>(in);
    }
    if (model.formatVersion >= 5) {
        model.ablationMode = readString(in);
    } else if (model.basisMethod == "steklov-schur") {
        model.ablationMode = "steklov-schur";
    }
    if (model.formatVersion >= 6) {
        model.sourceMode = readString(in);
        model.methodDescription = readString(in);
        model.traceSourceFingerprint = readScalar<std::uint64_t>(in);
        model.generalizedInputFingerprint =
            readScalar<std::uint64_t>(in);
        model.generalizedBoundaryFingerprint =
            readScalar<std::uint64_t>(in);
        model.generalizedHistoryFingerprint =
            readScalar<std::uint64_t>(in);
    }
    if (model.formatVersion >= 7) {
        model.historyCompressionMethod = readString(in);
        model.historyCompressionRank = readScalar<int>(in);
        model.historyCompressionTolerance = readScalar<double>(in);
    } else {
        model.historyCompressionMethod = "none";
        model.historyCompressionRank = 0;
        model.historyCompressionTolerance = 0.0;
    }
    model.fullInterfaceDofs = readScalar<int>(in);
    model.reducedInterfaceDofs = readScalar<int>(in);
    model.interfaceGlobalDofs = readVector<int>(in);
    const std::uint64_t portCount = readScalar<std::uint64_t>(in);
    model.ports.reserve(static_cast<std::size_t>(portCount));
    int rank = 0;
    for (std::uint64_t index = 0; index < portCount; ++index) {
        LocalPortBasis port;
        port.interfaceId = readScalar<int>(in);
        port.leftSubdomain = readScalar<int>(in);
        port.rightSubdomain = readScalar<int>(in);
        port.rows = readScalar<int>(in);
        port.rank = readScalar<int>(in);
        port.candidateColumns = readScalar<int>(in);
        port.acceptedColumns = readScalar<int>(in);
        port.templateId = readScalar<int>(in);
        port.templateReused = readScalar<std::uint8_t>(in) != 0;
        port.fingerprint = readScalar<std::uint64_t>(in);
        port.retainedEnergy = readScalar<double>(in);
        port.orthogonalityError = readScalar<double>(in);
        if (model.formatVersion >= 2) {
            port.mandatoryModes = readScalar<int>(in);
            port.spectralModes = readScalar<int>(in);
            port.spectralValues = readVector<double>(in);
            port.spectralResiduals = readVector<double>(in);
        }
        if (model.formatVersion >= 3) {
            port.targetFingerprint = readScalar<std::uint64_t>(in);
            port.sourceFingerprint = readScalar<std::uint64_t>(in);
            port.transferIndicator = readScalar<double>(in);
            port.patchSubdomains = readVector<int>(in);
            port.sourceIndices = readVector<int>(in);
        }
        if (model.formatVersion >= 6) {
            port.traceSourceFingerprint =
                readScalar<std::uint64_t>(in);
            port.inputSourceFingerprint =
                readScalar<std::uint64_t>(in);
            port.boundarySourceFingerprint =
                readScalar<std::uint64_t>(in);
            port.historySourceFingerprint =
                readScalar<std::uint64_t>(in);
            port.requestedTransferRank = readScalar<int>(in);
            port.traceSourceRows = readScalar<int>(in);
            port.inputSourceRows = readScalar<int>(in);
            port.boundarySourceRows = readScalar<int>(in);
            port.historySourceRows = readScalar<int>(in);
        }
        port.interfaceIndices = readVector<int>(in);
        port.basis = readVector<double>(in);
        if (port.rows != static_cast<int>(port.interfaceIndices.size())
            || port.basis.size() != static_cast<std::size_t>(port.rows * port.rank)
            || port.mandatoryModes < 0
            || port.spectralModes < 0
            || port.requestedTransferRank < 0
            || port.traceSourceRows < 0
            || port.inputSourceRows < 0
            || port.boundarySourceRows < 0
            || port.historySourceRows < 0
            || port.mandatoryModes + port.spectralModes > port.rank
            || port.spectralValues.size()
                != static_cast<std::size_t>(port.spectralModes)
            || port.spectralResiduals.size()
                != static_cast<std::size_t>(port.spectralModes)) {
            throw std::runtime_error("[Local port] Serialized basis dimensions are invalid.");
        }
        rank += port.rank;
        model.modelBytes += port.interfaceIndices.capacity() * sizeof(int)
            + port.sourceIndices.capacity() * sizeof(int)
            + port.patchSubdomains.capacity() * sizeof(int)
            + port.basis.capacity() * sizeof(double)
            + port.spectralValues.capacity() * sizeof(double)
            + port.spectralResiduals.capacity() * sizeof(double);
        model.ports.push_back(std::move(port));
    }
    if (rank != model.reducedInterfaceDofs
        || model.fullInterfaceDofs != static_cast<int>(model.interfaceGlobalDofs.size())) {
        throw std::runtime_error("[Local port] Serialized model dimensions are inconsistent.");
    }
    model.modelBytes += model.interfaceGlobalDofs.capacity() * sizeof(int);
    return model;
}

LocalPortReducedSchurSolver::LocalPortReducedSchurSolver(
    const local::Model& dynamicModel,
    const LocalPortModel& portModel)
    : dynamicModel_(dynamicModel), portModel_(portModel)
{
    if (dynamicModel.interfaceDofs != portModel.fullInterfaceDofs
        || dynamicModel.interfaceGlobalDofs != portModel.interfaceGlobalDofs
        || portModel.reducedInterfaceDofs <= 0) {
        throw std::runtime_error("[Local port] Dynamic/port model ordering mismatch.");
    }
    const auto assemblyStart = Clock::now();
    localFactors_.reserve(dynamicModel.subdomains.size());
    for (const local::SubdomainModel& subdomain : dynamicModel.subdomains) {
        localFactors_.push_back(local::factorDenseSymmetric(
            subdomain.reducedInterior, subdomain.rank));
    }
    std::vector<local::InterfaceEntry> entries = dynamicModel.interfaceEntries;
    std::sort(entries.begin(), entries.end(), [](const auto& left, const auto& right) {
        return left.row < right.row
            || (left.row == right.row && left.column < right.column);
    });
    rowPtr_.assign(static_cast<std::size_t>(dynamicModel.interfaceDofs + 1), 0);
    for (std::size_t begin = 0; begin < entries.size();) {
        std::size_t end = begin + 1;
        double value = entries[begin].value;
        while (end < entries.size() && entries[end].row == entries[begin].row
               && entries[end].column == entries[begin].column) {
            value += entries[end].value;
            ++end;
        }
        if (value != 0.0) {
            columnIndices_.push_back(entries[begin].column);
            interfaceValues_.push_back(value);
            ++rowPtr_[static_cast<std::size_t>(entries[begin].row + 1)];
        }
        begin = end;
    }
    std::partial_sum(rowPtr_.begin(), rowPtr_.end(), rowPtr_.begin());
    const int reduced = portModel.reducedInterfaceDofs;
    std::vector<double> matrix(static_cast<std::size_t>(reduced * reduced), 0.0);
    std::vector<double> unit(static_cast<std::size_t>(reduced), 0.0);
    for (int column = 0; column < reduced; ++column) {
        unit[static_cast<std::size_t>(column)] = 1.0;
        std::vector<double> full;
        lift(unit, full);
        std::vector<double> image;
        applyFullInterface(full, image);
        std::vector<double> restricted;
        restrict(image, restricted);
        for (int row = 0; row < reduced; ++row) {
            matrix[static_cast<std::size_t>(row * reduced + column)] =
                restricted[static_cast<std::size_t>(row)];
        }
        unit[static_cast<std::size_t>(column)] = 0.0;
    }
    long double asymmetrySquared = 0.0L;
    long double matrixSquared = 0.0L;
    for (int row = 0; row < reduced; ++row) {
        for (int column = 0; column < reduced; ++column) {
            const double value = matrix[static_cast<std::size_t>(
                row * reduced + column)];
            const double transpose = matrix[static_cast<std::size_t>(
                column * reduced + row)];
            const double difference = value - transpose;
            asymmetrySquared += static_cast<long double>(difference) * difference;
            matrixSquared += static_cast<long double>(value) * value;
        }
    }
    relativeAsymmetry_ = std::sqrt(static_cast<double>(asymmetrySquared))
        / std::max(1.0e-300, std::sqrt(static_cast<double>(matrixSquared)));
    assemblySeconds_ = secondsSince(assemblyStart);
    const auto factorStart = Clock::now();
#ifdef USE_MKL_PARDISO
    static_assert(sizeof(lapack_int) == sizeof(int),
        "LP64 LAPACK integer size is required for serialized port pivots.");
    portLu_ = std::move(matrix);
    portPivots_.assign(static_cast<std::size_t>(reduced), 0);
    const lapack_int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR,
        reduced, reduced, portLu_.data(), reduced,
        reinterpret_cast<lapack_int*>(portPivots_.data()));
    if (info != 0) {
        throw std::runtime_error(
            "[Local port] Projected Dynamic Schur LU factorization failed with info="
            + std::to_string(info));
    }
#else
    for (int row = 0; row < reduced; ++row) {
        for (int column = 0; column < row; ++column) {
            const double average = 0.5 * (
                matrix[static_cast<std::size_t>(row * reduced + column)]
                + matrix[static_cast<std::size_t>(column * reduced + row)]);
            matrix[static_cast<std::size_t>(row * reduced + column)] = average;
            matrix[static_cast<std::size_t>(column * reduced + row)] = average;
        }
    }
    portFactor_ = local::factorDenseSymmetric(matrix, reduced);
#endif
    factorizationSeconds_ = secondsSince(factorStart);
}

void LocalPortReducedSchurSolver::solveLocalCoordinates(
    std::vector<double>& values) const
{
    if (values.size()
        != static_cast<std::size_t>(portModel_.reducedInterfaceDofs)) {
        throw std::runtime_error(
            "[Local port] Local coordinate solve size mismatch.");
    }
#ifdef USE_MKL_PARDISO
    const lapack_int solveInfo = LAPACKE_dgetrs(
        LAPACK_ROW_MAJOR, 'N',
        portModel_.reducedInterfaceDofs, 1,
        portLu_.data(), portModel_.reducedInterfaceDofs,
        reinterpret_cast<const lapack_int*>(portPivots_.data()),
        values.data(), 1);
    if (solveInfo != 0) {
        throw std::runtime_error(
            "[Local port] Projected Dynamic Schur LU solve failed with info="
            + std::to_string(solveInfo));
    }
#else
    local::solveDenseSymmetric(portFactor_, values);
#endif
}

int LocalPortReducedSchurSolver::coarseDimension() const
{
    return coarseModel_ ? coarseModel_->rank : 0;
}

void LocalPortReducedSchurSolver::attachGlobalCoarse(
    const GlobalInterfaceCoarseModel& coarseModel,
    bool includeLocalBasis)
{
    if (coarseModel_ != nullptr) {
        throw std::runtime_error(
            "[Global coarse] A coarse model is already attached.");
    }
    if (coarseModel.fullInterfaceDofs != portModel_.fullInterfaceDofs
        || coarseModel.interfaceGlobalDofs
            != portModel_.interfaceGlobalDofs
        || coarseModel.rank <= 0
        || coarseModel.basis.size() != static_cast<std::size_t>(
            coarseModel.fullInterfaceDofs * coarseModel.rank)) {
        throw std::runtime_error(
            "[Global coarse] Coarse/local interface ordering mismatch.");
    }
    const auto assemblyStart = Clock::now();
    includeLocalBasis_ = includeLocalBasis;
    const int localRank = includeLocalBasis_
        ? portModel_.reducedInterfaceDofs : 0;
    const int coarseRank = coarseModel.rank;
    localCoarse_.assign(
        static_cast<std::size_t>(localRank * coarseRank), 0.0);
    std::vector<double> coarseMatrix(
        static_cast<std::size_t>(coarseRank * coarseRank), 0.0);
    std::vector<std::vector<double>> images(
        static_cast<std::size_t>(coarseRank));
    for (int column = 0; column < coarseRank; ++column) {
        if (coarseModel.schurImages.size()
            == coarseModel.basis.size()) {
            images[static_cast<std::size_t>(column)].assign(
                coarseModel.schurImages.begin()
                    + static_cast<std::ptrdiff_t>(
                        column * coarseModel.fullInterfaceDofs),
                coarseModel.schurImages.begin()
                    + static_cast<std::ptrdiff_t>(
                        (column + 1)
                        * coarseModel.fullInterfaceDofs));
        } else {
            std::vector<double> basisColumn(
                coarseModel.basis.begin()
                    + static_cast<std::ptrdiff_t>(
                        column * coarseModel.fullInterfaceDofs),
                coarseModel.basis.begin()
                    + static_cast<std::ptrdiff_t>(
                        (column + 1)
                        * coarseModel.fullInterfaceDofs));
            applyFullInterface(
                basisColumn, images[static_cast<std::size_t>(column)]);
        }
        if (includeLocalBasis_) {
            std::vector<double> restricted;
            restrict(images[static_cast<std::size_t>(column)], restricted);
            for (int row = 0; row < localRank; ++row) {
                localCoarse_[static_cast<std::size_t>(
                    row * coarseRank + column)] =
                    restricted[static_cast<std::size_t>(row)];
            }
        }
    }
    for (int row = 0; row < coarseRank; ++row) {
        const double* basis = coarseModel.basis.data()
            + static_cast<std::size_t>(
                row * coarseModel.fullInterfaceDofs);
        for (int column = 0; column < coarseRank; ++column) {
            coarseMatrix[static_cast<std::size_t>(
                row * coarseRank + column)] =
                dot(basis,
                    images[static_cast<std::size_t>(column)].data(),
                    coarseModel.fullInterfaceDofs);
        }
    }
    long double asymmetrySquared = 0.0L;
    long double matrixSquared = 0.0L;
    for (int row = 0; row < coarseRank; ++row) {
        for (int column = 0; column < coarseRank; ++column) {
            const double value = coarseMatrix[static_cast<std::size_t>(
                row * coarseRank + column)];
            const double transpose =
                coarseMatrix[static_cast<std::size_t>(
                    column * coarseRank + row)];
            const double difference = value - transpose;
            asymmetrySquared +=
                static_cast<long double>(difference) * difference;
            matrixSquared +=
                static_cast<long double>(value) * value;
        }
    }
    augmentedRelativeAsymmetry_ =
        std::sqrt(static_cast<double>(asymmetrySquared))
        / std::max(
            1.0e-300,
            std::sqrt(static_cast<double>(matrixSquared)));

    std::vector<double> solvedCross(
        static_cast<std::size_t>(localRank * coarseRank), 0.0);
    if (includeLocalBasis_) {
        for (int column = 0; column < coarseRank; ++column) {
            std::vector<double> values(
                static_cast<std::size_t>(localRank), 0.0);
            for (int row = 0; row < localRank; ++row) {
                values[static_cast<std::size_t>(row)] =
                    localCoarse_[static_cast<std::size_t>(
                        row * coarseRank + column)];
            }
            solveLocalCoordinates(values);
            for (int row = 0; row < localRank; ++row) {
                solvedCross[static_cast<std::size_t>(
                    row * coarseRank + column)] =
                    values[static_cast<std::size_t>(row)];
            }
        }
    }
    for (int row = 0; row < coarseRank; ++row) {
        for (int column = 0; column < coarseRank; ++column) {
            long double correction = 0.0L;
            for (int localMode = 0;
                 localMode < localRank; ++localMode) {
                correction += static_cast<long double>(
                    localCoarse_[static_cast<std::size_t>(
                        localMode * coarseRank + row)])
                    * solvedCross[static_cast<std::size_t>(
                        localMode * coarseRank + column)];
            }
            coarseMatrix[static_cast<std::size_t>(
                row * coarseRank + column)] -=
                static_cast<double>(correction);
        }
    }
    assemblySeconds_ += secondsSince(assemblyStart);
    const auto factorStart = Clock::now();
    coarseBorderFactor_ =
        local::factorDenseSymmetric(coarseMatrix, coarseRank);
    factorizationSeconds_ += secondsSince(factorStart);
    coarseModel_ = &coarseModel;
}

void LocalPortReducedSchurSolver::lift(const std::vector<double>& reduced,
                                       std::vector<double>& full) const
{
    if (reduced.size() != static_cast<std::size_t>(portModel_.reducedInterfaceDofs)) {
        throw std::runtime_error("[Local port] Reduced interface vector size mismatch.");
    }
    full.assign(static_cast<std::size_t>(portModel_.fullInterfaceDofs), 0.0);
    int offset = 0;
    for (const LocalPortBasis& port : portModel_.ports) {
        for (int mode = 0; mode < port.rank; ++mode) {
            const double coefficient = reduced[static_cast<std::size_t>(offset + mode)];
            for (int row = 0; row < port.rows; ++row) {
                full[static_cast<std::size_t>(port.interfaceIndices[static_cast<std::size_t>(row)])]
                    += coefficient * port.basis[static_cast<std::size_t>(mode * port.rows + row)];
            }
        }
        offset += port.rank;
    }
}

void LocalPortReducedSchurSolver::restrict(const std::vector<double>& full,
                                           std::vector<double>& reduced) const
{
    if (full.size() != static_cast<std::size_t>(portModel_.fullInterfaceDofs)) {
        throw std::runtime_error("[Local port] Full interface vector size mismatch.");
    }
    reduced.assign(static_cast<std::size_t>(portModel_.reducedInterfaceDofs), 0.0);
    int offset = 0;
    for (const LocalPortBasis& port : portModel_.ports) {
        for (int mode = 0; mode < port.rank; ++mode) {
            for (int row = 0; row < port.rows; ++row) {
                reduced[static_cast<std::size_t>(offset + mode)] +=
                    port.basis[static_cast<std::size_t>(mode * port.rows + row)]
                    * full[static_cast<std::size_t>(port.interfaceIndices[static_cast<std::size_t>(row)])];
            }
        }
        offset += port.rank;
    }
}

void LocalPortReducedSchurSolver::restrictLocal(
    const std::vector<double>& full,
    std::vector<double>& reduced) const
{
    restrict(full, reduced);
}

void LocalPortReducedSchurSolver::localGalerkinInterfaceResponse(
    const std::vector<double>& interfaceRhs,
    std::vector<double>& response) const
{
    std::vector<double> coordinates;
    restrict(interfaceRhs, coordinates);
    solveLocalCoordinates(coordinates);
    lift(coordinates, response);
}

void LocalPortReducedSchurSolver::projectSchurEnergyComplement(
    const std::vector<double>& input,
    std::vector<double>& output) const
{
    std::vector<double> image;
    applyFullInterface(input, image);
    std::vector<double> projection;
    localGalerkinInterfaceResponse(image, projection);
    output.resize(input.size());
    for (std::size_t row = 0; row < input.size(); ++row) {
        output[row] = input[row] - projection[row];
    }
}

void LocalPortReducedSchurSolver::liftCoarse(
    const std::vector<double>& reduced,
    std::vector<double>& full) const
{
    if (!coarseModel_
        || reduced.size()
            != static_cast<std::size_t>(coarseModel_->rank)) {
        throw std::runtime_error(
            "[Global coarse] Coarse lift dimensions are invalid.");
    }
    full.assign(
        static_cast<std::size_t>(coarseModel_->fullInterfaceDofs),
        0.0);
    for (int mode = 0; mode < coarseModel_->rank; ++mode) {
        const double coefficient =
            reduced[static_cast<std::size_t>(mode)];
        const double* basis = coarseModel_->basis.data()
            + static_cast<std::size_t>(
                mode * coarseModel_->fullInterfaceDofs);
        for (int row = 0;
             row < coarseModel_->fullInterfaceDofs; ++row) {
            full[static_cast<std::size_t>(row)] +=
                coefficient * basis[row];
        }
    }
}

void LocalPortReducedSchurSolver::restrictCoarse(
    const std::vector<double>& full,
    std::vector<double>& reduced) const
{
    if (!coarseModel_
        || full.size() != static_cast<std::size_t>(
            coarseModel_->fullInterfaceDofs)) {
        throw std::runtime_error(
            "[Global coarse] Coarse restriction dimensions are invalid.");
    }
    reduced.assign(
        static_cast<std::size_t>(coarseModel_->rank), 0.0);
    for (int mode = 0; mode < coarseModel_->rank; ++mode) {
        reduced[static_cast<std::size_t>(mode)] =
            dot(coarseModel_->basis.data()
                    + static_cast<std::size_t>(
                        mode * coarseModel_->fullInterfaceDofs),
                full.data(), coarseModel_->fullInterfaceDofs);
    }
}

double LocalPortReducedSchurSolver::relativeProjectionError(
    const std::vector<double>& full) const
{
    if (full.size() != static_cast<std::size_t>(portModel_.fullInterfaceDofs)) {
        throw std::runtime_error("[Local port] Projection vector size mismatch.");
    }
    long double errorSquared = 0.0L;
    long double referenceSquared = 0.0L;
    for (const LocalPortBasis& port : portModel_.ports) {
        std::vector<double> coefficients(static_cast<std::size_t>(port.rank), 0.0);
        for (int mode = 0; mode < port.rank; ++mode) {
            for (int row = 0; row < port.rows; ++row) {
                coefficients[static_cast<std::size_t>(mode)] +=
                    port.basis[static_cast<std::size_t>(mode * port.rows + row)]
                    * full[static_cast<std::size_t>(
                        port.interfaceIndices[static_cast<std::size_t>(row)])];
            }
        }
        for (int row = 0; row < port.rows; ++row) {
            double projected = 0.0;
            for (int mode = 0; mode < port.rank; ++mode) {
                projected += port.basis[static_cast<std::size_t>(
                    mode * port.rows + row)]
                    * coefficients[static_cast<std::size_t>(mode)];
            }
            const double reference = full[static_cast<std::size_t>(
                port.interfaceIndices[static_cast<std::size_t>(row)])];
            const double error = projected - reference;
            errorSquared += static_cast<long double>(error) * error;
            referenceSquared += static_cast<long double>(reference) * reference;
        }
    }
    return std::sqrt(static_cast<double>(errorSquared))
        / std::max(1.0e-300, std::sqrt(static_cast<double>(referenceSquared)));
}

void LocalPortReducedSchurSolver::applyFullInterface(
    const std::vector<double>& input,
    std::vector<double>& output) const
{
    if (input.size() != static_cast<std::size_t>(dynamicModel_.interfaceDofs)) {
        throw std::runtime_error("[Local port] Full Schur input size mismatch.");
    }
    output.assign(input.size(), 0.0);
    for (int row = 0; row < dynamicModel_.interfaceDofs; ++row) {
        for (int entry = rowPtr_[static_cast<std::size_t>(row)];
             entry < rowPtr_[static_cast<std::size_t>(row + 1)]; ++entry) {
            output[static_cast<std::size_t>(row)] +=
                interfaceValues_[static_cast<std::size_t>(entry)]
                * input[static_cast<std::size_t>(columnIndices_[static_cast<std::size_t>(entry)])];
        }
    }
    for (std::size_t slot = 0; slot < dynamicModel_.subdomains.size(); ++slot) {
        const local::SubdomainModel& local = dynamicModel_.subdomains[slot];
        std::vector<double> eliminated(static_cast<std::size_t>(local.rank), 0.0);
        for (int mode = 0; mode < local.rank; ++mode) {
            for (std::size_t gamma = 0; gamma < local.interfaceIndices.size(); ++gamma) {
                eliminated[static_cast<std::size_t>(mode)] +=
                    local.reducedInteriorInterface[static_cast<std::size_t>(
                        mode * local.localInterfaceDofs) + gamma]
                    * input[static_cast<std::size_t>(local.interfaceIndices[gamma])];
            }
        }
        local::solveDenseSymmetric(localFactors_[slot], eliminated);
        for (std::size_t gamma = 0; gamma < local.interfaceIndices.size(); ++gamma) {
            double correction = 0.0;
            for (int mode = 0; mode < local.rank; ++mode) {
                correction += local.reducedInterfaceInterior[static_cast<std::size_t>(
                    gamma * static_cast<std::size_t>(local.rank) + mode)]
                    * eliminated[static_cast<std::size_t>(mode)];
            }
            output[static_cast<std::size_t>(local.interfaceIndices[gamma])] -= correction;
        }
    }
}

LocalPortSolveResult LocalPortReducedSchurSolver::solve(
    const std::vector<double>& globalRhs) const
{
    if (globalRhs.size() != static_cast<std::size_t>(dynamicModel_.globalDofs)) {
        throw std::runtime_error("[Local port] Global RHS size mismatch.");
    }
    const auto totalStart = Clock::now();
    LocalPortSolveResult result;
    const auto condensedStart = Clock::now();
    std::vector<double> condensed(static_cast<std::size_t>(dynamicModel_.interfaceDofs), 0.0);
    for (int gamma = 0; gamma < dynamicModel_.interfaceDofs; ++gamma) {
        condensed[static_cast<std::size_t>(gamma)] = globalRhs[static_cast<std::size_t>(
            dynamicModel_.interfaceGlobalDofs[static_cast<std::size_t>(gamma)])];
    }
    std::vector<std::vector<double>> projected(dynamicModel_.subdomains.size());
    for (std::size_t slot = 0; slot < dynamicModel_.subdomains.size(); ++slot) {
        const local::SubdomainModel& local = dynamicModel_.subdomains[slot];
        projected[slot].assign(static_cast<std::size_t>(local.rank), 0.0);
        for (int mode = 0; mode < local.rank; ++mode) {
            for (int row = 0; row < local.interiorDofs; ++row) {
                projected[slot][static_cast<std::size_t>(mode)] +=
                    local.basis[static_cast<std::size_t>(mode * local.interiorDofs + row)]
                    * globalRhs[static_cast<std::size_t>(
                        local.interiorGlobalDofs[static_cast<std::size_t>(row)])];
            }
        }
        std::vector<double> eliminated = projected[slot];
        local::solveDenseSymmetric(localFactors_[slot], eliminated);
        for (std::size_t localGamma = 0;
             localGamma < local.interfaceIndices.size(); ++localGamma) {
            const int gamma = local.interfaceIndices[localGamma];
            for (int mode = 0; mode < local.rank; ++mode) {
                condensed[static_cast<std::size_t>(gamma)] -=
                    local.reducedInterfaceInterior[static_cast<std::size_t>(
                        localGamma * static_cast<std::size_t>(local.rank) + mode)]
                    * eliminated[static_cast<std::size_t>(mode)];
            }
        }
    }
    result.condensedRhsSeconds = secondsSince(condensedStart);
    std::vector<double> portRhs;
    if (!coarseModel_ || includeLocalBasis_) {
        restrict(condensed, portRhs);
    }
    const auto solveStart = Clock::now();
    std::vector<double> coordinates = portRhs;
    std::vector<double> coarseCoordinates;
    if (!coarseModel_) {
        solveLocalCoordinates(coordinates);
    } else if (includeLocalBasis_) {
        solveLocalCoordinates(coordinates);
        std::vector<double> coarseRhs;
        restrictCoarse(condensed, coarseRhs);
        for (int coarse = 0;
             coarse < coarseModel_->rank; ++coarse) {
            for (int localMode = 0;
                 localMode < portModel_.reducedInterfaceDofs;
                 ++localMode) {
                coarseRhs[static_cast<std::size_t>(coarse)] -=
                    localCoarse_[static_cast<std::size_t>(
                        localMode * coarseModel_->rank + coarse)]
                    * coordinates[static_cast<std::size_t>(localMode)];
            }
        }
        coarseCoordinates = coarseRhs;
        local::solveDenseSymmetric(
            coarseBorderFactor_, coarseCoordinates);
        coordinates = portRhs;
        for (int localMode = 0;
             localMode < portModel_.reducedInterfaceDofs;
             ++localMode) {
            for (int coarse = 0;
                 coarse < coarseModel_->rank; ++coarse) {
                coordinates[static_cast<std::size_t>(localMode)] -=
                    localCoarse_[static_cast<std::size_t>(
                        localMode * coarseModel_->rank + coarse)]
                    * coarseCoordinates[static_cast<std::size_t>(coarse)];
            }
        }
        solveLocalCoordinates(coordinates);
    } else {
        restrictCoarse(condensed, coarseCoordinates);
        local::solveDenseSymmetric(
            coarseBorderFactor_, coarseCoordinates);
    }
    result.reducedSolveSeconds = secondsSince(solveStart);
    if (!coarseModel_ || includeLocalBasis_) {
        lift(coordinates, result.solution.interfaceTemperature);
    } else {
        result.solution.interfaceTemperature.assign(
            condensed.size(), 0.0);
    }
    if (coarseModel_) {
        std::vector<double> coarseTemperature;
        liftCoarse(coarseCoordinates, coarseTemperature);
        for (std::size_t row = 0;
             row < coarseTemperature.size(); ++row) {
            result.solution.interfaceTemperature[row] +=
                coarseTemperature[row];
        }
    }
    std::vector<double> image;
    applyFullInterface(result.solution.interfaceTemperature, image);
    long double residualSquared = 0.0L;
    long double rhsSquared = 0.0L;
    for (std::size_t row = 0; row < condensed.size(); ++row) {
        const double residual = image[row] - condensed[row];
        residualSquared += static_cast<long double>(residual) * residual;
        rhsSquared += static_cast<long double>(condensed[row]) * condensed[row];
    }
    result.portRelativeResidual = std::sqrt(static_cast<double>(residualSquared))
        / std::max(1.0e-300, std::sqrt(static_cast<double>(rhsSquared)));
    std::vector<double> fullResidual(condensed.size(), 0.0);
    for (std::size_t row = 0; row < condensed.size(); ++row) {
        fullResidual[row] = image[row] - condensed[row];
    }
    std::vector<double> reducedResidual;
    if (!coarseModel_ || includeLocalBasis_) {
        restrict(fullResidual, reducedResidual);
    }
    std::vector<double> coarseResidual;
    if (coarseModel_) {
        restrictCoarse(fullResidual, coarseResidual);
    }
    long double reducedResidualSquared = 0.0L;
    long double reducedRhsSquared = 0.0L;
    for (std::size_t row = 0; row < reducedResidual.size(); ++row) {
        reducedResidualSquared += static_cast<long double>(reducedResidual[row])
            * reducedResidual[row];
        reducedRhsSquared += static_cast<long double>(portRhs[row]) * portRhs[row];
    }
    if (coarseModel_) {
        std::vector<double> coarseRhs;
        restrictCoarse(condensed, coarseRhs);
        for (std::size_t row = 0;
             row < coarseResidual.size(); ++row) {
            reducedResidualSquared +=
                static_cast<long double>(coarseResidual[row])
                * coarseResidual[row];
            reducedRhsSquared +=
                static_cast<long double>(coarseRhs[row])
                * coarseRhs[row];
        }
    }
    result.reducedRelativeResidual =
        std::sqrt(static_cast<double>(reducedResidualSquared))
        / std::max(1.0e-300, std::sqrt(static_cast<double>(reducedRhsSquared)));
    result.solution.timing.interfaceRelativeResidual = result.portRelativeResidual;
    result.solution.timing.interfaceSolveSeconds = result.reducedSolveSeconds;

    const auto recoveryStart = Clock::now();
    result.solution.temperature.assign(static_cast<std::size_t>(dynamicModel_.globalDofs), 0.0);
    for (int gamma = 0; gamma < dynamicModel_.interfaceDofs; ++gamma) {
        result.solution.temperature[static_cast<std::size_t>(
            dynamicModel_.interfaceGlobalDofs[static_cast<std::size_t>(gamma)])] =
            result.solution.interfaceTemperature[static_cast<std::size_t>(gamma)];
    }
    for (std::size_t slot = 0; slot < dynamicModel_.subdomains.size(); ++slot) {
        const local::SubdomainModel& local = dynamicModel_.subdomains[slot];
        std::vector<double> state = projected[slot];
        for (int mode = 0; mode < local.rank; ++mode) {
            for (std::size_t localGamma = 0;
                 localGamma < local.interfaceIndices.size(); ++localGamma) {
                state[static_cast<std::size_t>(mode)] -=
                    local.reducedInteriorInterface[static_cast<std::size_t>(
                        mode * local.localInterfaceDofs) + localGamma]
                    * result.solution.interfaceTemperature[static_cast<std::size_t>(
                        local.interfaceIndices[localGamma])];
            }
        }
        local::solveDenseSymmetric(localFactors_[slot], state);
        for (int row = 0; row < local.interiorDofs; ++row) {
            for (int mode = 0; mode < local.rank; ++mode) {
                result.solution.temperature[static_cast<std::size_t>(
                    local.interiorGlobalDofs[static_cast<std::size_t>(row)])] +=
                    local.basis[static_cast<std::size_t>(mode * local.interiorDofs + row)]
                    * state[static_cast<std::size_t>(mode)];
            }
        }
    }
    result.solution.timing.localRecoverySeconds = secondsSince(recoveryStart);
    result.solution.timing.totalSeconds = secondsSince(totalStart);
    result.solution.status = "success";
    return result;
}

std::size_t LocalPortReducedSchurSolver::factorMemoryBytes() const
{
    std::size_t bytes = rowPtr_.capacity() * sizeof(int)
        + columnIndices_.capacity() * sizeof(int)
        + interfaceValues_.capacity() * sizeof(double)
        + portFactor_.lower.capacity() * sizeof(double)
        + portFactor_.diagonal.capacity() * sizeof(double)
        + portLu_.capacity() * sizeof(double)
        + portPivots_.capacity() * sizeof(int)
        + localCoarse_.capacity() * sizeof(double)
        + coarseBorderFactor_.lower.capacity() * sizeof(double)
        + coarseBorderFactor_.diagonal.capacity() * sizeof(double);
    for (const auto& factor : localFactors_) {
        bytes += factor.lower.capacity() * sizeof(double)
            + factor.diagonal.capacity() * sizeof(double);
    }
    return bytes;
}

} // namespace mor::transient
