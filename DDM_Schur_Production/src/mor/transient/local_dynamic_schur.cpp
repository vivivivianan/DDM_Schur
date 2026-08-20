// Validated transient Dynamic Schur workflow.
//
// Production-reachable stages are: descriptor load/build, global-FOM trace
// construction, M1 local Block-Arnoldi, identical-subdomain template reuse,
// local-model cache I/O, augmented-direct setup, native reduced-history time
// stepping, residual/timing summaries, and optional summary-only FOM audit.
//
// This source is intentionally retained as one upstream unit to preserve cache
// fingerprints, diagnostics columns, and floating-point operation ordering.
// Historical port/coarse experiment branches remain syntactically present but
// are fixed off by main.cpp; their implementations are absent and protected by
// internal/removed_research_methods.cpp. Reaching one is a hard error.

#include "local_dynamic_schur.hpp"

#include "block_arnoldi.hpp"
#include "ddm_schur/interface_operator.hpp"
#include "ddm_schur/schur_fgmres.hpp"
#include "linear_solvers.hpp"
#include "global_interface_coarse.hpp"
#include "global_randomized_schur.hpp"
#include "local_port_reduced_schur.hpp"
#include "projection_diagnosis.hpp"
#include "mor/local/local_reduced_schur.hpp"
#include "mor/mor_diagnostics.hpp"
#include "optimal_port_space.hpp"
#include "port/randomized_transfer_port.hpp"
#include "thermal_descriptor_system.hpp"
#include "transient_input_waveform.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

namespace mor::transient {
namespace {

using Clock = std::chrono::steady_clock;

double elapsed(const Clock::time_point& start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

std::vector<std::string> parseCsvLine(const std::string& line)
{
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (std::size_t index = 0; index < line.size(); ++index) {
        const char value = line[index];
        if (value == '"') {
            if (quoted && index + 1 < line.size()
                && line[index + 1] == '"') {
                field.push_back('"');
                ++index;
            } else {
                quoted = !quoted;
            }
        } else if (value == ',' && !quoted) {
            fields.push_back(std::move(field));
            field.clear();
        } else {
            field.push_back(value);
        }
    }
    if (quoted) {
        throw std::runtime_error(
            "[Optimal port] Unterminated quoted topology-audit CSV field.");
    }
    fields.push_back(std::move(field));
    return fields;
}

struct ScalabilityAuditRow {
    int interfaceId = -1;
    int leftSubdomain = -1;
    int rightSubdomain = -1;
    int targetDofs = 0;
    int sourceDofs = 0;
};

std::vector<ScalabilityAuditRow> readScalabilityAudit(
    const std::filesystem::path& path,
    const std::string& caseName)
{
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error(
            "[Optimal port] Cannot read topology audit CSV: "
            + path.string());
    }
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error(
            "[Optimal port] Topology audit CSV is empty.");
    }
    std::vector<std::string> header = parseCsvLine(line);
    if (!header.empty() && header.front().size() >= 3
        && static_cast<unsigned char>(header.front()[0]) == 0xef
        && static_cast<unsigned char>(header.front()[1]) == 0xbb
        && static_cast<unsigned char>(header.front()[2]) == 0xbf) {
        header.front().erase(0, 3);
    }
    std::map<std::string, std::size_t> columns;
    for (std::size_t column = 0; column < header.size(); ++column) {
        columns.emplace(header[column], column);
    }
    const std::vector<std::string> required{
        "case", "interface_id", "left_subdomain", "right_subdomain",
        "target_dofs", "source_dofs"};
    for (const std::string& name : required) {
        if (columns.count(name) == 0) {
            throw std::runtime_error(
                "[Optimal port] Topology audit CSV lacks column " + name + '.');
        }
    }
    const auto cell = [&](const std::vector<std::string>& fields,
                          const std::string& name) -> const std::string& {
        const std::size_t column = columns.at(name);
        if (column >= fields.size()) {
            throw std::runtime_error(
                "[Optimal port] Short topology-audit CSV row.");
        }
        return fields[column];
    };
    std::vector<ScalabilityAuditRow> rows;
    std::set<int> identifiers;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> fields = parseCsvLine(line);
        if (cell(fields, "case") != caseName) continue;
        ScalabilityAuditRow row;
        row.interfaceId = std::stoi(cell(fields, "interface_id"));
        row.leftSubdomain = std::stoi(cell(fields, "left_subdomain"));
        row.rightSubdomain = std::stoi(cell(fields, "right_subdomain"));
        row.targetDofs = std::stoi(cell(fields, "target_dofs"));
        row.sourceDofs = std::stoi(cell(fields, "source_dofs"));
        if (!identifiers.insert(row.interfaceId).second) {
            throw std::runtime_error(
                "[Optimal port] Duplicate interface in topology audit CSV.");
        }
        rows.push_back(row);
    }
    if (rows.empty()) {
        throw std::runtime_error(
            "[Optimal port] Topology audit CSV has no rows for case "
            + caseName + '.');
    }
    return rows;
}

void writeCsvString(std::ostream& output, const std::string& value)
{
    output << '"';
    for (char character : value) {
        if (character == '"') output << '"';
        output << character;
    }
    output << '"';
}

double normSquared(const std::vector<double>& values)
{
    long double result = 0.0L;
    for (double value : values) result += static_cast<long double>(value) * value;
    return static_cast<double>(result);
}

double relativeDifference(const std::vector<double>& left,
                          const std::vector<double>& right)
{
    if (left.size() != right.size()) {
        return std::numeric_limits<double>::infinity();
    }
    double differenceScale = 0.0;
    double differenceSum = 1.0;
    double referenceScale = 0.0;
    double referenceSum = 1.0;
    const auto accumulateScaled = [](double value,
                                     double& scale,
                                     double& sum) {
        const double magnitude = std::abs(value);
        if (magnitude == 0.0) return;
        if (!std::isfinite(magnitude)) {
            scale = std::numeric_limits<double>::infinity();
            sum = 1.0;
        } else if (scale < magnitude) {
            const double ratio = scale / magnitude;
            sum = 1.0 + sum * ratio * ratio;
            scale = magnitude;
        } else {
            const double ratio = magnitude / scale;
            sum += ratio * ratio;
        }
    };
    for (std::size_t row = 0; row < left.size(); ++row) {
        const double delta = left[row] - right[row];
        accumulateScaled(delta, differenceScale, differenceSum);
        accumulateScaled(right[row], referenceScale, referenceSum);
    }
    if (differenceScale == 0.0) return 0.0;
    if (!std::isfinite(differenceScale)
        || !std::isfinite(referenceScale)) {
        return std::numeric_limits<double>::infinity();
    }
    if (referenceScale == 0.0) {
        return differenceScale * std::sqrt(differenceSum)
            / 1.0e-300;
    }
    return (differenceScale / referenceScale)
        * std::sqrt(differenceSum / referenceSum);
}

std::vector<double> euclideanOrthonormalBasis(
    const LocalPortBasis& port)
{
    std::vector<double> basis;
    basis.reserve(port.basis.size());
    int accepted = 0;
    for (int column = 0; column < port.rank; ++column) {
        std::vector<double> candidate(
            port.basis.begin()
                + static_cast<std::ptrdiff_t>(column * port.rows),
            port.basis.begin()
                + static_cast<std::ptrdiff_t>((column + 1) * port.rows));
        const double original =
            std::sqrt(std::max(0.0, normSquared(candidate)));
        for (int pass = 0; pass < 2; ++pass) {
            for (int prior = 0; prior < accepted; ++prior) {
                const double* vector = basis.data()
                    + static_cast<std::size_t>(prior * port.rows);
                long double coefficient = 0.0L;
                for (int row = 0; row < port.rows; ++row) {
                    coefficient += static_cast<long double>(vector[row])
                        * candidate[static_cast<std::size_t>(row)];
                }
                for (int row = 0; row < port.rows; ++row) {
                    candidate[static_cast<std::size_t>(row)] -=
                        static_cast<double>(coefficient) * vector[row];
                }
            }
        }
        const double magnitude =
            std::sqrt(std::max(0.0, normSquared(candidate)));
        if (!(magnitude > 1.0e-12 * std::max(1.0, original))) {
            continue;
        }
        for (double& value : candidate) value /= magnitude;
        basis.insert(basis.end(), candidate.begin(), candidate.end());
        ++accepted;
    }
    return basis;
}

double subspaceProjectorDifference(const LocalPortBasis& candidate,
                                   const LocalPortBasis& reference)
{
    if (candidate.rows != reference.rows) {
        return std::numeric_limits<double>::infinity();
    }
    const std::vector<double> left =
        euclideanOrthonormalBasis(candidate);
    const std::vector<double> right =
        euclideanOrthonormalBasis(reference);
    const int leftRank = candidate.rows > 0
        ? static_cast<int>(left.size()) / candidate.rows : 0;
    const int rightRank = reference.rows > 0
        ? static_cast<int>(right.size()) / reference.rows : 0;
    const auto projectionResidualSquared = [&](const std::vector<double>& from,
                                               int fromRank,
                                               const std::vector<double>& onto,
                                               int ontoRank) {
        long double total = 0.0L;
        std::vector<long double> residual(
            static_cast<std::size_t>(candidate.rows), 0.0L);
        for (int column = 0; column < fromRank; ++column) {
            const double* vector = from.data()
                + static_cast<std::size_t>(column * candidate.rows);
            for (int row = 0; row < candidate.rows; ++row) {
                residual[static_cast<std::size_t>(row)] = vector[row];
            }
            for (int ontoColumn = 0; ontoColumn < ontoRank; ++ontoColumn) {
                const double* ontoVector = onto.data()
                    + static_cast<std::size_t>(ontoColumn * candidate.rows);
                long double coefficient = 0.0L;
                for (int row = 0; row < candidate.rows; ++row) {
                    coefficient += residual[static_cast<std::size_t>(row)]
                        * ontoVector[row];
                }
                for (int row = 0; row < candidate.rows; ++row) {
                    residual[static_cast<std::size_t>(row)] -=
                        coefficient * ontoVector[row];
                }
            }
            for (long double value : residual) total += value * value;
        }
        return total;
    };
    // Form the two projection residuals explicitly. The algebraically
    // equivalent rankL + rankR - 2 ||Q_L^T Q_R||_F^2 expression loses
    // roughly sqrt(epsilon) accuracy when the two subspaces coincide.
    const long double squared = projectionResidualSquared(
        left, leftRank, right, rightRank)
        + projectionResidualSquared(right, rightRank, left, leftRank);
    return std::sqrt(static_cast<double>(std::max(0.0L, squared)))
        / std::sqrt(static_cast<double>(std::max(1, rightRank)));
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

template <typename T>
void hashVector(std::uint64_t& hash, const std::vector<T>& values)
{
    hashValue(hash, values.size());
    for (const T& value : values) hashValue(hash, value);
}

template <typename T>
std::uint64_t fingerprintVector(const std::vector<T>& values)
{
    std::uint64_t fingerprint = UINT64_C(1469598103934665603);
    hashVector(fingerprint, values);
    return fingerprint;
}

void hashSparse(std::uint64_t& hash, const SparseMatrix& matrix)
{
    hashValue(hash, matrix.n);
    hashVector(hash, matrix.rowPtr);
    hashVector(hash, matrix.colInd);
    hashVector(hash, matrix.values);
}

void hashCoupling(std::uint64_t& hash,
                  const std::vector<std::vector<std::pair<int, double>>>& rows)
{
    hashValue(hash, rows.size());
    for (const auto& row : rows) {
        hashValue(hash, row.size());
        for (const auto& entry : row) {
            hashValue(hash, entry.first);
            hashValue(hash, entry.second);
        }
    }
}

std::uint64_t fingerprintFile(const std::filesystem::path& path)
{
    if (path.empty()) return 0;
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error(
            "[Optimal port] Cannot fingerprint rank file: "
            + path.string());
    }
    std::uint64_t fingerprint = UINT64_C(1469598103934665603);
    std::uint64_t size = 0;
    char buffer[4096];
    while (input) {
        input.read(buffer, sizeof(buffer));
        const std::streamsize count = input.gcount();
        for (std::streamsize index = 0; index < count; ++index) {
            fingerprint ^= static_cast<unsigned char>(buffer[index]);
            fingerprint *= UINT64_C(1099511628211);
        }
        size += static_cast<std::uint64_t>(count);
    }
    hashValue(fingerprint, size);
    return fingerprint;
}

struct AdaptiveProductionRank {
    int historyRank = 64;
    int randomizedRank = 16;
    int residualRank = 4;
    bool adaptive = false;
};

AdaptiveProductionRank milestone8AdaptiveProductionRank(int interfaceId)
{
    AdaptiveProductionRank rank;
    if (interfaceId == 23 || interfaceId == 13
        || interfaceId == 18 || interfaceId == 4) {
        rank.historyRank = 256;
        rank.randomizedRank = 32;
        rank.residualRank = 8;
        rank.adaptive = true;
    } else if (interfaceId == 10) {
        rank.historyRank = 128;
        rank.randomizedRank = 32;
        rank.residualRank = 8;
        rank.adaptive = true;
    }
    return rank;
}

std::uint64_t milestone8AdaptiveProductionFingerprint()
{
    std::uint64_t fingerprint = UINT64_C(1469598103934665603);
    constexpr int policyVersion = 1;
    hashValue(fingerprint, policyVersion);
    for (int interfaceId = 0; interfaceId < 25; ++interfaceId) {
        const AdaptiveProductionRank rank =
            milestone8AdaptiveProductionRank(interfaceId);
        hashValue(fingerprint, interfaceId);
        hashValue(fingerprint, rank.historyRank);
        hashValue(fingerprint, rank.randomizedRank);
        hashValue(fingerprint, rank.residualRank);
    }
    return fingerprint;
}

std::uint64_t fingerprintDynamicPortOperator(const local::Model& model)
{
    std::uint64_t fingerprint = UINT64_C(1469598103934665603);
    hashValue(fingerprint, model.globalDofs);
    hashValue(fingerprint, model.interfaceDofs);
    hashVector(fingerprint, model.interfaceGlobalDofs);
    for (const local::InterfaceEntry& entry : model.interfaceEntries) {
        hashValue(fingerprint, entry.row);
        hashValue(fingerprint, entry.column);
        hashValue(fingerprint, entry.value);
    }
    for (const local::SubdomainModel& subdomain : model.subdomains) {
        hashValue(fingerprint, subdomain.subdomain);
        hashValue(fingerprint, subdomain.rank);
        hashVector(fingerprint, subdomain.interiorGlobalDofs);
        hashVector(fingerprint, subdomain.interfaceIndices);
        hashVector(fingerprint, subdomain.basis);
        hashVector(fingerprint, subdomain.reducedInterior);
        hashVector(fingerprint, subdomain.reducedInteriorInterface);
        hashVector(fingerprint, subdomain.reducedInterfaceInterior);
    }
    return fingerprint;
}

void writeOptimalPortTopologyAudit(
    const std::string& caseName,
    const std::vector<PortTopologyAudit>& rows,
    int globalDofs,
    int subdomains,
    int requestedRank,
    const std::string& requestedInnerSolver,
    const std::filesystem::path& outputDirectory)
{
    std::ofstream topology(
        outputDirectory / "optimal_port_topology_audit.csv");
    topology << "case,global_dofs,subdomains,physical_interfaces,"
        "interface_id,left_subdomain,right_subdomain,target_dofs,"
        "source_dofs,source_empty,mandatory_mode_count,"
        "heat_source_channel_count,external_boundary_channel_count,"
        "requested_rank,requested_inner_solver,estimated_workspace_bytes\n";
    std::ofstream memory(
        outputDirectory / "optimal_port_topology_memory_estimate.csv");
    memory << "case,interface_id,target_dofs,source_dofs,"
        "eigensolver_workspace_bytes,inner_solver_workspace_bytes,"
        "transfer_workspace_bytes,basis_storage_estimate_bytes,"
        "estimated_workspace_bytes\n";
    for (const PortTopologyAudit& row : rows) {
        topology << caseName << ',' << globalDofs << ',' << subdomains
            << ',' << rows.size() << ',' << row.interfaceId << ','
            << row.leftSubdomain << ',' << row.rightSubdomain << ','
            << row.targetDofs << ',' << row.sourceDofs << ','
            << (row.sourceEmpty ? 1 : 0) << ','
            << row.mandatoryModeCount << ','
            << row.heatSourceChannelCount << ','
            << row.externalBoundaryChannelCount << ','
            << requestedRank << ',' << requestedInnerSolver << ','
            << row.estimatedWorkspaceBytes << '\n';
        memory << caseName << ',' << row.interfaceId << ','
            << row.targetDofs << ',' << row.sourceDofs << ','
            << row.eigensolverWorkspaceBytes << ','
            << row.innerSolverWorkspaceBytes << ','
            << row.transferWorkspaceBytes << ','
            << row.basisStorageEstimateBytes << ','
            << row.estimatedWorkspaceBytes << '\n';
    }
}

bool sameSparse(const SparseMatrix& left, const SparseMatrix& right)
{
    return left.n == right.n && left.rowPtr == right.rowPtr
        && left.colInd == right.colInd && left.values == right.values;
}

bool sameCoupling(
    const std::vector<std::vector<std::pair<int, double>>>& left,
    const std::vector<std::vector<std::pair<int, double>>>& right)
{
    if (left.size() != right.size()) return false;
    for (std::size_t row = 0; row < left.size(); ++row) {
        if (left[row] != right[row]) return false;
    }
    return true;
}

double dot(const double* left, const double* right, int rows)
{
    long double result = 0.0L;
    for (int row = 0; row < rows; ++row) {
        result += static_cast<long double>(left[row]) * right[row];
    }
    return static_cast<double>(result);
}

struct MatrixPartition {
    std::vector<SparseMatrix> interior;
    std::vector<std::vector<std::vector<std::pair<int, double>>>> interiorInterface;
    std::vector<std::vector<std::vector<std::pair<int, double>>>> interfaceInterior;
    SparseMatrix interfaceBlock;
};

MatrixPartition partitionMatrix(
    const Mesh& mesh,
    const SparseMatrix& matrix,
    const ddm_schur::InterfacePartition& partition)
{
    MatrixPartition result;
    result.interfaceBlock = SparseMatrix(
        static_cast<int>(partition.interfaceGlobalDofs.size()));
    result.interior.reserve(partition.domains.size());
    result.interiorInterface.resize(partition.domains.size());
    result.interfaceInterior.resize(partition.domains.size());
    std::map<int, int> domainSlot;
    std::vector<int> globalToInterior(static_cast<std::size_t>(matrix.size()), -1);
    std::vector<int> globalToLocalGamma(static_cast<std::size_t>(matrix.size()), -1);
    for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
        const auto& domain = partition.domains[slot];
        domainSlot[domain.domainId] = static_cast<int>(slot);
        result.interior.emplace_back(static_cast<int>(domain.interiorGlobalDofs.size()));
        result.interiorInterface[slot].resize(domain.interiorGlobalDofs.size());
        result.interfaceInterior[slot].resize(domain.interfaceGlobalDofs.size());
        for (std::size_t local = 0; local < domain.interiorGlobalDofs.size(); ++local) {
            globalToInterior[static_cast<std::size_t>(
                domain.interiorGlobalDofs[local])] = static_cast<int>(local);
        }
        for (std::size_t local = 0; local < domain.interfaceGlobalDofs.size(); ++local) {
            globalToLocalGamma[static_cast<std::size_t>(
                domain.interfaceGlobalDofs[local])] = static_cast<int>(local);
        }
    }
    matrix.forEachEntry([&](int row, int column, double value) {
        const int gammaRow = partition.globalToInterface[static_cast<std::size_t>(row)];
        const int gammaColumn = partition.globalToInterface[static_cast<std::size_t>(column)];
        if (gammaRow >= 0 && gammaColumn >= 0) {
            result.interfaceBlock.add(gammaRow, gammaColumn, value);
            return;
        }
        if (gammaRow < 0 && gammaColumn < 0) {
            const int slot = domainSlot.at(mesh.nodes[static_cast<std::size_t>(row)].subdomain);
            if (slot != domainSlot.at(mesh.nodes[static_cast<std::size_t>(column)].subdomain)) {
                throw std::runtime_error("[Local transient] Cross-domain interior entry detected.");
            }
            result.interior[static_cast<std::size_t>(slot)].add(
                globalToInterior[static_cast<std::size_t>(row)],
                globalToInterior[static_cast<std::size_t>(column)], value);
            return;
        }
        if (gammaRow < 0) {
            const int slot = domainSlot.at(mesh.nodes[static_cast<std::size_t>(row)].subdomain);
            const int localGamma = globalToLocalGamma[static_cast<std::size_t>(column)];
            if (localGamma < 0 || mesh.nodes[static_cast<std::size_t>(row)].subdomain
                    != mesh.nodes[static_cast<std::size_t>(column)].subdomain) {
                throw std::runtime_error("[Local transient] Invalid interior-interface ownership.");
            }
            result.interiorInterface[static_cast<std::size_t>(slot)]
                [static_cast<std::size_t>(globalToInterior[static_cast<std::size_t>(row)])]
                .push_back({localGamma, value});
            return;
        }
        const int slot = domainSlot.at(mesh.nodes[static_cast<std::size_t>(column)].subdomain);
        const int localGamma = globalToLocalGamma[static_cast<std::size_t>(row)];
        if (localGamma < 0 || mesh.nodes[static_cast<std::size_t>(row)].subdomain
                != mesh.nodes[static_cast<std::size_t>(column)].subdomain) {
            throw std::runtime_error("[Local transient] Invalid interface-interior ownership.");
        }
        result.interfaceInterior[static_cast<std::size_t>(slot)]
            [static_cast<std::size_t>(localGamma)]
            .push_back({globalToInterior[static_cast<std::size_t>(column)], value});
    });
    for (SparseMatrix& block : result.interior) block.finalizeCsr();
    result.interfaceBlock.finalizeCsr();
    return result;
}

std::vector<double> projectSparse(const SparseMatrix& matrix,
                                  const std::vector<double>& basis,
                                  int rows,
                                  int rank)
{
    // Sparse-dense product Y=A*V followed by V^T*Y.  AII is never
    // densified and the matrix traversal is independent of the rank.
    std::vector<double> image(static_cast<std::size_t>(rows * rank), 0.0);
#ifdef USE_MKL_PARDISO
    parallelFor(static_cast<std::size_t>(rank), [&](std::size_t modeIndex) {
        const int mode = static_cast<int>(modeIndex);
        const double* vector = basis.data()
            + static_cast<std::size_t>(mode * rows);
        double* output = image.data() + static_cast<std::size_t>(mode * rows);
        for (int row = 0; row < rows; ++row) {
            double value = 0.0;
            for (int offset = matrix.rowPtr[static_cast<std::size_t>(row)];
                 offset < matrix.rowPtr[static_cast<std::size_t>(row + 1)];
                 ++offset) {
                value += matrix.values[static_cast<std::size_t>(offset)]
                    * vector[static_cast<std::size_t>(
                        matrix.colInd[static_cast<std::size_t>(offset)])];
            }
            output[static_cast<std::size_t>(row)] = value;
        }
    });
    std::vector<double> projected(
        static_cast<std::size_t>(rank * rank), 0.0);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
        rank, rank, rows, 1.0, basis.data(), rows,
        image.data(), rows, 0.0, projected.data(), rank);
    std::vector<double> result(static_cast<std::size_t>(rank * rank), 0.0);
    for (int left = 0; left < rank; ++left) {
        for (int right = 0; right < rank; ++right) {
            result[static_cast<std::size_t>(left * rank + right)] =
                projected[static_cast<std::size_t>(right * rank + left)];
        }
    }
    return result;
#else
    matrix.forEachEntry([&](int row, int column, double value) {
        for (int mode = 0; mode < rank; ++mode) {
            image[static_cast<std::size_t>(mode * rows + row)] += value
                * basis[static_cast<std::size_t>(mode * rows + column)];
        }
    });
    std::vector<double> result(static_cast<std::size_t>(rank * rank), 0.0);
    for (int left = 0; left < rank; ++left) {
        for (int right = 0; right < rank; ++right) {
            result[static_cast<std::size_t>(left * rank + right)] = dot(
                basis.data() + static_cast<std::size_t>(left * rows),
                image.data() + static_cast<std::size_t>(right * rows), rows);
        }
    }
    return result;
#endif
}

std::vector<double> projectInteriorInterface(
    const std::vector<std::vector<std::pair<int, double>>>& rows,
    const std::vector<double>& basis,
    int interiorRows,
    int rank,
    int gammaRows)
{
    std::vector<double> result(static_cast<std::size_t>(rank * gammaRows), 0.0);
    parallelFor(static_cast<std::size_t>(rank), [&](std::size_t modeIndex) {
        const int mode = static_cast<int>(modeIndex);
        double* output = result.data()
            + static_cast<std::size_t>(mode * gammaRows);
        const double* vector = basis.data()
            + static_cast<std::size_t>(mode * interiorRows);
        for (int row = 0; row < interiorRows; ++row) {
            for (const auto& entry : rows[static_cast<std::size_t>(row)]) {
                output[static_cast<std::size_t>(entry.first)] +=
                    vector[static_cast<std::size_t>(row)] * entry.second;
            }
        }
    });
    return result;
}

std::vector<double> projectInterfaceInterior(
    const std::vector<std::vector<std::pair<int, double>>>& rows,
    const std::vector<double>& basis,
    int interiorRows,
    int rank)
{
    std::vector<double> result(static_cast<std::size_t>(rows.size() * rank), 0.0);
    parallelFor(rows.size(), [&](std::size_t gamma) {
        for (const auto& entry : rows[gamma]) {
            for (int mode = 0; mode < rank; ++mode) {
                result[gamma * static_cast<std::size_t>(rank) + mode] +=
                    entry.second * basis[static_cast<std::size_t>(
                        mode * interiorRows + entry.first)];
            }
        }
    });
    return result;
}

std::vector<double> orthonormalTraceBasis(
    int rows,
    std::vector<double> candidates,
    int columns,
    double tolerance,
    int maximumRank)
{
#ifdef USE_MKL_PARDISO
    if (rows <= 0) return {};
    const int requestedRank = maximumRank > 0
        ? maximumRank : columns + 1;
    std::vector<double> basis(static_cast<std::size_t>(rows),
        1.0 / std::sqrt(static_cast<double>(rows)));
    if (requestedRank <= 1 || columns <= 0) return basis;
    if (candidates.size() != static_cast<std::size_t>(rows * columns)) {
        throw std::runtime_error(
            "[Local ROM] Interface trace candidate dimensions are invalid.");
    }
    std::vector<double> coefficients(static_cast<std::size_t>(columns), 0.0);
    for (int pass = 0; pass < 2; ++pass) {
        cblas_dgemv(CblasColMajor, CblasTrans, rows, columns, 1.0,
            candidates.data(), rows, basis.data(), 1, 0.0,
            coefficients.data(), 1);
        cblas_dger(CblasColMajor, rows, columns, -1.0,
            basis.data(), 1, coefficients.data(), 1,
            candidates.data(), rows);
    }
    std::vector<lapack_int> pivots(static_cast<std::size_t>(columns), 0);
    std::vector<double> tau(static_cast<std::size_t>(
        std::min(rows, columns)), 0.0);
    const lapack_int factorInfo = LAPACKE_dgeqp3(
        LAPACK_COL_MAJOR, static_cast<lapack_int>(rows),
        static_cast<lapack_int>(columns), candidates.data(),
        static_cast<lapack_int>(rows), pivots.data(), tau.data());
    if (factorInfo != 0) {
        throw std::runtime_error(
            "[Local ROM] Interface trace rank-revealing QR failed.");
    }
    const int maximumAccepted = std::min({
        requestedRank - 1, rows, columns});
    int accepted = 0;
    while (accepted < maximumAccepted) {
        const double diagonal = std::abs(candidates[static_cast<std::size_t>(
            accepted * rows + accepted)]);
        if (!(diagonal > tolerance) || !std::isfinite(diagonal)) break;
        ++accepted;
    }
    if (accepted > 0) {
        const lapack_int basisInfo = LAPACKE_dorgqr(
            LAPACK_COL_MAJOR, static_cast<lapack_int>(rows),
            static_cast<lapack_int>(accepted),
            static_cast<lapack_int>(accepted), candidates.data(),
            static_cast<lapack_int>(rows), tau.data());
        if (basisInfo != 0) {
            throw std::runtime_error(
                "[Local ROM] Interface trace basis generation failed.");
        }
        basis.insert(basis.end(), candidates.begin(), candidates.begin()
            + static_cast<std::ptrdiff_t>(rows * accepted));
    }
    return basis;
#else
    std::vector<double> basis;
    auto append = [&](std::vector<double> candidate) {
        const int columns = static_cast<int>(basis.size()) / std::max(1, rows);
        if (maximumRank > 0 && columns >= maximumRank) return;
        for (int pass = 0; pass < 2; ++pass) {
            const int currentColumns = static_cast<int>(basis.size()) / std::max(1, rows);
            for (int mode = 0; mode < currentColumns; ++mode) {
                const double* q = basis.data() + static_cast<std::size_t>(mode * rows);
                const double coefficient = dot(q, candidate.data(), rows);
                for (int row = 0; row < rows; ++row) candidate[static_cast<std::size_t>(row)] -=
                    coefficient * q[row];
            }
        }
        const double magnitude = std::sqrt(normSquared(candidate));
        if (!(magnitude > tolerance)) return;
        for (double& value : candidate) value /= magnitude;
        basis.insert(basis.end(), candidate.begin(), candidate.end());
    };
    append(std::vector<double>(static_cast<std::size_t>(rows), 1.0));
    for (int column = 0; column < columns; ++column) {
        const auto begin = candidates.begin() + static_cast<std::ptrdiff_t>(
            static_cast<std::size_t>(column) * rows);
        append(std::vector<double>(begin, begin + rows));
    }
    return basis;
#endif
}

std::vector<double> traceBasis(const ddm_schur::DomainBlocks& domain,
                               const std::vector<std::vector<double>>& traces,
                               double tolerance,
                               int maximumRank)
{
    const int rows = static_cast<int>(domain.interfaceGlobalDofs.size());
    const int columns = static_cast<int>(traces.size());
    std::vector<double> candidates(
        static_cast<std::size_t>(rows) * columns, 0.0);
    parallelFor(static_cast<std::size_t>(columns), [&](std::size_t columnIndex) {
        const auto& globalTrace = traces[columnIndex];
        double* local = candidates.data() + columnIndex
            * static_cast<std::size_t>(rows);
        for (int row = 0; row < rows; ++row) {
            local[static_cast<std::size_t>(row)] = globalTrace[
                static_cast<std::size_t>(domain.interfaceGlobalDofs[
                    static_cast<std::size_t>(row)])];
        }
    });
    return orthonormalTraceBasis(
        rows, std::move(candidates), columns, tolerance, maximumRank);
}

struct OperatorCoarseTraceBuild {
    std::vector<std::vector<double>> traces;
    int aggregates = 0;
    int krylovIterations = 0;
    int krylovMaximumIterations = 0;
    double krylovMaximumRelativeResidual = 0.0;
    double setupSeconds = 0.0;
};

std::vector<double> operatorCoarseTraceBasis(
    const Mesh& mesh,
    const ddm_schur::DomainBlocks& domain,
    const std::vector<std::vector<double>>& coarseTraces,
    double tolerance,
    int maximumRank);

OperatorCoarseTraceBuild buildOperatorCoarseTraces(
    const Mesh& mesh,
    const ThermalDescriptorSystem& descriptor,
    const ddm_schur::InterfacePartition& partition,
    const MatrixPartition& k,
    const Options& options)
{
    const auto start = Clock::now();
    OperatorCoarseTraceBuild result;
    const int interfaceDofs = static_cast<int>(
        partition.interfaceGlobalDofs.size());
    const int localCoarseRank = options.interfaceExcitationRank > 0
        ? std::min(12, options.interfaceExcitationRank) : 12;
    std::vector<std::vector<double>> localTraceBases(
        partition.domains.size());
    std::vector<int> rankOffsets(partition.domains.size() + 1, 0);
    for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
        localTraceBases[slot] = operatorCoarseTraceBasis(
            mesh, partition.domains[slot], {}, options.rankTolerance,
            localCoarseRank);
        const int localRows = static_cast<int>(
            partition.domains[slot].interfaceGlobalDofs.size());
        const int localRank = localRows > 0
            ? static_cast<int>(localTraceBases[slot].size()) / localRows : 0;
        rankOffsets[slot + 1] = rankOffsets[slot] + localRank;
    }
    result.aggregates = rankOffsets.back();
    if (result.aggregates <= 0) {
        throw std::runtime_error(
            "[Operator-coarse traces] Empty Schur coarse space.");
    }
    std::vector<double> z(static_cast<std::size_t>(
        interfaceDofs * result.aggregates), 0.0);
    for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
        const auto& domain = partition.domains[slot];
        const int localRows = static_cast<int>(domain.interfaceGlobalDofs.size());
        const int localRank = rankOffsets[slot + 1] - rankOffsets[slot];
        for (int mode = 0; mode < localRank; ++mode) {
            double* globalColumn = z.data() + static_cast<std::size_t>(
                (rankOffsets[slot] + mode) * interfaceDofs);
            const double* localColumn = localTraceBases[slot].data()
                + static_cast<std::size_t>(mode * localRows);
            for (int row = 0; row < localRows; ++row) {
                const int gamma = partition.globalToInterface[
                    static_cast<std::size_t>(
                        domain.interfaceGlobalDofs[static_cast<std::size_t>(row)])];
                globalColumn[static_cast<std::size_t>(gamma)] =
                    localColumn[static_cast<std::size_t>(row)];
            }
        }
    }

    const int requestedWorkers = std::max(1, options.localSolveThreads);
    const int workers = std::max(1, std::min({
        requestedWorkers, static_cast<int>(partition.domains.size()),
        static_cast<int>(solverParallelWorkers())}));
    const int factorThreads = std::max(
        1, options.localPardisoThreads / workers);
    std::vector<std::unique_ptr<SubdomainDirectSolver>> factors(
        partition.domains.size());
    std::vector<std::unique_ptr<SubdomainDirectSolver>> blockFactors(
        partition.domains.size());
    std::vector<std::exception_ptr> factorErrors(partition.domains.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(workers) if(workers > 1)
#endif
    for (int slotIndex = 0;
         slotIndex < static_cast<int>(partition.domains.size()); ++slotIndex) {
        try {
            ScopedDirectSolverMklThreads threads(factorThreads);
            factors[static_cast<std::size_t>(slotIndex)] =
                std::make_unique<SubdomainDirectSolver>(
                    static_cast<int>(partition.domains[
                        static_cast<std::size_t>(slotIndex)]
                            .interiorGlobalDofs.size()),
                    sparseMatrixEntries(k.interior[
                        static_cast<std::size_t>(slotIndex)]));
            const auto& domain = partition.domains[
                static_cast<std::size_t>(slotIndex)];
            std::vector<MatrixEntry> blockEntries;
            blockEntries.reserve(domain.fullBlockEntries.size());
            for (const auto& entry : domain.fullBlockEntries) {
                blockEntries.push_back(
                    {entry.row, entry.col, entry.value});
            }
            blockFactors[static_cast<std::size_t>(slotIndex)] =
                std::make_unique<SubdomainDirectSolver>(
                    static_cast<int>(domain.interiorGlobalDofs.size()
                        + domain.interfaceGlobalDofs.size()),
                    blockEntries);
        } catch (...) {
            factorErrors[static_cast<std::size_t>(slotIndex)] =
                std::current_exception();
        }
    }
    for (const auto& error : factorErrors) {
        if (error) std::rethrow_exception(error);
    }

    const auto applySchur = [&](const std::vector<double>& vectors,
                                int columns) {
        if (vectors.size() != static_cast<std::size_t>(
                interfaceDofs * columns)) {
            throw std::runtime_error(
                "[Operator-coarse traces] Schur input dimensions are invalid.");
        }
        std::vector<double> image(static_cast<std::size_t>(
            interfaceDofs * columns), 0.0);
        parallelForCoarse(static_cast<std::size_t>(columns),
            [&](std::size_t column) {
                std::vector<double> vector(
                    vectors.begin() + static_cast<std::ptrdiff_t>(
                        column * interfaceDofs),
                    vectors.begin() + static_cast<std::ptrdiff_t>(
                        (column + 1) * interfaceDofs));
                const std::vector<double> blockImage =
                    k.interfaceBlock.multiply(vector);
                std::copy(blockImage.begin(), blockImage.end(),
                    image.begin() + static_cast<std::ptrdiff_t>(
                        column * interfaceDofs));
            });
        std::vector<std::exception_ptr> applyErrors(
            partition.domains.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(workers) if(workers > 1)
#endif
        for (int slotIndex = 0;
             slotIndex < static_cast<int>(partition.domains.size()); ++slotIndex) {
            try {
                const std::size_t slot = static_cast<std::size_t>(slotIndex);
                const auto& domain = partition.domains[slot];
                const int interiorRows = static_cast<int>(
                    domain.interiorGlobalDofs.size());
                constexpr int blockColumns = 24;
                for (int first = 0; first < columns; first += blockColumns) {
                    const int activeColumns = std::min(
                        blockColumns, columns - first);
                    std::vector<double> rhs(static_cast<std::size_t>(
                        interiorRows * activeColumns), 0.0);
                    for (int localRow = 0; localRow < interiorRows; ++localRow) {
                        for (const auto& [localGamma, value] :
                             k.interiorInterface[slot][
                                 static_cast<std::size_t>(localRow)]) {
                            const int gamma = partition.globalToInterface[
                                static_cast<std::size_t>(
                                    domain.interfaceGlobalDofs[
                                        static_cast<std::size_t>(localGamma)])];
                            for (int column = 0; column < activeColumns; ++column) {
                                rhs[static_cast<std::size_t>(
                                    column * interiorRows + localRow)] -= value
                                    * vectors[static_cast<std::size_t>(
                                        (first + column) * interfaceDofs + gamma)];
                            }
                        }
                    }
                    std::vector<double> interior;
                    factors[slot]->solveMultiple(
                        rhs, activeColumns, interior, factorThreads);
                    for (int localGamma = 0;
                         localGamma < static_cast<int>(
                             domain.interfaceGlobalDofs.size());
                         ++localGamma) {
                        const int gamma = partition.globalToInterface[
                            static_cast<std::size_t>(domain.interfaceGlobalDofs[
                                static_cast<std::size_t>(localGamma)])];
                        for (int column = 0; column < activeColumns; ++column) {
                            long double value = 0.0L;
                            for (const auto& [interiorColumn, coefficient] :
                                 k.interfaceInterior[slot][
                                     static_cast<std::size_t>(localGamma)]) {
                                value += static_cast<long double>(coefficient)
                                    * interior[static_cast<std::size_t>(
                                        column * interiorRows
                                        + interiorColumn)];
                            }
                            image[static_cast<std::size_t>(
                                (first + column) * interfaceDofs + gamma)] +=
                                static_cast<double>(value);
                        }
                    }
                }
            } catch (...) {
                applyErrors[static_cast<std::size_t>(slotIndex)] =
                    std::current_exception();
            }
        }
        for (const auto& error : applyErrors) {
            if (error) std::rethrow_exception(error);
        }
        return image;
    };

    const std::vector<double> sz = applySchur(z, result.aggregates);

    std::vector<double> coarseSchur(static_cast<std::size_t>(
        result.aggregates * result.aggregates), 0.0);
#ifdef USE_MKL_PARDISO
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
        result.aggregates, result.aggregates, interfaceDofs, 1.0,
        z.data(), interfaceDofs, sz.data(), interfaceDofs, 0.0,
        coarseSchur.data(), result.aggregates);
#else
    for (int column = 0; column < result.aggregates; ++column) {
        for (int row = 0; row < result.aggregates; ++row) {
            coarseSchur[static_cast<std::size_t>(
                column * result.aggregates + row)] = dot(
                    z.data() + static_cast<std::size_t>(row * interfaceDofs),
                    sz.data() + static_cast<std::size_t>(column * interfaceDofs),
                    interfaceDofs);
        }
    }
#endif
    std::vector<MatrixEntry> coarseEntries;
    for (int row = 0; row < result.aggregates; ++row) {
        for (int column = row; column < result.aggregates; ++column) {
            const double value = 0.5 * (
                coarseSchur[static_cast<std::size_t>(
                    column * result.aggregates + row)]
                + coarseSchur[static_cast<std::size_t>(
                    row * result.aggregates + column)]);
            if (value != 0.0 || row == column) {
                coarseEntries.push_back({row, column, value});
            }
        }
    }
    ScopedDirectSolverMklThreads coarseThreads(1);
    SubdomainDirectSolver coarseFactor(result.aggregates, coarseEntries);

    const auto condensedRhs = [&](const std::vector<double>& globalRhs,
                                  int channels) {
        std::vector<double> condensed(static_cast<std::size_t>(
            interfaceDofs * channels), 0.0);
        for (int channel = 0; channel < channels; ++channel) {
            for (int gamma = 0; gamma < interfaceDofs; ++gamma) {
                condensed[static_cast<std::size_t>(
                    channel * interfaceDofs + gamma)] = globalRhs[
                        static_cast<std::size_t>(channel * descriptor.dofs
                            + partition.interfaceGlobalDofs[
                                static_cast<std::size_t>(gamma)])];
            }
        }
        std::vector<std::exception_ptr> errors(partition.domains.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(workers) if(workers > 1)
#endif
        for (int slotIndex = 0;
             slotIndex < static_cast<int>(partition.domains.size()); ++slotIndex) {
            try {
                const std::size_t slot = static_cast<std::size_t>(slotIndex);
                const auto& domain = partition.domains[slot];
                const int interiorRows = static_cast<int>(
                    domain.interiorGlobalDofs.size());
                std::vector<double> rhs(static_cast<std::size_t>(
                    interiorRows * channels), 0.0);
                for (int channel = 0; channel < channels; ++channel) {
                    for (int row = 0; row < interiorRows; ++row) {
                        rhs[static_cast<std::size_t>(
                            channel * interiorRows + row)] = globalRhs[
                                static_cast<std::size_t>(channel * descriptor.dofs
                                    + domain.interiorGlobalDofs[
                                        static_cast<std::size_t>(row)])];
                    }
                }
                std::vector<double> interior;
                factors[slot]->solveMultiple(
                    rhs, channels, interior, factorThreads);
                for (int localGamma = 0;
                     localGamma < static_cast<int>(domain.interfaceGlobalDofs.size());
                     ++localGamma) {
                    const int gamma = partition.globalToInterface[
                        static_cast<std::size_t>(domain.interfaceGlobalDofs[
                            static_cast<std::size_t>(localGamma)])];
                    for (int channel = 0; channel < channels; ++channel) {
                        long double correction = 0.0L;
                        for (const auto& [interiorColumn, coefficient] :
                             k.interfaceInterior[slot][
                                 static_cast<std::size_t>(localGamma)]) {
                            correction += static_cast<long double>(coefficient)
                                * interior[static_cast<std::size_t>(
                                    channel * interiorRows + interiorColumn)];
                        }
                        condensed[static_cast<std::size_t>(
                            channel * interfaceDofs + gamma)] -=
                            static_cast<double>(correction);
                    }
                }
            } catch (...) {
                errors[static_cast<std::size_t>(slotIndex)] =
                    std::current_exception();
            }
        }
        for (const auto& error : errors) {
            if (error) std::rethrow_exception(error);
        }
        return condensed;
    };
    const auto solveCoarse = [&](const std::vector<double>& condensed,
                                 int channels) {
        std::vector<double> projected(static_cast<std::size_t>(
            result.aggregates * channels), 0.0);
#ifdef USE_MKL_PARDISO
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
            result.aggregates, channels, interfaceDofs, 1.0,
            z.data(), interfaceDofs, condensed.data(), interfaceDofs,
            0.0, projected.data(), result.aggregates);
#else
        for (int channel = 0; channel < channels; ++channel) {
            for (int mode = 0; mode < result.aggregates; ++mode) {
                projected[static_cast<std::size_t>(
                    channel * result.aggregates + mode)] = dot(
                        z.data() + static_cast<std::size_t>(mode * interfaceDofs),
                        condensed.data() + static_cast<std::size_t>(
                            channel * interfaceDofs), interfaceDofs);
            }
        }
#endif
        std::vector<double> coordinates;
        coarseFactor.solveMultiple(projected, channels, coordinates, 1);
        std::vector<double> gamma(static_cast<std::size_t>(
            interfaceDofs * channels), 0.0);
#ifdef USE_MKL_PARDISO
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
            interfaceDofs, channels, result.aggregates, 1.0,
            z.data(), interfaceDofs, coordinates.data(), result.aggregates,
            0.0, gamma.data(), interfaceDofs);
#else
        for (int channel = 0; channel < channels; ++channel) {
            for (int mode = 0; mode < result.aggregates; ++mode) {
                const double coefficient = coordinates[static_cast<std::size_t>(
                    channel * result.aggregates + mode)];
                for (int row = 0; row < interfaceDofs; ++row) {
                    gamma[static_cast<std::size_t>(
                        channel * interfaceDofs + row)] += coefficient
                        * z[static_cast<std::size_t>(mode * interfaceDofs + row)];
                }
            }
        }
#endif
        return gamma;
    };
    const auto coarseCoordinates = [&](const std::vector<double>& input,
                                       const std::vector<double>& leftBasis,
                                       int channels) {
        std::vector<double> projected(static_cast<std::size_t>(
            result.aggregates * channels), 0.0);
#ifdef USE_MKL_PARDISO
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
            result.aggregates, channels, interfaceDofs, 1.0,
            leftBasis.data(), interfaceDofs, input.data(), interfaceDofs,
            0.0, projected.data(), result.aggregates);
#else
        for (int channel = 0; channel < channels; ++channel) {
            for (int mode = 0; mode < result.aggregates; ++mode) {
                projected[static_cast<std::size_t>(
                    channel * result.aggregates + mode)] = dot(
                        leftBasis.data() + static_cast<std::size_t>(
                            mode * interfaceDofs),
                        input.data() + static_cast<std::size_t>(
                            channel * interfaceDofs), interfaceDofs);
            }
        }
#endif
        std::vector<double> coordinates;
        coarseFactor.solveMultiple(projected, channels, coordinates, 1);
        return coordinates;
    };
    const auto expandCoarse = [&](const std::vector<double>& coordinates,
                                  const std::vector<double>& basis,
                                  int channels) {
        std::vector<double> expanded(static_cast<std::size_t>(
            interfaceDofs * channels), 0.0);
#ifdef USE_MKL_PARDISO
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
            interfaceDofs, channels, result.aggregates, 1.0,
            basis.data(), interfaceDofs, coordinates.data(),
            result.aggregates, 0.0, expanded.data(), interfaceDofs);
#else
        for (int channel = 0; channel < channels; ++channel) {
            for (int mode = 0; mode < result.aggregates; ++mode) {
                const double coefficient = coordinates[static_cast<std::size_t>(
                    channel * result.aggregates + mode)];
                for (int row = 0; row < interfaceDofs; ++row) {
                    expanded[static_cast<std::size_t>(
                        channel * interfaceDofs + row)] += coefficient
                        * basis[static_cast<std::size_t>(
                            mode * interfaceDofs + row)];
                }
            }
        }
#endif
        return expanded;
    };
    const auto applyPreconditioner = [&](const std::vector<double>& residual,
                                         int channels) {
        const std::vector<double> rightCoordinates =
            coarseCoordinates(residual, z, channels);
        const std::vector<double> coarseCorrection =
            expandCoarse(rightCoordinates, z, channels);
        const std::vector<double> coarseImage =
            expandCoarse(rightCoordinates, sz, channels);
        std::vector<double> projectedResidual(residual.size(), 0.0);
        for (std::size_t index = 0; index < residual.size(); ++index) {
            projectedResidual[index] = residual[index] - coarseImage[index];
        }
        std::vector<double> localCorrection(residual.size(), 0.0);
        std::vector<std::exception_ptr> errors(partition.domains.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(workers) if(workers > 1)
#endif
        for (int slotIndex = 0;
             slotIndex < static_cast<int>(partition.domains.size());
             ++slotIndex) {
            try {
                const std::size_t slot = static_cast<std::size_t>(slotIndex);
                const auto& domain = partition.domains[slot];
                const int interiorRows = static_cast<int>(
                    domain.interiorGlobalDofs.size());
                const int gammaRows = static_cast<int>(
                    domain.interfaceGlobalDofs.size());
                const int localRows = interiorRows + gammaRows;
                std::vector<double> localRhs(static_cast<std::size_t>(
                    localRows * channels), 0.0);
                for (int channel = 0; channel < channels; ++channel) {
                    for (int localGamma = 0; localGamma < gammaRows;
                         ++localGamma) {
                        const int gamma = partition.globalToInterface[
                            static_cast<std::size_t>(
                                domain.interfaceGlobalDofs[
                                    static_cast<std::size_t>(localGamma)])];
                        localRhs[static_cast<std::size_t>(
                            channel * localRows + interiorRows + localGamma)] =
                            projectedResidual[static_cast<std::size_t>(
                                channel * interfaceDofs + gamma)];
                    }
                }
                std::vector<double> localSolution;
                blockFactors[slot]->solveMultiple(
                    localRhs, channels, localSolution, factorThreads);
                for (int channel = 0; channel < channels; ++channel) {
                    for (int localGamma = 0; localGamma < gammaRows;
                         ++localGamma) {
                        const int gamma = partition.globalToInterface[
                            static_cast<std::size_t>(
                                domain.interfaceGlobalDofs[
                                    static_cast<std::size_t>(localGamma)])];
                        localCorrection[static_cast<std::size_t>(
                            channel * interfaceDofs + gamma)] =
                            localSolution[static_cast<std::size_t>(
                                channel * localRows + interiorRows
                                + localGamma)];
                    }
                }
            } catch (...) {
                errors[static_cast<std::size_t>(slotIndex)] =
                    std::current_exception();
            }
        }
        for (const auto& error : errors) {
            if (error) std::rethrow_exception(error);
        }
        const std::vector<double> leftCoordinates =
            coarseCoordinates(localCorrection, sz, channels);
        const std::vector<double> leftCorrection =
            expandCoarse(leftCoordinates, z, channels);
        std::vector<double> preconditioned(residual.size(), 0.0);
        for (std::size_t index = 0; index < residual.size(); ++index) {
            preconditioned[index] = coarseCorrection[index]
                + localCorrection[index] - leftCorrection[index];
        }
        return preconditioned;
    };
    const auto solveSchur = [&](const std::vector<double>& condensed,
                                int channels) {
        if (channels <= 0 || condensed.size() != static_cast<std::size_t>(
                interfaceDofs * channels)) {
            throw std::runtime_error(
                "[Operator-coarse traces] Condensed RHS dimensions are invalid.");
        }
        std::vector<double> solution = solveCoarse(condensed, channels);
        const std::vector<double> initialImage =
            applySchur(solution, channels);
        std::vector<double> residual(condensed.size(), 0.0);
        std::vector<double> rightNorm(static_cast<std::size_t>(channels), 0.0);
        std::vector<double> relativeResidual(
            static_cast<std::size_t>(channels), 0.0);
        std::vector<int> channelIterations(
            static_cast<std::size_t>(channels), 0);
        std::vector<unsigned char> active(
            static_cast<std::size_t>(channels), 0);
        int activeChannels = 0;
        const double tolerance = options.interfaceTolerance;
        for (int channel = 0; channel < channels; ++channel) {
            double* r = residual.data()
                + static_cast<std::size_t>(channel * interfaceDofs);
            const double* rhs = condensed.data()
                + static_cast<std::size_t>(channel * interfaceDofs);
            const double* image = initialImage.data()
                + static_cast<std::size_t>(channel * interfaceDofs);
            for (int row = 0; row < interfaceDofs; ++row) {
                r[row] = rhs[row] - image[row];
            }
            rightNorm[static_cast<std::size_t>(channel)] =
                std::sqrt(std::max(0.0,
                    dot(rhs, rhs, interfaceDofs)));
            relativeResidual[static_cast<std::size_t>(channel)] =
                std::sqrt(std::max(0.0, dot(r, r, interfaceDofs)))
                / std::max(1.0e-300,
                    rightNorm[static_cast<std::size_t>(channel)]);
            if (relativeResidual[static_cast<std::size_t>(channel)]
                    > tolerance) {
                active[static_cast<std::size_t>(channel)] = 1;
                ++activeChannels;
            }
        }
        std::vector<double> preconditioned =
            applyPreconditioner(residual, channels);
        std::vector<double> direction = preconditioned;
        std::vector<double> rho(static_cast<std::size_t>(channels), 0.0);
        for (int channel = 0; channel < channels; ++channel) {
            if (!active[static_cast<std::size_t>(channel)]) {
                std::fill_n(direction.begin() + static_cast<std::ptrdiff_t>(
                    channel * interfaceDofs), interfaceDofs, 0.0);
                continue;
            }
            rho[static_cast<std::size_t>(channel)] = dot(
                residual.data() + static_cast<std::size_t>(
                    channel * interfaceDofs),
                preconditioned.data() + static_cast<std::size_t>(
                    channel * interfaceDofs), interfaceDofs);
        }
        for (int iteration = 1;
             activeChannels > 0 && iteration <= options.interfaceMaxIterations;
             ++iteration) {
            const std::vector<double> image =
                applySchur(direction, channels);
            for (int channel = 0; channel < channels; ++channel) {
                if (!active[static_cast<std::size_t>(channel)]) continue;
                double* x = solution.data() + static_cast<std::size_t>(
                    channel * interfaceDofs);
                double* r = residual.data() + static_cast<std::size_t>(
                    channel * interfaceDofs);
                const double* p = direction.data() + static_cast<std::size_t>(
                    channel * interfaceDofs);
                const double* ap = image.data() + static_cast<std::size_t>(
                    channel * interfaceDofs);
                const double denominator = dot(p, ap, interfaceDofs);
                const double numerator = rho[static_cast<std::size_t>(channel)];
                if (!(denominator > 0.0) || !(numerator > 0.0)
                    || !std::isfinite(denominator)
                    || !std::isfinite(numerator)) {
                    throw std::runtime_error(
                        "[Operator-coarse traces] Construction PCG lost SPD.");
                }
                const double alpha = numerator / denominator;
                for (int row = 0; row < interfaceDofs; ++row) {
                    x[row] += alpha * p[row];
                    r[row] -= alpha * ap[row];
                }
                relativeResidual[static_cast<std::size_t>(channel)] =
                    std::sqrt(std::max(0.0, dot(r, r, interfaceDofs)))
                    / std::max(1.0e-300,
                        rightNorm[static_cast<std::size_t>(channel)]);
                channelIterations[static_cast<std::size_t>(channel)] =
                    iteration;
                if (relativeResidual[static_cast<std::size_t>(channel)]
                        <= tolerance) {
                    active[static_cast<std::size_t>(channel)] = 0;
                    --activeChannels;
                }
            }
            if (activeChannels == 0) break;
            std::vector<double> nextPreconditioned =
                applyPreconditioner(residual, channels);
            for (int channel = 0; channel < channels; ++channel) {
                double* p = direction.data() + static_cast<std::size_t>(
                    channel * interfaceDofs);
                if (!active[static_cast<std::size_t>(channel)]) {
                    std::fill_n(p, interfaceDofs, 0.0);
                    continue;
                }
                const double* r = residual.data() + static_cast<std::size_t>(
                    channel * interfaceDofs);
                const double* next = nextPreconditioned.data()
                    + static_cast<std::size_t>(channel * interfaceDofs);
                const double nextRho = dot(r, next, interfaceDofs);
                const double previousRho = rho[static_cast<std::size_t>(channel)];
                if (!(nextRho > 0.0) || !(previousRho > 0.0)
                    || !std::isfinite(nextRho)) {
                    throw std::runtime_error(
                        "[Operator-coarse traces] Construction preconditioner "
                        "is not positive definite.");
                }
                const double beta = nextRho / previousRho;
                for (int row = 0; row < interfaceDofs; ++row) {
                    p[row] = next[row] + beta * p[row];
                }
                rho[static_cast<std::size_t>(channel)] = nextRho;
            }
            preconditioned = std::move(nextPreconditioned);
        }
        if (activeChannels > 0) {
            const double worst = *std::max_element(
                relativeResidual.begin(), relativeResidual.end());
            throw std::runtime_error(
                "[Operator-coarse traces] Construction PCG did not converge; "
                "maximum relative residual=" + std::to_string(worst) + '.');
        }
        for (int channel = 0; channel < channels; ++channel) {
            result.krylovIterations +=
                channelIterations[static_cast<std::size_t>(channel)];
            result.krylovMaximumIterations = std::max(
                result.krylovMaximumIterations,
                channelIterations[static_cast<std::size_t>(channel)]);
            result.krylovMaximumRelativeResidual = std::max(
                result.krylovMaximumRelativeResidual,
                relativeResidual[static_cast<std::size_t>(channel)]);
        }
        return solution;
    };
    const auto recover = [&](const std::vector<double>& globalRhs,
                             const std::vector<double>& gamma,
                             int channels) {
        std::vector<double> full(static_cast<std::size_t>(
            descriptor.dofs * channels), 0.0);
        for (int channel = 0; channel < channels; ++channel) {
            for (int row = 0; row < interfaceDofs; ++row) {
                full[static_cast<std::size_t>(channel * descriptor.dofs
                    + partition.interfaceGlobalDofs[
                        static_cast<std::size_t>(row)])] = gamma[
                            static_cast<std::size_t>(channel * interfaceDofs + row)];
            }
        }
        std::vector<std::exception_ptr> errors(partition.domains.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(workers) if(workers > 1)
#endif
        for (int slotIndex = 0;
             slotIndex < static_cast<int>(partition.domains.size()); ++slotIndex) {
            try {
                const std::size_t slot = static_cast<std::size_t>(slotIndex);
                const auto& domain = partition.domains[slot];
                const int interiorRows = static_cast<int>(
                    domain.interiorGlobalDofs.size());
                std::vector<double> rhs(static_cast<std::size_t>(
                    interiorRows * channels), 0.0);
                for (int channel = 0; channel < channels; ++channel) {
                    for (int row = 0; row < interiorRows; ++row) {
                        double value = globalRhs[static_cast<std::size_t>(
                            channel * descriptor.dofs
                            + domain.interiorGlobalDofs[
                                static_cast<std::size_t>(row)])];
                        for (const auto& [localGamma, coefficient] :
                             k.interiorInterface[slot][
                                 static_cast<std::size_t>(row)]) {
                            const int globalGamma = partition.globalToInterface[
                                static_cast<std::size_t>(
                                    domain.interfaceGlobalDofs[
                                        static_cast<std::size_t>(localGamma)])];
                            value -= coefficient * gamma[static_cast<std::size_t>(
                                channel * interfaceDofs + globalGamma)];
                        }
                        rhs[static_cast<std::size_t>(
                            channel * interiorRows + row)] = value;
                    }
                }
                std::vector<double> interior;
                factors[slot]->solveMultiple(
                    rhs, channels, interior, factorThreads);
                for (int channel = 0; channel < channels; ++channel) {
                    for (int row = 0; row < interiorRows; ++row) {
                        full[static_cast<std::size_t>(channel * descriptor.dofs
                            + domain.interiorGlobalDofs[
                                static_cast<std::size_t>(row)])] = interior[
                                    static_cast<std::size_t>(
                                        channel * interiorRows + row)];
                    }
                }
            } catch (...) {
                errors[static_cast<std::size_t>(slotIndex)] =
                    std::current_exception();
            }
        }
        for (const auto& error : errors) {
            if (error) std::rethrow_exception(error);
        }
        return full;
    };

    const int channels = descriptor.sourceChannels;
    const std::vector<double> staticCondensed =
        condensedRhs(descriptor.input, channels);
    const std::vector<double> staticGamma =
        solveSchur(staticCondensed, channels);
    const std::vector<double> staticFull =
        recover(descriptor.input, staticGamma, channels);
    std::vector<double> momentRhs(staticFull.size(), 0.0);
    parallelForCoarse(static_cast<std::size_t>(channels),
        [&](std::size_t channelIndex) {
            std::vector<double> state(
                staticFull.begin() + static_cast<std::ptrdiff_t>(
                    channelIndex * descriptor.dofs),
                staticFull.begin() + static_cast<std::ptrdiff_t>(
                    (channelIndex + 1) * descriptor.dofs));
            std::vector<double> image = descriptor.capacity.multiply(state);
            for (int row = 0; row < descriptor.dofs; ++row) {
                momentRhs[static_cast<std::size_t>(
                    channelIndex * descriptor.dofs + row)] =
                    -image[static_cast<std::size_t>(row)];
            }
        });
    const std::vector<double> dynamicCondensed =
        condensedRhs(momentRhs, channels);
    const std::vector<double> dynamicGamma =
        solveSchur(dynamicCondensed, channels);

    result.traces.reserve(static_cast<std::size_t>(2 * channels));
    for (int channel = 0; channel < channels; ++channel) {
        std::vector<double> staticTrace(
            static_cast<std::size_t>(descriptor.dofs), 0.0);
        std::vector<double> dynamicTrace(
            static_cast<std::size_t>(descriptor.dofs), 0.0);
        for (int gamma = 0; gamma < interfaceDofs; ++gamma) {
            const int global = partition.interfaceGlobalDofs[
                static_cast<std::size_t>(gamma)];
            staticTrace[static_cast<std::size_t>(global)] = staticGamma[
                static_cast<std::size_t>(channel * interfaceDofs + gamma)];
            dynamicTrace[static_cast<std::size_t>(global)] = dynamicGamma[
                static_cast<std::size_t>(channel * interfaceDofs + gamma)];
        }
        result.traces.push_back(std::move(staticTrace));
        result.traces.push_back(std::move(dynamicTrace));
    }
    result.setupSeconds = elapsed(start);
    return result;
}

std::vector<double> operatorCoarseTraceBasis(
    const Mesh& mesh,
    const ddm_schur::DomainBlocks& domain,
    const std::vector<std::vector<double>>& coarseTraces,
    double tolerance,
    int maximumRank)
{
    const int rows = static_cast<int>(domain.interfaceGlobalDofs.size());
    std::vector<double> candidates;
    const auto append = [&](std::vector<double> column) {
        candidates.insert(candidates.end(), column.begin(), column.end());
    };
    for (const auto& globalTrace : coarseTraces) {
        std::vector<double> local(static_cast<std::size_t>(rows), 0.0);
        for (int row = 0; row < rows; ++row) {
            local[static_cast<std::size_t>(row)] = globalTrace[
                static_cast<std::size_t>(domain.interfaceGlobalDofs[
                    static_cast<std::size_t>(row)])];
        }
        append(std::move(local));
    }

    const auto quantize = [](double value) {
        constexpr double scale = 1.0e12;
        return std::round(value * scale) / scale;
    };
    for (const auto& [neighbor, globalDofs] :
         domain.interfaceGlobalDofsByNeighbor) {
        (void)neighbor;
        std::vector<int> localRows;
        localRows.reserve(globalDofs.size());
        std::array<double, 3> minimum{
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity()};
        std::array<double, 3> maximum{
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity()};
        for (int global : globalDofs) {
            const auto found = std::lower_bound(
                domain.interfaceGlobalDofs.begin(),
                domain.interfaceGlobalDofs.end(), global);
            if (found == domain.interfaceGlobalDofs.end() || *found != global) {
                throw std::runtime_error(
                    "[Operator-coarse traces] Port/local ordering mismatch.");
            }
            const int local = static_cast<int>(
                found - domain.interfaceGlobalDofs.begin());
            localRows.push_back(local);
            const Vec3& point = mesh.nodes[static_cast<std::size_t>(global)].p;
            const std::array<double, 3> coordinates{point.x, point.y, point.z};
            for (int axis = 0; axis < 3; ++axis) {
                minimum[static_cast<std::size_t>(axis)] = std::min(
                    minimum[static_cast<std::size_t>(axis)],
                    coordinates[static_cast<std::size_t>(axis)]);
                maximum[static_cast<std::size_t>(axis)] = std::max(
                    maximum[static_cast<std::size_t>(axis)],
                    coordinates[static_cast<std::size_t>(axis)]);
            }
        }
        std::vector<double> indicator(static_cast<std::size_t>(rows), 0.0);
        for (int row : localRows) indicator[static_cast<std::size_t>(row)] = 1.0;
        append(std::move(indicator));

        std::array<int, 3> axes{0, 1, 2};
        std::sort(axes.begin(), axes.end(), [&](int left, int right) {
            return maximum[static_cast<std::size_t>(left)]
                    - minimum[static_cast<std::size_t>(left)]
                > maximum[static_cast<std::size_t>(right)]
                    - minimum[static_cast<std::size_t>(right)];
        });
        std::vector<std::vector<double>> coordinates;
        for (int axis : axes) {
            const double range = maximum[static_cast<std::size_t>(axis)]
                - minimum[static_cast<std::size_t>(axis)];
            const double coordinateScale = std::max({
                1.0, std::abs(minimum[static_cast<std::size_t>(axis)]),
                std::abs(maximum[static_cast<std::size_t>(axis)])});
            if (!(range > 1.0e-12 * coordinateScale)) continue;
            const double center = 0.5 * (
                minimum[static_cast<std::size_t>(axis)]
                + maximum[static_cast<std::size_t>(axis)]);
            std::vector<double> mode(static_cast<std::size_t>(rows), 0.0);
            for (int row : localRows) {
                const Vec3& point = mesh.nodes[static_cast<std::size_t>(
                    domain.interfaceGlobalDofs[static_cast<std::size_t>(row)])].p;
                const std::array<double, 3> values{point.x, point.y, point.z};
                mode[static_cast<std::size_t>(row)] = quantize(
                    2.0 * (values[static_cast<std::size_t>(axis)] - center) / range);
            }
            coordinates.push_back(std::move(mode));
            if (coordinates.size() == 2) break;
        }
        if (coordinates.empty()) continue;
        const std::vector<double>& u = coordinates[0];
        const std::vector<double>* v = coordinates.size() > 1
            ? &coordinates[1] : nullptr;
        append(u);
        if (v != nullptr) append(*v);
        const auto polynomial = [&](int up, int vp) {
            std::vector<double> mode(static_cast<std::size_t>(rows), 0.0);
            for (int row : localRows) {
                const double a = u[static_cast<std::size_t>(row)];
                const double b = v != nullptr
                    ? (*v)[static_cast<std::size_t>(row)] : 1.0;
                mode[static_cast<std::size_t>(row)] = quantize(
                    std::pow(a, up) * std::pow(b, vp));
            }
            append(std::move(mode));
        };
        polynomial(2, 0);
        polynomial(3, 0);
        polynomial(4, 0);
        if (v != nullptr) {
            polynomial(0, 2);
            polynomial(1, 1);
            polynomial(0, 3);
            polynomial(2, 1);
            polynomial(1, 2);
            polynomial(0, 4);
        }
    }
    const int columns = rows > 0
        ? static_cast<int>(candidates.size()) / rows : 0;
    return orthonormalTraceBasis(
        rows, std::move(candidates), columns, tolerance, maximumRank);
}

void couplingTimesBasis(
    const std::vector<std::vector<std::pair<int, double>>>& rows,
    const std::vector<double>& gammaBasis,
    int gammaRows,
    int gammaRank,
    double scale,
    std::vector<double>& columns)
{
    const int interiorRows = static_cast<int>(rows.size());
    const std::size_t begin = columns.size();
    columns.resize(begin + static_cast<std::size_t>(interiorRows)
        * static_cast<std::size_t>(gammaRank), 0.0);
    parallelFor(static_cast<std::size_t>(gammaRank), [&](std::size_t modeIndex) {
        const int mode = static_cast<int>(modeIndex);
        double* output = columns.data() + begin
            + static_cast<std::size_t>(mode * interiorRows);
        for (int row = 0; row < interiorRows; ++row) {
            for (const auto& entry : rows[static_cast<std::size_t>(row)]) {
                output[static_cast<std::size_t>(row)] += scale * entry.second
                    * gammaBasis[static_cast<std::size_t>(mode * gammaRows + entry.first)];
            }
        }
    });
}

struct LocalModel {
    int domainId = -1;
    int interiorDofs = 0;
    int gammaDofs = 0;
    int physicalChannels = 0;
    int excitationRank = 0;
    int initialBlockRank = 0;
    int rank = 0;
    int deflated = 0;
    int templateId = -1;
    bool templateReused = false;
    std::uint64_t templateFingerprint = 0;
    std::vector<int> interiorGlobal;
    std::vector<int> gammaIndices;
    std::vector<double> basis;
    std::vector<ArnoldiHistoryRow> history;
    BlockArnoldiTiming arnoldiTiming;
    std::vector<double> cii;
    std::vector<double> kii;
    std::vector<double> ciGamma;
    std::vector<double> kiGamma;
    std::vector<double> cGammaI;
    std::vector<double> kGammaI;
    std::vector<double> reducedInput;
    std::vector<double> reducedBoundary;
    std::vector<double> referenceInterior;
};

constexpr std::uint64_t localDynamicCacheMagic =
    UINT64_C(0x4c444d4f52434143); // "LDMORCAC"
constexpr int localDynamicCacheVersion = 4;
constexpr int constructionTraceLocalDynamicCacheVersion = 3;
constexpr int secondMomentLocalDynamicCacheVersion = 2;
constexpr int legacyLocalDynamicCacheVersion = 1;
constexpr std::uint64_t localDynamicReferenceCacheMagic =
    UINT64_C(0x4c444d4f52524546); // "LDMORREF"
constexpr int localDynamicReferenceCacheVersion = 1;

int constructionTraceModeCode(const std::string& mode)
{
    if (mode == "global-fom") return 0;
    if (mode == "operator-coarse") return 1;
    throw std::runtime_error(
        "[Local dynamic cache] Unsupported construction trace mode.");
}

template <typename T>
void writeCacheScalar(std::ofstream& output, const T& value)
{
    static_assert(std::is_trivially_copyable_v<T>);
    output.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template <typename T>
T readCacheScalar(std::ifstream& input)
{
    static_assert(std::is_trivially_copyable_v<T>);
    T value{};
    input.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!input) {
        throw std::runtime_error(
            "[Local dynamic cache] Truncated cache header.");
    }
    return value;
}

template <typename T>
void writeCacheVector(std::ofstream& output,
                      const std::vector<T>& values)
{
    static_assert(std::is_trivially_copyable_v<T>);
    writeCacheScalar(output, static_cast<std::uint64_t>(values.size()));
    if (!values.empty()) {
        output.write(reinterpret_cast<const char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(T)));
    }
}

template <typename T>
std::vector<T> readCacheVector(std::ifstream& input)
{
    static_assert(std::is_trivially_copyable_v<T>);
    const std::uint64_t count = readCacheScalar<std::uint64_t>(input);
    if (count > (UINT64_C(1) << 38)) {
        throw std::runtime_error(
            "[Local dynamic cache] Invalid vector length.");
    }
    std::vector<T> values(static_cast<std::size_t>(count));
    if (!values.empty()) {
        input.read(reinterpret_cast<char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(T)));
        if (!input) {
            throw std::runtime_error(
                "[Local dynamic cache] Truncated vector payload.");
        }
    }
    return values;
}

std::filesystem::path localDynamicCachePath(
    std::filesystem::path path)
{
    if (!path.empty() && path.extension().empty()) {
        path /= "local_dynamic_interior_model.bin";
    }
    return path;
}

std::filesystem::path localDynamicReferenceCachePath(
    std::filesystem::path path)
{
    if (path.empty()) return {};
    if (path.extension().empty()) {
        return path / "local_dynamic_reference.bin";
    }
    return path.parent_path()
        / (path.stem().string() + "_reference.bin");
}

std::filesystem::path localDynamicDescriptorCachePath(
    std::filesystem::path path)
{
    if (path.empty()) return {};
    if (path.extension().empty()) {
        return path / "thermal_descriptor.bin";
    }
    return path.parent_path()
        / (path.stem().string() + "_thermal_descriptor.bin");
}

struct LocalDynamicReferenceCache {
    std::vector<double> reference;
    std::vector<double> boundaryOffset;
};

void saveLocalDynamicReference(
    const std::filesystem::path& path,
    const ThermalDescriptorSystem& descriptor,
    const LocalDynamicReferenceCache& cache)
{
    if (path.empty()) return;
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path, std::ios::binary);
    if (!output) {
        throw std::runtime_error(
            "[Local dynamic reference cache] Cannot create cache file.");
    }
    writeCacheScalar(output, localDynamicReferenceCacheMagic);
    writeCacheScalar(output, localDynamicReferenceCacheVersion);
    writeCacheScalar(output, descriptor.dofs);
    writeCacheScalar(output, descriptor.fingerprints.mesh);
    writeCacheScalar(output, descriptor.fingerprints.conductivity);
    writeCacheScalar(output, descriptor.fingerprints.boundary);
    writeCacheVector(output, cache.reference);
    writeCacheVector(output, cache.boundaryOffset);
    if (!output) {
        throw std::runtime_error(
            "[Local dynamic reference cache] Failed while writing payload.");
    }
}

LocalDynamicReferenceCache loadLocalDynamicReference(
    const std::filesystem::path& path,
    const ThermalDescriptorSystem& descriptor)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error(
            "[Local dynamic reference cache] Cannot open cache file.");
    }
    const bool headerValid =
        readCacheScalar<std::uint64_t>(input)
            == localDynamicReferenceCacheMagic
        && readCacheScalar<int>(input)
            == localDynamicReferenceCacheVersion
        && readCacheScalar<int>(input) == descriptor.dofs
        && readCacheScalar<std::uint64_t>(input)
            == descriptor.fingerprints.mesh
        && readCacheScalar<std::uint64_t>(input)
            == descriptor.fingerprints.conductivity
        && readCacheScalar<std::uint64_t>(input)
            == descriptor.fingerprints.boundary;
    if (!headerValid) {
        throw std::runtime_error(
            "[Local dynamic reference cache] Fingerprint mismatch; reuse refused.");
    }
    LocalDynamicReferenceCache cache;
    cache.reference = readCacheVector<double>(input);
    cache.boundaryOffset = readCacheVector<double>(input);
    if (cache.reference.size()
            != static_cast<std::size_t>(descriptor.dofs)
        || cache.boundaryOffset.size() != cache.reference.size()) {
        throw std::runtime_error(
            "[Local dynamic reference cache] Dimension mismatch; reuse refused.");
    }
    return cache;
}

void writeLocalDynamicCacheHeader(
    std::ofstream& output,
    const ThermalDescriptorSystem& descriptor,
    const Options& options)
{
    writeCacheScalar(output, localDynamicCacheMagic);
    writeCacheScalar(output, localDynamicCacheVersion);
    writeCacheScalar(output, descriptor.dofs);
    writeCacheScalar(output, descriptor.sourceChannels);
    writeCacheScalar(output, descriptor.fingerprints.mesh);
    writeCacheScalar(output, descriptor.fingerprints.capacity);
    writeCacheScalar(output, descriptor.fingerprints.conductivity);
    writeCacheScalar(output, descriptor.fingerprints.input);
    writeCacheScalar(output, descriptor.fingerprints.boundary);
    writeCacheScalar(output, descriptor.fingerprints.sources);
    writeCacheScalar(output, options.moments);
    writeCacheScalar(output, options.expansionPoint);
    writeCacheScalar(output, options.rankTolerance);
    writeCacheScalar(output, options.secondMomentEnergy);
    writeCacheScalar(output, options.secondMomentMaximumColumns);
    writeCacheScalar(output,
        constructionTraceModeCode(options.constructionTraceMode));
    writeCacheScalar(output, options.interfaceExcitationRank);
    writeCacheScalar(output, options.reuseIdenticalSubdomains);
}

void validateLocalDynamicCacheHeader(
    std::ifstream& input,
    const ThermalDescriptorSystem& descriptor,
    const Options& options)
{
    const std::uint64_t magic = readCacheScalar<std::uint64_t>(input);
    const int version = readCacheScalar<int>(input);
    const int dofs = readCacheScalar<int>(input);
    const int sourceChannels = readCacheScalar<int>(input);
    const std::uint64_t mesh = readCacheScalar<std::uint64_t>(input);
    const std::uint64_t capacity = readCacheScalar<std::uint64_t>(input);
    const std::uint64_t conductivity = readCacheScalar<std::uint64_t>(input);
    const std::uint64_t inputFingerprint = readCacheScalar<std::uint64_t>(input);
    const std::uint64_t boundary = readCacheScalar<std::uint64_t>(input);
    const std::uint64_t sources = readCacheScalar<std::uint64_t>(input);
    const int moments = readCacheScalar<int>(input);
    const double expansionPoint = readCacheScalar<double>(input);
    const double rankTolerance = readCacheScalar<double>(input);
    double secondMomentEnergy = 1.0;
    int secondMomentMaximumColumns = 0;
    if (version == localDynamicCacheVersion
        || version == constructionTraceLocalDynamicCacheVersion
        || version == secondMomentLocalDynamicCacheVersion) {
        secondMomentEnergy = readCacheScalar<double>(input);
        secondMomentMaximumColumns = readCacheScalar<int>(input);
    }
    int constructionTraceMode = 0;
    if (version == localDynamicCacheVersion
        || version == constructionTraceLocalDynamicCacheVersion) {
        constructionTraceMode = readCacheScalar<int>(input);
    }
    const int interfaceExcitationRank = readCacheScalar<int>(input);
    const bool reuseIdenticalSubdomains = readCacheScalar<bool>(input);
    const bool valid = magic == localDynamicCacheMagic
        && (version == localDynamicCacheVersion
            || version == constructionTraceLocalDynamicCacheVersion
            || version == secondMomentLocalDynamicCacheVersion
            || version == legacyLocalDynamicCacheVersion)
        && dofs == descriptor.dofs
        && sourceChannels == descriptor.sourceChannels
        && mesh == descriptor.fingerprints.mesh
        && capacity == descriptor.fingerprints.capacity
        && conductivity == descriptor.fingerprints.conductivity
        && inputFingerprint == descriptor.fingerprints.input
        && boundary == descriptor.fingerprints.boundary
        && sources == descriptor.fingerprints.sources
        && moments == options.moments
        && expansionPoint == options.expansionPoint
        && rankTolerance == options.rankTolerance
        && secondMomentEnergy == options.secondMomentEnergy
        && secondMomentMaximumColumns == options.secondMomentMaximumColumns
        && constructionTraceMode
            == constructionTraceModeCode(options.constructionTraceMode)
        && (version != constructionTraceLocalDynamicCacheVersion
            || constructionTraceMode == 0)
        && interfaceExcitationRank == options.interfaceExcitationRank
        && reuseIdenticalSubdomains == options.reuseIdenticalSubdomains;
    if (!valid) {
        throw std::runtime_error(
            "[Local dynamic cache] Fingerprint or basis parameter mismatch; reuse refused.");
    }
}

void saveLocalDynamicModels(
    const std::filesystem::path& path,
    const ThermalDescriptorSystem& descriptor,
    const Options& options,
    const std::vector<LocalModel>& models)
{
    if (path.empty()) return;
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::ofstream output(path, std::ios::binary);
    if (!output) {
        throw std::runtime_error(
            "[Local dynamic cache] Cannot create cache file.");
    }
    writeLocalDynamicCacheHeader(output, descriptor, options);
    writeCacheScalar(output, static_cast<std::uint64_t>(models.size()));
    for (const LocalModel& model : models) {
        writeCacheScalar(output, model.domainId);
        writeCacheScalar(output, model.interiorDofs);
        writeCacheScalar(output, model.gammaDofs);
        writeCacheScalar(output, model.physicalChannels);
        writeCacheScalar(output, model.excitationRank);
        writeCacheScalar(output, model.initialBlockRank);
        writeCacheScalar(output, model.rank);
        writeCacheScalar(output, model.deflated);
        writeCacheScalar(output, model.templateId);
        writeCacheScalar(output, model.templateReused);
        writeCacheScalar(output, model.templateFingerprint);
        writeCacheVector(output, model.interiorGlobal);
        writeCacheVector(output, model.gammaIndices);
        writeCacheVector(output, model.basis);
        writeCacheVector(output, model.history);
        writeCacheScalar(output, model.arnoldiTiming);
        writeCacheVector(output, model.cii);
        writeCacheVector(output, model.kii);
        writeCacheVector(output, model.ciGamma);
        writeCacheVector(output, model.kiGamma);
        writeCacheVector(output, model.cGammaI);
        writeCacheVector(output, model.kGammaI);
        writeCacheVector(output, model.reducedInput);
        writeCacheVector(output, model.reducedBoundary);
        writeCacheVector(output, model.referenceInterior);
    }
    if (!output) {
        throw std::runtime_error(
            "[Local dynamic cache] Failed while writing cache payload.");
    }
}

std::vector<LocalModel> loadLocalDynamicModels(
    const std::filesystem::path& path,
    const ThermalDescriptorSystem& descriptor,
    const Options& options,
    const ddm_schur::InterfacePartition& partition)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error(
            "[Local dynamic cache] Cannot open cache file.");
    }
    validateLocalDynamicCacheHeader(input, descriptor, options);
    const std::uint64_t count = readCacheScalar<std::uint64_t>(input);
    if (count != partition.domains.size()) {
        throw std::runtime_error(
            "[Local dynamic cache] Subdomain count mismatch; reuse refused.");
    }
    std::vector<LocalModel> models(static_cast<std::size_t>(count));
    for (std::size_t slot = 0; slot < models.size(); ++slot) {
        LocalModel& model = models[slot];
        model.domainId = readCacheScalar<int>(input);
        model.interiorDofs = readCacheScalar<int>(input);
        model.gammaDofs = readCacheScalar<int>(input);
        model.physicalChannels = readCacheScalar<int>(input);
        model.excitationRank = readCacheScalar<int>(input);
        model.initialBlockRank = readCacheScalar<int>(input);
        model.rank = readCacheScalar<int>(input);
        model.deflated = readCacheScalar<int>(input);
        model.templateId = readCacheScalar<int>(input);
        model.templateReused = readCacheScalar<bool>(input);
        model.templateFingerprint = readCacheScalar<std::uint64_t>(input);
        model.interiorGlobal = readCacheVector<int>(input);
        model.gammaIndices = readCacheVector<int>(input);
        model.basis = readCacheVector<double>(input);
        model.history = readCacheVector<ArnoldiHistoryRow>(input);
        model.arnoldiTiming = readCacheScalar<BlockArnoldiTiming>(input);
        model.cii = readCacheVector<double>(input);
        model.kii = readCacheVector<double>(input);
        model.ciGamma = readCacheVector<double>(input);
        model.kiGamma = readCacheVector<double>(input);
        model.cGammaI = readCacheVector<double>(input);
        model.kGammaI = readCacheVector<double>(input);
        model.reducedInput = readCacheVector<double>(input);
        model.reducedBoundary = readCacheVector<double>(input);
        model.referenceInterior = readCacheVector<double>(input);

        const auto& domain = partition.domains[slot];
        std::vector<int> expectedGamma;
        expectedGamma.reserve(domain.interfaceGlobalDofs.size());
        for (int global : domain.interfaceGlobalDofs) {
            expectedGamma.push_back(partition.globalToInterface[
                static_cast<std::size_t>(global)]);
        }
        const bool dimensionsValid = model.domainId == domain.domainId
            && model.interiorDofs
                == static_cast<int>(domain.interiorGlobalDofs.size())
            && model.gammaDofs
                == static_cast<int>(domain.interfaceGlobalDofs.size())
            && model.interiorGlobal == domain.interiorGlobalDofs
            && model.gammaIndices == expectedGamma
            && model.rank > 0
            && model.basis.size() == static_cast<std::size_t>(
                model.interiorDofs) * model.rank
            && model.cii.size() == static_cast<std::size_t>(
                model.rank) * model.rank
            && model.kii.size() == model.cii.size()
            && model.ciGamma.size() == static_cast<std::size_t>(
                model.rank) * model.gammaDofs
            && model.kiGamma.size() == model.ciGamma.size()
            && model.cGammaI.size() == static_cast<std::size_t>(
                model.gammaDofs) * model.rank
            && model.kGammaI.size() == model.cGammaI.size();
        if (!dimensionsValid) {
            throw std::runtime_error(
                "[Local dynamic cache] Partition/order/dimension mismatch; reuse refused.");
        }
    }
    return models;
}

std::vector<double> projectInput(const ThermalDescriptorSystem& descriptor,
                                 const LocalModel& local)
{
    std::vector<double> result(static_cast<std::size_t>(local.rank)
        * descriptor.sourceChannels, 0.0);
    for (int channel = 0; channel < descriptor.sourceChannels; ++channel) {
        for (int mode = 0; mode < local.rank; ++mode) {
            long double value = 0.0L;
            for (int row = 0; row < local.interiorDofs; ++row) {
                value += static_cast<long double>(local.basis[static_cast<std::size_t>(
                    mode * local.interiorDofs + row)]) * descriptor.input[static_cast<std::size_t>(
                    channel * descriptor.dofs + local.interiorGlobal[static_cast<std::size_t>(row)])];
            }
            result[static_cast<std::size_t>(channel * local.rank + mode)] =
                static_cast<double>(value);
        }
    }
    return result;
}

void addInput(const ThermalDescriptorSystem& descriptor,
              const std::vector<double>& powers,
              std::vector<double>& rhs)
{
    for (int channel = 0; channel < descriptor.sourceChannels; ++channel) {
        const double coefficient = powers[static_cast<std::size_t>(channel)];
        const double* input = descriptor.input.data()
            + static_cast<std::size_t>(channel * descriptor.dofs);
        for (int row = 0; row < descriptor.dofs; ++row) {
            rhs[static_cast<std::size_t>(row)] += coefficient * input[row];
        }
    }
}

void addReducedInput(const LocalModel& local,
                     const std::vector<double>& powers,
                     std::vector<double>& rhs)
{
    for (std::size_t channel = 0; channel < powers.size(); ++channel) {
        for (int row = 0; row < local.rank; ++row) {
            rhs[static_cast<std::size_t>(row)] += powers[channel]
                * local.reducedInput[channel * static_cast<std::size_t>(local.rank)
                    + static_cast<std::size_t>(row)];
        }
    }
}

void denseMatvec(const std::vector<double>& matrix,
                 int rows,
                 int columns,
                 const std::vector<double>& x,
                 std::vector<double>& y,
                 double alpha = 1.0,
                 double beta = 0.0)
{
    if (y.size() != static_cast<std::size_t>(rows)) y.assign(rows, 0.0);
    for (int row = 0; row < rows; ++row) {
        long double value = 0.0L;
        for (int column = 0; column < columns; ++column) {
            value += static_cast<long double>(matrix[static_cast<std::size_t>(
                row * columns + column)]) * x[static_cast<std::size_t>(column)];
        }
        y[static_cast<std::size_t>(row)] = alpha * static_cast<double>(value)
            + beta * y[static_cast<std::size_t>(row)];
    }
}

struct NativeReducedHistoryRhs {
    std::vector<std::vector<double>> interior;
    std::vector<double> interfaceRhs;
};

NativeReducedHistoryRhs buildNativeReducedHistoryRhs(
    const ThermalDescriptorSystem& descriptor,
    const MatrixPartition& capacity,
    const ddm_schur::InterfacePartition& partition,
    const std::vector<LocalModel>& locals,
    const std::vector<std::vector<double>>& localStates,
    const std::vector<double>& interfaceState,
    const std::vector<double>& powers,
    const std::vector<double>& boundaryOffset,
    double timeStep)
{
    if (!(timeStep > 0.0) || locals.size() != localStates.size()
        || interfaceState.size() != partition.interfaceGlobalDofs.size()) {
        throw std::runtime_error(
            "[Local transient] Native reduced history state is inconsistent.");
    }
    NativeReducedHistoryRhs result;
    result.interior.resize(locals.size());
    const std::vector<double> interfaceHistory =
        capacity.interfaceBlock.multiply(interfaceState);
    result.interfaceRhs.resize(interfaceState.size(), 0.0);
    for (std::size_t gamma = 0; gamma < result.interfaceRhs.size(); ++gamma) {
        const int global = partition.interfaceGlobalDofs[gamma];
        double value = interfaceHistory[gamma] / timeStep
            + boundaryOffset[static_cast<std::size_t>(global)];
        for (int channel = 0; channel < descriptor.sourceChannels; ++channel) {
            value += powers[static_cast<std::size_t>(channel)]
                * descriptor.input[static_cast<std::size_t>(
                    channel * descriptor.dofs + global)];
        }
        result.interfaceRhs[gamma] = value;
    }
    for (std::size_t slot = 0; slot < locals.size(); ++slot) {
        const LocalModel& local = locals[slot];
        if (localStates[slot].size() != static_cast<std::size_t>(local.rank)) {
            throw std::runtime_error(
                "[Local transient] Native reduced interior state has the wrong rank.");
        }
        std::vector<double>& localRhs = result.interior[slot];
        denseMatvec(local.cii, local.rank, local.rank, localStates[slot],
            localRhs, 1.0 / timeStep, 0.0);
        std::vector<double> localGamma(
            static_cast<std::size_t>(local.gammaDofs), 0.0);
        for (int row = 0; row < local.gammaDofs; ++row) {
            localGamma[static_cast<std::size_t>(row)] = interfaceState[
                static_cast<std::size_t>(local.gammaIndices[
                    static_cast<std::size_t>(row)])];
        }
        denseMatvec(local.ciGamma, local.rank, local.gammaDofs, localGamma,
            localRhs, 1.0 / timeStep, 1.0);
        addReducedInput(local, powers, localRhs);
        for (int mode = 0; mode < local.rank; ++mode) {
            localRhs[static_cast<std::size_t>(mode)] +=
                local.reducedBoundary[static_cast<std::size_t>(mode)];
        }
        std::vector<double> interfaceContribution;
        denseMatvec(local.cGammaI, local.gammaDofs, local.rank,
            localStates[slot], interfaceContribution, 1.0 / timeStep, 0.0);
        for (int row = 0; row < local.gammaDofs; ++row) {
            result.interfaceRhs[static_cast<std::size_t>(local.gammaIndices[
                static_cast<std::size_t>(row)])] +=
                interfaceContribution[static_cast<std::size_t>(row)];
        }
    }
    return result;
}

struct DynamicLocal {
    const LocalModel* model = nullptr;
    std::vector<double> aii;
    std::vector<double> aiGamma;
    std::vector<double> aGammaI;
    local::DenseSymmetricFactor factor;
    std::vector<double> solvedCoupling;
};

local::Model makeDynamicReducedModel(
    int globalDofs,
    const ddm_schur::InterfacePartition& partition,
    const MatrixPartition& k,
    const MatrixPartition& c,
    const std::vector<LocalModel>& locals,
    double timeStep)
{
    local::Model model;
    model.globalDofs = globalDofs;
    model.interfaceDofs = static_cast<int>(partition.interfaceGlobalDofs.size());
    model.interfaceGlobalDofs = partition.interfaceGlobalDofs;
    const auto appendInterface = [&](const SparseMatrix& matrix, double scale) {
        matrix.forEachEntry([&](int row, int column, double value) {
            model.interfaceEntries.push_back({row, column, scale * value});
        });
    };
    appendInterface(k.interfaceBlock, 1.0);
    appendInterface(c.interfaceBlock, 1.0 / timeStep);
    model.subdomains.reserve(locals.size());
    for (std::size_t slot = 0; slot < locals.size(); ++slot) {
        const LocalModel& source = locals[slot];
        local::SubdomainModel reduced;
        reduced.subdomain = source.domainId;
        reduced.interiorDofs = source.interiorDofs;
        reduced.localInterfaceDofs = source.gammaDofs;
        reduced.rank = source.rank;
        reduced.numericalRank = source.rank;
        reduced.boundaryInterfaceFingerprint = source.templateFingerprint;
        reduced.interiorGlobalDofs = source.interiorGlobal;
        reduced.interfaceGlobalDofs = partition.domains[slot].interfaceGlobalDofs;
        reduced.interfaceIndices = source.gammaIndices;
        reduced.referenceInterior.assign(
            static_cast<std::size_t>(source.interiorDofs), 0.0);
        reduced.basis = source.basis;
        reduced.reducedInterior.resize(source.kii.size());
        for (std::size_t entry = 0; entry < reduced.reducedInterior.size(); ++entry) {
            reduced.reducedInterior[entry] = source.kii[entry]
                + source.cii[entry] / timeStep;
        }
        reduced.reducedInteriorInterface.resize(source.kiGamma.size());
        for (std::size_t entry = 0;
             entry < reduced.reducedInteriorInterface.size(); ++entry) {
            reduced.reducedInteriorInterface[entry] = source.kiGamma[entry]
                + source.ciGamma[entry] / timeStep;
        }
        reduced.reducedInterfaceInterior.resize(source.kGammaI.size());
        for (std::size_t entry = 0;
             entry < reduced.reducedInterfaceInterior.size(); ++entry) {
            reduced.reducedInterfaceInterior[entry] = source.kGammaI[entry]
                + source.cGammaI[entry] / timeStep;
        }
        reduced.interiorReferenceImage.assign(
            static_cast<std::size_t>(source.interiorDofs), 0.0);
        reduced.interfaceReferenceImage.assign(
            static_cast<std::size_t>(source.gammaDofs), 0.0);
        model.totalLocalRank += source.rank;
        model.subdomains.push_back(std::move(reduced));
    }
    model.modelBytes = local::estimateModelBytes(model);
    return model;
}

struct ReducedHistorySources {
    int channels = 0;
    std::vector<double> condensedColumns;
};

ReducedHistorySources buildReducedHistorySources(
    const local::Model& dynamicModel,
    const std::vector<LocalModel>& locals,
    double timeStep)
{
    if (!(timeStep > 0.0)
        || dynamicModel.subdomains.size() != locals.size()) {
        throw std::runtime_error(
            "[Optimal port] Invalid Local Block Arnoldi history source.");
    }
    ReducedHistorySources result;
    for (const LocalModel& localData : locals) {
        result.channels += localData.rank;
    }
    result.condensedColumns.assign(static_cast<std::size_t>(
        dynamicModel.interfaceDofs * result.channels), 0.0);
    int channel = 0;
    for (std::size_t slot = 0; slot < locals.size(); ++slot) {
        const LocalModel& localData = locals[slot];
        const local::SubdomainModel& dynamic =
            dynamicModel.subdomains[slot];
        const local::DenseSymmetricFactor factor =
            local::factorDenseSymmetric(
                dynamic.reducedInterior, dynamic.rank);
        for (int previousMode = 0;
             previousMode < localData.rank; ++previousMode, ++channel) {
            std::vector<double> reducedInterior(
                static_cast<std::size_t>(localData.rank), 0.0);
            for (int row = 0; row < localData.rank; ++row) {
                reducedInterior[static_cast<std::size_t>(row)] =
                    localData.cii[static_cast<std::size_t>(
                        row * localData.rank + previousMode)] / timeStep;
            }
            local::solveDenseSymmetric(factor, reducedInterior);
            for (int gamma = 0; gamma < localData.gammaDofs; ++gamma) {
                double value = localData.cGammaI[static_cast<std::size_t>(
                    gamma * localData.rank + previousMode)] / timeStep;
                for (int mode = 0; mode < localData.rank; ++mode) {
                    value -= dynamic.reducedInterfaceInterior[
                        static_cast<std::size_t>(
                            gamma * localData.rank + mode)]
                        * reducedInterior[static_cast<std::size_t>(mode)];
                }
                const int globalGamma =
                    localData.gammaIndices[static_cast<std::size_t>(gamma)];
                result.condensedColumns[static_cast<std::size_t>(
                    channel * dynamicModel.interfaceDofs
                    + globalGamma)] += value;
            }
        }
    }
    return result;
}

void hashZeroDoubles(std::uint64_t& hash, std::size_t count)
{
    // FNV-1a with a zero byte is multiplication by the FNV prime.  Apply the
    // repeated multiplication by exponentiation so a sparse, logically
    // dense history operator can be fingerprinted without allocating or
    // scanning its zero entries.
    std::uint64_t power = UINT64_C(1099511628211);
    std::uint64_t multiplier = 1;
    std::size_t exponent = count * sizeof(double);
    while (exponent != 0) {
        if ((exponent & 1U) != 0) multiplier *= power;
        power *= power;
        exponent >>= 1U;
    }
    hash *= multiplier;
}

std::uint64_t fingerprintReducedHistorySourcesStreaming(
    const local::Model& dynamicModel,
    const std::vector<LocalModel>& locals,
    double timeStep)
{
    if (!(timeStep > 0.0)
        || dynamicModel.subdomains.size() != locals.size()) {
        throw std::runtime_error(
            "[Global randomized] Invalid streaming history fingerprint.");
    }
    int channels = 0;
    for (const LocalModel& localData : locals) {
        channels += localData.rank;
    }
    const std::size_t logicalSize =
        static_cast<std::size_t>(dynamicModel.interfaceDofs)
        * static_cast<std::size_t>(channels);
    std::uint64_t fingerprint = UINT64_C(1469598103934665603);
    hashValue(fingerprint, logicalSize);
    for (std::size_t slot = 0; slot < locals.size(); ++slot) {
        const LocalModel& localData = locals[slot];
        const local::SubdomainModel& dynamic =
            dynamicModel.subdomains[slot];
        const local::DenseSymmetricFactor factor =
            local::factorDenseSymmetric(
                dynamic.reducedInterior, dynamic.rank);
        for (int previousMode = 0;
             previousMode < localData.rank; ++previousMode) {
            std::vector<double> reducedInterior(
                static_cast<std::size_t>(localData.rank), 0.0);
            for (int row = 0; row < localData.rank; ++row) {
                reducedInterior[static_cast<std::size_t>(row)] =
                    localData.cii[static_cast<std::size_t>(
                        row * localData.rank + previousMode)]
                    / timeStep;
            }
            local::solveDenseSymmetric(factor, reducedInterior);
            std::map<int, double> sparseColumn;
            for (int gamma = 0;
                 gamma < localData.gammaDofs; ++gamma) {
                double value = localData.cGammaI[
                    static_cast<std::size_t>(
                        gamma * localData.rank + previousMode)]
                    / timeStep;
                for (int mode = 0; mode < localData.rank; ++mode) {
                    value -= dynamic.reducedInterfaceInterior[
                        static_cast<std::size_t>(
                            gamma * localData.rank + mode)]
                        * reducedInterior[
                            static_cast<std::size_t>(mode)];
                }
                sparseColumn[
                    localData.gammaIndices[
                        static_cast<std::size_t>(gamma)]] += value;
            }
            std::size_t cursor = 0;
            for (const auto& [row, value] : sparseColumn) {
                if (row < 0 || row >= dynamicModel.interfaceDofs
                    || static_cast<std::size_t>(row) < cursor) {
                    throw std::runtime_error(
                        "[Global randomized] Streaming history "
                        "fingerprint has an invalid interface row.");
                }
                hashZeroDoubles(
                    fingerprint,
                    static_cast<std::size_t>(row) - cursor);
                hashValue(fingerprint, value);
                cursor = static_cast<std::size_t>(row) + 1;
            }
            hashZeroDoubles(
                fingerprint,
                static_cast<std::size_t>(
                    dynamicModel.interfaceDofs) - cursor);
        }
    }
    return fingerprint;
}

struct Accuracy {
    double spaceTimeRelativeL2 = 0.0;
    double maximumAbsolute = 0.0;
    double maximumTemperatureError = 0.0;
    double maximumFullResidual = 0.0;
    double maximumFullResidualBeforeGate = 0.0;
    bool residualGateAllPassed = true;
    int residualToleranceViolationSteps = 0;
    int residualFallbackSteps = 0;
    double residualFallbackSolveSeconds = 0.0;
    double maximumReducedResidual = 0.0;
    double maximumJump = 0.0;
    double maximumFluxImbalance = 0.0;
    double maximumFluxRelativeL2 = 0.0;
    double localCoreSeconds = 0.0;
    double interfaceSolveSeconds = 0.0;
    double recoverySeconds = 0.0;
    double proxySolveSeconds = 0.0;
    double coarseSolveSeconds = 0.0;
    double portForwardSolveSeconds = 0.0;
    double portCoreSolveSeconds = 0.0;
    double portBackSubstitutionSeconds = 0.0;
    double interfaceOperatorSeconds = 0.0;
    double interfacePreconditionerSeconds = 0.0;
    double interfaceOrthogonalizationSeconds = 0.0;
    double interfaceVectorUpdateSeconds = 0.0;
    double interfacePredictorSeconds = 0.0;
    int interfacePredictorAppliedSteps = 0;
    int interfacePredictorAcceptedSteps = 0;
    double nativeReducedRhsSeconds = 0.0;
    int nativeReducedHistorySteps = 0;
    double stepRhsSeconds = 0.0;
    double fullResidualSeconds = 0.0;
    int adaptiveInterfaceRetrySteps = 0;
    int adaptiveInterfaceRetryIterations = 0;
    double adaptiveInterfaceRetrySeconds = 0.0;
    int interfaceIterationsTotal = 0;
    int interfaceIterationsMaximum = 0;
    int interfaceMatvecs = 0;
    int interfaceKrylovFallbackSteps = 0;
    std::string interfaceKrylovActual = "not_run";
    double maximumInterfaceResidual = 0.0;
    double maximumPortProjectionError = 0.0;
    double maximumPortReducedResidual = 0.0;
    double fomFactorSeconds = 0.0;
    double fomSolveSeconds = 0.0;
    std::size_t factorBytes = 0;
    std::size_t dynamicSchurFactorBytes = 0;
    std::size_t fomFactorBytes = 0;
    double correctedRelativeL2 = 0.0;
    double correctedMaximumAbsolute = 0.0;
    double correctedMaximumResidual = 0.0;
    double correctedSolveSeconds = 0.0;
};

std::vector<double> reconstruct(const std::vector<double>& reference,
                                const ddm_schur::InterfacePartition& partition,
                                const std::vector<LocalModel>& locals,
                                const std::vector<std::vector<double>>& states,
                                const std::vector<double>& gamma)
{
    std::vector<double> result = reference;
    for (std::size_t slot = 0; slot < locals.size(); ++slot) {
        const LocalModel& local = locals[slot];
        for (int row = 0; row < local.interiorDofs; ++row) {
            long double value = 0.0L;
            for (int mode = 0; mode < local.rank; ++mode) {
                value += static_cast<long double>(local.basis[static_cast<std::size_t>(
                    mode * local.interiorDofs + row)])
                    * states[slot][static_cast<std::size_t>(mode)];
            }
            result[static_cast<std::size_t>(local.interiorGlobal[static_cast<std::size_t>(row)])]
                += static_cast<double>(value);
        }
    }
    for (std::size_t row = 0; row < partition.interfaceGlobalDofs.size(); ++row) {
        result[static_cast<std::size_t>(partition.interfaceGlobalDofs[row])] += gamma[row];
    }
    return result;
}

std::vector<double> projectInitial(const LocalModel& local,
                                   const SparseMatrix& cii,
                                   const std::vector<double>& theta)
{
    std::vector<double> interior(static_cast<std::size_t>(local.interiorDofs), 0.0);
    for (int row = 0; row < local.interiorDofs; ++row) {
        interior[static_cast<std::size_t>(row)] = theta[static_cast<std::size_t>(
            local.interiorGlobal[static_cast<std::size_t>(row)])];
    }
    const std::vector<double> weighted = cii.multiply(interior);
    std::vector<double> rhs(static_cast<std::size_t>(local.rank), 0.0);
    for (int mode = 0; mode < local.rank; ++mode) {
        rhs[static_cast<std::size_t>(mode)] = dot(
            local.basis.data() + static_cast<std::size_t>(mode * local.interiorDofs),
            weighted.data(), local.interiorDofs);
    }
    local::DenseSymmetricFactor factor = local::factorDenseSymmetric(local.cii, local.rank);
    local::solveDenseSymmetric(factor, rhs);
    return rhs;
}

std::vector<double> initialTemperature(const ThermalDescriptorSystem& descriptor,
                                       const Options& options,
                                       const std::vector<double>& reference,
                                       SubdomainDirectSolver* steadyFactor)
{
    if (options.initialMode == "steady") {
        if (steadyFactor == nullptr) {
            throw std::runtime_error(
                "[Local dynamic] Steady initial mode requires a steady factor.");
        }
        std::vector<double> rhs = descriptorInputRhs(
            descriptor, descriptor.nominalPowersW);
        std::vector<double> result;
        steadyFactor->solve(rhs, result);
        return result;
    }
    if (options.initialMode == "uniform" || std::isfinite(options.initialTemperature)) {
        const double value = std::isfinite(options.initialTemperature)
            ? options.initialTemperature : 300.0;
        return std::vector<double>(static_cast<std::size_t>(descriptor.dofs), value);
    }
    return reference;
}

struct EnrichmentResult {
    LocalPortSnapshotFamilies snapshots;
    int addedRank = 0;
    double factorizationSeconds = 0.0;
    double solveSeconds = 0.0;
    double orthogonalizationSeconds = 0.0;
    double totalSeconds = 0.0;
};

void writePortRankDiagnostics(const LocalPortModel& model,
                              const std::filesystem::path& outputDirectory)
{
    std::ofstream output(outputDirectory / "local_port_rank_by_interface.csv");
    if (model.basisMethod == "port-pod") {
        output << "interface_id,left_subdomain,right_subdomain,full_interface_rows,"
            "candidate_columns,accepted_columns,selected_rank,retained_energy,"
            "orthogonality_error,template_id,template_reused,fingerprint\n";
    } else {
        output << "method,interface_id,left_subdomain,right_subdomain,"
            "full_interface_rows,candidate_columns,accepted_columns,total_port_rank,"
            "mandatory_rank,requested_transfer_rank,converged_transfer_rank,"
            "retained_energy,orthogonality_error,"
            "template_id,template_reused,fingerprint\n";
    }
    output << std::setprecision(17);
    for (const LocalPortBasis& port : model.ports) {
        if (model.basisMethod != "port-pod") {
            output << model.basisMethod << ',';
        }
        output << port.interfaceId << ',' << port.leftSubdomain << ','
            << port.rightSubdomain << ',' << port.rows << ','
            << port.candidateColumns << ',' << port.acceptedColumns << ','
            << port.rank << ',';
        if (model.basisMethod != "port-pod") {
            output << port.mandatoryModes << ','
                << port.requestedTransferRank << ','
                << port.spectralModes << ',';
        }
        output << port.retainedEnergy << ',' << port.orthogonalityError << ','
            << port.templateId << ',' << (port.templateReused ? 1 : 0)
            << ',' << port.fingerprint << '\n';
    }
}

void writeOptimalPortDiagnostics(
    const LocalPortModel& model,
    const OptimalPortBuildResult& build,
    const std::filesystem::path& outputDirectory)
{
    std::ofstream rank(
        outputDirectory / "optimal_port_rank_by_interface.csv");
    rank << "interface_id,left_subdomain,right_subdomain,target_rows,"
        "source_rows,mandatory_rank,requested_transfer_rank,"
        "converged_transfer_rank,total_port_rank,"
        "trace_source_rows,input_source_rows,boundary_source_rows,"
        "history_source_rows,"
        "transfer_indicator,orthogonality_error,target_fingerprint,"
        "trace_source_fingerprint,input_source_fingerprint,"
        "boundary_source_fingerprint,history_source_fingerprint,"
        "basis_fingerprint\n"
        << std::setprecision(17);
    for (const LocalPortBasis& port : model.ports) {
        rank << port.interfaceId << ',' << port.leftSubdomain << ','
            << port.rightSubdomain << ',' << port.rows << ','
            << (port.traceSourceRows + port.inputSourceRows
                + port.boundarySourceRows + port.historySourceRows) << ','
            << port.mandatoryModes << ','
            << port.requestedTransferRank << ','
            << port.spectralModes << ',' << port.rank << ','
            << port.traceSourceRows << ',' << port.inputSourceRows << ','
            << port.boundarySourceRows << ',' << port.historySourceRows << ','
            << port.transferIndicator << ',' << port.orthogonalityError
            << ',' << port.targetFingerprint << ','
            << port.traceSourceFingerprint << ','
            << port.inputSourceFingerprint << ','
            << port.boundarySourceFingerprint << ','
            << port.historySourceFingerprint << ','
            << port.fingerprint << '\n';
    }

    std::ofstream decay(
        outputDirectory / "optimal_port_eigenvalue_decay.csv");
    decay << "interface_id,mode,eigenvalue,relative_sqrt_decay\n"
        << std::setprecision(17);
    std::ofstream residual(
        outputDirectory / "optimal_port_eigenpair_residual.csv");
    residual << "interface_id,mode,eigenvalue,relative_residual\n"
        << std::setprecision(17);
    for (const LocalPortBasis& port : model.ports) {
        const double leading = port.spectralValues.empty()
            ? 0.0 : port.spectralValues.front();
        for (std::size_t mode = 0;
             mode < port.spectralValues.size(); ++mode) {
            const double relative = leading > 0.0
                ? std::sqrt(std::max(
                    0.0, port.spectralValues[mode] / leading))
                : 0.0;
            decay << port.interfaceId << ',' << mode << ','
                << port.spectralValues[mode] << ',' << relative << '\n';
            residual << port.interfaceId << ',' << mode << ','
                << port.spectralValues[mode] << ','
                << port.spectralResiduals[mode] << '\n';
        }
    }

    std::ofstream operators(
        outputDirectory / "optimal_port_operator_diagnostics.csv");
    operators << "interface_id,target_rows,source_rows,gram_minimum,"
        "gram_maximum,gram_regularization,schur_relative_asymmetry,"
        "adjoint_relative_error,explicit_column_reference_error,"
        "eigen_status,eigen_converged,eigen_iterations,"
        "eigen_operator_applies,eigenpair_residual,"
        "mandatory_rank,requested_transfer_rank,"
        "converged_transfer_rank,total_port_rank,"
        "peak_workspace_bytes,transfer_indicator\n"
        << std::setprecision(17);
    std::ofstream inner(
        outputDirectory / "optimal_port_inner_solver.csv");
    inner << "interface_id,requested_solver,actual_solver,solver_path,status,"
        "fallback_triggered,fallback_reason,setup_applies,solve_calls,"
        "setup_seconds,total_solve_seconds,mean_solve_seconds,"
        "max_solve_seconds,total_iterations,mean_iterations,"
        "max_iterations,final_relative_residual,max_relative_residual,"
        "relative_asymmetry,"
        "diagonal_shift,factor_bytes,a_tt_dimension,a_tt_nnz,"
        "a_tt_assembly_time_s,a_tt_factorization_time_s,"
        "a_tt_factor_bytes,reduced_correction_rank,w_setup_time_s,"
        "w_bytes,q_dimension,q_assembly_time_s,"
        "q_factorization_time_s,q_bytes,transfer_apply_time_s,"
        "transpose_apply_time_s,reference_target_action_error,"
        "reference_cross_action_error,"
        "reference_cross_transpose_error,"
        "reference_action_check_time_s,"
        "peak_incremental_memory_bytes,"
        "A_tt_solve_relative_residual,"
        "Q_solve_pre_refinement_residual,"
        "Q_solve_relative_residual,Q_refinement_iterations,"
        "woodbury_pre_refinement_residual,"
        "woodbury_post_refinement_residual,refinement_iterations,"
        "refinement_residual_0,refinement_residual_1,"
        "refinement_residual_2,refinement_residual_3,"
        "refinement_correction_relative_norm,"
        "refinement_reduction_factor,refinement_converged,"
        "refinement_triggered_solve_calls,"
        "refinement_time_s,pardiso_internal_refinement_steps,"
        "Q_min_abs_factor_diagonal,Q_max_abs_factor_diagonal,"
        "Q_factor_diagonal_ratio,"
        "woodbury_cancellation_factor\n"
        << std::setprecision(17);
    for (const auto& port : build.interfaces) {
        operators << port.interfaceId << ',' << port.targetRows << ','
            << port.sourceRows << ',' << port.gramMinimum << ','
            << port.gramMaximum << ',' << port.gramRegularization << ','
            << port.innerSolver.relativeAsymmetry << ','
            << port.adjointRelativeError << ','
            << port.explicitColumnReferenceError << ','
            << port.eigenStatus << ',' << (port.eigenConverged ? 1 : 0)
            << ',' << port.eigenIterations << ','
            << port.eigenOperatorApplies << ','
            << port.maximumEigenpairResidual << ','
            << port.mandatoryRank << ','
            << port.requestedTransferRank << ','
            << port.convergedTransferRank << ','
            << port.totalPortRank << ','
            << port.peakWorkspaceBytes << ','
            << port.transferIndicator << '\n';
        const double meanSolve = port.innerSolver.solveCalls > 0
            ? port.innerSolver.totalSolveSeconds
                / port.innerSolver.solveCalls
            : 0.0;
        const double meanIterations = port.innerSolver.solveCalls > 0
            ? static_cast<double>(port.innerSolver.totalIterations)
                / port.innerSolver.solveCalls
            : 0.0;
        inner << port.interfaceId << ','
            << port.innerSolver.requestedSolver << ','
            << port.innerSolver.actualSolver << ','
            << port.innerSolver.solver << ','
            << port.innerSolver.status << ','
            << (port.innerSolver.fallbackReason.empty() ? 0 : 1) << ','
            << port.innerSolver.fallbackReason << ','
            << port.innerSolver.setupApplies << ','
            << port.innerSolver.solveCalls << ','
            << port.innerSolver.setupSeconds << ','
            << port.innerSolver.totalSolveSeconds << ',' << meanSolve
            << ',' << port.innerSolver.maximumSolveSeconds << ','
            << port.innerSolver.totalIterations << ',' << meanIterations
            << ',' << port.innerSolver.maximumIterations << ','
            << port.innerSolver.finalRelativeResidual << ','
            << port.innerSolver.maximumRelativeResidual << ','
            << port.innerSolver.relativeAsymmetry << ','
            << port.innerSolver.diagonalShift << ','
            << port.innerSolver.factorBytes << ','
            << port.innerSolver.aTtDimension << ','
            << port.innerSolver.aTtNonzeros << ','
            << port.innerSolver.aTtAssemblySeconds << ','
            << port.innerSolver.aTtFactorizationSeconds << ','
            << port.innerSolver.aTtFactorBytes << ','
            << port.innerSolver.reducedCorrectionRank << ','
            << port.innerSolver.wSetupSeconds << ','
            << port.innerSolver.wBytes << ','
            << port.innerSolver.qDimension << ','
            << port.innerSolver.qAssemblySeconds << ','
            << port.innerSolver.qFactorizationSeconds << ','
            << port.innerSolver.qBytes << ','
            << port.innerSolver.transferApplySeconds << ','
            << port.innerSolver.transposeApplySeconds << ','
            << port.innerSolver.referenceTargetActionError << ','
            << port.innerSolver.referenceCrossActionError << ','
            << port.innerSolver.referenceCrossTransposeError << ','
            << port.innerSolver.referenceActionCheckSeconds << ','
            << port.innerSolver.peakIncrementalMemoryBytes << ','
            << port.innerSolver.aTtSolveRelativeResidual << ','
            << port.innerSolver.qSolvePreRefinementResidual << ','
            << port.innerSolver.qSolveRelativeResidual << ','
            << port.innerSolver.qRefinementIterations << ','
            << port.innerSolver.woodburyPreRefinementResidual << ','
            << port.innerSolver.woodburyPostRefinementResidual << ','
            << port.innerSolver.refinementIterations << ','
            << port.innerSolver.refinementResidual0 << ','
            << port.innerSolver.refinementResidual1 << ','
            << port.innerSolver.refinementResidual2 << ','
            << port.innerSolver.refinementResidual3 << ','
            << port.innerSolver.refinementCorrectionRelativeNorm << ','
            << port.innerSolver.refinementReductionFactor << ','
            << (port.innerSolver.refinementConverged ? 1 : 0) << ','
            << port.innerSolver.refinementTriggeredSolveCalls << ','
            << port.innerSolver.refinementSeconds << ','
            << port.innerSolver.pardisoInternalRefinementSteps << ','
            << port.innerSolver.qMinimumAbsoluteFactorDiagonal << ','
            << port.innerSolver.qMaximumAbsoluteFactorDiagonal << ','
            << port.innerSolver.qFactorDiagonalRatio << ','
            << port.innerSolver.woodburyCancellationFactor << '\n';
    }

    std::ofstream timing(
        outputDirectory / "optimal_port_timing.csv");
    timing << "patch_setup_seconds,gram_assembly_seconds,"
        "mandatory_mode_seconds,eigen_solve_seconds,"
        "orthogonalization_seconds,total_port_offline_seconds\n"
        << std::setprecision(17)
        << build.patchSetupSeconds << ',' << build.gramAssemblySeconds
        << ',' << build.mandatoryModeSeconds << ','
        << build.eigenSolveSeconds << ','
        << build.orthogonalizationSeconds << ','
        << build.totalSeconds << '\n';
}

void writeResidualKrylovDiagnostics(
    const ResidualKrylovBuildResult& build,
    const std::filesystem::path& outputDirectory)
{
    std::ofstream output(
        outputDirectory
        / "residual_krylov_interface_diagnostics.csv");
    output << "interface_id,target_dofs,source_dofs,constant_rank,"
        "geometry_rank,input_rank,boundary_rank,history_rank,"
        "raw_history_channels,active_history_channels,"
        "requested_history_rank,compressed_history_rank,"
        "deflated_history_channels,history_target_rhs,"
        "history_compression_method,"
        "history_compression_relative_error,"
        "history_compression_time_s,"
        "history_compression_workspace_bytes,"
        "history_compression_fingerprint,"
        "mandatory_rank_total,requested_randomized_rank,"
        "accepted_randomized_rank,raw_probe_columns,"
        "independent_probe_columns,deflated_probe_columns,"
        "probe_block_rank,requested_enrichment_rank,"
        "accepted_enrichment_rank,enrichment_sweeps,"
        "initial_max_probe_residual,final_max_probe_residual,"
        "residual_reduction_factor,target_solve_count,"
        "target_solve_phase33_calls,schur_apply_count,"
        "target_residual,weighted_adjoint_error,"
        "setup_time_s,peak_incremental_memory_bytes,"
        "deflation_limited,status,corrected,snapshot_used,"
        "pod_used,svd_used\n" << std::setprecision(17);
    for (const auto& row : build.interfaces) {
        output << row.interfaceId << ',' << row.targetRows << ','
            << row.sourceRows << ',' << row.constantRank << ','
            << row.geometryRank << ',' << row.inputRank << ','
            << row.boundaryRank << ',' << row.historyRank << ','
            << row.rawHistoryChannels << ','
            << row.activeHistoryChannels << ','
            << row.requestedHistoryRank << ','
            << row.compressedHistoryRank << ','
            << row.deflatedHistoryChannels << ','
            << row.historyTargetRightHandSides << ','
            << row.historyCompressionMethod << ','
            << row.historyCompressionRelativeError << ','
            << row.historyCompressionSeconds << ','
            << row.historyCompressionWorkspaceBytes << ','
            << row.historyCompressionFingerprint << ','
            << row.mandatoryRankTotal << ','
            << row.requestedRandomizedRank << ','
            << row.acceptedRandomizedRank << ','
            << row.rawProbeColumns << ','
            << row.independentProbeColumns << ','
            << row.deflatedProbeColumns << ','
            << row.probeBlockRank << ','
            << row.requestedEnrichmentRank << ','
            << row.acceptedEnrichmentRank << ','
            << row.enrichmentSweeps << ','
            << row.initialMaximumProbeResidual << ','
            << row.finalMaximumProbeResidual << ','
            << row.residualReductionFactor << ','
            << row.targetSolveCount << ','
            << row.innerSolver.solveCalls << ','
            << row.schurApplyCount << ','
            << row.innerSolver.maximumRelativeResidual << ','
            << row.weightedAdjointError << ','
            << row.totalSeconds << ','
            << row.peakIncrementalMemoryBytes << ','
            << (row.deflationLimited ? 1 : 0) << ','
            << row.status << ",0,0,0,0\n";
    }
    std::ofstream timing(
        outputDirectory / "residual_krylov_build_timing.csv");
    timing << "mandatory_mode_seconds,probe_setup_seconds,"
        "enrichment_seconds,orthogonalization_seconds,"
        "total_basis_seconds\n" << std::setprecision(17)
        << build.mandatoryModeSeconds << ','
        << build.probeSetupSeconds << ','
        << build.enrichmentSeconds << ','
        << build.orthogonalizationSeconds << ','
        << build.totalSeconds << '\n';
}

void writeRandomizedTransferDiagnostics(
    const port::RandomizedTransferBuildResult& build,
    const std::filesystem::path& outputDirectory)
{
    std::ofstream interfaceOutput(
        outputDirectory
        / "randomized_transfer_interface_diagnostics.csv");
    interfaceOutput
        << "method,interface_id,target_dofs,source_dofs,"
        "requested_rank,oversampling,probe_columns,accepted_rank,"
        "power_iterations,seed,apply_count,transpose_apply_count,"
        "target_solve_count,target_solve_phase33_calls,"
        "basis_error_indicator,orthogonality_error,"
        "weighted_adjoint_error,target_residual,"
        "probe_generation_time_s,transfer_apply_time_s,"
        "transpose_apply_time_s,qr_time_s,basis_build_time_s,"
        "probe_matrix_bytes,Y_matrix_bytes,qr_workspace_bytes,"
        "final_basis_bytes,peak_incremental_memory_bytes,"
        "inner_solver_requested,inner_solver_actual,"
        "source_mode,status,snapshot_used,fom_used_for_basis\n"
        << std::setprecision(17);
    for (const port::PortBasisResult& row : build.interfaces) {
        interfaceOutput
            << row.methodName << ',' << row.physicalInterfaceId << ','
            << row.targetDofs << ',' << row.sourceDofs << ','
            << row.requestedRank << ',' << row.oversampling << ','
            << row.probeColumns << ',' << row.acceptedRank << ','
            << row.powerIterations << ',' << row.seed << ','
            << row.applyCount << ',' << row.transposeApplyCount << ','
            << row.targetSolveCount << ','
            << row.targetSolvePhase33Calls << ','
            << row.basisErrorIndicator << ','
            << row.orthogonalityError << ','
            << row.weightedAdjointError << ',' << row.residual << ','
            << row.probeGenerationSeconds << ','
            << row.transferApplySeconds << ','
            << row.transposeApplySeconds << ',' << row.qrSeconds << ','
            << row.basisBuildTime << ',' << row.probeMatrixBytes << ','
            << row.sampleMatrixBytes << ',' << row.qrWorkspaceBytes
            << ',' << row.finalBasisBytes << ',' << row.memoryPeak << ','
            << row.innerSolver.requestedSolver << ','
            << row.innerSolver.actualSolver << ',' << row.sourceMode
            << ',' << row.status << ','
            << (row.snapshotUsed ? 1 : 0) << ','
            << (row.fomUsedForBasis ? 1 : 0) << '\n';
    }
    std::ofstream timing(
        outputDirectory / "randomized_transfer_timing.csv");
    timing
        << "patch_setup_time_s,source_setup_time_s,"
        "probe_generation_time_s,transfer_apply_time_s,"
        "transpose_apply_time_s,qr_time_s,"
        "basis_serialization_time_s,total_basis_time_s\n"
        << std::setprecision(17)
        << build.patchSetupSeconds << ',' << build.sourceSetupSeconds
        << ',' << build.probeGenerationSeconds << ','
        << build.transferApplySeconds << ','
        << build.transposeApplySeconds << ',' << build.qrSeconds << ','
        << build.serializationSeconds << ',' << build.totalSeconds
        << '\n';
    std::ofstream memory(
        outputDirectory / "randomized_transfer_memory.csv");
    memory
        << "interface_id,probe_matrix_bytes,Y_matrix_bytes,"
        "qr_workspace_bytes,final_basis_bytes,"
        "woodbury_peak_bytes,peak_incremental_memory_bytes\n"
        << std::setprecision(17);
    for (const port::PortBasisResult& row : build.interfaces) {
        memory << row.physicalInterfaceId << ','
            << row.probeMatrixBytes << ',' << row.sampleMatrixBytes << ','
            << row.qrWorkspaceBytes << ',' << row.finalBasisBytes << ','
            << row.innerSolver.peakIncrementalMemoryBytes << ','
            << row.memoryPeak << '\n';
    }
}

double basisOrthogonalityError(const LocalModel& local)
{
#ifdef USE_MKL_PARDISO
    if (local.rank <= 0) return 0.0;
    std::vector<double> gram(static_cast<std::size_t>(
        local.rank * local.rank), 0.0);
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
        local.rank, local.rank, local.interiorDofs, 1.0,
        local.basis.data(), local.interiorDofs,
        local.basis.data(), local.interiorDofs,
        0.0, gram.data(), local.rank);
    double maximum = 0.0;
    for (int column = 0; column < local.rank; ++column) {
        for (int row = 0; row <= column; ++row) {
            maximum = std::max(maximum, std::abs(
                gram[static_cast<std::size_t>(column * local.rank + row)]
                - (row == column ? 1.0 : 0.0)));
        }
    }
    return maximum;
#else
    double maximum = 0.0;
    for (int left = 0; left < local.rank; ++left) {
        for (int right = 0; right <= left; ++right) {
            const double product = dot(
                local.basis.data() + static_cast<std::size_t>(left * local.interiorDofs),
                local.basis.data() + static_cast<std::size_t>(right * local.interiorDofs),
                local.interiorDofs);
            maximum = std::max(maximum,
                std::abs(product - (left == right ? 1.0 : 0.0)));
        }
    }
    return maximum;
#endif
}

void reorthogonalizeExistingBasis(LocalModel& local)
{
    if (local.rank <= 0) return;
    std::vector<double> corrected;
    corrected.reserve(local.basis.size());
    for (int mode = 0; mode < local.rank; ++mode) {
        std::vector<double> column(
            local.basis.begin() + static_cast<std::ptrdiff_t>(
                mode * local.interiorDofs),
            local.basis.begin() + static_cast<std::ptrdiff_t>(
                (mode + 1) * local.interiorDofs));
        for (int pass = 0; pass < 2; ++pass) {
            for (int previous = 0; previous < mode; ++previous) {
                const double coefficient = dot(
                    corrected.data() + static_cast<std::size_t>(
                        previous * local.interiorDofs),
                    column.data(), local.interiorDofs);
                for (int row = 0; row < local.interiorDofs; ++row) {
                    column[static_cast<std::size_t>(row)] -= coefficient
                        * corrected[static_cast<std::size_t>(
                            previous * local.interiorDofs + row)];
                }
            }
        }
        const double magnitude = std::sqrt(normSquared(column));
        if (!(magnitude > 1.0e-14) || !std::isfinite(magnitude)) {
            throw std::runtime_error(
                "Existing Local Block Arnoldi basis became rank deficient during reorthogonalization.");
        }
        for (double& value : column) value /= magnitude;
        corrected.insert(corrected.end(), column.begin(), column.end());
    }
    local.basis = std::move(corrected);
}

struct ReducedSchurValidationReport {
    double maximumProjectionError = 0.0;
    double averageProjectionError = 0.0;
    double maximumOperatorError = 0.0;
    double averageOperatorError = 0.0;
    std::vector<double> randomAverageBySubdomain;
    std::vector<double> randomMaximumBySubdomain;
};

ReducedSchurValidationReport validateExistingInteriorBasis(
    const MatrixPartition& k,
    const MatrixPartition& c,
    const ddm_schur::InterfacePartition& partition,
    const std::vector<LocalModel>& locals,
    double timeStep,
    std::uint64_t seed,
    const std::filesystem::path& outputDirectory)
{
    if (!(timeStep > 0.0)) {
        throw std::runtime_error("Reduced Schur validation requires positive dt.");
    }
    const std::size_t gammaSize = partition.interfaceGlobalDofs.size();
    std::vector<SubdomainDirectSolver> localFactors;
    localFactors.reserve(locals.size());
    for (std::size_t slot = 0; slot < locals.size(); ++slot) {
        SparseMatrix localStep(locals[slot].interiorDofs);
        localStep.appendScaledEntries(c.interior[slot], 1.0 / timeStep);
        localStep.appendScaledEntries(k.interior[slot], 1.0);
        localStep.finalizeCsr();
        localFactors.emplace_back(locals[slot].interiorDofs,
            sparseMatrixEntries(localStep));
    }

    const auto applyInteriorGamma = [&](std::size_t slot,
                                        const std::vector<double>& x,
                                        std::vector<double>& y) {
        const auto& domain = partition.domains[slot];
        y.assign(domain.interiorGlobalDofs.size(), 0.0);
        for (std::size_t row = 0; row < y.size(); ++row) {
            for (const auto& entry : k.interiorInterface[slot][row]) {
                y[row] += entry.second * x[static_cast<std::size_t>(
                    locals[slot].gammaIndices[static_cast<std::size_t>(entry.first)])];
            }
            for (const auto& entry : c.interiorInterface[slot][row]) {
                y[row] += entry.second / timeStep * x[static_cast<std::size_t>(
                    locals[slot].gammaIndices[static_cast<std::size_t>(entry.first)])];
            }
        }
    };
    const auto applyGammaInterior = [&](std::size_t slot,
                                        const std::vector<double>& u,
                                        std::vector<double>& y) {
        const auto& domain = partition.domains[slot];
        y.assign(domain.interfaceGlobalDofs.size(), 0.0);
        for (std::size_t row = 0; row < y.size(); ++row) {
            for (const auto& entry : k.interfaceInterior[slot][row]) {
                y[row] += entry.second * u[static_cast<std::size_t>(entry.first)];
            }
            for (const auto& entry : c.interfaceInterior[slot][row]) {
                y[row] += entry.second / timeStep * u[static_cast<std::size_t>(entry.first)];
            }
        }
    };

    std::ofstream basisOut(outputDirectory / "local_dynamic_schur_basis_response_validation.csv");
    basisOut << "subdomain,rank,interior_dofs,interface_dofs,orthogonality_error,"
        << "random_vectors,average_projection_error,maximum_projection_error\n"
        << std::setprecision(17);
    std::mt19937_64 generator(seed);
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    ReducedSchurValidationReport report;
    report.randomAverageBySubdomain.assign(locals.size(), 0.0);
    report.randomMaximumBySubdomain.assign(locals.size(), 0.0);
    long double projectionSum = 0.0L;
    std::size_t projectionCount = 0;
    for (std::size_t slot = 0; slot < locals.size(); ++slot) {
        const LocalModel& local = locals[slot];
        double localMaximum = 0.0;
        long double localSum = 0.0L;
        constexpr int probes = 100;
        for (int probe = 0; probe < probes; ++probe) {
            std::vector<double> x(gammaSize, 0.0);
            for (double& value : x) value = distribution(generator);
            std::vector<double> rhs;
            applyInteriorGamma(slot, x, rhs);
            std::vector<double> solution;
            localFactors[slot].solve(rhs, solution);
            std::vector<double> coefficients(static_cast<std::size_t>(local.rank), 0.0);
            for (int mode = 0; mode < local.rank; ++mode) {
                coefficients[static_cast<std::size_t>(mode)] = dot(
                    local.basis.data() + static_cast<std::size_t>(mode * local.interiorDofs),
                    solution.data(), local.interiorDofs);
            }
            long double projectedNorm = 0.0L;
            long double errorNorm = 0.0L;
            for (int row = 0; row < local.interiorDofs; ++row) {
                long double projected = 0.0L;
                for (int mode = 0; mode < local.rank; ++mode) {
                    projected += static_cast<long double>(local.basis[static_cast<std::size_t>(
                        mode * local.interiorDofs + row)])
                        * coefficients[static_cast<std::size_t>(mode)];
                }
                projectedNorm += projected * projected;
                const long double error = static_cast<long double>(solution[static_cast<std::size_t>(row)]) - projected;
                errorNorm += error * error;
            }
            const double relative = std::sqrt(static_cast<double>(errorNorm))
                / std::max(1.0e-300, std::sqrt(static_cast<double>(errorNorm + projectedNorm)));
            localMaximum = std::max(localMaximum, relative);
            localSum += relative;
            projectionSum += relative;
            ++projectionCount;
        }
        const double average = static_cast<double>(localSum / probes);
        report.randomAverageBySubdomain[slot] = average;
        report.randomMaximumBySubdomain[slot] = localMaximum;
        report.maximumProjectionError = std::max(report.maximumProjectionError, localMaximum);
        basisOut << local.domainId << ',' << local.rank << ',' << local.interiorDofs << ','
            << local.gammaDofs << ',' << basisOrthogonalityError(local) << ',' << probes << ','
            << average << ',' << localMaximum << '\n';
    }
    report.averageProjectionError = projectionCount == 0 ? 0.0
        : static_cast<double>(projectionSum / projectionCount);

    SparseMatrix interfaceStep(static_cast<int>(gammaSize));
    interfaceStep.appendScaledEntries(c.interfaceBlock, 1.0 / timeStep);
    interfaceStep.appendScaledEntries(k.interfaceBlock, 1.0);
    interfaceStep.finalizeCsr();
    std::vector<local::DenseSymmetricFactor> factors;
    factors.reserve(locals.size());
    for (const LocalModel& local : locals) {
        std::vector<double> ar(local.kii.size(), 0.0);
        for (std::size_t entry = 0; entry < ar.size(); ++entry) {
            ar[entry] = local.kii[entry] + local.cii[entry] / timeStep;
        }
        factors.push_back(local::factorDenseSymmetric(ar, local.rank));
    }
    long double operatorSum = 0.0L;
    constexpr int operatorProbes = 8;
    for (int probe = 0; probe < operatorProbes; ++probe) {
        std::vector<double> x(gammaSize, 0.0);
        for (double& value : x) value = distribution(generator);
        const std::vector<double> base = interfaceStep.multiply(x);
        std::vector<double> fom = base;
        std::vector<double> rom = base;
        for (std::size_t slot = 0; slot < locals.size(); ++slot) {
            std::vector<double> rhs;
            applyInteriorGamma(slot, x, rhs);
            std::vector<double> u;
            localFactors[slot].solve(rhs, u);
            std::vector<double> image;
            applyGammaInterior(slot, u, image);
            for (std::size_t row = 0; row < image.size(); ++row) {
                fom[static_cast<std::size_t>(locals[slot].gammaIndices[row])] -= image[row];
            }
            const LocalModel& local = locals[slot];
            std::vector<double> reducedRhs(static_cast<std::size_t>(local.rank), 0.0);
            for (int mode = 0; mode < local.rank; ++mode) {
                for (int row = 0; row < local.gammaDofs; ++row) {
                    reducedRhs[static_cast<std::size_t>(mode)] +=
                        (local.kiGamma[static_cast<std::size_t>(mode * local.gammaDofs + row)]
                         + local.ciGamma[static_cast<std::size_t>(mode * local.gammaDofs + row)] / timeStep)
                        * x[static_cast<std::size_t>(local.gammaIndices[static_cast<std::size_t>(row)])];
                }
            }
            solveDenseSymmetric(factors[slot], reducedRhs);
            for (int row = 0; row < local.gammaDofs; ++row) {
                long double image = 0.0L;
                for (int mode = 0; mode < local.rank; ++mode) {
                    image += static_cast<long double>(
                        local.kGammaI[static_cast<std::size_t>(row * local.rank + mode)]
                        + local.cGammaI[static_cast<std::size_t>(row * local.rank + mode)] / timeStep)
                        * reducedRhs[static_cast<std::size_t>(mode)];
                }
                rom[static_cast<std::size_t>(local.gammaIndices[static_cast<std::size_t>(row)])]
                    -= static_cast<double>(image);
            }
        }
        long double error = 0.0L;
        long double reference = 0.0L;
        for (std::size_t row = 0; row < gammaSize; ++row) {
            const long double difference = static_cast<long double>(rom[row]) - fom[row];
            error += difference * difference;
            reference += static_cast<long double>(fom[row]) * fom[row];
        }
        const double relative = std::sqrt(static_cast<double>(error))
            / std::max(1.0e-300, std::sqrt(static_cast<double>(reference)));
        report.maximumOperatorError = std::max(report.maximumOperatorError, relative);
        operatorSum += relative;
    }
    report.averageOperatorError = static_cast<double>(operatorSum / operatorProbes);
    std::ofstream operatorOut(outputDirectory / "local_dynamic_schur_operator_validation.csv");
    operatorOut << "operator_probes,average_relative_error,maximum_relative_error,gate_tolerance,status\n"
        << std::setprecision(17) << operatorProbes << ',' << report.averageOperatorError << ','
        << report.maximumOperatorError << ",1e-4,"
        << (report.maximumOperatorError < 1.0e-4 ? "pass" : "fail") << '\n';
    return report;
}

struct SourceAlignedValidationReport {
    double averagePhysicalOperatorError = 0.0;
    double maximumPhysicalOperatorError = 0.0;
};

SourceAlignedValidationReport runSourceAlignedInterfaceValidation(
    const Mesh& mesh,
    const CaseConfig& physics,
    const ThermalDescriptorSystem& descriptor,
    const MatrixPartition& k,
    const MatrixPartition& c,
    const ddm_schur::InterfacePartition& partition,
    const std::vector<LocalModel>& locals,
    const std::vector<double>& thetaInitial,
    const std::vector<double>& boundaryOffset,
    const Options& options,
    const SparseMatrix& fullStep,
    const ReducedSchurValidationReport& randomReport,
    const std::filesystem::path& outputDirectory)
{
    const int steps = 10;
    const std::size_t gammaSize = partition.interfaceGlobalDofs.size();
    ddm_schur::Options schurOptions;
    schurOptions.maxIterations = options.interfaceMaxIterations;
    schurOptions.restart = options.interfaceRestart;
    schurOptions.relativeTolerance = std::min(options.interfaceTolerance, 1.0e-10);
    schurOptions.coarseLinearXY = options.coarseLinearXY;
    schurOptions.coarseLinearZ = options.coarseLinearZ;
    schurOptions.proxyEnabled = false;
    ddm_schur::DdmSchurSolver fomSchur(mesh, fullStep, physics, schurOptions);

    std::vector<SubdomainDirectSolver> localFactors;
    localFactors.reserve(locals.size());
    for (std::size_t slot = 0; slot < locals.size(); ++slot) {
        SparseMatrix localStep(locals[slot].interiorDofs);
        localStep.appendScaledEntries(c.interior[slot], 1.0 / options.timeStep);
        localStep.appendScaledEntries(k.interior[slot], 1.0);
        localStep.finalizeCsr();
        localFactors.emplace_back(locals[slot].interiorDofs,
            sparseMatrixEntries(localStep));
    }
    SparseMatrix interfaceStep(static_cast<int>(gammaSize));
    interfaceStep.appendScaledEntries(c.interfaceBlock, 1.0 / options.timeStep);
    interfaceStep.appendScaledEntries(k.interfaceBlock, 1.0);
    interfaceStep.finalizeCsr();
    std::vector<local::DenseSymmetricFactor> reducedFactors;
    reducedFactors.reserve(locals.size());
    for (const LocalModel& local : locals) {
        std::vector<double> reducedMatrix(local.kii.size(), 0.0);
        for (std::size_t entry = 0; entry < reducedMatrix.size(); ++entry) {
            reducedMatrix[entry] = local.kii[entry] + local.cii[entry] / options.timeStep;
        }
        reducedFactors.push_back(local::factorDenseSymmetric(
            reducedMatrix, local.rank));
    }
    const auto projectError = [](const LocalModel& local,
                                 const std::vector<double>& vector) {
        long double projectionSquared = 0.0L;
        long double errorSquared = 0.0L;
        for (int mode = 0; mode < local.rank; ++mode) {
            const double coefficient = dot(
                local.basis.data() + static_cast<std::size_t>(
                    mode * local.interiorDofs), vector.data(), local.interiorDofs);
            for (int row = 0; row < local.interiorDofs; ++row) {
                const double value = local.basis[static_cast<std::size_t>(
                    mode * local.interiorDofs + row)] * coefficient;
                projectionSquared += static_cast<long double>(value) * value;
            }
        }
        for (int row = 0; row < local.interiorDofs; ++row) {
            long double projected = 0.0L;
            for (int mode = 0; mode < local.rank; ++mode) {
                const double coefficient = dot(
                    local.basis.data() + static_cast<std::size_t>(
                        mode * local.interiorDofs), vector.data(), local.interiorDofs);
                projected += static_cast<long double>(local.basis[static_cast<std::size_t>(
                    mode * local.interiorDofs + row)]) * coefficient;
            }
            const long double error = static_cast<long double>(vector[static_cast<std::size_t>(row)])
                - projected;
            errorSquared += error * error;
        }
        return std::sqrt(static_cast<double>(errorSquared))
            / std::max(1.0e-300, std::sqrt(static_cast<double>(
                errorSquared + projectionSquared)));
    };
    const auto applyInteriorGamma = [&](std::size_t slot,
                                        const std::vector<double>& gamma,
                                        std::vector<double>& rhs) {
        rhs.assign(locals[slot].interiorDofs, 0.0);
        for (int row = 0; row < locals[slot].interiorDofs; ++row) {
            for (const auto& entry : k.interiorInterface[slot][static_cast<std::size_t>(row)]) {
                rhs[static_cast<std::size_t>(row)] += entry.second * gamma[static_cast<std::size_t>(
                    locals[slot].gammaIndices[static_cast<std::size_t>(entry.first)])];
            }
            for (const auto& entry : c.interiorInterface[slot][static_cast<std::size_t>(row)]) {
                rhs[static_cast<std::size_t>(row)] += entry.second / options.timeStep
                    * gamma[static_cast<std::size_t>(
                        locals[slot].gammaIndices[static_cast<std::size_t>(entry.first)])];
            }
        }
    };
    std::ofstream trajectoryOut(
        outputDirectory / "local_dynamic_schur_source_aligned_interface_trajectory.csv");
    trajectoryOut << "step,time_s,subdomain,local_interface_index,global_interface_index,value\n"
        << std::setprecision(17);
    std::vector<long double> transientSum(locals.size(), 0.0L);
    std::vector<long double> physicalSum(locals.size(), 0.0L);
    std::vector<double> transientMaximum(locals.size(), 0.0);
    std::vector<double> physicalMaximum(locals.size(), 0.0);
    std::ofstream operatorOut(
        outputDirectory / "local_dynamic_schur_physical_operator_validation.csv");
    operatorOut << "step,time_s,relative_operator_error,gate_tolerance,status\n"
        << std::setprecision(17);
    std::ofstream fomTimingOut(
        outputDirectory / "local_dynamic_schur_fom_timing.csv");
    fomTimingOut << "step,time_s,iterations,local_solve_seconds,interface_solve_seconds,"
        "fgmres_seconds,total_solve_seconds,total_seconds,interface_relative_residual,status\n"
        << std::setprecision(17);
    long double operatorErrorSum = 0.0L;
    double operatorErrorMaximum = 0.0;
    std::vector<double> previous = thetaInitial;
    const PowerWaveform waveform = makeBuiltinWaveform(
        options.waveform, descriptor.nominalPowersW, options.timeStep, steps, options.seed);
    for (int step = 1; step <= steps; ++step) {
        std::vector<double> rhs = descriptor.capacity.multiply(previous);
        for (double& value : rhs) value /= options.timeStep;
        addInput(descriptor, waveform.sample(step * options.timeStep), rhs);
        for (std::size_t row = 0; row < rhs.size(); ++row) rhs[row] += boundaryOffset[row];
        const ddm_schur::SolveResult fom = fomSchur.solve(rhs);
        if (fom.report.status != "success") {
            throw std::runtime_error(
                "source_aligned_validation_failed: FOM Dynamic Schur did not converge.");
        }
        fomTimingOut << step << ',' << step * options.timeStep << ','
            << fom.report.iterations << ',' << fom.report.localSolveSeconds << ','
            << fom.report.interfaceSolveSeconds << ',' << fom.report.fgmresSeconds << ','
            << fom.report.totalSolveSeconds << ',' << fom.report.totalSeconds << ','
            << fom.report.interfaceRelativeResidual << ',' << fom.report.status << '\n';
        previous = fom.temperature;
        std::vector<double> exactSchurImage;
        fomSchur.applyExactSchur(fom.interfaceSolution, exactSchurImage);
        std::vector<double> reducedSchurImage = interfaceStep.multiply(
            fom.interfaceSolution);
        for (std::size_t slot = 0; slot < locals.size(); ++slot) {
            const LocalModel& local = locals[slot];
            std::vector<double> reducedRhs(static_cast<std::size_t>(local.rank), 0.0);
            for (int mode = 0; mode < local.rank; ++mode) {
                for (int row = 0; row < local.gammaDofs; ++row) {
                    reducedRhs[static_cast<std::size_t>(mode)] +=
                        (local.kiGamma[static_cast<std::size_t>(
                            mode * local.gammaDofs + row)]
                         + local.ciGamma[static_cast<std::size_t>(
                            mode * local.gammaDofs + row)] / options.timeStep)
                        * fom.interfaceSolution[static_cast<std::size_t>(
                            local.gammaIndices[static_cast<std::size_t>(row)])];
                }
            }
            local::solveDenseSymmetric(reducedFactors[slot], reducedRhs);
            for (int row = 0; row < local.gammaDofs; ++row) {
                long double image = 0.0L;
                for (int mode = 0; mode < local.rank; ++mode) {
                    image += static_cast<long double>(
                        local.kGammaI[static_cast<std::size_t>(row * local.rank + mode)]
                        + local.cGammaI[static_cast<std::size_t>(row * local.rank + mode)]
                            / options.timeStep)
                        * reducedRhs[static_cast<std::size_t>(mode)];
                }
                reducedSchurImage[static_cast<std::size_t>(
                    local.gammaIndices[static_cast<std::size_t>(row)])] -=
                    static_cast<double>(image);
            }
        }
        long double operatorErrorSquared = 0.0L;
        long double operatorReferenceSquared = 0.0L;
        for (std::size_t row = 0; row < exactSchurImage.size(); ++row) {
            const long double error = static_cast<long double>(
                reducedSchurImage[row]) - exactSchurImage[row];
            operatorErrorSquared += error * error;
            operatorReferenceSquared += static_cast<long double>(
                exactSchurImage[row]) * exactSchurImage[row];
        }
        const double operatorError = std::sqrt(static_cast<double>(
            operatorErrorSquared)) / std::max(1.0e-300, std::sqrt(static_cast<double>(
                operatorReferenceSquared)));
        operatorErrorSum += operatorError;
        operatorErrorMaximum = std::max(operatorErrorMaximum, operatorError);
        operatorOut << step << ',' << step * options.timeStep << ',' << operatorError
            << ",1e-4," << (operatorError < 1.0e-4 ? "pass" : "fail") << '\n';
        for (std::size_t slot = 0; slot < locals.size(); ++slot) {
            const LocalModel& local = locals[slot];
            std::vector<double> interior(static_cast<std::size_t>(local.interiorDofs), 0.0);
            for (int row = 0; row < local.interiorDofs; ++row) {
                interior[static_cast<std::size_t>(row)] = fom.temperature[static_cast<std::size_t>(
                    local.interiorGlobal[static_cast<std::size_t>(row)])];
            }
            const double transientError = projectError(local, interior);
            transientSum[slot] += transientError;
            transientMaximum[slot] = std::max(transientMaximum[slot], transientError);
            std::vector<double> gamma(local.gammaDofs, 0.0);
            for (int row = 0; row < local.gammaDofs; ++row) {
                const int gammaIndex = local.gammaIndices[static_cast<std::size_t>(row)];
                gamma[static_cast<std::size_t>(row)] =
                    fom.interfaceSolution[static_cast<std::size_t>(gammaIndex)];
                trajectoryOut << step << ',' << step * options.timeStep << ','
                    << local.domainId << ',' << row << ','
                    << partition.domains[slot].interfaceGlobalDofs[static_cast<std::size_t>(row)]
                    << ',' << gamma[static_cast<std::size_t>(row)] << '\n';
            }
            std::vector<double> response;
            applyInteriorGamma(slot, fom.interfaceSolution, response);
            std::vector<double> responseSolution;
            localFactors[slot].solve(response, responseSolution);
            const double physicalError = projectError(local, responseSolution);
            physicalSum[slot] += physicalError;
            physicalMaximum[slot] = std::max(physicalMaximum[slot], physicalError);
        }
    }
    std::ofstream comparison(
        outputDirectory / "local_dynamic_schur_source_aligned_comparison.csv");
    comparison << "subdomain,rank,steps,transient_interior_average_error,"
        "transient_interior_maximum_error,random_interface_average_error,"
        "random_interface_maximum_error,physical_interface_average_error,"
        "physical_interface_maximum_error\n" << std::setprecision(17);
    for (std::size_t slot = 0; slot < locals.size(); ++slot) {
        comparison << locals[slot].domainId << ',' << locals[slot].rank << ',' << steps << ','
            << static_cast<double>(transientSum[slot] / steps) << ','
            << transientMaximum[slot] << ','
            << randomReport.randomAverageBySubdomain[slot] << ','
            << randomReport.randomMaximumBySubdomain[slot] << ','
            << static_cast<double>(physicalSum[slot] / steps) << ','
            << physicalMaximum[slot] << '\n';
    }
    operatorOut << "# average," << static_cast<double>(operatorErrorSum / steps)
        << "\n# maximum," << operatorErrorMaximum << '\n';
    std::cout << "[Source-aligned validation] completed " << steps
        << " FOM Dynamic Schur steps; trajectory and physical-response diagnostics written.\n";
    return {static_cast<double>(operatorErrorSum / steps), operatorErrorMaximum};
}

bool appendEnrichmentDirection(LocalModel& local,
                               std::vector<double> candidate,
                               double tolerance)
{
    if (local.rank >= local.interiorDofs) return false;
    const double originalNorm = std::sqrt(normSquared(candidate));
    if (!(originalNorm > 0.0)) return false;
    for (int pass = 0; pass < 2; ++pass) {
        for (int mode = 0; mode < local.rank; ++mode) {
            const double* basis = local.basis.data()
                + static_cast<std::size_t>(mode * local.interiorDofs);
            const double coefficient = dot(
                basis, candidate.data(), local.interiorDofs);
            for (int row = 0; row < local.interiorDofs; ++row) {
                candidate[static_cast<std::size_t>(row)] -= coefficient * basis[row];
            }
        }
    }
    const double magnitude = std::sqrt(normSquared(candidate));
    if (!(magnitude > tolerance * originalNorm)) return false;
    for (double& value : candidate) value /= magnitude;
    local.basis.insert(local.basis.end(), candidate.begin(), candidate.end());
    ++local.rank;
    return true;
}

EnrichmentResult enrichLocalDynamicBases(
    const ThermalDescriptorSystem& descriptor,
    const MatrixPartition& k,
    const MatrixPartition& c,
    const ddm_schur::InterfacePartition& partition,
    const std::vector<double>& boundaryOffset,
    const std::vector<double>& thetaInitial,
    const Options& options,
    int steps,
    std::vector<LocalModel>& locals)
{
    EnrichmentResult report;
    if (options.localPortEnrichmentRounds <= 0) return report;
    const auto totalStart = Clock::now();
    SparseMatrix fullStep(descriptor.dofs);
    fullStep.appendScaledEntries(descriptor.capacity, 1.0 / options.timeStep);
    fullStep.appendScaledEntries(descriptor.conductivity, 1.0);
    fullStep.finalizeCsr();
    const auto factorStart = Clock::now();
    SubdomainDirectSolver fomFactor(descriptor.dofs, sparseMatrixEntries(fullStep));
    std::vector<SubdomainDirectSolver> localFactors;
    localFactors.reserve(locals.size());
    for (std::size_t slot = 0; slot < locals.size(); ++slot) {
        SparseMatrix localStep(locals[slot].interiorDofs);
        localStep.appendScaledEntries(c.interior[slot], 1.0 / options.timeStep);
        localStep.appendScaledEntries(k.interior[slot], 1.0);
        localStep.finalizeCsr();
        localFactors.emplace_back(locals[slot].interiorDofs,
            sparseMatrixEntries(localStep));
    }
    report.factorizationSeconds = elapsed(factorStart);
    struct Candidate {
        double norm = 0.0;
        std::vector<double> values;
    };
    std::vector<std::vector<Candidate>> candidates(locals.size());
    std::vector<std::string> waveformNames{
        options.waveform,
        "single_step",
        "rectangular_pulse",
        "piecewise_multilevel",
        "variable_duty_cycle",
        "mixed_frequency",
        "asynchronous_hotspots",
        "unseen_waveform"};
    std::sort(waveformNames.begin(), waveformNames.end());
    waveformNames.erase(std::unique(waveformNames.begin(), waveformNames.end()),
        waveformNames.end());
    const int trainingSteps = std::min(steps,
        std::max(4, 2 + 4 * options.localPortEnrichmentRounds));
    for (std::size_t waveformIndex = 0;
         waveformIndex < waveformNames.size(); ++waveformIndex) {
        const PowerWaveform training = makeBuiltinWaveform(
            waveformNames[waveformIndex], descriptor.nominalPowersW,
            options.timeStep, std::max(steps, trainingSteps),
            waveformNames[waveformIndex] == options.waveform
                ? options.seed
                : options.seed + 101U + static_cast<std::uint64_t>(waveformIndex));
        std::vector<double> previous = thetaInitial;
        for (int step = 1; step <= trainingSteps; ++step) {
            std::vector<double> rhs = descriptor.capacity.multiply(previous);
            for (double& value : rhs) value /= options.timeStep;
            addInput(descriptor, training.sample(step * options.timeStep), rhs);
            for (std::size_t row = 0; row < rhs.size(); ++row) {
                rhs[row] += boundaryOffset[row];
            }
            const auto solveStart = Clock::now();
            std::vector<double> exact;
            fomFactor.solve(rhs, exact);
            report.solveSeconds += elapsed(solveStart);
            report.snapshots.temperature.push_back(exact);
            std::vector<double> increment(exact.size(), 0.0);
            for (std::size_t row = 0; row < exact.size(); ++row) {
                increment[row] = exact[row] - previous[row];
            }
            if (waveformNames[waveformIndex] == options.waveform) {
                report.snapshots.mandatoryTemperature.push_back(exact);
                report.snapshots.mandatoryTemperature.push_back(increment);
            }
            report.snapshots.temperature.push_back(std::move(increment));
            report.snapshots.flux.push_back(
                descriptor.conductivity.multiply(exact));

            std::vector<double> approximate(static_cast<std::size_t>(descriptor.dofs), 0.0);
            for (std::size_t gamma = 0;
                 gamma < partition.interfaceGlobalDofs.size(); ++gamma) {
                const int global = partition.interfaceGlobalDofs[gamma];
                approximate[static_cast<std::size_t>(global)] =
                    exact[static_cast<std::size_t>(global)];
            }
            for (std::size_t slot = 0; slot < locals.size(); ++slot) {
                const LocalModel& local = locals[slot];
                std::vector<double> coordinates(static_cast<std::size_t>(local.rank), 0.0);
                for (int mode = 0; mode < local.rank; ++mode) {
                    for (int row = 0; row < local.interiorDofs; ++row) {
                        coordinates[static_cast<std::size_t>(mode)] +=
                            local.basis[static_cast<std::size_t>(
                                mode * local.interiorDofs + row)]
                            * exact[static_cast<std::size_t>(
                                local.interiorGlobal[static_cast<std::size_t>(row)])];
                    }
                }
                for (int row = 0; row < local.interiorDofs; ++row) {
                    for (int mode = 0; mode < local.rank; ++mode) {
                        approximate[static_cast<std::size_t>(
                            local.interiorGlobal[static_cast<std::size_t>(row)])] +=
                            local.basis[static_cast<std::size_t>(
                                mode * local.interiorDofs + row)]
                            * coordinates[static_cast<std::size_t>(mode)];
                    }
                }
            }
            std::vector<double> fullResidual = fullStep.multiply(approximate);
            for (std::size_t row = 0; row < fullResidual.size(); ++row) {
                fullResidual[row] -= rhs[row];
            }
            report.snapshots.residual.push_back(fullResidual);
            for (std::size_t slot = 0; slot < locals.size(); ++slot) {
                std::vector<double> localResidual(
                    static_cast<std::size_t>(locals[slot].interiorDofs), 0.0);
                for (int row = 0; row < locals[slot].interiorDofs; ++row) {
                    localResidual[static_cast<std::size_t>(row)] =
                        -fullResidual[static_cast<std::size_t>(
                            locals[slot].interiorGlobal[static_cast<std::size_t>(row)])];
                }
                std::vector<double> correction;
                const auto localSolveStart = Clock::now();
                localFactors[slot].solve(localResidual, correction);
                report.solveSeconds += elapsed(localSolveStart);
                candidates[slot].push_back(
                    {std::sqrt(normSquared(correction)), std::move(correction)});
            }
            previous = std::move(exact);
        }
    }

    // The solves above are grouped by waveform so each trajectory can reuse
    // its own previous state.  Reorder the port snapshots by time and then by
    // waveform before the fixed-rank QR.  This keeps early columns from being
    // monopolized by one waveform/power-channel combination.
    const auto interleaveWaveforms = [&](std::vector<std::vector<double>>& family,
                                         int columnsPerStep) {
        if (family.empty()) return;
        const int columnsPerWaveform = trainingSteps * columnsPerStep;
        if (family.size() != waveformNames.size()
                * static_cast<std::size_t>(columnsPerWaveform)) {
            throw std::runtime_error(
                "[Local port] Unexpected enrichment snapshot-family size.");
        }
        std::vector<std::vector<double>> interleaved;
        interleaved.reserve(family.size());
        for (int step = 0; step < trainingSteps; ++step) {
            for (std::size_t waveformIndex = 0;
                 waveformIndex < waveformNames.size(); ++waveformIndex) {
                const std::size_t base = waveformIndex
                    * static_cast<std::size_t>(columnsPerWaveform)
                    + static_cast<std::size_t>(step * columnsPerStep);
                for (int column = 0; column < columnsPerStep; ++column) {
                    interleaved.push_back(std::move(
                        family[base + static_cast<std::size_t>(column)]));
                }
            }
        }
        family = std::move(interleaved);
    };
    interleaveWaveforms(report.snapshots.temperature, 2);
    interleaveWaveforms(report.snapshots.flux, 1);
    interleaveWaveforms(report.snapshots.residual, 1);

    const auto orthStart = Clock::now();
    const int maximumAdded = 4 * options.localPortEnrichmentRounds;
    for (std::size_t slot = 0; slot < locals.size(); ++slot) {
        LocalModel& local = locals[slot];
        std::sort(candidates[slot].begin(), candidates[slot].end(),
            [](const Candidate& left, const Candidate& right) {
                return left.norm > right.norm;
            });
        const int before = local.rank;
        int considered = 0;
        for (Candidate& candidate : candidates[slot]) {
            if (local.rank - before >= maximumAdded) break;
            ++considered;
            appendEnrichmentDirection(local, std::move(candidate.values),
                options.rankTolerance);
        }
        const int added = local.rank - before;
        ArnoldiHistoryRow history;
        history.moment = options.moments + 1;
        history.inputColumns = considered;
        history.addedRank = added;
        history.cumulativeRank = local.rank;
        history.deflatedColumns = considered - added;
        history.orthogonalityError = basisOrthogonalityError(local);
        history.basisBytes = local.basis.size() * sizeof(double);
        local.history.push_back(history);
        local.deflated += history.deflatedColumns;
        report.addedRank += added;

        local.cii = projectSparse(c.interior[slot], local.basis,
            local.interiorDofs, local.rank);
        local.kii = projectSparse(k.interior[slot], local.basis,
            local.interiorDofs, local.rank);
        local.ciGamma = projectInteriorInterface(c.interiorInterface[slot],
            local.basis, local.interiorDofs, local.rank, local.gammaDofs);
        local.kiGamma = projectInteriorInterface(k.interiorInterface[slot],
            local.basis, local.interiorDofs, local.rank, local.gammaDofs);
        local.cGammaI = projectInterfaceInterior(c.interfaceInterior[slot],
            local.basis, local.interiorDofs, local.rank);
        local.kGammaI = projectInterfaceInterior(k.interfaceInterior[slot],
            local.basis, local.interiorDofs, local.rank);
        local.reducedInput = projectInput(descriptor, local);
        local.reducedBoundary.assign(static_cast<std::size_t>(local.rank), 0.0);
        for (int row = 0; row < local.interiorDofs; ++row) {
            const int global = local.interiorGlobal[static_cast<std::size_t>(row)];
            for (int mode = 0; mode < local.rank; ++mode) {
                local.reducedBoundary[static_cast<std::size_t>(mode)] +=
                    local.basis[static_cast<std::size_t>(
                        mode * local.interiorDofs + row)]
                    * boundaryOffset[static_cast<std::size_t>(global)];
            }
        }
    }
    report.orthogonalizationSeconds = elapsed(orthStart);
    report.totalSeconds = elapsed(totalStart);
    return report;
}

} // namespace

void runLocalBlockArnoldiDynamicSchurWorkflow(
    const Mesh& mesh,
    const CaseConfig& physics,
    const Options& options,
    const std::filesystem::path& outputDirectory)
{
    const auto totalStart = Clock::now();
    if (options.integrator != "backward-euler") {
        throw std::runtime_error(
            "Local Block Arnoldi Milestone 4 supports fixed-dt backward Euler.");
    }
    if (!(options.timeStep > 0.0) || !(options.endTime > 0.0)) {
        throw std::runtime_error("Local dynamic Schur requires positive dt and end time.");
    }
    if (!(options.fullResidualTolerance > 0.0)
        || !std::isfinite(options.fullResidualTolerance)) {
        throw std::runtime_error(
            "--mor-full-residual-tolerance must be finite and positive.");
    }
    if (options.constructionTraceMode != "global-fom"
        && options.constructionTraceMode != "operator-coarse") {
        throw std::runtime_error(
            "[Local dynamic] Unsupported construction trace mode.");
    }
    const bool portPod = options.portReduction
        && options.portBasisMethod == "port-pod";
    const bool steklovSchur = options.portReduction
        && options.portBasisMethod == "steklov-schur";
    const bool optimalTransfer = options.portReduction
        && options.portBasisMethod == "optimal-transfer";
    const bool randomizedTransfer = options.portReduction
        && options.portBasisMethod == "randomized-transfer";
    const bool hybridRandomized = options.portReduction
        && options.portBasisMethod == "hybrid-randomized";
    const bool residualKrylov = options.portReduction
        && options.portBasisMethod == "residual-krylov";
    const bool mandatoryOnly = options.portReduction
        && options.portBasisMethod == "mandatory-only";
    const bool operatorInformedPort =
        residualKrylov || mandatoryOnly;
    const bool reducedMethod =
        options.localInteriorArnoldiReducedSchurValidation;
    const bool reducedValidation =
        reducedMethod && !options.skipReducedSchurValidation;
    if (options.sourceAlignedInterfaceValidation && !reducedMethod) {
        throw std::runtime_error(
            "--mor-source-aligned-validation requires "
            "--mor-transient-method local-interior-arnoldi-reduced-schur.");
    }
    std::string effectiveOptimalPortSourceMode =
        options.optimalPortSourceMode;
    if (optimalTransfer
        && (options.optimalPortAblation == "trace-transfer-only"
        || options.optimalPortAblation == "constant-geometry-trace"
        || options.optimalPortAblation == "original-mandatory-trace")) {
        effectiveOptimalPortSourceMode = "trace-only";
    } else if (optimalTransfer
        && (options.optimalPortAblation
                   == "generalized-transfer-only"
               || options.optimalPortAblation
                   == "constant-geometry-generalized")) {
        effectiveOptimalPortSourceMode = "generalized-dynamic";
    }
    if (options.portReduction && !portPod && !steklovSchur
        && !optimalTransfer && !randomizedTransfer
        && !hybridRandomized && !operatorInformedPort) {
        throw std::runtime_error("[Optimal port] Unsupported interface basis method.");
    }
    std::filesystem::create_directories(outputDirectory);
    const int steps = std::max(1, static_cast<int>(std::llround(
        options.endTime / options.timeStep)));
    const bool collectInterfaceGram =
        (steklovSchur || optimalTransfer || randomizedTransfer
            || hybridRandomized || operatorInformedPort)
        && !options.optimalPortTopologyAudit;
    const std::filesystem::path descriptorLoadPath =
        !options.portReduction
        ? localDynamicDescriptorCachePath(options.loadPath)
        : std::filesystem::path{};
    const std::filesystem::path descriptorSavePath =
        !options.portReduction
        ? localDynamicDescriptorCachePath(options.savePath)
        : std::filesystem::path{};
    const std::filesystem::path descriptorReadPath =
        !descriptorLoadPath.empty() ? descriptorLoadPath : descriptorSavePath;
    const auto descriptorFingerprintStart = Clock::now();
    const std::uint64_t descriptorInputFingerprint =
        thermalDescriptorInputFingerprint(
            mesh, physics, options.massType, collectInterfaceGram);
    const double descriptorFingerprintSeconds =
        elapsed(descriptorFingerprintStart);
    bool descriptorCacheHit = false;
    bool descriptorCacheSaved = false;
    double descriptorCacheLoadSeconds = 0.0;
    double descriptorCacheSaveSeconds = 0.0;
    double descriptorAssemblySeconds = 0.0;
    std::uintmax_t descriptorCacheBytes = 0;
    std::error_code descriptorExistsError;
    const bool descriptorCacheExists = !descriptorReadPath.empty()
        && std::filesystem::exists(
            descriptorReadPath, descriptorExistsError);
    if (descriptorExistsError) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Cannot inspect cache path.");
    }
    ThermalDescriptorSystem descriptor;
    if (descriptorCacheExists) {
        const auto loadStart = Clock::now();
        descriptor = loadThermalDescriptorCache(
            descriptorReadPath, descriptorInputFingerprint,
            static_cast<int>(mesh.nodes.size()));
        descriptorCacheLoadSeconds = elapsed(loadStart);
        descriptorCacheHit = true;
        if (!descriptorSavePath.empty()
            && descriptorSavePath != descriptorReadPath) {
            const auto saveStart = Clock::now();
            saveThermalDescriptorCache(
                descriptorSavePath, descriptorInputFingerprint, descriptor);
            descriptorCacheSaveSeconds = elapsed(saveStart);
            descriptorCacheSaved = true;
        }
    } else {
        const auto assemblyStart = Clock::now();
        descriptor = assembleThermalDescriptorSystem(
            mesh, physics, options.massType, collectInterfaceGram);
        descriptorAssemblySeconds = elapsed(assemblyStart);
        const std::filesystem::path writePath =
            !descriptorSavePath.empty()
            ? descriptorSavePath : descriptorLoadPath;
        if (!writePath.empty()) {
            const auto saveStart = Clock::now();
            saveThermalDescriptorCache(
                writePath, descriptorInputFingerprint, descriptor);
            descriptorCacheSaveSeconds = elapsed(saveStart);
            descriptorCacheSaved = true;
        }
    }
    const std::filesystem::path descriptorSizePath = descriptorCacheSaved
        ? descriptorSavePath
        : (descriptorCacheHit
        ? descriptorReadPath
        : (!descriptorSavePath.empty()
            ? descriptorSavePath : descriptorLoadPath));
    if (!descriptorSizePath.empty()) {
        std::error_code sizeError;
        descriptorCacheBytes =
            std::filesystem::file_size(descriptorSizePath, sizeError);
        if (sizeError) descriptorCacheBytes = 0;
    }
    std::cout << "[Thermal descriptor cache] "
              << (descriptorCacheHit ? "hit" : "miss")
              << ", bytes=" << descriptorCacheBytes
              << ", fingerprint=" << descriptorFingerprintSeconds
              << " s, assembly=" << descriptorAssemblySeconds
              << " s, load=" << descriptorCacheLoadSeconds
              << " s, save=" << descriptorCacheSaveSeconds << " s\n";
    const auto interfacePartitionStart = Clock::now();
    const ddm_schur::InterfacePartition partition =
        ddm_schur::buildInterfacePartition(mesh, descriptor.conductivity);
    const double interfacePartitionSeconds =
        elapsed(interfacePartitionStart);
    if (options.optimalPortTopologyAudit) {
        if (!optimalTransfer) {
            throw std::runtime_error(
                "[Optimal port] Topology audit requires --port-basis-method optimal-transfer.");
        }
        const std::vector<PortTopologyAudit> audit =
            auditOptimalPortTopology(
                mesh, physics, partition, descriptor.sourceSubdomains,
                options.optimalPortOversamplingLayers,
                options.optimalPortRank,
                options.optimalPortInnerSolver);
        writeOptimalPortTopologyAudit(
            physics.name, audit, descriptor.dofs,
            static_cast<int>(partition.domains.size()),
            options.optimalPortRank,
            options.optimalPortInnerSolver,
            outputDirectory);
        std::cout << "[Optimal port] topology audit only: case="
                  << physics.name << ", physical interfaces="
                  << audit.size()
                  << "; eigensolve/transient skipped.\n";
        return;
    }
    const auto matrixPartitionStart = Clock::now();
    MatrixPartition k = partitionMatrix(mesh, descriptor.conductivity, partition);
    MatrixPartition c = partitionMatrix(mesh, descriptor.capacity, partition);
    const double matrixPartitionSeconds = elapsed(matrixPartitionStart);
    const double operatorPreparationSeconds = elapsed(totalStart);

    const std::filesystem::path dynamicLoadPath =
        !options.portReduction
        ? localDynamicCachePath(options.loadPath)
        : std::filesystem::path{};
    const std::filesystem::path dynamicSavePath =
        !options.portReduction
        ? localDynamicCachePath(options.savePath)
        : std::filesystem::path{};
    const std::filesystem::path referenceLoadPath =
        !options.portReduction
        ? localDynamicReferenceCachePath(options.loadPath)
        : std::filesystem::path{};
    const std::filesystem::path referenceSavePath =
        !options.portReduction
        ? localDynamicReferenceCachePath(options.savePath)
        : std::filesystem::path{};

    const auto referenceStart = Clock::now();
    std::unique_ptr<SubdomainDirectSolver> steadyFactor;
    bool globalConstructionFactorUsed = false;
    double constructionGlobalFactorSeconds = 0.0;
    double constructionGlobalSolveSeconds = 0.0;
    const int constructionPardisoThreads = std::max(1, std::min(
        static_cast<int>(solverParallelWorkers()),
        std::max(options.localSolveThreads, options.localPardisoThreads)));
    auto ensureSteadyFactor = [&]() {
        if (!steadyFactor) {
            // Keep the global K factor lazy: cache hits and an analytic thermal
            // equilibrium must not pay this cold-construction cost.
            globalConstructionFactorUsed = true;
            const auto factorStart = Clock::now();
            ScopedDirectSolverMklThreads threads(constructionPardisoThreads);
            steadyFactor = std::make_unique<SubdomainDirectSolver>(
                descriptor.dofs,
                sparseMatrixEntries(descriptor.conductivity));
            constructionGlobalFactorSeconds += elapsed(factorStart);
        }
    };
    std::vector<double> reference;
    std::vector<double> boundaryOffset;
    bool referenceCacheHit = false;
    bool referenceCacheSaved = false;
    double referenceCacheLoadSeconds = 0.0;
    double referenceCacheSaveSeconds = 0.0;
    bool analyticReferenceUsed = false;
    double analyticReferenceRelativeResidual =
        std::numeric_limits<double>::quiet_NaN();
    std::uintmax_t referenceCacheBytes = 0;
    std::error_code referenceExistsError;
    const bool referenceCacheExists = !referenceLoadPath.empty()
        && std::filesystem::exists(referenceLoadPath, referenceExistsError);
    if (referenceExistsError) {
        throw std::runtime_error(
            "[Local dynamic reference cache] Cannot inspect cache path.");
    }
    if (referenceCacheExists) {
        const auto loadStart = Clock::now();
        LocalDynamicReferenceCache cached =
            loadLocalDynamicReference(referenceLoadPath, descriptor);
        reference = std::move(cached.reference);
        boundaryOffset = std::move(cached.boundaryOffset);
        referenceCacheLoadSeconds = elapsed(loadStart);
        referenceCacheHit = true;
        if (!referenceSavePath.empty()
            && referenceSavePath != referenceLoadPath) {
            const auto saveStart = Clock::now();
            saveLocalDynamicReference(
                referenceSavePath, descriptor, {reference, boundaryOffset});
            referenceCacheSaveSeconds = elapsed(saveStart);
            referenceCacheSaved = true;
        }
    } else {
        {
            // A constant temperature is an exact equilibrium when every active
            // Dirichlet value and convection ambient agrees. The residual test
            // below protects this shortcut against incompatible assembly data.
            std::optional<double> equilibrium;
            const auto acceptTemperature = [&](double value) {
                if (!equilibrium.has_value()) {
                    equilibrium = value;
                    return true;
                }
                const double scale = std::max({
                    1.0, std::abs(*equilibrium), std::abs(value)});
                return std::abs(value - *equilibrium) <= 1.0e-12 * scale;
            };
            bool compatible = true;
            for (const auto& condition : physics.dirichletConditions) {
                compatible = compatible
                    && acceptTemperature(condition.temperature);
            }
            for (const auto& condition : physics.convectionConditions) {
                if (condition.coefficient > 0.0) {
                    compatible = compatible
                        && acceptTemperature(condition.ambientTemperature);
                }
            }
            if (compatible && equilibrium.has_value()) {
                reference.assign(
                    static_cast<std::size_t>(descriptor.dofs), *equilibrium);
                const std::vector<double> kReference =
                    descriptor.conductivity.multiply(reference);
                boundaryOffset = descriptor.boundaryRhs;
                long double residualSquared = 0.0L;
                long double scaleSquared = 0.0L;
                for (int row = 0; row < descriptor.dofs; ++row) {
                    boundaryOffset[static_cast<std::size_t>(row)] -=
                        kReference[static_cast<std::size_t>(row)];
                    const long double residual = boundaryOffset[
                        static_cast<std::size_t>(row)];
                    residualSquared += residual * residual;
                    const long double rhs = descriptor.boundaryRhs[
                        static_cast<std::size_t>(row)];
                    const long double image = kReference[
                        static_cast<std::size_t>(row)];
                    scaleSquared += std::max(rhs * rhs, image * image);
                }
                analyticReferenceRelativeResidual = std::sqrt(
                    static_cast<double>(residualSquared
                        / std::max(1.0e-300L, scaleSquared)));
                analyticReferenceUsed =
                    analyticReferenceRelativeResidual <= 1.0e-10;
            }
        }
        if (!analyticReferenceUsed) {
            ensureSteadyFactor();
            const auto solveStart = Clock::now();
            steadyFactor->solveMultiple(
                descriptor.boundaryRhs, 1, reference,
                constructionPardisoThreads);
            constructionGlobalSolveSeconds += elapsed(solveStart);
            const std::vector<double> kReference =
                descriptor.conductivity.multiply(reference);
            boundaryOffset = descriptor.boundaryRhs;
            for (int row = 0; row < descriptor.dofs; ++row) {
                boundaryOffset[static_cast<std::size_t>(row)] -=
                    kReference[static_cast<std::size_t>(row)];
            }
        }
        const std::filesystem::path writePath = !referenceSavePath.empty()
            ? referenceSavePath : referenceLoadPath;
        if (!writePath.empty()) {
            const auto saveStart = Clock::now();
            saveLocalDynamicReference(
                writePath, descriptor, {reference, boundaryOffset});
            referenceCacheSaveSeconds = elapsed(saveStart);
            referenceCacheSaved = true;
        }
    }
    const double referenceSeconds = elapsed(referenceStart);
    const std::filesystem::path referenceSizePath = referenceCacheSaved
        ? referenceSavePath
        : (referenceCacheHit
        ? referenceLoadPath
        : (!referenceSavePath.empty() ? referenceSavePath : referenceLoadPath));
    if (!referenceSizePath.empty()) {
        std::error_code sizeError;
        referenceCacheBytes =
            std::filesystem::file_size(referenceSizePath, sizeError);
        if (sizeError) referenceCacheBytes = 0;
    }
    std::cout << "[Local dynamic reference cache] "
              << (referenceCacheHit ? "hit" : "miss")
              << ", bytes=" << referenceCacheBytes
              << ", setup=" << referenceSeconds << " s\n";

    // Global FOM solves are used only to generate compressed local interface
    // trace excitations.  A validated local-dynamic cache already contains the
    // resulting bases, so cache-load runs must not repeat these global solves.
    // No global state basis is retained or sliced.
    std::vector<double> staticResponses;
    std::vector<double> momentRhs;
    std::vector<double> dynamicResponses;
    std::vector<std::vector<double>> traceCandidates;
    double constructionTraceSetupSeconds = 0.0;
    int operatorCoarseAggregates = 0;
    int operatorTraceKrylovIterations = 0;
    int operatorTraceKrylovMaximumIterations = 0;
    double operatorTraceKrylovMaximumRelativeResidual = 0.0;
    const bool requireTraceCandidates = options.portReduction
        || options.loadPath.empty();
    if (requireTraceCandidates) {
        const auto traceSetupStart = Clock::now();
        if (options.constructionTraceMode == "operator-coarse") {
            OperatorCoarseTraceBuild coarse = buildOperatorCoarseTraces(
                mesh, descriptor, partition, k, options);
            traceCandidates = std::move(coarse.traces);
            operatorCoarseAggregates = coarse.aggregates;
            operatorTraceKrylovIterations = coarse.krylovIterations;
            operatorTraceKrylovMaximumIterations =
                coarse.krylovMaximumIterations;
            operatorTraceKrylovMaximumRelativeResidual =
                coarse.krylovMaximumRelativeResidual;
        } else {
            ensureSteadyFactor();
            // At s0=0 the first two source transfer moments are
            // X0=K^-1 B and X1=-K^-1 C X0. Solve each complete source block as
            // one multi-RHS operation so symbolic/numerical work is shared.
            const auto staticSolveStart = Clock::now();
            steadyFactor->solveMultiple(
                descriptor.input, descriptor.sourceChannels, staticResponses,
                constructionPardisoThreads);
            constructionGlobalSolveSeconds += elapsed(staticSolveStart);
            momentRhs.assign(staticResponses.size(), 0.0);
            for (int channel = 0; channel < descriptor.sourceChannels; ++channel) {
                std::vector<double> response(
                    static_cast<std::size_t>(descriptor.dofs), 0.0);
                std::copy_n(staticResponses.begin() + static_cast<std::ptrdiff_t>(
                    channel * descriptor.dofs), descriptor.dofs, response.begin());
                std::vector<double> image = descriptor.capacity.multiply(response);
                for (int row = 0; row < descriptor.dofs; ++row) {
                    momentRhs[static_cast<std::size_t>(
                        channel * descriptor.dofs + row)] =
                        -image[static_cast<std::size_t>(row)];
                }
            }
            const auto dynamicSolveStart = Clock::now();
            steadyFactor->solveMultiple(
                momentRhs, descriptor.sourceChannels, dynamicResponses,
                constructionPardisoThreads);
            constructionGlobalSolveSeconds += elapsed(dynamicSolveStart);
            traceCandidates.reserve(
                static_cast<std::size_t>(2 * descriptor.sourceChannels));
            for (int channel = 0; channel < descriptor.sourceChannels; ++channel) {
                traceCandidates.emplace_back(
                    staticResponses.begin() + static_cast<std::ptrdiff_t>(
                        channel * descriptor.dofs),
                    staticResponses.begin() + static_cast<std::ptrdiff_t>(
                        (channel + 1) * descriptor.dofs));
                traceCandidates.emplace_back(
                    dynamicResponses.begin() + static_cast<std::ptrdiff_t>(
                        channel * descriptor.dofs),
                    dynamicResponses.begin() + static_cast<std::ptrdiff_t>(
                        (channel + 1) * descriptor.dofs));
            }
            std::vector<double>().swap(staticResponses);
            std::vector<double>().swap(momentRhs);
            std::vector<double>().swap(dynamicResponses);
        }
        constructionTraceSetupSeconds = elapsed(traceSetupStart);
    }

    const auto localBuildStart = Clock::now();
    std::vector<LocalModel> locals;
    locals.reserve(partition.domains.size());
    double localSymbolic = 0.0;
    double localNumerical = 0.0;
    double localSolve = 0.0;
    double localOrthogonalization = 0.0;
    double localTraceBasisSeconds = 0.0;
    double localInputSetupSeconds = 0.0;
    double localFingerprintSeconds = 0.0;
    double localOrthogonalityAuditSeconds = 0.0;
    double localProjectionSeconds = 0.0;
    int totalRank = 0;
    int uniqueTemplates = 0;
    int reusedInstances = 0;
    std::vector<std::pair<double, double>> basisOrthogonalityRows;
    basisOrthogonalityRows.reserve(partition.domains.size());
    std::map<std::uint64_t, std::vector<std::size_t>> templateSlots;
    bool localDynamicCacheHit = false;
    double localDynamicCacheLoadSeconds = 0.0;
    double localDynamicCacheSaveSeconds = 0.0;
    std::uintmax_t localDynamicCacheBytes = 0;
    if (!dynamicLoadPath.empty()) {
        const auto loadStart = Clock::now();
        locals = loadLocalDynamicModels(
            dynamicLoadPath, descriptor, options, partition);
        localDynamicCacheLoadSeconds = elapsed(loadStart);
        localDynamicCacheHit = true;
        std::set<int> templateIdentifiers;
        for (const LocalModel& localModel : locals) {
            totalRank += localModel.rank;
            templateIdentifiers.insert(localModel.templateId);
            reusedInstances += localModel.templateReused ? 1 : 0;
            basisOrthogonalityRows.emplace_back(0.0, 0.0);
        }
        uniqueTemplates = static_cast<int>(templateIdentifiers.size());
        std::error_code sizeError;
        localDynamicCacheBytes =
            std::filesystem::file_size(dynamicLoadPath, sizeError);
        if (sizeError) localDynamicCacheBytes = 0;
        std::cout << "[Local dynamic cache] loaded "
                  << dynamicLoadPath << ", bytes="
                  << localDynamicCacheBytes << ", seconds="
                  << localDynamicCacheLoadSeconds << '\n';
    }
    if (!localDynamicCacheHit) {
        struct PreparedLocalBuild {
            LocalModel model;
            std::vector<double> excitation;
            std::vector<double> effectiveInput;
            std::vector<double> nominalPowers;
            std::size_t prototypeSlot =
                std::numeric_limits<std::size_t>::max();
            double traceBasisSeconds = 0.0;
            double inputSetupSeconds = 0.0;
            double fingerprintSeconds = 0.0;
            double orthogonalityAuditSeconds = 0.0;
            double projectionSeconds = 0.0;
        };
        const std::size_t domainCount = partition.domains.size();
        std::vector<PreparedLocalBuild> prepared(domainCount);
        for (std::size_t slot = 0; slot < domainCount; ++slot) {
            const auto& domain = partition.domains[slot];
            PreparedLocalBuild& item = prepared[slot];
            LocalModel& localModel = item.model;
            localModel.domainId = domain.domainId;
            localModel.interiorDofs = static_cast<int>(
                domain.interiorGlobalDofs.size());
            localModel.gammaDofs = static_cast<int>(
                domain.interfaceGlobalDofs.size());
            localModel.interiorGlobal = domain.interiorGlobalDofs;
            localModel.gammaIndices.reserve(domain.interfaceGlobalDofs.size());
            for (int global : domain.interfaceGlobalDofs) {
                localModel.gammaIndices.push_back(partition.globalToInterface[
                    static_cast<std::size_t>(global)]);
            }
            const auto traceBasisStart = Clock::now();
            item.excitation = traceBasis(
                domain, traceCandidates, options.rankTolerance,
                options.interfaceExcitationRank);
            item.traceBasisSeconds = elapsed(traceBasisStart);
            localModel.excitationRank = static_cast<int>(item.excitation.size())
                / std::max(1, localModel.gammaDofs);

            const auto inputSetupStart = Clock::now();
            for (int channel = 0; channel < descriptor.sourceChannels;
                 ++channel) {
                std::vector<double> column(
                    static_cast<std::size_t>(localModel.interiorDofs), 0.0);
                double magnitude = 0.0;
                for (int row = 0; row < localModel.interiorDofs; ++row) {
                    column[static_cast<std::size_t>(row)] = descriptor.input[
                        static_cast<std::size_t>(channel * descriptor.dofs
                            + domain.interiorGlobalDofs[
                                static_cast<std::size_t>(row)])];
                    magnitude = std::max(magnitude,
                        std::abs(column[static_cast<std::size_t>(row)]));
                }
                if (magnitude > 0.0) {
                    item.effectiveInput.insert(item.effectiveInput.end(),
                        column.begin(), column.end());
                    item.nominalPowers.push_back(
                        descriptor.nominalPowersW[
                            static_cast<std::size_t>(channel)]);
                    ++localModel.physicalChannels;
                }
            }
            couplingTimesBasis(k.interiorInterface[slot], item.excitation,
                localModel.gammaDofs, localModel.excitationRank, -1.0,
                item.effectiveInput);
            couplingTimesBasis(c.interiorInterface[slot], item.excitation,
                localModel.gammaDofs, localModel.excitationRank, -1.0,
                item.effectiveInput);
            std::vector<double> ones(
                static_cast<std::size_t>(localModel.interiorDofs), 1.0);
            const std::vector<double> uniformImage =
                k.interior[slot].multiply(ones);
            item.effectiveInput.insert(item.effectiveInput.end(),
                uniformImage.begin(), uniformImage.end());
            localModel.initialBlockRank = static_cast<int>(
                item.effectiveInput.size())
                / std::max(1, localModel.interiorDofs);
            item.inputSetupSeconds = elapsed(inputSetupStart);

            const auto fingerprintStart = Clock::now();
            std::uint64_t fingerprint = UINT64_C(1469598103934665603);
            hashValue(fingerprint, localModel.interiorDofs);
            hashValue(fingerprint, localModel.gammaDofs);
            hashValue(fingerprint, localModel.physicalChannels);
            hashValue(fingerprint, localModel.excitationRank);
            hashSparse(fingerprint, k.interior[slot]);
            hashSparse(fingerprint, c.interior[slot]);
            hashCoupling(fingerprint, k.interiorInterface[slot]);
            hashCoupling(fingerprint, k.interfaceInterior[slot]);
            hashCoupling(fingerprint, c.interiorInterface[slot]);
            hashCoupling(fingerprint, c.interfaceInterior[slot]);
            hashVector(fingerprint, item.nominalPowers);
            hashVector(fingerprint, item.excitation);
            hashVector(fingerprint, item.effectiveInput);
            localModel.templateFingerprint = fingerprint;
            item.fingerprintSeconds = elapsed(fingerprintStart);

            if (options.reuseIdenticalSubdomains) {
                const auto candidates = templateSlots.find(fingerprint);
                if (candidates != templateSlots.end()) {
                    for (std::size_t candidate : candidates->second) {
                        const PreparedLocalBuild& prototype = prepared[candidate];
                        if (sameSparse(k.interior[slot], k.interior[candidate])
                            && sameSparse(c.interior[slot], c.interior[candidate])
                            && sameCoupling(k.interiorInterface[slot],
                                k.interiorInterface[candidate])
                            && sameCoupling(k.interfaceInterior[slot],
                                k.interfaceInterior[candidate])
                            && sameCoupling(c.interiorInterface[slot],
                                c.interiorInterface[candidate])
                            && sameCoupling(c.interfaceInterior[slot],
                                c.interfaceInterior[candidate])
                            && prototype.nominalPowers == item.nominalPowers
                            && prototype.excitation == item.excitation
                            && prototype.effectiveInput == item.effectiveInput) {
                            item.prototypeSlot = candidate;
                            break;
                        }
                    }
                }
            }
            if (item.prototypeSlot
                    == std::numeric_limits<std::size_t>::max()) {
                localModel.templateId = uniqueTemplates++;
                templateSlots[fingerprint].push_back(slot);
            } else {
                localModel.templateReused = true;
                ++reusedInstances;
            }
        }

        const int localBuildWorkers = std::max(1, std::min({
            options.localSolveThreads, static_cast<int>(domainCount),
            static_cast<int>(solverParallelWorkers())}));
        // This is intentional two-level parallelism: OpenMP owns independent
        // subdomains, while each local PARDISO factor stays at the separately
        // bounded localPardisoThreads count to avoid nested oversubscription.
        std::vector<std::exception_ptr> buildErrors(domainCount);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(localBuildWorkers) if(localBuildWorkers > 1)
#endif
        for (int slotIndex = 0;
             slotIndex < static_cast<int>(domainCount); ++slotIndex) {
            const std::size_t slot = static_cast<std::size_t>(slotIndex);
            PreparedLocalBuild& item = prepared[slot];
            if (item.prototypeSlot
                    != std::numeric_limits<std::size_t>::max()) {
                continue;
            }
            try {
                LocalModel& localModel = item.model;
                ThermalDescriptorSystem localDescriptor;
                localDescriptor.dofs = localModel.interiorDofs;
                localDescriptor.sourceChannels = localModel.initialBlockRank;
                localDescriptor.capacity = c.interior[slot];
                localDescriptor.conductivity = k.interior[slot];
                localDescriptor.input = std::move(item.effectiveInput);
                localDescriptor.boundaryRhs.assign(
                    static_cast<std::size_t>(localModel.interiorDofs), 0.0);
                BlockArnoldiResult arnoldi = buildBlockArnoldiBasis(
                    localDescriptor, options.moments, options.expansionPoint,
                    options.rankTolerance, options.secondMomentEnergy,
                    options.secondMomentMaximumColumns,
                    std::max(1, options.localPardisoThreads));
                item.effectiveInput = std::move(localDescriptor.input);
                localModel.rank = arnoldi.rank;
                localModel.basis = std::move(arnoldi.basis);
                localModel.history = std::move(arnoldi.history);
                localModel.arnoldiTiming = arnoldi.timing;
            } catch (...) {
                buildErrors[slot] = std::current_exception();
            }
        }
        for (const auto& error : buildErrors) {
            if (error) std::rethrow_exception(error);
        }
        for (std::size_t slot = 0; slot < domainCount; ++slot) {
            PreparedLocalBuild& item = prepared[slot];
            if (item.prototypeSlot
                    == std::numeric_limits<std::size_t>::max()) {
                continue;
            }
            const LocalModel& prototype =
                prepared[item.prototypeSlot].model;
            item.model.rank = prototype.rank;
            item.model.basis = prototype.basis;
            item.model.history = prototype.history;
            item.model.templateId = prototype.templateId;
        }

        locals.resize(domainCount);
        basisOrthogonalityRows.resize(domainCount);
        std::vector<std::exception_ptr> projectionErrors(domainCount);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(localBuildWorkers) if(localBuildWorkers > 1)
#endif
        for (int slotIndex = 0;
             slotIndex < static_cast<int>(domainCount); ++slotIndex) {
            const std::size_t slot = static_cast<std::size_t>(slotIndex);
            try {
                PreparedLocalBuild& item = prepared[slot];
                LocalModel& localModel = item.model;
                const auto auditStart = Clock::now();
                const double orthogonalityBefore =
                    basisOrthogonalityError(localModel);
                if (orthogonalityBefore > 1.0e-12) {
                    reorthogonalizeExistingBasis(localModel);
                }
                basisOrthogonalityRows[slot] = {
                    orthogonalityBefore,
                    basisOrthogonalityError(localModel)};
                item.orthogonalityAuditSeconds = elapsed(auditStart);
                for (const ArnoldiHistoryRow& row : localModel.history) {
                    localModel.deflated += row.deflatedColumns;
                }

                const auto projectionStart = Clock::now();
                localModel.cii = projectSparse(c.interior[slot],
                    localModel.basis, localModel.interiorDofs, localModel.rank);
                localModel.kii = projectSparse(k.interior[slot],
                    localModel.basis, localModel.interiorDofs, localModel.rank);
                localModel.ciGamma = projectInteriorInterface(
                    c.interiorInterface[slot], localModel.basis,
                    localModel.interiorDofs, localModel.rank,
                    localModel.gammaDofs);
                localModel.kiGamma = projectInteriorInterface(
                    k.interiorInterface[slot], localModel.basis,
                    localModel.interiorDofs, localModel.rank,
                    localModel.gammaDofs);
                localModel.cGammaI = projectInterfaceInterior(
                    c.interfaceInterior[slot], localModel.basis,
                    localModel.interiorDofs, localModel.rank);
                localModel.kGammaI = projectInterfaceInterior(
                    k.interfaceInterior[slot], localModel.basis,
                    localModel.interiorDofs, localModel.rank);
                localModel.reducedInput = projectInput(descriptor, localModel);
                localModel.reducedBoundary.assign(
                    static_cast<std::size_t>(localModel.rank), 0.0);
                localModel.referenceInterior.resize(
                    static_cast<std::size_t>(localModel.interiorDofs));
                for (int row = 0; row < localModel.interiorDofs; ++row) {
                    const int global = localModel.interiorGlobal[
                        static_cast<std::size_t>(row)];
                    localModel.referenceInterior[
                        static_cast<std::size_t>(row)] = reference[
                            static_cast<std::size_t>(global)];
                    for (int mode = 0; mode < localModel.rank; ++mode) {
                        localModel.reducedBoundary[
                            static_cast<std::size_t>(mode)] +=
                            localModel.basis[static_cast<std::size_t>(
                                mode * localModel.interiorDofs + row)]
                            * boundaryOffset[static_cast<std::size_t>(global)];
                    }
                }
                item.projectionSeconds = elapsed(projectionStart);
                locals[slot] = std::move(localModel);
            } catch (...) {
                projectionErrors[slot] = std::current_exception();
            }
        }
        for (const auto& error : projectionErrors) {
            if (error) std::rethrow_exception(error);
        }
        for (const PreparedLocalBuild& item : prepared) {
            localTraceBasisSeconds += item.traceBasisSeconds;
            localInputSetupSeconds += item.inputSetupSeconds;
            localFingerprintSeconds += item.fingerprintSeconds;
            localOrthogonalityAuditSeconds +=
                item.orthogonalityAuditSeconds;
            localProjectionSeconds += item.projectionSeconds;
        }
        for (const LocalModel& localModel : locals) {
            localSymbolic +=
                localModel.arnoldiTiming.symbolicAnalysisSeconds;
            localNumerical +=
                localModel.arnoldiTiming.numericalFactorizationSeconds;
            localSolve += localModel.arnoldiTiming.multiRhsSolveSeconds;
            localOrthogonalization +=
                localModel.arnoldiTiming.orthogonalizationSeconds;
            totalRank += localModel.rank;
        }
    }
    const double localBuildSeconds = elapsed(localBuildStart);
    const auto postLocalSetupStart = Clock::now();

    ReducedSchurValidationReport reducedValidationReport;
    if (reducedValidation) {
        std::ofstream orthogonalityOut(
            outputDirectory / "local_dynamic_schur_basis_orthogonality.csv");
        orthogonalityOut << "subdomain,rank,basis_storage,basis_layout,"
            "orthogonality_before,orthogonality_after,reorthogonalized,"
            "gate_tolerance,status\n" << std::setprecision(17);
        for (std::size_t slot = 0; slot < locals.size(); ++slot) {
            const auto row = basisOrthogonalityRows[slot];
            orthogonalityOut << locals[slot].domainId << ',' << locals[slot].rank
                << ",LocalModel::basis,in-memory-column-major,"
                << row.first << ',' << row.second << ','
                << (row.first > 1.0e-12 ? 1 : 0) << ",1e-12,"
                << (row.second <= 1.0e-12 ? "pass" : "fail") << '\n';
        }
        reducedValidationReport = validateExistingInteriorBasis(
            k, c, partition, locals, options.timeStep, options.seed, outputDirectory);
        const ReducedSchurValidationReport& validation = reducedValidationReport;
        std::cout << "[Reduced Dynamic Schur validation] basis response average/max="
            << validation.averageProjectionError << '/' << validation.maximumProjectionError
            << ", operator average/max=" << validation.averageOperatorError << '/'
            << validation.maximumOperatorError << "\n";
        std::ofstream validationSummary(
            outputDirectory / "local_dynamic_schur_validation_summary.csv");
        validationSummary
            << "status,subdomains,total_local_rank,full_interface_dofs,"
            "average_projection_error,maximum_projection_error,"
            "average_operator_error,maximum_operator_error,operator_tolerance\n"
            << std::setprecision(17)
            << (validation.maximumOperatorError < 1.0e-4
                ? "random_diagnostic_passed" : "random_diagnostic_failed") << ','
            << locals.size() << ',' << totalRank << ','
            << partition.interfaceGlobalDofs.size() << ','
            << validation.averageProjectionError << ','
            << validation.maximumProjectionError << ','
            << validation.averageOperatorError << ','
            << validation.maximumOperatorError << ",1e-4\n";
        if (!(validation.maximumOperatorError < 1.0e-4)) {
            std::cout << "[Reduced Dynamic Schur validation] random-interface "
                "diagnostic failed; it is not a physical trajectory gate.\n";
        }
    }

    std::unique_ptr<LocalPortModel> portModel;
    std::unique_ptr<GlobalInterfaceCoarseModel>
        globalCoarseModel;
    std::unique_ptr<GlobalRandomizedSchurDiagnostics>
        globalRandomizedDiagnostics;
    std::unique_ptr<LocalPortReducedSchurSolver> portSolver;
    double portSnapshotSeconds = 0.0;
    double portBasisSeconds = 0.0;
    double portLocalPilotSeconds = 0.0;
    std::unique_ptr<OptimalPortBuildResult> optimalPortBuild;
    std::unique_ptr<ResidualKrylovBuildResult>
        residualKrylovBuild;
    std::unique_ptr<port::RandomizedTransferBuildResult>
        randomizedTransferBuild;
    if (portPod) {
        const auto portStart = Clock::now();
        std::filesystem::path portPath = options.savePath;
        if (!portPath.empty() && portPath.extension().empty()) {
            portPath /= "local_port_basis.bin";
        }
        std::filesystem::path loadPath = options.loadPath;
        if (!loadPath.empty() && loadPath.extension().empty()) {
            loadPath /= "local_port_basis.bin";
        }
        if (!loadPath.empty()) {
            portModel = std::make_unique<LocalPortModel>(
                loadLocalPortModel(loadPath));
            if (portModel->interfaceGlobalDofs != partition.interfaceGlobalDofs) {
                throw std::runtime_error(
                    "[Local port] Loaded basis does not match the current interface ordering.");
            }
        } else {
            LocalPortSnapshotFamilies families;
            const int requested = options.localPortRank > 0
                ? options.localPortRank : 32;
            const int maximumTraces = std::max(16, 2 * requested + 8);
            const int stride = std::max(1, static_cast<int>(traceCandidates.size())
                / maximumTraces);
            for (int trace = 0;
                 trace < static_cast<int>(traceCandidates.size())
                    && static_cast<int>(families.temperature.size()) < maximumTraces;
                 trace += stride) {
                families.temperature.push_back(
                    traceCandidates[static_cast<std::size_t>(trace)]);
                families.flux.push_back(descriptor.conductivity.multiply(
                    traceCandidates[static_cast<std::size_t>(trace)]));
                std::vector<double> residual = descriptor.capacity.multiply(
                    traceCandidates[static_cast<std::size_t>(trace)]);
                const std::vector<double> stiffness = descriptor.conductivity.multiply(
                    traceCandidates[static_cast<std::size_t>(trace)]);
                for (std::size_t row = 0; row < residual.size(); ++row) {
                    residual[row] = residual[row] / options.timeStep + stiffness[row];
                }
                families.residual.push_back(std::move(residual));
            }
            LocalPortOptions portOptions;
            portOptions.requestedRank = options.localPortRank;
            portOptions.discardedEnergyTolerance =
                options.localPortEnergyTolerance;
            portOptions.relativeTolerance = options.rankTolerance;
            portOptions.rankFile = options.localPortRankFile;
            portOptions.temperatureWeight = options.localPortTemperatureWeight;
            portOptions.fluxWeight = options.localPortFluxWeight;
            portOptions.residualWeight = options.localPortResidualWeight;
            portModel = std::make_unique<LocalPortModel>(buildLocalPortModel(
                mesh, physics, partition, families, portOptions));
            if (portPath.empty()) {
                portPath = outputDirectory / "local_port_basis.bin";
            }
            saveLocalPortModel(*portModel, portPath);
        }
        portSnapshotSeconds = portModel->snapshotSeconds;
        portBasisSeconds = portModel->basisSeconds;
        writePortRankDiagnostics(*portModel, outputDirectory);
        std::cout << "[Local port] physical interfaces=" << portModel->ports.size()
                  << ", full/reduced interface=" << portModel->fullInterfaceDofs
                  << '/' << portModel->reducedInterfaceDofs
                  << ", setup=" << elapsed(portStart) << " s\n";
    }

    const PowerWaveform waveform = options.inputPath.empty()
        ? makeBuiltinWaveform(options.waveform, descriptor.nominalPowersW,
            options.timeStep, steps, options.seed)
        : loadPowerWaveformCsv(options.inputPath, descriptor.sourceChannels);
    if (options.initialMode == "steady") {
        ensureSteadyFactor();
    }
    const std::vector<double> initial = initialTemperature(
        descriptor, options, reference, steadyFactor.get());
    std::vector<double> thetaInitial(initial.size(), 0.0);
    for (std::size_t row = 0; row < initial.size(); ++row) {
        thetaInitial[row] = initial[row] - reference[row];
    }
    EnrichmentResult enrichment;
    if (portPod || !options.portReduction) {
        enrichment = enrichLocalDynamicBases(
            descriptor, k, c, partition, boundaryOffset, thetaInitial,
            options, steps, locals);
    }
    if (enrichment.addedRank > 0) {
        totalRank = 0;
        for (LocalModel& local : locals) {
            totalRank += local.rank;
            if (local.history.back().addedRank > 0) {
                local.templateReused = false;
                local.templateId = local.domainId;
            }
        }
        localSolve += enrichment.solveSeconds;
        localOrthogonalization += enrichment.orthogonalizationSeconds;
        if (portPod) {
            const auto pilotStart = Clock::now();
            local::Model pilotModel = makeDynamicReducedModel(
                descriptor.dofs, partition, k, c, locals, options.timeStep);
            ddm_schur::Options pilotSchurOptions;
            pilotSchurOptions.maxIterations = options.interfaceMaxIterations;
            pilotSchurOptions.restart = options.interfaceRestart;
            pilotSchurOptions.relativeTolerance = options.interfaceTolerance;
            pilotSchurOptions.coarseLinearXY = options.coarseLinearXY;
            pilotSchurOptions.coarseLinearZ = options.coarseLinearZ;
            pilotSchurOptions.proxyEnabled = options.proxyEnabled;
            pilotSchurOptions.proxyDisableCoarse = options.proxyDisableCoarse;
            pilotSchurOptions.proxyHighConductivityThreshold =
                options.proxyHighConductivityThreshold;
            pilotSchurOptions.proxyUseMaterialConnectivity =
                options.proxyUseMaterialConnectivity;
            pilotSchurOptions.proxyRing = options.proxyRing;
            pilotSchurOptions.proxyProbeColumns = options.proxyProbeColumns;
            pilotSchurOptions.proxyBlockSize = options.proxyBlockSize;
            pilotSchurOptions.proxyValidateBlockEquivalence =
                options.proxyValidateBlockEquivalence;
            const std::filesystem::path pilotOutput =
                outputDirectory / "local_port_full_interface_pilot";
            std::filesystem::create_directories(pilotOutput);
            pilotSchurOptions.proxyOutputDirectory = pilotOutput.string();
            local::LocalReducedSchurSolver pilotSolver(
                pilotModel, mesh, physics, partition, pilotSchurOptions,
                options.matrixFreeInterfaceThreshold, pilotOutput);
            const int pilotSteps = std::min(steps,
                std::max(4, 2 + 4 * options.localPortEnrichmentRounds));
            const PowerWaveform pilotWaveform = options.inputPath.empty()
                ? makeBuiltinWaveform(options.waveform, descriptor.nominalPowersW,
                    options.timeStep, steps, options.seed)
                : loadPowerWaveformCsv(
                    options.inputPath, descriptor.sourceChannels);
            std::vector<double> previousPilot = thetaInitial;
            std::vector<double> previousPilotInterface(
                partition.interfaceGlobalDofs.size(), 0.0);
            for (std::size_t gammaIndex = 0;
                 gammaIndex < partition.interfaceGlobalDofs.size(); ++gammaIndex) {
                previousPilotInterface[gammaIndex] = thetaInitial[
                    static_cast<std::size_t>(
                        partition.interfaceGlobalDofs[gammaIndex])];
            }
            for (int step = 1; step <= pilotSteps; ++step) {
                std::vector<double> pilotRhs =
                    descriptor.capacity.multiply(previousPilot);
                for (double& value : pilotRhs) value /= options.timeStep;
                addInput(descriptor,
                    pilotWaveform.sample(step * options.timeStep), pilotRhs);
                for (std::size_t row = 0; row < pilotRhs.size(); ++row) {
                    pilotRhs[row] += boundaryOffset[row];
                }
                local::SolveResult pilotSolve = pilotSolver.solve(pilotRhs);
                if (pilotSolve.status != "success") {
                    throw std::runtime_error(
                        "[Local port] Full-interface local-ROM pilot did not converge.");
                }
                std::vector<double> trace(
                    static_cast<std::size_t>(descriptor.dofs), 0.0);
                std::vector<double> increment(
                    static_cast<std::size_t>(descriptor.dofs), 0.0);
                for (std::size_t gammaIndex = 0;
                     gammaIndex < partition.interfaceGlobalDofs.size(); ++gammaIndex) {
                    const int global = partition.interfaceGlobalDofs[gammaIndex];
                    const double value =
                        pilotSolve.interfaceTemperature[gammaIndex];
                    trace[static_cast<std::size_t>(global)] = value;
                    increment[static_cast<std::size_t>(global)] =
                        value - previousPilotInterface[gammaIndex];
                    previousPilotInterface[gammaIndex] = value;
                }
                enrichment.snapshots.mandatoryTemperature.push_back(
                    std::move(trace));
                enrichment.snapshots.mandatoryTemperature.push_back(
                    std::move(increment));
                previousPilot = std::move(pilotSolve.temperature);
            }
            portLocalPilotSeconds = elapsed(pilotStart);
            LocalPortOptions portOptions;
            portOptions.requestedRank = options.localPortRank;
            portOptions.discardedEnergyTolerance =
                options.localPortEnergyTolerance;
            portOptions.relativeTolerance = options.rankTolerance;
            portOptions.rankFile = options.localPortRankFile;
            portOptions.temperatureWeight = options.localPortTemperatureWeight;
            portOptions.fluxWeight = options.localPortFluxWeight;
            portOptions.residualWeight = options.localPortResidualWeight;
            portModel = std::make_unique<LocalPortModel>(buildLocalPortModel(
                mesh, physics, partition, enrichment.snapshots, portOptions));
            portSnapshotSeconds += portModel->snapshotSeconds;
            portBasisSeconds += portModel->basisSeconds;
            std::filesystem::path portPath = options.savePath;
            if (portPath.empty()) {
                portPath = outputDirectory / "local_port_basis.bin";
            } else if (portPath.extension().empty()) {
                portPath /= "local_port_basis.bin";
            }
            saveLocalPortModel(*portModel, portPath);
            writePortRankDiagnostics(*portModel, outputDirectory);
        }
        std::cout << "[Local port] residual enrichment added "
                  << enrichment.addedRank << " local modes in "
                  << enrichment.totalSeconds << " s\n";
    }
    std::vector<double> gamma(static_cast<std::size_t>(partition.interfaceGlobalDofs.size()), 0.0);
    for (std::size_t row = 0; row < gamma.size(); ++row) {
        gamma[row] = thetaInitial[static_cast<std::size_t>(partition.interfaceGlobalDofs[row])];
    }
    std::vector<std::vector<double>> localStates;
    localStates.reserve(locals.size());
    for (std::size_t slot = 0; slot < locals.size(); ++slot) {
        localStates.push_back(projectInitial(locals[slot], c.interior[slot], thetaInitial));
    }
    // The steady PARDISO factors and full trace snapshots are offline-only.
    // Release them before constructing the fixed-dt dynamic factors so the
    // RRAM26 workflow does not retain two full monolithic factorizations.
    steadyFactor.reset();
    staticResponses.clear();
    staticResponses.shrink_to_fit();
    momentRhs.clear();
    momentRhs.shrink_to_fit();
    dynamicResponses.clear();
    dynamicResponses.shrink_to_fit();
    traceCandidates.clear();
    traceCandidates.shrink_to_fit();

    std::unique_ptr<local::Model> preparedPortDynamicModel;
    if (steklovSchur || optimalTransfer || randomizedTransfer
        || hybridRandomized
        || operatorInformedPort) {
        const auto portStart = Clock::now();
        const double workflowBeforePortSeconds = elapsed(totalStart);
        const auto dynamicModelStart = Clock::now();
        preparedPortDynamicModel = std::make_unique<local::Model>(
            makeDynamicReducedModel(
                descriptor.dofs, partition, k, c, locals, options.timeStep));
        const double dynamicModelBuildSeconds =
            elapsed(dynamicModelStart);
        std::filesystem::path loadPath = options.loadPath;
        if (!loadPath.empty() && loadPath.extension().empty()) {
            loadPath /= "local_port_basis.bin";
        }
        const bool generalizedHistoryRequested =
            operatorInformedPort || hybridRandomized
            || ((optimalTransfer || randomizedTransfer)
                && effectiveOptimalPortSourceMode
                    == "generalized-dynamic");
        const bool streamHistoryFingerprintOnly =
            (options.globalRandomizedSchur
                || options.projectionDiagnosis)
            && !loadPath.empty();
        ReducedHistorySources reducedHistory =
            generalizedHistoryRequested
                && !streamHistoryFingerprintOnly
            ? buildReducedHistorySources(
                *preparedPortDynamicModel, locals, options.timeStep)
            : ReducedHistorySources{};
        const std::uint64_t generalizedHistoryFingerprint =
            generalizedHistoryRequested
            ? (streamHistoryFingerprintOnly
                ? fingerprintReducedHistorySourcesStreaming(
                    *preparedPortDynamicModel, locals,
                    options.timeStep)
                : fingerprintVector(
                    reducedHistory.condensedColumns))
            : fingerprintVector(
                reducedHistory.condensedColumns);
        if (hybridRandomized) {
            std::cout << "[Hybrid port] reduced history channels="
                << reducedHistory.channels << ", bytes="
                << reducedHistory.condensedColumns.size()
                    * sizeof(double)
                << '\n' << std::flush;
        }
        const std::uint64_t dynamicOperatorFingerprint =
            fingerprintDynamicPortOperator(*preparedPortDynamicModel);
        std::uint64_t penaltyFingerprint =
            UINT64_C(1469598103934665603);
        hashVector(
            penaltyFingerprint,
            descriptor.interfacePenaltyMassDiagonal);
        const std::uint64_t rankFileFingerprint =
            options.milestone8AdaptiveProduction
            ? milestone8AdaptiveProductionFingerprint()
            : fingerprintFile(options.optimalPortRankFile);
        if (options.optimalPortBasisPilot && !loadPath.empty()) {
            throw std::runtime_error(
                "[Optimal port] Basis pilot must build, not load, its interfaces.");
        }
        if (options.optimalPortRefinementValidation) {
            const std::vector<PortPatch> patches =
                buildOptimalPortPatches(
                    mesh, partition,
                    options.optimalPortOversamplingLayers);
            ReducedDynamicSchurOperator sharedRefinementSchur(
                *preparedPortDynamicModel, true);
            auto refinementOptions =
                [&](int refinementIterations) {
                    OptimalTransferPortOptions value;
                    value.requestedRank = options.optimalPortRank;
                    value.rankMode = options.optimalPortRankMode;
                    value.rankFile = options.optimalPortRankFile;
                    value.eigenvalueTolerance =
                        options.optimalPortEigenvalueTolerance;
                    value.minimumRank = options.optimalPortMinimumRank;
                    value.maximumRank = options.optimalPortMaximumRank;
                    value.innerProduct =
                        options.optimalPortInnerProduct;
                    value.oversamplingLayers =
                        options.optimalPortOversamplingLayers;
                    value.eigensolverMaximumIterations =
                        options.optimalPortEigenMaximumIterations;
                    value.eigensolverTolerance =
                        options.optimalPortEigenTolerance;
                    value.relativeDeflationTolerance =
                        options.rankTolerance;
                    value.ablationMode = options.optimalPortAblation;
                    value.sourceMode =
                        effectiveOptimalPortSourceMode;
                    value.innerSolver.innerSolver =
                        "woodbury-exact";
                    value.innerSolver.relativeTolerance =
                        options.optimalPortInnerTolerance;
                    value.innerSolver.maximumIterations =
                        options.optimalPortInnerMaximumIterations;
                    value.innerSolver.refinementMaximumIterations =
                        refinementIterations;
                    value.innerSolver.refinementTolerance =
                        options.optimalPortInnerRefinementTolerance;
                    return value;
                };
            const std::vector<int> refinementIterations{
                0,
                options.optimalPortInnerRefinementMaximumIterations};
            std::vector<OptimalPortBuildResult> builds;
            builds.reserve(refinementIterations.size());
            for (int iterations : refinementIterations) {
                builds.push_back(buildOptimalTransferPortModel(
                    mesh, partition, *preparedPortDynamicModel,
                    descriptor.interfaceTraceMassDiagonal,
                    descriptor.interfacePenaltyMassDiagonal,
                    descriptor.input, descriptor.sourceChannels,
                    descriptor.boundaryRhs,
                    reducedHistory.condensedColumns,
                    reducedHistory.channels,
                    refinementOptions(iterations),
                    &sharedRefinementSchur));
            }

            struct RefinementValidationRow {
                double solutionDifference = 0.0;
                double targetResidual = 0.0;
                double preResidual = 0.0;
                double postResidual = 0.0;
                int refinementIterations = 0;
                bool refinementConverged = true;
                double eigenvalueDifference = 0.0;
                double subspaceDifference = 0.0;
                double adjointError = 0.0;
                double totalBasisSeconds = 0.0;
                double targetSolveSeconds = 0.0;
                std::string status = "pass";
            };
            std::vector<RefinementValidationRow> rows(
                refinementIterations.size());
            for (std::size_t candidate = 0;
                 candidate < builds.size(); ++candidate) {
                rows[candidate].totalBasisSeconds =
                    builds[candidate].totalSeconds;
                for (const auto& diagnostics :
                     builds[candidate].interfaces) {
                    rows[candidate].targetResidual = std::max(
                        rows[candidate].targetResidual,
                        diagnostics.innerSolver.maximumRelativeResidual);
                    rows[candidate].preResidual = std::max(
                        rows[candidate].preResidual,
                        diagnostics.innerSolver
                            .woodburyPreRefinementResidual);
                    rows[candidate].postResidual = std::max(
                        rows[candidate].postResidual,
                        diagnostics.innerSolver
                            .woodburyPostRefinementResidual);
                    rows[candidate].refinementIterations = std::max(
                        rows[candidate].refinementIterations,
                        diagnostics.innerSolver.refinementIterations);
                    rows[candidate].refinementConverged =
                        rows[candidate].refinementConverged
                        && diagnostics.innerSolver.refinementConverged;
                    rows[candidate].adjointError = std::max(
                        rows[candidate].adjointError,
                        diagnostics.adjointRelativeError);
                    rows[candidate].targetSolveSeconds +=
                        diagnostics.innerSolver.totalSolveSeconds;
                }
            }
            for (const PortPatch& patch : patches) {
                std::vector<std::unique_ptr<PatchTransferOperator>>
                    targetSolvers;
                for (int iterations : refinementIterations) {
                    targetSolvers.push_back(
                        std::make_unique<PatchTransferOperator>(
                            sharedRefinementSchur, patch,
                            refinementOptions(iterations).innerSolver));
                }
                for (int probe = 0; probe < 2; ++probe) {
                    std::vector<double> rightHandSide(
                        patch.target.size(), 0.0);
                    for (std::size_t row = 0;
                         row < rightHandSide.size(); ++row) {
                        rightHandSide[row] =
                            std::sin((0.231 + 0.117 * probe)
                                * static_cast<double>(row + 1))
                            + 0.1 * std::cos(
                                0.419
                                * static_cast<double>(row + 1));
                    }
                    std::vector<double> disabledSolution;
                    targetSolvers.front()->solveTargetResponse(
                        rightHandSide, disabledSolution);
                    for (std::size_t candidate = 1;
                         candidate < targetSolvers.size();
                         ++candidate) {
                        std::vector<double> solution;
                        targetSolvers[candidate]->solveTargetResponse(
                            rightHandSide, solution);
                        rows[candidate].solutionDifference = std::max(
                            rows[candidate].solutionDifference,
                            relativeDifference(
                                solution, disabledSolution));
                    }
                }
                for (std::size_t candidate = 0;
                     candidate < targetSolvers.size();
                     ++candidate) {
                    const PatchInnerSolverStatistics& solver =
                        targetSolvers[candidate]->statistics();
                    rows[candidate].targetResidual = std::max(
                        rows[candidate].targetResidual,
                        solver.maximumRelativeResidual);
                    rows[candidate].preResidual = std::max(
                        rows[candidate].preResidual,
                        solver.woodburyPreRefinementResidual);
                    rows[candidate].postResidual = std::max(
                        rows[candidate].postResidual,
                        solver.woodburyPostRefinementResidual);
                    rows[candidate].refinementIterations = std::max(
                        rows[candidate].refinementIterations,
                        solver.refinementIterations);
                    rows[candidate].refinementConverged =
                        rows[candidate].refinementConverged
                        && solver.refinementConverged;
                    rows[candidate].targetSolveSeconds +=
                        solver.totalSolveSeconds;
                }
            }
            if (builds.front().model.ports.size()
                != builds.back().model.ports.size()) {
                rows.back().status = "basis_count_mismatch";
            } else {
                for (std::size_t port = 0;
                     port < builds.front().model.ports.size(); ++port) {
                    const LocalPortBasis& referencePort =
                        builds.front().model.ports[port];
                    const LocalPortBasis& candidatePort =
                        builds.back().model.ports[port];
                    if (candidatePort.spectralValues.size()
                        != referencePort.spectralValues.size()) {
                        rows.back().eigenvalueDifference =
                            std::numeric_limits<double>::infinity();
                    } else if (!candidatePort.spectralValues.empty()) {
                        rows.back().eigenvalueDifference = std::max(
                            rows.back().eigenvalueDifference,
                            relativeDifference(
                                candidatePort.spectralValues,
                                referencePort.spectralValues));
                    }
                    rows.back().subspaceDifference = std::max(
                        rows.back().subspaceDifference,
                        subspaceProjectorDifference(
                            candidatePort, referencePort));
                }
            }
            if (rows.back().solutionDifference > 1.0e-10
                || rows.back().targetResidual > 1.0e-9
                || rows.back().adjointError > 1.0e-10
                || rows.back().eigenvalueDifference > 1.0e-8
                || rows.back().subspaceDifference > 1.0e-8) {
                rows.back().status = "failed";
            }
            if (rows.front().targetResidual > 1.0e-9
                || rows.front().adjointError > 1.0e-10) {
                rows.front().status = "failed";
            }

            std::ofstream validation(
                outputDirectory
                / "milestone8_woodbury_refinement_validation.csv");
            validation
                << "case,refinement_enabled,"
                "refinement_max_iterations,interfaces,"
                "relative_solution_difference_vs_disabled,"
                "target_solve_residual,"
                "woodbury_pre_refinement_residual,"
                "woodbury_post_refinement_residual,"
                "refinement_iterations,refinement_converged,"
                "eigenvalue_relative_difference_vs_disabled,"
                "port_subspace_projector_difference_vs_disabled,"
                "weighted_adjoint_error,total_basis_time_s,"
                "target_solve_time_s,status\n"
                << std::setprecision(17);
            for (std::size_t candidate = 0;
                 candidate < rows.size(); ++candidate) {
                const RefinementValidationRow& row =
                    rows[candidate];
                validation << physics.name << ','
                    << (refinementIterations[candidate] > 0 ? 1 : 0)
                    << ',' << refinementIterations[candidate] << ','
                    << patches.size() << ','
                    << row.solutionDifference << ','
                    << row.targetResidual << ','
                    << row.preResidual << ','
                    << row.postResidual << ','
                    << row.refinementIterations << ','
                    << (row.refinementConverged ? 1 : 0) << ','
                    << row.eigenvalueDifference << ','
                    << row.subspaceDifference << ','
                    << row.adjointError << ','
                    << row.totalBasisSeconds << ','
                    << row.targetSolveSeconds << ','
                    << row.status << '\n';
            }
            std::cout
                << "[Optimal port] Woodbury refinement validation "
                << "complete; transient skipped.\n";
            return;
        } else if (options.optimalPortTargetSolverComparison) {
            const std::vector<std::string> solvers{
                "iterative-schur",
                "assembled-dense",
                "woodbury-exact"};
            const std::vector<PortPatch> patches =
                buildOptimalPortPatches(
                    mesh, partition,
                    options.optimalPortOversamplingLayers);
            ReducedDynamicSchurOperator sharedComparisonSchur(
                *preparedPortDynamicModel, true);
            auto comparisonOptions =
                [&](const std::string& solver) {
                    OptimalTransferPortOptions value;
                    value.requestedRank = options.optimalPortRank;
                    value.rankMode = options.optimalPortRankMode;
                    value.rankFile = options.optimalPortRankFile;
                    value.eigenvalueTolerance =
                        options.optimalPortEigenvalueTolerance;
                    value.minimumRank = options.optimalPortMinimumRank;
                    value.maximumRank = options.optimalPortMaximumRank;
                    value.innerProduct =
                        options.optimalPortInnerProduct;
                    value.oversamplingLayers =
                        options.optimalPortOversamplingLayers;
                    value.eigensolverMaximumIterations =
                        options.optimalPortEigenMaximumIterations;
                    value.eigensolverTolerance =
                        options.optimalPortEigenTolerance;
                    value.relativeDeflationTolerance =
                        options.rankTolerance;
                    value.ablationMode = options.optimalPortAblation;
                    value.sourceMode =
                        effectiveOptimalPortSourceMode;
                    value.innerSolver.innerSolver = solver;
                    value.innerSolver.relativeTolerance =
                        options.optimalPortInnerTolerance;
                    value.innerSolver.maximumIterations =
                        options.optimalPortInnerMaximumIterations;
                    value.innerSolver.refinementMaximumIterations =
                        options.optimalPortInnerRefinementMaximumIterations;
                    value.innerSolver.refinementTolerance =
                        options.optimalPortInnerRefinementTolerance;
                    return value;
                };
            std::vector<OptimalPortBuildResult> builds;
            builds.reserve(solvers.size());
            for (const std::string& solver : solvers) {
                builds.push_back(buildOptimalTransferPortModel(
                    mesh, partition, *preparedPortDynamicModel,
                    descriptor.interfaceTraceMassDiagonal,
                    descriptor.interfacePenaltyMassDiagonal,
                    descriptor.input, descriptor.sourceChannels,
                    descriptor.boundaryRhs,
                    reducedHistory.condensedColumns,
                    reducedHistory.channels,
                    comparisonOptions(solver),
                    &sharedComparisonSchur));
            }
            struct ComparisonRow {
                double solutionDifference = 0.0;
                double targetResidual = 0.0;
                double adjointError = 0.0;
                double eigenvalueDifference = 0.0;
                double subspaceDifference = 0.0;
                double setupSeconds = 0.0;
                std::string status = "pass";
            };
            std::vector<ComparisonRow> rows(solvers.size());
            for (std::size_t solver = 0;
                 solver < solvers.size(); ++solver) {
                rows[solver].setupSeconds =
                    builds[solver].totalSeconds;
                for (const auto& diagnostics :
                     builds[solver].interfaces) {
                    rows[solver].targetResidual = std::max(
                        rows[solver].targetResidual,
                        diagnostics.innerSolver.maximumRelativeResidual);
                    rows[solver].adjointError = std::max(
                        rows[solver].adjointError,
                        diagnostics.adjointRelativeError);
                }
            }
            for (const PortPatch& patch : patches) {
                std::vector<std::unique_ptr<PatchTransferOperator>>
                    targetSolvers;
                for (const std::string& solver : solvers) {
                    targetSolvers.push_back(
                        std::make_unique<PatchTransferOperator>(
                            sharedComparisonSchur, patch,
                            comparisonOptions(solver).innerSolver));
                }
                std::vector<std::vector<std::vector<double>>> solutions(
                    solvers.size());
                for (int probe = 0; probe < 2; ++probe) {
                    std::vector<double> rightHandSide(
                        patch.target.size(), 0.0);
                    for (std::size_t row = 0;
                         row < rightHandSide.size(); ++row) {
                        rightHandSide[row] =
                            std::sin((0.231 + 0.117 * probe)
                                * static_cast<double>(row + 1))
                            + 0.1 * std::cos(
                                0.419 * static_cast<double>(row + 1));
                    }
                    for (std::size_t solver = 0;
                         solver < solvers.size(); ++solver) {
                        std::vector<double> solution;
                        targetSolvers[solver]->solveTargetResponse(
                            rightHandSide, solution);
                        solutions[solver].push_back(
                            std::move(solution));
                    }
                }
                for (std::size_t solver = 0;
                     solver < solvers.size(); ++solver) {
                    for (int probe = 0; probe < 2; ++probe) {
                        rows[solver].solutionDifference = std::max(
                            rows[solver].solutionDifference,
                            relativeDifference(
                                solutions[solver][
                                    static_cast<std::size_t>(probe)],
                                solutions[1][
                                    static_cast<std::size_t>(probe)]));
                    }
                    rows[solver].targetResidual = std::max(
                        rows[solver].targetResidual,
                        targetSolvers[solver]->statistics()
                            .maximumRelativeResidual);
                }
            }
            for (std::size_t solver = 0;
                 solver < solvers.size(); ++solver) {
                if (builds[solver].model.ports.size()
                    != builds[1].model.ports.size()) {
                    rows[solver].status = "basis_count_mismatch";
                    continue;
                }
                for (std::size_t port = 0;
                     port < builds[solver].model.ports.size(); ++port) {
                    const LocalPortBasis& candidate =
                        builds[solver].model.ports[port];
                    const LocalPortBasis& referencePort =
                        builds[1].model.ports[port];
                    if (candidate.spectralValues.size()
                        != referencePort.spectralValues.size()) {
                        rows[solver].eigenvalueDifference =
                            std::numeric_limits<double>::infinity();
                    } else if (!candidate.spectralValues.empty()) {
                        rows[solver].eigenvalueDifference = std::max(
                            rows[solver].eigenvalueDifference,
                            relativeDifference(
                                candidate.spectralValues,
                                referencePort.spectralValues));
                    }
                    rows[solver].subspaceDifference = std::max(
                        rows[solver].subspaceDifference,
                        subspaceProjectorDifference(
                            candidate, referencePort));
                }
                if (rows[solver].solutionDifference > 1.0e-10
                    || rows[solver].targetResidual > 1.0e-9
                    || rows[solver].adjointError > 1.0e-10
                    || rows[solver].eigenvalueDifference > 1.0e-8
                    || rows[solver].subspaceDifference > 1.0e-8) {
                    rows[solver].status = "failed";
                }
            }
            std::ofstream comparison(
                outputDirectory
                / "milestone8_target_solver_comparison.csv");
            comparison
                << "case,solver,interfaces,relative_solution_difference,"
                "target_solve_residual,weighted_adjoint_error,"
                "eigenvalue_relative_difference,"
                "basis_projector_relative_difference,setup_time_s,status\n"
                << std::setprecision(17);
            for (std::size_t solver = 0;
                 solver < solvers.size(); ++solver) {
                comparison << physics.name << ',' << solvers[solver]
                    << ',' << patches.size() << ','
                    << rows[solver].solutionDifference << ','
                    << rows[solver].targetResidual << ','
                    << rows[solver].adjointError << ','
                    << rows[solver].eigenvalueDifference << ','
                    << rows[solver].subspaceDifference << ','
                    << rows[solver].setupSeconds << ','
                    << rows[solver].status << '\n';
            }
            std::cout
                << "[Optimal port] target-solver comparison complete; "
                << "transient skipped.\n";
            return;
        } else if (options.optimalPortWoodburyPilot) {
            OptimalTransferPortOptions pilotOptions;
            pilotOptions.requestedRank = 8;
            pilotOptions.requestedTransferRank = 4;
            pilotOptions.rankMode = "fixed";
            pilotOptions.eigenvalueTolerance =
                options.optimalPortEigenvalueTolerance;
            pilotOptions.minimumRank = 1;
            pilotOptions.maximumRank = 8;
            pilotOptions.innerProduct =
                options.optimalPortInnerProduct;
            pilotOptions.oversamplingLayers =
                options.optimalPortOversamplingLayers;
            pilotOptions.eigensolverMaximumIterations =
                options.optimalPortEigenMaximumIterations;
            pilotOptions.eigensolverTolerance =
                options.optimalPortEigenTolerance;
            pilotOptions.relativeDeflationTolerance =
                options.rankTolerance;
            pilotOptions.ablationMode =
                "constant-geometry-trace";
            pilotOptions.sourceMode = "trace-only";
            pilotOptions.selectedInterfaceIds = {16};
            pilotOptions.targetSolverPilotPreflight = true;
            pilotOptions.innerSolver.innerSolver =
                "woodbury-exact";
            pilotOptions.innerSolver.relativeTolerance =
                options.optimalPortInnerTolerance;
            pilotOptions.innerSolver.maximumIterations =
                options.optimalPortInnerMaximumIterations;
            pilotOptions.innerSolver.refinementMaximumIterations =
                options.optimalPortInnerRefinementMaximumIterations;
            pilotOptions.innerSolver.refinementTolerance =
                options.optimalPortInnerRefinementTolerance;
            const auto totalPilotStart = Clock::now();
            ReducedDynamicSchurOperator sharedPilotSchur(
                *preparedPortDynamicModel, false);
            pilotOptions.maximumPilotSeconds = std::max(
                1.0e-12, 600.0 - elapsed(totalPilotStart));
            OptimalPortBuildResult pilotBuild =
                buildOptimalTransferPortModel(
                    mesh, partition, *preparedPortDynamicModel,
                    descriptor.interfaceTraceMassDiagonal,
                    descriptor.interfacePenaltyMassDiagonal,
                    descriptor.input, descriptor.sourceChannels,
                    descriptor.boundaryRhs,
                    reducedHistory.condensedColumns,
                    reducedHistory.channels,
                    pilotOptions, &sharedPilotSchur);
            const double totalPilotSeconds =
                elapsed(totalPilotStart);
            if (pilotBuild.interfaces.size() != 1
                || pilotBuild.model.ports.size() != 1) {
                throw std::runtime_error(
                    "[Optimal port] Woodbury pilot did not return interface 16.");
            }
            const OptimalPortInterfaceDiagnostics& diagnostics =
                pilotBuild.interfaces.front();
            if (diagnostics.interfaceId != 16
                || diagnostics.targetRows != 3299
                || diagnostics.traceSourceRows != 4169) {
                throw std::runtime_error(
                    "[Optimal port] RRAM26 interface-16 topology changed.");
            }
            const PatchInnerSolverStatistics& solver =
                diagnostics.innerSolver;
            const double meanSolve = solver.solveCalls > 0
                ? solver.totalSolveSeconds / solver.solveCalls : 0.0;
            std::string pilotStatus =
                diagnostics.pilotStatus == "passed_preflight"
                ? "passed" : diagnostics.pilotStatus;
            if (pilotStatus == "not_requested"
                || pilotStatus == "running") {
                pilotStatus = "failed_internal_pilot_status";
            } else if (pilotStatus != "passed") {
                // The build stopped at the first failed preflight gate.
            } else if (solver.setupSeconds > 300.0) {
                pilotStatus = "failed_setup_time";
            } else if (meanSolve > 5.0) {
                pilotStatus = "failed_target_solve_time";
            } else if (diagnostics.operatorCheckSeconds > 60.0) {
                pilotStatus = "failed_operator_adjoint_time";
            } else if (totalPilotSeconds > 600.0) {
                pilotStatus = "failed_total_time";
            } else if (solver.peakIncrementalMemoryBytes
                    > UINT64_C(2147483648)) {
                pilotStatus = "failed_incremental_workspace";
            } else if (solver.maximumRelativeResidual > 1.0e-9) {
                pilotStatus = "failed_target_residual";
            } else if (diagnostics.adjointRelativeError > 1.0e-8) {
                pilotStatus = "failed_weighted_adjoint";
            } else if (!diagnostics.eigenConverged
                       || diagnostics.convergedTransferRank != 4) {
                pilotStatus = "failed_transfer_eigensolver";
            }

            std::ofstream pilot(
                outputDirectory
                / "milestone8_rram26_woodbury_pilot.csv");
            pilot
                << "interface_id,target_dofs,source_dofs,"
                "mandatory_rank,requested_transfer_rank,"
                "converged_transfer_rank,a_tt_dimension,a_tt_nnz,"
                "a_tt_assembly_time_s,a_tt_factorization_time_s,"
                "a_tt_factor_bytes,reduced_correction_rank,"
                "w_setup_time_s,w_bytes,q_dimension,"
                "q_assembly_time_s,q_factorization_time_s,q_bytes,"
                "target_solve_calls,mean_target_solve_time_s,"
                "max_target_solve_time_s,target_solve_residual,"
                "transfer_apply_time_s,transpose_apply_time_s,"
                "operator_adjoint_check_time_s,weighted_adjoint_error,"
                "eigensolver_iterations,eigenpair_residual,"
                "peak_incremental_memory_bytes,"
                "process_peak_memory_bytes,total_pilot_time_s,"
                "eigenvalues,pilot_status\n"
                << std::setprecision(17)
                << diagnostics.interfaceId << ','
                << diagnostics.targetRows << ','
                << diagnostics.traceSourceRows << ','
                << diagnostics.mandatoryRank << ','
                << diagnostics.requestedTransferRank << ','
                << diagnostics.convergedTransferRank << ','
                << solver.aTtDimension << ','
                << solver.aTtNonzeros << ','
                << solver.aTtAssemblySeconds << ','
                << solver.aTtFactorizationSeconds << ','
                << solver.aTtFactorBytes << ','
                << solver.reducedCorrectionRank << ','
                << solver.wSetupSeconds << ','
                << solver.wBytes << ','
                << solver.qDimension << ','
                << solver.qAssemblySeconds << ','
                << solver.qFactorizationSeconds << ','
                << solver.qBytes << ','
                << solver.solveCalls << ',' << meanSolve << ','
                << solver.maximumSolveSeconds << ','
                << solver.maximumRelativeResidual << ','
                << diagnostics.transferApplySeconds << ','
                << diagnostics.transposeApplySeconds << ','
                << diagnostics.operatorCheckSeconds << ','
                << diagnostics.adjointRelativeError << ','
                << diagnostics.eigenIterations << ','
                << diagnostics.maximumEigenpairResidual << ','
                << solver.peakIncrementalMemoryBytes << ','
                << peakWorkingSetBytes() << ','
                << totalPilotSeconds << ",\"";
            for (std::size_t eigen = 0;
                 eigen < diagnostics.eigenvalues.size(); ++eigen) {
                if (eigen != 0) pilot << ';';
                pilot << diagnostics.eigenvalues[eigen];
            }
            pilot << "\"," << pilotStatus << '\n';

            std::ofstream memory(
                outputDirectory
                / "milestone8_target_solver_memory.csv");
            memory
                << "case,interface_id,solver,a_tt_factor_bytes,"
                "w_bytes,q_bytes,peak_incremental_memory_bytes,"
                "process_peak_memory_bytes\n"
                << physics.name << ',' << diagnostics.interfaceId
                << ",woodbury-exact,"
                << solver.aTtFactorBytes << ','
                << solver.wBytes << ',' << solver.qBytes << ','
                << solver.peakIncrementalMemoryBytes << ','
                << peakWorkingSetBytes() << '\n';
            std::cout
                << "[Optimal port] RRAM26 interface-16 Woodbury pilot "
                << pilotStatus
                << "; all other interfaces and transient skipped.\n";
            return;
        } else if (options.randomizedPortRepresentativePilot
                   || options.historyCompressionMaximumInterfacePilot
                   || options.adaptivePortLocalPilot) {
            const bool historyCompressionPilot =
                options.historyCompressionMaximumInterfacePilot;
            const bool adaptivePilot =
                options.adaptivePortLocalPilot;
            std::vector<ScalabilityAuditRow> ordered =
                readScalabilityAudit(
                    options.optimalPortTopologyAuditCsv,
                    physics.name);
            std::sort(
                ordered.begin(), ordered.end(),
                [](const ScalabilityAuditRow& left,
                   const ScalabilityAuditRow& right) {
                    return left.targetDofs < right.targetDofs
                        || (left.targetDofs == right.targetDofs
                            && left.interfaceId < right.interfaceId);
                });
            std::vector<ScalabilityAuditRow> selected;
            std::vector<std::string> labels;
            if (adaptivePilot) {
                for (const int interfaceId :
                     options.adaptivePortInterfaceIds) {
                    const auto selectedRow = std::find_if(
                        ordered.begin(), ordered.end(),
                        [&](const ScalabilityAuditRow& row) {
                            return row.interfaceId == interfaceId;
                        });
                    if (selectedRow == ordered.end()) {
                        throw std::runtime_error(
                            "[Adaptive port] Requested interface is "
                            "absent from the topology audit.");
                    }
                    selected.push_back(*selectedRow);
                    labels.push_back(
                        "interface-"
                        + std::to_string(selectedRow->interfaceId));
                }
            } else if (historyCompressionPilot) {
                selected = {ordered.back()};
                labels = {"maximum"};
            } else {
                selected = {
                    ordered.front(),
                    ordered[ordered.size() / 2],
                    ordered.back()};
                labels = {"minimum", "median", "maximum"};
            }
            std::ofstream output(
                outputDirectory
                / (adaptivePilot
                    ? "milestone8_adaptive_port_local_pilot.csv"
                    : (historyCompressionPilot
                    ? "milestone8_rram26_history_compression_max_interface.csv"
                    : (hybridRandomized
                    ? "milestone8_hybrid_rram26_representative.csv"
                    : "milestone8_rram26_randomized_representative.csv"))));
            output
                << "case,selection,interface_id,target_dofs,"
                "source_dofs,requested_rank,oversampling,"
                "probe_columns,accepted_rank,power_iterations,seed,"
                "apply_count,transpose_apply_count,target_solve_count,"
                "target_solve_phase33_calls,basis_error_indicator,"
                "orthogonality_error,weighted_adjoint_error,"
                "target_residual,basis_build_time_s,"
                "probe_generation_time_s,transfer_apply_time_s,"
                "transpose_apply_time_s,qr_time_s,"
                "probe_matrix_bytes,Y_matrix_bytes,"
                "qr_workspace_bytes,final_basis_bytes,"
                "incremental_memory_bytes,process_peak_memory_bytes,"
                "inner_solver_requested,inner_solver_actual,status,"
                "transient_advanced,full_field_read,snapshot_used,"
                "fom_used_for_basis,mandatory_rank,"
                "requested_randomized_rank,accepted_randomized_rank,"
                "requested_enrichment_rank,accepted_enrichment_rank,"
                "initial_max_probe_residual,final_max_probe_residual,"
                "total_port_rank,randomized_target_solve_count,"
                "residual_target_solve_count,total_target_solve_count,"
                "randomized_build_time_s,residual_build_time_s,"
                "total_basis_build_time_s,raw_history_channels,"
                "active_history_channels,requested_history_rank,"
                "compressed_history_rank,deflated_history_channels,"
                "history_target_rhs,history_compression_method,"
                "history_compression_relative_error,"
                "history_compression_time_s,"
                "history_compression_workspace_bytes,"
                "history_compression_fingerprint\n"
                << std::setprecision(17);
            ReducedDynamicSchurOperator sharedRandomizedSchur(
                *preparedPortDynamicModel, true);
            bool passed = true;
            for (std::size_t index = 0;
                 index < selected.size(); ++index) {
                port::RandomizedTransferPortOptions pilotOptions;
                pilotOptions.fluxAware = options.fluxAwarePort;
                pilotOptions.fluxType = options.fluxAwareFluxType;
                pilotOptions.requestedRank =
                    options.randomizedPortRank;
                pilotOptions.oversampling =
                    options.randomizedPortOversampling;
                pilotOptions.powerIterations =
                    options.randomizedPortPowerIterations;
                pilotOptions.seed = options.randomizedPortSeed;
                pilotOptions.innerProduct =
                    options.optimalPortInnerProduct;
                pilotOptions.sourceMode = "trace-only";
                pilotOptions.oversamplingLayers =
                    options.optimalPortOversamplingLayers;
                pilotOptions.relativeDeflationTolerance =
                    options.rankTolerance;
                pilotOptions.selectedInterfaceIds = {
                    selected[index].interfaceId};
                pilotOptions.innerSolver.innerSolver =
                    "woodbury-exact";
                pilotOptions.innerSolver.relativeTolerance =
                    options.optimalPortInnerTolerance;
                pilotOptions.innerSolver.maximumIterations =
                    options.optimalPortInnerMaximumIterations;
                pilotOptions.innerSolver
                    .refinementMaximumIterations =
                    options
                        .optimalPortInnerRefinementMaximumIterations;
                pilotOptions.innerSolver.refinementTolerance =
                    options.optimalPortInnerRefinementTolerance;
                const port::RandomizedTransferBuildResult pilot =
                    port::buildRandomizedTransferPortModel(
                        mesh, physics, partition, *preparedPortDynamicModel,
                        descriptor.interfaceTraceMassDiagonal,
                        descriptor.interfacePenaltyMassDiagonal,
                        descriptor.input, descriptor.sourceChannels,
                        descriptor.boundaryRhs,
                        reducedHistory.condensedColumns,
                        reducedHistory.channels,
                        pilotOptions,
                        &sharedRandomizedSchur);
                if (hybridRandomized) {
                    std::cout << "[Hybrid port pilot] randomized "
                        << labels[index] << " interface complete\n"
                        << std::flush;
                }
                if (pilot.interfaces.size() != 1) {
                    throw std::runtime_error(
                        "[Randomized port] Representative pilot did "
                        "not return exactly one selected interface.");
                }
                const port::PortBasisResult& row =
                    pilot.interfaces.front();
                std::unique_ptr<ResidualKrylovBuildResult>
                    hybridResidualBuild;
                const ResidualKrylovInterfaceDiagnostics*
                    hybridResidual = nullptr;
                if (hybridRandomized) {
                    ResidualKrylovPortOptions residualOptions;
                    residualOptions.basisMethod =
                        "hybrid-randomized";
                    residualOptions.innerProduct =
                        options.optimalPortInnerProduct;
                    residualOptions.oversamplingLayers =
                        options.optimalPortOversamplingLayers;
                    residualOptions.maximumEnrichmentRank =
                        options.residualKrylovMaximumRank;
                    residualOptions.maximumSweeps =
                        options.residualKrylovMaximumSweeps;
                    residualOptions.blockSize =
                        options.residualKrylovBlockSize;
                    residualOptions.residualTolerance =
                        options.residualKrylovTolerance;
                    residualOptions.relativeDeflationTolerance =
                        options.rankTolerance;
                    residualOptions.probeMode =
                        options.residualKrylovProbeMode;
                    residualOptions.historyCompressionMethod =
                        options.historyCompressionMethod;
                    residualOptions.historyCompressionRank =
                        options.historyCompressionRank;
                    residualOptions.historyCompressionTolerance =
                        options.historyCompressionTolerance;
                    residualOptions.selectedInterfaceIds = {
                        selected[index].interfaceId};
                    residualOptions.innerSolver.innerSolver =
                        options.residualKrylovInnerSolver;
                    residualOptions.innerSolver.relativeTolerance =
                        options.optimalPortInnerTolerance;
                    residualOptions.innerSolver.maximumIterations =
                        options.optimalPortInnerMaximumIterations;
                    residualOptions.innerSolver
                        .refinementMaximumIterations =
                        options
                            .optimalPortInnerRefinementMaximumIterations;
                    residualOptions.innerSolver.refinementTolerance =
                        options.optimalPortInnerRefinementTolerance;
                    hybridResidualBuild = std::make_unique<
                        ResidualKrylovBuildResult>(
                        buildResidualKrylovPortModel(
                            mesh, partition,
                            *preparedPortDynamicModel,
                            descriptor.interfaceTraceMassDiagonal,
                            descriptor.interfacePenaltyMassDiagonal,
                            descriptor.input,
                            descriptor.sourceChannels,
                            descriptor.boundaryRhs,
                            reducedHistory.condensedColumns,
                            reducedHistory.channels,
                            residualOptions,
                            &sharedRandomizedSchur,
                            &pilot.model));
                    if (hybridResidualBuild->interfaces.size()
                            != 1
                        || hybridResidualBuild->model.ports.size()
                            != 1) {
                        throw std::runtime_error(
                            "[Hybrid port] Representative pilot did "
                            "not return exactly one selected interface.");
                    }
                    hybridResidual =
                        &hybridResidualBuild->interfaces.front();
                }
                const double combinedBuildSeconds =
                    row.basisBuildTime
                    + (hybridResidual == nullptr
                        ? 0.0 : hybridResidual->totalSeconds);
                const std::size_t combinedMemory =
                    std::max(
                        row.memoryPeak,
                        hybridResidual == nullptr
                        ? std::size_t{0}
                        : hybridResidual
                            ->peakIncrementalMemoryBytes);
                const int residualTargetSolves =
                    hybridResidual == nullptr
                    ? 0 : hybridResidual->targetSolveCount;
                const int totalTargetSolves =
                    row.targetSolveCount + residualTargetSolves;
                const int totalPhase33Calls =
                    row.targetSolvePhase33Calls
                    + (hybridResidual == nullptr
                        ? 0 : hybridResidual
                            ->innerSolver.solveCalls);
                const double combinedTargetResidual = std::max(
                    row.residual,
                    hybridResidual == nullptr
                    ? 0.0 : hybridResidual->innerSolver
                        .maximumRelativeResidual);
                const double combinedAdjointError = std::max(
                    row.weightedAdjointError,
                    hybridResidual == nullptr
                    ? 0.0
                    : hybridResidual->weightedAdjointError);
                std::string status = row.status;
                if (status.rfind("success", 0) == 0
                    && hybridResidual != nullptr
                    && hybridResidual->status != "success") {
                    status = hybridResidual->status;
                }
                if (status.rfind("success", 0) == 0
                    && row.targetDofs != selected[index].targetDofs) {
                    status = "audit_target_dimension_mismatch";
                }
                if (status.rfind("success", 0) == 0
                    && row.sourceDofs != selected[index].sourceDofs) {
                    status = "audit_source_dimension_mismatch";
                }
                if (status.rfind("success", 0) == 0
                    && (historyCompressionPilot || adaptivePilot
                        ? combinedBuildSeconds >= 120.0
                        : combinedBuildSeconds > 300.0)) {
                    status = "basis_build_time_gate_failed";
                }
                if (status.rfind("success", 0) == 0
                    && (historyCompressionPilot || adaptivePilot
                        ? combinedMemory >= UINT64_C(1073741824)
                        : combinedMemory > UINT64_C(1073741824))) {
                    status = "incremental_memory_gate_failed";
                }
                if (status.rfind("success", 0) == 0
                    && combinedTargetResidual > 1.0e-9) {
                    status = "target_residual_gate_failed";
                }
                if (status.rfind("success", 0) == 0
                    && combinedAdjointError > 1.0e-8) {
                    status = "weighted_adjoint_gate_failed";
                }
                output << physics.name << ',' << labels[index] << ','
                    << row.physicalInterfaceId << ',' << row.targetDofs
                    << ',' << row.sourceDofs << ',' << row.requestedRank
                    << ',' << row.oversampling << ',' << row.probeColumns
                    << ',' << row.acceptedRank << ','
                    << row.powerIterations << ',' << row.seed << ','
                    << row.applyCount << ','
                    << row.transposeApplyCount << ','
                    << totalTargetSolves << ','
                    << totalPhase33Calls << ','
                    << row.basisErrorIndicator << ','
                    << row.orthogonalityError << ','
                    << combinedAdjointError << ','
                    << combinedTargetResidual
                    << ',' << combinedBuildSeconds << ','
                    << row.probeGenerationSeconds << ','
                    << row.transferApplySeconds << ','
                    << row.transposeApplySeconds << ','
                    << row.qrSeconds << ',' << row.probeMatrixBytes
                    << ',' << row.sampleMatrixBytes << ','
                    << row.qrWorkspaceBytes << ','
                    << row.finalBasisBytes << ',' << combinedMemory
                    << ',' << peakWorkingSetBytes() << ','
                    << row.innerSolver.requestedSolver << ','
                    << row.innerSolver.actualSolver << ',' << status
                    << ",0,0,0,0,"
                    << (hybridResidual == nullptr
                        ? 0 : hybridResidual->mandatoryRankTotal)
                    << ',' << row.requestedRank << ','
                    << (hybridResidual == nullptr
                        ? row.acceptedRank
                        : hybridResidual
                            ->acceptedRandomizedRank)
                    << ',' << (hybridResidual == nullptr
                        ? 0 : hybridResidual
                            ->requestedEnrichmentRank)
                    << ',' << (hybridResidual == nullptr
                        ? 0 : hybridResidual
                            ->acceptedEnrichmentRank)
                    << ',' << (hybridResidual == nullptr
                        ? 0.0 : hybridResidual
                            ->initialMaximumProbeResidual)
                    << ',' << (hybridResidual == nullptr
                        ? 0.0 : hybridResidual
                            ->finalMaximumProbeResidual)
                    << ',' << (hybridResidual == nullptr
                        ? row.acceptedRank
                        : hybridResidualBuild->model.ports.front().rank)
                    << ',' << row.targetSolveCount
                    << ',' << residualTargetSolves
                    << ',' << totalTargetSolves
                    << ',' << row.basisBuildTime
                    << ',' << (hybridResidual == nullptr
                        ? 0.0 : hybridResidual->totalSeconds)
                    << ',' << combinedBuildSeconds
                    << ',' << (hybridResidual == nullptr
                        ? 0 : hybridResidual->rawHistoryChannels)
                    << ',' << (hybridResidual == nullptr
                        ? 0 : hybridResidual->activeHistoryChannels)
                    << ',' << (hybridResidual == nullptr
                        ? 0 : hybridResidual->requestedHistoryRank)
                    << ',' << (hybridResidual == nullptr
                        ? 0 : hybridResidual->compressedHistoryRank)
                    << ',' << (hybridResidual == nullptr
                        ? 0 : hybridResidual->deflatedHistoryChannels)
                    << ',' << (hybridResidual == nullptr
                        ? 0 : hybridResidual
                            ->historyTargetRightHandSides)
                    << ',' << (hybridResidual == nullptr
                        ? "none" : hybridResidual
                            ->historyCompressionMethod)
                    << ',' << (hybridResidual == nullptr
                        ? 0.0 : hybridResidual
                            ->historyCompressionRelativeError)
                    << ',' << (hybridResidual == nullptr
                        ? 0.0 : hybridResidual
                            ->historyCompressionSeconds)
                    << ',' << (hybridResidual == nullptr
                        ? std::size_t{0} : hybridResidual
                            ->historyCompressionWorkspaceBytes)
                    << ',' << (hybridResidual == nullptr
                        ? UINT64_C(0) : hybridResidual
                            ->historyCompressionFingerprint)
                    << '\n';
                std::cout
                    << "[Randomized port pilot] "
                    << labels[index] << " interface="
                    << row.physicalInterfaceId << ", rank="
                    << row.acceptedRank << ", target RHS="
                    << row.targetSolveCount << ", setup="
                    << combinedBuildSeconds << " s, status="
                    << status << '\n';
                if (status.rfind("success", 0) != 0) {
                    passed = false;
                    break;
                }
            }
            std::ofstream stop(
                outputDirectory
                / (adaptivePilot
                    ? "milestone8_adaptive_port_local_pilot_stop.csv"
                    : (historyCompressionPilot
                    ? "milestone8_history_compression_max_pilot_stop.csv"
                    : (hybridRandomized
                    ? "milestone8_hybrid_rram26_pilot_stop.csv"
                    : "milestone8_rram26_randomized_pilot_stop.csv"))));
            stop
                << "status,requested_rank,requested_history_rank,"
                "requested_residual_rank,interfaces_requested,"
                "transient_advanced,full_field_read,snapshot_used\n"
                << (passed
                    ? (adaptivePilot
                        ? "adaptive_local_pilot_passed"
                        : (historyCompressionPilot
                        ? "maximum_interface_passed"
                        : "representative_passed"))
                    : (adaptivePilot
                        ? "adaptive_local_pilot_failed"
                        : (historyCompressionPilot
                        ? "maximum_interface_gate_failed"
                        : "representative_gate_failed")))
                << ',' << options.randomizedPortRank
                << ',' << options.historyCompressionRank
                << ',' << options.residualKrylovMaximumRank
                << ',' << selected.size() << ",0,0,0\n";
            std::cout
                << (adaptivePilot
                    ? "[Adaptive port pilot] "
                    : (hybridRandomized
                    ? "[Hybrid port pilot] "
                    : "[Randomized port pilot] "))
                << "representative basis-only "
                << "run complete; transient skipped.\n";
            return;
        } else if (options.residualKrylovRepresentativePilot
                   || options.residualKrylovAllInterfaceBasis) {
            const bool allInterfaces =
                options.residualKrylovAllInterfaceBasis;
            const std::vector<ScalabilityAuditRow> auditRows =
                readScalabilityAudit(
                    options.optimalPortTopologyAuditCsv,
                    physics.name);
            std::vector<ScalabilityAuditRow> ordered = auditRows;
            std::sort(
                ordered.begin(), ordered.end(),
                [](const ScalabilityAuditRow& left,
                   const ScalabilityAuditRow& right) {
                    return left.targetDofs < right.targetDofs
                        || (left.targetDofs == right.targetDofs
                            && left.interfaceId < right.interfaceId);
                });
            std::vector<ScalabilityAuditRow> selected;
            std::vector<std::string> labels;
            if (allInterfaces) {
                const std::filesystem::path gate =
                    options.optimalPortTopologyAuditCsv.parent_path()
                    / "milestone8_rram26_residual_krylov_stage_b_gate.csv";
                std::ifstream gateInput(gate);
                std::string gateText(
                    (std::istreambuf_iterator<char>(gateInput)),
                    std::istreambuf_iterator<char>());
                if (!gateInput || gateText.find("stage_b_passed")
                        == std::string::npos) {
                    throw std::runtime_error(
                        "[Residual Krylov] Stage C is locked until "
                        "the representative-interface gate passes.");
                }
                selected = auditRows;
                std::sort(
                    selected.begin(), selected.end(),
                    [](const ScalabilityAuditRow& left,
                       const ScalabilityAuditRow& right) {
                        return left.interfaceId < right.interfaceId;
                    });
                labels.assign(selected.size(), "all");
            } else {
                selected = {
                    ordered.front(),
                    ordered[ordered.size() / 2],
                    ordered.back()};
                labels = {"minimum", "median", "maximum"};
            }
            std::vector<int> selectedIds;
            for (const auto& row : selected) {
                selectedIds.push_back(row.interfaceId);
            }
            ResidualKrylovPortOptions residualOptions;
            residualOptions.basisMethod = options.portBasisMethod;
            residualOptions.innerProduct =
                options.optimalPortInnerProduct;
            residualOptions.oversamplingLayers =
                options.optimalPortOversamplingLayers;
            residualOptions.maximumEnrichmentRank =
                mandatoryOnly ? 0
                    : options.residualKrylovMaximumRank;
            residualOptions.maximumSweeps =
                options.residualKrylovMaximumSweeps;
            residualOptions.blockSize =
                options.residualKrylovBlockSize;
            residualOptions.residualTolerance =
                options.residualKrylovTolerance;
            residualOptions.relativeDeflationTolerance =
                options.rankTolerance;
            residualOptions.probeMode =
                options.residualKrylovProbeMode;
            residualOptions.historyCompressionMethod =
                options.historyCompressionMethod;
            residualOptions.historyCompressionRank =
                options.historyCompressionRank;
            residualOptions.historyCompressionTolerance =
                options.historyCompressionTolerance;
            residualOptions.selectedInterfaceIds =
                selectedIds;
            residualOptions.innerSolver.innerSolver =
                options.residualKrylovInnerSolver;
            residualOptions.innerSolver.relativeTolerance =
                options.optimalPortInnerTolerance;
            residualOptions.innerSolver.maximumIterations =
                options.optimalPortInnerMaximumIterations;
            residualOptions.innerSolver.refinementMaximumIterations =
                options.optimalPortInnerRefinementMaximumIterations;
            residualOptions.innerSolver.refinementTolerance =
                options.optimalPortInnerRefinementTolerance;
            ReducedDynamicSchurOperator sharedResidualSchur(
                *preparedPortDynamicModel, true);
            ResidualKrylovBuildResult build =
                buildResidualKrylovPortModel(
                    mesh, partition, *preparedPortDynamicModel,
                    descriptor.interfaceTraceMassDiagonal,
                    descriptor.interfacePenaltyMassDiagonal,
                    descriptor.input, descriptor.sourceChannels,
                    descriptor.boundaryRhs,
                    reducedHistory.condensedColumns,
                    reducedHistory.channels,
                    residualOptions, &sharedResidualSchur);
            const std::filesystem::path csv = outputDirectory /
                (allInterfaces
                    ? "milestone8_rram26_residual_krylov_all_interfaces.csv"
                    : "milestone8_rram26_residual_krylov_representative.csv");
            std::ofstream output(csv);
            output << "case,selection,interface_id,target_dofs,"
                "source_dofs,constant_rank,geometry_rank,input_rank,"
                "boundary_rank,history_rank,mandatory_rank_total,"
                "raw_probe_columns,independent_probe_columns,"
                "deflated_probe_columns,requested_enrichment_rank,"
                "accepted_enrichment_rank,enrichment_sweeps,"
                "initial_max_probe_residual,final_max_probe_residual,"
                "residual_reduction_factor,target_solve_count,"
                "schur_apply_count,target_residual,"
                "weighted_adjoint_error,setup_time_s,"
                "peak_incremental_memory_bytes,deflation_limited,"
                "status,transient_advanced,full_field_read,"
                "snapshot_basis_used\n" << std::setprecision(17);
            bool passed = true;
            for (std::size_t index = 0;
                 index < build.interfaces.size(); ++index) {
                const auto& row = build.interfaces[index];
                const auto audit = std::find_if(
                    selected.begin(), selected.end(),
                    [&](const ScalabilityAuditRow& value) {
                        return value.interfaceId == row.interfaceId;
                    });
                if (audit == selected.end()
                    || audit->targetDofs != row.targetRows
                    || audit->sourceDofs != row.sourceRows) {
                    throw std::runtime_error(
                        "[Residual Krylov] Current topology differs "
                        "from the audited representative interface.");
                }
                const auto labelPosition = static_cast<std::size_t>(
                    std::distance(selected.begin(), audit));
                const bool rowPassed =
                    row.status == "success"
                    && row.totalSeconds <= 120.0
                    && row.innerSolver.maximumRelativeResidual
                        <= 1.0e-9
                    && row.weightedAdjointError <= 1.0e-8
                    && row.peakIncrementalMemoryBytes
                        <= UINT64_C(1073741824)
                    && (row.acceptedEnrichmentRank
                            >= row.requestedEnrichmentRank
                        || row.deflationLimited
                        || row.finalMaximumProbeResidual
                            <= options.residualKrylovTolerance)
                    && (row.requestedEnrichmentRank == 0
                        || row.finalMaximumProbeResidual
                            < row.initialMaximumProbeResidual);
                const std::string gateStatus =
                    row.innerSolver.maximumRelativeResidual
                            > 1.0e-9
                        ? "target_residual_gate_failed"
                        : (row.weightedAdjointError > 1.0e-8
                            ? "weighted_adjoint_gate_failed"
                            : (row.totalSeconds > 120.0
                                ? "single_interface_time_gate_failed"
                                : (row.peakIncrementalMemoryBytes
                                        > UINT64_C(1073741824)
                                    ? "incremental_memory_gate_failed"
                                    : (row.acceptedEnrichmentRank
                                            < row.requestedEnrichmentRank
                                        && !row.deflationLimited
                                        && row.finalMaximumProbeResidual
                                            > options.residualKrylovTolerance
                                        ? "enrichment_rank_gate_failed"
                                        : (row.requestedEnrichmentRank > 0
                                            && row.finalMaximumProbeResidual
                                                >= row.initialMaximumProbeResidual
                                            ? "probe_residual_gate_failed"
                                            : row.status)))));
                passed = passed && rowPassed;
                output << physics.name << ','
                    << labels[labelPosition] << ','
                    << row.interfaceId << ',' << row.targetRows << ','
                    << row.sourceRows << ',' << row.constantRank << ','
                    << row.geometryRank << ',' << row.inputRank << ','
                    << row.boundaryRank << ',' << row.historyRank << ','
                    << row.mandatoryRankTotal << ','
                    << row.rawProbeColumns << ','
                    << row.independentProbeColumns << ','
                    << row.deflatedProbeColumns << ','
                    << row.requestedEnrichmentRank << ','
                    << row.acceptedEnrichmentRank << ','
                    << row.enrichmentSweeps << ','
                    << row.initialMaximumProbeResidual << ','
                    << row.finalMaximumProbeResidual << ','
                    << row.residualReductionFactor << ','
                    << row.targetSolveCount << ','
                    << row.schurApplyCount << ','
                    << row.innerSolver.maximumRelativeResidual << ','
                    << row.weightedAdjointError << ','
                    << row.totalSeconds << ','
                    << row.peakIncrementalMemoryBytes << ','
                    << (row.deflationLimited ? 1 : 0) << ','
                    << (rowPassed ? "success" : gateStatus)
                    << ",0,0,0\n";
                if (!rowPassed) {
                    std::ofstream failure(
                        outputDirectory
                        / "milestone8_residual_krylov_failures.csv",
                        std::ios::app);
                    if (failure.tellp() == 0) {
                        failure << "case,interface_id,"
                            "requested_enrichment_rank,status\n";
                    }
                    failure << physics.name << ','
                        << row.interfaceId << ','
                        << row.requestedEnrichmentRank << ','
                        << gateStatus << '\n';
                }
            }
            writeResidualKrylovDiagnostics(build, outputDirectory);
            std::ofstream stop(
                outputDirectory
                / "milestone8_residual_krylov_stop.csv");
            stop << "status,interfaces_built,"
                "transient_advanced,full_field_read,"
                "snapshot_basis_used,pod_used,svd_used\n"
                << (passed
                    ? (allInterfaces
                        ? "stage_c_passed"
                        : "stage_b_configuration_passed")
                    : (allInterfaces
                        ? "stage_c_failed"
                        : "stage_b_configuration_failed"))
                << ',' << build.interfaces.size()
                << ",0,0,0,0,0\n";
            std::cout
                << "[Residual Krylov] "
                << (allInterfaces
                    ? "all-interface" : "representative-interface")
                << " basis-only build complete; transient skipped.\n";
            return;
        } else if (options.optimalPortRepresentativeInterfacePilot
                   || options.optimalPortMaximumInterfaceRefinementPilot
                   || options.optimalPortAllInterfaceBasis) {
            const bool runAllInterfaces =
                options.optimalPortAllInterfaceBasis;
            const bool runMaximumRefinement =
                options.optimalPortMaximumInterfaceRefinementPilot;
            const auto scalabilityStart = Clock::now();
            const std::vector<ScalabilityAuditRow> auditRows =
                readScalabilityAudit(
                    options.optimalPortTopologyAuditCsv,
                    physics.name);
            const std::vector<PortPatch> allPatches =
                buildOptimalPortPatches(
                    mesh, partition,
                    options.optimalPortOversamplingLayers);
            if (auditRows.size() != 25 || allPatches.size() != 25) {
                throw std::runtime_error(
                    "[Optimal port] RRAM26 scalability validation "
                    "requires exactly 25 physical interfaces.");
            }
            std::map<int, const PortPatch*> patchesById;
            for (const PortPatch& patch : allPatches) {
                patchesById.emplace(patch.interfaceId, &patch);
            }
            std::map<int, ScalabilityAuditRow> auditById;
            for (const ScalabilityAuditRow& audit : auditRows) {
                auditById.emplace(audit.interfaceId, audit);
                const auto found = patchesById.find(audit.interfaceId);
                if (found == patchesById.end()) {
                    throw std::runtime_error(
                        "[Optimal port] Topology-audit interface is "
                        "absent from the current owner map.");
                }
                const PortPatch& patch = *found->second;
                const bool hasLeft = std::find(
                    patch.patchSubdomains.begin(),
                    patch.patchSubdomains.end(),
                    patch.leftSubdomain)
                    != patch.patchSubdomains.end();
                const bool hasRight = std::find(
                    patch.patchSubdomains.begin(),
                    patch.patchSubdomains.end(),
                    patch.rightSubdomain)
                    != patch.patchSubdomains.end();
                std::set<int> targetUnique(
                    patch.target.begin(), patch.target.end());
                std::set<int> sourceUnique(
                    patch.source.begin(), patch.source.end());
                std::vector<int> overlap;
                std::set_intersection(
                    targetUnique.begin(), targetUnique.end(),
                    sourceUnique.begin(), sourceUnique.end(),
                    std::back_inserter(overlap));
                if (audit.leftSubdomain != patch.leftSubdomain
                    || audit.rightSubdomain != patch.rightSubdomain
                    || audit.targetDofs
                        != static_cast<int>(patch.target.size())
                    || audit.sourceDofs
                        != static_cast<int>(patch.source.size())
                    || !hasLeft || !hasRight
                    || targetUnique.size() != patch.target.size()
                    || sourceUnique.size() != patch.source.size()
                    || !overlap.empty()) {
                    throw std::runtime_error(
                        "[Optimal port] Topology/source/owner mapping "
                        "does not match the audited interface "
                        + std::to_string(audit.interfaceId) + '.');
                }
            }

            std::vector<ScalabilityAuditRow> orderedAudit = auditRows;
            std::sort(
                orderedAudit.begin(), orderedAudit.end(),
                [](const ScalabilityAuditRow& left,
                   const ScalabilityAuditRow& right) {
                    return left.targetDofs < right.targetDofs
                        || (left.targetDofs == right.targetDofs
                            && left.interfaceId < right.interfaceId);
                });
            std::vector<ScalabilityAuditRow> selectedAudit;
            std::vector<std::string> selectionLabels;
            std::vector<int> transferRanks;
            if (runAllInterfaces) {
                selectedAudit = auditRows;
                std::sort(
                    selectedAudit.begin(), selectedAudit.end(),
                    [](const ScalabilityAuditRow& left,
                       const ScalabilityAuditRow& right) {
                        return left.interfaceId < right.interfaceId;
                    });
                selectionLabels.assign(
                    selectedAudit.size(), "all");
                transferRanks = {4};
            } else if (runMaximumRefinement) {
                selectedAudit = {orderedAudit.back()};
                selectionLabels = {"maximum"};
                transferRanks = {4, 8};
            } else {
                selectedAudit = {
                    orderedAudit.front(),
                    orderedAudit[orderedAudit.size() / 2],
                    orderedAudit.back()};
                selectionLabels = {
                    "minimum", "median", "maximum"};
                transferRanks = {4, 8};
            }

            const auto writeDetailedHeader = [](std::ostream& output) {
                output
                    << "selection,interface_id,adjacent_subdomains,"
                    "target_dofs,source_dofs,mandatory_rank,"
                    "requested_transfer_rank,converged_transfer_rank,"
                    "total_port_rank,A_tt_nnz,woodbury_setup_time_s,"
                    "W_setup_time_s,Q_setup_time_s,"
                    "mean_target_solve_time_s,max_target_solve_time_s,"
                    "target_solve_residual,weighted_adjoint_error,"
                    "column_consistency_error,eigensolver_iterations,"
                    "eigensolver_residual,transfer_eigenvalues,"
                    "total_interface_basis_time_s,"
                    "incremental_workspace_bytes,"
                    "process_peak_memory_bytes,status,failure_reason,"
                    "inner_solver_actual,fallback_triggered,"
                    "fallback_reason\n";
            };
            std::ofstream representative(
                outputDirectory
                / "milestone8_rram26_representative_interface_pilot.csv");
            std::ofstream allInterface(
                outputDirectory
                / "milestone8_rram26_all_interface_basis.csv");
            std::ofstream summary(
                outputDirectory
                / "milestone8_rram26_all_interface_summary.csv");
            std::ofstream eigenvalues(
                outputDirectory
                / "milestone8_rram26_interface_eigenvalues.csv");
            std::ofstream timing(
                outputDirectory
                / "milestone8_rram26_basis_timing.csv");
            std::ofstream memory(
                outputDirectory
                / "milestone8_rram26_basis_memory.csv");
            std::ofstream failures(
                outputDirectory
                / "milestone8_rram26_basis_failures.csv");
            std::ofstream residualBreakdown(
                outputDirectory
                / "milestone8_woodbury_residual_breakdown.csv");
            if (!representative || !allInterface || !summary
                || !eigenvalues || !timing || !memory || !failures
                || !residualBreakdown) {
                throw std::runtime_error(
                    "[Optimal port] Cannot create scalability CSV outputs.");
            }
            writeDetailedHeader(representative);
            writeDetailedHeader(allInterface);
            eigenvalues
                << "scope,selection,interface_id,"
                "requested_transfer_rank,mode_index,eigenvalue,"
                "eigenpair_residual\n";
            timing
                << "scope,selection,interface_id,"
                "requested_transfer_rank,phase,time_s\n";
            memory
                << "scope,selection,interface_id,"
                "requested_transfer_rank,"
                "process_peak_memory_bytes,"
                "baseline_memory_before_port_build_bytes,"
                "incremental_port_build_memory_bytes,"
                "A_tt_factor_bytes,W_bytes,Q_bytes,"
                "eigensolver_workspace_bytes,"
                "final_port_basis_storage_bytes,"
                "serialized_model_bytes\n";
            failures
                << "scope,selection,interface_id,"
                "requested_transfer_rank,stage,status,failure_reason\n";
            residualBreakdown
                << "case,selection,interface_id,"
                "requested_transfer_rank,"
                "A_tt_solve_relative_residual,"
                "Q_solve_pre_refinement_residual,"
                "Q_solve_relative_residual,"
                "Q_refinement_iterations,"
                "woodbury_pre_refinement_residual,"
                "woodbury_post_refinement_residual,"
                "refinement_iterations,refinement_residual_0,"
                "refinement_residual_1,refinement_residual_2,"
                "refinement_residual_3,"
                "refinement_correction_relative_norm,"
                "refinement_reduction_factor,"
                "refinement_converged,"
                "refinement_triggered_solve_calls,"
                "pardiso_internal_refinement_steps,"
                "Q_min_abs_factor_diagonal,"
                "Q_max_abs_factor_diagonal,"
                "Q_factor_diagonal_ratio,"
                "woodbury_cancellation_factor,"
                "final_target_solve_residual,status,"
                "failure_reason\n";
            summary
                << "physical_interfaces,successful_interfaces,"
                "failed_interfaces,full_interface_dofs,"
                "mandatory_rank_total,requested_transfer_rank_total,"
                "converged_transfer_rank_total,total_port_rank,"
                "basis_build_total_time_s,basis_build_mean_time_s,"
                "basis_build_max_time_s,woodbury_setup_total_time_s,"
                "eigensolver_total_time_s,"
                "maximum_incremental_workspace_bytes,"
                "process_peak_memory_bytes,serialized_port_basis_bytes,"
                "snapshot_used,fom_used_for_basis,status\n";
            representative << std::setprecision(17);
            allInterface << std::setprecision(17);
            eigenvalues << std::setprecision(17);
            timing << std::setprecision(17);
            memory << std::setprecision(17);
            summary << std::setprecision(17);
            residualBreakdown << std::setprecision(17);

            const std::size_t baselineMemoryBytes =
                currentWorkingSetBytes();
            const auto sharedSchurStart = Clock::now();
            ReducedDynamicSchurOperator sharedScalabilitySchur(
                *preparedPortDynamicModel, false);
            const double sharedSchurSetupSeconds =
                elapsed(sharedSchurStart);
            const double meshModelPreparationSeconds = std::max(
                0.0, workflowBeforePortSeconds
                    - referenceSeconds - localBuildSeconds);
            timing
                << (runAllInterfaces ? "all-interface" : "representative")
                << ",all,-1,0,mesh_model_preparation,"
                << meshModelPreparationSeconds << '\n'
                << (runAllInterfaces ? "all-interface" : "representative")
                << ",all,-1,0,local_block_arnoldi_build,"
                << localBuildSeconds << '\n'
                << (runAllInterfaces ? "all-interface" : "representative")
                << ",all,-1,0,dynamic_reduced_schur_model_build,"
                << dynamicModelBuildSeconds << '\n'
                << (runAllInterfaces ? "all-interface" : "representative")
                << ",all,-1,0,reduced_schur_operator_setup,"
                << sharedSchurSetupSeconds << '\n';

            LocalPortModel aggregateModel;
            bool aggregateInitialized = false;
            bool stopped = false;
            int successfulInterfaces = 0;
            int failedInterfaces = 0;
            int mandatoryRankTotal = 0;
            int requestedTransferRankTotal = 0;
            int convergedTransferRankTotal = 0;
            int totalPortRank = 0;
            double basisBuildTotalSeconds = 0.0;
            double basisBuildMaximumSeconds = 0.0;
            double woodburySetupTotalSeconds = 0.0;
            double eigenSolveTotalSeconds = 0.0;
            std::size_t maximumIncrementalBytes = 0;
            std::size_t serializedBytes = 0;

            for (std::size_t selectedIndex = 0;
                 selectedIndex < selectedAudit.size() && !stopped;
                 ++selectedIndex) {
                const ScalabilityAuditRow& audit =
                    selectedAudit[selectedIndex];
                const std::string& label =
                    selectionLabels[selectedIndex];
                for (int requestedTransferRank : transferRanks) {
                    OptimalTransferPortOptions trialOptions;
                    trialOptions.requestedRank =
                        4 + requestedTransferRank;
                    trialOptions.requestedTransferRank =
                        requestedTransferRank;
                    trialOptions.rankMode = "fixed";
                    trialOptions.eigenvalueTolerance =
                        options.optimalPortEigenvalueTolerance;
                    trialOptions.minimumRank = 1;
                    trialOptions.maximumRank = 12;
                    trialOptions.innerProduct =
                        options.optimalPortInnerProduct;
                    trialOptions.oversamplingLayers =
                        options.optimalPortOversamplingLayers;
                    trialOptions.eigensolverMaximumIterations =
                        std::min(
                            500,
                            options.optimalPortEigenMaximumIterations);
                    trialOptions.eigensolverTolerance =
                        options.optimalPortEigenTolerance;
                    trialOptions.relativeDeflationTolerance =
                        options.rankTolerance;
                    trialOptions.ablationMode =
                        "constant-geometry-trace";
                    trialOptions.sourceMode = "trace-only";
                    trialOptions.selectedInterfaceIds = {
                        audit.interfaceId};
                    trialOptions.targetSolverPilotPreflight = true;
                    trialOptions.columnConsistencyCheck = true;
                    trialOptions.maximumPilotSeconds = 600.0;
                    trialOptions.maximumTargetSetupSeconds = 30.0;
                    trialOptions.maximumMeanTargetSolveSeconds = 0.1;
                    if (runMaximumRefinement) {
                        // M8.5b-A1 qualifies on the unchanged residual,
                        // eigensolver, total-time, and memory gates. Mean
                        // solve time is recorded but is not an A1 gate.
                        trialOptions.maximumMeanTargetSolveSeconds =
                            std::numeric_limits<double>::infinity();
                    }
                    trialOptions.maximumOperatorCheckSeconds = 60.0;
                    trialOptions.maximumTargetResidual = 1.0e-9;
                    trialOptions.maximumWeightedAdjointError = 1.0e-8;
                    trialOptions.maximumIncrementalWorkspaceBytes =
                        UINT64_C(1073741824);
                    trialOptions.innerSolver.innerSolver =
                        "woodbury-exact";
                    trialOptions.innerSolver.relativeTolerance =
                        options.optimalPortInnerTolerance;
                    trialOptions.innerSolver.maximumIterations =
                        options.optimalPortInnerMaximumIterations;
                    trialOptions.innerSolver.refinementMaximumIterations =
                        options.optimalPortInnerRefinementMaximumIterations;
                    trialOptions.innerSolver.refinementTolerance =
                        options.optimalPortInnerRefinementTolerance;

                    OptimalPortBuildResult trial;
                    std::string exceptionReason;
                    try {
                        trial = buildOptimalTransferPortModel(
                            mesh, partition,
                            *preparedPortDynamicModel,
                            descriptor.interfaceTraceMassDiagonal,
                            descriptor.interfacePenaltyMassDiagonal,
                            descriptor.input,
                            descriptor.sourceChannels,
                            descriptor.boundaryRhs,
                            reducedHistory.condensedColumns,
                            reducedHistory.channels,
                            trialOptions,
                            &sharedScalabilitySchur);
                    } catch (const std::exception& error) {
                        exceptionReason = error.what();
                    }
                    if (!exceptionReason.empty()
                        || trial.interfaces.size() != 1
                        || trial.model.ports.size() != 1) {
                        const std::string reason =
                            !exceptionReason.empty()
                            ? exceptionReason
                            : "basis_build_did_not_return_one_interface";
                        std::ostream& detail = runAllInterfaces
                            ? static_cast<std::ostream&>(allInterface)
                            : static_cast<std::ostream&>(representative);
                        detail << label << ',' << audit.interfaceId
                            << ",\"" << audit.leftSubdomain << ';'
                            << audit.rightSubdomain << "\","
                            << audit.targetDofs << ','
                            << audit.sourceDofs
                            << ",0," << requestedTransferRank
                            << ",0,0,0,0,0,0,0,0,0,0,0,0,0,\"\",0,0,"
                            << peakWorkingSetBytes()
                            << ",failed,";
                        writeCsvString(detail, reason);
                        detail << ",\"\",0,\"\"\n";
                        failures
                            << (runAllInterfaces
                                ? "all-interface" : "representative")
                            << ',' << label << ','
                            << audit.interfaceId << ','
                            << requestedTransferRank
                            << ",basis_build_exception,failed,";
                        writeCsvString(failures, reason);
                        failures << '\n';
                        residualBreakdown
                            << physics.name << ',' << label << ','
                            << audit.interfaceId << ','
                            << requestedTransferRank
                            << ",0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,"
                            << "failed,";
                        writeCsvString(
                            residualBreakdown, reason);
                        residualBreakdown << '\n';
                        ++failedInterfaces;
                        stopped = true;
                        break;
                    }

                    const OptimalPortInterfaceDiagnostics& diagnostics =
                        trial.interfaces.front();
                    const LocalPortBasis& port =
                        trial.model.ports.front();
                    const PatchInnerSolverStatistics& solver =
                        diagnostics.innerSolver;
                    const double meanSolve = solver.solveCalls > 0
                        ? solver.totalSolveSeconds / solver.solveCalls
                        : 0.0;
                    const double qSetupSeconds =
                        solver.qAssemblySeconds
                        + solver.qFactorizationSeconds;
                    std::string status = "passed";
                    std::string failureReason;
                    std::string failureStage;
                    const auto fail = [&](const std::string& stage,
                                          const std::string& reason) {
                        if (status == "passed") {
                            status = "failed";
                            failureStage = stage;
                            failureReason = reason;
                        }
                    };
                    if (diagnostics.interfaceId != audit.interfaceId
                        || diagnostics.targetRows != audit.targetDofs
                        || diagnostics.traceSourceRows
                            != audit.sourceDofs) {
                        fail("topology_owner_mapping",
                             "runtime_mapping_differs_from_audit");
                    }
                    if (diagnostics.pilotStatus
                            != "passed_preflight") {
                        if (runMaximumRefinement
                            && diagnostics.pilotStatus
                                == "failed_target_residual") {
                            fail("target_solve",
                                 "woodbury_target_solve_accuracy_gate_failed");
                        } else {
                            fail("target_solver_preflight",
                                 diagnostics.pilotStatus);
                        }
                    }
                    if (solver.setupSeconds > 30.0) {
                        fail("woodbury_setup",
                             "woodbury_setup_exceeded_30_s");
                    }
                    if (!runMaximumRefinement && meanSolve > 0.1) {
                        fail("target_solve",
                             "mean_target_solve_exceeded_0.1_s");
                    }
                    if (solver.maximumRelativeResidual > 1.0e-9) {
                        fail("target_solve",
                             runMaximumRefinement
                                ? "woodbury_target_solve_accuracy_gate_failed"
                                : "target_residual_exceeded_1e-9");
                    }
                    if (diagnostics.adjointRelativeError > 1.0e-8) {
                        fail("weighted_adjoint",
                             "weighted_adjoint_error_exceeded_1e-8");
                    }
                    if (diagnostics.explicitColumnReferenceError
                            > 1.0e-10) {
                        fail("column_consistency",
                             "column_consistency_error_exceeded_1e-10");
                    }
                    if (!diagnostics.eigenConverged
                        || diagnostics.convergedTransferRank
                            != requestedTransferRank) {
                        fail("transfer_eigensolver",
                             "requested_transfer_rank_not_converged");
                    }
                    if (diagnostics.maximumEigenpairResidual > 1.0e-8) {
                        fail("transfer_eigensolver",
                             "eigensolver_residual_exceeded_1e-8");
                    }
                    if (diagnostics.eigenIterations > 500) {
                        fail("transfer_eigensolver",
                             "eigensolver_iterations_exceeded_500");
                    }
                    if (trial.totalSeconds > 600.0) {
                        fail("interface_basis",
                             "single_interface_time_exceeded_600_s");
                    }
                    if (diagnostics.peakWorkspaceBytes
                            > UINT64_C(1073741824)) {
                        fail("memory",
                             "incremental_workspace_exceeded_1_gib");
                    }

                    std::ostream& detail = runAllInterfaces
                        ? static_cast<std::ostream&>(allInterface)
                        : static_cast<std::ostream&>(representative);
                    detail << label << ',' << diagnostics.interfaceId
                        << ",\"" << audit.leftSubdomain << ';'
                        << audit.rightSubdomain << "\","
                        << diagnostics.targetRows << ','
                        << diagnostics.traceSourceRows << ','
                        << diagnostics.mandatoryRank << ','
                        << diagnostics.requestedTransferRank << ','
                        << diagnostics.convergedTransferRank << ','
                        << diagnostics.totalPortRank << ','
                        << solver.aTtNonzeros << ','
                        << solver.setupSeconds << ','
                        << solver.wSetupSeconds << ','
                        << qSetupSeconds << ','
                        << meanSolve << ','
                        << solver.maximumSolveSeconds << ','
                        << solver.maximumRelativeResidual << ','
                        << diagnostics.adjointRelativeError << ','
                        << diagnostics.explicitColumnReferenceError << ','
                        << diagnostics.eigenIterations << ','
                        << diagnostics.maximumEigenpairResidual << ",\"";
                    for (std::size_t eigen = 0;
                         eigen < diagnostics.eigenvalues.size(); ++eigen) {
                        if (eigen != 0) detail << ';';
                        detail << diagnostics.eigenvalues[eigen];
                    }
                    detail << "\"," << trial.totalSeconds << ','
                        << diagnostics.peakWorkspaceBytes << ','
                        << peakWorkingSetBytes() << ','
                        << status << ',';
                    writeCsvString(detail, failureReason);
                    detail << ',';
                    writeCsvString(detail, solver.actualSolver);
                    detail << ','
                        << (solver.fallbackReason.empty() ? 0 : 1)
                        << ',';
                    writeCsvString(detail, solver.fallbackReason);
                    detail << '\n';
                    detail.flush();

                    residualBreakdown
                        << physics.name << ',' << label << ','
                        << diagnostics.interfaceId << ','
                        << requestedTransferRank << ','
                        << solver.aTtSolveRelativeResidual << ','
                        << solver.qSolvePreRefinementResidual << ','
                        << solver.qSolveRelativeResidual << ','
                        << solver.qRefinementIterations << ','
                        << solver.woodburyPreRefinementResidual << ','
                        << solver.woodburyPostRefinementResidual << ','
                        << solver.refinementIterations << ','
                        << solver.refinementResidual0 << ','
                        << solver.refinementResidual1 << ','
                        << solver.refinementResidual2 << ','
                        << solver.refinementResidual3 << ','
                        << solver.refinementCorrectionRelativeNorm << ','
                        << solver.refinementReductionFactor << ','
                        << (solver.refinementConverged ? 1 : 0) << ','
                        << solver.refinementTriggeredSolveCalls << ','
                        << solver.pardisoInternalRefinementSteps << ','
                        << solver.qMinimumAbsoluteFactorDiagonal << ','
                        << solver.qMaximumAbsoluteFactorDiagonal << ','
                        << solver.qFactorDiagonalRatio << ','
                        << solver.woodburyCancellationFactor << ','
                        << solver.maximumRelativeResidual << ','
                        << status << ',';
                    writeCsvString(
                        residualBreakdown, failureReason);
                    residualBreakdown << '\n';
                    residualBreakdown.flush();

                    for (std::size_t eigen = 0;
                         eigen < diagnostics.eigenvalues.size(); ++eigen) {
                        eigenvalues
                            << (runAllInterfaces
                                ? "all-interface" : "representative")
                            << ',' << label << ','
                            << diagnostics.interfaceId << ','
                            << requestedTransferRank << ','
                            << eigen << ','
                            << diagnostics.eigenvalues[eigen] << ','
                            << (eigen < port.spectralResiduals.size()
                                ? port.spectralResiduals[eigen] : 0.0)
                            << '\n';
                    }
                    const std::string scope = runAllInterfaces
                        ? "all-interface" : "representative";
                    timing << scope << ',' << label << ','
                        << diagnostics.interfaceId << ','
                        << requestedTransferRank
                        << ",topology_owner_mapping,"
                        << trial.patchSetupSeconds << '\n'
                        << scope << ',' << label << ','
                        << diagnostics.interfaceId << ','
                        << requestedTransferRank
                        << ",target_solver_setup,"
                        << solver.setupSeconds << '\n'
                        << scope << ',' << label << ','
                        << diagnostics.interfaceId << ','
                        << requestedTransferRank << ",eigensolver,"
                        << trial.eigenSolveSeconds << '\n'
                        << scope << ',' << label << ','
                        << diagnostics.interfaceId << ','
                        << requestedTransferRank
                        << ",orthogonalization,"
                        << trial.orthogonalizationSeconds << '\n'
                        << scope << ',' << label << ','
                        << diagnostics.interfaceId << ','
                        << requestedTransferRank << ",serialization,0\n"
                        << scope << ',' << label << ','
                        << diagnostics.interfaceId << ','
                        << requestedTransferRank
                        << ",pure_port_basis_build,"
                        << trial.totalSeconds << '\n';
                    memory << scope << ',' << label << ','
                        << diagnostics.interfaceId << ','
                        << requestedTransferRank << ','
                        << peakWorkingSetBytes() << ','
                        << baselineMemoryBytes << ','
                        << diagnostics.peakWorkspaceBytes << ','
                        << solver.aTtFactorBytes << ','
                        << solver.wBytes << ','
                        << solver.qBytes << ','
                        << diagnostics.eigensolverWorkspaceBytes << ','
                        << diagnostics.finalBasisStorageBytes
                        << ",0\n";
                    timing.flush();
                    memory.flush();
                    eigenvalues.flush();

                    basisBuildTotalSeconds += trial.totalSeconds;
                    basisBuildMaximumSeconds = std::max(
                        basisBuildMaximumSeconds,
                        trial.totalSeconds);
                    woodburySetupTotalSeconds +=
                        solver.setupSeconds;
                    eigenSolveTotalSeconds +=
                        trial.eigenSolveSeconds;
                    maximumIncrementalBytes = std::max(
                        maximumIncrementalBytes,
                        diagnostics.peakWorkspaceBytes);
                    if (status != "passed") {
                        failures << scope << ',' << label << ','
                            << diagnostics.interfaceId << ','
                            << requestedTransferRank << ','
                            << failureStage << ",failed,";
                        writeCsvString(failures, failureReason);
                        failures << '\n';
                        failures.flush();
                        ++failedInterfaces;
                        stopped = true;
                        break;
                    }
                    ++successfulInterfaces;

                    if (runAllInterfaces) {
                        if (!aggregateInitialized) {
                            aggregateModel =
                                std::move(trial.model);
                            aggregateInitialized = true;
                        } else {
                            LocalPortBasis retained =
                                std::move(trial.model.ports.front());
                            aggregateModel.reducedInterfaceDofs +=
                                retained.rank;
                            aggregateModel.modelBytes +=
                                diagnostics.finalBasisStorageBytes;
                            aggregateModel.ports.push_back(
                                std::move(retained));
                        }
                        mandatoryRankTotal +=
                            diagnostics.mandatoryRank;
                        requestedTransferRankTotal +=
                            diagnostics.requestedTransferRank;
                        convergedTransferRankTotal +=
                            diagnostics.convergedTransferRank;
                        totalPortRank += diagnostics.totalPortRank;
                    }
                    std::cout
                        << "[Optimal port scalability] selection="
                        << label << ", interface="
                        << diagnostics.interfaceId << ", transfer_rank="
                        << requestedTransferRank << ", target="
                        << diagnostics.targetRows << ", status="
                        << status << ", time=" << trial.totalSeconds
                        << " s\n";
                }
            }

            double serializationSeconds = 0.0;
            if (runAllInterfaces && !stopped
                && successfulInterfaces == 25
                && aggregateInitialized) {
                std::uint64_t traceFingerprint =
                    UINT64_C(1469598103934665603);
                for (const PortPatch& patch : allPatches) {
                    hashValue(traceFingerprint, patch.interfaceId);
                    hashValue(
                        traceFingerprint, patch.sourceFingerprint);
                }
                aggregateModel.traceSourceFingerprint =
                    traceFingerprint;
                aggregateModel.timeStep = options.timeStep;
                aggregateModel.meshFingerprint =
                    descriptor.fingerprints.mesh;
                aggregateModel.materialFingerprint =
                    descriptor.fingerprints.conductivity;
                aggregateModel.operatorFingerprint =
                    dynamicOperatorFingerprint;
                aggregateModel.inputFingerprint =
                    descriptor.fingerprints.input;
                aggregateModel.penaltyFingerprint =
                    penaltyFingerprint;
                aggregateModel.boundaryFingerprint =
                    descriptor.fingerprints.boundary;
                aggregateModel.sourceMetadataFingerprint =
                    descriptor.fingerprints.sources;
                aggregateModel.rankFileFingerprint =
                    rankFileFingerprint;
                aggregateModel.requestedRank = 8;
                aggregateModel.minimumRank = 1;
                aggregateModel.maximumRank = 12;
                aggregateModel.eigensolverMaximumIterations =
                    std::min(
                        500,
                        options.optimalPortEigenMaximumIterations);
                aggregateModel.innerSolver = "woodbury-exact";
                aggregateModel.basisSeconds =
                    basisBuildTotalSeconds;
                aggregateModel.snapshotSeconds = 0.0;
                aggregateModel.buildTimestamp = std::to_string(
                    std::chrono::duration_cast<std::chrono::seconds>(
                        std::chrono::system_clock::now()
                            .time_since_epoch()).count());
                aggregateModel.sourceCommit =
                    DDM_SCHUR_SOURCE_COMMIT;
                const std::filesystem::path basisPath =
                    outputDirectory
                    / "milestone8_rram26_rank4_port_basis.bin";
                const auto serializationStart = Clock::now();
                saveLocalPortModel(aggregateModel, basisPath);
                serializationSeconds =
                    elapsed(serializationStart);
                serializedBytes =
                    static_cast<std::size_t>(
                        std::filesystem::file_size(basisPath));
                timing << "all-interface,all,-1,4,serialization,"
                    << serializationSeconds << '\n';
                memory << "all-interface,all,-1,4,"
                    << peakWorkingSetBytes() << ','
                    << baselineMemoryBytes << ','
                    << maximumIncrementalBytes
                    << ",0,0,0,0,"
                    << aggregateModel.modelBytes << ','
                    << serializedBytes << '\n';
            }

            const int basisCases = successfulInterfaces
                + failedInterfaces;
            const std::string allStatus = runAllInterfaces
                ? (!stopped && successfulInterfaces == 25
                    ? "passed" : "failed_stopped_at_first_interface")
                : "not_run_stage_a_stop";
            summary << 25 << ','
                << (runAllInterfaces ? successfulInterfaces : 0)
                << ',' << (runAllInterfaces ? failedInterfaces : 0)
                << ',' << preparedPortDynamicModel->interfaceDofs
                << ',' << (runAllInterfaces
                    ? mandatoryRankTotal : 0)
                << ',' << (runAllInterfaces
                    ? requestedTransferRankTotal : 0)
                << ',' << (runAllInterfaces
                    ? convergedTransferRankTotal : 0)
                << ',' << (runAllInterfaces ? totalPortRank : 0)
                << ',' << (runAllInterfaces
                    ? basisBuildTotalSeconds : 0.0)
                << ',' << (runAllInterfaces && basisCases > 0
                    ? basisBuildTotalSeconds / basisCases : 0.0)
                << ',' << (runAllInterfaces
                    ? basisBuildMaximumSeconds : 0.0)
                << ',' << (runAllInterfaces
                    ? woodburySetupTotalSeconds : 0.0)
                << ',' << (runAllInterfaces
                    ? eigenSolveTotalSeconds : 0.0)
                << ',' << (runAllInterfaces
                    ? maximumIncrementalBytes : 0)
                << ',' << peakWorkingSetBytes()
                << ',' << serializedBytes
                << ",0,0," << allStatus << '\n';
            timing
                << (runAllInterfaces ? "all-interface" : "representative")
                << ",all,-1,0,workflow_wall_before_return,"
                << elapsed(totalStart) << '\n'
                << (runAllInterfaces ? "all-interface" : "representative")
                << ",all,-1,0,scalability_stage_wall,"
                << elapsed(scalabilityStart) << '\n';
            timing.flush();
            memory.flush();
            summary.flush();
            failures.flush();

            if (runAllInterfaces) {
                std::cout
                    << "[Optimal port scalability] all-interface rank-4 "
                    << allStatus
                    << "; transient intentionally skipped.\n";
            } else {
                std::cout
                    << "[Optimal port scalability] representative Stage A "
                    << (stopped ? "failed" : "passed")
                    << "; all-interface basis and transient intentionally "
                    << "skipped.\n";
            }
            return;
        } else if (options.optimalPortBasisPilot) {
            const std::vector<PortPatch> allPatches =
                buildOptimalPortPatches(
                    mesh, partition,
                    options.optimalPortOversamplingLayers);
            if (allPatches.size() < 3) {
                throw std::runtime_error(
                    "[Optimal port] Formal basis pilot needs at least three interfaces.");
            }
            std::vector<const PortPatch*> ordered;
            ordered.reserve(allPatches.size());
            for (const PortPatch& patch : allPatches) {
                ordered.push_back(&patch);
            }
            std::sort(
                ordered.begin(), ordered.end(),
                [](const PortPatch* left, const PortPatch* right) {
                    return left->target.size() < right->target.size()
                        || (left->target.size() == right->target.size()
                            && left->interfaceId < right->interfaceId);
                });
            const std::vector<const PortPatch*> selected{
                ordered.front(),
                ordered[ordered.size() / 2],
                ordered.back()};
            const std::vector<int> budgets{8, 12, 16};
            std::ofstream pilot(
                outputDirectory
                / "milestone8_rram26_basis_pilot.csv");
            pilot << "selection,total_rank_budget,interface_id,"
                "target_dofs,source_dofs,mandatory_rank,"
                "requested_transfer_rank,converged_transfer_rank,"
                "total_port_rank,eigensolver_iterations,"
                "eigensolver_status,eigensolver_converged,"
                "eigenpair_residual,weighted_adjoint_error,"
                "inner_solver_requested,inner_solver_actual,"
                "fallback_triggered,fallback_reason,setup_time_s,"
                "shared_schur_setup_time_s,peak_workspace_bytes,"
                "preconditioner_diagonal,eigenvalues\n"
                << std::setprecision(17);
            const std::vector<std::string> labels{
                "minimum", "median", "maximum"};
            // The dynamic matrix, local factors, and selected Schur-Jacobi
            // diagonal do not change with interface selection or rank budget.
            // Construct them once for the full pilot.
            const auto pilotStart = Clock::now();
            const auto sharedSchurStart = Clock::now();
            ReducedDynamicSchurOperator sharedPilotSchur(
                *preparedPortDynamicModel,
                false);
            const double sharedSchurSetupSeconds =
                elapsed(sharedSchurStart);
            bool allConverged = true;
            for (std::size_t selectedIndex = 0;
                 selectedIndex < selected.size(); ++selectedIndex) {
                for (int budget : budgets) {
                    OptimalTransferPortOptions pilotOptions;
                    pilotOptions.requestedRank = budget;
                    pilotOptions.rankMode = "fixed";
                    pilotOptions.eigenvalueTolerance =
                        options.optimalPortEigenvalueTolerance;
                    pilotOptions.minimumRank = 1;
                    pilotOptions.maximumRank = 16;
                    pilotOptions.innerProduct =
                        options.optimalPortInnerProduct;
                    pilotOptions.oversamplingLayers =
                        options.optimalPortOversamplingLayers;
                    pilotOptions.eigensolverMaximumIterations =
                        options.optimalPortEigenMaximumIterations;
                    pilotOptions.eigensolverTolerance =
                        options.optimalPortEigenTolerance;
                    pilotOptions.relativeDeflationTolerance =
                        options.rankTolerance;
                    pilotOptions.ablationMode =
                        options.optimalPortAblation;
                    pilotOptions.sourceMode =
                        effectiveOptimalPortSourceMode;
                    pilotOptions.selectedInterfaceIds = {
                        selected[selectedIndex]->interfaceId};
                    pilotOptions.innerSolver.innerSolver =
                        options.optimalPortInnerSolver;
                    pilotOptions.innerSolver.relativeTolerance =
                        options.optimalPortInnerTolerance;
                    pilotOptions.innerSolver.maximumIterations =
                        options.optimalPortInnerMaximumIterations;
                    pilotOptions.innerSolver.refinementMaximumIterations =
                        options.optimalPortInnerRefinementMaximumIterations;
                    pilotOptions.innerSolver.refinementTolerance =
                        options.optimalPortInnerRefinementTolerance;
                    OptimalPortBuildResult trial =
                        buildOptimalTransferPortModel(
                            mesh, partition,
                            *preparedPortDynamicModel,
                            descriptor.interfaceTraceMassDiagonal,
                            descriptor.interfacePenaltyMassDiagonal,
                            descriptor.input,
                            descriptor.sourceChannels,
                            descriptor.boundaryRhs,
                            reducedHistory.condensedColumns,
                            reducedHistory.channels,
                            pilotOptions,
                            &sharedPilotSchur);
                    if (trial.interfaces.size() != 1
                        || trial.model.ports.size() != 1) {
                        throw std::runtime_error(
                            "[Optimal port] Pilot did not return one interface.");
                    }
                    const OptimalPortInterfaceDiagnostics& row =
                        trial.interfaces.front();
                    pilot << labels[selectedIndex] << ',' << budget
                        << ',' << row.interfaceId
                        << ',' << row.targetRows
                        << ',' << row.sourceRows
                        << ',' << row.mandatoryRank
                        << ',' << row.requestedTransferRank
                        << ',' << row.convergedTransferRank
                        << ',' << row.totalPortRank
                        << ',' << row.eigenIterations
                        << ',' << row.eigenStatus
                        << ',' << (row.eigenConverged ? 1 : 0)
                        << ',' << row.maximumEigenpairResidual
                        << ',' << row.adjointRelativeError
                        << ',' << row.innerSolver.requestedSolver
                        << ',' << row.innerSolver.actualSolver
                        << ','
                        << (row.innerSolver.fallbackReason.empty()
                            ? 0 : 1)
                        << ',' << row.innerSolver.fallbackReason
                        << ',' << trial.totalSeconds
                        << ',' << sharedSchurSetupSeconds
                        << ',' << row.peakWorkspaceBytes
                        << ",assembled-interface,";
                    for (std::size_t eigen = 0;
                         eigen < row.eigenvalues.size(); ++eigen) {
                        if (eigen != 0) pilot << ';';
                        pilot << row.eigenvalues[eigen];
                    }
                    pilot << '\n';
                    pilot.flush();
                    allConverged =
                        allConverged && row.eigenConverged;
                    std::cout
                        << "[Optimal port pilot] selection="
                        << labels[selectedIndex]
                        << ", interface=" << row.interfaceId
                        << ", budget=" << budget
                        << ", mandatory=" << row.mandatoryRank
                        << ", transfer="
                        << row.convergedTransferRank
                        << ", converged="
                        << (row.eigenConverged ? 1 : 0)
                        << ", setup=" << trial.totalSeconds
                        << " s\n";
                }
            }
            std::ofstream stop(
                outputDirectory
                / "milestone8_rram26_basis_pilot_stop.csv");
            stop << "status,interfaces_total,interfaces_piloted,"
                "rank_budgets,transient_advanced,rank24_used,"
                "full_field_read,snapshot_basis_used,"
                "shared_schur_setup_time_s,pilot_wall_time_s,"
                "preconditioner_diagonal,peak_working_set_bytes\n"
                << (allConverged
                    ? "pilot_passed" : "pilot_has_nonconverged_cases")
                << ',' << allPatches.size()
                << ",3,\"8;12;16\",0,0,0,0,"
                << sharedSchurSetupSeconds << ','
                << elapsed(pilotStart)
                << ",assembled-interface,"
                << peakWorkingSetBytes() << '\n';
            std::cout
                << "[Optimal port pilot] complete; transient and all-interface "
                << "basis build intentionally skipped.\n";
            return;
        } else if (!loadPath.empty() && !options.fluxAwarePort) {
            portModel = std::make_unique<LocalPortModel>(
                loadLocalPortModel(loadPath));
            if (portModel->basisMethod != options.portBasisMethod) {
                throw std::runtime_error(
                    "[Optimal port] Loaded basis method does not match the requested method.");
            }
            if (portModel->interfaceGlobalDofs
                != partition.interfaceGlobalDofs) {
                throw std::runtime_error(
                    "[Optimal port] Loaded basis does not match the interface ordering.");
            }
            if ((portModel->basisMethod == "optimal-transfer"
                    || portModel->basisMethod
                        == "randomized-transfer"
                    || portModel->basisMethod
                        == "hybrid-randomized")
                && portModel->formatVersion < 6) {
                throw std::runtime_error(
                    "[Optimal port] Cached basis predates generalized-source provenance; rebuild it.");
            }
            if (portModel->formatVersion < 4) {
                throw std::runtime_error(
                    "[Optimal port] Cached M8 basis lacks complete version-4 provenance; rebuild it.");
            }
            if (operatorInformedPort || hybridRandomized) {
                if (options.historyCompressionMethod != "none"
                    && portModel->formatVersion < 7) {
                    throw std::runtime_error(
                        "[History compression] Cached basis predates "
                        "version-7 compression provenance; rebuild it.");
                }
                if (portModel->formatVersion >= 7
                    && (portModel->historyCompressionMethod
                            != options.historyCompressionMethod
                        || portModel->historyCompressionRank
                            != options.historyCompressionRank
                        || portModel->historyCompressionTolerance
                            != options.historyCompressionTolerance)) {
                    throw std::runtime_error(
                        "[History compression] Loaded basis compression "
                        "configuration does not match the request.");
                }
            }
            std::vector<std::string> provenanceMismatches;
            auto requireProvenance =
                [&](bool matches, const char* field) {
                    if (!matches) {
                        provenanceMismatches.emplace_back(field);
                    }
                };
            requireProvenance(
                portModel->timeStep == options.timeStep, "time_step");
            requireProvenance(
                portModel->meshFingerprint
                    == descriptor.fingerprints.mesh,
                "mesh");
            requireProvenance(
                portModel->materialFingerprint
                    == descriptor.fingerprints.conductivity,
                "material");
            requireProvenance(
                portModel->operatorFingerprint
                    == dynamicOperatorFingerprint,
                "dynamic_operator");
            requireProvenance(
                portModel->inputFingerprint
                    == descriptor.fingerprints.input,
                "input");
            requireProvenance(
                portModel->penaltyFingerprint == penaltyFingerprint,
                "penalty");
            requireProvenance(
                portModel->boundaryFingerprint
                    == descriptor.fingerprints.boundary,
                "boundary");
            requireProvenance(
                portModel->sourceMetadataFingerprint
                    == descriptor.fingerprints.sources,
                "source_metadata");
            if (portModel->formatVersion >= 6) {
                requireProvenance(
                    portModel->generalizedInputFingerprint
                        == fingerprintVector(descriptor.input),
                    "generalized_input");
                requireProvenance(
                    portModel->generalizedBoundaryFingerprint
                        == fingerprintVector(descriptor.boundaryRhs),
                    "generalized_boundary");
                requireProvenance(
                    portModel->generalizedHistoryFingerprint
                        == generalizedHistoryFingerprint,
                    "generalized_history");
            }
            requireProvenance(
                portModel->rankFileFingerprint
                    == rankFileFingerprint,
                "rank_file");
            const bool sourceOnlyGlobalPrototype =
                options.globalInterfaceCoarsePrototype
                || options.globalRandomizedSchur
                || options.projectionDiagnosis;
            requireProvenance(
                sourceOnlyGlobalPrototype
                    || portModel->sourceCommit
                        == DDM_SCHUR_SOURCE_COMMIT,
                "source_commit");
            if (sourceOnlyGlobalPrototype
                && portModel->sourceCommit
                    != DDM_SCHUR_SOURCE_COMMIT) {
                std::cout
                    << "[Global interface prototype] Reusing the frozen M8.9 local "
                    << "basis across a source-only prototype commit; "
                    << "all mesh/operator/source fingerprints passed.\n";
            }
            if (!provenanceMismatches.empty()) {
                std::ostringstream message;
                for (std::size_t index = 0;
                     index < provenanceMismatches.size(); ++index) {
                    if (index != 0) message << ';';
                    message << provenanceMismatches[index];
                }
                throw std::runtime_error(
                    "[Optimal port] Loaded basis provenance does not "
                    "match the current operator: " + message.str());
            }
            const std::string expectedInnerProduct = steklovSchur
                ? "penalty-weighted-mass"
                : options.optimalPortInnerProduct;
            const std::string expectedRankMode = steklovSchur
                ? "fixed"
                : (hybridRandomized
                    ? "hybrid-randomized-residual"
                    : (operatorInformedPort
                    ? "residual-tolerance"
                    : (randomizedTransfer
                        ? "fixed" : options.optimalPortRankMode)));
            const int expectedRank = steklovSchur
                ? std::max(1, options.localPortRank)
                : (hybridRandomized
                    ? options.randomizedPortRank
                    : (operatorInformedPort
                    ? (mandatoryOnly
                        ? 0 : options.residualKrylovMaximumRank)
                    : (randomizedTransfer
                        ? options.randomizedPortRank
                        : options.optimalPortRank)));
            const int expectedMinimumRank = steklovSchur
                ? expectedRank
                : (hybridRandomized
                    ? 0
                    : (operatorInformedPort
                    ? 0 : (randomizedTransfer
                        ? expectedRank
                        : options.optimalPortMinimumRank)));
            const int expectedMaximumRank = steklovSchur
                ? expectedRank
                : (hybridRandomized
                    ? options.residualKrylovMaximumRank
                    : (operatorInformedPort
                    ? expectedRank : (randomizedTransfer
                        ? expectedRank
                        : options.optimalPortMaximumRank)));
            const double expectedEigenvalueTolerance = steklovSchur
                ? 0.0
                : ((operatorInformedPort || hybridRandomized)
                    ? options.residualKrylovTolerance
                    : (randomizedTransfer
                        ? 0.0
                        : options.optimalPortEigenvalueTolerance));
            const double expectedEigensolverTolerance = steklovSchur
                ? options.rankTolerance
                : ((operatorInformedPort || hybridRandomized)
                    ? 0.0 : (randomizedTransfer
                        ? 0.0 : options.optimalPortEigenTolerance));
            const int expectedEigensolverIterations = steklovSchur
                ? 24
                : ((operatorInformedPort || hybridRandomized)
                    ? options.residualKrylovMaximumSweeps
                    : (randomizedTransfer
                        ? options.randomizedPortPowerIterations
                        : options.optimalPortEigenMaximumIterations));
            const std::string expectedInnerSolver = steklovSchur
                ? "direct"
                : ((operatorInformedPort || hybridRandomized)
                    ? options.residualKrylovInnerSolver
                    : (randomizedTransfer
                        ? "woodbury-exact"
                        : options.optimalPortInnerSolver));
            const double expectedInnerTolerance = steklovSchur
                ? 0.0 : options.optimalPortInnerTolerance;
            const int expectedInnerIterations = steklovSchur
                ? 0 : options.optimalPortInnerMaximumIterations;
            if (portModel->innerProduct != expectedInnerProduct
                || portModel->rankMode != expectedRankMode
                || portModel->ablationMode
                    != (steklovSchur
                        ? "steklov-schur"
                        : (hybridRandomized
                            ? "hybrid-randomized"
                            : (operatorInformedPort
                            ? options.portBasisMethod
                            : (randomizedTransfer
                                ? "randomized-range"
                                : options.optimalPortAblation))))
                || ((optimalTransfer || randomizedTransfer
                        || hybridRandomized || operatorInformedPort)
                    && portModel->sourceMode
                        != (hybridRandomized
                            ? effectiveOptimalPortSourceMode
                            : (operatorInformedPort
                            ? options.residualKrylovProbeMode
                            : effectiveOptimalPortSourceMode)))
                || portModel->oversamplingLayers
                    != (steklovSchur
                        ? 0 : options.optimalPortOversamplingLayers)
                || portModel->requestedRank != expectedRank
                || portModel->minimumRank != expectedMinimumRank
                || portModel->maximumRank != expectedMaximumRank
                || portModel->eigenvalueTolerance
                    != expectedEigenvalueTolerance
                || portModel->eigensolverTolerance
                    != expectedEigensolverTolerance
                || portModel->eigensolverMaximumIterations
                    != expectedEigensolverIterations
                || portModel->relativeDeflationTolerance
                    != options.rankTolerance
                || portModel->innerSolver != expectedInnerSolver
                || portModel->innerSolverTolerance
                    != expectedInnerTolerance
                || portModel->innerSolverMaximumIterations
                    != expectedInnerIterations) {
                throw std::runtime_error(
                    "[Optimal port] Loaded basis build options do not match the requested configuration.");
            }
            const std::vector<PortPatch> expectedPatches =
                buildOptimalPortPatches(
                    mesh, partition,
                    steklovSchur
                        ? 0 : options.optimalPortOversamplingLayers);
            if ((optimalTransfer || randomizedTransfer
                    || hybridRandomized)
                && portModel->formatVersion >= 6) {
                std::uint64_t traceFingerprint =
                    UINT64_C(1469598103934665603);
                for (const PortPatch& patch : expectedPatches) {
                    hashValue(traceFingerprint, patch.interfaceId);
                    hashValue(
                        traceFingerprint, patch.sourceFingerprint);
                }
                if (portModel->traceSourceFingerprint
                    != traceFingerprint) {
                    throw std::runtime_error(
                        "[Optimal port] Cached trace-source fingerprint does not match topology.");
                }
            }
            if (portModel->ports.size() != expectedPatches.size()) {
                throw std::runtime_error(
                    "[Optimal port] Loaded physical-port count does not match the current owner map.");
            }
            for (std::size_t index = 0;
                 index < expectedPatches.size(); ++index) {
                const LocalPortBasis& port = portModel->ports[index];
                const PortPatch& patch = expectedPatches[index];
                if (port.interfaceId != patch.interfaceId
                    || port.leftSubdomain != patch.leftSubdomain
                    || port.rightSubdomain != patch.rightSubdomain
                    || port.interfaceIndices != patch.target
                    || port.sourceIndices != patch.source
                    || port.patchSubdomains != patch.patchSubdomains
                    || port.targetFingerprint != patch.targetFingerprint
                    || port.sourceFingerprint != patch.sourceFingerprint) {
                    throw std::runtime_error(
                        "[Optimal port] Loaded target/source patch mapping does not match the current topology.");
                }
            }
        } else if (steklovSchur) {
            SteklovPortOptions steklovOptions;
            steklovOptions.requestedRank = std::max(1, options.localPortRank);
            steklovOptions.relativeTolerance = options.rankTolerance;
            portModel = std::make_unique<LocalPortModel>(
                buildSteklovSchurPortModel(
                    mesh, partition, *preparedPortDynamicModel,
                    descriptor.interfaceTraceMassDiagonal,
                    descriptor.interfacePenaltyMassDiagonal,
                    steklovOptions));
        } else if (hybridRandomized) {
            if (options.milestone8AdaptiveProduction) {
                const std::vector<PortPatch> adaptivePatches =
                    buildOptimalPortPatches(
                        mesh, partition,
                        options.optimalPortOversamplingLayers);
                if (adaptivePatches.size() != 25) {
                    throw std::runtime_error(
                        "[M8.9 adaptive] The frozen RRAM26 policy "
                        "requires exactly 25 physical interfaces.");
                }
                std::map<std::tuple<int, int, int>, std::vector<int>>
                    rankGroups;
                std::set<int> policyInterfaces;
                for (const PortPatch& patch : adaptivePatches) {
                    if (patch.interfaceId < 0
                        || patch.interfaceId >= 25
                        || !policyInterfaces.insert(
                            patch.interfaceId).second) {
                        throw std::runtime_error(
                            "[M8.9 adaptive] Physical-interface ids "
                            "must be the unique range 0..24.");
                    }
                    const AdaptiveProductionRank rank =
                        milestone8AdaptiveProductionRank(
                            patch.interfaceId);
                    rankGroups[std::make_tuple(
                        rank.historyRank,
                        rank.randomizedRank,
                        rank.residualRank)].push_back(
                            patch.interfaceId);
                }
                if (policyInterfaces.size() != 25
                    || *policyInterfaces.begin() != 0
                    || *policyInterfaces.rbegin() != 24) {
                    throw std::runtime_error(
                        "[M8.9 adaptive] Frozen policy coverage is "
                        "incomplete.");
                }

                ReducedDynamicSchurOperator sharedHybridSchur(
                    *preparedPortDynamicModel, true);
                randomizedTransferBuild = std::make_unique<
                    port::RandomizedTransferBuildResult>();
                residualKrylovBuild = std::make_unique<
                    ResidualKrylovBuildResult>();
                bool aggregateInitialized = false;
                for (const auto& group : rankGroups) {
                    const int historyRank =
                        std::get<0>(group.first);
                    const int randomizedRank =
                        std::get<1>(group.first);
                    const int residualRank =
                        std::get<2>(group.first);
                    const std::vector<int>& interfaceIds =
                        group.second;
                    std::cout
                        << "[M8.9 adaptive] group history/randomized/"
                        << "residual=" << historyRank << '/'
                        << randomizedRank << '/' << residualRank
                        << ", interfaces=" << interfaceIds.size()
                        << '\n' << std::flush;

                    port::RandomizedTransferPortOptions randomizedOptions;
                    randomizedOptions.fluxAware =
                        options.fluxAwarePort;
                    randomizedOptions.fluxType =
                        options.fluxAwareFluxType;
                    randomizedOptions.requestedRank = randomizedRank;
                    randomizedOptions.oversampling =
                        options.randomizedPortOversampling;
                    randomizedOptions.powerIterations =
                        options.randomizedPortPowerIterations;
                    randomizedOptions.seed = options.randomizedPortSeed;
                    randomizedOptions.innerProduct =
                        options.optimalPortInnerProduct;
                    randomizedOptions.sourceMode =
                        effectiveOptimalPortSourceMode;
                    randomizedOptions.oversamplingLayers =
                        options.optimalPortOversamplingLayers;
                    randomizedOptions.relativeDeflationTolerance =
                        options.rankTolerance;
                    randomizedOptions.selectedInterfaceIds =
                        interfaceIds;
                    randomizedOptions.innerSolver.innerSolver =
                        "woodbury-exact";
                    randomizedOptions.innerSolver.relativeTolerance =
                        options.optimalPortInnerTolerance;
                    randomizedOptions.innerSolver.maximumIterations =
                        options.optimalPortInnerMaximumIterations;
                    randomizedOptions.innerSolver
                        .refinementMaximumIterations =
                        options
                            .optimalPortInnerRefinementMaximumIterations;
                    randomizedOptions.innerSolver.refinementTolerance =
                        options.optimalPortInnerRefinementTolerance;
                    port::RandomizedTransferBuildResult randomizedGroup =
                        port::buildRandomizedTransferPortModel(
                            mesh, physics, partition, *preparedPortDynamicModel,
                            descriptor.interfaceTraceMassDiagonal,
                            descriptor.interfacePenaltyMassDiagonal,
                            descriptor.input, descriptor.sourceChannels,
                            descriptor.boundaryRhs,
                            reducedHistory.condensedColumns,
                            reducedHistory.channels,
                            randomizedOptions, &sharedHybridSchur);
                    for (const auto& interfaceResult :
                         randomizedGroup.interfaces) {
                        if (interfaceResult.status.rfind(
                                "success", 0) != 0) {
                            throw std::runtime_error(
                                "[M8.9 adaptive] Randomized interface "
                                + std::to_string(
                                    interfaceResult
                                        .physicalInterfaceId)
                                + " failed: "
                                + interfaceResult.status);
                        }
                    }

                    ResidualKrylovPortOptions residualOptions;
                    residualOptions.basisMethod = "hybrid-randomized";
                    residualOptions.innerProduct =
                        options.optimalPortInnerProduct;
                    residualOptions.oversamplingLayers =
                        options.optimalPortOversamplingLayers;
                    residualOptions.maximumEnrichmentRank =
                        residualRank;
                    residualOptions.maximumSweeps =
                        options.residualKrylovMaximumSweeps;
                    residualOptions.blockSize =
                        options.residualKrylovBlockSize;
                    residualOptions.residualTolerance =
                        options.residualKrylovTolerance;
                    residualOptions.relativeDeflationTolerance =
                        options.rankTolerance;
                    residualOptions.probeMode =
                        options.residualKrylovProbeMode;
                    residualOptions.historyCompressionMethod =
                        options.historyCompressionMethod;
                    residualOptions.historyCompressionRank =
                        historyRank;
                    residualOptions.historyCompressionTolerance =
                        options.historyCompressionTolerance;
                    residualOptions.selectedInterfaceIds =
                        interfaceIds;
                    residualOptions.innerSolver.innerSolver =
                        options.residualKrylovInnerSolver;
                    residualOptions.innerSolver.relativeTolerance =
                        options.optimalPortInnerTolerance;
                    residualOptions.innerSolver.maximumIterations =
                        options.optimalPortInnerMaximumIterations;
                    residualOptions.innerSolver
                        .refinementMaximumIterations =
                        options
                            .optimalPortInnerRefinementMaximumIterations;
                    residualOptions.innerSolver.refinementTolerance =
                        options.optimalPortInnerRefinementTolerance;
                    ResidualKrylovBuildResult residualGroup =
                        buildResidualKrylovPortModel(
                            mesh, partition,
                            *preparedPortDynamicModel,
                            descriptor.interfaceTraceMassDiagonal,
                            descriptor.interfacePenaltyMassDiagonal,
                            descriptor.input,
                            descriptor.sourceChannels,
                            descriptor.boundaryRhs,
                            reducedHistory.condensedColumns,
                            reducedHistory.channels,
                            residualOptions, &sharedHybridSchur,
                            &randomizedGroup.model);
                    for (const auto& interfaceResult :
                         residualGroup.interfaces) {
                        if (interfaceResult.status != "success") {
                            throw std::runtime_error(
                                "[M8.9 adaptive] Residual interface "
                                + std::to_string(
                                    interfaceResult.interfaceId)
                                + " failed: "
                                + interfaceResult.status);
                        }
                    }

                    if (!aggregateInitialized) {
                        *randomizedTransferBuild =
                            std::move(randomizedGroup);
                        *residualKrylovBuild =
                            std::move(residualGroup);
                        aggregateInitialized = true;
                    } else {
                        randomizedTransferBuild->patchSetupSeconds +=
                            randomizedGroup.patchSetupSeconds;
                        randomizedTransferBuild->sourceSetupSeconds +=
                            randomizedGroup.sourceSetupSeconds;
                        randomizedTransferBuild
                            ->probeGenerationSeconds +=
                            randomizedGroup.probeGenerationSeconds;
                        randomizedTransferBuild->transferApplySeconds +=
                            randomizedGroup.transferApplySeconds;
                        randomizedTransferBuild
                            ->transposeApplySeconds +=
                            randomizedGroup.transposeApplySeconds;
                        randomizedTransferBuild->qrSeconds +=
                            randomizedGroup.qrSeconds;
                        randomizedTransferBuild->totalSeconds +=
                            randomizedGroup.totalSeconds;
                        randomizedTransferBuild->model
                            .reducedInterfaceDofs +=
                            randomizedGroup.model
                                .reducedInterfaceDofs;
                        randomizedTransferBuild->model.modelBytes +=
                            randomizedGroup.model.modelBytes;
                        for (LocalPortBasis& port :
                             randomizedGroup.model.ports) {
                            randomizedTransferBuild->model.ports
                                .push_back(std::move(port));
                        }
                        for (auto& diagnostics :
                             randomizedGroup.interfaces) {
                            randomizedTransferBuild->interfaces
                                .push_back(
                                    std::move(diagnostics));
                        }

                        residualKrylovBuild->mandatoryModeSeconds +=
                            residualGroup.mandatoryModeSeconds;
                        residualKrylovBuild->probeSetupSeconds +=
                            residualGroup.probeSetupSeconds;
                        residualKrylovBuild->enrichmentSeconds +=
                            residualGroup.enrichmentSeconds;
                        residualKrylovBuild
                            ->orthogonalizationSeconds +=
                            residualGroup.orthogonalizationSeconds;
                        residualKrylovBuild->totalSeconds +=
                            residualGroup.totalSeconds;
                        residualKrylovBuild->model
                            .reducedInterfaceDofs +=
                            residualGroup.model.reducedInterfaceDofs;
                        residualKrylovBuild->model.modelBytes +=
                            residualGroup.model.modelBytes;
                        for (LocalPortBasis& port :
                             residualGroup.model.ports) {
                            residualKrylovBuild->model.ports
                                .push_back(std::move(port));
                        }
                        for (auto& diagnostics :
                             residualGroup.interfaces) {
                            residualKrylovBuild->interfaces
                                .push_back(
                                    std::move(diagnostics));
                        }
                    }
                }
                if (!aggregateInitialized
                    || randomizedTransferBuild->interfaces.size() != 25
                    || residualKrylovBuild->interfaces.size() != 25
                    || residualKrylovBuild->model.ports.size() != 25) {
                    throw std::runtime_error(
                        "[M8.9 adaptive] Grouped build did not cover "
                        "all 25 interfaces.");
                }
                const auto portOrder = [](const LocalPortBasis& left,
                                          const LocalPortBasis& right) {
                    return left.interfaceId < right.interfaceId;
                };
                std::sort(
                    randomizedTransferBuild->model.ports.begin(),
                    randomizedTransferBuild->model.ports.end(),
                    portOrder);
                std::sort(
                    residualKrylovBuild->model.ports.begin(),
                    residualKrylovBuild->model.ports.end(),
                    portOrder);
                std::sort(
                    randomizedTransferBuild->interfaces.begin(),
                    randomizedTransferBuild->interfaces.end(),
                    [](const port::PortBasisResult& left,
                       const port::PortBasisResult& right) {
                        return left.physicalInterfaceId
                            < right.physicalInterfaceId;
                    });
                std::sort(
                    residualKrylovBuild->interfaces.begin(),
                    residualKrylovBuild->interfaces.end(),
                    [](const ResidualKrylovInterfaceDiagnostics& left,
                       const ResidualKrylovInterfaceDiagnostics& right) {
                        return left.interfaceId < right.interfaceId;
                    });
                residualKrylovBuild->model.historyCompressionRank =
                    options.historyCompressionRank;
                residualKrylovBuild->model.requestedRank =
                    options.randomizedPortRank;
                residualKrylovBuild->model.maximumRank =
                    options.residualKrylovMaximumRank;
                residualKrylovBuild->model.methodDescription =
                    "Local Block Arnoldi with Adaptive Operator-Informed "
                    "Port Space, Randomized Transfer Enrichment, and "
                    "Schur-Residual Enrichment";
                residualKrylovBuild->model.basisSeconds =
                    randomizedTransferBuild->totalSeconds
                    + residualKrylovBuild->totalSeconds;
                portModel = std::make_unique<LocalPortModel>(
                    std::move(residualKrylovBuild->model));
            } else {
            port::RandomizedTransferPortOptions randomizedOptions;
            randomizedOptions.fluxAware = options.fluxAwarePort;
            randomizedOptions.fluxType =
                options.fluxAwareFluxType;
            if (options.fluxAwarePort
                && options.projectionDiagnosis) {
                randomizedOptions.selectedInterfaceIds =
                    options.projectionInterfaceIds;
            }
            randomizedOptions.requestedRank =
                options.randomizedPortRank;
            randomizedOptions.oversampling =
                options.randomizedPortOversampling;
            randomizedOptions.powerIterations =
                options.randomizedPortPowerIterations;
            randomizedOptions.seed = options.randomizedPortSeed;
            randomizedOptions.innerProduct =
                options.optimalPortInnerProduct;
            randomizedOptions.sourceMode =
                effectiveOptimalPortSourceMode;
            randomizedOptions.oversamplingLayers =
                options.optimalPortOversamplingLayers;
            randomizedOptions.relativeDeflationTolerance =
                options.rankTolerance;
            randomizedOptions.innerSolver.innerSolver =
                "woodbury-exact";
            randomizedOptions.innerSolver.relativeTolerance =
                options.optimalPortInnerTolerance;
            randomizedOptions.innerSolver.maximumIterations =
                options.optimalPortInnerMaximumIterations;
            randomizedOptions.innerSolver
                .refinementMaximumIterations =
                options.optimalPortInnerRefinementMaximumIterations;
            randomizedOptions.innerSolver.refinementTolerance =
                options.optimalPortInnerRefinementTolerance;
            ReducedDynamicSchurOperator sharedHybridSchur(
                *preparedPortDynamicModel, true);
            randomizedTransferBuild =
                std::make_unique<
                    port::RandomizedTransferBuildResult>(
                    port::buildRandomizedTransferPortModel(
                        mesh, physics, partition, *preparedPortDynamicModel,
                        descriptor.interfaceTraceMassDiagonal,
                        descriptor.interfacePenaltyMassDiagonal,
                        descriptor.input, descriptor.sourceChannels,
                        descriptor.boundaryRhs,
                        reducedHistory.condensedColumns,
                        reducedHistory.channels,
                        randomizedOptions, &sharedHybridSchur));
            for (const port::PortBasisResult& interfaceResult :
                 randomizedTransferBuild->interfaces) {
                if (interfaceResult.status.rfind("success", 0)
                    != 0) {
                    throw std::runtime_error(
                        "[Hybrid port] Randomized interface "
                        + std::to_string(
                            interfaceResult.physicalInterfaceId)
                        + " failed: " + interfaceResult.status);
                }
            }
            std::unique_ptr<LocalPortModel> fluxAugmentedInitial;
            if (options.fluxAwarePort && !loadPath.empty()) {
                LocalPortModel frozen = loadLocalPortModel(loadPath);
                fluxAugmentedInitial =
                    std::make_unique<LocalPortModel>(frozen);
                fluxAugmentedInitial->ports.clear();
                fluxAugmentedInitial->reducedInterfaceDofs = 0;
                for (const LocalPortBasis& enrichment :
                     randomizedTransferBuild->model.ports) {
                    const auto found = std::find_if(
                        frozen.ports.begin(), frozen.ports.end(),
                        [&](const LocalPortBasis& candidate) {
                            return candidate.interfaceId
                                == enrichment.interfaceId;
                        });
                    if (found == frozen.ports.end()
                        || found->rows != enrichment.rows
                        || found->interfaceIndices
                            != enrichment.interfaceIndices) {
                        throw std::runtime_error(
                            "[Flux-aware port] Frozen M8.9 basis "
                            "does not match the selected interface.");
                    }
                    LocalPortBasis combined = *found;
                    combined.basis.insert(
                        combined.basis.end(),
                        enrichment.basis.begin(),
                        enrichment.basis.end());
                    combined.rank += enrichment.rank;
                    combined.acceptedColumns = combined.rank;
                    combined.requestedTransferRank +=
                        enrichment.requestedTransferRank;
                    fluxAugmentedInitial->reducedInterfaceDofs +=
                        combined.rank;
                    fluxAugmentedInitial->ports.push_back(
                        std::move(combined));
                }
                // The residual builder accepts an operator-side initial
                // transfer span under this compatibility name.  The span
                // itself is frozen-M8.9 plus flux-aware randomized columns.
                fluxAugmentedInitial->basisMethod =
                    "randomized-transfer";
            }
            ResidualKrylovPortOptions residualOptions;
            residualOptions.basisMethod = "hybrid-randomized";
            residualOptions.innerProduct =
                options.optimalPortInnerProduct;
            residualOptions.oversamplingLayers =
                options.optimalPortOversamplingLayers;
            residualOptions.maximumEnrichmentRank =
                options.residualKrylovMaximumRank;
            residualOptions.maximumSweeps =
                options.residualKrylovMaximumSweeps;
            residualOptions.blockSize =
                options.residualKrylovBlockSize;
            residualOptions.residualTolerance =
                options.residualKrylovTolerance;
            residualOptions.relativeDeflationTolerance =
                options.rankTolerance;
            residualOptions.probeMode =
                options.residualKrylovProbeMode;
            residualOptions.historyCompressionMethod =
                options.historyCompressionMethod;
            residualOptions.historyCompressionRank =
                options.historyCompressionRank;
            residualOptions.historyCompressionTolerance =
                options.historyCompressionTolerance;
            residualOptions.innerSolver.innerSolver =
                options.residualKrylovInnerSolver;
            residualOptions.innerSolver.relativeTolerance =
                options.optimalPortInnerTolerance;
            residualOptions.innerSolver.maximumIterations =
                options.optimalPortInnerMaximumIterations;
            residualOptions.innerSolver.refinementMaximumIterations =
                options.optimalPortInnerRefinementMaximumIterations;
            residualOptions.innerSolver.refinementTolerance =
                options.optimalPortInnerRefinementTolerance;
            if (options.fluxAwarePort
                && options.projectionDiagnosis) {
                residualOptions.selectedInterfaceIds =
                    options.projectionInterfaceIds;
            }
            residualKrylovBuild =
                std::make_unique<ResidualKrylovBuildResult>(
                    buildResidualKrylovPortModel(
                        mesh, partition, *preparedPortDynamicModel,
                        descriptor.interfaceTraceMassDiagonal,
                        descriptor.interfacePenaltyMassDiagonal,
                        descriptor.input, descriptor.sourceChannels,
                        descriptor.boundaryRhs,
                        reducedHistory.condensedColumns,
                        reducedHistory.channels,
                        residualOptions, &sharedHybridSchur,
                        fluxAugmentedInitial
                            ? fluxAugmentedInitial.get()
                            : &randomizedTransferBuild->model));
            residualKrylovBuild->model.basisSeconds +=
                randomizedTransferBuild->totalSeconds;
            portModel = std::make_unique<LocalPortModel>(
                std::move(residualKrylovBuild->model));
            }
        } else if (operatorInformedPort) {
            ResidualKrylovPortOptions residualOptions;
            residualOptions.basisMethod = options.portBasisMethod;
            residualOptions.innerProduct =
                options.optimalPortInnerProduct;
            residualOptions.oversamplingLayers =
                options.optimalPortOversamplingLayers;
            residualOptions.maximumEnrichmentRank =
                mandatoryOnly ? 0
                    : options.residualKrylovMaximumRank;
            residualOptions.maximumSweeps =
                options.residualKrylovMaximumSweeps;
            residualOptions.blockSize =
                options.residualKrylovBlockSize;
            residualOptions.residualTolerance =
                options.residualKrylovTolerance;
            residualOptions.relativeDeflationTolerance =
                options.rankTolerance;
            residualOptions.probeMode =
                options.residualKrylovProbeMode;
            residualOptions.historyCompressionMethod =
                options.historyCompressionMethod;
            residualOptions.historyCompressionRank =
                options.historyCompressionRank;
            residualOptions.historyCompressionTolerance =
                options.historyCompressionTolerance;
            residualOptions.innerSolver.innerSolver =
                options.residualKrylovInnerSolver;
            residualOptions.innerSolver.relativeTolerance =
                options.optimalPortInnerTolerance;
            residualOptions.innerSolver.maximumIterations =
                options.optimalPortInnerMaximumIterations;
            residualOptions.innerSolver.refinementMaximumIterations =
                options.optimalPortInnerRefinementMaximumIterations;
            residualOptions.innerSolver.refinementTolerance =
                options.optimalPortInnerRefinementTolerance;
            residualKrylovBuild =
                std::make_unique<ResidualKrylovBuildResult>(
                    buildResidualKrylovPortModel(
                        mesh, partition, *preparedPortDynamicModel,
                        descriptor.interfaceTraceMassDiagonal,
                        descriptor.interfacePenaltyMassDiagonal,
                        descriptor.input, descriptor.sourceChannels,
                        descriptor.boundaryRhs,
                        reducedHistory.condensedColumns,
                        reducedHistory.channels,
                        residualOptions));
            portModel = std::make_unique<LocalPortModel>(
                std::move(residualKrylovBuild->model));
        } else if (randomizedTransfer) {
            port::RandomizedTransferPortOptions randomizedOptions;
            randomizedOptions.fluxAware = options.fluxAwarePort;
            randomizedOptions.fluxType =
                options.fluxAwareFluxType;
            randomizedOptions.requestedRank =
                options.randomizedPortRank;
            randomizedOptions.oversampling =
                options.randomizedPortOversampling;
            randomizedOptions.powerIterations =
                options.randomizedPortPowerIterations;
            randomizedOptions.seed = options.randomizedPortSeed;
            randomizedOptions.innerProduct =
                options.optimalPortInnerProduct;
            randomizedOptions.sourceMode =
                effectiveOptimalPortSourceMode;
            randomizedOptions.oversamplingLayers =
                options.optimalPortOversamplingLayers;
            randomizedOptions.relativeDeflationTolerance =
                options.rankTolerance;
            randomizedOptions.innerSolver.innerSolver =
                "woodbury-exact";
            randomizedOptions.innerSolver.relativeTolerance =
                options.optimalPortInnerTolerance;
            randomizedOptions.innerSolver.maximumIterations =
                options.optimalPortInnerMaximumIterations;
            randomizedOptions.innerSolver
                .refinementMaximumIterations =
                options.optimalPortInnerRefinementMaximumIterations;
            randomizedOptions.innerSolver.refinementTolerance =
                options.optimalPortInnerRefinementTolerance;
            ReducedDynamicSchurOperator sharedRandomizedSchur(
                *preparedPortDynamicModel, true);
            randomizedTransferBuild =
                std::make_unique<
                    port::RandomizedTransferBuildResult>(
                    port::buildRandomizedTransferPortModel(
                        mesh, physics, partition, *preparedPortDynamicModel,
                        descriptor.interfaceTraceMassDiagonal,
                        descriptor.interfacePenaltyMassDiagonal,
                        descriptor.input, descriptor.sourceChannels,
                        descriptor.boundaryRhs,
                        reducedHistory.condensedColumns,
                        reducedHistory.channels,
                        randomizedOptions,
                        &sharedRandomizedSchur));
            for (const port::PortBasisResult& interfaceResult :
                 randomizedTransferBuild->interfaces) {
                if (interfaceResult.status.rfind("success", 0)
                    != 0) {
                    throw std::runtime_error(
                        "[Randomized port] Interface "
                        + std::to_string(
                            interfaceResult.physicalInterfaceId)
                        + " failed: " + interfaceResult.status);
                }
            }
            if (options.randomizedPortCompareOptimal) {
                OptimalTransferPortOptions referenceOptions;
                referenceOptions.requestedRank =
                    options.randomizedPortRank;
                referenceOptions.rankMode = "fixed";
                referenceOptions.minimumRank =
                    options.randomizedPortRank;
                referenceOptions.maximumRank =
                    options.randomizedPortRank;
                referenceOptions.innerProduct =
                    options.optimalPortInnerProduct;
                referenceOptions.oversamplingLayers =
                    options.optimalPortOversamplingLayers;
                referenceOptions.eigensolverMaximumIterations =
                    std::max(
                        200,
                        options.optimalPortEigenMaximumIterations);
                referenceOptions.eigensolverTolerance =
                    options.optimalPortEigenTolerance;
                referenceOptions.relativeDeflationTolerance =
                    options.rankTolerance;
                referenceOptions.sourceMode =
                    effectiveOptimalPortSourceMode;
                referenceOptions.ablationMode =
                    effectiveOptimalPortSourceMode
                            == "trace-only"
                    ? "trace-transfer-only"
                    : (effectiveOptimalPortSourceMode
                                == "generalized-dynamic"
                        ? "generalized-transfer-only"
                        : "transfer-only");
                referenceOptions.innerSolver.innerSolver =
                    "woodbury-exact";
                referenceOptions.innerSolver.relativeTolerance =
                    options.optimalPortInnerTolerance;
                referenceOptions.innerSolver.maximumIterations =
                    options.optimalPortInnerMaximumIterations;
                referenceOptions.innerSolver
                    .refinementMaximumIterations =
                    options
                        .optimalPortInnerRefinementMaximumIterations;
                referenceOptions.innerSolver.refinementTolerance =
                    options.optimalPortInnerRefinementTolerance;
                const OptimalPortBuildResult referenceBuild =
                    buildOptimalTransferPortModel(
                        mesh, partition, *preparedPortDynamicModel,
                        descriptor.interfaceTraceMassDiagonal,
                        descriptor.interfacePenaltyMassDiagonal,
                        descriptor.input, descriptor.sourceChannels,
                        descriptor.boundaryRhs,
                        reducedHistory.condensedColumns,
                        reducedHistory.channels,
                        referenceOptions,
                        &sharedRandomizedSchur);
                std::ofstream comparison(
                    outputDirectory
                    / "randomized_port_subspace_comparison.csv");
                comparison
                    << "interface_id,randomized_rank,optimal_rank,"
                    "maximum_principal_angle_radians,"
                    "projector_difference,optimal_target_solve_count,"
                    "optimal_basis_build_time_s,source_mode,status\n"
                    << std::setprecision(17);
                for (const LocalPortBasis& randomizedPort :
                     randomizedTransferBuild->model.ports) {
                    const auto referencePort = std::find_if(
                        referenceBuild.model.ports.begin(),
                        referenceBuild.model.ports.end(),
                        [&](const LocalPortBasis& value) {
                            return value.interfaceId
                                == randomizedPort.interfaceId;
                        });
                    const auto referenceDiagnostic = std::find_if(
                        referenceBuild.interfaces.begin(),
                        referenceBuild.interfaces.end(),
                        [&](const OptimalPortInterfaceDiagnostics& value) {
                            return value.interfaceId
                                == randomizedPort.interfaceId;
                        });
                    if (referencePort
                            == referenceBuild.model.ports.end()
                        || referenceDiagnostic
                            == referenceBuild.interfaces.end()) {
                        throw std::runtime_error(
                            "[Randomized port] Missing same-source "
                            "Optimal reference interface.");
                    }
                    const port::WeightedPortSubspaceComparison metrics =
                        port::compareWeightedPortSubspaces(
                            randomizedPort, *referencePort,
                            *preparedPortDynamicModel,
                            descriptor.interfaceTraceMassDiagonal,
                            descriptor.interfacePenaltyMassDiagonal,
                            options.optimalPortInnerProduct);
                    comparison << randomizedPort.interfaceId << ','
                        << randomizedPort.rank << ','
                        << referencePort->rank << ','
                        << metrics.maximumPrincipalAngleRadians << ','
                        << metrics.projectorDifference << ','
                        << referenceDiagnostic->innerSolver
                            .solveRightHandSides
                        << ',' << referenceBuild.totalSeconds << ','
                        << effectiveOptimalPortSourceMode << ','
                        << (referenceDiagnostic->eigenConverged
                            ? "success"
                            : referenceDiagnostic->eigenStatus)
                        << '\n';
                }
            }
            portModel = std::make_unique<LocalPortModel>(
                std::move(randomizedTransferBuild->model));
        } else {
            OptimalTransferPortOptions optimalOptions;
            optimalOptions.requestedRank = options.optimalPortRank;
            optimalOptions.rankMode = options.optimalPortRankMode;
            optimalOptions.rankFile = options.optimalPortRankFile;
            optimalOptions.eigenvalueTolerance =
                options.optimalPortEigenvalueTolerance;
            optimalOptions.minimumRank = options.optimalPortMinimumRank;
            optimalOptions.maximumRank = options.optimalPortMaximumRank;
            optimalOptions.innerProduct =
                options.optimalPortInnerProduct;
            optimalOptions.oversamplingLayers =
                options.optimalPortOversamplingLayers;
            optimalOptions.eigensolverMaximumIterations =
                options.optimalPortEigenMaximumIterations;
            optimalOptions.eigensolverTolerance =
                options.optimalPortEigenTolerance;
            optimalOptions.relativeDeflationTolerance =
                options.rankTolerance;
            optimalOptions.ablationMode = options.optimalPortAblation;
            optimalOptions.sourceMode =
                effectiveOptimalPortSourceMode;
            optimalOptions.innerSolver.innerSolver =
                options.optimalPortInnerSolver;
            optimalOptions.innerSolver.relativeTolerance =
                options.optimalPortInnerTolerance;
            optimalOptions.innerSolver.maximumIterations =
                options.optimalPortInnerMaximumIterations;
            optimalOptions.innerSolver.refinementMaximumIterations =
                options.optimalPortInnerRefinementMaximumIterations;
            optimalOptions.innerSolver.refinementTolerance =
                options.optimalPortInnerRefinementTolerance;
            optimalPortBuild =
                std::make_unique<OptimalPortBuildResult>(
                    buildOptimalTransferPortModel(
                        mesh, partition, *preparedPortDynamicModel,
                        descriptor.interfaceTraceMassDiagonal,
                        descriptor.interfacePenaltyMassDiagonal,
                        descriptor.input, descriptor.sourceChannels,
                        descriptor.boundaryRhs,
                        reducedHistory.condensedColumns,
                        reducedHistory.channels,
                        optimalOptions));
            portModel = std::make_unique<LocalPortModel>(
                std::move(optimalPortBuild->model));
        }
        if (loadPath.empty()) {
            portModel->timeStep = options.timeStep;
            portModel->meshFingerprint = descriptor.fingerprints.mesh;
            portModel->materialFingerprint =
                descriptor.fingerprints.conductivity;
            portModel->operatorFingerprint = dynamicOperatorFingerprint;
            portModel->inputFingerprint = descriptor.fingerprints.input;
            portModel->penaltyFingerprint = penaltyFingerprint;
            portModel->boundaryFingerprint =
                descriptor.fingerprints.boundary;
            portModel->sourceMetadataFingerprint =
                descriptor.fingerprints.sources;
            portModel->rankFileFingerprint = rankFileFingerprint;
            portModel->buildTimestamp = std::to_string(
                std::chrono::duration_cast<std::chrono::seconds>(
                    std::chrono::system_clock::now()
                        .time_since_epoch()).count());
            portModel->sourceCommit = DDM_SCHUR_SOURCE_COMMIT;
            std::filesystem::path savePath = options.savePath;
            if (savePath.empty()) {
                savePath = outputDirectory / "local_port_basis.bin";
            } else if (savePath.extension().empty()) {
                savePath /= "local_port_basis.bin";
            }
            const auto serializationStart = Clock::now();
            saveLocalPortModel(*portModel, savePath);
            if (randomizedTransferBuild) {
                randomizedTransferBuild->serializationSeconds =
                    elapsed(serializationStart);
                for (port::PortBasisResult& interfaceResult :
                     randomizedTransferBuild->interfaces) {
                    interfaceResult.serializationSeconds =
                        randomizedTransferBuild
                            ->serializationSeconds;
                }
            }
        }
        portSnapshotSeconds = 0.0;
        portBasisSeconds = portModel->basisSeconds;
        writePortRankDiagnostics(*portModel, outputDirectory);
        if (optimalPortBuild) {
            writeOptimalPortDiagnostics(
                *portModel, *optimalPortBuild, outputDirectory);
        }
        if (residualKrylovBuild) {
            writeResidualKrylovDiagnostics(
                *residualKrylovBuild, outputDirectory);
        }
        if (randomizedTransferBuild) {
            writeRandomizedTransferDiagnostics(
                *randomizedTransferBuild, outputDirectory);
        }
        std::cout << "[Optimal port] method=" << portModel->basisMethod
                  << ", physical interfaces=" << portModel->ports.size()
                  << ", full/reduced interface="
                  << portModel->fullInterfaceDofs << '/'
                  << portModel->reducedInterfaceDofs
                  << ", setup=" << elapsed(portStart) << " s\n";
        if (options.projectionDiagnosis) {
            if (!preparedPortDynamicModel || !portModel) {
                throw std::runtime_error(
                    "[Projection diagnosis] Dynamic and port models "
                    "are required.");
            }
            ProjectionDiagnosisOptions projectionOptions;
            projectionOptions.interfaceIds =
                options.projectionInterfaceIds;
            projectionOptions.timeStep = options.timeStep;
            projectionOptions.time = options.timeStep;
            projectionOptions.powers =
                waveform.sample(options.timeStep);
            projectionOptions.outputDirectory = outputDirectory;
            projectionOptions.fluxOperatorAudit =
                options.fluxOperatorAudit;
            projectionOptions.schurOptions.maxIterations =
                options.interfaceMaxIterations;
            projectionOptions.schurOptions.restart =
                options.interfaceRestart;
            projectionOptions.schurOptions.relativeTolerance =
                options.interfaceTolerance;
            projectionOptions.schurOptions.coarseLinearXY =
                options.coarseLinearXY;
            projectionOptions.schurOptions.coarseLinearZ =
                options.coarseLinearZ;
            projectionOptions.schurOptions.proxyEnabled =
                options.proxyEnabled;
            projectionOptions.schurOptions.proxyDisableCoarse =
                options.proxyDisableCoarse;
            projectionOptions.schurOptions
                .proxyHighConductivityThreshold =
                options.proxyHighConductivityThreshold;
            projectionOptions.schurOptions
                .proxyUseMaterialConnectivity =
                options.proxyUseMaterialConnectivity;
            projectionOptions.schurOptions.proxyRing =
                options.proxyRing;
            projectionOptions.schurOptions.proxyProbeColumns =
                options.proxyProbeColumns;
            projectionOptions.schurOptions.proxyBlockSize =
                options.proxyBlockSize;
            projectionOptions.schurOptions
                .proxyValidateBlockEquivalence =
                options.proxyValidateBlockEquivalence;
            const ProjectionDiagnosisSummary diagnosis =
                runFullInterfaceProjectionDiagnosis(
                    mesh, physics, partition, descriptor,
                    *preparedPortDynamicModel, *portModel,
                    thetaInitial, reference, boundaryOffset,
                    projectionOptions);
            std::cout
                << "[Projection diagnosis] completed="
                << diagnosis.completedInterfaces << '/'
                << diagnosis.requestedInterfaces
                << ", full-interface setup/solve="
                << diagnosis.fullInterfaceSetupSeconds << '/'
                << diagnosis.fullInterfaceSolveSeconds
                << " s, projection="
                << diagnosis.projectionSeconds
                << " s, status=" << diagnosis.status << '\n';
            return;
        }
        if (options.milestone8ProductionBasisOnly) {
            if (!hybridRandomized || !randomizedTransferBuild
                || !residualKrylovBuild
                || randomizedTransferBuild->interfaces.size()
                    != portModel->ports.size()
                || residualKrylovBuild->interfaces.size()
                    != portModel->ports.size()) {
                throw std::runtime_error(
                    "[M8 production] Complete hybrid interface "
                    "diagnostics are required for the basis-only stop.");
            }
            std::map<int, const port::PortBasisResult*> randomizedById;
            for (const auto& row : randomizedTransferBuild->interfaces) {
                randomizedById.emplace(
                    row.physicalInterfaceId, &row);
            }
            std::map<int, const ResidualKrylovInterfaceDiagnostics*>
                residualById;
            for (const auto& row : residualKrylovBuild->interfaces) {
                residualById.emplace(row.interfaceId, &row);
            }
            int totalPortRank = 0;
            int targetSolveCount = 0;
            int targetSolvePhase33Calls = 0;
            double maximumInterfaceSeconds = 0.0;
            std::size_t peakIncrementalBytes = 0;
            double maximumTargetResidual = 0.0;
            double maximumWeightedAdjointError = 0.0;
            int adaptiveInterfaceCount = 0;
            int defaultInterfaceCount = 0;
            bool passed = true;
            std::ofstream adaptiveDiagnostics;
            if (options.milestone8AdaptiveProduction) {
                adaptiveDiagnostics.open(
                    outputDirectory
                    / "milestone8_adaptive_port_interface_diagnostics.csv");
                adaptiveDiagnostics
                    << "interface_id,history_rank,"
                    "randomized_transfer_rank,residual_rank,"
                    "mandatory_rank,total_port_rank,basis_build_time_s,"
                    "target_solve_count,target_rhs_count,"
                    "peak_incremental_memory_bytes,target_residual,"
                    "weighted_adjoint_error,probe_residual_before,"
                    "probe_residual_after,adaptive,status\n"
                    << std::setprecision(17);
            }
            for (const LocalPortBasis& port : portModel->ports) {
                const auto randomized =
                    randomizedById.find(port.interfaceId);
                const auto residual =
                    residualById.find(port.interfaceId);
                if (randomized == randomizedById.end()
                    || residual == residualById.end()) {
                    throw std::runtime_error(
                        "[M8 production] Interface diagnostics do not "
                        "cover the serialized physical-port partition.");
                }
                totalPortRank += port.rank;
                targetSolveCount +=
                    randomized->second->targetSolveCount
                    + residual->second->targetSolveCount;
                targetSolvePhase33Calls +=
                    randomized->second->targetSolvePhase33Calls
                    + residual->second->innerSolver.solveCalls;
                maximumInterfaceSeconds = std::max(
                    maximumInterfaceSeconds,
                    randomized->second->basisBuildTime
                        + residual->second->totalSeconds);
                peakIncrementalBytes = std::max({
                    peakIncrementalBytes,
                    randomized->second->memoryPeak,
                    residual->second->peakIncrementalMemoryBytes});
                const double interfaceTargetResidual = std::max(
                    randomized->second->residual,
                    residual->second->innerSolver
                        .maximumRelativeResidual);
                const double interfaceAdjointError = std::max(
                    randomized->second->weightedAdjointError,
                    residual->second->weightedAdjointError);
                maximumTargetResidual = std::max(
                    maximumTargetResidual,
                    interfaceTargetResidual);
                maximumWeightedAdjointError = std::max(
                    maximumWeightedAdjointError,
                    interfaceAdjointError);
                const double interfaceSeconds =
                    randomized->second->basisBuildTime
                    + residual->second->totalSeconds;
                const std::size_t interfaceMemory = std::max(
                    randomized->second->memoryPeak,
                    residual->second->peakIncrementalMemoryBytes);
                passed = passed
                    && randomized->second->status.rfind("success", 0) == 0
                    && residual->second->status == "success";
                if (options.milestone8AdaptiveProduction) {
                    const AdaptiveProductionRank policy =
                        milestone8AdaptiveProductionRank(
                            port.interfaceId);
                    if (policy.adaptive) {
                        ++adaptiveInterfaceCount;
                    } else {
                        ++defaultInterfaceCount;
                    }
                    const bool interfacePassed =
                        randomized->second->status.rfind(
                            "success", 0) == 0
                        && residual->second->status == "success"
                        && interfaceSeconds < 180.0
                        && interfaceMemory
                            < UINT64_C(1024) * 1024 * 1024
                        && interfaceTargetResidual < 1.0e-9
                        && interfaceAdjointError < 1.0e-8;
                    passed = passed && interfacePassed;
                    adaptiveDiagnostics
                        << port.interfaceId << ','
                        << policy.historyRank << ','
                        << policy.randomizedRank << ','
                        << policy.residualRank << ','
                        << residual->second->mandatoryRankTotal << ','
                        << port.rank << ','
                        << interfaceSeconds << ','
                        << (randomized->second
                                ->targetSolvePhase33Calls
                            + residual->second
                                ->innerSolver.solveCalls)
                        << ','
                        << (randomized->second->targetSolveCount
                            + residual->second->targetSolveCount)
                        << ',' << interfaceMemory << ','
                        << interfaceTargetResidual << ','
                        << interfaceAdjointError << ','
                        << residual->second
                            ->initialMaximumProbeResidual
                        << ','
                        << residual->second
                            ->finalMaximumProbeResidual
                        << ',' << (policy.adaptive ? 1 : 0)
                        << ','
                        << (interfacePassed
                            ? "passed" : "gate_failed")
                        << '\n';
                }
            }
            std::filesystem::path serializedPath = options.savePath;
            if (serializedPath.empty()) {
                serializedPath =
                    outputDirectory / "local_port_basis.bin";
            } else if (serializedPath.extension().empty()) {
                serializedPath /= "local_port_basis.bin";
            }
            std::error_code fileError;
            const std::uintmax_t serializedBytes =
                std::filesystem::file_size(
                    serializedPath, fileError);
            if (fileError) {
                throw std::runtime_error(
                    "[M8 production] Serialized basis is missing after "
                    "the basis-only build.");
            }
            std::ofstream summary(
                outputDirectory
                / "milestone8_production_basis_summary.csv");
            summary
                << "case,physical_interfaces,history_rank,"
                "randomized_rank,residual_rank,total_port_rank,"
                "total_basis_build_time_s,max_interface_time_s,"
                "peak_incremental_memory_bytes,"
                "process_peak_memory_bytes,target_solve_count,"
                "target_solve_phase33_calls,serialized_model_bytes,"
                "snapshot_used,fom_used_for_basis,pod_used,svd_used,"
                "full_field_read,transient_advanced,status\n"
                << std::setprecision(17)
                << physics.name << ',' << portModel->ports.size()
                << ',' << options.historyCompressionRank
                << ',' << options.randomizedPortRank
                << ',' << options.residualKrylovMaximumRank
                << ',' << totalPortRank
                << ',' << portModel->basisSeconds
                << ',' << maximumInterfaceSeconds
                << ',' << peakIncrementalBytes
                << ',' << peakWorkingSetBytes()
                << ',' << targetSolveCount
                << ',' << targetSolvePhase33Calls
                << ',' << serializedBytes
                << ",0,0,0,0,0,0,"
                << (passed ? "success" : "basis_failed")
                << '\n';
            if (options.milestone8AdaptiveProduction) {
                std::ofstream adaptiveSummary(
                    outputDirectory
                    / "milestone8_adaptive_port_basis_summary.csv");
                adaptiveSummary
                    << "case,physical_interfaces,adaptive_interfaces,"
                    "default_interfaces,total_port_rank,"
                    "total_basis_build_time_s,max_interface_time_s,"
                    "peak_incremental_memory_bytes,"
                    "process_peak_memory_bytes,max_target_residual,"
                    "max_weighted_adjoint_error,target_solve_count,"
                    "target_rhs_count,serialized_model_bytes,"
                    "snapshot_used,fom_used_for_basis,pod_used,"
                    "svd_used,full_field_read,transient_advanced,"
                    "corrected,status\n"
                    << std::setprecision(17)
                    << physics.name << ',' << portModel->ports.size()
                    << ',' << adaptiveInterfaceCount
                    << ',' << defaultInterfaceCount
                    << ',' << totalPortRank
                    << ',' << portModel->basisSeconds
                    << ',' << maximumInterfaceSeconds
                    << ',' << peakIncrementalBytes
                    << ',' << peakWorkingSetBytes()
                    << ',' << maximumTargetResidual
                    << ',' << maximumWeightedAdjointError
                    << ',' << targetSolvePhase33Calls
                    << ',' << targetSolveCount
                    << ',' << serializedBytes
                    << ",0,0,0,0,0,0,0,"
                    << (passed ? "passed" : "gate_failed")
                    << '\n';
            }
            if (!passed) {
                throw std::runtime_error(
                    "[M8 production] At least one physical interface "
                    "failed its hybrid basis diagnostics.");
            }
            std::cout
                << "[M8 production] serialized all-interface basis; "
                << "Schur assembly and transient skipped.\n";
            return;
        }
        if (options.globalInterfaceCoarsePrototype) {
            if (!preparedPortDynamicModel || !portModel) {
                throw std::runtime_error(
                    "[Global coarse] Dynamic and local-port models "
                    "are required.");
            }
            portSolver =
                std::make_unique<LocalPortReducedSchurSolver>(
                    *preparedPortDynamicModel, *portModel);
            ReducedDynamicSchurOperator exactGlobalSchur(
                *preparedPortDynamicModel,
                options.globalInterfaceCoarseInverseMode
                    == "exact-pcg");
            const std::size_t baselineWorkingSet =
                globalCoarseCurrentWorkingSetBytes();
            GeneralizedTransferSourceBlocks sources =
                buildGeneralizedTransferSourceBlocks(
                    *preparedPortDynamicModel,
                    descriptor.input,
                    descriptor.sourceChannels,
                    descriptor.boundaryRhs,
                    {}, 0);
            // The RRAM26 history operator is multi-gigabyte. Transfer its
            // ownership into the prototype source block after all frozen
            // M8.9 provenance checks instead of duplicating it.
            sources.historyChannels =
                reducedHistory.channels;
            sources.history =
                std::move(reducedHistory.condensedColumns);
            GlobalInterfaceCoarseOptions coarseOptions;
            coarseOptions.inverseMode =
                options.globalInterfaceCoarseInverseMode;
            coarseOptions.explicitReference =
                options.globalInterfaceCoarseExplicitReference;
            coarseOptions.requestedRank =
                options.globalInterfaceCoarseRank;
            coarseOptions.candidateDimension =
                options.globalInterfaceCoarseCandidateDimension;
            coarseOptions.maximumIterations =
                options.globalInterfaceCoarseMaximumIterations;
            coarseOptions.innerMaximumIterations =
                options
                    .globalInterfaceCoarseInnerMaximumIterations;
            coarseOptions.historyRank = 64;
            coarseOptions.krylovSweeps =
                options.globalInterfaceCoarseKrylovSweeps;
            coarseOptions.ritzTolerance =
                options.globalInterfaceCoarseTolerance;
            coarseOptions.innerTolerance =
                options.globalInterfaceCoarseInnerTolerance;
            coarseOptions.orthogonalityTolerance = 1.0e-10;
            coarseOptions.deflationTolerance =
                options.rankTolerance;
            coarseOptions.selectedInterfaceIds =
                options.globalInterfaceCoarseInterfaceIds;
            std::vector<double> coarseInterfaceMetric(
                partition.interfaceGlobalDofs.size(), 0.0);
            for (std::size_t interfaceRow = 0;
                 interfaceRow < partition.interfaceGlobalDofs.size();
                 ++interfaceRow) {
                coarseInterfaceMetric[interfaceRow] =
                    descriptor.interfacePenaltyMassDiagonal[
                        static_cast<std::size_t>(
                            partition.interfaceGlobalDofs[
                                interfaceRow])];
            }
            globalCoarseModel =
                std::make_unique<GlobalInterfaceCoarseModel>(
                    buildGlobalInterfaceCoarsePrototype(
                        *portModel, *portSolver,
                        exactGlobalSchur,
                        buildOptimalPortPatches(
                            mesh, partition, 0),
                        coarseInterfaceMetric,
                        sources, coarseOptions));
            auto& coarseDiagnostics =
                globalCoarseModel->diagnostics;
            coarseDiagnostics.baselineWorkingSetBytes =
                baselineWorkingSet;
            coarseDiagnostics.peakWorkingSetBytes = std::max(
                coarseDiagnostics.peakWorkingSetBytes,
                globalCoarseCurrentWorkingSetBytes());
            coarseDiagnostics.peakIncrementalMemoryBytes =
                coarseDiagnostics.peakWorkingSetBytes
                    > baselineWorkingSet
                ? coarseDiagnostics.peakWorkingSetBytes
                    - baselineWorkingSet
                : 0;
            const bool finalGate =
                globalCoarseModel->rank
                    == options.globalInterfaceCoarseRank
                && coarseDiagnostics.innerSolveConverged
                && coarseDiagnostics.maximumInnerSolveResidual
                    <= options.globalInterfaceCoarseInnerTolerance
                && coarseDiagnostics.symmetryError < 1.0e-12
                && coarseDiagnostics.maximumRitzResidual
                    < options.globalInterfaceCoarseTolerance
                && coarseDiagnostics
                    .maximumLocalCoarseOrthogonality < 1.0e-10
                && coarseDiagnostics.coarseSchurGramError
                    < 1.0e-10
                && coarseDiagnostics.basisSeconds < 300.0
                && coarseDiagnostics
                    .peakIncrementalMemoryBytes
                    < UINT64_C(1024) * 1024 * 1024;
            if (finalGate) {
                coarseDiagnostics.status = "passed";
            } else if (coarseDiagnostics.status == "passed") {
                coarseDiagnostics.status =
                    "global_coarse_final_gate_failed";
            }
            writeGlobalInterfaceCoarseDiagnostics(
                *globalCoarseModel, physics.name,
                outputDirectory);
            if (!finalGate) {
                std::cout
                    << "[Global coarse] Rank-4 algebra/resource gate "
                    << "failed with status="
                    << coarseDiagnostics.status
                    << "; transient validation skipped.\n";
                return;
            }
            portSolver->attachGlobalCoarse(
                *globalCoarseModel);
            std::cout
                << "[Global coarse] local/coarse/augmented="
                << portModel->reducedInterfaceDofs << '/'
                << globalCoarseModel->rank << '/'
                << portSolver->portDimension()
                << ", ritz residual="
                << coarseDiagnostics.maximumRitzResidual
                << ", orthogonality="
                << coarseDiagnostics
                    .maximumLocalCoarseOrthogonality
                << ", setup=" << coarseDiagnostics.basisSeconds
                << " s\n";
        }
        if (options.globalRandomizedSchur) {
            if (!preparedPortDynamicModel || !portModel) {
                throw std::runtime_error(
                    "[Global randomized] Dynamic and local-port models "
                    "are required.");
            }
            if (!portSolver) {
                portSolver =
                    std::make_unique<LocalPortReducedSchurSolver>(
                        *preparedPortDynamicModel, *portModel);
            }
            const auto globalBuildStart = Clock::now();
            const std::size_t baselineWorkingSet =
                globalCoarseCurrentWorkingSetBytes();
            ReducedDynamicSchurOperator exactGlobalSchur(
                *preparedPortDynamicModel, true);
            std::vector<double> globalInterfaceMetric(
                partition.interfaceGlobalDofs.size(), 0.0);
            for (std::size_t interfaceRow = 0;
                 interfaceRow < partition.interfaceGlobalDofs.size();
                 ++interfaceRow) {
                globalInterfaceMetric[interfaceRow] =
                    descriptor.interfacePenaltyMassDiagonal[
                        static_cast<std::size_t>(
                            partition.interfaceGlobalDofs[
                                interfaceRow])];
            }
            GlobalRandomizedSchurOptions randomizedOptions;
            randomizedOptions.requestedRank =
                options.globalRandomizedRank;
            randomizedOptions.seed =
                options.globalRandomizedSeed;
            randomizedOptions.innerMaximumIterations =
                options.globalRandomizedInnerMaximumIterations;
            randomizedOptions.innerTolerance =
                options.globalRandomizedInnerTolerance;
            randomizedOptions.orthogonalityTolerance = 1.0e-10;
            randomizedOptions.deflationTolerance =
                options.rankTolerance;
            randomizedOptions.maximumBasisTimeSeconds = 600.0;
            randomizedOptions.composition =
                options.globalRandomizedComposition;
            GlobalRandomizedSchurResult randomized =
                buildGlobalRandomizedSchurPortSpace(
                    *portModel, *portSolver, exactGlobalSchur,
                    globalInterfaceMetric, randomizedOptions);
            randomized.diagnostics.basisBuildTimeSeconds =
                elapsed(globalBuildStart);
            randomized.diagnostics.baselineWorkingSetBytes =
                baselineWorkingSet;
            randomized.diagnostics.peakWorkingSetBytes = std::max(
                randomized.diagnostics.peakWorkingSetBytes,
                globalCoarseCurrentWorkingSetBytes());
            randomized.diagnostics.peakIncrementalMemoryBytes =
                randomized.diagnostics.peakWorkingSetBytes
                    > baselineWorkingSet
                ? randomized.diagnostics.peakWorkingSetBytes
                    - baselineWorkingSet
                : 0;
            const bool finalGate =
                randomized.model.rank
                    == options.globalRandomizedRank
                && randomized.diagnostics.orthogonalityError
                    <= 1.0e-10
                && randomized.diagnostics.schurResidual
                    <= 1.0e-8
                && randomized.diagnostics
                    .maximumTargetSolveResidual
                    <= options.globalRandomizedInnerTolerance
                && randomized.diagnostics
                    .peakIncrementalMemoryBytes
                    < UINT64_C(1024) * 1024 * 1024
                && randomized.diagnostics
                    .basisBuildTimeSeconds < 600.0;
            if (finalGate) {
                randomized.diagnostics.status = "passed";
            } else if (randomized.diagnostics.status == "passed") {
                randomized.diagnostics.status =
                    "global_randomized_port_failed";
                if (randomized.diagnostics
                        .basisBuildTimeSeconds >= 600.0) {
                    randomized.diagnostics.failureReason =
                        "basis_build_time_gate_failed";
                } else if (randomized.diagnostics
                               .peakIncrementalMemoryBytes
                           >= UINT64_C(1024) * 1024 * 1024) {
                    randomized.diagnostics.failureReason =
                        "incremental_memory_gate_failed";
                } else {
                    randomized.diagnostics.failureReason =
                        "combined_algebra_resource_gate_failed";
                }
            }
            writeGlobalRandomizedSchurDiagnostics(
                randomized, physics.name, outputDirectory);
            if (!finalGate) {
                std::cout
                    << "[Global randomized] Algebra/resource gate "
                    << "failed with status="
                    << randomized.diagnostics.status
                    << "; single-step validation skipped.\n";
                return;
            }
            globalRandomizedDiagnostics =
                std::make_unique<
                    GlobalRandomizedSchurDiagnostics>(
                        std::move(randomized.diagnostics));
            globalCoarseModel =
                std::make_unique<GlobalInterfaceCoarseModel>(
                    std::move(randomized.model));
            portSolver->attachGlobalCoarse(
                *globalCoarseModel,
                options.globalRandomizedComposition
                    == "augment-local");
            std::cout
                << "[Global randomized] composition="
                << options.globalRandomizedComposition
                << ", local/global/active="
                << portModel->reducedInterfaceDofs << '/'
                << globalCoarseModel->rank << '/'
                << portSolver->portDimension()
                << ", orthogonality="
                << globalRandomizedDiagnostics
                    ->orthogonalityError
                << ", solve residual="
                << globalRandomizedDiagnostics
                    ->maximumTargetSolveResidual
                << ", setup="
                << globalRandomizedDiagnostics
                    ->basisBuildTimeSeconds << " s\n";
        }
    }

    const double postLocalSetupSeconds = elapsed(postLocalSetupStart);
    const auto schurStart = Clock::now();
    const bool matrixFreeDynamicSchur = !options.portReduction &&
        static_cast<int>(partition.interfaceGlobalDofs.size())
            > options.matrixFreeInterfaceThreshold;
    bool nativeReducedHistoryAligned = matrixFreeDynamicSchur
        && options.nativeReducedHistory;
    std::vector<DynamicLocal> dynamic;
    std::vector<MatrixEntry> schurMatrixEntries;
    std::unique_ptr<SubdomainDirectSolver> schurFactor;
    std::unique_ptr<local::Model> matrixFreeModel =
        std::move(preparedPortDynamicModel);
    std::unique_ptr<local::LocalReducedSchurSolver> matrixFreeSolver;
    double schurFactorSeconds = 0.0;
    if (options.portReduction) {
        if (!matrixFreeModel) {
            matrixFreeModel = std::make_unique<local::Model>(
                makeDynamicReducedModel(
                    descriptor.dofs, partition, k, c, locals,
                    options.timeStep));
        }
        if (!portModel) {
            throw std::runtime_error(
                "[Local port] Interface basis was not constructed.");
        }
        if (!portSolver) {
            portSolver =
                std::make_unique<LocalPortReducedSchurSolver>(
                    *matrixFreeModel, *portModel);
        }
        schurFactorSeconds = portSolver->factorizationSeconds();
    } else if (matrixFreeDynamicSchur) {
        matrixFreeModel = std::make_unique<local::Model>(makeDynamicReducedModel(
            descriptor.dofs, partition, k, c, locals, options.timeStep));
        ddm_schur::Options schurOptions;
        schurOptions.maxIterations = options.interfaceMaxIterations;
        schurOptions.restart = options.interfaceRestart;
        schurOptions.relativeTolerance = options.interfaceTolerance;
        schurOptions.coarseLinearXY = options.coarseLinearXY;
        schurOptions.coarseLinearZ = options.coarseLinearZ;
        schurOptions.proxyEnabled = options.proxyEnabled;
        schurOptions.proxyDisableCoarse = options.proxyDisableCoarse;
        schurOptions.proxyHighConductivityThreshold =
            options.proxyHighConductivityThreshold;
        schurOptions.proxyUseMaterialConnectivity =
            options.proxyUseMaterialConnectivity;
        schurOptions.proxyRing = options.proxyRing;
        schurOptions.proxyProbeColumns = options.proxyProbeColumns;
        schurOptions.proxyBlockSize = options.proxyBlockSize;
        schurOptions.proxyValidateBlockEquivalence =
            options.proxyValidateBlockEquivalence;
        schurOptions.localSolveThreads = options.localSolveThreads;
        // Local Arnoldi factors benefit from one MKL thread per concurrently
        // built subdomain.  The augmented interface system is one global
        // sparse factor, so let it use the full outer thread budget instead.
        schurOptions.localPardisoThreads =
            options.interfaceKrylov == "augmented-direct"
            ? std::max(options.localPardisoThreads, options.localSolveThreads)
            : options.localPardisoThreads;
        schurOptions.interfaceKrylov = options.interfaceKrylov;
        schurOptions.portCoreCacheEnabled = options.portCoreCacheEnabled;
        schurOptions.portCoreCachePath =
            options.portCoreCachePath.string();
        schurOptions.proxyCacheEnabled = options.proxyCacheEnabled;
        schurOptions.proxyCachePath = options.proxyCachePath.string();
        schurOptions.interfaceOperatorCoarseRank =
            options.interfaceOperatorCoarseRank;
        schurOptions.interfaceOperatorCoarseSweeps =
            options.interfaceOperatorCoarseSweeps;
        schurOptions.interfaceOperatorCoarsePredictor =
            options.interfaceOperatorCoarsePredictor;
        schurOptions.interfaceOperatorCoarseCachePath =
            options.interfaceOperatorCoarseCachePath.string();
        schurOptions.proxyOutputDirectory = outputDirectory.string();
        matrixFreeSolver = std::make_unique<local::LocalReducedSchurSolver>(
            *matrixFreeModel, mesh, physics, partition, schurOptions,
            options.matrixFreeInterfaceThreshold, outputDirectory);
        schurFactorSeconds = matrixFreeSolver->factorizationSeconds();
    } else {
        dynamic.reserve(locals.size());
        std::map<std::pair<int, int>, double> schurEntries;
        k.interfaceBlock.forEachEntry([&](int row, int column, double value) {
            const auto key = std::minmax(row, column);
            schurEntries[{key.first, key.second}] += 0.5 * value;
        });
        c.interfaceBlock.forEachEntry([&](int row, int column, double value) {
            const auto key = std::minmax(row, column);
            schurEntries[{key.first, key.second}] += 0.5 * value / options.timeStep;
        });
        // The 0.5 accumulation sees both halves. Restore the half diagonal.
        k.interfaceBlock.forEachEntry([&](int row, int column, double value) {
            if (row == column) schurEntries[{row, column}] += 0.5 * value;
        });
        c.interfaceBlock.forEachEntry([&](int row, int column, double value) {
            if (row == column) schurEntries[{row, column}] += 0.5 * value / options.timeStep;
        });
        for (const LocalModel& localModel : locals) {
            DynamicLocal localDynamic;
            localDynamic.model = &localModel;
            localDynamic.aii.resize(localModel.kii.size());
            for (std::size_t entry = 0; entry < localDynamic.aii.size(); ++entry) {
                localDynamic.aii[entry] = localModel.kii[entry]
                    + localModel.cii[entry] / options.timeStep;
            }
            localDynamic.aiGamma.resize(localModel.kiGamma.size());
            for (std::size_t entry = 0; entry < localDynamic.aiGamma.size(); ++entry) {
                localDynamic.aiGamma[entry] = localModel.kiGamma[entry]
                    + localModel.ciGamma[entry] / options.timeStep;
            }
            localDynamic.aGammaI.resize(localModel.kGammaI.size());
            for (std::size_t entry = 0; entry < localDynamic.aGammaI.size(); ++entry) {
                localDynamic.aGammaI[entry] = localModel.kGammaI[entry]
                    + localModel.cGammaI[entry] / options.timeStep;
            }
            localDynamic.factor = local::factorDenseSymmetric(
                localDynamic.aii, localModel.rank);
            localDynamic.solvedCoupling.assign(static_cast<std::size_t>(localModel.rank)
                * localModel.gammaDofs, 0.0);
            for (int column = 0; column < localModel.gammaDofs; ++column) {
                std::vector<double> rhs(static_cast<std::size_t>(localModel.rank), 0.0);
                for (int mode = 0; mode < localModel.rank; ++mode) {
                    rhs[static_cast<std::size_t>(mode)] = localDynamic.aiGamma[
                        static_cast<std::size_t>(mode * localModel.gammaDofs + column)];
                }
                local::solveDenseSymmetric(localDynamic.factor, rhs);
                for (int mode = 0; mode < localModel.rank; ++mode) {
                    localDynamic.solvedCoupling[static_cast<std::size_t>(
                        mode * localModel.gammaDofs + column)] = rhs[static_cast<std::size_t>(mode)];
                }
            }
            for (int localRow = 0; localRow < localModel.gammaDofs; ++localRow) {
                for (int localColumn = localRow; localColumn < localModel.gammaDofs; ++localColumn) {
                    long double update = 0.0L;
                    for (int mode = 0; mode < localModel.rank; ++mode) {
                        update += static_cast<long double>(localDynamic.aGammaI[
                            static_cast<std::size_t>(localRow * localModel.rank + mode)])
                            * localDynamic.solvedCoupling[static_cast<std::size_t>(
                                mode * localModel.gammaDofs + localColumn)];
                    }
                    const int globalRow = localModel.gammaIndices[static_cast<std::size_t>(localRow)];
                    const int globalColumn = localModel.gammaIndices[static_cast<std::size_t>(localColumn)];
                    const auto key = std::minmax(globalRow, globalColumn);
                    schurEntries[{key.first, key.second}] -= static_cast<double>(update);
                }
            }
            dynamic.push_back(std::move(localDynamic));
        }
        schurMatrixEntries.reserve(schurEntries.size());
        for (const auto& entry : schurEntries) {
            if (std::abs(entry.second) > 0.0) {
                schurMatrixEntries.push_back(
                    {entry.first.first, entry.first.second, entry.second});
            }
        }
        const auto schurFactorStart = Clock::now();
        schurFactor = std::make_unique<SubdomainDirectSolver>(
            static_cast<int>(partition.interfaceGlobalDofs.size()), schurMatrixEntries);
        schurFactorSeconds = elapsed(schurFactorStart);
    }
    const double schurSetupSeconds = elapsed(schurStart);
    const auto postSchurSetupStart = Clock::now();

    std::unique_ptr<SparseMatrix> fullStep;
    double fullStepAssemblySeconds = 0.0;
    auto ensureFullStep = [&]() -> SparseMatrix& {
        if (!fullStep) {
            const auto assemblyStart = Clock::now();
            fullStep = std::make_unique<SparseMatrix>(descriptor.dofs);
            fullStep->appendScaledEntries(
                descriptor.capacity, 1.0 / options.timeStep);
            fullStep->appendScaledEntries(descriptor.conductivity, 1.0);
            fullStep->finalizeCsr();
            fullStepAssemblySeconds += elapsed(assemblyStart);
        }
        return *fullStep;
    };
    auto applyFullStep = [&](const std::vector<double>& state) {
        if (fullStep) {
            return fullStep->multiply(state);
        }
        std::vector<double> result = descriptor.capacity.multiply(state);
        const std::vector<double> conductivityImage =
            descriptor.conductivity.multiply(state);
        for (std::size_t row = 0; row < result.size(); ++row) {
            result[row] = result[row] / options.timeStep
                + conductivityImage[row];
        }
        return result;
    };
    if (reducedValidation) {
        const SourceAlignedValidationReport physicalValidation =
            runSourceAlignedInterfaceValidation(
            mesh, physics, descriptor, k, c, partition, locals, thetaInitial,
            boundaryOffset, options, ensureFullStep(), reducedValidationReport,
            outputDirectory);
        if (!(physicalValidation.maximumPhysicalOperatorError < 1.0e-4)) {
            throw std::runtime_error(
                "Reduced Dynamic Schur physical operator validation failed: "
                "FOM transient interface trajectory error exceeds 1e-4.");
        }
        if (options.sourceAlignedInterfaceValidation) return;
    }
    Accuracy accuracy;
    accuracy.dynamicSchurFactorBytes = options.portReduction
        ? portSolver->factorMemoryBytes()
        : (matrixFreeDynamicSchur
            ? matrixFreeSolver->factorMemoryBytes() : schurFactor->memoryBytes());
    std::unique_ptr<SubdomainDirectSolver> fomFactor;
    const auto ensureFomFactor = [&]() -> SubdomainDirectSolver& {
        if (!fomFactor) {
            const auto factorStart = Clock::now();
            fomFactor = std::make_unique<SubdomainDirectSolver>(
                descriptor.dofs, sparseMatrixEntries(ensureFullStep()));
            accuracy.fomFactorSeconds += elapsed(factorStart);
            accuracy.fomFactorBytes = fomFactor->memoryBytes();
            accuracy.factorBytes = accuracy.fomFactorBytes
                + accuracy.dynamicSchurFactorBytes;
        }
        return *fomFactor;
    };
    accuracy.factorBytes = accuracy.dynamicSchurFactorBytes;
    if (options.compareFom || options.localPortCorrected) {
        ensureFomFactor();
    }
    std::vector<double> fomTheta = thetaInitial;
    std::vector<double> previousRomTheta = thetaInitial;
    std::vector<double> matrixFreeRomTheta = thetaInitial;
    long double totalErrorSquared = 0.0L;
    long double totalReferenceSquared = 0.0L;
    long double correctedErrorSquared = 0.0L;
    long double correctedReferenceSquared = 0.0L;
    std::ofstream byTime(outputDirectory / "local_dynamic_schur_accuracy_by_time.csv");
    byTime << "waveform,initial_mode,step,time_s,relative_l2,maximum_absolute_k,"
        "fom_maximum_k,local_maximum_k,maximum_temperature_error_k,full_residual,"
        "full_residual_before_gate,full_residual_tolerance,"
        "residual_fallback_enabled,residual_gate_passed,"
        "residual_fallback_used,full_residual_after_fallback,"
        "reduced_residual,temperature_jump_rms_k,relative_flux_imbalance,"
        "fom_rom_flux_relative_l2,worst_physical_interface,worst_face_pair,"
        "worst_integration_triangle,interface_iterations,interface_matvecs,"
        "interface_true_relative_residual,port_projection_relative_error\n"
        << std::setprecision(17);
    std::ofstream reducedTiming;
    if (reducedMethod || matrixFreeDynamicSchur || options.portReduction) {
        reducedTiming.open(outputDirectory / "local_dynamic_schur_reduced_timing.csv");
        reducedTiming << "step,time_s,interface_initial_guess,"
            "interface_krylov_actual,interface_krylov_fallback,"
            "interface_krylov_fallback_reason,"
            "adaptive_tolerance_used,adaptive_retry,adaptive_retry_seconds,"
            "interface_initial_relative_residual,interface_predictor_applied,"
            "interface_predictor_accepted,interface_predictor_initial_relative_residual,"
            "interface_predictor_relative_residual,interface_predictor_seconds,"
            "local_reduced_solve_seconds,"
            "native_reduced_rhs_seconds,step_rhs_seconds,"
            "interface_solve_seconds,proxy_solve_seconds,"
            "coarse_solve_seconds,port_forward_solve_seconds,"
            "port_core_solve_seconds,port_back_substitution_seconds,"
            "schur_apply_seconds,preconditioner_seconds,"
            "orthogonalization_seconds,vector_update_seconds,"
            "full_residual_seconds,state_change_relative,"
            "total_timestep_seconds,interface_iterations,"
            "full_residual\n"
            << std::setprecision(17);
    }
    std::ofstream fluxOut;
    if (options.compareFom && !options.compareFomSummaryOnly) {
        fluxOut.open(outputDirectory / "local_dynamic_schur_interface_flux.csv");
        fluxOut << "waveform,initial_mode,step,time_s,physical_interface_id,face_pair_id,"
            "integration_triangle_id,left_subdomain,right_subdomain,left_boundary_entity,"
            "right_boundary_entity,area_m2,rom_temperature_jump_rms_k,"
            "rom_left_physical_normal_flux_w_m2,rom_right_physical_normal_flux_w_m2,"
            "rom_sipg_numerical_flux_w_m2,rom_flux_imbalance_l2_w_m2,"
            "rom_relative_flux_imbalance,fom_temperature_jump_rms_k,"
            "fom_left_physical_normal_flux_w_m2,fom_right_physical_normal_flux_w_m2,"
            "fom_sipg_numerical_flux_w_m2,fom_flux_imbalance_l2_w_m2,"
            "fom_relative_flux_imbalance,sipg_flux_error_w_m2,sipg_flux_relative_error\n"
            << std::setprecision(17);
    }
    std::ofstream adaptivePhysicalDiagnostics;
    if (options.milestone8AdaptiveProduction) {
        adaptivePhysicalDiagnostics.open(
            outputDirectory
            / "milestone8_adaptive_port_physical_diagnostics.csv");
        adaptivePhysicalDiagnostics
            << "step,time_s,interface_id,history_rank,"
            "randomized_transfer_rank,residual_rank,total_port_rank,"
            "target_dofs,temperature_relative_l2,max_nodal_error_k,"
            "local_full_residual,flux_relative_l2,"
            "temperature_jump_rms_k,max_relative_flux_imbalance,"
            "corrected\n" << std::setprecision(17);
    }
    std::map<std::pair<int, int>, int>
        adaptivePhysicalInterfaceBySubdomains;
    if (options.milestone8AdaptiveProduction) {
        for (const LocalPortBasis& port : portModel->ports) {
            const auto adjacent = std::minmax(
                port.leftSubdomain, port.rightSubdomain);
            const auto inserted =
                adaptivePhysicalInterfaceBySubdomains.emplace(
                    std::make_pair(adjacent.first, adjacent.second),
                    port.interfaceId);
            if (!inserted.second
                && inserted.first->second != port.interfaceId) {
                throw std::runtime_error(
                    "[M8.9 adaptive] Adjacent subdomain pair has "
                    "multiple physical-port ids.");
            }
        }
    }
    std::vector<double> finalRom;
    std::vector<double> gammaTwoStepsAgo = gamma;
    struct ExactSteadyStateReuse {
        bool active = false;
        std::vector<double> rom;
        double romMaximum = 0.0;
        double fullResidual = 0.0;
        double fullResidualBeforeGate = 0.0;
        double fullResidualAfterFallback = 0.0;
        double reducedResidual = 0.0;
        InterfacePhysicsMetrics interfaceMetrics;
    } steadyReuse;
    int steadyStateDetectedStep = -1;
    int steadyStateReusedSteps = 0;
    double lastStateChangeRelative =
        std::numeric_limits<double>::quiet_NaN();
    const double postSchurSetupSeconds = elapsed(postSchurSetupStart);
    const auto timeSteppingStart = Clock::now();

    for (int step = 0; step <= steps; ++step) {
        const auto stepStart = Clock::now();
        double stepLocalSolveSeconds = 0.0;
        double stepInterfaceSolveSeconds = 0.0;
        double stepProxySolveSeconds = 0.0;
        double stepCoarseSolveSeconds = 0.0;
        double stepPortForwardSolveSeconds = 0.0;
        double stepPortCoreSolveSeconds = 0.0;
        double stepPortBackSubstitutionSeconds = 0.0;
        double stepSchurApplySeconds = 0.0;
        double stepPreconditionerSeconds = 0.0;
        double stepOrthogonalizationSeconds = 0.0;
        double stepVectorUpdateSeconds = 0.0;
        bool stepAdaptiveToleranceUsed = false;
        bool stepAdaptiveRetry = false;
        double stepAdaptiveRetrySeconds = 0.0;
        struct AdaptivePhysicalInterfaceRow {
            double temperatureRelativeL2 = 0.0;
            double maximumNodalError = 0.0;
            double localFullResidual = 0.0;
            long double fluxErrorSquared = 0.0L;
            long double fluxReferenceSquared = 0.0L;
            long double jumpAreaWeightedSquared = 0.0L;
            long double area = 0.0L;
            double maximumRelativeFluxImbalance = 0.0;
        };
        std::map<int, AdaptivePhysicalInterfaceRow>
            adaptivePhysicalRows;
        const double time = step * options.timeStep;
        int stepInterfaceIterations = 0;
        int stepInterfaceMatvecs = 0;
        double stepInterfaceInitialResidual = 0.0;
        double stepInterfaceResidual = 0.0;
        bool stepInterfacePredictorApplied = false;
        bool stepInterfacePredictorAccepted = false;
        double stepInterfacePredictorInitialResidual = 0.0;
        double stepInterfacePredictorResidual = 0.0;
        double stepInterfacePredictorSeconds = 0.0;
        std::string stepInterfaceKrylovActual = "not_applicable";
        bool stepInterfaceKrylovFallback = false;
        std::string stepInterfaceKrylovFallbackReason;
        double stepPortProjectionError = 0.0;
        double stepNativeReducedRhsSeconds = 0.0;
        bool stepUsedNativeReducedHistory = false;
        double stepRhsSeconds = 0.0;
        double stepResidualSeconds = 0.0;
        double stepStateChangeRelative =
            std::numeric_limits<double>::quiet_NaN();
        bool stepStateExactlyUnchanged = false;
        std::vector<double> stepPowers;
        std::vector<double> stepRhs;
        std::vector<std::vector<double>> nativeInteriorRhs;
        std::vector<double> nativeInterfaceRhs;
        if (step > 0 && steadyReuse.active) {
            ++steadyStateReusedSteps;
            const double unavailable =
                std::numeric_limits<double>::quiet_NaN();
            byTime << waveform.name << ',' << options.initialMode << ','
                << step << ',' << time << ',' << unavailable << ','
                << unavailable << ',' << unavailable << ','
                << steadyReuse.romMaximum << ',' << unavailable << ','
                << steadyReuse.fullResidual << ','
                << steadyReuse.fullResidualBeforeGate << ','
                << options.fullResidualTolerance << ','
                << (options.fullResidualFallback ? 1 : 0) << ",1,0,"
                << steadyReuse.fullResidualAfterFallback << ','
                << steadyReuse.reducedResidual << ','
                << steadyReuse.interfaceMetrics.temperatureJumpRms << ','
                << steadyReuse.interfaceMetrics.relativeFluxImbalance << ','
                << unavailable << ",-1,-1,-1,0,0,0,0\n";
            if (reducedMethod || matrixFreeDynamicSchur
                || options.portReduction) {
                reducedTiming << step << ',' << time
                    << ",steady_exact_reuse,not_run,0,none,0,0,0,0,0,0,0,0,0,"
                    << "0,0,0,0,0,0,0,0,0,0,0,"
                    << "0," << elapsed(stepStart) << ",0,"
                    << steadyReuse.fullResidual << '\n';
            }
            lastStateChangeRelative = 0.0;
            if (step == steps) finalRom = steadyReuse.rom;
            continue;
        }
        if (step > 0) {
            stepPowers = waveform.sample(time);
            if (matrixFreeDynamicSchur && options.nativeReducedHistory
                && nativeReducedHistoryAligned) {
                const auto nativeRhsStart = Clock::now();
                NativeReducedHistoryRhs nativeRhs = buildNativeReducedHistoryRhs(
                    descriptor, c, partition, locals, localStates, gamma,
                    stepPowers, boundaryOffset, options.timeStep);
                nativeInteriorRhs = std::move(nativeRhs.interior);
                nativeInterfaceRhs = std::move(nativeRhs.interfaceRhs);
                stepNativeReducedRhsSeconds = elapsed(nativeRhsStart);
                stepUsedNativeReducedHistory = true;
                accuracy.nativeReducedRhsSeconds += stepNativeReducedRhsSeconds;
                ++accuracy.nativeReducedHistorySteps;
            }
            const auto rhsStart = Clock::now();
            stepRhs = descriptor.capacity.multiply(previousRomTheta);
            for (double& value : stepRhs) value /= options.timeStep;
            addInput(descriptor, stepPowers, stepRhs);
            for (std::size_t row = 0; row < stepRhs.size(); ++row) {
                stepRhs[row] += boundaryOffset[row];
            }
            stepRhsSeconds = elapsed(rhsStart);
            if (options.compareFom) {
                const auto fomSolveStart = Clock::now();
                std::vector<double> fomRhs = descriptor.capacity.multiply(fomTheta);
                for (double& value : fomRhs) value /= options.timeStep;
                addInput(descriptor, stepPowers, fomRhs);
                for (int row = 0; row < descriptor.dofs; ++row) {
                    fomRhs[static_cast<std::size_t>(row)] +=
                        boundaryOffset[static_cast<std::size_t>(row)];
                }
                std::vector<double> fomNext;
                ensureFomFactor().solve(fomRhs, fomNext);
                fomTheta = std::move(fomNext);
                accuracy.fomSolveSeconds += elapsed(fomSolveStart);
            }
            if (options.compareFom && options.portReduction) {
                std::vector<double> fomInterface(
                    partition.interfaceGlobalDofs.size(), 0.0);
                for (std::size_t gammaIndex = 0;
                     gammaIndex < partition.interfaceGlobalDofs.size(); ++gammaIndex) {
                    fomInterface[gammaIndex] = fomTheta[static_cast<std::size_t>(
                        partition.interfaceGlobalDofs[gammaIndex])];
                }
                stepPortProjectionError =
                    portSolver->relativeProjectionError(fomInterface);
                accuracy.maximumPortProjectionError = std::max(
                    accuracy.maximumPortProjectionError,
                    stepPortProjectionError);
            }

            const auto localCoreStart = Clock::now();
            if (matrixFreeDynamicSchur || options.portReduction) {
                local::SolveResult solve;
                double portCondensedSeconds = 0.0;
                if (options.portReduction) {
                    LocalPortSolveResult portSolve = portSolver->solve(stepRhs);
                    portCondensedSeconds = portSolve.condensedRhsSeconds;
                    accuracy.maximumPortReducedResidual = std::max(
                        accuracy.maximumPortReducedResidual,
                        portSolve.reducedRelativeResidual);
                    solve = std::move(portSolve.solution);
                } else {
                    const std::vector<double> gammaBeforeSolve = gamma;
                    std::vector<double> extrapolatedGamma;
                    const std::vector<double>* initialGuess = nullptr;
                    if (options.interfaceInitialGuess == "previous") {
                        initialGuess = &gamma;
                    } else if (options.interfaceInitialGuess == "extrapolated") {
                        if (step > 1) {
                            extrapolatedGamma.resize(gamma.size(), 0.0);
                            for (std::size_t row = 0; row < gamma.size(); ++row) {
                                extrapolatedGamma[row] = 2.0 * gamma[row]
                                    - gammaTwoStepsAgo[row];
                            }
                            initialGuess = &extrapolatedGamma;
                        } else {
                            initialGuess = &gamma;
                        }
                    }
                    stepAdaptiveToleranceUsed = !options.compareFom
                        && options.interfaceKrylov == "fgmres"
                        && options.adaptiveInterfaceTolerance
                            > options.interfaceTolerance;
                    solve = stepUsedNativeReducedHistory
                        ? matrixFreeSolver->solveReducedRhs(
                            nativeInteriorRhs, nativeInterfaceRhs, initialGuess,
                            stepAdaptiveToleranceUsed
                                ? options.adaptiveInterfaceTolerance : 0.0)
                        : matrixFreeSolver->solve(
                            stepRhs, initialGuess,
                            stepAdaptiveToleranceUsed
                                ? options.adaptiveInterfaceTolerance : 0.0);
                    gammaTwoStepsAgo = gammaBeforeSolve;
                }
                if (solve.status != "success") {
                    throw std::runtime_error(
                        "[Local transient] Dynamic Schur interface solve failed: "
                        + solve.status);
                }
                matrixFreeRomTheta = std::move(solve.temperature);
                gamma = std::move(solve.interfaceTemperature);
                if (matrixFreeDynamicSchur && options.nativeReducedHistory) {
                    if (solve.localReducedCoordinates.size() == locals.size()) {
                        localStates = std::move(solve.localReducedCoordinates);
                        nativeReducedHistoryAligned = true;
                    } else {
                        nativeReducedHistoryAligned = false;
                    }
                }
                accuracy.interfaceSolveSeconds += solve.timing.interfaceSolveSeconds;
                accuracy.interfaceSolveSeconds += portCondensedSeconds;
                accuracy.recoverySeconds += solve.timing.localRecoverySeconds
                    + solve.timing.fullFieldReconstructionSeconds;
                accuracy.proxySolveSeconds += solve.timing.proxySolveSeconds;
                accuracy.coarseSolveSeconds += solve.timing.coarseSolveSeconds;
                accuracy.portForwardSolveSeconds +=
                    solve.timing.portForwardSolveSeconds;
                accuracy.portCoreSolveSeconds +=
                    solve.timing.portCoreSolveSeconds;
                accuracy.portBackSubstitutionSeconds +=
                    solve.timing.portBackSubstitutionSeconds;
                accuracy.interfaceOperatorSeconds +=
                    solve.timing.interfaceOperatorSeconds;
                accuracy.interfacePreconditionerSeconds +=
                    solve.timing.interfacePreconditionerSeconds;
                accuracy.interfaceOrthogonalizationSeconds +=
                    solve.timing.interfaceOrthogonalizationSeconds;
                accuracy.interfaceVectorUpdateSeconds +=
                    solve.timing.interfaceVectorUpdateSeconds;
                accuracy.interfacePredictorSeconds +=
                    solve.timing.interfacePredictorSeconds;
                if (solve.timing.interfacePredictorApplied) {
                    ++accuracy.interfacePredictorAppliedSteps;
                }
                if (solve.timing.interfacePredictorAccepted) {
                    ++accuracy.interfacePredictorAcceptedSteps;
                }
                accuracy.interfaceIterationsTotal += solve.timing.interfaceIterations;
                accuracy.interfaceIterationsMaximum = std::max(
                    accuracy.interfaceIterationsMaximum,
                    solve.timing.interfaceIterations);
                accuracy.interfaceMatvecs += solve.timing.interfaceMatvecs;
                accuracy.interfaceKrylovActual =
                    solve.timing.interfaceKrylovActual;
                if (solve.timing.interfaceKrylovFallback) {
                    ++accuracy.interfaceKrylovFallbackSteps;
                }
                accuracy.maximumInterfaceResidual = std::max(
                    accuracy.maximumInterfaceResidual,
                    solve.timing.interfaceRelativeResidual);
                stepInterfaceIterations = solve.timing.interfaceIterations;
                stepInterfaceMatvecs = solve.timing.interfaceMatvecs;
                stepInterfaceInitialResidual =
                    solve.timing.interfaceInitialRelativeResidual;
                stepInterfaceResidual = solve.timing.interfaceRelativeResidual;
                stepInterfacePredictorApplied =
                    solve.timing.interfacePredictorApplied;
                stepInterfacePredictorAccepted =
                    solve.timing.interfacePredictorAccepted;
                stepInterfacePredictorInitialResidual =
                    solve.timing.interfacePredictorInitialRelativeResidual;
                stepInterfacePredictorResidual =
                    solve.timing.interfacePredictorRelativeResidual;
                stepInterfacePredictorSeconds =
                    solve.timing.interfacePredictorSeconds;
                stepInterfaceKrylovActual =
                    solve.timing.interfaceKrylovActual;
                stepInterfaceKrylovFallback =
                    solve.timing.interfaceKrylovFallback;
                stepInterfaceKrylovFallbackReason =
                    solve.timing.interfaceKrylovFallbackReason;
                stepLocalSolveSeconds = solve.timing.localReducedAssemblySeconds;
                stepInterfaceSolveSeconds = solve.timing.interfaceSolveSeconds
                    + portCondensedSeconds;
                stepProxySolveSeconds = solve.timing.proxySolveSeconds;
                stepCoarseSolveSeconds = solve.timing.coarseSolveSeconds;
                stepPortForwardSolveSeconds =
                    solve.timing.portForwardSolveSeconds;
                stepPortCoreSolveSeconds = solve.timing.portCoreSolveSeconds;
                stepPortBackSubstitutionSeconds =
                    solve.timing.portBackSubstitutionSeconds;
                stepSchurApplySeconds = solve.timing.interfaceOperatorSeconds;
                stepPreconditionerSeconds =
                    solve.timing.interfacePreconditionerSeconds;
                stepOrthogonalizationSeconds =
                    solve.timing.interfaceOrthogonalizationSeconds;
                stepVectorUpdateSeconds =
                    solve.timing.interfaceVectorUpdateSeconds;
            } else {
            std::vector<std::vector<double>> localRhs(locals.size());
            std::vector<double> gammaRhs(static_cast<std::size_t>(
                partition.interfaceGlobalDofs.size()), 0.0);
            const std::vector<double> cGamma = c.interfaceBlock.multiply(gamma);
            const auto reducedLocalStart = Clock::now();
            for (std::size_t row = 0; row < gammaRhs.size(); ++row) {
                const int global = partition.interfaceGlobalDofs[row];
                gammaRhs[row] = cGamma[row] / options.timeStep
                    + boundaryOffset[static_cast<std::size_t>(global)];
                for (int channel = 0; channel < descriptor.sourceChannels; ++channel) {
                    gammaRhs[row] += stepPowers[static_cast<std::size_t>(channel)]
                        * descriptor.input[static_cast<std::size_t>(channel * descriptor.dofs + global)];
                }
            }
            for (std::size_t slot = 0; slot < locals.size(); ++slot) {
                const LocalModel& localModel = locals[slot];
                std::vector<double> rhs;
                denseMatvec(localModel.cii, localModel.rank, localModel.rank,
                    localStates[slot], rhs, 1.0 / options.timeStep, 0.0);
                std::vector<double> localGamma(static_cast<std::size_t>(localModel.gammaDofs), 0.0);
                for (int row = 0; row < localModel.gammaDofs; ++row) {
                    localGamma[static_cast<std::size_t>(row)] = gamma[static_cast<std::size_t>(
                        localModel.gammaIndices[static_cast<std::size_t>(row)])];
                }
                denseMatvec(localModel.ciGamma, localModel.rank,
                    localModel.gammaDofs, localGamma, rhs,
                    1.0 / options.timeStep, 1.0);
                addReducedInput(localModel, stepPowers, rhs);
                for (int mode = 0; mode < localModel.rank; ++mode) {
                    rhs[static_cast<std::size_t>(mode)] +=
                        localModel.reducedBoundary[static_cast<std::size_t>(mode)];
                }
                std::vector<double> gammaHistory;
                denseMatvec(localModel.cGammaI, localModel.gammaDofs,
                    localModel.rank, localStates[slot], gammaHistory,
                    1.0 / options.timeStep, 0.0);
                for (int row = 0; row < localModel.gammaDofs; ++row) {
                    gammaRhs[static_cast<std::size_t>(localModel.gammaIndices[
                        static_cast<std::size_t>(row)])] += gammaHistory[static_cast<std::size_t>(row)];
                }
                localRhs[slot] = std::move(rhs);
            }
            for (std::size_t slot = 0; slot < dynamic.size(); ++slot) {
                std::vector<double> solved = localRhs[slot];
                local::solveDenseSymmetric(dynamic[slot].factor, solved);
                std::vector<double> image;
                denseMatvec(dynamic[slot].aGammaI, locals[slot].gammaDofs,
                    locals[slot].rank, solved, image);
                for (int row = 0; row < locals[slot].gammaDofs; ++row) {
                    gammaRhs[static_cast<std::size_t>(locals[slot].gammaIndices[
                        static_cast<std::size_t>(row)])] -= image[static_cast<std::size_t>(row)];
                }
            }
            stepLocalSolveSeconds = elapsed(reducedLocalStart);
            const auto interfaceStart = Clock::now();
            std::vector<double> nextGamma;
            schurFactor->solve(gammaRhs, nextGamma);
            stepInterfaceSolveSeconds = elapsed(interfaceStart);
            accuracy.interfaceSolveSeconds += stepInterfaceSolveSeconds;
            const auto recoveryStart = Clock::now();
            std::vector<std::vector<double>> nextStates(locals.size());
            for (std::size_t slot = 0; slot < locals.size(); ++slot) {
                std::vector<double> localGamma(static_cast<std::size_t>(locals[slot].gammaDofs), 0.0);
                for (int row = 0; row < locals[slot].gammaDofs; ++row) {
                    localGamma[static_cast<std::size_t>(row)] = nextGamma[static_cast<std::size_t>(
                        locals[slot].gammaIndices[static_cast<std::size_t>(row)])];
                }
                std::vector<double> image;
                denseMatvec(dynamic[slot].aiGamma, locals[slot].rank,
                    locals[slot].gammaDofs, localGamma, image);
                nextStates[slot] = localRhs[slot];
                for (int mode = 0; mode < locals[slot].rank; ++mode) {
                    nextStates[slot][static_cast<std::size_t>(mode)] -= image[static_cast<std::size_t>(mode)];
                }
                local::solveDenseSymmetric(dynamic[slot].factor, nextStates[slot]);
            }
            gamma = std::move(nextGamma);
            localStates = std::move(nextStates);
            accuracy.recoverySeconds += elapsed(recoveryStart);
            }
            accuracy.localCoreSeconds += elapsed(localCoreStart);
        }

        std::vector<double> rom;
        if (matrixFreeDynamicSchur || options.portReduction) {
            rom = reference;
            for (std::size_t row = 0; row < rom.size(); ++row) {
                rom[row] += matrixFreeRomTheta[row];
            }
        } else {
            rom = reconstruct(reference, partition, locals, localStates, gamma);
        }

        double fullResidual = 0.0;
        double fullResidualBeforeGate = 0.0;
        double fullResidualAfterFallback = 0.0;
        double residualScale = 1.0;
        bool residualGatePassed = true;
        bool residualFallbackUsed = false;
        std::vector<double> acceptedTheta;
        std::vector<double> residual;
        if (step > 0) {
            acceptedTheta = rom;
            for (std::size_t row = 0; row < acceptedTheta.size(); ++row) {
                acceptedTheta[row] -= reference[row];
            }
            const auto residualStart = Clock::now();
            residual = applyFullStep(acceptedTheta);
            for (std::size_t row = 0; row < residual.size(); ++row) {
                residual[row] -= stepRhs[row];
            }
            stepResidualSeconds = elapsed(residualStart);
            accuracy.stepRhsSeconds += stepRhsSeconds;
            accuracy.fullResidualSeconds += stepResidualSeconds;
            residualScale = std::max({
                std::sqrt(normSquared(stepRhs)),
                1.0e-8 * std::sqrt(normSquared(descriptor.boundaryRhs)),
                1.0e-300});
            fullResidualBeforeGate = std::sqrt(normSquared(residual))
                / residualScale;
            if (!std::isfinite(fullResidualBeforeGate)) {
                accuracy.residualGateAllPassed = false;
                throw std::runtime_error(
                    "residual_gate_failed: non-finite full discrete residual.");
            }

            fullResidualAfterFallback = fullResidualBeforeGate;
            if (fullResidualAfterFallback > options.fullResidualTolerance
                && stepAdaptiveToleranceUsed) {
                stepAdaptiveRetry = true;
                ++accuracy.adaptiveInterfaceRetrySteps;
                const auto adaptiveRetryStart = Clock::now();
                local::SolveResult strictSolve = stepUsedNativeReducedHistory
                    ? matrixFreeSolver->solveReducedRhs(
                        nativeInteriorRhs, nativeInterfaceRhs, &gamma,
                        options.interfaceTolerance)
                    : matrixFreeSolver->solve(
                        stepRhs, &gamma, options.interfaceTolerance);
                if (strictSolve.status != "success") {
                    throw std::runtime_error(
                        "[Local transient] Strict adaptive Dynamic Schur "
                        "retry did not converge.");
                }
                matrixFreeRomTheta = std::move(strictSolve.temperature);
                gamma = std::move(strictSolve.interfaceTemperature);
                if (options.nativeReducedHistory) {
                    if (strictSolve.localReducedCoordinates.size()
                        == locals.size()) {
                        localStates = std::move(strictSolve.localReducedCoordinates);
                        nativeReducedHistoryAligned = true;
                    } else {
                        nativeReducedHistoryAligned = false;
                    }
                }
                stepAdaptiveRetrySeconds = elapsed(adaptiveRetryStart);
                accuracy.adaptiveInterfaceRetrySeconds +=
                    stepAdaptiveRetrySeconds;
                accuracy.adaptiveInterfaceRetryIterations +=
                    strictSolve.timing.interfaceIterations;
                accuracy.interfaceSolveSeconds +=
                    strictSolve.timing.interfaceSolveSeconds;
                accuracy.recoverySeconds +=
                    strictSolve.timing.localRecoverySeconds
                    + strictSolve.timing.fullFieldReconstructionSeconds;
                accuracy.proxySolveSeconds +=
                    strictSolve.timing.proxySolveSeconds;
                accuracy.coarseSolveSeconds +=
                    strictSolve.timing.coarseSolveSeconds;
                accuracy.portForwardSolveSeconds +=
                    strictSolve.timing.portForwardSolveSeconds;
                accuracy.portCoreSolveSeconds +=
                    strictSolve.timing.portCoreSolveSeconds;
                accuracy.portBackSubstitutionSeconds +=
                    strictSolve.timing.portBackSubstitutionSeconds;
                accuracy.interfaceOperatorSeconds +=
                    strictSolve.timing.interfaceOperatorSeconds;
                accuracy.interfacePreconditionerSeconds +=
                    strictSolve.timing.interfacePreconditionerSeconds;
                accuracy.interfaceOrthogonalizationSeconds +=
                    strictSolve.timing.interfaceOrthogonalizationSeconds;
                accuracy.interfaceVectorUpdateSeconds +=
                    strictSolve.timing.interfaceVectorUpdateSeconds;
                accuracy.interfaceIterationsTotal +=
                    strictSolve.timing.interfaceIterations;
                accuracy.interfaceMatvecs += strictSolve.timing.interfaceMatvecs;
                accuracy.maximumInterfaceResidual = std::max(
                    accuracy.maximumInterfaceResidual,
                    strictSolve.timing.interfaceRelativeResidual);
                stepInterfaceIterations +=
                    strictSolve.timing.interfaceIterations;
                stepInterfaceMatvecs += strictSolve.timing.interfaceMatvecs;
                stepInterfaceResidual =
                    strictSolve.timing.interfaceRelativeResidual;
                stepLocalSolveSeconds +=
                    strictSolve.timing.localReducedAssemblySeconds;
                stepInterfaceSolveSeconds +=
                    strictSolve.timing.interfaceSolveSeconds;
                stepProxySolveSeconds += strictSolve.timing.proxySolveSeconds;
                stepCoarseSolveSeconds += strictSolve.timing.coarseSolveSeconds;
                stepPortForwardSolveSeconds +=
                    strictSolve.timing.portForwardSolveSeconds;
                stepPortCoreSolveSeconds +=
                    strictSolve.timing.portCoreSolveSeconds;
                stepPortBackSubstitutionSeconds +=
                    strictSolve.timing.portBackSubstitutionSeconds;
                stepSchurApplySeconds +=
                    strictSolve.timing.interfaceOperatorSeconds;
                stepPreconditionerSeconds +=
                    strictSolve.timing.interfacePreconditionerSeconds;
                stepOrthogonalizationSeconds +=
                    strictSolve.timing.interfaceOrthogonalizationSeconds;
                stepVectorUpdateSeconds +=
                    strictSolve.timing.interfaceVectorUpdateSeconds;
                accuracy.interfaceIterationsMaximum = std::max(
                    accuracy.interfaceIterationsMaximum,
                    stepInterfaceIterations);
                accuracy.localCoreSeconds += stepAdaptiveRetrySeconds;

                rom = reference;
                for (std::size_t row = 0; row < rom.size(); ++row) {
                    rom[row] += matrixFreeRomTheta[row];
                }
                acceptedTheta = matrixFreeRomTheta;
                const auto retryResidualStart = Clock::now();
                residual = applyFullStep(acceptedTheta);
                for (std::size_t row = 0; row < residual.size(); ++row) {
                    residual[row] -= stepRhs[row];
                }
                const double retryResidualSeconds =
                    elapsed(retryResidualStart);
                stepResidualSeconds += retryResidualSeconds;
                accuracy.fullResidualSeconds += retryResidualSeconds;
                fullResidualAfterFallback = std::sqrt(normSquared(residual))
                    / residualScale;
                if (!std::isfinite(fullResidualAfterFallback)) {
                    fullResidualAfterFallback =
                        std::numeric_limits<double>::infinity();
                }
            }
            if (fullResidualAfterFallback > options.fullResidualTolerance) {
                ++accuracy.residualToleranceViolationSteps;
                if (!options.fullResidualFallback) {
                    residualGatePassed = false;
                    accuracy.residualGateAllPassed = false;
                } else {
                    residualFallbackUsed = true;
                    ++accuracy.residualFallbackSteps;
                    const auto fallbackStart = Clock::now();
                    // The factor is for the exact current-step operator C/dt + K.
                    ensureFomFactor().solve(stepRhs, acceptedTheta);
                    accuracy.residualFallbackSolveSeconds += elapsed(fallbackStart);
                    for (std::size_t row = 0; row < rom.size(); ++row) {
                        rom[row] = reference[row] + acceptedTheta[row];
                    }
                    residual = applyFullStep(acceptedTheta);
                    for (std::size_t row = 0; row < residual.size(); ++row) {
                        residual[row] -= stepRhs[row];
                    }
                    fullResidualAfterFallback = std::sqrt(normSquared(residual))
                        / residualScale;
                    if (!std::isfinite(fullResidualAfterFallback)
                        || fullResidualAfterFallback
                            > options.fullResidualTolerance) {
                        residualGatePassed = false;
                        accuracy.residualGateAllPassed = false;
                        throw std::runtime_error(
                            "residual_gate_failed: full-order fallback did not satisfy "
                            "--mor-full-residual-tolerance.");
                    }

                    // Synchronize every history representation with the
                    // accepted fallback state before advancing the next step.
                    if (matrixFreeDynamicSchur || options.portReduction) {
                        matrixFreeRomTheta = acceptedTheta;
                        for (std::size_t gammaIndex = 0;
                             gammaIndex < partition.interfaceGlobalDofs.size();
                             ++gammaIndex) {
                            gamma[gammaIndex] = acceptedTheta[
                                static_cast<std::size_t>(
                                    partition.interfaceGlobalDofs[gammaIndex])];
                        }
                        if (matrixFreeDynamicSchur
                            && options.nativeReducedHistory) {
                            nativeReducedHistoryAligned = false;
                        }
                    } else {
                        for (std::size_t gammaIndex = 0;
                             gammaIndex < partition.interfaceGlobalDofs.size();
                             ++gammaIndex) {
                            gamma[gammaIndex] = acceptedTheta[
                                static_cast<std::size_t>(
                                    partition.interfaceGlobalDofs[gammaIndex])];
                        }
                        for (std::size_t slot = 0; slot < locals.size(); ++slot) {
                            localStates[slot] = projectInitial(
                                locals[slot], c.interior[slot], acceptedTheta);
                        }
                    }
                }
            }
            fullResidual = fullResidualAfterFallback;
            accuracy.maximumFullResidualBeforeGate = std::max(
                accuracy.maximumFullResidualBeforeGate,
                fullResidualBeforeGate);
            accuracy.maximumFullResidual = std::max(
                accuracy.maximumFullResidual, fullResidual);
            long double stateChangeSquared = 0.0L;
            long double previousStateSquared = 0.0L;
            stepStateExactlyUnchanged =
                acceptedTheta.size() == previousRomTheta.size();
            for (std::size_t row = 0; row < acceptedTheta.size(); ++row) {
                const double difference = acceptedTheta[row]
                    - previousRomTheta[row];
                stateChangeSquared +=
                    static_cast<long double>(difference) * difference;
                previousStateSquared +=
                    static_cast<long double>(previousRomTheta[row])
                    * previousRomTheta[row];
                stepStateExactlyUnchanged = stepStateExactlyUnchanged
                    && acceptedTheta[row] == previousRomTheta[row];
            }
            stepStateChangeRelative =
                std::sqrt(static_cast<double>(stateChangeSquared))
                / std::max(1.0e-300,
                    std::sqrt(static_cast<double>(previousStateSquared)));
            lastStateChangeRelative = stepStateChangeRelative;
            previousRomTheta = acceptedTheta;
        }
        const double unavailable = std::numeric_limits<double>::quiet_NaN();
        std::vector<double> fom = rom;
        long double errorSquared = 0.0L;
        long double referenceSquared = 0.0L;
        double maximum = unavailable;
        double fomMax = unavailable;
        const double romMax = *std::max_element(rom.begin(), rom.end());
        double stepRelativeL2 = unavailable;
        double stepMaximumTemperatureError = unavailable;
        if (options.compareFom) {
            fom = reference;
            for (std::size_t row = 0; row < fom.size(); ++row) {
                fom[row] += fomTheta[row];
            }
            maximum = 0.0;
            for (std::size_t row = 0; row < rom.size(); ++row) {
                const double error = rom[row] - fom[row];
                errorSquared += static_cast<long double>(error) * error;
                referenceSquared += static_cast<long double>(fom[row]) * fom[row];
                maximum = std::max(maximum, std::abs(error));
            }
            stepRelativeL2 = std::sqrt(static_cast<double>(errorSquared))
                / std::max(1.0e-300,
                    std::sqrt(static_cast<double>(referenceSquared)));
            fomMax = *std::max_element(fom.begin(), fom.end());
            stepMaximumTemperatureError = std::abs(fomMax - romMax);
            totalErrorSquared += errorSquared;
            totalReferenceSquared += referenceSquared;
            accuracy.maximumAbsolute = std::max(accuracy.maximumAbsolute, maximum);
            accuracy.maximumTemperatureError = std::max(
                accuracy.maximumTemperatureError,
                stepMaximumTemperatureError);
        }

        double reducedResidual = 0.0;
        std::vector<double> corrected = rom;
        if (step > 0) {
            if (options.compareFom && options.milestone8AdaptiveProduction) {
                for (const LocalPortBasis& port :
                     portModel->ports) {
                    long double localErrorSquared = 0.0L;
                    long double localReferenceSquared = 0.0L;
                    long double localResidualSquared = 0.0L;
                    double localMaximumError = 0.0;
                    for (int gammaIndex :
                         port.interfaceIndices) {
                        const int global =
                            partition.interfaceGlobalDofs[
                                static_cast<std::size_t>(
                                    gammaIndex)];
                        const double error =
                            rom[static_cast<std::size_t>(global)]
                            - fom[static_cast<std::size_t>(global)];
                        localErrorSquared +=
                            static_cast<long double>(error) * error;
                        localReferenceSquared +=
                            static_cast<long double>(
                                fom[static_cast<std::size_t>(global)])
                            * fom[static_cast<std::size_t>(global)];
                        localResidualSquared +=
                            static_cast<long double>(
                                residual[
                                    static_cast<std::size_t>(global)])
                            * residual[
                                static_cast<std::size_t>(global)];
                        localMaximumError = std::max(
                            localMaximumError,
                            std::abs(error));
                    }
                    AdaptivePhysicalInterfaceRow& row =
                        adaptivePhysicalRows[port.interfaceId];
                    row.temperatureRelativeL2 =
                        std::sqrt(static_cast<double>(
                            localErrorSquared))
                        / std::max(
                            1.0e-300,
                            std::sqrt(static_cast<double>(
                                localReferenceSquared)));
                    row.maximumNodalError =
                        localMaximumError;
                    row.localFullResidual =
                        std::sqrt(static_cast<double>(
                            localResidualSquared))
                        / residualScale;
                }
            }
            std::vector<double> reducedPieces;
            for (const LocalModel& localModel : locals) {
                for (int mode = 0; mode < localModel.rank; ++mode) {
                    long double value = 0.0L;
                    for (int row = 0; row < localModel.interiorDofs; ++row) {
                        value += static_cast<long double>(localModel.basis[static_cast<std::size_t>(
                            mode * localModel.interiorDofs + row)])
                            * residual[static_cast<std::size_t>(localModel.interiorGlobal[
                                static_cast<std::size_t>(row)])];
                    }
                    reducedPieces.push_back(static_cast<double>(value));
                }
            }
            for (int global : partition.interfaceGlobalDofs) {
                reducedPieces.push_back(residual[static_cast<std::size_t>(global)]);
            }
            reducedResidual = std::sqrt(normSquared(reducedPieces))
                / residualScale;
            accuracy.maximumReducedResidual = std::max(
                accuracy.maximumReducedResidual, reducedResidual);
            if (options.portReduction && options.localPortCorrected) {
                std::vector<double> correctionRhs(residual.size(), 0.0);
                for (std::size_t row = 0; row < residual.size(); ++row) {
                    correctionRhs[row] = -residual[row];
                }
                const auto correctedStart = Clock::now();
                std::vector<double> correction;
                ensureFomFactor().solve(correctionRhs, correction);
                accuracy.correctedSolveSeconds += elapsed(correctedStart);
                std::vector<double> correctedTheta = acceptedTheta;
                for (std::size_t row = 0; row < correctedTheta.size(); ++row) {
                    correctedTheta[row] += correction[row];
                    corrected[row] = reference[row] + correctedTheta[row];
                }
                std::vector<double> correctedResidual =
                    applyFullStep(correctedTheta);
                for (std::size_t row = 0; row < correctedResidual.size(); ++row) {
                    correctedResidual[row] -= stepRhs[row];
                }
                accuracy.correctedMaximumResidual = std::max(
                    accuracy.correctedMaximumResidual,
                    std::sqrt(normSquared(correctedResidual)) / residualScale);
            }
        }
        if (options.compareFom
            && options.portReduction && options.localPortCorrected) {
            for (std::size_t row = 0; row < corrected.size(); ++row) {
                const double error = corrected[row] - fom[row];
                correctedErrorSquared += static_cast<long double>(error) * error;
                correctedReferenceSquared +=
                    static_cast<long double>(fom[row]) * fom[row];
                accuracy.correctedMaximumAbsolute = std::max(
                    accuracy.correctedMaximumAbsolute, std::abs(error));
            }
        }
        const bool writeDetailedFlux = options.compareFom
            && !options.compareFomSummaryOnly
            && (partition.interfaceGlobalDofs.size() <= 200000U || step == steps);
        InterfacePhysicsMetrics aggregateInterfaceMetrics;
        double fluxRelativeL2 = std::numeric_limits<double>::quiet_NaN();
        int worstInterfaceId = -1;
        int worstFacePairId = -1;
        int worstIntegrationTriangle = -1;
        if (writeDetailedFlux) {
            const DetailedInterfacePhysicsMetrics interfaceMetrics =
                calculateDetailedInterfacePhysicsMetrics(mesh, physics, rom);
            const DetailedInterfacePhysicsMetrics fomInterfaceMetrics =
                calculateDetailedInterfacePhysicsMetrics(mesh, physics, fom);
            if (interfaceMetrics.triangles.size()
                    != fomInterfaceMetrics.triangles.size()) {
                throw std::runtime_error(
                    "[Local transient] FOM/ROM interface diagnostic ordering mismatch.");
            }
            aggregateInterfaceMetrics = interfaceMetrics.aggregate;
            worstInterfaceId = interfaceMetrics.worstInterfaceId;
            worstFacePairId = interfaceMetrics.worstFacePairId;
            worstIntegrationTriangle = interfaceMetrics.worstIntegrationTriangle;
        long double fluxErrorSquared = 0.0L;
        long double fluxReferenceSquared = 0.0L;
        for (std::size_t triangle = 0;
             triangle < interfaceMetrics.triangles.size(); ++triangle) {
            const InterfaceTriangleFluxRecord& reducedFlux =
                interfaceMetrics.triangles[triangle];
            const InterfaceTriangleFluxRecord& fomFlux =
                fomInterfaceMetrics.triangles[triangle];
            const double fluxError = reducedFlux.sipgNumericalFlux
                - fomFlux.sipgNumericalFlux;
            fluxErrorSquared += static_cast<long double>(reducedFlux.area)
                * fluxError * fluxError;
            fluxReferenceSquared += static_cast<long double>(fomFlux.area)
                * fomFlux.sipgNumericalFlux * fomFlux.sipgNumericalFlux;
            const double relativeError = std::abs(fluxError) / std::max({
                std::abs(fomFlux.sipgNumericalFlux),
                fomInterfaceMetrics.relativeFluxFloor, 1.0e-300});
            if (options.milestone8AdaptiveProduction
                && step > 0) {
                const auto adjacent = std::minmax(
                    reducedFlux.leftSubdomain,
                    reducedFlux.rightSubdomain);
                const auto physicalInterface =
                    adaptivePhysicalInterfaceBySubdomains.find(
                        std::make_pair(
                            adjacent.first, adjacent.second));
                if (physicalInterface
                    == adaptivePhysicalInterfaceBySubdomains.end()) {
                    throw std::runtime_error(
                        "[M8.9 adaptive] Flux triangle adjacent "
                        "subdomains do not map to a physical port.");
                }
                AdaptivePhysicalInterfaceRow& row =
                    adaptivePhysicalRows[
                        physicalInterface->second];
                row.fluxErrorSquared +=
                    static_cast<long double>(reducedFlux.area)
                    * fluxError * fluxError;
                row.fluxReferenceSquared +=
                    static_cast<long double>(fomFlux.area)
                    * fomFlux.sipgNumericalFlux
                    * fomFlux.sipgNumericalFlux;
                row.jumpAreaWeightedSquared +=
                    static_cast<long double>(reducedFlux.area)
                    * reducedFlux.temperatureJumpRms
                    * reducedFlux.temperatureJumpRms;
                row.area += reducedFlux.area;
                row.maximumRelativeFluxImbalance = std::max(
                    row.maximumRelativeFluxImbalance,
                    reducedFlux.relativeFluxImbalance);
            }
            fluxOut << waveform.name << ',' << options.initialMode << ',' << step << ','
                << time << ',' << reducedFlux.interfaceId << ','
                << reducedFlux.facePairId << ',' << reducedFlux.integrationTriangleId << ','
                << reducedFlux.leftSubdomain << ',' << reducedFlux.rightSubdomain << ','
                << reducedFlux.leftBoundaryEntity << ','
                << reducedFlux.rightBoundaryEntity << ',' << reducedFlux.area << ','
                << reducedFlux.temperatureJumpRms << ','
                << reducedFlux.leftPhysicalNormalFlux << ','
                << reducedFlux.rightPhysicalNormalFlux << ','
                << reducedFlux.sipgNumericalFlux << ','
                << reducedFlux.fluxImbalanceL2 << ','
                << reducedFlux.relativeFluxImbalance << ','
                << fomFlux.temperatureJumpRms << ','
                << fomFlux.leftPhysicalNormalFlux << ','
                << fomFlux.rightPhysicalNormalFlux << ','
                << fomFlux.sipgNumericalFlux << ',' << fomFlux.fluxImbalanceL2 << ','
                << fomFlux.relativeFluxImbalance << ',' << fluxError << ','
                << relativeError << '\n';
        }
            fluxRelativeL2 = std::sqrt(static_cast<double>(fluxErrorSquared))
            / std::max(1.0e-300,
                std::sqrt(static_cast<double>(fluxReferenceSquared)));
            accuracy.maximumFluxRelativeL2 = std::max(
                accuracy.maximumFluxRelativeL2, fluxRelativeL2);
        } else {
            aggregateInterfaceMetrics = calculateInterfacePhysicsMetrics(
                mesh, physics, rom);
        }
        accuracy.maximumJump = std::max(
            accuracy.maximumJump, aggregateInterfaceMetrics.temperatureJumpRms);
        accuracy.maximumFluxImbalance = std::max(
            accuracy.maximumFluxImbalance,
            aggregateInterfaceMetrics.relativeFluxImbalance);
        if (options.compareFom && options.milestone8AdaptiveProduction
            && step > 0) {
            for (const LocalPortBasis& port :
                 portModel->ports) {
                const AdaptiveProductionRank policy =
                    milestone8AdaptiveProductionRank(
                        port.interfaceId);
                const AdaptivePhysicalInterfaceRow& row =
                    adaptivePhysicalRows[port.interfaceId];
                adaptivePhysicalDiagnostics
                    << step << ',' << time << ','
                    << port.interfaceId << ','
                    << policy.historyRank << ','
                    << policy.randomizedRank << ','
                    << policy.residualRank << ','
                    << port.rank << ',' << port.rows << ','
                    << row.temperatureRelativeL2 << ','
                    << row.maximumNodalError << ','
                    << row.localFullResidual << ','
                    << (row.fluxReferenceSquared > 0.0L
                        ? std::sqrt(static_cast<double>(
                            row.fluxErrorSquared
                            / row.fluxReferenceSquared))
                        : std::numeric_limits<double>::quiet_NaN())
                    << ','
                    << (row.area > 0.0L
                        ? std::sqrt(static_cast<double>(
                            row.jumpAreaWeightedSquared
                            / row.area))
                        : 0.0)
                    << ','
                    << row.maximumRelativeFluxImbalance
                    << ",0\n";
            }
            adaptivePhysicalDiagnostics.flush();
        }
        byTime << waveform.name << ',' << options.initialMode << ',' << step << ',' << time << ','
            << stepRelativeL2 << ',' << maximum << ',' << fomMax << ','
            << romMax << ',' << stepMaximumTemperatureError
            << ',' << fullResidual << ',' << fullResidualBeforeGate << ','
            << options.fullResidualTolerance << ','
            << (options.fullResidualFallback ? 1 : 0) << ','
            << (residualGatePassed ? 1 : 0) << ','
            << (residualFallbackUsed ? 1 : 0) << ','
            << fullResidualAfterFallback << ',' << reducedResidual << ','
            << aggregateInterfaceMetrics.temperatureJumpRms << ','
            << aggregateInterfaceMetrics.relativeFluxImbalance << ','
            << fluxRelativeL2 << ',' << worstInterfaceId << ','
            << worstFacePairId << ',' << worstIntegrationTriangle << ','
            << stepInterfaceIterations << ',' << stepInterfaceMatvecs << ','
            << stepInterfaceResidual << ',' << stepPortProjectionError << '\n';
        if ((reducedMethod || matrixFreeDynamicSchur || options.portReduction)
            && step > 0) {
            reducedTiming << step << ',' << time << ','
                << (options.portReduction ? "port_reduced_zero"
                    : options.interfaceInitialGuess)
                << ',' << stepInterfaceKrylovActual << ','
                << (stepInterfaceKrylovFallback ? 1 : 0) << ','
                << (stepInterfaceKrylovFallbackReason.empty()
                    ? "none" : stepInterfaceKrylovFallbackReason) << ','
                << (stepAdaptiveToleranceUsed
                    ? options.adaptiveInterfaceTolerance
                    : options.interfaceTolerance) << ','
                << (stepAdaptiveRetry ? 1 : 0) << ','
                << stepAdaptiveRetrySeconds << ','
                << stepInterfaceInitialResidual << ','
                << (stepInterfacePredictorApplied ? 1 : 0) << ','
                << (stepInterfacePredictorAccepted ? 1 : 0) << ','
                << stepInterfacePredictorInitialResidual << ','
                << stepInterfacePredictorResidual << ','
                << stepInterfacePredictorSeconds << ','
                << stepLocalSolveSeconds << ',' << stepNativeReducedRhsSeconds
                << ',' << stepRhsSeconds << ','
                << stepInterfaceSolveSeconds << ','
                << stepProxySolveSeconds << ',' << stepCoarseSolveSeconds << ','
                << stepPortForwardSolveSeconds << ','
                << stepPortCoreSolveSeconds << ','
                << stepPortBackSubstitutionSeconds << ','
                << stepSchurApplySeconds << ','
                << stepPreconditionerSeconds << ','
                << stepOrthogonalizationSeconds << ','
                << stepVectorUpdateSeconds << ','
                << stepResidualSeconds << ',' << stepStateChangeRelative << ','
                << elapsed(stepStart) << ','
                << stepInterfaceIterations << ',' << fullResidual << '\n';
        }
        if (step > 0 && step < steps && stepStateExactlyUnchanged
            && residualGatePassed && !options.compareFom
            && !options.portReduction) {
            bool futureInputUnchanged = true;
            for (int futureStep = step + 1;
                 futureStep <= steps && futureInputUnchanged;
                 ++futureStep) {
                futureInputUnchanged = waveform.sample(
                    futureStep * options.timeStep) == stepPowers;
            }
            if (futureInputUnchanged) {
                steadyReuse.active = true;
                steadyReuse.rom = rom;
                steadyReuse.romMaximum = romMax;
                steadyReuse.fullResidual = fullResidual;
                steadyReuse.fullResidualBeforeGate = fullResidualBeforeGate;
                steadyReuse.fullResidualAfterFallback =
                    fullResidualAfterFallback;
                steadyReuse.reducedResidual = reducedResidual;
                steadyReuse.interfaceMetrics = aggregateInterfaceMetrics;
                steadyStateDetectedStep = step;
            }
        }
        if (step == steps) finalRom = rom;
    }
    const double timeSteppingSeconds = elapsed(timeSteppingStart);
    const auto postTimeSteppingStart = Clock::now();
    if (options.compareFom) {
        accuracy.spaceTimeRelativeL2 =
            std::sqrt(static_cast<double>(totalErrorSquared))
            / std::max(1.0e-300,
                std::sqrt(static_cast<double>(totalReferenceSquared)));
    } else {
        const double unavailable = std::numeric_limits<double>::quiet_NaN();
        accuracy.spaceTimeRelativeL2 = unavailable;
        accuracy.maximumAbsolute = unavailable;
        accuracy.maximumTemperatureError = unavailable;
        accuracy.maximumFluxRelativeL2 = unavailable;
    }
    if (options.portReduction && options.localPortCorrected) {
        accuracy.correctedRelativeL2 =
            std::sqrt(static_cast<double>(correctedErrorSquared))
            / std::max(1.0e-300,
                std::sqrt(static_cast<double>(correctedReferenceSquared)));
        std::ofstream correctedOut(
            outputDirectory / "local_port_corrected_accuracy.csv");
        correctedOut << "relative_l2,maximum_absolute_k,maximum_full_residual,"
            "correction_solve_seconds,status\n" << std::setprecision(17)
            << accuracy.correctedRelativeL2 << ','
            << accuracy.correctedMaximumAbsolute << ','
            << accuracy.correctedMaximumResidual << ','
            << accuracy.correctedSolveSeconds << ','
            << ((accuracy.correctedRelativeL2 < 1.0e-4
                    && accuracy.correctedMaximumAbsolute < 0.1
                    && accuracy.correctedMaximumResidual < 1.0e-8)
                    ? "success" : "accuracy_failed") << '\n';
    }

    if (options.deploymentRhsCount > 1) {
        if (!matrixFreeDynamicSchur && !options.portReduction) {
            throw std::runtime_error(
                "[Local transient] Multi-waveform deployment benchmark requires matrix-free Dynamic Schur.");
        }
        const auto deploymentStart = Clock::now();
        double onlineSeconds = 0.0;
        long long deploymentIterations = 0;
        int maximumDeploymentIterations = 0;
        double maximumDeploymentResidual = 0.0;
        double outputChecksum = 0.0;
        for (int deployment = 0;
             deployment < options.deploymentRhsCount; ++deployment) {
            const PowerWaveform deploymentWaveform = makeBuiltinWaveform(
                "unseen_waveform", descriptor.nominalPowersW,
                options.timeStep, 1, options.seed
                    + static_cast<std::uint64_t>(deployment + 1));
            const std::vector<double> powers = deploymentWaveform.sample(
                options.timeStep);
            const auto onlineStart = Clock::now();
            std::vector<double> globalRhs =
                descriptor.capacity.multiply(thetaInitial);
            for (double& value : globalRhs) value /= options.timeStep;
            addInput(descriptor, powers, globalRhs);
            for (std::size_t row = 0; row < globalRhs.size(); ++row) {
                globalRhs[row] += boundaryOffset[row];
            }
            local::SolveResult solve;
            if (options.portReduction) {
                solve = portSolver->solve(globalRhs).solution;
            } else {
                solve = matrixFreeSolver->solve(globalRhs);
            }
            if (solve.status != "success") {
                throw std::runtime_error(
                    "[Local transient] Multi-waveform deployment FGMRES did not converge.");
            }
            onlineSeconds += elapsed(onlineStart);
            deploymentIterations += solve.timing.interfaceIterations;
            maximumDeploymentIterations = std::max(
                maximumDeploymentIterations, solve.timing.interfaceIterations);
            maximumDeploymentResidual = std::max(
                maximumDeploymentResidual,
                solve.timing.interfaceRelativeResidual);
            outputChecksum += *std::max_element(
                solve.temperature.begin(), solve.temperature.end());
        }
        const double deploymentWallSeconds = elapsed(deploymentStart);
        std::ofstream deploymentOut(
            outputDirectory / "local_dynamic_schur_deployment_timing.csv");
        deploymentOut << "waveforms,steps_per_waveform,setup_reused,proxy_reused,"
            "port_factor_reused,local_factors_reused,total_online_seconds,wall_seconds,"
            "average_online_seconds_per_waveform,average_interface_iterations,"
            "maximum_interface_iterations,maximum_interface_relative_residual,"
            "maximum_temperature_checksum_k,peak_working_set_bytes\n"
            << std::setprecision(17) << options.deploymentRhsCount << ",1,1,"
            << (options.portReduction ? 0 : 1) << ','
            << (options.portReduction ? 1 : 0) << ",1,"
            << onlineSeconds << ',' << deploymentWallSeconds << ','
            << onlineSeconds / options.deploymentRhsCount << ','
            << static_cast<double>(deploymentIterations)
                / options.deploymentRhsCount << ','
            << maximumDeploymentIterations << ',' << maximumDeploymentResidual << ','
            << outputChecksum << ',' << peakWorkingSetBytes() << '\n';
    }

    if ((options.compareFom && !options.compareFomSummaryOnly)
        || options.outputMode == "full-field") {
        std::ofstream finalTemperature(
            outputDirectory / "local_dynamic_schur_final_temperature.csv");
        finalTemperature << "global_dof,x_m,y_m,z_m,subdomain,temperature_k\n"
            << std::setprecision(17);
        for (int row = 0; row < descriptor.dofs; ++row) {
            const DeploymentDof& dof = descriptor.deploymentDofs[
                static_cast<std::size_t>(row)];
            finalTemperature << row << ',' << dof.x << ',' << dof.y << ','
                << dof.z << ',' << dof.subdomain << ','
                << finalRom[static_cast<std::size_t>(row)] << '\n';
        }
    }

    std::ofstream rankOut(outputDirectory / "local_block_arnoldi_rank.csv");
    rankOut << "subdomain,template_id,template_reused,template_fingerprint,"
        "interior_dofs,interface_dofs,physical_power_channels,"
        "interface_excitation_rank,initial_block_rank,moment,input_columns,added_rank,cumulative_rank,"
        "deflated_columns,orthogonality_error,arnoldi_residual,symbolic_calls,"
        "numerical_calls,symbolic_seconds,numerical_seconds,multi_rhs_solve_seconds,"
        "orthogonalization_seconds\n" << std::setprecision(17);
    for (const LocalModel& localModel : locals) {
        for (const ArnoldiHistoryRow& row : localModel.history) {
            rankOut << localModel.domainId << ',' << localModel.templateId << ','
                << (localModel.templateReused ? 1 : 0) << ','
                << localModel.templateFingerprint << ',' << localModel.interiorDofs << ','
                << localModel.gammaDofs << ',' << localModel.physicalChannels << ','
                << localModel.excitationRank << ',' << localModel.initialBlockRank << ','
                << row.moment << ',' << row.inputColumns << ',' << row.addedRank << ','
                << row.cumulativeRank << ','
                << row.deflatedColumns << ',' << row.orthogonalityError << ','
                << row.arnoldiResidual << ','
                << localModel.arnoldiTiming.symbolicAnalysisCalls << ','
                << localModel.arnoldiTiming.numericalFactorizationCalls << ','
                << localModel.arnoldiTiming.symbolicAnalysisSeconds << ','
                << localModel.arnoldiTiming.numericalFactorizationSeconds << ','
                << localModel.arnoldiTiming.multiRhsSolveSeconds << ','
                << localModel.arnoldiTiming.orthogonalizationSeconds << '\n';
        }
    }

    if (options.portReduction) {
        std::ofstream enrichmentOut(
            outputDirectory / "local_port_enrichment_history.csv");
        enrichmentOut << "subdomain,round,input_columns,added_rank,cumulative_rank,"
            "deflated_columns,orthogonality_error,factorization_seconds,"
            "solve_seconds,orthogonalization_seconds,total_seconds\n"
            << std::setprecision(17);
        for (const LocalModel& localModel : locals) {
            for (const ArnoldiHistoryRow& row : localModel.history) {
                if (row.moment <= options.moments) continue;
                enrichmentOut << localModel.domainId << ','
                    << (row.moment - options.moments) << ','
                    << row.inputColumns << ',' << row.addedRank << ','
                    << row.cumulativeRank << ',' << row.deflatedColumns << ','
                    << row.orthogonalityError << ','
                    << enrichment.factorizationSeconds << ','
                    << enrichment.solveSeconds << ','
                    << enrichment.orthogonalizationSeconds << ','
                    << enrichment.totalSeconds << '\n';
            }
        }
    }

    std::size_t modelBytes = partition.interfaceGlobalDofs.size() * sizeof(int)
        + schurMatrixEntries.size() * sizeof(MatrixEntry);
    for (const LocalModel& localModel : locals) {
        modelBytes += localModel.interiorGlobal.size() * sizeof(int)
            + localModel.gammaIndices.size() * sizeof(int)
            + localModel.basis.size() * sizeof(double)
            + localModel.cii.size() * sizeof(double)
            + localModel.kii.size() * sizeof(double)
            + localModel.ciGamma.size() * sizeof(double)
            + localModel.kiGamma.size() * sizeof(double)
            + localModel.cGammaI.size() * sizeof(double)
            + localModel.kGammaI.size() * sizeof(double)
            + localModel.reducedInput.size() * sizeof(double)
            + localModel.reducedBoundary.size() * sizeof(double)
            + localModel.referenceInterior.size() * sizeof(double);
    }
    if (matrixFreeModel) {
        modelBytes += matrixFreeModel->interfaceEntries.size()
            * sizeof(local::InterfaceEntry);
    }
    if (portModel) modelBytes += portModel->modelBytes;
    if (globalCoarseModel) {
        modelBytes +=
            globalCoarseModel->basis.capacity() * sizeof(double)
            + globalCoarseModel->schurImages.capacity()
                * sizeof(double)
            + globalCoarseModel->eigenvalues.capacity()
                * sizeof(double)
            + globalCoarseModel->ritzResiduals.capacity()
                * sizeof(double);
    }
    const std::string interfaceSolver = options.portReduction
        ? "port-reduced-dense"
        : (matrixFreeDynamicSchur
            ? matrixFreeSolver->interfaceSolver() : "sparse-pardiso");
    const int schurSymbolicCalls = options.portReduction ? 0 : (matrixFreeDynamicSchur
        ? matrixFreeSolver->symbolicAnalysisCalls()
        : schurFactor->symbolicAnalysisCalls());
    const int schurNumericalCalls = options.portReduction ? 1 : (matrixFreeDynamicSchur
        ? matrixFreeSolver->numericalFactorizationCalls()
        : schurFactor->numericalFactorizationCalls());
    const char* finalStatus = options.fullResidualFallback
        && !accuracy.residualGateAllPassed
        ? "residual_gate_failed"
        : (!options.compareFom
            ? "success"
            : ((accuracy.spaceTimeRelativeL2 < 1.0e-4
                && accuracy.maximumAbsolute < 0.1
                && accuracy.maximumTemperatureError < 0.01)
            ? "success" : "accuracy_failed"));
    // Enrichment changes each basis and every projected local block. Publish
    // the cache only after reduced factorization, time stepping, and all
    // requested residual/accuracy gates accept the completed workflow.
    if (std::string(finalStatus) == "success" && !dynamicSavePath.empty()
        && (!localDynamicCacheHit || dynamicSavePath != dynamicLoadPath
            || enrichment.addedRank > 0)) {
        const auto saveStart = Clock::now();
        saveLocalDynamicModels(
            dynamicSavePath, descriptor, options, locals);
        localDynamicCacheSaveSeconds = elapsed(saveStart);
        std::error_code sizeError;
        localDynamicCacheBytes =
            std::filesystem::file_size(dynamicSavePath, sizeError);
        if (sizeError) localDynamicCacheBytes = 0;
        std::cout << "[Local dynamic cache] saved "
                  << dynamicSavePath << ", bytes="
                  << localDynamicCacheBytes << ", seconds="
                  << localDynamicCacheSaveSeconds << '\n';
    }
    const std::size_t peakBytes = peakWorkingSetBytes();
    const double postTimeSteppingSeconds = elapsed(postTimeSteppingStart);
    const double totalSeconds = elapsed(totalStart);
    std::ofstream summary(outputDirectory / "local_dynamic_schur_summary.csv");
    summary << "status,method,subdomains,global_dofs,full_interface_dofs,total_local_rank,"
        "unique_templates,reused_instances,reuse_requested,waveform,initial_mode,dt_s,steps,"
        "steady_state_detected_step,steady_state_reused_steps,"
        "last_state_change_relative,space_time_relative_l2,maximum_absolute_k,"
        "maximum_temperature_error_k,maximum_full_residual,"
        "maximum_full_residual_before_gate,full_residual_tolerance,"
        "residual_fallback_enabled,residual_gate_passed,"
        "residual_tolerance_violation_steps,residual_fallback_steps,"
        "residual_fallback_solve_seconds,maximum_reduced_residual,"
        "maximum_temperature_jump_rms_k,maximum_relative_flux_imbalance,"
        "maximum_fom_rom_flux_relative_l2,reference_setup_seconds,"
        "reference_cache_hit,reference_cache_load_seconds,"
        "reference_cache_save_seconds,reference_cache_bytes,"
        "local_basis_setup_seconds,arnoldi_second_moment_energy,"
        "arnoldi_second_moment_max_columns,local_model_cache_hit,"
        "local_model_cache_load_seconds,local_model_cache_save_seconds,"
        "local_model_cache_bytes,local_symbolic_seconds,local_numerical_seconds,"
        "local_multi_rhs_seconds,local_orthogonalization_seconds,"
        "local_trace_basis_seconds,local_input_setup_seconds,"
        "local_fingerprint_seconds,local_orthogonality_audit_seconds,"
        "local_projection_seconds,dynamic_schur_setup_seconds,"
        "dynamic_schur_factor_seconds,dynamic_schur_symbolic_calls,dynamic_schur_numerical_calls,"
        "port_core_cache_enabled,port_core_cache_hit,"
        "port_core_cache_load_seconds,port_core_cache_save_seconds,"
        "port_core_cache_bytes,port_core_partition_seconds,"
        "port_core_coupling_assembly_seconds,port_core_leaf_csr_seconds,"
        "port_core_leaf_factor_seconds,port_core_elimination_seconds,"
        "port_core_multi_rhs_port_seconds,port_core_schur_product_port_seconds,"
        "port_core_core_accumulation_seconds,port_core_core_csr_seconds,"
        "port_core_core_factor_seconds,"
        "interface_solver,interface_krylov_requested,adaptive_interface_tolerance,"
        "interface_krylov_actual,interface_krylov_fallback_steps,"
        "adaptive_interface_retry_steps,adaptive_interface_retry_iterations,"
        "adaptive_interface_retry_seconds,interface_iterations_total,"
        "interface_iterations_maximum,"
        "interface_matvecs,maximum_interface_relative_residual,coarse_dimension,"
        "geometric_coarse_dimension,operator_coarse_dimension,"
        "operator_coarse_cache_hit,operator_coarse_setup_seconds,"
        "operator_coarse_cache_load_seconds,operator_coarse_cache_save_seconds,"
        "interface_predictor_applied_steps,interface_predictor_accepted_steps,"
        "interface_predictor_seconds,"
        "proxy_colors,proxy_probing_applies,proxy_nonzeros,proxy_setup_seconds,"
        "proxy_symbolic_seconds,proxy_numerical_seconds,proxy_solve_seconds,"
        "coarse_solve_seconds,port_forward_solve_seconds,"
        "port_core_solve_seconds,port_back_substitution_seconds,"
        "interface_solve_seconds,"
        "interface_operator_seconds,interface_preconditioner_seconds,"
        "interface_orthogonalization_seconds,interface_vector_update_seconds,"
        "native_reduced_history_enabled,native_reduced_history_steps,"
        "native_reduced_rhs_seconds,step_rhs_seconds,full_residual_seconds,"
        "local_recovery_seconds,"
        "local_online_core_seconds,"
        "fom_factor_seconds,fom_solve_seconds,dynamic_schur_factor_memory_bytes,"
        "fom_factor_memory_bytes,factor_memory_bytes,model_bytes,peak_working_set_bytes,"
        "operator_preparation_seconds,descriptor_fingerprint_seconds,"
        "descriptor_cache_hit,descriptor_cache_load_seconds,"
        "descriptor_cache_save_seconds,descriptor_cache_bytes,"
        "descriptor_assembly_seconds,"
        "interface_partition_seconds,matrix_partition_seconds,"
        "post_local_setup_seconds,"
        "post_schur_setup_seconds,full_step_assembly_seconds,"
        "time_stepping_seconds,"
        "post_time_stepping_seconds,total_seconds,"
        "port_reduction,port_dimension,port_snapshot_seconds,port_basis_seconds,"
        "port_schur_assembly_seconds,port_schur_factor_seconds,"
        "port_local_full_interface_pilot_seconds,"
        "enrichment_rounds,enrichment_added_rank,enrichment_factorization_seconds,"
        "enrichment_solve_seconds,enrichment_orthogonalization_seconds,"
        "enrichment_total_seconds,corrected_relative_l2,corrected_maximum_absolute_k,"
        "corrected_maximum_full_residual,corrected_solve_seconds,"
        "maximum_port_projection_relative_error,"
        "maximum_port_reduced_relative_residual,port_schur_relative_asymmetry,"
        "port_basis_method,port_space_method,optimal_port_source_mode,"
        "optimal_port_ablation,port_snapshot_used,"
        "port_fom_used_for_basis,fom_comparison_enabled,"
        "local_solve_threads,local_pardiso_threads,interface_factor_threads,"
        "proxy_cache_enabled,"
        "proxy_matrix_cache_hit,proxy_factor_cache_hit,"
        "construction_trace_mode,construction_trace_setup_seconds,"
        "construction_pardiso_threads,"
        "construction_global_factor_seconds,"
        "construction_global_solve_seconds,"
        "operator_coarse_trace_dimension,"
        "operator_trace_krylov_iterations,"
        "operator_trace_krylov_maximum_iterations,"
        "operator_trace_krylov_maximum_relative_residual,"
        "analytic_reference_used,"
        "analytic_reference_relative_residual,"
        "global_construction_factor_used\n"
        << std::setprecision(17)
        << finalStatus
        << ',' << (reducedMethod
                ? "Local Interior Arnoldi + Reduced Dynamic Schur"
                : (optimalTransfer
                ? "Local Block Arnoldi + Operator-Informed Mandatory Port Modes + Transfer Spectral Enrichment + Reduced Dynamic Schur"
                : (options.globalRandomizedSchur
                    ? "Local Block Arnoldi + Global Randomized Schur Port Space + Reduced Dynamic Schur"
                : (operatorInformedPort
                    ? "Local Block Arnoldi with Operator-Informed Port Space and Schur-Residual Krylov Enrichment"
                : (options.portReduction
                ? "Local Block Arnoldi + Port-Reduced Dynamic Schur"
                : "Local Block Arnoldi + Dynamic Schur"))
                )))
        << ',' << locals.size() << ','
        << descriptor.dofs << ',' << partition.interfaceGlobalDofs.size() << ','
        << totalRank << ',' << uniqueTemplates << ',' << reusedInstances << ','
        << (options.reuseIdenticalSubdomains ? 1 : 0) << ','
        << waveform.name << ',' << options.initialMode << ','
        << options.timeStep << ',' << steps << ','
        << steadyStateDetectedStep << ',' << steadyStateReusedSteps << ','
        << lastStateChangeRelative << ',' << accuracy.spaceTimeRelativeL2 << ','
        << accuracy.maximumAbsolute << ',' << accuracy.maximumTemperatureError << ','
        << accuracy.maximumFullResidual << ','
        << accuracy.maximumFullResidualBeforeGate << ','
        << options.fullResidualTolerance << ','
        << (options.fullResidualFallback ? 1 : 0) << ','
        << (accuracy.residualGateAllPassed ? 1 : 0) << ','
        << accuracy.residualToleranceViolationSteps << ','
        << accuracy.residualFallbackSteps << ','
        << accuracy.residualFallbackSolveSeconds << ','
        << accuracy.maximumReducedResidual << ','
        << accuracy.maximumJump << ',' << accuracy.maximumFluxImbalance << ','
        << accuracy.maximumFluxRelativeL2 << ','
        << referenceSeconds << ','
        << (referenceCacheHit ? 1 : 0) << ','
        << referenceCacheLoadSeconds << ','
        << referenceCacheSaveSeconds << ','
        << referenceCacheBytes << ','
        << localBuildSeconds << ',' << options.secondMomentEnergy << ','
        << options.secondMomentMaximumColumns << ','
        << (localDynamicCacheHit ? 1 : 0) << ','
        << localDynamicCacheLoadSeconds << ','
        << localDynamicCacheSaveSeconds << ','
        << localDynamicCacheBytes << ',' << localSymbolic << ','
        << localNumerical << ',' << localSolve << ',' << localOrthogonalization << ','
        << localTraceBasisSeconds << ',' << localInputSetupSeconds << ','
        << localFingerprintSeconds << ',' << localOrthogonalityAuditSeconds << ','
        << localProjectionSeconds << ','
        << schurSetupSeconds << ',' << schurFactorSeconds << ','
        << schurSymbolicCalls << ',' << schurNumericalCalls << ','
        << (options.portCoreCacheEnabled ? 1 : 0) << ','
        << (matrixFreeDynamicSchur
                && matrixFreeSolver->portCoreCacheHit() ? 1 : 0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->portCoreCacheLoadSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->portCoreCacheSaveSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->portCoreCacheBytes() : 0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->portCorePartitionSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->portCoreCouplingAssemblySeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->portCoreLeafCsrSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->portCoreLeafFactorSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->portCoreEliminationSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->portCoreMultiRhsPortSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->portCoreSchurProductPortSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->portCoreCoreAccumulationSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->portCoreCoreCsrSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->portCoreCoreFactorSeconds() : 0.0) << ','
        << interfaceSolver << ',' << options.interfaceKrylov << ','
        << options.adaptiveInterfaceTolerance << ','
        << accuracy.interfaceKrylovActual << ','
        << accuracy.interfaceKrylovFallbackSteps << ','
        << accuracy.adaptiveInterfaceRetrySteps << ','
        << accuracy.adaptiveInterfaceRetryIterations << ','
        << accuracy.adaptiveInterfaceRetrySeconds << ','
        << accuracy.interfaceIterationsTotal << ','
        << accuracy.interfaceIterationsMaximum << ',' << accuracy.interfaceMatvecs << ','
        << accuracy.maximumInterfaceResidual << ','
        << (globalCoarseModel
                ? globalCoarseModel->rank
                : (matrixFreeDynamicSchur
                    ? matrixFreeSolver->coarseDimension() : 0))
        << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->geometricCoarseDimension() : 0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->operatorCoarseDimension() : 0) << ','
        << (matrixFreeDynamicSchur
                && matrixFreeSolver->operatorCoarseCacheHit() ? 1 : 0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->operatorCoarseSetupSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->operatorCoarseCacheLoadSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur
                ? matrixFreeSolver->operatorCoarseCacheSaveSeconds() : 0.0) << ','
        << accuracy.interfacePredictorAppliedSteps << ','
        << accuracy.interfacePredictorAcceptedSteps << ','
        << accuracy.interfacePredictorSeconds << ','
        << (matrixFreeDynamicSchur ? matrixFreeSolver->proxyColors() : 0) << ','
        << (matrixFreeDynamicSchur ? matrixFreeSolver->proxyProbingApplies() : 0) << ','
        << (matrixFreeDynamicSchur ? matrixFreeSolver->proxyNonzeros() : 0) << ','
        << (matrixFreeDynamicSchur ? matrixFreeSolver->proxySetupSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur ? matrixFreeSolver->proxySymbolicSeconds() : 0.0) << ','
        << (matrixFreeDynamicSchur ? matrixFreeSolver->proxyNumericalSeconds() : 0.0) << ','
        << accuracy.proxySolveSeconds << ',' << accuracy.coarseSolveSeconds << ','
        << accuracy.portForwardSolveSeconds << ','
        << accuracy.portCoreSolveSeconds << ','
        << accuracy.portBackSubstitutionSeconds << ','
        << accuracy.interfaceSolveSeconds << ','
        << accuracy.interfaceOperatorSeconds << ','
        << accuracy.interfacePreconditionerSeconds << ','
        << accuracy.interfaceOrthogonalizationSeconds << ','
        << accuracy.interfaceVectorUpdateSeconds << ','
        << (matrixFreeDynamicSchur && options.nativeReducedHistory ? 1 : 0) << ','
        << accuracy.nativeReducedHistorySteps << ','
        << accuracy.nativeReducedRhsSeconds << ','
        << accuracy.stepRhsSeconds << ',' << accuracy.fullResidualSeconds << ','
        << accuracy.recoverySeconds << ','
        << accuracy.localCoreSeconds << ',' << accuracy.fomFactorSeconds << ','
        << accuracy.fomSolveSeconds << ',' << accuracy.dynamicSchurFactorBytes << ','
        << accuracy.fomFactorBytes << ',' << accuracy.factorBytes << ',' << modelBytes << ','
        << peakBytes << ',' << operatorPreparationSeconds << ','
        << descriptorFingerprintSeconds << ','
        << (descriptorCacheHit ? 1 : 0) << ','
        << descriptorCacheLoadSeconds << ','
        << descriptorCacheSaveSeconds << ','
        << descriptorCacheBytes << ','
        << descriptorAssemblySeconds << ',' << interfacePartitionSeconds << ','
        << matrixPartitionSeconds << ',' << postLocalSetupSeconds << ','
        << postSchurSetupSeconds << ',' << fullStepAssemblySeconds << ','
        << timeSteppingSeconds << ',' << postTimeSteppingSeconds << ','
        << totalSeconds << ','
        << (options.portReduction ? 1 : 0) << ','
        << (portSolver
                ? portSolver->portDimension()
                : (portModel
                    ? portModel->reducedInterfaceDofs : 0))
        << ','
        << portSnapshotSeconds << ',' << portBasisSeconds << ','
        << (portSolver ? portSolver->assemblySeconds() : 0.0) << ','
        << (portSolver ? portSolver->factorizationSeconds() : 0.0) << ','
        << portLocalPilotSeconds << ','
        << options.localPortEnrichmentRounds << ',' << enrichment.addedRank << ','
        << enrichment.factorizationSeconds << ',' << enrichment.solveSeconds << ','
        << enrichment.orthogonalizationSeconds << ',' << enrichment.totalSeconds << ','
        << accuracy.correctedRelativeL2 << ','
        << accuracy.correctedMaximumAbsolute << ','
        << accuracy.correctedMaximumResidual << ','
        << accuracy.correctedSolveSeconds << ','
        << accuracy.maximumPortProjectionError << ','
        << accuracy.maximumPortReducedResidual << ','
        << (portSolver ? portSolver->relativeAsymmetry() : 0.0) << ','
        << options.portBasisMethod << ','
        << (options.globalRandomizedSchur
            ? "Global Randomized Schur Port Space"
            : (optimalTransfer
            ? "Operator-Informed Port Space with Transfer Spectral Enrichment"
            : (hybridRandomized
                ? "Operator-Informed Mandatory + Randomized Transfer + Schur-Residual Enrichment"
                : (randomizedTransfer
                ? "Randomized Transfer Range Finder Port Space"
                : (operatorInformedPort
                ? "Operator-Informed Mandatory Port Basis with Schur-Residual Krylov Enrichment"
                : "not_applicable")))))
        << ','
        << ((optimalTransfer || randomizedTransfer
                || hybridRandomized)
            ? effectiveOptimalPortSourceMode
            : (operatorInformedPort
                ? options.residualKrylovProbeMode
                : "not_applicable"))
        << ','
        << (optimalTransfer
            ? options.optimalPortAblation
            : (hybridRandomized
                ? "hybrid-randomized"
                : (randomizedTransfer
                ? "randomized-range"
                : (operatorInformedPort
                ? options.portBasisMethod
                : (steklovSchur
                    ? "steklov-schur" : "not_applicable")))))
        << ','
        << (portPod ? 1 : 0) << ','
        << (portPod ? 1 : 0) << ','
        << (options.compareFom ? 1 : 0) << ','
        << options.localSolveThreads << ','
        << options.localPardisoThreads << ','
        << ((matrixFreeDynamicSchur
                && options.interfaceKrylov == "augmented-direct")
            ? std::max(options.localPardisoThreads,
                options.localSolveThreads)
            : options.localPardisoThreads) << ','
        << (options.proxyCacheEnabled ? 1 : 0) << ','
        << (matrixFreeDynamicSchur
            && matrixFreeSolver->proxyMatrixCacheHit() ? 1 : 0) << ','
        << (matrixFreeDynamicSchur
            && matrixFreeSolver->proxyFactorCacheHit() ? 1 : 0) << ','
        << options.constructionTraceMode << ','
        << constructionTraceSetupSeconds << ','
        << constructionPardisoThreads << ','
        << constructionGlobalFactorSeconds << ','
        << constructionGlobalSolveSeconds << ','
        << operatorCoarseAggregates << ','
        << operatorTraceKrylovIterations << ','
        << operatorTraceKrylovMaximumIterations << ','
        << operatorTraceKrylovMaximumRelativeResidual << ','
        << (analyticReferenceUsed ? 1 : 0) << ','
        << analyticReferenceRelativeResidual << ','
        << (globalConstructionFactorUsed ? 1 : 0) << '\n';

    if (globalCoarseModel && !options.globalRandomizedSchur) {
        std::ofstream coarseAccuracy(
            outputDirectory
            / "milestone8_global_coarse_single_step.csv");
        coarseAccuracy
            << "case,method,scope,local_port_rank,coarse_rank,"
            "augmented_rank,temperature_relative_l2,"
            "max_error_k,max_temperature_error_k,flux_relative_l2,"
            "interface_residual,full_residual,basis_time_s,"
            "peak_incremental_memory_bytes,corrected,"
            "snapshot_used,fom_used_for_basis,pod_used,svd_used,"
            "status\n"
            << std::setprecision(17)
            << physics.name
            << ",M8.10.1 Local + Global Schur Spectral Coarse,"
            << "one_step,"
            << portModel->reducedInterfaceDofs << ','
            << globalCoarseModel->rank << ','
            << portSolver->portDimension() << ','
            << accuracy.spaceTimeRelativeL2 << ','
            << accuracy.maximumAbsolute << ','
            << accuracy.maximumTemperatureError << ','
            << accuracy.maximumFluxRelativeL2 << ','
            << accuracy.maximumInterfaceResidual << ','
            << accuracy.maximumFullResidual << ','
            << globalCoarseModel->diagnostics.basisSeconds << ','
            << globalCoarseModel->diagnostics
                .peakIncrementalMemoryBytes
            << ",0,0,0,0,0,"
            << ((accuracy.spaceTimeRelativeL2 < 1.0e-4
                    && accuracy.maximumAbsolute < 0.1)
                ? "success" : "accuracy_failed")
            << '\n';
    }
    if (globalCoarseModel && options.globalRandomizedSchur
        && globalRandomizedDiagnostics) {
        std::ofstream randomizedAccuracy(
            outputDirectory
            / "milestone8_global_randomized_single_step.csv");
        randomizedAccuracy
            << "case,method,scope,composition,full_interface_dofs,"
            "local_port_rank,global_port_rank,active_port_rank,"
            "compression_ratio,temperature_relative_l2,max_error_k,"
            "max_temperature_error_k,flux_relative_l2,"
            "interface_residual,full_residual,basis_time_s,"
            "global_schur_apply_count,pardiso_phase33_calls,"
            "global_rhs_count,global_inner_iterations,"
            "mean_solve_time_s,peak_incremental_memory_bytes,"
            "corrected,snapshot_used,fom_used_for_basis,pod_used,"
            "svd_used,training_waveform_used,status\n"
            << std::setprecision(17)
            << physics.name
            << ",M8.11 Global Randomized Schur Port,"
            << "one_step,"
            << options.globalRandomizedComposition << ','
            << partition.interfaceGlobalDofs.size() << ','
            << portModel->reducedInterfaceDofs << ','
            << globalCoarseModel->rank << ','
            << portSolver->portDimension() << ','
            << globalRandomizedDiagnostics->compressionRatio << ','
            << accuracy.spaceTimeRelativeL2 << ','
            << accuracy.maximumAbsolute << ','
            << accuracy.maximumTemperatureError << ','
            << accuracy.maximumFluxRelativeL2 << ','
            << accuracy.maximumInterfaceResidual << ','
            << accuracy.maximumFullResidual << ','
            << globalRandomizedDiagnostics
                ->basisBuildTimeSeconds << ','
            << globalRandomizedDiagnostics
                ->globalSchurApplyCount << ','
            << globalRandomizedDiagnostics
                ->pardisoPhase33Calls << ','
            << globalRandomizedDiagnostics
                ->globalRhsCount << ','
            << globalRandomizedDiagnostics
                ->globalInnerIterations << ','
            << globalRandomizedDiagnostics
                ->meanSolveTimeSeconds << ','
            << globalRandomizedDiagnostics
                ->peakIncrementalMemoryBytes
            << ",0,0,0,0,0,0,"
            << ((accuracy.spaceTimeRelativeL2 < 1.0e-4
                    && accuracy.maximumAbsolute < 0.1)
                ? "success" : "accuracy_failed")
            << '\n';
    }

    std::cout << "Local Block Arnoldi + Dynamic Schur: subdomains=" << locals.size()
              << ", interface=" << partition.interfaceGlobalDofs.size()
              << ", total local rank=" << totalRank
              << ", interface solver=" << interfaceSolver
              << ", max interface iterations=" << accuracy.interfaceIterationsMaximum
              << ", space-time L2=" << accuracy.spaceTimeRelativeL2
              << ", max error=" << accuracy.maximumAbsolute << " K\n";
}

} // namespace mor::transient
