// Assemble and serialize the transient descriptor (K, C, source channels,
// boundary offsets, and fingerprints). Cache loading validates dimensions and
// every physics/mesh fingerprint before exposing the operator to Arnoldi or
// time integration; a warm mismatch fails instead of silently rebuilding.

#include "thermal_descriptor_system.hpp"

#include "config_io.hpp"
#include "mesh_loader.hpp"
#include "fem_assembly.hpp"
#include "linear_solvers.hpp"
#include "mor/model_io.hpp"
#include "mor/source_parameterization.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace mor::transient {
namespace {

using Clock = std::chrono::steady_clock;

template <typename T>
void hashAdd(std::uint64_t& hash, const T& value)
{
    const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
    for (std::size_t i = 0; i < sizeof(T); ++i) {
        hash ^= bytes[i];
        hash *= UINT64_C(1099511628211);
    }
}

void hashString(std::uint64_t& hash, const std::string& value)
{
    hashAdd(hash, static_cast<std::uint64_t>(value.size()));
    for (const unsigned char byte : value) hashAdd(hash, byte);
}

void hashVec3(std::uint64_t& hash, const Vec3& value)
{
    hashAdd(hash, value.x);
    hashAdd(hash, value.y);
    hashAdd(hash, value.z);
}

void hashMaterial(std::uint64_t& hash, const Material& material)
{
    hashString(hash, material.name);
    hashAdd(hash, material.conductivity);
    hashAdd(hash, material.conductivityX);
    hashAdd(hash, material.conductivityY);
    hashAdd(hash, material.conductivityZ);
    hashAdd(hash, material.density);
    hashAdd(hash, material.heatCapacity);
}

std::uint64_t hashMatrix(const SparseMatrix& matrix)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    hashAdd(hash, matrix.n);
    matrix.forEachEntry([&](int row, int column, double value) {
        hashAdd(hash, row);
        hashAdd(hash, column);
        hashAdd(hash, value);
    });
    return hash;
}

std::uint64_t hashVector(const std::vector<double>& vector)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    const std::uint64_t size = static_cast<std::uint64_t>(vector.size());
    hashAdd(hash, size);
    for (double value : vector) hashAdd(hash, value);
    return hash;
}

double symmetryError(const SparseMatrix& matrix)
{
    double numerator = 0.0;
    double denominator = 0.0;
    for (int row = 0; row < matrix.size(); ++row) {
        for (int offset = matrix.rowPtr[static_cast<std::size_t>(row)];
             offset < matrix.rowPtr[static_cast<std::size_t>(row + 1)]; ++offset) {
            const int column = matrix.colInd[static_cast<std::size_t>(offset)];
            const double value = matrix.values[static_cast<std::size_t>(offset)];
            const int begin = matrix.rowPtr[static_cast<std::size_t>(column)];
            const int end = matrix.rowPtr[static_cast<std::size_t>(column + 1)];
            const auto found = std::lower_bound(
                matrix.colInd.begin() + begin, matrix.colInd.begin() + end, row);
            const double transpose = found != matrix.colInd.begin() + end && *found == row
                ? matrix.values[static_cast<std::size_t>(
                    std::distance(matrix.colInd.begin(), found))] : 0.0;
            const double delta = value - transpose;
            numerator += delta * delta;
            denominator += value * value;
        }
    }
    return std::sqrt(numerator) /
        std::max(1.0e-300, std::sqrt(denominator));
}

MatrixDiagnostic diagnose(const SparseMatrix& matrix,
                          const Mesh& mesh,
                          const std::string& name,
                          const std::string& units,
                          bool capacity)
{
    MatrixDiagnostic result;
    result.name = name;
    result.units = units;
    result.rows = matrix.size();
    result.nonzeros = matrix.values.size();
    result.symmetryError = symmetryError(matrix);
    result.minimumDiagonal = std::numeric_limits<double>::infinity();
    result.maximumDiagonal = -std::numeric_limits<double>::infinity();
    result.positiveOnFreeDofs = true;
    for (int row = 0; row < matrix.size(); ++row) {
        const double diagonal = matrix.diagonal(row);
        result.minimumDiagonal = std::min(result.minimumDiagonal, diagonal);
        result.maximumDiagonal = std::max(result.maximumDiagonal, diagonal);
        if (std::abs(diagonal) <= 1.0e-30) ++result.zeroDiagonalCount;
        const bool constrainedCapacity = capacity
            && mesh.nodes[static_cast<std::size_t>(row)].dirichlet;
        if (!constrainedCapacity && !(diagonal > 0.0)) {
            result.positiveOnFreeDofs = false;
        }
    }
    return result;
}

void removeDirichletCapacity(const Mesh& mesh, SparseMatrix& capacity)
{
    std::vector<MatrixEntry> filtered;
    filtered.reserve(capacity.triplets.size());
    capacity.forEachEntry([&](int row, int column, double value) {
        if (!mesh.nodes[static_cast<std::size_t>(row)].dirichlet
            && !mesh.nodes[static_cast<std::size_t>(column)].dirichlet) {
            filtered.push_back({row, column, value});
        }
    });
    capacity.triplets = std::move(filtered);
    capacity.rowPtr.clear();
    capacity.colInd.clear();
    capacity.values.clear();
    capacity.csrReady = false;
}

void lumpCapacity(SparseMatrix& capacity)
{
    std::vector<double> rowSum(static_cast<std::size_t>(capacity.size()), 0.0);
    capacity.forEachEntry([&](int row, int, double value) {
        rowSum[static_cast<std::size_t>(row)] += value;
    });
    capacity.triplets.clear();
    capacity.rowPtr.clear();
    capacity.colInd.clear();
    capacity.values.clear();
    capacity.csrReady = false;
    for (int row = 0; row < capacity.size(); ++row) {
        if (std::abs(rowSum[static_cast<std::size_t>(row)]) > 0.0) {
            capacity.add(row, row, rowSum[static_cast<std::size_t>(row)]);
        }
    }
}

constexpr std::uint64_t descriptorCacheMagic =
    UINT64_C(0x5448444553434143); // "THDESCAC"
constexpr std::uint32_t descriptorCacheVersion = 1;
// Bump whenever the mathematical descriptor assembly changes in a way not
// represented by the hashed inputs below.
constexpr std::uint32_t descriptorAssemblyAlgorithmVersion = 1;

template <typename T>
void writeScalar(std::ofstream& output, const T& value)
{
    static_assert(std::is_trivially_copyable_v<T>);
    output.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template <typename T>
T readScalar(std::ifstream& input, const char* field)
{
    static_assert(std::is_trivially_copyable_v<T>);
    T value{};
    input.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!input) {
        throw std::runtime_error(
            std::string("[Thermal descriptor cache] Truncated field: ") + field);
    }
    return value;
}

void writeString(std::ofstream& output, const std::string& value)
{
    writeScalar(output, static_cast<std::uint64_t>(value.size()));
    if (!value.empty()) {
        output.write(value.data(), static_cast<std::streamsize>(value.size()));
    }
}

std::string readString(std::ifstream& input, const char* field)
{
    const std::uint64_t count = readScalar<std::uint64_t>(input, field);
    if (count > 4096) {
        throw std::runtime_error(
            std::string("[Thermal descriptor cache] Invalid string length: ") + field);
    }
    std::string value(static_cast<std::size_t>(count), '\0');
    if (!value.empty()) {
        input.read(value.data(), static_cast<std::streamsize>(value.size()));
        if (!input) {
            throw std::runtime_error(
                std::string("[Thermal descriptor cache] Truncated string: ") + field);
        }
    }
    return value;
}

template <typename T>
void writeVector(std::ofstream& output, const std::vector<T>& values)
{
    static_assert(std::is_trivially_copyable_v<T>);
    writeScalar(output, static_cast<std::uint64_t>(values.size()));
    if (!values.empty()) {
        output.write(reinterpret_cast<const char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(T)));
    }
}

template <typename T>
std::vector<T> readVectorLimited(std::ifstream& input,
                                 std::uint64_t maximumCount,
                                 const char* field)
{
    static_assert(std::is_trivially_copyable_v<T>);
    const std::uint64_t count = readScalar<std::uint64_t>(input, field);
    if (count > maximumCount
        || count > static_cast<std::uint64_t>(
            std::numeric_limits<std::size_t>::max() / sizeof(T))) {
        throw std::runtime_error(
            std::string("[Thermal descriptor cache] Invalid vector length: ") + field);
    }
    std::vector<T> values(static_cast<std::size_t>(count));
    if (!values.empty()) {
        input.read(reinterpret_cast<char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(T)));
        if (!input) {
            throw std::runtime_error(
                std::string("[Thermal descriptor cache] Truncated vector: ") + field);
        }
    }
    return values;
}

template <typename T>
std::vector<T> readVectorExact(std::ifstream& input,
                               std::uint64_t expectedCount,
                               const char* field)
{
    std::vector<T> values = readVectorLimited<T>(input, expectedCount, field);
    if (values.size() != static_cast<std::size_t>(expectedCount)) {
        throw std::runtime_error(
            std::string("[Thermal descriptor cache] Dimension mismatch: ") + field);
    }
    return values;
}

void writeMatrix(std::ofstream& output, const SparseMatrix& matrix)
{
    if (!matrix.csrReady) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Matrix must be finalized before save.");
    }
    writeScalar(output, matrix.n);
    writeVector(output, matrix.rowPtr);
    writeVector(output, matrix.colInd);
    writeVector(output, matrix.values);
}

SparseMatrix readMatrix(std::ifstream& input,
                        int expectedDofs,
                        const char* field)
{
    const int rows = readScalar<int>(input, field);
    if (rows != expectedDofs || rows <= 0) {
        throw std::runtime_error(
            std::string("[Thermal descriptor cache] Matrix row mismatch: ") + field);
    }
    SparseMatrix matrix(rows);
    matrix.rowPtr = readVectorExact<int>(
        input, static_cast<std::uint64_t>(rows) + 1, "matrix row pointer");
    matrix.colInd = readVectorLimited<int>(input,
        static_cast<std::uint64_t>(std::numeric_limits<int>::max()),
        "matrix column indices");
    matrix.values = readVectorExact<double>(input,
        static_cast<std::uint64_t>(matrix.colInd.size()), "matrix values");
    if (matrix.rowPtr.empty() || matrix.rowPtr.front() != 0
        || matrix.rowPtr.back() != static_cast<int>(matrix.colInd.size())) {
        throw std::runtime_error(
            std::string("[Thermal descriptor cache] Invalid CSR row pointer: ") + field);
    }
    for (int row = 0; row < rows; ++row) {
        const int begin = matrix.rowPtr[static_cast<std::size_t>(row)];
        const int end = matrix.rowPtr[static_cast<std::size_t>(row + 1)];
        if (begin < 0 || begin > end
            || end > static_cast<int>(matrix.colInd.size())) {
            throw std::runtime_error(
                std::string("[Thermal descriptor cache] Invalid CSR offsets: ") + field);
        }
        int previous = -1;
        for (int offset = begin; offset < end; ++offset) {
            const int column = matrix.colInd[static_cast<std::size_t>(offset)];
            const double value = matrix.values[static_cast<std::size_t>(offset)];
            if (column < 0 || column >= rows || column <= previous
                || !std::isfinite(value)) {
                throw std::runtime_error(
                    std::string("[Thermal descriptor cache] Invalid CSR payload: ") + field);
            }
            previous = column;
        }
    }
    matrix.csrReady = true;
    return matrix;
}

void writeDiagnostic(std::ofstream& output, const MatrixDiagnostic& value)
{
    writeString(output, value.name);
    writeString(output, value.units);
    writeScalar(output, value.rows);
    writeScalar(output, static_cast<std::uint64_t>(value.nonzeros));
    writeScalar(output, value.symmetryError);
    writeScalar(output, value.minimumDiagonal);
    writeScalar(output, value.maximumDiagonal);
    writeScalar(output, value.zeroDiagonalCount);
    writeScalar(output, static_cast<std::uint8_t>(value.positiveOnFreeDofs));
}

MatrixDiagnostic readDiagnostic(std::ifstream& input, int expectedDofs)
{
    MatrixDiagnostic value;
    value.name = readString(input, "diagnostic name");
    value.units = readString(input, "diagnostic units");
    value.rows = readScalar<int>(input, "diagnostic rows");
    value.nonzeros = static_cast<std::size_t>(
        readScalar<std::uint64_t>(input, "diagnostic nonzeros"));
    value.symmetryError = readScalar<double>(input, "diagnostic symmetry");
    value.minimumDiagonal = readScalar<double>(input, "diagnostic minimum diagonal");
    value.maximumDiagonal = readScalar<double>(input, "diagnostic maximum diagonal");
    value.zeroDiagonalCount = readScalar<int>(input, "diagnostic zero diagonal count");
    value.positiveOnFreeDofs =
        readScalar<std::uint8_t>(input, "diagnostic positivity") != 0;
    if (value.rows != expectedDofs || value.nonzeros >
        static_cast<std::size_t>(std::numeric_limits<int>::max())
        || !std::isfinite(value.symmetryError)
        || !std::isfinite(value.minimumDiagonal)
        || !std::isfinite(value.maximumDiagonal)) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Invalid matrix diagnostic payload.");
    }
    return value;
}

} // namespace

std::uint64_t thermalDescriptorInputFingerprint(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::string& massType,
    bool collectInterfaceGram)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    hashAdd(hash, descriptorAssemblyAlgorithmVersion);
    hashString(hash, massType);
    hashAdd(hash, static_cast<std::uint8_t>(collectInterfaceGram));

    hashAdd(hash, static_cast<std::uint64_t>(mesh.nodes.size()));
    for (const Node& node : mesh.nodes) {
        hashVec3(hash, node.p);
        hashAdd(hash, node.subdomain);
        hashAdd(hash, node.sourceVertex);
        hashAdd(hash, static_cast<std::uint8_t>(node.dirichlet));
        hashAdd(hash, node.dirichletValue);
    }
    hashAdd(hash, static_cast<std::uint64_t>(mesh.tets.size()));
    for (const Tet& tet : mesh.tets) {
        for (int value : tet.v) hashAdd(hash, value);
        for (int value : tet.source) hashAdd(hash, value);
        for (int value : tet.dof) hashAdd(hash, value);
        hashAdd(hash, tet.subdomain);
        hashAdd(hash, tet.domainEntity);
    }
    hashAdd(hash, static_cast<std::uint64_t>(mesh.boundaryFaces.size()));
    for (const BoundaryFace& face : mesh.boundaryFaces) {
        hashAdd(hash, face.tet);
        hashAdd(hash, face.subdomain);
        hashAdd(hash, face.boundaryEntity);
        for (int value : face.local) hashAdd(hash, value);
        for (const Vec3& point : face.points) hashVec3(hash, point);
        hashVec3(hash, face.normal);
    }
    hashAdd(hash, static_cast<std::uint64_t>(mesh.interfaceFaces.size()));
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        hashAdd(hash, face.leftTet);
        hashAdd(hash, face.rightTet);
        hashAdd(hash, face.leftFaceId);
        hashAdd(hash, face.rightFaceId);
        hashAdd(hash, face.leftBoundaryEntity);
        hashAdd(hash, face.rightBoundaryEntity);
        for (int value : face.leftLocal) hashAdd(hash, value);
        for (int value : face.rightLocal) hashAdd(hash, value);
        hashVec3(hash, face.leftNormal);
        hashVec3(hash, face.rightNormal);
        hashAdd(hash, static_cast<std::uint64_t>(face.integrationTriangles.size()));
        for (const auto& triangle : face.integrationTriangles) {
            for (const Vec3& point : triangle) hashVec3(hash, point);
        }
        hashAdd(hash, face.overlapPolygonVertices);
        hashAdd(hash, face.overlapArea);
    }

    hashAdd(hash, physics.coordinateScale);
    hashString(hash, physics.dirichletMethod);
    hashAdd(hash, physics.nitschePenaltyFactor);
    hashString(hash, physics.interfaceScheme);
    hashString(hash, physics.penaltyMode);
    hashString(hash, physics.penaltyScaling);
    hashAdd(hash, physics.penaltyFactor);
    hashAdd(hash, physics.thermalSourceScale);
    hashAdd(hash, static_cast<std::uint64_t>(physics.domains.size()));
    for (const DomainConfig& domain : physics.domains) {
        hashMaterial(hash, domain.material);
        hashVec3(hash, domain.translationMeters);
        hashAdd(hash,
            static_cast<std::uint64_t>(domain.materialsByDomainEntity.size()));
        for (const auto& item : domain.materialsByDomainEntity) {
            hashAdd(hash, item.first);
            hashMaterial(hash, item.second);
        }
    }
    hashAdd(hash,
        static_cast<std::uint64_t>(physics.dirichletConditions.size()));
    for (const BoundaryCondition& condition : physics.dirichletConditions) {
        hashAdd(hash, condition.subdomain);
        hashAdd(hash, condition.boundaryEntity);
        hashAdd(hash, condition.temperature);
    }
    hashAdd(hash,
        static_cast<std::uint64_t>(physics.convectionConditions.size()));
    for (const ConvectionCondition& condition : physics.convectionConditions) {
        hashAdd(hash, condition.subdomain);
        hashAdd(hash, condition.boundaryEntity);
        hashAdd(hash, condition.coefficient);
        hashAdd(hash, condition.ambientTemperature);
    }
    hashAdd(hash, static_cast<std::uint64_t>(physics.heatSources.size()));
    for (const HeatSource& source : physics.heatSources) {
        hashAdd(hash, source.subdomain);
        hashAdd(hash, source.domainEntity);
        hashAdd(hash, source.heatRateW);
    }
    return hash;
}

ThermalDescriptorSystem assembleThermalDescriptorSystem(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::string& massType,
    bool collectInterfaceGram)
{
    const auto start = Clock::now();
    if (massType != "consistent" && massType != "lumped") {
        throw std::runtime_error(
            "Transient MOR mass type must be consistent or lumped.");
    }
    ThermalDescriptorSystem result;
    result.dofs = static_cast<int>(mesh.nodes.size());
    result.massType = massType;
    result.capacity = SparseMatrix(result.dofs);
    result.conductivity = SparseMatrix(result.dofs);
    std::vector<double> assembledSource(
        static_cast<std::size_t>(result.dofs), 0.0);
    assembleVolume(mesh, physics, &result.capacity,
                   result.conductivity, assembledSource);
    if (physics.thermalSourceScale != 1.0) {
        for (double& value : assembledSource) value *= physics.thermalSourceScale;
    }
    const std::vector<double> heatOnlySource = assembledSource;
    assembleConvectionBoundaries(
        mesh, physics, result.conductivity, assembledSource, false, nullptr);
    assembleSipgInterface(mesh, physics, result.conductivity, false, true, nullptr);
    if (collectInterfaceGram) {
        AssemblyDiagnostics interfaceDiagnostics;
        interfaceDiagnostics.interfaceConsistencyDiag.assign(
            static_cast<std::size_t>(result.dofs), 0.0);
        interfaceDiagnostics.interfacePenaltyDiag.assign(
            static_cast<std::size_t>(result.dofs), 0.0);
        interfaceDiagnostics.interfaceTraceMassDiag.assign(
            static_cast<std::size_t>(result.dofs), 0.0);
        assembleSipgInterface(mesh, physics, result.conductivity, true, false,
                              &interfaceDiagnostics);
        result.interfaceTraceMassDiagonal =
            std::move(interfaceDiagnostics.interfaceTraceMassDiag);
        result.interfacePenaltyMassDiagonal =
            std::move(interfaceDiagnostics.interfacePenaltyDiag);
    } else {
        assembleSipgInterface(
            mesh, physics, result.conductivity, true, false, nullptr);
    }

    std::vector<double> fixedAdjust;
    if (physics.dirichletMethod == "nitsche") {
        assembleNitscheDirichletBoundaries(
            mesh, physics, result.conductivity, assembledSource, nullptr);
    } else if (physics.dirichletMethod == "strong") {
        fixedAdjust = makeDirichletAdjustedSystem(mesh, result.conductivity);
        removeDirichletCapacity(mesh, result.capacity);
    } else {
        throw std::runtime_error(
            "Transient MOR requires the existing strong or Nitsche Dirichlet path.");
    }
    if (massType == "lumped") lumpCapacity(result.capacity);
    result.capacity.finalizeCsr();
    result.conductivity.finalizeCsr();

    const SourceParameterization sources = buildSourceParameterization(
        mesh, physics, assembledSource, heatOnlySource, fixedAdjust);
    result.sourceChannels = static_cast<int>(sources.channels.size());
    result.boundaryRhs = sources.referenceRhs;
    result.input.assign(static_cast<std::size_t>(result.dofs)
        * static_cast<std::size_t>(result.sourceChannels), 0.0);
    for (int channel = 0; channel < result.sourceChannels; ++channel) {
        const SourceChannel& source = sources.channels[static_cast<std::size_t>(channel)];
        std::copy(source.rhsPerWatt.begin(), source.rhsPerWatt.end(),
                  result.input.begin() + static_cast<std::ptrdiff_t>(
                      static_cast<std::size_t>(channel) * result.dofs));
        result.nominalPowersW.push_back(source.nominalPowerW);
        result.minimumPowersW.push_back(source.minimumPowerW);
        result.maximumPowersW.push_back(source.maximumPowerW);
        result.sourceSubdomains.push_back(source.subdomain);
        result.sourceDomainEntities.push_back(source.domainEntity);
    }
    result.deploymentDofs.reserve(mesh.nodes.size());
    for (const Node& node : mesh.nodes) {
        result.deploymentDofs.push_back({
            node.p.x, node.p.y, node.p.z, node.subdomain, node.sourceVertex});
    }

    const Fingerprints common = computeFingerprints(
        mesh, result.conductivity, physics, {});
    result.fingerprints.mesh = common.mesh;
    result.fingerprints.sources = common.sources;
    result.fingerprints.capacity = hashMatrix(result.capacity);
    result.fingerprints.conductivity = hashMatrix(result.conductivity);
    result.fingerprints.input = hashVector(result.input);
    result.fingerprints.boundary = hashVector(result.boundaryRhs);
    result.capacityDiagnostic = diagnose(
        result.capacity, mesh, "capacity", "J/K", true);
    result.conductivityDiagnostic = diagnose(
        result.conductivity, mesh, "conductivity", "W/K", false);
    result.assemblySeconds = std::chrono::duration<double>(
        Clock::now() - start).count();
    return result;
}

void saveThermalDescriptorCache(
    const std::filesystem::path& path,
    std::uint64_t inputFingerprint,
    const ThermalDescriptorSystem& system)
{
    if (path.empty()) return;
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    std::filesystem::path temporary = path;
    temporary += ".tmp";
    std::ofstream output(temporary, std::ios::binary | std::ios::trunc);
    if (!output) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Cannot create temporary cache file.");
    }
    writeScalar(output, descriptorCacheMagic);
    writeScalar(output, descriptorCacheVersion);
    writeScalar(output, inputFingerprint);
    writeScalar(output, system.dofs);
    writeScalar(output, system.sourceChannels);
    writeString(output, system.massType);
    writeMatrix(output, system.capacity);
    writeMatrix(output, system.conductivity);
    writeVector(output, system.interfaceTraceMassDiagonal);
    writeVector(output, system.interfacePenaltyMassDiagonal);
    writeVector(output, system.input);
    writeVector(output, system.boundaryRhs);
    writeVector(output, system.nominalPowersW);
    writeVector(output, system.minimumPowersW);
    writeVector(output, system.maximumPowersW);
    writeVector(output, system.sourceSubdomains);
    writeVector(output, system.sourceDomainEntities);
    writeScalar(output,
        static_cast<std::uint64_t>(system.deploymentDofs.size()));
    for (const DeploymentDof& dof : system.deploymentDofs) {
        writeScalar(output, dof.x);
        writeScalar(output, dof.y);
        writeScalar(output, dof.z);
        writeScalar(output, dof.subdomain);
        writeScalar(output, dof.sourceVertex);
    }
    writeScalar(output, system.fingerprints.mesh);
    writeScalar(output, system.fingerprints.capacity);
    writeScalar(output, system.fingerprints.conductivity);
    writeScalar(output, system.fingerprints.input);
    writeScalar(output, system.fingerprints.boundary);
    writeScalar(output, system.fingerprints.sources);
    writeDiagnostic(output, system.capacityDiagnostic);
    writeDiagnostic(output, system.conductivityDiagnostic);
    writeScalar(output, system.assemblySeconds);
    output.close();
    if (!output) {
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        throw std::runtime_error(
            "[Thermal descriptor cache] Failed while writing cache payload.");
    }
    std::error_code removeError;
    std::filesystem::remove(path, removeError);
    if (removeError) {
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        throw std::runtime_error(
            "[Thermal descriptor cache] Cannot replace existing cache file: "
            + removeError.message());
    }
    std::error_code renameError;
    std::filesystem::rename(temporary, path, renameError);
    if (renameError) {
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        throw std::runtime_error(
            "[Thermal descriptor cache] Cannot publish cache file: "
            + renameError.message());
    }
}

ThermalDescriptorSystem loadThermalDescriptorCache(
    const std::filesystem::path& path,
    std::uint64_t expectedInputFingerprint,
    int expectedDofs)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Cannot open cache file.");
    }
    const std::uint64_t magic =
        readScalar<std::uint64_t>(input, "magic");
    const std::uint32_t version =
        readScalar<std::uint32_t>(input, "version");
    const std::uint64_t fingerprint =
        readScalar<std::uint64_t>(input, "input fingerprint");
    if (magic != descriptorCacheMagic || version != descriptorCacheVersion) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Unsupported cache format; reuse refused.");
    }
    if (fingerprint != expectedInputFingerprint) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Assembly-input fingerprint mismatch; reuse refused.");
    }
    ThermalDescriptorSystem system;
    system.dofs = readScalar<int>(input, "DOF count");
    system.sourceChannels = readScalar<int>(input, "source-channel count");
    if (system.dofs != expectedDofs || system.dofs <= 0
        || system.sourceChannels < 0
        || system.sourceChannels > 1000000) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Header dimension mismatch; reuse refused.");
    }
    system.massType = readString(input, "mass type");
    if (system.massType != "consistent" && system.massType != "lumped") {
        throw std::runtime_error(
            "[Thermal descriptor cache] Invalid mass type.");
    }
    system.capacity = readMatrix(input, expectedDofs, "capacity");
    system.conductivity = readMatrix(input, expectedDofs, "conductivity");
    system.interfaceTraceMassDiagonal = readVectorLimited<double>(
        input, static_cast<std::uint64_t>(expectedDofs), "interface trace mass");
    system.interfacePenaltyMassDiagonal = readVectorLimited<double>(
        input, static_cast<std::uint64_t>(expectedDofs), "interface penalty mass");
    if ((!system.interfaceTraceMassDiagonal.empty()
            && system.interfaceTraceMassDiagonal.size()
                != static_cast<std::size_t>(expectedDofs))
        || (!system.interfacePenaltyMassDiagonal.empty()
            && system.interfacePenaltyMassDiagonal.size()
                != static_cast<std::size_t>(expectedDofs))) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Interface Gram dimension mismatch.");
    }
    const std::uint64_t inputCount =
        static_cast<std::uint64_t>(expectedDofs)
        * static_cast<std::uint64_t>(system.sourceChannels);
    if (system.sourceChannels > 0
        && inputCount / static_cast<std::uint64_t>(system.sourceChannels)
            != static_cast<std::uint64_t>(expectedDofs)) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Input dimension overflow.");
    }
    system.input = readVectorExact<double>(input, inputCount, "input operator");
    system.boundaryRhs = readVectorExact<double>(
        input, static_cast<std::uint64_t>(expectedDofs), "boundary RHS");
    const std::uint64_t channels =
        static_cast<std::uint64_t>(system.sourceChannels);
    system.nominalPowersW = readVectorExact<double>(
        input, channels, "nominal powers");
    system.minimumPowersW = readVectorExact<double>(
        input, channels, "minimum powers");
    system.maximumPowersW = readVectorExact<double>(
        input, channels, "maximum powers");
    system.sourceSubdomains = readVectorExact<int>(
        input, channels, "source subdomains");
    system.sourceDomainEntities = readVectorExact<int>(
        input, channels, "source domain entities");
    const std::uint64_t deploymentCount =
        readScalar<std::uint64_t>(input, "deployment DOF count");
    if (deploymentCount != static_cast<std::uint64_t>(expectedDofs)) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Deployment DOF dimension mismatch.");
    }
    system.deploymentDofs.resize(static_cast<std::size_t>(deploymentCount));
    for (DeploymentDof& dof : system.deploymentDofs) {
        dof.x = readScalar<double>(input, "deployment x");
        dof.y = readScalar<double>(input, "deployment y");
        dof.z = readScalar<double>(input, "deployment z");
        dof.subdomain = readScalar<int>(input, "deployment subdomain");
        dof.sourceVertex = readScalar<int>(input, "deployment source vertex");
        if (!std::isfinite(dof.x) || !std::isfinite(dof.y)
            || !std::isfinite(dof.z)) {
            throw std::runtime_error(
                "[Thermal descriptor cache] Non-finite deployment coordinate.");
        }
    }
    system.fingerprints.mesh =
        readScalar<std::uint64_t>(input, "mesh fingerprint");
    system.fingerprints.capacity =
        readScalar<std::uint64_t>(input, "capacity fingerprint");
    system.fingerprints.conductivity =
        readScalar<std::uint64_t>(input, "conductivity fingerprint");
    system.fingerprints.input =
        readScalar<std::uint64_t>(input, "input fingerprint");
    system.fingerprints.boundary =
        readScalar<std::uint64_t>(input, "boundary fingerprint");
    system.fingerprints.sources =
        readScalar<std::uint64_t>(input, "source fingerprint");
    if (system.fingerprints.capacity != hashMatrix(system.capacity)
        || system.fingerprints.conductivity
            != hashMatrix(system.conductivity)
        || system.fingerprints.input != hashVector(system.input)
        || system.fingerprints.boundary != hashVector(system.boundaryRhs)) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Payload fingerprint mismatch; reuse refused.");
    }
    const auto finiteVector = [](const std::vector<double>& values) {
        return std::all_of(values.begin(), values.end(),
            [](double value) { return std::isfinite(value); });
    };
    if (!finiteVector(system.interfaceTraceMassDiagonal)
        || !finiteVector(system.interfacePenaltyMassDiagonal)
        || !finiteVector(system.input)
        || !finiteVector(system.boundaryRhs)
        || !finiteVector(system.nominalPowersW)
        || !finiteVector(system.minimumPowersW)
        || !finiteVector(system.maximumPowersW)) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Non-finite vector payload; reuse refused.");
    }
    system.capacityDiagnostic = readDiagnostic(input, expectedDofs);
    system.conductivityDiagnostic = readDiagnostic(input, expectedDofs);
    system.assemblySeconds =
        readScalar<double>(input, "original assembly seconds");
    if (!std::isfinite(system.assemblySeconds)
        || system.capacityDiagnostic.nonzeros != system.capacity.values.size()
        || system.conductivityDiagnostic.nonzeros
            != system.conductivity.values.size()) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Diagnostic/matrix mismatch.");
    }
    const int trailing = input.peek();
    if (trailing != std::char_traits<char>::eof()) {
        throw std::runtime_error(
            "[Thermal descriptor cache] Unexpected trailing payload.");
    }
    return system;
}

std::vector<double> descriptorInputRhs(
    const ThermalDescriptorSystem& system,
    const std::vector<double>& powersW)
{
    if (powersW.size() != static_cast<std::size_t>(system.sourceChannels)) {
        throw std::runtime_error("Transient power vector dimension mismatch.");
    }
    std::vector<double> rhs = system.boundaryRhs;
    for (int channel = 0; channel < system.sourceChannels; ++channel) {
        const double power = powersW[static_cast<std::size_t>(channel)];
        const double* column = system.input.data()
            + static_cast<std::size_t>(channel) * system.dofs;
        for (int row = 0; row < system.dofs; ++row) {
            rhs[static_cast<std::size_t>(row)] += power * column[row];
        }
    }
    return rhs;
}

void writeTransientMatrixDiagnostics(
    const ThermalDescriptorSystem& system,
    const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "matrix,rows,nonzeros,units,symmetry_error,minimum_diagonal,"
        << "maximum_diagonal,zero_diagonal_count,positive_on_free_dofs,mass_type,"
        << "source_channels\n" << std::setprecision(17);
    for (const MatrixDiagnostic* diagnostic : {
            &system.capacityDiagnostic, &system.conductivityDiagnostic}) {
        out << diagnostic->name << ',' << diagnostic->rows << ','
            << diagnostic->nonzeros << ',' << diagnostic->units << ','
            << diagnostic->symmetryError << ',' << diagnostic->minimumDiagonal << ','
            << diagnostic->maximumDiagonal << ',' << diagnostic->zeroDiagonalCount << ','
            << (diagnostic->positiveOnFreeDofs ? 1 : 0) << ','
            << system.massType << ',' << system.sourceChannels << '\n';
    }
}

bool sameFingerprints(const TransientFingerprints& left,
                      const TransientFingerprints& right)
{
    return left.mesh == right.mesh
        && left.capacity == right.capacity
        && left.conductivity == right.conductivity
        && left.input == right.input
        && left.boundary == right.boundary
        && left.sources == right.sources;
}

} // namespace mor::transient
