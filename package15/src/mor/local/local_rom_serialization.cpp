#include "local_rom_serialization.hpp"

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <type_traits>

namespace mor::local {
namespace {

constexpr std::uint64_t magic = 0x4c4f43414c524f4dULL;

template <typename T>
void writeScalar(std::ofstream& out, const T& value)
{
    static_assert(std::is_trivially_copyable_v<T>);
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template <typename T>
T readScalar(std::ifstream& in)
{
    static_assert(std::is_trivially_copyable_v<T>);
    T value{};
    in.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!in) throw std::runtime_error("[Local ROM] Truncated model file.");
    return value;
}

template <typename T>
void writeVector(std::ofstream& out, const std::vector<T>& values)
{
    const std::uint64_t size = static_cast<std::uint64_t>(values.size());
    writeScalar(out, size);
    if (size > 0) {
        out.write(reinterpret_cast<const char*>(values.data()),
            static_cast<std::streamsize>(size * sizeof(T)));
    }
}

template <typename T>
std::vector<T> readVector(std::ifstream& in)
{
    const std::uint64_t size = readScalar<std::uint64_t>(in);
    if (size > (std::uint64_t{1} << 38)) {
        throw std::runtime_error("[Local ROM] Invalid serialized vector size.");
    }
    std::vector<T> values(static_cast<std::size_t>(size));
    if (size > 0) {
        in.read(reinterpret_cast<char*>(values.data()),
            static_cast<std::streamsize>(size * sizeof(T)));
        if (!in) throw std::runtime_error("[Local ROM] Truncated model vector.");
    }
    return values;
}

void writeFingerprints(std::ofstream& out, const Fingerprints& fingerprints)
{
    writeScalar(out, fingerprints.mesh);
    writeScalar(out, fingerprints.system);
    writeScalar(out, fingerprints.interfaceOrdering);
    writeScalar(out, fingerprints.sources);
}

Fingerprints readFingerprints(std::ifstream& in)
{
    Fingerprints result;
    result.mesh = readScalar<std::uint64_t>(in);
    result.system = readScalar<std::uint64_t>(in);
    result.interfaceOrdering = readScalar<std::uint64_t>(in);
    result.sources = readScalar<std::uint64_t>(in);
    return result;
}

bool equalFingerprints(const Fingerprints& left, const Fingerprints& right)
{
    return left.mesh == right.mesh && left.system == right.system
        && left.interfaceOrdering == right.interfaceOrdering
        && left.sources == right.sources;
}

} // namespace

void saveLocalRomModel(const Model& model, const std::filesystem::path& directory)
{
    std::filesystem::create_directories(directory);
    std::ofstream out(directory / "local_rom_model.bin", std::ios::binary);
    if (!out) throw std::runtime_error("[Local ROM] Cannot create model file.");
    writeScalar(out, magic);
    writeScalar(out, model.formatVersion);
    writeScalar(out, model.globalDofs);
    writeScalar(out, model.interfaceDofs);
    writeScalar(out, model.totalLocalRank);
    writeFingerprints(out, model.fingerprints);
    writeVector(out, model.interfaceGlobalDofs);
    writeVector(out, model.interfaceEntries);
    writeScalar(out, static_cast<std::uint64_t>(model.subdomains.size()));
    for (const SubdomainModel& local : model.subdomains) {
        writeScalar(out, local.subdomain);
        writeScalar(out, local.interiorDofs);
        writeScalar(out, local.localInterfaceDofs);
        writeScalar(out, local.snapshots);
        writeScalar(out, local.numericalRank);
        writeScalar(out, local.rank);
        writeScalar(out, local.meshFingerprint);
        writeScalar(out, local.materialFingerprint);
        writeScalar(out, local.boundaryInterfaceFingerprint);
        writeScalar(out, local.templateId);
        writeScalar(out, local.templateReused);
        writeScalar(out, local.templateConsistencyDifference);
        writeScalar(out, local.orthogonalityError);
        writeScalar(out, local.snapshotExtractionSeconds);
        writeScalar(out, local.podSeconds);
        writeScalar(out, local.projectionSeconds);
        writeScalar(out, local.reducedInteriorSymmetryError);
        writeScalar(out, local.reducedInteriorMinimumEigenvalue);
        writeScalar(out, local.reducedInteriorMaximumEigenvalue);
        writeScalar(out, local.reducedInteriorConditionEstimate);
        writeScalar(out, local.couplingSymmetryError);
        writeScalar(out, local.localSchurSymmetryError);
        writeVector(out, local.interiorGlobalDofs);
        writeVector(out, local.interfaceGlobalDofs);
        writeVector(out, local.interfaceIndices);
        writeVector(out, local.referenceInterior);
        writeVector(out, local.basis);
        writeVector(out, local.singularValues);
        writeVector(out, local.retainedEnergy);
        writeVector(out, local.reducedInterior);
        writeVector(out, local.reducedInteriorInterface);
        writeVector(out, local.reducedInterfaceInterior);
        writeVector(out, local.interiorReferenceImage);
        writeVector(out, local.interfaceReferenceImage);
    }
    if (!out) throw std::runtime_error("[Local ROM] Failed while writing model data.");

    std::ofstream metadata(directory / "local_rom_model.json");
    metadata << "{\n"
        << "  \"method\": \"Local-POD-Schur-ROM\",\n"
        << "  \"basis_scope\": \"independent_per_subdomain\",\n"
        << "  \"interface_scope\": \"full_physical_SIPG_interface\",\n"
        << "  \"format_version\": " << model.formatVersion << ",\n"
        << "  \"global_dofs\": " << model.globalDofs << ",\n"
        << "  \"interface_dofs\": " << model.interfaceDofs << ",\n"
        << "  \"subdomains\": " << model.subdomains.size() << ",\n"
        << "  \"total_local_rank\": " << model.totalLocalRank << "\n"
        << "}\n";
}

Model loadLocalRomModel(const std::filesystem::path& directory,
                        const Fingerprints& expectedFingerprints,
                        int expectedGlobalDofs,
                        int expectedInterfaceDofs)
{
    std::ifstream in(directory / "local_rom_model.bin", std::ios::binary);
    if (!in) throw std::runtime_error("[Local ROM] Cannot open serialized model.");
    if (readScalar<std::uint64_t>(in) != magic) {
        throw std::runtime_error("[Local ROM] Invalid model magic.");
    }
    Model model;
    model.formatVersion = readScalar<int>(in);
    model.globalDofs = readScalar<int>(in);
    model.interfaceDofs = readScalar<int>(in);
    model.totalLocalRank = readScalar<int>(in);
    model.fingerprints = readFingerprints(in);
    if (model.formatVersion != 6 || model.globalDofs != expectedGlobalDofs
        || model.interfaceDofs != expectedInterfaceDofs
        || !equalFingerprints(model.fingerprints, expectedFingerprints)) {
        throw std::runtime_error(
            "[Local ROM] Model fingerprint mismatch: mesh/matrix/interface/source ordering changed.");
    }
    model.interfaceGlobalDofs = readVector<int>(in);
    model.interfaceEntries = readVector<InterfaceEntry>(in);
    const std::uint64_t subdomains = readScalar<std::uint64_t>(in);
    model.subdomains.resize(static_cast<std::size_t>(subdomains));
    for (SubdomainModel& local : model.subdomains) {
        local.subdomain = readScalar<int>(in);
        local.interiorDofs = readScalar<int>(in);
        local.localInterfaceDofs = readScalar<int>(in);
        local.snapshots = readScalar<int>(in);
        local.numericalRank = readScalar<int>(in);
        local.rank = readScalar<int>(in);
        local.meshFingerprint = readScalar<std::uint64_t>(in);
        local.materialFingerprint = readScalar<std::uint64_t>(in);
        local.boundaryInterfaceFingerprint = readScalar<std::uint64_t>(in);
        local.templateId = readScalar<int>(in);
        local.templateReused = readScalar<bool>(in);
        local.templateConsistencyDifference = readScalar<double>(in);
        local.orthogonalityError = readScalar<double>(in);
        local.snapshotExtractionSeconds = readScalar<double>(in);
        local.podSeconds = readScalar<double>(in);
        local.projectionSeconds = readScalar<double>(in);
        local.reducedInteriorSymmetryError = readScalar<double>(in);
        local.reducedInteriorMinimumEigenvalue = readScalar<double>(in);
        local.reducedInteriorMaximumEigenvalue = readScalar<double>(in);
        local.reducedInteriorConditionEstimate = readScalar<double>(in);
        local.couplingSymmetryError = readScalar<double>(in);
        local.localSchurSymmetryError = readScalar<double>(in);
        local.interiorGlobalDofs = readVector<int>(in);
        local.interfaceGlobalDofs = readVector<int>(in);
        local.interfaceIndices = readVector<int>(in);
        local.referenceInterior = readVector<double>(in);
        local.basis = readVector<double>(in);
        local.singularValues = readVector<double>(in);
        local.retainedEnergy = readVector<double>(in);
        local.reducedInterior = readVector<double>(in);
        local.reducedInteriorInterface = readVector<double>(in);
        local.reducedInterfaceInterior = readVector<double>(in);
        local.interiorReferenceImage = readVector<double>(in);
        local.interfaceReferenceImage = readVector<double>(in);
    }
    model.modelBytes = estimateModelBytes(model);
    return model;
}

} // namespace mor::local
