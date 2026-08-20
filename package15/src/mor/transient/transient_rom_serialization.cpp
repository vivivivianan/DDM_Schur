#include "transient_rom_serialization.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace mor::transient {
namespace {

template <typename T>
void writeVector(const std::filesystem::path& path,
                 const std::vector<T>& values)
{
    static_assert(std::is_trivially_copyable<T>::value,
                  "Transient model vectors must be trivially copyable.");
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("Cannot write transient model: " + path.string());
    const std::uint64_t count = static_cast<std::uint64_t>(values.size());
    out.write(reinterpret_cast<const char*>(&count), sizeof(count));
    if (!values.empty()) {
        out.write(reinterpret_cast<const char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(T)));
    }
    out.flush();
    if (!out) {
        throw std::runtime_error(
            "Transient model write failed or disk is full: " + path.string());
    }
}

template <typename T>
std::vector<T> readVector(const std::filesystem::path& path,
                          std::size_t expected)
{
    static_assert(std::is_trivially_copyable<T>::value,
                  "Transient model vectors must be trivially copyable.");
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot read transient model: " + path.string());
    std::uint64_t count = 0;
    in.read(reinterpret_cast<char*>(&count), sizeof(count));
    if (count != expected) {
        throw std::runtime_error("Transient model array dimension mismatch: " + path.string());
    }
    std::vector<T> values(static_cast<std::size_t>(count));
    if (!values.empty()) {
        in.read(reinterpret_cast<char*>(values.data()),
            static_cast<std::streamsize>(values.size() * sizeof(T)));
    }
    if (!in) throw std::runtime_error("Transient model array is truncated: " + path.string());
    return values;
}

std::string readText(const std::filesystem::path& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read transient metadata: " + path.string());
    return std::string(std::istreambuf_iterator<char>(in),
                       std::istreambuf_iterator<char>());
}

std::string jsonToken(const std::string& text, const std::string& key)
{
    const std::string marker = "\"" + key + "\"";
    const std::size_t keyPosition = text.find(marker);
    if (keyPosition == std::string::npos) {
        throw std::runtime_error("Missing transient metadata key: " + key);
    }
    const std::size_t colon = text.find(':', keyPosition + marker.size());
    std::size_t begin = text.find_first_not_of(" \t\r\n", colon + 1);
    if (text[begin] == '"') {
        const std::size_t end = text.find('"', begin + 1);
        return text.substr(begin + 1, end - begin - 1);
    }
    const std::size_t end = text.find_first_of(",}\r\n", begin);
    return text.substr(begin, end - begin);
}

std::uint64_t jsonUnsigned(const std::string& text, const std::string& key)
{
    const std::string value = jsonToken(text, key);
    return std::stoull(value, nullptr, 10);
}

std::uint64_t jsonHex(const std::string& text, const std::string& key)
{
    return std::stoull(jsonToken(text, key), nullptr, 16);
}

double jsonDouble(const std::string& text, const std::string& key)
{
    return std::stod(jsonToken(text, key));
}

std::string hex(std::uint64_t value)
{
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << value;
    return out.str();
}

} // namespace

void saveTransientReducedModel(
    const TransientReducedModel& model,
    const std::filesystem::path& directory)
{
    std::filesystem::create_directories(directory);
    std::ofstream metadata(directory / "metadata.json");
    metadata << std::setprecision(17)
        << "{\n"
        << "  \"format_version\": " << model.formatVersion << ",\n"
        << "  \"model_type\": \"transient-block-arnoldi\",\n"
        << "  \"global_dofs\": " << model.globalDofs << ",\n"
        << "  \"rank\": " << model.rank << ",\n"
        << "  \"block_size\": " << model.blockSize << ",\n"
        << "  \"moments\": " << model.moments << ",\n"
        << "  \"source_channels\": " << model.sourceChannels << ",\n"
        << "  \"expansion_point\": " << model.expansionPoint << ",\n"
        << "  \"rank_tolerance\": " << model.rankTolerance << ",\n"
        << "  \"mass_type\": \"" << model.massType << "\",\n"
        << "  \"dt_policy\": \"arbitrary positive fixed dt; refactor reduced matrix per dt\",\n"
        << "  \"basis_orthogonality_error\": "
            << model.basisOrthogonalityError << ",\n"
        << "  \"reference_residual\": " << model.referenceResidual << ",\n"
        << "  \"capacity_symmetry_error\": "
            << model.capacityDiagnostic.symmetryError << ",\n"
        << "  \"capacity_min_eigenvalue\": "
            << model.capacityDiagnostic.minimumEigenvalue << ",\n"
        << "  \"capacity_max_eigenvalue\": "
            << model.capacityDiagnostic.maximumEigenvalue << ",\n"
        << "  \"capacity_cholesky\": "
            << (model.capacityDiagnostic.choleskySucceeded ? 1 : 0) << ",\n"
        << "  \"conductivity_symmetry_error\": "
            << model.conductivityDiagnostic.symmetryError << ",\n"
        << "  \"conductivity_min_eigenvalue\": "
            << model.conductivityDiagnostic.minimumEigenvalue << ",\n"
        << "  \"conductivity_max_eigenvalue\": "
            << model.conductivityDiagnostic.maximumEigenvalue << ",\n"
        << "  \"conductivity_cholesky\": "
            << (model.conductivityDiagnostic.choleskySucceeded ? 1 : 0) << "\n"
        << "}\n";

    writeVector(directory / "basis.bin", model.basis);
    writeVector(directory / "reduced_capacity.bin", model.reducedCapacity);
    writeVector(directory / "reduced_conductivity.bin", model.reducedConductivity);
    writeVector(directory / "reduced_input.bin", model.reducedInput);
    writeVector(directory / "reduced_boundary.bin", model.reducedBoundary);
    writeVector(directory / "reference_temperature.bin", model.referenceTemperature);
    writeVector(directory / "nominal_powers.bin", model.nominalPowersW);
    writeVector(directory / "minimum_powers.bin", model.minimumPowersW);
    writeVector(directory / "maximum_powers.bin", model.maximumPowersW);
    writeVector(directory / "source_subdomains.bin", model.sourceSubdomains);
    writeVector(directory / "source_entities.bin", model.sourceDomainEntities);
    writeVector(directory / "dofs.bin", model.deploymentDofs);

    std::ofstream sources(directory / "source_channels.csv");
    sources << "channel,subdomain,domain_entity,nominal_power_w,minimum_power_w,maximum_power_w\n"
        << std::setprecision(17);
    for (int channel = 0; channel < model.sourceChannels; ++channel) {
        sources << channel << ','
            << model.sourceSubdomains[static_cast<std::size_t>(channel)] << ','
            << model.sourceDomainEntities[static_cast<std::size_t>(channel)] << ','
            << model.nominalPowersW[static_cast<std::size_t>(channel)] << ','
            << model.minimumPowersW[static_cast<std::size_t>(channel)] << ','
            << model.maximumPowersW[static_cast<std::size_t>(channel)] << '\n';
    }
    std::ofstream history(directory / "arnoldi_history.csv");
    history << "moment,input_columns,added_rank,cumulative_rank,deflated_columns,"
        << "orthogonality_error,arnoldi_residual,solve_seconds,orthogonalization_seconds,basis_bytes\n"
        << std::setprecision(17);
    for (const ArnoldiHistoryRow& row : model.arnoldiHistory) {
        history << row.moment << ',' << row.inputColumns << ',' << row.addedRank << ','
            << row.cumulativeRank << ',' << row.deflatedColumns << ','
            << row.orthogonalityError << ',' << row.arnoldiResidual << ','
            << row.solveSeconds << ',' << row.orthogonalizationSeconds << ','
            << row.basisBytes << '\n';
    }
    std::ofstream expansion(directory / "expansion_points.csv");
    expansion << "index,expansion_point\n0," << std::setprecision(17)
        << model.expansionPoint << '\n';
    std::ofstream fingerprints(directory / "fingerprints.json");
    fingerprints << "{\n"
        << "  \"mesh\": \"" << hex(model.fingerprints.mesh) << "\",\n"
        << "  \"capacity\": \"" << hex(model.fingerprints.capacity) << "\",\n"
        << "  \"conductivity\": \"" << hex(model.fingerprints.conductivity) << "\",\n"
        << "  \"input\": \"" << hex(model.fingerprints.input) << "\",\n"
        << "  \"boundary\": \"" << hex(model.fingerprints.boundary) << "\",\n"
        << "  \"sources\": \"" << hex(model.fingerprints.sources) << "\"\n"
        << "}\n";
    metadata.flush();
    sources.flush();
    history.flush();
    expansion.flush();
    fingerprints.flush();
    if (!metadata || !sources || !history || !expansion || !fingerprints) {
        throw std::runtime_error(
            "Transient model metadata write failed or disk is full: "
            + directory.string());
    }
}

TransientReducedModel loadTransientReducedModel(
    const std::filesystem::path& directory,
    const TransientFingerprints* expectedFingerprints)
{
    const std::string metadata = readText(directory / "metadata.json");
    const std::string fingerprints = readText(directory / "fingerprints.json");
    TransientReducedModel model;
    model.formatVersion = static_cast<int>(jsonUnsigned(metadata, "format_version"));
    if (model.formatVersion != 1
        || jsonToken(metadata, "model_type") != "transient-block-arnoldi") {
        throw std::runtime_error("Unsupported transient Block Arnoldi model format.");
    }
    model.globalDofs = static_cast<int>(jsonUnsigned(metadata, "global_dofs"));
    model.rank = static_cast<int>(jsonUnsigned(metadata, "rank"));
    model.blockSize = static_cast<int>(jsonUnsigned(metadata, "block_size"));
    model.moments = static_cast<int>(jsonUnsigned(metadata, "moments"));
    model.sourceChannels = static_cast<int>(jsonUnsigned(metadata, "source_channels"));
    model.expansionPoint = jsonDouble(metadata, "expansion_point");
    model.rankTolerance = jsonDouble(metadata, "rank_tolerance");
    model.massType = jsonToken(metadata, "mass_type");
    model.basisOrthogonalityError = jsonDouble(metadata, "basis_orthogonality_error");
    model.referenceResidual = jsonDouble(metadata, "reference_residual");
    model.capacityDiagnostic.symmetryError = jsonDouble(metadata, "capacity_symmetry_error");
    model.capacityDiagnostic.minimumEigenvalue = jsonDouble(metadata, "capacity_min_eigenvalue");
    model.capacityDiagnostic.maximumEigenvalue = jsonDouble(metadata, "capacity_max_eigenvalue");
    model.capacityDiagnostic.choleskySucceeded =
        jsonUnsigned(metadata, "capacity_cholesky") != 0;
    model.conductivityDiagnostic.symmetryError =
        jsonDouble(metadata, "conductivity_symmetry_error");
    model.conductivityDiagnostic.minimumEigenvalue =
        jsonDouble(metadata, "conductivity_min_eigenvalue");
    model.conductivityDiagnostic.maximumEigenvalue =
        jsonDouble(metadata, "conductivity_max_eigenvalue");
    model.conductivityDiagnostic.choleskySucceeded =
        jsonUnsigned(metadata, "conductivity_cholesky") != 0;
    model.fingerprints.mesh = jsonHex(fingerprints, "mesh");
    model.fingerprints.capacity = jsonHex(fingerprints, "capacity");
    model.fingerprints.conductivity = jsonHex(fingerprints, "conductivity");
    model.fingerprints.input = jsonHex(fingerprints, "input");
    model.fingerprints.boundary = jsonHex(fingerprints, "boundary");
    model.fingerprints.sources = jsonHex(fingerprints, "sources");
    if (expectedFingerprints != nullptr
        && !sameFingerprints(model.fingerprints, *expectedFingerprints)) {
        throw std::runtime_error(
            "Transient model fingerprint rejection: mesh/C/K/B/boundary/source definition changed.");
    }
    model.basis = readVector<double>(directory / "basis.bin",
        static_cast<std::size_t>(model.globalDofs) * model.rank);
    model.reducedCapacity = readVector<double>(directory / "reduced_capacity.bin",
        static_cast<std::size_t>(model.rank) * model.rank);
    model.reducedConductivity = readVector<double>(directory / "reduced_conductivity.bin",
        static_cast<std::size_t>(model.rank) * model.rank);
    model.reducedInput = readVector<double>(directory / "reduced_input.bin",
        static_cast<std::size_t>(model.rank) * model.sourceChannels);
    model.reducedBoundary = readVector<double>(directory / "reduced_boundary.bin",
        static_cast<std::size_t>(model.rank));
    model.referenceTemperature = readVector<double>(directory / "reference_temperature.bin",
        static_cast<std::size_t>(model.globalDofs));
    model.nominalPowersW = readVector<double>(directory / "nominal_powers.bin",
        static_cast<std::size_t>(model.sourceChannels));
    model.minimumPowersW = readVector<double>(directory / "minimum_powers.bin",
        static_cast<std::size_t>(model.sourceChannels));
    model.maximumPowersW = readVector<double>(directory / "maximum_powers.bin",
        static_cast<std::size_t>(model.sourceChannels));
    model.sourceSubdomains = readVector<int>(directory / "source_subdomains.bin",
        static_cast<std::size_t>(model.sourceChannels));
    model.sourceDomainEntities = readVector<int>(directory / "source_entities.bin",
        static_cast<std::size_t>(model.sourceChannels));
    model.deploymentDofs = readVector<DeploymentDof>(directory / "dofs.bin",
        static_cast<std::size_t>(model.globalDofs));
    return model;
}

std::uint64_t transientModelFileBytes(const std::filesystem::path& directory)
{
    std::uint64_t bytes = 0;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file()) bytes += entry.file_size();
    }
    return bytes;
}

} // namespace mor::transient
