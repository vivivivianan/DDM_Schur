#include "interior_model_serialization.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace mor {
namespace {

std::string hexValue(std::uint64_t value)
{
    std::ostringstream stream;
    stream << std::hex << std::setw(16) << std::setfill('0') << value;
    return stream.str();
}

std::string readText(const std::filesystem::path& path)
{
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("[ROM Deployment] Cannot read metadata: " + path.string());
    }
    return std::string(std::istreambuf_iterator<char>(in),
                       std::istreambuf_iterator<char>());
}

std::string jsonToken(const std::string& text, const std::string& key)
{
    const std::string marker = "\"" + key + "\"";
    std::size_t position = text.find(marker);
    if (position == std::string::npos) {
        throw std::runtime_error("[ROM Deployment] Missing metadata key: " + key);
    }
    position = text.find(':', position + marker.size());
    if (position == std::string::npos) {
        throw std::runtime_error("[ROM Deployment] Malformed metadata key: " + key);
    }
    ++position;
    while (position < text.size() && std::isspace(static_cast<unsigned char>(text[position]))) {
        ++position;
    }
    if (position < text.size() && text[position] == '"') {
        const std::size_t end = text.find('"', position + 1);
        if (end == std::string::npos) {
            throw std::runtime_error("[ROM Deployment] Malformed string metadata: " + key);
        }
        return text.substr(position + 1, end - position - 1);
    }
    std::size_t end = position;
    while (end < text.size() && text[end] != ',' && text[end] != '\n' && text[end] != '}') {
        ++end;
    }
    return text.substr(position, end - position);
}

std::uint64_t jsonUnsigned(const std::string& text, const std::string& key)
{
    const std::string token = jsonToken(text, key);
    const bool hexadecimal = token.size() == 16
        && token.find_first_not_of("0123456789abcdefABCDEF") == std::string::npos;
    return std::stoull(token, nullptr, hexadecimal ? 16 : 10);
}

void writeIntegers(const std::filesystem::path& path, const std::vector<int>& values)
{
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("[ROM Deployment] Cannot write integer array: " + path.string());
    }
    for (int value : values) {
        const std::int32_t fixed = static_cast<std::int32_t>(value);
        out.write(reinterpret_cast<const char*>(&fixed), sizeof(fixed));
    }
}

std::vector<int> readIntegers(const std::filesystem::path& path, std::size_t count)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("[ROM Deployment] Cannot read integer array: " + path.string());
    }
    std::vector<int> result(count, 0);
    for (std::size_t i = 0; i < count; ++i) {
        std::int32_t fixed = 0;
        in.read(reinterpret_cast<char*>(&fixed), sizeof(fixed));
        if (!in) {
            throw std::runtime_error("[ROM Deployment] Truncated integer array: " + path.string());
        }
        result[i] = static_cast<int>(fixed);
    }
    return result;
}

void writeScalars(const std::filesystem::path& path,
                  const std::vector<double>& values,
                  const std::string& precision)
{
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("[ROM Deployment] Cannot write response array: " + path.string());
    }
    if (precision == "float64") {
        out.write(reinterpret_cast<const char*>(values.data()),
                  static_cast<std::streamsize>(values.size() * sizeof(double)));
    } else if (precision == "float32") {
        for (double value : values) {
            const float converted = static_cast<float>(value);
            out.write(reinterpret_cast<const char*>(&converted), sizeof(converted));
        }
    } else {
        throw std::runtime_error("--mor-storage-precision must be float64 or float32.");
    }
}

std::vector<double> readScalars(const std::filesystem::path& path,
                                std::size_t count,
                                const std::string& precision)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("[ROM Deployment] Cannot read response array: " + path.string());
    }
    std::vector<double> result(count, 0.0);
    if (precision == "float64") {
        in.read(reinterpret_cast<char*>(result.data()),
                static_cast<std::streamsize>(count * sizeof(double)));
    } else if (precision == "float32") {
        for (double& value : result) {
            float converted = 0.0F;
            in.read(reinterpret_cast<char*>(&converted), sizeof(converted));
            value = static_cast<double>(converted);
        }
    } else {
        throw std::runtime_error("[ROM Deployment] Unsupported storage precision.");
    }
    if (!in) {
        throw std::runtime_error("[ROM Deployment] Truncated response array: " + path.string());
    }
    return result;
}

std::filesystem::path domainPath(const std::filesystem::path& root, int domain)
{
    std::ostringstream name;
    name << "subdomain_" << std::setw(3) << std::setfill('0') << domain;
    return root / "interior" / name.str();
}

std::uint64_t hashDofs(const std::vector<int>& dofs)
{
    std::uint64_t value = 1469598103934665603ULL;
    auto addBytes = [&](const void* data, std::size_t count) {
        const auto* bytes = static_cast<const unsigned char*>(data);
        for (std::size_t i = 0; i < count; ++i) {
            value ^= static_cast<std::uint64_t>(bytes[i]);
            value *= 1099511628211ULL;
        }
    };
    const std::uint64_t count = static_cast<std::uint64_t>(dofs.size());
    addBytes(&count, sizeof(count));
    for (int dof : dofs) {
        addBytes(&dof, sizeof(dof));
    }
    return value;
}

std::uint64_t hashGlobalDofs(const std::vector<DeploymentDof>& dofs)
{
    std::uint64_t value = 1469598103934665603ULL;
    auto addBytes = [&](const void* data, std::size_t count) {
        const auto* bytes = static_cast<const unsigned char*>(data);
        for (std::size_t i = 0; i < count; ++i) {
            value ^= static_cast<std::uint64_t>(bytes[i]);
            value *= 1099511628211ULL;
        }
    };
    const std::uint64_t count = static_cast<std::uint64_t>(dofs.size());
    addBytes(&count, sizeof(count));
    for (std::size_t index = 0; index < dofs.size(); ++index) {
        const DeploymentDof& dof = dofs[index];
        addBytes(&index, sizeof(index));
        addBytes(&dof.subdomain, sizeof(dof.subdomain));
        addBytes(&dof.sourceVertex, sizeof(dof.sourceVertex));
        addBytes(&dof.x, sizeof(dof.x));
        addBytes(&dof.y, sizeof(dof.y));
        addBytes(&dof.z, sizeof(dof.z));
    }
    return value;
}

void writeMetadata(const DeploymentResponseModel& model,
                   const std::filesystem::path& deployment,
                   std::uint64_t deploymentBytes,
                   std::uint64_t completeBytes)
{
    std::ofstream out(deployment / "full_response_metadata.json");
    if (!out) {
        throw std::runtime_error("[ROM Deployment] Cannot write deployment metadata.");
    }
    out << "{\n"
        << "  \"format_version\": 2,\n"
        << "  \"model_type\": \"fixed-matrix-input-to-temperature-response\",\n"
        << "  \"interior_mode\": \"" << model.interiorMode << "\",\n"
        << "  \"storage_precision\": \"" << model.storagePrecision << "\",\n"
        << "  \"layout\": \"column-major\",\n"
        << "  \"global_dofs\": " << model.globalDofs << ",\n"
        << "  \"interface_dofs\": " << model.interfaceDofs << ",\n"
        << "  \"interface_rank\": " << model.interfaceRank << ",\n"
        << "  \"source_channels\": " << model.sourceChannels << ",\n"
        << "  \"subdomains\": " << model.interiors.size() << ",\n"
        << "  \"mesh_fingerprint\": \"" << hexValue(model.fingerprints.mesh) << "\",\n"
        << "  \"system_fingerprint\": \"" << hexValue(model.fingerprints.system) << "\",\n"
        << "  \"interface_fingerprint\": \""
        << hexValue(model.fingerprints.interfaceOrdering) << "\",\n"
        << "  \"source_fingerprint\": \"" << hexValue(model.fingerprints.sources) << "\",\n"
        << "  \"global_dof_ordering_hash\": \""
        << hexValue(model.globalDofOrderingHash) << "\",\n"
        << "  \"deployment_file_bytes\": " << deploymentBytes << ",\n"
        << "  \"complete_model_file_bytes\": " << completeBytes << "\n"
        << "}\n";
}

} // namespace

std::uint64_t recursiveFileBytes(const std::filesystem::path& directory)
{
    std::uint64_t result = 0;
    if (!std::filesystem::exists(directory)) {
        return 0;
    }
    for (const auto& entry : std::filesystem::recursive_directory_iterator(directory)) {
        if (entry.is_regular_file()) {
            result += static_cast<std::uint64_t>(entry.file_size());
        }
    }
    return result;
}

void saveDeploymentResponseModel(DeploymentResponseModel& model,
                                 const std::filesystem::path& modelDirectory,
                                 const std::string& interiorMode,
                                 const std::string& storagePrecision)
{
    const auto start = std::chrono::steady_clock::now();
    if (interiorMode != "exact-response" && interiorMode != "compressed-rb") {
        throw std::runtime_error("[MOR Interior] Cannot serialize unsupported interior mode.");
    }
    if (storagePrecision != "float64" && storagePrecision != "float32") {
        throw std::runtime_error("--mor-storage-precision must be float64 or float32.");
    }
    model.interiorMode = interiorMode;
    model.storagePrecision = storagePrecision;
    const std::filesystem::path deployment = modelDirectory / "deployment";
    std::error_code removeError;
    std::filesystem::remove_all(deployment, removeError);
    if (removeError) {
        throw std::runtime_error("[ROM Deployment] Cannot replace prior deployment directory.");
    }
    std::filesystem::create_directories(deployment / "interior");

    std::vector<double> coordinates;
    std::vector<int> subdomains;
    std::vector<int> sourceVertices;
    coordinates.reserve(model.dofs.size() * 3);
    subdomains.reserve(model.dofs.size());
    sourceVertices.reserve(model.dofs.size());
    for (const DeploymentDof& dof : model.dofs) {
        coordinates.push_back(dof.x);
        coordinates.push_back(dof.y);
        coordinates.push_back(dof.z);
        subdomains.push_back(dof.subdomain);
        sourceVertices.push_back(dof.sourceVertex);
    }
    writeScalars(deployment / "global_dof_coordinates.bin", coordinates, "float64");
    writeIntegers(deployment / "global_dof_subdomains.bin", subdomains);
    writeIntegers(deployment / "global_dof_source_vertices.bin", sourceVertices);
    writeIntegers(deployment / "interface_global_dofs.bin", model.interfaceGlobalDofs);
    writeScalars(deployment / "interface_reference.bin", model.interfaceReference,
                 storagePrecision);
    writeScalars(deployment / "interface_power_response.bin",
                 model.interfacePowerResponse, storagePrecision);
    writeScalars(deployment / "reduced_input_operator.bin",
                 model.reducedInputOperator, storagePrecision);
    writeScalars(deployment / "nominal_powers.bin", model.nominalPowers, "float64");
    writeScalars(deployment / "minimum_powers.bin", model.minimumPowers, "float64");
    writeScalars(deployment / "maximum_powers.bin", model.maximumPowers, "float64");
    {
        std::ofstream channels(deployment / "power_channels.csv");
        channels << "channel,subdomain,domain_entity,nominal_power_w,minimum_power_w,maximum_power_w\n"
                 << std::setprecision(17);
        for (int channel = 0; channel < model.sourceChannels; ++channel) {
            const int subdomain = static_cast<std::size_t>(channel) < model.sourceSubdomains.size()
                ? model.sourceSubdomains[static_cast<std::size_t>(channel)] : -1;
            const int entity = static_cast<std::size_t>(channel) < model.sourceDomainEntities.size()
                ? model.sourceDomainEntities[static_cast<std::size_t>(channel)] : -1;
            channels << channel << ',' << subdomain << ',' << entity << ','
                     << model.nominalPowers[static_cast<std::size_t>(channel)] << ','
                     << model.minimumPowers[static_cast<std::size_t>(channel)] << ','
                     << model.maximumPowers[static_cast<std::size_t>(channel)] << '\n';
        }
    }

    for (const InteriorResponseBlock& block : model.interiors) {
        const std::filesystem::path path = domainPath(deployment, block.subdomain);
        std::filesystem::create_directories(path);
        writeIntegers(path / "global_dofs.bin", block.globalDofs);
        writeScalars(path / "reference_temperature.bin",
                     block.referenceTemperature, storagePrecision);
        std::ofstream channels(path / "power_channel_map.csv");
        channels << "local_channel,global_power_channel\n";
        for (std::size_t local = 0; local < block.directPowerChannels.size(); ++local) {
            channels << local << ',' << block.directPowerChannels[local] << '\n';
        }
        if (interiorMode == "exact-response") {
            std::filesystem::remove(path / "local_basis.bin");
            std::filesystem::remove(path / "local_coordinate_map.bin");
            std::filesystem::remove(path / "singular_values.csv");
            writeScalars(path / "exact_response.bin", block.exactResponse,
                         storagePrecision);
        } else {
            if (block.rank <= 0) {
                throw std::runtime_error("[MOR Interior] Compressed block has zero rank.");
            }
            std::filesystem::remove(path / "exact_response.bin");
            writeScalars(path / "local_basis.bin", block.localBasis,
                         storagePrecision);
            writeScalars(path / "local_coordinate_map.bin",
                         block.localCoordinateMap, storagePrecision);
            std::ofstream singular(path / "singular_values.csv");
            singular << "index,singular_value,retained_energy\n" << std::setprecision(17);
            for (std::size_t i = 0; i < block.singularValues.size(); ++i) {
                const double energy = i < block.retainedEnergy.size()
                    ? block.retainedEnergy[i] : 1.0;
                singular << i << ',' << block.singularValues[i] << ',' << energy << '\n';
            }
        }
        std::ofstream metadata(path / "metadata.json");
        metadata << "{\n"
            << "  \"subdomain\": " << block.subdomain << ",\n"
            << "  \"interior_dofs\": " << block.globalDofs.size() << ",\n"
            << "  \"power_channels\": " << model.sourceChannels << ",\n"
            << "  \"direct_power_channels\": " << block.directPowerChannels.size() << ",\n"
            << "  \"rank\": " << (interiorMode == "compressed-rb" ? block.rank : 0) << ",\n"
            << "  \"dof_ordering_hash\": \"" << hexValue(block.dofOrderingHash) << "\"\n"
            << "}\n";
    }
    writeMetadata(model, deployment, 0, 0);
    for (int pass = 0; pass < 3; ++pass) {
        model.deploymentFileBytes = recursiveFileBytes(deployment);
        model.completeModelFileBytes = recursiveFileBytes(modelDirectory);
        writeMetadata(model, deployment, model.deploymentFileBytes,
                      model.completeModelFileBytes);
    }
    model.deploymentFileBytes = recursiveFileBytes(deployment);
    model.completeModelFileBytes = recursiveFileBytes(modelDirectory);
    model.serializationSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
}

DeploymentResponseModel loadDeploymentResponseModel(
    const std::filesystem::path& modelDirectory,
    const std::string& requestedInteriorMode,
    const Fingerprints* expectedFingerprints,
    int expectedGlobalDofs,
    int expectedInterfaceDofs,
    int expectedInterfaceRank,
    int expectedSourceChannels)
{
    const auto start = std::chrono::steady_clock::now();
    const std::filesystem::path deployment = modelDirectory / "deployment";
    const std::string metadata = readText(deployment / "full_response_metadata.json");
    DeploymentResponseModel model;
    model.formatVersion = static_cast<int>(jsonUnsigned(metadata, "format_version"));
    if (model.formatVersion != 2) {
        throw std::runtime_error("[ROM Deployment] Unsupported deployment model format.");
    }
    model.interiorMode = jsonToken(metadata, "interior_mode");
    model.storagePrecision = jsonToken(metadata, "storage_precision");
    model.globalDofs = static_cast<int>(jsonUnsigned(metadata, "global_dofs"));
    model.interfaceDofs = static_cast<int>(jsonUnsigned(metadata, "interface_dofs"));
    model.interfaceRank = static_cast<int>(jsonUnsigned(metadata, "interface_rank"));
    model.sourceChannels = static_cast<int>(jsonUnsigned(metadata, "source_channels"));
    const int subdomains = static_cast<int>(jsonUnsigned(metadata, "subdomains"));
    model.fingerprints.mesh = jsonUnsigned(metadata, "mesh_fingerprint");
    model.fingerprints.system = jsonUnsigned(metadata, "system_fingerprint");
    model.fingerprints.interfaceOrdering = jsonUnsigned(metadata, "interface_fingerprint");
    model.fingerprints.sources = jsonUnsigned(metadata, "source_fingerprint");
    model.globalDofOrderingHash = jsonUnsigned(metadata, "global_dof_ordering_hash");
    if (!requestedInteriorMode.empty() && requestedInteriorMode != "pardiso"
        && requestedInteriorMode != model.interiorMode) {
        throw std::runtime_error("[ROM Deployment] Requested interior representation is not stored in this model.");
    }
    if (expectedFingerprints
        && (model.fingerprints.mesh != expectedFingerprints->mesh
            || model.fingerprints.system != expectedFingerprints->system
            || model.fingerprints.interfaceOrdering != expectedFingerprints->interfaceOrdering
            || model.fingerprints.sources != expectedFingerprints->sources)) {
        throw std::runtime_error("[ROM Deployment] Mesh/matrix/interface/source fingerprint mismatch.");
    }
    if ((expectedGlobalDofs > 0 && expectedGlobalDofs != model.globalDofs)
        || (expectedInterfaceDofs > 0 && expectedInterfaceDofs != model.interfaceDofs)
        || (expectedInterfaceRank > 0 && expectedInterfaceRank != model.interfaceRank)
        || (expectedSourceChannels > 0 && expectedSourceChannels != model.sourceChannels)) {
        throw std::runtime_error("[ROM Deployment] Model dimension or Stage 2A rank mismatch.");
    }

    const std::vector<double> coordinates = readScalars(
        deployment / "global_dof_coordinates.bin",
        static_cast<std::size_t>(3 * model.globalDofs), "float64");
    const std::vector<int> domainIds = readIntegers(
        deployment / "global_dof_subdomains.bin",
        static_cast<std::size_t>(model.globalDofs));
    const std::vector<int> sourceVertices = readIntegers(
        deployment / "global_dof_source_vertices.bin",
        static_cast<std::size_t>(model.globalDofs));
    model.dofs.resize(static_cast<std::size_t>(model.globalDofs));
    for (int dof = 0; dof < model.globalDofs; ++dof) {
        DeploymentDof& descriptor = model.dofs[static_cast<std::size_t>(dof)];
        descriptor.x = coordinates[static_cast<std::size_t>(3 * dof)];
        descriptor.y = coordinates[static_cast<std::size_t>(3 * dof + 1)];
        descriptor.z = coordinates[static_cast<std::size_t>(3 * dof + 2)];
        descriptor.subdomain = domainIds[static_cast<std::size_t>(dof)];
        descriptor.sourceVertex = sourceVertices[static_cast<std::size_t>(dof)];
    }
    if (hashGlobalDofs(model.dofs) != model.globalDofOrderingHash) {
        throw std::runtime_error("[ROM Deployment] Global DOF mapping hash mismatch.");
    }
    model.interfaceGlobalDofs = readIntegers(
        deployment / "interface_global_dofs.bin",
        static_cast<std::size_t>(model.interfaceDofs));
    model.interfaceReference = readScalars(deployment / "interface_reference.bin",
        static_cast<std::size_t>(model.interfaceDofs), model.storagePrecision);
    model.interfacePowerResponse = readScalars(
        deployment / "interface_power_response.bin",
        static_cast<std::size_t>(model.interfaceDofs * model.sourceChannels),
        model.storagePrecision);
    model.reducedInputOperator = readScalars(
        deployment / "reduced_input_operator.bin",
        static_cast<std::size_t>(model.interfaceRank * model.sourceChannels),
        model.storagePrecision);
    model.nominalPowers = readScalars(deployment / "nominal_powers.bin",
        static_cast<std::size_t>(model.sourceChannels), "float64");
    model.minimumPowers = readScalars(deployment / "minimum_powers.bin",
        static_cast<std::size_t>(model.sourceChannels), "float64");
    model.maximumPowers = readScalars(deployment / "maximum_powers.bin",
        static_cast<std::size_t>(model.sourceChannels), "float64");
    {
        std::filesystem::path channelPath = deployment / "power_channels.csv";
        if (!std::filesystem::exists(channelPath)) {
            channelPath = modelDirectory / "source_channels.csv";
        }
        model.sourceSubdomains.assign(static_cast<std::size_t>(model.sourceChannels), -1);
        model.sourceDomainEntities.assign(static_cast<std::size_t>(model.sourceChannels), -1);
        if (std::filesystem::exists(channelPath)) {
            std::ifstream channels(channelPath);
            std::string line;
            std::getline(channels, line);
            int rows = 0;
            while (std::getline(channels, line)) {
                std::replace(line.begin(), line.end(), ',', ' ');
                std::istringstream row(line);
                int channel = -1;
                int subdomain = -1;
                int entity = -1;
                double nominal = 0.0;
                double minimum = 0.0;
                double maximum = 0.0;
                if (!(row >> channel >> subdomain >> entity >> nominal >> minimum >> maximum)
                    || channel != rows || channel >= model.sourceChannels) {
                    throw std::runtime_error("[ROM Deployment] Invalid power-channel definition.");
                }
                const auto close = [](double left, double right) {
                    return std::abs(left - right) <= 32.0 * std::numeric_limits<double>::epsilon()
                        * std::max({1.0, std::abs(left), std::abs(right)});
                };
                if (!close(nominal, model.nominalPowers[static_cast<std::size_t>(channel)])
                    || !close(minimum, model.minimumPowers[static_cast<std::size_t>(channel)])
                    || !close(maximum, model.maximumPowers[static_cast<std::size_t>(channel)])) {
                    throw std::runtime_error("[ROM Deployment] Power-channel limits mismatch.");
                }
                model.sourceSubdomains[static_cast<std::size_t>(channel)] = subdomain;
                model.sourceDomainEntities[static_cast<std::size_t>(channel)] = entity;
                ++rows;
            }
            if (rows != model.sourceChannels) {
                throw std::runtime_error("[ROM Deployment] Power-channel count mismatch.");
            }
        }
    }

    std::vector<char> covered(static_cast<std::size_t>(model.globalDofs), 0);
    for (int dof : model.interfaceGlobalDofs) {
        if (dof < 0 || dof >= model.globalDofs
            || covered[static_cast<std::size_t>(dof)]) {
            throw std::runtime_error("[ROM Deployment] Invalid interface DOF map.");
        }
        covered[static_cast<std::size_t>(dof)] = 1;
    }
    model.interiors.resize(static_cast<std::size_t>(subdomains));
    for (int domain = 0; domain < subdomains; ++domain) {
        const std::filesystem::path path = domainPath(deployment, domain);
        const std::string blockMetadata = readText(path / "metadata.json");
        InteriorResponseBlock& block = model.interiors[static_cast<std::size_t>(domain)];
        block.subdomain = static_cast<int>(jsonUnsigned(blockMetadata, "subdomain"));
        const int rows = static_cast<int>(jsonUnsigned(blockMetadata, "interior_dofs"));
        block.rank = static_cast<int>(jsonUnsigned(blockMetadata, "rank"));
        block.dofOrderingHash = jsonUnsigned(blockMetadata, "dof_ordering_hash");
        block.globalDofs = readIntegers(path / "global_dofs.bin",
                                        static_cast<std::size_t>(rows));
        {
            std::ifstream channels(path / "power_channel_map.csv");
            std::string line;
            std::getline(channels, line);
            while (std::getline(channels, line)) {
                const std::size_t comma = line.find(',');
                if (comma == std::string::npos) {
                    continue;
                }
                block.directPowerChannels.push_back(std::stoi(line.substr(comma + 1)));
            }
        }
        if (block.subdomain != domain || hashDofs(block.globalDofs) != block.dofOrderingHash) {
            throw std::runtime_error("[ROM Deployment] Interior DOF ordering hash mismatch.");
        }
        for (int dof : block.globalDofs) {
            if (dof < 0 || dof >= model.globalDofs
                || covered[static_cast<std::size_t>(dof)]
                || model.dofs[static_cast<std::size_t>(dof)].subdomain != domain) {
                throw std::runtime_error("[ROM Deployment] Invalid interior DOF map.");
            }
            covered[static_cast<std::size_t>(dof)] = 1;
        }
        block.referenceTemperature = readScalars(path / "reference_temperature.bin",
            static_cast<std::size_t>(rows), model.storagePrecision);
        if (model.interiorMode == "exact-response") {
            block.exactResponse = readScalars(path / "exact_response.bin",
                static_cast<std::size_t>(rows * model.sourceChannels),
                model.storagePrecision);
        } else {
            if (block.rank <= 0) {
                throw std::runtime_error("[ROM Deployment] Stored compressed rank is invalid.");
            }
            block.localBasis = readScalars(path / "local_basis.bin",
                static_cast<std::size_t>(rows * block.rank), model.storagePrecision);
            block.localCoordinateMap = readScalars(path / "local_coordinate_map.bin",
                static_cast<std::size_t>(block.rank * model.sourceChannels),
                model.storagePrecision);
            std::ifstream singular(path / "singular_values.csv");
            std::string line;
            std::getline(singular, line);
            while (std::getline(singular, line)) {
                std::replace(line.begin(), line.end(), ',', ' ');
                std::istringstream row(line);
                int index = 0;
                double sigma = 0.0;
                double retained = 0.0;
                if (row >> index >> sigma >> retained) {
                    block.singularValues.push_back(sigma);
                    block.retainedEnergy.push_back(retained);
                }
            }
        }
    }
    if (std::find(covered.begin(), covered.end(), 0) != covered.end()) {
        throw std::runtime_error("[ROM Deployment] Deployment DOF maps do not cover the full solution.");
    }
    model.deploymentFileBytes = recursiveFileBytes(deployment);
    model.completeModelFileBytes = recursiveFileBytes(modelDirectory);
    model.loadSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    return model;
}

} // namespace mor
