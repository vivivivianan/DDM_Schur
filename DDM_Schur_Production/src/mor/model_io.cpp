// Versioned binary I/O for local-ROM caches. Readers use bounded counts and
// fingerprints before allocation; writers use temporary files followed by an
// atomic rename so interrupted cold starts cannot masquerade as valid caches.

#include "model_io.hpp"

#include "sipg_core.hpp"

#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace mor {
namespace {

class Fnv1a {
public:
    template <typename T>
    void add(const T& value)
    {
        const unsigned char* bytes = reinterpret_cast<const unsigned char*>(&value);
        for (std::size_t i = 0; i < sizeof(T); ++i) {
            value_ ^= static_cast<std::uint64_t>(bytes[i]);
            value_ *= 1099511628211ULL;
        }
    }

    void addString(const std::string& value)
    {
        for (unsigned char byte : value) {
            value_ ^= static_cast<std::uint64_t>(byte);
            value_ *= 1099511628211ULL;
        }
        const std::uint64_t length = static_cast<std::uint64_t>(value.size());
        add(length);
    }

    std::uint64_t value() const { return value_; }

private:
    std::uint64_t value_ = 1469598103934665603ULL;
};

void writeVector(const std::filesystem::path& path,
                 const std::vector<double>& values)
{
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("[MOR] Cannot write model array: " + path.string());
    }
    const std::uint64_t count = static_cast<std::uint64_t>(values.size());
    out.write(reinterpret_cast<const char*>(&count), sizeof(count));
    if (!values.empty()) {
        out.write(reinterpret_cast<const char*>(values.data()),
                  static_cast<std::streamsize>(values.size() * sizeof(double)));
    }
}

std::vector<double> readVector(const std::filesystem::path& path,
                               std::size_t expectedCount)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("[MOR] Cannot read model array: " + path.string());
    }
    std::uint64_t count = 0;
    in.read(reinterpret_cast<char*>(&count), sizeof(count));
    if (count != expectedCount) {
        throw std::runtime_error("[MOR] Model array size mismatch: " + path.string());
    }
    std::vector<double> values(static_cast<std::size_t>(count));
    if (count > 0) {
        in.read(reinterpret_cast<char*>(values.data()),
                static_cast<std::streamsize>(values.size() * sizeof(double)));
    }
    if (!in) {
        throw std::runtime_error("[MOR] Model array is truncated: " + path.string());
    }
    return values;
}

std::string readText(const std::filesystem::path& path)
{
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("[MOR] Cannot read model metadata: " + path.string());
    }
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

std::uint64_t jsonUnsigned(const std::string& text, const std::string& key)
{
    const std::string marker = "\"" + key + "\"";
    const std::size_t keyPosition = text.find(marker);
    if (keyPosition == std::string::npos) {
        throw std::runtime_error("[MOR] Missing metadata key: " + key);
    }
    const std::size_t colon = text.find(':', keyPosition + marker.size());
    const std::size_t quote = text.find('"', colon + 1);
    if (quote != std::string::npos && quote < text.find_first_of(",}\n", colon + 1)) {
        const std::size_t endQuote = text.find('"', quote + 1);
        return std::stoull(text.substr(quote + 1, endQuote - quote - 1), nullptr, 16);
    }
    const std::size_t begin = text.find_first_of("0123456789", colon + 1);
    const std::size_t end = text.find_first_not_of("0123456789", begin);
    return std::stoull(text.substr(begin, end - begin));
}

double jsonDouble(const std::string& text, const std::string& key)
{
    const std::string marker = "\"" + key + "\"";
    const std::size_t keyPosition = text.find(marker);
    if (keyPosition == std::string::npos) {
        throw std::runtime_error("[MOR] Missing metadata key: " + key);
    }
    const std::size_t colon = text.find(':', keyPosition + marker.size());
    const std::size_t begin = text.find_first_of("-+0123456789.", colon + 1);
    const std::size_t end = text.find_first_not_of("-+0123456789.eE", begin);
    return std::stod(text.substr(begin, end - begin));
}

std::string hexValue(std::uint64_t value)
{
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << value;
    return out.str();
}

} // namespace

Fingerprints computeFingerprints(const Mesh& mesh,
                                 const SparseMatrix& system,
                                 const CaseConfig& physics,
                                 const std::vector<int>& interfaceGlobalDofs)
{
    Fingerprints result;
    Fnv1a meshHash;
    const std::uint64_t nodeCount = static_cast<std::uint64_t>(mesh.nodes.size());
    const std::uint64_t tetCount = static_cast<std::uint64_t>(mesh.tets.size());
    meshHash.add(nodeCount);
    meshHash.add(tetCount);
    for (const Node& node : mesh.nodes) {
        meshHash.add(node.p.x);
        meshHash.add(node.p.y);
        meshHash.add(node.p.z);
        meshHash.add(node.subdomain);
        meshHash.add(node.dirichlet);
        meshHash.add(node.dirichletValue);
    }
    for (const Tet& tet : mesh.tets) {
        for (int dof : tet.dof) {
            meshHash.add(dof);
        }
        meshHash.add(tet.subdomain);
        meshHash.add(tet.domainEntity);
    }
    result.mesh = meshHash.value();

    Fnv1a systemHash;
    systemHash.add(system.n);
    system.forEachEntry([&](int row, int column, double value) {
        systemHash.add(row);
        systemHash.add(column);
        systemHash.add(value);
    });
    result.system = systemHash.value();

    Fnv1a interfaceHash;
    const std::uint64_t interfaceCount =
        static_cast<std::uint64_t>(interfaceGlobalDofs.size());
    interfaceHash.add(interfaceCount);
    for (int dof : interfaceGlobalDofs) {
        interfaceHash.add(dof);
    }
    result.interfaceOrdering = interfaceHash.value();

    Fnv1a sourceHash;
    sourceHash.add(physics.thermalSourceScale);
    const std::uint64_t sourceCount =
        static_cast<std::uint64_t>(physics.heatSources.size());
    sourceHash.add(sourceCount);
    for (const HeatSource& source : physics.heatSources) {
        sourceHash.add(source.subdomain);
        sourceHash.add(source.domainEntity);
        sourceHash.add(source.heatRateW);
    }
    const std::uint64_t dirichletCount =
        static_cast<std::uint64_t>(physics.dirichletConditions.size());
    sourceHash.add(dirichletCount);
    for (const BoundaryCondition& condition : physics.dirichletConditions) {
        sourceHash.add(condition.subdomain);
        sourceHash.add(condition.boundaryEntity);
        sourceHash.add(condition.temperature);
    }
    const std::uint64_t convectionCount =
        static_cast<std::uint64_t>(physics.convectionConditions.size());
    sourceHash.add(convectionCount);
    for (const ConvectionCondition& condition : physics.convectionConditions) {
        sourceHash.add(condition.subdomain);
        sourceHash.add(condition.boundaryEntity);
        sourceHash.add(condition.coefficient);
        sourceHash.add(condition.ambientTemperature);
    }
    result.sources = sourceHash.value();
    return result;
}

void saveModel(const ReducedSchurModel& model,
               const std::filesystem::path& directory)
{
    std::filesystem::create_directories(directory);
    std::ofstream metadata(directory / "metadata.json");
    if (!metadata) {
        throw std::runtime_error("[MOR] Cannot write model metadata.");
    }
    metadata << std::setprecision(17)
        << "{\n"
        << "  \"format_version\": 1,\n"
        << "  \"model_type\": \"exact-interface-reduced-schur\",\n"
        << "  \"interface_dofs\": " << model.interfaceDofs << ",\n"
        << "  \"rank\": " << model.rank << ",\n"
        << "  \"source_channels\": " << model.sourceChannels << ",\n"
        << "  \"mesh_fingerprint\": \"" << hexValue(model.fingerprints.mesh) << "\",\n"
        << "  \"system_fingerprint\": \"" << hexValue(model.fingerprints.system) << "\",\n"
        << "  \"interface_fingerprint\": \"" << hexValue(model.fingerprints.interfaceOrdering) << "\",\n"
        << "  \"source_fingerprint\": \"" << hexValue(model.fingerprints.sources) << "\",\n"
        << "  \"reduced_symmetry_error\": " << model.symmetryError << ",\n"
        << "  \"reduced_min_symmetric_eigenvalue\": " << model.minimumSymmetricEigenvalue << ",\n"
        << "  \"reduced_max_symmetric_eigenvalue\": " << model.maximumSymmetricEigenvalue << ",\n"
        << "  \"reduced_condition_estimate\": " << model.conditionEstimate << ",\n"
        << "  \"exact_schur_apply_seconds\": " << model.exactSchurApplySeconds << ",\n"
        << "  \"reduced_assembly_seconds\": " << model.reducedAssemblySeconds << "\n"
        << "}\n";
    writeVector(directory / "reference_interface.bin", model.referenceInterface);
    writeVector(directory / "basis.bin", model.basis);
    writeVector(directory / "schur_basis.bin", model.schurBasis);
    writeVector(directory / "reduced_operator.bin", model.reducedOperator);
    writeVector(directory / "singular_values.bin", model.singularValues);
}

ReducedSchurModel loadModel(const std::filesystem::path& directory,
                            const Fingerprints& expected,
                            int expectedInterfaceDofs,
                            int expectedSourceChannels)
{
    const std::string metadata = readText(directory / "metadata.json");
    if (jsonUnsigned(metadata, "format_version") != 1) {
        throw std::runtime_error("[MOR] Unsupported reduced model format.");
    }
    ReducedSchurModel model;
    model.interfaceDofs = static_cast<int>(jsonUnsigned(metadata, "interface_dofs"));
    model.rank = static_cast<int>(jsonUnsigned(metadata, "rank"));
    model.sourceChannels = static_cast<int>(jsonUnsigned(metadata, "source_channels"));
    model.fingerprints.mesh = jsonUnsigned(metadata, "mesh_fingerprint");
    model.fingerprints.system = jsonUnsigned(metadata, "system_fingerprint");
    model.fingerprints.interfaceOrdering = jsonUnsigned(metadata, "interface_fingerprint");
    model.fingerprints.sources = jsonUnsigned(metadata, "source_fingerprint");
    model.symmetryError = jsonDouble(metadata, "reduced_symmetry_error");
    model.minimumSymmetricEigenvalue = jsonDouble(metadata, "reduced_min_symmetric_eigenvalue");
    model.maximumSymmetricEigenvalue = jsonDouble(metadata, "reduced_max_symmetric_eigenvalue");
    model.conditionEstimate = jsonDouble(metadata, "reduced_condition_estimate");
    if (model.interfaceDofs != expectedInterfaceDofs
        || model.sourceChannels != expectedSourceChannels
        || model.fingerprints.mesh != expected.mesh
        || model.fingerprints.system != expected.system
        || model.fingerprints.interfaceOrdering != expected.interfaceOrdering
        || model.fingerprints.sources != expected.sources) {
        throw std::runtime_error(
            "[MOR] Reduced model fingerprint mismatch: mesh/matrix/interface/source ordering changed.");
    }
    model.referenceInterface = readVector(directory / "reference_interface.bin",
        static_cast<std::size_t>(model.interfaceDofs));
    model.basis = readVector(directory / "basis.bin",
        static_cast<std::size_t>(model.interfaceDofs * model.rank));
    model.schurBasis = readVector(directory / "schur_basis.bin",
        static_cast<std::size_t>(model.interfaceDofs * model.rank));
    model.reducedOperator = readVector(directory / "reduced_operator.bin",
        static_cast<std::size_t>(model.rank * model.rank));
    model.singularValues = readVector(directory / "singular_values.bin",
        static_cast<std::size_t>(model.rank));
    return model;
}

} // namespace mor
