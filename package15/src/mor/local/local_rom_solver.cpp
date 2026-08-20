#include "local_rom_solver.hpp"

#include "sipg_core.hpp"
#include "config_io.hpp"
#include "mesh_loader.hpp"
#include "fem_assembly.hpp"
#include "linear_solvers.hpp"
#include "diagnostics_io.hpp"
#include "ddm_schur/interface_operator.hpp"
#include "ddm_schur/schur_fgmres.hpp"
#include "mor/mor_diagnostics.hpp"
#include "mor/model_io.hpp"
#include "mor/source_parameterization.hpp"
#include "local_reduced_blocks.hpp"
#include "local_reduced_schur.hpp"
#include "local_rom_serialization.hpp"
#include "local_snapshot_basis.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <set>
#include <stdexcept>
#include <tuple>
#include <type_traits>

namespace mor::local {
namespace {

struct EvaluationRow {
    int caseIndex = -1;
    std::string split;
    std::string family;
    SolveResult solve;
    ErrorMetrics pureError;
    double temperatureJump = 0.0;
    double relativeFluxImbalance = 0.0;
    double fluxImbalanceL2 = 0.0;
    double maximumFluxImbalance = 0.0;
    int worstInterfaceId = -1;
    int worstFacePairId = -1;
    int worstIntegrationTriangle = -1;
    double fomSeconds = 0.0;
    double sourceProjectionSeconds = 0.0;
    double diagnosticSeconds = 0.0;
};

struct MultiRhsTimingRow {
    int rhsCount = 0;
    double sourceProjectionSeconds = 0.0;
    double reducedAssemblySeconds = 0.0;
    double interfaceSolveSeconds = 0.0;
    double localBackSubstitutionSeconds = 0.0;
    double reconstructionSeconds = 0.0;
    double totalOnlineSeconds = 0.0;
    double checksum = 0.0;
};

double norm2(const std::vector<double>& values)
{
    long double sum = 0.0L;
    for (double value : values) {
        sum += static_cast<long double>(value) * value;
    }
    return std::sqrt(static_cast<double>(sum));
}

class LocalHash {
public:
    template <typename T>
    void add(const T& value)
    {
        static_assert(std::is_trivially_copyable_v<T>);
        const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
        for (std::size_t i = 0; i < sizeof(T); ++i) {
            value_ ^= static_cast<std::uint64_t>(bytes[i]);
            value_ *= 1099511628211ULL;
        }
    }

    void addString(const std::string& value)
    {
        add(value.size());
        for (unsigned char byte : value) {
            value_ ^= static_cast<std::uint64_t>(byte);
            value_ *= 1099511628211ULL;
        }
    }

    std::uint64_t value() const { return value_; }

private:
    std::uint64_t value_ = 1469598103934665603ULL;
};

std::string hexadecimal(std::uint64_t value)
{
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << value;
    return out.str();
}

void addMaterial(LocalHash& hash, const Material& material)
{
    hash.addString(material.name);
    hash.add(material.conductivityX);
    hash.add(material.conductivityY);
    hash.add(material.conductivityZ);
    hash.add(material.density);
    hash.add(material.heatCapacity);
}

void assignLocalFingerprintsAndTemplates(
    const Mesh& mesh,
    const CaseConfig& physics,
    Model& model,
    bool reuseTemplates)
{
    std::map<std::tuple<std::uint64_t, std::uint64_t, std::uint64_t, int>, int> ids;
    std::map<int, std::size_t> representative;
    int nextTemplate = 0;
    for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
        SubdomainModel& local = model.subdomains[slot];
        const int domain = local.subdomain;
        std::vector<int> localNode(mesh.nodes.size(), -1);
        std::vector<int> nodes;
        Vec3 minimum{std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max(),
                     std::numeric_limits<double>::max()};
        for (std::size_t global = 0; global < mesh.nodes.size(); ++global) {
            if (mesh.nodes[global].subdomain != domain) continue;
            localNode[global] = static_cast<int>(nodes.size());
            nodes.push_back(static_cast<int>(global));
            minimum.x = std::min(minimum.x, mesh.nodes[global].p.x);
            minimum.y = std::min(minimum.y, mesh.nodes[global].p.y);
            minimum.z = std::min(minimum.z, mesh.nodes[global].p.z);
        }
        LocalHash meshHash;
        meshHash.add(nodes.size());
        for (int global : nodes) {
            const Node& node = mesh.nodes[static_cast<std::size_t>(global)];
            // Canonicalize translated copies at 1e-14 m. This removes only
            // floating-point offset roundoff; a physical mesh change remains distinct.
            constexpr double meshFingerprintScale = 1.0e14;
            meshHash.add(static_cast<std::int64_t>(std::llround(
                (node.p.x - minimum.x) * meshFingerprintScale)));
            meshHash.add(static_cast<std::int64_t>(std::llround(
                (node.p.y - minimum.y) * meshFingerprintScale)));
            meshHash.add(static_cast<std::int64_t>(std::llround(
                (node.p.z - minimum.z) * meshFingerprintScale)));
            meshHash.add(node.sourceVertex);
        }
        std::size_t localTets = 0;
        for (const Tet& tet : mesh.tets) {
            if (tet.subdomain != domain) continue;
            ++localTets;
            meshHash.add(tet.domainEntity);
            for (int global : tet.dof) {
                meshHash.add(localNode[static_cast<std::size_t>(global)]);
            }
        }
        meshHash.add(localTets);
        local.meshFingerprint = meshHash.value();

        LocalHash materialHash;
        const DomainConfig& domainConfig = physics.domains[static_cast<std::size_t>(domain)];
        addMaterial(materialHash, domainConfig.material);
        materialHash.add(domainConfig.materialsByDomainEntity.size());
        for (const auto& material : domainConfig.materialsByDomainEntity) {
            materialHash.add(material.first);
            addMaterial(materialHash, material.second);
        }
        local.materialFingerprint = materialHash.value();

        LocalHash boundaryHash;
        boundaryHash.addString(physics.interfaceScheme);
        boundaryHash.addString(physics.penaltyMode);
        boundaryHash.addString(physics.penaltyScaling);
        boundaryHash.add(physics.penaltyFactor);
        for (const BoundaryCondition& condition : physics.dirichletConditions) {
            if (condition.subdomain == -1 || condition.subdomain == domain) {
                boundaryHash.add(1);
                boundaryHash.add(condition.boundaryEntity);
                boundaryHash.add(condition.temperature);
            }
        }
        for (const ConvectionCondition& condition : physics.convectionConditions) {
            if (condition.subdomain == -1 || condition.subdomain == domain) {
                boundaryHash.add(2);
                boundaryHash.add(condition.boundaryEntity);
                boundaryHash.add(condition.coefficient);
                boundaryHash.add(condition.ambientTemperature);
            }
        }
        for (const HeatSource& source : physics.heatSources) {
            if (source.subdomain == domain) {
                boundaryHash.add(3);
                boundaryHash.add(source.domainEntity);
            }
        }
        std::vector<std::tuple<int, std::vector<int>, std::vector<int>>> interfaces;
        for (const InterfaceConfig& interfaceConfig : physics.interfaces) {
            if (interfaceConfig.leftSubdomain == domain) {
                interfaces.emplace_back(1, interfaceConfig.leftBoundaryEntities,
                    interfaceConfig.rightBoundaryEntities);
            } else if (interfaceConfig.rightSubdomain == domain) {
                interfaces.emplace_back(-1, interfaceConfig.rightBoundaryEntities,
                    interfaceConfig.leftBoundaryEntities);
            }
        }
        std::sort(interfaces.begin(), interfaces.end());
        boundaryHash.add(interfaces.size());
        for (const auto& interfaceSignature : interfaces) {
            boundaryHash.add(std::get<0>(interfaceSignature));
            for (int boundary : std::get<1>(interfaceSignature)) boundaryHash.add(boundary);
            boundaryHash.add(-7);
            for (int boundary : std::get<2>(interfaceSignature)) boundaryHash.add(boundary);
            boundaryHash.add(-11);
        }
        local.boundaryInterfaceFingerprint = boundaryHash.value();

        const auto key = std::make_tuple(local.meshFingerprint,
            local.materialFingerprint, local.boundaryInterfaceFingerprint, local.rank);
        auto found = ids.find(key);
        if (found == ids.end()) {
            local.templateId = nextTemplate;
            ids.emplace(key, nextTemplate);
            representative.emplace(nextTemplate, slot);
            ++nextTemplate;
        } else {
            local.templateId = found->second;
            if (reuseTemplates) {
                const SubdomainModel& source = model.subdomains[
                    representative.at(local.templateId)];
                if (source.interiorDofs != local.interiorDofs
                    || source.localInterfaceDofs != local.localInterfaceDofs
                    || source.rank != local.rank
                    || source.basis.size() != local.basis.size()) {
                    throw std::runtime_error(
                        "[Local ROM] Structurally identical template has incompatible dimensions.");
                }
                local.basis = source.basis;
                local.orthogonalityError = source.orthogonalityError;
                local.templateReused = true;
            }
        }
    }
}

void writeStructureDiagnostics(const Mesh& mesh,
                               const CaseConfig& physics,
                               const ddm_schur::InterfacePartition& partition,
                               const Model& model,
                               const std::filesystem::path& outputDirectory)
{
    std::map<int, int> traceMultiplicity;
    std::size_t traceReferences = 0;
    for (const SubdomainModel& local : model.subdomains) {
        for (int global : local.interfaceGlobalDofs) {
            ++traceMultiplicity[global];
            ++traceReferences;
        }
    }
    std::size_t duplicateReferences = 0;
    for (const auto& value : traceMultiplicity) {
        if (value.second > 1) {
            duplicateReferences += static_cast<std::size_t>(value.second - 1);
        }
    }
    std::ofstream domains(outputDirectory / "local_rom_subdomain_structure.csv");
    domains << "subdomain_id,interior_dofs,local_interface_dofs,global_unique_interface_dofs,"
        << "all_local_trace_references,unique_trace_union,duplicate_trace_references,"
        << "mesh_fingerprint,material_fingerprint,boundary_interface_fingerprint,"
        << "template_id,template_reused\n";
    for (const SubdomainModel& local : model.subdomains) {
        domains << local.subdomain << ',' << local.interiorDofs << ','
            << local.localInterfaceDofs << ',' << partition.interfaceGlobalDofs.size() << ','
            << traceReferences << ',' << traceMultiplicity.size() << ','
            << duplicateReferences << ',' << hexadecimal(local.meshFingerprint) << ','
            << hexadecimal(local.materialFingerprint) << ','
            << hexadecimal(local.boundaryInterfaceFingerprint) << ','
            << local.templateId << ',' << (local.templateReused ? 1 : 0) << '\n';
    }

    std::ofstream interfaces(outputDirectory / "local_rom_interface_structure.csv");
    interfaces << "interface_pair_id,master_subdomain,slave_subdomain,"
        << "master_boundary_entity_count,slave_boundary_entity_count,face_pair_count,"
        << "triangle_overlap_count,master_area_m2,slave_area_m2,matched_overlap_area_m2\n"
        << std::setprecision(17);
    for (std::size_t index = 0; index < mesh.interfaceSummaries.size(); ++index) {
        const InterfaceBuildSummary& summary = mesh.interfaceSummaries[index];
        interfaces << index << ',' << summary.leftSubdomain << ','
            << summary.rightSubdomain << ',' << summary.leftBoundaryEntityCount << ','
            << summary.rightBoundaryEntityCount << ',' << summary.facePairCount << ','
            << summary.integrationTriangleCount << ',' << summary.leftArea << ','
            << summary.rightArea << ',' << summary.matchedOverlapArea << '\n';
    }

    std::ofstream definition(outputDirectory / "local_rom_partition_definition.csv");
    definition << "subdomain_id,physical_role,geometric_domain_count,geometric_domain_ids,"
        << "material_region_count,material_names,power_channel_count,power_domain_entities,"
        << "neighbor_count,neighbor_subdomains,interior_dofs,interface_dofs\n";
    for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
        const ddm_schur::DomainBlocks& domain = partition.domains[slot];
        std::set<int> geometricDomains;
        std::set<std::string> materials;
        for (const Tet& tet : mesh.tets) {
            if (tet.subdomain == domain.domainId) {
                geometricDomains.insert(tet.domainEntity);
                materials.insert(materialForTet(physics, tet).name);
            }
        }
        std::vector<int> powerEntities;
        for (const HeatSource& source : physics.heatSources) {
            if (source.subdomain == domain.domainId) {
                powerEntities.push_back(source.domainEntity);
            }
        }
        std::string role = "physical_ddm_module";
        std::string joinedMaterials;
        for (const std::string& material : materials) {
            if (!joinedMaterials.empty()) joinedMaterials += ';';
            joinedMaterials += material;
        }
        std::string lowerMaterials = joinedMaterials;
        std::transform(lowerMaterials.begin(), lowerMaterials.end(),
            lowerMaterials.begin(), [](unsigned char value) {
                return static_cast<char>(std::tolower(value));
            });
        if (lowerMaterials.find("heatsink") != std::string::npos
            || lowerMaterials.find("aluminum") != std::string::npos) {
            role = "heatsink_lid_stack";
        } else if (powerEntities.size() > 1
                   && lowerMaterials.find("silicon") != std::string::npos) {
            role = "multi_chiplet_die_interconnect_assembly";
        } else if (!powerEntities.empty()) {
            role = "powered_die_or_active_layer";
        }
        std::vector<int> neighbors;
        for (const auto& neighbor : domain.interfaceGlobalDofsByNeighbor) {
            neighbors.push_back(neighbor.first);
        }
        const auto writeIntList = [&](const auto& values) {
            bool first = true;
            for (const auto& value : values) {
                if (!first) definition << ';';
                definition << value;
                first = false;
            }
        };
        definition << domain.domainId << ',' << role << ','
            << geometricDomains.size() << ",\"";
        writeIntList(geometricDomains);
        definition << "\"," << materials.size() << ",\"" << joinedMaterials
            << "\"," << powerEntities.size() << ",\"";
        writeIntList(powerEntities);
        definition << "\"," << neighbors.size() << ",\"";
        writeIntList(neighbors);
        definition << "\"," << domain.interiorGlobalDofs.size() << ','
            << domain.interfaceGlobalDofs.size() << '\n';
    }
}

struct TemplateStorageSummary {
    int uniqueTemplates = 0;
    int reusedInstances = 0;
    std::size_t withoutReuseBytes = 0;
    std::size_t withReuseBytes = 0;
};

std::size_t reusablePayloadBytes(const SubdomainModel& local)
{
    const std::size_t interior = static_cast<std::size_t>(local.interiorDofs);
    const std::size_t trace = static_cast<std::size_t>(local.localInterfaceDofs);
    const std::size_t rank = static_cast<std::size_t>(local.rank);
    return sizeof(double) * (interior * rank + rank * rank
        + 2 * rank * trace + 2 * interior + trace);
}

TemplateStorageSummary writeTemplateStorage(
    const Model& model,
    bool reuseEnabled,
    const std::filesystem::path& outputDirectory)
{
    TemplateStorageSummary result;
    std::set<int> templates;
    std::set<int> storedTemplates;
    for (const SubdomainModel& local : model.subdomains) {
        templates.insert(local.templateId);
        result.withoutReuseBytes += reusablePayloadBytes(local);
        if (local.templateReused) {
            ++result.reusedInstances;
        }
        if (storedTemplates.insert(local.templateId).second) {
            result.withReuseBytes += reusablePayloadBytes(local);
        }
    }
    result.uniqueTemplates = static_cast<int>(templates.size());
    std::ofstream out(outputDirectory / "local_rom_template_reuse.csv");
    out << "reuse_enabled,unique_template_count,subdomain_instance_count,"
        << "reused_instance_count,basis_storage_without_reuse_bytes,"
        << "basis_storage_with_reuse_bytes,model_size_reduction,"
        << "online_relative_l2_difference,maximum_temperature_difference_k\n"
        << (reuseEnabled ? 1 : 0) << ',' << result.uniqueTemplates << ','
        << model.subdomains.size() << ',' << result.reusedInstances << ','
        << result.withoutReuseBytes << ',' << result.withReuseBytes << ','
        << (1.0 - static_cast<double>(result.withReuseBytes)
            / std::max(1.0, static_cast<double>(result.withoutReuseBytes)))
        << ",nan,nan\n";
    return result;
}

ErrorMetrics errors(const SparseMatrix& system,
                    const std::vector<double>& rhs,
                    const std::vector<double>& approximation,
                    const std::vector<double>* truth)
{
    ErrorMetrics result;
    std::vector<double> residual = system.multiply(approximation);
    for (std::size_t row = 0; row < residual.size(); ++row) {
        residual[row] -= rhs[row];
    }
    result.globalRelativeResidual = norm2(residual)
        / std::max(1.0e-300, norm2(rhs));
    if (truth != nullptr) {
        long double errorSquared = 0.0L;
        long double truthSquared = 0.0L;
        double maximum = 0.0;
        for (std::size_t row = 0; row < approximation.size(); ++row) {
            const double difference = approximation[row] - (*truth)[row];
            errorSquared += static_cast<long double>(difference) * difference;
            truthSquared += static_cast<long double>((*truth)[row]) * (*truth)[row];
            maximum = std::max(maximum, std::abs(difference));
        }
        result.relativeL2 = std::sqrt(static_cast<double>(errorSquared))
            / std::max(1.0e-300, std::sqrt(static_cast<double>(truthSquared)));
        result.maximumAbsolute = maximum;
        const auto approximateMaximum = std::max_element(
            approximation.begin(), approximation.end());
        const auto truthMaximum = std::max_element(truth->begin(), truth->end());
        result.maximumTemperatureError = std::abs(*approximateMaximum - *truthMaximum);
    }
    return result;
}

double exactInterfaceResidual(ddm_schur::DdmSchurSolver& solver,
                              const std::vector<double>& rhs,
                              const std::vector<double>& interfaceTemperature)
{
    const std::vector<double> condensed = solver.condensedRhs(rhs);
    std::vector<double> image;
    solver.applyExactSchur(interfaceTemperature, image);
    for (std::size_t row = 0; row < image.size(); ++row) {
        image[row] = condensed[row] - image[row];
    }
    return norm2(image) / std::max(1.0e-300, norm2(condensed));
}

double physicalInterfaceResidual(const SparseMatrix& system,
                                 const std::vector<double>& rhs,
                                 const std::vector<double>& temperature,
                                 const std::vector<int>& interfaceGlobalDofs)
{
    const std::vector<double> image = system.multiply(temperature);
    long double residualSquared = 0.0L;
    for (int global : interfaceGlobalDofs) {
        const double value = image[static_cast<std::size_t>(global)]
            - rhs[static_cast<std::size_t>(global)];
        residualSquared += static_cast<long double>(value) * value;
    }
    return std::sqrt(static_cast<double>(residualSquared))
        / std::max(1.0e-300, norm2(rhs));
}

std::vector<SnapshotDatabase> solveLocalTrainingBatch(
    SubdomainDirectSolver& solver,
    const SourceParameterization& sources,
    const std::vector<ParameterCase>& cases,
    const ddm_schur::InterfacePartition& partition,
    const std::vector<double>& referenceTemperature)
{
    const int globalDofs = partition.totalDofs;
    if (referenceTemperature.size() != static_cast<std::size_t>(globalDofs)) {
        throw std::runtime_error(
            "[Local ROM] Training reference has the wrong size.");
    }
    std::vector<SnapshotDatabase> result(partition.domains.size());
    for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
        SnapshotDatabase& local = result[slot];
        local.rows = static_cast<int>(
            partition.domains[slot].interiorGlobalDofs.size());
        local.cases = cases;
        local.values.assign(static_cast<std::size_t>(local.rows)
            * cases.size(), 0.0);
    }

    const int blockSize = globalDofs > 500000
        ? 4 : std::max(1, static_cast<int>(cases.size()));
    for (std::size_t first = 0; first < cases.size();
         first += static_cast<std::size_t>(blockSize)) {
        const int active = std::min(
            blockSize, static_cast<int>(cases.size() - first));
        std::vector<double> batch(static_cast<std::size_t>(globalDofs)
            * static_cast<std::size_t>(active), 0.0);
        for (int rhs = 0; rhs < active; ++rhs) {
            const std::vector<double> source = composeRhs(
                sources, cases[first + static_cast<std::size_t>(rhs)].powersW);
            std::copy(source.begin(), source.end(),
                batch.begin() + static_cast<std::ptrdiff_t>(
                    static_cast<std::size_t>(rhs)
                    * static_cast<std::size_t>(globalDofs)));
        }
        std::vector<double> solutions;
        solver.solveMultiple(batch, active, solutions);
        for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
            const ddm_schur::DomainBlocks& domain = partition.domains[slot];
            SnapshotDatabase& local = result[slot];
            for (int rhs = 0; rhs < active; ++rhs) {
                const std::size_t sample = first + static_cast<std::size_t>(rhs);
                const std::size_t solutionOffset = static_cast<std::size_t>(rhs)
                    * static_cast<std::size_t>(globalDofs);
                const std::size_t snapshotOffset = sample
                    * static_cast<std::size_t>(local.rows);
                for (int row = 0; row < local.rows; ++row) {
                    const int global = domain.interiorGlobalDofs[
                        static_cast<std::size_t>(row)];
                    local.values[snapshotOffset + static_cast<std::size_t>(row)] =
                        solutions[solutionOffset + static_cast<std::size_t>(global)]
                        - referenceTemperature[static_cast<std::size_t>(global)];
                }
            }
        }
    }
    return result;
}

std::vector<ParameterCase> makeLocalTrainingCases(
    const SourceParameterization& sources,
    int requestedCount,
    std::uint64_t seed)
{
    const int channels = static_cast<int>(sources.channels.size());
    const int required = channels + (channels > 1 ? 7 : 1);
    std::vector<ParameterCase> cases = makeParameterCases(
        sources, "training", std::max(requestedCount, required), seed, true);
    auto replaceCase = [&](int offset, const std::string& family,
                           const std::vector<double>& powers) {
        const std::size_t slot = static_cast<std::size_t>(channels + offset);
        if (slot >= cases.size()) return;
        cases[slot].family = family;
        cases[slot].powersW = powers;
    };
    if (channels > 1) {
        std::vector<double> powers(static_cast<std::size_t>(channels), 0.0);
        const int middle = channels / 2;
        powers[static_cast<std::size_t>(std::max(0, middle - 1))] =
            sources.channels[static_cast<std::size_t>(std::max(0, middle - 1))].nominalPowerW;
        powers[static_cast<std::size_t>(middle)] =
            sources.channels[static_cast<std::size_t>(middle)].nominalPowerW;
        replaceCase(0, "adjacent_pair", powers);

        powers.assign(static_cast<std::size_t>(channels), 0.0);
        powers.front() = sources.channels.front().nominalPowerW;
        powers.back() = sources.channels.back().nominalPowerW;
        replaceCase(1, "first_last_pair", powers);

        powers.assign(static_cast<std::size_t>(channels), 0.0);
        powers[static_cast<std::size_t>(middle)] =
            sources.channels[static_cast<std::size_t>(middle)].nominalPowerW;
        replaceCase(2, "middle_source", powers);

        powers.assign(static_cast<std::size_t>(channels), 0.0);
        for (int channel = 0; channel < channels; ++channel) {
            if (channel % 3 == 0) {
                powers[static_cast<std::size_t>(channel)] =
                    sources.channels[static_cast<std::size_t>(channel)].nominalPowerW;
            }
        }
        replaceCase(3, "multiple_hotspots", powers);

        powers.assign(static_cast<std::size_t>(channels), 0.0);
        for (int channel = 0; channel < channels; ++channel) {
            powers[static_cast<std::size_t>(channel)] =
                sources.channels[static_cast<std::size_t>(channel)].nominalPowerW;
        }
        replaceCase(4, "all_uniform", powers);

        powers.assign(static_cast<std::size_t>(channels), 0.0);
        double nominalTotal = 0.0;
        for (const SourceChannel& source : sources.channels) {
            nominalTotal += source.nominalPowerW;
        }
        powers.front() = 0.75 * nominalTotal;
        powers.back() = 0.25 * nominalTotal;
        replaceCase(5, "same_total_endpoint_skew", powers);

        powers.assign(static_cast<std::size_t>(channels), 0.0);
        powers.front() = 1.3 * sources.channels.front().nominalPowerW;
        powers[1] = 0.7 * sources.channels[1].nominalPowerW;
        replaceCase(6, "remote_interface_driven", powers);
    }
    return cases;
}

std::vector<int> readRankFile(const std::filesystem::path& path,
                              const ddm_schur::InterfacePartition& partition)
{
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("[Local ROM] Cannot open local rank file: " + path.string());
    }
    std::map<int, int> bySubdomain;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::replace(line.begin(), line.end(), ';', ',');
        std::istringstream stream(line);
        std::string first;
        std::string second;
        if (!std::getline(stream, first, ',') || !std::getline(stream, second, ',')) continue;
        try {
            bySubdomain[std::stoi(first)] = std::stoi(second);
        } catch (const std::exception&) {
            // Allow a header such as subdomain,rank.
        }
    }
    std::vector<int> ranks;
    ranks.reserve(partition.domains.size());
    for (const ddm_schur::DomainBlocks& domain : partition.domains) {
        const auto found = bySubdomain.find(domain.domainId);
        if (found == bySubdomain.end() || found->second <= 0) {
            throw std::runtime_error(
                "[Local ROM] Rank file must provide a positive rank for every subdomain.");
        }
        ranks.push_back(found->second);
    }
    return ranks;
}

void writeRankTable(const Model& model,
                    const std::vector<double>& worstRelative,
                    const std::vector<double>& worstMaximum,
                    const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "subdomain,subdomain_id,interior_dofs,interface_dofs,snapshots,snapshot_count,"
        << "numerical_rank,numerical_snapshot_rank,selected_rank,orthogonality_error,"
        << "basis_orthogonality,largest_singular_value,smallest_retained_singular_value,"
        << "retained_energy,template_id,template_reused,template_consistency_difference,"
        << "snapshot_extraction_seconds,"
        << "pod_seconds,projection_seconds,reduced_interior_symmetry_error,"
        << "reduced_interior_min_eigenvalue,reduced_interior_max_eigenvalue,"
        << "reduced_interior_condition_estimate,coupling_symmetry_error,local_schur_symmetry_error,"
        << "worst_validation_relative_l2,worst_validation_max_node_k,singular_values\n"
        << std::setprecision(17);
    for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
        const SubdomainModel& local = model.subdomains[slot];
        const double energy = local.rank > 0
            && static_cast<std::size_t>(local.rank) <= local.retainedEnergy.size()
            ? local.retainedEnergy[static_cast<std::size_t>(local.rank - 1)] : 0.0;
        const double largestSingular = local.singularValues.empty()
            ? 0.0 : local.singularValues.front();
        const double smallestRetained = local.rank > 0
            && static_cast<std::size_t>(local.rank) <= local.singularValues.size()
            ? local.singularValues[static_cast<std::size_t>(local.rank - 1)] : 0.0;
        out << local.subdomain << ',' << local.subdomain << ',' << local.interiorDofs << ','
            << local.localInterfaceDofs << ',' << local.snapshots << ',' << local.snapshots << ','
            << local.numericalRank << ',' << local.numericalRank << ',' << local.rank << ','
            << local.orthogonalityError << ',' << local.orthogonalityError << ','
            << largestSingular << ',' << smallestRetained << ',' << energy << ','
            << local.templateId << ',' << (local.templateReused ? 1 : 0) << ','
            << local.templateConsistencyDifference << ','
            << local.snapshotExtractionSeconds << ',' << local.podSeconds << ','
            << local.projectionSeconds << ',' << local.reducedInteriorSymmetryError << ','
            << local.reducedInteriorMinimumEigenvalue << ','
            << local.reducedInteriorMaximumEigenvalue << ','
            << local.reducedInteriorConditionEstimate << ','
            << local.couplingSymmetryError << ',' << local.localSchurSymmetryError << ','
            << worstRelative[slot] << ',' << worstMaximum[slot] << ",\"";
        for (std::size_t mode = 0; mode < local.singularValues.size(); ++mode) {
            if (mode > 0) out << ';';
            out << local.singularValues[mode];
        }
        out << "\"\n";
    }
}

void writeAccuracy(const std::vector<EvaluationRow>& rows,
                   const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "case,split,family,mode,interface_fgmres_iterations,interface_fgmres_matvecs,"
        << "interface_fgmres_true_residual,zero_guess_iterations,local_rom_guess_iterations,"
        << "interface_relative_residual,global_relative_residual,relative_l2,"
        << "max_node_error_k,max_temperature_error_k,interface_temperature_jump_rms_k,"
        << "relative_flux_imbalance,area_weighted_flux_imbalance_l2_w_m2,"
        << "maximum_flux_imbalance_w_m2,worst_interface_id,worst_face_pair_id,worst_integration_triangle,"
        << "relative_flux_floor_w_m2,initial_true_residual,final_true_residual,"
        << "fom_seconds,online_seconds,status\n" << std::setprecision(17);
    for (const EvaluationRow& row : rows) {
        out << row.caseIndex << ',' << row.split << ',' << row.family << ','
            << (row.solve.correctedIterations >= 0 ? "corrected" : "pure") << ','
            << row.solve.timing.interfaceIterations << ','
            << row.solve.timing.interfaceMatvecs << ','
            << row.solve.timing.interfaceRelativeResidual << ','
            << row.solve.zeroGuessIterations << ',' << row.solve.correctedIterations << ','
            << row.solve.error.interfaceRelativeResidual << ','
            << row.solve.error.globalRelativeResidual << ','
            << row.solve.error.relativeL2 << ',' << row.solve.error.maximumAbsolute << ','
            << row.solve.error.maximumTemperatureError << ',' << row.temperatureJump << ','
            << row.relativeFluxImbalance << ',' << row.fluxImbalanceL2 << ','
            << row.maximumFluxImbalance << ',' << row.worstInterfaceId << ','
            << row.worstFacePairId << ',' << row.worstIntegrationTriangle << ",1e-12,"
            << row.pureError.interfaceRelativeResidual << ','
            << row.solve.correctedTrueResidual << ',' << row.fomSeconds << ','
            << row.solve.timing.totalSeconds << ',' << row.solve.status << '\n';
    }
}

void writeOnlineTiming(const std::vector<EvaluationRow>& rows,
                       const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "case,split,source_projection_seconds,local_reduced_assembly_seconds,"
        << "interface_fgmres_iterations,interface_fgmres_matvecs,interface_fgmres_true_residual,"
        << "interface_solve_seconds,proxy_solve_seconds,coarse_solve_seconds,"
        << "local_recovery_seconds,full_field_reconstruction_seconds,"
        << "correction_seconds,diagnostic_seconds,core_online_seconds,total_with_source_and_diagnostics_seconds\n"
        << std::setprecision(17);
    for (const EvaluationRow& row : rows) {
        const OnlineTiming& timing = row.solve.timing;
        out << row.caseIndex << ',' << row.split << ',' << row.sourceProjectionSeconds << ','
            << timing.localReducedAssemblySeconds << ','
            << timing.interfaceIterations << ',' << timing.interfaceMatvecs << ','
            << timing.interfaceRelativeResidual << ','
            << timing.interfaceSolveSeconds << ',' << timing.proxySolveSeconds << ','
            << timing.coarseSolveSeconds << ',' << timing.localRecoverySeconds << ','
            << timing.fullFieldReconstructionSeconds << ','
            << timing.correctionSeconds << ',' << row.diagnosticSeconds << ','
            << timing.totalSeconds << ','
            << (row.sourceProjectionSeconds + timing.totalSeconds + row.diagnosticSeconds) << '\n';
    }
}

void writeMultiRhsTiming(const std::vector<MultiRhsTimingRow>& rows,
                         const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "rhs_count,source_projection_total_seconds,source_projection_average_seconds,"
        << "local_reduced_assembly_total_seconds,interface_solve_total_seconds,"
        << "local_back_substitution_total_seconds,full_field_reconstruction_total_seconds,"
        << "online_total_seconds,average_online_seconds,checksum\n" << std::setprecision(17);
    for (const MultiRhsTimingRow& row : rows) {
        const double inverse = 1.0 / std::max(1, row.rhsCount);
        out << row.rhsCount << ',' << row.sourceProjectionSeconds << ','
            << row.sourceProjectionSeconds * inverse << ','
            << row.reducedAssemblySeconds << ',' << row.interfaceSolveSeconds << ','
            << row.localBackSubstitutionSeconds << ',' << row.reconstructionSeconds << ','
            << row.totalOnlineSeconds << ',' << row.totalOnlineSeconds * inverse << ','
            << row.checksum << '\n';
    }
}

} // namespace

WorkflowResult runLocalRomSchurWorkflow(
    const Mesh& mesh,
    const SparseMatrix& system,
    const CaseConfig& physics,
    const std::vector<double>& assembledSource,
    const std::vector<double>& heatOnlySource,
    const std::vector<double>& fixedAdjust,
    const ddm_schur::Options& schurOptions,
    const Options& options,
    const std::filesystem::path& outputDirectory)
{
    if (!options.generate && options.loadPath.empty()) {
        throw std::runtime_error("[Local ROM] Generate a model or provide --local-mor-load.");
    }
    if (options.method != "pod" && options.method != "rb") {
        throw std::runtime_error("[Local ROM] --local-mor-method must be pod or rb.");
    }
    if (options.mode != "pure" && options.mode != "corrected") {
        throw std::runtime_error("[Local ROM] --local-mor-mode must be pure or corrected.");
    }
    if (options.interfaceMode != "full") {
        throw std::runtime_error("[Local ROM] Steady Local-ROM requires --local-interface-mode full.");
    }
    std::filesystem::create_directories(outputDirectory);
    const auto offlineStart = std::chrono::steady_clock::now();

    const SourceParameterization sources = buildSourceParameterization(
        mesh, physics, assembledSource, heatOnlySource, fixedAdjust);
    if (sources.channels.empty()) {
        throw std::runtime_error("[Local ROM] At least one independent physical power channel is required.");
    }
    const auto partitionStart = std::chrono::steady_clock::now();
    const ddm_schur::InterfacePartition partition =
        ddm_schur::buildInterfacePartition(mesh, system);
    const double partitionSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - partitionStart).count();
    const Fingerprints fingerprints = computeFingerprints(
        mesh, system, physics, partition.interfaceGlobalDofs);

    double fomFactorWallSeconds = 0.0;
    std::unique_ptr<SubdomainDirectSolver> fomSolver;
    if (options.generate || options.compareFom) {
        const auto factorStart = std::chrono::steady_clock::now();
        fomSolver = std::make_unique<SubdomainDirectSolver>(
            system.size(), sparseMatrixEntries(system));
        fomFactorWallSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - factorStart).count();
    }
    Model model;
    double modelLoadSeconds = 0.0;
    int trainingSnapshotCount = 0;
    std::vector<ParameterCase> trainingCases;
    if (!options.loadPath.empty() && !options.generate) {
        const auto loadStart = std::chrono::steady_clock::now();
        model = loadLocalRomModel(
            options.loadPath, fingerprints, system.size(),
            static_cast<int>(partition.interfaceGlobalDofs.size()));
        modelLoadSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - loadStart).count();
        trainingSnapshotCount = model.subdomains.empty()
            ? 0 : model.subdomains.front().snapshots;
    } else {
        std::vector<double> referenceTemperature;
        fomSolver->solve(sources.referenceRhs, referenceTemperature);
        trainingCases = makeLocalTrainingCases(
            sources, options.trainingCases, options.seed);
        trainingSnapshotCount = static_cast<int>(trainingCases.size());
        writeSourceChannels(sources, outputDirectory / "local_rom_source_channels.csv");
        writeParameterCases(trainingCases,
            outputDirectory / "local_rom_training_cases.csv");
        const auto snapshotStart = std::chrono::steady_clock::now();
        std::vector<SnapshotDatabase> localTrainingSnapshots = solveLocalTrainingBatch(
            *fomSolver, sources, trainingCases, partition, referenceTemperature);
        const double snapshotSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - snapshotStart).count();
        Options basisOptions = options;
        if (!options.rankFile.empty()) {
            basisOptions.rankPerSubdomain = readRankFile(options.rankFile, partition);
        }
        model = buildIndependentSnapshotBasesFromLocalDatabases(
            partition, referenceTemperature,
            std::move(localTrainingSnapshots), basisOptions);
        model.fingerprints = fingerprints;
        model.snapshotSeconds = snapshotSeconds;
        assignLocalFingerprintsAndTemplates(
            mesh, physics, model, options.reuseIdenticalSubdomains);
        projectLocalReducedBlocks(model, partition);
    }
    LocalReducedSchurSolver online(
        model, mesh, physics, partition, schurOptions,
        options.matrixFreeInterfaceThreshold, outputDirectory);
    model.schurConstructionSeconds = online.assemblySeconds();
    model.modelBytes = estimateModelBytes(model);
    double serializationSeconds = 0.0;
    if (!options.savePath.empty()) {
        const auto serializationStart = std::chrono::steady_clock::now();
        saveLocalRomModel(model, options.savePath);
        serializationSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - serializationStart).count();
    }
    model.offlineSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - offlineStart).count();
    writeStructureDiagnostics(mesh, physics, partition, model, outputDirectory);
    const TemplateStorageSummary templateStorage = writeTemplateStorage(
        model, options.reuseIdenticalSubdomains, outputDirectory);

    std::unique_ptr<ddm_schur::DdmSchurSolver> exactSolver;
    double exactSchurSetupSeconds = 0.0;
    if (options.mode == "corrected") {
        ddm_schur::Options exactOptions = schurOptions;
        exactOptions.proxyOutputDirectory = outputDirectory.string();
        exactSolver = std::make_unique<ddm_schur::DdmSchurSolver>(
            mesh, system, physics, exactOptions);
    }

    std::vector<ParameterCase> evaluationCases;
    ParameterCase nominal;
    nominal.index = 0;
    nominal.split = "nominal";
    nominal.family = "nominal_power";
    for (const SourceChannel& channel : sources.channels) {
        nominal.powersW.push_back(channel.nominalPowerW);
    }
    evaluationCases.push_back(nominal);
    std::vector<ParameterCase> validation = makeParameterCases(
        sources, "validation", options.validationCases, options.seed + 101ULL, false);
    std::vector<ParameterCase> test = makeParameterCases(
        sources, "test", options.testCases, options.seed + 1009ULL, false);
    evaluationCases.insert(evaluationCases.end(), validation.begin(), validation.end());
    evaluationCases.insert(evaluationCases.end(), test.begin(), test.end());
    writeParameterCases(validation, outputDirectory / "local_rom_validation_cases.csv");
    writeParameterCases(test, outputDirectory / "local_rom_test_cases.csv");

    std::vector<EvaluationRow> rows;
    std::vector<double> worstLocalRelative(model.subdomains.size(), 0.0);
    std::vector<double> worstLocalMaximum(model.subdomains.size(), 0.0);
    std::ofstream fluxOut(outputDirectory / "local_rom_interface_flux.csv");
    fluxOut << "case,split,family,power_vector_w,physical_interface_id,face_pair_id,integration_triangle_id,left_subdomain,right_subdomain,left_boundary_entity,"
        << "right_boundary_entity,area_m2,temperature_jump_rms_k,"
        << "left_physical_normal_flux_w_m2,right_physical_normal_flux_w_m2,"
        << "sipg_numerical_flux_w_m2,flux_imbalance_l2_w_m2,"
        << "relative_flux_imbalance,relative_flux_floor_w_m2,"
        << "maximum_flux_imbalance_w_m2,worst_integration_triangle\n"
        << std::setprecision(17);

    bool accepted = true;
    for (std::size_t caseSlot = 0; caseSlot < evaluationCases.size(); ++caseSlot) {
        const ParameterCase& parameterCase = evaluationCases[caseSlot];
        const auto sourceStart = std::chrono::steady_clock::now();
        const std::vector<double> rhs = composeRhs(sources, parameterCase.powersW);
        const double sourceSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - sourceStart).count();
        SolveResult reduced = online.solve(rhs);

        std::vector<double> truth;
        double fomSeconds = 0.0;
        if (options.compareFom) {
            const auto fomStart = std::chrono::steady_clock::now();
            fomSolver->solve(rhs, truth);
            fomSeconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - fomStart).count();
        }
        reduced.error = errors(system, rhs, reduced.temperature,
            truth.empty() ? nullptr : &truth);
        reduced.error.interfaceRelativeResidual = exactSolver
            ? exactInterfaceResidual(*exactSolver, rhs, reduced.interfaceTemperature)
            : physicalInterfaceResidual(
                system, rhs, reduced.temperature, partition.interfaceGlobalDofs);
        const ErrorMetrics pureError = reduced.error;

        if (options.mode == "corrected") {
            const ddm_schur::SolveResult zero = exactSolver->solve(rhs);
            const auto correctionStart = std::chrono::steady_clock::now();
            ddm_schur::SolveResult corrected = exactSolver->solveCorrection(
                rhs, reduced.interfaceTemperature);
            reduced.timing.correctionSeconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - correctionStart).count();
            reduced.timing.totalSeconds += reduced.timing.correctionSeconds;
            reduced.zeroGuessIterations = zero.report.iterations;
            reduced.correctedIterations = corrected.report.iterations;
            reduced.correctedTrueResidual = corrected.report.interfaceRelativeResidual;
            reduced.interfaceTemperature = std::move(corrected.interfaceSolution);
            reduced.temperature = std::move(corrected.temperature);
            reduced.status = corrected.report.status;
            reduced.error = errors(system, rhs, reduced.temperature,
                truth.empty() ? nullptr : &truth);
            reduced.error.interfaceRelativeResidual = corrected.report.interfaceRelativeResidual;
        }

        const auto diagnosticStart = std::chrono::steady_clock::now();
        const DetailedInterfacePhysicsMetrics interfaceMetrics =
            calculateDetailedInterfacePhysicsMetrics(mesh, physics, reduced.temperature);
        const double diagnosticSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - diagnosticStart).count();
        std::ostringstream powerVector;
        for (std::size_t power = 0; power < parameterCase.powersW.size(); ++power) {
            if (power > 0) powerVector << ';';
            powerVector << parameterCase.powersW[power];
        }
        for (const InterfaceTriangleFluxRecord& flux : interfaceMetrics.triangles) {
            fluxOut << caseSlot << ',' << parameterCase.split << ','
                << parameterCase.family << ",\"" << powerVector.str() << "\","
                << flux.interfaceId << ',' << flux.facePairId << ','
                << flux.integrationTriangleId << ','
                << flux.leftSubdomain << ',' << flux.rightSubdomain << ','
                << flux.leftBoundaryEntity << ',' << flux.rightBoundaryEntity << ','
                << flux.area << ',' << flux.temperatureJumpRms << ','
                << flux.leftPhysicalNormalFlux << ',' << flux.rightPhysicalNormalFlux << ','
                << flux.sipgNumericalFlux << ',' << flux.fluxImbalanceL2 << ','
                << flux.relativeFluxImbalance << ',' << interfaceMetrics.relativeFluxFloor << ','
                << flux.maximumFluxImbalance << ','
                << interfaceMetrics.worstIntegrationTriangle << '\n';
        }

        if (!truth.empty()) {
            for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
                const SubdomainModel& local = model.subdomains[slot];
                long double localError = 0.0L;
                long double localTruth = 0.0L;
                double localMaximum = 0.0;
                for (int global : local.interiorGlobalDofs) {
                    const double difference = reduced.temperature[static_cast<std::size_t>(global)]
                        - truth[static_cast<std::size_t>(global)];
                    localError += static_cast<long double>(difference) * difference;
                    localTruth += static_cast<long double>(truth[static_cast<std::size_t>(global)])
                        * truth[static_cast<std::size_t>(global)];
                    localMaximum = std::max(localMaximum, std::abs(difference));
                }
                worstLocalRelative[slot] = std::max(worstLocalRelative[slot],
                    std::sqrt(static_cast<double>(localError))
                    / std::max(1.0e-300, std::sqrt(static_cast<double>(localTruth))));
                worstLocalMaximum[slot] = std::max(
                    worstLocalMaximum[slot], localMaximum);
            }
        }

        EvaluationRow row;
        row.caseIndex = static_cast<int>(caseSlot);
        row.split = parameterCase.split;
        row.family = parameterCase.family;
        row.pureError = pureError;
        row.solve = std::move(reduced);
        row.temperatureJump = interfaceMetrics.aggregate.temperatureJumpRms;
        row.relativeFluxImbalance = interfaceMetrics.aggregate.relativeFluxImbalance;
        row.fluxImbalanceL2 = interfaceMetrics.areaWeightedFluxImbalanceL2;
        row.maximumFluxImbalance = interfaceMetrics.maximumFluxImbalance;
        row.worstInterfaceId = interfaceMetrics.worstInterfaceId;
        row.worstFacePairId = interfaceMetrics.worstFacePairId;
        row.worstIntegrationTriangle = interfaceMetrics.worstIntegrationTriangle;
        row.fomSeconds = fomSeconds;
        row.sourceProjectionSeconds = sourceSeconds;
        row.diagnosticSeconds = diagnosticSeconds;
        if (options.compareFom) {
            accepted = accepted
                && row.solve.error.relativeL2 < 1.0e-5
                && row.solve.error.maximumAbsolute < 0.1
                && row.solve.error.maximumTemperatureError < 0.01;
        }
        accepted = accepted && row.solve.status == "success";
        rows.push_back(std::move(row));
    }

    if (exactSolver) {
        exactSchurSetupSeconds = exactSolver->setupReport().setupSeconds;
        model.offlineSeconds += exactSchurSetupSeconds;
    }

    writeRankTable(model, worstLocalRelative, worstLocalMaximum,
        outputDirectory / "local_rom_rank_by_subdomain.csv");
    writeAccuracy(rows, outputDirectory / "local_rom_accuracy_by_case.csv");
    writeOnlineTiming(rows, outputDirectory / "local_rom_online_timing.csv");

    const std::vector<ParameterCase> deploymentCases = makeParameterCases(
        sources, "deployment", 100, options.seed + 10007ULL, false);
    std::vector<MultiRhsTimingRow> multiRhsRows;
    const std::vector<int> deploymentBatchSizes = options.mode == "pure"
        ? std::vector<int>{1, 10, 100}
        : std::vector<int>{};
    for (int rhsCount : deploymentBatchSizes) {
        MultiRhsTimingRow aggregate;
        aggregate.rhsCount = rhsCount;
        for (int sample = 0; sample < rhsCount; ++sample) {
            const auto sourceStart = std::chrono::steady_clock::now();
            const std::vector<double> rhs = composeRhs(
                sources, deploymentCases[static_cast<std::size_t>(sample)].powersW);
            aggregate.sourceProjectionSeconds += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - sourceStart).count();
            const SolveResult batchResult = online.solve(rhs);
            aggregate.reducedAssemblySeconds +=
                batchResult.timing.localReducedAssemblySeconds;
            aggregate.interfaceSolveSeconds += batchResult.timing.interfaceSolveSeconds;
            aggregate.localBackSubstitutionSeconds += batchResult.timing.localRecoverySeconds;
            aggregate.reconstructionSeconds +=
                batchResult.timing.fullFieldReconstructionSeconds;
            aggregate.totalOnlineSeconds += batchResult.timing.totalSeconds;
            if (!batchResult.temperature.empty()) {
                aggregate.checksum += batchResult.temperature.front()
                    + batchResult.temperature.back();
            }
        }
        multiRhsRows.push_back(aggregate);
    }
    writeMultiRhsTiming(
        multiRhsRows, outputDirectory / "local_rom_multi_rhs_timing.csv");

    std::ofstream offline(outputDirectory / "local_rom_offline_timing.csv");
    offline << "model_load_seconds,interface_partition_seconds,fom_symbolic_analysis_seconds,"
        << "fom_numerical_factorization_seconds,fom_factorization_wall_seconds,"
        << "snapshot_solve_seconds,local_basis_seconds,local_projection_seconds,"
        << "reduced_schur_construction_seconds,interface_solver,interface_nnz,"
        << "coarse_dimension,proxy_colors,proxy_probing_applies,proxy_setup_seconds,"
        << "interface_symbolic_analysis_seconds,interface_numerical_factorization_seconds,"
        << "reduced_factorization_seconds,exact_schur_setup_seconds,serialization_seconds,"
        << "total_offline_seconds,peak_working_set_bytes\n" << std::setprecision(17)
        << modelLoadSeconds << ',' << partitionSeconds << ','
        << (fomSolver ? fomSolver->symbolicAnalysisSeconds() : 0.0) << ','
        << (fomSolver ? fomSolver->numericalFactorizationSeconds() : 0.0) << ','
        << fomFactorWallSeconds << ','
        << model.snapshotSeconds << ',' << model.basisSeconds << ','
        << model.projectionSeconds << ',' << model.schurConstructionSeconds << ','
        << online.interfaceSolver() << ',' << online.interfaceNonzeros() << ','
        << online.coarseDimension() << ',' << online.proxyColors() << ','
        << online.proxyProbingApplies() << ',' << online.proxySetupSeconds() << ','
        << online.symbolicAnalysisSeconds() << ','
        << online.numericalFactorizationSeconds() << ','
        << online.factorizationSeconds() << ',' << exactSchurSetupSeconds << ','
        << serializationSeconds << ','
        << model.offlineSeconds << ','
        << peakWorkingSetBytes() << '\n';

    std::ofstream memory(outputDirectory / "local_rom_memory.csv");
    memory << "model_bytes,fom_factor_bytes,interface_factor_bytes,proxy_memory_bytes,peak_working_set_bytes,global_dofs,"
        << "interface_dofs,total_local_rank,unique_template_count,reused_instance_count,"
        << "basis_storage_without_reuse_bytes,basis_storage_with_reuse_bytes\n"
        << model.modelBytes << ',' << (fomSolver ? fomSolver->memoryBytes() : 0) << ','
        << online.factorMemoryBytes() << ',' << online.proxyMemoryBytes() << ','
        << peakWorkingSetBytes() << ',' << model.globalDofs << ','
        << model.interfaceDofs << ',' << model.totalLocalRank << ','
        << templateStorage.uniqueTemplates << ',' << templateStorage.reusedInstances << ','
        << templateStorage.withoutReuseBytes << ',' << templateStorage.withReuseBytes << '\n';

    const EvaluationRow& nominalRow = rows.front();
    const bool milestone3 = model.globalDofs > 100000;
    const bool milestone2 = !milestone3 && model.subdomains.size() > 2;
    const std::string milestoneName = milestone3
        ? "Milestone 3" : (milestone2 ? "Milestone 2" : "Milestone 1");
    const std::string milestoneToken = milestone3
        ? "milestone3" : (milestone2 ? "milestone2" : "milestone1");
    std::ofstream summary(outputDirectory / "local_rom_schur_summary.csv");
    summary << "method,global_dofs,subdomains,interface_dofs,total_local_rank,interface_solver,interface_nnz,"
        << "coarse_dimension,proxy_colors,proxy_probing_applies,proxy_setup_seconds,"
        << "unique_template_count,reused_instance_count,"
        << "training_snapshots,mode,offline_seconds,nominal_online_seconds,"
        << "relative_l2,max_node_error_k,max_temperature_error_k,"
        << "interface_temperature_jump_rms_k,relative_flux_imbalance,"
        << "worst_interface_id,worst_face_pair_id,worst_integration_triangle,relative_flux_floor_w_m2,"
        << "interface_fgmres_iterations,interface_fgmres_matvecs,interface_fgmres_true_residual,"
        << "zero_guess_iterations,local_rom_guess_iterations,initial_true_residual,"
        << "final_true_residual,correction_seconds,"
        << "model_bytes,status\n" << std::setprecision(17)
        << "Local-POD-Schur-ROM," << model.globalDofs << ',' << model.subdomains.size() << ','
        << model.interfaceDofs << ',' << model.totalLocalRank << ','
        << online.interfaceSolver() << ',' << online.interfaceNonzeros() << ','
        << online.coarseDimension() << ',' << online.proxyColors() << ','
        << online.proxyProbingApplies() << ',' << online.proxySetupSeconds() << ','
        << templateStorage.uniqueTemplates << ',' << templateStorage.reusedInstances << ','
        << trainingSnapshotCount << ',' << options.mode << ',' << model.offlineSeconds << ','
        << nominalRow.solve.timing.totalSeconds << ',' << nominalRow.solve.error.relativeL2 << ','
        << nominalRow.solve.error.maximumAbsolute << ','
        << nominalRow.solve.error.maximumTemperatureError << ','
        << nominalRow.temperatureJump << ',' << nominalRow.relativeFluxImbalance << ','
        << nominalRow.worstInterfaceId << ',' << nominalRow.worstFacePairId << ','
        << nominalRow.worstIntegrationTriangle << ",1e-12,"
        << nominalRow.solve.timing.interfaceIterations << ','
        << nominalRow.solve.timing.interfaceMatvecs << ','
        << nominalRow.solve.timing.interfaceRelativeResidual << ','
        << nominalRow.solve.zeroGuessIterations << ','
        << nominalRow.solve.correctedIterations << ','
        << nominalRow.pureError.interfaceRelativeResidual << ','
        << nominalRow.solve.correctedTrueResidual << ','
        << nominalRow.solve.timing.correctionSeconds << ',' << model.modelBytes << ','
        << (accepted ? "success" : "accuracy_gate_failed") << '\n';

    std::ofstream comparison(outputDirectory / "local_rom_vs_global_rom.csv");
    comparison << "method,basis_scope,interface_scope,status\n"
        << "Local-POD-Schur-ROM,independent_per_subdomain,full_physical_interface,"
        << (accepted ? "success" : "accuracy_gate_failed") << '\n'
        << "Global-Reduced-Schur-ROM,global_interface_basis,reduced,"
        << "preserved_benchmark_not_executed_in_" << milestoneToken << "\n"
        << "Global-Block-Arnoldi-ROM,global_volume_basis,global,"
        << "preserved_transient_benchmark_not_executed_in_" << milestoneToken << "\n";

    std::ofstream report(outputDirectory / "local_rom_schur_report.md");
    report << "# Local POD/RB + full-interface Schur ROM (" << milestoneName << ")\n\n"
        << "This implementation uses one independently trained interior basis per physical "
        << "subdomain. It does not slice a global basis. The original SIPG matrix supplies "
        << "all A_II, A_IΓ, A_ΓI, and A_ΓΓ blocks; interface DOFs remain full order.\n\n"
        << "- Global DOFs: " << model.globalDofs << "\n"
        << "- Subdomains: " << model.subdomains.size() << "\n"
        << "- Full interface DOFs: " << model.interfaceDofs << "\n"
        << "- Total local rank: " << model.totalLocalRank << "\n"
        << "- Offline time: " << model.offlineSeconds << " s\n"
        << "- Nominal online time: " << nominalRow.solve.timing.totalSeconds << " s\n"
        << "- Nominal relative L2: " << nominalRow.solve.error.relativeL2 << "\n"
        << "- Nominal max-node error: " << nominalRow.solve.error.maximumAbsolute << " K\n"
        << "- Status: " << (accepted ? "success" : "accuracy_gate_failed") << "\n";

    double worstEvaluationL2 = 0.0;
    double worstEvaluationMaximum = 0.0;
    double worstEvaluationHotspot = 0.0;
    double worstEvaluationJump = 0.0;
    double worstEvaluationFlux = 0.0;
    int worstEvaluationCase = -1;
    std::string worstEvaluationSplit;
    std::string worstEvaluationFamily;
    for (const EvaluationRow& row : rows) {
        if (row.solve.error.relativeL2 >= worstEvaluationL2) {
            worstEvaluationL2 = row.solve.error.relativeL2;
            worstEvaluationCase = row.caseIndex;
            worstEvaluationSplit = row.split;
            worstEvaluationFamily = row.family;
        }
        worstEvaluationMaximum = std::max(
            worstEvaluationMaximum, row.solve.error.maximumAbsolute);
        worstEvaluationHotspot = std::max(
            worstEvaluationHotspot, row.solve.error.maximumTemperatureError);
        worstEvaluationJump = std::max(worstEvaluationJump, row.temperatureJump);
        worstEvaluationFlux = std::max(
            worstEvaluationFlux, row.maximumFluxImbalance);
    }

    // Replace the short console-style summary above with the complete milestone
    // technical report. Keep this generated so every reproduction run is auditable.
    report.close();
    std::ofstream detailedReport(outputDirectory / "local_rom_schur_report.md");
    detailedReport << std::setprecision(12)
        << "# Local POD/RB + full-interface Schur ROM (" << milestoneName << ")\n\n"
        << "This implementation uses one independently trained interior basis per physical "
        << "subdomain. It does not slice a global basis. The original SIPG matrix supplies "
        << "all A_II, A_I_Gamma, A_Gamma_I, and A_Gamma_Gamma blocks; interface DOFs "
        << "remain full order. Reduced local Schur contributions are assembled in the "
        << "original global interface ordering.\n\n"
        << "## Configuration\n\n"
        << "| global DOFs | subdomains | full interface DOFs | total local rank | unique templates | reused instances | interface solver | mode |\n"
        << "|---:|---:|---:|---:|---:|---:|---|---|\n"
        << "| " << model.globalDofs << " | " << model.subdomains.size() << " | "
        << model.interfaceDofs << " | " << model.totalLocalRank << " | "
        << templateStorage.uniqueTemplates << " | " << templateStorage.reusedInstances << " | "
        << online.interfaceSolver() << " | "
        << options.mode << " |\n\n"
        << "| subdomain | interior DOFs | local interface DOFs | snapshots | numerical rank | selected rank | orthogonality error | retained energy |\n"
        << "|---:|---:|---:|---:|---:|---:|---:|---:|\n";
    for (const SubdomainModel& local : model.subdomains) {
        const double energy = local.rank > 0
            && static_cast<std::size_t>(local.rank) <= local.retainedEnergy.size()
            ? local.retainedEnergy[static_cast<std::size_t>(local.rank - 1)] : 0.0;
        detailedReport << "| " << local.subdomain << " | " << local.interiorDofs << " | "
            << local.localInterfaceDofs << " | " << local.snapshots << " | "
            << local.numericalRank << " | " << local.rank << " | "
            << local.orthogonalityError << " | " << energy << " |\n";
    }
    detailedReport << "\n## Accuracy and interface physics\n\n"
        << "| relative L2 | max-node error (K) | max-temperature error (K) | interface jump RMS (K) | area-weighted flux imbalance L2 (W/m2) | max flux imbalance (W/m2) | interface residual |\n"
        << "|---:|---:|---:|---:|---:|---:|---:|\n"
        << "| " << nominalRow.solve.error.relativeL2 << " | "
        << nominalRow.solve.error.maximumAbsolute << " | "
        << nominalRow.solve.error.maximumTemperatureError << " | "
        << nominalRow.temperatureJump << " | " << nominalRow.fluxImbalanceL2 << " | "
        << nominalRow.maximumFluxImbalance << " | "
        << nominalRow.solve.error.interfaceRelativeResidual << " |\n\n"
        << "Worst over all nominal/validation/test cases: relative L2 "
        << worstEvaluationL2 << " (case " << worstEvaluationCase << ", "
        << worstEvaluationSplit << ", " << worstEvaluationFamily << "), max-node error "
        << worstEvaluationMaximum << " K, max-temperature error "
        << worstEvaluationHotspot << " K, interface jump RMS "
        << worstEvaluationJump << " K, and maximum flux imbalance "
        << worstEvaluationFlux << " W/m2. Nominal worst physical-interface/face-pair/triangle IDs are "
        << nominalRow.worstInterfaceId << "/" << nominalRow.worstFacePairId << "/"
        << nominalRow.worstIntegrationTriangle << ".\n\n"
        << "The interface diagnostics integrate physical normal flux and the original SIPG "
        << "numerical flux on the existing nonmatching triangle-overlap/BVH faces. No nodal "
        << "one-to-one interface was introduced. The relative flux denominator is "
        << "max(||q_left||, ||q_right||, 1e-12 W/m2 * sqrt(interface area)); absolute L2 "
        << "and maximum imbalance remain explicit acceptance diagnostics.\n\n"
        << "## Timing and memory\n\n"
        << "| FOM factor wall (s) | snapshot solves (s) | local bases (s) | projection (s) | reduced Schur build (s) | reduced factorization (s) | total offline (s) | nominal online (s) |\n"
        << "|---:|---:|---:|---:|---:|---:|---:|---:|\n"
        << "| " << fomFactorWallSeconds << " | " << model.snapshotSeconds << " | "
        << model.basisSeconds << " | " << model.projectionSeconds << " | "
        << model.schurConstructionSeconds << " | " << online.factorizationSeconds()
        << " | " << model.offlineSeconds << " | "
        << nominalRow.solve.timing.totalSeconds << " |\n\n"
        << "| model bytes | FOM factor bytes | interface factor bytes | peak working-set bytes |\n"
        << "|---:|---:|---:|---:|\n"
        << "| " << model.modelBytes << " | "
        << (fomSolver ? fomSolver->memoryBytes() : 0) << " | "
        << online.factorMemoryBytes() << " | "
        << peakWorkingSetBytes() << " |\n\n"
        << "### Fixed-matrix multi-RHS deployment\n\n"
        << "| RHS count | source projection total (s) | core online total (s) | average online (s) |\n"
        << "|---:|---:|---:|---:|\n";
    for (const MultiRhsTimingRow& timing : multiRhsRows) {
        detailedReport << "| " << timing.rhsCount << " | "
            << timing.sourceProjectionSeconds << " | " << timing.totalOnlineSeconds
            << " | " << timing.totalOnlineSeconds / std::max(1, timing.rhsCount) << " |\n";
    }
    detailedReport << "\n"
        << "## Milestone decision\n\n"
        << "- Accuracy gate: full-field relative L2 < 1e-5, max-node error < 0.1 K, "
        << "and max-temperature error < 0.01 K.\n"
        << "- Result: **" << (accepted ? "pass" : "fail") << "**.\n"
        << (milestone3
            ? "- Scope: this run is the RRAM26/Chiplet steady Milestone 3 Local-ROM + full-interface Schur path. Local Block Arnoldi and dynamic Schur remain out of scope until the next milestone.\n"
            : (milestone2
                ? "- Scope: this commit completes ten-cube steady Milestone 2 with scalable local bases and the full physical Schur interface. Local Block Arnoldi and dynamic Schur are explicitly out of scope.\n"
                : "- Scope: this run is the preserved two-cube steady Milestone 1 regression.\n"))
        << "- The existing Global Reduced Schur and Global-Block-Arnoldi-ROM remain unchanged "
        << "benchmark/regression paths, not the Local-ROM-Schur mainline.\n";

    WorkflowResult workflow;
    workflow.nominalTemperature = nominalRow.solve.temperature;
    workflow.subdomains = static_cast<int>(model.subdomains.size());
    workflow.interfaceDofs = model.interfaceDofs;
    workflow.totalLocalRank = model.totalLocalRank;
    workflow.correctionIterations = nominalRow.solve.correctedIterations >= 0
        ? nominalRow.solve.correctedIterations
        : nominalRow.solve.timing.interfaceIterations;
    workflow.coarseDimension = online.coarseDimension();
    workflow.setupSeconds = model.offlineSeconds;
    workflow.onlineSeconds = nominalRow.solve.timing.totalSeconds;
    workflow.coarseSolveSeconds = nominalRow.solve.timing.coarseSolveSeconds;
    workflow.globalResidual = nominalRow.solve.error.globalRelativeResidual;
    workflow.interfaceResidual = nominalRow.solve.error.interfaceRelativeResidual;
    workflow.status = accepted ? "success" : "accuracy_gate_failed";

    std::cout << std::setprecision(12)
        << "[Local ROM] subdomains/interface DOFs/total rank: "
        << workflow.subdomains << '/' << workflow.interfaceDofs << '/'
        << workflow.totalLocalRank << '\n'
        << "[Local ROM] offline/nominal online: " << workflow.setupSeconds
        << " / " << workflow.onlineSeconds << " s\n"
        << "[Local ROM] nominal relative L2/max error: "
        << nominalRow.solve.error.relativeL2 << " / "
        << nominalRow.solve.error.maximumAbsolute << " K\n";
    return workflow;
}

} // namespace mor::local
