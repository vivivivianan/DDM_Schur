#include "local_snapshot_basis.hpp"

#include "ddm_schur/interface_operator.hpp"
#include "mor/pod_basis.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace mor::local {
namespace {

double orthogonalityError(const std::vector<double>& basis, int rows, int rank)
{
    double maximum = 0.0;
    for (int left = 0; left < rank; ++left) {
        for (int right = 0; right < rank; ++right) {
            double product = 0.0;
            for (int row = 0; row < rows; ++row) {
                product += basis[static_cast<std::size_t>(left * rows + row)]
                    * basis[static_cast<std::size_t>(right * rows + row)];
            }
            maximum = std::max(maximum,
                std::abs(product - (left == right ? 1.0 : 0.0)));
        }
    }
    return maximum;
}

int requestedRankFor(const Options& options, std::size_t slot, std::size_t domains)
{
    if (!options.rankPerSubdomain.empty()) {
        if (options.rankPerSubdomain.size() != domains) {
            throw std::runtime_error(
                "[Local ROM] --local-mor-rank-per-subdomain must contain one rank per subdomain.");
        }
        return options.rankPerSubdomain[slot];
    }
    return options.rank;
}

} // namespace

Model buildIndependentSnapshotBases(
    const ddm_schur::InterfacePartition& partition,
    const std::vector<double>& referenceTemperature,
    const std::vector<std::vector<double>>& trainingTemperatures,
    const Options& options)
{
    if (referenceTemperature.size() != static_cast<std::size_t>(partition.totalDofs)
        || trainingTemperatures.empty()) {
        throw std::runtime_error("[Local ROM] Invalid local snapshot dimensions.");
    }
    for (const auto& temperature : trainingTemperatures) {
        if (temperature.size() != referenceTemperature.size()) {
            throw std::runtime_error("[Local ROM] Training snapshots have inconsistent sizes.");
        }
    }

    std::vector<SnapshotDatabase> localSnapshots(partition.domains.size());
    for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
        const ddm_schur::DomainBlocks& domain = partition.domains[slot];
        SnapshotDatabase& snapshots = localSnapshots[slot];
        snapshots.rows = static_cast<int>(domain.interiorGlobalDofs.size());
        snapshots.cases.resize(trainingTemperatures.size());
        snapshots.values.assign(static_cast<std::size_t>(snapshots.rows)
            * trainingTemperatures.size(), 0.0);
        for (std::size_t sample = 0; sample < trainingTemperatures.size(); ++sample) {
            snapshots.cases[sample].index = static_cast<int>(sample);
            snapshots.cases[sample].split = "train";
            snapshots.cases[sample].family = "global_power_and_interface_response";
            for (int localRow = 0; localRow < snapshots.rows; ++localRow) {
                const int global = domain.interiorGlobalDofs[
                    static_cast<std::size_t>(localRow)];
                snapshots.values[sample * static_cast<std::size_t>(snapshots.rows)
                    + static_cast<std::size_t>(localRow)] =
                    trainingTemperatures[sample][static_cast<std::size_t>(global)]
                    - referenceTemperature[static_cast<std::size_t>(global)];
            }
        }
    }
    return buildIndependentSnapshotBasesFromLocalDatabases(
        partition, referenceTemperature, std::move(localSnapshots), options);
}

Model buildIndependentSnapshotBasesFromLocalDatabases(
    const ddm_schur::InterfacePartition& partition,
    const std::vector<double>& referenceTemperature,
    std::vector<SnapshotDatabase> localSnapshots,
    const Options& options)
{
    if (referenceTemperature.size() != static_cast<std::size_t>(partition.totalDofs)
        || localSnapshots.size() != partition.domains.size()) {
        throw std::runtime_error(
            "[Local ROM] Invalid local snapshot database dimensions.");
    }

    const auto start = std::chrono::steady_clock::now();
    Model model;
    model.globalDofs = partition.totalDofs;
    model.interfaceDofs = static_cast<int>(partition.interfaceGlobalDofs.size());
    model.interfaceGlobalDofs = partition.interfaceGlobalDofs;
    model.subdomains.reserve(partition.domains.size());

    for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
        const ddm_schur::DomainBlocks& domain = partition.domains[slot];
        const int rows = static_cast<int>(domain.interiorGlobalDofs.size());
        if (rows <= 0) {
            throw std::runtime_error("[Local ROM] A subdomain has no interior DOFs.");
        }

        SnapshotDatabase& snapshots = localSnapshots[slot];
        if (snapshots.rows != rows || snapshots.cases.empty()
            || snapshots.values.size() != static_cast<std::size_t>(rows)
                * snapshots.cases.size()) {
            throw std::runtime_error(
                "[Local ROM] A local snapshot database is inconsistent.");
        }
        const auto podStart = std::chrono::steady_clock::now();
        const PodResult pod = buildGramPod(
            snapshots,
            requestedRankFor(options, slot, partition.domains.size()),
            options.energyTolerance,
            options.singularValueTolerance);
        const double podSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - podStart).count();

        SubdomainModel local;
        local.subdomain = domain.domainId;
        local.interiorDofs = rows;
        local.localInterfaceDofs = static_cast<int>(domain.interfaceGlobalDofs.size());
        local.snapshots = static_cast<int>(snapshots.cases.size());
        local.numericalRank = pod.numericalRank;
        local.rank = pod.selectedRank;
        local.interiorGlobalDofs = domain.interiorGlobalDofs;
        local.interfaceGlobalDofs = domain.interfaceGlobalDofs;
        local.interfaceIndices.reserve(domain.interfaceGlobalDofs.size());
        for (int global : domain.interfaceGlobalDofs) {
            local.interfaceIndices.push_back(
                partition.globalToInterface[static_cast<std::size_t>(global)]);
        }
        local.referenceInterior.resize(static_cast<std::size_t>(rows), 0.0);
        for (int row = 0; row < rows; ++row) {
            local.referenceInterior[static_cast<std::size_t>(row)] =
                referenceTemperature[static_cast<std::size_t>(
                    domain.interiorGlobalDofs[static_cast<std::size_t>(row)])];
        }
        local.basis = pod.basis;
        local.singularValues = pod.singularValues;
        local.retainedEnergy = pod.retainedEnergy;
        local.orthogonalityError = orthogonalityError(local.basis, rows, local.rank);
        local.snapshotExtractionSeconds = 0.0;
        local.podSeconds = podSeconds;
        model.totalLocalRank += local.rank;
        model.subdomains.push_back(std::move(local));
        snapshots.values.clear();
        snapshots.values.shrink_to_fit();
    }

    model.basisSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    return model;
}

} // namespace mor::local
