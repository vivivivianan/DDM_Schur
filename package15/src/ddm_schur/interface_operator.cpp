#include "../sipg_core.hpp"
#include "interface_operator.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>

namespace ddm_schur {

InterfacePartition buildInterfacePartition(const Mesh& mesh, const SparseMatrix& system)
{
    if (system.size() != static_cast<int>(mesh.nodes.size())) {
        throw std::runtime_error("[Schur] Matrix size and mesh DOF count do not match.");
    }

    InterfacePartition result;
    result.totalDofs = system.size();
    std::set<int> domainIds;
    for (const Node& node : mesh.nodes) {
        domainIds.insert(node.subdomain);
    }
    std::map<int, int> domainSlot;
    for (int id : domainIds) {
        domainSlot[id] = static_cast<int>(result.domains.size());
        DomainBlocks domain;
        domain.domainId = id;
        result.domains.push_back(std::move(domain));
    }

    std::vector<char> isInterface(static_cast<std::size_t>(system.size()), 0);
    system.forEachEntry([&](int row, int col, double) {
        if (row < 0 || row >= system.size() || col < 0 || col >= system.size()) {
            return;
        }
        if (mesh.nodes[static_cast<std::size_t>(row)].subdomain
            != mesh.nodes[static_cast<std::size_t>(col)].subdomain) {
            isInterface[static_cast<std::size_t>(row)] = 1;
            isInterface[static_cast<std::size_t>(col)] = 1;
        }
    });

    result.globalToInterface.assign(static_cast<std::size_t>(system.size()), -1);
    std::vector<int> globalToInterior(static_cast<std::size_t>(system.size()), -1);
    std::vector<int> globalToDomainInterface(static_cast<std::size_t>(system.size()), -1);
    for (int dof = 0; dof < system.size(); ++dof) {
        const int id = mesh.nodes[static_cast<std::size_t>(dof)].subdomain;
        DomainBlocks& domain = result.domains[static_cast<std::size_t>(domainSlot.at(id))];
        if (isInterface[static_cast<std::size_t>(dof)]) {
            result.globalToInterface[static_cast<std::size_t>(dof)] =
                static_cast<int>(result.interfaceGlobalDofs.size());
            result.interfaceGlobalDofs.push_back(dof);
            globalToDomainInterface[static_cast<std::size_t>(dof)] =
                static_cast<int>(domain.interfaceGlobalDofs.size());
            domain.interfaceGlobalDofs.push_back(dof);
        } else {
            globalToInterior[static_cast<std::size_t>(dof)] =
                static_cast<int>(domain.interiorGlobalDofs.size());
            domain.interiorGlobalDofs.push_back(dof);
        }
    }

    for (DomainBlocks& domain : result.domains) {
        domain.interiorInterfaceRows.resize(domain.interiorGlobalDofs.size());
        domain.interfaceInteriorRows.resize(domain.interfaceGlobalDofs.size());
    }

    system.forEachEntry([&](int row, int col, double value) {
        if (row < 0 || row >= system.size() || col < 0 || col >= system.size()) {
            return;
        }
        const bool rowGamma = isInterface[static_cast<std::size_t>(row)] != 0;
        const bool colGamma = isInterface[static_cast<std::size_t>(col)] != 0;
        const int rowId = mesh.nodes[static_cast<std::size_t>(row)].subdomain;
        const int colId = mesh.nodes[static_cast<std::size_t>(col)].subdomain;

        if (!rowGamma && !colGamma) {
            if (rowId != colId) {
                throw std::runtime_error("[Schur] Cross-domain coupling escaped interface detection.");
            }
            DomainBlocks& domain = result.domains[static_cast<std::size_t>(domainSlot.at(rowId))];
            const int localRow = globalToInterior[static_cast<std::size_t>(row)];
            const int localCol = globalToInterior[static_cast<std::size_t>(col)];
            domain.interiorEntries.push_back({localRow, localCol, value});
            domain.fullBlockEntries.push_back({localRow, localCol, value});
            return;
        }

        if (!rowGamma && colGamma) {
            DomainBlocks& domain = result.domains[static_cast<std::size_t>(domainSlot.at(rowId))];
            if (rowId != colId) {
                throw std::runtime_error("[Schur] Interior DOF is coupled to a foreign interface owner.");
            }
            const int localRow = globalToInterior[static_cast<std::size_t>(row)];
            const int gammaCol = result.globalToInterface[static_cast<std::size_t>(col)];
            const int localGammaCol = globalToDomainInterface[static_cast<std::size_t>(col)];
            domain.interiorInterfaceRows[static_cast<std::size_t>(localRow)].push_back({gammaCol, value});
            domain.fullBlockEntries.push_back(
                {localRow, static_cast<int>(domain.interiorGlobalDofs.size()) + localGammaCol, value});
            return;
        }

        if (rowGamma && !colGamma) {
            DomainBlocks& domain = result.domains[static_cast<std::size_t>(domainSlot.at(colId))];
            if (rowId != colId) {
                throw std::runtime_error("[Schur] Foreign interface row is coupled to an interior DOF.");
            }
            const int localGammaRow = globalToDomainInterface[static_cast<std::size_t>(row)];
            const int localCol = globalToInterior[static_cast<std::size_t>(col)];
            domain.interfaceInteriorRows[static_cast<std::size_t>(localGammaRow)].push_back({localCol, value});
            domain.fullBlockEntries.push_back(
                {static_cast<int>(domain.interiorGlobalDofs.size()) + localGammaRow, localCol, value});
            return;
        }

        const int gammaRow = result.globalToInterface[static_cast<std::size_t>(row)];
        const int gammaCol = result.globalToInterface[static_cast<std::size_t>(col)];
        result.interfaceEntries.push_back({gammaRow, gammaCol, value});
        if (rowId != colId) {
            result.domains[static_cast<std::size_t>(domainSlot.at(rowId))]
                .interfaceGlobalDofsByNeighbor[colId].push_back(row);
            result.domains[static_cast<std::size_t>(domainSlot.at(colId))]
                .interfaceGlobalDofsByNeighbor[rowId].push_back(col);
        }
        if (rowId == colId) {
            DomainBlocks& domain = result.domains[static_cast<std::size_t>(domainSlot.at(rowId))];
            const int offset = static_cast<int>(domain.interiorGlobalDofs.size());
            domain.fullBlockEntries.push_back(
                {offset + globalToDomainInterface[static_cast<std::size_t>(row)],
                 offset + globalToDomainInterface[static_cast<std::size_t>(col)],
                 value});
        }
    });

    for (DomainBlocks& domain : result.domains) {
        for (auto& neighborDofs : domain.interfaceGlobalDofsByNeighbor) {
            std::vector<int>& dofs = neighborDofs.second;
            std::sort(dofs.begin(), dofs.end());
            dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
        }
    }

    return result;
}

std::size_t partitionMemoryBytes(const InterfacePartition& partition)
{
    std::size_t bytes = partition.globalToInterface.size() * sizeof(int)
        + partition.interfaceGlobalDofs.size() * sizeof(int)
        + partition.interfaceEntries.size() * sizeof(Entry);
    for (const DomainBlocks& domain : partition.domains) {
        bytes += domain.interiorGlobalDofs.size() * sizeof(int);
        bytes += domain.interfaceGlobalDofs.size() * sizeof(int);
        for (const auto& neighborDofs : domain.interfaceGlobalDofsByNeighbor) {
            bytes += sizeof(neighborDofs.first) + neighborDofs.second.size() * sizeof(int);
        }
        bytes += domain.interiorEntries.size() * sizeof(Entry);
        bytes += domain.fullBlockEntries.size() * sizeof(Entry);
        for (const auto& row : domain.interiorInterfaceRows) {
            bytes += row.size() * sizeof(std::pair<int, double>);
        }
        for (const auto& row : domain.interfaceInteriorRows) {
            bytes += row.size() * sizeof(std::pair<int, double>);
        }
    }
    return bytes;
}

} // namespace ddm_schur
