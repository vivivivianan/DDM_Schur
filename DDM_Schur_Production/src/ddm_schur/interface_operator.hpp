#pragma once

// Domain-local matrix blocks plus the one-to-one maps between global and
// physical-interface DOFs. The production solver assembles local contributions
// in parallel but uses this common interface ordering for the coupled solve.

#include "common.hpp"

#include <map>
#include <utility>
#include <vector>

namespace ddm_schur {

struct Entry {
    int row = 0;
    int col = 0;
    double value = 0.0;
};

struct DomainBlocks {
    int domainId = -1;
    std::vector<int> interiorGlobalDofs;
    std::vector<int> interfaceGlobalDofs;
    std::map<int, std::vector<int>> interfaceGlobalDofsByNeighbor;
    std::vector<Entry> interiorEntries;
    std::vector<Entry> fullBlockEntries;
    std::vector<std::vector<std::pair<int, double>>> interiorInterfaceRows;
    std::vector<std::vector<std::pair<int, double>>> interfaceInteriorRows;
};

struct InterfacePartition {
    int totalDofs = 0;
    std::vector<int> globalToInterface;
    std::vector<int> interfaceGlobalDofs;
    std::vector<Entry> interfaceEntries;
    std::vector<DomainBlocks> domains;
};

InterfacePartition buildInterfacePartition(const Mesh& mesh, const SparseMatrix& system);
std::size_t partitionMemoryBytes(const InterfacePartition& partition);

} // namespace ddm_schur
