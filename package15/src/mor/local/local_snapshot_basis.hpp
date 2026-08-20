#pragma once

#include "local_subdomain_model.hpp"

#include <vector>

namespace ddm_schur { struct InterfacePartition; }

namespace mor::local {

Model buildIndependentSnapshotBases(
    const ddm_schur::InterfacePartition& partition,
    const std::vector<double>& referenceTemperature,
    const std::vector<std::vector<double>>& trainingTemperatures,
    const Options& options);

Model buildIndependentSnapshotBasesFromLocalDatabases(
    const ddm_schur::InterfacePartition& partition,
    const std::vector<double>& referenceTemperature,
    std::vector<::mor::SnapshotDatabase> localSnapshots,
    const Options& options);

} // namespace mor::local
