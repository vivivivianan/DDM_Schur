#pragma once

#include "types.hpp"

namespace mor {

PodResult buildGramPod(const SnapshotDatabase& snapshots,
                       int requestedRank,
                       double discardedEnergyTolerance,
                       double relativeSingularTolerance);

std::vector<double> symmetricEigenvalues(std::vector<double> matrix, int n);

} // namespace mor
