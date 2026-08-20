#pragma once

#include "types.hpp"

namespace mor {

void buildCompressedInteriorBases(DeploymentResponseModel& model,
                                  int requestedRank,
                                  double discardedEnergyTolerance,
                                  double relativeSingularTolerance);

} // namespace mor
