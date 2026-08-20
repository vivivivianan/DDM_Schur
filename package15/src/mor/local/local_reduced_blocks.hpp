#pragma once

#include "local_subdomain_model.hpp"

namespace ddm_schur { struct InterfacePartition; }

namespace mor::local {

void projectLocalReducedBlocks(Model& model,
                               const ddm_schur::InterfacePartition& partition);

} // namespace mor::local
