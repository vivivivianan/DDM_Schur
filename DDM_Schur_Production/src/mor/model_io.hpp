#pragma once

// Cache serialization API for production local models and fingerprints.

#include "types.hpp"

#include <filesystem>
#include <vector>

struct CaseConfig;
struct Mesh;
struct SparseMatrix;

namespace mor {

Fingerprints computeFingerprints(const Mesh& mesh,
                                 const SparseMatrix& system,
                                 const CaseConfig& physics,
                                 const std::vector<int>& interfaceGlobalDofs);

void saveModel(const ReducedSchurModel& model,
               const std::filesystem::path& directory);

ReducedSchurModel loadModel(const std::filesystem::path& directory,
                            const Fingerprints& expected,
                            int expectedInterfaceDofs,
                            int expectedSourceChannels);

} // namespace mor
