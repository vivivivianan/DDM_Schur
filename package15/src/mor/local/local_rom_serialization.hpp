#pragma once

#include "local_subdomain_model.hpp"

#include <filesystem>

namespace mor::local {

void saveLocalRomModel(const Model& model, const std::filesystem::path& directory);
Model loadLocalRomModel(const std::filesystem::path& directory,
                        const Fingerprints& expectedFingerprints,
                        int expectedGlobalDofs,
                        int expectedInterfaceDofs);

} // namespace mor::local
