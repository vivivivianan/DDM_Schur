#pragma once

#include "types.hpp"

#include <filesystem>
#include <string>

namespace mor {

void saveDeploymentResponseModel(DeploymentResponseModel& model,
                                 const std::filesystem::path& modelDirectory,
                                 const std::string& interiorMode,
                                 const std::string& storagePrecision);

DeploymentResponseModel loadDeploymentResponseModel(
    const std::filesystem::path& modelDirectory,
    const std::string& requestedInteriorMode,
    const Fingerprints* expectedFingerprints = nullptr,
    int expectedGlobalDofs = 0,
    int expectedInterfaceDofs = 0,
    int expectedInterfaceRank = 0,
    int expectedSourceChannels = 0);

std::uint64_t recursiveFileBytes(const std::filesystem::path& directory);

} // namespace mor
