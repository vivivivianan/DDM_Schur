#pragma once

#include "transient_reduced_model.hpp"

#include <filesystem>

namespace mor::transient {

void saveTransientReducedModel(
    const TransientReducedModel& model,
    const std::filesystem::path& directory);

TransientReducedModel loadTransientReducedModel(
    const std::filesystem::path& directory,
    const TransientFingerprints* expectedFingerprints = nullptr);

std::uint64_t transientModelFileBytes(
    const std::filesystem::path& directory);

} // namespace mor::transient
