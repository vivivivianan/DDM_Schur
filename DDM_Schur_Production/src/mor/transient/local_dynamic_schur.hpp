#pragma once

// Single production workflow entry point. Inputs are an assembled case mesh,
// physical configuration, fixed options, and an empty output directory.

#include "transient_workflow.hpp"

#include <filesystem>

struct CaseConfig;
struct Mesh;

namespace mor::transient {

// Stage 2 local transient mainline.  The interface remains in the original
// physical ordering; only each subdomain interior is reduced by an
// independently generated Block Arnoldi basis.
void runLocalBlockArnoldiDynamicSchurWorkflow(
    const Mesh& mesh,
    const CaseConfig& physics,
    const Options& options,
    const std::filesystem::path& outputDirectory);

} // namespace mor::transient
