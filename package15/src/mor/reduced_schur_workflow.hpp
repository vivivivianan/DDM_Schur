#pragma once

#include "types.hpp"

#include <filesystem>
#include <vector>

struct CaseConfig;
struct Mesh;
struct SparseMatrix;
namespace ddm_schur { struct Options; }

namespace mor {

WorkflowResult runReducedSchurWorkflow(
    const Mesh& mesh,
    const SparseMatrix& system,
    const CaseConfig& physics,
    const std::vector<double>& assembledSource,
    const std::vector<double>& heatOnlySource,
    const std::vector<double>& fixedAdjust,
    const ddm_schur::Options& schurOptions,
    const Options& options,
    const std::filesystem::path& outputDirectory);

} // namespace mor
