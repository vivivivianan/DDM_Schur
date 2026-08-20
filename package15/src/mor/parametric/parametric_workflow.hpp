#pragma once

#include "mor/types.hpp"

#include <filesystem>
#include <vector>

struct CaseConfig;
struct Mesh;
struct SparseMatrix;
namespace ddm_schur { struct Options; }

namespace mor::parametric {

WorkflowResult runParametricReducedSchurWorkflow(
    const Mesh& mesh,
    const SparseMatrix& referenceSystem,
    const CaseConfig& physics,
    const std::vector<double>& assembledSource,
    const std::vector<double>& heatOnlySource,
    const std::vector<double>& fixedAdjust,
    const ddm_schur::Options& schurOptions,
    const mor::Options& options,
    const std::filesystem::path& outputDirectory);

void runParametricDeploymentOnly(
    const std::filesystem::path& modelDirectory,
    const std::filesystem::path& outputDirectory,
    const mor::Options& options);

} // namespace mor::parametric
