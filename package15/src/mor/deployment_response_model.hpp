#pragma once

#include "types.hpp"

#include <vector>

struct Mesh;
namespace ddm_schur { class DdmSchurSolver; }

namespace mor {

struct SourceParameterization;
class ReducedSchurOnlineSolver;

DeploymentResponseModel buildDeploymentResponseModel(
    const Mesh& mesh,
    const SourceParameterization& sources,
    const std::vector<int>& interfaceGlobalDofs,
    const ReducedSchurModel& interfaceModel,
    ReducedSchurOnlineSolver& onlineSolver);

SolveResult evaluateDeploymentResponse(
    const DeploymentResponseModel& model,
    const std::vector<double>& powersW,
    const std::string& interiorMode,
    int rankOverride = 0);

DeploymentRunResult runDeploymentOnly(
    const std::filesystem::path& modelDirectory,
    const std::filesystem::path& outputDirectory,
    const Options& options);

} // namespace mor
