#pragma once

// Compatibility declaration for a removed full-interface projection study.

#include "ddm_schur/interface_operator.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

struct CaseConfig;
struct Mesh;

namespace ddm_schur { struct InterfacePartition; }
namespace mor::local { struct Model; }

namespace mor::transient {

struct LocalPortModel;
struct ThermalDescriptorSystem;

struct ProjectionDiagnosisOptions {
    std::vector<int> interfaceIds;
    double timeStep = 0.0;
    double time = 0.0;
    std::vector<double> powers;
    ddm_schur::Options schurOptions;
    std::filesystem::path outputDirectory;
    bool fluxOperatorAudit = false;
};

struct ProjectionDiagnosisSummary {
    int requestedInterfaces = 0;
    int completedInterfaces = 0;
    double fullInterfaceSetupSeconds = 0.0;
    double fullInterfaceSolveSeconds = 0.0;
    double projectionSeconds = 0.0;
    double maximumTargetSolveResidual = 0.0;
    std::size_t peakWorkingSetBytes = 0;
    bool fullFieldRead = false;
    bool snapshotUsed = false;
    bool fomUsedForBasis = false;
    bool podUsed = false;
    bool svdUsed = false;
    std::string status = "not_run";
};

ProjectionDiagnosisSummary runFullInterfaceProjectionDiagnosis(
    const Mesh& mesh,
    const CaseConfig& physics,
    const ddm_schur::InterfacePartition& partition,
    const ThermalDescriptorSystem& descriptor,
    const local::Model& dynamicModel,
    const LocalPortModel& portModel,
    const std::vector<double>& previousTheta,
    const std::vector<double>& referenceTemperature,
    const std::vector<double>& boundaryOffset,
    const ProjectionDiagnosisOptions& options);

} // namespace mor::transient
