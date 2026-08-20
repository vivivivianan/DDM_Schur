// Link boundary for research methods intentionally removed from this export.
//
// local_dynamic_schur.cpp is retained from the validated implementation so
// descriptor construction, cache fingerprints, timing columns, and numerical
// ordering stay byte-for-byte compatible with prior runs. That upstream file
// also contains dormant calls to experimental port/coarse/Krylov methods.
// Their large implementations are not shipped here. These small definitions
// satisfy the static C++ references and fail loudly if a future code change
// accidentally makes one reachable from the production CLI.
//
// This is preferable to silent no-op behavior: the executable either follows
// the documented augmented-direct path or terminates with a precise message.

#include "ddm_schur/schur_fgmres.hpp"
#include "mor/transient/global_interface_coarse.hpp"
#include "mor/transient/global_randomized_schur.hpp"
#include "mor/transient/optimal_port_space.hpp"
#include "mor/transient/projection_diagnosis.hpp"
#include "mor/transient/port/randomized_transfer_port.hpp"

#include <stdexcept>
#include <utility>

namespace {

[[noreturn]] void removed(const char* method)
{
    throw std::runtime_error(
        std::string("[Production build] Removed research method reached: ")
        + method);
}

} // namespace

namespace ddm_schur {

// Complete the private type so the deliberately empty unique_ptr can be
// destroyed in this translation unit.
struct DdmSchurSolver::Impl {};

DdmSchurSolver::DdmSchurSolver(
    const Mesh&, const SparseMatrix&, const CaseConfig&, Options)
{
    removed("Stage-1 iterative DDM Schur");
}

DdmSchurSolver::~DdmSchurSolver() = default;

SolveResult DdmSchurSolver::solve(const std::vector<double>&)
{
    removed("Stage-1 iterative DDM Schur solve");
}

void DdmSchurSolver::applyExactSchur(
    const std::vector<double>&, std::vector<double>&)
{
    removed("Stage-1 exact Schur action");
}

} // namespace ddm_schur

namespace mor::transient {

LocalPortModel buildLocalPortModel(
    const Mesh&, const CaseConfig&, const ddm_schur::InterfacePartition&,
    const LocalPortSnapshotFamilies&, const LocalPortOptions&)
{
    removed("snapshot/POD port basis");
}

void saveLocalPortModel(const LocalPortModel&,
                        const std::filesystem::path&)
{
    removed("port-model serialization");
}

LocalPortModel loadLocalPortModel(const std::filesystem::path&)
{
    removed("port-model deserialization");
}

LocalPortReducedSchurSolver::LocalPortReducedSchurSolver(
    const local::Model& dynamicModel,
    const LocalPortModel& portModel)
    : dynamicModel_(dynamicModel), portModel_(portModel)
{
    removed("port-reduced Schur solver");
}

LocalPortSolveResult LocalPortReducedSchurSolver::solve(
    const std::vector<double>&) const
{
    removed("port-reduced Schur solve");
}

void LocalPortReducedSchurSolver::attachGlobalCoarse(
    const GlobalInterfaceCoarseModel&, bool)
{
    removed("global interface coarse attachment");
}

double LocalPortReducedSchurSolver::relativeProjectionError(
    const std::vector<double>&) const
{
    removed("port projection diagnostic");
}

int LocalPortReducedSchurSolver::coarseDimension() const
{
    return 0;
}

std::size_t LocalPortReducedSchurSolver::factorMemoryBytes() const
{
    return 0;
}

std::vector<PortPatch> buildOptimalPortPatches(
    const Mesh&, const ddm_schur::InterfacePartition&, int)
{
    removed("optimal-port patch construction");
}

std::vector<PortTopologyAudit> auditOptimalPortTopology(
    const Mesh&, const CaseConfig&,
    const ddm_schur::InterfacePartition&, const std::vector<int>&,
    int, int, const std::string&, int)
{
    removed("optimal-port topology audit");
}

ReducedDynamicSchurOperator::ReducedDynamicSchurOperator(
    const local::Model& model, bool)
    : model_(model)
{
    removed("experimental reduced Schur operator");
}

GeneralizedTransferSourceBlocks buildGeneralizedTransferSourceBlocks(
    const local::Model&, const std::vector<double>&, int,
    const std::vector<double>&, const std::vector<double>&, int)
{
    removed("operator-derived particular trace sources");
}

// Complete the private implementation types needed by unique_ptr destruction.
struct PatchTransferOperator::PatchAlgebra {};
struct PatchTransferOperator::SparseTargetFactor {};

PatchTransferOperator::PatchTransferOperator(
    const ReducedDynamicSchurOperator& schur,
    PortPatch patch,
    const PatchTransferOptions& options)
    : schur_(schur), patch_(std::move(patch)), options_(options)
{
    removed("local transfer operator");
}

PatchTransferOperator::~PatchTransferOperator() = default;

void PatchTransferOperator::solveTargetResponse(
    const std::vector<double>&, std::vector<double>&)
{
    removed("local transfer target solve");
}

GlobalInterfaceCoarseModel buildGlobalInterfaceCoarsePrototype(
    const LocalPortModel&, LocalPortReducedSchurSolver&,
    const ReducedDynamicSchurOperator&, const std::vector<PortPatch>&,
    const std::vector<double>&, const GeneralizedTransferSourceBlocks&,
    const GlobalInterfaceCoarseOptions&)
{
    removed("global interface spectral coarse prototype");
}

void writeGlobalInterfaceCoarseDiagnostics(
    const GlobalInterfaceCoarseModel&, const std::string&,
    const std::filesystem::path&)
{
    removed("global interface coarse diagnostics");
}

std::size_t globalCoarseCurrentWorkingSetBytes()
{
    return 0;
}

GlobalRandomizedSchurResult buildGlobalRandomizedSchurPortSpace(
    const LocalPortModel&, LocalPortReducedSchurSolver&,
    const ReducedDynamicSchurOperator&, const std::vector<double>&,
    const GlobalRandomizedSchurOptions&)
{
    removed("global randomized Schur space");
}

void writeGlobalRandomizedSchurDiagnostics(
    const GlobalRandomizedSchurResult&, const std::string&,
    const std::filesystem::path&)
{
    removed("global randomized Schur diagnostics");
}

LocalPortModel buildSteklovSchurPortModel(
    const Mesh&, const ddm_schur::InterfacePartition&,
    const local::Model&, const std::vector<double>&,
    const std::vector<double>&, const SteklovPortOptions&)
{
    removed("local spectral Steklov port basis");
}

OptimalPortBuildResult buildOptimalTransferPortModel(
    const Mesh&, const ddm_schur::InterfacePartition&,
    const local::Model&, const std::vector<double>&,
    const std::vector<double>&, const std::vector<double>&, int,
    const std::vector<double>&, const std::vector<double>&, int,
    const OptimalTransferPortOptions&,
    const ReducedDynamicSchurOperator*)
{
    removed("optimal transfer port basis");
}

ResidualKrylovBuildResult buildResidualKrylovPortModel(
    const Mesh&, const ddm_schur::InterfacePartition&,
    const local::Model&, const std::vector<double>&,
    const std::vector<double>&, const std::vector<double>&, int,
    const std::vector<double>&, const std::vector<double>&, int,
    const ResidualKrylovPortOptions&,
    const ReducedDynamicSchurOperator*, const LocalPortModel*)
{
    removed("residual Krylov port enrichment");
}

ProjectionDiagnosisSummary runFullInterfaceProjectionDiagnosis(
    const Mesh&, const CaseConfig&, const ddm_schur::InterfacePartition&,
    const ThermalDescriptorSystem&, const local::Model&,
    const LocalPortModel&, const std::vector<double>&,
    const std::vector<double>&, const std::vector<double>&,
    const ProjectionDiagnosisOptions&)
{
    removed("port projection diagnosis");
}

namespace port {

RandomizedTransferBuildResult buildRandomizedTransferPortModel(
    const Mesh&, const CaseConfig&, const ddm_schur::InterfacePartition&,
    const local::Model&, const std::vector<double>&,
    const std::vector<double>&, const std::vector<double>&, int,
    const std::vector<double>&, const std::vector<double>&, int,
    const RandomizedTransferPortOptions&,
    const ReducedDynamicSchurOperator*)
{
    removed("randomized transfer port basis");
}

WeightedPortSubspaceComparison compareWeightedPortSubspaces(
    const LocalPortBasis&, const LocalPortBasis&, const local::Model&,
    const std::vector<double>&, const std::vector<double>&,
    const std::string&)
{
    removed("weighted port subspace comparison");
}

} // namespace port
} // namespace mor::transient
