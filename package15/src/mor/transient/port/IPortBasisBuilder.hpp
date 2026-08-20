#pragma once

#include "../local_port_reduced_schur.hpp"
#include "../transfer_operator.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

struct Mesh;
struct CaseConfig;
namespace ddm_schur { struct InterfacePartition; }

namespace mor::transient::port {

struct PortBasisBuildContext {
    const Mesh& mesh;
    const CaseConfig& physics;
    const ddm_schur::InterfacePartition& partition;
    const local::Model& dynamicModel;
    const std::vector<double>& traceMassDiagonal;
    const std::vector<double>& penaltyMassDiagonal;
    const GeneralizedTransferSourceBlocks& sourceBlocks;
    const ReducedDynamicSchurOperator& schur;
};

struct PortBasisResult {
    std::string methodName;
    int physicalInterfaceId = -1;
    int targetDofs = 0;
    int sourceDofs = 0;
    int mandatoryRank = 0;
    int portRank = 0;
    int basisDimension = 0;
    double basisBuildTime = 0.0;
    std::size_t memoryPeak = 0;
    int solverCalls = 0;
    int eigenIterations = 0;
    double residual = 0.0;
    double orthogonalityError = 0.0;
    int serializationVersion = 0;
    bool snapshotUsed = false;
    bool fomUsedForBasis = false;

    int requestedRank = 0;
    int oversampling = 0;
    int probeColumns = 0;
    int acceptedRank = 0;
    int powerIterations = 0;
    std::uint64_t seed = 0;
    int applyCount = 0;
    int transposeApplyCount = 0;
    int targetSolveCount = 0;
    int targetSolvePhase33Calls = 0;
    double basisErrorIndicator = 0.0;
    double weightedAdjointError = 0.0;
    double probeGenerationSeconds = 0.0;
    double transferApplySeconds = 0.0;
    double transposeApplySeconds = 0.0;
    double qrSeconds = 0.0;
    double serializationSeconds = 0.0;
    std::size_t probeMatrixBytes = 0;
    std::size_t sampleMatrixBytes = 0;
    std::size_t qrWorkspaceBytes = 0;
    std::size_t finalBasisBytes = 0;
    std::string sourceMode;
    std::string status = "not_run";
    PatchInnerSolverStatistics innerSolver;
    LocalPortBasis basis;
};

class IPortBasisBuilder {
public:
    virtual ~IPortBasisBuilder() = default;
    virtual std::string methodName() const = 0;
    virtual PortBasisResult build(
        const PortBasisBuildContext& context,
        const PortPatch& patch) = 0;
};

} // namespace mor::transient::port
