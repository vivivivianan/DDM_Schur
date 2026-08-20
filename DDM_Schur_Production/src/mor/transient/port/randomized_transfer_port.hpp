#pragma once

// Compatibility declarations for the removed randomized transfer-port method.

#include "IPortBasisBuilder.hpp"

#include <memory>
#include <string>
#include <vector>

namespace mor::transient::port {

struct RandomizedTransferPortOptions {
    int requestedRank = 8;
    int oversampling = 5;
    int powerIterations = 1;
    std::uint64_t seed = 12345;
    std::string innerProduct = "penalty-weighted-mass";
    std::string sourceMode = "trace-only";
    int oversamplingLayers = 0;
    double relativeDeflationTolerance = 1.0e-12;
    bool fluxAware = false;
    std::string fluxType = "both";
    std::vector<int> selectedInterfaceIds;
    PatchTransferOptions innerSolver;
};

struct RandomizedTransferBuildResult {
    LocalPortModel model;
    std::vector<PortBasisResult> interfaces;
    double patchSetupSeconds = 0.0;
    double sourceSetupSeconds = 0.0;
    double probeGenerationSeconds = 0.0;
    double transferApplySeconds = 0.0;
    double transposeApplySeconds = 0.0;
    double qrSeconds = 0.0;
    double serializationSeconds = 0.0;
    double totalSeconds = 0.0;
};

class RandomizedTransferBasisBuilder final : public IPortBasisBuilder {
public:
    explicit RandomizedTransferBasisBuilder(
        RandomizedTransferPortOptions options);

    std::string methodName() const override;
    PortBasisResult build(
        const PortBasisBuildContext& context,
        const PortPatch& patch) override;

private:
    RandomizedTransferPortOptions options_;
};

RandomizedTransferBuildResult buildRandomizedTransferPortModel(
    const Mesh& mesh,
    const CaseConfig& physics,
    const ddm_schur::InterfacePartition& partition,
    const local::Model& dynamicModel,
    const std::vector<double>& traceMassDiagonal,
    const std::vector<double>& penaltyMassDiagonal,
    const std::vector<double>& input,
    int sourceChannels,
    const std::vector<double>& boundaryLoad,
    const std::vector<double>& condensedHistory,
    int historyChannels,
    const RandomizedTransferPortOptions& options,
    const ReducedDynamicSchurOperator* sharedSchur = nullptr);

struct WeightedPortSubspaceComparison {
    double maximumPrincipalAngleRadians = 0.0;
    double projectorDifference = 0.0;
};

WeightedPortSubspaceComparison compareWeightedPortSubspaces(
    const LocalPortBasis& candidate,
    const LocalPortBasis& reference,
    const local::Model& dynamicModel,
    const std::vector<double>& traceMassDiagonal,
    const std::vector<double>& penaltyMassDiagonal,
    const std::string& innerProduct);

} // namespace mor::transient::port
