#pragma once

#include "mor/local/local_reduced_schur.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

struct Mesh;
struct CaseConfig;
namespace ddm_schur { struct InterfacePartition; }

namespace mor::transient {

struct PortPatch {
    int interfaceId = -1;
    int leftSubdomain = -1;
    int rightSubdomain = -1;
    std::vector<int> patchSubdomains;
    std::vector<int> target;
    std::vector<int> source;
    std::uint64_t targetFingerprint = 0;
    std::uint64_t sourceFingerprint = 0;
};

std::vector<PortPatch> buildOptimalPortPatches(
    const Mesh& mesh,
    const ddm_schur::InterfacePartition& partition,
    int oversamplingLayers);

struct PortTopologyAudit {
    int interfaceId = -1;
    int leftSubdomain = -1;
    int rightSubdomain = -1;
    int targetDofs = 0;
    int sourceDofs = 0;
    bool sourceEmpty = false;
    int mandatoryModeCount = 0;
    int heatSourceChannelCount = 0;
    int externalBoundaryChannelCount = 0;
    std::size_t eigensolverWorkspaceBytes = 0;
    std::size_t innerSolverWorkspaceBytes = 0;
    std::size_t transferWorkspaceBytes = 0;
    std::size_t basisStorageEstimateBytes = 0;
    std::size_t estimatedWorkspaceBytes = 0;
};

std::vector<PortTopologyAudit> auditOptimalPortTopology(
    const Mesh& mesh,
    const CaseConfig& physics,
    const ddm_schur::InterfacePartition& partition,
    const std::vector<int>& sourceSubdomains,
    int oversamplingLayers,
    int requestedRank,
    const std::string& requestedInnerSolver,
    int directRowLimit = 2048);

class ReducedDynamicSchurOperator {
public:
    explicit ReducedDynamicSchurOperator(
        const local::Model& model,
        // This affects only the iterative target-solver Jacobi
        // preconditioner.  apply/applyTranspose remain the exact same Schur
        // actions.  Large pilots may use the assembled-interface diagonal
        // after the exact diagonal setup has been shown to be intractable.
        bool exactJacobiDiagonal = true);

    int size() const { return model_.interfaceDofs; }
    const std::vector<double>& diagonal() const { return diagonal_; }
    const local::Model& model() const { return model_; }

    void apply(const std::vector<double>& input,
               std::vector<double>& output) const;
    void applyTranspose(const std::vector<double>& input,
                        std::vector<double>& output) const;
    void solveReducedInterior(std::size_t subdomainSlot,
                              std::vector<double>& rightHandSide) const;
    void solveReducedInteriorMultiple(
        std::size_t subdomainSlot,
        std::vector<double>& rightHandSides,
        int rightHandSideCount) const;

private:
    const local::Model& model_;
    std::vector<local::DenseSymmetricFactor> factors_;
    std::vector<int> rowPtr_;
    std::vector<int> columns_;
    std::vector<double> values_;
    std::vector<double> diagonal_;
};

struct PatchTransferOptions {
    std::string innerSolver = "auto";
    int directRowLimit = 2048;
    double relativeTolerance = 1.0e-10;
    int maximumIterations = 1000;
    int refinementMaximumIterations = 3;
    double refinementTolerance = 1.0e-10;
};

struct PatchInnerSolverStatistics {
    std::string requestedSolver;
    std::string actualSolver;
    std::string solver;
    std::string status = "not_run";
    std::string fallbackReason;
    int setupApplies = 0;
    int solveCalls = 0;
    int solveRightHandSides = 0;
    int totalIterations = 0;
    int maximumIterations = 0;
    double setupSeconds = 0.0;
    double totalSolveSeconds = 0.0;
    double maximumSolveSeconds = 0.0;
    double finalRelativeResidual = 0.0;
    double maximumRelativeResidual = 0.0;
    double relativeAsymmetry = 0.0;
    double diagonalShift = 0.0;
    std::size_t factorBytes = 0;
    int aTtDimension = 0;
    std::size_t aTtNonzeros = 0;
    double aTtAssemblySeconds = 0.0;
    double aTtFactorizationSeconds = 0.0;
    std::size_t aTtFactorBytes = 0;
    int reducedCorrectionRank = 0;
    double wSetupSeconds = 0.0;
    std::size_t wBytes = 0;
    int qDimension = 0;
    double qAssemblySeconds = 0.0;
    double qFactorizationSeconds = 0.0;
    std::size_t qBytes = 0;
    double transferApplySeconds = 0.0;
    double transposeApplySeconds = 0.0;
    double referenceTargetActionError = 0.0;
    double referenceCrossActionError = 0.0;
    double referenceCrossTransposeError = 0.0;
    double referenceActionCheckSeconds = 0.0;
    std::size_t peakIncrementalMemoryBytes = 0;
    double aTtSolveRelativeResidual = 0.0;
    double qSolvePreRefinementResidual = 0.0;
    double qSolveRelativeResidual = 0.0;
    int qRefinementIterations = 0;
    double woodburyPreRefinementResidual = 0.0;
    double woodburyPostRefinementResidual = 0.0;
    int refinementIterations = 0;
    double refinementResidual0 = 0.0;
    double refinementResidual1 = 0.0;
    double refinementResidual2 = 0.0;
    double refinementResidual3 = 0.0;
    double refinementCorrectionRelativeNorm = 0.0;
    double refinementReductionFactor = 1.0;
    bool refinementConverged = true;
    int refinementTriggeredSolveCalls = 0;
    double refinementSeconds = 0.0;
    int pardisoInternalRefinementSteps = 0;
    double qMinimumAbsoluteFactorDiagonal = 0.0;
    double qMaximumAbsoluteFactorDiagonal = 0.0;
    double qFactorDiagonalRatio = 0.0;
    double woodburyCancellationFactor = 0.0;
};

struct GeneralizedTransferSourceBlocks {
    // Column-major interfaceDofs x channels.  Every column is an already
    // condensed Schur-interface right-hand side.  The trace block remains
    // matrix-free and is represented by PortPatch::source.
    int interfaceDofs = 0;
    int inputChannels = 0;
    int boundaryChannels = 0;
    int historyChannels = 0;
    std::vector<double> input;
    std::vector<double> boundary;
    std::vector<double> history;
    std::uint64_t inputFingerprint = 0;
    std::uint64_t boundaryFingerprint = 0;
    std::uint64_t historyFingerprint = 0;
};

GeneralizedTransferSourceBlocks buildGeneralizedTransferSourceBlocks(
    const local::Model& model,
    const std::vector<double>& input,
    int sourceChannels,
    const std::vector<double>& boundaryLoad,
    const std::vector<double>& condensedHistory,
    int historyChannels);

class PatchTransferOperator {
public:
    PatchTransferOperator(const ReducedDynamicSchurOperator& schur,
                          PortPatch patch,
                          const PatchTransferOptions& options);
    ~PatchTransferOperator();

    const PortPatch& patch() const { return patch_; }
    int targetRows() const { return static_cast<int>(patch_.target.size()); }
    int sourceRows() const { return static_cast<int>(patch_.source.size()); }

    void apply(const std::vector<double>& source,
               std::vector<double>& target);
    void applyTranspose(const std::vector<double>& target,
                        std::vector<double>& source);
    void applyMultiple(const std::vector<double>& sources,
                       int rightHandSideCount,
                       std::vector<double>& targets);
    void applyTransposeMultiple(const std::vector<double>& targets,
                                int rightHandSideCount,
                                std::vector<double>& sources);
    void solveTargetResponse(const std::vector<double>& rightHandSide,
                             std::vector<double>& target);
    // Column-major targetRows x rightHandSideCount block solve.  The
    // Woodbury path uses one PARDISO phase-33 multi-RHS call per block.
    void solveTargetResponses(
        const std::vector<double>& rightHandSides,
        int rightHandSideCount,
        std::vector<double>& targets);
    void solveTargetResponsesTranspose(
        const std::vector<double>& rightHandSides,
        int rightHandSideCount,
        std::vector<double>& targets);
    void solveTargetResponseTranspose(
        const std::vector<double>& rightHandSide,
        std::vector<double>& target);
    // Exact matrix-free target Schur actions used by residual-driven port
    // enrichment.  These do not assemble or alter S_tt.
    void applyTargetAction(const std::vector<double>& input,
                           std::vector<double>& output) const;
    void applyTargetActionTranspose(const std::vector<double>& input,
                                    std::vector<double>& output) const;
    // Exact reduced-block S_to actions without a target solve.  These keep
    // the trace source matrix-free while avoiding a global-interface vector.
    void formTraceRightHandSide(const std::vector<double>& source,
                                std::vector<double>& target) const;
    void applyTraceRightHandSideTranspose(
        const std::vector<double>& target,
        std::vector<double>& source) const;

    std::vector<double> explicitMatrix();
    const PatchInnerSolverStatistics& statistics() const { return statistics_; }

private:
    struct PatchAlgebra;
    struct SparseTargetFactor;

    void applyTargetBlock(const std::vector<double>& input,
                          std::vector<double>& output,
                          bool transpose) const;
    void applyCrossBlock(const std::vector<double>& input,
                         std::vector<double>& output,
                         bool transpose) const;
    void initializePatchAlgebra();
    void initializeAssembledDense();
    void initializeWoodbury();
    void initializeIterativePreconditioner();
    bool solveFgmres(const std::vector<double>& rightHandSide,
                     bool transpose,
                     std::vector<double>& solution,
                     int& iterations,
                     double& relativeResidual);
    void solveTarget(std::vector<double>& rightHandSide, bool transpose);
    void solveWoodburyOnce(const std::vector<double>& rightHandSide,
                           bool transpose,
                           std::vector<double>& solution);
    void solveWoodburyMultipleOnce(
        const std::vector<double>& rightHandSides,
        int rightHandSideCount,
        bool transpose,
        std::vector<double>& solutions);
    void solveTargetMultiple(
        const std::vector<double>& rightHandSides,
        int rightHandSideCount,
        bool transpose,
        std::vector<double>& solutions);
    double exactTargetResidual(
        const std::vector<double>& rightHandSide,
        const std::vector<double>& solution,
        bool transpose,
        std::vector<double>* residual = nullptr) const;
    double sparseAResidual(
        const std::vector<double>& rightHandSide,
        const std::vector<double>& solution,
        bool transpose) const;
    double qResidual(
        const std::vector<double>& rightHandSide,
        const std::vector<double>& solution,
        bool transpose,
        std::vector<double>* residual = nullptr) const;
    void solveSparseA(const std::vector<double>& rightHandSide,
                      bool transpose,
                      std::vector<double>& solution) const;
    void solveQ(std::vector<double>& rightHandSide, bool transpose) const;

    const ReducedDynamicSchurOperator& schur_;
    PortPatch patch_;
    PatchTransferOptions options_;
    bool direct_ = false;
    bool assembledDense_ = false;
    bool woodbury_ = false;
    bool forceFgmres_ = false;
    local::DenseSymmetricFactor factor_;
    local::DenseSymmetricFactor qFactor_;
    std::vector<double> preconditionerDiagonal_;
    std::vector<double> woodburyU_;
    std::vector<double> woodburyV_;
    std::vector<double> woodburyW_;
    std::vector<double> woodburyWTranspose_;
    std::vector<double> woodburyQ_;
    std::unique_ptr<PatchAlgebra> algebra_;
    std::unique_ptr<SparseTargetFactor> sparseFactor_;
    std::unique_ptr<SparseTargetFactor> qSparseFactor_;
    PatchInnerSolverStatistics statistics_;
};

// Matrix-free transfer from
//   xi = [outer trace; physical inputs; boundary loads; previous local state]
// to the target trace.  S_tt^{-1} is never formed; the existing target solver
// is reused by apply and transpose apply.
class GeneralizedPatchTransferOperator {
public:
    GeneralizedPatchTransferOperator(
        const ReducedDynamicSchurOperator& schur,
        PortPatch patch,
        GeneralizedTransferSourceBlocks blocks,
        bool includeTrace,
        bool includeInput,
        bool includeBoundary,
        bool includeHistory,
        const PatchTransferOptions& options);

    const PortPatch& patch() const { return targetSolver_.patch(); }
    int targetRows() const { return targetSolver_.targetRows(); }
    int sourceRows() const;
    int traceRows() const { return includeTrace_ ? targetSolver_.sourceRows() : 0; }
    int inputRows() const { return includeInput_ ? blocks_.inputChannels : 0; }
    int boundaryRows() const {
        return includeBoundary_ ? blocks_.boundaryChannels : 0;
    }
    int historyRows() const {
        return includeHistory_ ? blocks_.historyChannels : 0;
    }
    std::uint64_t traceFingerprint() const {
        return targetSolver_.patch().sourceFingerprint;
    }
    std::uint64_t inputFingerprint() const {
        return includeInput_ ? blocks_.inputFingerprint : 0;
    }
    std::uint64_t boundaryFingerprint() const {
        return includeBoundary_ ? blocks_.boundaryFingerprint : 0;
    }
    std::uint64_t historyFingerprint() const {
        return includeHistory_ ? blocks_.historyFingerprint : 0;
    }

    void apply(const std::vector<double>& source,
               std::vector<double>& target);
    void applyTranspose(const std::vector<double>& target,
                        std::vector<double>& source);
    void applyMultiple(const std::vector<double>& sources,
                       int rightHandSideCount,
                       std::vector<double>& targets);
    void applyTransposeMultiple(const std::vector<double>& targets,
                                int rightHandSideCount,
                                std::vector<double>& sources);
    void solveTargetResponse(const std::vector<double>& rightHandSide,
                             std::vector<double>& target);
    PatchInnerSolverStatistics statistics() const;

private:
    void addTraceRightHandSide(const double* source,
                               std::vector<double>& target) const;
    void addColumns(const std::vector<double>& columns,
                    int channels,
                    const double* coefficients,
                    std::vector<double>& target) const;
    void transposeColumns(const std::vector<double>& columns,
                          int channels,
                          const std::vector<double>& target,
                          double* coefficients) const;

    const ReducedDynamicSchurOperator& schur_;
    GeneralizedTransferSourceBlocks blocks_;
    bool includeTrace_ = false;
    bool includeInput_ = false;
    bool includeBoundary_ = false;
    bool includeHistory_ = false;
    double transferApplySeconds_ = 0.0;
    double transposeApplySeconds_ = 0.0;
    PatchTransferOperator targetSolver_;
};

} // namespace mor::transient
