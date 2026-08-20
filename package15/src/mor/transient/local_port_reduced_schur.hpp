#pragma once

#include "mor/local/local_reduced_schur.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

struct CaseConfig;
struct Mesh;
namespace ddm_schur { struct InterfacePartition; }

namespace mor::transient {

struct GlobalInterfaceCoarseModel;

struct LocalPortOptions {
    // M7 remains the default.  M8 builders are dispatched explicitly rather
    // than replacing the snapshot/POD path in place.
    std::string basisMethod = "port-pod";
    int requestedRank = 16;
    double discardedEnergyTolerance = 1.0e-8;
    double relativeTolerance = 1.0e-12;
    std::filesystem::path rankFile;
    double temperatureWeight = 1.0;
    double fluxWeight = 1.0;
    double residualWeight = 1.0;
};

struct LocalPortBasis {
    int interfaceId = -1;
    int leftSubdomain = -1;
    int rightSubdomain = -1;
    int rows = 0;
    int rank = 0;
    int candidateColumns = 0;
    int acceptedColumns = 0;
    int templateId = -1;
    bool templateReused = false;
    std::uint64_t fingerprint = 0;
    std::uint64_t targetFingerprint = 0;
    std::uint64_t sourceFingerprint = 0;
    std::uint64_t traceSourceFingerprint = 0;
    std::uint64_t inputSourceFingerprint = 0;
    std::uint64_t boundarySourceFingerprint = 0;
    std::uint64_t historySourceFingerprint = 0;
    double retainedEnergy = 0.0;
    double orthogonalityError = 0.0;
    double transferIndicator = 0.0;
    int mandatoryModes = 0;
    int spectralModes = 0;
    int requestedTransferRank = 0;
    int traceSourceRows = 0;
    int inputSourceRows = 0;
    int boundarySourceRows = 0;
    int historySourceRows = 0;
    std::vector<double> spectralValues;
    std::vector<double> spectralResiduals;
    std::vector<int> patchSubdomains;
    std::vector<int> sourceIndices;
    std::vector<int> interfaceIndices;
    // Column-major rows x rank.  A column contains both nonmatching sides of
    // one physical SIPG interface and therefore has one shared port
    // coordinate.
    std::vector<double> basis;
};

struct LocalPortModel {
    // Version 7 adds operator-side history-channel compression metadata.
    // The reader keeps versions 1-6, including the exact M7
    // LPORT001/version-1 path, available.
    int formatVersion = 7;
    std::string basisMethod = "port-pod";
    std::string ablationMode = "mandatory-transfer";
    std::string sourceMode = "trace-only";
    std::string methodDescription;
    std::string historyCompressionMethod = "none";
    int historyCompressionRank = 0;
    double historyCompressionTolerance = 0.0;
    std::string innerProduct = "euclidean";
    std::string rankMode = "fixed";
    double timeStep = 0.0;
    int oversamplingLayers = 0;
    std::uint64_t meshFingerprint = 0;
    std::uint64_t materialFingerprint = 0;
    std::uint64_t operatorFingerprint = 0;
    std::uint64_t inputFingerprint = 0;
    std::uint64_t penaltyFingerprint = 0;
    std::uint64_t boundaryFingerprint = 0;
    std::uint64_t sourceMetadataFingerprint = 0;
    std::uint64_t traceSourceFingerprint = 0;
    std::uint64_t generalizedInputFingerprint = 0;
    std::uint64_t generalizedBoundaryFingerprint = 0;
    std::uint64_t generalizedHistoryFingerprint = 0;
    std::uint64_t rankFileFingerprint = 0;
    int requestedRank = 0;
    int minimumRank = 0;
    int maximumRank = 0;
    double eigenvalueTolerance = 0.0;
    double eigensolverTolerance = 0.0;
    int eigensolverMaximumIterations = 0;
    double relativeDeflationTolerance = 0.0;
    std::string innerSolver;
    double innerSolverTolerance = 0.0;
    int innerSolverMaximumIterations = 0;
    std::string buildTimestamp;
    std::string sourceCommit;
    int fullInterfaceDofs = 0;
    int reducedInterfaceDofs = 0;
    std::vector<int> interfaceGlobalDofs;
    std::vector<LocalPortBasis> ports;
    double snapshotSeconds = 0.0;
    double basisSeconds = 0.0;
    std::size_t modelBytes = 0;
};

struct LocalPortSnapshotFamilies {
    // Target deployment trajectory modes are inserted exactly before POD so
    // highly conditioned Schur systems cannot amplify a tiny POD truncation.
    std::vector<std::vector<double>> mandatoryTemperature;
    std::vector<std::vector<double>> temperature;
    std::vector<std::vector<double>> flux;
    std::vector<std::vector<double>> residual;
};

struct LocalPortSolveResult {
    local::SolveResult solution;
    double portRelativeResidual = 0.0;
    double reducedRelativeResidual = 0.0;
    double condensedRhsSeconds = 0.0;
    double reducedSolveSeconds = 0.0;
};

LocalPortModel buildLocalPortModel(
    const Mesh& mesh,
    const CaseConfig& physics,
    const ddm_schur::InterfacePartition& partition,
    const LocalPortSnapshotFamilies& snapshots,
    const LocalPortOptions& options);

void saveLocalPortModel(const LocalPortModel& model,
                        const std::filesystem::path& path);
LocalPortModel loadLocalPortModel(const std::filesystem::path& path);

class LocalPortReducedSchurSolver {
public:
    LocalPortReducedSchurSolver(const local::Model& dynamicModel,
                                const LocalPortModel& portModel);

    LocalPortSolveResult solve(const std::vector<double>& globalRhs) const;
    void applyFullInterface(const std::vector<double>& input,
                            std::vector<double>& output) const;
    void attachGlobalCoarse(
        const GlobalInterfaceCoarseModel& coarseModel,
        bool includeLocalBasis = true);
    void localGalerkinInterfaceResponse(
        const std::vector<double>& interfaceRhs,
        std::vector<double>& response) const;
    void projectSchurEnergyComplement(
        const std::vector<double>& input,
        std::vector<double>& output) const;
    void restrictLocal(const std::vector<double>& full,
                       std::vector<double>& reduced) const;
    double relativeProjectionError(const std::vector<double>& full) const;
    int portDimension() const {
        return ((coarseModel_ && !includeLocalBasis_)
                    ? 0 : portModel_.reducedInterfaceDofs)
            + coarseDimension();
    }
    int localPortDimension() const {
        return portModel_.reducedInterfaceDofs;
    }
    int coarseDimension() const;
    double assemblySeconds() const { return assemblySeconds_; }
    double factorizationSeconds() const { return factorizationSeconds_; }
    double relativeAsymmetry() const { return relativeAsymmetry_; }
    double augmentedRelativeAsymmetry() const {
        return augmentedRelativeAsymmetry_;
    }
    std::size_t factorMemoryBytes() const;

private:
    void lift(const std::vector<double>& reduced,
              std::vector<double>& full) const;
    void restrict(const std::vector<double>& full,
                  std::vector<double>& reduced) const;
    void solveLocalCoordinates(std::vector<double>& values) const;
    void liftCoarse(const std::vector<double>& reduced,
                    std::vector<double>& full) const;
    void restrictCoarse(const std::vector<double>& full,
                        std::vector<double>& reduced) const;

    const local::Model& dynamicModel_;
    const LocalPortModel& portModel_;
    std::vector<local::DenseSymmetricFactor> localFactors_;
    std::vector<int> rowPtr_;
    std::vector<int> columnIndices_;
    std::vector<double> interfaceValues_;
    local::DenseSymmetricFactor portFactor_;
    std::vector<double> portLu_;
    std::vector<int> portPivots_;
    const GlobalInterfaceCoarseModel* coarseModel_ = nullptr;
    // Row-major local-rank by coarse-rank cross block.
    std::vector<double> localCoarse_;
    local::DenseSymmetricFactor coarseBorderFactor_;
    bool includeLocalBasis_ = true;
    double augmentedRelativeAsymmetry_ = 0.0;
    double assemblySeconds_ = 0.0;
    double factorizationSeconds_ = 0.0;
    double relativeAsymmetry_ = 0.0;
};

} // namespace mor::transient
