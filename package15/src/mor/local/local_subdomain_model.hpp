#pragma once

#include "mor/types.hpp"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace mor::local {

struct Options {
    bool generate = false;
    std::filesystem::path loadPath;
    std::filesystem::path savePath;
    std::string method = "pod";
    std::string mode = "pure";
    std::string interfaceMode = "full";
    int rank = 0;
    std::vector<int> rankPerSubdomain;
    std::filesystem::path rankFile;
    double energyTolerance = 1.0e-10;
    double singularValueTolerance = 1.0e-12;
    int trainingCases = 12;
    int validationCases = 8;
    int testCases = 8;
    std::uint64_t seed = 20260722ULL;
    bool compareFom = true;
    bool reuseIdenticalSubdomains = false;
    int matrixFreeInterfaceThreshold = 20000;
};

struct SubdomainModel {
    int subdomain = -1;
    int interiorDofs = 0;
    int localInterfaceDofs = 0;
    int snapshots = 0;
    int numericalRank = 0;
    int rank = 0;
    std::uint64_t meshFingerprint = 0;
    std::uint64_t materialFingerprint = 0;
    std::uint64_t boundaryInterfaceFingerprint = 0;
    int templateId = -1;
    bool templateReused = false;
    double templateConsistencyDifference = 0.0;
    std::vector<int> interiorGlobalDofs;
    std::vector<int> interfaceGlobalDofs;
    std::vector<int> interfaceIndices;
    std::vector<double> referenceInterior;
    // Column-major interiorDofs x rank.
    std::vector<double> basis;
    // Optional column-major Petrov test basis W.  An empty vector denotes
    // the Galerkin path W=V.  The trial basis is always `basis`.
    std::vector<double> testBasis;
    std::vector<double> singularValues;
    std::vector<double> retainedEnergy;
    double orthogonalityError = 0.0;
    double snapshotExtractionSeconds = 0.0;
    double podSeconds = 0.0;
    double projectionSeconds = 0.0;
    double reducedInteriorSymmetryError = 0.0;
    double reducedInteriorMinimumEigenvalue = 0.0;
    double reducedInteriorMaximumEigenvalue = 0.0;
    double reducedInteriorConditionEstimate = 0.0;
    double couplingSymmetryError = 0.0;
    double localSchurSymmetryError = 0.0;

    // Dense row-major reduced blocks.
    std::vector<double> reducedInterior;          // rank x rank
    std::vector<double> reducedInteriorInterface; // rank x local interface DOFs
    std::vector<double> reducedInterfaceInterior; // local interface DOFs x rank
    std::vector<double> interiorReferenceImage;   // interior DOFs
    std::vector<double> interfaceReferenceImage;  // local interface DOFs
};

struct InterfaceEntry {
    int row = 0;
    int column = 0;
    double value = 0.0;
};

struct Model {
    int formatVersion = 6;
    int globalDofs = 0;
    int interfaceDofs = 0;
    int totalLocalRank = 0;
    // CN-consistent Petrov-LSPG has distinct interior test and trial spaces,
    // hence its augmented interface system is generally nonsymmetric.
    bool usesPetrovTestSpace = false;
    Fingerprints fingerprints;
    std::vector<int> interfaceGlobalDofs;
    // Sparse full-order A_GammaGamma entries. The exact reduced Schur update
    // remains in factored low-rank local form and is assembled without drops.
    std::vector<InterfaceEntry> interfaceEntries;
    std::vector<SubdomainModel> subdomains;
    double snapshotSeconds = 0.0;
    double basisSeconds = 0.0;
    double projectionSeconds = 0.0;
    double schurConstructionSeconds = 0.0;
    double offlineSeconds = 0.0;
    std::size_t modelBytes = 0;
};

struct ErrorMetrics {
    double globalRelativeResidual = 0.0;
    double interfaceRelativeResidual = 0.0;
    double relativeL2 = 0.0;
    double maximumAbsolute = 0.0;
    double maximumTemperatureError = 0.0;
};

struct OnlineTiming {
    double localReducedAssemblySeconds = 0.0;
    double interfaceSolveSeconds = 0.0;
    double proxySolveSeconds = 0.0;
    double coarseSolveSeconds = 0.0;
    double portForwardSolveSeconds = 0.0;
    double portCoreSolveSeconds = 0.0;
    double portBackSubstitutionSeconds = 0.0;
    double interfaceOperatorSeconds = 0.0;
    double interfacePreconditionerSeconds = 0.0;
    double interfaceOrthogonalizationSeconds = 0.0;
    double interfaceVectorUpdateSeconds = 0.0;
    double interfacePredictorSeconds = 0.0;
    double localRecoverySeconds = 0.0;
    double fullFieldReconstructionSeconds = 0.0;
    double correctionSeconds = 0.0;
    double totalSeconds = 0.0;
    int interfaceIterations = 0;
    int interfaceMatvecs = 0;
    double interfaceInitialRelativeResidual = 0.0;
    double interfaceRelativeResidual = 0.0;
    double interfacePredictorInitialRelativeResidual = 0.0;
    double interfacePredictorRelativeResidual = 0.0;
    bool interfacePredictorApplied = false;
    bool interfacePredictorAccepted = false;
    std::string interfaceKrylovActual = "not_run";
    bool interfaceKrylovFallback = false;
    std::string interfaceKrylovFallbackReason;
};

struct SolveResult {
    std::vector<double> interfaceTemperature;
    std::vector<double> temperature;
    // Interior coordinates from the accepted local recovery.  The transient
    // workflow uses these to form the next reduced RHS without projecting a
    // reconstructed global RHS back onto every local basis.
    std::vector<std::vector<double>> localReducedCoordinates;
    OnlineTiming timing;
    ErrorMetrics error;
    int zeroGuessIterations = -1;
    int correctedIterations = -1;
    double correctedTrueResidual = 0.0;
    std::string status = "not_run";
};

struct WorkflowResult {
    std::vector<double> nominalTemperature;
    int subdomains = 0;
    int interfaceDofs = 0;
    int totalLocalRank = 0;
    int correctionIterations = 0;
    int coarseDimension = 0;
    double setupSeconds = 0.0;
    double onlineSeconds = 0.0;
    double coarseSolveSeconds = 0.0;
    double globalResidual = 0.0;
    double interfaceResidual = 0.0;
    std::string status = "not_run";
};

std::size_t estimateModelBytes(const Model& model);

} // namespace mor::local
