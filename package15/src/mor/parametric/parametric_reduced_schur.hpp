#pragma once

#include "affine_fem_components.hpp"
#include "ddm_schur/interface_operator.hpp"
#include "mor/types.hpp"

#include <filesystem>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace mor::parametric {

struct ParametricLocalBlock {
    int domainId = -1;
    int rank = 0;
    std::uint64_t dofOrderingHash = 0;
    std::vector<int> globalDofs;
    std::vector<double> basis;
    std::vector<double> singularValues;
    std::vector<double> aiiConstant;
    std::vector<double> aiiLinear;
    std::vector<double> aiiHarmonic;
    std::vector<double> aiGammaConstant;
    std::vector<double> aiGammaLinear;
    std::vector<double> aiGammaHarmonic;
    std::vector<double> aGammaIConstant;
    std::vector<double> aGammaILinear;
    std::vector<double> aGammaIHarmonic;
    std::vector<double> rhsConstant;
    std::vector<double> rhsLinear;
    std::vector<double> rhsHarmonic;
    std::vector<double> sourceChannels;
};

struct ParametricReducedModel {
    int formatVersion = 3;
    int globalDofs = 0;
    int interfaceDofs = 0;
    int interfaceRank = 0;
    int sourceChannels = 0;
    AffineParameter parameter;
    Fingerprints fingerprints;
    std::uint64_t affineConstantHash = 0;
    std::uint64_t affineLinearHash = 0;
    std::vector<std::uint64_t> affineHarmonicHashes;
    std::uint64_t sourceDefinitionHash = 0;
    std::vector<DeploymentDof> dofs;
    std::vector<int> interfaceGlobalDofs;
    std::vector<double> referenceTemperature;
    std::vector<double> interfaceBasis;
    std::vector<double> interfaceSingularValues;
    std::vector<double> aGammaGammaConstant;
    std::vector<double> aGammaGammaLinear;
    std::vector<double> aGammaGammaHarmonic;
    std::vector<double> rhsGammaConstant;
    std::vector<double> rhsGammaLinear;
    std::vector<double> rhsGammaHarmonic;
    std::vector<double> sourceGamma;
    std::vector<double> nominalPowersW;
    std::vector<double> minimumPowersW;
    std::vector<double> maximumPowersW;
    std::vector<int> sourceSubdomains;
    std::vector<int> sourceDomainEntities;
    std::vector<ParametricLocalBlock> locals;
    double basisSeconds = 0.0;
    double projectionPreparationSeconds = 0.0;
    double interfaceProjectionSeconds = 0.0;
    double localProjectionSeconds = 0.0;
    double projectionSeconds = 0.0;
    double serializationSeconds = 0.0;
    std::uint64_t fileBytes = 0;
};

struct ParametricTrainingSnapshots {
    std::vector<double> referenceTemperature;
    SnapshotDatabase interfaceSnapshots;
    std::vector<SnapshotDatabase> interiors;
};

struct ParametricOnlineTiming {
    double coefficientSeconds = 0.0;
    double localAssemblySeconds = 0.0;
    double localFactorSeconds = 0.0;
    double reducedSchurAssemblySeconds = 0.0;
    double reducedSolveSeconds = 0.0;
    double reconstructionSeconds = 0.0;
    double outputSeconds = 0.0;
    double totalSeconds = 0.0;
};

struct ParametricOnlineResult {
    double parameterValue = 0.0;
    std::vector<double> powersW;
    std::vector<double> temperature;
    std::vector<double> interfaceTemperature;
    ParametricOnlineTiming timing;
    bool extrapolated = false;
    std::string status = "not_run";
};

// All mutable storage needed by one pure-online solve.  The workspace is sized
// once from a serialized model and then reused for every parameter/power case;
// the reconstructed fields always describe the most recently completed case.
struct ParametricDenseFactorWorkspace {
    std::vector<double> values;
    std::vector<int> pivots;
    bool cholesky = true;
};

struct ParametricLocalWorkspace {
    int rank = 0;
    ParametricDenseFactorWorkspace factor;
    std::vector<double> eliminatedCoupling;
    std::vector<double> rhsSolution;
    std::vector<double> aGammaI;
};

struct ParametricRomWorkspace {
    int globalDofs = 0;
    int interfaceDofs = 0;
    int interfaceRank = 0;
    int sourceChannels = 0;
    std::vector<double> parameterCoefficients;
    std::vector<double> powersW;
    std::vector<double> factorBackup;
    ParametricDenseFactorWorkspace schurFactor;
    std::vector<double> interfaceCoordinates;
    std::vector<ParametricLocalWorkspace> locals;
    std::vector<double> fullTemperature;
    std::vector<double> interfaceTemperature;
    std::vector<double> errorScratch;
    std::vector<char> outputScratch;

    void initialize(const ParametricReducedModel& model);
    void resetForNextCase();
    std::size_t workspaceBytes() const;
    std::size_t reconstructionBytes() const;
};

struct ParametricModelMemoryBreakdown {
    std::size_t interfaceBasisBytes = 0;
    std::size_t localBasisBytes = 0;
    std::size_t reducedBlockBytes = 0;
    std::size_t sourceBytes = 0;
    std::size_t metadataBytes = 0;

    std::size_t totalBytes() const;
};

ParametricReducedModel buildParametricReducedModel(
    const Mesh& mesh,
    const AffineFemComponents& components,
    const ddm_schur::InterfacePartition& partition,
    const ParametricTrainingSnapshots& snapshots,
    const mor::Options& options);

ParametricOnlineResult solveParametricRom(
    const ParametricReducedModel& model,
    double parameterValue,
    const std::vector<double>& powersW,
    bool allowExtrapolation);

ParametricOnlineResult solveParametricRomInWorkspace(
    const ParametricReducedModel& model,
    double parameterValue,
    const std::vector<double>& powersW,
    bool allowExtrapolation,
    ParametricRomWorkspace& workspace);

ParametricModelMemoryBreakdown parametricModelMemoryBreakdown(
    const ParametricReducedModel& model);

void saveParametricModel(const ParametricReducedModel& model,
                         const std::filesystem::path& directory);

ParametricReducedModel loadParametricModel(
    const std::filesystem::path& directory,
    const Fingerprints* expectedFingerprints = nullptr,
    const AffineFemComponents* expectedAffine = nullptr);

void truncateParametricModel(ParametricReducedModel& model,
                             int interfaceRank,
                             int localRank);

double modelOrthogonalityError(const std::vector<double>& basis,
                               int rows,
                               int columns);

} // namespace mor::parametric
