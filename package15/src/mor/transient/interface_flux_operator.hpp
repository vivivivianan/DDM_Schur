#pragma once

#include <cstddef>
#include <string>
#include <vector>

struct CaseConfig;
struct Mesh;

namespace mor::local { struct Model; }

namespace mor::transient {

class ReducedDynamicSchurOperator;

struct InterfaceFluxResponse {
    int interfaceId = -1;
    int leftSubdomain = -1;
    int rightSubdomain = -1;
    std::vector<int> facePairIds;
    std::vector<int> triangleIds;
    std::vector<double> areas;
    std::vector<double> penalties;
    std::vector<double> leftPhysicalOutward;
    std::vector<double> rightPhysicalOutward;
    std::vector<double> numerical;
};

struct InterfaceFluxOrientationAudit {
    int interfaceId = -1;
    std::string ownerSide = "temperature_trace_unique_owner";
    int masterSubdomain = -1;
    int slaveSubdomain = -1;
    double ownerNormalX = 0.0;
    double ownerNormalY = 0.0;
    double ownerNormalZ = 0.0;
    double masterNormalX = 0.0;
    double masterNormalY = 0.0;
    double masterNormalZ = 0.0;
    double slaveNormalX = 0.0;
    double slaveNormalY = 0.0;
    double slaveNormalZ = 0.0;
    double normalDotProduct = 0.0;
    double minimumNormalDotProduct = 0.0;
    double maximumNormalDotProduct = 0.0;
    int targetDofCount = 0;
    int fluxTriangleCount = 0;
    int fluxQuadratureCount = 0;
    int junctionDuplicateCount = 0;
    bool deterministicOwnerOrdering = false;
    bool oppositeNormals = false;
    bool projectionUsesAbsoluteValue = false;
    bool dofQuadratureOneToOneAssumed = false;
    std::string physicalFluxSignConvention;
    std::string numericalFluxSignConvention;
    std::string mappingConvention;
    std::string status = "not_run";
};

class InterfaceFluxOperator {
public:
    InterfaceFluxOperator(
        const Mesh& mesh,
        const CaseConfig& physics,
        const local::Model& dynamicModel,
        const ReducedDynamicSchurOperator& schur,
        int interfaceId,
        int leftSubdomain,
        int rightSubdomain);

    InterfaceFluxResponse evaluateRecoveredField(
        const std::vector<double>& temperature) const;

    InterfaceFluxResponse applyGlobalTrace(
        const std::vector<double>& interfaceTemperature) const;

    InterfaceFluxResponse applyTargetTrace(
        const std::vector<int>& targetIndices,
        const std::vector<double>& targetTemperature) const;

    InterfaceFluxResponse evaluateAffine(
        const std::vector<double>& globalRhs,
        const std::vector<double>& referenceTemperature) const;

    InterfaceFluxResponse evaluateTotal(
        const std::vector<double>& globalRhs,
        const std::vector<double>& interfaceTemperature,
        const std::vector<double>& referenceTemperature) const;

    InterfaceFluxOrientationAudit auditOrientation(
        int targetDofCount,
        int junctionDuplicateCount) const;

    std::vector<double> recover(
        const std::vector<double>& globalRhs,
        const std::vector<double>& interfaceTemperature) const;

private:
    const Mesh& mesh_;
    const CaseConfig& physics_;
    const local::Model& model_;
    const ReducedDynamicSchurOperator& schur_;
    int interfaceId_ = -1;
    int leftSubdomain_ = -1;
    int rightSubdomain_ = -1;
};

double interfaceFluxRelativeDifference(
    const InterfaceFluxResponse& candidate,
    const InterfaceFluxResponse& reference,
    const std::string& fluxType);

double interfaceFluxRelativeLinearityError(
    const InterfaceFluxResponse& combined,
    const InterfaceFluxResponse& first,
    double firstScale,
    const InterfaceFluxResponse& second,
    double secondScale,
    const std::string& fluxType);

double interfaceFluxNorm(
    const InterfaceFluxResponse& response,
    const std::string& fluxType);

InterfaceFluxResponse combineInterfaceFlux(
    const InterfaceFluxResponse& first,
    double firstScale,
    const InterfaceFluxResponse& second,
    double secondScale);

} // namespace mor::transient
