#include "interface_flux_operator.hpp"

#include "sipg_core.hpp"
#include "mor/local/local_reduced_schur.hpp"
#include "mor/mor_diagnostics.hpp"
#include "mor/transient/transfer_operator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace mor::transient {
namespace {

void requireSameLayout(
    const InterfaceFluxResponse& left,
    const InterfaceFluxResponse& right)
{
    if (left.interfaceId != right.interfaceId
        || left.facePairIds != right.facePairIds
        || left.triangleIds != right.triangleIds
        || left.areas.size() != right.areas.size()
        || left.leftPhysicalOutward.size() != right.areas.size()
        || left.rightPhysicalOutward.size() != right.areas.size()
        || left.numerical.size() != right.areas.size()) {
        throw std::runtime_error(
            "[Flux operator] Interface flux layouts do not match.");
    }
}

double squaredDifference(
    const InterfaceFluxResponse& candidate,
    const InterfaceFluxResponse& reference,
    const std::string& fluxType,
    bool referenceNorm)
{
    requireSameLayout(candidate, reference);
    long double value = 0.0L;
    for (std::size_t row = 0; row < reference.areas.size(); ++row) {
        const long double area = reference.areas[row];
        if (fluxType == "numerical") {
            const double delta = referenceNorm
                ? reference.numerical[row]
                : candidate.numerical[row]
                    - reference.numerical[row];
            value += area * delta * delta;
        } else if (fluxType == "physical") {
            const double left = referenceNorm
                ? reference.leftPhysicalOutward[row]
                : candidate.leftPhysicalOutward[row]
                    - reference.leftPhysicalOutward[row];
            const double right = referenceNorm
                ? reference.rightPhysicalOutward[row]
                : candidate.rightPhysicalOutward[row]
                    - reference.rightPhysicalOutward[row];
            value += area * (left * left + right * right);
        } else {
            throw std::runtime_error(
                "[Flux operator] Flux type must be physical or numerical.");
        }
    }
    return static_cast<double>(value);
}

} // namespace

InterfaceFluxOperator::InterfaceFluxOperator(
    const Mesh& mesh,
    const CaseConfig& physics,
    const local::Model& dynamicModel,
    const ReducedDynamicSchurOperator& schur,
    int interfaceId,
    int leftSubdomain,
    int rightSubdomain)
    : mesh_(mesh),
      physics_(physics),
      model_(dynamicModel),
      schur_(schur),
      interfaceId_(interfaceId),
      leftSubdomain_(leftSubdomain),
      rightSubdomain_(rightSubdomain)
{
    if (interfaceId < 0 || leftSubdomain < 0 || rightSubdomain < 0
        || leftSubdomain == rightSubdomain
        || schur.size() != dynamicModel.interfaceDofs) {
        throw std::runtime_error(
            "[Flux operator] Invalid interface operator construction.");
    }
}

std::vector<double> InterfaceFluxOperator::recover(
    const std::vector<double>& globalRhs,
    const std::vector<double>& interfaceTemperature) const
{
    if (globalRhs.size()
            != static_cast<std::size_t>(model_.globalDofs)
        || interfaceTemperature.size()
            != static_cast<std::size_t>(model_.interfaceDofs)) {
        throw std::runtime_error(
            "[Flux operator] Recovery dimensions are invalid.");
    }
    std::vector<double> recovered(
        static_cast<std::size_t>(model_.globalDofs), 0.0);
    for (int gamma = 0; gamma < model_.interfaceDofs; ++gamma) {
        recovered[static_cast<std::size_t>(
            model_.interfaceGlobalDofs[
                static_cast<std::size_t>(gamma)])]
            = interfaceTemperature[static_cast<std::size_t>(gamma)];
    }
    for (std::size_t slot = 0;
         slot < model_.subdomains.size(); ++slot) {
        const local::SubdomainModel& local =
            model_.subdomains[slot];
        std::vector<double> reduced(
            static_cast<std::size_t>(local.rank), 0.0);
        for (int mode = 0; mode < local.rank; ++mode) {
            const double* basis = local.basis.data()
                + static_cast<std::size_t>(
                    mode * local.interiorDofs);
            long double value = 0.0L;
            for (int row = 0; row < local.interiorDofs; ++row) {
                value += static_cast<long double>(basis[row])
                    * globalRhs[static_cast<std::size_t>(
                        local.interiorGlobalDofs[
                            static_cast<std::size_t>(row)])];
            }
            for (int gamma = 0;
                 gamma < local.localInterfaceDofs; ++gamma) {
                value -= static_cast<long double>(
                    local.reducedInteriorInterface[
                        static_cast<std::size_t>(
                            mode * local.localInterfaceDofs
                            + gamma)])
                    * interfaceTemperature[
                        static_cast<std::size_t>(
                            local.interfaceIndices[
                                static_cast<std::size_t>(gamma)])];
            }
            reduced[static_cast<std::size_t>(mode)] =
                static_cast<double>(value);
        }
        schur_.solveReducedInterior(slot, reduced);
        for (int row = 0; row < local.interiorDofs; ++row) {
            long double value = 0.0L;
            for (int mode = 0; mode < local.rank; ++mode) {
                value += static_cast<long double>(
                    local.basis[static_cast<std::size_t>(
                        mode * local.interiorDofs + row)])
                    * reduced[static_cast<std::size_t>(mode)];
            }
            recovered[static_cast<std::size_t>(
                local.interiorGlobalDofs[
                    static_cast<std::size_t>(row)])]
                = static_cast<double>(value);
        }
    }
    return recovered;
}

InterfaceFluxResponse InterfaceFluxOperator::evaluateRecoveredField(
    const std::vector<double>& temperature) const
{
    const DetailedInterfacePhysicsMetrics metrics =
        calculateDetailedInterfacePhysicsMetricsForSubdomains(
            mesh_, physics_, temperature,
            leftSubdomain_, rightSubdomain_);
    InterfaceFluxResponse response;
    response.interfaceId = interfaceId_;
    response.leftSubdomain = leftSubdomain_;
    response.rightSubdomain = rightSubdomain_;
    response.facePairIds.reserve(metrics.triangles.size());
    response.triangleIds.reserve(metrics.triangles.size());
    response.areas.reserve(metrics.triangles.size());
    response.penalties.reserve(metrics.triangles.size());
    response.leftPhysicalOutward.reserve(metrics.triangles.size());
    response.rightPhysicalOutward.reserve(metrics.triangles.size());
    response.numerical.reserve(metrics.triangles.size());
    for (const InterfaceTriangleFluxRecord& row : metrics.triangles) {
        response.facePairIds.push_back(row.facePairId);
        response.triangleIds.push_back(row.integrationTriangleId);
        response.areas.push_back(row.area);
        response.penalties.push_back(row.sipgPenalty);
        response.leftPhysicalOutward.push_back(
            row.leftPhysicalNormalFlux);
        response.rightPhysicalOutward.push_back(
            row.rightPhysicalNormalFlux);
        response.numerical.push_back(row.sipgNumericalFlux);
    }
    if (response.areas.empty()) {
        throw std::runtime_error(
            "[Flux operator] Selected physical interface has no "
            "quadrature triangles.");
    }
    return response;
}

InterfaceFluxResponse InterfaceFluxOperator::applyGlobalTrace(
    const std::vector<double>& interfaceTemperature) const
{
    const std::vector<double> zeroRhs(
        static_cast<std::size_t>(model_.globalDofs), 0.0);
    return evaluateRecoveredField(
        recover(zeroRhs, interfaceTemperature));
}

InterfaceFluxResponse InterfaceFluxOperator::applyTargetTrace(
    const std::vector<int>& targetIndices,
    const std::vector<double>& targetTemperature) const
{
    if (targetIndices.size() != targetTemperature.size()) {
        throw std::runtime_error(
            "[Flux operator] Target trace dimensions are invalid.");
    }
    std::vector<double> full(
        static_cast<std::size_t>(model_.interfaceDofs), 0.0);
    for (std::size_t row = 0; row < targetIndices.size(); ++row) {
        const int gamma = targetIndices[row];
        if (gamma < 0 || gamma >= model_.interfaceDofs) {
            throw std::runtime_error(
                "[Flux operator] Target trace index is out of range.");
        }
        full[static_cast<std::size_t>(gamma)] =
            targetTemperature[row];
    }
    return applyGlobalTrace(full);
}

InterfaceFluxResponse InterfaceFluxOperator::evaluateAffine(
    const std::vector<double>& globalRhs,
    const std::vector<double>& referenceTemperature) const
{
    std::vector<double> zeroGamma(
        static_cast<std::size_t>(model_.interfaceDofs), 0.0);
    std::vector<double> field = recover(globalRhs, zeroGamma);
    if (field.size() != referenceTemperature.size()) {
        throw std::runtime_error(
            "[Flux operator] Reference temperature size mismatch.");
    }
    for (std::size_t row = 0; row < field.size(); ++row) {
        field[row] += referenceTemperature[row];
    }
    return evaluateRecoveredField(field);
}

InterfaceFluxResponse InterfaceFluxOperator::evaluateTotal(
    const std::vector<double>& globalRhs,
    const std::vector<double>& interfaceTemperature,
    const std::vector<double>& referenceTemperature) const
{
    std::vector<double> field =
        recover(globalRhs, interfaceTemperature);
    if (field.size() != referenceTemperature.size()) {
        throw std::runtime_error(
            "[Flux operator] Reference temperature size mismatch.");
    }
    for (std::size_t row = 0; row < field.size(); ++row) {
        field[row] += referenceTemperature[row];
    }
    return evaluateRecoveredField(field);
}

InterfaceFluxOrientationAudit
InterfaceFluxOperator::auditOrientation(
    int targetDofCount,
    int junctionDuplicateCount) const
{
    InterfaceFluxOrientationAudit audit;
    audit.interfaceId = interfaceId_;
    audit.targetDofCount = targetDofCount;
    audit.junctionDuplicateCount = junctionDuplicateCount;
    audit.deterministicOwnerOrdering =
        junctionDuplicateCount == 0;
    audit.physicalFluxSignConvention =
        "left/right outward; conservation is q_left+q_right=0";
    audit.numericalFluxSignConvention =
        "owner-left normal; -average(k grad T dot n)+tau jump";
    audit.mappingConvention =
        "temperature DOFs map through recovered FE field; flux uses "
        "face quadrature, never node-to-quadrature indexing";
    audit.projectionUsesAbsoluteValue = false;
    audit.dofQuadratureOneToOneAssumed = false;
    bool initialized = false;
    audit.minimumNormalDotProduct =
        std::numeric_limits<double>::max();
    audit.maximumNormalDotProduct =
        -std::numeric_limits<double>::max();
    // The production SIPG face evaluator uses the fixed 4x4 Duffy-rule
    // quadrature from makeTriangleQuadrature().
    constexpr int quadraturePerTriangle = 16;
    for (const InterfaceFace& face : mesh_.interfaceFaces) {
        const Tet& left =
            mesh_.tets[static_cast<std::size_t>(face.leftTet)];
        const Tet& right =
            mesh_.tets[static_cast<std::size_t>(face.rightTet)];
        const std::pair<int, int> adjacent{
            std::min(left.subdomain, right.subdomain),
            std::max(left.subdomain, right.subdomain)};
        const std::pair<int, int> selected{
            std::min(leftSubdomain_, rightSubdomain_),
            std::max(leftSubdomain_, rightSubdomain_)};
        if (adjacent != selected) {
            continue;
        }
        const Vec3 master = normalized(face.leftNormal);
        const Vec3 slave = normalized(face.rightNormal);
        const double normalDot = dot(master, slave);
        if (!initialized) {
            audit.masterSubdomain = left.subdomain;
            audit.slaveSubdomain = right.subdomain;
            audit.ownerNormalX = master.x;
            audit.ownerNormalY = master.y;
            audit.ownerNormalZ = master.z;
            audit.masterNormalX = master.x;
            audit.masterNormalY = master.y;
            audit.masterNormalZ = master.z;
            audit.slaveNormalX = slave.x;
            audit.slaveNormalY = slave.y;
            audit.slaveNormalZ = slave.z;
            audit.normalDotProduct = normalDot;
            initialized = true;
        }
        audit.minimumNormalDotProduct =
            std::min(audit.minimumNormalDotProduct, normalDot);
        audit.maximumNormalDotProduct =
            std::max(audit.maximumNormalDotProduct, normalDot);
        audit.fluxTriangleCount +=
            static_cast<int>(face.integrationTriangles.size());
        audit.fluxQuadratureCount += quadraturePerTriangle
            * static_cast<int>(face.integrationTriangles.size());
    }
    audit.oppositeNormals = initialized
        && audit.maximumNormalDotProduct <= -1.0 + 1.0e-10;
    audit.status = initialized
        && audit.oppositeNormals
        && audit.targetDofCount > 0
        && audit.fluxQuadratureCount > 0
        && audit.junctionDuplicateCount == 0
        ? "passed" : "orientation_gate_failed";
    return audit;
}

double interfaceFluxNorm(
    const InterfaceFluxResponse& response,
    const std::string& fluxType)
{
    return std::sqrt(std::max(0.0,
        squaredDifference(
            response, response, fluxType, true)));
}

double interfaceFluxRelativeDifference(
    const InterfaceFluxResponse& candidate,
    const InterfaceFluxResponse& reference,
    const std::string& fluxType)
{
    const double error = std::sqrt(std::max(
        0.0, squaredDifference(
            candidate, reference, fluxType, false)));
    return error / std::max(
        1.0e-300, interfaceFluxNorm(reference, fluxType));
}

InterfaceFluxResponse combineInterfaceFlux(
    const InterfaceFluxResponse& first,
    double firstScale,
    const InterfaceFluxResponse& second,
    double secondScale)
{
    requireSameLayout(first, second);
    InterfaceFluxResponse result = first;
    for (std::size_t row = 0; row < result.areas.size(); ++row) {
        result.leftPhysicalOutward[row] =
            firstScale * first.leftPhysicalOutward[row]
            + secondScale * second.leftPhysicalOutward[row];
        result.rightPhysicalOutward[row] =
            firstScale * first.rightPhysicalOutward[row]
            + secondScale * second.rightPhysicalOutward[row];
        result.numerical[row] =
            firstScale * first.numerical[row]
            + secondScale * second.numerical[row];
    }
    return result;
}

double interfaceFluxRelativeLinearityError(
    const InterfaceFluxResponse& combined,
    const InterfaceFluxResponse& first,
    double firstScale,
    const InterfaceFluxResponse& second,
    double secondScale,
    const std::string& fluxType)
{
    return interfaceFluxRelativeDifference(
        combined,
        combineInterfaceFlux(
            first, firstScale, second, secondScale),
        fluxType);
}

} // namespace mor::transient
