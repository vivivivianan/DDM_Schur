// Compact numerical diagnostics used during model construction: basis
// orthogonality, projected symmetry, eigenvalue bounds, and residual summaries.
// These routines never construct a basis from validation results.

#include "mor_diagnostics.hpp"

#include "sipg_core.hpp"
#include "config_io.hpp"
#include "mesh_loader.hpp"
#include "fem_assembly.hpp"
#include "linear_solvers.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace mor {

namespace {

SipgInterfaceScales sipgScales(const Mesh& mesh,
                               const CaseConfig& physics,
                               const InterfaceFace& face,
                               const Tet& left,
                               const Tet& right,
                               const Material& leftMaterial,
                               const Material& rightMaterial)
{
    const std::array<Vec3, 3> leftPhysical{{
        mesh.nodes[static_cast<std::size_t>(left.v[static_cast<std::size_t>(face.leftLocal[0])])].p,
        mesh.nodes[static_cast<std::size_t>(left.v[static_cast<std::size_t>(face.leftLocal[1])])].p,
        mesh.nodes[static_cast<std::size_t>(left.v[static_cast<std::size_t>(face.leftLocal[2])])].p}};
    const std::array<Vec3, 3> rightPhysical{{
        mesh.nodes[static_cast<std::size_t>(right.v[static_cast<std::size_t>(face.rightLocal[0])])].p,
        mesh.nodes[static_cast<std::size_t>(right.v[static_cast<std::size_t>(face.rightLocal[1])])].p,
        mesh.nodes[static_cast<std::size_t>(right.v[static_cast<std::size_t>(face.rightLocal[2])])].p}};
    const double leftArea = triangleArea(leftPhysical);
    const double rightArea = triangleArea(rightPhysical);
    const double leftVolume = std::abs(elementGeometry(mesh, left).detJ) / 6.0;
    const double rightVolume = std::abs(elementGeometry(mesh, right).detJ) / 6.0;
    const double hLeft = 3.0 * leftVolume / std::max(1.0e-30, leftArea);
    const double hRight = 3.0 * rightVolume / std::max(1.0e-30, rightArea);
    const Vec3 normal = norm(face.leftNormal) > 1.0e-30
        ? face.leftNormal : normalized(face.rightNormal * -1.0);
    const double kLeft = std::max(1.0e-30, normalConductivity(leftMaterial, normal));
    const double kRight = std::max(1.0e-30, normalConductivity(rightMaterial, normal));
    return sipgInterfaceScales(physics, kLeft, kRight, hLeft, hRight);
}

} // namespace

InterfacePhysicsMetrics calculateInterfacePhysicsMetrics(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::vector<double>& temperature)
{
    if (temperature.size() != mesh.nodes.size()) {
        throw std::runtime_error("[MOR] Interface diagnostic temperature has the wrong size.");
    }
    const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
    double area = 0.0;
    double jumpSquared = 0.0;
    double fluxMismatchSquared = 0.0;
    double fluxScaleSquared = 0.0;
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        const Tet& left = mesh.tets[static_cast<std::size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<std::size_t>(face.rightTet)];
        const ElementGeometry leftGeometry = elementGeometry(mesh, left);
        const ElementGeometry rightGeometry = elementGeometry(mesh, right);
        const Material& leftMaterial = materialForTet(physics, left);
        const Material& rightMaterial = materialForTet(physics, right);
        for (const auto& triangle : face.integrationTriangles) {
            const double jacobian = norm(cross(
                triangle[1] - triangle[0], triangle[2] - triangle[0]));
            for (const TriangleQuadraturePoint& point : quadrature) {
                const Vec3 position = (1.0 - point.a - point.b) * triangle[0]
                    + point.a * triangle[1] + point.b * triangle[2];
                const double weight = point.weight * jacobian;
                const auto leftLambda = lambdaOnTetFace(
                    position, left, face.leftLocal, mesh);
                const auto rightLambda = lambdaOnTetFace(
                    position, right, face.rightLocal, mesh);
                const auto leftShape = shapeP2(leftLambda);
                const auto rightShape = shapeP2(rightLambda);
                const auto leftGradientShape = gradShapeP2(leftLambda, leftGeometry);
                const auto rightGradientShape = gradShapeP2(rightLambda, rightGeometry);
                double leftTemperature = 0.0;
                double rightTemperature = 0.0;
                Vec3 leftGradient{};
                Vec3 rightGradient{};
                for (int local = 0; local < 10; ++local) {
                    const double leftValue = temperature[static_cast<std::size_t>(
                        left.dof[static_cast<std::size_t>(local)])];
                    const double rightValue = temperature[static_cast<std::size_t>(
                        right.dof[static_cast<std::size_t>(local)])];
                    leftTemperature += leftShape[static_cast<std::size_t>(local)] * leftValue;
                    rightTemperature += rightShape[static_cast<std::size_t>(local)] * rightValue;
                    leftGradient = leftGradient
                        + leftValue * leftGradientShape[static_cast<std::size_t>(local)];
                    rightGradient = rightGradient
                        + rightValue * rightGradientShape[static_cast<std::size_t>(local)];
                }
                const double leftOutwardFlux = -(
                    leftMaterial.conductivityX * leftGradient.x * face.leftNormal.x
                    + leftMaterial.conductivityY * leftGradient.y * face.leftNormal.y
                    + leftMaterial.conductivityZ * leftGradient.z * face.leftNormal.z);
                const double rightOutwardFlux = -(
                    rightMaterial.conductivityX * rightGradient.x * face.rightNormal.x
                    + rightMaterial.conductivityY * rightGradient.y * face.rightNormal.y
                    + rightMaterial.conductivityZ * rightGradient.z * face.rightNormal.z);
                const double mismatch = leftOutwardFlux + rightOutwardFlux;
                area += weight;
                const double jump = leftTemperature - rightTemperature;
                jumpSquared += jump * jump * weight;
                fluxMismatchSquared += mismatch * mismatch * weight;
                fluxScaleSquared += 0.5 * (leftOutwardFlux * leftOutwardFlux
                    + rightOutwardFlux * rightOutwardFlux) * weight;
            }
        }
    }
    InterfacePhysicsMetrics result;
    result.temperatureJumpRms = std::sqrt(jumpSquared / std::max(1.0e-300, area));
    result.relativeFluxImbalance = std::sqrt(fluxMismatchSquared)
        / std::max(1.0e-300, std::sqrt(fluxScaleSquared));
    return result;
}

namespace {

DetailedInterfacePhysicsMetrics calculateDetailedInterfacePhysicsMetricsImpl(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::vector<double>& temperature,
    const std::pair<int, int>* selectedSubdomains)
{
    if (temperature.size() != mesh.nodes.size()) {
        throw std::runtime_error("[MOR] Interface diagnostic temperature has the wrong size.");
    }
    DetailedInterfacePhysicsMetrics result;
    const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
    double totalArea = 0.0;
    double totalJumpSquared = 0.0;
    double totalMismatchSquared = 0.0;
    double totalLeftFluxSquared = 0.0;
    double totalRightFluxSquared = 0.0;
    result.interfaces.reserve(mesh.interfaceFaces.size());
    for (std::size_t faceIndex = 0; faceIndex < mesh.interfaceFaces.size(); ++faceIndex) {
        const InterfaceFace& face = mesh.interfaceFaces[faceIndex];
        const Tet& left = mesh.tets[static_cast<std::size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<std::size_t>(face.rightTet)];
        if (selectedSubdomains != nullptr) {
            const std::pair<int, int> adjacent{
                std::min(left.subdomain, right.subdomain),
                std::max(left.subdomain, right.subdomain)};
            if (adjacent != *selectedSubdomains) continue;
        }
        const ElementGeometry leftGeometry = elementGeometry(mesh, left);
        const ElementGeometry rightGeometry = elementGeometry(mesh, right);
        const Material& leftMaterial = materialForTet(physics, left);
        const Material& rightMaterial = materialForTet(physics, right);
        const SipgInterfaceScales scales = sipgScales(
            mesh, physics, face, left, right, leftMaterial, rightMaterial);
        const double penalty = scales.penalty;
        const Vec3 normal = norm(face.leftNormal) > 1.0e-30
            ? face.leftNormal : normalized(face.rightNormal * -1.0);

        InterfaceFluxRecord record;
        record.facePairId = static_cast<int>(faceIndex);
        record.leftSubdomain = left.subdomain;
        record.rightSubdomain = right.subdomain;
        for (std::size_t interfaceIndex = 0;
             interfaceIndex < mesh.interfaceSummaries.size(); ++interfaceIndex) {
            const InterfaceBuildSummary& summary = mesh.interfaceSummaries[interfaceIndex];
            if ((summary.leftSubdomain == record.leftSubdomain
                    && summary.rightSubdomain == record.rightSubdomain)
                || (summary.leftSubdomain == record.rightSubdomain
                    && summary.rightSubdomain == record.leftSubdomain)) {
                record.interfaceId = static_cast<int>(interfaceIndex);
                break;
            }
        }
        record.leftBoundaryEntity = face.leftBoundaryEntity;
        record.rightBoundaryEntity = face.rightBoundaryEntity;
        double jumpSquared = 0.0;
        double mismatchSquared = 0.0;
        double leftFluxSquared = 0.0;
        double rightFluxSquared = 0.0;
        double leftFluxIntegral = 0.0;
        double rightFluxIntegral = 0.0;
        double numericalFluxIntegral = 0.0;
        for (std::size_t triangleIndex = 0;
             triangleIndex < face.integrationTriangles.size(); ++triangleIndex) {
            const auto& triangle = face.integrationTriangles[triangleIndex];
            const double jacobian = norm(cross(
                triangle[1] - triangle[0], triangle[2] - triangle[0]));
            double triangleArea = 0.0;
            double triangleJumpSquared = 0.0;
            double triangleMismatchSquared = 0.0;
            double triangleLeftFluxSquared = 0.0;
            double triangleRightFluxSquared = 0.0;
            double triangleLeftFluxIntegral = 0.0;
            double triangleRightFluxIntegral = 0.0;
            double triangleNumericalFluxIntegral = 0.0;
            double triangleMaximumMismatch = 0.0;
            for (const TriangleQuadraturePoint& point : quadrature) {
                const Vec3 position = (1.0 - point.a - point.b) * triangle[0]
                    + point.a * triangle[1] + point.b * triangle[2];
                const double weight = point.weight * jacobian;
                const auto leftLambda = lambdaOnTetFace(
                    position, left, face.leftLocal, mesh);
                const auto rightLambda = lambdaOnTetFace(
                    position, right, face.rightLocal, mesh);
                const auto leftShape = shapeP2(leftLambda);
                const auto rightShape = shapeP2(rightLambda);
                const auto leftGradientShape = gradShapeP2(leftLambda, leftGeometry);
                const auto rightGradientShape = gradShapeP2(rightLambda, rightGeometry);
                double leftTemperature = 0.0;
                double rightTemperature = 0.0;
                Vec3 leftGradient{};
                Vec3 rightGradient{};
                for (int local = 0; local < 10; ++local) {
                    const double leftValue = temperature[static_cast<std::size_t>(
                        left.dof[static_cast<std::size_t>(local)])];
                    const double rightValue = temperature[static_cast<std::size_t>(
                        right.dof[static_cast<std::size_t>(local)])];
                    leftTemperature += leftShape[static_cast<std::size_t>(local)] * leftValue;
                    rightTemperature += rightShape[static_cast<std::size_t>(local)] * rightValue;
                    leftGradient = leftGradient
                        + leftValue * leftGradientShape[static_cast<std::size_t>(local)];
                    rightGradient = rightGradient
                        + rightValue * rightGradientShape[static_cast<std::size_t>(local)];
                }
                const double leftDirectional =
                    leftMaterial.conductivityX * leftGradient.x * normal.x
                    + leftMaterial.conductivityY * leftGradient.y * normal.y
                    + leftMaterial.conductivityZ * leftGradient.z * normal.z;
                const double rightDirectional =
                    rightMaterial.conductivityX * rightGradient.x * normal.x
                    + rightMaterial.conductivityY * rightGradient.y * normal.y
                    + rightMaterial.conductivityZ * rightGradient.z * normal.z;
                const double leftOutwardFlux = -leftDirectional;
                const double rightOutwardFlux = rightDirectional;
                const double jump = leftTemperature - rightTemperature;
                const double mismatch = leftOutwardFlux + rightOutwardFlux;
                const double numericalFlux =
                    -(scales.weightLeft * leftDirectional
                      + scales.weightRight * rightDirectional)
                    + penalty * jump;
                record.area += weight;
                jumpSquared += jump * jump * weight;
                mismatchSquared += mismatch * mismatch * weight;
                leftFluxSquared += leftOutwardFlux * leftOutwardFlux * weight;
                rightFluxSquared += rightOutwardFlux * rightOutwardFlux * weight;
                leftFluxIntegral += leftOutwardFlux * weight;
                rightFluxIntegral += rightOutwardFlux * weight;
                numericalFluxIntegral += numericalFlux * weight;
                triangleArea += weight;
                triangleJumpSquared += jump * jump * weight;
                triangleMismatchSquared += mismatch * mismatch * weight;
                triangleLeftFluxSquared += leftOutwardFlux * leftOutwardFlux * weight;
                triangleRightFluxSquared += rightOutwardFlux * rightOutwardFlux * weight;
                triangleLeftFluxIntegral += leftOutwardFlux * weight;
                triangleRightFluxIntegral += rightOutwardFlux * weight;
                triangleNumericalFluxIntegral += numericalFlux * weight;
                triangleMaximumMismatch = std::max(
                    triangleMaximumMismatch, std::abs(mismatch));
                if (std::abs(mismatch) > record.maximumFluxImbalance) {
                    record.maximumFluxImbalance = std::abs(mismatch);
                    record.worstIntegrationTriangle = static_cast<int>(triangleIndex);
                }
            }
            InterfaceTriangleFluxRecord triangleRecord;
            triangleRecord.interfaceId = record.interfaceId;
            triangleRecord.facePairId = record.facePairId;
            triangleRecord.integrationTriangleId = static_cast<int>(triangleIndex);
            triangleRecord.leftSubdomain = record.leftSubdomain;
            triangleRecord.rightSubdomain = record.rightSubdomain;
            triangleRecord.leftBoundaryEntity = record.leftBoundaryEntity;
            triangleRecord.rightBoundaryEntity = record.rightBoundaryEntity;
            triangleRecord.area = triangleArea;
            const double inverseTriangleArea =
                1.0 / std::max(1.0e-300, triangleArea);
            triangleRecord.temperatureJumpRms =
                std::sqrt(triangleJumpSquared * inverseTriangleArea);
            triangleRecord.leftPhysicalNormalFlux =
                triangleLeftFluxIntegral * inverseTriangleArea;
            triangleRecord.rightPhysicalNormalFlux =
                triangleRightFluxIntegral * inverseTriangleArea;
            triangleRecord.sipgNumericalFlux =
                triangleNumericalFluxIntegral * inverseTriangleArea;
            triangleRecord.sipgPenalty = penalty;
            triangleRecord.fluxImbalanceL2 =
                std::sqrt(triangleMismatchSquared * inverseTriangleArea);
            const double triangleFluxScale = std::max({
                std::sqrt(triangleLeftFluxSquared),
                std::sqrt(triangleRightFluxSquared),
                result.relativeFluxFloor * std::sqrt(triangleArea)});
            triangleRecord.relativeFluxImbalance =
                std::sqrt(triangleMismatchSquared)
                / std::max(1.0e-300, triangleFluxScale);
            triangleRecord.maximumFluxImbalance = triangleMaximumMismatch;
            result.triangles.push_back(std::move(triangleRecord));
        }
        const double inverseArea = 1.0 / std::max(1.0e-300, record.area);
        record.temperatureJumpRms = std::sqrt(jumpSquared * inverseArea);
        record.leftPhysicalNormalFlux = leftFluxIntegral * inverseArea;
        record.rightPhysicalNormalFlux = rightFluxIntegral * inverseArea;
        record.sipgNumericalFlux = numericalFluxIntegral * inverseArea;
        record.fluxImbalanceL2 = std::sqrt(mismatchSquared * inverseArea);
        const double localFluxScale = std::max({
            std::sqrt(leftFluxSquared), std::sqrt(rightFluxSquared),
            result.relativeFluxFloor * std::sqrt(record.area)});
        record.relativeFluxImbalance = std::sqrt(mismatchSquared)
            / std::max(1.0e-300, localFluxScale);
        totalArea += record.area;
        totalJumpSquared += jumpSquared;
        totalMismatchSquared += mismatchSquared;
        totalLeftFluxSquared += leftFluxSquared;
        totalRightFluxSquared += rightFluxSquared;
        if (record.maximumFluxImbalance > result.maximumFluxImbalance) {
            result.maximumFluxImbalance = record.maximumFluxImbalance;
            result.worstInterfaceId = record.interfaceId;
            result.worstFacePairId = record.facePairId;
            result.worstIntegrationTriangle = record.worstIntegrationTriangle;
        }
        result.interfaces.push_back(record);
    }
    result.aggregate.temperatureJumpRms = std::sqrt(
        totalJumpSquared / std::max(1.0e-300, totalArea));
    const double globalFluxScale = std::max({
        std::sqrt(totalLeftFluxSquared), std::sqrt(totalRightFluxSquared),
        result.relativeFluxFloor * std::sqrt(totalArea)});
    result.aggregate.relativeFluxImbalance = std::sqrt(totalMismatchSquared)
        / std::max(1.0e-300, globalFluxScale);
    result.areaWeightedFluxImbalanceL2 = std::sqrt(
        totalMismatchSquared / std::max(1.0e-300, totalArea));
    return result;
}

} // namespace

DetailedInterfacePhysicsMetrics calculateDetailedInterfacePhysicsMetrics(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::vector<double>& temperature)
{
    return calculateDetailedInterfacePhysicsMetricsImpl(
        mesh, physics, temperature, nullptr);
}

DetailedInterfacePhysicsMetrics
calculateDetailedInterfacePhysicsMetricsForSubdomains(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::vector<double>& temperature,
    int firstSubdomain,
    int secondSubdomain)
{
    const std::pair<int, int> selected =
        std::minmax(firstSubdomain, secondSubdomain);
    return calculateDetailedInterfacePhysicsMetricsImpl(
        mesh, physics, temperature, &selected);
}

} // namespace mor
