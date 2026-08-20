#pragma once

// Read-only quality checks for local bases and projected operators.

#include <vector>

struct CaseConfig;
struct Mesh;

namespace mor {

struct InterfacePhysicsMetrics {
    double temperatureJumpRms = 0.0;
    double relativeFluxImbalance = 0.0;
};

struct InterfaceFluxRecord {
    int interfaceId = -1;
    int facePairId = -1;
    int leftSubdomain = -1;
    int rightSubdomain = -1;
    int leftBoundaryEntity = -1;
    int rightBoundaryEntity = -1;
    double area = 0.0;
    double temperatureJumpRms = 0.0;
    double leftPhysicalNormalFlux = 0.0;
    double rightPhysicalNormalFlux = 0.0;
    double sipgNumericalFlux = 0.0;
    double fluxImbalanceL2 = 0.0;
    double relativeFluxImbalance = 0.0;
    double maximumFluxImbalance = 0.0;
    int worstIntegrationTriangle = -1;
};

struct InterfaceTriangleFluxRecord {
    int interfaceId = -1;
    int facePairId = -1;
    int integrationTriangleId = -1;
    int leftSubdomain = -1;
    int rightSubdomain = -1;
    int leftBoundaryEntity = -1;
    int rightBoundaryEntity = -1;
    double area = 0.0;
    double temperatureJumpRms = 0.0;
    double leftPhysicalNormalFlux = 0.0;
    double rightPhysicalNormalFlux = 0.0;
    double sipgNumericalFlux = 0.0;
    double sipgPenalty = 0.0;
    double fluxImbalanceL2 = 0.0;
    double relativeFluxImbalance = 0.0;
    double maximumFluxImbalance = 0.0;
};

struct DetailedInterfacePhysicsMetrics {
    InterfacePhysicsMetrics aggregate;
    double areaWeightedFluxImbalanceL2 = 0.0;
    double maximumFluxImbalance = 0.0;
    int worstInterfaceId = -1;
    int worstFacePairId = -1;
    int worstIntegrationTriangle = -1;
    double relativeFluxFloor = 1.0e-12;
    std::vector<InterfaceFluxRecord> interfaces;
    std::vector<InterfaceTriangleFluxRecord> triangles;
};

InterfacePhysicsMetrics calculateInterfacePhysicsMetrics(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::vector<double>& temperature);

DetailedInterfacePhysicsMetrics calculateDetailedInterfacePhysicsMetrics(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::vector<double>& temperature);

// Evaluate the identical SIPG face-quadrature path, restricted to one
// unordered adjacent-subdomain pair.  This is used by operator audits that
// must not pay for every physical interface on each linearity probe.
DetailedInterfacePhysicsMetrics
calculateDetailedInterfacePhysicsMetricsForSubdomains(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::vector<double>& temperature,
    int firstSubdomain,
    int secondSubdomain);

} // namespace mor
