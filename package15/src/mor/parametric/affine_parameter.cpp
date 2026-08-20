#include "affine_parameter.hpp"

#include "sipg_core.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <set>
#include <stdexcept>

namespace mor::parametric {
namespace {

class Hash64 {
public:
    template <typename T>
    void add(const T& value)
    {
        const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
        for (std::size_t i = 0; i < sizeof(T); ++i) {
            value_ ^= static_cast<std::uint64_t>(bytes[i]);
            value_ *= UINT64_C(1099511628211);
        }
    }

    void add(const std::string& value)
    {
        for (unsigned char byte : value) {
            value_ ^= byte;
            value_ *= UINT64_C(1099511628211);
        }
    }

    std::uint64_t value() const { return value_; }

private:
    std::uint64_t value_ = UINT64_C(1469598103934665603);
};

bool closeValue(double left, double right)
{
    return std::abs(left - right)
        <= 1.0e-12 * std::max({1.0, std::abs(left), std::abs(right)});
}

} // namespace

bool tetSelected(const AffineParameter& parameter,
                 int subdomain,
                 int domainEntity)
{
    return parameter.name == "material-k"
        && (parameter.subdomain < 0 || parameter.subdomain == subdomain)
        && parameter.regionId == domainEntity;
}

bool boundarySelected(const AffineParameter& parameter,
                      int subdomain,
                      int boundaryEntity)
{
    return parameter.name == "convection-h"
        && (parameter.subdomain < 0 || parameter.subdomain == subdomain)
        && parameter.regionId == boundaryEntity;
}

double harmonicTheta(const AffineParameter& parameter,
                     double value,
                     std::size_t group)
{
    if (group >= parameter.harmonicNeighborConductivities.size()) {
        return 0.0;
    }
    const double fixed = parameter.harmonicNeighborConductivities[group];
    return 2.0 * value * fixed / (value + fixed);
}

AffineParameter resolveAffineParameter(const Mesh& mesh,
                                       const CaseConfig& physics,
                                       const mor::Options& options)
{
    if (options.matrixParameter != "convection-h"
        && options.matrixParameter != "material-k") {
        throw std::runtime_error(
            "--mor-matrix-parameter must be convection-h or material-k.");
    }
    if (options.parameterRegionId < 0) {
        throw std::runtime_error(
            "Stage 2B.1 requires --mor-parameter-region-id <id>.");
    }
    if (!std::isfinite(options.parameterMinimum)
        || !std::isfinite(options.parameterMaximum)
        || !std::isfinite(options.parameterReference)) {
        throw std::runtime_error(
            "Stage 2B.1 requires explicit --mor-parameter-min, --mor-parameter-max, and --mor-parameter-reference.");
    }
    if (!(options.parameterMinimum > 0.0)
        || !(options.parameterMaximum > options.parameterMinimum)
        || options.parameterReference < options.parameterMinimum
        || options.parameterReference > options.parameterMaximum) {
        throw std::runtime_error("The Stage 2B.1 physical parameter interval is invalid.");
    }

    AffineParameter result;
    result.name = options.matrixParameter;
    result.units = result.name == "convection-h" ? "W/(m^2 K)" : "W/(m K)";
    result.subdomain = options.parameterSubdomain;
    result.regionId = options.parameterRegionId;
    result.minimum = options.parameterMinimum;
    result.maximum = options.parameterMaximum;
    result.reference = options.parameterReference;

    if (result.name == "convection-h") {
        for (const ConvectionCondition& condition : physics.convectionConditions) {
            if (!boundarySelected(result, condition.subdomain, condition.boundaryEntity)) {
                continue;
            }
            if (result.selectedBoundaryFaceCount == 0) {
                result.ambientTemperature = condition.ambientTemperature;
            } else if (!closeValue(result.ambientTemperature, condition.ambientTemperature)) {
                throw std::runtime_error(
                    "Selected convection conditions must share one ambient temperature.");
            }
            if (!closeValue(condition.coefficient, result.reference)) {
                throw std::runtime_error(
                    "--mor-parameter-reference must match the selected convection coefficient in the case configuration.");
            }
            ++result.selectedBoundaryFaceCount;
        }
        if (result.selectedBoundaryFaceCount == 0) {
            throw std::runtime_error(
                "No configured convection condition matches the selected subdomain/boundary ID.");
        }
        int meshFaces = 0;
        for (const BoundaryFace& face : mesh.boundaryFaces) {
            if (boundarySelected(result, face.subdomain, face.boundaryEntity)) {
                ++meshFaces;
            }
        }
        result.selectedBoundaryFaceCount = meshFaces;
        if (meshFaces == 0) {
            throw std::runtime_error("The selected convection boundary has no mesh faces.");
        }
    } else {
        std::set<int> selectedTetIndices;
        double configuredReference = -1.0;
        for (std::size_t index = 0; index < mesh.tets.size(); ++index) {
            const Tet& tet = mesh.tets[index];
            if (!tetSelected(result, tet.subdomain, tet.domainEntity)) {
                continue;
            }
            const Material& material = materialForTet(physics, tet);
            if (!closeValue(material.conductivityX, material.conductivityY)
                || !closeValue(material.conductivityX, material.conductivityZ)) {
                throw std::runtime_error(
                    "Stage 2B.1 material-k currently requires an isotropic selected material.");
            }
            if (configuredReference < 0.0) {
                configuredReference = material.conductivityX;
            } else if (!closeValue(configuredReference, material.conductivityX)) {
                throw std::runtime_error(
                    "Selected material regions do not share one reference conductivity.");
            }
            selectedTetIndices.insert(static_cast<int>(index));
        }
        result.selectedTetCount = static_cast<int>(selectedTetIndices.size());
        if (result.selectedTetCount == 0) {
            throw std::runtime_error(
                "No tetrahedra match the selected material subdomain/region ID.");
        }
        if (!closeValue(configuredReference, result.reference)) {
            throw std::runtime_error(
                "--mor-parameter-reference must match the selected material conductivity in the case configuration.");
        }
        std::vector<double> neighborConductivities;
        for (const InterfaceFace& face : mesh.interfaceFaces) {
            const Tet& left = mesh.tets[static_cast<std::size_t>(face.leftTet)];
            const Tet& right = mesh.tets[static_cast<std::size_t>(face.rightTet)];
            const bool selectedLeft = selectedTetIndices.count(face.leftTet) != 0;
            const bool selectedRight = selectedTetIndices.count(face.rightTet) != 0;
            if (!selectedLeft && !selectedRight) {
                continue;
            }
            result.touchesInterface = true;
            if (selectedLeft && selectedRight) {
                continue;
            }
            const Material& fixed = materialForTet(
                physics, selectedLeft ? right : left);
            if (!closeValue(fixed.conductivityX, fixed.conductivityY)
                || !closeValue(fixed.conductivityX, fixed.conductivityZ)) {
                throw std::runtime_error(
                    "Material-k SIPG affine decomposition requires an isotropic fixed neighbor material.");
            }
            bool found = false;
            for (double existing : neighborConductivities) {
                if (closeValue(existing, fixed.conductivityX)) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                neighborConductivities.push_back(fixed.conductivityX);
            }
        }
        std::sort(neighborConductivities.begin(), neighborConductivities.end());
        result.harmonicNeighborConductivities = std::move(neighborConductivities);
        if (result.touchesInterface && physics.penaltyMode != "harmonic") {
            throw std::runtime_error(
                "Material-k Stage 2B.1 interface parameterization requires harmonic SIPG penalty mode.");
        }
    }

    Hash64 hash;
    hash.add(result.name);
    hash.add(result.units);
    hash.add(result.subdomain);
    hash.add(result.regionId);
    hash.add(result.minimum);
    hash.add(result.maximum);
    hash.add(result.reference);
    hash.add(result.ambientTemperature);
    for (double conductivity : result.harmonicNeighborConductivities) {
        hash.add(conductivity);
    }
    hash.add(result.touchesInterface);
    hash.add(result.selectedTetCount);
    hash.add(result.selectedBoundaryFaceCount);
    result.definitionHash = hash.value();
    return result;
}

CaseConfig physicsAtParameter(const CaseConfig& physics,
                              const AffineParameter& parameter,
                              double value)
{
    if (!(value > 0.0) || !std::isfinite(value)) {
        throw std::runtime_error("The Stage 2B.1 physical parameter must be positive and finite.");
    }
    CaseConfig result = physics;
    if (parameter.name == "convection-h") {
        int matches = 0;
        for (ConvectionCondition& condition : result.convectionConditions) {
            if (boundarySelected(parameter, condition.subdomain, condition.boundaryEntity)) {
                condition.coefficient = value;
                ++matches;
            }
        }
        if (matches == 0) {
            throw std::runtime_error("Selected convection parameter disappeared from the configuration.");
        }
        return result;
    }

    for (std::size_t domain = 0; domain < result.domains.size(); ++domain) {
        if (parameter.subdomain >= 0
            && parameter.subdomain != static_cast<int>(domain)) {
            continue;
        }
        auto existing = result.domains[domain].materialsByDomainEntity.find(
            parameter.regionId);
        Material selected = existing != result.domains[domain].materialsByDomainEntity.end()
            ? existing->second : result.domains[domain].material;
        setIsotropicConductivity(selected, value);
        result.domains[domain].materialsByDomainEntity[parameter.regionId] = selected;
    }
    return result;
}

} // namespace mor::parametric
