#include "affine_fem_components.hpp"

#include "config_io.hpp"
#include "mesh_loader.hpp"
#include "fem_assembly.hpp"
#include "linear_solvers.hpp"
#include "diagnostics_io.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

namespace mor::parametric {
namespace {

void assembleUnitMaterialVolume(const Mesh& mesh,
                                const AffineParameter& parameter,
                                SparseMatrix& component)
{
    const std::vector<TetQuadraturePoint> quadrature = makeTetQuadrature();
    for (const Tet& tet : mesh.tets) {
        if (!tetSelected(parameter, tet.subdomain, tet.domainEntity)) {
            continue;
        }
        const ElementGeometry geometry = elementGeometry(mesh, tet);
        std::array<double, 100> local{};
        for (const TetQuadraturePoint& point : quadrature) {
            const double weight = point.weight * geometry.detJ;
            const auto gradient = physicalGradP2(point, geometry);
            for (int row = 0; row < 10; ++row) {
                for (int column = 0; column < 10; ++column) {
                    local[static_cast<std::size_t>(row * 10 + column)] +=
                        dot(gradient[static_cast<std::size_t>(row)],
                            gradient[static_cast<std::size_t>(column)]) * weight;
                }
            }
        }
        for (int row = 0; row < 10; ++row) {
            for (int column = 0; column < 10; ++column) {
                component.add(tet.dof[static_cast<std::size_t>(row)],
                              tet.dof[static_cast<std::size_t>(column)],
                              local[static_cast<std::size_t>(row * 10 + column)]);
            }
        }
    }
}

void zeroAllConductivities(CaseConfig& physics)
{
    for (DomainConfig& domain : physics.domains) {
        setIsotropicConductivity(domain.material, 0.0);
        for (auto& entry : domain.materialsByDomainEntity) {
            setIsotropicConductivity(entry.second, 0.0);
        }
    }
}

void assembleUnitMaterialConsistency(const Mesh& mesh,
                                     const CaseConfig& physics,
                                     const AffineParameter& parameter,
                                     SparseMatrix& component)
{
    CaseConfig unit = physics;
    zeroAllConductivities(unit);
    for (std::size_t domain = 0; domain < unit.domains.size(); ++domain) {
        if (parameter.subdomain >= 0
            && parameter.subdomain != static_cast<int>(domain)) {
            continue;
        }
        Material material;
        material.name = "Stage2B1_unit_k";
        setIsotropicConductivity(material, 1.0);
        unit.domains[domain].materialsByDomainEntity[parameter.regionId] = material;
    }
    assembleSipgInterface(mesh, unit, component, false, true, nullptr);
}

void assembleUnitMaterialPenalty(const Mesh& mesh,
                                 const CaseConfig& physics,
                                 const AffineParameter& parameter,
                                 SparseMatrix& linearComponent,
                                 std::vector<SparseMatrix>& harmonicComponents)
{
    const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
    constexpr double polynomialOrder = 2.0;
    const double pScale = physics.penaltyScaling == "p1_squared"
        ? (polynomialOrder + 1.0) * (polynomialOrder + 1.0)
        : polynomialOrder * (polynomialOrder + 1.0);
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        const Tet& left = mesh.tets[static_cast<std::size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<std::size_t>(face.rightTet)];
        const bool selectedLeft = tetSelected(
            parameter, left.subdomain, left.domainEntity);
        const bool selectedRight = tetSelected(
            parameter, right.subdomain, right.domainEntity);
        if (!selectedLeft && !selectedRight) {
            continue;
        }
        SparseMatrix* target = &linearComponent;
        if (!(selectedLeft && selectedRight)) {
            const Material& fixed = materialForTet(
                physics, selectedLeft ? right : left);
            std::size_t group = 0;
            for (; group < parameter.harmonicNeighborConductivities.size(); ++group) {
                const double candidate = parameter.harmonicNeighborConductivities[group];
                if (std::abs(candidate - fixed.conductivityX)
                    <= 1.0e-12 * std::max({1.0, std::abs(candidate),
                                           std::abs(fixed.conductivityX)})) {
                    break;
                }
            }
            if (group >= harmonicComponents.size()) {
                throw std::runtime_error("SIPG harmonic affine group lookup failed.");
            }
            target = &harmonicComponents[group];
        }
        const ElementGeometry leftGeometry = elementGeometry(mesh, left);
        const ElementGeometry rightGeometry = elementGeometry(mesh, right);
        const std::array<Vec3, 3> leftPhysicalFace{{
            mesh.nodes[static_cast<std::size_t>(left.v[static_cast<std::size_t>(face.leftLocal[0])])].p,
            mesh.nodes[static_cast<std::size_t>(left.v[static_cast<std::size_t>(face.leftLocal[1])])].p,
            mesh.nodes[static_cast<std::size_t>(left.v[static_cast<std::size_t>(face.leftLocal[2])])].p}};
        const std::array<Vec3, 3> rightPhysicalFace{{
            mesh.nodes[static_cast<std::size_t>(right.v[static_cast<std::size_t>(face.rightLocal[0])])].p,
            mesh.nodes[static_cast<std::size_t>(right.v[static_cast<std::size_t>(face.rightLocal[1])])].p,
            mesh.nodes[static_cast<std::size_t>(right.v[static_cast<std::size_t>(face.rightLocal[2])])].p}};
        const double areaSum = triangleArea(leftPhysicalFace)
            + triangleArea(rightPhysicalFace);
        const double volumeSum = std::abs(leftGeometry.detJ) / 6.0
            + std::abs(rightGeometry.detJ) / 6.0;
        const double hFace = volumeSum / std::max(1.0e-30, areaSum);
        const double unitPenalty = physics.penaltyFactor * pScale
            / std::max(1.0e-30, hFace);
        std::array<double, 400> local{};
        for (const auto& triangle : face.integrationTriangles) {
            const double jacobian = norm(cross(
                triangle[1] - triangle[0], triangle[2] - triangle[0]));
            for (const TriangleQuadraturePoint& point : quadrature) {
                const double weight = point.weight * jacobian;
                const Vec3 position = (1.0 - point.a - point.b) * triangle[0]
                    + point.a * triangle[1] + point.b * triangle[2];
                const auto leftShape = shapeP2(lambdaOnTetFace(
                    position, left, face.leftLocal, mesh));
                const auto rightShape = shapeP2(lambdaOnTetFace(
                    position, right, face.rightLocal, mesh));
                for (int test = 0; test < 10; ++test) {
                    for (int trial = 0; trial < 10; ++trial) {
                        const double vLeft = leftShape[static_cast<std::size_t>(test)];
                        const double vRight = rightShape[static_cast<std::size_t>(test)];
                        const double uLeft = leftShape[static_cast<std::size_t>(trial)];
                        const double uRight = rightShape[static_cast<std::size_t>(trial)];
                        local[static_cast<std::size_t>(test * 20 + trial)] +=
                            unitPenalty * uLeft * vLeft * weight;
                        local[static_cast<std::size_t>(test * 20 + 10 + trial)] -=
                            unitPenalty * uRight * vLeft * weight;
                        local[static_cast<std::size_t>((10 + test) * 20 + trial)] -=
                            unitPenalty * uLeft * vRight * weight;
                        local[static_cast<std::size_t>((10 + test) * 20 + 10 + trial)] +=
                            unitPenalty * uRight * vRight * weight;
                    }
                }
            }
        }
        std::array<int, 20> dofs{};
        for (int localDof = 0; localDof < 10; ++localDof) {
            dofs[static_cast<std::size_t>(localDof)] =
                left.dof[static_cast<std::size_t>(localDof)];
            dofs[static_cast<std::size_t>(10 + localDof)] =
                right.dof[static_cast<std::size_t>(localDof)];
        }
        for (int row = 0; row < 20; ++row) {
            for (int column = 0; column < 20; ++column) {
                target->add(dofs[static_cast<std::size_t>(row)],
                            dofs[static_cast<std::size_t>(column)],
                            local[static_cast<std::size_t>(row * 20 + column)]);
            }
        }
    }
}

void makeHomogeneousDirichletComponent(const Mesh& mesh,
                                       SparseMatrix& matrix,
                                       std::vector<double>& rhs)
{
    const std::vector<double> fixed = makeDirichletAdjustedSystem(mesh, matrix);
    std::vector<MatrixEntry> freeEntries;
    matrix.forEachEntry([&](int row, int column, double value) {
        if (!mesh.nodes[static_cast<std::size_t>(row)].dirichlet
            && !mesh.nodes[static_cast<std::size_t>(column)].dirichlet) {
            freeEntries.push_back({row, column, value});
        }
    });
    matrix = SparseMatrix(static_cast<int>(mesh.nodes.size()));
    matrix.appendEntries(freeEntries);
    for (std::size_t row = 0; row < rhs.size(); ++row) {
        if (mesh.nodes[row].dirichlet) {
            rhs[row] = 0.0;
        } else {
            rhs[row] -= fixed[row];
        }
    }
    matrix.finalizeCsr();
}

std::uint64_t hashValues(const std::vector<double>& values)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    for (double value : values) {
        const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
        for (std::size_t i = 0; i < sizeof(value); ++i) {
            hash ^= bytes[i];
            hash *= UINT64_C(1099511628211);
        }
    }
    return hash;
}

std::uint64_t hashMatrix(const SparseMatrix& matrix)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    matrix.forEachEntry([&](int row, int column, double value) {
        const auto add = [&](const auto& item) {
            const auto* bytes = reinterpret_cast<const unsigned char*>(&item);
            for (std::size_t i = 0; i < sizeof(item); ++i) {
                hash ^= bytes[i];
                hash *= UINT64_C(1099511628211);
            }
        };
        add(row);
        add(column);
        add(value);
    });
    return hash;
}

double vectorRelativeDifference(const std::vector<double>& left,
                                const std::vector<double>& right)
{
    double difference = 0.0;
    double reference = 0.0;
    for (std::size_t i = 0; i < left.size(); ++i) {
        const double delta = left[i] - right[i];
        difference += delta * delta;
        reference += left[i] * left[i];
    }
    return std::sqrt(difference) / std::max(1.0e-300, std::sqrt(reference));
}

double matrixRelativeDifference(const SparseMatrix& left,
                                const SparseMatrix& right)
{
    if (!left.csrReady || !right.csrReady || left.size() != right.size()) {
        throw std::runtime_error("Affine validation requires finalized matching matrices.");
    }
    double difference = 0.0;
    double reference = 0.0;
    for (int row = 0; row < left.size(); ++row) {
        int li = left.rowPtr[static_cast<std::size_t>(row)];
        int ri = right.rowPtr[static_cast<std::size_t>(row)];
        const int lend = left.rowPtr[static_cast<std::size_t>(row + 1)];
        const int rend = right.rowPtr[static_cast<std::size_t>(row + 1)];
        while (li < lend || ri < rend) {
            const int lc = li < lend ? left.colInd[static_cast<std::size_t>(li)]
                                     : std::numeric_limits<int>::max();
            const int rc = ri < rend ? right.colInd[static_cast<std::size_t>(ri)]
                                     : std::numeric_limits<int>::max();
            double lv = 0.0;
            double rv = 0.0;
            if (lc <= rc) {
                lv = left.values[static_cast<std::size_t>(li++)];
            }
            if (rc <= lc) {
                rv = right.values[static_cast<std::size_t>(ri++)];
            }
            const double delta = lv - rv;
            difference += delta * delta;
            reference += lv * lv;
        }
    }
    return std::sqrt(difference) / std::max(1.0e-300, std::sqrt(reference));
}

double symmetryError(const SparseMatrix& matrix)
{
    double numerator = 0.0;
    double denominator = 0.0;
    for (int row = 0; row < matrix.size(); ++row) {
        for (int offset = matrix.rowPtr[static_cast<std::size_t>(row)];
             offset < matrix.rowPtr[static_cast<std::size_t>(row + 1)]; ++offset) {
            const int column = matrix.colInd[static_cast<std::size_t>(offset)];
            const double value = matrix.values[static_cast<std::size_t>(offset)];
            double transpose = 0.0;
            const int begin = matrix.rowPtr[static_cast<std::size_t>(column)];
            const int end = matrix.rowPtr[static_cast<std::size_t>(column + 1)];
            const auto found = std::lower_bound(
                matrix.colInd.begin() + begin, matrix.colInd.begin() + end, row);
            if (found != matrix.colInd.begin() + end && *found == row) {
                transpose = matrix.values[static_cast<std::size_t>(
                    std::distance(matrix.colInd.begin(), found))];
            }
            const double delta = value - transpose;
            numerator += delta * delta;
            denominator += value * value;
        }
    }
    return std::sqrt(numerator) / std::max(1.0e-300, std::sqrt(denominator));
}

} // namespace

DirectParametricSystem assembleDirectParametricSystem(
    const Mesh& mesh,
    const CaseConfig& referencePhysics,
    const AffineParameter& parameter,
    double value)
{
    CaseConfig physics = physicsAtParameter(referencePhysics, parameter, value);
    if (physics.dirichletMethod != "strong") {
        throw std::runtime_error(
            "Stage 2B.1 currently requires the existing strong Dirichlet path.");
    }
    DirectParametricSystem result;
    const int size = static_cast<int>(mesh.nodes.size());
    result.matrix = SparseMatrix(size);
    result.rawSource.assign(static_cast<std::size_t>(size), 0.0);
    assembleVolume(mesh, physics, nullptr, result.matrix, result.rawSource);
    if (physics.thermalSourceScale != 1.0) {
        for (double& item : result.rawSource) {
            item *= physics.thermalSourceScale;
        }
    }
    result.heatOnlySource = result.rawSource;
    assembleConvectionBoundaries(
        mesh, physics, result.matrix, result.rawSource, false, nullptr);
    assembleSipgInterface(mesh, physics, result.matrix, false, true, nullptr);
    assembleSipgInterface(mesh, physics, result.matrix, true, false, nullptr);
    result.fixedAdjust = makeDirichletAdjustedSystem(mesh, result.matrix);
    result.rhs = result.rawSource;
    applyDirichletRhs(mesh, result.fixedAdjust, result.rhs);
    result.matrix.finalizeCsr();
    return result;
}

AffineFemComponents buildAffineFemComponents(
    const Mesh& mesh,
    const CaseConfig& physics,
    const SparseMatrix& referenceSystem,
    const std::vector<double>& assembledSource,
    const std::vector<double>& heatOnlySource,
    const std::vector<double>& fixedAdjust,
    const mor::Options& options)
{
    const auto start = std::chrono::steady_clock::now();
    if (physics.dirichletMethod != "strong") {
        throw std::runtime_error(
            "Stage 2B.1 exact affine components currently support strong Dirichlet conditions.");
    }
    AffineFemComponents result;
    result.parameter = resolveAffineParameter(mesh, physics, options);
    const int size = referenceSystem.size();
    result.matrixLinear = SparseMatrix(size);
    result.matrixHarmonic.assign(
        result.parameter.harmonicNeighborConductivities.size(), SparseMatrix(size));
    result.rhsLinear.assign(static_cast<std::size_t>(size), 0.0);
    result.rhsHarmonic.assign(result.matrixHarmonic.size(),
        std::vector<double>(static_cast<std::size_t>(size), 0.0));

    if (result.parameter.name == "material-k") {
        assembleUnitMaterialVolume(mesh, result.parameter, result.matrixLinear);
        if (result.parameter.touchesInterface) {
            assembleUnitMaterialConsistency(
                mesh, physics, result.parameter, result.matrixLinear);
            assembleUnitMaterialPenalty(
                mesh, physics, result.parameter, result.matrixLinear,
                result.matrixHarmonic);
        }
    } else {
        CaseConfig unit = physics;
        unit.convectionConditions.clear();
        for (const ConvectionCondition& condition : physics.convectionConditions) {
            if (!boundarySelected(result.parameter,
                                  condition.subdomain,
                                  condition.boundaryEntity)) {
                continue;
            }
            ConvectionCondition selected = condition;
            selected.coefficient = 1.0;
            unit.convectionConditions.push_back(selected);
        }
        assembleConvectionBoundaries(
            mesh, unit, result.matrixLinear, result.rhsLinear, false, nullptr);
    }
    makeHomogeneousDirichletComponent(
        mesh, result.matrixLinear, result.rhsLinear);
    for (std::size_t group = 0; group < result.matrixHarmonic.size(); ++group) {
        makeHomogeneousDirichletComponent(
            mesh, result.matrixHarmonic[group], result.rhsHarmonic[group]);
    }

    result.matrixConstant = SparseMatrix(size);
    result.matrixConstant.appendScaledEntries(referenceSystem, 1.0);
    result.matrixConstant.appendScaledEntries(
        result.matrixLinear, -result.parameter.reference);
    for (std::size_t group = 0; group < result.matrixHarmonic.size(); ++group) {
        result.matrixConstant.appendScaledEntries(
            result.matrixHarmonic[group],
            -harmonicTheta(result.parameter, result.parameter.reference, group));
    }
    result.matrixConstant.finalizeCsr();

    result.sources = buildSourceParameterization(
        mesh, physics, assembledSource, heatOnlySource, fixedAdjust);
    result.rhsConstant = result.sources.referenceRhs;
    for (std::size_t row = 0; row < result.rhsConstant.size(); ++row) {
        result.rhsConstant[row] -= result.parameter.reference
            * result.rhsLinear[row];
        for (std::size_t group = 0; group < result.rhsHarmonic.size(); ++group) {
            result.rhsConstant[row] -= harmonicTheta(
                result.parameter, result.parameter.reference, group)
                * result.rhsHarmonic[group][row];
        }
    }

    result.constantMatrixHash = hashMatrix(result.matrixConstant);
    result.linearMatrixHash = hashMatrix(result.matrixLinear);
    for (const SparseMatrix& matrix : result.matrixHarmonic) {
        result.harmonicMatrixHashes.push_back(hashMatrix(matrix));
    }
    result.constantRhsHash = hashValues(result.rhsConstant);
    result.linearRhsHash = hashValues(result.rhsLinear);
    for (const std::vector<double>& rhs : result.rhsHarmonic) {
        result.harmonicRhsHashes.push_back(hashValues(rhs));
    }
    result.assemblySeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    return result;
}

SparseMatrix composeMatrix(const AffineFemComponents& components, double value)
{
    SparseMatrix result(components.matrixConstant.size());
    result.appendScaledEntries(components.matrixConstant, 1.0);
    result.appendScaledEntries(components.matrixLinear, value);
    for (std::size_t group = 0; group < components.matrixHarmonic.size(); ++group) {
        result.appendScaledEntries(components.matrixHarmonic[group],
            harmonicTheta(components.parameter, value, group));
    }
    result.finalizeCsr();
    return result;
}

std::vector<double> composeRhs(const AffineFemComponents& components,
                               double value,
                               const std::vector<double>& powersW)
{
    if (powersW.size() != components.sources.channels.size()) {
        throw std::runtime_error("Stage 2B.1 power vector dimension mismatch.");
    }
    std::vector<double> result = components.rhsConstant;
    for (std::size_t row = 0; row < result.size(); ++row) {
        result[row] += value * components.rhsLinear[row];
        for (std::size_t group = 0; group < components.rhsHarmonic.size(); ++group) {
            result[row] += harmonicTheta(components.parameter, value, group)
                * components.rhsHarmonic[group][row];
        }
    }
    for (std::size_t channel = 0; channel < powersW.size(); ++channel) {
        const std::vector<double>& load = components.sources.channels[channel].rhsPerWatt;
        for (std::size_t row = 0; row < result.size(); ++row) {
            result[row] += powersW[channel] * load[row];
        }
    }
    return result;
}

AffineValidationRow validateAffineOperator(
    const Mesh& mesh,
    const CaseConfig& physics,
    const AffineFemComponents& components,
    double value,
    const std::vector<double>& powersW)
{
    DirectParametricSystem direct = assembleDirectParametricSystem(
        mesh, physics, components.parameter, value);
    SparseMatrix affine = composeMatrix(components, value);
    const std::vector<double> affineRhs = composeRhs(components, value, powersW);

    std::vector<double> directRhs = direct.rhs;
    const std::vector<double> nominal = [&]() {
        std::vector<double> values;
        for (const SourceChannel& channel : components.sources.channels) {
            values.push_back(channel.nominalPowerW);
        }
        return values;
    }();
    for (std::size_t channel = 0; channel < powersW.size(); ++channel) {
        const double delta = powersW[channel] - nominal[channel];
        const std::vector<double>& load = components.sources.channels[channel].rhsPerWatt;
        for (std::size_t row = 0; row < directRhs.size(); ++row) {
            directRhs[row] += delta * load[row];
        }
    }

    AffineValidationRow row;
    row.parameterValue = value;
    row.matrixRelativeDifference = matrixRelativeDifference(direct.matrix, affine);
    row.rhsRelativeDifference = vectorRelativeDifference(directRhs, affineRhs);
    row.symmetryError = symmetryError(affine);
    row.minimumDiagonal = std::numeric_limits<double>::infinity();
    for (int index = 0; index < affine.size(); ++index) {
        row.minimumDiagonal = std::min(row.minimumDiagonal, affine.diagonal(index));
    }
    row.nonzeros = affine.values.size();
    return row;
}

void writeAffineValidation(const std::vector<AffineValidationRow>& rows,
                           const std::filesystem::path& path)
{
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot write affine operator validation table.");
    }
    out << "parameter_value,matrix_relative_difference,rhs_relative_difference,"
        << "symmetry_error,minimum_diagonal,nnz\n" << std::setprecision(17);
    for (const AffineValidationRow& row : rows) {
        out << row.parameterValue << ',' << row.matrixRelativeDifference << ','
            << row.rhsRelativeDifference << ',' << row.symmetryError << ','
            << row.minimumDiagonal << ',' << row.nonzeros << '\n';
    }
}

} // namespace mor::parametric
