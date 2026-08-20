#include "../sipg_core.hpp"
#include "interface_operator.hpp"
#include "local_solver.hpp"
#include "schur_operator.hpp"

#include <algorithm>
#include <chrono>
#include <map>
#include <stdexcept>
#include <tuple>

namespace ddm_schur {

namespace {

struct CoarseVector {
    std::vector<int> indices;
    std::vector<double> values;
};

void appendInterfaceStencilDofs(const Tet& tet, std::vector<int>& dofs)
{
    dofs.insert(dofs.end(), tet.dof.begin(), tet.dof.end());
}

} // namespace

struct SchurOperator::Impl {
    InterfacePartition partition;
    SparseMatrix interfaceMatrix;
    std::vector<LocalSolver> interiorSolvers;
    std::vector<LocalSolver> blockSolvers;
    LocalSolver coarseSolver;
    std::vector<CoarseVector> coarseBasis;
    std::vector<CoarseVector> coarseImages;
    double factorSeconds = 0.0;
    double symbolicSeconds = 0.0;
    double numericalSeconds = 0.0;
    double localSeconds = 0.0;
    double coarseSeconds = 0.0;
    int solveCalls = 0;
    int symbolicCalls = 0;
    int numericalCalls = 0;
    int matvecs = 0;
    int patchCount = 0;
    std::size_t bytes = 0;

    Impl(const Mesh& mesh,
         const SparseMatrix& system,
         bool coarseLinearXY,
         bool coarseLinearZ,
         bool coarseGlobalQuadraticZ,
         bool coarseInterfacePatches,
         bool coarseInterfaceLinearXY)
        : partition(buildInterfacePartition(mesh, system)),
          interfaceMatrix(static_cast<int>(partition.interfaceGlobalDofs.size()))
    {
        for (const Entry& entry : partition.interfaceEntries) {
            interfaceMatrix.add(entry.row, entry.col, entry.value);
        }
        interfaceMatrix.finalizeCsr();

        interiorSolvers.reserve(partition.domains.size());
        blockSolvers.reserve(partition.domains.size());
        for (const DomainBlocks& domain : partition.domains) {
            interiorSolvers.emplace_back(
                static_cast<int>(domain.interiorGlobalDofs.size()), domain.interiorEntries);
            blockSolvers.emplace_back(
                static_cast<int>(domain.interiorGlobalDofs.size() + domain.interfaceGlobalDofs.size()),
                domain.fullBlockEntries);
        }
        buildCoarseCorrection(mesh,
                              coarseLinearXY,
                              coarseLinearZ,
                              coarseGlobalQuadraticZ,
                              coarseInterfacePatches,
                              coarseInterfaceLinearXY);
        for (const LocalSolver& solver : interiorSolvers) {
            symbolicSeconds += solver.symbolicAnalysisSeconds();
            numericalSeconds += solver.numericalFactorizationSeconds();
            symbolicCalls += solver.symbolicAnalysisCalls();
            numericalCalls += solver.numericalFactorizationCalls();
        }
        for (const LocalSolver& solver : blockSolvers) {
            symbolicSeconds += solver.symbolicAnalysisSeconds();
            numericalSeconds += solver.numericalFactorizationSeconds();
            symbolicCalls += solver.symbolicAnalysisCalls();
            numericalCalls += solver.numericalFactorizationCalls();
        }
        factorSeconds = symbolicSeconds + numericalSeconds;

        bytes = partitionMemoryBytes(partition)
            + interfaceMatrix.rowPtr.size() * sizeof(int)
            + interfaceMatrix.colInd.size() * sizeof(int)
            + interfaceMatrix.values.size() * sizeof(double);
        for (const LocalSolver& solver : interiorSolvers) {
            bytes += solver.memoryBytes();
        }
        for (const LocalSolver& solver : blockSolvers) {
            bytes += solver.memoryBytes();
        }
        bytes += coarseSolver.memoryBytes()
            + coarseBasis.size() * sizeof(CoarseVector);
        for (const CoarseVector& basis : coarseBasis) {
            bytes += basis.indices.size() * sizeof(int) + basis.values.size() * sizeof(double);
        }
        for (const CoarseVector& image : coarseImages) {
            bytes += image.indices.size() * sizeof(int) + image.values.size() * sizeof(double);
        }
        solveCalls = 0;
        matvecs = 0;
        localSeconds = 0.0;
        coarseSeconds = 0.0;
    }

    void solveLocal(LocalSolver& solver,
                    const std::vector<double>& rhs,
                    std::vector<double>& solution)
    {
        const auto start = std::chrono::steady_clock::now();
        solver.solve(rhs, solution);
        localSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        ++solveCalls;
    }

    void solveCoarse(const std::vector<double>& rhs, std::vector<double>& solution)
    {
        const auto start = std::chrono::steady_clock::now();
        coarseSolver.solve(rhs, solution);
        coarseSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }

    void applyRaw(const std::vector<double>& interfaceVector,
                  std::vector<double>& result,
                  bool countMatvec)
    {
        result = interfaceMatrix.multiply(interfaceVector);
        for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
            const DomainBlocks& domain = partition.domains[slot];
            std::vector<double> interiorProduct(domain.interiorGlobalDofs.size(), 0.0);
            for (std::size_t row = 0; row < domain.interiorInterfaceRows.size(); ++row) {
                for (const auto& entry : domain.interiorInterfaceRows[row]) {
                    interiorProduct[row] += entry.second * interfaceVector[static_cast<std::size_t>(entry.first)];
                }
            }
            std::vector<double> eliminated;
            solveLocal(interiorSolvers[slot], interiorProduct, eliminated);
            for (std::size_t localGamma = 0; localGamma < domain.interfaceGlobalDofs.size(); ++localGamma) {
                double correction = 0.0;
                for (const auto& entry : domain.interfaceInteriorRows[localGamma]) {
                    correction += entry.second * eliminated[static_cast<std::size_t>(entry.first)];
                }
                const int gamma = partition.globalToInterface[
                    static_cast<std::size_t>(domain.interfaceGlobalDofs[localGamma])];
                result[static_cast<std::size_t>(gamma)] -= correction;
            }
        }
        if (countMatvec) {
            ++matvecs;
        }
    }

    void buildCoarseCorrection(const Mesh& mesh,
                               bool coarseLinearXY,
                               bool coarseLinearZ,
                               bool coarseGlobalQuadraticZ,
                               bool coarseInterfacePatches,
                               bool coarseInterfaceLinearXY)
    {
        const bool addLinearXY = coarseInterfacePatches
            ? coarseInterfaceLinearXY
            : coarseLinearXY;
        const bool addLinearZ = !coarseInterfacePatches && coarseLinearZ;
        const bool addGlobalQuadraticZ = !coarseInterfacePatches
            && coarseGlobalQuadraticZ;
        std::vector<double> globalQuadraticZ(partition.interfaceGlobalDofs.size(), 0.0);

        if (coarseInterfacePatches) {
            using PatchKey = std::tuple<int, int, int, int>;
            std::map<PatchKey, std::vector<int>> patchDofs;
            for (const InterfaceFace& face : mesh.interfaceFaces) {
                if (face.leftTet < 0 || face.rightTet < 0) {
                    continue;
                }
                const Tet& left = mesh.tets[static_cast<std::size_t>(face.leftTet)];
                const Tet& right = mesh.tets[static_cast<std::size_t>(face.rightTet)];
                // A physical interface contributes one patch on each
                // subdomain boundary.  Keeping the two trace sides separate
                // lets the coarse space represent both mean and jump errors.
                const PatchKey leftKey{left.subdomain,
                                       right.subdomain,
                                       face.leftBoundaryEntity,
                                       face.rightBoundaryEntity};
                const PatchKey rightKey{right.subdomain,
                                        left.subdomain,
                                        face.rightBoundaryEntity,
                                        face.leftBoundaryEntity};
                // SIPG normal-flux couplings make every P2 DOF in the
                // boundary-tet stencil part of the algebraic Schur interface,
                // not only the six DOFs with nonzero geometric face trace.
                appendInterfaceStencilDofs(left, patchDofs[leftKey]);
                appendInterfaceStencilDofs(right, patchDofs[rightKey]);
            }

            // Make physical patches a disjoint partition of interface DOFs.
            // Shared edge/corner DOFs are assigned to the first deterministic
            // (domain pair, entity pair) key, avoiding dependent patch modes.
            std::vector<int> patchOwner(partition.interfaceGlobalDofs.size(), -1);
            for (auto& patch : patchDofs) {
                std::vector<int>& dofs = patch.second;
                std::sort(dofs.begin(), dofs.end());
                dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
                std::vector<int> ownedDofs;
                ownedDofs.reserve(dofs.size());
                for (int globalDof : dofs) {
                    const int gamma = partition.globalToInterface[static_cast<std::size_t>(globalDof)];
                    if (gamma >= 0 && patchOwner[static_cast<std::size_t>(gamma)] < 0) {
                        patchOwner[static_cast<std::size_t>(gamma)] = patchCount;
                        ownedDofs.push_back(globalDof);
                    }
                }
                if (ownedDofs.empty()) {
                    continue;
                }

                std::vector<int> interfaceIndices;
                interfaceIndices.reserve(ownedDofs.size());
                for (int globalDof : ownedDofs) {
                    interfaceIndices.push_back(
                        partition.globalToInterface[static_cast<std::size_t>(globalDof)]);
                }
                const std::size_t count = ownedDofs.size();
                CoarseVector constant;
                constant.indices = interfaceIndices;
                constant.values.assign(count, 1.0 / std::sqrt(static_cast<double>(count)));
                coarseBasis.push_back(std::move(constant));

                if (coarseInterfaceLinearXY) {
                    std::vector<std::vector<double>> localLinearModes;
                    for (int axis = 0; axis < 2; ++axis) {
                        double mean = 0.0;
                        for (int globalDof : ownedDofs) {
                            const Vec3& point = mesh.nodes[static_cast<std::size_t>(globalDof)].p;
                            mean += axis == 0 ? point.x : point.y;
                        }
                        mean /= static_cast<double>(count);
                        std::vector<double> values;
                        values.reserve(count);
                        for (int globalDof : ownedDofs) {
                            const Vec3& point = mesh.nodes[static_cast<std::size_t>(globalDof)].p;
                            values.push_back((axis == 0 ? point.x : point.y) - mean);
                        }
                        for (const std::vector<double>& previous : localLinearModes) {
                            double projection = 0.0;
                            for (std::size_t i = 0; i < count; ++i) {
                                projection += values[i] * previous[i];
                            }
                            for (std::size_t i = 0; i < count; ++i) {
                                values[i] -= projection * previous[i];
                            }
                        }
                        double normSquared = 0.0;
                        for (double value : values) {
                            normSquared += value * value;
                        }
                        if (normSquared > 1.0e-24) {
                            const double inverseNorm = 1.0 / std::sqrt(normSquared);
                            for (double& value : values) {
                                value *= inverseNorm;
                            }
                            CoarseVector linear;
                            linear.indices = interfaceIndices;
                            linear.values = values;
                            coarseBasis.push_back(std::move(linear));
                            localLinearModes.push_back(std::move(values));
                        }
                    }
                }
                ++patchCount;
            }
        }

        for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
            const DomainBlocks& domain = partition.domains[slot];
            if (domain.interfaceGlobalDofs.empty()) {
                continue;
            }
            if (coarseInterfacePatches) {
                continue;
            }

            const auto appendSubdomainModes = [&](const std::vector<int>& subdomainDofs) {
                const std::size_t count = subdomainDofs.size();
                if (count == 0) {
                    return;
                }
                std::vector<int> interfaceIndices;
                interfaceIndices.reserve(count);
                for (int globalDof : subdomainDofs) {
                    interfaceIndices.push_back(
                        partition.globalToInterface[static_cast<std::size_t>(globalDof)]);
                }

                CoarseVector constant;
                constant.indices = interfaceIndices;
                constant.values.assign(count, 1.0 / std::sqrt(static_cast<double>(count)));
                coarseBasis.push_back(std::move(constant));

                // Build a numerically stable local basis for span{1,x,y,z}.
                // Centering removes the constant component.  Modified
                // Gram-Schmidt removes correlations among coordinate modes on
                // thin or tilted interfaces without changing their span.
                std::vector<std::vector<double>> localLinearModes;
                const auto coordinate = [&](int globalDof, int axis) {
                    const Vec3& point = mesh.nodes[static_cast<std::size_t>(globalDof)].p;
                    return axis == 0 ? point.x : (axis == 1 ? point.y : point.z);
                };
                const auto appendLinearMode = [&](int axis) {
                    double mean = 0.0;
                    for (int globalDof : subdomainDofs) {
                        mean += coordinate(globalDof, axis);
                    }
                    mean /= static_cast<double>(count);

                    std::vector<double> values;
                    values.reserve(count);
                    for (int globalDof : subdomainDofs) {
                        values.push_back(coordinate(globalDof, axis) - mean);
                    }

                    for (const std::vector<double>& previous : localLinearModes) {
                        double projection = 0.0;
                        for (std::size_t i = 0; i < count; ++i) {
                            projection += values[i] * previous[i];
                        }
                        for (std::size_t i = 0; i < count; ++i) {
                            values[i] -= projection * previous[i];
                        }
                    }

                    double normSquared = 0.0;
                    for (double value : values) {
                        normSquared += value * value;
                    }
                    if (normSquared > 1.0e-24) {
                        const double inverseNorm = 1.0 / std::sqrt(normSquared);
                        for (double& value : values) {
                            value *= inverseNorm;
                        }
                        CoarseVector linear;
                        linear.indices = interfaceIndices;
                        linear.values = values;
                        coarseBasis.push_back(std::move(linear));
                        localLinearModes.push_back(std::move(values));
                        return true;
                    }
                    return false;
                };

                if (addLinearXY) {
                    appendLinearMode(0);
                    appendLinearMode(1);
                }
                const bool zModeAdded = addLinearZ && appendLinearMode(2);

                if (addGlobalQuadraticZ && zModeAdded) {
                    double zMean = 0.0;
                    for (int globalDof : subdomainDofs) {
                        zMean += coordinate(globalDof, 2);
                    }
                    zMean /= static_cast<double>(count);
                    std::vector<double> quadraticValues;
                    quadraticValues.reserve(count);
                    double quadraticMean = 0.0;
                    for (int globalDof : subdomainDofs) {
                        const double delta = coordinate(globalDof, 2) - zMean;
                        quadraticValues.push_back(delta * delta);
                        quadraticMean += delta * delta;
                    }
                    quadraticMean /= static_cast<double>(count);
                    for (double& value : quadraticValues) {
                        value -= quadraticMean;
                    }
                    for (const std::vector<double>& linearMode : localLinearModes) {
                        double projection = 0.0;
                        for (std::size_t i = 0; i < count; ++i) {
                            projection += quadraticValues[i] * linearMode[i];
                        }
                        for (std::size_t i = 0; i < count; ++i) {
                            quadraticValues[i] -= projection * linearMode[i];
                        }
                    }

                    double quadraticNormSquared = 0.0;
                    for (double value : quadraticValues) {
                        quadraticNormSquared += value * value;
                    }
                    if (quadraticNormSquared > 1.0e-24) {
                        for (std::size_t i = 0; i < count; ++i) {
                            globalQuadraticZ[static_cast<std::size_t>(interfaceIndices[i])] +=
                                quadraticValues[i];
                        }
                    }
                }
            };

            appendSubdomainModes(domain.interfaceGlobalDofs);
        }

        double globalQuadraticNormSquared = 0.0;
        for (double value : globalQuadraticZ) {
            globalQuadraticNormSquared += value * value;
        }
        if (addGlobalQuadraticZ && globalQuadraticNormSquared > 1.0e-24) {
            const double inverseNorm = 1.0 / std::sqrt(globalQuadraticNormSquared);
            CoarseVector quadratic;
            for (std::size_t i = 0; i < globalQuadraticZ.size(); ++i) {
                if (globalQuadraticZ[i] != 0.0) {
                    quadratic.indices.push_back(static_cast<int>(i));
                    quadratic.values.push_back(globalQuadraticZ[i] * inverseNorm);
                }
            }
            coarseBasis.push_back(std::move(quadratic));
        }

        const int coarseDim = static_cast<int>(coarseBasis.size());
        if (coarseDim == 0) {
            return;
        }

        std::vector<Entry> coarseEntries;
        std::vector<std::vector<std::pair<int, double>>> basisByDof(
            partition.interfaceGlobalDofs.size());
        for (int coarse = 0; coarse < coarseDim; ++coarse) {
            const CoarseVector& basis = coarseBasis[static_cast<std::size_t>(coarse)];
            for (std::size_t i = 0; i < basis.indices.size(); ++i) {
                basisByDof[static_cast<std::size_t>(basis.indices[i])].push_back(
                    {coarse, basis.values[i]});
            }
        }
        std::vector<double> coarseColumn(static_cast<std::size_t>(coarseDim), 0.0);
        std::vector<int> touchedRows;
        touchedRows.reserve(static_cast<std::size_t>(coarseDim));
        std::vector<unsigned char> rowTouched(static_cast<std::size_t>(coarseDim), 0);
        coarseImages.reserve(static_cast<std::size_t>(coarseDim));
        for (int col = 0; col < coarseDim; ++col) {
            std::vector<double> basis(partition.interfaceGlobalDofs.size(), 0.0);
            const CoarseVector& sparseBasis = coarseBasis[static_cast<std::size_t>(col)];
            for (std::size_t i = 0; i < sparseBasis.indices.size(); ++i) {
                basis[static_cast<std::size_t>(sparseBasis.indices[i])] = sparseBasis.values[i];
            }
            std::vector<double> image;
            applyRaw(basis, image, false);
            CoarseVector sparseImage;
            for (std::size_t i = 0; i < image.size(); ++i) {
                if (image[i] != 0.0) {
                    sparseImage.indices.push_back(static_cast<int>(i));
                    sparseImage.values.push_back(image[i]);
                }
            }
            for (std::size_t i = 0; i < sparseImage.indices.size(); ++i) {
                const int dof = sparseImage.indices[i];
                const double imageValue = sparseImage.values[i];
                for (const auto& incidence : basisByDof[static_cast<std::size_t>(dof)]) {
                    const int row = incidence.first;
                    if (rowTouched[static_cast<std::size_t>(row)] == 0) {
                        rowTouched[static_cast<std::size_t>(row)] = 1;
                        touchedRows.push_back(row);
                    }
                    coarseColumn[static_cast<std::size_t>(row)] +=
                        incidence.second * imageValue;
                }
            }
            std::sort(touchedRows.begin(), touchedRows.end());
            for (int row : touchedRows) {
                const double value = coarseColumn[static_cast<std::size_t>(row)];
                if (value != 0.0) {
                    coarseEntries.push_back({row, col, value});
                }
                coarseColumn[static_cast<std::size_t>(row)] = 0.0;
                rowTouched[static_cast<std::size_t>(row)] = 0;
            }
            touchedRows.clear();
            coarseImages.push_back(std::move(sparseImage));
        }
        coarseSolver = LocalSolver(coarseDim, coarseEntries);
    }

    void applyCoarse(const std::vector<double>& residual,
                     std::vector<double>& correction,
                     std::vector<double>* image)
    {
        correction.assign(residual.size(), 0.0);
        if (image != nullptr) {
            image->assign(residual.size(), 0.0);
        }
        if (coarseBasis.empty()) {
            return;
        }
        std::vector<double> coarseRhs(coarseBasis.size(), 0.0);
        for (std::size_t coarse = 0; coarse < coarseBasis.size(); ++coarse) {
            for (std::size_t i = 0; i < coarseBasis[coarse].indices.size(); ++i) {
                coarseRhs[coarse] += coarseBasis[coarse].values[i]
                    * residual[static_cast<std::size_t>(coarseBasis[coarse].indices[i])];
            }
        }
        std::vector<double> coefficients;
        solveCoarse(coarseRhs, coefficients);
        for (std::size_t coarse = 0; coarse < coarseBasis.size(); ++coarse) {
            for (std::size_t i = 0; i < coarseBasis[coarse].indices.size(); ++i) {
                correction[static_cast<std::size_t>(coarseBasis[coarse].indices[i])] +=
                    coarseBasis[coarse].values[i] * coefficients[coarse];
            }
            if (image != nullptr) {
                const CoarseVector& coarseImage = coarseImages[coarse];
                for (std::size_t i = 0; i < coarseImage.indices.size(); ++i) {
                    (*image)[static_cast<std::size_t>(coarseImage.indices[i])] +=
                        coarseImage.values[i] * coefficients[coarse];
                }
            }
        }
    }

    void applyLeftCoarseProjection(const std::vector<double>& vector,
                                   std::vector<double>& correction)
    {
        correction.assign(vector.size(), 0.0);
        if (coarseBasis.empty()) {
            return;
        }

        // coarseImages[j] = S * Z_j.  The assembled Schur complement is
        // symmetric, so (S * Z_j)^T * vector = Z_j^T * S * vector.  Reusing
        // these setup-time images applies Q*S without another Schur matvec.
        std::vector<double> coarseRhs(coarseBasis.size(), 0.0);
        for (std::size_t coarse = 0; coarse < coarseBasis.size(); ++coarse) {
            const CoarseVector& image = coarseImages[coarse];
            for (std::size_t i = 0; i < image.indices.size(); ++i) {
                coarseRhs[coarse] += image.values[i]
                    * vector[static_cast<std::size_t>(image.indices[i])];
            }
        }

        std::vector<double> coefficients;
        solveCoarse(coarseRhs, coefficients);
        for (std::size_t coarse = 0; coarse < coarseBasis.size(); ++coarse) {
            for (std::size_t i = 0; i < coarseBasis[coarse].indices.size(); ++i) {
                correction[static_cast<std::size_t>(coarseBasis[coarse].indices[i])] +=
                    coarseBasis[coarse].values[i] * coefficients[coarse];
            }
        }
    }
};

SchurOperator::SchurOperator(const Mesh& mesh,
                             const SparseMatrix& system,
                             bool coarseLinearXY,
                             bool coarseLinearZ,
                             bool coarseGlobalQuadraticZ,
                             bool coarseInterfacePatches,
                             bool coarseInterfaceLinearXY)
    : impl_(std::make_unique<Impl>(mesh,
                                   system,
                                   coarseLinearXY,
                                   coarseLinearZ,
                                   coarseGlobalQuadraticZ,
                                   coarseInterfacePatches,
                                   coarseInterfaceLinearXY)) {}

SchurOperator::~SchurOperator() = default;
SchurOperator::SchurOperator(SchurOperator&&) noexcept = default;
SchurOperator& SchurOperator::operator=(SchurOperator&&) noexcept = default;

int SchurOperator::domains() const { return static_cast<int>(impl_->partition.domains.size()); }
int SchurOperator::totalDofs() const { return impl_->partition.totalDofs; }
int SchurOperator::interfaceDofs() const { return static_cast<int>(impl_->partition.interfaceGlobalDofs.size()); }
int SchurOperator::interiorDofs() const { return totalDofs() - interfaceDofs(); }
int SchurOperator::coarseDimension() const { return static_cast<int>(impl_->coarseBasis.size()); }
int SchurOperator::interfacePatchCount() const { return impl_->patchCount; }
double SchurOperator::localFactorizationSeconds() const { return impl_->factorSeconds; }
double SchurOperator::localSymbolicAnalysisSeconds() const { return impl_->symbolicSeconds; }
double SchurOperator::localNumericalFactorizationSeconds() const { return impl_->numericalSeconds; }
double SchurOperator::localSolveSeconds() const { return impl_->localSeconds; }
double SchurOperator::coarseSolveSeconds() const { return impl_->coarseSeconds; }
std::size_t SchurOperator::memoryBytes() const { return impl_->bytes; }
int SchurOperator::localSolveCalls() const { return impl_->solveCalls; }
int SchurOperator::localSymbolicAnalysisCalls() const { return impl_->symbolicCalls; }
int SchurOperator::localNumericalFactorizationCalls() const { return impl_->numericalCalls; }
int SchurOperator::matvecCalls() const { return impl_->matvecs; }

std::vector<double> SchurOperator::condensedRhs(const std::vector<double>& globalRhs)
{
    if (static_cast<int>(globalRhs.size()) != totalDofs()) {
        throw std::runtime_error("[Schur] Global right-hand side has the wrong size.");
    }
    std::vector<double> result(static_cast<std::size_t>(interfaceDofs()), 0.0);
    for (int gamma = 0; gamma < interfaceDofs(); ++gamma) {
        result[static_cast<std::size_t>(gamma)] =
            globalRhs[static_cast<std::size_t>(impl_->partition.interfaceGlobalDofs[static_cast<std::size_t>(gamma)])];
    }
    for (std::size_t slot = 0; slot < impl_->partition.domains.size(); ++slot) {
        const DomainBlocks& domain = impl_->partition.domains[slot];
        std::vector<double> localRhs(domain.interiorGlobalDofs.size(), 0.0);
        for (std::size_t i = 0; i < domain.interiorGlobalDofs.size(); ++i) {
            localRhs[i] = globalRhs[static_cast<std::size_t>(domain.interiorGlobalDofs[i])];
        }
        std::vector<double> eliminated;
        impl_->solveLocal(impl_->interiorSolvers[slot], localRhs, eliminated);
        for (std::size_t localGamma = 0; localGamma < domain.interfaceGlobalDofs.size(); ++localGamma) {
            double correction = 0.0;
            for (const auto& entry : domain.interfaceInteriorRows[localGamma]) {
                correction += entry.second * eliminated[static_cast<std::size_t>(entry.first)];
            }
            const int gamma = impl_->partition.globalToInterface[
                static_cast<std::size_t>(domain.interfaceGlobalDofs[localGamma])];
            result[static_cast<std::size_t>(gamma)] -= correction;
        }
    }
    return result;
}

void SchurOperator::apply(const std::vector<double>& interfaceVector, std::vector<double>& result)
{
    if (static_cast<int>(interfaceVector.size()) != interfaceDofs()) {
        throw std::runtime_error("[Schur] Interface vector has the wrong size.");
    }
    impl_->applyRaw(interfaceVector, result, true);
}

void SchurOperator::applyBlockPreconditioner(const std::vector<double>& residual,
                                             std::vector<double>& result)
{
    if (static_cast<int>(residual.size()) != interfaceDofs()) {
        throw std::runtime_error("[Schur] Preconditioner vector has the wrong size.");
    }
    std::vector<double> coarseCorrection;
    std::vector<double> coarseImage;
    impl_->applyCoarse(residual, coarseCorrection, &coarseImage);
    std::vector<double> projectedResidual(residual.size(), 0.0);
    for (std::size_t i = 0; i < residual.size(); ++i) {
        projectedResidual[i] = residual[i] - coarseImage[i];
    }

    std::vector<double> localCorrection(residual.size(), 0.0);
    for (std::size_t slot = 0; slot < impl_->partition.domains.size(); ++slot) {
        const DomainBlocks& domain = impl_->partition.domains[slot];
        const std::size_t interiorCount = domain.interiorGlobalDofs.size();
        std::vector<double> localRhs(interiorCount + domain.interfaceGlobalDofs.size(), 0.0);
        for (std::size_t localGamma = 0; localGamma < domain.interfaceGlobalDofs.size(); ++localGamma) {
            const int gamma = impl_->partition.globalToInterface[
                static_cast<std::size_t>(domain.interfaceGlobalDofs[localGamma])];
            localRhs[interiorCount + localGamma] = projectedResidual[static_cast<std::size_t>(gamma)];
        }
        std::vector<double> localSolution;
        impl_->solveLocal(impl_->blockSolvers[slot], localRhs, localSolution);
        for (std::size_t localGamma = 0; localGamma < domain.interfaceGlobalDofs.size(); ++localGamma) {
            const int gamma = impl_->partition.globalToInterface[
                static_cast<std::size_t>(domain.interfaceGlobalDofs[localGamma])];
            localCorrection[static_cast<std::size_t>(gamma)] = localSolution[interiorCount + localGamma];
        }
    }

    std::vector<double> leftCoarseProjection;
    impl_->applyLeftCoarseProjection(localCorrection, leftCoarseProjection);
    result.resize(residual.size());
    for (std::size_t i = 0; i < residual.size(); ++i) {
        result[i] = coarseCorrection[i] + localCorrection[i] - leftCoarseProjection[i];
    }
}

std::vector<double> SchurOperator::recover(const std::vector<double>& globalRhs,
                                           const std::vector<double>& interfaceSolution)
{
    if (static_cast<int>(globalRhs.size()) != totalDofs()
        || static_cast<int>(interfaceSolution.size()) != interfaceDofs()) {
        throw std::runtime_error("[Schur] Recovery vector has the wrong size.");
    }
    std::vector<double> solution(static_cast<std::size_t>(totalDofs()), 0.0);
    for (int gamma = 0; gamma < interfaceDofs(); ++gamma) {
        solution[static_cast<std::size_t>(impl_->partition.interfaceGlobalDofs[static_cast<std::size_t>(gamma)])] =
            interfaceSolution[static_cast<std::size_t>(gamma)];
    }
    for (std::size_t slot = 0; slot < impl_->partition.domains.size(); ++slot) {
        const DomainBlocks& domain = impl_->partition.domains[slot];
        std::vector<double> localRhs(domain.interiorGlobalDofs.size(), 0.0);
        for (std::size_t row = 0; row < domain.interiorGlobalDofs.size(); ++row) {
            localRhs[row] = globalRhs[static_cast<std::size_t>(domain.interiorGlobalDofs[row])];
            for (const auto& entry : domain.interiorInterfaceRows[row]) {
                localRhs[row] -= entry.second * interfaceSolution[static_cast<std::size_t>(entry.first)];
            }
        }
        std::vector<double> localSolution;
        impl_->solveLocal(impl_->interiorSolvers[slot], localRhs, localSolution);
        for (std::size_t i = 0; i < domain.interiorGlobalDofs.size(); ++i) {
            solution[static_cast<std::size_t>(domain.interiorGlobalDofs[i])] = localSolution[i];
        }
    }
    return solution;
}

} // namespace ddm_schur
