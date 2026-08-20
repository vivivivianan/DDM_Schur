#pragma once

// Monolithic-equivalent Schwarz / block Gauss-Seidel solver built from the
// already assembled SIPG matrix. This keeps the original SIPG discretization
// intact and only changes the linear solve path.

struct SchwarzStepStats {
    int timeStep = 0;
    double time = 0.0;
    int iterations = 0;
    double relUpdate = std::numeric_limits<double>::quiet_NaN();
    double relResidual = std::numeric_limits<double>::quiet_NaN();
    double globalRelResidual = std::numeric_limits<double>::quiet_NaN();
    double freeRelResidual = std::numeric_limits<double>::quiet_NaN();
    double freeAbsResidual = std::numeric_limits<double>::quiet_NaN();
    double interfaceJump = std::numeric_limits<double>::quiet_NaN();
    double interfaceJumpL2 = std::numeric_limits<double>::quiet_NaN();
    double fluxBalanceL2 = std::numeric_limits<double>::quiet_NaN();
    double assemblyTime = 0.0;
    double solveTime = 0.0;
    double couplingTime = 0.0;
    double totalTime = 0.0;
    bool converged = false;
};

struct SchwarzSubdomainSolveStats {
    int timeStep = 0;
    double time = 0.0;
    int subdomain = -1;
    int dofs = 0;
    double solveTime = 0.0;
};

struct SchwarzInterfaceFluxRow {
    int timeStep = 0;
    int schwarzIter = 0;
    int interfaceId = -1;
    int leftSubdomain = -1;
    int rightSubdomain = -1;
    double meanTLeft = std::numeric_limits<double>::quiet_NaN();
    double meanTRight = std::numeric_limits<double>::quiet_NaN();
    double maxJumpT = std::numeric_limits<double>::quiet_NaN();
    double l2JumpT = std::numeric_limits<double>::quiet_NaN();
    double meanFluxLeft = std::numeric_limits<double>::quiet_NaN();
    double meanFluxRight = std::numeric_limits<double>::quiet_NaN();
    double fluxBalanceL2 = std::numeric_limits<double>::quiet_NaN();
    double fluxBalanceMax = std::numeric_limits<double>::quiet_NaN();
};

struct SchwarzTracePointData {
    double tLeft = 0.0;
    double tRight = 0.0;
    double qLeft = 0.0;
    double qRight = 0.0;
};

enum class SchwarzInterfaceCondition {
    None,
    Dirichlet,
    Neumann,
    Robin
};

struct SchwarzInterfaceDirichletConstraint {
    int globalDof = -1;
    std::array<int, 10> sourceDofs{};
    std::array<double, 10> sourceShape{};
};

static SchwarzInterfaceCondition schwarzConditionForSide(const std::string& standaloneMode,
                                                         const std::string& orientation,
                                                         bool sideIsLeft)
{
    if (standaloneMode == "dirichlet_dirichlet") {
        return SchwarzInterfaceCondition::Dirichlet;
    }
    const bool forward = orientation != "reverse";
    const bool dirichletSide = forward ? sideIsLeft : !sideIsLeft;
    if (standaloneMode == "dirichlet_neumann") {
        return dirichletSide ? SchwarzInterfaceCondition::Dirichlet
                             : SchwarzInterfaceCondition::Neumann;
    }
    if (standaloneMode == "dirichlet_robin") {
        return dirichletSide ? SchwarzInterfaceCondition::Dirichlet
                             : SchwarzInterfaceCondition::Robin;
    }
    if (standaloneMode == "robin") {
        return SchwarzInterfaceCondition::Robin;
    }
    return SchwarzInterfaceCondition::None;
}

static bool schwarzModeHasCondition(const std::string& standaloneMode,
                                    const std::string& orientation,
                                    SchwarzInterfaceCondition condition)
{
    return schwarzConditionForSide(standaloneMode, orientation, true) == condition
        || schwarzConditionForSide(standaloneMode, orientation, false) == condition;
}

static std::array<int, 6> p2FaceDofIndices(const std::array<int, 3>& faceLocal)
{
    const auto edgeDof = [](int a, int b) {
        if (a > b) {
            std::swap(a, b);
        }
        if (a == 0 && b == 1) { return 4; }
        if (a == 0 && b == 2) { return 5; }
        if (a == 0 && b == 3) { return 6; }
        if (a == 1 && b == 2) { return 7; }
        if (a == 1 && b == 3) { return 8; }
        if (a == 2 && b == 3) { return 9; }
        throw std::runtime_error("Invalid P2 tetrahedral face edge.");
    };
    return {{
        faceLocal[0],
        faceLocal[1],
        faceLocal[2],
        edgeDof(faceLocal[0], faceLocal[1]),
        edgeDof(faceLocal[0], faceLocal[2]),
        edgeDof(faceLocal[1], faceLocal[2])
    }};
}

static Vec3 schwarzInterfaceNormal(const Mesh& mesh,
                                   const InterfaceFace& face,
                                   const Tet& left,
                                   const Tet& right)
{
    return norm(face.leftNormal) > 1.0e-30
        ? face.leftNormal
        : normalized(subdomainCenter(mesh, right.subdomain) - subdomainCenter(mesh, left.subdomain));
}

static Vec3 schwarzRightInterfaceNormal(const Vec3& leftNormal,
                                        const InterfaceFace& face)
{
    return norm(face.rightNormal) > 1.0e-30 ? face.rightNormal : (-1.0 * leftNormal);
}

static double schwarzRobinAlphaForFace(const Mesh& mesh,
                                       const CaseConfig& config,
                                       const InterfaceFace& face,
                                       const Tet& left,
                                       const Tet& right,
                                       const ElementGeometry& gLeft,
                                       const ElementGeometry& gRight,
                                       const Vec3& normal)
{
    const Material& leftMaterial = materialForTet(config, left);
    const Material& rightMaterial = materialForTet(config, right);
    const std::array<Vec3, 3> leftPhysicalFace{{
        mesh.nodes[static_cast<size_t>(left.v[static_cast<size_t>(face.leftLocal[0])])].p,
        mesh.nodes[static_cast<size_t>(left.v[static_cast<size_t>(face.leftLocal[1])])].p,
        mesh.nodes[static_cast<size_t>(left.v[static_cast<size_t>(face.leftLocal[2])])].p
    }};
    const std::array<Vec3, 3> rightPhysicalFace{{
        mesh.nodes[static_cast<size_t>(right.v[static_cast<size_t>(face.rightLocal[0])])].p,
        mesh.nodes[static_cast<size_t>(right.v[static_cast<size_t>(face.rightLocal[1])])].p,
        mesh.nodes[static_cast<size_t>(right.v[static_cast<size_t>(face.rightLocal[2])])].p
    }};
    const double leftArea = triangleArea(leftPhysicalFace);
    const double rightArea = triangleArea(rightPhysicalFace);
    const double leftVolume = gLeft.detJ / 6.0;
    const double rightVolume = gRight.detJ / 6.0;
    const double hFace = (leftVolume + rightVolume) / std::max(1.0e-30, leftArea + rightArea);
    const double kLeft = std::max(1.0e-30, normalConductivity(leftMaterial, normal));
    const double kRight = std::max(1.0e-30, normalConductivity(rightMaterial, normal));
    const double kHarm = 2.0 * kLeft * kRight / std::max(1.0e-30, kLeft + kRight);
    return config.schwarz.robinAlphaFactor * kHarm / std::max(1.0e-30, hFace);
}

static SchwarzTracePointData schwarzTracePointData(const CaseConfig& config,
                                                   const Tet& left,
                                                   const Tet& right,
                                                   const ElementGeometry& gLeft,
                                                   const ElementGeometry& gRight,
                                                   const std::array<double, 4>& lLeft,
                                                   const std::array<double, 4>& lRight,
                                                   const Vec3& normal,
                                                   const Vec3& rightNormal,
                                                   double robinAlpha,
                                                   const std::vector<double>& temperature)
{
    const Material& leftMaterial = materialForTet(config, left);
    const Material& rightMaterial = materialForTet(config, right);
    const auto nLeft = shapeP2(lLeft);
    const auto nRight = shapeP2(lRight);
    const auto gradLeft = gradShapeP2(lLeft, gLeft);
    const auto gradRight = gradShapeP2(lRight, gRight);

    double tLeft = 0.0;
    double tRight = 0.0;
    Vec3 gradTLeft{};
    Vec3 gradTRight{};
    for (int i = 0; i < 10; ++i) {
        const double valueLeft = temperature[static_cast<size_t>(left.dof[static_cast<size_t>(i)])];
        const double valueRight = temperature[static_cast<size_t>(right.dof[static_cast<size_t>(i)])];
        tLeft += nLeft[static_cast<size_t>(i)] * valueLeft;
        tRight += nRight[static_cast<size_t>(i)] * valueRight;
        gradTLeft = gradTLeft + valueLeft * gradLeft[static_cast<size_t>(i)];
        gradTRight = gradTRight + valueRight * gradRight[static_cast<size_t>(i)];
    }

    const double dLeft = leftMaterial.conductivityX * gradTLeft.x * normal.x
                       + leftMaterial.conductivityY * gradTLeft.y * normal.y
                       + leftMaterial.conductivityZ * gradTLeft.z * normal.z;
    const double dRightAlongLeftNormal = rightMaterial.conductivityX * gradTRight.x * normal.x
                                       + rightMaterial.conductivityY * gradTRight.y * normal.y
                                       + rightMaterial.conductivityZ * gradTRight.z * normal.z;
    const double dRight = rightMaterial.conductivityX * gradTRight.x * rightNormal.x
                        + rightMaterial.conductivityY * gradTRight.y * rightNormal.y
                        + rightMaterial.conductivityZ * gradTRight.z * rightNormal.z;

    double qLeft = -dLeft;
    double qRight = -dRight;
    if (config.schwarz.fluxEval == "sipg_numeric") {
        const double averageNormalDerivative = 0.5 * (dLeft + dRightAlongLeftNormal);
        const double jump = tLeft - tRight;
        qLeft = -averageNormalDerivative + robinAlpha * jump;
        qRight = averageNormalDerivative - robinAlpha * jump;
    }
    return {tLeft, tRight, qLeft, qRight};
}

static std::vector<SchwarzInterfaceFluxRow> collectSchwarzInterfaceFluxRows(
    const Mesh& mesh,
    const CaseConfig& config,
    const std::vector<double>& temperature,
    int timeStep,
    int schwarzIter)
{
    const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
    std::vector<SchwarzInterfaceFluxRow> rows;
    rows.reserve(mesh.interfaceFaces.size());
    for (size_t faceIndex = 0; faceIndex < mesh.interfaceFaces.size(); ++faceIndex) {
        const InterfaceFace& face = mesh.interfaceFaces[faceIndex];
        const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
        const ElementGeometry gLeft = elementGeometry(mesh, left);
        const ElementGeometry gRight = elementGeometry(mesh, right);
        const Vec3 normal = schwarzInterfaceNormal(mesh, face, left, right);
        const Vec3 rightNormal = schwarzRightInterfaceNormal(normal, face);
        const double robinAlpha = schwarzRobinAlphaForFace(mesh, config, face, left, right, gLeft, gRight, normal);

        double area = 0.0;
        double sumTLeft = 0.0;
        double sumTRight = 0.0;
        double sumFluxLeft = 0.0;
        double sumFluxRight = 0.0;
        double jumpSquared = 0.0;
        double balanceSquared = 0.0;
        double maxJump = 0.0;
        double maxBalance = 0.0;

        for (const auto& integrationTriangle : face.integrationTriangles) {
            const Vec3 p0 = integrationTriangle[0];
            const Vec3 p1 = integrationTriangle[1];
            const Vec3 p2 = integrationTriangle[2];
            const double jacFace = norm(cross(p1 - p0, p2 - p0));
            for (const TriangleQuadraturePoint& qp : quadrature) {
                const double weight = qp.weight * jacFace;
                const Vec3 q = (1.0 - qp.a - qp.b) * p0 + qp.a * p1 + qp.b * p2;
                const auto lLeft = lambdaOnTetFace(q, left, face.leftLocal, mesh);
                const auto lRight = lambdaOnTetFace(q, right, face.rightLocal, mesh);
                const SchwarzTracePointData point =
                    schwarzTracePointData(config, left, right, gLeft, gRight,
                                          lLeft, lRight, normal, rightNormal, robinAlpha, temperature);
                const double jump = point.tLeft - point.tRight;
                const double balance = point.qLeft + point.qRight;
                area += weight;
                sumTLeft += point.tLeft * weight;
                sumTRight += point.tRight * weight;
                sumFluxLeft += point.qLeft * weight;
                sumFluxRight += point.qRight * weight;
                jumpSquared += jump * jump * weight;
                balanceSquared += balance * balance * weight;
                maxJump = std::max(maxJump, std::abs(jump));
                maxBalance = std::max(maxBalance, std::abs(balance));
            }
        }

        const double invArea = 1.0 / std::max(1.0e-300, area);
        rows.push_back({timeStep,
                        schwarzIter,
                        static_cast<int>(faceIndex),
                        left.subdomain,
                        right.subdomain,
                        sumTLeft * invArea,
                        sumTRight * invArea,
                        maxJump,
                        std::sqrt(jumpSquared * invArea),
                        sumFluxLeft * invArea,
                        sumFluxRight * invArea,
                        std::sqrt(balanceSquared * invArea),
                        maxBalance});
    }
    return rows;
}

struct SchwarzStepResult {
    std::vector<double> temperature;
    SchwarzStepStats stats;
    std::vector<SchwarzSubdomainSolveStats> subdomainStats;
    std::vector<SchwarzInterfaceFluxRow> interfaceFluxRows;
};

class SchwarzBlockSolver {
public:
    SchwarzBlockSolver(const Mesh& mesh,
                       const SparseMatrix& a,
                       const SchwarzOptions& options,
                       const CaseConfig* config = nullptr,
                       const SparseMatrix* robinBaseMatrix = nullptr,
                       bool factorLocalSolvers = true)
        : mesh_(&mesh),
          a_(a),
          options_(options),
          config_(config),
          robinBaseMatrix_(robinBaseMatrix),
          factorLocalSolvers_(factorLocalSolvers)
    {
        build(mesh, a);
    }

    SchwarzStepResult solveStep(const Mesh& mesh,
                                const std::vector<double>& rhs,
                                const std::vector<double>& initialGuess,
                                int timeStep,
                                double time,
                                double assemblyTimeForStep)
    {
        if (static_cast<int>(rhs.size()) != a_.size()) {
            throw std::runtime_error("Schwarz RHS size does not match the system size.");
        }
        if (static_cast<int>(initialGuess.size()) != a_.size()) {
            throw std::runtime_error("Schwarz initial guess size does not match the system size.");
        }
        if (usesPhysicalTransmission() && config_ == nullptr) {
            throw std::runtime_error("Physical Schwarz standalone modes require CaseConfig for interface transmission data.");
        }

        SchwarzStepResult result;
        result.temperature = initialGuess;
        result.stats.timeStep = timeStep;
        result.stats.time = time;
        result.stats.assemblyTime = assemblyTimeForStep;
        result.subdomainStats.reserve(blockDofs_.size());

        const auto stepStart = std::chrono::steady_clock::now();
        const int maxIters = std::max(1, options_.maxIters);
        const double relaxation = options_.relaxation;
        if (!(relaxation > 0.0) || !std::isfinite(relaxation)) {
            throw std::runtime_error("Schwarz relaxation must be a positive finite value.");
        }

        double solveSeconds = 0.0;
        double couplingSeconds = 0.0;
        const double rhsNorm = std::max(1.0e-300, l2Norm(rhs));

        for (int iter = 1; iter <= maxIters; ++iter) {
            const std::vector<double> iterationStart = result.temperature;
            std::vector<double> additiveCandidate;
            if (options_.type == "additive") {
                additiveCandidate = result.temperature;
            }

            for (size_t blockIndex = 0; blockIndex < blockDofs_.size(); ++blockIndex) {
                const std::vector<int>& coreDofs = blockDofs_[blockIndex];
                if (coreDofs.empty() || localDofs_[blockIndex].empty()) {
                    continue;
                }

                const std::vector<double>& neighborTemperature =
                    options_.type == "additive" ? iterationStart : result.temperature;

                const auto couplingStart = std::chrono::steady_clock::now();
                std::vector<double> localRhs = usesPhysicalTransmission()
                    ? buildPhysicalTransmissionLocalRhs(blockIndex, rhs, neighborTemperature)
                    : buildLocalRhs(blockIndex, rhs, neighborTemperature);
                const auto couplingEnd = std::chrono::steady_clock::now();
                couplingSeconds += std::chrono::duration<double>(couplingEnd - couplingStart).count();

                std::vector<double> localSolution;
                const auto solveStart = std::chrono::steady_clock::now();
                {
                    ScopedMklSingleThread mklThreads;
                    if (usesPhysicalTransmission()) {
                        robinSolvers_[blockIndex].solve(localRhs, localSolution);
                    } else {
                        solvers_[blockIndex].solve(localRhs, localSolution);
                    }
                }
                const auto solveEnd = std::chrono::steady_clock::now();
                const double localSolveSeconds = std::chrono::duration<double>(solveEnd - solveStart).count();
                solveSeconds += localSolveSeconds;

                if (iter == 1) {
                    result.subdomainStats.push_back({
                        timeStep,
                        time,
                        static_cast<int>(blockIndex),
                        static_cast<int>(localDofs_[blockIndex].size()),
                        localSolveSeconds
                    });
                } else if (blockIndex < result.subdomainStats.size()) {
                    result.subdomainStats[blockIndex].solveTime += localSolveSeconds;
                }

                std::vector<double>& target = options_.type == "additive"
                    ? additiveCandidate
                    : result.temperature;
                for (int globalDof : coreDofs) {
                    const int local = localIndex_[blockIndex].at(globalDof);
                    const double oldValue = iterationStart[static_cast<size_t>(globalDof)];
                    target[static_cast<size_t>(globalDof)] =
                        (1.0 - relaxation) * oldValue + relaxation * localSolution[static_cast<size_t>(local)];
                }
            }

            if (options_.type == "additive") {
                result.temperature = std::move(additiveCandidate);
            }

            result.stats.iterations = iter;
            result.stats.relUpdate = relativeUpdateNorm(result.temperature, iterationStart);
            result.stats.globalRelResidual = relativeSystemResidual(result.temperature, rhs, rhsNorm);
            result.stats.freeRelResidual =
                relativeFreeSystemResidual(result.temperature, rhs, result.stats.freeAbsResidual);
            result.stats.relResidual = result.stats.globalRelResidual;
            if (options_.writeInterfaceFlux && config_ != nullptr) {
                result.interfaceFluxRows =
                    collectSchwarzInterfaceFluxRows(mesh, *config_, result.temperature, timeStep, iter);
                double maxJumpL2 = 0.0;
                double maxFluxBalanceL2 = 0.0;
                for (const SchwarzInterfaceFluxRow& row : result.interfaceFluxRows) {
                    maxJumpL2 = std::max(maxJumpL2, row.l2JumpT);
                    maxFluxBalanceL2 = std::max(maxFluxBalanceL2, row.fluxBalanceL2);
                }
                result.stats.interfaceJumpL2 = maxJumpL2;
                result.stats.fluxBalanceL2 = maxFluxBalanceL2;
            }
            const double convergenceResidual = std::isfinite(result.stats.freeRelResidual)
                ? result.stats.freeRelResidual
                : result.stats.relResidual;
            if (result.stats.relUpdate <= options_.tolRelUpdate
                && convergenceResidual <= options_.tolRelResidual) {
                result.stats.converged = true;
                break;
            }
        }

        result.stats.solveTime = solveSeconds;
        result.stats.couplingTime = couplingSeconds;
        result.stats.interfaceJump = options_.checkInterfaceJump
            ? interfaceAverageJump(mesh, result.temperature)
            : std::numeric_limits<double>::quiet_NaN();
        result.stats.totalTime = assemblyTimeForStep
            + std::chrono::duration<double>(std::chrono::steady_clock::now() - stepStart).count();
        return result;
    }

    void apply(const std::vector<double>& r, std::vector<double>& z)
    {
        if (static_cast<int>(r.size()) != a_.size()) {
            throw std::runtime_error("Schwarz preconditioner input size does not match the system size.");
        }
        ++applyCallCount_;
        std::vector<std::vector<double>> localRhs(localDofs_.size());
        std::vector<std::vector<double>> localSolutions(localDofs_.size());

        const auto haloUpdateStart = std::chrono::steady_clock::now();
        parallelForCoarse(blockDofs_.size(), [&](size_t blockIndex) {
            const std::vector<int>& dofs = localDofs_[blockIndex];
            if (dofs.empty()) {
                return;
            }
            std::vector<double> rhsBlock(dofs.size(), 0.0);
            for (size_t local = 0; local < dofs.size(); ++local) {
                rhsBlock[local] = r[static_cast<size_t>(dofs[local])];
            }
            localRhs[blockIndex] = std::move(rhsBlock);
        });
        const auto haloUpdateEnd = std::chrono::steady_clock::now();
        communicationOrHaloUpdateSeconds_ +=
            std::chrono::duration<double>(haloUpdateEnd - haloUpdateStart).count();

        const auto localSolveStart = std::chrono::steady_clock::now();
        parallelForCoarse(blockDofs_.size(), [&](size_t blockIndex) {
            if (localRhs[blockIndex].empty()) {
                return;
            }
            std::vector<double> localSolution;
            {
                ScopedMklSingleThread mklThreads;
                solvers_[blockIndex].solve(localRhs[blockIndex], localSolution);
            }
            localSolutions[blockIndex] = std::move(localSolution);
        });
        const auto localSolveEnd = std::chrono::steady_clock::now();
        localSolveApplySeconds_ +=
            std::chrono::duration<double>(localSolveEnd - localSolveStart).count();

        const auto restrictionStart = std::chrono::steady_clock::now();
        z.assign(r.size(), 0.0);
        if (options_.overlapMode == "halo" && options_.overlapLayers > 0) {
            std::vector<int> weights(r.size(), 0);
            for (size_t blockIndex = 0; blockIndex < localDofs_.size(); ++blockIndex) {
                const std::vector<int>& dofs = localDofs_[blockIndex];
                const std::vector<double>& localSolution = localSolutions[blockIndex];
                for (size_t local = 0; local < dofs.size() && local < localSolution.size(); ++local) {
                    const int globalDof = dofs[local];
                    z[static_cast<size_t>(globalDof)] += localSolution[local];
                    ++weights[static_cast<size_t>(globalDof)];
                }
            }
            for (size_t i = 0; i < z.size(); ++i) {
                if (weights[i] > 1) {
                    z[i] /= static_cast<double>(weights[i]);
                }
            }
            const auto restrictionEnd = std::chrono::steady_clock::now();
            restrictionSeconds_ +=
                std::chrono::duration<double>(restrictionEnd - restrictionStart).count();
            return;
        }

        for (size_t blockIndex = 0; blockIndex < blockDofs_.size(); ++blockIndex) {
            const std::vector<double>& localSolution = localSolutions[blockIndex];
            for (int globalDof : blockDofs_[blockIndex]) {
                const int local = localIndex_[blockIndex].at(globalDof);
                if (local >= 0 && local < static_cast<int>(localSolution.size())) {
                    z[static_cast<size_t>(globalDof)] = localSolution[static_cast<size_t>(local)];
                }
            }
        }
        const auto restrictionEnd = std::chrono::steady_clock::now();
        restrictionSeconds_ +=
            std::chrono::duration<double>(restrictionEnd - restrictionStart).count();
    }

    size_t memoryBytes() const
    {
        size_t bytes = globalToBlock_.size() * sizeof(int)
            + tetToBlock_.size() * sizeof(int)
            + globalToLocal_.size() * sizeof(int);
        for (const auto& dofs : blockDofs_) {
            bytes += dofs.size() * sizeof(int);
        }
        for (const auto& dofs : localDofs_) {
            bytes += dofs.size() * sizeof(int);
        }
        for (const auto& map : localIndex_) {
            bytes += map.size() * (sizeof(int) * 2 + sizeof(size_t));
        }
        for (const auto& entries : offBlockEntries_) {
            bytes += entries.size() * sizeof(MatrixEntry);
        }
        for (const auto& solver : solvers_) {
            bytes += solver.memoryBytes();
        }
        for (const auto& solver : robinSolvers_) {
            bytes += solver.memoryBytes();
        }
        return bytes;
    }

    size_t blockCount() const { return blockDofs_.size(); }
    size_t overlapLayers() const { return static_cast<size_t>(std::max(0, options_.overlapLayers)); }
    std::string overlapMode() const { return options_.overlapMode; }
    std::string partitionMode() const { return options_.partitionMode; }
    size_t coreDofs(size_t blockIndex) const { return blockIndex < blockDofs_.size() ? blockDofs_[blockIndex].size() : 0; }
    size_t localDofs(size_t blockIndex) const { return blockIndex < localDofs_.size() ? localDofs_[blockIndex].size() : 0; }
    size_t haloDofs(size_t blockIndex) const
    {
        return localDofs(blockIndex) >= coreDofs(blockIndex)
            ? localDofs(blockIndex) - coreDofs(blockIndex)
            : 0;
    }
    size_t ownedElements(size_t blockIndex) const { return blockIndex < ownedElementCounts_.size() ? ownedElementCounts_[blockIndex] : 0; }
    size_t localElements(size_t blockIndex) const { return blockIndex < localElementCounts_.size() ? localElementCounts_[blockIndex] : 0; }
    size_t haloElements(size_t blockIndex) const
    {
        return localElements(blockIndex) >= ownedElements(blockIndex)
            ? localElements(blockIndex) - ownedElements(blockIndex)
            : 0;
    }
    size_t diagonalBlockNonzeros() const { return diagonalBlockNonzeros_; }
    size_t couplingBlockNonzeros() const { return couplingBlockNonzeros_; }
    size_t globalNonzeros() const { return globalNonzeros_; }
    size_t partitionedNonzeros() const { return diagonalBlockNonzeros_ + couplingBlockNonzeros_; }
    size_t invalidMatrixEntries() const { return invalidMatrixEntries_; }
    double reconstructionMaxAbsError() const { return reconstructionMaxAbsError_; }
    double reconstructionL1Error() const { return reconstructionL1Error_; }
    double reconstructionFrobeniusError() const { return reconstructionFrobeniusError_; }
    double globalEntryL1Norm() const { return globalEntryL1Norm_; }
    double partitionedEntryL1Norm() const { return partitionedEntryL1Norm_; }
    double globalEntryFrobeniusNorm() const { return std::sqrt(globalEntryFrobeniusSquared_); }
    double partitionedEntryFrobeniusNorm() const { return std::sqrt(partitionedEntryFrobeniusSquared_); }
    int dofBlock(int globalDof) const
    {
        return globalDof >= 0 && globalDof < static_cast<int>(globalToBlock_.size())
            ? globalToBlock_[static_cast<size_t>(globalDof)]
            : -1;
    }
    int tetBlock(int tetId) const
    {
        return tetId >= 0 && tetId < static_cast<int>(tetToBlock_.size())
            ? tetToBlock_[static_cast<size_t>(tetId)]
            : -1;
    }
    void copyTimingTo(SolverStatistics& stats) const
    {
        stats.rasHaloBuildSeconds = haloBuildSeconds_;
        stats.rasLocalMatrixAssemblySeconds = localMatrixAssemblySeconds_;
        stats.rasLocalFactorizationSeconds = localFactorizationSeconds_;
        stats.rasLocalSolveApplySeconds = localSolveApplySeconds_;
        stats.rasRestrictionSeconds = restrictionSeconds_;
        stats.rasCommunicationOrHaloUpdateSeconds = communicationOrHaloUpdateSeconds_;
        stats.rasSetupCount = setupCount_;
        stats.rasFactorizationReuse = setupCount_ == 1 && applyCallCount_ > 0;
    }

private:
    const Mesh* mesh_ = nullptr;
    const SparseMatrix& a_;
    SchwarzOptions options_;
    const CaseConfig* config_ = nullptr;
    const SparseMatrix* robinBaseMatrix_ = nullptr;
    bool factorLocalSolvers_ = true;
    std::vector<std::vector<int>> blockDofs_;
    std::vector<std::vector<int>> localDofs_;
    std::vector<std::unordered_map<int, int>> localIndex_;
    std::vector<size_t> ownedElementCounts_;
    std::vector<size_t> localElementCounts_;
    std::vector<int> globalToBlock_;
    std::vector<int> tetToBlock_;
    std::vector<int> globalToLocal_;
    std::vector<std::vector<MatrixEntry>> offBlockEntries_;
    std::vector<std::vector<MatrixEntry>> localDirichletColumnEntries_;
    std::vector<std::vector<int>> interfaceDirichletLocalRows_;
    std::vector<std::vector<char>> interfaceDirichletMask_;
    std::vector<std::vector<SchwarzInterfaceDirichletConstraint>> interfaceDirichletConstraints_;
    std::vector<SubdomainDirectSolver> solvers_;
    std::vector<SubdomainDirectSolver> robinSolvers_;
    size_t globalNonzeros_ = 0;
    size_t diagonalBlockNonzeros_ = 0;
    size_t couplingBlockNonzeros_ = 0;
    size_t invalidMatrixEntries_ = 0;
    double globalEntryL1Norm_ = 0.0;
    double partitionedEntryL1Norm_ = 0.0;
    double globalEntryFrobeniusSquared_ = 0.0;
    double partitionedEntryFrobeniusSquared_ = 0.0;
    double reconstructionMaxAbsError_ = 0.0;
    double reconstructionL1Error_ = 0.0;
    double reconstructionFrobeniusError_ = 0.0;
    double haloBuildSeconds_ = 0.0;
    double localMatrixAssemblySeconds_ = 0.0;
    double localFactorizationSeconds_ = 0.0;
    double localSolveApplySeconds_ = 0.0;
    double restrictionSeconds_ = 0.0;
    double communicationOrHaloUpdateSeconds_ = 0.0;
    int setupCount_ = 0;
    int applyCallCount_ = 0;

    bool usesPhysicalTransmission() const
    {
        return isPhysicalSchwarzStandaloneMode(options_.standaloneMode);
    }

    static Vec3 tetCentroid(const Mesh& mesh, const Tet& tet)
    {
        Vec3 c{};
        for (int vertex : tet.v) {
            c = c + mesh.nodes[static_cast<size_t>(vertex)].p;
        }
        return c * 0.25;
    }

    static std::string materialPartitionKey(const Material& material)
    {
        std::ostringstream out;
        out << material.name << ":"
            << std::setprecision(8) << material.conductivityX << ":"
            << material.conductivityY << ":"
            << material.conductivityZ;
        return out.str();
    }

    int materialBlockForTet(const Tet& tet, std::map<std::string, int>& blockByMaterial) const
    {
        if (config_ == nullptr) {
            return tet.subdomain;
        }
        const Material& material = materialForTet(*config_, tet);
        const std::string key = materialPartitionKey(material);
        auto found = blockByMaterial.find(key);
        if (found != blockByMaterial.end()) {
            return found->second;
        }
        const int block = static_cast<int>(blockByMaterial.size());
        blockByMaterial[key] = block;
        return block;
    }

    bool tetHasHeatSource(const Tet& tet) const
    {
        if (config_ == nullptr) {
            return false;
        }
        for (const HeatSource& source : config_->heatSources) {
            if (tetMatchesHeatSource(tet, source)) {
                return true;
            }
        }
        return false;
    }

    int verticalColumnBlock(const Vec3& p, const Vec3& lo, const Vec3& hi) const
    {
        constexpr int binsX = 5;
        constexpr int binsY = 5;
        const double dx = std::max(1.0e-300, hi.x - lo.x);
        const double dy = std::max(1.0e-300, hi.y - lo.y);
        int ix = static_cast<int>(std::floor((p.x - lo.x) / dx * binsX));
        int iy = static_cast<int>(std::floor((p.y - lo.y) / dy * binsY));
        ix = std::max(0, std::min(binsX - 1, ix));
        iy = std::max(0, std::min(binsY - 1, iy));
        return iy * binsX + ix;
    }

    void finalizePartitionMaps(const Mesh& mesh,
                               std::vector<int> nodeToBlock,
                               std::vector<int> tetToBlock)
    {
        if (nodeToBlock.size() != mesh.nodes.size()) {
            throw std::runtime_error("Schwarz partition node map has invalid size.");
        }
        if (tetToBlock.size() != mesh.tets.size()) {
            throw std::runtime_error("Schwarz partition tet map has invalid size.");
        }

        int maxBlock = 0;
        for (int& block : nodeToBlock) {
            if (block < 0) {
                block = 0;
            }
            maxBlock = std::max(maxBlock, block);
        }
        for (size_t tetId = 0; tetId < mesh.tets.size(); ++tetId) {
            int& block = tetToBlock[tetId];
            if (block < 0) {
                block = nodeToBlock[static_cast<size_t>(mesh.tets[tetId].dof[0])];
            }
            maxBlock = std::max(maxBlock, block);
        }

        blockDofs_.assign(static_cast<size_t>(maxBlock + 1), {});
        globalToBlock_.assign(mesh.nodes.size(), -1);
        globalToLocal_.assign(mesh.nodes.size(), -1);
        tetToBlock_ = std::move(tetToBlock);

        for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
            const int block = nodeToBlock[static_cast<size_t>(i)];
            if (block < 0 || block > maxBlock) {
                throw std::runtime_error("Invalid Schwarz partition block id.");
            }
            globalToBlock_[static_cast<size_t>(i)] = block;
            globalToLocal_[static_cast<size_t>(i)] =
                static_cast<int>(blockDofs_[static_cast<size_t>(block)].size());
            blockDofs_[static_cast<size_t>(block)].push_back(i);
        }
    }

    void buildCurrentPartitionMaps(const Mesh& mesh)
    {
        std::vector<int> nodeToBlock(mesh.nodes.size(), 0);
        std::vector<int> tetToBlock(mesh.tets.size(), 0);
        for (size_t i = 0; i < mesh.nodes.size(); ++i) {
            nodeToBlock[i] = mesh.nodes[i].subdomain;
        }
        for (size_t tetId = 0; tetId < mesh.tets.size(); ++tetId) {
            tetToBlock[tetId] = mesh.tets[tetId].subdomain;
        }
        finalizePartitionMaps(mesh, std::move(nodeToBlock), std::move(tetToBlock));
    }

    void buildMaterialAlignedPartitionMaps(const Mesh& mesh)
    {
        std::map<std::string, int> blockByMaterial;
        std::vector<int> tetToBlock(mesh.tets.size(), 0);
        std::vector<std::unordered_map<int, int>> nodeVotes(mesh.nodes.size());
        for (size_t tetId = 0; tetId < mesh.tets.size(); ++tetId) {
            const int block = materialBlockForTet(mesh.tets[tetId], blockByMaterial);
            tetToBlock[tetId] = block;
            for (int dof : mesh.tets[tetId].dof) {
                ++nodeVotes[static_cast<size_t>(dof)][block];
            }
        }

        std::vector<int> nodeToBlock(mesh.nodes.size(), 0);
        for (size_t i = 0; i < nodeVotes.size(); ++i) {
            int bestBlock = 0;
            int bestVotes = -1;
            for (const auto& vote : nodeVotes[i]) {
                if (vote.second > bestVotes || (vote.second == bestVotes && vote.first < bestBlock)) {
                    bestBlock = vote.first;
                    bestVotes = vote.second;
                }
            }
            nodeToBlock[i] = bestBlock;
        }
        finalizePartitionMaps(mesh, std::move(nodeToBlock), std::move(tetToBlock));
    }

    void buildHotspotContainedPartitionMaps(const Mesh& mesh)
    {
        std::vector<int> nodeToBlock(mesh.nodes.size(), 0);
        std::vector<int> tetToBlock(mesh.tets.size(), 0);
        for (size_t i = 0; i < mesh.nodes.size(); ++i) {
            nodeToBlock[i] = mesh.nodes[i].subdomain + 1;
        }
        for (size_t tetId = 0; tetId < mesh.tets.size(); ++tetId) {
            const Tet& tet = mesh.tets[tetId];
            const bool hotspot = tetHasHeatSource(tet);
            tetToBlock[tetId] = hotspot ? 0 : tet.subdomain + 1;
            if (hotspot) {
                for (int dof : tet.dof) {
                    nodeToBlock[static_cast<size_t>(dof)] = 0;
                }
            }
        }
        finalizePartitionMaps(mesh, std::move(nodeToBlock), std::move(tetToBlock));
    }

    void buildVerticalHeatFlowAlignedPartitionMaps(const Mesh& mesh)
    {
        Vec3 lo{std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max()};
        Vec3 hi{-std::numeric_limits<double>::max(),
                -std::numeric_limits<double>::max(),
                -std::numeric_limits<double>::max()};
        for (const Node& node : mesh.nodes) {
            lo.x = std::min(lo.x, node.p.x);
            lo.y = std::min(lo.y, node.p.y);
            lo.z = std::min(lo.z, node.p.z);
            hi.x = std::max(hi.x, node.p.x);
            hi.y = std::max(hi.y, node.p.y);
            hi.z = std::max(hi.z, node.p.z);
        }

        std::vector<int> nodeToBlock(mesh.nodes.size(), 0);
        for (size_t i = 0; i < mesh.nodes.size(); ++i) {
            nodeToBlock[i] = verticalColumnBlock(mesh.nodes[i].p, lo, hi);
        }
        std::vector<int> tetToBlock(mesh.tets.size(), 0);
        for (size_t tetId = 0; tetId < mesh.tets.size(); ++tetId) {
            tetToBlock[tetId] = verticalColumnBlock(tetCentroid(mesh, mesh.tets[tetId]), lo, hi);
        }
        finalizePartitionMaps(mesh, std::move(nodeToBlock), std::move(tetToBlock));
    }

    void buildPartitionMaps(const Mesh& mesh)
    {
        if (options_.partitionMode == "current") {
            buildCurrentPartitionMaps(mesh);
        } else if (options_.partitionMode == "material_aligned") {
            buildMaterialAlignedPartitionMaps(mesh);
        } else if (options_.partitionMode == "hotspot_contained") {
            buildHotspotContainedPartitionMaps(mesh);
        } else if (options_.partitionMode == "vertical_heat_flow_aligned") {
            buildVerticalHeatFlowAlignedPartitionMaps(mesh);
        } else {
            throw std::runtime_error("Unknown Schwarz partition mode: " + options_.partitionMode);
        }
    }

    std::vector<std::vector<int>> buildTetAdjacency(const Mesh& mesh) const
    {
        std::vector<std::set<int>> adjacencySets(mesh.tets.size());
        std::unordered_map<int, std::vector<int>> tetsByDof;
        tetsByDof.reserve(mesh.nodes.size());
        for (int tetId = 0; tetId < static_cast<int>(mesh.tets.size()); ++tetId) {
            const Tet& tet = mesh.tets[static_cast<size_t>(tetId)];
            for (int dof : tet.dof) {
                tetsByDof[dof].push_back(tetId);
            }
        }
        for (const auto& item : tetsByDof) {
            const std::vector<int>& tets = item.second;
            for (int a : tets) {
                for (int b : tets) {
                    if (a != b) {
                        adjacencySets[static_cast<size_t>(a)].insert(b);
                    }
                }
            }
        }
        for (const InterfaceFace& face : mesh.interfaceFaces) {
            if (face.leftTet >= 0 && face.rightTet >= 0
                && face.leftTet < static_cast<int>(mesh.tets.size())
                && face.rightTet < static_cast<int>(mesh.tets.size())) {
                adjacencySets[static_cast<size_t>(face.leftTet)].insert(face.rightTet);
                adjacencySets[static_cast<size_t>(face.rightTet)].insert(face.leftTet);
            }
        }
        std::vector<std::vector<int>> adjacency(adjacencySets.size());
        for (size_t i = 0; i < adjacencySets.size(); ++i) {
            adjacency[i].assign(adjacencySets[i].begin(), adjacencySets[i].end());
        }
        return adjacency;
    }

    std::vector<std::vector<int>> buildOverlapDofsFromHaloElements(const Mesh& mesh, int requestedLayers)
    {
        const int layers = options_.overlapMode == "none" ? 0 : std::max(0, std::min(3, requestedLayers));

        std::vector<std::vector<int>> ownedTets(blockDofs_.size());
        for (int tetId = 0; tetId < static_cast<int>(mesh.tets.size()); ++tetId) {
            const int block = tetBlock(tetId);
            if (block >= 0 && block < static_cast<int>(ownedTets.size())) {
                ownedTets[static_cast<size_t>(block)].push_back(tetId);
            }
        }
        ownedElementCounts_.assign(blockDofs_.size(), 0);
        localElementCounts_.assign(blockDofs_.size(), 0);
        for (size_t blockIndex = 0; blockIndex < ownedTets.size(); ++blockIndex) {
            ownedElementCounts_[blockIndex] = ownedTets[blockIndex].size();
            localElementCounts_[blockIndex] = ownedTets[blockIndex].size();
        }
        if (layers == 0) {
            return blockDofs_;
        }

        const std::vector<std::vector<int>> adjacency = buildTetAdjacency(mesh);
        std::vector<std::vector<int>> result(blockDofs_.size());
        for (size_t blockIndex = 0; blockIndex < blockDofs_.size(); ++blockIndex) {
            std::vector<char> inTet(mesh.tets.size(), 0);
            std::vector<int> localTets = ownedTets[blockIndex];
            std::vector<int> frontier = ownedTets[blockIndex];
            for (int tet : localTets) {
                inTet[static_cast<size_t>(tet)] = 1;
            }
            for (int layer = 0; layer < layers; ++layer) {
                std::vector<int> nextFrontier;
                for (int tet : frontier) {
                    for (int neighbor : adjacency[static_cast<size_t>(tet)]) {
                        if (!inTet[static_cast<size_t>(neighbor)]) {
                            inTet[static_cast<size_t>(neighbor)] = 1;
                            localTets.push_back(neighbor);
                            nextFrontier.push_back(neighbor);
                        }
                    }
                }
                frontier = std::move(nextFrontier);
                if (frontier.empty()) {
                    break;
                }
            }

            std::vector<char> inDof(mesh.nodes.size(), 0);
            std::vector<int> dofs;
            dofs.reserve(blockDofs_[blockIndex].size());
            localElementCounts_[blockIndex] = localTets.size();
            for (int dof : blockDofs_[blockIndex]) {
                if (!inDof[static_cast<size_t>(dof)]) {
                    inDof[static_cast<size_t>(dof)] = 1;
                    dofs.push_back(dof);
                }
            }
            for (int tetId : localTets) {
                const Tet& tet = mesh.tets[static_cast<size_t>(tetId)];
                for (int dof : tet.dof) {
                    if (!inDof[static_cast<size_t>(dof)]) {
                        inDof[static_cast<size_t>(dof)] = 1;
                        dofs.push_back(dof);
                    }
                }
            }
            std::sort(dofs.begin(), dofs.end());
            dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
            result[blockIndex] = std::move(dofs);
        }
        return result;
    }

    void build(const Mesh& mesh, const SparseMatrix& a)
    {
        haloBuildSeconds_ = 0.0;
        localMatrixAssemblySeconds_ = 0.0;
        localFactorizationSeconds_ = 0.0;
        localSolveApplySeconds_ = 0.0;
        restrictionSeconds_ = 0.0;
        communicationOrHaloUpdateSeconds_ = 0.0;
        applyCallCount_ = 0;
        setupCount_ = 1;
        if (!a.csrReady) {
            throw std::runtime_error("SchwarzBlockSolver requires CSR-finalized system matrix.");
        }
        if (options_.type != "multiplicative" && options_.type != "additive") {
            throw std::runtime_error("Schwarz type must be multiplicative or additive.");
        }
        options_.standaloneMode = normalizeSchwarzTransmissionModeName(options_.standaloneMode);
        options_.transmission = schwarzTransmissionForStandaloneMode(options_.standaloneMode);
        options_.transmissionOrientation =
            normalizeSchwarzTransmissionOrientationName(options_.transmissionOrientation);
        if (!isKnownSchwarzStandaloneMode(options_.standaloneMode)) {
            throw std::runtime_error(
                "Schwarz standalone_mode must be algebraic, dirichlet_neumann, dirichlet_dirichlet, dirichlet_robin, or robin.");
        }
        if (!isKnownSchwarzTransmissionOrientation(options_.transmissionOrientation)) {
            throw std::runtime_error("Schwarz transmission_orientation must be forward or reverse.");
        }
        if (usesPhysicalTransmission()) {
            if (config_ == nullptr) {
                throw std::runtime_error("Physical Schwarz standalone modes require CaseConfig.");
            }
            if (robinBaseMatrix_ == nullptr) {
                throw std::runtime_error("Physical Schwarz standalone modes require the interface-free transmission matrix.");
            }
        }
        if (options_.fluxEval != "physical_gradient" && options_.fluxEval != "sipg_numeric") {
            throw std::runtime_error("Schwarz flux_eval must be physical_gradient or sipg_numeric.");
        }
        options_.partitionMode = normalizeSchwarzPartitionModeName(options_.partitionMode);
        if (!isKnownSchwarzPartitionMode(options_.partitionMode)) {
            throw std::runtime_error(
                "Schwarz partition mode must be current, material_aligned, hotspot_contained, or vertical_heat_flow_aligned.");
        }
        if (mesh.nodes.empty()) {
            throw std::runtime_error("Cannot build Schwarz blocks for an empty mesh.");
        }

        buildPartitionMaps(mesh);

        const auto haloBuildStart = std::chrono::steady_clock::now();
        localDofs_ = buildOverlapDofsFromHaloElements(mesh, std::max(0, options_.overlapLayers));
        localIndex_.assign(localDofs_.size(), {});
        for (size_t blockIndex = 0; blockIndex < localDofs_.size(); ++blockIndex) {
            auto& map = localIndex_[blockIndex];
            map.reserve(localDofs_[blockIndex].size());
            for (int local = 0; local < static_cast<int>(localDofs_[blockIndex].size()); ++local) {
                map[localDofs_[blockIndex][static_cast<size_t>(local)]] = local;
            }
        }
        haloBuildSeconds_ =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - haloBuildStart).count();

        const auto localMatrixAssemblyStart = std::chrono::steady_clock::now();
        offBlockEntries_.assign(blockDofs_.size(), {});
        globalNonzeros_ = a.values.size();
        diagonalBlockNonzeros_ = 0;
        couplingBlockNonzeros_ = 0;
        invalidMatrixEntries_ = 0;
        globalEntryL1Norm_ = 0.0;
        partitionedEntryL1Norm_ = 0.0;
        globalEntryFrobeniusSquared_ = 0.0;
        partitionedEntryFrobeniusSquared_ = 0.0;
        std::vector<std::vector<MatrixEntry>> diagonalEntries(blockDofs_.size());
        for (size_t blockIndex = 0; blockIndex < blockDofs_.size(); ++blockIndex) {
            diagonalEntries[blockIndex].reserve(localDofs_[blockIndex].size() * 32);
            offBlockEntries_[blockIndex].reserve(localDofs_[blockIndex].size() * 8);
        }

        for (size_t blockIndex = 0; blockIndex < blockDofs_.size(); ++blockIndex) {
            const int block = static_cast<int>(blockIndex);
            for (int globalRow : blockDofs_[blockIndex]) {
                for (int k = a.rowPtr[static_cast<size_t>(globalRow)];
                     k < a.rowPtr[static_cast<size_t>(globalRow + 1)];
                     ++k) {
                    const int globalCol = a.colInd[static_cast<size_t>(k)];
                    const double value = a.values[static_cast<size_t>(k)];
                    globalEntryL1Norm_ += std::abs(value);
                    globalEntryFrobeniusSquared_ += value * value;
                    if (globalCol < 0 || globalCol >= static_cast<int>(globalToBlock_.size())) {
                        ++invalidMatrixEntries_;
                        continue;
                    }
                    const int colBlock = globalToBlock_[static_cast<size_t>(globalCol)];
                    partitionedEntryL1Norm_ += std::abs(value);
                    partitionedEntryFrobeniusSquared_ += value * value;
                    if (colBlock == block) {
                        ++diagonalBlockNonzeros_;
                    } else {
                        ++couplingBlockNonzeros_;
                    }
                }
            }
        }
        for (size_t blockIndex = 0; blockIndex < blockDofs_.size(); ++blockIndex) {
            const auto& map = localIndex_[blockIndex];
            for (int globalRow : localDofs_[blockIndex]) {
                const int localRow = map.at(globalRow);
                for (int k = a.rowPtr[static_cast<size_t>(globalRow)];
                     k < a.rowPtr[static_cast<size_t>(globalRow + 1)];
                     ++k) {
                    const int globalCol = a.colInd[static_cast<size_t>(k)];
                    const auto found = map.find(globalCol);
                    if (found != map.end()) {
                        diagonalEntries[blockIndex].push_back({localRow,
                                                               found->second,
                                                               a.values[static_cast<size_t>(k)]});
                    } else {
                        offBlockEntries_[blockIndex].push_back({localRow,
                                                                globalCol,
                                                                a.values[static_cast<size_t>(k)]});
                    }
                }
            }
        }
        const double droppedL1 = std::abs(globalEntryL1Norm_ - partitionedEntryL1Norm_);
        const double droppedFrobeniusSquared =
            std::max(0.0, globalEntryFrobeniusSquared_ - partitionedEntryFrobeniusSquared_);
        reconstructionMaxAbsError_ = invalidMatrixEntries_ == 0 ? 0.0 : std::numeric_limits<double>::quiet_NaN();
        reconstructionL1Error_ = droppedL1;
        reconstructionFrobeniusError_ = std::sqrt(droppedFrobeniusSquared);

        std::vector<std::vector<MatrixEntry>> transmissionEntries = diagonalEntries;
        localDirichletColumnEntries_.clear();
        interfaceDirichletLocalRows_.clear();
        interfaceDirichletMask_.clear();
        interfaceDirichletConstraints_.clear();
        if (usesPhysicalTransmission()) {
            transmissionEntries = makeDiagonalBlockEntries(*robinBaseMatrix_);
            if (schwarzModeHasCondition(options_.standaloneMode,
                                        options_.transmissionOrientation,
                                        SchwarzInterfaceCondition::Robin)) {
                appendRobinMatrixEntries(mesh, *config_, transmissionEntries);
            }
            buildInterfaceDirichletConstraints(mesh);
            applyLocalDirichletConstraintsToMatrix(transmissionEntries);
        }
        localMatrixAssemblySeconds_ =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - localMatrixAssemblyStart).count();

        if (!factorLocalSolvers_) {
            solvers_.clear();
            robinSolvers_.clear();
            localFactorizationSeconds_ = 0.0;
            return;
        }

        solvers_.clear();
        solvers_.resize(blockDofs_.size());
        robinSolvers_.clear();
        if (usesPhysicalTransmission()) {
            robinSolvers_.resize(blockDofs_.size());
        }
        std::vector<std::string> factorErrors(blockDofs_.size());
        const auto localFactorizationStart = std::chrono::steady_clock::now();
        parallelForCoarse(blockDofs_.size(), [&](size_t blockIndex) {
            try {
                if (localDofs_[blockIndex].empty()) {
                    return;
                }
                ScopedMklSingleThread mklThreads;
                solvers_[blockIndex] =
                    SubdomainDirectSolver(static_cast<int>(localDofs_[blockIndex].size()),
                                          diagonalEntries[blockIndex]);
                if (usesPhysicalTransmission()) {
                    robinSolvers_[blockIndex] =
                        SubdomainDirectSolver(static_cast<int>(localDofs_[blockIndex].size()),
                                              transmissionEntries[blockIndex]);
                }
            } catch (const std::exception& err) {
                factorErrors[blockIndex] = err.what();
            }
        });
        localFactorizationSeconds_ =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - localFactorizationStart).count();
        for (size_t blockIndex = 0; blockIndex < factorErrors.size(); ++blockIndex) {
            if (!factorErrors[blockIndex].empty()) {
                throw std::runtime_error("Schwarz Aii block " + std::to_string(blockIndex)
                    + " factorization failed: " + factorErrors[blockIndex]);
            }
        }
    }

    std::vector<std::vector<MatrixEntry>> makeDiagonalBlockEntries(const SparseMatrix& matrix) const
    {
        if (!matrix.csrReady) {
            throw std::runtime_error("Schwarz robin base matrix must be CSR-finalized.");
        }
        if (matrix.size() != a_.size()) {
            throw std::runtime_error("Schwarz robin base matrix size does not match the SIPG system.");
        }
        std::vector<std::vector<MatrixEntry>> entries(blockDofs_.size());
        for (size_t blockIndex = 0; blockIndex < blockDofs_.size(); ++blockIndex) {
            const auto& map = localIndex_[blockIndex];
            entries[blockIndex].reserve(localDofs_[blockIndex].size() * 32);
            for (int globalRow : localDofs_[blockIndex]) {
                const int localRow = map.at(globalRow);
                for (int k = matrix.rowPtr[static_cast<size_t>(globalRow)];
                     k < matrix.rowPtr[static_cast<size_t>(globalRow + 1)];
                     ++k) {
                    const int globalCol = matrix.colInd[static_cast<size_t>(k)];
                    const auto found = map.find(globalCol);
                    if (found == map.end()) { continue; }
                    entries[blockIndex].push_back({
                        localRow,
                        found->second,
                        matrix.values[static_cast<size_t>(k)]
                    });
                }
            }
        }
        return entries;
    }

    void addDirichletConstraintsForSide(const Mesh& mesh,
                                        int targetBlock,
                                        const Tet& target,
                                        const std::array<int, 3>& targetFaceLocal,
                                        const Tet& source,
                                        const std::array<int, 3>& sourceFaceLocal)
    {
        if (targetBlock < 0 || targetBlock >= static_cast<int>(interfaceDirichletConstraints_.size())) {
            return;
        }
        const std::array<int, 6> faceDofs = p2FaceDofIndices(targetFaceLocal);
        for (int targetLocalDof : faceDofs) {
            SchwarzInterfaceDirichletConstraint constraint;
            constraint.globalDof = target.dof[static_cast<size_t>(targetLocalDof)];
            if (constraint.globalDof < 0
                || constraint.globalDof >= static_cast<int>(mesh.nodes.size())) {
                continue;
            }
            const Vec3 p = mesh.nodes[static_cast<size_t>(constraint.globalDof)].p;
            const auto sourceLambda = lambdaOnTetFace(p, source, sourceFaceLocal, mesh);
            constraint.sourceShape = shapeP2(sourceLambda);
            constraint.sourceDofs = source.dof;
            interfaceDirichletConstraints_[static_cast<size_t>(targetBlock)].push_back(constraint);
        }
    }

    void buildInterfaceDirichletConstraints(const Mesh& mesh)
    {
        interfaceDirichletConstraints_.assign(blockDofs_.size(), {});
        interfaceDirichletLocalRows_.assign(blockDofs_.size(), {});
        interfaceDirichletMask_.assign(blockDofs_.size(), {});
        if (!schwarzModeHasCondition(options_.standaloneMode,
                                     options_.transmissionOrientation,
                                     SchwarzInterfaceCondition::Dirichlet)) {
            return;
        }

        for (const InterfaceFace& face : mesh.interfaceFaces) {
            const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
            const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
            const int leftBlock = tetBlock(face.leftTet);
            const int rightBlock = tetBlock(face.rightTet);
            if (schwarzConditionForSide(options_.standaloneMode,
                                        options_.transmissionOrientation,
                                        true)
                == SchwarzInterfaceCondition::Dirichlet) {
                addDirichletConstraintsForSide(mesh,
                                               leftBlock,
                                               left,
                                               face.leftLocal,
                                               right,
                                               face.rightLocal);
            }
            if (schwarzConditionForSide(options_.standaloneMode,
                                        options_.transmissionOrientation,
                                        false)
                == SchwarzInterfaceCondition::Dirichlet) {
                addDirichletConstraintsForSide(mesh,
                                               rightBlock,
                                               right,
                                               face.rightLocal,
                                               left,
                                               face.leftLocal);
            }
        }

        for (size_t blockIndex = 0; blockIndex < blockDofs_.size(); ++blockIndex) {
            std::vector<char> mask(localDofs_[blockIndex].size(), 0);
            std::vector<int> rows;
            rows.reserve(interfaceDirichletConstraints_[blockIndex].size());
            for (const SchwarzInterfaceDirichletConstraint& constraint :
                 interfaceDirichletConstraints_[blockIndex]) {
                const auto found = localIndex_[blockIndex].find(constraint.globalDof);
                if (found == localIndex_[blockIndex].end()) {
                    continue;
                }
                const int local = found->second;
                if (local >= 0 && local < static_cast<int>(mask.size()) && !mask[static_cast<size_t>(local)]) {
                    mask[static_cast<size_t>(local)] = 1;
                    rows.push_back(local);
                }
            }
            std::sort(rows.begin(), rows.end());
            interfaceDirichletMask_[blockIndex] = std::move(mask);
            interfaceDirichletLocalRows_[blockIndex] = std::move(rows);
        }
    }

    void applyLocalDirichletConstraintsToMatrix(std::vector<std::vector<MatrixEntry>>& localEntries)
    {
        localDirichletColumnEntries_.assign(localEntries.size(), {});
        for (size_t blockIndex = 0; blockIndex < localEntries.size(); ++blockIndex) {
            if (blockIndex >= interfaceDirichletMask_.size()
                || interfaceDirichletLocalRows_[blockIndex].empty()) {
                continue;
            }
            const std::vector<char>& mask = interfaceDirichletMask_[blockIndex];
            std::vector<MatrixEntry> filtered;
            filtered.reserve(localEntries[blockIndex].size());
            localDirichletColumnEntries_[blockIndex].reserve(interfaceDirichletLocalRows_[blockIndex].size() * 16);
            for (const MatrixEntry& entry : localEntries[blockIndex]) {
                const bool rowDirichlet = entry.row >= 0
                    && entry.row < static_cast<int>(mask.size())
                    && mask[static_cast<size_t>(entry.row)];
                const bool colDirichlet = entry.col >= 0
                    && entry.col < static_cast<int>(mask.size())
                    && mask[static_cast<size_t>(entry.col)];
                if (!rowDirichlet && colDirichlet) {
                    localDirichletColumnEntries_[blockIndex].push_back(entry);
                }
                if (rowDirichlet || colDirichlet) {
                    continue;
                }
                filtered.push_back(entry);
            }
            for (int localRow : interfaceDirichletLocalRows_[blockIndex]) {
                filtered.push_back({localRow, localRow, 1.0});
            }
            localEntries[blockIndex] = std::move(filtered);
        }
    }

    void appendRobinMatrixEntries(const Mesh& mesh,
                                  const CaseConfig& config,
                                  std::vector<std::vector<MatrixEntry>>& localEntries) const
    {
        const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
        for (const InterfaceFace& face : mesh.interfaceFaces) {
            const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
            const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
            const ElementGeometry gLeft = elementGeometry(mesh, left);
            const ElementGeometry gRight = elementGeometry(mesh, right);
            const Vec3 normal = schwarzInterfaceNormal(mesh, face, left, right);
            const double alpha = schwarzRobinAlphaForFace(mesh, config, face, left, right, gLeft, gRight, normal);
            for (const auto& integrationTriangle : face.integrationTriangles) {
                const Vec3 p0 = integrationTriangle[0];
                const Vec3 p1 = integrationTriangle[1];
                const Vec3 p2 = integrationTriangle[2];
                const double jacFace = norm(cross(p1 - p0, p2 - p0));
                for (const TriangleQuadraturePoint& qp : quadrature) {
                    const double weight = qp.weight * jacFace;
                    const Vec3 q = (1.0 - qp.a - qp.b) * p0 + qp.a * p1 + qp.b * p2;
                    const auto lLeft = lambdaOnTetFace(q, left, face.leftLocal, mesh);
                    const auto lRight = lambdaOnTetFace(q, right, face.rightLocal, mesh);
                    const auto nLeft = shapeP2(lLeft);
                    const auto nRight = shapeP2(lRight);
                    const int leftBlock = tetBlock(face.leftTet);
                    const int rightBlock = tetBlock(face.rightTet);
                    if (leftBlock < 0 || rightBlock < 0
                        || leftBlock >= static_cast<int>(localIndex_.size())
                        || rightBlock >= static_cast<int>(localIndex_.size())) {
                        continue;
                    }
                    for (int a = 0; a < 10; ++a) {
                        const auto leftRowIt = localIndex_[static_cast<size_t>(leftBlock)].find(left.dof[static_cast<size_t>(a)]);
                        const auto rightRowIt = localIndex_[static_cast<size_t>(rightBlock)].find(right.dof[static_cast<size_t>(a)]);
                        if (leftRowIt == localIndex_[static_cast<size_t>(leftBlock)].end()
                            || rightRowIt == localIndex_[static_cast<size_t>(rightBlock)].end()) {
                            continue;
                        }
                        const int leftRow = leftRowIt->second;
                        const int rightRow = rightRowIt->second;
                        for (int b = 0; b < 10; ++b) {
                            const auto leftColIt = localIndex_[static_cast<size_t>(leftBlock)].find(left.dof[static_cast<size_t>(b)]);
                            const auto rightColIt = localIndex_[static_cast<size_t>(rightBlock)].find(right.dof[static_cast<size_t>(b)]);
                            if (leftColIt == localIndex_[static_cast<size_t>(leftBlock)].end()
                                || rightColIt == localIndex_[static_cast<size_t>(rightBlock)].end()) {
                                continue;
                            }
                            const int leftCol = leftColIt->second;
                            const int rightCol = rightColIt->second;
                            const double leftValue = alpha * nLeft[static_cast<size_t>(a)]
                                * nLeft[static_cast<size_t>(b)] * weight;
                            const double rightValue = alpha * nRight[static_cast<size_t>(a)]
                                * nRight[static_cast<size_t>(b)] * weight;
                            if (std::abs(leftValue) > 0.0) {
                                localEntries[static_cast<size_t>(leftBlock)].push_back({leftRow, leftCol, leftValue});
                            }
                            if (std::abs(rightValue) > 0.0) {
                                localEntries[static_cast<size_t>(rightBlock)].push_back({rightRow, rightCol, rightValue});
                            }
                        }
                    }
                }
            }
        }
    }

    std::vector<double> buildLocalRhs(size_t blockIndex,
                                      const std::vector<double>& rhs,
                                      const std::vector<double>& neighborTemperature) const
    {
        const std::vector<int>& dofs = localDofs_[blockIndex];
        std::vector<double> localRhs(dofs.size(), 0.0);
        for (size_t local = 0; local < dofs.size(); ++local) {
            localRhs[local] = rhs[static_cast<size_t>(dofs[local])];
        }
        for (const MatrixEntry& entry : offBlockEntries_[blockIndex]) {
            localRhs[static_cast<size_t>(entry.row)] -=
                entry.value * neighborTemperature[static_cast<size_t>(entry.col)];
        }
        return localRhs;
    }

    std::vector<double> interfaceDirichletValues(size_t blockIndex,
                                                 const std::vector<double>& neighborTemperature) const
    {
        const size_t localSize = blockIndex < localDofs_.size() ? localDofs_[blockIndex].size() : 0;
        std::vector<double> values(localSize, std::numeric_limits<double>::quiet_NaN());
        if (blockIndex >= interfaceDirichletConstraints_.size()
            || interfaceDirichletConstraints_[blockIndex].empty()) {
            return values;
        }
        std::vector<double> sums(localSize, 0.0);
        std::vector<int> counts(localSize, 0);
        for (const SchwarzInterfaceDirichletConstraint& constraint :
             interfaceDirichletConstraints_[blockIndex]) {
            const auto found = localIndex_[blockIndex].find(constraint.globalDof);
            if (found == localIndex_[blockIndex].end()) {
                continue;
            }
            double value = 0.0;
            for (int i = 0; i < 10; ++i) {
                const int sourceDof = constraint.sourceDofs[static_cast<size_t>(i)];
                if (sourceDof >= 0 && sourceDof < static_cast<int>(neighborTemperature.size())) {
                    value += constraint.sourceShape[static_cast<size_t>(i)]
                           * neighborTemperature[static_cast<size_t>(sourceDof)];
                }
            }
            const int local = found->second;
            if (local >= 0 && local < static_cast<int>(localSize)) {
                sums[static_cast<size_t>(local)] += value;
                ++counts[static_cast<size_t>(local)];
            }
        }
        for (size_t local = 0; local < localSize; ++local) {
            if (counts[local] > 0) {
                values[local] = sums[local] / static_cast<double>(counts[local]);
            }
        }
        return values;
    }

    bool isInterfaceDirichletLocalRow(size_t blockIndex, int localRow) const
    {
        return blockIndex < interfaceDirichletMask_.size()
            && localRow >= 0
            && localRow < static_cast<int>(interfaceDirichletMask_[blockIndex].size())
            && interfaceDirichletMask_[blockIndex][static_cast<size_t>(localRow)] != 0;
    }

    void addPhysicalFluxTransmissionRhs(size_t blockIndex,
                                        const std::vector<double>& neighborTemperature,
                                        std::vector<double>& localRhs) const
    {
        if (mesh_ == nullptr || config_ == nullptr) {
            throw std::runtime_error("Physical Schwarz RHS requires mesh and case config.");
        }
        const Mesh& mesh = *mesh_;
        const CaseConfig& config = *config_;
        const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
        for (const InterfaceFace& face : mesh.interfaceFaces) {
            const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
            const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
            const int leftBlock = tetBlock(face.leftTet);
            const int rightBlock = tetBlock(face.rightTet);
            const bool useLeft = static_cast<int>(blockIndex) == leftBlock;
            const bool useRight = static_cast<int>(blockIndex) == rightBlock;
            const SchwarzInterfaceCondition leftCondition =
                schwarzConditionForSide(options_.standaloneMode,
                                        options_.transmissionOrientation,
                                        true);
            const SchwarzInterfaceCondition rightCondition =
                schwarzConditionForSide(options_.standaloneMode,
                                        options_.transmissionOrientation,
                                        false);
            if ((!useLeft || (leftCondition != SchwarzInterfaceCondition::Neumann
                              && leftCondition != SchwarzInterfaceCondition::Robin))
                && (!useRight || (rightCondition != SchwarzInterfaceCondition::Neumann
                                  && rightCondition != SchwarzInterfaceCondition::Robin))) {
                continue;
            }
            const ElementGeometry gLeft = elementGeometry(mesh, left);
            const ElementGeometry gRight = elementGeometry(mesh, right);
            const Vec3 normal = schwarzInterfaceNormal(mesh, face, left, right);
            const Vec3 rightNormal = schwarzRightInterfaceNormal(normal, face);
            const double alpha = schwarzRobinAlphaForFace(mesh, config, face, left, right, gLeft, gRight, normal);
            for (const auto& integrationTriangle : face.integrationTriangles) {
                const Vec3 p0 = integrationTriangle[0];
                const Vec3 p1 = integrationTriangle[1];
                const Vec3 p2 = integrationTriangle[2];
                const double jacFace = norm(cross(p1 - p0, p2 - p0));
                for (const TriangleQuadraturePoint& qp : quadrature) {
                    const double weight = qp.weight * jacFace;
                    const Vec3 q = (1.0 - qp.a - qp.b) * p0 + qp.a * p1 + qp.b * p2;
                    const auto lLeft = lambdaOnTetFace(q, left, face.leftLocal, mesh);
                    const auto lRight = lambdaOnTetFace(q, right, face.rightLocal, mesh);
                    const auto nLeft = shapeP2(lLeft);
                    const auto nRight = shapeP2(lRight);
                    const SchwarzTracePointData point =
                        schwarzTracePointData(config, left, right, gLeft, gRight,
                                              lLeft, lRight, normal, rightNormal, alpha, neighborTemperature);
                    if (useLeft
                        && (leftCondition == SchwarzInterfaceCondition::Neumann
                            || leftCondition == SchwarzInterfaceCondition::Robin)) {
                        const double gTransmission = leftCondition == SchwarzInterfaceCondition::Robin
                            ? point.qRight + alpha * point.tRight
                            : point.qRight;
                        for (int a = 0; a < 10; ++a) {
                            const auto rowIt = localIndex_[blockIndex].find(left.dof[static_cast<size_t>(a)]);
                            if (rowIt == localIndex_[blockIndex].end()) { continue; }
                            const int localRow = rowIt->second;
                            if (isInterfaceDirichletLocalRow(blockIndex, localRow)) { continue; }
                            localRhs[static_cast<size_t>(localRow)] +=
                                gTransmission * nLeft[static_cast<size_t>(a)] * weight;
                        }
                    }
                    if (useRight
                        && (rightCondition == SchwarzInterfaceCondition::Neumann
                            || rightCondition == SchwarzInterfaceCondition::Robin)) {
                        const double gTransmission = rightCondition == SchwarzInterfaceCondition::Robin
                            ? point.qLeft + alpha * point.tLeft
                            : point.qLeft;
                        for (int a = 0; a < 10; ++a) {
                            const auto rowIt = localIndex_[blockIndex].find(right.dof[static_cast<size_t>(a)]);
                            if (rowIt == localIndex_[blockIndex].end()) { continue; }
                            const int localRow = rowIt->second;
                            if (isInterfaceDirichletLocalRow(blockIndex, localRow)) { continue; }
                            localRhs[static_cast<size_t>(localRow)] +=
                                gTransmission * nRight[static_cast<size_t>(a)] * weight;
                        }
                    }
                }
            }
        }
    }

    std::vector<double> buildPhysicalTransmissionLocalRhs(
        size_t blockIndex,
        const std::vector<double>& rhs,
        const std::vector<double>& neighborTemperature) const
    {
        const std::vector<int>& dofs = localDofs_[blockIndex];
        std::vector<double> localRhs(dofs.size(), 0.0);
        for (size_t local = 0; local < dofs.size(); ++local) {
            localRhs[local] = rhs[static_cast<size_t>(dofs[local])];
        }

        const std::vector<double> dirichletValues =
            interfaceDirichletValues(blockIndex, neighborTemperature);
        if (blockIndex < localDirichletColumnEntries_.size()) {
            for (const MatrixEntry& entry : localDirichletColumnEntries_[blockIndex]) {
                if (entry.row < 0 || entry.row >= static_cast<int>(localRhs.size())
                    || entry.col < 0 || entry.col >= static_cast<int>(dirichletValues.size())) {
                    continue;
                }
                const double value = dirichletValues[static_cast<size_t>(entry.col)];
                if (std::isfinite(value)) {
                    localRhs[static_cast<size_t>(entry.row)] -= entry.value * value;
                }
            }
        }

        addPhysicalFluxTransmissionRhs(blockIndex, neighborTemperature, localRhs);

        if (blockIndex < interfaceDirichletLocalRows_.size()) {
            for (int localRow : interfaceDirichletLocalRows_[blockIndex]) {
                if (localRow < 0 || localRow >= static_cast<int>(localRhs.size())) {
                    continue;
                }
                const double value = localRow < static_cast<int>(dirichletValues.size())
                    ? dirichletValues[static_cast<size_t>(localRow)]
                    : std::numeric_limits<double>::quiet_NaN();
                localRhs[static_cast<size_t>(localRow)] =
                    std::isfinite(value) ? value : neighborTemperature[static_cast<size_t>(dofs[static_cast<size_t>(localRow)])];
            }
        }

        return localRhs;
    }

    std::vector<double> buildRobinLocalRhs(size_t blockIndex,
                                          const std::vector<double>& rhs,
                                          const std::vector<double>& neighborTemperature) const
    {
        if (mesh_ == nullptr || config_ == nullptr) {
            throw std::runtime_error("Robin Schwarz RHS requires mesh and case config.");
        }
        const Mesh& mesh = *mesh_;
        const CaseConfig& config = *config_;
        const std::vector<int>& dofs = localDofs_[blockIndex];
        std::vector<double> localRhs(dofs.size(), 0.0);
        for (size_t local = 0; local < dofs.size(); ++local) {
            localRhs[local] = rhs[static_cast<size_t>(dofs[local])];
        }

        const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
        for (const InterfaceFace& face : mesh.interfaceFaces) {
            const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
            const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
            const bool useLeft = static_cast<int>(blockIndex) == left.subdomain;
            const bool useRight = static_cast<int>(blockIndex) == right.subdomain;
            if (!useLeft && !useRight) {
                continue;
            }
            const ElementGeometry gLeft = elementGeometry(mesh, left);
            const ElementGeometry gRight = elementGeometry(mesh, right);
            const Vec3 normal = schwarzInterfaceNormal(mesh, face, left, right);
            const Vec3 rightNormal = schwarzRightInterfaceNormal(normal, face);
            const double alpha = schwarzRobinAlphaForFace(mesh, config, face, left, right, gLeft, gRight, normal);
            for (const auto& integrationTriangle : face.integrationTriangles) {
                const Vec3 p0 = integrationTriangle[0];
                const Vec3 p1 = integrationTriangle[1];
                const Vec3 p2 = integrationTriangle[2];
                const double jacFace = norm(cross(p1 - p0, p2 - p0));
                for (const TriangleQuadraturePoint& qp : quadrature) {
                    const double weight = qp.weight * jacFace;
                    const Vec3 q = (1.0 - qp.a - qp.b) * p0 + qp.a * p1 + qp.b * p2;
                    const auto lLeft = lambdaOnTetFace(q, left, face.leftLocal, mesh);
                    const auto lRight = lambdaOnTetFace(q, right, face.rightLocal, mesh);
                    const auto nLeft = shapeP2(lLeft);
                    const auto nRight = shapeP2(lRight);
                    const SchwarzTracePointData point =
                        schwarzTracePointData(config, left, right, gLeft, gRight,
                                              lLeft, lRight, normal, rightNormal, alpha, neighborTemperature);
                    if (useLeft) {
                        const double gTransmission = point.qRight + alpha * point.tRight;
                        for (int a = 0; a < 10; ++a) {
                            const auto rowIt = localIndex_[blockIndex].find(left.dof[static_cast<size_t>(a)]);
                            if (rowIt == localIndex_[blockIndex].end()) { continue; }
                            const int localRow = rowIt->second;
                            localRhs[static_cast<size_t>(localRow)] +=
                                gTransmission * nLeft[static_cast<size_t>(a)] * weight;
                        }
                    }
                    if (useRight) {
                        const double gTransmission = point.qLeft + alpha * point.tLeft;
                        for (int a = 0; a < 10; ++a) {
                            const auto rowIt = localIndex_[blockIndex].find(right.dof[static_cast<size_t>(a)]);
                            if (rowIt == localIndex_[blockIndex].end()) { continue; }
                            const int localRow = rowIt->second;
                            localRhs[static_cast<size_t>(localRow)] +=
                                gTransmission * nRight[static_cast<size_t>(a)] * weight;
                        }
                    }
                }
            }
        }
        return localRhs;
    }

    static double relativeUpdateNorm(const std::vector<double>& current,
                                     const std::vector<double>& previous)
    {
        double numerator = 0.0;
        double denominator = 0.0;
        for (size_t i = 0; i < current.size(); ++i) {
            const double diff = current[i] - previous[i];
            numerator += diff * diff;
            denominator += current[i] * current[i];
        }
        return std::sqrt(numerator) / std::sqrt(std::max(1.0e-300, denominator));
    }

    double relativeSystemResidual(const std::vector<double>& x,
                                  const std::vector<double>& rhs,
                                  double rhsNorm) const
    {
        const std::vector<double> ax = a_.multiply(x);
        double sum = 0.0;
        for (size_t i = 0; i < rhs.size(); ++i) {
            const double r = ax[i] - rhs[i];
            sum += r * r;
        }
        return std::sqrt(sum) / rhsNorm;
    }

    double relativeFreeSystemResidual(const std::vector<double>& x,
                                      const std::vector<double>& rhs,
                                      double& absResidual) const
    {
        const std::vector<double> ax = a_.multiply(x);
        double sum = 0.0;
        double rhsSum = 0.0;
        for (size_t i = 0; i < rhs.size(); ++i) {
            if (mesh_ != nullptr
                && i < mesh_->nodes.size()
                && mesh_->nodes[i].dirichlet) {
                continue;
            }
            const double r = ax[i] - rhs[i];
            sum += r * r;
            rhsSum += rhs[i] * rhs[i];
        }
        absResidual = std::sqrt(sum);
        return absResidual / std::sqrt(std::max(1.0e-300, rhsSum));
    }
};

class TwoLevelSchwarzPreconditioner {
public:
    TwoLevelSchwarzPreconditioner(SchwarzBlockSolver& localPreconditioner,
                                  const Mesh& mesh,
                                  const SparseMatrix& a,
                                  SubdomainConstantCoarseSpace coarseSpace,
                                  std::string correctionMode)
        : local_(localPreconditioner),
          a_(a),
          coarseSpace_(std::move(coarseSpace)),
          correctionMode_(std::move(correctionMode))
    {
        (void)mesh;
        const auto setupStart = std::chrono::steady_clock::now();
        buildCoarseMatrix();
        coarseSetupSeconds_ =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - setupStart).count();
    }

    TwoLevelSchwarzPreconditioner(SchwarzBlockSolver& localPreconditioner,
                                  const Mesh& mesh,
                                  const SparseMatrix& a)
        : TwoLevelSchwarzPreconditioner(localPreconditioner,
                                        mesh,
                                        a,
                                        buildSubdomainConstantCoarseSpace(mesh),
                                        "additive")
    {
    }

    void apply(const std::vector<double>& r, std::vector<double>& z)
    {
        if (coarseSpace_.dofsByMode.empty()) {
            local_.apply(r, z);
            return;
        }

        std::vector<double> coarseCorrection;
        if (!coarseCorrectionForResidual(r, coarseCorrection)) {
            local_.apply(r, z);
            return;
        }

        if (correctionMode_ == "multiplicative") {
            const std::vector<double> azc = a_.multiply(coarseCorrection);
            std::vector<double> residualAfterCoarse(r.size(), 0.0);
            parallelFor(r.size(), [&](size_t i) {
                residualAfterCoarse[i] = r[i] - azc[i];
            });
            std::vector<double> localCorrection;
            local_.apply(residualAfterCoarse, localCorrection);
            z.resize(r.size());
            parallelFor(r.size(), [&](size_t i) {
                z[i] = coarseCorrection[i] + localCorrection[i];
            });
            accumulateCorrectionNorms(coarseCorrection, localCorrection);
        } else {
            local_.apply(r, z);
            accumulateCorrectionNorms(coarseCorrection, z);
            parallelFor(z.size(), [&](size_t i) {
                z[i] += coarseCorrection[i];
            });
        }
    }

    size_t memoryBytes() const
    {
        size_t bytes = local_.memoryBytes() + coarseSpace_.memoryBytes();
        for (const auto& row : coarseMatrix_) {
            bytes += row.size() * sizeof(double);
        }
        return bytes;
    }

    void copyTimingTo(SolverStatistics& stats) const
    {
        local_.copyTimingTo(stats);
        stats.coarseEnabled = coarseDim() > 0;
        stats.coarseDim = coarseDim();
        stats.coarseSetupSeconds = coarseSetupSeconds_;
        stats.coarseSolveSeconds = coarseSolveSeconds_;
        stats.coarseMatrixNnz = coarseMatrixNnz_;
        stats.coarseResidualNorm = coarseResidualNorm_;
        stats.coarseRhsNorm = coarseRhsNorm_;
        stats.coarseSolutionNorm = coarseSolutionNorm_;
        stats.coarseCorrectionNorm = coarseCorrectionNorm_;
        stats.localCorrectionNorm = localCorrectionNorm_;
        stats.coarseToLocalNormRatio = coarseToLocalNormRatio_;
    }

    int coarseDim() const { return coarseSpace_.dim(); }
    double coarseSetupSeconds() const { return coarseSetupSeconds_; }
    double coarseSolveSeconds() const { return coarseSolveSeconds_; }
    size_t coarseMatrixNnz() const { return coarseMatrixNnz_; }
    double coarseResidualNorm() const { return coarseResidualNorm_; }

private:
    SchwarzBlockSolver& local_;
    const SparseMatrix& a_;
    SubdomainConstantCoarseSpace coarseSpace_;
    std::string correctionMode_ = "additive";
    std::vector<std::vector<double>> coarseMatrix_;
    size_t coarseMatrixNnz_ = 0;
    double coarseSetupSeconds_ = 0.0;
    double coarseSolveSeconds_ = 0.0;
    double coarseResidualNorm_ = std::numeric_limits<double>::quiet_NaN();
    double coarseRhsNorm_ = std::numeric_limits<double>::quiet_NaN();
    double coarseSolutionNorm_ = std::numeric_limits<double>::quiet_NaN();
    double coarseCorrectionNorm_ = std::numeric_limits<double>::quiet_NaN();
    double localCorrectionNorm_ = std::numeric_limits<double>::quiet_NaN();
    double coarseToLocalNormRatio_ = std::numeric_limits<double>::quiet_NaN();
    double coarseCorrectionNormSum_ = 0.0;
    double localCorrectionNormSum_ = 0.0;

    void buildCoarseMatrix()
    {
        const size_t dim = coarseSpace_.dofsByMode.size();
        coarseMatrix_.assign(dim, std::vector<double>(dim, 0.0));
        coarseMatrixNnz_ = 0;
        if (dim == 0) {
            return;
        }
        for (int row = 0; row < a_.n; ++row) {
            const int rowMode = coarseSpace_.globalToMode[static_cast<size_t>(row)];
            if (rowMode < 0) {
                continue;
            }
            const double rowScale = coarseSpace_.invNormByMode[static_cast<size_t>(rowMode)];
            for (int k = a_.rowPtr[static_cast<size_t>(row)];
                 k < a_.rowPtr[static_cast<size_t>(row + 1)];
                 ++k) {
                const int col = a_.colInd[static_cast<size_t>(k)];
                if (col < 0 || col >= static_cast<int>(coarseSpace_.globalToMode.size())) {
                    continue;
                }
                const int colMode = coarseSpace_.globalToMode[static_cast<size_t>(col)];
                if (colMode < 0) {
                    continue;
                }
                const double colScale = coarseSpace_.invNormByMode[static_cast<size_t>(colMode)];
                coarseMatrix_[static_cast<size_t>(rowMode)][static_cast<size_t>(colMode)] +=
                    rowScale * a_.values[static_cast<size_t>(k)] * colScale;
            }
        }
        for (const auto& row : coarseMatrix_) {
            for (double value : row) {
                if (std::abs(value) > 0.0) { ++coarseMatrixNnz_; }
            }
        }
    }

    bool coarseCorrectionForResidual(const std::vector<double>& r,
                                     std::vector<double>& correction)
    {
        const auto coarseSolveStart = std::chrono::steady_clock::now();
        std::vector<double> rhs(coarseSpace_.dofsByMode.size(), 0.0);
        for (size_t mode = 0; mode < coarseSpace_.dofsByMode.size(); ++mode) {
            double sum = 0.0;
            for (int dof : coarseSpace_.dofsByMode[mode]) {
                sum += r[static_cast<size_t>(dof)];
            }
            rhs[mode] = coarseSpace_.invNormByMode[mode] * sum;
        }
        std::vector<double> coeffs;
        if (!solveSmallDenseSystem(coarseMatrix_, rhs, coeffs)) {
            coarseResidualNorm_ = std::numeric_limits<double>::infinity();
            coarseSolveSeconds_ +=
                std::chrono::duration<double>(std::chrono::steady_clock::now() - coarseSolveStart).count();
            return false;
        }
        correction.assign(r.size(), 0.0);
        for (size_t mode = 0; mode < coarseSpace_.dofsByMode.size(); ++mode) {
            const double alpha = coeffs[mode] * coarseSpace_.invNormByMode[mode];
            for (int dof : coarseSpace_.dofsByMode[mode]) {
                correction[static_cast<size_t>(dof)] += alpha;
            }
        }
        coarseRhsNorm_ = vectorNorm(rhs);
        coarseSolutionNorm_ = vectorNorm(coeffs);
        coarseResidualNorm_ = denseResidualNorm(coarseMatrix_, coeffs, rhs);
        coarseSolveSeconds_ +=
            std::chrono::duration<double>(std::chrono::steady_clock::now() - coarseSolveStart).count();
        return true;
    }

    void accumulateCorrectionNorms(const std::vector<double>& coarseCorrection,
                                   const std::vector<double>& localCorrection)
    {
        coarseCorrectionNorm_ = vectorNorm(coarseCorrection);
        localCorrectionNorm_ = vectorNorm(localCorrection);
        coarseCorrectionNormSum_ += coarseCorrectionNorm_;
        localCorrectionNormSum_ += localCorrectionNorm_;
        coarseToLocalNormRatio_ =
            coarseCorrectionNormSum_ / std::max(1.0e-300, localCorrectionNormSum_);
    }

    static double denseResidualNorm(const std::vector<std::vector<double>>& matrix,
                                    const std::vector<double>& x,
                                    const std::vector<double>& rhs)
    {
        double sum = 0.0;
        for (size_t row = 0; row < matrix.size(); ++row) {
            double ax = 0.0;
            for (size_t col = 0; col < x.size(); ++col) {
                ax += matrix[row][col] * x[col];
            }
            const double residual = ax - rhs[row];
            sum += residual * residual;
        }
        return std::sqrt(sum);
    }
};

struct KrylovIterationRow {
    std::string solver;
    int timeStep = 0;
    double time = 0.0;
    int iteration = 0;
    int restartCycle = 0;
    int innerIteration = 0;
    double absoluteResidual = std::numeric_limits<double>::quiet_NaN();
    double relativeResidual = std::numeric_limits<double>::quiet_NaN();
    double trueRelativeResidual = std::numeric_limits<double>::quiet_NaN();
    double preconditionedResidual = std::numeric_limits<double>::quiet_NaN();
    double preconditionerTime = 0.0;
    double globalMatvecTime = 0.0;
    double totalSolveTime = 0.0;
};

static std::vector<double> trueResidualVectorTimed(const SparseMatrix& a,
                                                   const std::vector<double>& x,
                                                   const std::vector<double>& b,
                                                   double& matvecSeconds)
{
    const auto matvecStart = std::chrono::steady_clock::now();
    const std::vector<double> ax = a.multiply(x);
    const auto matvecEnd = std::chrono::steady_clock::now();
    matvecSeconds += std::chrono::duration<double>(matvecEnd - matvecStart).count();
    std::vector<double> r(b.size(), 0.0);
    parallelFor(b.size(), [&](size_t i) {
        r[i] = b[i] - ax[i];
    });
    return r;
}

static std::vector<double> fgmresCandidate(const std::vector<double>& base,
                                           const std::vector<std::vector<double>>& z,
                                           const std::vector<std::vector<double>>& h,
                                           const std::vector<double>& g,
                                           int activeColumns)
{
    std::vector<double> y(static_cast<size_t>(activeColumns), 0.0);
    for (int i = activeColumns - 1; i >= 0; --i) {
        double sum = g[static_cast<size_t>(i)];
        for (int j = i + 1; j < activeColumns; ++j) {
            sum -= h[static_cast<size_t>(i)][static_cast<size_t>(j)] * y[static_cast<size_t>(j)];
        }
        const double diag = h[static_cast<size_t>(i)][static_cast<size_t>(i)];
        if (!(std::abs(diag) > 0.0) || !std::isfinite(diag)) {
            throw std::runtime_error("Schwarz-preconditioned FGMRES singular least-squares system.");
        }
        y[static_cast<size_t>(i)] = sum / diag;
    }

    std::vector<double> candidate = base;
    for (int i = 0; i < activeColumns; ++i) {
        axpy(y[static_cast<size_t>(i)], z[static_cast<size_t>(i)], candidate);
    }
    return candidate;
}

template <typename Preconditioner>
static std::vector<double> flexibleGmresWithKrylovLog(const SparseMatrix& a,
                                                      const std::vector<double>& b,
                                                      std::vector<double> x,
                                                      Preconditioner& preconditioner,
                                                      int& iterations,
                                                      std::vector<KrylovIterationRow>& logRows,
                                                      int timeStep,
                                                      double time,
                                                      SolverStatistics* stats,
                                                      const std::string& solverLabel,
                                                      int maxIterations,
                                                      double relativeTolerance,
                                                      int restart)
{
    const int n = a.size();
    const int restartLength = std::max(1, restart);
    const double rhsNorm = std::max(1.0e-300, vectorNorm(b));
    iterations = 0;
    const auto solveStart = std::chrono::steady_clock::now();

    double initialMatvecSeconds = 0.0;
    std::vector<double> r = trueResidualVectorTimed(a, x, b, initialMatvecSeconds);
    double beta = vectorNorm(r);
    double relativeResidual = beta / rhsNorm;
    if (stats != nullptr) {
        stats->matvecSeconds += initialMatvecSeconds;
        recordInitialResidualDiagnostics(stats, a, b, x, relativeTolerance);
    }
    logRows.push_back({solverLabel,
                       timeStep,
                       time,
                       0,
                       0,
                       0,
                       beta,
                       relativeResidual,
                       relativeResidual,
                       std::numeric_limits<double>::quiet_NaN(),
                       0.0,
                       initialMatvecSeconds,
                       std::chrono::duration<double>(std::chrono::steady_clock::now() - solveStart).count()});
    if (stats != nullptr) {
        stats->pcgResidualNorm = beta;
        stats->finalRelativeResidual = relativeResidual;
    }
    if (relativeResidual <= relativeTolerance) {
        if (stats != nullptr && stats->status != "failed") {
            stats->status = "success";
        }
        return x;
    }

    int restartCycle = 0;
    while (iterations < maxIterations) {
        ++restartCycle;
        const std::vector<double> xBase = x;
        std::vector<std::vector<double>> v(static_cast<size_t>(restartLength + 1),
                                           std::vector<double>(static_cast<size_t>(n), 0.0));
        std::vector<std::vector<double>> z(static_cast<size_t>(restartLength),
                                           std::vector<double>(static_cast<size_t>(n), 0.0));
        std::vector<std::vector<double>> h(static_cast<size_t>(restartLength + 1),
                                           std::vector<double>(static_cast<size_t>(restartLength), 0.0));
        std::vector<double> cs(static_cast<size_t>(restartLength), 0.0);
        std::vector<double> sn(static_cast<size_t>(restartLength), 0.0);
        std::vector<double> g(static_cast<size_t>(restartLength + 1), 0.0);

        parallelFor(v[0].size(), [&](size_t i) {
            v[0][i] = r[i] / beta;
        });
        g[0] = beta;

        int innerCount = 0;
        std::vector<double> latestCandidate = x;
        std::vector<double> latestResidual = r;
        double latestResidualNorm = beta;
        bool haveLatest = false;

        for (; innerCount < restartLength && iterations < maxIterations; ++innerCount) {
            const auto preconditionerStart = std::chrono::steady_clock::now();
            preconditioner.apply(v[static_cast<size_t>(innerCount)], z[static_cast<size_t>(innerCount)]);
            const auto preconditionerEnd = std::chrono::steady_clock::now();
            const double preconditionerSeconds =
                std::chrono::duration<double>(preconditionerEnd - preconditionerStart).count();
            const double preconditionedResidual =
                vectorNorm(z[static_cast<size_t>(innerCount)]) / rhsNorm;
            if (stats != nullptr) {
                stats->preconditionerApplySeconds += preconditionerSeconds;
                ++stats->preconditionerApplyCalls;
            }
            if (vectorHasNonFinite(z[static_cast<size_t>(innerCount)])) {
                if (stats != nullptr) {
                    stats->status = "failed";
                    stats->failureReason = solverLabel + " breakdown: preconditioner produced NaN/Inf";
                    stats->breakdownIteration = iterations;
                }
                return x;
            }

            double iterationMatvecSeconds = 0.0;
            const auto matvecStart = std::chrono::steady_clock::now();
            std::vector<double> w = a.multiply(z[static_cast<size_t>(innerCount)]);
            const auto matvecEnd = std::chrono::steady_clock::now();
            iterationMatvecSeconds += std::chrono::duration<double>(matvecEnd - matvecStart).count();
            if (stats != nullptr) {
                stats->matvecSeconds += iterationMatvecSeconds;
            }
            if (vectorHasNonFinite(w)) {
                if (stats != nullptr) {
                    stats->status = "failed";
                    stats->failureReason = solverLabel + " breakdown: A*z produced NaN/Inf";
                    stats->breakdownIteration = iterations;
                }
                return x;
            }

            for (int i = 0; i <= innerCount; ++i) {
                h[static_cast<size_t>(i)][static_cast<size_t>(innerCount)] =
                    vectorDot(w, v[static_cast<size_t>(i)]);
                axpy(-h[static_cast<size_t>(i)][static_cast<size_t>(innerCount)],
                     v[static_cast<size_t>(i)],
                     w);
            }
            h[static_cast<size_t>(innerCount + 1)][static_cast<size_t>(innerCount)] = vectorNorm(w);
            if (h[static_cast<size_t>(innerCount + 1)][static_cast<size_t>(innerCount)] > 0.0) {
                parallelFor(w.size(), [&](size_t i) {
                    v[static_cast<size_t>(innerCount + 1)][i] =
                        w[i] / h[static_cast<size_t>(innerCount + 1)][static_cast<size_t>(innerCount)];
                });
            }

            for (int i = 0; i < innerCount; ++i) {
                const double hi = h[static_cast<size_t>(i)][static_cast<size_t>(innerCount)];
                const double hip1 = h[static_cast<size_t>(i + 1)][static_cast<size_t>(innerCount)];
                h[static_cast<size_t>(i)][static_cast<size_t>(innerCount)] =
                    cs[static_cast<size_t>(i)] * hi + sn[static_cast<size_t>(i)] * hip1;
                h[static_cast<size_t>(i + 1)][static_cast<size_t>(innerCount)] =
                    -sn[static_cast<size_t>(i)] * hi + cs[static_cast<size_t>(i)] * hip1;
            }

            const double hii = h[static_cast<size_t>(innerCount)][static_cast<size_t>(innerCount)];
            const double hnext = h[static_cast<size_t>(innerCount + 1)][static_cast<size_t>(innerCount)];
            const double denom = std::hypot(hii, hnext);
            if (!(denom > 0.0) || !std::isfinite(denom)) {
                if (stats != nullptr) {
                    stats->status = "failed";
                    stats->failureReason = solverLabel + " breakdown: Arnoldi Hessenberg pivot is zero/nonfinite";
                    stats->breakdownIteration = iterations;
                }
                return x;
            }
            cs[static_cast<size_t>(innerCount)] = hii / denom;
            sn[static_cast<size_t>(innerCount)] = hnext / denom;
            h[static_cast<size_t>(innerCount)][static_cast<size_t>(innerCount)] = denom;
            h[static_cast<size_t>(innerCount + 1)][static_cast<size_t>(innerCount)] = 0.0;

            const double gi = g[static_cast<size_t>(innerCount)];
            g[static_cast<size_t>(innerCount)] = cs[static_cast<size_t>(innerCount)] * gi;
            g[static_cast<size_t>(innerCount + 1)] = -sn[static_cast<size_t>(innerCount)] * gi;

            ++iterations;
            try {
                latestCandidate = fgmresCandidate(xBase, z, h, g, innerCount + 1);
            } catch (const std::exception& err) {
                if (stats != nullptr) {
                    stats->status = "failed";
                    stats->failureReason = err.what();
                    stats->breakdownIteration = iterations;
                }
                return x;
            }
            latestResidual = trueResidualVectorTimed(a, latestCandidate, b, iterationMatvecSeconds);
            latestResidualNorm = vectorNorm(latestResidual);
            relativeResidual = latestResidualNorm / rhsNorm;
            haveLatest = true;

            if (stats != nullptr) {
                stats->breakdownIteration = iterations;
                stats->pcgResidualNorm = latestResidualNorm;
                stats->finalRelativeResidual = relativeResidual;
            }
            logRows.push_back({solverLabel,
                               timeStep,
                               time,
                               iterations,
                               restartCycle,
                               innerCount + 1,
                               latestResidualNorm,
                               relativeResidual,
                               relativeResidual,
                               preconditionedResidual,
                               preconditionerSeconds,
                               iterationMatvecSeconds,
                               std::chrono::duration<double>(std::chrono::steady_clock::now() - solveStart).count()});
            if (relativeResidual <= relativeTolerance) {
                if (stats != nullptr && stats->status != "failed") {
                    stats->status = "success";
                }
                return latestCandidate;
            }
        }

        if (haveLatest) {
            x = std::move(latestCandidate);
            r = std::move(latestResidual);
            beta = latestResidualNorm;
        } else {
            double matvecSeconds = 0.0;
            r = trueResidualVectorTimed(a, x, b, matvecSeconds);
            beta = vectorNorm(r);
        }
        if (!(beta > 0.0) || !std::isfinite(beta)) {
            if (stats != nullptr) {
                stats->status = "failed";
                stats->failureReason = solverLabel + " breakdown: restarted residual is zero/nonfinite";
                stats->breakdownIteration = iterations;
            }
            return x;
        }
    }

    if (stats != nullptr) {
        stats->finalRelativeResidual = relativeResidualNorm(a, x, b);
        if (stats->status != "failed") {
            stats->status = stats->finalRelativeResidual <= relativeTolerance ? "success" : "failed";
        }
        if (stats->status == "failed" && stats->failureReason.empty()) {
            stats->failureReason = solverLabel + " did not converge";
        }
    }
    return x;
}

struct SchwarzRunResult {
    std::vector<double> temperature;
    std::vector<SchwarzStepStats> stepStats;
    std::vector<SchwarzSubdomainSolveStats> subdomainStats;
    std::vector<SchwarzInterfaceFluxRow> interfaceFluxRows;
};

static SchwarzRunResult runSchwarzAnalysis(bool steady,
                                           const std::string& solverName,
                                           const Mesh& mesh,
                                           const CaseConfig& physics,
                                           const SparseMatrix& mass,
                                           const SparseMatrix& system,
                                           const std::vector<double>& source,
                                           const std::vector<double>& fixedAdjust,
                                           SchwarzBlockSolver& schwarzSolver,
                                           double blockAssemblySeconds,
                                           SolverStatistics& stats)
{
    stats.name = solverName;
    stats.totalIterations = 0;
    stats.maxIterations = 0;
    stats.parallelWorkers = 1;
    stats.initialGuessType = physics.initialGuessType;
    (void)system;
    if (stats.workingSetBeforeBytes == 0) {
        stats.workingSetBeforeBytes = currentWorkingSetBytes();
    }
    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());

    SchwarzRunResult run;
    run.temperature = initialTemperatureVector(mesh, physics);
    const auto solveStart = std::chrono::steady_clock::now();
    bool allConverged = true;

    const int stepBegin = steady ? 0 : 1;
    const int stepEnd = steady ? 0 : physics.timeSteps;
    for (int step = stepBegin; step <= stepEnd; ++step) {
        std::vector<double> rhs = steady
            ? source
            : makeTransientRhs(mass, run.temperature, source, physics.timeStep);
        applyDirichletRhs(mesh, fixedAdjust, rhs);
        recordRhsDiagnostics(&stats, rhs, physics, fixedAdjust);
        const std::vector<double> x0 = (!steady && physics.disableWarmStart)
            ? initialTemperatureVector(mesh, physics)
            : run.temperature;

        const double time = steady ? physics.startTime : step * physics.timeStep;
        const double stepAssemblySeconds = (step == stepBegin) ? blockAssemblySeconds : 0.0;
        SchwarzStepResult stepResult = schwarzSolver.solveStep(mesh,
                                                               rhs,
                                                               x0,
                                                               step,
                                                               time,
                                                               stepAssemblySeconds);
        run.temperature = std::move(stepResult.temperature);
        stats.totalIterations += stepResult.stats.iterations;
        stats.maxIterations = std::max(stats.maxIterations, stepResult.stats.iterations);
        stats.finalRelativeResidual = std::isfinite(stepResult.stats.freeRelResidual)
            ? stepResult.stats.freeRelResidual
            : stepResult.stats.relResidual;
        stats.trueRelativeResidual = std::isfinite(stepResult.stats.globalRelResidual)
            ? stepResult.stats.globalRelResidual
            : stepResult.stats.relResidual;
        allConverged = allConverged && stepResult.stats.converged;
        run.subdomainStats.insert(run.subdomainStats.end(),
                                  stepResult.subdomainStats.begin(),
                                  stepResult.subdomainStats.end());
        run.interfaceFluxRows.insert(run.interfaceFluxRows.end(),
                                     stepResult.interfaceFluxRows.begin(),
                                     stepResult.interfaceFluxRows.end());
        run.stepStats.push_back(stepResult.stats);

        if (steady || step == 1 || step % 10 == 0 || step == physics.timeSteps) {
            const auto minmax = std::minmax_element(run.temperature.begin(), run.temperature.end());
            const double avg = std::accumulate(run.temperature.begin(), run.temperature.end(), 0.0)
                / static_cast<double>(run.temperature.size());
            std::cout << "[" << solverName << "] "
                      << (steady ? "steady" : ("step " + std::to_string(step)))
                      << "  time=" << std::setw(10) << time
                      << " s  Tmin=" << std::setw(12) << *minmax.first
                      << "  Tmax=" << std::setw(12) << *minmax.second
                      << "  Tavg=" << std::setw(12) << avg
                      << "  interface_avg_jump=" << stepResult.stats.interfaceJump
                      << "  iterations=" << stepResult.stats.iterations
                      << "  rel_update=" << stepResult.stats.relUpdate
                      << "  rel_residual=" << stepResult.stats.relResidual
                      << "  free_rel_residual=" << stepResult.stats.freeRelResidual << "\n";
        }
    }

    stats.solveSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - solveStart).count();
    stats.workingSetAfterBytes = currentWorkingSetBytes();
    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
    if (stats.status == "not_run") {
        if (vectorHasNonFinite(run.temperature)) {
            stats.status = "failed";
            stats.failureReason = "solution contains NaN/Inf";
        } else if (!allConverged) {
            stats.status = "failed";
            stats.failureReason = "Schwarz iteration did not satisfy both convergence criteria on every step";
        } else {
            stats.status = "success";
        }
    }
    return run;
}

struct SchwarzPreconditionedKrylovRunResult {
    std::vector<double> temperature;
    std::vector<double> finalRhs;
    std::vector<KrylovIterationRow> iterationRows;
};

template <typename Preconditioner>
static SchwarzPreconditionedKrylovRunResult runSchwarzPreconditionedFgmresAnalysis(
    bool steady,
    const std::string& solverName,
    const Mesh& mesh,
    const CaseConfig& physics,
    const SparseMatrix& mass,
    const SparseMatrix& system,
    const std::vector<double>& source,
    const std::vector<double>& fixedAdjust,
    Preconditioner& preconditioner,
    SolverStatistics& stats,
    int maxIterations,
    double relativeTolerance,
    int restart)
{
    stats.name = solverName;
    stats.totalIterations = 0;
    stats.maxIterations = 0;
    stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
    stats.initialGuessType = physics.initialGuessType;
    if (stats.workingSetBeforeBytes == 0) {
        stats.workingSetBeforeBytes = currentWorkingSetBytes();
    }
    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());

    SchwarzPreconditionedKrylovRunResult run;
    run.temperature = initialTemperatureVector(mesh, physics);
    const auto solveStart = std::chrono::steady_clock::now();
    bool allStepsConverged = true;

    const int stepBegin = steady ? 0 : 1;
    const int stepEnd = steady ? 0 : physics.timeSteps;
    for (int step = stepBegin; step <= stepEnd; ++step) {
        std::vector<double> rhs = steady
            ? source
            : makeTransientRhs(mass, run.temperature, source, physics.timeStep);
        applyDirichletRhs(mesh, fixedAdjust, rhs);
        run.finalRhs = rhs;
        recordRhsDiagnostics(&stats, rhs, physics, fixedAdjust);

        const double physicalTime = steady ? physics.startTime : step * physics.timeStep;
        int iterations = 0;
        std::vector<double> x0 = (!steady && physics.disableWarmStart)
            ? initialTemperatureVector(mesh, physics)
            : run.temperature;
        run.temperature = flexibleGmresWithKrylovLog(system,
                                                     rhs,
                                                     std::move(x0),
                                                     preconditioner,
                                                     iterations,
                                                     run.iterationRows,
                                                     step,
                                                     physicalTime,
                                                     &stats,
                                                     solverName,
                                                     maxIterations,
                                                     relativeTolerance,
                                                     restart);
        stats.totalIterations += iterations;
        stats.maxIterations = std::max(stats.maxIterations, iterations);
        const double stepRelResidual = relativeResidualNorm(system, run.temperature, rhs);
        allStepsConverged = allStepsConverged
            && std::isfinite(stepRelResidual)
            && stepRelResidual <= relativeTolerance;
        stats.finalRelativeResidual = stepRelResidual;
        stats.trueRelativeResidual = stepRelResidual;

        if (steady || step == 1 || step % 10 == 0 || step == physics.timeSteps) {
            const auto minmax = std::minmax_element(run.temperature.begin(), run.temperature.end());
            const double avg = std::accumulate(run.temperature.begin(), run.temperature.end(), 0.0)
                / static_cast<double>(run.temperature.size());
            std::cout << "[" << solverName << "] "
                      << (steady ? "steady" : ("step " + std::to_string(step)))
                      << "  time=" << std::setw(10) << physicalTime
                      << " s  Tmin=" << std::setw(12) << *minmax.first
                      << "  Tmax=" << std::setw(12) << *minmax.second
                      << "  Tavg=" << std::setw(12) << avg
                      << "  interface_avg_jump=" << interfaceAverageJump(mesh, run.temperature)
                      << "  iterations=" << iterations
                      << "  rel_residual=" << stepRelResidual << "\n";
        }
    }

    stats.solveSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - solveStart).count();
    stats.workingSetAfterBytes = currentWorkingSetBytes();
    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
    if (vectorHasNonFinite(run.temperature)) {
        stats.status = "failed";
        stats.failureReason = "solution contains NaN/Inf";
    } else if (!allStepsConverged) {
        stats.status = "failed";
        if (stats.failureReason.empty()) {
            stats.failureReason = "Schwarz-preconditioned FGMRES did not converge on every step";
        }
    } else if (stats.status != "failed") {
        stats.status = "success";
    }
    return run;
}

static void writeSchwarzIterationStats(const std::vector<SchwarzStepStats>& rows,
                                       const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "time_step,time,iterations,rel_update,rel_residual,global_rel_residual,"
        << "free_rel_residual,free_abs_residual,interface_jump,"
        << "interface_jump_l2,flux_balance_l2,assembly_time,solve_time,coupling_time,total_time\n";
    out << std::setprecision(16);
    for (const SchwarzStepStats& row : rows) {
        out << row.timeStep << ','
            << row.time << ','
            << row.iterations << ','
            << row.relUpdate << ','
            << row.relResidual << ','
            << row.globalRelResidual << ','
            << row.freeRelResidual << ','
            << row.freeAbsResidual << ','
            << row.interfaceJump << ','
            << row.interfaceJumpL2 << ','
            << row.fluxBalanceL2 << ','
            << row.assemblyTime << ','
            << row.solveTime << ','
            << row.couplingTime << ','
            << row.totalTime << '\n';
    }
}

static void writeSchwarzInterfaceFluxRows(const std::vector<SchwarzInterfaceFluxRow>& rows,
                                          const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "time_step,schwarz_iter,interface_id,left_sd,right_sd,mean_T_left,mean_T_right,"
        << "max_jump_T,l2_jump_T,mean_flux_left,mean_flux_right,flux_balance_l2,flux_balance_max\n";
    out << std::setprecision(16);
    for (const SchwarzInterfaceFluxRow& row : rows) {
        out << row.timeStep << ','
            << row.schwarzIter << ','
            << row.interfaceId << ','
            << row.leftSubdomain << ','
            << row.rightSubdomain << ','
            << row.meanTLeft << ','
            << row.meanTRight << ','
            << row.maxJumpT << ','
            << row.l2JumpT << ','
            << row.meanFluxLeft << ','
            << row.meanFluxRight << ','
            << row.fluxBalanceL2 << ','
            << row.fluxBalanceMax << '\n';
    }
}

static void writeKrylovIterationLog(const std::vector<KrylovIterationRow>& rows,
                                    const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "solver,time_step,time,iteration,restart_cycle,inner_iteration,"
        << "absolute_residual,relative_residual,true_relative_residual,preconditioned_residual,preconditioner_time,"
        << "global_matvec_time,total_solve_time\n";
    out << std::setprecision(16);
    for (const KrylovIterationRow& row : rows) {
        out << csvEscape(row.solver) << ','
            << row.timeStep << ','
            << row.time << ','
            << row.iteration << ','
            << row.restartCycle << ','
            << row.innerIteration << ','
            << row.absoluteResidual << ','
            << row.relativeResidual << ','
            << row.trueRelativeResidual << ','
            << row.preconditionedResidual << ','
            << row.preconditionerTime << ','
            << row.globalMatvecTime << ','
            << row.totalSolveTime << '\n';
    }
}

static void writeSchwarzSubdomainSolveStats(const std::vector<SchwarzSubdomainSolveStats>& rows,
                                            const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "time_step,time,subdomain,dofs,solve_time\n";
    out << std::setprecision(16);
    for (const SchwarzSubdomainSolveStats& row : rows) {
        out << row.timeStep << ','
            << row.time << ','
            << row.subdomain << ','
            << row.dofs << ','
            << row.solveTime << '\n';
    }
}

static void writeSchwarzBlockMatrixCheck(const SchwarzBlockSolver& solver,
                                         const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "metric,value\n";
    out << "subdomain_count," << solver.blockCount() << '\n';
    out << "overlap_mode," << csvEscape(solver.overlapMode()) << '\n';
    out << "overlap_layers," << solver.overlapLayers() << '\n';
    out << "partition_mode," << csvEscape(solver.partitionMode()) << '\n';
    for (size_t block = 0; block < solver.blockCount(); ++block) {
        out << "block_" << block << "_owned_elements," << solver.ownedElements(block) << '\n';
        out << "block_" << block << "_local_elements," << solver.localElements(block) << '\n';
        out << "block_" << block << "_halo_elements," << solver.haloElements(block) << '\n';
        out << "block_" << block << "_owned_dofs," << solver.coreDofs(block) << '\n';
        out << "block_" << block << "_local_dofs," << solver.localDofs(block) << '\n';
        out << "block_" << block << "_halo_dofs," << solver.haloDofs(block) << '\n';
    }
    out << "global_nnz," << solver.globalNonzeros() << '\n';
    out << "partitioned_nnz," << solver.partitionedNonzeros() << '\n';
    out << "diagonal_block_nnz," << solver.diagonalBlockNonzeros() << '\n';
    out << "coupling_block_nnz," << solver.couplingBlockNonzeros() << '\n';
    out << "invalid_matrix_entries," << solver.invalidMatrixEntries() << '\n';
    out << "global_entry_l1_norm," << std::setprecision(16) << solver.globalEntryL1Norm() << '\n';
    out << "partitioned_entry_l1_norm," << std::setprecision(16) << solver.partitionedEntryL1Norm() << '\n';
    out << "global_entry_frobenius_norm," << std::setprecision(16) << solver.globalEntryFrobeniusNorm() << '\n';
    out << "partitioned_entry_frobenius_norm," << std::setprecision(16) << solver.partitionedEntryFrobeniusNorm() << '\n';
    out << "reconstruction_max_abs_error," << std::setprecision(16) << solver.reconstructionMaxAbsError() << '\n';
    out << "reconstruction_l1_error," << std::setprecision(16) << solver.reconstructionL1Error() << '\n';
    out << "reconstruction_frobenius_error," << std::setprecision(16) << solver.reconstructionFrobeniusError() << '\n';
    out << "reconstruction_note,"
        << csvEscape("blocks are partitioned directly from the finalized global SIPG matrix entries")
        << '\n';
}

struct RramTetDiagnostic {
    double sourcePower = 0.0;
    double sourceDensity = 0.0;
    double gradientMagnitude = std::numeric_limits<double>::quiet_NaN();
    double lengthScale = 0.0;
};

struct RramPhysicalOverlapStats {
    int layer = 0;
    int boundaryPairCount = 0;
    double minThickness = std::numeric_limits<double>::quiet_NaN();
    double meanThickness = std::numeric_limits<double>::quiet_NaN();
    double maxThickness = std::numeric_limits<double>::quiet_NaN();
};

struct RramInterfaceDiagnosticRow {
    int interfaceId = -1;
    int leftSubdomain = -1;
    int rightSubdomain = -1;
    std::string leftMaterial;
    std::string rightMaterial;
    double kLeft = std::numeric_limits<double>::quiet_NaN();
    double kRight = std::numeric_limits<double>::quiet_NaN();
    double kRatio = std::numeric_limits<double>::quiet_NaN();
    double interfaceArea = 0.0;
    double maxHeatSourceNearInterface = 0.0;
    double maxGradTNearInterface = std::numeric_limits<double>::quiet_NaN();
    double overlapPhysicalThicknessLayer1 = std::numeric_limits<double>::quiet_NaN();
    double overlapPhysicalThicknessLayer2 = std::numeric_limits<double>::quiet_NaN();
    bool cutsFilamentOrActiveOxide = false;
    bool cutsMaxHeatSourceRegion = false;
    bool cutsMaxTemperatureGradientRegion = false;
    bool cutsElectrodeOxideInterface = false;
};

static double rramScalarConductivity(const Material& material)
{
    return (material.conductivityX + material.conductivityY + material.conductivityZ) / 3.0;
}

static std::string rramMaterialId(const Material& material, int domainEntity)
{
    return material.name + "(" + std::to_string(domainEntity) + ")";
}

static std::string joinStrings(const std::set<std::string>& values, const std::string& delimiter)
{
    std::ostringstream out;
    bool first = true;
    for (const std::string& value : values) {
        if (!first) {
            out << delimiter;
        }
        first = false;
        out << value;
    }
    return out.str();
}

static bool rramMaterialIsFilamentOrActiveOxide(const std::string& name)
{
    const std::string value = lowerString(name);
    return value == "cf"
        || value.find("filament") != std::string::npos
        || value.find("hfo") != std::string::npos
        || value.find("tio2") != std::string::npos
        || value.find("oxide") != std::string::npos;
}

static bool rramMaterialIsElectrodeLike(const std::string& name)
{
    const std::string value = lowerString(name);
    return value.find("wl_bl") != std::string::npos
        || value.find("electrode") != std::string::npos
        || value.find("nickel") != std::string::npos
        || value == "ni"
        || value == "ti"
        || value == "pt";
}

static double rramTetGradientMagnitude(const Mesh& mesh,
                                       const Tet& tet,
                                       const std::vector<double>& temperature)
{
    if (temperature.size() != mesh.nodes.size()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const std::array<double, 4> lambda{{0.25, 0.25, 0.25, 0.25}};
    const ElementGeometry geo = elementGeometry(mesh, tet);
    const auto grad = gradShapeP2(lambda, geo);
    Vec3 gradT{};
    for (int i = 0; i < 10; ++i) {
        const int dof = tet.dof[static_cast<size_t>(i)];
        gradT = gradT + temperature[static_cast<size_t>(dof)] * grad[static_cast<size_t>(i)];
    }
    return norm(gradT);
}

static std::vector<RramTetDiagnostic> collectRramTetDiagnostics(const Mesh& mesh,
                                                                const CaseConfig& physics,
                                                                const std::vector<double>& temperature)
{
    std::vector<RramTetDiagnostic> diagnostics(mesh.tets.size());
    std::vector<double> volumes(mesh.tets.size(), 0.0);
    for (size_t tetId = 0; tetId < mesh.tets.size(); ++tetId) {
        const double volume = std::max(0.0, tetVolume(mesh, mesh.tets[tetId]));
        volumes[tetId] = volume;
        diagnostics[tetId].lengthScale = std::pow(std::max(1.0e-300, 6.0 * volume), 1.0 / 3.0);
        diagnostics[tetId].gradientMagnitude =
            rramTetGradientMagnitude(mesh, mesh.tets[tetId], temperature);
    }

    for (const HeatSource& source : physics.heatSources) {
        const double sourceVolume = selectedSourceVolume(mesh, source);
        if (!(sourceVolume > 0.0) || !std::isfinite(sourceVolume)) {
            continue;
        }
        const double density = source.heatRateW / sourceVolume;
        for (size_t tetId = 0; tetId < mesh.tets.size(); ++tetId) {
            if (!tetMatchesHeatSource(mesh.tets[tetId], source)) {
                continue;
            }
            diagnostics[tetId].sourceDensity += density;
            diagnostics[tetId].sourcePower += density * volumes[tetId];
        }
    }
    return diagnostics;
}

static std::vector<std::vector<int>> buildRramTetAdjacency(const Mesh& mesh)
{
    std::vector<std::set<int>> adjacencySets(mesh.tets.size());
    std::unordered_map<int, std::vector<int>> tetsByDof;
    tetsByDof.reserve(mesh.nodes.size());
    for (int tetId = 0; tetId < static_cast<int>(mesh.tets.size()); ++tetId) {
        for (int dof : mesh.tets[static_cast<size_t>(tetId)].dof) {
            tetsByDof[dof].push_back(tetId);
        }
    }
    for (const auto& item : tetsByDof) {
        const std::vector<int>& incident = item.second;
        for (int a : incident) {
            for (int b : incident) {
                if (a != b) {
                    adjacencySets[static_cast<size_t>(a)].insert(b);
                }
            }
        }
    }
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        if (face.leftTet >= 0 && face.rightTet >= 0
            && face.leftTet < static_cast<int>(mesh.tets.size())
            && face.rightTet < static_cast<int>(mesh.tets.size())) {
            adjacencySets[static_cast<size_t>(face.leftTet)].insert(face.rightTet);
            adjacencySets[static_cast<size_t>(face.rightTet)].insert(face.leftTet);
        }
    }

    std::vector<std::vector<int>> adjacency(adjacencySets.size());
    for (size_t i = 0; i < adjacencySets.size(); ++i) {
        adjacency[i].assign(adjacencySets[i].begin(), adjacencySets[i].end());
    }
    return adjacency;
}

static RramPhysicalOverlapStats computeRramPhysicalOverlapStats(const Mesh& mesh,
                                                                const SchwarzBlockSolver& solver,
                                                                int layer,
                                                                const std::vector<RramTetDiagnostic>& tetDiagnostics)
{
    RramPhysicalOverlapStats stats;
    stats.layer = layer;
    if (layer <= 0 || mesh.tets.empty()) {
        stats.minThickness = 0.0;
        stats.meanThickness = 0.0;
        stats.maxThickness = 0.0;
        return stats;
    }

    const std::vector<std::vector<int>> adjacency = buildRramTetAdjacency(mesh);
    double sum = 0.0;
    for (int tetId = 0; tetId < static_cast<int>(adjacency.size()); ++tetId) {
        const int block = solver.tetBlock(tetId);
        for (int neighbor : adjacency[static_cast<size_t>(tetId)]) {
            if (neighbor <= tetId) {
                continue;
            }
            if (block == solver.tetBlock(neighbor)) {
                continue;
            }
            const double h = 0.5 * (tetDiagnostics[static_cast<size_t>(tetId)].lengthScale
                                  + tetDiagnostics[static_cast<size_t>(neighbor)].lengthScale);
            const double thickness = static_cast<double>(layer) * h;
            if (stats.boundaryPairCount == 0) {
                stats.minThickness = thickness;
                stats.maxThickness = thickness;
            } else {
                stats.minThickness = std::min(stats.minThickness, thickness);
                stats.maxThickness = std::max(stats.maxThickness, thickness);
            }
            sum += thickness;
            ++stats.boundaryPairCount;
        }
    }
    if (stats.boundaryPairCount > 0) {
        stats.meanThickness = sum / static_cast<double>(stats.boundaryPairCount);
    }
    return stats;
}

static void writeRramPhysicalOverlapStats(const Mesh& mesh,
                                          const SchwarzBlockSolver& solver,
                                          const std::vector<RramTetDiagnostic>& tetDiagnostics,
                                          const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "partition_mode,overlap_layers,min_overlap_thickness,mean_overlap_thickness,"
        << "max_overlap_thickness,boundary_pair_count\n";
    out << std::setprecision(16);
    for (int layer = 1; layer <= 3; ++layer) {
        const RramPhysicalOverlapStats stats =
            computeRramPhysicalOverlapStats(mesh, solver, layer, tetDiagnostics);
        out << csvEscape(solver.partitionMode()) << ','
            << layer << ','
            << stats.minThickness << ','
            << stats.meanThickness << ','
            << stats.maxThickness << ','
            << stats.boundaryPairCount << '\n';
    }
}

static void writeRramPartitionDiagnostics(const Mesh& mesh,
                                          const CaseConfig& physics,
                                          const SchwarzBlockSolver& solver,
                                          const std::vector<RramTetDiagnostic>& tetDiagnostics,
                                          const std::filesystem::path& path)
{
    struct BlockAggregate {
        std::set<std::string> materialIds;
        double minK = std::numeric_limits<double>::max();
        double maxK = -std::numeric_limits<double>::max();
        double totalHeatSource = 0.0;
        double maxHeatSource = 0.0;
        double maxGradient = std::numeric_limits<double>::quiet_NaN();
    };

    std::vector<BlockAggregate> blocks(solver.blockCount());
    for (int tetId = 0; tetId < static_cast<int>(mesh.tets.size()); ++tetId) {
        const int block = solver.tetBlock(tetId);
        if (block < 0 || block >= static_cast<int>(blocks.size())) {
            continue;
        }
        const Tet& tet = mesh.tets[static_cast<size_t>(tetId)];
        const Material& material = materialForTet(physics, tet);
        const double k = rramScalarConductivity(material);
        BlockAggregate& aggregate = blocks[static_cast<size_t>(block)];
        aggregate.materialIds.insert(rramMaterialId(material, tet.domainEntity));
        aggregate.minK = std::min(aggregate.minK, k);
        aggregate.maxK = std::max(aggregate.maxK, k);
        aggregate.totalHeatSource += tetDiagnostics[static_cast<size_t>(tetId)].sourcePower;
        aggregate.maxHeatSource = std::max(aggregate.maxHeatSource,
                                           tetDiagnostics[static_cast<size_t>(tetId)].sourceDensity);
        const double gradient = tetDiagnostics[static_cast<size_t>(tetId)].gradientMagnitude;
        if (std::isfinite(gradient)) {
            aggregate.maxGradient = std::isfinite(aggregate.maxGradient)
                ? std::max(aggregate.maxGradient, gradient)
                : gradient;
        }
    }

    std::ofstream out(path);
    out << "subdomain_id,owned_elements,halo_elements,owned_dofs,halo_dofs,material_ids,"
        << "min_k,max_k,k_contrast,total_heat_source,max_heat_source,max_temperature_gradient\n";
    out << std::setprecision(16);
    for (size_t block = 0; block < blocks.size(); ++block) {
        const BlockAggregate& aggregate = blocks[block];
        const double minK = aggregate.minK == std::numeric_limits<double>::max()
            ? std::numeric_limits<double>::quiet_NaN()
            : aggregate.minK;
        const double maxK = aggregate.maxK == -std::numeric_limits<double>::max()
            ? std::numeric_limits<double>::quiet_NaN()
            : aggregate.maxK;
        const double contrast = std::isfinite(minK) && minK > 0.0 && std::isfinite(maxK)
            ? maxK / minK
            : std::numeric_limits<double>::quiet_NaN();
        out << block << ','
            << solver.ownedElements(block) << ','
            << solver.haloElements(block) << ','
            << solver.coreDofs(block) << ','
            << solver.haloDofs(block) << ','
            << csvEscape(joinStrings(aggregate.materialIds, "|")) << ','
            << minK << ','
            << maxK << ','
            << contrast << ','
            << aggregate.totalHeatSource << ','
            << aggregate.maxHeatSource << ','
            << aggregate.maxGradient << '\n';
    }
}

static std::string dominantWeightedKey(const std::map<std::string, double>& weights)
{
    std::string best;
    double bestWeight = -1.0;
    for (const auto& item : weights) {
        if (item.second > bestWeight) {
            best = item.first;
            bestWeight = item.second;
        }
    }
    return best;
}

static std::vector<RramInterfaceDiagnosticRow> collectRramInterfaceDiagnostics(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::vector<RramTetDiagnostic>& tetDiagnostics)
{
    struct Aggregate {
        int left = -1;
        int right = -1;
        double area = 0.0;
        double kLeftAreaSum = 0.0;
        double kRightAreaSum = 0.0;
        double hSum = 0.0;
        int hCount = 0;
        double maxSource = 0.0;
        double maxGradient = std::numeric_limits<double>::quiet_NaN();
        std::map<std::string, double> leftMaterialWeights;
        std::map<std::string, double> rightMaterialWeights;
        bool cutsFilamentOrActiveOxide = false;
        bool cutsElectrodeOxideInterface = false;
    };

    std::map<std::pair<int, int>, Aggregate> aggregates;
    double globalMaxSource = 0.0;
    double globalMaxGradient = 0.0;
    for (const RramTetDiagnostic& diag : tetDiagnostics) {
        globalMaxSource = std::max(globalMaxSource, diag.sourceDensity);
        if (std::isfinite(diag.gradientMagnitude)) {
            globalMaxGradient = std::max(globalMaxGradient, diag.gradientMagnitude);
        }
    }

    for (const InterfaceFace& face : mesh.interfaceFaces) {
        if (face.leftTet < 0 || face.rightTet < 0
            || face.leftTet >= static_cast<int>(mesh.tets.size())
            || face.rightTet >= static_cast<int>(mesh.tets.size())) {
            continue;
        }
        const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
        const Material& leftMaterial = materialForTet(physics, left);
        const Material& rightMaterial = materialForTet(physics, right);
        const Vec3 normal = schwarzInterfaceNormal(mesh, face, left, right);
        double area = integrationArea(face.integrationTriangles);
        if (!(area > 0.0)) {
            area = face.overlapArea;
        }
        if (!(area > 0.0)) {
            area = 1.0;
        }

        const std::pair<int, int> key{left.subdomain, right.subdomain};
        Aggregate& aggregate = aggregates[key];
        aggregate.left = left.subdomain;
        aggregate.right = right.subdomain;
        aggregate.area += area;
        aggregate.kLeftAreaSum += normalConductivity(leftMaterial, normal) * area;
        aggregate.kRightAreaSum += normalConductivity(rightMaterial, normal) * area;
        aggregate.leftMaterialWeights[rramMaterialId(leftMaterial, left.domainEntity)] += area;
        aggregate.rightMaterialWeights[rramMaterialId(rightMaterial, right.domainEntity)] += area;
        const double h = 0.5 * (tetDiagnostics[static_cast<size_t>(face.leftTet)].lengthScale
                              + tetDiagnostics[static_cast<size_t>(face.rightTet)].lengthScale);
        aggregate.hSum += h;
        ++aggregate.hCount;
        aggregate.maxSource = std::max(aggregate.maxSource,
            std::max(tetDiagnostics[static_cast<size_t>(face.leftTet)].sourceDensity,
                     tetDiagnostics[static_cast<size_t>(face.rightTet)].sourceDensity));
        const double leftGradient = tetDiagnostics[static_cast<size_t>(face.leftTet)].gradientMagnitude;
        const double rightGradient = tetDiagnostics[static_cast<size_t>(face.rightTet)].gradientMagnitude;
        const double faceGradient = std::max(std::isfinite(leftGradient) ? leftGradient : 0.0,
                                             std::isfinite(rightGradient) ? rightGradient : 0.0);
        aggregate.maxGradient = std::isfinite(aggregate.maxGradient)
            ? std::max(aggregate.maxGradient, faceGradient)
            : faceGradient;
        const bool leftActive = rramMaterialIsFilamentOrActiveOxide(leftMaterial.name);
        const bool rightActive = rramMaterialIsFilamentOrActiveOxide(rightMaterial.name);
        const bool leftElectrode = rramMaterialIsElectrodeLike(leftMaterial.name);
        const bool rightElectrode = rramMaterialIsElectrodeLike(rightMaterial.name);
        aggregate.cutsFilamentOrActiveOxide =
            aggregate.cutsFilamentOrActiveOxide || leftActive || rightActive;
        aggregate.cutsElectrodeOxideInterface =
            aggregate.cutsElectrodeOxideInterface
            || (leftActive && rightElectrode)
            || (rightActive && leftElectrode);
    }

    std::vector<RramInterfaceDiagnosticRow> rows;
    rows.reserve(aggregates.size());
    int interfaceId = 0;
    for (const auto& item : aggregates) {
        const Aggregate& aggregate = item.second;
        RramInterfaceDiagnosticRow row;
        row.interfaceId = interfaceId++;
        row.leftSubdomain = aggregate.left;
        row.rightSubdomain = aggregate.right;
        row.leftMaterial = dominantWeightedKey(aggregate.leftMaterialWeights);
        row.rightMaterial = dominantWeightedKey(aggregate.rightMaterialWeights);
        row.interfaceArea = aggregate.area;
        row.kLeft = aggregate.kLeftAreaSum / std::max(1.0e-300, aggregate.area);
        row.kRight = aggregate.kRightAreaSum / std::max(1.0e-300, aggregate.area);
        row.kRatio = std::max(row.kLeft, row.kRight) / std::max(1.0e-300, std::min(row.kLeft, row.kRight));
        row.maxHeatSourceNearInterface = aggregate.maxSource;
        row.maxGradTNearInterface = aggregate.maxGradient;
        const double layer1 = aggregate.hCount > 0
            ? aggregate.hSum / static_cast<double>(aggregate.hCount)
            : std::numeric_limits<double>::quiet_NaN();
        row.overlapPhysicalThicknessLayer1 = layer1;
        row.overlapPhysicalThicknessLayer2 = std::isfinite(layer1) ? 2.0 * layer1
            : std::numeric_limits<double>::quiet_NaN();
        row.cutsFilamentOrActiveOxide = aggregate.cutsFilamentOrActiveOxide;
        row.cutsMaxHeatSourceRegion = globalMaxSource > 0.0
            && row.maxHeatSourceNearInterface >= 0.9 * globalMaxSource;
        row.cutsMaxTemperatureGradientRegion = globalMaxGradient > 0.0
            && std::isfinite(row.maxGradTNearInterface)
            && row.maxGradTNearInterface >= 0.9 * globalMaxGradient;
        row.cutsElectrodeOxideInterface = aggregate.cutsElectrodeOxideInterface;
        rows.push_back(std::move(row));
    }
    return rows;
}

static void writeRramInterfaceDiagnostics(const std::vector<RramInterfaceDiagnosticRow>& rows,
                                          const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "interface_id,left_subdomain,right_subdomain,left_material,right_material,"
        << "k_left,k_right,k_ratio,interface_area,max_heat_source_near_interface,"
        << "max_grad_T_near_interface,overlap_physical_thickness_layer1,"
        << "overlap_physical_thickness_layer2\n";
    out << std::setprecision(16);
    for (const RramInterfaceDiagnosticRow& row : rows) {
        out << row.interfaceId << ','
            << row.leftSubdomain << ','
            << row.rightSubdomain << ','
            << csvEscape(row.leftMaterial) << ','
            << csvEscape(row.rightMaterial) << ','
            << row.kLeft << ','
            << row.kRight << ','
            << row.kRatio << ','
            << row.interfaceArea << ','
            << row.maxHeatSourceNearInterface << ','
            << row.maxGradTNearInterface << ','
            << row.overlapPhysicalThicknessLayer1 << ','
            << row.overlapPhysicalThicknessLayer2 << '\n';
    }
}

static void writeRramInterfaceCutChecks(const std::vector<RramInterfaceDiagnosticRow>& rows,
                                        const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "interface_id,left_subdomain,right_subdomain,cuts_filament_or_active_oxide,"
        << "cuts_max_heat_source_region,cuts_max_temperature_gradient_region,"
        << "cuts_electrode_oxide_interface,notes\n";
    for (const RramInterfaceDiagnosticRow& row : rows) {
        std::vector<std::string> notes;
        if (row.cutsFilamentOrActiveOxide) {
            notes.push_back("filament_or_active_oxide_adjacent");
        }
        if (row.cutsMaxHeatSourceRegion) {
            notes.push_back("near_global_max_heat_source");
        }
        if (row.cutsMaxTemperatureGradientRegion) {
            notes.push_back("near_global_max_temperature_gradient");
        }
        if (row.cutsElectrodeOxideInterface) {
            notes.push_back("electrode_oxide_material_pair");
        }
        std::set<std::string> noteSet(notes.begin(), notes.end());
        out << row.interfaceId << ','
            << row.leftSubdomain << ','
            << row.rightSubdomain << ','
            << (row.cutsFilamentOrActiveOxide ? 1 : 0) << ','
            << (row.cutsMaxHeatSourceRegion ? 1 : 0) << ','
            << (row.cutsMaxTemperatureGradientRegion ? 1 : 0) << ','
            << (row.cutsElectrodeOxideInterface ? 1 : 0) << ','
            << csvEscape(joinStrings(noteSet, "|")) << '\n';
    }
}

static double rramMinimumLayerThickness(const Mesh& mesh)
{
    double minThickness = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < mesh.subdomainMin.size() && i < mesh.subdomainMax.size(); ++i) {
        const double dz = mesh.subdomainMax[i].z - mesh.subdomainMin[i].z;
        if (dz > 0.0 && std::isfinite(dz)) {
            minThickness = std::min(minThickness, dz);
        }
    }
    return std::isfinite(minThickness) ? minThickness : std::numeric_limits<double>::quiet_NaN();
}

static double rramMaxThermalDiffusionLength(const Mesh& mesh, const CaseConfig& physics)
{
    double maxAlpha = 0.0;
    for (const Tet& tet : mesh.tets) {
        const Material& material = materialForTet(physics, tet);
        const double alpha = rramScalarConductivity(material)
            / std::max(1.0e-300, material.density * material.heatCapacity);
        if (std::isfinite(alpha)) {
            maxAlpha = std::max(maxAlpha, alpha);
        }
    }
    return std::sqrt(std::max(0.0, maxAlpha) * std::max(0.0, physics.timeStep));
}

static void writeRramDiagnosticsBundle(const Mesh& mesh,
                                       const CaseConfig& physics,
                                       const SchwarzBlockSolver& solver,
                                       const std::vector<double>& temperature,
                                       const std::filesystem::path& outputDir)
{
    const std::vector<RramTetDiagnostic> tetDiagnostics =
        collectRramTetDiagnostics(mesh, physics, temperature);
    writeRramPartitionDiagnostics(mesh,
                                  physics,
                                  solver,
                                  tetDiagnostics,
                                  outputDir / "rram_partition_diagnostics.csv");
    const std::vector<RramInterfaceDiagnosticRow> interfaceRows =
        collectRramInterfaceDiagnostics(mesh, physics, tetDiagnostics);
    writeRramInterfaceDiagnostics(interfaceRows,
                                  outputDir / "interface_diagnostics.csv");
    writeRramInterfaceCutChecks(interfaceRows,
                                outputDir / "rram_interface_cut_checks.csv");
    writeRramPhysicalOverlapStats(mesh,
                                  solver,
                                  tetDiagnostics,
                                  outputDir / "rram_physical_overlap_stats.csv");
}

template <typename SolverRunResult>
static void writeRramSolverDiagnosticsSummary(const std::vector<SolverRunResult>& results,
                                              const Mesh& mesh,
                                              const CaseConfig& physics,
                                              const SchwarzBlockSolver& solver,
                                              const std::filesystem::path& path)
{
    const std::vector<RramTetDiagnostic> tetDiagnostics =
        collectRramTetDiagnostics(mesh, physics, {});
    const std::vector<RramInterfaceDiagnosticRow> interfaceRows =
        collectRramInterfaceDiagnostics(mesh, physics, tetDiagnostics);
    double maxKRatio = 0.0;
    bool cutsCriticalRegion = false;
    for (const RramInterfaceDiagnosticRow& row : interfaceRows) {
        if (std::isfinite(row.kRatio)) {
            maxKRatio = std::max(maxKRatio, row.kRatio);
        }
        cutsCriticalRegion = cutsCriticalRegion
            || row.cutsFilamentOrActiveOxide
            || row.cutsMaxHeatSourceRegion
            || row.cutsMaxTemperatureGradientRegion
            || row.cutsElectrodeOxideInterface;
    }

    std::string bestSolver = "not_run";
    int bestIterations = std::numeric_limits<int>::max();
    double bestTotalTime = std::numeric_limits<double>::infinity();
    bool anyFailure = false;
    for (const SolverRunResult& result : results) {
        anyFailure = anyFailure || result.stats.status != "success";
        const bool candidate = result.stats.status == "success"
            && std::isfinite(result.stats.finalRelativeResidual);
        if (candidate
            && (result.stats.totalIterations < bestIterations
                || (result.stats.totalIterations == bestIterations
                    && result.stats.setupSeconds + result.stats.solveSeconds < bestTotalTime))) {
            bestSolver = result.stats.name;
            bestIterations = result.stats.totalIterations;
            bestTotalTime = result.stats.setupSeconds + result.stats.solveSeconds;
        }
    }

    const RramPhysicalOverlapStats layer2Stats =
        computeRramPhysicalOverlapStats(mesh, solver, 2, tetDiagnostics);
    const double minLayerThickness = rramMinimumLayerThickness(mesh);
    const double diffusionLength = rramMaxThermalDiffusionLength(mesh, physics);
    const bool layer2TooThin =
        std::isfinite(layer2Stats.meanThickness)
        && ((std::isfinite(minLayerThickness) && layer2Stats.meanThickness < 0.5 * minLayerThickness)
            || (std::isfinite(diffusionLength) && layer2Stats.meanThickness < 0.5 * diffusionLength));

    const bool highKContrast = maxKRatio >= 10.0;
    const bool needsAmgOrCoarse = highKContrast || anyFailure || layer2TooThin;
    const bool recommendRepartition = cutsCriticalRegion || layer2TooThin || anyFailure;

    std::set<std::string> notes;
    notes.insert("RAS_FGMRES_remains_mainline");
    notes.insert("constant_two_level_coarse_correction_remains_opt_in");
    if (highKContrast) {
        notes.insert("high_k_contrast_interface_detected");
    }
    if (cutsCriticalRegion) {
        notes.insert("current_interfaces_touch_rram_critical_regions");
    }
    if (layer2TooThin) {
        notes.insert("layer2_overlap_is_physically_thin_consider_length_based_overlap");
    }
    if (anyFailure) {
        notes.insert("at_least_one_solver_failed_or_did_not_converge");
    }

    std::ofstream out(path);
    const double reportedBestTotalTime = std::isfinite(bestTotalTime)
        ? bestTotalTime
        : std::numeric_limits<double>::quiet_NaN();
    out << "partition_mode,best_partition,overlap_effective,has_high_k_contrast_interface,"
        << "needs_amg_or_coarse_correction,recommend_repartition,best_solver,best_iterations,"
        << "best_total_time,max_k_ratio,layer2_mean_overlap_thickness,min_layer_thickness,"
        << "thermal_diffusion_length,notes\n";
    out << std::setprecision(16)
        << csvEscape(solver.partitionMode()) << ','
        << csvEscape(physics.schwarz.partitionMode) << ','
        << csvEscape("single_run_needs_partition_sweep") << ','
        << (highKContrast ? 1 : 0) << ','
        << (needsAmgOrCoarse ? 1 : 0) << ','
        << (recommendRepartition ? 1 : 0) << ','
        << csvEscape(bestSolver) << ','
        << (bestIterations == std::numeric_limits<int>::max() ? -1 : bestIterations) << ','
        << reportedBestTotalTime << ','
        << maxKRatio << ','
        << layer2Stats.meanThickness << ','
        << minLayerThickness << ','
        << diffusionLength << ','
        << csvEscape(joinStrings(notes, "|")) << '\n';
}

struct ResidualNormPair {
    double absolute = std::numeric_limits<double>::quiet_NaN();
    double relative = std::numeric_limits<double>::quiet_NaN();
};

template <typename SolverRunResult>
static ResidualNormPair finalResidualNormPair(const SparseMatrix& system,
                                              const SolverRunResult* result)
{
    if (result == nullptr
        || result->temperature.empty()
        || result->finalRhs.empty()
        || result->temperature.size() != result->finalRhs.size()
        || result->temperature.size() != static_cast<size_t>(system.size())) {
        return {};
    }

    std::vector<double> residual = trueResidualVector(system, result->temperature, result->finalRhs);
    const double absolute = l2Norm(residual);
    return {absolute, absolute / std::max(1.0e-300, l2Norm(result->finalRhs))};
}

template <typename SolverRunResult>
static void writeMonolithicVsSchwarzComparison(const std::vector<SolverRunResult>& results,
                                               const CaseConfig& physics,
                                               bool steady,
                                               const SparseMatrix& system,
                                               const std::filesystem::path& path)
{
    const SolverRunResult* monolithic = nullptr;
    const SolverRunResult* schwarzPrecondFgmres = nullptr;
    bool hasSchwarzFamily = false;
    for (const SolverRunResult& result : results) {
        if (result.stats.name.find("Schwarz-SIPG") != std::string::npos) {
            hasSchwarzFamily = true;
        }
        if (result.stats.name.find("Schwarz-Precond-FGMRES") != std::string::npos) {
            hasSchwarzFamily = true;
            schwarzPrecondFgmres = &result;
        }
        if (result.stats.name == "Global-PARDISO-General-Direct") {
            monolithic = &result;
        }
    }
    if (monolithic == nullptr) {
        for (const SolverRunResult& result : results) {
            if (result.stats.name.find("Global-PARDISO") != std::string::npos) {
                monolithic = &result;
                break;
            }
        }
    }
    if (!hasSchwarzFamily
        && physics.solverMethod != "schwarz"
        && physics.solverMethod != "schwarz_precond_fgmres"
        && physics.solverMethod != "schwarz_precond_fgmres_two_level") {
        return;
    }

    std::ofstream out(path);
    out << "solver,case,time_step,time,rel_L2_error,max_abs_error,Tmax_mono,Tmax_solver,"
        << "Tmax_diff,total_time,iteration_count,status,direct_abs_residual,direct_rel_residual,"
        << "fgmres_abs_residual,fgmres_rel_residual\n";
    out << std::setprecision(16);

    const ResidualNormPair directResidual = finalResidualNormPair(system, monolithic);
    const ResidualNormPair fgmresResidual = finalResidualNormPair(system, schwarzPrecondFgmres);

    const int outputStep = steady ? 0 : physics.timeSteps;
    const double outputTime = steady ? physics.startTime : physics.timeStep * physics.timeSteps;
    if (monolithic == nullptr || monolithic->temperature.empty()) {
        out << "missing_monolithic_reference," << csvEscape(physics.name) << ','
            << outputStep << ','
            << outputTime
            << ",nan,nan,nan,nan,nan,nan,nan,failed,nan,nan,nan,nan\n";
        return;
    }

    const auto monoMinMax = std::minmax_element(monolithic->temperature.begin(), monolithic->temperature.end());
    for (const SolverRunResult& result : results) {
        const double totalTime = result.stats.setupSeconds + result.stats.solveSeconds;
        if (result.temperature.empty()
            || result.temperature.size() != monolithic->temperature.size()) {
            out << csvEscape(result.stats.name) << ','
                << csvEscape(physics.name) << ','
                << outputStep << ','
                << outputTime
                << ",nan,nan," << *monoMinMax.second
                << ",nan,nan," << totalTime << ','
                << result.stats.totalIterations << ','
                << csvEscape(result.stats.status) << ','
                << directResidual.absolute << ','
                << directResidual.relative << ','
                << fgmresResidual.absolute << ','
                << fgmresResidual.relative << '\n';
            continue;
        }
        const auto solverMinMax = std::minmax_element(result.temperature.begin(), result.temperature.end());
        out << csvEscape(result.stats.name) << ','
            << csvEscape(physics.name) << ','
            << outputStep << ','
            << outputTime << ','
            << relativeL2Difference(result.temperature, monolithic->temperature) << ','
            << maxAbsDifference(result.temperature, monolithic->temperature) << ','
            << *monoMinMax.second << ','
            << *solverMinMax.second << ','
            << (*solverMinMax.second - *monoMinMax.second) << ','
            << totalTime << ','
            << result.stats.totalIterations << ','
            << csvEscape(result.stats.status) << ','
            << directResidual.absolute << ','
            << directResidual.relative << ','
            << fgmresResidual.absolute << ','
            << fgmresResidual.relative << '\n';
    }
}

template <typename SolverRunResult>
static void writeOverlapDiagnostics(const std::vector<SolverRunResult>& results,
                                    const std::vector<double>& maxDiffs,
                                    const std::vector<double>& relativeL2Diffs,
                                    const CaseConfig& physics,
                                    const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "mode,overlap_layers,partition_mode,iterations,final_rel_residual,rel_L2_error,max_abs_error,"
        << "total_time,preconditioner_apply_time\n";
    out << std::setprecision(16);
    for (size_t i = 0; i < results.size(); ++i) {
        const SolverRunResult& result = results[i];
        const bool isSchwarz = result.stats.name.find("Schwarz") != std::string::npos
            || result.stats.name.find("RAS") != std::string::npos;
        if (!isSchwarz && result.stats.name.find("Global-PARDISO") == std::string::npos) {
            continue;
        }
        out << csvEscape(result.stats.name) << ','
            << physics.schwarz.overlapLayers << ','
            << csvEscape(physics.schwarz.partitionMode) << ','
            << result.stats.totalIterations << ','
            << result.stats.finalRelativeResidual << ','
            << (i < relativeL2Diffs.size() ? relativeL2Diffs[i] : std::numeric_limits<double>::quiet_NaN()) << ','
            << (i < maxDiffs.size() ? maxDiffs[i] : std::numeric_limits<double>::quiet_NaN()) << ','
            << (result.stats.setupSeconds + result.stats.solveSeconds) << ','
            << result.stats.preconditionerApplySeconds << '\n';
    }
}

template <typename SolverRunResult>
static void writeLargeScaleSolverSummary(const std::vector<SolverRunResult>& results,
                                         const std::vector<double>& maxDiffs,
                                         const std::vector<double>& relativeL2Diffs,
                                         const Mesh& mesh,
                                         const CaseConfig& physics,
                                         bool steady,
                                         const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "solver,status,failure_reason,case,time_step,time,overlap_mode,overlap_layers,partition_mode,"
        << "solver_method,preconditioner,ras_type,local_solver,initial_guess_type,"
        << "iterations,max_step_iterations,final_rel_residual,total_time,setup_time,solve_time,"
        << "true_relative_residual,rhs_norm,rhs_min,rhs_max,rhs_nonzero_count,rhs_near_zero_count,"
        << "heat_source_total,dirichlet_rhs_norm,nitsche_rhs_norm,robin_rhs_norm,neumann_rhs_norm,"
        << "initial_x_norm,initial_ax_norm,initial_residual_norm,initial_relative_residual,"
        << "solver_tolerance,zero_iteration_due_to_initial_residual,"
        << "Tmin,Tmax,Tavg,interface_avg_jump,rel_L2_error,max_abs_error,"
        << "preconditioner_apply_time,preconditioner_apply_calls,preconditioner_apply_time_per_call,"
        << "preconditioner_memory_mb,working_set_before_mb,working_set_after_mb,peak_working_set_mb,"
        << "local_diag_scaling,diag_scaling_eps,local_ic_shift_mode,local_ic_shift_used_max,"
        << "ras_halo_build_time,ras_local_matrix_assembly_time,"
        << "ras_local_factorization_time,ras_local_solve_apply_time,ras_restriction_time,"
        << "ras_communication_or_halo_update_time,coarse_dim,coarse_setup_time,"
        << "coarse_solve_time,coarse_matrix_nnz,coarse_residual_norm,"
        << "coarse_enabled,coarse_space,coarse_correction,coarse_rhs_norm,coarse_solution_norm,"
        << "coarse_correction_norm,local_correction_norm,coarse_to_local_norm_ratio,"
        << "deflation_enabled,deflation_dim,deflation_setup_time,deflation_apply_time,"
        << "deflation_correction_norm,matvec_time\n";
    out << std::setprecision(16);

    const int outputStep = steady ? 0 : physics.timeSteps;
    const double outputTime = steady ? physics.startTime : physics.timeStep * physics.timeSteps;
    for (size_t i = 0; i < results.size(); ++i) {
        const SolverRunResult& result = results[i];
        const double totalTime = result.stats.setupSeconds + result.stats.solveSeconds;
        const double applyTimePerCall = result.stats.preconditionerApplyCalls > 0
            ? result.stats.preconditionerApplySeconds / static_cast<double>(result.stats.preconditionerApplyCalls)
            : 0.0;
        const double jump = result.temperature.empty()
            ? std::numeric_limits<double>::quiet_NaN()
            : interfaceAverageJump(mesh, result.temperature);
        const std::string solverMethod = result.stats.solverMethod.empty()
            ? result.stats.name
            : result.stats.solverMethod;
        const std::string preconditioner = result.stats.preconditioner.empty()
            ? result.stats.name
            : result.stats.preconditioner;
        const std::string rasType = result.stats.rasType.empty()
            ? (physics.schwarz.overlapMode == "halo" ? "additive" : physics.schwarz.overlapMode)
            : result.stats.rasType;
        std::string localSolver = result.stats.localSolver;
        if (localSolver.empty()) {
            if (result.stats.name.find("ILUT") != std::string::npos) {
                localSolver = "ilut";
            } else if (result.stats.name.find("IC") != std::string::npos) {
                localSolver = "ic";
            } else if (result.stats.name.find("PARDISO") != std::string::npos
                       || result.stats.name.find("Schwarz") != std::string::npos) {
                localSolver = "direct";
            }
        }
        const bool coarseEnabled = result.stats.coarseEnabled || result.stats.coarseDim > 0;
        const bool deflationEnabled = result.stats.deflationEnabled
            || result.stats.name.find("Deflated") != std::string::npos
            || result.stats.name.find("Eig") != std::string::npos;
        out << csvEscape(result.stats.name) << ','
            << csvEscape(result.stats.status) << ','
            << csvEscape(result.stats.failureReason) << ','
            << csvEscape(physics.name) << ','
            << outputStep << ','
            << outputTime << ','
            << csvEscape(physics.schwarz.overlapMode) << ','
            << physics.schwarz.overlapLayers << ','
            << csvEscape(physics.schwarz.partitionMode) << ','
            << csvEscape(solverMethod) << ','
            << csvEscape(preconditioner) << ','
            << csvEscape(rasType) << ','
            << csvEscape(localSolver) << ','
            << csvEscape(result.stats.initialGuessType) << ','
            << result.stats.totalIterations << ','
            << result.stats.maxIterations << ','
            << result.stats.finalRelativeResidual << ','
            << totalTime << ','
            << result.stats.setupSeconds << ','
            << result.stats.solveSeconds << ','
            << result.stats.trueRelativeResidual << ','
            << result.stats.rhsNorm << ','
            << result.stats.rhsMin << ','
            << result.stats.rhsMax << ','
            << result.stats.rhsNonzeroCount << ','
            << result.stats.rhsNearZeroCount << ','
            << result.stats.heatSourceTotal << ','
            << result.stats.dirichletRhsContributionNorm << ','
            << result.stats.nitscheRhsContributionNorm << ','
            << result.stats.robinRhsContributionNorm << ','
            << result.stats.neumannRhsContributionNorm << ','
            << result.stats.initialXNorm << ','
            << result.stats.initialAxNorm << ','
            << result.stats.initialResidualNorm << ','
            << result.stats.initialRelativeResidual << ','
            << result.stats.solverTolerance << ','
            << (result.stats.zeroIterationDueToInitialResidual ? 1 : 0) << ','
            << result.stats.temperatureMin << ','
            << result.stats.temperatureMax << ','
            << result.stats.temperatureAverage << ','
            << jump << ','
            << (i < relativeL2Diffs.size() ? relativeL2Diffs[i] : std::numeric_limits<double>::quiet_NaN()) << ','
            << (i < maxDiffs.size() ? maxDiffs[i] : std::numeric_limits<double>::quiet_NaN()) << ','
            << result.stats.preconditionerApplySeconds << ','
            << result.stats.preconditionerApplyCalls << ','
            << applyTimePerCall << ','
            << megabytes(result.stats.preconditionerBytes) << ','
            << megabytes(result.stats.workingSetBeforeBytes) << ','
            << megabytes(result.stats.workingSetAfterBytes) << ','
            << megabytes(result.stats.peakWorkingSetBytes) << ','
            << (result.stats.localDiagScaling ? 1 : 0) << ','
            << result.stats.diagScalingEps << ','
            << csvEscape(result.stats.localIcShiftMode) << ','
            << result.stats.localIcShiftUsedMax << ','
            << result.stats.rasHaloBuildSeconds << ','
            << result.stats.rasLocalMatrixAssemblySeconds << ','
            << result.stats.rasLocalFactorizationSeconds << ','
            << result.stats.rasLocalSolveApplySeconds << ','
            << result.stats.rasRestrictionSeconds << ','
            << result.stats.rasCommunicationOrHaloUpdateSeconds << ','
            << result.stats.coarseDim << ','
            << result.stats.coarseSetupSeconds << ','
            << result.stats.coarseSolveSeconds << ','
            << result.stats.coarseMatrixNnz << ','
            << result.stats.coarseResidualNorm << ','
            << (coarseEnabled ? 1 : 0) << ','
            << csvEscape(result.stats.coarseSpace) << ','
            << csvEscape(result.stats.coarseCorrection) << ','
            << result.stats.coarseRhsNorm << ','
            << result.stats.coarseSolutionNorm << ','
            << result.stats.coarseCorrectionNorm << ','
            << result.stats.localCorrectionNorm << ','
            << result.stats.coarseToLocalNormRatio << ','
            << (deflationEnabled ? 1 : 0) << ','
            << result.stats.deflationDim << ','
            << result.stats.deflationSetupSeconds << ','
            << result.stats.deflationApplySeconds << ','
            << result.stats.deflationCorrectionNorm << ','
            << result.stats.matvecSeconds << '\n';
    }
}
