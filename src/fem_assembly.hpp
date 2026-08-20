#pragma once

// Finite-element geometry, quadrature, volume/boundary/interface assembly, and interface diagnostics.
// This file is intentionally included from main.cpp after the preceding SIPG modules.

static ElementGeometry elementGeometry(const Mesh& mesh, const Tet& tet)
{
    const Vec3 p0 = mesh.nodes[static_cast<size_t>(tet.v[0])].p;
    const Vec3 p1 = mesh.nodes[static_cast<size_t>(tet.v[1])].p;
    const Vec3 p2 = mesh.nodes[static_cast<size_t>(tet.v[2])].p;
    const Vec3 p3 = mesh.nodes[static_cast<size_t>(tet.v[3])].p;
    const Vec3 a = p1 - p0;
    const Vec3 b = p2 - p0;
    const Vec3 c = p3 - p0;
    const double detJ = dot(a, cross(b, c));
    if (std::abs(detJ) < 1.0e-30) {
        throw std::runtime_error("Degenerate tetrahedron detected.");
    }

    const std::array<Vec3, 4> gradRef{{{-1.0, -1.0, -1.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    const Vec3 r0 = cross(b, c) * (1.0 / detJ);
    const Vec3 r1 = cross(c, a) * (1.0 / detJ);
    const Vec3 r2 = cross(a, b) * (1.0 / detJ);

    ElementGeometry geo;
    geo.detJ = std::abs(detJ);
    for (int i = 0; i < 4; ++i) {
        const Vec3 g = gradRef[static_cast<size_t>(i)];
        geo.gradLambda[static_cast<size_t>(i)] = {r0.x * g.x + r1.x * g.y + r2.x * g.z,
                                                  r0.y * g.x + r1.y * g.y + r2.y * g.z,
                                                  r0.z * g.x + r1.z * g.y + r2.z * g.z};
    }
    return geo;
}

static double tetVolume(const Mesh& mesh, const Tet& tet)
{
    return elementGeometry(mesh, tet).detJ / 6.0;
}

static std::array<double, 10> shapeP2(const std::array<double, 4>& l)
{
    return {{
        l[0] * (2.0 * l[0] - 1.0), l[1] * (2.0 * l[1] - 1.0),
        l[2] * (2.0 * l[2] - 1.0), l[3] * (2.0 * l[3] - 1.0),
        4.0 * l[0] * l[1], 4.0 * l[0] * l[2], 4.0 * l[0] * l[3],
        4.0 * l[1] * l[2], 4.0 * l[1] * l[3], 4.0 * l[2] * l[3]
    }};
}

static std::array<Vec3, 10> gradShapeP2(const std::array<double, 4>& l, const ElementGeometry& geo)
{
    std::array<Vec3, 10> g{};
    for (int i = 0; i < 4; ++i) {
        g[static_cast<size_t>(i)] = (4.0 * l[static_cast<size_t>(i)] - 1.0) * geo.gradLambda[static_cast<size_t>(i)];
    }
    const std::array<std::array<int, 2>, 6> edges{{
        {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}, {{1, 3}}, {{2, 3}}
    }};
    for (int e = 0; e < 6; ++e) {
        const int a = edges[static_cast<size_t>(e)][0];
        const int b = edges[static_cast<size_t>(e)][1];
        g[static_cast<size_t>(4 + e)] =
            4.0 * (l[static_cast<size_t>(a)] * geo.gradLambda[static_cast<size_t>(b)]
                 + l[static_cast<size_t>(b)] * geo.gradLambda[static_cast<size_t>(a)]);
    }
    return g;
}

struct TetQuadraturePoint {
    double weight = 0.0;
    std::array<double, 10> shape{};
    std::array<std::array<double, 4>, 10> gradCoeff{};
};

struct TriangleQuadraturePoint {
    double a = 0.0;
    double b = 0.0;
    double weight = 0.0;
};

static std::array<std::array<double, 4>, 10> gradShapeP2Coefficients(const std::array<double, 4>& l)
{
    std::array<std::array<double, 4>, 10> coeff{};
    for (int i = 0; i < 4; ++i) {
        coeff[static_cast<size_t>(i)][static_cast<size_t>(i)] = 4.0 * l[static_cast<size_t>(i)] - 1.0;
    }
    const std::array<std::array<int, 2>, 6> edges{{
        {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}, {{1, 3}}, {{2, 3}}
    }};
    for (int e = 0; e < 6; ++e) {
        const int a = edges[static_cast<size_t>(e)][0];
        const int b = edges[static_cast<size_t>(e)][1];
        coeff[static_cast<size_t>(4 + e)][static_cast<size_t>(a)] = 4.0 * l[static_cast<size_t>(b)];
        coeff[static_cast<size_t>(4 + e)][static_cast<size_t>(b)] = 4.0 * l[static_cast<size_t>(a)];
    }
    return coeff;
}

static std::array<Vec3, 10> physicalGradP2(const TetQuadraturePoint& qp, const ElementGeometry& geo)
{
    std::array<Vec3, 10> grad{};
    for (int a = 0; a < 10; ++a) {
        Vec3 value{};
        for (int i = 0; i < 4; ++i) {
            value = value + qp.gradCoeff[static_cast<size_t>(a)][static_cast<size_t>(i)]
                          * geo.gradLambda[static_cast<size_t>(i)];
        }
        grad[static_cast<size_t>(a)] = value;
    }
    return grad;
}

static std::vector<TetQuadraturePoint> makeTetQuadrature()
{
    const std::array<double, 4> gp{{0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262}};
    const std::array<double, 4> gw{{0.17392742256872692, 0.32607257743127307, 0.32607257743127307, 0.17392742256872692}};
    std::vector<TetQuadraturePoint> points;
    points.reserve(64);
    for (int iu = 0; iu < 4; ++iu) {
        for (int iv = 0; iv < 4; ++iv) {
            for (int iw = 0; iw < 4; ++iw) {
                const double u = gp[static_cast<size_t>(iu)];
                const double v = gp[static_cast<size_t>(iv)];
                const double w = gp[static_cast<size_t>(iw)];
                const double r = u;
                const double s = (1.0 - u) * v;
                const double t = (1.0 - u) * (1.0 - v) * w;
                const std::array<double, 4> lambda{{1.0 - r - s - t, r, s, t}};
                TetQuadraturePoint qp;
                qp.weight = gw[static_cast<size_t>(iu)] * gw[static_cast<size_t>(iv)] * gw[static_cast<size_t>(iw)]
                    * (1.0 - u) * (1.0 - u) * (1.0 - v);
                qp.shape = shapeP2(lambda);
                qp.gradCoeff = gradShapeP2Coefficients(lambda);
                points.push_back(qp);
            }
        }
    }
    return points;
}

static std::vector<TriangleQuadraturePoint> makeTriangleQuadrature()
{
    const std::array<double, 4> gp{{0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262}};
    const std::array<double, 4> gw{{0.17392742256872692, 0.32607257743127307, 0.32607257743127307, 0.17392742256872692}};
    std::vector<TriangleQuadraturePoint> points;
    points.reserve(16);
    for (int iu = 0; iu < 4; ++iu) {
        for (int iv = 0; iv < 4; ++iv) {
            const double u = gp[static_cast<size_t>(iu)];
            const double v = gp[static_cast<size_t>(iv)];
            TriangleQuadraturePoint qp;
            qp.a = u;
            qp.b = (1.0 - u) * v;
            qp.weight = gw[static_cast<size_t>(iu)] * gw[static_cast<size_t>(iv)] * (1.0 - u);
            points.push_back(qp);
        }
    }
    return points;
}

static std::array<double, 3> barycentricOnTriangle3D(const Vec3& p, const std::array<Vec3, 3>& tri)
{
    const Vec3 v0 = tri[1] - tri[0];
    const Vec3 v1 = tri[2] - tri[0];
    const Vec3 v2 = p - tri[0];
    const double d00 = dot(v0, v0);
    const double d01 = dot(v0, v1);
    const double d11 = dot(v1, v1);
    const double d20 = dot(v2, v0);
    const double d21 = dot(v2, v1);
    const double denom = d00 * d11 - d01 * d01;
    const double scale = std::max(d00 * d11, std::numeric_limits<double>::min());
    const double tolerance = std::numeric_limits<double>::epsilon() * scale * 64.0;
    if (std::abs(denom) < tolerance) {
        throw std::runtime_error("Degenerate boundary triangle detected.");
    }
    const double l1 = (d11 * d20 - d01 * d21) / denom;
    const double l2 = (d00 * d21 - d01 * d20) / denom;
    return {{1.0 - l1 - l2, l1, l2}};
}

static std::array<double, 4> lambdaOnTetFace(const Vec3& p,
                                             const Tet& tet,
                                             const std::array<int, 3>& locals,
                                             const Mesh& mesh)
{
    std::array<Vec3, 3> tri{{
        mesh.nodes[static_cast<size_t>(tet.v[static_cast<size_t>(locals[0])])].p,
        mesh.nodes[static_cast<size_t>(tet.v[static_cast<size_t>(locals[1])])].p,
        mesh.nodes[static_cast<size_t>(tet.v[static_cast<size_t>(locals[2])])].p
    }};
    const std::array<double, 3> faceLambda = barycentricOnTriangle3D(p, tri);
    std::array<double, 4> lambda{{0.0, 0.0, 0.0, 0.0}};
    for (int i = 0; i < 3; ++i) {
        lambda[static_cast<size_t>(locals[static_cast<size_t>(i)])] = faceLambda[static_cast<size_t>(i)];
    }
    return lambda;
}

static bool barycentricInTet(const Mesh& mesh, const Tet& tet, const ElementGeometry& geo, const Vec3& p, std::array<double, 4>& lambda)
{
    const Vec3 p0 = mesh.nodes[static_cast<size_t>(tet.v[0])].p;
    lambda[1] = dot(geo.gradLambda[1], p - p0);
    lambda[2] = dot(geo.gradLambda[2], p - p0);
    lambda[3] = dot(geo.gradLambda[3], p - p0);
    lambda[0] = 1.0 - lambda[1] - lambda[2] - lambda[3];
    return std::all_of(lambda.begin(), lambda.end(), [](double value) {
        return value >= -1.0e-9 && value <= 1.0 + 1.0e-9;
    });
}

static double selectedSourceVolume(const Mesh& mesh, const HeatSource& source)
{
    double volume = 0.0;
    for (const Tet& tet : mesh.tets) {
        if (tetMatchesHeatSource(tet, source)) {
            volume += tetVolume(mesh, tet);
        }
    }
    return volume;
}

static std::vector<HeatSourceAssemblyDiagnostic> collectHeatSourceDiagnostics(const Mesh& mesh,
                                                                              const CaseConfig& config)
{
    const std::vector<TetQuadraturePoint> quadrature = makeTetQuadrature();
    std::vector<HeatSourceAssemblyDiagnostic> diagnostics(config.heatSources.size());
    for (size_t i = 0; i < config.heatSources.size(); ++i) {
        HeatSourceAssemblyDiagnostic& diag = diagnostics[i];
        const HeatSource& source = config.heatSources[i];
        diag.index = static_cast<int>(i);
        diag.subdomain = source.subdomain;
        diag.domainEntity = source.domainEntity;
        diag.configuredPowerW = source.heatRateW;
        diag.sourceVolumeUsedForDensity = selectedSourceVolume(mesh, source);
        if (diag.sourceVolumeUsedForDensity > 0.0) {
            diag.volumetricQ = source.heatRateW / diag.sourceVolumeUsedForDensity;
        }
    }

    for (const Tet& tet : mesh.tets) {
        const ElementGeometry geo = elementGeometry(mesh, tet);
        const double volume = tetVolume(mesh, tet);
        for (size_t i = 0; i < config.heatSources.size(); ++i) {
            HeatSourceAssemblyDiagnostic& diag = diagnostics[i];
            if (!tetMatchesHeatSource(tet, config.heatSources[i])) {
                continue;
            }
            ++diag.tetCount;
            diag.sourceVolumeAssembled += volume;
            for (const TetQuadraturePoint& qp : quadrature) {
                double shapeSum = 0.0;
                for (double value : qp.shape) {
                    shapeSum += value;
                }
                diag.quadraturePowerW += diag.volumetricQ * shapeSum * qp.weight * geo.detJ;
            }
        }
    }

    for (HeatSourceAssemblyDiagnostic& diag : diagnostics) {
        diag.expectedPowerW = diag.volumetricQ * diag.sourceVolumeAssembled;
    }
    return diagnostics;
}

static void mergeSourceEntries(std::vector<std::vector<VectorEntry>>& localSources, std::vector<double>& source)
{
    std::vector<VectorEntry> entries;
    size_t count = 0;
    for (const auto& local : localSources) {
        count += local.size();
    }
    entries.reserve(count);
    for (auto& local : localSources) {
        entries.insert(entries.end(),
                       std::make_move_iterator(local.begin()),
                       std::make_move_iterator(local.end()));
    }
    std::sort(entries.begin(), entries.end(), [](const VectorEntry& a, const VectorEntry& b) {
        return a.index < b.index;
    });

    size_t i = 0;
    while (i < entries.size()) {
        const int index = entries[i].index;
        double sum = 0.0;
        while (i < entries.size() && entries[i].index == index) {
            sum += entries[i].value;
            ++i;
        }
        if (index >= 0 && index < static_cast<int>(source.size())) {
            source[static_cast<size_t>(index)] += sum;
        }
    }
}

static void assembleVolume(const Mesh& mesh,
                           const CaseConfig& config,
                           SparseMatrix* mass,
                           SparseMatrix& stiffness,
                           std::vector<double>& source)
{
    const std::vector<TetQuadraturePoint> quadrature = makeTetQuadrature();
    std::vector<double> volumetricHeat(config.heatSources.size(), 0.0);
    for (size_t i = 0; i < config.heatSources.size(); ++i) {
        const double sourceVolume = selectedSourceVolume(mesh, config.heatSources[i]);
        if (sourceVolume > 0.0) {
            volumetricHeat[i] = config.heatSources[i].heatRateW / sourceVolume;
        }
    }

    const bool assembleMass = mass != nullptr;
    const unsigned int requestedThreads = solverParallelWorkers();
    const int threadCount = static_cast<int>(std::max(1u, std::min<unsigned int>(
        requestedThreads, static_cast<unsigned int>(std::max<size_t>(1, mesh.tets.size())))));
    std::vector<std::vector<MatrixEntry>> localStiffness(static_cast<size_t>(threadCount));
    std::vector<std::vector<MatrixEntry>> localMass(static_cast<size_t>(assembleMass ? threadCount : 0));
    std::vector<std::vector<VectorEntry>> localSources(static_cast<size_t>(threadCount));
    const size_t tetsPerThread = (mesh.tets.size() + static_cast<size_t>(threadCount) - 1) / static_cast<size_t>(threadCount);
    for (int t = 0; t < threadCount; ++t) {
        localStiffness[static_cast<size_t>(t)].reserve(tetsPerThread * 100);
        if (assembleMass) {
            localMass[static_cast<size_t>(t)].reserve(tetsPerThread * 100);
        }
        if (!config.heatSources.empty()) {
            localSources[static_cast<size_t>(t)].reserve(tetsPerThread * 10);
        }
    }

    const auto assembleTet = [&](size_t tetIndex,
                                 std::vector<MatrixEntry>& massEntries,
                                 std::vector<MatrixEntry>& stiffnessEntries,
                                 std::vector<VectorEntry>& sourceEntries) {
        const Tet& tet = mesh.tets[tetIndex];
        const ElementGeometry geo = elementGeometry(mesh, tet);
        const Material& mat = materialForTet(config, tet);
        double localVolumetricHeat = 0.0;
        for (size_t i = 0; i < config.heatSources.size(); ++i) {
            if (tetMatchesHeatSource(tet, config.heatSources[i])) {
                localVolumetricHeat += volumetricHeat[i];
            }
        }

        std::array<double, 100> localStiffnessMatrix{};
        std::array<double, 100> localMassMatrix{};
        std::array<double, 10> localSource{};
        const double heatCapacity = mat.density * mat.heatCapacity;
        for (const TetQuadraturePoint& qp : quadrature) {
            const double weight = qp.weight * geo.detJ;
            const auto grad = physicalGradP2(qp, geo);
            for (int a = 0; a < 10; ++a) {
                if (std::abs(localVolumetricHeat) > 0.0) {
                    localSource[static_cast<size_t>(a)] += localVolumetricHeat * qp.shape[static_cast<size_t>(a)] * weight;
                }
                for (int b = 0; b < 10; ++b) {
                    const size_t offset = static_cast<size_t>(a * 10 + b);
                    localStiffnessMatrix[offset] +=
                        conductivityGradientDot(mat, grad[static_cast<size_t>(a)], grad[static_cast<size_t>(b)]) * weight;
                    if (assembleMass) {
                        localMassMatrix[offset] += heatCapacity
                            * qp.shape[static_cast<size_t>(a)] * qp.shape[static_cast<size_t>(b)] * weight;
                    }
                }
            }
        }

        for (int a = 0; a < 10; ++a) {
            const int row = tet.dof[static_cast<size_t>(a)];
            if (std::abs(localSource[static_cast<size_t>(a)]) > 0.0) {
                sourceEntries.push_back({row, localSource[static_cast<size_t>(a)]});
            }
            for (int b = 0; b < 10; ++b) {
                const int col = tet.dof[static_cast<size_t>(b)];
                const size_t offset = static_cast<size_t>(a * 10 + b);
                if (std::abs(localStiffnessMatrix[offset]) > 0.0) {
                    stiffnessEntries.push_back({row, col, localStiffnessMatrix[offset]});
                }
                if (assembleMass && std::abs(localMassMatrix[offset]) > 0.0) {
                    massEntries.push_back({row, col, localMassMatrix[offset]});
                }
            }
        }
    };

#ifdef _OPENMP
#pragma omp parallel num_threads(threadCount)
    {
        const int tid = omp_get_thread_num();
        auto& stiffnessEntries = localStiffness[static_cast<size_t>(tid)];
        auto& sourceEntries = localSources[static_cast<size_t>(tid)];
        std::vector<MatrixEntry> unusedMassEntries;
        auto& massEntries = assembleMass ? localMass[static_cast<size_t>(tid)] : unusedMassEntries;
#pragma omp for schedule(static)
        for (long long i = 0; i < static_cast<long long>(mesh.tets.size()); ++i) {
            assembleTet(static_cast<size_t>(i), massEntries, stiffnessEntries, sourceEntries);
        }
    }
#else
    for (size_t i = 0; i < mesh.tets.size(); ++i) {
        std::vector<MatrixEntry> unusedMassEntries;
        assembleTet(i,
                    assembleMass ? localMass[0] : unusedMassEntries,
                    localStiffness[0],
                    localSources[0]);
    }
#endif

    for (auto& entries : localStiffness) {
        stiffness.appendEntries(entries);
    }
    if (assembleMass) {
        for (auto& entries : localMass) {
            mass->appendEntries(entries);
        }
    }
    mergeSourceEntries(localSources, source);
}

static int assembleConvectionBoundaries(const Mesh& mesh,
                                        const CaseConfig& config,
                                        SparseMatrix& stiffness,
                                        std::vector<double>& source,
                                        bool disableLhs = false,
                                        AssemblyDiagnostics* diagnostics = nullptr)
{
    if (config.convectionConditions.empty()) {
        return 0;
    }

    const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
    std::vector<MatrixEntry> stiffnessEntries;
    std::vector<VectorEntry> sourceEntries;
    int matchedFaces = 0;

    for (const BoundaryFace& face : mesh.boundaryFaces) {
        for (const ConvectionCondition& condition : config.convectionConditions) {
            if (!boundaryMatches(face, condition)) {
                continue;
            }
            ++matchedFaces;
            const Tet& tet = mesh.tets[static_cast<size_t>(face.tet)];
            std::array<double, 100> localStiffness{};
            std::array<double, 10> localSource{};

            const double jacFace = norm(cross(face.points[1] - face.points[0], face.points[2] - face.points[0]));
            for (const TriangleQuadraturePoint& qp : quadrature) {
                const double weight = qp.weight * jacFace;
                const Vec3 q = (1.0 - qp.a - qp.b) * face.points[0]
                             + qp.a * face.points[1]
                             + qp.b * face.points[2];
                const auto lambda = lambdaOnTetFace(q, tet, face.local, mesh);
                const auto shape = shapeP2(lambda);
                for (int a = 0; a < 10; ++a) {
                    localSource[static_cast<size_t>(a)] +=
                        condition.coefficient * condition.ambientTemperature
                        * shape[static_cast<size_t>(a)] * weight;
                    for (int b = 0; b < 10; ++b) {
                        localStiffness[static_cast<size_t>(a * 10 + b)] +=
                            condition.coefficient
                            * shape[static_cast<size_t>(a)]
                            * shape[static_cast<size_t>(b)]
                            * weight;
                    }
                }
            }

            for (int a = 0; a < 10; ++a) {
                const int row = tet.dof[static_cast<size_t>(a)];
                if (std::abs(localSource[static_cast<size_t>(a)]) > 0.0) {
                    sourceEntries.push_back({row, localSource[static_cast<size_t>(a)]});
                }
                if (!disableLhs) {
                    for (int b = 0; b < 10; ++b) {
                        const double value = localStiffness[static_cast<size_t>(a * 10 + b)];
                        if (std::abs(value) > 0.0) {
                            stiffnessEntries.push_back({row, tet.dof[static_cast<size_t>(b)], value});
                        }
                        if (diagnostics != nullptr) {
                            addDiagonalEntry(diagnostics->robinDiag,
                                             row,
                                             tet.dof[static_cast<size_t>(b)],
                                             value);
                        }
                    }
                }
            }
        }
    }

    stiffness.appendEntries(stiffnessEntries);
    std::vector<std::vector<VectorEntry>> localSources;
    localSources.push_back(std::move(sourceEntries));
    mergeSourceEntries(localSources, source);
    return matchedFaces;
}

static int assembleHeatFluxBoundaries(const Mesh& mesh,
                                      const CaseConfig& config,
                                      std::vector<double>& source)
{
    if (config.heatFluxConditions.empty()) return 0;
    const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
    std::vector<VectorEntry> entries;
    int matchedFaces = 0;
    for (const BoundaryFace& face : mesh.boundaryFaces) {
        for (const HeatFluxCondition& condition : config.heatFluxConditions) {
            if (!boundaryMatches(face, condition)) continue;
            ++matchedFaces;
            const Tet& tet = mesh.tets[static_cast<size_t>(face.tet)];
            const double jacFace = norm(cross(face.points[1] - face.points[0],
                                              face.points[2] - face.points[0]));
            for (const TriangleQuadraturePoint& qp : quadrature) {
                const double weight = qp.weight * jacFace;
                const Vec3 q = (1.0 - qp.a - qp.b) * face.points[0] + qp.a * face.points[1] + qp.b * face.points[2];
                const auto shape = shapeP2(lambdaOnTetFace(q, tet, face.local, mesh));
                for (int a = 0; a < 10; ++a) {
                    entries.push_back({tet.dof[static_cast<size_t>(a)],
                                       condition.inwardFluxWm2 * shape[static_cast<size_t>(a)] * weight});
                }
            }
        }
    }
    std::vector<std::vector<VectorEntry>> grouped;
    grouped.push_back(std::move(entries));
    mergeSourceEntries(grouped, source);
    return matchedFaces;
}

static int assembleNitscheDirichletBoundaries(const Mesh& mesh,
                                              const CaseConfig& config,
                                              SparseMatrix& system,
                                              std::vector<double>& source,
                                              AssemblyDiagnostics* diagnostics = nullptr)
{
    if (config.dirichletConditions.empty()) {
        return 0;
    }

    const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
    std::vector<MatrixEntry> stiffnessEntries;
    std::vector<VectorEntry> sourceEntries;
    int matchedFaces = 0;

    for (const BoundaryFace& face : mesh.boundaryFaces) {
        for (const BoundaryCondition& condition : config.dirichletConditions) {
            if (!boundaryMatches(face, condition)) {
                continue;
            }
            ++matchedFaces;
            const Tet& tet = mesh.tets[static_cast<size_t>(face.tet)];
            const ElementGeometry geo = elementGeometry(mesh, tet);
            const Material& material = materialForTet(config, tet);
            const Vec3 normal = norm(face.normal) > 1.0e-30
                ? face.normal
                : normalized(cross(face.points[1] - face.points[0], face.points[2] - face.points[0]));

            const double faceArea = boundaryFaceArea(face);
            const double volume = std::abs(geo.detJ) / 6.0;
            if (!(faceArea > 0.0) || !(volume > 0.0)) {
                throw std::runtime_error("Invalid boundary face area or cell volume in Nitsche Dirichlet assembly.");
            }
            const double hFace = 3.0 * volume / std::max(1.0e-30, faceArea);
            constexpr double polynomialOrder = 2.0;
            const double pScale = config.penaltyScaling == "p1_squared"
                ? (polynomialOrder + 1.0) * (polynomialOrder + 1.0)
                : polynomialOrder * (polynomialOrder + 1.0);
            const double kNormal = std::max(1.0e-30, normalConductivity(material, normal));
            const double penalty = config.nitschePenaltyFactor * pScale * kNormal / std::max(1.0e-30, hFace);

            std::array<double, 100> localStiffness{};
            std::array<double, 10> localSource{};
            const double jacFace = norm(cross(face.points[1] - face.points[0], face.points[2] - face.points[0]));
            for (const TriangleQuadraturePoint& qp : quadrature) {
                const double weight = qp.weight * jacFace;
                const Vec3 q = (1.0 - qp.a - qp.b) * face.points[0]
                             + qp.a * face.points[1]
                             + qp.b * face.points[2];
                const auto lambda = lambdaOnTetFace(q, tet, face.local, mesh);
                const auto shape = shapeP2(lambda);
                const auto grad = gradShapeP2(lambda, geo);

                for (int aTest = 0; aTest < 10; ++aTest) {
                    const double v = shape[static_cast<size_t>(aTest)];
                    const double dV = material.conductivityX * grad[static_cast<size_t>(aTest)].x * normal.x
                                    + material.conductivityY * grad[static_cast<size_t>(aTest)].y * normal.y
                                    + material.conductivityZ * grad[static_cast<size_t>(aTest)].z * normal.z;
                    localSource[static_cast<size_t>(aTest)] +=
                        (-dV * condition.temperature + penalty * condition.temperature * v) * weight;

                    for (int bTrial = 0; bTrial < 10; ++bTrial) {
                        const double u = shape[static_cast<size_t>(bTrial)];
                        const double dU = material.conductivityX * grad[static_cast<size_t>(bTrial)].x * normal.x
                                        + material.conductivityY * grad[static_cast<size_t>(bTrial)].y * normal.y
                                        + material.conductivityZ * grad[static_cast<size_t>(bTrial)].z * normal.z;
                        localStiffness[static_cast<size_t>(aTest * 10 + bTrial)] +=
                            (-dU * v - dV * u + penalty * u * v) * weight;
                    }
                }
            }

            for (int a = 0; a < 10; ++a) {
                const int row = tet.dof[static_cast<size_t>(a)];
                if (std::abs(localSource[static_cast<size_t>(a)]) > 0.0) {
                    sourceEntries.push_back({row, localSource[static_cast<size_t>(a)]});
                }
                for (int b = 0; b < 10; ++b) {
                    const double value = localStiffness[static_cast<size_t>(a * 10 + b)];
                    if (std::abs(value) > 0.0) {
                        stiffnessEntries.push_back({row, tet.dof[static_cast<size_t>(b)], value});
                    }
                    if (diagnostics != nullptr) {
                        addDiagonalEntry(diagnostics->dirichletDiag,
                                         row,
                                         tet.dof[static_cast<size_t>(b)],
                                         value);
                    }
                }
            }
        }
    }

    system.appendEntries(stiffnessEntries);
    std::vector<std::vector<VectorEntry>> localSources;
    localSources.push_back(std::move(sourceEntries));
    mergeSourceEntries(localSources, source);
    return matchedFaces;
}

static Vec3 subdomainCenter(const Mesh& mesh, int subdomain)
{
    return 0.5 * (mesh.subdomainMin[static_cast<size_t>(subdomain)]
                + mesh.subdomainMax[static_cast<size_t>(subdomain)]);
}

static void assembleSipgInterface(const Mesh& mesh,
                                  const CaseConfig& config,
                                  SparseMatrix& stiffness,
                                  bool disableConsistency,
                                  bool disablePenalty,
                                  AssemblyDiagnostics* diagnostics,
                                  const std::vector<char>* disabledConsistencyFaces = nullptr)
{
    const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
    const unsigned int requestedThreads = diagnostics == nullptr ? solverParallelWorkers() : 1u;
    const int threadCount = static_cast<int>(std::max(1u, std::min<unsigned int>(
        requestedThreads, static_cast<unsigned int>(std::max<size_t>(1, mesh.interfaceFaces.size())))));
    std::vector<std::vector<MatrixEntry>> localStiffness(static_cast<size_t>(threadCount));
    const size_t facesPerThread = (mesh.interfaceFaces.size() + static_cast<size_t>(threadCount) - 1) / static_cast<size_t>(threadCount);
    for (auto& entries : localStiffness) {
        entries.reserve(facesPerThread * 400);
    }

    const auto assembleFace = [&](size_t faceIndex, std::vector<MatrixEntry>& entries) {
        const InterfaceFace& face = mesh.interfaceFaces[faceIndex];
        const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
        const ElementGeometry gLeft = elementGeometry(mesh, left);
        const ElementGeometry gRight = elementGeometry(mesh, right);
        const Material& leftMaterial = materialForTet(config, left);
        const Material& rightMaterial = materialForTet(config, right);
        const Vec3 normal = norm(face.leftNormal) > 1.0e-30
            ? face.leftNormal
            : normalized(subdomainCenter(mesh, right.subdomain) - subdomainCenter(mesh, left.subdomain));

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
        const double leftFaceArea = triangleArea(leftPhysicalFace);
        const double rightFaceArea = triangleArea(rightPhysicalFace);

        // ��� detJ �п���Ϊ��������ȡ abs
        const double leftVolume = std::abs(gLeft.detJ) / 6.0;
        const double rightVolume = std::abs(gRight.detJ) / 6.0;

        // ���෨��߶ȣ�max ģʽ��Ȼʹ��
        const double hLeft = 3.0 * leftVolume / std::max(1.0e-30, leftFaceArea);
        const double hRight = 3.0 * rightVolume / std::max(1.0e-30, rightFaceArea);

        // ����Ҫ�Ľ����������ȣ�
        // h_char = (V_L + V_R) / (A_L + A_R)
        const double areaSum = leftFaceArea + rightFaceArea;
        const double volumeSum = leftVolume + rightVolume;

        if (!(areaSum > 0.0) || !(volumeSum > 0.0)) {
            throw std::runtime_error("Invalid face area or cell volume in SIPG penalty.");
        }

        const double hChar = volumeSum / areaSum;

        // ���� hFace ���������������뱨���
        // ���� hFace �����㶨����������������ȡ�
        const double hFace = hChar;

        constexpr double polynomialOrder = 2.0;

        // max ģʽ����ԭ���� pScale
        const double pScale = config.penaltyScaling == "p1_squared"
            ? (polynomialOrder + 1.0) * (polynomialOrder + 1.0)
            : polynomialOrder * (polynomialOrder + 1.0);

        // �����ȵ���
        const double kLeftNormal = std::max(1.0e-30, normalConductivity(leftMaterial, normal));
        const double kRightNormal = std::max(1.0e-30, normalConductivity(rightMaterial, normal));

        const double alphaLeft = kLeftNormal / std::max(1.0e-30, hLeft);
        const double alphaRight = kRightNormal / std::max(1.0e-30, hRight);

        double penalty = 0.0;

        if (config.penaltyMode == "max") {
            // max ģʽ��
            // sigma = C_pen * pScale * max(k_L / h_L, k_R / h_R)
            penalty = config.penaltyFactor * pScale * std::max(alphaLeft, alphaRight);

        }
        else {
            // harmonic ģʽ������Ĺ�ʽ��
            //
            // sigma = C_pen * p^2 * k_harm / h_char
            //
            // k_harm = 2 k_L k_R / (k_L + k_R)
            // h_char = (V_L + V_R) / (A_L + A_R)

            const double kHarm = 2.0 * kLeftNormal * kRightNormal
                / std::max(1.0e-30, kLeftNormal + kRightNormal);

            const double pSquared = polynomialOrder * polynomialOrder;

            penalty = config.penaltyFactor * pSquared * kHarm / std::max(1.0e-30, hFace);
        }
        const double adjointSign = config.interfaceScheme == "nipg" ? 1.0 : -1.0;
        std::array<double, 400> consistencyMatrix{};
        std::array<double, 400> penaltyMatrix{};

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
                const auto gradLeft = gradShapeP2(lLeft, gLeft);
                const auto gradRight = gradShapeP2(lRight, gRight);

                for (int aTest = 0; aTest < 10; ++aTest) {
                    const double vL = nLeft[static_cast<size_t>(aTest)];
                    const double vR = nRight[static_cast<size_t>(aTest)];
                    const double dVL = leftMaterial.conductivityX * gradLeft[static_cast<size_t>(aTest)].x * normal.x
                                      + leftMaterial.conductivityY * gradLeft[static_cast<size_t>(aTest)].y * normal.y
                                      + leftMaterial.conductivityZ * gradLeft[static_cast<size_t>(aTest)].z * normal.z;
                    const double dVR = rightMaterial.conductivityX * gradRight[static_cast<size_t>(aTest)].x * normal.x
                                      + rightMaterial.conductivityY * gradRight[static_cast<size_t>(aTest)].y * normal.y
                                      + rightMaterial.conductivityZ * gradRight[static_cast<size_t>(aTest)].z * normal.z;

                    for (int bTrial = 0; bTrial < 10; ++bTrial) {
                        const double uL = nLeft[static_cast<size_t>(bTrial)];
                        const double uR = nRight[static_cast<size_t>(bTrial)];
                        const double dUL = leftMaterial.conductivityX * gradLeft[static_cast<size_t>(bTrial)].x * normal.x
                                          + leftMaterial.conductivityY * gradLeft[static_cast<size_t>(bTrial)].y * normal.y
                                          + leftMaterial.conductivityZ * gradLeft[static_cast<size_t>(bTrial)].z * normal.z;
                        const double dUR = rightMaterial.conductivityX * gradRight[static_cast<size_t>(bTrial)].x * normal.x
                                          + rightMaterial.conductivityY * gradRight[static_cast<size_t>(bTrial)].y * normal.y
                                          + rightMaterial.conductivityZ * gradRight[static_cast<size_t>(bTrial)].z * normal.z;

                        consistencyMatrix[static_cast<size_t>(aTest * 20 + bTrial)] +=
                            (-0.5 * dUL * vL + adjointSign * 0.5 * dVL * uL) * weight;
                        consistencyMatrix[static_cast<size_t>(aTest * 20 + 10 + bTrial)] +=
                            (-0.5 * dUR * vL - adjointSign * 0.5 * dVL * uR) * weight;
                        consistencyMatrix[static_cast<size_t>((10 + aTest) * 20 + bTrial)] +=
                            (0.5 * dUL * vR + adjointSign * 0.5 * dVR * uL) * weight;
                        consistencyMatrix[static_cast<size_t>((10 + aTest) * 20 + 10 + bTrial)] +=
                            (0.5 * dUR * vR - adjointSign * 0.5 * dVR * uR) * weight;

                        penaltyMatrix[static_cast<size_t>(aTest * 20 + bTrial)] +=
                            (penalty * uL * vL) * weight;
                        penaltyMatrix[static_cast<size_t>(aTest * 20 + 10 + bTrial)] +=
                            (-penalty * uR * vL) * weight;
                        penaltyMatrix[static_cast<size_t>((10 + aTest) * 20 + bTrial)] +=
                            (-penalty * uL * vR) * weight;
                        penaltyMatrix[static_cast<size_t>((10 + aTest) * 20 + 10 + bTrial)] +=
                            (penalty * uR * vR) * weight;
                    }
                }
            }
        }

        std::array<int, 20> dofs{};
        for (int i = 0; i < 10; ++i) {
            dofs[static_cast<size_t>(i)] = left.dof[static_cast<size_t>(i)];
            dofs[static_cast<size_t>(10 + i)] = right.dof[static_cast<size_t>(i)];
        }
        for (int row = 0; row < 20; ++row) {
            for (int col = 0; col < 20; ++col) {
                const bool faceConsistencyDisabled = disabledConsistencyFaces != nullptr
                    && faceIndex < disabledConsistencyFaces->size()
                    && (*disabledConsistencyFaces)[faceIndex] != 0;
                const double consistencyValue = (disableConsistency || faceConsistencyDisabled)
                    ? 0.0
                    : consistencyMatrix[static_cast<size_t>(row * 20 + col)];
                const double penaltyValue = disablePenalty ? 0.0 : penaltyMatrix[static_cast<size_t>(row * 20 + col)];
                const double value = consistencyValue + penaltyValue;
                if (std::abs(value) > 0.0) {
                    entries.push_back({dofs[static_cast<size_t>(row)], dofs[static_cast<size_t>(col)], value});
                }
                if (diagnostics != nullptr) {
                    addDiagonalEntry(diagnostics->interfaceConsistencyDiag,
                                     dofs[static_cast<size_t>(row)],
                                     dofs[static_cast<size_t>(col)],
                                     consistencyValue);
                    addDiagonalEntry(diagnostics->interfacePenaltyDiag,
                                     dofs[static_cast<size_t>(row)],
                                     dofs[static_cast<size_t>(col)],
                                     penaltyValue);
                }
            }
        }

        if (diagnostics != nullptr && !disablePenalty) {
            InterfacePenaltyStats& stats = diagnostics->interfacePenaltyStats;
            ++stats.facePairCount;
            stats.etaMin = std::min(stats.etaMin, penalty);
            stats.etaMax = std::max(stats.etaMax, penalty);
            stats.etaSum += penalty;
            stats.hMin = std::min(stats.hMin, hFace);
            stats.hMax = std::max(stats.hMax, hFace);
            stats.hSum += hFace;
            const double kRatio = kLeftNormal / std::max(1.0e-30, kRightNormal);
            stats.kRatioMin = std::min(stats.kRatioMin, kRatio);
            stats.kRatioMax = std::max(stats.kRatioMax, kRatio);
        }
    };

#ifdef _OPENMP
#pragma omp parallel num_threads(threadCount)
    {
        const int tid = omp_get_thread_num();
        auto& entries = localStiffness[static_cast<size_t>(tid)];
#pragma omp for schedule(static)
        for (long long i = 0; i < static_cast<long long>(mesh.interfaceFaces.size()); ++i) {
            assembleFace(static_cast<size_t>(i), entries);
        }
    }
#else
    for (size_t i = 0; i < mesh.interfaceFaces.size(); ++i) {
        assembleFace(i, localStiffness[0]);
    }
#endif

    for (auto& entries : localStiffness) {
        stiffness.appendEntries(entries);
    }
}

static double coordinateFieldValue(const Vec3& p, int field)
{
    switch (field) {
    case 1: return p.x;
    case 2: return p.y;
    case 3: return p.z;
    default: return 1.0;
    }
}

static const char* coordinateFieldName(int field)
{
    switch (field) {
    case 1: return "x";
    case 2: return "y";
    case 3: return "z";
    default: return "one";
    }
}

static double localFieldValue(const Mesh& mesh, const Tet& tet, const std::array<double, 10>& shape, int field)
{
    double value = 0.0;
    for (int i = 0; i < 10; ++i) {
        const Node& node = mesh.nodes[static_cast<size_t>(tet.dof[static_cast<size_t>(i)])];
        value += shape[static_cast<size_t>(i)] * coordinateFieldValue(node.p, field);
    }
    return value;
}

static Vec3 coordinateFieldGradient(int field)
{
    switch (field) {
    case 1: return {1.0, 0.0, 0.0};
    case 2: return {0.0, 1.0, 0.0};
    case 3: return {0.0, 0.0, 1.0};
    default: return {0.0, 0.0, 0.0};
    }
}

struct InterfaceEnergyDiagnostic {
    std::string field;
    double jumpL2 = 0.0;
    double penaltyEnergy = 0.0;
    double consistencyEnergy = 0.0;
    double totalEnergy = 0.0;
};

struct InterfaceScaleDiagnostic {
    int facePairCount = 0;
    double hLeftMin = std::numeric_limits<double>::max();
    double hLeftMax = 0.0;
    double hLeftSum = 0.0;
    double hRightMin = std::numeric_limits<double>::max();
    double hRightMax = 0.0;
    double hRightSum = 0.0;
    double alphaLeftMin = std::numeric_limits<double>::max();
    double alphaLeftMax = 0.0;
    double alphaLeftSum = 0.0;
    double alphaRightMin = std::numeric_limits<double>::max();
    double alphaRightMax = 0.0;
    double alphaRightSum = 0.0;
    double etaMin = std::numeric_limits<double>::max();
    double etaMax = 0.0;
    double etaSum = 0.0;
    double faceAreaMin = std::numeric_limits<double>::max();
    double faceAreaMax = 0.0;
    double faceAreaSum = 0.0;
};

struct InterfaceLocalDiagnostic {
    int facePairId = 0;
    int leftElementId = -1;
    int rightElementId = -1;
    int leftBoundaryId = -1;
    int rightBoundaryId = -1;
    Vec3 normalLeft{};
    Vec3 normalRight{};
    double normalDot = 0.0;
    double faceArea = 0.0;
    double hLeft = 0.0;
    double hRight = 0.0;
    double kLeftNormal = 0.0;
    double kRightNormal = 0.0;
    double alphaLeft = 0.0;
    double alphaRight = 0.0;
    double eta = 0.0;
    double consistencyDiagMin = 0.0;
    double consistencyDiagMax = 0.0;
    double penaltyDiagMin = 0.0;
    double penaltyDiagMax = 0.0;
    double localSymmetryError = 0.0;
    std::array<double, 4> energyOneXyz{};
};

struct InterfacePointDiagnostic {
    int facePairId = 0;
    int pointId = 0;
    Vec3 x{};
    std::array<double, 4> lambdaLeft{};
    std::array<double, 4> lambdaRight{};
    double lambdaLeftSum = 0.0;
    double lambdaRightSum = 0.0;
    double lambdaLeftMin = 0.0;
    double lambdaRightMin = 0.0;
    double shapeLeftSum = 0.0;
    double shapeRightSum = 0.0;
};

static void collectInterfaceDiagnostics(const Mesh& mesh,
                                        const CaseConfig& config,
                                        std::vector<InterfaceEnergyDiagnostic>& energyDiagnostics,
                                        InterfaceScaleDiagnostic& scaleDiagnostics,
                                        std::vector<InterfaceLocalDiagnostic>& localDiagnostics,
                                        std::vector<InterfacePointDiagnostic>& pointDiagnostics)
{
    energyDiagnostics.clear();
    for (int field = 0; field < 4; ++field) {
        InterfaceEnergyDiagnostic diag;
        diag.field = coordinateFieldName(field);
        energyDiagnostics.push_back(diag);
    }
    scaleDiagnostics = InterfaceScaleDiagnostic{};
    localDiagnostics.clear();
    pointDiagnostics.clear();

    const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
    int pointId = 0;
    for (size_t faceIndex = 0; faceIndex < mesh.interfaceFaces.size(); ++faceIndex) {
        const InterfaceFace& face = mesh.interfaceFaces[faceIndex];
        const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
        const ElementGeometry gLeft = elementGeometry(mesh, left);
        const ElementGeometry gRight = elementGeometry(mesh, right);
        const Material& leftMaterial = materialForTet(config, left);
        const Material& rightMaterial = materialForTet(config, right);
        const Vec3 normal = norm(face.leftNormal) > 1.0e-30
            ? face.leftNormal
            : normalized(subdomainCenter(mesh, right.subdomain) - subdomainCenter(mesh, left.subdomain));

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
        const double leftFaceArea = triangleArea(leftPhysicalFace);
        const double rightFaceArea = triangleArea(rightPhysicalFace);
        const double leftVolume = gLeft.detJ / 6.0;
        const double rightVolume = gRight.detJ / 6.0;
        const double hLeft = 3.0 * leftVolume / std::max(1.0e-30, leftFaceArea);
        const double hRight = 3.0 * rightVolume / std::max(1.0e-30, rightFaceArea);
        constexpr double polynomialOrder = 2.0;
        const double pScale = (polynomialOrder + 1.0) * (polynomialOrder + 1.0);
        const double kLeftNormal = std::max(1.0e-30, normalConductivity(leftMaterial, normal));
        const double kRightNormal = std::max(1.0e-30, normalConductivity(rightMaterial, normal));
        const double alphaLeft = kLeftNormal / std::max(1.0e-30, hLeft);
        const double alphaRight = kRightNormal / std::max(1.0e-30, hRight);
        const double eta = config.penaltyMode == "max"
            ? config.penaltyFactor * pScale * std::max(alphaLeft, alphaRight)
            : config.penaltyFactor * pScale
                * (2.0 * alphaLeft * alphaRight / std::max(1.0e-30, alphaLeft + alphaRight));
        const double adjointSign = config.interfaceScheme == "nipg" ? 1.0 : -1.0;

        ++scaleDiagnostics.facePairCount;
        scaleDiagnostics.hLeftMin = std::min(scaleDiagnostics.hLeftMin, hLeft);
        scaleDiagnostics.hLeftMax = std::max(scaleDiagnostics.hLeftMax, hLeft);
        scaleDiagnostics.hLeftSum += hLeft;
        scaleDiagnostics.hRightMin = std::min(scaleDiagnostics.hRightMin, hRight);
        scaleDiagnostics.hRightMax = std::max(scaleDiagnostics.hRightMax, hRight);
        scaleDiagnostics.hRightSum += hRight;
        scaleDiagnostics.alphaLeftMin = std::min(scaleDiagnostics.alphaLeftMin, alphaLeft);
        scaleDiagnostics.alphaLeftMax = std::max(scaleDiagnostics.alphaLeftMax, alphaLeft);
        scaleDiagnostics.alphaLeftSum += alphaLeft;
        scaleDiagnostics.alphaRightMin = std::min(scaleDiagnostics.alphaRightMin, alphaRight);
        scaleDiagnostics.alphaRightMax = std::max(scaleDiagnostics.alphaRightMax, alphaRight);
        scaleDiagnostics.alphaRightSum += alphaRight;
        scaleDiagnostics.etaMin = std::min(scaleDiagnostics.etaMin, eta);
        scaleDiagnostics.etaMax = std::max(scaleDiagnostics.etaMax, eta);
        scaleDiagnostics.etaSum += eta;
        const double faceArea = integrationArea(face.integrationTriangles);
        scaleDiagnostics.faceAreaMin = std::min(scaleDiagnostics.faceAreaMin, faceArea);
        scaleDiagnostics.faceAreaMax = std::max(scaleDiagnostics.faceAreaMax, faceArea);
        scaleDiagnostics.faceAreaSum += faceArea;

        std::array<double, 400> consistencyMatrix{};
        std::array<double, 400> penaltyMatrix{};
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
                const auto gradLeft = gradShapeP2(lLeft, gLeft);
                const auto gradRight = gradShapeP2(lRight, gRight);

                if (pointDiagnostics.size() < 20) {
                    InterfacePointDiagnostic point;
                    point.facePairId = static_cast<int>(faceIndex);
                    point.pointId = pointId++;
                    point.x = q;
                    std::array<double, 4> baryLeft{};
                    std::array<double, 4> baryRight{};
                    barycentricInTet(mesh, left, gLeft, q, baryLeft);
                    barycentricInTet(mesh, right, gRight, q, baryRight);
                    point.lambdaLeft = baryLeft;
                    point.lambdaRight = baryRight;
                    point.lambdaLeftSum = std::accumulate(baryLeft.begin(), baryLeft.end(), 0.0);
                    point.lambdaRightSum = std::accumulate(baryRight.begin(), baryRight.end(), 0.0);
                    point.lambdaLeftMin = *std::min_element(baryLeft.begin(), baryLeft.end());
                    point.lambdaRightMin = *std::min_element(baryRight.begin(), baryRight.end());
                    point.shapeLeftSum = std::accumulate(nLeft.begin(), nLeft.end(), 0.0);
                    point.shapeRightSum = std::accumulate(nRight.begin(), nRight.end(), 0.0);
                    pointDiagnostics.push_back(point);
                }

                for (int field = 0; field < 4; ++field) {
                    const double uLeft = localFieldValue(mesh, left, nLeft, field);
                    const double uRight = localFieldValue(mesh, right, nRight, field);
                    const double jump = uLeft - uRight;
                    const Vec3 gradField = coordinateFieldGradient(field);
                    const double dLeft = leftMaterial.conductivityX * gradField.x * normal.x
                                       + leftMaterial.conductivityY * gradField.y * normal.y
                                       + leftMaterial.conductivityZ * gradField.z * normal.z;
                    const double dRight = rightMaterial.conductivityX * gradField.x * normal.x
                                        + rightMaterial.conductivityY * gradField.y * normal.y
                                        + rightMaterial.conductivityZ * gradField.z * normal.z;
                    energyDiagnostics[static_cast<size_t>(field)].jumpL2 += jump * jump * weight;
                    energyDiagnostics[static_cast<size_t>(field)].penaltyEnergy += eta * jump * jump * weight;
                    energyDiagnostics[static_cast<size_t>(field)].consistencyEnergy +=
                        -(dLeft + dRight) * jump * weight;
                }

                for (int aTest = 0; aTest < 10; ++aTest) {
                    const double vL = nLeft[static_cast<size_t>(aTest)];
                    const double vR = nRight[static_cast<size_t>(aTest)];
                    const double dVL = leftMaterial.conductivityX * gradLeft[static_cast<size_t>(aTest)].x * normal.x
                                      + leftMaterial.conductivityY * gradLeft[static_cast<size_t>(aTest)].y * normal.y
                                      + leftMaterial.conductivityZ * gradLeft[static_cast<size_t>(aTest)].z * normal.z;
                    const double dVR = rightMaterial.conductivityX * gradRight[static_cast<size_t>(aTest)].x * normal.x
                                      + rightMaterial.conductivityY * gradRight[static_cast<size_t>(aTest)].y * normal.y
                                      + rightMaterial.conductivityZ * gradRight[static_cast<size_t>(aTest)].z * normal.z;
                    for (int bTrial = 0; bTrial < 10; ++bTrial) {
                        const double uL = nLeft[static_cast<size_t>(bTrial)];
                        const double uR = nRight[static_cast<size_t>(bTrial)];
                        const double dUL = leftMaterial.conductivityX * gradLeft[static_cast<size_t>(bTrial)].x * normal.x
                                          + leftMaterial.conductivityY * gradLeft[static_cast<size_t>(bTrial)].y * normal.y
                                          + leftMaterial.conductivityZ * gradLeft[static_cast<size_t>(bTrial)].z * normal.z;
                        const double dUR = rightMaterial.conductivityX * gradRight[static_cast<size_t>(bTrial)].x * normal.x
                                          + rightMaterial.conductivityY * gradRight[static_cast<size_t>(bTrial)].y * normal.y
                                          + rightMaterial.conductivityZ * gradRight[static_cast<size_t>(bTrial)].z * normal.z;
                        consistencyMatrix[static_cast<size_t>(aTest * 20 + bTrial)] +=
                            (-0.5 * dUL * vL + adjointSign * 0.5 * dVL * uL) * weight;
                        consistencyMatrix[static_cast<size_t>(aTest * 20 + 10 + bTrial)] +=
                            (-0.5 * dUR * vL - adjointSign * 0.5 * dVL * uR) * weight;
                        consistencyMatrix[static_cast<size_t>((10 + aTest) * 20 + bTrial)] +=
                            (0.5 * dUL * vR + adjointSign * 0.5 * dVR * uL) * weight;
                        consistencyMatrix[static_cast<size_t>((10 + aTest) * 20 + 10 + bTrial)] +=
                            (0.5 * dUR * vR - adjointSign * 0.5 * dVR * uR) * weight;
                        penaltyMatrix[static_cast<size_t>(aTest * 20 + bTrial)] += (eta * uL * vL) * weight;
                        penaltyMatrix[static_cast<size_t>(aTest * 20 + 10 + bTrial)] += (-eta * uR * vL) * weight;
                        penaltyMatrix[static_cast<size_t>((10 + aTest) * 20 + bTrial)] += (-eta * uL * vR) * weight;
                        penaltyMatrix[static_cast<size_t>((10 + aTest) * 20 + 10 + bTrial)] += (eta * uR * vR) * weight;
                    }
                }
            }
        }
        for (InterfaceEnergyDiagnostic& diag : energyDiagnostics) {
            diag.totalEnergy = diag.consistencyEnergy + diag.penaltyEnergy;
        }

        if (localDiagnostics.size() < 20) {
            InterfaceLocalDiagnostic local;
            local.facePairId = static_cast<int>(faceIndex);
            local.leftElementId = face.leftTet;
            local.rightElementId = face.rightTet;
            local.leftBoundaryId = face.leftBoundaryEntity;
            local.rightBoundaryId = face.rightBoundaryEntity;
            local.normalLeft = face.leftNormal;
            local.normalRight = face.rightNormal;
            local.normalDot = dot(face.leftNormal, face.rightNormal);
            local.faceArea = faceArea;
            local.hLeft = hLeft;
            local.hRight = hRight;
            local.kLeftNormal = kLeftNormal;
            local.kRightNormal = kRightNormal;
            local.alphaLeft = alphaLeft;
            local.alphaRight = alphaRight;
            local.eta = eta;
            local.consistencyDiagMin = std::numeric_limits<double>::max();
            local.consistencyDiagMax = -std::numeric_limits<double>::max();
            local.penaltyDiagMin = std::numeric_limits<double>::max();
            local.penaltyDiagMax = -std::numeric_limits<double>::max();
            double normSq = 0.0;
            double asymSq = 0.0;
            for (int i = 0; i < 20; ++i) {
                local.consistencyDiagMin = std::min(local.consistencyDiagMin, consistencyMatrix[static_cast<size_t>(i * 20 + i)]);
                local.consistencyDiagMax = std::max(local.consistencyDiagMax, consistencyMatrix[static_cast<size_t>(i * 20 + i)]);
                local.penaltyDiagMin = std::min(local.penaltyDiagMin, penaltyMatrix[static_cast<size_t>(i * 20 + i)]);
                local.penaltyDiagMax = std::max(local.penaltyDiagMax, penaltyMatrix[static_cast<size_t>(i * 20 + i)]);
                for (int j = 0; j < 20; ++j) {
                    const double value = consistencyMatrix[static_cast<size_t>(i * 20 + j)]
                                       + penaltyMatrix[static_cast<size_t>(i * 20 + j)];
                    const double tvalue = consistencyMatrix[static_cast<size_t>(j * 20 + i)]
                                        + penaltyMatrix[static_cast<size_t>(j * 20 + i)];
                    normSq += value * value;
                    const double diff = value - tvalue;
                    asymSq += diff * diff;
                }
            }
            local.localSymmetryError = std::sqrt(asymSq) / std::sqrt(std::max(1.0e-300, normSq));
            std::array<double, 20> fieldValues{};
            for (int field = 0; field < 4; ++field) {
                for (int i = 0; i < 10; ++i) {
                    fieldValues[static_cast<size_t>(i)] =
                        coordinateFieldValue(mesh.nodes[static_cast<size_t>(left.dof[static_cast<size_t>(i)])].p, field);
                    fieldValues[static_cast<size_t>(10 + i)] =
                        coordinateFieldValue(mesh.nodes[static_cast<size_t>(right.dof[static_cast<size_t>(i)])].p, field);
                }
                double energy = 0.0;
                for (int i = 0; i < 20; ++i) {
                    for (int j = 0; j < 20; ++j) {
                        energy += fieldValues[static_cast<size_t>(i)]
                            * (consistencyMatrix[static_cast<size_t>(i * 20 + j)]
                               + penaltyMatrix[static_cast<size_t>(i * 20 + j)])
                            * fieldValues[static_cast<size_t>(j)];
                    }
                }
                local.energyOneXyz[static_cast<size_t>(field)] = energy;
            }
            localDiagnostics.push_back(local);
        }
    }
    for (InterfaceScaleDiagnostic* stats : {&scaleDiagnostics}) {
        if (stats->facePairCount == 0) {
            stats->hLeftMin = stats->hRightMin = stats->alphaLeftMin = stats->alphaRightMin
                = stats->etaMin = stats->faceAreaMin = 0.0;
        }
    }
}

static const char* interfaceFieldLabel(const std::string& field)
{
    if (field == "one") {
        return "T=1";
    }
    if (field == "x") {
        return "T=x";
    }
    if (field == "y") {
        return "T=y";
    }
    if (field == "z") {
        return "T=z";
    }
    return "T=?";
}

static bool barycentricOutOfBounds(const std::array<double, 4>& lambda)
{
    constexpr double tolerance = 1.0e-8;
    const double sum = std::accumulate(lambda.begin(), lambda.end(), 0.0);
    const double minLambda = *std::min_element(lambda.begin(), lambda.end());
    const double maxLambda = *std::max_element(lambda.begin(), lambda.end());
    return minLambda < -tolerance || maxLambda > 1.0 + tolerance || std::abs(sum - 1.0) > tolerance;
}

static void printInterfaceContinuousFieldDiagnostics(const std::vector<InterfaceEnergyDiagnostic>& diagnostics)
{
    std::cout << "interface continuous field test\n";
    for (const InterfaceEnergyDiagnostic& diag : diagnostics) {
        const double jumpL2 = std::sqrt(std::max(0.0, diag.jumpL2));
        std::cout << "  " << interfaceFieldLabel(diag.field) << ":\n"
                  << "    jump_L2_on_interface=" << jumpL2 << "\n"
                  << "    consistency_energy=" << diag.consistencyEnergy << "\n"
                  << "    penalty_energy=" << diag.penaltyEnergy << "\n"
                  << "    total_interface_energy=" << diag.totalEnergy << "\n";
    }
}

static void writeInterfaceContinuousFieldDiagnostics(const std::vector<InterfaceEnergyDiagnostic>& diagnostics,
                                                     const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "field,jump_L2_on_interface,consistency_energy,penalty_energy,total_interface_energy\n";
    for (const InterfaceEnergyDiagnostic& diag : diagnostics) {
        const double jumpL2 = std::sqrt(std::max(0.0, diag.jumpL2));
        out << interfaceFieldLabel(diag.field) << ','
            << jumpL2 << ','
            << diag.consistencyEnergy << ','
            << diag.penaltyEnergy << ','
            << diag.totalEnergy << '\n';
    }
}

static void writeLocalInterfacePairDiagnostics(const std::vector<InterfaceLocalDiagnostic>& diagnostics,
                                               const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "face_pair_id,left_element_id,right_element_id,left_boundary_id,right_boundary_id,"
        << "normal_L_x,normal_L_y,normal_L_z,normal_R_x,normal_R_y,normal_R_z,"
        << "dot_normal_L_normal_R,face_area,h_L,h_R,k_L_normal,k_R_normal,alpha_L,alpha_R,eta,"
        << "local_consistency_diagonal_min,local_consistency_diagonal_max,"
        << "local_penalty_diagonal_min,local_penalty_diagonal_max,local_symmetry_error,"
        << "energy_T_1,energy_T_x,energy_T_y,energy_T_z\n";
    for (const InterfaceLocalDiagnostic& diag : diagnostics) {
        out << diag.facePairId << ','
            << diag.leftElementId << ','
            << diag.rightElementId << ','
            << diag.leftBoundaryId << ','
            << diag.rightBoundaryId << ','
            << diag.normalLeft.x << ',' << diag.normalLeft.y << ',' << diag.normalLeft.z << ','
            << diag.normalRight.x << ',' << diag.normalRight.y << ',' << diag.normalRight.z << ','
            << diag.normalDot << ','
            << diag.faceArea << ','
            << diag.hLeft << ','
            << diag.hRight << ','
            << diag.kLeftNormal << ','
            << diag.kRightNormal << ','
            << diag.alphaLeft << ','
            << diag.alphaRight << ','
            << diag.eta << ','
            << diag.consistencyDiagMin << ','
            << diag.consistencyDiagMax << ','
            << diag.penaltyDiagMin << ','
            << diag.penaltyDiagMax << ','
            << diag.localSymmetryError << ','
            << diag.energyOneXyz[0] << ','
            << diag.energyOneXyz[1] << ','
            << diag.energyOneXyz[2] << ','
            << diag.energyOneXyz[3] << '\n';
    }
}

static void writeBarycentricDiagnostics(const std::vector<InterfacePointDiagnostic>& diagnostics,
                                        const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "face_pair_id,point_id,x_q,y_q,z_q,"
        << "lambda_L_0,lambda_L_1,lambda_L_2,lambda_L_3,"
        << "lambda_R_0,lambda_R_1,lambda_R_2,lambda_R_3,"
        << "lambda_L_sum,lambda_R_sum,N_L_sum,N_R_sum,left_out_of_bounds,right_out_of_bounds\n";
    for (const InterfacePointDiagnostic& diag : diagnostics) {
        out << diag.facePairId << ','
            << diag.pointId << ','
            << diag.x.x << ',' << diag.x.y << ',' << diag.x.z << ','
            << diag.lambdaLeft[0] << ',' << diag.lambdaLeft[1] << ','
            << diag.lambdaLeft[2] << ',' << diag.lambdaLeft[3] << ','
            << diag.lambdaRight[0] << ',' << diag.lambdaRight[1] << ','
            << diag.lambdaRight[2] << ',' << diag.lambdaRight[3] << ','
            << diag.lambdaLeftSum << ','
            << diag.lambdaRightSum << ','
            << diag.shapeLeftSum << ','
            << diag.shapeRightSum << ','
            << (barycentricOutOfBounds(diag.lambdaLeft) ? 1 : 0) << ','
            << (barycentricOutOfBounds(diag.lambdaRight) ? 1 : 0) << '\n';
    }
}

static Vec3 triangleCenter(const std::array<Vec3, 3>& tri)
{
    return (tri[0] + tri[1] + tri[2]) / 3.0;
}

static std::array<Vec3, 3> boundaryTriangleFromTetFace(const Mesh& mesh,
                                                       const Tet& tet,
                                                       const std::array<int, 3>& local)
{
    return {{
        mesh.nodes[static_cast<size_t>(tet.v[static_cast<size_t>(local[0])])].p,
        mesh.nodes[static_cast<size_t>(tet.v[static_cast<size_t>(local[1])])].p,
        mesh.nodes[static_cast<size_t>(tet.v[static_cast<size_t>(local[2])])].p
    }};
}

static bool sameLocalFace(const std::array<int, 3>& a, const std::array<int, 3>& b)
{
    return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

static int boundaryFaceIndexForInterfaceSide(const Mesh& mesh,
                                             int tet,
                                             int boundaryEntity,
                                             const std::array<int, 3>& local)
{
    for (size_t i = 0; i < mesh.boundaryFaces.size(); ++i) {
        const BoundaryFace& face = mesh.boundaryFaces[i];
        if (face.tet == tet
            && face.boundaryEntity == boundaryEntity
            && sameLocalFace(face.local, local)) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

struct InterfaceFacePairDiagnosticRow {
    int pairId = -1;
    int leftDomainId = -1;
    int rightDomainId = -1;
    int leftFaceId = -1;
    int rightFaceId = -1;
    Vec3 leftCenter{};
    Vec3 rightCenter{};
    double centerDistance = 0.0;
    double leftArea = 0.0;
    double rightArea = 0.0;
    double areaRatio = 0.0;
    Vec3 leftNormal{};
    Vec3 rightNormal{};
    double normalDot = 0.0;
    double kLeft = 0.0;
    double kRight = 0.0;
    double hLeft = 0.0;
    double hRight = 0.0;
    double sigma = 0.0;
    double penaltyDiagMin = 0.0;
    double penaltyDiagMax = 0.0;
    double consistencyDiagMin = 0.0;
    double consistencyDiagMax = 0.0;
    double consistencyPenaltyRatio = 0.0;
    double consistencyPenaltyFrobeniusRatio = 0.0;
    int suspicious = 0;
    std::string suspiciousReason;
};

static std::vector<InterfaceFacePairDiagnosticRow> collectInterfaceFacePairDiagnostics(const Mesh& mesh,
                                                                                       const CaseConfig& config)
{
    const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
    std::vector<InterfaceFacePairDiagnosticRow> rows;
    rows.reserve(mesh.interfaceFaces.size());
    for (size_t faceIndex = 0; faceIndex < mesh.interfaceFaces.size(); ++faceIndex) {
        const InterfaceFace& face = mesh.interfaceFaces[faceIndex];
        const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
        const ElementGeometry gLeft = elementGeometry(mesh, left);
        const ElementGeometry gRight = elementGeometry(mesh, right);
        const Material& leftMaterial = materialForTet(config, left);
        const Material& rightMaterial = materialForTet(config, right);
        const Vec3 normal = norm(face.leftNormal) > 1.0e-30
            ? face.leftNormal
            : normalized(subdomainCenter(mesh, right.subdomain) - subdomainCenter(mesh, left.subdomain));
        const std::array<Vec3, 3> leftPhysicalFace =
            boundaryTriangleFromTetFace(mesh, left, face.leftLocal);
        const std::array<Vec3, 3> rightPhysicalFace =
            boundaryTriangleFromTetFace(mesh, right, face.rightLocal);
        const double leftFaceArea = triangleArea(leftPhysicalFace);
        const double rightFaceArea = triangleArea(rightPhysicalFace);
        const double leftVolume = gLeft.detJ / 6.0;
        const double rightVolume = gRight.detJ / 6.0;
        const double hLeft = 3.0 * leftVolume / std::max(1.0e-30, leftFaceArea);
        const double hRight = 3.0 * rightVolume / std::max(1.0e-30, rightFaceArea);
        const double hFace = std::min(hLeft, hRight);
        constexpr double polynomialOrder = 2.0;
        const double pScale = config.penaltyScaling == "p1_squared"
            ? (polynomialOrder + 1.0) * (polynomialOrder + 1.0)
            : polynomialOrder * (polynomialOrder + 1.0);
        const double kLeftNormal = std::max(1.0e-30, normalConductivity(leftMaterial, normal));
        const double kRightNormal = std::max(1.0e-30, normalConductivity(rightMaterial, normal));
        const double alphaLeft = kLeftNormal / std::max(1.0e-30, hLeft);
        const double alphaRight = kRightNormal / std::max(1.0e-30, hRight);
        const double penalty = config.penaltyMode == "max"
            ? config.penaltyFactor * pScale * std::max(alphaLeft, alphaRight)
            : config.penaltyFactor * pScale
                * (2.0 * alphaLeft * alphaRight / std::max(1.0e-30, alphaLeft + alphaRight));
        const double adjointSign = config.interfaceScheme == "nipg" ? 1.0 : -1.0;
        std::array<double, 400> consistencyMatrix{};
        std::array<double, 400> penaltyMatrix{};
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
                const auto gradLeft = gradShapeP2(lLeft, gLeft);
                const auto gradRight = gradShapeP2(lRight, gRight);
                for (int aTest = 0; aTest < 10; ++aTest) {
                    const double vL = nLeft[static_cast<size_t>(aTest)];
                    const double vR = nRight[static_cast<size_t>(aTest)];
                    const double dVL = leftMaterial.conductivityX * gradLeft[static_cast<size_t>(aTest)].x * normal.x
                                      + leftMaterial.conductivityY * gradLeft[static_cast<size_t>(aTest)].y * normal.y
                                      + leftMaterial.conductivityZ * gradLeft[static_cast<size_t>(aTest)].z * normal.z;
                    const double dVR = rightMaterial.conductivityX * gradRight[static_cast<size_t>(aTest)].x * normal.x
                                      + rightMaterial.conductivityY * gradRight[static_cast<size_t>(aTest)].y * normal.y
                                      + rightMaterial.conductivityZ * gradRight[static_cast<size_t>(aTest)].z * normal.z;
                    for (int bTrial = 0; bTrial < 10; ++bTrial) {
                        const double uL = nLeft[static_cast<size_t>(bTrial)];
                        const double uR = nRight[static_cast<size_t>(bTrial)];
                        const double dUL = leftMaterial.conductivityX * gradLeft[static_cast<size_t>(bTrial)].x * normal.x
                                          + leftMaterial.conductivityY * gradLeft[static_cast<size_t>(bTrial)].y * normal.y
                                          + leftMaterial.conductivityZ * gradLeft[static_cast<size_t>(bTrial)].z * normal.z;
                        const double dUR = rightMaterial.conductivityX * gradRight[static_cast<size_t>(bTrial)].x * normal.x
                                          + rightMaterial.conductivityY * gradRight[static_cast<size_t>(bTrial)].y * normal.y
                                          + rightMaterial.conductivityZ * gradRight[static_cast<size_t>(bTrial)].z * normal.z;
                        consistencyMatrix[static_cast<size_t>(aTest * 20 + bTrial)] +=
                            (-0.5 * dUL * vL + adjointSign * 0.5 * dVL * uL) * weight;
                        consistencyMatrix[static_cast<size_t>(aTest * 20 + 10 + bTrial)] +=
                            (-0.5 * dUR * vL - adjointSign * 0.5 * dVL * uR) * weight;
                        consistencyMatrix[static_cast<size_t>((10 + aTest) * 20 + bTrial)] +=
                            (0.5 * dUL * vR + adjointSign * 0.5 * dVR * uL) * weight;
                        consistencyMatrix[static_cast<size_t>((10 + aTest) * 20 + 10 + bTrial)] +=
                            (0.5 * dUR * vR - adjointSign * 0.5 * dVR * uR) * weight;
                        penaltyMatrix[static_cast<size_t>(aTest * 20 + bTrial)] +=
                            (penalty * uL * vL) * weight;
                        penaltyMatrix[static_cast<size_t>(aTest * 20 + 10 + bTrial)] +=
                            (-penalty * uR * vL) * weight;
                        penaltyMatrix[static_cast<size_t>((10 + aTest) * 20 + bTrial)] +=
                            (-penalty * uL * vR) * weight;
                        penaltyMatrix[static_cast<size_t>((10 + aTest) * 20 + 10 + bTrial)] +=
                            (penalty * uR * vR) * weight;
                    }
                }
            }
        }

        InterfaceFacePairDiagnosticRow row;
        row.pairId = static_cast<int>(faceIndex);
        row.leftDomainId = left.subdomain;
        row.rightDomainId = right.subdomain;
        row.leftFaceId = face.leftFaceId >= 0
            ? face.leftFaceId
            : boundaryFaceIndexForInterfaceSide(mesh, face.leftTet, face.leftBoundaryEntity, face.leftLocal);
        row.rightFaceId = face.rightFaceId >= 0
            ? face.rightFaceId
            : boundaryFaceIndexForInterfaceSide(mesh, face.rightTet, face.rightBoundaryEntity, face.rightLocal);
        row.leftCenter = triangleCenter(leftPhysicalFace);
        row.rightCenter = triangleCenter(rightPhysicalFace);
        row.centerDistance = norm(row.leftCenter - row.rightCenter);
        row.leftArea = leftFaceArea;
        row.rightArea = rightFaceArea;
        row.areaRatio = std::max(leftFaceArea, rightFaceArea) / std::max(1.0e-300, std::min(leftFaceArea, rightFaceArea));
        row.leftNormal = face.leftNormal;
        row.rightNormal = face.rightNormal;
        row.normalDot = dot(face.leftNormal, face.rightNormal);
        row.kLeft = kLeftNormal;
        row.kRight = kRightNormal;
        row.hLeft = hLeft;
        row.hRight = hRight;
        row.sigma = penalty;
        row.penaltyDiagMin = std::numeric_limits<double>::max();
        row.penaltyDiagMax = -std::numeric_limits<double>::max();
        row.consistencyDiagMin = std::numeric_limits<double>::max();
        row.consistencyDiagMax = -std::numeric_limits<double>::max();
        double maxAbsConsistencyDiag = 0.0;
        double maxAbsPenaltyDiag = 0.0;
        double consistencyFrobenius = 0.0;
        double penaltyFrobenius = 0.0;
        for (int i = 0; i < 20; ++i) {
            const double cDiag = consistencyMatrix[static_cast<size_t>(i * 20 + i)];
            const double pDiag = penaltyMatrix[static_cast<size_t>(i * 20 + i)];
            row.consistencyDiagMin = std::min(row.consistencyDiagMin, cDiag);
            row.consistencyDiagMax = std::max(row.consistencyDiagMax, cDiag);
            row.penaltyDiagMin = std::min(row.penaltyDiagMin, pDiag);
            row.penaltyDiagMax = std::max(row.penaltyDiagMax, pDiag);
            maxAbsConsistencyDiag = std::max(maxAbsConsistencyDiag, std::abs(cDiag));
            maxAbsPenaltyDiag = std::max(maxAbsPenaltyDiag, std::abs(pDiag));
        }
        for (int i = 0; i < 400; ++i) {
            consistencyFrobenius += consistencyMatrix[static_cast<size_t>(i)] * consistencyMatrix[static_cast<size_t>(i)];
            penaltyFrobenius += penaltyMatrix[static_cast<size_t>(i)] * penaltyMatrix[static_cast<size_t>(i)];
        }
        row.consistencyPenaltyRatio = maxAbsConsistencyDiag / std::max(1.0e-300, maxAbsPenaltyDiag);
        row.consistencyPenaltyFrobeniusRatio =
            std::sqrt(consistencyFrobenius) / std::sqrt(std::max(1.0e-300, penaltyFrobenius));
        std::vector<std::string> reasons;
        if (row.normalDot > -0.99) {
            reasons.push_back("normal_dot_not_near_minus_one");
        }
        if (row.centerDistance > 1.0e-8) {
            reasons.push_back("large_center_distance");
        }
        if (row.areaRatio > 10.0) {
            reasons.push_back("large_area_ratio");
        }
        if (!(row.sigma > 0.0) || !std::isfinite(row.sigma)) {
            reasons.push_back("nonpositive_or_nonfinite_sigma");
        }
        if (row.penaltyDiagMax <= 0.0 || row.penaltyDiagMin < -1.0e-20) {
            reasons.push_back("bad_penalty_diagonal");
        }
        if (row.consistencyPenaltyRatio > 1.0 || row.consistencyPenaltyFrobeniusRatio > 1.0) {
            reasons.push_back("consistency_large_vs_penalty");
        }
        row.suspicious = reasons.empty() ? 0 : 1;
        for (size_t i = 0; i < reasons.size(); ++i) {
            if (i > 0) {
                row.suspiciousReason += ';';
            }
            row.suspiciousReason += reasons[i];
        }
        rows.push_back(row);
    }
    return rows;
}

static void writeInterfaceFacePairDiagnostics(const std::vector<InterfaceFacePairDiagnosticRow>& rows,
                                              const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "pair_id,left_domain_id,right_domain_id,left_face_id,right_face_id,"
        << "left_center_x,left_center_y,left_center_z,right_center_x,right_center_y,right_center_z,"
        << "center_distance,left_area,right_area,area_ratio,"
        << "left_normal_x,left_normal_y,left_normal_z,right_normal_x,right_normal_y,right_normal_z,"
        << "dot_n_left_n_right,k_left,k_right,h_left,h_right,sigma,"
        << "penalty_diagonal_contribution_min,penalty_diagonal_contribution_max,"
        << "consistency_diagonal_contribution_min,consistency_diagonal_contribution_max,"
        << "consistency_to_penalty_ratio,consistency_to_penalty_frobenius_ratio,"
        << "suspicious_face,suspicious_reason\n";
    out << std::setprecision(16);
    for (const InterfaceFacePairDiagnosticRow& row : rows) {
        out << row.pairId << ','
            << row.leftDomainId << ','
            << row.rightDomainId << ','
            << row.leftFaceId << ','
            << row.rightFaceId << ','
            << row.leftCenter.x << ',' << row.leftCenter.y << ',' << row.leftCenter.z << ','
            << row.rightCenter.x << ',' << row.rightCenter.y << ',' << row.rightCenter.z << ','
            << row.centerDistance << ','
            << row.leftArea << ','
            << row.rightArea << ','
            << row.areaRatio << ','
            << row.leftNormal.x << ',' << row.leftNormal.y << ',' << row.leftNormal.z << ','
            << row.rightNormal.x << ',' << row.rightNormal.y << ',' << row.rightNormal.z << ','
            << row.normalDot << ','
            << row.kLeft << ','
            << row.kRight << ','
            << row.hLeft << ','
            << row.hRight << ','
            << row.sigma << ','
            << row.penaltyDiagMin << ','
            << row.penaltyDiagMax << ','
            << row.consistencyDiagMin << ','
            << row.consistencyDiagMax << ','
            << row.consistencyPenaltyRatio << ','
            << row.consistencyPenaltyFrobeniusRatio << ','
            << row.suspicious << ','
            << csvEscape(row.suspiciousReason) << '\n';
    }
}

static void writeNonmatchingInterfaceProjectionDiagnostics(const std::vector<InterfaceFacePairDiagnosticRow>& rows,
                                                           const std::filesystem::path& path)
{
    std::map<int, int> leftMultiplicity;
    std::map<int, int> rightMultiplicity;
    double areaRatioMin = std::numeric_limits<double>::infinity();
    double areaRatioMax = 0.0;
    double areaRatioSum = 0.0;
    double centerMin = std::numeric_limits<double>::infinity();
    double centerMax = 0.0;
    double centerSum = 0.0;
    int largeAreaRatioCount = 0;
    int largeCenterDistanceCount = 0;
    double leftAreaSum = 0.0;
    double rightAreaSum = 0.0;
    for (const InterfaceFacePairDiagnosticRow& row : rows) {
        ++leftMultiplicity[row.leftFaceId];
        ++rightMultiplicity[row.rightFaceId];
        areaRatioMin = std::min(areaRatioMin, row.areaRatio);
        areaRatioMax = std::max(areaRatioMax, row.areaRatio);
        areaRatioSum += row.areaRatio;
        centerMin = std::min(centerMin, row.centerDistance);
        centerMax = std::max(centerMax, row.centerDistance);
        centerSum += row.centerDistance;
        if (row.areaRatio > 3.0) {
            ++largeAreaRatioCount;
        }
        if (row.centerDistance > 1.0e-8) {
            ++largeCenterDistanceCount;
        }
        leftAreaSum += row.leftArea;
        rightAreaSum += row.rightArea;
    }
    const int oneLeftToManyRight = static_cast<int>(std::count_if(leftMultiplicity.begin(), leftMultiplicity.end(),
        [](const auto& entry) { return entry.second > 1; }));
    const int oneRightToManyLeft = static_cast<int>(std::count_if(rightMultiplicity.begin(), rightMultiplicity.end(),
        [](const auto& entry) { return entry.second > 1; }));
    const double count = static_cast<double>(std::max<size_t>(1, rows.size()));
    std::ofstream out(path);
    out << "metric,value\n";
    out << std::setprecision(16);
    out << "face_pair_count," << rows.size() << '\n';
    out << "unique_left_faces," << leftMultiplicity.size() << '\n';
    out << "unique_right_faces," << rightMultiplicity.size() << '\n';
    out << "left_faces_with_multiple_right_faces," << oneLeftToManyRight << '\n';
    out << "right_faces_with_multiple_left_faces," << oneRightToManyLeft << '\n';
    out << "area_ratio_min," << (std::isfinite(areaRatioMin) ? areaRatioMin : 0.0) << '\n';
    out << "area_ratio_avg," << areaRatioSum / count << '\n';
    out << "area_ratio_max," << areaRatioMax << '\n';
    out << "center_distance_min," << (std::isfinite(centerMin) ? centerMin : 0.0) << '\n';
    out << "center_distance_avg," << centerSum / count << '\n';
    out << "center_distance_max," << centerMax << '\n';
    out << "large_area_ratio_count_gt_3," << largeAreaRatioCount << '\n';
    out << "large_center_distance_count_gt_1e-8," << largeCenterDistanceCount << '\n';
    out << "summed_left_face_area_over_pairs," << leftAreaSum << '\n';
    out << "summed_right_face_area_over_pairs," << rightAreaSum << '\n';
    out << "uses_true_geometric_overlap_integration,1\n";
    out << "uses_mortar_projection,0\n";
    out << "potential_nonmatching_interface_integration_error,"
        << ((oneLeftToManyRight > 0 || oneRightToManyLeft > 0 || largeAreaRatioCount > 0 || largeCenterDistanceCount > 0) ? 1 : 0)
        << '\n';
}

static std::string joinInts(const std::set<int>& values)
{
    std::ostringstream out;
    bool first = true;
    for (int value : values) {
        if (!first) {
            out << ';';
        }
        out << value;
        first = false;
    }
    return out.str();
}

static long long quantizedCellKey(const Vec3& p, double cellSize)
{
    const long long ix = static_cast<long long>(std::floor(p.x / cellSize));
    const long long iy = static_cast<long long>(std::floor(p.y / cellSize));
    const long long iz = static_cast<long long>(std::floor(p.z / cellSize));
    return ((ix & 0x1fffffLL) << 42) ^ ((iy & 0x1fffffLL) << 21) ^ (iz & 0x1fffffLL);
}

static void writeMissingInterfaceDiagnostics(const Mesh& mesh,
                                             const CaseConfig& config,
                                             const std::vector<InterfaceFacePairDiagnosticRow>& pairRows,
                                             const std::filesystem::path& path)
{
    std::set<int> configuredLeftFaces;
    std::set<int> configuredRightFaces;
    std::set<int> detectedLeftFaces;
    std::set<int> detectedRightFaces;
    std::map<std::pair<int, int>, int> pairCounts;
    for (const InterfaceFacePairDiagnosticRow& row : pairRows) {
        if (row.leftFaceId >= 0) {
            configuredLeftFaces.insert(row.leftFaceId);
            detectedLeftFaces.insert(row.leftFaceId);
        }
        if (row.rightFaceId >= 0) {
            configuredRightFaces.insert(row.rightFaceId);
            detectedRightFaces.insert(row.rightFaceId);
        }
        ++pairCounts[{row.leftFaceId, row.rightFaceId}];
    }

    std::vector<int> rightCandidateFaces;
    rightCandidateFaces.reserve(mesh.boundaryFaces.size());
    std::map<long long, std::vector<int>> rightByCell;
    constexpr double contactCellSize = 5.0e-8;
    for (size_t j = 0; j < mesh.boundaryFaces.size(); ++j) {
        const BoundaryFace& right = mesh.boundaryFaces[j];
        if (right.subdomain != 1 || isDirichletBoundary(right, config)) {
            continue;
        }
        rightCandidateFaces.push_back(static_cast<int>(j));
        rightByCell[quantizedCellKey(triangleCenter(right.points), contactCellSize)].push_back(static_cast<int>(j));
    }

    for (size_t i = 0; i < mesh.boundaryFaces.size(); ++i) {
        const BoundaryFace& left = mesh.boundaryFaces[i];
        if (left.subdomain != 0 || isDirichletBoundary(left, config)) {
            continue;
        }
        const Vec3 center = triangleCenter(left.points);
        std::set<int> candidateIds;
        for (int dx = -2; dx <= 2; ++dx) {
            for (int dy = -2; dy <= 2; ++dy) {
                for (int dz = -2; dz <= 2; ++dz) {
                    const Vec3 shifted{
                        center.x + static_cast<double>(dx) * contactCellSize,
                        center.y + static_cast<double>(dy) * contactCellSize,
                        center.z + static_cast<double>(dz) * contactCellSize
                    };
                    const auto found = rightByCell.find(quantizedCellKey(shifted, contactCellSize));
                    if (found != rightByCell.end()) {
                        candidateIds.insert(found->second.begin(), found->second.end());
                    }
                }
            }
        }
        for (int rightId : candidateIds) {
            const BoundaryFace& right = mesh.boundaryFaces[static_cast<size_t>(rightId)];
            std::vector<std::array<Vec3, 3>> pieces = intersectTrianglesAsFan(left.points, right.points);
            if (!pieces.empty()) {
                detectedLeftFaces.insert(static_cast<int>(i));
                detectedRightFaces.insert(rightId);
            }
        }
    }

    std::set<int> missingLeft;
    std::set<int> missingRight;
    std::set_difference(detectedLeftFaces.begin(), detectedLeftFaces.end(),
                        configuredLeftFaces.begin(), configuredLeftFaces.end(),
                        std::inserter(missingLeft, missingLeft.begin()));
    std::set_difference(detectedRightFaces.begin(), detectedRightFaces.end(),
                        configuredRightFaces.begin(), configuredRightFaces.end(),
                        std::inserter(missingRight, missingRight.begin()));
    std::set<int> unmatchedLeft;
    std::set<int> unmatchedRight;
    std::set_difference(configuredLeftFaces.begin(), configuredLeftFaces.end(),
                        detectedLeftFaces.begin(), detectedLeftFaces.end(),
                        std::inserter(unmatchedLeft, unmatchedLeft.begin()));
    std::set_difference(configuredRightFaces.begin(), configuredRightFaces.end(),
                        detectedRightFaces.begin(), detectedRightFaces.end(),
                        std::inserter(unmatchedRight, unmatchedRight.begin()));
    int duplicatePairs = 0;
    for (const auto& entry : pairCounts) {
        if (entry.second > 1) {
            ++duplicatePairs;
        }
    }
    const int suspiciousPairs = static_cast<int>(std::count_if(pairRows.begin(), pairRows.end(),
        [](const InterfaceFacePairDiagnosticRow& row) { return row.suspicious != 0; }));

    std::ofstream report(path.parent_path() / "rram_interface_missing_faces_report.csv");
    report << "side,face_id,subdomain,boundary_entity,status,center_x,center_y,center_z,area,"
           << "normal_x,normal_y,normal_z\n";
    report << std::setprecision(16);
    auto writeFaceReport = [&](const std::string& side,
                               const std::set<int>& faceIds,
                               const std::set<int>& configured,
                               const std::set<int>& missing) {
        for (int faceId : faceIds) {
            const BoundaryFace& face = mesh.boundaryFaces[static_cast<size_t>(faceId)];
            std::string status = "detected";
            if (missing.find(faceId) != missing.end()) {
                status = "missing";
            } else if (configured.find(faceId) != configured.end()) {
                status = "configured";
            }
            const Vec3 center = triangleCenter(face.points);
            report << side << ','
                   << faceId << ','
                   << face.subdomain << ','
                   << face.boundaryEntity << ','
                   << status << ','
                   << center.x << ','
                   << center.y << ','
                   << center.z << ','
                   << boundaryFaceArea(face) << ','
                   << face.normal.x << ','
                   << face.normal.y << ','
                   << face.normal.z << '\n';
        }
    };
    writeFaceReport("left", detectedLeftFaces, configuredLeftFaces, missingLeft);
    writeFaceReport("right", detectedRightFaces, configuredRightFaces, missingRight);

    std::ofstream out(path);
    out << "detected_contact_faces_left,detected_contact_faces_right,"
        << "configured_interface_faces_left,configured_interface_faces_right,"
        << "missing_left_faces,missing_right_faces,unmatched_left_faces,unmatched_right_faces,"
        << "duplicate_pairs,suspicious_pairs,missing_left_face_ids,missing_right_face_ids,"
        << "unmatched_left_face_ids,unmatched_right_face_ids\n";
    out << detectedLeftFaces.size() << ','
        << detectedRightFaces.size() << ','
        << configuredLeftFaces.size() << ','
        << configuredRightFaces.size() << ','
        << missingLeft.size() << ','
        << missingRight.size() << ','
        << unmatchedLeft.size() << ','
        << unmatchedRight.size() << ','
        << duplicatePairs << ','
        << suspiciousPairs << ','
        << csvEscape(joinInts(missingLeft)) << ','
        << csvEscape(joinInts(missingRight)) << ','
        << csvEscape(joinInts(unmatchedLeft)) << ','
        << csvEscape(joinInts(unmatchedRight)) << '\n';
}
