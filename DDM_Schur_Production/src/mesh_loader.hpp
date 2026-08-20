#pragma once

// COMSOL text mesh import, deterministic P2 edge-node completion, boundary
// extraction, translated template instances, and nonmatching interface pairing.
// Global DOF order produced here is part of every descriptor/model fingerprint.

static int nextNumericLine(const std::vector<std::string>& lines, int start)
{
    for (int i = start; i < static_cast<int>(lines.size()); ++i) {
        if (!numericPart(lines[static_cast<size_t>(i)]).empty()) {
            return i;
        }
    }
    throw std::runtime_error("Unexpected end of COMSOL mphtxt file.");
}

static int findLineContaining(const std::vector<std::string>& lines, const std::string& text, int start = 0)
{
    for (int i = start; i < static_cast<int>(lines.size()); ++i) {
        if (lines[static_cast<size_t>(i)].find(text) != std::string::npos) {
            return i;
        }
    }
    throw std::runtime_error("Cannot find section: " + text);
}

static void parseElementSection(const std::vector<std::string>& lines,
                                const std::string& typeName,
                                int verticesPerElementExpected,
                                std::vector<std::vector<int>>& elements,
                                std::vector<int>& entityIds)
{
    const int typeLine = findLineContaining(lines, typeName + " # type name");
    const int nvLine = nextNumericLine(lines, typeLine + 1);
    const int countLine = nextNumericLine(lines, nvLine + 1);
    const int verticesPerElement = firstInt(lines[static_cast<size_t>(nvLine)]);
    const int elementCount = firstInt(lines[static_cast<size_t>(countLine)]);
    if (verticesPerElement != verticesPerElementExpected) {
        throw std::runtime_error("Unexpected vertex count in " + typeName + " section.");
    }

    const int elementsHeader = findLineContaining(lines, "# Elements", countLine + 1);
    elements.clear();
    elements.reserve(static_cast<size_t>(elementCount));
    int lineIndex = elementsHeader + 1;
    for (int i = 0; i < elementCount; ++i) {
        lineIndex = nextNumericLine(lines, lineIndex);
        std::vector<int> values = parseInts(lines[static_cast<size_t>(lineIndex)]);
        if (static_cast<int>(values.size()) < verticesPerElement) {
            throw std::runtime_error("Malformed element connectivity in " + typeName + " section.");
        }
        values.resize(static_cast<size_t>(verticesPerElement));
        elements.push_back(values);
        ++lineIndex;
    }

    const int entityHeader = findLineContaining(lines, "# Geometric entity indices", lineIndex);
    entityIds.clear();
    entityIds.reserve(static_cast<size_t>(elementCount));
    lineIndex = entityHeader + 1;
    for (int i = 0; i < elementCount; ++i) {
        lineIndex = nextNumericLine(lines, lineIndex);
        entityIds.push_back(firstInt(lines[static_cast<size_t>(lineIndex)]));
        ++lineIndex;
    }
}

static ImportedMesh readComsolMphtxt(const std::filesystem::path& path, double coordinateScale)
{
    const std::vector<std::string> lines = readLines(path);
    const int vertexCountLine = findLineContaining(lines, "# number of mesh vertices");
    const int vertexCount = firstInt(lines[static_cast<size_t>(vertexCountLine)]);
    const int coordHeader = findLineContaining(lines, "# Mesh vertex coordinates", vertexCountLine + 1);

    ImportedMesh mesh;
    mesh.vertices.reserve(static_cast<size_t>(vertexCount));
    int lineIndex = coordHeader + 1;
    for (int i = 0; i < vertexCount; ++i) {
        lineIndex = nextNumericLine(lines, lineIndex);
        const std::vector<double> values = parseDoubles(lines[static_cast<size_t>(lineIndex)]);
        if (values.size() < 3) {
            throw std::runtime_error("Malformed vertex coordinate line in " + path.string());
        }
        mesh.vertices.push_back({coordinateScale * values[0], coordinateScale * values[1], coordinateScale * values[2]});
        ++lineIndex;
    }

    std::vector<std::vector<int>> triElements;
    std::vector<int> triEntities;
    parseElementSection(lines, "tri", 3, triElements, triEntities);
    mesh.triangles.reserve(triElements.size());
    for (const std::vector<int>& tri : triElements) {
        mesh.triangles.push_back({tri[0], tri[1], tri[2]});
    }
    mesh.triangleEntities = std::move(triEntities);

    std::vector<std::vector<int>> tetElements;
    std::vector<int> tetEntities;
    parseElementSection(lines, "tet", 4, tetElements, tetEntities);
    mesh.tets.reserve(tetElements.size());
    for (const std::vector<int>& tet : tetElements) {
        mesh.tets.push_back({tet[0], tet[1], tet[2], tet[3]});
    }
    mesh.tetEntities = std::move(tetEntities);
    return mesh;
}

static std::filesystem::path findProjectRoot()
{
    std::filesystem::path p = std::filesystem::current_path();
    for (int i = 0; i < 6; ++i) {
        if (std::filesystem::exists(p / "case_config.txt")
            || std::filesystem::exists(p / "CMakeLists.txt")) {
            return p;
        }
        if (!p.has_parent_path() || p == p.parent_path()) {
            break;
        }
        p = p.parent_path();
    }
    throw std::runtime_error("Cannot locate project root from current directory.");
}

static std::uint64_t edgeKey(int a, int b)
{
    if (a > b) {
        std::swap(a, b);
    }
    return (static_cast<std::uint64_t>(static_cast<unsigned int>(a)) << 32)
        | static_cast<std::uint64_t>(static_cast<unsigned int>(b));
}

static std::array<int, 3> sortedFaceKey(std::array<int, 3> f)
{
    std::sort(f.begin(), f.end());
    return f;
}

static int localEdgeDof(int a, int b)
{
    if (a > b) {
        std::swap(a, b);
    }
    if (a == 0 && b == 1) return 4;
    if (a == 0 && b == 2) return 5;
    if (a == 0 && b == 3) return 6;
    if (a == 1 && b == 2) return 7;
    if (a == 1 && b == 3) return 8;
    if (a == 2 && b == 3) return 9;
    throw std::runtime_error("Invalid local edge.");
}

struct DomainBuildContext {
    std::vector<int> vertexToNode;
    std::unordered_map<std::uint64_t, int> edgeToNode;
    std::map<std::array<int, 3>, std::pair<int, std::array<int, 3>>> faceOwner;
};

static int addMidpointNode(Mesh& mesh,
                           DomainBuildContext& context,
                           int sourceA,
                           int sourceB,
                           int nodeA,
                           int nodeB,
                           int subdomain)
{
    const std::uint64_t key = edgeKey(sourceA, sourceB);
    const auto found = context.edgeToNode.find(key);
    if (found != context.edgeToNode.end()) {
        return found->second;
    }
    Node node;
    node.p = 0.5 * (mesh.nodes[static_cast<size_t>(nodeA)].p + mesh.nodes[static_cast<size_t>(nodeB)].p);
    node.subdomain = subdomain;
    node.sourceVertex = -1;
    const int id = static_cast<int>(mesh.nodes.size());
    mesh.nodes.push_back(node);
    context.edgeToNode.emplace(key, id);
    return id;
}

static void completeP2Dofs(Mesh& mesh, DomainBuildContext& context, Tet& tet)
{
    for (int i = 0; i < 4; ++i) {
        tet.dof[static_cast<size_t>(i)] = tet.v[static_cast<size_t>(i)];
    }
    const std::array<std::array<int, 2>, 6> edges{{
        {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}, {{1, 3}}, {{2, 3}}
    }};
    for (int e = 0; e < 6; ++e) {
        const int a = edges[static_cast<size_t>(e)][0];
        const int b = edges[static_cast<size_t>(e)][1];
        tet.dof[static_cast<size_t>(4 + e)] =
            addMidpointNode(mesh, context,
                            tet.source[static_cast<size_t>(a)], tet.source[static_cast<size_t>(b)],
                            tet.v[static_cast<size_t>(a)], tet.v[static_cast<size_t>(b)],
                            tet.subdomain);
    }
}

static void updateBounds(Vec3& lo, Vec3& hi, const Vec3& p)
{
    lo.x = std::min(lo.x, p.x);
    lo.y = std::min(lo.y, p.y);
    lo.z = std::min(lo.z, p.z);
    hi.x = std::max(hi.x, p.x);
    hi.y = std::max(hi.y, p.y);
    hi.z = std::max(hi.z, p.z);
}

static void appendImportedDomain(Mesh& mesh, const ImportedMesh& imported, int subdomain)
{
    DomainBuildContext context;
    context.vertexToNode.resize(imported.vertices.size());
    if (mesh.subdomainMin.size() <= static_cast<size_t>(subdomain)) {
        mesh.subdomainMin.resize(static_cast<size_t>(subdomain + 1));
        mesh.subdomainMax.resize(static_cast<size_t>(subdomain + 1));
    }
    Vec3 lo{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    Vec3 hi{-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};

    for (int i = 0; i < static_cast<int>(imported.vertices.size()); ++i) {
        Node node;
        node.p = imported.vertices[static_cast<size_t>(i)];
        node.subdomain = subdomain;
        node.sourceVertex = i;
        context.vertexToNode[static_cast<size_t>(i)] = static_cast<int>(mesh.nodes.size());
        mesh.nodes.push_back(node);
        updateBounds(lo, hi, node.p);
    }
    mesh.subdomainMin[static_cast<size_t>(subdomain)] = lo;
    mesh.subdomainMax[static_cast<size_t>(subdomain)] = hi;

    const std::array<std::array<int, 3>, 4> facePattern{{
        {{1, 2, 3}}, {{0, 3, 2}}, {{0, 1, 3}}, {{0, 2, 1}}
    }};

    for (int e = 0; e < static_cast<int>(imported.tets.size()); ++e) {
        Tet tet;
        tet.subdomain = subdomain;
        tet.domainEntity = imported.tetEntities[static_cast<size_t>(e)];
        for (int i = 0; i < 4; ++i) {
            tet.source[static_cast<size_t>(i)] = imported.tets[static_cast<size_t>(e)][static_cast<size_t>(i)];
            tet.v[static_cast<size_t>(i)] = context.vertexToNode[static_cast<size_t>(tet.source[static_cast<size_t>(i)])];
        }
        completeP2Dofs(mesh, context, tet);
        const int tetId = static_cast<int>(mesh.tets.size());
        mesh.tets.push_back(tet);

        for (const auto& face : facePattern) {
            std::array<int, 3> sourceFace{{
                tet.source[static_cast<size_t>(face[0])],
                tet.source[static_cast<size_t>(face[1])],
                tet.source[static_cast<size_t>(face[2])]
            }};
            context.faceOwner[sortedFaceKey(sourceFace)] = {tetId, face};
        }
    }

    for (int i = 0; i < static_cast<int>(imported.triangles.size()); ++i) {
        const std::array<int, 3> tri = imported.triangles[static_cast<size_t>(i)];
        const auto found = context.faceOwner.find(sortedFaceKey(tri));
        if (found == context.faceOwner.end()) {
            continue;
        }
        BoundaryFace face;
        face.tet = found->second.first;
        face.local = found->second.second;
        face.subdomain = subdomain;
        face.boundaryEntity = imported.triangleEntities[static_cast<size_t>(i)];
        for (int k = 0; k < 3; ++k) {
            face.points[static_cast<size_t>(k)] =
                mesh.nodes[static_cast<size_t>(context.vertexToNode[static_cast<size_t>(tri[static_cast<size_t>(k)])])].p;
        }
        Vec3 rawNormal = cross(face.points[1] - face.points[0], face.points[2] - face.points[0]);
        Vec3 centroid = (face.points[0] + face.points[1] + face.points[2]) / 3.0;
        int oppositeLocal = 0;
        for (; oppositeLocal < 4; ++oppositeLocal) {
            if (oppositeLocal != face.local[0]
                && oppositeLocal != face.local[1]
                && oppositeLocal != face.local[2]) {
                break;
            }
        }
        if (oppositeLocal < 4) {
            const Tet& ownerTet = mesh.tets[static_cast<size_t>(face.tet)];
            const Vec3 opposite = mesh.nodes[static_cast<size_t>(ownerTet.v[static_cast<size_t>(oppositeLocal)])].p;
            if (dot(rawNormal, centroid - opposite) < 0.0) {
                rawNormal = -1.0 * rawNormal;
            }
        }
        face.normal = normalized(rawNormal);
        mesh.boundaryFaces.push_back(face);
    }
}

static void markDirichletFaceDofs(Mesh& mesh, const BoundaryFace& face, double value)
{
    Tet& tet = mesh.tets[static_cast<size_t>(face.tet)];
    for (int i = 0; i < 3; ++i) {
        Node& node = mesh.nodes[static_cast<size_t>(tet.dof[static_cast<size_t>(face.local[static_cast<size_t>(i)])])];
        node.dirichlet = true;
        node.dirichletValue = value;
    }
    for (int i = 0; i < 3; ++i) {
        const int a = face.local[static_cast<size_t>(i)];
        const int b = face.local[static_cast<size_t>((i + 1) % 3)];
        Node& node = mesh.nodes[static_cast<size_t>(tet.dof[static_cast<size_t>(localEdgeDof(a, b))])];
        node.dirichlet = true;
        node.dirichletValue = value;
    }
}

static bool boundaryMatches(const BoundaryFace& face, const BoundaryCondition& condition)
{
    return condition.boundaryEntity == face.boundaryEntity
        && (condition.subdomain < 0 || condition.subdomain == face.subdomain);
}

static bool boundaryMatches(const BoundaryFace& face, const ConvectionCondition& condition)
{
    return condition.boundaryEntity == face.boundaryEntity
        && (condition.subdomain < 0 || condition.subdomain == face.subdomain);
}

static bool tetMatchesHeatSource(const Tet& tet, const HeatSource& source)
{
    return source.domainEntity == tet.domainEntity
        && (source.subdomain < 0 || source.subdomain == tet.subdomain);
}

static bool isDirichletBoundary(const BoundaryFace& face, const CaseConfig& config)
{
    for (const BoundaryCondition& condition : config.dirichletConditions) {
        if (boundaryMatches(face, condition)) {
            return true;
        }
    }
    return false;
}

static bool isConvectionBoundary(const BoundaryFace& face, const CaseConfig& config)
{
    for (const ConvectionCondition& condition : config.convectionConditions) {
        if (boundaryMatches(face, condition)) {
            return true;
        }
    }
    return false;
}

static double triangleArea(const std::array<Vec3, 3>& tri)
{
    return 0.5 * norm(cross(tri[1] - tri[0], tri[2] - tri[0]));
}

static double boundaryFaceArea(const BoundaryFace& face)
{
    return triangleArea(face.points);
}

static double integrationArea(const std::vector<std::array<Vec3, 3>>& triangles)
{
    double area = 0.0;
    for (const auto& tri : triangles) {
        area += triangleArea(tri);
    }
    return area;
}

struct Projector {
    int dropAxis = 2;
    double planeCoordinate = 0.0;

    Vec2 to2D(const Vec3& p) const
    {
        if (dropAxis == 0) return {p.y, p.z};
        if (dropAxis == 1) return {p.x, p.z};
        return {p.x, p.y};
    }

    Vec3 to3D(const Vec2& p) const
    {
        if (dropAxis == 0) return {planeCoordinate, p.a, p.b};
        if (dropAxis == 1) return {p.a, planeCoordinate, p.b};
        return {p.a, p.b, planeCoordinate};
    }
};

static double cross2D(const Vec2& a, const Vec2& b, const Vec2& c)
{
    return (b.a - a.a) * (c.b - a.b) - (b.b - a.b) * (c.a - a.a);
}

static double polygonArea(const std::vector<Vec2>& poly)
{
    double area = 0.0;
    for (size_t i = 0; i < poly.size(); ++i) {
        const Vec2& a = poly[i];
        const Vec2& b = poly[(i + 1) % poly.size()];
        area += a.a * b.b - a.b * b.a;
    }
    return 0.5 * area;
}

static Vec2 lineIntersection(const Vec2& p0, const Vec2& p1, const Vec2& q0, const Vec2& q1)
{
    const double a1 = p1.b - p0.b;
    const double b1 = p0.a - p1.a;
    const double c1 = a1 * p0.a + b1 * p0.b;
    const double a2 = q1.b - q0.b;
    const double b2 = q0.a - q1.a;
    const double c2 = a2 * q0.a + b2 * q0.b;
    const double det = a1 * b2 - a2 * b1;
    if (std::abs(det) < 1.0e-24) {
        return p1;
    }
    return {(b2 * c1 - b1 * c2) / det, (a1 * c2 - a2 * c1) / det};
}

static std::vector<Vec2> clipByEdge(const std::vector<Vec2>& subject, const Vec2& a, const Vec2& b, double orientation)
{
    std::vector<Vec2> output;
    if (subject.empty()) {
        return output;
    }
    const auto inside = [&](const Vec2& p) {
        return orientation * cross2D(a, b, p) >= -1.0e-15;
    };
    Vec2 previous = subject.back();
    bool previousInside = inside(previous);
    for (const Vec2& current : subject) {
        const bool currentInside = inside(current);
        if (currentInside != previousInside) {
            output.push_back(lineIntersection(previous, current, a, b));
        }
        if (currentInside) {
            output.push_back(current);
        }
        previous = current;
        previousInside = currentInside;
    }
    return output;
}

static bool bboxOverlap(const std::array<Vec3, 3>& a, const std::array<Vec3, 3>& b);

static std::vector<Vec2> clipByEdgeScaled(const std::vector<Vec2>& subject,
                                          const Vec2& a,
                                          const Vec2& b,
                                          double orientation,
                                          double tolerance)
{
    std::vector<Vec2> output;
    if (subject.empty()) {
        return output;
    }
    const auto inside = [&](const Vec2& p) {
        return orientation * cross2D(a, b, p) >= -tolerance;
    };
    Vec2 previous = subject.back();
    bool previousInside = inside(previous);
    for (const Vec2& current : subject) {
        const bool currentInside = inside(current);
        if (currentInside != previousInside) {
            output.push_back(lineIntersection(previous, current, a, b));
        }
        if (currentInside) {
            output.push_back(current);
        }
        previous = current;
        previousInside = currentInside;
    }
    std::vector<Vec2> cleaned;
    cleaned.reserve(output.size());
    const double mergeTol2 = std::max(1.0e-60, tolerance);
    for (const Vec2& p : output) {
        if (!cleaned.empty()) {
            const Vec2 d = p - cleaned.back();
            if (dot2D(d, d) <= mergeTol2) {
                continue;
            }
        }
        cleaned.push_back(p);
    }
    if (cleaned.size() > 1) {
        const Vec2 d = cleaned.front() - cleaned.back();
        if (dot2D(d, d) <= mergeTol2) {
            cleaned.pop_back();
        }
    }
    return cleaned;
}

struct InterfacePlane2D {
    Vec3 origin{};
    Vec3 e1{};
    Vec3 e2{};
    Vec3 normal{};
};

static InterfacePlane2D makeInterfacePlane2D(const std::array<Vec3, 3>& left, const Vec3& leftNormal)
{
    InterfacePlane2D plane;
    plane.origin = left[0];
    plane.normal = norm(leftNormal) > 1.0e-30
        ? normalized(leftNormal)
        : normalized(cross(left[1] - left[0], left[2] - left[0]));
    plane.e1 = normalized(left[1] - left[0]);
    if (norm(plane.e1) < 1.0e-30) {
        plane.e1 = normalized(left[2] - left[0]);
    }
    plane.e2 = normalized(cross(plane.normal, plane.e1));
    if (norm(plane.e2) < 1.0e-30) {
        plane.e2 = normalized(cross(normalized(cross(left[1] - left[0], left[2] - left[0])), plane.e1));
    }
    return plane;
}

static Vec2 projectToInterfacePlane(const InterfacePlane2D& plane, const Vec3& p)
{
    const Vec3 d = p - plane.origin;
    return {dot(d, plane.e1), dot(d, plane.e2)};
}

static Vec3 liftFromInterfacePlane(const InterfacePlane2D& plane, const Vec2& p)
{
    return plane.origin + p.a * plane.e1 + p.b * plane.e2;
}

static double maxTriangleEdgeLength2D(const std::array<Vec2, 3>& tri)
{
    double maxLen = 0.0;
    for (int i = 0; i < 3; ++i) {
        const Vec2 d = tri[static_cast<size_t>((i + 1) % 3)] - tri[static_cast<size_t>(i)];
        maxLen = std::max(maxLen, std::sqrt(dot2D(d, d)));
    }
    return maxLen;
}

struct TriangleOverlapResult {
    std::vector<std::array<Vec3, 3>> triangles;
    int polygonVertices = 0;
    double area = 0.0;
};

static TriangleOverlapResult intersectTrianglesOnInterfacePlane(const std::array<Vec3, 3>& left,
                                                                const std::array<Vec3, 3>& right,
                                                                const Vec3& leftNormal,
                                                                const Vec3& rightNormal)
{
    TriangleOverlapResult result;
    if (!bboxOverlap(left, right)) {
        return result;
    }
    const Vec3 nLeft = norm(leftNormal) > 1.0e-30
        ? normalized(leftNormal)
        : normalized(cross(left[1] - left[0], left[2] - left[0]));
    const Vec3 nRight = norm(rightNormal) > 1.0e-30
        ? normalized(rightNormal)
        : normalized(cross(right[1] - right[0], right[2] - right[0]));
    if (std::abs(dot(nLeft, nRight)) < 1.0 - 1.0e-8) {
        return result;
    }

    const InterfacePlane2D plane = makeInterfacePlane2D(left, nLeft);
    const std::array<Vec2, 3> left2{{
        projectToInterfacePlane(plane, left[0]),
        projectToInterfacePlane(plane, left[1]),
        projectToInterfacePlane(plane, left[2])
    }};
    const std::array<Vec2, 3> right2{{
        projectToInterfacePlane(plane, right[0]),
        projectToInterfacePlane(plane, right[1]),
        projectToInterfacePlane(plane, right[2])
    }};
    const double maxEdge = std::max(maxTriangleEdgeLength2D(left2), maxTriangleEdgeLength2D(right2));
    const double clipTolerance = std::max(1.0e-30, maxEdge * maxEdge * 1.0e-10);
    const double areaTolerance = std::max(1.0e-30, maxEdge * maxEdge * 1.0e-12);
    std::vector<Vec2> polygon{left2[0], left2[1], left2[2]};
    const double clipOrientation = polygonArea({right2[0], right2[1], right2[2]}) >= 0.0 ? 1.0 : -1.0;
    for (int i = 0; i < 3; ++i) {
        polygon = clipByEdgeScaled(polygon,
                                   right2[static_cast<size_t>(i)],
                                   right2[static_cast<size_t>((i + 1) % 3)],
                                   clipOrientation,
                                   clipTolerance);
    }
    double signedArea = polygonArea(polygon);
    if (polygon.size() < 3 || std::abs(signedArea) <= areaTolerance) {
        return result;
    }
    if (signedArea < 0.0) {
        std::reverse(polygon.begin(), polygon.end());
        signedArea = -signedArea;
    }
    result.polygonVertices = static_cast<int>(polygon.size());
    result.area = signedArea;
    result.triangles.reserve(polygon.size() - 2);
    for (size_t i = 1; i + 1 < polygon.size(); ++i) {
        result.triangles.push_back({{
            liftFromInterfacePlane(plane, polygon[0]),
            liftFromInterfacePlane(plane, polygon[i]),
            liftFromInterfacePlane(plane, polygon[i + 1])
        }});
    }
    return result;
}

static bool bboxOverlap(const std::array<Vec3, 3>& a, const std::array<Vec3, 3>& b)
{
    Vec3 amin{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    Vec3 amax{-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};
    Vec3 bmin = amin;
    Vec3 bmax = amax;
    for (const Vec3& p : a) updateBounds(amin, amax, p);
    for (const Vec3& p : b) updateBounds(bmin, bmax, p);
    const double eps = 1.0e-12;
    return amin.x <= bmax.x + eps && amax.x + eps >= bmin.x
        && amin.y <= bmax.y + eps && amax.y + eps >= bmin.y
        && amin.z <= bmax.z + eps && amax.z + eps >= bmin.z;
}

static Projector makeProjector(const std::array<Vec3, 3>& tri)
{
    const Vec3 n = cross(tri[1] - tri[0], tri[2] - tri[0]);
    const double ax = std::abs(n.x);
    const double ay = std::abs(n.y);
    const double az = std::abs(n.z);
    Projector p;
    if (ax >= ay && ax >= az) {
        p.dropAxis = 0;
        p.planeCoordinate = (tri[0].x + tri[1].x + tri[2].x) / 3.0;
    } else if (ay >= ax && ay >= az) {
        p.dropAxis = 1;
        p.planeCoordinate = (tri[0].y + tri[1].y + tri[2].y) / 3.0;
    } else {
        p.dropAxis = 2;
        p.planeCoordinate = (tri[0].z + tri[1].z + tri[2].z) / 3.0;
    }
    return p;
}

static std::vector<std::array<Vec3, 3>> intersectTrianglesAsFan(const std::array<Vec3, 3>& left,
                                                                const std::array<Vec3, 3>& right)
{
    if (!bboxOverlap(left, right)) {
        return {};
    }
    const Vec3 leftNormalRaw = cross(left[1] - left[0], left[2] - left[0]);
    const Vec3 rightNormalRaw = cross(right[1] - right[0], right[2] - right[0]);
    const double leftNormalNorm = norm(leftNormalRaw);
    const double rightNormalNorm = norm(rightNormalRaw);
    if (leftNormalNorm < 1.0e-30 || rightNormalNorm < 1.0e-30) {
        return {};
    }
    const Vec3 leftNormal = leftNormalRaw / leftNormalNorm;
    const Vec3 rightNormal = rightNormalRaw / rightNormalNorm;
    if (std::abs(dot(leftNormal, rightNormal)) < 1.0 - 1.0e-10) {
        return {};
    }
    for (const Vec3& p : right) {
        if (std::abs(dot(leftNormal, p - left[0])) > 1.0e-10) {
            return {};
        }
    }
    const Projector projector = makeProjector(left);
    std::vector<Vec2> polygon{projector.to2D(left[0]), projector.to2D(left[1]), projector.to2D(left[2])};
    std::array<Vec2, 3> clip{{projector.to2D(right[0]), projector.to2D(right[1]), projector.to2D(right[2])}};
    const double orientation = polygonArea({clip[0], clip[1], clip[2]}) >= 0.0 ? 1.0 : -1.0;
    for (int i = 0; i < 3; ++i) {
        polygon = clipByEdge(polygon, clip[static_cast<size_t>(i)], clip[static_cast<size_t>((i + 1) % 3)], orientation);
    }
    if (polygon.size() < 3 || std::abs(polygonArea(polygon)) < 1.0e-20) {
        return {};
    }
    std::vector<std::array<Vec3, 3>> triangles;
    for (size_t i = 1; i + 1 < polygon.size(); ++i) {
        triangles.push_back({{
            projector.to3D(polygon[0]),
            projector.to3D(polygon[i]),
            projector.to3D(polygon[i + 1])
        }});
    }
    return triangles;
}

static Mesh loadComsolDomainDecompositionMesh(const CaseConfig& config)
{
    Mesh mesh;
    mesh.subdomainMin.resize(config.domains.size());
    mesh.subdomainMax.resize(config.domains.size());
    for (int subdomain = 0; subdomain < static_cast<int>(config.domains.size()); ++subdomain) {
        ImportedMesh imported = readComsolMphtxt(config.domains[static_cast<size_t>(subdomain)].meshPath,
                                                 config.coordinateScale);
        const Vec3 offset = config.domains[static_cast<std::size_t>(subdomain)].translationMeters;
        for (Vec3& vertex : imported.vertices) {
            vertex = vertex + offset;
        }
        appendImportedDomain(mesh, imported, subdomain);
    }

    for (const BoundaryFace& face : mesh.boundaryFaces) {
        for (const BoundaryCondition& condition : config.dirichletConditions) {
            if (boundaryMatches(face, condition)) {
                markDirichletFaceDofs(mesh, face, condition.temperature);
            }
        }
    }

    std::vector<InterfaceConfig> interfaces = config.interfaces;
    if (interfaces.empty() && config.autoInterfaces) {
        std::map<int, std::set<int>> boundaryEntitiesBySubdomain;
        for (const BoundaryFace& face : mesh.boundaryFaces) {
            if (!isDirichletBoundary(face, config)) {
                boundaryEntitiesBySubdomain[face.subdomain].insert(face.boundaryEntity);
            }
        }
        std::vector<int> subdomains;
        for (const auto& entry : boundaryEntitiesBySubdomain) {
            subdomains.push_back(entry.first);
        }
        for (size_t i = 0; i < subdomains.size(); ++i) {
            for (size_t j = i + 1; j < subdomains.size(); ++j) {
                const std::set<int>& leftSet = boundaryEntitiesBySubdomain[subdomains[i]];
                const std::set<int>& rightSet = boundaryEntitiesBySubdomain[subdomains[j]];
                interfaces.push_back({
                    subdomains[i],
                    subdomains[j],
                    std::vector<int>(leftSet.begin(), leftSet.end()),
                    std::vector<int>(rightSet.begin(), rightSet.end())
                });
            }
        }
    }

    for (const InterfaceConfig& interfaceConfig : interfaces) {
        std::vector<const BoundaryFace*> leftFaces;
        std::vector<const BoundaryFace*> rightFaces;
        InterfaceBuildSummary summary;
        summary.leftSubdomain = interfaceConfig.leftSubdomain;
        summary.rightSubdomain = interfaceConfig.rightSubdomain;
        summary.leftBoundaryEntityCount = static_cast<int>(interfaceConfig.leftBoundaryEntities.size());
        summary.rightBoundaryEntityCount = static_cast<int>(interfaceConfig.rightBoundaryEntities.size());
        for (const BoundaryFace& face : mesh.boundaryFaces) {
            if (isDirichletBoundary(face, config)) {
                continue;
            }
            if (face.subdomain == interfaceConfig.leftSubdomain) {
                if (containsEntity(interfaceConfig.leftBoundaryEntities, face.boundaryEntity)) {
                    leftFaces.push_back(&face);
                    summary.leftArea += boundaryFaceArea(face);
                }
            } else if (face.subdomain == interfaceConfig.rightSubdomain) {
                if (containsEntity(interfaceConfig.rightBoundaryEntities, face.boundaryEntity)) {
                    rightFaces.push_back(&face);
                    summary.rightArea += boundaryFaceArea(face);
                }
            }
        }
        summary.leftFaceCount = static_cast<int>(leftFaces.size());
        summary.rightFaceCount = static_cast<int>(rightFaces.size());

        for (const BoundaryFace* left : leftFaces) {
            for (const BoundaryFace* right : rightFaces) {
                TriangleOverlapResult overlap = intersectTrianglesOnInterfacePlane(left->points,
                                                                                   right->points,
                                                                                   left->normal,
                                                                                   right->normal);
                if (overlap.triangles.empty()) {
                    continue;
                }
                const double normalDot = dot(left->normal, right->normal);
                summary.normalDotMin = std::min(summary.normalDotMin, normalDot);
                summary.normalDotMax = std::max(summary.normalDotMax, normalDot);
                summary.normalDotSum += normalDot;
                summary.matchedOverlapArea += overlap.area;
                ++summary.facePairCount;
                summary.integrationTriangleCount += static_cast<int>(overlap.triangles.size());
                InterfaceFace pair;
                pair.leftTet = left->tet;
                pair.rightTet = right->tet;
                pair.leftFaceId = static_cast<int>(left - mesh.boundaryFaces.data());
                pair.rightFaceId = static_cast<int>(right - mesh.boundaryFaces.data());
                pair.leftLocal = left->local;
                pair.rightLocal = right->local;
                pair.leftNormal = left->normal;
                pair.rightNormal = right->normal;
                pair.leftBoundaryEntity = left->boundaryEntity;
                pair.rightBoundaryEntity = right->boundaryEntity;
                pair.integrationTriangles = std::move(overlap.triangles);
                pair.overlapPolygonVertices = overlap.polygonVertices;
                pair.overlapArea = overlap.area;
                mesh.interfaceFaces.push_back(std::move(pair));
            }
        }
        if (summary.facePairCount == 0) {
            summary.normalDotMin = 0.0;
            summary.normalDotMax = 0.0;
        }
        mesh.interfaceSummaries.push_back(summary);
    }

    if (!config.explicitInterfaceFacePairs.empty()) {
        InterfaceBuildSummary summary;
        summary.leftSubdomain = 0;
        summary.rightSubdomain = 1;
        std::set<int> leftFaceIds;
        std::set<int> rightFaceIds;
        std::set<int> leftBoundaryEntities;
        std::set<int> rightBoundaryEntities;
        for (const auto& facePair : config.explicitInterfaceFacePairs) {
            const int leftFaceId = facePair.first;
            const int rightFaceId = facePair.second;
            if (leftFaceId < 0 || leftFaceId >= static_cast<int>(mesh.boundaryFaces.size())
                || rightFaceId < 0 || rightFaceId >= static_cast<int>(mesh.boundaryFaces.size())) {
                continue;
            }
            const BoundaryFace& left = mesh.boundaryFaces[static_cast<size_t>(leftFaceId)];
            const BoundaryFace& right = mesh.boundaryFaces[static_cast<size_t>(rightFaceId)];
            if (isDirichletBoundary(left, config) || isDirichletBoundary(right, config)) {
                continue;
            }
            TriangleOverlapResult overlap = intersectTrianglesOnInterfacePlane(left.points,
                                                                               right.points,
                                                                               left.normal,
                                                                               right.normal);
            if (overlap.triangles.empty()) {
                continue;
            }
            const double normalDot = dot(left.normal, right.normal);
            summary.normalDotMin = std::min(summary.normalDotMin, normalDot);
            summary.normalDotMax = std::max(summary.normalDotMax, normalDot);
            summary.normalDotSum += normalDot;
            summary.matchedOverlapArea += overlap.area;
            ++summary.facePairCount;
            summary.integrationTriangleCount += static_cast<int>(overlap.triangles.size());
            leftFaceIds.insert(leftFaceId);
            rightFaceIds.insert(rightFaceId);
            leftBoundaryEntities.insert(left.boundaryEntity);
            rightBoundaryEntities.insert(right.boundaryEntity);
            InterfaceFace pair;
            pair.leftTet = left.tet;
            pair.rightTet = right.tet;
            pair.leftFaceId = leftFaceId;
            pair.rightFaceId = rightFaceId;
            pair.leftLocal = left.local;
            pair.rightLocal = right.local;
            pair.leftNormal = left.normal;
            pair.rightNormal = right.normal;
            pair.leftBoundaryEntity = left.boundaryEntity;
            pair.rightBoundaryEntity = right.boundaryEntity;
            pair.integrationTriangles = std::move(overlap.triangles);
            pair.overlapPolygonVertices = overlap.polygonVertices;
            pair.overlapArea = overlap.area;
            mesh.interfaceFaces.push_back(std::move(pair));
        }
        summary.leftFaceCount = static_cast<int>(leftFaceIds.size());
        summary.rightFaceCount = static_cast<int>(rightFaceIds.size());
        summary.leftBoundaryEntityCount = static_cast<int>(leftBoundaryEntities.size());
        summary.rightBoundaryEntityCount = static_cast<int>(rightBoundaryEntities.size());
        for (int id : leftFaceIds) {
            summary.leftArea += boundaryFaceArea(mesh.boundaryFaces[static_cast<size_t>(id)]);
        }
        for (int id : rightFaceIds) {
            summary.rightArea += boundaryFaceArea(mesh.boundaryFaces[static_cast<size_t>(id)]);
        }
        if (summary.facePairCount == 0) {
            summary.normalDotMin = 0.0;
            summary.normalDotMax = 0.0;
        }
        mesh.interfaceSummaries.push_back(summary);
    }

    return mesh;
}
