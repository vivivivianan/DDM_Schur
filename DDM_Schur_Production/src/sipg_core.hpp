#pragma once

// Core mesh/material/config types, sparse triplet storage, parallel utilities,
// memory counters, and platform/MKL includes. This legacy header is broad, but
// it is the common numerical vocabulary used by assembly and validation.

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <psapi.h>
#endif

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Vec2 {
    double a = 0.0;
    double b = 0.0;
};

static Vec3 operator+(const Vec3& u, const Vec3& v) { return {u.x + v.x, u.y + v.y, u.z + v.z}; }
static Vec3 operator-(const Vec3& u, const Vec3& v) { return {u.x - v.x, u.y - v.y, u.z - v.z}; }
static Vec3 operator*(double s, const Vec3& u) { return {s * u.x, s * u.y, s * u.z}; }
static Vec3 operator*(const Vec3& u, double s) { return s * u; }
static Vec3 operator/(const Vec3& u, double s) { return {u.x / s, u.y / s, u.z / s}; }
[[maybe_unused]] static Vec2 operator+(const Vec2& u, const Vec2& v) { return {u.a + v.a, u.b + v.b}; }
static Vec2 operator-(const Vec2& u, const Vec2& v) { return {u.a - v.a, u.b - v.b}; }
[[maybe_unused]] static Vec2 operator*(double s, const Vec2& u) { return {s * u.a, s * u.b}; }
static double dot(const Vec3& u, const Vec3& v) { return u.x * v.x + u.y * v.y + u.z * v.z; }
static double dot2D(const Vec2& u, const Vec2& v) { return u.a * v.a + u.b * v.b; }
static Vec3 cross(const Vec3& u, const Vec3& v)
{
    return {u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
}
static double norm(const Vec3& u) { return std::sqrt(dot(u, u)); }
static Vec3 normalized(const Vec3& u)
{
    const double n = norm(u);
    if (n < 1.0e-30) {
        return {0.0, 0.0, 1.0};
    }
    return u / n;
}

struct Node {
    Vec3 p;
    int subdomain = 0;
    int sourceVertex = -1;
    bool dirichlet = false;
    double dirichletValue = 293.15;
};

struct Tet {
    std::array<int, 4> v{};
    std::array<int, 4> source{};
    std::array<int, 10> dof{};
    int subdomain = 0;
    int domainEntity = -1;
};

struct BoundaryFace {
    int tet = -1;
    int subdomain = 0;
    int boundaryEntity = -1;
    std::array<int, 3> local{};
    std::array<Vec3, 3> points{};
    Vec3 normal{};
};

struct InterfaceFace {
    int leftTet = -1;
    int rightTet = -1;
    int leftFaceId = -1;
    int rightFaceId = -1;
    int leftBoundaryEntity = -1;
    int rightBoundaryEntity = -1;
    std::array<int, 3> leftLocal{};
    std::array<int, 3> rightLocal{};
    Vec3 leftNormal{};
    Vec3 rightNormal{};
    std::vector<std::array<Vec3, 3>> integrationTriangles;
    int overlapPolygonVertices = 0;
    double overlapArea = 0.0;
};

struct InterfaceBuildSummary {
    int leftSubdomain = -1;
    int rightSubdomain = -1;
    int leftBoundaryEntityCount = 0;
    int rightBoundaryEntityCount = 0;
    int leftFaceCount = 0;
    int rightFaceCount = 0;
    double leftArea = 0.0;
    double rightArea = 0.0;
    double matchedOverlapArea = 0.0;
    int facePairCount = 0;
    int integrationTriangleCount = 0;
    double normalDotMin = std::numeric_limits<double>::max();
    double normalDotMax = -std::numeric_limits<double>::max();
    double normalDotSum = 0.0;
};

struct Mesh {
    std::vector<Node> nodes;
    std::vector<Tet> tets;
    std::vector<BoundaryFace> boundaryFaces;
    std::vector<InterfaceFace> interfaceFaces;
    std::vector<InterfaceBuildSummary> interfaceSummaries;
    std::vector<Vec3> subdomainMin;
    std::vector<Vec3> subdomainMax;
};

struct Material {
    std::string name = "Copper";
    double conductivity = 400.0;
    double conductivityX = 400.0;
    double conductivityY = 400.0;
    double conductivityZ = 400.0;
    double density = 8960.0;
    double heatCapacity = 385.0;
};

struct DomainConfig {
    std::filesystem::path meshPath;
    Material material;
    std::map<int, Material> materialsByDomainEntity;
    Vec3 translationMeters{};
};

struct BoundaryCondition {
    int subdomain = -1;
    int boundaryEntity = -1;
    double temperature = 293.15;
};

struct HeatSource {
    int subdomain = -1;
    int domainEntity = -1;
    double heatRateW = 0.0;
};

struct HeatSourceAssemblyDiagnostic {
    int index = -1;
    int subdomain = -1;
    int domainEntity = -1;
    double configuredPowerW = 0.0;
    int tetCount = 0;
    double sourceVolumeUsedForDensity = 0.0;
    double sourceVolumeAssembled = 0.0;
    double volumetricQ = 0.0;
    double expectedPowerW = 0.0;
    double quadraturePowerW = 0.0;
};

struct ConvectionCondition {
    int subdomain = -1;
    int boundaryEntity = -1;
    double coefficient = 0.0;
    double ambientTemperature = 293.15;
};

struct InterfaceConfig {
    int leftSubdomain = -1;
    int rightSubdomain = -1;
    std::vector<int> leftBoundaryEntities;
    std::vector<int> rightBoundaryEntities;
};

struct SchwarzOptions {
    bool enabled = false;
    std::string type = "multiplicative";
    std::string standaloneMode = "algebraic";
    std::string transmission = "none";
    std::string transmissionOrientation = "forward";
    std::string fluxEval = "sipg_numeric";
    std::string overlapMode = "ras";
    std::string partitionMode = "current";
    int overlapLayers = 1;
    int maxIters = 1000;
    double tolRelUpdate = 1.0e-10;
    double tolRelResidual = 1.0e-10;
    double relaxation = 1.0;
    double robinAlphaFactor = 10.0;
    std::string initialGuess = "previous_time_step";
    bool checkInterfaceJump = true;
    bool outputIterationLog = true;
    bool writeInterfaceFlux = false;
    bool validateAgainstMonolithic = false;
};

struct CaseConfig {
    std::string name = "default";
    std::string solverMethod = "schwarz_precond_fgmres";
    std::string timeIntegrator = "backward_euler";
    SchwarzOptions schwarz;
    double coordinateScale = 1.0e-3;
    double initialTemperature = 293.15;
    double startTime = 0.0;
    double timeStep = 0.1*1e-3;
    int timeSteps = 100;
    double penaltyFactor = 15.0;
    std::string dirichletMethod = "strong";
    double nitschePenaltyFactor = 15.0;
    std::string penaltyMode = "harmonic";
    std::string interfaceScheme = "sipg";
    std::string penaltyScaling = "p_p1";
    bool autoInterfaces = true;
    bool disableWarmStart = false;
    bool forceNontrivialRhs = false;
    double thermalSourceScale = 1.0;
    std::string initialGuessType = "current";
    std::vector<std::pair<int, int>> explicitInterfaceFacePairs;
    std::filesystem::path outputDir;
    std::vector<DomainConfig> domains;
    std::vector<BoundaryCondition> dirichletConditions;
    std::vector<ConvectionCondition> convectionConditions;
    std::vector<HeatSource> heatSources;
    std::vector<InterfaceConfig> interfaces;
    std::filesystem::path comsolComparisonPath;
};

static void setIsotropicConductivity(Material& material, double conductivity)
{
    material.conductivity = conductivity;
    material.conductivityX = conductivity;
    material.conductivityY = conductivity;
    material.conductivityZ = conductivity;
}

static const Material& materialForTet(const CaseConfig& config, const Tet& tet)
{
    const DomainConfig& domain = config.domains[static_cast<size_t>(tet.subdomain)];
    const auto found = domain.materialsByDomainEntity.find(tet.domainEntity);
    if (found != domain.materialsByDomainEntity.end()) {
        return found->second;
    }
    return domain.material;
}

static double normalConductivity(const Material& material, const Vec3& normal)
{
    return material.conductivityX * normal.x * normal.x
         + material.conductivityY * normal.y * normal.y
         + material.conductivityZ * normal.z * normal.z;
}

static double conductivityGradientDot(const Material& material, const Vec3& a, const Vec3& b)
{
    return material.conductivityX * a.x * b.x
         + material.conductivityY * a.y * b.y
         + material.conductivityZ * a.z * b.z;
}

enum class AnalysisMode {
    Transient,
    Steady
};

static std::string analysisModeName(AnalysisMode mode)
{
    return mode == AnalysisMode::Steady ? "steady" : "transient";
}

static std::string normalizeSolverMethodName(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    std::replace(value.begin(), value.end(), '-', '_');
    return value;
}

static std::string normalizeSchwarzPartitionModeName(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    std::replace(value.begin(), value.end(), '-', '_');
    return value;
}

static std::string normalizeSchwarzTransmissionModeName(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    std::replace(value.begin(), value.end(), '-', '_');
    if (value == "none") {
        return "algebraic";
    }
    if (value == "dn") {
        return "dirichlet_neumann";
    }
    if (value == "dd") {
        return "dirichlet_dirichlet";
    }
    if (value == "dr") {
        return "dirichlet_robin";
    }
    if (value == "rr" || value == "robin_robin") {
        return "robin";
    }
    return value;
}

static bool isKnownSchwarzStandaloneMode(const std::string& value)
{
    return value == "algebraic"
        || value == "dirichlet_neumann"
        || value == "dirichlet_dirichlet"
        || value == "dirichlet_robin"
        || value == "robin";
}

static bool isKnownSchwarzTransmissionMode(const std::string& value)
{
    return value == "none" || isKnownSchwarzStandaloneMode(value);
}

static std::string schwarzTransmissionForStandaloneMode(const std::string& value)
{
    return value == "algebraic" ? "none" : value;
}

static bool isPhysicalSchwarzStandaloneMode(const std::string& value)
{
    return value != "algebraic";
}

static std::string normalizeSchwarzTransmissionOrientationName(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    std::replace(value.begin(), value.end(), '-', '_');
    if (value == "left_dirichlet" || value == "left_to_right" || value == "lr") {
        return "forward";
    }
    if (value == "right_dirichlet" || value == "right_to_left" || value == "rl") {
        return "reverse";
    }
    return value;
}

static bool isKnownSchwarzTransmissionOrientation(const std::string& value)
{
    return value == "forward" || value == "reverse";
}

static bool isKnownSchwarzPartitionMode(const std::string& value)
{
    return value == "current"
        || value == "material_aligned"
        || value == "hotspot_contained"
        || value == "vertical_heat_flow_aligned";
}

static bool isKnownSolverMethod(const std::string& value)
{
    return value == "monolithic"
        || value == "schwarz"
        || value == "schwarz_precond_fgmres"
        || value == "schwarz_precond_fgmres_two_level"
        || value == "schur-direct-exact"
        || value == "schur_direct_exact";
}

struct ElementGeometry {
    double detJ = 0.0;
    std::array<Vec3, 4> gradLambda{};
};

struct ImportedMesh {
    std::vector<Vec3> vertices;
    std::vector<std::array<int, 4>> tets;
    std::vector<int> tetEntities;
    std::vector<std::array<int, 3>> triangles;
    std::vector<int> triangleEntities;
};

struct ComparisonPoint {
    Vec3 p;
    double comsolTemperature = 0.0;
};

struct CachedTet {
    int tetId = -1;
    ElementGeometry geo;
    Vec3 lo;
    Vec3 hi;
};

struct SolverStatistics {
    std::string name;
    std::string status = "not_run";
    std::string failureReason;
    int totalIterations = 0;
    int maxIterations = 0;
    int breakdownIteration = -1;
    int parallelWorkers = 1;
    std::string solverMethod;
    std::string preconditioner;
    std::string rasType;
    std::string localSolver;
    std::string initialGuessType = "current";
    std::string localIcShiftMode = "auto";
    std::string coarseSpace;
    std::string coarseCorrection;
    double finalRelativeResidual = std::numeric_limits<double>::quiet_NaN();
    double trueRelativeResidual = std::numeric_limits<double>::quiet_NaN();
    double rhsNorm = std::numeric_limits<double>::quiet_NaN();
    double rhsMin = std::numeric_limits<double>::quiet_NaN();
    double rhsMax = std::numeric_limits<double>::quiet_NaN();
    int rhsNonzeroCount = 0;
    int rhsNearZeroCount = 0;
    double heatSourceTotal = std::numeric_limits<double>::quiet_NaN();
    double dirichletRhsContributionNorm = 0.0;
    double nitscheRhsContributionNorm = 0.0;
    double robinRhsContributionNorm = 0.0;
    double neumannRhsContributionNorm = 0.0;
    double initialXNorm = std::numeric_limits<double>::quiet_NaN();
    double initialAxNorm = std::numeric_limits<double>::quiet_NaN();
    double initialResidualNorm = std::numeric_limits<double>::quiet_NaN();
    double initialRelativeResidual = std::numeric_limits<double>::quiet_NaN();
    double solverTolerance = std::numeric_limits<double>::quiet_NaN();
    bool zeroIterationDueToInitialResidual = false;
    double pcgPAp = std::numeric_limits<double>::quiet_NaN();
    double pcgRTr = std::numeric_limits<double>::quiet_NaN();
    double pcgRMz = std::numeric_limits<double>::quiet_NaN();
    double pcgAlpha = std::numeric_limits<double>::quiet_NaN();
    double pcgBeta = std::numeric_limits<double>::quiet_NaN();
    double pcgResidualNorm = std::numeric_limits<double>::quiet_NaN();
    bool pcgPHasNonFinite = false;
    bool pcgApHasNonFinite = false;
    double preconditionerApplySeconds = 0.0;
    int preconditionerApplyCalls = 0;
    double rasHaloBuildSeconds = 0.0;
    double rasLocalMatrixAssemblySeconds = 0.0;
    double rasLocalFactorizationSeconds = 0.0;
    double rasLocalSolveApplySeconds = 0.0;
    double rasRestrictionSeconds = 0.0;
    double rasCommunicationOrHaloUpdateSeconds = 0.0;
    int rasSetupCount = 0;
    bool rasFactorizationReuse = false;
    bool localDiagScaling = false;
    double diagScalingEps = 1.0e-30;
    double localIcShiftUsedMax = std::numeric_limits<double>::quiet_NaN();
    bool coarseEnabled = false;
    int coarseDim = 0;
    double coarseSetupSeconds = 0.0;
    double coarseSolveSeconds = 0.0;
    size_t coarseMatrixNnz = 0;
    double coarseResidualNorm = std::numeric_limits<double>::quiet_NaN();
    double coarseRhsNorm = std::numeric_limits<double>::quiet_NaN();
    double coarseSolutionNorm = std::numeric_limits<double>::quiet_NaN();
    double coarseCorrectionNorm = std::numeric_limits<double>::quiet_NaN();
    double localCorrectionNorm = std::numeric_limits<double>::quiet_NaN();
    double coarseToLocalNormRatio = std::numeric_limits<double>::quiet_NaN();
    bool deflationEnabled = false;
    int deflationDim = 0;
    double deflationSetupSeconds = 0.0;
    double deflationApplySeconds = 0.0;
    double deflationCorrectionNorm = std::numeric_limits<double>::quiet_NaN();
    double matvecSeconds = 0.0;
    std::vector<double> acceptedIcShiftBySubdomain;
    std::vector<int> icPivotNonpositiveBySubdomain;
    std::vector<int> icPivotTinyBySubdomain;
    std::vector<int> icNonfiniteLBySubdomain;
    double setupSeconds = 0.0;
    double solveSeconds = 0.0;
    double temperatureMin = std::numeric_limits<double>::quiet_NaN();
    double temperatureMax = std::numeric_limits<double>::quiet_NaN();
    double temperatureAverage = std::numeric_limits<double>::quiet_NaN();
    size_t preconditionerBytes = 0;
    size_t workingSetBeforeBytes = 0;
    size_t workingSetAfterBytes = 0;
    size_t peakWorkingSetBytes = 0;
};

struct IcFactorDiagnostics {
    int subdomain = -1;
    int dofs = 0;
    std::string ordering = "natural";
    double appliedShift = 0.0;
    bool diagonalScaling = true;
    bool exactSpd = false;
    bool accepted = false;
    double pivotMin = std::numeric_limits<double>::max();
    double pivotMax = 0.0;
    int pivotNonpositiveCount = 0;
    int pivotTinyCount = 0;
    int firstBadPivotRow = -1;
    double firstBadPivotValue = 0.0;
    int nonFiniteLCount = 0;
    bool breakdown = false;
};

struct ProgramTiming {
    double preprocessingSeconds = 0.0;
    double volumeAssemblySeconds = 0.0;
    double interfaceAssemblySeconds = 0.0;
    double systemBuildSeconds = 0.0;
    double dirichletAssemblySeconds = 0.0;
    double csrFinalizeSeconds = 0.0;
    double assemblySeconds = 0.0;
    double postprocessingSeconds = 0.0;
    double totalSeconds = 0.0;
};

struct DiagonalStats {
    double minDiagonal = 0.0;
    double maxDiagonal = 0.0;
    double sumDiagonal = 0.0;
    int negativeEntries = 0;
    int zeroEntries = 0;
};

struct InterfacePenaltyStats {
    int facePairCount = 0;
    double etaMin = std::numeric_limits<double>::max();
    double etaMax = 0.0;
    double etaSum = 0.0;
    double hMin = std::numeric_limits<double>::max();
    double hMax = 0.0;
    double hSum = 0.0;
    double kRatioMin = std::numeric_limits<double>::max();
    double kRatioMax = 0.0;
};

struct AssemblyDiagnostics {
    std::vector<double> volumeDiag;
    std::vector<double> robinDiag;
    std::vector<double> interfaceConsistencyDiag;
    std::vector<double> interfacePenaltyDiag;
    // Positive trace norms assembled with the same nonmatching face
    // quadrature as SIPG.  The unweighted form is the interface L2 mass
    // diagonal; the penalty form is interfacePenaltyDiag.
    std::vector<double> interfaceTraceMassDiag;
    std::vector<double> preDirichletDiag;
    std::vector<double> dirichletDiag;
    std::vector<double> finalDiag;
    std::vector<char> interfaceDof;
    std::vector<std::pair<int, int>> interfaceBoundaryByDof;
    InterfacePenaltyStats interfacePenaltyStats;
};

struct MatrixEntry {
    int row = 0;
    int col = 0;
    double value = 0.0;
};

struct VectorEntry {
    int index = 0;
    double value = 0.0;
};

static unsigned int solverParallelWorkers()
{
    const unsigned int hw = std::thread::hardware_concurrency();
    const unsigned int hardware = hw == 0 ? 1u : hw;
    std::string workerValue;
#ifdef _WIN32
    char* envBuffer = nullptr;
    size_t envSize = 0;
    if (_dupenv_s(&envBuffer, &envSize, "SIPG_SOLVER_WORKERS") == 0 && envBuffer != nullptr) {
        workerValue = envBuffer;
        std::free(envBuffer);
    }
#else
    const char* env = std::getenv("SIPG_SOLVER_WORKERS");
    if (env != nullptr) {
        workerValue = env;
    }
#endif
    if (!workerValue.empty()) {
        try {
            const int requested = std::stoi(workerValue);
            if (requested > 0) {
                return std::max(1u, std::min(hardware, static_cast<unsigned int>(requested)));
            }
        } catch (...) {
            // Ignore malformed environment values and keep the conservative default.
        }
    }
    return std::max(1u, std::min(hardware, 4u));
}

// CSV is the production diagnostic interchange format. Keeping this tiny
// utility beside the shared core types avoids depending on the research-only
// diagnostics header merely to quote a field safely.
static std::string csvEscape(const std::string& value)
{
    if (value.find_first_of(",\"\n\r") == std::string::npos) {
        return value;
    }
    std::string escaped = "\"";
    for (const char character : value) {
        if (character == '\"') {
            escaped += "\"\"";
        } else {
            escaped += character;
        }
    }
    escaped += '\"';
    return escaped;
}

template <typename Work>
static void parallelFor(size_t count, Work&& work)
{
    const unsigned int workers = solverParallelWorkers();
    constexpr size_t minChunk = 1024;
    if (workers <= 1 || count < minChunk * 2) {
        for (size_t i = 0; i < count; ++i) {
            work(i);
        }
        return;
    }

#ifdef _OPENMP
    const int threadCount = static_cast<int>(std::min<size_t>(workers, (count + minChunk - 1) / minChunk));
#pragma omp parallel for schedule(static) num_threads(threadCount)
    for (long long i = 0; i < static_cast<long long>(count); ++i) {
        work(static_cast<size_t>(i));
    }
#else
    for (size_t i = 0; i < count; ++i) {
        work(i);
    }
#endif
}

template <typename Work>
static void parallelForCoarse(size_t count, Work&& work)
{
    if (count == 0) {
        return;
    }
    const unsigned int workers = std::min<unsigned int>(solverParallelWorkers(), static_cast<unsigned int>(count));
    if (workers <= 1) {
        for (size_t i = 0; i < count; ++i) {
            work(i);
        }
        return;
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(static_cast<int>(workers))
    for (long long i = 0; i < static_cast<long long>(count); ++i) {
        work(static_cast<size_t>(i));
    }
#else
    for (size_t i = 0; i < count; ++i) {
        work(i);
    }
#endif
}

struct SparseMatrix {
    int n = 0;
    std::vector<MatrixEntry> triplets;
    std::vector<int> rowPtr;
    std::vector<int> colInd;
    std::vector<double> values;
    bool csrReady = false;

    explicit SparseMatrix(int nIn = 0) : n(nIn) {}
    int size() const { return n; }

    void add(int i, int j, double value)
    {
        if (std::abs(value) > 0.0) {
            triplets.push_back({i, j, value});
            csrReady = false;
        }
    }

    void appendEntries(std::vector<MatrixEntry>& entries)
    {
        if (entries.empty()) {
            return;
        }
        triplets.insert(triplets.end(),
                        std::make_move_iterator(entries.begin()),
                        std::make_move_iterator(entries.end()));
        std::vector<MatrixEntry>().swap(entries);
        csrReady = false;
    }

    void appendScaledEntries(const SparseMatrix& source, double scale)
    {
        if (std::abs(scale) == 0.0) {
            return;
        }
        source.forEachEntry([&](int row, int col, double value) {
            add(row, col, scale * value);
        });
    }

    template <typename Fn>
    void forEachEntry(Fn&& fn) const
    {
        if (csrReady) {
            for (int row = 0; row < n; ++row) {
                for (int k = rowPtr[static_cast<size_t>(row)]; k < rowPtr[static_cast<size_t>(row + 1)]; ++k) {
                    fn(row, colInd[static_cast<size_t>(k)], values[static_cast<size_t>(k)]);
                }
            }
            return;
        }
        for (const MatrixEntry& entry : triplets) {
            fn(entry.row, entry.col, entry.value);
        }
    }

    void finalizeCsr()
    {
        std::sort(triplets.begin(), triplets.end(), [](const MatrixEntry& a, const MatrixEntry& b) {
            return std::tie(a.row, a.col) < std::tie(b.row, b.col);
        });

        rowPtr.assign(static_cast<size_t>(n) + 1, 0);
        colInd.clear();
        values.clear();
        colInd.reserve(triplets.size());
        values.reserve(triplets.size());

        size_t i = 0;
        while (i < triplets.size()) {
            const int row = triplets[i].row;
            const int col = triplets[i].col;
            double sum = 0.0;
            while (i < triplets.size() && triplets[i].row == row && triplets[i].col == col) {
                sum += triplets[i].value;
                ++i;
            }
            if (row < 0 || row >= n || col < 0 || col >= n || std::abs(sum) == 0.0) {
                continue;
            }
            if (values.size() >= static_cast<size_t>(std::numeric_limits<int>::max())) {
                throw std::runtime_error("Sparse matrix has too many nonzeros for 32-bit CSR indices.");
            }
            colInd.push_back(col);
            values.push_back(sum);
            ++rowPtr[static_cast<size_t>(row + 1)];
        }

        for (int row = 0; row < n; ++row) {
            rowPtr[static_cast<size_t>(row + 1)] += rowPtr[static_cast<size_t>(row)];
        }
        std::vector<MatrixEntry>().swap(triplets);
        csrReady = true;
    }

    double diagonal(int i) const
    {
        if (csrReady) {
            for (int k = rowPtr[static_cast<size_t>(i)]; k < rowPtr[static_cast<size_t>(i + 1)]; ++k) {
                if (colInd[static_cast<size_t>(k)] == i) {
                    return values[static_cast<size_t>(k)];
                }
            }
            return 0.0;
        }
        double diag = 0.0;
        for (const MatrixEntry& entry : triplets) {
            if (entry.row == i && entry.col == i) {
                diag += entry.value;
            }
        }
        return diag;
    }

    std::vector<double> multiply(const std::vector<double>& x) const
    {
        std::vector<double> y(static_cast<size_t>(size()), 0.0);
        if (csrReady) {
            parallelFor(y.size(), [&](size_t i) {
                double sum = 0.0;
                for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
                    sum += values[static_cast<size_t>(k)] * x[static_cast<size_t>(colInd[static_cast<size_t>(k)])];
                }
                y[i] = sum;
            });
            return y;
        }
        for (const MatrixEntry& entry : triplets) {
            if (entry.row >= 0 && entry.row < n && entry.col >= 0 && entry.col < n) {
                y[static_cast<size_t>(entry.row)] += entry.value * x[static_cast<size_t>(entry.col)];
            }
        }
        return y;
    }
};

static std::vector<double> matrixDiagonalVector(const SparseMatrix& matrix)
{
    std::vector<double> diag(static_cast<size_t>(matrix.size()), 0.0);
    matrix.forEachEntry([&](int row, int col, double value) {
        if (row == col && row >= 0 && row < matrix.size()) {
            diag[static_cast<size_t>(row)] += value;
        }
    });
    return diag;
}

static void addDiagonalEntry(std::vector<double>& diag, int row, int col, double value)
{
    if (row == col && row >= 0 && row < static_cast<int>(diag.size())) {
        diag[static_cast<size_t>(row)] += value;
    }
}

static DiagonalStats diagonalStats(const std::vector<double>& diagonal)
{
    DiagonalStats stats;
    if (diagonal.empty()) {
        return stats;
    }
    stats.minDiagonal = std::numeric_limits<double>::max();
    stats.maxDiagonal = -std::numeric_limits<double>::max();
    constexpr double zeroTolerance = 1.0e-30;
    for (double value : diagonal) {
        stats.minDiagonal = std::min(stats.minDiagonal, value);
        stats.maxDiagonal = std::max(stats.maxDiagonal, value);
        stats.sumDiagonal += value;
        if (value < 0.0) {
            ++stats.negativeEntries;
        }
        if (std::abs(value) <= zeroTolerance) {
            ++stats.zeroEntries;
        }
    }
    return stats;
}

static size_t currentWorkingSetBytes()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX counters{};
    counters.cb = sizeof(counters);
    if (GetProcessMemoryInfo(GetCurrentProcess(),
                             reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&counters),
                             sizeof(counters)) != 0) {
        return static_cast<size_t>(counters.WorkingSetSize);
    }
#endif
    return 0;
}

static size_t peakWorkingSetBytes()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX counters{};
    counters.cb = sizeof(counters);
    if (GetProcessMemoryInfo(GetCurrentProcess(),
                             reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&counters),
                             sizeof(counters)) != 0) {
        return static_cast<size_t>(counters.PeakWorkingSetSize);
    }
#endif
    return 0;
}

static double megabytes(size_t bytes)
{
    return static_cast<double>(bytes) / (1024.0 * 1024.0);
}

class ScopedMklSingleThread {
public:
    ScopedMklSingleThread()
    {
#ifdef USE_MKL_PARDISO
        previousThreads_ = mkl_get_max_threads();
        mkl_set_num_threads_local(1);
#endif
    }

    ~ScopedMklSingleThread()
    {
#ifdef USE_MKL_PARDISO
        mkl_set_num_threads_local(previousThreads_);
#endif
    }

private:
#ifdef USE_MKL_PARDISO
    int previousThreads_ = 1;
#endif
};

