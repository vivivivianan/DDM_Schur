#include "optimal_port_space.hpp"

#include "ddm_schur/interface_operator.hpp"
#include "mor/local/local_reduced_schur.hpp"
#include "mor/local/local_subdomain_model.hpp"
#include "port_eigensolver.hpp"
#include "sipg_core.hpp"
#include "transfer_operator.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

namespace mor::transient {
namespace {

using Clock = std::chrono::steady_clock;

double dot(const double* x, const double* y, int n)
{
    long double value = 0.0L;
    for (int i = 0; i < n; ++i) value += static_cast<long double>(x[i]) * y[i];
    return static_cast<double>(value);
}

double massDot(const double* x, const double* y,
               const std::vector<double>& mass)
{
    long double value = 0.0L;
    for (std::size_t i = 0; i < mass.size(); ++i) {
        value += static_cast<long double>(mass[i]) * x[i] * y[i];
    }
    return static_cast<double>(value);
}

template <typename T>
void hashValue(std::uint64_t& hash, const T& value)
{
    const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
    for (std::size_t i = 0; i < sizeof(T); ++i) {
        hash ^= bytes[i];
        hash *= UINT64_C(1099511628211);
    }
}

std::vector<double> traceWeights(
    const local::Model& model, const std::vector<int>& target,
    const std::vector<double>& traceMassDiagonal,
    const std::vector<double>& penaltyMassDiagonal,
    const std::string& innerProduct,
    double* rawMinimum = nullptr,
    double* rawMaximum = nullptr,
    double* regularization = nullptr)
{
    if (traceMassDiagonal.size() != penaltyMassDiagonal.size()
        || traceMassDiagonal.empty()) {
        throw std::runtime_error(
            "[Optimal port] SIPG trace Gram diagonals are unavailable.");
    }
    const std::vector<double>* diagonal = nullptr;
    if (innerProduct == "trace-mass") {
        diagonal = &traceMassDiagonal;
    } else if (innerProduct == "penalty-weighted-mass") {
        diagonal = &penaltyMassDiagonal;
    } else {
        throw std::runtime_error(
            "[Optimal port] Unsupported inner product: " + innerProduct);
    }
    std::vector<double> result(target.size(), 0.0);
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = 0.0;
    for (std::size_t row = 0; row < target.size(); ++row) {
        const int gamma = target[row];
        if (gamma < 0
            || gamma >= static_cast<int>(model.interfaceGlobalDofs.size())) {
            throw std::runtime_error(
                "[Optimal port] Target trace index is out of range.");
        }
        const int global =
            model.interfaceGlobalDofs[static_cast<std::size_t>(gamma)];
        if (global < 0
            || global >= static_cast<int>(diagonal->size())) {
            throw std::runtime_error(
                "[Optimal port] Target global DOF is out of Gram range.");
        }
        const double value = (*diagonal)[static_cast<std::size_t>(global)];
        if (!std::isfinite(value) || value < 0.0) {
            throw std::runtime_error(
                "[Optimal port] Trace Gram diagonal is not nonnegative finite.");
        }
        result[row] = value;
        minimum = std::min(minimum, value);
        maximum = std::max(maximum, value);
    }
    if (!(maximum > 0.0)) {
        throw std::runtime_error(
            "[Optimal port] Trace Gram is identically zero on a target port.");
    }
    // SIPG consistency couples a small set of normal-derivative DOFs whose
    // value trace is exactly zero.  The physical Gram is therefore positive
    // semidefinite.  A machine-scale diagonal floor makes the generalized
    // eigenproblem definite without replacing the assembled trace metric.
    const double floor = 128.0 * std::numeric_limits<double>::epsilon()
        * maximum;
    for (double& value : result) {
        value = std::max(value, floor);
    }
    if (rawMinimum != nullptr) *rawMinimum = minimum;
    if (rawMaximum != nullptr) *rawMaximum = maximum;
    if (regularization != nullptr) {
        *regularization = minimum < floor ? floor : 0.0;
    }
    return result;
}

std::map<int, int> readRankMap(const std::filesystem::path& path)
{
    std::map<int, int> result;
    if (path.empty()) return result;
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error(
            "[Optimal port] Cannot open rank file: " + path.string());
    }
    std::string header;
    std::getline(input, header);
    std::vector<std::string> names;
    for (std::size_t begin = 0; begin <= header.size();) {
        const std::size_t end = header.find(',', begin);
        names.push_back(header.substr(begin,
            end == std::string::npos ? std::string::npos : end - begin));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    int interfaceColumn = -1;
    int rankColumn = -1;
    for (int column = 0; column < static_cast<int>(names.size()); ++column) {
        if (names[static_cast<std::size_t>(column)] == "interface_id") {
            interfaceColumn = column;
        }
        if (names[static_cast<std::size_t>(column)] == "selected_rank"
            || names[static_cast<std::size_t>(column)] == "rank") {
            rankColumn = column;
        }
    }
    if (interfaceColumn < 0 || rankColumn < 0) {
        throw std::runtime_error(
            "[Optimal port] Rank file needs interface_id and selected_rank columns.");
    }
    std::string line;
    while (std::getline(input, line)) {
        std::vector<std::string> values;
        for (std::size_t begin = 0; begin <= line.size();) {
            const std::size_t end = line.find(',', begin);
            values.push_back(line.substr(begin,
                end == std::string::npos
                    ? std::string::npos : end - begin));
            if (end == std::string::npos) break;
            begin = end + 1;
        }
        if (interfaceColumn >= static_cast<int>(values.size())
            || rankColumn >= static_cast<int>(values.size())) {
            continue;
        }
        const int interfaceId =
            std::stoi(values[static_cast<std::size_t>(interfaceColumn)]);
        const int rank =
            std::stoi(values[static_cast<std::size_t>(rankColumn)]);
        if (rank <= 0) {
            throw std::runtime_error(
                "[Optimal port] Rank-file values must be positive.");
        }
        result[interfaceId] = rank;
    }
    return result;
}

std::uint64_t fingerprintVector(const std::vector<double>& values)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    hashValue(hash, values.size());
    for (double value : values) hashValue(hash, value);
    return hash;
}

double portOrthogonalityError(const LocalPortBasis& port,
                              const std::vector<double>& mass)
{
    double maximum = 0.0;
    for (int left = 0; left < port.rank; ++left) {
        for (int right = 0; right <= left; ++right) {
            const double product = massDot(
                port.basis.data()
                    + static_cast<std::size_t>(left * port.rows),
                port.basis.data()
                    + static_cast<std::size_t>(right * port.rows),
                mass);
            maximum = std::max(maximum,
                std::abs(product - (left == right ? 1.0 : 0.0)));
        }
    }
    return maximum;
}

void symmetricJacobi(std::vector<double>& matrix, int n,
                     std::vector<double>& values, std::vector<double>& vectors)
{
    vectors.assign(static_cast<std::size_t>(n * n), 0.0);
    for (int i = 0; i < n; ++i) vectors[static_cast<std::size_t>(i * n + i)] = 1.0;
    const int sweeps = std::max(12, 8 * n * n);
    for (int sweep = 0; sweep < sweeps; ++sweep) {
        int p = 0, q = 0;
        double maximum = 0.0;
        for (int row = 0; row < n; ++row) for (int column = row + 1; column < n; ++column) {
            const double value = std::abs(matrix[static_cast<std::size_t>(row * n + column)]);
            if (value > maximum) { maximum = value; p = row; q = column; }
        }
        if (maximum <= 64.0 * std::numeric_limits<double>::epsilon()) break;
        const double app = matrix[static_cast<std::size_t>(p * n + p)];
        const double aqq = matrix[static_cast<std::size_t>(q * n + q)];
        const double apq = matrix[static_cast<std::size_t>(p * n + q)];
        const double angle = 0.5 * std::atan2(2.0 * apq, aqq - app);
        const double c = std::cos(angle), s = std::sin(angle);
        for (int k = 0; k < n; ++k) {
            const double mkp = matrix[static_cast<std::size_t>(k * n + p)];
            const double mkq = matrix[static_cast<std::size_t>(k * n + q)];
            matrix[static_cast<std::size_t>(k * n + p)] = c * mkp - s * mkq;
            matrix[static_cast<std::size_t>(k * n + q)] = s * mkp + c * mkq;
        }
        for (int k = 0; k < n; ++k) {
            const double mpk = matrix[static_cast<std::size_t>(p * n + k)];
            const double mqk = matrix[static_cast<std::size_t>(q * n + k)];
            matrix[static_cast<std::size_t>(p * n + k)] = c * mpk - s * mqk;
            matrix[static_cast<std::size_t>(q * n + k)] = s * mpk + c * mqk;
        }
        for (int k = 0; k < n; ++k) {
            const double vkp = vectors[static_cast<std::size_t>(k * n + p)];
            const double vkq = vectors[static_cast<std::size_t>(k * n + q)];
            vectors[static_cast<std::size_t>(k * n + p)] = c * vkp - s * vkq;
            vectors[static_cast<std::size_t>(k * n + q)] = s * vkp + c * vkq;
        }
    }
    values.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) values[static_cast<std::size_t>(i)] = matrix[static_cast<std::size_t>(i * n + i)];
}

bool appendMassMode(LocalPortBasis& port, std::vector<double> candidate,
                    const std::vector<double>& mass, int maximumRank,
                    double tolerance)
{
    ++port.candidateColumns;
    if (port.rank >= maximumRank) return false;
    const double original = std::sqrt(std::max(0.0, massDot(candidate.data(), candidate.data(), mass)));
    if (!(original > 0.0)) return false;
    for (int pass = 0; pass < 2; ++pass) for (int mode = 0; mode < port.rank; ++mode) {
        const double* basis = port.basis.data() + static_cast<std::size_t>(mode * port.rows);
        const double coefficient = massDot(basis, candidate.data(), mass);
        for (int row = 0; row < port.rows; ++row) candidate[static_cast<std::size_t>(row)] -= coefficient * basis[row];
    }
    const double norm = std::sqrt(std::max(0.0, massDot(candidate.data(), candidate.data(), mass)));
    if (!(norm > tolerance * original)) return false;
    for (double& value : candidate) value /= norm;
    port.basis.insert(port.basis.end(), candidate.begin(), candidate.end());
    ++port.rank;
    ++port.acceptedColumns;
    return true;
}

} // namespace

LocalPortModel buildSteklovSchurPortModel(
    const Mesh& mesh, const ddm_schur::InterfacePartition& partition,
    const local::Model& dynamicModel,
    const std::vector<double>& traceMassDiagonal,
    const std::vector<double>& penaltyMassDiagonal,
    const SteklovPortOptions& options)
{
    if (options.requestedRank <= 0 || options.directRowLimit <= 0
        || options.inverseIterations <= 0 || !(options.relativeTolerance > 0.0)) {
        throw std::runtime_error("[Optimal port] Invalid Steklov options.");
    }
    const auto start = Clock::now();
    LocalPortModel model;
    model.formatVersion = 5;
    model.basisMethod = "steklov-schur";
    model.ablationMode = "steklov-schur";
    model.innerProduct = options.innerProduct;
    model.rankMode = "fixed";
    model.requestedRank = options.requestedRank;
    model.minimumRank = options.requestedRank;
    model.maximumRank = options.requestedRank;
    model.eigensolverTolerance = options.relativeTolerance;
    model.eigensolverMaximumIterations = options.inverseIterations;
    model.relativeDeflationTolerance = options.relativeTolerance;
    model.innerSolver = "direct";
    model.fullInterfaceDofs = dynamicModel.interfaceDofs;
    model.interfaceGlobalDofs = dynamicModel.interfaceGlobalDofs;
    ReducedDynamicSchurOperator action(dynamicModel);
    const auto patches = buildOptimalPortPatches(mesh, partition, 0);
    for (const PortPatch& patch : patches) {
        const int rows = static_cast<int>(patch.target.size());
        if (rows > options.directRowLimit) {
            throw std::runtime_error("[Optimal port] Steklov direct patch limit exceeded on interface "
                + std::to_string(patch.interfaceId)
                + "; use optimal-transfer matrix-free mode for this port.");
        }
        LocalPortBasis port;
        port.interfaceId = patch.interfaceId;
        port.leftSubdomain = patch.leftSubdomain;
        port.rightSubdomain = patch.rightSubdomain;
        port.rows = rows;
        port.interfaceIndices = patch.target;
        port.sourceIndices = patch.source;
        port.patchSubdomains = patch.patchSubdomains;
        port.targetFingerprint = patch.targetFingerprint;
        port.sourceFingerprint = patch.sourceFingerprint;
        port.templateId = port.interfaceId;
        const std::vector<double> mass = traceWeights(
            dynamicModel, port.interfaceIndices, traceMassDiagonal,
            penaltyMassDiagonal, options.innerProduct);
        std::vector<double> matrix(static_cast<std::size_t>(rows * rows), 0.0);
        std::vector<double> full(static_cast<std::size_t>(dynamicModel.interfaceDofs), 0.0), image;
        for (int column = 0; column < rows; ++column) {
            full[static_cast<std::size_t>(port.interfaceIndices[static_cast<std::size_t>(column)])] = 1.0;
            action.apply(full, image);
            for (int row = 0; row < rows; ++row) matrix[static_cast<std::size_t>(row * rows + column)] =
                image[static_cast<std::size_t>(port.interfaceIndices[static_cast<std::size_t>(row)])];
            full[static_cast<std::size_t>(port.interfaceIndices[static_cast<std::size_t>(column)])] = 0.0;
        }
        for (int row = 0; row < rows; ++row) for (int column = 0; column < row; ++column) {
            const double average = 0.5 * (matrix[static_cast<std::size_t>(row * rows + column)]
                + matrix[static_cast<std::size_t>(column * rows + row)]);
            matrix[static_cast<std::size_t>(row * rows + column)] = average;
            matrix[static_cast<std::size_t>(column * rows + row)] = average;
        }
        const local::DenseSymmetricFactor factor = local::factorDenseSymmetric(matrix, rows);
        if (!factor.cholesky) {
            throw std::runtime_error(
                "[Optimal port] Steklov S_tt is not SPD on interface "
                + std::to_string(patch.interfaceId) + '.');
        }
        const int block = std::min(rows, options.requestedRank);
        std::vector<double> q(static_cast<std::size_t>(rows * block), 0.0);
        for (int column = 0; column < block; ++column) for (int row = 0; row < rows; ++row) {
            q[static_cast<std::size_t>(column * rows + row)] = std::sin(
                0.6180339887498949 * static_cast<double>((row + 1) * (column + 3)));
        }
        for (int iteration = 0; iteration < options.inverseIterations; ++iteration) {
            for (int column = 0; column < block; ++column) {
                std::vector<double> rhs(static_cast<std::size_t>(rows));
                for (int row = 0; row < rows; ++row) rhs[static_cast<std::size_t>(row)] = mass[static_cast<std::size_t>(row)]
                    * q[static_cast<std::size_t>(column * rows + row)];
                local::solveDenseSymmetric(factor, rhs);
                std::copy(rhs.begin(), rhs.end(), q.begin() + static_cast<std::ptrdiff_t>(column * rows));
                for (int pass = 0; pass < 2; ++pass) for (int prior = 0; prior < column; ++prior) {
                    double* current = q.data() + static_cast<std::size_t>(column * rows);
                    const double* previous = q.data() + static_cast<std::size_t>(prior * rows);
                    const double coefficient = massDot(previous, current, mass);
                    for (int row = 0; row < rows; ++row) current[row] -= coefficient * previous[row];
                }
                double* current = q.data() + static_cast<std::size_t>(column * rows);
                const double norm = std::sqrt(std::max(0.0, massDot(current, current, mass)));
                if (!(norm > options.relativeTolerance)) throw std::runtime_error("[Optimal port] Steklov inverse iteration deflated.");
                for (int row = 0; row < rows; ++row) current[row] /= norm;
            }
        }
        std::vector<double> reduced(static_cast<std::size_t>(block * block), 0.0);
        for (int column = 0; column < block; ++column) {
            for (int row = 0; row < rows; ++row) {
                double value = 0.0;
                for (int k = 0; k < rows; ++k) value += matrix[static_cast<std::size_t>(row * rows + k)]
                    * q[static_cast<std::size_t>(column * rows + k)];
                full[static_cast<std::size_t>(row)] = value;
            }
            for (int left = 0; left < block; ++left) reduced[static_cast<std::size_t>(left * block + column)] =
                dot(q.data() + static_cast<std::size_t>(left * rows), full.data(), rows);
        }
        std::vector<double> eigenvalues, eigenvectors;
        symmetricJacobi(reduced, block, eigenvalues, eigenvectors);
        std::vector<int> order(static_cast<std::size_t>(block));
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](int a, int b) { return eigenvalues[static_cast<std::size_t>(a)] < eigenvalues[static_cast<std::size_t>(b)]; });
        for (int coordinate = 0; coordinate < 4; ++coordinate) {
            std::vector<double> candidate(static_cast<std::size_t>(rows), 1.0);
            for (int row = 0; row < rows; ++row) {
                const Vec3& point = mesh.nodes[static_cast<std::size_t>(dynamicModel.interfaceGlobalDofs[static_cast<std::size_t>(port.interfaceIndices[static_cast<std::size_t>(row)])])].p;
                if (coordinate == 1) candidate[static_cast<std::size_t>(row)] = point.x;
                if (coordinate == 2) candidate[static_cast<std::size_t>(row)] = point.y;
                if (coordinate == 3) candidate[static_cast<std::size_t>(row)] = point.z;
            }
            if (appendMassMode(port, std::move(candidate), mass, options.requestedRank, options.relativeTolerance)) ++port.mandatoryModes;
        }
        for (int selected : order) {
            std::vector<double> candidate(static_cast<std::size_t>(rows), 0.0);
            for (int column = 0; column < block; ++column) {
                const double coefficient = eigenvectors[static_cast<std::size_t>(column * block + selected)];
                for (int row = 0; row < rows; ++row) candidate[static_cast<std::size_t>(row)] += coefficient
                    * q[static_cast<std::size_t>(column * rows + row)];
            }
            std::vector<double> residual(static_cast<std::size_t>(rows), 0.0);
            for (int row = 0; row < rows; ++row) for (int column = 0; column < rows; ++column) residual[static_cast<std::size_t>(row)] +=
                matrix[static_cast<std::size_t>(row * rows + column)] * candidate[static_cast<std::size_t>(column)];
            for (int row = 0; row < rows; ++row) residual[static_cast<std::size_t>(row)] -= eigenvalues[static_cast<std::size_t>(selected)]
                * mass[static_cast<std::size_t>(row)] * candidate[static_cast<std::size_t>(row)];
            const double denominator = std::max(1.0e-300, std::sqrt(dot(candidate.data(), candidate.data(), rows)));
            const double residualNorm = std::sqrt(dot(residual.data(), residual.data(), rows)) / denominator;
            if (appendMassMode(port, std::move(candidate), mass, options.requestedRank, options.relativeTolerance)) {
                ++port.spectralModes;
                port.spectralValues.push_back(eigenvalues[static_cast<std::size_t>(selected)]);
                port.spectralResiduals.push_back(residualNorm);
            }
        }
        for (int left = 0; left < port.rank; ++left) for (int right = 0; right <= left; ++right) {
            const double product = massDot(port.basis.data() + static_cast<std::size_t>(left * rows),
                port.basis.data() + static_cast<std::size_t>(right * rows), mass);
            port.orthogonalityError = std::max(port.orthogonalityError,
                std::abs(product - (left == right ? 1.0 : 0.0)));
        }
        std::uint64_t fingerprint = UINT64_C(1469598103934665603);
        hashValue(fingerprint, port.interfaceId); hashValue(fingerprint, port.rank);
        for (int index : port.interfaceIndices) hashValue(fingerprint, index);
        for (double value : port.basis) hashValue(fingerprint, value);
        port.fingerprint = fingerprint;
        model.reducedInterfaceDofs += port.rank;
        model.ports.push_back(std::move(port));
    }
    model.basisSeconds = std::chrono::duration<double>(Clock::now() - start).count();
    model.modelBytes = model.interfaceGlobalDofs.capacity() * sizeof(int);
    for (const auto& port : model.ports) model.modelBytes += port.interfaceIndices.capacity() * sizeof(int)
        + port.basis.capacity() * sizeof(double) + port.spectralValues.capacity() * sizeof(double)
        + port.spectralResiduals.capacity() * sizeof(double);
    return model;
}

OptimalPortBuildResult buildOptimalTransferPortModel(
    const Mesh& mesh,
    const ddm_schur::InterfacePartition& partition,
    const local::Model& dynamicModel,
    const std::vector<double>& traceMassDiagonal,
    const std::vector<double>& penaltyMassDiagonal,
    const std::vector<double>& input,
    int sourceChannels,
    const std::vector<double>& boundaryLoad,
    const std::vector<double>& condensedHistory,
    int historyChannels,
    const OptimalTransferPortOptions& options,
    const ReducedDynamicSchurOperator* sharedSchur)
{
    if (options.requestedRank <= 0
        || options.minimumRank <= 0
        || options.maximumRank < options.minimumRank
        || !(options.eigenvalueTolerance > 0.0)
        || !(options.eigensolverTolerance > 0.0)
        || options.eigensolverMaximumIterations <= 0
        || !(options.relativeDeflationTolerance > 0.0)) {
        throw std::runtime_error(
            "[Optimal port] Invalid optimal-transfer options.");
    }
    if (options.rankMode != "fixed"
        && options.rankMode != "match-m7"
        && options.rankMode != "eigenvalue-tolerance") {
        throw std::runtime_error(
            "[Optimal port] Unsupported rank mode.");
    }
    if (options.sourceMode != "trace-only"
        && options.sourceMode != "trace-plus-input"
        && options.sourceMode != "generalized-dynamic") {
        throw std::runtime_error(
            "[Optimal port] Unsupported generalized source mode.");
    }
    if (options.ablationMode != "mandatory-transfer"
        && options.ablationMode != "mandatory-only"
        && options.ablationMode != "transfer-only"
        && options.ablationMode != "constant-only"
        && options.ablationMode != "geometry-particular"
        && options.ablationMode != "constant-geometry"
        && options.ablationMode != "input-particular-only"
        && options.ablationMode != "trace-transfer-only"
        && options.ablationMode != "generalized-transfer-only"
        && options.ablationMode != "constant-geometry-trace"
        && options.ablationMode != "constant-geometry-generalized"
        && options.ablationMode != "original-mandatory-trace") {
        throw std::runtime_error(
            "[Optimal port] Unsupported ablation mode.");
    }
    if (historyChannels < 0
        || condensedHistory.size() != static_cast<std::size_t>(
            dynamicModel.interfaceDofs * historyChannels)) {
        throw std::runtime_error(
            "[Optimal port] Condensed history source dimensions are invalid.");
    }
    const auto totalStart = Clock::now();
    const bool includeConstant =
        options.ablationMode == "mandatory-transfer"
        || options.ablationMode == "mandatory-only"
        || options.ablationMode == "constant-only"
        || options.ablationMode == "constant-geometry"
        || options.ablationMode == "constant-geometry-trace"
        || options.ablationMode == "constant-geometry-generalized"
        || options.ablationMode == "original-mandatory-trace";
    const bool includeGeometry =
        options.ablationMode == "mandatory-transfer"
        || options.ablationMode == "mandatory-only"
        || options.ablationMode == "geometry-particular"
        || options.ablationMode == "constant-geometry"
        || options.ablationMode == "constant-geometry-trace"
        || options.ablationMode == "constant-geometry-generalized"
        || options.ablationMode == "original-mandatory-trace";
    const bool includeParticular =
        options.ablationMode == "mandatory-transfer"
        || options.ablationMode == "mandatory-only"
        || options.ablationMode == "geometry-particular"
        || options.ablationMode == "input-particular-only"
        || options.ablationMode == "original-mandatory-trace";
    const bool includeTransfer =
        options.ablationMode == "mandatory-transfer"
        || options.ablationMode == "transfer-only"
        || options.ablationMode == "trace-transfer-only"
        || options.ablationMode == "generalized-transfer-only"
        || options.ablationMode == "constant-geometry-trace"
        || options.ablationMode == "constant-geometry-generalized"
        || options.ablationMode == "original-mandatory-trace";
    std::string sourceMode = options.sourceMode;
    if (options.ablationMode == "trace-transfer-only"
        || options.ablationMode == "constant-geometry-trace"
        || options.ablationMode == "original-mandatory-trace") {
        sourceMode = "trace-only";
    } else if (options.ablationMode == "generalized-transfer-only"
               || options.ablationMode
                    == "constant-geometry-generalized") {
        sourceMode = "generalized-dynamic";
    }
    const bool includeTraceSource = includeTransfer;
    const bool includeInputSource = includeTransfer
        && (sourceMode == "trace-plus-input"
            || sourceMode == "generalized-dynamic");
    const bool includeBoundarySource = includeTransfer
        && sourceMode == "generalized-dynamic";
    const bool includeHistorySource = includeTransfer
        && sourceMode == "generalized-dynamic";
    OptimalPortBuildResult result;
    result.model.formatVersion = 7;
    result.model.basisMethod = "optimal-transfer";
    result.model.methodDescription =
        "Operator-Informed Port Space with Transfer Spectral Enrichment";
    result.model.ablationMode = options.ablationMode;
    result.model.sourceMode = sourceMode;
    result.model.innerProduct = options.innerProduct;
    result.model.rankMode = options.rankMode;
    result.model.oversamplingLayers = options.oversamplingLayers;
    result.model.requestedRank = options.requestedRank;
    result.model.minimumRank = options.minimumRank;
    result.model.maximumRank = options.maximumRank;
    result.model.eigenvalueTolerance = options.eigenvalueTolerance;
    result.model.eigensolverTolerance = options.eigensolverTolerance;
    result.model.eigensolverMaximumIterations =
        options.eigensolverMaximumIterations;
    result.model.relativeDeflationTolerance =
        options.relativeDeflationTolerance;
    result.model.innerSolver = options.innerSolver.innerSolver;
    result.model.innerSolverTolerance =
        options.innerSolver.relativeTolerance;
    result.model.innerSolverMaximumIterations =
        options.innerSolver.maximumIterations;
    result.model.fullInterfaceDofs = dynamicModel.interfaceDofs;
    result.model.interfaceGlobalDofs = dynamicModel.interfaceGlobalDofs;

    const auto patchStart = Clock::now();
    std::vector<PortPatch> patches = buildOptimalPortPatches(
        mesh, partition, options.oversamplingLayers);
    if (!options.selectedInterfaceIds.empty()) {
        const std::set<int> selected(
            options.selectedInterfaceIds.begin(),
            options.selectedInterfaceIds.end());
        patches.erase(std::remove_if(
            patches.begin(), patches.end(),
            [&](const PortPatch& patch) {
                return selected.count(patch.interfaceId) == 0;
            }), patches.end());
        if (patches.size() != selected.size()) {
            throw std::runtime_error(
                "[Optimal port] A selected pilot interface is absent.");
        }
    }
    result.patchSetupSeconds =
        std::chrono::duration<double>(Clock::now() - patchStart).count();
    const std::map<int, int> rankMap = readRankMap(options.rankFile);
    if (options.rankMode == "match-m7" && rankMap.empty()) {
        throw std::runtime_error(
            "[Optimal port] match-m7 needs a nonempty rank file.");
    }
    GeneralizedTransferSourceBlocks sourceBlocks =
        buildGeneralizedTransferSourceBlocks(
            dynamicModel, input, sourceChannels,
            boundaryLoad, condensedHistory, historyChannels);
    const std::vector<double>& condensedInputs = sourceBlocks.input;
    const int boundaryChannels = sourceBlocks.boundaryChannels;
    result.model.generalizedInputFingerprint =
        sourceBlocks.inputFingerprint;
    result.model.generalizedBoundaryFingerprint =
        sourceBlocks.boundaryFingerprint;
    result.model.generalizedHistoryFingerprint =
        sourceBlocks.historyFingerprint;
    std::unique_ptr<ReducedDynamicSchurOperator> ownedSchur;
    if (sharedSchur == nullptr) {
        ownedSchur =
            std::make_unique<ReducedDynamicSchurOperator>(dynamicModel);
        sharedSchur = ownedSchur.get();
    }
    const ReducedDynamicSchurOperator& schur = *sharedSchur;
    std::uint64_t traceSourceFingerprint =
        UINT64_C(1469598103934665603);
    for (const PortPatch& patch : patches) {
        hashValue(traceSourceFingerprint, patch.interfaceId);
        hashValue(traceSourceFingerprint, patch.sourceFingerprint);
    }
    result.model.traceSourceFingerprint = traceSourceFingerprint;

    result.model.ports.reserve(patches.size());
    result.interfaces.reserve(patches.size());
    for (const PortPatch& patch : patches) {
        OptimalPortInterfaceDiagnostics diagnostics;
        diagnostics.interfaceId = patch.interfaceId;
        diagnostics.targetRows = static_cast<int>(patch.target.size());

        const auto gramStart = Clock::now();
        std::vector<double> targetMass = traceWeights(
            dynamicModel, patch.target, traceMassDiagonal,
            penaltyMassDiagonal, options.innerProduct,
            &diagnostics.gramMinimum, &diagnostics.gramMaximum,
            &diagnostics.gramRegularization);
        std::vector<double> sourceMass;
        if (includeTraceSource && !patch.source.empty()) {
            sourceMass = traceWeights(
                dynamicModel, patch.source, traceMassDiagonal,
                penaltyMassDiagonal, options.innerProduct);
        }
        if (includeInputSource) {
            sourceMass.insert(
                sourceMass.end(),
                static_cast<std::size_t>(sourceChannels), 1.0);
        }
        if (includeBoundarySource) {
            sourceMass.insert(
                sourceMass.end(),
                static_cast<std::size_t>(boundaryChannels), 1.0);
        }
        if (includeHistorySource) {
            sourceMass.insert(
                sourceMass.end(),
                static_cast<std::size_t>(historyChannels), 1.0);
        }
        result.gramAssemblySeconds +=
            std::chrono::duration<double>(
                Clock::now() - gramStart).count();

        GeneralizedPatchTransferOperator transfer(
            schur, patch, sourceBlocks,
            includeTraceSource, includeInputSource,
            includeBoundarySource, includeHistorySource,
            options.innerSolver);
        diagnostics.sourceRows = transfer.sourceRows();
        diagnostics.traceSourceRows = transfer.traceRows();
        diagnostics.inputSourceRows = transfer.inputRows();
        diagnostics.boundarySourceRows = transfer.boundaryRows();
        diagnostics.historySourceRows = transfer.historyRows();
        bool pilotBlocked = false;
        auto blockPilot = [&](const std::string& status) {
            pilotBlocked = true;
            diagnostics.pilotStatus = status;
        };
        if (options.targetSolverPilotPreflight) {
            diagnostics.pilotStatus = "running";
            PatchInnerSolverStatistics current =
                transfer.statistics();
            if (current.setupSeconds
                    > options.maximumTargetSetupSeconds) {
                blockPilot("failed_setup_time");
            } else if (current.peakIncrementalMemoryBytes
                    > options.maximumIncrementalWorkspaceBytes) {
                blockPilot("failed_incremental_workspace");
            }
            if (!pilotBlocked) {
                const auto preflightStart = Clock::now();
                for (int probe = 0; probe < 2; ++probe) {
                    std::vector<double> rightHandSide(
                        patch.target.size(), 0.0);
                    for (std::size_t row = 0;
                         row < rightHandSide.size(); ++row) {
                        rightHandSide[row] =
                            std::sin((0.271 + 0.113 * probe)
                                * static_cast<double>(row + 1))
                            + 0.1 * std::cos(
                                0.397 * static_cast<double>(row + 1));
                    }
                    std::vector<double> solution;
                    transfer.solveTargetResponse(
                        rightHandSide, solution);
                }
                diagnostics.targetSolvePreflightSeconds =
                    std::chrono::duration<double>(
                        Clock::now() - preflightStart).count();
                current = transfer.statistics();
                const double meanSolve = current.solveCalls > 0
                    ? current.totalSolveSeconds / current.solveCalls : 0.0;
                if (meanSolve
                        > options.maximumMeanTargetSolveSeconds) {
                    blockPilot("failed_target_solve_time");
                } else if (current.maximumRelativeResidual
                        > options.maximumTargetResidual) {
                    blockPilot("failed_target_residual");
                }
            }
        }
        const auto operatorCheckStart = Clock::now();
        if (!pilotBlocked && transfer.sourceRows() > 0) {
            const int referenceColumns =
                std::min(4, transfer.sourceRows());
            std::vector<double> sourceProbe(
                static_cast<std::size_t>(transfer.sourceRows()), 0.0);
            for (int column = 0; column < referenceColumns; ++column) {
                sourceProbe[static_cast<std::size_t>(column)] =
                    0.25 * static_cast<double>(column + 1);
            }
            std::vector<double> matrixFreeImage;
            const auto transferApplyStart = Clock::now();
            transfer.apply(sourceProbe, matrixFreeImage);
            diagnostics.transferApplySeconds =
                std::chrono::duration<double>(
                    Clock::now() - transferApplyStart).count();
            if (!options.targetSolverPilotPreflight
                || options.columnConsistencyCheck) {
                std::vector<double> explicitImage(
                    patch.target.size(), 0.0);
                for (int column = 0; column < referenceColumns; ++column) {
                    std::vector<double> unit(
                        static_cast<std::size_t>(transfer.sourceRows()), 0.0);
                    unit[static_cast<std::size_t>(column)] = 1.0;
                    std::vector<double> columnImage;
                    transfer.apply(unit, columnImage);
                    for (std::size_t row = 0;
                         row < explicitImage.size(); ++row) {
                        explicitImage[row] +=
                            sourceProbe[static_cast<std::size_t>(column)]
                            * columnImage[row];
                    }
                }
                std::vector<double> difference = matrixFreeImage;
                for (std::size_t row = 0; row < difference.size(); ++row) {
                    difference[row] -= explicitImage[row];
                }
                diagnostics.explicitColumnReferenceError =
                    std::sqrt(std::max(0.0,
                        dot(difference.data(), difference.data(),
                            static_cast<int>(difference.size()))))
                    / std::max(1.0e-300,
                        std::sqrt(std::max(0.0,
                            dot(explicitImage.data(), explicitImage.data(),
                                static_cast<int>(explicitImage.size())))));
            }

            std::vector<double> targetProbe(
                patch.target.size(), 0.0);
            for (std::size_t row = 0; row < targetProbe.size(); ++row) {
                targetProbe[row] = std::sin(
                    0.3819660112501051 * static_cast<double>(row + 1));
            }
            std::vector<double> weightedTarget = targetProbe;
            for (std::size_t row = 0; row < weightedTarget.size(); ++row) {
                weightedTarget[row] *= targetMass[row];
            }
            std::vector<double> transposeImage;
            const auto transposeApplyStart = Clock::now();
            transfer.applyTranspose(
                weightedTarget, transposeImage);
            diagnostics.transposeApplySeconds =
                std::chrono::duration<double>(
                    Clock::now() - transposeApplyStart).count();
            std::vector<double> weightedAdjoint = transposeImage;
            for (std::size_t row = 0;
                 row < weightedAdjoint.size(); ++row) {
                weightedAdjoint[row] /= sourceMass[row];
            }
            const double left = massDot(
                targetProbe.data(), matrixFreeImage.data(), targetMass);
            const double right = massDot(
                weightedAdjoint.data(), sourceProbe.data(), sourceMass);
            diagnostics.adjointRelativeError =
                std::abs(left - right)
                / std::max({1.0e-300, std::abs(left), std::abs(right)});
        }
        diagnostics.operatorCheckSeconds =
            std::chrono::duration<double>(
                Clock::now() - operatorCheckStart).count();
        if (options.targetSolverPilotPreflight && !pilotBlocked) {
            const PatchInnerSolverStatistics current =
                transfer.statistics();
            if (diagnostics.operatorCheckSeconds
                    > options.maximumOperatorCheckSeconds) {
                blockPilot("failed_operator_adjoint_time");
            } else if (current.maximumRelativeResidual
                    > options.maximumTargetResidual) {
                blockPilot("failed_target_residual");
            } else if (diagnostics.adjointRelativeError
                    > options.maximumWeightedAdjointError) {
                blockPilot("failed_weighted_adjoint");
            } else if (options.columnConsistencyCheck
                    && diagnostics.explicitColumnReferenceError
                        > 1.0e-10) {
                blockPilot("failed_column_consistency");
            }
        }
        LocalPortBasis port;
        port.interfaceId = patch.interfaceId;
        port.leftSubdomain = patch.leftSubdomain;
        port.rightSubdomain = patch.rightSubdomain;
        port.rows = static_cast<int>(patch.target.size());
        port.interfaceIndices = patch.target;
        port.sourceIndices = patch.source;
        port.patchSubdomains = patch.patchSubdomains;
        port.targetFingerprint = patch.targetFingerprint;
        port.sourceFingerprint = patch.sourceFingerprint;
        port.traceSourceFingerprint = transfer.traceFingerprint();
        port.inputSourceFingerprint = transfer.inputFingerprint();
        port.boundarySourceFingerprint = transfer.boundaryFingerprint();
        port.historySourceFingerprint = transfer.historyFingerprint();
        port.traceSourceRows = transfer.traceRows();
        port.inputSourceRows = transfer.inputRows();
        port.boundarySourceRows = transfer.boundaryRows();
        port.historySourceRows = transfer.historyRows();
        port.templateId = patch.interfaceId;

        int targetRank = options.requestedRank;
        if (options.rankMode == "match-m7") {
            const auto found = rankMap.find(patch.interfaceId);
            if (found == rankMap.end()) {
                throw std::runtime_error(
                    "[Optimal port] M7 rank is missing for interface "
                    + std::to_string(patch.interfaceId) + '.');
            }
            targetRank = found->second;
        } else if (options.rankMode == "eigenvalue-tolerance") {
            targetRank = options.maximumRank;
        }
        targetRank = std::min(targetRank, port.rows);
        if (targetRank <= 0) {
            throw std::runtime_error(
                "[Optimal port] Selected target rank is zero.");
        }

        const auto mandatoryStart = Clock::now();
        const int firstCoordinate = includeConstant ? 0 : 1;
        const int lastCoordinate = includeGeometry ? 3 : 0;
        for (int coordinate = firstCoordinate;
             coordinate <= lastCoordinate; ++coordinate) {
            std::vector<double> candidate(
                static_cast<std::size_t>(port.rows), 1.0);
            for (int row = 0; row < port.rows; ++row) {
                const int gamma =
                    port.interfaceIndices[static_cast<std::size_t>(row)];
                const int global = dynamicModel.interfaceGlobalDofs[
                    static_cast<std::size_t>(gamma)];
                const Vec3& point =
                    mesh.nodes[static_cast<std::size_t>(global)].p;
                if (coordinate == 1) {
                    candidate[static_cast<std::size_t>(row)] = point.x;
                } else if (coordinate == 2) {
                    candidate[static_cast<std::size_t>(row)] = point.y;
                } else if (coordinate == 3) {
                    candidate[static_cast<std::size_t>(row)] = point.z;
                }
            }
            if (appendMassMode(
                    port, std::move(candidate), targetMass, targetRank,
                    options.relativeDeflationTolerance)) {
                ++port.mandatoryModes;
            }
        }
        for (int channel = 0;
             includeParticular && channel < sourceChannels
                 && port.rank < targetRank;
             ++channel) {
            std::vector<double> rightHandSide(
                static_cast<std::size_t>(port.rows), 0.0);
            for (int row = 0; row < port.rows; ++row) {
                rightHandSide[static_cast<std::size_t>(row)] =
                    condensedInputs[static_cast<std::size_t>(
                        channel * dynamicModel.interfaceDofs
                        + port.interfaceIndices[
                            static_cast<std::size_t>(row)])];
            }
            std::vector<double> candidate;
            transfer.solveTargetResponse(rightHandSide, candidate);
            if (appendMassMode(
                    port, std::move(candidate), targetMass, targetRank,
                    options.relativeDeflationTolerance)) {
                ++port.mandatoryModes;
            }
        }
        result.mandatoryModeSeconds +=
            std::chrono::duration<double>(
                Clock::now() - mandatoryStart).count();

        if (options.requestedTransferRank >= 0) {
            targetRank = std::min(
                port.rows,
                port.rank + options.requestedTransferRank);
        }
        if (options.targetSolverPilotPreflight
            && options.maximumPilotSeconds > 0.0
            && std::chrono::duration<double>(
                Clock::now() - totalStart).count()
                >= options.maximumPilotSeconds) {
            blockPilot("failed_total_time");
        }
        MatrixFreeEigenResult eigen;
        if (!pilotBlocked && includeTransfer
            && transfer.sourceRows() > 0 && port.rank < targetRank) {
            MatrixFreeEigenOptions eigenOptions;
            // One look-ahead Ritz vector materially improves subspace
            // convergence and supplies the fixed-rank truncation indicator.
            // It is never counted as a requested/converged Transfer mode.
            constexpr int lookahead = 1;
            const int maximumTransfer = std::max(
                1, std::min(port.rows - port.rank,
                    targetRank - port.rank + lookahead));
            eigenOptions.requestedEigenpairs = maximumTransfer;
            eigenOptions.oversamplingVectors = 4;
            eigenOptions.maximumIterations =
                options.eigensolverMaximumIterations;
            eigenOptions.relativeTolerance =
                options.eigensolverTolerance;
            eigenOptions.deflationTolerance =
                options.relativeDeflationTolerance;
            if (options.maximumPilotSeconds > 0.0) {
                const double usedSeconds =
                    std::chrono::duration<double>(
                        Clock::now() - totalStart).count();
                eigenOptions.maximumWallSeconds =
                    std::max(1.0e-12,
                        options.maximumPilotSeconds - usedSeconds);
            }
            const auto eigenStart = Clock::now();
            eigen = solveLargestSymmetricEigenpairs(
                port.rows,
                [&](const std::vector<double>& inputVector,
                    std::vector<double>& outputVector) {
                    std::vector<double> weightedTarget(
                        inputVector.size(), 0.0);
                    for (std::size_t row = 0;
                         row < inputVector.size(); ++row) {
                        weightedTarget[row] =
                            std::sqrt(targetMass[row])
                            * inputVector[row];
                    }
                    std::vector<double> sourceVector;
                    transfer.applyTranspose(
                        weightedTarget, sourceVector);
                    for (std::size_t row = 0;
                         row < sourceVector.size(); ++row) {
                        sourceVector[row] /= sourceMass[row];
                    }
                    std::vector<double> targetVector;
                    transfer.apply(sourceVector, targetVector);
                    outputVector.resize(targetVector.size());
                    for (std::size_t row = 0;
                         row < targetVector.size(); ++row) {
                        outputVector[row] =
                            std::sqrt(targetMass[row])
                            * targetVector[row];
                    }
                },
                eigenOptions);
            result.eigenSolveSeconds +=
                std::chrono::duration<double>(
                    Clock::now() - eigenStart).count();
        } else {
            eigen.dimension = port.rows;
            eigen.converged = !pilotBlocked;
            eigen.status = pilotBlocked
                ? diagnostics.pilotStatus
                : (!includeTransfer
                ? "transfer_disabled_by_ablation"
                : (transfer.sourceRows() == 0
                    ? "empty_source" : "rank_filled_by_mandatory"));
        }

        int transferLimit = std::max(0, targetRank - port.rank);
        port.requestedTransferRank = transferLimit;
        if (options.rankMode == "eigenvalue-tolerance"
            && !eigen.eigenvalues.empty()) {
            transferLimit = 0;
            const double leading = std::max(
                std::numeric_limits<double>::min(),
                eigen.eigenvalues.front());
            for (int retained = 0;
                 retained < static_cast<int>(eigen.eigenvalues.size())
                    && port.rank + retained < options.maximumRank;
                 ++retained) {
                transferLimit = retained + 1;
                const int next = retained + 1;
                const double indicator =
                    next < static_cast<int>(eigen.eigenvalues.size())
                    ? std::sqrt(std::max(0.0,
                        eigen.eigenvalues[
                            static_cast<std::size_t>(next)] / leading))
                    : 0.0;
                if (port.rank + transferLimit >= options.minimumRank
                    && indicator <= options.eigenvalueTolerance) {
                    break;
                }
            }
        }

        const auto orthogonalizationStart = Clock::now();
        for (int selected = 0;
             selected < static_cast<int>(eigen.eigenvalues.size())
                && port.spectralModes < transferLimit
                && port.rank < targetRank;
             ++selected) {
            if (selected >= static_cast<int>(eigen.residuals.size())
                || eigen.residuals[static_cast<std::size_t>(selected)]
                    > options.eigensolverTolerance) {
                continue;
            }
            std::vector<double> candidate(
                static_cast<std::size_t>(port.rows), 0.0);
            for (int row = 0; row < port.rows; ++row) {
                candidate[static_cast<std::size_t>(row)] =
                    eigen.eigenvectors[static_cast<std::size_t>(
                        selected * port.rows + row)]
                    / std::sqrt(targetMass[static_cast<std::size_t>(row)]);
            }
            if (appendMassMode(
                    port, std::move(candidate), targetMass, targetRank,
                    options.relativeDeflationTolerance)) {
                ++port.spectralModes;
                port.spectralValues.push_back(
                    eigen.eigenvalues[static_cast<std::size_t>(selected)]);
                port.spectralResiduals.push_back(
                    eigen.residuals[static_cast<std::size_t>(selected)]);
            }
        }
        result.orthogonalizationSeconds +=
            std::chrono::duration<double>(
                Clock::now() - orthogonalizationStart).count();

        if (!eigen.eigenvalues.empty()) {
            const double leading = std::max(
                std::numeric_limits<double>::min(),
                eigen.eigenvalues.front());
            const int next = port.spectralModes;
            port.transferIndicator =
                next < static_cast<int>(eigen.eigenvalues.size())
                ? std::sqrt(std::max(0.0,
                    eigen.eigenvalues[static_cast<std::size_t>(next)]
                        / leading))
                : 0.0;
        }
        port.orthogonalityError =
            portOrthogonalityError(port, targetMass);
        std::uint64_t basisFingerprint =
            UINT64_C(1469598103934665603);
        hashValue(basisFingerprint, port.interfaceId);
        hashValue(basisFingerprint, port.rank);
        hashValue(basisFingerprint, port.targetFingerprint);
        hashValue(basisFingerprint, port.sourceFingerprint);
        hashValue(basisFingerprint, port.inputSourceFingerprint);
        hashValue(basisFingerprint, port.boundarySourceFingerprint);
        hashValue(basisFingerprint, port.historySourceFingerprint);
        for (double value : port.basis) {
            hashValue(basisFingerprint, value);
        }
        port.fingerprint = basisFingerprint;

        diagnostics.mandatoryRank = port.mandatoryModes;
        diagnostics.requestedTransferRank =
            port.requestedTransferRank;
        diagnostics.convergedTransferRank = port.spectralModes;
        diagnostics.totalPortRank = port.rank;
        diagnostics.eigenIterations = eigen.iterations;
        diagnostics.eigenOperatorApplies = eigen.operatorApplies;
        diagnostics.eigenConverged =
            port.spectralModes == port.requestedTransferRank;
        diagnostics.eigenStatus = pilotBlocked
            ? diagnostics.pilotStatus
            : (diagnostics.eigenConverged
                ? "success"
                : (port.spectralModes > 0
                    ? "partial_convergence"
                    : eigen.status));
        diagnostics.transferIndicator = port.transferIndicator;
        diagnostics.eigenvalues = port.spectralValues;
        const int residualCount = std::min(
            port.requestedTransferRank,
            static_cast<int>(eigen.residuals.size()));
        for (int selected = 0; selected < residualCount; ++selected) {
            diagnostics.maximumEigenpairResidual = std::max(
                diagnostics.maximumEigenpairResidual,
                eigen.residuals[static_cast<std::size_t>(selected)]);
        }
        diagnostics.innerSolver = transfer.statistics();
        if (options.targetSolverPilotPreflight
            && diagnostics.pilotStatus == "running") {
            diagnostics.pilotStatus =
                eigen.status == "time_limit"
                ? "failed_total_time"
                : "passed_preflight";
        }
        const std::size_t target = static_cast<std::size_t>(port.rows);
        const std::size_t block = static_cast<std::size_t>(
            std::min(port.rows,
                std::max(1, port.requestedTransferRank + 4)));
        diagnostics.peakWorkspaceBytes =
            std::max(
                diagnostics.innerSolver.factorBytes,
                diagnostics.innerSolver.peakIncrementalMemoryBytes)
            + 4 * target * block * sizeof(double)
            + static_cast<std::size_t>(
                5 * transfer.targetRows()
                + 3 * transfer.sourceRows()) * sizeof(double);
        diagnostics.eigensolverWorkspaceBytes =
            4 * target * block * sizeof(double)
            + static_cast<std::size_t>(
                5 * transfer.targetRows()
                + 3 * transfer.sourceRows()) * sizeof(double);
        diagnostics.finalBasisStorageBytes =
            port.basis.capacity() * sizeof(double)
            + port.interfaceIndices.capacity() * sizeof(int)
            + port.sourceIndices.capacity() * sizeof(int)
            + port.patchSubdomains.capacity() * sizeof(int)
            + port.spectralValues.capacity() * sizeof(double)
            + port.spectralResiduals.capacity() * sizeof(double);

        result.model.reducedInterfaceDofs += port.rank;
        result.model.modelBytes +=
            port.interfaceIndices.capacity() * sizeof(int)
            + port.sourceIndices.capacity() * sizeof(int)
            + port.patchSubdomains.capacity() * sizeof(int)
            + port.basis.capacity() * sizeof(double)
            + port.spectralValues.capacity() * sizeof(double)
            + port.spectralResiduals.capacity() * sizeof(double);
        result.model.ports.push_back(std::move(port));
        result.interfaces.push_back(std::move(diagnostics));
    }
    result.model.snapshotSeconds = 0.0;
    result.model.basisSeconds =
        std::chrono::duration<double>(Clock::now() - totalStart).count();
    result.model.modelBytes +=
        result.model.interfaceGlobalDofs.capacity() * sizeof(int);
    result.totalSeconds = result.model.basisSeconds;
    return result;
}

namespace {

struct HistoryCompressionSelection {
    std::vector<int> selectedColumns;
    int activeColumns = 0;
    double relativeError = 0.0;
    double seconds = 0.0;
    std::size_t workspaceBytes = 0;
    std::uint64_t fingerprint = UINT64_C(1469598103934665603);
};

HistoryCompressionSelection selectHistoryColumns(
    const std::vector<double>& globalColumns,
    int globalRows,
    int channelCount,
    const std::vector<int>& target,
    const std::vector<double>& targetMass,
    const ResidualKrylovPortOptions& options)
{
    const auto start = Clock::now();
    HistoryCompressionSelection result;
    if (channelCount == 0) {
        result.seconds =
            std::chrono::duration<double>(Clock::now() - start).count();
        return result;
    }
    if (globalRows <= 0
        || globalColumns.size() != static_cast<std::size_t>(
            globalRows) * static_cast<std::size_t>(channelCount)
        || target.size() != targetMass.size()) {
        throw std::runtime_error(
            "[History compression] Operator block dimensions are invalid.");
    }

    std::vector<double> normsSquared(
        static_cast<std::size_t>(channelCount), 0.0);
    double maximumNormSquared = 0.0;
    for (int column = 0; column < channelCount; ++column) {
        long double norm = 0.0L;
        const double* values = globalColumns.data()
            + static_cast<std::size_t>(column)
                * static_cast<std::size_t>(globalRows);
        for (std::size_t row = 0; row < target.size(); ++row) {
            const double value =
                values[target[row]];
            norm += static_cast<long double>(value) * value
                / targetMass[row];
        }
        normsSquared[static_cast<std::size_t>(column)] =
            std::max(0.0, static_cast<double>(norm));
        maximumNormSquared = std::max(
            maximumNormSquared,
            normsSquared[static_cast<std::size_t>(column)]);
    }
    const double activeThreshold =
        options.historyCompressionTolerance
        * options.historyCompressionTolerance
        * maximumNormSquared;
    std::vector<int> active;
    active.reserve(static_cast<std::size_t>(channelCount));
    for (int column = 0; column < channelCount; ++column) {
        if (normsSquared[static_cast<std::size_t>(column)]
                > activeThreshold
            && normsSquared[static_cast<std::size_t>(column)] > 0.0) {
            active.push_back(column);
        }
    }
    result.activeColumns = static_cast<int>(active.size());

    if (options.historyCompressionMethod == "none") {
        result.selectedColumns.resize(
            static_cast<std::size_t>(channelCount));
        std::iota(
            result.selectedColumns.begin(),
            result.selectedColumns.end(), 0);
        for (int column : result.selectedColumns) {
            hashValue(result.fingerprint, column);
        }
        result.workspaceBytes =
            normsSquared.capacity() * sizeof(double)
            + active.capacity() * sizeof(int)
            + result.selectedColumns.capacity() * sizeof(int);
        result.seconds =
            std::chrono::duration<double>(Clock::now() - start).count();
        return result;
    }
    if (active.empty()) {
        result.relativeError = 0.0;
        result.workspaceBytes =
            normsSquared.capacity() * sizeof(double)
            + active.capacity() * sizeof(int);
        result.seconds =
            std::chrono::duration<double>(Clock::now() - start).count();
        return result;
    }

    const int activeCount = static_cast<int>(active.size());
    const int requested = std::min(
        options.historyCompressionRank, activeCount);
    std::vector<double> gram(
        static_cast<std::size_t>(activeCount)
            * static_cast<std::size_t>(activeCount), 0.0);
    constexpr int rowBlockSize = 128;
    std::vector<double> scaledBlock(static_cast<std::size_t>(
        rowBlockSize * activeCount), 0.0);
    for (int begin = 0;
         begin < static_cast<int>(target.size());
         begin += rowBlockSize) {
        const int rows = std::min(
            rowBlockSize,
            static_cast<int>(target.size()) - begin);
        for (int column = 0; column < activeCount; ++column) {
            const int sourceColumn =
                active[static_cast<std::size_t>(column)];
            const double* values = globalColumns.data()
                + static_cast<std::size_t>(sourceColumn)
                    * static_cast<std::size_t>(globalRows);
            double* scaled = scaledBlock.data()
                + static_cast<std::size_t>(column * rowBlockSize);
            for (int localRow = 0; localRow < rows; ++localRow) {
                const int targetRow = begin + localRow;
                scaled[localRow] =
                    values[target[static_cast<std::size_t>(targetRow)]]
                    / std::sqrt(
                        targetMass[static_cast<std::size_t>(targetRow)]);
            }
        }
#ifdef USE_MKL_PARDISO
        cblas_dsyrk(
            CblasColMajor, CblasUpper, CblasTrans,
            activeCount, rows, 1.0, scaledBlock.data(),
            rowBlockSize, 1.0, gram.data(), activeCount);
#else
        for (int right = 0; right < activeCount; ++right) {
            const double* rightValues = scaledBlock.data()
                + static_cast<std::size_t>(right * rowBlockSize);
            for (int left = 0; left <= right; ++left) {
                const double* leftValues = scaledBlock.data()
                    + static_cast<std::size_t>(left * rowBlockSize);
                long double value = 0.0L;
                for (int row = 0; row < rows; ++row) {
                    value += static_cast<long double>(leftValues[row])
                        * rightValues[row];
                }
                gram[static_cast<std::size_t>(
                    left + right * activeCount)] +=
                    static_cast<double>(value);
            }
        }
#endif
    }
    for (int column = 0; column < activeCount; ++column) {
        for (int row = column + 1; row < activeCount; ++row) {
            gram[static_cast<std::size_t>(
                row + column * activeCount)] =
                gram[static_cast<std::size_t>(
                    column + row * activeCount)];
        }
    }

    std::vector<double> residualDiagonal(
        static_cast<std::size_t>(activeCount), 0.0);
    double initialMaximum = 0.0;
    for (int column = 0; column < activeCount; ++column) {
        residualDiagonal[static_cast<std::size_t>(column)] =
            std::max(0.0, gram[static_cast<std::size_t>(
                column + column * activeCount)]);
        initialMaximum = std::max(
            initialMaximum,
            residualDiagonal[static_cast<std::size_t>(column)]);
    }
    std::vector<double> cholesky(
        static_cast<std::size_t>(activeCount)
            * static_cast<std::size_t>(requested), 0.0);
    std::vector<unsigned char> selected(
        static_cast<std::size_t>(activeCount), 0);
    // The weighted Gram matrix is formed in double precision.  Do not let a
    // user tolerance below the reliable Gram-diagonal resolution turn
    // roundoff into artificial rank.
    const double pivotThreshold = std::max(
        options.historyCompressionTolerance
            * options.historyCompressionTolerance,
        64.0 * std::numeric_limits<double>::epsilon())
        * initialMaximum;
    result.selectedColumns.reserve(static_cast<std::size_t>(requested));
    for (int mode = 0; mode < requested; ++mode) {
        int pivot = -1;
        double maximumResidual = -1.0;
        for (int candidate = 0;
             candidate < activeCount; ++candidate) {
            if (selected[static_cast<std::size_t>(candidate)] != 0) {
                continue;
            }
            const double value =
                residualDiagonal[static_cast<std::size_t>(candidate)];
            if (value > maximumResidual
                || (value == maximumResidual && pivot >= 0
                    && active[static_cast<std::size_t>(candidate)]
                        < active[static_cast<std::size_t>(pivot)])) {
                maximumResidual = value;
                pivot = candidate;
            }
        }
        if (pivot < 0 || !(maximumResidual > pivotThreshold)) {
            break;
        }
        selected[static_cast<std::size_t>(pivot)] = 1;
        result.selectedColumns.push_back(
            active[static_cast<std::size_t>(pivot)]);
        const double diagonal = std::sqrt(maximumResidual);
        cholesky[static_cast<std::size_t>(
            pivot + mode * activeCount)] = diagonal;
        for (int row = 0; row < activeCount; ++row) {
            if (selected[static_cast<std::size_t>(row)] != 0) {
                continue;
            }
            double value = gram[static_cast<std::size_t>(
                row + pivot * activeCount)];
            for (int previous = 0; previous < mode; ++previous) {
                value -= cholesky[static_cast<std::size_t>(
                    row + previous * activeCount)]
                    * cholesky[static_cast<std::size_t>(
                        pivot + previous * activeCount)];
            }
            const double factor = value / diagonal;
            cholesky[static_cast<std::size_t>(
                row + mode * activeCount)] = factor;
            residualDiagonal[static_cast<std::size_t>(row)] =
                std::max(
                    0.0,
                    residualDiagonal[static_cast<std::size_t>(row)]
                        - factor * factor);
        }
        residualDiagonal[static_cast<std::size_t>(pivot)] = 0.0;
    }
    double maximumRemaining = 0.0;
    for (double value : residualDiagonal) {
        maximumRemaining = std::max(maximumRemaining, value);
    }
    result.relativeError = initialMaximum > 0.0
        ? std::sqrt(maximumRemaining / initialMaximum) : 0.0;
    hashValue(result.fingerprint, channelCount);
    hashValue(result.fingerprint, result.activeColumns);
    hashValue(result.fingerprint, options.historyCompressionRank);
    for (int column : result.selectedColumns) {
        hashValue(result.fingerprint, column);
    }
    result.workspaceBytes =
        normsSquared.capacity() * sizeof(double)
        + active.capacity() * sizeof(int)
        + gram.capacity() * sizeof(double)
        + scaledBlock.capacity() * sizeof(double)
        + residualDiagonal.capacity() * sizeof(double)
        + cholesky.capacity() * sizeof(double)
        + selected.capacity() * sizeof(unsigned char)
        + result.selectedColumns.capacity() * sizeof(int);
    result.seconds =
        std::chrono::duration<double>(Clock::now() - start).count();
    return result;
}

} // namespace

ResidualKrylovBuildResult buildResidualKrylovPortModel(
    const Mesh& mesh,
    const ddm_schur::InterfacePartition& partition,
    const local::Model& dynamicModel,
    const std::vector<double>& traceMassDiagonal,
    const std::vector<double>& penaltyMassDiagonal,
    const std::vector<double>& input,
    int sourceChannels,
    const std::vector<double>& boundaryLoad,
    const std::vector<double>& condensedHistory,
    int historyChannels,
    const ResidualKrylovPortOptions& options,
    const ReducedDynamicSchurOperator* sharedSchur,
    const LocalPortModel* initialTransferModel)
{
    const auto totalStart = Clock::now();
    if ((options.basisMethod != "mandatory-only"
            && options.basisMethod != "residual-krylov"
            && options.basisMethod != "hybrid-randomized")
        || options.maximumEnrichmentRank < 0
        || options.maximumSweeps <= 0 || options.blockSize <= 0
        || !(options.residualTolerance > 0.0)
        || !(options.relativeDeflationTolerance > 0.0)
        || (options.historyCompressionMethod != "none"
            && options.historyCompressionMethod
                != "deterministic-rrqr")
        || (options.historyCompressionMethod
                == "deterministic-rrqr"
            && (options.historyCompressionRank <= 0
                || !(options.historyCompressionTolerance > 0.0)))
        || (options.probeMode != "operator-geometry"
            && options.probeMode != "particular-only")
        || options.innerSolver.innerSolver != "woodbury-exact") {
        throw std::runtime_error(
            "[Residual Krylov] Invalid port-space options.");
    }
    if (historyChannels < 0
        || condensedHistory.size() != static_cast<std::size_t>(
            dynamicModel.interfaceDofs * historyChannels)) {
        throw std::runtime_error(
            "[Residual Krylov] Condensed history dimensions are invalid.");
    }

    ResidualKrylovBuildResult result;
    result.model.formatVersion = 7;
    result.model.basisMethod = options.basisMethod;
    result.model.methodDescription =
        options.basisMethod == "hybrid-randomized"
        ? "Local Block Arnoldi with Operator-Informed Mandatory, "
          "Randomized Transfer, and Schur-Residual Port Enrichment"
        : "Local Block Arnoldi with Operator-Informed Port Space and "
          "Schur-Residual Krylov Enrichment";
    result.model.ablationMode = options.basisMethod;
    result.model.sourceMode =
        initialTransferModel == nullptr
        ? options.probeMode : initialTransferModel->sourceMode;
    result.model.innerProduct = options.innerProduct;
    result.model.rankMode =
        initialTransferModel == nullptr
        ? "residual-tolerance" : "hybrid-randomized-residual";
    result.model.oversamplingLayers = options.oversamplingLayers;
    result.model.requestedRank =
        initialTransferModel == nullptr
        ? options.maximumEnrichmentRank
        : initialTransferModel->requestedRank;
    result.model.minimumRank = 0;
    result.model.maximumRank = options.maximumEnrichmentRank;
    result.model.eigenvalueTolerance = options.residualTolerance;
    result.model.eigensolverTolerance = 0.0;
    result.model.eigensolverMaximumIterations = options.maximumSweeps;
    result.model.relativeDeflationTolerance =
        options.relativeDeflationTolerance;
    result.model.innerSolver = options.innerSolver.innerSolver;
    result.model.innerSolverTolerance =
        options.innerSolver.relativeTolerance;
    result.model.innerSolverMaximumIterations =
        options.innerSolver.maximumIterations;
    result.model.historyCompressionMethod =
        options.historyCompressionMethod;
    result.model.historyCompressionRank =
        options.historyCompressionRank;
    result.model.historyCompressionTolerance =
        options.historyCompressionTolerance;
    result.model.fullInterfaceDofs = dynamicModel.interfaceDofs;
    result.model.interfaceGlobalDofs = dynamicModel.interfaceGlobalDofs;
    result.model.generalizedInputFingerprint = fingerprintVector(input);
    result.model.generalizedBoundaryFingerprint =
        fingerprintVector(boundaryLoad);
    result.model.generalizedHistoryFingerprint =
        fingerprintVector(condensedHistory);
    if (initialTransferModel != nullptr) {
        if (options.basisMethod != "hybrid-randomized"
            || initialTransferModel->basisMethod
                != "randomized-transfer"
            || initialTransferModel->innerProduct
                != options.innerProduct
            || initialTransferModel->oversamplingLayers
                != options.oversamplingLayers
            || initialTransferModel->fullInterfaceDofs
                != dynamicModel.interfaceDofs
            || initialTransferModel->interfaceGlobalDofs
                != dynamicModel.interfaceGlobalDofs) {
            throw std::runtime_error(
                "[Hybrid port] Initial randomized-transfer model "
                "does not match the target operator or metric.");
        }
    } else if (options.basisMethod == "hybrid-randomized") {
        throw std::runtime_error(
            "[Hybrid port] A randomized-transfer model is required.");
    }

    std::vector<PortPatch> allPatches = buildOptimalPortPatches(
        mesh, partition, options.oversamplingLayers);
    std::map<int, int> physicalOwner;
    for (const PortPatch& physical : allPatches) {
        for (int gamma : physical.target) {
            const auto insertion =
                physicalOwner.emplace(gamma, physical.interfaceId);
            if (!insertion.second
                && insertion.first->second != physical.interfaceId) {
                throw std::runtime_error(
                    "[Residual Krylov] Target DOF has multiple physical-port owners.");
            }
        }
    }
    std::vector<PortPatch> patches = allPatches;
    if (!options.selectedInterfaceIds.empty()) {
        const std::set<int> selected(
            options.selectedInterfaceIds.begin(),
            options.selectedInterfaceIds.end());
        patches.erase(std::remove_if(
            patches.begin(), patches.end(),
            [&](const PortPatch& patch) {
                return selected.count(patch.interfaceId) == 0;
            }), patches.end());
        if (patches.size() != selected.size()) {
            throw std::runtime_error(
                "[Residual Krylov] A selected interface is absent.");
        }
    }

    GeneralizedTransferSourceBlocks residualSourceBlocks =
        buildGeneralizedTransferSourceBlocks(
            dynamicModel, input, sourceChannels,
            boundaryLoad, {}, 0);
    const std::vector<double>& condensedInputs =
        residualSourceBlocks.input;
    const int boundaryChannels =
        residualSourceBlocks.boundaryChannels;
    const std::vector<double>& condensedBoundary =
        residualSourceBlocks.boundary;
    std::unique_ptr<ReducedDynamicSchurOperator> ownedSchur;
    if (sharedSchur == nullptr) {
        ownedSchur =
            std::make_unique<ReducedDynamicSchurOperator>(dynamicModel);
        sharedSchur = ownedSchur.get();
    }
    std::uint64_t traceFingerprint =
        UINT64_C(1469598103934665603);
    for (const PortPatch& patch : allPatches) {
        hashValue(traceFingerprint, patch.interfaceId);
        hashValue(traceFingerprint, patch.sourceFingerprint);
    }
    result.model.traceSourceFingerprint = traceFingerprint;

    auto inverseMassNorm = [](const std::vector<double>& vector,
                              const std::vector<double>& mass) {
        long double value = 0.0L;
        for (std::size_t row = 0; row < vector.size(); ++row) {
            value += static_cast<long double>(vector[row]) * vector[row]
                / mass[row];
        }
        return std::sqrt(std::max(0.0, static_cast<double>(value)));
    };
    auto inverseMassDot = [](const std::vector<double>& left,
                             const std::vector<double>& right,
                             const std::vector<double>& mass) {
        long double value = 0.0L;
        for (std::size_t row = 0; row < left.size(); ++row) {
            value += static_cast<long double>(left[row]) * right[row]
                / mass[row];
        }
        return static_cast<double>(value);
    };
    result.model.ports.reserve(patches.size());
    result.interfaces.reserve(patches.size());
    for (const PortPatch& patch : patches) {
        const auto interfaceStart = Clock::now();
        ResidualKrylovInterfaceDiagnostics diagnostics;
        diagnostics.interfaceId = patch.interfaceId;
        diagnostics.targetRows = static_cast<int>(patch.target.size());
        diagnostics.sourceRows = static_cast<int>(patch.source.size());
        diagnostics.requestedEnrichmentRank =
            options.basisMethod == "mandatory-only"
            ? 0 : options.maximumEnrichmentRank;

        const std::vector<double> targetMass = traceWeights(
            dynamicModel, patch.target, traceMassDiagonal,
            penaltyMassDiagonal, options.innerProduct);
        const HistoryCompressionSelection historySelection =
            selectHistoryColumns(
                condensedHistory, dynamicModel.interfaceDofs,
                historyChannels, patch.target, targetMass, options);
        diagnostics.rawHistoryChannels = historyChannels;
        diagnostics.activeHistoryChannels =
            historySelection.activeColumns;
        diagnostics.requestedHistoryRank =
            options.historyCompressionMethod == "none"
            ? historyChannels : options.historyCompressionRank;
        diagnostics.compressedHistoryRank =
            static_cast<int>(historySelection.selectedColumns.size());
        diagnostics.deflatedHistoryChannels =
            historyChannels - diagnostics.compressedHistoryRank;
        diagnostics.historyTargetRightHandSides =
            diagnostics.compressedHistoryRank;
        diagnostics.historyCompressionRelativeError =
            historySelection.relativeError;
        diagnostics.historyCompressionSeconds =
            historySelection.seconds;
        diagnostics.historyCompressionWorkspaceBytes =
            historySelection.workspaceBytes;
        diagnostics.historyCompressionFingerprint =
            historySelection.fingerprint;
        diagnostics.historyCompressionMethod =
            options.historyCompressionMethod;
        PatchTransferOperator targetSolver(
            *sharedSchur, patch, options.innerSolver);
        if (options.basisMethod == "hybrid-randomized") {
            std::cout << "[Hybrid port] interface="
                << patch.interfaceId << ", target="
                << patch.target.size() << ", source="
                << patch.source.size()
                << ", target solver ready\n" << std::flush;
        }
        LocalPortBasis port;
        port.interfaceId = patch.interfaceId;
        port.leftSubdomain = patch.leftSubdomain;
        port.rightSubdomain = patch.rightSubdomain;
        port.rows = diagnostics.targetRows;
        port.interfaceIndices = patch.target;
        port.sourceIndices = patch.source;
        port.patchSubdomains = patch.patchSubdomains;
        port.targetFingerprint = patch.targetFingerprint;
        port.sourceFingerprint = patch.sourceFingerprint;
        port.traceSourceFingerprint = patch.sourceFingerprint;
        port.inputSourceFingerprint =
            result.model.generalizedInputFingerprint;
        port.boundarySourceFingerprint =
            result.model.generalizedBoundaryFingerprint;
        port.historySourceFingerprint =
            result.model.generalizedHistoryFingerprint;
        hashValue(
            port.historySourceFingerprint,
            historySelection.fingerprint);
        port.traceSourceRows = diagnostics.sourceRows;
        port.inputSourceRows = sourceChannels;
        port.boundarySourceRows = boundaryChannels;
        port.historySourceRows =
            diagnostics.compressedHistoryRank;
        port.templateId = patch.interfaceId;
        port.requestedTransferRank =
            diagnostics.requestedEnrichmentRank;
        const int initialTransferRank =
            initialTransferModel == nullptr
            ? 0 : initialTransferModel->requestedRank;
        const int expectedPortColumns = std::min(
            port.rows,
            targetSolver.statistics().reducedCorrectionRank
                + 4 + initialTransferRank
                + diagnostics.requestedEnrichmentRank);
        port.basis.reserve(static_cast<std::size_t>(
            port.rows * expectedPortColumns));

        const auto mandatoryStart = Clock::now();
        for (int coordinate = 0; coordinate <= 3; ++coordinate) {
            std::vector<double> candidate(patch.target.size(), 1.0);
            for (std::size_t row = 0; row < patch.target.size(); ++row) {
                const int global = dynamicModel.interfaceGlobalDofs[
                    static_cast<std::size_t>(patch.target[row])];
                const Vec3& point =
                    mesh.nodes[static_cast<std::size_t>(global)].p;
                if (coordinate == 1) candidate[row] = point.x;
                if (coordinate == 2) candidate[row] = point.y;
                if (coordinate == 3) candidate[row] = point.z;
            }
            if (appendMassMode(
                    port, std::move(candidate), targetMass, port.rows,
                    options.relativeDeflationTolerance)) {
                if (coordinate == 0) {
                    ++diagnostics.constantRank;
                } else {
                    ++diagnostics.geometryRank;
                }
            }
        }
        std::size_t mandatoryTransientWorkspacePeak = 0;
        auto appendParticularFamily =
            [&](const std::vector<double>& globalColumns,
                int channelCount, int& accepted,
                const std::vector<int>* channelMap = nullptr) {
                if (channelMap != nullptr
                    && channelMap->size()
                        != static_cast<std::size_t>(channelCount)) {
                    throw std::runtime_error(
                        "[History compression] Selected channel map "
                        "has the wrong dimension.");
                }
                constexpr int particularBlockSize = 64;
                for (int begin = 0;
                     begin < channelCount;
                     begin += particularBlockSize) {
                    const int count = std::min(
                        particularBlockSize,
                        channelCount - begin);
                    std::vector<double> block(
                        static_cast<std::size_t>(
                            count * port.rows), 0.0);
                    for (int column = 0; column < count; ++column) {
                        const int sourceColumn = channelMap == nullptr
                            ? begin + column
                            : (*channelMap)[static_cast<std::size_t>(
                                begin + column)];
                        for (int row = 0; row < port.rows; ++row) {
                            block[static_cast<std::size_t>(
                                column * port.rows + row)] =
                                globalColumns[static_cast<std::size_t>(
                                    sourceColumn)
                                        * static_cast<std::size_t>(
                                            dynamicModel.interfaceDofs)
                                    + static_cast<std::size_t>(
                                        patch.target[
                                            static_cast<std::size_t>(
                                                row)])];
                        }
                    }
                    std::vector<double> particulars;
                    targetSolver.solveTargetResponses(
                        block, count, particulars);
                    mandatoryTransientWorkspacePeak = std::max(
                        mandatoryTransientWorkspacePeak,
                        (block.capacity() + particulars.capacity())
                            * sizeof(double));
#ifdef USE_MKL_PARDISO
                    // Weighted two-pass block MGS against the existing port.
                    // This is algebraically the same deterministic
                    // orthogonalization as appendMassMode, but turns the
                    // O(rows*rank*block) work into BLAS-3 operations.
                    const int baseRank = port.rank;
                    std::vector<double> originalNorms(
                        static_cast<std::size_t>(count), 0.0);
                    for (int column = 0; column < count; ++column) {
                        originalNorms[static_cast<std::size_t>(column)] =
                            std::sqrt(std::max(
                                0.0, massDot(
                                    particulars.data()
                                        + static_cast<std::size_t>(
                                            column * port.rows),
                                    particulars.data()
                                        + static_cast<std::size_t>(
                                            column * port.rows),
                                    targetMass)));
                    }
                    if (baseRank > 0) {
                        std::vector<double> weighted =
                            particulars;
                        std::vector<double> coefficients(
                            static_cast<std::size_t>(
                                baseRank * count), 0.0);
                        mandatoryTransientWorkspacePeak = std::max(
                            mandatoryTransientWorkspacePeak,
                            (block.capacity()
                                + particulars.capacity()
                                + weighted.capacity()
                                + coefficients.capacity()
                                + originalNorms.capacity())
                                * sizeof(double));
                        for (int pass = 0; pass < 2; ++pass) {
                            weighted = particulars;
                            for (int column = 0;
                                 column < count; ++column) {
                                double* values = weighted.data()
                                    + static_cast<std::size_t>(
                                        column * port.rows);
                                for (int row = 0;
                                     row < port.rows; ++row) {
                                    values[row] *= targetMass[
                                        static_cast<std::size_t>(row)];
                                }
                            }
                            cblas_dgemm(
                                CblasColMajor, CblasTrans,
                                CblasNoTrans, baseRank, count,
                                port.rows, 1.0, port.basis.data(),
                                port.rows, weighted.data(), port.rows,
                                0.0, coefficients.data(), baseRank);
                            cblas_dgemm(
                                CblasColMajor, CblasNoTrans,
                                CblasNoTrans, port.rows, count,
                                baseRank, -1.0, port.basis.data(),
                                port.rows, coefficients.data(),
                                baseRank, 1.0, particulars.data(),
                                port.rows);
                        }
                    }
                    port.candidateColumns += count;
                    for (int column = 0; column < count; ++column) {
                        std::vector<double> candidate(
                            particulars.begin()
                                + static_cast<std::ptrdiff_t>(
                                    column * port.rows),
                            particulars.begin()
                                + static_cast<std::ptrdiff_t>(
                                    (column + 1) * port.rows));
                        mandatoryTransientWorkspacePeak = std::max(
                            mandatoryTransientWorkspacePeak,
                            (block.capacity()
                                + particulars.capacity()
                                + originalNorms.capacity()
                                + candidate.capacity())
                                * sizeof(double));
                        for (int pass = 0; pass < 2; ++pass) {
                            for (int mode = baseRank;
                                 mode < port.rank; ++mode) {
                                const double* basis =
                                    port.basis.data()
                                    + static_cast<std::size_t>(
                                        mode * port.rows);
                                const double coefficient = massDot(
                                    basis, candidate.data(),
                                    targetMass);
                                for (int row = 0;
                                     row < port.rows; ++row) {
                                    candidate[
                                        static_cast<std::size_t>(row)]
                                        -= coefficient * basis[row];
                                }
                            }
                        }
                        const double norm = std::sqrt(std::max(
                            0.0, massDot(
                                candidate.data(), candidate.data(),
                                targetMass)));
                        if (!(originalNorms[
                                    static_cast<std::size_t>(column)]
                                > 0.0)
                            || !(norm
                                > options.relativeDeflationTolerance
                                    * originalNorms[
                                        static_cast<std::size_t>(
                                            column)])) {
                            continue;
                        }
                        for (double& value : candidate) value /= norm;
                        port.basis.insert(
                            port.basis.end(),
                            candidate.begin(), candidate.end());
                        ++port.rank;
                        ++port.acceptedColumns;
                        ++accepted;
                    }
#else
                    for (int column = 0;
                         column < count; ++column) {
                        std::vector<double> particular(
                            particulars.begin()
                                + static_cast<std::ptrdiff_t>(
                                    column * port.rows),
                            particulars.begin()
                                + static_cast<std::ptrdiff_t>(
                                    (column + 1) * port.rows));
                        mandatoryTransientWorkspacePeak = std::max(
                            mandatoryTransientWorkspacePeak,
                            (block.capacity()
                                + particulars.capacity()
                                + particular.capacity())
                                * sizeof(double));
                        if (appendMassMode(
                                port, std::move(particular),
                                targetMass, port.rows,
                                options.relativeDeflationTolerance)) {
                            ++accepted;
                        }
                    }
#endif
                }
            };
        appendParticularFamily(
            condensedInputs, sourceChannels, diagnostics.inputRank);
        appendParticularFamily(
            condensedBoundary, boundaryChannels,
            diagnostics.boundaryRank);
        appendParticularFamily(
            condensedHistory,
            diagnostics.compressedHistoryRank,
            diagnostics.historyRank,
            &historySelection.selectedColumns);
        if (options.basisMethod == "hybrid-randomized") {
            std::cout << "[Hybrid port] interface="
                << patch.interfaceId << ", mandatory rank="
                << port.rank << " (constant="
                << diagnostics.constantRank << ", geometry="
                << diagnostics.geometryRank << ", input="
                << diagnostics.inputRank << ", boundary="
                << diagnostics.boundaryRank << ", history="
                << diagnostics.historyRank << ")\n" << std::flush;
        }
        diagnostics.mandatoryRankTotal = port.rank;
        port.mandatoryModes = port.rank;
        result.mandatoryModeSeconds +=
            std::chrono::duration<double>(
                Clock::now() - mandatoryStart).count();

        if (initialTransferModel != nullptr) {
            const auto transfer = std::find_if(
                initialTransferModel->ports.begin(),
                initialTransferModel->ports.end(),
                [&](const LocalPortBasis& value) {
                    return value.interfaceId == patch.interfaceId;
                });
            if (transfer == initialTransferModel->ports.end()) {
                throw std::runtime_error(
                    "[Hybrid port] Missing randomized basis for "
                    "interface " + std::to_string(patch.interfaceId));
            }
            if (transfer->rows != port.rows
                || transfer->interfaceIndices != port.interfaceIndices
                || transfer->targetFingerprint
                    != port.targetFingerprint) {
                throw std::runtime_error(
                    "[Hybrid port] Randomized target mapping differs "
                    "from the mandatory port.");
            }
            diagnostics.requestedRandomizedRank =
                transfer->requestedTransferRank;
            port.requestedTransferRank =
                transfer->requestedTransferRank;
            const auto randomizedStart = Clock::now();
            for (int mode = 0; mode < transfer->rank; ++mode) {
                std::vector<double> candidate(
                    transfer->basis.begin()
                        + static_cast<std::ptrdiff_t>(
                            mode * transfer->rows),
                    transfer->basis.begin()
                        + static_cast<std::ptrdiff_t>(
                            (mode + 1) * transfer->rows));
                if (appendMassMode(
                        port, std::move(candidate), targetMass,
                        port.rows,
                        options.relativeDeflationTolerance)) {
                    ++diagnostics.acceptedRandomizedRank;
                    ++port.spectralModes;
                    port.spectralValues.push_back(
                        mode < static_cast<int>(
                            transfer->spectralValues.size())
                        ? transfer->spectralValues[
                            static_cast<std::size_t>(mode)]
                        : 0.0);
                    port.spectralResiduals.push_back(
                        mode < static_cast<int>(
                            transfer->spectralResiduals.size())
                        ? transfer->spectralResiduals[
                            static_cast<std::size_t>(mode)]
                        : transfer->transferIndicator);
                }
            }
            result.orthogonalizationSeconds +=
                std::chrono::duration<double>(
                    Clock::now() - randomizedStart).count();
            std::cout << "[Hybrid port] interface="
                << patch.interfaceId << ", randomized accepted="
                << diagnostics.acceptedRandomizedRank << '/'
                << diagnostics.requestedRandomizedRank
                << '\n' << std::flush;
        }

        const auto probeStart = Clock::now();
        std::vector<std::vector<double>> rawProbes;
        // Input, boundary, and history channels already have their exact
        // S_tt^{-1} particular responses in Phi_mandatory. Their Schur
        // residual is therefore structurally zero (including a channel
        // deflated as dependent on earlier particulars), so retaining full
        // target-sized copies here only duplicates memory and work.
        const int structurallyRepresentedProbes =
            sourceChannels + boundaryChannels
            + diagnostics.compressedHistoryRank;
        if (options.probeMode == "operator-geometry"
            && !patch.source.empty()) {
            const std::vector<double> sourceMass = traceWeights(
                dynamicModel, patch.source, traceMassDiagonal,
                penaltyMassDiagonal, options.innerProduct);
            std::map<int, std::vector<int>> sourceGroups;
            for (int row = 0;
                 row < static_cast<int>(patch.source.size()); ++row) {
                const auto owner =
                    physicalOwner.find(patch.source[
                        static_cast<std::size_t>(row)]);
                const int group = owner == physicalOwner.end()
                    ? -1 : owner->second;
                sourceGroups[group].push_back(row);
            }
            std::vector<std::vector<double>> sourceGeometry;
            for (const auto& group : sourceGroups) {
                for (int coordinate = 0; coordinate <= 3; ++coordinate) {
                    std::vector<double> candidate(
                        patch.source.size(), 0.0);
                    for (int row : group.second) {
                        const int gamma = patch.source[
                            static_cast<std::size_t>(row)];
                        const int global =
                            dynamicModel.interfaceGlobalDofs[
                                static_cast<std::size_t>(gamma)];
                        const Vec3& point =
                            mesh.nodes[
                                static_cast<std::size_t>(global)].p;
                        candidate[static_cast<std::size_t>(row)] =
                            coordinate == 0 ? 1.0
                            : (coordinate == 1 ? point.x
                            : (coordinate == 2 ? point.y : point.z));
                    }
                    const double original = std::sqrt(std::max(
                        0.0, massDot(candidate.data(),
                            candidate.data(), sourceMass)));
                    if (!(original > 0.0)) continue;
                    for (int pass = 0; pass < 2; ++pass) {
                        for (const auto& basis : sourceGeometry) {
                            const double coefficient = massDot(
                                basis.data(), candidate.data(),
                                sourceMass);
                            for (std::size_t row = 0;
                                 row < candidate.size(); ++row) {
                                candidate[row] -=
                                    coefficient * basis[row];
                            }
                        }
                    }
                    const double norm = std::sqrt(std::max(
                        0.0, massDot(candidate.data(),
                            candidate.data(), sourceMass)));
                    if (!(norm
                            > options.relativeDeflationTolerance
                                * original)) {
                        continue;
                    }
                    for (double& value : candidate) value /= norm;
                    sourceGeometry.push_back(candidate);
                    std::vector<double> targetRhs;
                    targetSolver.formTraceRightHandSide(
                        candidate, targetRhs);
                    rawProbes.push_back(std::move(targetRhs));
                }
            }
        }
        diagnostics.rawProbeColumns =
            structurallyRepresentedProbes
            + static_cast<int>(rawProbes.size());

        // Deterministic weighted pivoted QR in the G_t^{-1} metric.
        std::vector<std::vector<double>> probes;
        std::vector<bool> used(rawProbes.size(), false);
        for (;;) {
            int pivot = -1;
            double pivotNorm = 0.0;
            std::vector<double> pivotVector;
            for (int candidateIndex = 0;
                 candidateIndex < static_cast<int>(rawProbes.size());
                 ++candidateIndex) {
                if (used[static_cast<std::size_t>(candidateIndex)]) {
                    continue;
                }
                std::vector<double> candidate =
                    rawProbes[static_cast<std::size_t>(candidateIndex)];
                const double original =
                    inverseMassNorm(candidate, targetMass);
                if (!(original > 0.0)) continue;
                for (double& value : candidate) value /= original;
                for (int pass = 0; pass < 2; ++pass) {
                    for (const auto& basis : probes) {
                        const double coefficient = inverseMassDot(
                            basis, candidate, targetMass);
                        for (std::size_t row = 0;
                             row < candidate.size(); ++row) {
                            candidate[row] -=
                                coefficient * basis[row];
                        }
                    }
                }
                const double norm =
                    inverseMassNorm(candidate, targetMass);
                if (norm > pivotNorm) {
                    pivotNorm = norm;
                    pivot = candidateIndex;
                    pivotVector = std::move(candidate);
                }
            }
            if (pivot < 0
                || !(pivotNorm
                    > options.relativeDeflationTolerance)) {
                break;
            }
            used[static_cast<std::size_t>(pivot)] = true;
            for (double& value : pivotVector) value /= pivotNorm;
            probes.push_back(std::move(pivotVector));
        }
        diagnostics.independentProbeColumns =
            static_cast<int>(probes.size());
        diagnostics.deflatedProbeColumns =
            diagnostics.rawProbeColumns
            - diagnostics.independentProbeColumns;
        diagnostics.probeBlockRank =
            diagnostics.independentProbeColumns;
        if (options.basisMethod == "hybrid-randomized") {
            std::cout << "[Hybrid port] interface="
                << patch.interfaceId << ", probes raw/independent="
                << diagnostics.rawProbeColumns << '/'
                << diagnostics.independentProbeColumns
                << '\n' << std::flush;
        }
        result.probeSetupSeconds +=
            std::chrono::duration<double>(
                Clock::now() - probeStart).count();

        std::vector<std::vector<double>> schurImages;
        auto updateSchurImages = [&]() {
            while (schurImages.size()
                   < static_cast<std::size_t>(port.rank)) {
                const int mode =
                    static_cast<int>(schurImages.size());
                std::vector<double> basisColumn(
                    static_cast<std::size_t>(port.rows), 0.0);
                std::copy_n(
                    port.basis.data()
                        + static_cast<std::size_t>(
                            mode * port.rows),
                    port.rows, basisColumn.data());
                std::vector<double> image;
                targetSolver.applyTargetAction(
                    basisColumn, image);
                ++diagnostics.schurApplyCount;
                schurImages.push_back(std::move(image));
            }
        };
        auto residuals = [&]() {
            updateSchurImages();
            std::vector<std::vector<double>> values;
            values.reserve(probes.size());
            if (port.rank == 0) return values;
            std::vector<double> projected(
                static_cast<std::size_t>(
                    port.rank * port.rank), 0.0);
            for (int rowMode = 0; rowMode < port.rank; ++rowMode) {
                const double* rowBasis = port.basis.data()
                    + static_cast<std::size_t>(
                        rowMode * port.rows);
                for (int columnMode = 0;
                     columnMode < port.rank; ++columnMode) {
                    projected[static_cast<std::size_t>(
                        rowMode * port.rank + columnMode)] =
                        dot(rowBasis,
                            schurImages[
                                static_cast<std::size_t>(
                                    columnMode)].data(),
                            port.rows);
                }
            }
            for (int row = 0; row < port.rank; ++row) {
                for (int column = row + 1;
                     column < port.rank; ++column) {
                    const double symmetric = 0.5 * (
                        projected[static_cast<std::size_t>(
                            row * port.rank + column)]
                        + projected[static_cast<std::size_t>(
                            column * port.rank + row)]);
                    projected[static_cast<std::size_t>(
                        row * port.rank + column)] = symmetric;
                    projected[static_cast<std::size_t>(
                        column * port.rank + row)] = symmetric;
                }
            }
#ifdef USE_MKL_PARDISO
            std::vector<double> projectedLu = projected;
            std::vector<lapack_int> projectedPivots(
                static_cast<std::size_t>(port.rank), 0);
            const lapack_int factorInfo = LAPACKE_dgetrf(
                LAPACK_ROW_MAJOR, port.rank, port.rank,
                projectedLu.data(), port.rank,
                projectedPivots.data());
            if (factorInfo != 0) {
                throw std::runtime_error(
                    "[Residual Krylov] Projected target Schur "
                    "LU factorization failed on interface "
                    + std::to_string(port.interfaceId)
                    + " at rank " + std::to_string(port.rank)
                    + " with info=" + std::to_string(factorInfo));
            }
#else
            local::DenseSymmetricFactor factor;
            try {
                factor =
                    local::factorDenseSymmetric(projected, port.rank);
            } catch (const std::exception& error) {
                throw std::runtime_error(
                    "[Residual Krylov] Projected target Schur "
                    "factorization failed on interface "
                    + std::to_string(port.interfaceId)
                    + " at rank " + std::to_string(port.rank)
                    + ": " + error.what());
            }
#endif
            for (const auto& probe : probes) {
                std::vector<double> coefficients(
                    static_cast<std::size_t>(port.rank), 0.0);
                for (int mode = 0; mode < port.rank; ++mode) {
                    coefficients[static_cast<std::size_t>(mode)] =
                        dot(port.basis.data()
                                + static_cast<std::size_t>(
                                    mode * port.rows),
                            probe.data(), port.rows);
                }
#ifdef USE_MKL_PARDISO
                const lapack_int solveInfo = LAPACKE_dgetrs(
                    LAPACK_ROW_MAJOR, 'N', port.rank, 1,
                    projectedLu.data(), port.rank,
                    projectedPivots.data(),
                    coefficients.data(), 1);
                if (solveInfo != 0) {
                    throw std::runtime_error(
                        "[Residual Krylov] Projected target Schur "
                        "solve failed with info="
                        + std::to_string(solveInfo));
                }
#else
                local::solveDenseSymmetric(factor, coefficients);
#endif
                std::vector<double> residual = probe;
                for (int mode = 0; mode < port.rank; ++mode) {
                    const double coefficient =
                        coefficients[
                            static_cast<std::size_t>(mode)];
                    for (int row = 0; row < port.rows; ++row) {
                        residual[static_cast<std::size_t>(row)]
                            -= coefficient
                            * schurImages[
                                static_cast<std::size_t>(mode)]
                                [static_cast<std::size_t>(row)];
                    }
                }
                values.push_back(std::move(residual));
            }
            return values;
        };
        auto indicators = [&](const std::vector<std::vector<double>>& values) {
            std::vector<double> eta(values.size(), 0.0);
            for (std::size_t column = 0;
                 column < values.size(); ++column) {
                eta[column] = inverseMassNorm(
                    values[column], targetMass)
                    / std::max(1.0e-300,
                        inverseMassNorm(
                            probes[column], targetMass));
            }
            return eta;
        };
        auto residualBlock = residuals();
        auto eta = indicators(residualBlock);
        diagnostics.initialMaximumProbeResidual = eta.empty()
            ? 0.0 : *std::max_element(eta.begin(), eta.end());

        const auto enrichmentStart = Clock::now();
        while (diagnostics.acceptedEnrichmentRank
                    < diagnostics.requestedEnrichmentRank
            && diagnostics.enrichmentSweeps
                    < options.maximumSweeps
            && !eta.empty()
            && *std::max_element(eta.begin(), eta.end())
                    > options.residualTolerance) {
            ++diagnostics.enrichmentSweeps;
            std::vector<int> order(eta.size());
            std::iota(order.begin(), order.end(), 0);
            std::stable_sort(
                order.begin(), order.end(),
                [&](int left, int right) {
                    return eta[static_cast<std::size_t>(left)]
                        > eta[static_cast<std::size_t>(right)];
                });
            const int remaining =
                diagnostics.requestedEnrichmentRank
                - diagnostics.acceptedEnrichmentRank;
            const int selectedCount = std::min({
                remaining, options.blockSize,
                static_cast<int>(order.size())});
            int acceptedThisSweep = 0;
            for (int selected = 0;
                 selected < selectedCount; ++selected) {
                const int probeIndex =
                    order[static_cast<std::size_t>(selected)];
                if (eta[static_cast<std::size_t>(probeIndex)]
                        <= options.residualTolerance) {
                    break;
                }
                std::vector<double> correction;
                targetSolver.solveTargetResponse(
                    residualBlock[
                        static_cast<std::size_t>(probeIndex)],
                    correction);
                const auto orthStart = Clock::now();
                if (appendMassMode(
                        port, std::move(correction), targetMass,
                        port.rows,
                        options.relativeDeflationTolerance)) {
                    ++diagnostics.acceptedEnrichmentRank;
                    ++port.spectralModes;
                    port.spectralValues.push_back(0.0);
                    port.spectralResiduals.push_back(
                        eta[static_cast<std::size_t>(
                            probeIndex)]);
                    ++acceptedThisSweep;
                } else {
                    diagnostics.deflationLimited = true;
                }
                result.orthogonalizationSeconds +=
                    std::chrono::duration<double>(
                        Clock::now() - orthStart).count();
            }
            if (acceptedThisSweep == 0) break;
            residualBlock = residuals();
            eta = indicators(residualBlock);
        }
        result.enrichmentSeconds +=
            std::chrono::duration<double>(
                Clock::now() - enrichmentStart).count();
        diagnostics.finalMaximumProbeResidual = eta.empty()
            ? 0.0 : *std::max_element(eta.begin(), eta.end());
        diagnostics.residualReductionFactor =
            diagnostics.initialMaximumProbeResidual > 0.0
            ? diagnostics.finalMaximumProbeResidual
                / diagnostics.initialMaximumProbeResidual
            : 1.0;

        // Weighted target-solve adjoint check.  It uses the exact transpose
        // path and the same G_t metric as the basis construction.
        std::vector<double> adjointInput(
            patch.target.size(), 0.0);
        std::vector<double> adjointTest(
            patch.target.size(), 0.0);
        for (std::size_t row = 0;
             row < patch.target.size(); ++row) {
            adjointInput[row] =
                std::sin(0.271 * static_cast<double>(row + 1));
            adjointTest[row] =
                std::cos(0.419 * static_cast<double>(row + 1));
        }
        std::vector<double> forward;
        targetSolver.solveTargetResponse(
            adjointInput, forward);
        std::vector<double> weightedTest = adjointTest;
        for (std::size_t row = 0;
             row < weightedTest.size(); ++row) {
            weightedTest[row] *= targetMass[row];
        }
        std::vector<double> transpose;
        targetSolver.solveTargetResponseTranspose(
            weightedTest, transpose);
        const double adjointLeft = massDot(
            adjointTest.data(), forward.data(), targetMass);
        const double adjointRight = dot(
            transpose.data(), adjointInput.data(),
            static_cast<int>(transpose.size()));
        diagnostics.weightedAdjointError =
            std::abs(adjointLeft - adjointRight)
            / std::max({
                1.0e-300, std::abs(adjointLeft),
                std::abs(adjointRight)});

        diagnostics.innerSolver = targetSolver.statistics();
        diagnostics.targetSolveCount =
            diagnostics.innerSolver.solveRightHandSides;
        auto nestedVectorBytes =
            [](const std::vector<std::vector<double>>& vectors) {
                std::size_t bytes = 0;
                for (const auto& vector : vectors) {
                    bytes += vector.capacity() * sizeof(double);
                }
                return bytes;
            };
        const std::size_t actualAlgebraWorkspaceBytes =
            port.basis.capacity() * sizeof(double)
            + targetMass.capacity() * sizeof(double)
            + historySelection.selectedColumns.capacity() * sizeof(int)
            + nestedVectorBytes(rawProbes)
            + nestedVectorBytes(probes)
            + nestedVectorBytes(schurImages)
            + nestedVectorBytes(residualBlock)
            + eta.capacity() * sizeof(double);
        diagnostics.peakIncrementalMemoryBytes = std::max(
            diagnostics.historyCompressionWorkspaceBytes,
            diagnostics.innerSolver.peakIncrementalMemoryBytes
                + actualAlgebraWorkspaceBytes
                + mandatoryTransientWorkspacePeak);
        diagnostics.totalSeconds =
            std::chrono::duration<double>(
                Clock::now() - interfaceStart).count();
        diagnostics.status =
            diagnostics.innerSolver.maximumRelativeResidual
                    > 1.0e-9
                ? "target_residual_failed"
                : (diagnostics.weightedAdjointError > 1.0e-8
                    ? "weighted_adjoint_failed"
                    : (diagnostics.acceptedEnrichmentRank
                            < diagnostics.requestedEnrichmentRank
                        && diagnostics.finalMaximumProbeResidual
                            > options.residualTolerance
                        && !diagnostics.deflationLimited
                        ? "enrichment_rank_failed"
                        : (diagnostics.requestedEnrichmentRank > 0
                            && diagnostics.finalMaximumProbeResidual
                                > options.residualTolerance
                            && diagnostics.finalMaximumProbeResidual
                                >= diagnostics.initialMaximumProbeResidual
                            ? "probe_residual_not_reduced"
                            : "success")));

        diagnostics.mandatoryRankTotal = port.mandatoryModes;
        port.orthogonalityError =
            portOrthogonalityError(port, targetMass);
        std::uint64_t basisFingerprint =
            UINT64_C(1469598103934665603);
        hashValue(basisFingerprint, port.interfaceId);
        hashValue(basisFingerprint, port.targetFingerprint);
        hashValue(basisFingerprint, port.sourceFingerprint);
        for (double value : port.basis) {
            hashValue(basisFingerprint, value);
        }
        port.fingerprint = basisFingerprint;
        result.model.reducedInterfaceDofs += port.rank;
        result.model.modelBytes +=
            port.interfaceIndices.capacity() * sizeof(int)
            + port.sourceIndices.capacity() * sizeof(int)
            + port.patchSubdomains.capacity() * sizeof(int)
            + port.basis.capacity() * sizeof(double);
        result.model.ports.push_back(std::move(port));
        result.interfaces.push_back(std::move(diagnostics));
    }
    result.model.snapshotSeconds = 0.0;
    result.model.basisSeconds =
        std::chrono::duration<double>(
            Clock::now() - totalStart).count();
    result.model.modelBytes +=
        result.model.interfaceGlobalDofs.capacity() * sizeof(int);
    result.totalSeconds = result.model.basisSeconds;
    return result;
}

} // namespace mor::transient
