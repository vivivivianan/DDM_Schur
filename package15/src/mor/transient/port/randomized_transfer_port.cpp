#include "randomized_transfer_port.hpp"

#include "port_basis_factory.hpp"
#include "mor/transient/interface_flux_operator.hpp"

#include "ddm_schur/interface_operator.hpp"
#include "sipg_core.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

namespace mor::transient::port {
namespace {

using Clock = std::chrono::steady_clock;

double elapsed(const Clock::time_point& start)
{
    return std::chrono::duration<double>(
        Clock::now() - start).count();
}

template <typename T>
void hashValue(std::uint64_t& hash, const T& value)
{
    const auto* bytes =
        reinterpret_cast<const unsigned char*>(&value);
    for (std::size_t byte = 0; byte < sizeof(T); ++byte) {
        hash ^= bytes[byte];
        hash *= UINT64_C(1099511628211);
    }
}

std::uint64_t fingerprintVector(
    const std::vector<double>& values)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    hashValue(hash, values.size());
    for (double value : values) hashValue(hash, value);
    return hash;
}

std::uint64_t splitmix64(std::uint64_t value)
{
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30))
        * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27))
        * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

double weightedDot(
    const double* left,
    const double* right,
    const std::vector<double>& mass)
{
    long double result = 0.0L;
    for (std::size_t row = 0; row < mass.size(); ++row) {
        result += static_cast<long double>(mass[row])
            * left[row] * right[row];
    }
    return static_cast<double>(result);
}

double weightedNorm(
    const double* vector,
    const std::vector<double>& mass)
{
    return std::sqrt(std::max(
        0.0, weightedDot(vector, vector, mass)));
}

std::vector<double> metricWeights(
    const local::Model& model,
    const std::vector<int>& indices,
    const std::vector<double>& traceMassDiagonal,
    const std::vector<double>& penaltyMassDiagonal,
    const std::string& innerProduct)
{
    if (traceMassDiagonal.size() != penaltyMassDiagonal.size()
        || traceMassDiagonal.empty()) {
        throw std::runtime_error(
            "[Randomized port] SIPG trace metrics are unavailable.");
    }
    const std::vector<double>* diagonal = nullptr;
    if (innerProduct == "trace-mass") {
        diagonal = &traceMassDiagonal;
    } else if (innerProduct == "penalty-weighted-mass") {
        diagonal = &penaltyMassDiagonal;
    } else {
        throw std::runtime_error(
            "[Randomized port] Unsupported inner product.");
    }
    std::vector<double> weights(indices.size(), 0.0);
    double maximum = 0.0;
    for (std::size_t row = 0; row < indices.size(); ++row) {
        const int gamma = indices[row];
        if (gamma < 0
            || gamma
                >= static_cast<int>(
                    model.interfaceGlobalDofs.size())) {
            throw std::runtime_error(
                "[Randomized port] Interface index is out of range.");
        }
        const int global = model.interfaceGlobalDofs[
            static_cast<std::size_t>(gamma)];
        if (global < 0
            || global >= static_cast<int>(diagonal->size())) {
            throw std::runtime_error(
                "[Randomized port] Metric index is out of range.");
        }
        const double value =
            (*diagonal)[static_cast<std::size_t>(global)];
        if (!std::isfinite(value) || value < 0.0) {
            throw std::runtime_error(
                "[Randomized port] Metric is not nonnegative finite.");
        }
        weights[row] = value;
        maximum = std::max(maximum, value);
    }
    if (!(maximum > 0.0)) {
        throw std::runtime_error(
            "[Randomized port] Metric is zero on an interface.");
    }
    const double floor =
        128.0 * std::numeric_limits<double>::epsilon()
        * maximum;
    for (double& value : weights) {
        value = std::max(value, floor);
    }
    return weights;
}

struct WeightedQrResult {
    std::vector<double> basis;
    int rank = 0;
    double orthogonalityError = 0.0;
    double sampledProjectionError = 0.0;
};

void swapColumns(
    std::vector<double>& matrix,
    int rows,
    int left,
    int right);

WeightedQrResult mixedFluxPivotedMgs(
    std::vector<double> candidates,
    int rows,
    int columns,
    int requestedRank,
    const std::vector<double>& targetMass,
    const PortPatch& patch,
    const InterfaceFluxOperator& fluxOperator,
    const std::string& fluxType,
    double relativeTolerance)
{
    WeightedQrResult result;
    if (rows <= 0 || columns <= 0
        || candidates.size()
            != static_cast<std::size_t>(rows * columns)
        || targetMass.size() != static_cast<std::size_t>(rows)) {
        return result;
    }
    const bool physical =
        fluxType == "physical" || fluxType == "both";
    const bool numerical =
        fluxType == "numerical" || fluxType == "both";
    std::vector<InterfaceFluxResponse> images;
    images.reserve(static_cast<std::size_t>(columns));
    for (int column = 0; column < columns; ++column) {
        std::vector<double> trace(
            candidates.begin()
                + static_cast<std::ptrdiff_t>(column * rows),
            candidates.begin()
                + static_cast<std::ptrdiff_t>((column + 1) * rows));
        images.push_back(
            fluxOperator.applyTargetTrace(patch.target, trace));
    }
    const int triangles =
        static_cast<int>(images.front().areas.size());
    const int augmentedRows = rows
        + (physical ? 2 * triangles : 0)
        + (numerical ? triangles : 0);
    std::vector<double> augmented(
        static_cast<std::size_t>(augmentedRows * columns), 0.0);
    for (int column = 0; column < columns; ++column) {
        double* out = augmented.data()
            + static_cast<std::size_t>(column * augmentedRows);
        const double* temperature = candidates.data()
            + static_cast<std::size_t>(column * rows);
        for (int row = 0; row < rows; ++row) {
            out[row] = std::sqrt(targetMass[
                static_cast<std::size_t>(row)]) * temperature[row];
        }
        int offset = rows;
        for (int triangle = 0; triangle < triangles; ++triangle) {
            const std::size_t index =
                static_cast<std::size_t>(triangle);
            const double scale =
                std::sqrt(images[column].areas[index])
                / std::max(1.0e-300,
                    std::abs(images[column].penalties[index]));
            if (physical) {
                out[offset++] = scale
                    * images[column].leftPhysicalOutward[index];
                out[offset++] = scale
                    * images[column].rightPhysicalOutward[index];
            }
            if (numerical) {
                out[offset++] =
                    scale * images[column].numerical[index];
            }
        }
    }
    auto columnNorm = [&](int column) {
        const double* vector = augmented.data()
            + static_cast<std::size_t>(column * augmentedRows);
        long double value = 0.0L;
        for (int row = 0; row < augmentedRows; ++row) {
            value += static_cast<long double>(vector[row])
                * vector[row];
        }
        return std::sqrt(static_cast<double>(value));
    };
    std::vector<double> originalNorms(
        static_cast<std::size_t>(columns), 0.0);
    double largestOriginal = 0.0;
    for (int column = 0; column < columns; ++column) {
        originalNorms[static_cast<std::size_t>(column)] =
            columnNorm(column);
        largestOriginal = std::max(
            largestOriginal,
            originalNorms[static_cast<std::size_t>(column)]);
    }
    const double threshold = relativeTolerance
        * std::max(1.0e-300, largestOriginal);
    const int maximumRank =
        std::min({requestedRank, rows, columns});
    for (int selected = 0; selected < maximumRank; ++selected) {
        int pivot = selected;
        double pivotNorm = -1.0;
        for (int column = selected; column < columns; ++column) {
            const double value = columnNorm(column);
            if (value > pivotNorm) {
                pivotNorm = value;
                pivot = column;
            }
        }
        if (!(pivotNorm > threshold) || !std::isfinite(pivotNorm)) {
            break;
        }
        swapColumns(candidates, rows, selected, pivot);
        swapColumns(augmented, augmentedRows, selected, pivot);
        std::swap(originalNorms[static_cast<std::size_t>(selected)],
                  originalNorms[static_cast<std::size_t>(pivot)]);
        double* basis = candidates.data()
            + static_cast<std::size_t>(selected * rows);
        double* augmentedBasis = augmented.data()
            + static_cast<std::size_t>(selected * augmentedRows);
        for (int pass = 0; pass < 2; ++pass) {
            for (int prior = 0; prior < selected; ++prior) {
                const double* priorAugmented = augmented.data()
                    + static_cast<std::size_t>(
                        prior * augmentedRows);
                long double coefficient = 0.0L;
                for (int row = 0; row < augmentedRows; ++row) {
                    coefficient +=
                        static_cast<long double>(
                            priorAugmented[row])
                        * augmentedBasis[row];
                }
                const double alpha =
                    static_cast<double>(coefficient);
                const double* priorBasis = candidates.data()
                    + static_cast<std::size_t>(prior * rows);
                for (int row = 0; row < rows; ++row) {
                    basis[row] -= alpha * priorBasis[row];
                }
                for (int row = 0; row < augmentedRows; ++row) {
                    augmentedBasis[row]
                        -= alpha * priorAugmented[row];
                }
            }
        }
        pivotNorm = columnNorm(selected);
        if (!(pivotNorm > threshold) || !std::isfinite(pivotNorm)) {
            break;
        }
        for (int row = 0; row < rows; ++row) {
            basis[row] /= pivotNorm;
        }
        for (int row = 0; row < augmentedRows; ++row) {
            augmentedBasis[row] /= pivotNorm;
        }
        ++result.rank;
        for (int column = selected + 1;
             column < columns; ++column) {
            double* candidateAugmented = augmented.data()
                + static_cast<std::size_t>(
                    column * augmentedRows);
            double* candidate = candidates.data()
                + static_cast<std::size_t>(column * rows);
            for (int pass = 0; pass < 2; ++pass) {
                long double coefficient = 0.0L;
                for (int row = 0; row < augmentedRows; ++row) {
                    coefficient +=
                        static_cast<long double>(
                            augmentedBasis[row])
                        * candidateAugmented[row];
                }
                const double alpha =
                    static_cast<double>(coefficient);
                for (int row = 0; row < rows; ++row) {
                    candidate[row] -= alpha * basis[row];
                }
                for (int row = 0; row < augmentedRows; ++row) {
                    candidateAugmented[row]
                        -= alpha * augmentedBasis[row];
                }
            }
        }
    }
    result.basis.assign(
        candidates.begin(),
        candidates.begin()
            + static_cast<std::ptrdiff_t>(result.rank * rows));
    for (int column = result.rank; column < columns; ++column) {
        result.sampledProjectionError = std::max(
            result.sampledProjectionError,
            columnNorm(column) / std::max(
                1.0e-300,
                originalNorms[static_cast<std::size_t>(column)]));
    }
    for (int left = 0; left < result.rank; ++left) {
        const double* leftColumn = augmented.data()
            + static_cast<std::size_t>(left * augmentedRows);
        for (int right = 0; right <= left; ++right) {
            const double* rightColumn = augmented.data()
                + static_cast<std::size_t>(right * augmentedRows);
            long double product = 0.0L;
            for (int row = 0; row < augmentedRows; ++row) {
                product += static_cast<long double>(
                    leftColumn[row]) * rightColumn[row];
            }
            result.orthogonalityError = std::max(
                result.orthogonalityError,
                std::abs(static_cast<double>(product)
                    - (left == right ? 1.0 : 0.0)));
        }
    }
    return result;
}

void swapColumns(
    std::vector<double>& matrix,
    int rows,
    int left,
    int right)
{
    if (left == right) return;
    for (int row = 0; row < rows; ++row) {
        std::swap(
            matrix[static_cast<std::size_t>(left * rows + row)],
            matrix[static_cast<std::size_t>(right * rows + row)]);
    }
}

WeightedQrResult weightedPivotedMgs(
    std::vector<double> candidates,
    int rows,
    int columns,
    int requestedRank,
    const std::vector<double>& mass,
    double relativeTolerance)
{
    WeightedQrResult result;
    if (rows <= 0 || columns <= 0 || requestedRank <= 0
        || candidates.size() != static_cast<std::size_t>(
            rows * columns)
        || mass.size() != static_cast<std::size_t>(rows)) {
        return result;
    }
    std::vector<double> originalNorms(
        static_cast<std::size_t>(columns), 0.0);
    double largestOriginal = 0.0;
    for (int column = 0; column < columns; ++column) {
        originalNorms[static_cast<std::size_t>(column)] =
            weightedNorm(
                candidates.data()
                    + static_cast<std::size_t>(column * rows),
                mass);
        largestOriginal = std::max(
            largestOriginal,
            originalNorms[static_cast<std::size_t>(column)]);
    }
    const double threshold =
        relativeTolerance * std::max(1.0e-300, largestOriginal);
    const int maximumRank =
        std::min({requestedRank, rows, columns});
    for (int selected = 0;
         selected < maximumRank; ++selected) {
        int pivot = selected;
        double pivotNorm = -1.0;
        for (int column = selected; column < columns; ++column) {
            const double norm = weightedNorm(
                candidates.data()
                    + static_cast<std::size_t>(column * rows),
                mass);
            if (norm > pivotNorm) {
                pivotNorm = norm;
                pivot = column;
            }
        }
        if (!(pivotNorm > threshold) || !std::isfinite(pivotNorm)) {
            break;
        }
        swapColumns(candidates, rows, selected, pivot);
        std::swap(
            originalNorms[static_cast<std::size_t>(selected)],
            originalNorms[static_cast<std::size_t>(pivot)]);
        double* basis = candidates.data()
            + static_cast<std::size_t>(selected * rows);
        for (int pass = 0; pass < 2; ++pass) {
            for (int prior = 0; prior < selected; ++prior) {
                const double* previous = candidates.data()
                    + static_cast<std::size_t>(prior * rows);
                const double coefficient =
                    weightedDot(previous, basis, mass);
                for (int row = 0; row < rows; ++row) {
                    basis[row] -= coefficient * previous[row];
                }
            }
        }
        pivotNorm = weightedNorm(basis, mass);
        if (!(pivotNorm > threshold) || !std::isfinite(pivotNorm)) {
            break;
        }
        for (int row = 0; row < rows; ++row) {
            basis[row] /= pivotNorm;
        }
        int signRow = 0;
        for (int row = 1; row < rows; ++row) {
            if (std::abs(basis[row])
                > std::abs(basis[signRow])) {
                signRow = row;
            }
        }
        if (basis[signRow] < 0.0) {
            for (int row = 0; row < rows; ++row) {
                basis[row] = -basis[row];
            }
        }
        ++result.rank;
        for (int column = selected + 1;
             column < columns; ++column) {
            double* candidate = candidates.data()
                + static_cast<std::size_t>(column * rows);
            for (int pass = 0; pass < 2; ++pass) {
                const double coefficient =
                    weightedDot(basis, candidate, mass);
                for (int row = 0; row < rows; ++row) {
                    candidate[row] -= coefficient * basis[row];
                }
            }
        }
    }
    result.basis.assign(
        candidates.begin(),
        candidates.begin()
            + static_cast<std::ptrdiff_t>(result.rank * rows));
    for (int column = result.rank;
         column < columns; ++column) {
        const double residual = weightedNorm(
            candidates.data()
                + static_cast<std::size_t>(column * rows),
            mass);
        result.sampledProjectionError = std::max(
            result.sampledProjectionError,
            residual / std::max(
                1.0e-300,
                originalNorms[static_cast<std::size_t>(column)]));
    }
    for (int left = 0; left < result.rank; ++left) {
        for (int right = 0; right <= left; ++right) {
            const double product = weightedDot(
                result.basis.data()
                    + static_cast<std::size_t>(left * rows),
                result.basis.data()
                    + static_cast<std::size_t>(right * rows),
                mass);
            result.orthogonalityError = std::max(
                result.orthogonalityError,
                std::abs(
                    product - (left == right ? 1.0 : 0.0)));
        }
    }
    return result;
}

std::vector<double> deterministicProbes(
    int rows,
    int columns,
    const std::vector<double>& sourceMass,
    std::uint64_t seed,
    int interfaceId,
    std::uint64_t sourceFingerprint)
{
    std::vector<double> probes(
        static_cast<std::size_t>(rows * columns), 0.0);
    for (int column = 0; column < columns; ++column) {
        for (int row = 0; row < rows; ++row) {
            std::uint64_t key = seed;
            key ^= splitmix64(
                static_cast<std::uint64_t>(
                    static_cast<std::uint32_t>(interfaceId)));
            key ^= splitmix64(sourceFingerprint);
            key ^= splitmix64(
                static_cast<std::uint64_t>(column)
                    * UINT64_C(0x9e3779b97f4a7c15)
                + static_cast<std::uint64_t>(row));
            const double sign =
                (splitmix64(key) & UINT64_C(1)) != 0
                ? 1.0 : -1.0;
            probes[static_cast<std::size_t>(
                column * rows + row)] =
                sign / std::sqrt(sourceMass[
                    static_cast<std::size_t>(row)]);
        }
    }
    return probes;
}

std::uint64_t basisFingerprint(const LocalPortBasis& basis)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    hashValue(hash, basis.interfaceId);
    hashValue(hash, basis.targetFingerprint);
    hashValue(hash, basis.sourceFingerprint);
    hashValue(hash, basis.rank);
    for (double value : basis.basis) hashValue(hash, value);
    return hash;
}

void jacobiEigenvalues(std::vector<double> matrix,
                       int n,
                       std::vector<double>& values)
{
    const int maximumSweeps = std::max(16, 12 * n * n);
    for (int sweep = 0; sweep < maximumSweeps; ++sweep) {
        int p = 0;
        int q = 0;
        double maximum = 0.0;
        for (int row = 0; row < n; ++row) {
            for (int column = row + 1; column < n; ++column) {
                const double magnitude = std::abs(
                    matrix[static_cast<std::size_t>(
                        row * n + column)]);
                if (magnitude > maximum) {
                    maximum = magnitude;
                    p = row;
                    q = column;
                }
            }
        }
        if (maximum <= 128.0
                * std::numeric_limits<double>::epsilon()) {
            break;
        }
        const double app = matrix[static_cast<std::size_t>(
            p * n + p)];
        const double aqq = matrix[static_cast<std::size_t>(
            q * n + q)];
        const double apq = matrix[static_cast<std::size_t>(
            p * n + q)];
        const double angle =
            0.5 * std::atan2(2.0 * apq, aqq - app);
        const double c = std::cos(angle);
        const double s = std::sin(angle);
        for (int k = 0; k < n; ++k) {
            const double mkp = matrix[static_cast<std::size_t>(
                k * n + p)];
            const double mkq = matrix[static_cast<std::size_t>(
                k * n + q)];
            matrix[static_cast<std::size_t>(k * n + p)] =
                c * mkp - s * mkq;
            matrix[static_cast<std::size_t>(k * n + q)] =
                s * mkp + c * mkq;
        }
        for (int k = 0; k < n; ++k) {
            const double mpk = matrix[static_cast<std::size_t>(
                p * n + k)];
            const double mqk = matrix[static_cast<std::size_t>(
                q * n + k)];
            matrix[static_cast<std::size_t>(p * n + k)] =
                c * mpk - s * mqk;
            matrix[static_cast<std::size_t>(q * n + k)] =
                s * mpk + c * mqk;
        }
    }
    values.resize(static_cast<std::size_t>(n), 0.0);
    for (int index = 0; index < n; ++index) {
        values[static_cast<std::size_t>(index)] =
            matrix[static_cast<std::size_t>(index * n + index)];
    }
}

} // namespace

RandomizedTransferBasisBuilder::RandomizedTransferBasisBuilder(
    RandomizedTransferPortOptions options)
    : options_(std::move(options))
{
    if (options_.requestedRank <= 0
        || options_.oversampling < 0
        || (options_.powerIterations != 0
            && options_.powerIterations != 1)
        || !(options_.relativeDeflationTolerance > 0.0)
        || (options_.sourceMode != "trace-only"
            && options_.sourceMode != "trace-plus-input"
            && options_.sourceMode != "generalized-dynamic")
        || options_.innerSolver.innerSolver
            != "woodbury-exact"
        || (options_.fluxType != "numerical"
            && options_.fluxType != "physical"
            && options_.fluxType != "both")) {
        throw std::runtime_error(
            "[Randomized port] Invalid builder options.");
    }
}

std::string RandomizedTransferBasisBuilder::methodName() const
{
    return "randomized-transfer";
}

PortBasisResult RandomizedTransferBasisBuilder::build(
    const PortBasisBuildContext& context,
    const PortPatch& patch)
{
    const auto totalStart = Clock::now();
    PortBasisResult result;
    result.methodName = methodName();
    result.physicalInterfaceId = patch.interfaceId;
    result.targetDofs = static_cast<int>(patch.target.size());
    result.requestedRank = options_.requestedRank;
    result.oversampling = options_.oversampling;
    result.powerIterations = options_.powerIterations;
    result.seed = options_.seed;
    result.sourceMode = options_.sourceMode;
    result.serializationVersion = 6;
    result.snapshotUsed = false;
    result.fomUsedForBasis = false;

    const bool includeTrace = true;
    const bool includeInput =
        options_.sourceMode == "trace-plus-input"
        || options_.sourceMode == "generalized-dynamic";
    const bool includeBoundary =
        options_.sourceMode == "generalized-dynamic";
    const bool includeHistory =
        options_.sourceMode == "generalized-dynamic";
    GeneralizedPatchTransferOperator transfer(
        context.schur, patch, context.sourceBlocks,
        includeTrace, includeInput, includeBoundary,
        includeHistory, options_.innerSolver);
    result.sourceDofs = transfer.sourceRows();
    if (result.sourceDofs <= 0) {
        result.status = "empty_transfer_source";
        result.innerSolver = transfer.statistics();
        result.basisBuildTime = elapsed(totalStart);
        return result;
    }
    const int probeColumns = std::min(
        result.sourceDofs,
        options_.requestedRank + options_.oversampling);
    result.probeColumns = probeColumns;
    if (probeColumns < options_.requestedRank) {
        result.status = "source_rank_below_requested_rank";
        result.innerSolver = transfer.statistics();
        result.basisBuildTime = elapsed(totalStart);
        return result;
    }

    const std::vector<double> targetMass = metricWeights(
        context.dynamicModel, patch.target,
        context.traceMassDiagonal,
        context.penaltyMassDiagonal,
        options_.innerProduct);
    std::vector<double> sourceMass;
    if (includeTrace && !patch.source.empty()) {
        sourceMass = metricWeights(
            context.dynamicModel, patch.source,
            context.traceMassDiagonal,
            context.penaltyMassDiagonal,
            options_.innerProduct);
    }
    if (includeInput) {
        sourceMass.insert(
            sourceMass.end(),
            static_cast<std::size_t>(
                context.sourceBlocks.inputChannels), 1.0);
    }
    if (includeBoundary) {
        sourceMass.insert(
            sourceMass.end(),
            static_cast<std::size_t>(
                context.sourceBlocks.boundaryChannels), 1.0);
    }
    if (includeHistory) {
        sourceMass.insert(
            sourceMass.end(),
            static_cast<std::size_t>(
                context.sourceBlocks.historyChannels), 1.0);
    }
    if (sourceMass.size()
        != static_cast<std::size_t>(result.sourceDofs)) {
        throw std::runtime_error(
            "[Randomized port] Source metric dimensions are invalid.");
    }

    const auto probeStart = Clock::now();
    std::vector<double> probes = deterministicProbes(
        result.sourceDofs, probeColumns, sourceMass,
        options_.seed, patch.interfaceId,
        patch.sourceFingerprint);
    result.probeGenerationSeconds = elapsed(probeStart);
    result.probeMatrixBytes =
        probes.capacity() * sizeof(double);
    std::vector<double> firstProbe(
        probes.begin(),
        probes.begin() + result.sourceDofs);

    auto applyStart = Clock::now();
    std::vector<double> samples;
    transfer.applyMultiple(
        probes, probeColumns, samples);
    result.transferApplySeconds += elapsed(applyStart);
    result.applyCount += probeColumns;
    std::vector<double> firstImage(
        samples.begin(),
        samples.begin() + result.targetDofs);
    result.sampleMatrixBytes =
        samples.capacity() * sizeof(double);
    int activeColumns = probeColumns;

    for (int power = 0;
         power < options_.powerIterations; ++power) {
        // Intermediate weighted orthogonalization is a numerically stable
        // change of coordinates inside range(T Omega); it does not form T
        // and prevents the unscaled T^T T product from losing its smaller
        // requested directions before the final range extraction.
        WeightedQrResult targetRange = weightedPivotedMgs(
            std::move(samples), result.targetDofs,
            activeColumns, activeColumns, targetMass,
            std::min(
                options_.relativeDeflationTolerance, 1.0e-14));
        activeColumns = targetRange.rank;
        if (activeColumns <= 0) {
            result.status = "range_deflated_to_zero";
            result.innerSolver = transfer.statistics();
            result.basisBuildTime = elapsed(totalStart);
            return result;
        }
        const auto transposeStart = Clock::now();
        std::vector<double> sourceImage;
        transfer.applyTransposeMultiple(
            targetRange.basis, activeColumns, sourceImage);
        result.transposeApplySeconds += elapsed(transposeStart);
        result.transposeApplyCount += activeColumns;
        WeightedQrResult sourceRange = weightedPivotedMgs(
            std::move(sourceImage), result.sourceDofs,
            activeColumns, activeColumns, sourceMass,
            std::min(
                options_.relativeDeflationTolerance, 1.0e-14));
        activeColumns = sourceRange.rank;
        if (activeColumns <= 0) {
            result.status = "range_deflated_to_zero";
            result.innerSolver = transfer.statistics();
            result.basisBuildTime = elapsed(totalStart);
            return result;
        }
        applyStart = Clock::now();
        transfer.applyMultiple(
            sourceRange.basis, activeColumns, samples);
        result.transferApplySeconds += elapsed(applyStart);
        result.applyCount += activeColumns;
        result.sampleMatrixBytes = std::max(
            result.sampleMatrixBytes,
            samples.capacity() * sizeof(double));
    }
    probes.clear();
    probes.shrink_to_fit();

    std::vector<double> targetTest(
        static_cast<std::size_t>(result.targetDofs), 0.0);
    for (int row = 0; row < result.targetDofs; ++row) {
        targetTest[static_cast<std::size_t>(row)] =
            std::sin(0.3819660112501051
                * static_cast<double>(row + 1));
    }
    std::vector<double> weightedTarget = targetTest;
    for (int row = 0; row < result.targetDofs; ++row) {
        weightedTarget[static_cast<std::size_t>(row)] *=
            targetMass[static_cast<std::size_t>(row)];
    }
    const auto transposeStart = Clock::now();
    std::vector<double> transposeImage;
    transfer.applyTranspose(
        weightedTarget, transposeImage);
    result.transposeApplySeconds += elapsed(transposeStart);
    ++result.transposeApplyCount;
    const double adjointLeft = weightedDot(
        targetTest.data(), firstImage.data(), targetMass);
    long double adjointRight = 0.0L;
    for (int row = 0; row < result.sourceDofs; ++row) {
        adjointRight += static_cast<long double>(
            transposeImage[static_cast<std::size_t>(row)])
            * firstProbe[static_cast<std::size_t>(row)];
    }
    result.weightedAdjointError =
        std::abs(adjointLeft - static_cast<double>(adjointRight))
        / std::max({
            1.0e-300, std::abs(adjointLeft),
            std::abs(static_cast<double>(adjointRight))});

    const auto qrStart = Clock::now();
    WeightedQrResult qr;
    if (options_.fluxAware) {
        InterfaceFluxOperator fluxOperator(
            context.mesh, context.physics,
            context.dynamicModel, context.schur,
            patch.interfaceId, patch.leftSubdomain,
            patch.rightSubdomain);
        qr = mixedFluxPivotedMgs(
            std::move(samples), result.targetDofs,
            activeColumns, options_.requestedRank,
            targetMass, patch, fluxOperator,
            options_.fluxType,
            options_.relativeDeflationTolerance);
    } else {
        qr = weightedPivotedMgs(
            std::move(samples), result.targetDofs,
            activeColumns, options_.requestedRank,
            targetMass, options_.relativeDeflationTolerance);
    }
    result.qrSeconds = elapsed(qrStart);
    result.qrWorkspaceBytes = static_cast<std::size_t>(
        result.targetDofs * activeColumns)
        * sizeof(double) * (options_.fluxAware ? 4U : 1U);
    result.acceptedRank = qr.rank;
    result.portRank = qr.rank;
    result.basisDimension = qr.rank;
    result.orthogonalityError = qr.orthogonalityError;
    result.basisErrorIndicator =
        qr.sampledProjectionError;

    LocalPortBasis basis;
    basis.interfaceId = patch.interfaceId;
    basis.leftSubdomain = patch.leftSubdomain;
    basis.rightSubdomain = patch.rightSubdomain;
    basis.rows = result.targetDofs;
    basis.rank = qr.rank;
    basis.candidateColumns = probeColumns;
    basis.acceptedColumns = qr.rank;
    basis.templateId = patch.interfaceId;
    basis.targetFingerprint = patch.targetFingerprint;
    basis.sourceFingerprint = patch.sourceFingerprint;
    basis.traceSourceFingerprint =
        patch.sourceFingerprint;
    basis.inputSourceFingerprint =
        includeInput
        ? context.sourceBlocks.inputFingerprint : 0;
    basis.boundarySourceFingerprint =
        includeBoundary
        ? context.sourceBlocks.boundaryFingerprint : 0;
    basis.historySourceFingerprint =
        includeHistory
        ? context.sourceBlocks.historyFingerprint : 0;
    basis.traceSourceRows =
        includeTrace
        ? static_cast<int>(patch.source.size()) : 0;
    basis.inputSourceRows =
        includeInput
        ? context.sourceBlocks.inputChannels : 0;
    basis.boundarySourceRows =
        includeBoundary
        ? context.sourceBlocks.boundaryChannels : 0;
    basis.historySourceRows =
        includeHistory
        ? context.sourceBlocks.historyChannels : 0;
    basis.requestedTransferRank =
        options_.requestedRank;
    basis.spectralModes = qr.rank;
    basis.spectralValues.assign(
        static_cast<std::size_t>(qr.rank), 0.0);
    basis.spectralResiduals.assign(
        static_cast<std::size_t>(qr.rank),
        qr.sampledProjectionError);
    basis.retainedEnergy = std::max(
        0.0, 1.0 - qr.sampledProjectionError
            * qr.sampledProjectionError);
    basis.orthogonalityError =
        qr.orthogonalityError;
    basis.transferIndicator =
        qr.sampledProjectionError;
    basis.patchSubdomains = patch.patchSubdomains;
    basis.sourceIndices = patch.source;
    basis.interfaceIndices = patch.target;
    basis.basis = std::move(qr.basis);
    basis.fingerprint = basisFingerprint(basis);
    result.finalBasisBytes =
        basis.basis.capacity() * sizeof(double);
    result.basis = std::move(basis);

    result.innerSolver = transfer.statistics();
    result.solverCalls =
        result.innerSolver.solveRightHandSides;
    result.targetSolveCount =
        result.innerSolver.solveRightHandSides;
    result.targetSolvePhase33Calls =
        result.innerSolver.solveCalls;
    result.residual =
        result.innerSolver.maximumRelativeResidual;
    const std::size_t algebraBytes = std::max({
        result.probeMatrixBytes + result.sampleMatrixBytes,
        2 * result.sampleMatrixBytes,
        result.qrWorkspaceBytes + result.finalBasisBytes});
    result.memoryPeak =
        result.innerSolver.peakIncrementalMemoryBytes
        + algebraBytes;
    result.basisBuildTime = elapsed(totalStart);
    result.status =
        result.residual > 1.0e-9
            ? "target_residual_failed"
            : (result.weightedAdjointError > 1.0e-8
                ? "weighted_adjoint_failed"
                : (result.orthogonalityError > 1.0e-10
                    ? "orthogonality_failed"
                    : (result.acceptedRank
                            < options_.requestedRank
                        ? "success_deflated"
                        : "success")));
    return result;
}

RandomizedTransferBuildResult buildRandomizedTransferPortModel(
    const Mesh& mesh,
    const CaseConfig& physics,
    const ddm_schur::InterfacePartition& partition,
    const local::Model& dynamicModel,
    const std::vector<double>& traceMassDiagonal,
    const std::vector<double>& penaltyMassDiagonal,
    const std::vector<double>& input,
    int sourceChannels,
    const std::vector<double>& boundaryLoad,
    const std::vector<double>& condensedHistory,
    int historyChannels,
    const RandomizedTransferPortOptions& options,
    const ReducedDynamicSchurOperator* sharedSchur)
{
    const auto totalStart = Clock::now();
    RandomizedTransferBuildResult result;
    result.model.formatVersion = 6;
    result.model.basisMethod = "randomized-transfer";
    result.model.methodDescription =
        "Local Block Arnoldi with Randomized Transfer Port Space "
        "and Reduced Dynamic Schur";
    result.model.ablationMode = "randomized-range";
    result.model.sourceMode = options.sourceMode;
    result.model.innerProduct = options.innerProduct;
    result.model.rankMode = "fixed";
    result.model.oversamplingLayers =
        options.oversamplingLayers;
    result.model.requestedRank = options.requestedRank;
    result.model.minimumRank = options.requestedRank;
    result.model.maximumRank = options.requestedRank;
    result.model.eigenvalueTolerance = 0.0;
    result.model.eigensolverTolerance = 0.0;
    result.model.eigensolverMaximumIterations =
        options.powerIterations;
    result.model.relativeDeflationTolerance =
        options.relativeDeflationTolerance;
    result.model.innerSolver =
        options.innerSolver.innerSolver;
    result.model.innerSolverTolerance =
        options.innerSolver.relativeTolerance;
    result.model.innerSolverMaximumIterations =
        options.innerSolver.maximumIterations;
    result.model.fullInterfaceDofs =
        dynamicModel.interfaceDofs;
    result.model.interfaceGlobalDofs =
        dynamicModel.interfaceGlobalDofs;

    const auto patchStart = Clock::now();
    std::vector<PortPatch> patches =
        buildOptimalPortPatches(
            mesh, partition, options.oversamplingLayers);
    if (!options.selectedInterfaceIds.empty()) {
        const std::set<int> selected(
            options.selectedInterfaceIds.begin(),
            options.selectedInterfaceIds.end());
        patches.erase(
            std::remove_if(
                patches.begin(), patches.end(),
                [&](const PortPatch& patch) {
                    return selected.count(
                        patch.interfaceId) == 0;
                }),
            patches.end());
        if (patches.size() != selected.size()) {
            throw std::runtime_error(
                "[Randomized port] Selected interface is absent.");
        }
    }
    result.patchSetupSeconds = elapsed(patchStart);

    const auto sourceStart = Clock::now();
    GeneralizedTransferSourceBlocks sourceBlocks;
    if (options.sourceMode == "trace-only") {
        // The trace-only range finder never consumes generalized source
        // columns. Avoid copying the potentially very large Local Block
        // Arnoldi history block merely to compute provenance.
        sourceBlocks.interfaceDofs = dynamicModel.interfaceDofs;
        sourceBlocks.inputFingerprint = fingerprintVector(input);
        sourceBlocks.boundaryFingerprint =
            fingerprintVector(boundaryLoad);
        sourceBlocks.historyFingerprint =
            fingerprintVector(condensedHistory);
    } else {
        sourceBlocks = buildGeneralizedTransferSourceBlocks(
            dynamicModel, input, sourceChannels,
            boundaryLoad, condensedHistory, historyChannels);
    }
    result.sourceSetupSeconds = elapsed(sourceStart);
    result.model.generalizedInputFingerprint =
        sourceBlocks.inputFingerprint;
    result.model.generalizedBoundaryFingerprint =
        sourceBlocks.boundaryFingerprint;
    result.model.generalizedHistoryFingerprint =
        sourceBlocks.historyFingerprint;
    std::uint64_t traceFingerprint =
        UINT64_C(1469598103934665603);
    for (const PortPatch& patch : patches) {
        hashValue(traceFingerprint, patch.interfaceId);
        hashValue(
            traceFingerprint, patch.sourceFingerprint);
    }
    result.model.traceSourceFingerprint =
        traceFingerprint;

    std::unique_ptr<ReducedDynamicSchurOperator> ownedSchur;
    if (sharedSchur == nullptr) {
        ownedSchur =
            std::make_unique<ReducedDynamicSchurOperator>(
                dynamicModel);
        sharedSchur = ownedSchur.get();
    }
    PortBasisBuildContext context{
        mesh, physics, partition, dynamicModel,
        traceMassDiagonal, penaltyMassDiagonal,
        sourceBlocks, *sharedSchur};
    std::unique_ptr<IPortBasisBuilder> builder =
        makePortBasisBuilder(options);
    result.model.ports.reserve(patches.size());
    result.interfaces.reserve(patches.size());
    for (const PortPatch& patch : patches) {
        PortBasisResult interfaceResult =
            builder->build(context, patch);
        result.probeGenerationSeconds +=
            interfaceResult.probeGenerationSeconds;
        result.transferApplySeconds +=
            interfaceResult.transferApplySeconds;
        result.transposeApplySeconds +=
            interfaceResult.transposeApplySeconds;
        result.qrSeconds += interfaceResult.qrSeconds;
        if (interfaceResult.status.rfind("success", 0) != 0) {
            result.interfaces.push_back(
                std::move(interfaceResult));
            continue;
        }
        result.model.reducedInterfaceDofs +=
            interfaceResult.basis.rank;
        result.model.modelBytes +=
            interfaceResult.basis.interfaceIndices.capacity()
                * sizeof(int)
            + interfaceResult.basis.sourceIndices.capacity()
                * sizeof(int)
            + interfaceResult.basis.patchSubdomains.capacity()
                * sizeof(int)
            + interfaceResult.basis.basis.capacity()
                * sizeof(double)
            + interfaceResult.basis.spectralValues.capacity()
                * sizeof(double)
            + interfaceResult.basis.spectralResiduals.capacity()
                * sizeof(double);
        result.model.ports.push_back(
            interfaceResult.basis);
        result.interfaces.push_back(
            std::move(interfaceResult));
    }
    result.model.modelBytes +=
        result.model.interfaceGlobalDofs.capacity()
        * sizeof(int);
    result.totalSeconds = elapsed(totalStart);
    result.model.basisSeconds = result.totalSeconds;
    result.model.snapshotSeconds = 0.0;
    return result;
}

WeightedPortSubspaceComparison compareWeightedPortSubspaces(
    const LocalPortBasis& candidate,
    const LocalPortBasis& reference,
    const local::Model& dynamicModel,
    const std::vector<double>& traceMassDiagonal,
    const std::vector<double>& penaltyMassDiagonal,
    const std::string& innerProduct)
{
    WeightedPortSubspaceComparison result;
    if (candidate.rows != reference.rows
        || candidate.rank <= 0 || reference.rank <= 0
        || candidate.interfaceIndices
            != reference.interfaceIndices) {
        result.maximumPrincipalAngleRadians =
            std::numeric_limits<double>::infinity();
        result.projectorDifference =
            std::numeric_limits<double>::infinity();
        return result;
    }
    const std::vector<double> mass = metricWeights(
        dynamicModel, candidate.interfaceIndices,
        traceMassDiagonal, penaltyMassDiagonal,
        innerProduct);
    std::vector<double> cross(static_cast<std::size_t>(
        candidate.rank * reference.rank), 0.0);
    long double overlap = 0.0L;
    for (int left = 0; left < candidate.rank; ++left) {
        for (int right = 0; right < reference.rank; ++right) {
            const double value = weightedDot(
                candidate.basis.data()
                    + static_cast<std::size_t>(
                        left * candidate.rows),
                reference.basis.data()
                    + static_cast<std::size_t>(
                        right * reference.rows),
                mass);
            cross[static_cast<std::size_t>(
                left * reference.rank + right)] = value;
            overlap += static_cast<long double>(value) * value;
        }
    }
    const int commonRank =
        std::min(candidate.rank, reference.rank);
    std::vector<double> gram(static_cast<std::size_t>(
        commonRank * commonRank), 0.0);
    if (reference.rank <= candidate.rank) {
        for (int row = 0; row < reference.rank; ++row) {
            for (int column = 0;
                 column < reference.rank; ++column) {
                long double value = 0.0L;
                for (int mode = 0;
                     mode < candidate.rank; ++mode) {
                    value += static_cast<long double>(
                        cross[static_cast<std::size_t>(
                            mode * reference.rank + row)])
                        * cross[static_cast<std::size_t>(
                            mode * reference.rank + column)];
                }
                gram[static_cast<std::size_t>(
                    row * commonRank + column)] =
                    static_cast<double>(value);
            }
        }
    } else {
        for (int row = 0; row < candidate.rank; ++row) {
            for (int column = 0;
                 column < candidate.rank; ++column) {
                long double value = 0.0L;
                for (int mode = 0;
                     mode < reference.rank; ++mode) {
                    value += static_cast<long double>(
                        cross[static_cast<std::size_t>(
                            row * reference.rank + mode)])
                        * cross[static_cast<std::size_t>(
                            column * reference.rank + mode)];
                }
                gram[static_cast<std::size_t>(
                    row * commonRank + column)] =
                    static_cast<double>(value);
            }
        }
    }
    std::vector<double> eigenvalues;
    jacobiEigenvalues(
        std::move(gram), commonRank, eigenvalues);
    const double minimum = std::clamp(
        *std::min_element(
            eigenvalues.begin(), eigenvalues.end()),
        0.0, 1.0);
    result.maximumPrincipalAngleRadians =
        std::acos(std::sqrt(minimum));
    const long double squared = std::max(
        0.0L,
        static_cast<long double>(
            candidate.rank + reference.rank)
            - 2.0L * overlap);
    result.projectorDifference =
        std::sqrt(static_cast<double>(squared))
        / std::sqrt(static_cast<double>(
            std::max(1, reference.rank)));
    return result;
}

} // namespace mor::transient::port
