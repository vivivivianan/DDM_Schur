#include "global_interface_coarse.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#endif

namespace mor::transient {
namespace {

using Clock = std::chrono::steady_clock;

double secondsSince(const Clock::time_point& start)
{
    return std::chrono::duration<double>(
        Clock::now() - start).count();
}

double dot(const double* left, const double* right, int size)
{
    long double value = 0.0L;
    for (int row = 0; row < size; ++row) {
        value += static_cast<long double>(left[row]) * right[row];
    }
    return static_cast<double>(value);
}

double norm(const std::vector<double>& values)
{
    return std::sqrt(std::max(
        0.0, dot(values.data(), values.data(),
            static_cast<int>(values.size()))));
}

double relativeDifference(const std::vector<double>& left,
                          const std::vector<double>& right)
{
    if (left.size() != right.size()) {
        throw std::runtime_error(
            "[Global coarse] Relative-difference size mismatch.");
    }
    long double difference = 0.0L;
    long double magnitude = 0.0L;
    for (std::size_t row = 0; row < left.size(); ++row) {
        const double delta = left[row] - right[row];
        difference += static_cast<long double>(delta) * delta;
        magnitude += static_cast<long double>(right[row]) * right[row];
    }
    return std::sqrt(static_cast<double>(difference))
        / std::max(
            std::numeric_limits<double>::epsilon(),
            std::sqrt(static_cast<double>(magnitude)));
}

std::size_t workingSetBytes()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX counters{};
    counters.cb = sizeof(counters);
    if (GetProcessMemoryInfo(
            GetCurrentProcess(),
            reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&counters),
            sizeof(counters))) {
        return static_cast<std::size_t>(counters.WorkingSetSize);
    }
#endif
    return 0;
}

template <typename T>
void hashValue(std::uint64_t& hash, const T& value)
{
    const auto* bytes =
        reinterpret_cast<const unsigned char*>(&value);
    for (std::size_t index = 0; index < sizeof(T); ++index) {
        hash ^= bytes[index];
        hash *= UINT64_C(1099511628211);
    }
}

std::uint64_t localBasisFingerprint(
    const LocalPortModel& model)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    hashValue(hash, model.fullInterfaceDofs);
    hashValue(hash, model.reducedInterfaceDofs);
    for (const LocalPortBasis& port : model.ports) {
        hashValue(hash, port.interfaceId);
        hashValue(hash, port.rank);
        for (int index : port.interfaceIndices) {
            hashValue(hash, index);
        }
        for (double value : port.basis) {
            hashValue(hash, value);
        }
    }
    return hash;
}

struct CompactSeed {
    std::vector<double> values;
    int family = 0;
    int sourceColumn = -1;
};

double inverseMetricDot(
    const std::vector<double>& left,
    const std::vector<double>& right,
    const std::vector<int>& target,
    const std::vector<double>& metric)
{
    long double value = 0.0L;
    for (std::size_t row = 0; row < target.size(); ++row) {
        value += static_cast<long double>(left[row]) * right[row]
            / metric[static_cast<std::size_t>(target[row])];
    }
    return static_cast<double>(value);
}

double inverseMetricNorm(
    const std::vector<double>& value,
    const std::vector<int>& target,
    const std::vector<double>& metric)
{
    return std::sqrt(std::max(
        0.0, inverseMetricDot(value, value, target, metric)));
}

std::vector<int> topHistoryColumns(
    const GeneralizedTransferSourceBlocks& sources,
    const std::vector<int>& target,
    const std::vector<double>& metric,
    int requested)
{
    std::vector<std::pair<double, int>> norms;
    norms.reserve(static_cast<std::size_t>(sources.historyChannels));
    for (int column = 0;
         column < sources.historyChannels; ++column) {
        const double* values = sources.history.data()
            + static_cast<std::size_t>(
                column * sources.interfaceDofs);
        long double squared = 0.0L;
        for (int row : target) {
            squared += static_cast<long double>(values[row])
                * values[row]
                / metric[static_cast<std::size_t>(row)];
        }
        if (squared > 0.0L) {
            norms.emplace_back(
                static_cast<double>(squared), column);
        }
    }
    std::sort(norms.begin(), norms.end(),
        [](const auto& left, const auto& right) {
            return left.first > right.first
                || (left.first == right.first
                    && left.second < right.second);
        });
    const int retained = std::min(
        requested, static_cast<int>(norms.size()));
    std::vector<int> selected;
    selected.reserve(static_cast<std::size_t>(retained));
    for (int index = 0; index < retained; ++index) {
        selected.push_back(
            norms[static_cast<std::size_t>(index)].second);
    }
    return selected;
}

void appendSourceColumns(
    std::vector<CompactSeed>& seeds,
    const std::vector<double>& source,
    int rows,
    const std::vector<int>& columns,
    const std::vector<int>& target,
    int family)
{
    for (int column : columns) {
        CompactSeed seed;
        seed.family = family;
        seed.sourceColumn = column;
        seed.values.resize(target.size(), 0.0);
        const double* values = source.data()
            + static_cast<std::size_t>(column * rows);
        for (std::size_t row = 0; row < target.size(); ++row) {
            seed.values[row] =
                values[target[row]];
        }
        seeds.push_back(std::move(seed));
    }
}

std::vector<CompactSeed> pivotedSeedQr(
    const std::vector<CompactSeed>& raw,
    const std::vector<int>& target,
    const std::vector<double>& metric,
    int maximumColumns,
    double tolerance)
{
    std::vector<CompactSeed> accepted;
    std::vector<unsigned char> used(raw.size(), 0);
    accepted.reserve(static_cast<std::size_t>(maximumColumns));
    while (static_cast<int>(accepted.size()) < maximumColumns) {
        int pivot = -1;
        double pivotNorm = 0.0;
        CompactSeed pivotSeed;
        for (int candidate = 0;
             candidate < static_cast<int>(raw.size()); ++candidate) {
            if (used[static_cast<std::size_t>(candidate)] != 0) {
                continue;
            }
            CompactSeed value =
                raw[static_cast<std::size_t>(candidate)];
            const double original = inverseMetricNorm(
                value.values, target, metric);
            if (!(original > 0.0)) continue;
            for (double& entry : value.values) {
                entry /= original;
            }
            for (int pass = 0; pass < 2; ++pass) {
                for (const CompactSeed& basis : accepted) {
                    const double coefficient =
                        inverseMetricDot(
                            basis.values, value.values,
                            target, metric);
                    for (std::size_t row = 0;
                         row < value.values.size(); ++row) {
                        value.values[row] -= coefficient
                            * basis.values[row];
                    }
                }
            }
            const double candidateNorm = inverseMetricNorm(
                value.values, target, metric);
            if (candidateNorm > pivotNorm) {
                pivotNorm = candidateNorm;
                pivot = candidate;
                pivotSeed = std::move(value);
            }
        }
        if (pivot < 0 || !(pivotNorm > tolerance)) break;
        used[static_cast<std::size_t>(pivot)] = 1;
        for (double& entry : pivotSeed.values) {
            entry /= pivotNorm;
        }
        accepted.push_back(std::move(pivotSeed));
    }
    return accepted;
}

struct CandidateSpace {
    int rows = 0;
    std::vector<double> vectors;
    std::vector<double> images;
};

bool appendEnergyCandidate(
    CandidateSpace& space,
    std::vector<double> candidate,
    LocalPortReducedSchurSolver& localSolver,
    const ReducedDynamicSchurOperator& schur,
    double tolerance,
    int& schurApplyCount)
{
    for (int projectionPass = 0;
         projectionPass < 2; ++projectionPass) {
        std::vector<double> projected;
        localSolver.projectSchurEnergyComplement(
            candidate, projected);
        ++schurApplyCount;
        candidate = std::move(projected);
    }
    std::vector<double> image;
    schur.apply(candidate, image);
    ++schurApplyCount;
    const double originalEnergy = dot(
        candidate.data(), image.data(), space.rows);
    if (!(originalEnergy > 0.0)
        || !std::isfinite(originalEnergy)) {
        return false;
    }
    for (int pass = 0; pass < 2; ++pass) {
        const int rank = static_cast<int>(
            space.vectors.size()
            / static_cast<std::size_t>(space.rows));
        for (int mode = 0; mode < rank; ++mode) {
            const double* basis = space.vectors.data()
                + static_cast<std::size_t>(mode * space.rows);
            const double* basisImage = space.images.data()
                + static_cast<std::size_t>(mode * space.rows);
            const double coefficient =
                dot(basis, image.data(), space.rows);
            for (int row = 0; row < space.rows; ++row) {
                candidate[static_cast<std::size_t>(row)] -=
                    coefficient * basis[row];
                image[static_cast<std::size_t>(row)] -=
                    coefficient * basisImage[row];
            }
        }
    }
    const double energy = dot(
        candidate.data(), image.data(), space.rows);
    if (!(energy > tolerance * tolerance * originalEnergy)
        || !std::isfinite(energy)) {
        return false;
    }
    const double inverse = 1.0 / std::sqrt(energy);
    for (double& value : candidate) value *= inverse;
    for (double& value : image) value *= inverse;
    space.vectors.insert(
        space.vectors.end(), candidate.begin(), candidate.end());
    space.images.insert(
        space.images.end(), image.begin(), image.end());
    return true;
}

struct ExactInnerSolveResult {
    std::vector<double> solution;
    int iterations = 0;
    double relativeResidual =
        std::numeric_limits<double>::infinity();
    bool converged = false;
};

ExactInnerSolveResult solveExactSchurPcg(
    const std::vector<double>& rightHandSide,
    LocalPortReducedSchurSolver& localSolver,
    const ReducedDynamicSchurOperator& schur,
    int maximumIterations,
    double tolerance,
    int& schurApplyCount,
    int& preconditionerApplyCount)
{
    const int rows = schur.size();
    ExactInnerSolveResult result;
    result.solution.assign(
        static_cast<std::size_t>(rows), 0.0);
    const double rightHandSideNorm = norm(rightHandSide);
    if (!(rightHandSideNorm > 0.0)) {
        result.relativeResidual = 0.0;
        result.converged = true;
        return result;
    }
    const std::vector<double>& diagonal = schur.diagonal();
    double diagonalScale = 0.0;
    for (double value : diagonal) {
        diagonalScale = std::max(
            diagonalScale, std::abs(value));
    }
    const double diagonalFloor =
        1024.0 * std::numeric_limits<double>::epsilon()
        * std::max(
            std::numeric_limits<double>::min(),
            diagonalScale);

    const auto applyTwoLevelPreconditioner =
        [&](const std::vector<double>& residual,
            std::vector<double>& output) {
            // Symmetric two-level deflated Jacobi:
            // Z E^-1 Z^T r + P D^-1 P^T r,
            // P = I - Z E^-1 Z^T S.
            std::vector<double> coarse;
            localSolver.localGalerkinInterfaceResponse(
                residual, coarse);
            std::vector<double> coarseImage;
            schur.apply(coarse, coarseImage);
            ++schurApplyCount;
            std::vector<double> dualComplement(
                residual.size(), 0.0);
            for (std::size_t row = 0;
                 row < residual.size(); ++row) {
                dualComplement[row] =
                    residual[row] - coarseImage[row];
            }
            std::vector<double> jacobi(
                residual.size(), 0.0);
            for (std::size_t row = 0;
                 row < residual.size(); ++row) {
                if (!(std::abs(diagonal[row])
                          > diagonalFloor)
                    || !std::isfinite(diagonal[row])) {
                    throw std::runtime_error(
                        "[Global coarse] Exact-PCG Jacobi "
                        "preconditioner has a tiny/nonfinite diagonal.");
                }
                jacobi[row] =
                    dualComplement[row] / diagonal[row];
            }
            std::vector<double> primalComplement;
            localSolver.projectSchurEnergyComplement(
                jacobi, primalComplement);
            ++schurApplyCount;
            output.resize(residual.size());
            for (std::size_t row = 0;
                 row < residual.size(); ++row) {
                output[row] =
                    coarse[row] + primalComplement[row];
            }
            ++preconditionerApplyCount;
        };

    std::vector<double> residual = rightHandSide;
    std::vector<double> preconditioned;
    applyTwoLevelPreconditioner(
        residual, preconditioned);
    std::vector<double> direction = preconditioned;
    double residualPreconditioned = dot(
        residual.data(), preconditioned.data(), rows);
    if (!(residualPreconditioned > 0.0)
        || !std::isfinite(residualPreconditioned)) {
        return result;
    }
    for (int iteration = 0;
         iteration < maximumIterations; ++iteration) {
        std::vector<double> image;
        schur.apply(direction, image);
        ++schurApplyCount;
        const double denominator = dot(
            direction.data(), image.data(), rows);
        if (!(denominator > 0.0)
            || !std::isfinite(denominator)) {
            break;
        }
        const double alpha =
            residualPreconditioned / denominator;
        for (int row = 0; row < rows; ++row) {
            result.solution[static_cast<std::size_t>(row)]
                += alpha
                * direction[static_cast<std::size_t>(row)];
            residual[static_cast<std::size_t>(row)]
                -= alpha * image[static_cast<std::size_t>(row)];
        }
        result.iterations = iteration + 1;
        result.relativeResidual =
            norm(residual) / rightHandSideNorm;
        if (result.relativeResidual <= tolerance) {
            result.converged = true;
            break;
        }
        std::vector<double> nextPreconditioned;
        applyTwoLevelPreconditioner(
            residual, nextPreconditioned);
        const double nextResidualPreconditioned = dot(
            residual.data(),
            nextPreconditioned.data(), rows);
        if (!(nextResidualPreconditioned > 0.0)
            || !std::isfinite(
                nextResidualPreconditioned)) {
            break;
        }
        const double beta =
            nextResidualPreconditioned
            / residualPreconditioned;
        for (int row = 0; row < rows; ++row) {
            direction[static_cast<std::size_t>(row)] =
                nextPreconditioned[
                    static_cast<std::size_t>(row)]
                + beta
                * direction[static_cast<std::size_t>(row)];
        }
        preconditioned =
            std::move(nextPreconditioned);
        residualPreconditioned =
            nextResidualPreconditioned;
    }
    // Always recompute the true, unpreconditioned Schur residual.
    std::vector<double> finalImage;
    schur.apply(result.solution, finalImage);
    ++schurApplyCount;
    for (std::size_t row = 0;
         row < finalImage.size(); ++row) {
        finalImage[row] -= rightHandSide[row];
    }
    result.relativeResidual =
        norm(finalImage) / rightHandSideNorm;
    result.converged =
        result.relativeResidual <= tolerance;
    return result;
}

struct RitzData {
    int rank = 0;
    std::vector<double> values;
    std::vector<double> residuals;
    std::vector<double> projectedResiduals;
    std::vector<double> inverseMetricResiduals;
    std::vector<double> vectors;
    std::vector<double> images;
    std::vector<std::vector<double>> residualVectors;
};

RitzData computeRitz(
    const CandidateSpace& space,
    const std::vector<double>& metric,
    int requested)
{
    const int candidateRank = static_cast<int>(
        space.vectors.size()
        / static_cast<std::size_t>(space.rows));
    if (candidateRank <= 0) {
        return {};
    }
    std::vector<double> stiffness(
        static_cast<std::size_t>(
            candidateRank * candidateRank), 0.0);
    std::vector<double> mass(stiffness.size(), 0.0);
    for (int column = 0; column < candidateRank; ++column) {
        const double* right = space.vectors.data()
            + static_cast<std::size_t>(column * space.rows);
        const double* rightImage = space.images.data()
            + static_cast<std::size_t>(column * space.rows);
        for (int rowMode = 0;
             rowMode < candidateRank; ++rowMode) {
            const double* left = space.vectors.data()
                + static_cast<std::size_t>(
                    rowMode * space.rows);
            stiffness[static_cast<std::size_t>(
                rowMode + column * candidateRank)] =
                dot(left, rightImage, space.rows);
            long double massValue = 0.0L;
            for (int row = 0; row < space.rows; ++row) {
                massValue += static_cast<long double>(left[row])
                    * metric[static_cast<std::size_t>(row)]
                    * right[row];
            }
            mass[static_cast<std::size_t>(
                rowMode + column * candidateRank)] =
                static_cast<double>(massValue);
        }
    }
#ifdef USE_MKL_PARDISO
    std::vector<double> eigenvalues(
        static_cast<std::size_t>(candidateRank), 0.0);
    const lapack_int info = LAPACKE_dsygvd(
        LAPACK_COL_MAJOR, 1, 'V', 'U',
        candidateRank, stiffness.data(), candidateRank,
        mass.data(), candidateRank, eigenvalues.data());
    if (info != 0) {
        throw std::runtime_error(
            "[Global coarse] Projected generalized eigensolve failed with info="
            + std::to_string(info));
    }
#else
    throw std::runtime_error(
        "[Global coarse] Prototype requires MKL generalized eigensolver.");
#endif
    RitzData result;
    result.rank = std::min(requested, candidateRank);
    result.values.resize(static_cast<std::size_t>(result.rank));
    result.residuals.resize(static_cast<std::size_t>(result.rank));
    result.projectedResiduals.resize(
        static_cast<std::size_t>(result.rank), 0.0);
    result.inverseMetricResiduals.resize(
        static_cast<std::size_t>(result.rank), 0.0);
    result.vectors.assign(static_cast<std::size_t>(
        space.rows * result.rank), 0.0);
    result.images.assign(result.vectors.size(), 0.0);
    result.residualVectors.resize(
        static_cast<std::size_t>(result.rank));
    for (int selected = 0; selected < result.rank; ++selected) {
        const double eigenvalue =
            eigenvalues[static_cast<std::size_t>(selected)];
        result.values[static_cast<std::size_t>(selected)] =
            eigenvalue;
        double* vector = result.vectors.data()
            + static_cast<std::size_t>(selected * space.rows);
        double* image = result.images.data()
            + static_cast<std::size_t>(selected * space.rows);
        for (int source = 0; source < candidateRank; ++source) {
            const double coefficient =
                stiffness[static_cast<std::size_t>(
                    source + selected * candidateRank)];
            const double* sourceVector = space.vectors.data()
                + static_cast<std::size_t>(source * space.rows);
            const double* sourceImage = space.images.data()
                + static_cast<std::size_t>(source * space.rows);
            for (int row = 0; row < space.rows; ++row) {
                vector[row] += coefficient * sourceVector[row];
                image[row] += coefficient * sourceImage[row];
            }
        }
        const double energy =
            dot(vector, image, space.rows);
        if (!(energy > 0.0) || !std::isfinite(energy)) {
            throw std::runtime_error(
                "[Global coarse] Ritz vector has nonpositive Schur energy.");
        }
        const double inverseEnergyNorm =
            1.0 / std::sqrt(energy);
        for (int row = 0; row < space.rows; ++row) {
            vector[row] *= inverseEnergyNorm;
            image[row] *= inverseEnergyNorm;
        }
        std::vector<double>& residual =
            result.residualVectors[static_cast<std::size_t>(selected)];
        residual.resize(static_cast<std::size_t>(space.rows));
        long double residualSquared = 0.0L;
        long double imageSquared = 0.0L;
        for (int row = 0; row < space.rows; ++row) {
            residual[static_cast<std::size_t>(row)] =
                image[row] - eigenvalue
                    * metric[static_cast<std::size_t>(row)]
                    * vector[row];
            residualSquared += static_cast<long double>(
                residual[static_cast<std::size_t>(row)])
                * residual[static_cast<std::size_t>(row)];
            imageSquared +=
                static_cast<long double>(image[row]) * image[row];
        }
        result.residuals[static_cast<std::size_t>(selected)] =
            std::sqrt(static_cast<double>(residualSquared))
            / std::max(
                std::numeric_limits<double>::epsilon(),
                std::sqrt(static_cast<double>(imageSquared)));
    }
    return result;
}

struct OrthogonalityData {
    double absolute = 0.0;
    double relative = 0.0;
};

OrthogonalityData localCoarseOverlap(
    LocalPortReducedSchurSolver& localSolver,
    const RitzData& ritz)
{
    OrthogonalityData result;
    if (ritz.rank <= 0) return result;
    const int rows = static_cast<int>(
        ritz.images.size()
        / static_cast<std::size_t>(ritz.rank));
    for (int mode = 0; mode < ritz.rank; ++mode) {
        std::vector<double> image(
            ritz.images.begin()
                + static_cast<std::ptrdiff_t>(
                    mode * rows),
            ritz.images.begin()
                + static_cast<std::ptrdiff_t>(
                    (mode + 1) * rows));
        std::vector<double> restricted;
        localSolver.restrictLocal(image, restricted);
        for (double value : restricted) {
            result.absolute = std::max(
                result.absolute, std::abs(value));
        }
        std::vector<double> localProjection;
        localSolver.localGalerkinInterfaceResponse(
            image, localProjection);
        const double projectionEnergy =
            dot(localProjection.data(), image.data(), rows);
        const double* vector = ritz.vectors.data()
            + static_cast<std::size_t>(mode * rows);
        const double modeEnergy =
            dot(vector, image.data(), rows);
        result.relative = std::max(
            result.relative,
            std::sqrt(
                std::max(0.0, std::abs(projectionEnergy))
                / std::max(
                    std::numeric_limits<double>::epsilon(),
                    std::abs(modeEnergy))));
    }
    return result;
}

void finalizeRitzModes(
    RitzData& ritz,
    const std::vector<double>& metric,
    LocalPortReducedSchurSolver& localSolver,
    const ReducedDynamicSchurOperator& schur,
    GlobalInterfaceCoarseDiagnostics& diagnostics)
{
    if (ritz.rank <= 0) return;
    const int rows = schur.size();
    std::vector<double> finalizedVectors;
    std::vector<double> finalizedImages;
    finalizedVectors.reserve(ritz.vectors.size());
    finalizedImages.reserve(ritz.images.size());
    diagnostics.maximumRitzResidual = 0.0;
    diagnostics.maximumProjectedRitzResidual = 0.0;
    diagnostics.maximumInverseMetricRitzResidual = 0.0;
    for (int mode = 0; mode < ritz.rank; ++mode) {
        std::vector<double> vector(
            ritz.vectors.begin()
                + static_cast<std::ptrdiff_t>(mode * rows),
            ritz.vectors.begin()
                + static_cast<std::ptrdiff_t>((mode + 1) * rows));
        for (int pass = 0; pass < 2; ++pass) {
            std::vector<double> projected;
            localSolver.projectSchurEnergyComplement(
                vector, projected);
            ++diagnostics.schurApplyCount;
            vector = std::move(projected);
        }
        std::vector<double> image;
        schur.apply(vector, image);
        ++diagnostics.schurApplyCount;
        for (int pass = 0; pass < 2; ++pass) {
            for (int previous = 0;
                 previous < mode; ++previous) {
                const double* previousVector =
                    finalizedVectors.data()
                    + static_cast<std::size_t>(
                        previous * rows);
                const double* previousImage =
                    finalizedImages.data()
                    + static_cast<std::size_t>(
                        previous * rows);
                const double coefficient = dot(
                    previousVector, image.data(), rows);
                for (int row = 0; row < rows; ++row) {
                    vector[static_cast<std::size_t>(row)]
                        -= coefficient * previousVector[row];
                    image[static_cast<std::size_t>(row)]
                        -= coefficient * previousImage[row];
                }
            }
        }
        const double energy =
            dot(vector.data(), image.data(), rows);
        if (!(energy > 0.0) || !std::isfinite(energy)) {
            throw std::runtime_error(
                "[Global coarse] Final Schur-MGS produced "
                "nonpositive energy.");
        }
        const double inverseEnergy =
            1.0 / std::sqrt(energy);
        for (int row = 0; row < rows; ++row) {
            vector[static_cast<std::size_t>(row)]
                *= inverseEnergy;
            image[static_cast<std::size_t>(row)]
                *= inverseEnergy;
        }
        long double massEnergy = 0.0L;
        for (int row = 0; row < rows; ++row) {
            massEnergy += static_cast<long double>(
                vector[static_cast<std::size_t>(row)])
                * metric[static_cast<std::size_t>(row)]
                * vector[static_cast<std::size_t>(row)];
        }
        if (!(massEnergy > 0.0L)) {
            throw std::runtime_error(
                "[Global coarse] Final mode has "
                "nonpositive trace-metric energy.");
        }
        const double eigenvalue =
            dot(vector.data(), image.data(), rows)
            / static_cast<double>(massEnergy);
        ritz.values[static_cast<std::size_t>(mode)] =
            eigenvalue;
        std::vector<double> residual(
            static_cast<std::size_t>(rows), 0.0);
        long double residualSquared = 0.0L;
        long double imageSquared = 0.0L;
        long double massImageSquared = 0.0L;
        long double inverseMetricResidualSquared = 0.0L;
        long double inverseMetricImageSquared = 0.0L;
        long double inverseMetricMassImageSquared = 0.0L;
        for (int row = 0; row < rows; ++row) {
            const double massImage =
                metric[static_cast<std::size_t>(row)]
                * vector[static_cast<std::size_t>(row)];
            const double value =
                image[static_cast<std::size_t>(row)]
                - eigenvalue * massImage;
            residual[static_cast<std::size_t>(row)] = value;
            residualSquared +=
                static_cast<long double>(value) * value;
            imageSquared += static_cast<long double>(
                image[static_cast<std::size_t>(row)])
                * image[static_cast<std::size_t>(row)];
            massImageSquared +=
                static_cast<long double>(massImage) * massImage;
            inverseMetricResidualSquared +=
                static_cast<long double>(value) * value
                / metric[static_cast<std::size_t>(row)];
            inverseMetricImageSquared +=
                static_cast<long double>(
                    image[static_cast<std::size_t>(row)])
                * image[static_cast<std::size_t>(row)]
                / metric[static_cast<std::size_t>(row)];
            inverseMetricMassImageSquared +=
                static_cast<long double>(massImage) * massImage
                / metric[static_cast<std::size_t>(row)];
        }
        const double residualDenominator = std::max({
            std::numeric_limits<double>::epsilon(),
            std::sqrt(static_cast<double>(imageSquared)),
            std::abs(eigenvalue)
                * std::sqrt(
                    static_cast<double>(massImageSquared))});
        const double relativeResidual =
            std::sqrt(static_cast<double>(residualSquared))
            / residualDenominator;
        const double inverseMetricDenominator = std::max({
            std::numeric_limits<double>::epsilon(),
            std::sqrt(static_cast<double>(
                inverseMetricImageSquared)),
            std::abs(eigenvalue)
                * std::sqrt(static_cast<double>(
                    inverseMetricMassImageSquared))});
        const double inverseMetricResidual =
            std::sqrt(static_cast<double>(
                inverseMetricResidualSquared))
            / inverseMetricDenominator;
        std::vector<double> localResidualResponse;
        localSolver.localGalerkinInterfaceResponse(
            residual, localResidualResponse);
        std::vector<double> localResidualImage;
        schur.apply(
            localResidualResponse, localResidualImage);
        ++diagnostics.schurApplyCount;
        std::vector<double> projectedResidual = residual;
        for (int row = 0; row < rows; ++row) {
            projectedResidual[static_cast<std::size_t>(row)]
                -= localResidualImage[static_cast<std::size_t>(row)];
        }
        const double projectedRelativeResidual =
            norm(projectedResidual) / residualDenominator;
        ritz.residuals[static_cast<std::size_t>(mode)] =
            relativeResidual;
        ritz.projectedResiduals[
            static_cast<std::size_t>(mode)] =
            projectedRelativeResidual;
        ritz.inverseMetricResiduals[
            static_cast<std::size_t>(mode)] =
            inverseMetricResidual;
        ritz.residualVectors[static_cast<std::size_t>(mode)] =
            std::move(residual);
        diagnostics.maximumRitzResidual = std::max(
            diagnostics.maximumRitzResidual,
            relativeResidual);
        diagnostics.maximumProjectedRitzResidual =
            std::max(
                diagnostics.maximumProjectedRitzResidual,
                projectedRelativeResidual);
        diagnostics.maximumInverseMetricRitzResidual =
            std::max(
                diagnostics.maximumInverseMetricRitzResidual,
                inverseMetricResidual);
        finalizedVectors.insert(
            finalizedVectors.end(),
            vector.begin(), vector.end());
        finalizedImages.insert(
            finalizedImages.end(),
            image.begin(), image.end());
    }
    ritz.vectors = std::move(finalizedVectors);
    ritz.images = std::move(finalizedImages);

    diagnostics.coarseSchurGramError = 0.0;
    diagnostics.coarseMassGramError = 0.0;
    std::vector<double> massDiagonal(
        static_cast<std::size_t>(ritz.rank), 0.0);
    for (int left = 0; left < ritz.rank; ++left) {
        const double* leftVector = ritz.vectors.data()
            + static_cast<std::size_t>(left * rows);
        for (int right = 0; right < ritz.rank; ++right) {
            const double* rightVector = ritz.vectors.data()
                + static_cast<std::size_t>(right * rows);
            const double* rightImage = ritz.images.data()
                + static_cast<std::size_t>(right * rows);
            const double schurGram =
                dot(leftVector, rightImage, rows);
            diagnostics.coarseSchurGramError = std::max(
                diagnostics.coarseSchurGramError,
                std::abs(
                    schurGram - (left == right ? 1.0 : 0.0)));
            long double massGram = 0.0L;
            for (int row = 0; row < rows; ++row) {
                massGram +=
                    static_cast<long double>(leftVector[row])
                    * metric[static_cast<std::size_t>(row)]
                    * rightVector[row];
            }
            if (left == right) {
                massDiagonal[static_cast<std::size_t>(left)] =
                    static_cast<double>(massGram);
            }
        }
    }
    for (int left = 0; left < ritz.rank; ++left) {
        const double* leftVector = ritz.vectors.data()
            + static_cast<std::size_t>(left * rows);
        for (int right = left + 1;
             right < ritz.rank; ++right) {
            const double* rightVector = ritz.vectors.data()
                + static_cast<std::size_t>(right * rows);
            long double massGram = 0.0L;
            for (int row = 0; row < rows; ++row) {
                massGram +=
                    static_cast<long double>(leftVector[row])
                    * metric[static_cast<std::size_t>(row)]
                    * rightVector[row];
            }
            const double correlation =
                std::abs(static_cast<double>(massGram))
                / std::sqrt(std::max(
                    std::numeric_limits<double>::epsilon(),
                    massDiagonal[static_cast<std::size_t>(left)]
                    * massDiagonal[static_cast<std::size_t>(right)]));
            diagnostics.coarseMassGramError = std::max(
                diagnostics.coarseMassGramError,
                correlation);
        }
    }
}

void runExplicitProjectedReference(
    const LocalPortModel& localModel,
    LocalPortReducedSchurSolver& localSolver,
    const ReducedDynamicSchurOperator& schur,
    const std::vector<double>& metric,
    const GlobalInterfaceCoarseModel& model,
    GlobalInterfaceCoarseDiagnostics& diagnostics)
{
    diagnostics.explicitReferenceRequested = true;
    const int rows = schur.size();
    const int localRank = localModel.reducedInterfaceDofs;
    if (rows > 1024) {
        diagnostics.explicitReferenceStatus =
            "skipped_dimension_above_dense_reference_limit";
        return;
    }
    if (localRank <= 0 || localRank >= rows
        || model.rank <= 0) {
        diagnostics.explicitReferenceStatus =
            "invalid_reference_dimensions";
        return;
    }
#ifndef USE_MKL_PARDISO
    diagnostics.explicitReferenceStatus =
        "mkl_lapack_required";
    return;
#else
    std::vector<double> denseSchur(
        static_cast<std::size_t>(rows * rows), 0.0);
    for (int column = 0; column < rows; ++column) {
        std::vector<double> unit(
            static_cast<std::size_t>(rows), 0.0);
        unit[static_cast<std::size_t>(column)] = 1.0;
        std::vector<double> image;
        schur.apply(unit, image);
        std::copy(
            image.begin(), image.end(),
            denseSchur.begin()
                + static_cast<std::ptrdiff_t>(
                    column * rows));
    }
    std::vector<double> localBasis(
        static_cast<std::size_t>(rows * localRank), 0.0);
    int offset = 0;
    for (const LocalPortBasis& port : localModel.ports) {
        for (int mode = 0; mode < port.rank; ++mode) {
            for (int localRow = 0;
                 localRow < port.rows; ++localRow) {
                localBasis[static_cast<std::size_t>(
                    port.interfaceIndices[
                        static_cast<std::size_t>(localRow)]
                    + (offset + mode) * rows)] =
                    port.basis[static_cast<std::size_t>(
                        localRow + mode * port.rows)];
            }
        }
        offset += port.rank;
    }
    std::vector<double> schurLocalBasis(
        localBasis.size(), 0.0);
    cblas_dgemm(
        CblasColMajor, CblasNoTrans, CblasNoTrans,
        rows, localRank, rows, 1.0,
        denseSchur.data(), rows,
        localBasis.data(), rows, 0.0,
        schurLocalBasis.data(), rows);

    // C^T phi=0 is exactly V_L^T S phi=0. The complete right singular
    // vectors of C^T provide an explicit Euclidean nullspace reference.
    std::vector<double> constraintTranspose(
        static_cast<std::size_t>(localRank * rows), 0.0);
    for (int column = 0; column < rows; ++column) {
        for (int row = 0; row < localRank; ++row) {
            constraintTranspose[static_cast<std::size_t>(
                row + column * localRank)] =
                schurLocalBasis[static_cast<std::size_t>(
                    column + row * rows)];
        }
    }
    std::vector<double> singularValues(
        static_cast<std::size_t>(localRank), 0.0);
    std::vector<double> rightSingularVectors(
        static_cast<std::size_t>(rows * rows), 0.0);
    std::vector<double> superb(
        static_cast<std::size_t>(
            std::max(1, localRank - 1)), 0.0);
    double unusedLeft = 0.0;
    const lapack_int svdInfo = LAPACKE_dgesvd(
        LAPACK_COL_MAJOR, 'N', 'A',
        localRank, rows,
        constraintTranspose.data(), localRank,
        singularValues.data(), &unusedLeft, 1,
        rightSingularVectors.data(), rows,
        superb.data());
    if (svdInfo != 0) {
        diagnostics.explicitReferenceStatus =
            "constraint_nullspace_svd_failed";
        return;
    }
    const double singularScale = singularValues.empty()
        ? 0.0 : singularValues.front();
    int constraintRank = 0;
    for (double value : singularValues) {
        if (value > 1.0e-12
                * std::max(
                    std::numeric_limits<double>::min(),
                    singularScale)) {
            ++constraintRank;
        }
    }
    const int complementRank = rows - constraintRank;
    if (complementRank < model.rank) {
        diagnostics.explicitReferenceStatus =
            "explicit_complement_too_small";
        return;
    }
    diagnostics.explicitReferenceDimension = complementRank;
    std::vector<double> complement(
        static_cast<std::size_t>(
            rows * complementRank), 0.0);
    for (int mode = 0; mode < complementRank; ++mode) {
        for (int row = 0; row < rows; ++row) {
            complement[static_cast<std::size_t>(
                row + mode * rows)] =
                rightSingularVectors[static_cast<std::size_t>(
                    constraintRank + mode + row * rows)];
        }
    }
    std::vector<double> schurComplement(
        complement.size(), 0.0);
    cblas_dgemm(
        CblasColMajor, CblasNoTrans, CblasNoTrans,
        rows, complementRank, rows, 1.0,
        denseSchur.data(), rows,
        complement.data(), rows, 0.0,
        schurComplement.data(), rows);
    std::vector<double> massComplement = complement;
    for (int mode = 0; mode < complementRank; ++mode) {
        for (int row = 0; row < rows; ++row) {
            massComplement[static_cast<std::size_t>(
                row + mode * rows)] *=
                metric[static_cast<std::size_t>(row)];
        }
    }
    std::vector<double> stiffness(
        static_cast<std::size_t>(
            complementRank * complementRank), 0.0);
    std::vector<double> mass(stiffness.size(), 0.0);
    cblas_dgemm(
        CblasColMajor, CblasTrans, CblasNoTrans,
        complementRank, complementRank, rows, 1.0,
        complement.data(), rows,
        schurComplement.data(), rows, 0.0,
        stiffness.data(), complementRank);
    cblas_dgemm(
        CblasColMajor, CblasTrans, CblasNoTrans,
        complementRank, complementRank, rows, 1.0,
        complement.data(), rows,
        massComplement.data(), rows, 0.0,
        mass.data(), complementRank);
    std::vector<double> eigenvalues(
        static_cast<std::size_t>(complementRank), 0.0);
    const lapack_int eigenInfo = LAPACKE_dsygvd(
        LAPACK_COL_MAJOR, 1, 'V', 'U',
        complementRank, stiffness.data(), complementRank,
        mass.data(), complementRank, eigenvalues.data());
    if (eigenInfo != 0) {
        diagnostics.explicitReferenceStatus =
            "explicit_projected_eigensolve_failed";
        return;
    }
    std::vector<double> referenceModes(
        static_cast<std::size_t>(rows * model.rank), 0.0);
    cblas_dgemm(
        CblasColMajor, CblasNoTrans, CblasNoTrans,
        rows, model.rank, complementRank, 1.0,
        complement.data(), rows,
        stiffness.data(), complementRank, 0.0,
        referenceModes.data(), rows);
    std::vector<double> referenceImages(
        referenceModes.size(), 0.0);
    cblas_dgemm(
        CblasColMajor, CblasNoTrans, CblasNoTrans,
        rows, model.rank, rows, 1.0,
        denseSchur.data(), rows,
        referenceModes.data(), rows, 0.0,
        referenceImages.data(), rows);

    diagnostics.explicitReferenceFullResidual = 0.0;
    diagnostics.explicitReferenceProjectedResidual = 0.0;
    diagnostics.explicitEigenvalueRelativeError = 0.0;
    for (int mode = 0; mode < model.rank; ++mode) {
        double* vector = referenceModes.data()
            + static_cast<std::size_t>(mode * rows);
        double* image = referenceImages.data()
            + static_cast<std::size_t>(mode * rows);
        const double energy = dot(vector, image, rows);
        const double scale = 1.0 / std::sqrt(energy);
        for (int row = 0; row < rows; ++row) {
            vector[row] *= scale;
            image[row] *= scale;
        }
        long double massEnergy = 0.0L;
        for (int row = 0; row < rows; ++row) {
            massEnergy += static_cast<long double>(vector[row])
                * metric[static_cast<std::size_t>(row)]
                * vector[row];
        }
        const double lambda =
            dot(vector, image, rows)
            / static_cast<double>(massEnergy);
        diagnostics.explicitEigenvalueRelativeError = std::max(
            diagnostics.explicitEigenvalueRelativeError,
            std::abs(
                model.eigenvalues[static_cast<std::size_t>(mode)]
                - lambda)
            / std::max(
                std::numeric_limits<double>::epsilon(),
                std::abs(lambda)));
        std::vector<double> residual(
            static_cast<std::size_t>(rows), 0.0);
        long double residualSquared = 0.0L;
        long double imageSquared = 0.0L;
        long double massImageSquared = 0.0L;
        for (int row = 0; row < rows; ++row) {
            const double massImage =
                metric[static_cast<std::size_t>(row)]
                * vector[row];
            residual[static_cast<std::size_t>(row)] =
                image[row] - lambda * massImage;
            residualSquared += static_cast<long double>(
                residual[static_cast<std::size_t>(row)])
                * residual[static_cast<std::size_t>(row)];
            imageSquared +=
                static_cast<long double>(image[row]) * image[row];
            massImageSquared +=
                static_cast<long double>(massImage) * massImage;
        }
        const double denominator = std::max({
            std::numeric_limits<double>::epsilon(),
            std::sqrt(static_cast<double>(imageSquared)),
            std::abs(lambda)
                * std::sqrt(
                    static_cast<double>(massImageSquared))});
        diagnostics.explicitReferenceFullResidual = std::max(
            diagnostics.explicitReferenceFullResidual,
            std::sqrt(static_cast<double>(residualSquared))
                / denominator);
        std::vector<double> localResponse;
        localSolver.localGalerkinInterfaceResponse(
            residual, localResponse);
        std::vector<double> localImage;
        schur.apply(localResponse, localImage);
        for (int row = 0; row < rows; ++row) {
            residual[static_cast<std::size_t>(row)]
                -= localImage[static_cast<std::size_t>(row)];
        }
        diagnostics.explicitReferenceProjectedResidual =
            std::max(
                diagnostics.explicitReferenceProjectedResidual,
                norm(residual) / denominator);
    }
    std::vector<double> cross(
        static_cast<std::size_t>(model.rank * model.rank),
        0.0);
    cblas_dgemm(
        CblasColMajor, CblasTrans, CblasNoTrans,
        model.rank, model.rank, rows, 1.0,
        model.basis.data(), rows,
        referenceImages.data(), rows, 0.0,
        cross.data(), model.rank);
    std::vector<double> principalSingularValues(
        static_cast<std::size_t>(model.rank), 0.0);
    std::vector<double> principalSuperb(
        static_cast<std::size_t>(
            std::max(1, model.rank - 1)), 0.0);
    double unusedPrincipalLeft = 0.0;
    double unusedPrincipalRight = 0.0;
    const lapack_int principalInfo = LAPACKE_dgesvd(
        LAPACK_COL_MAJOR, 'N', 'N',
        model.rank, model.rank, cross.data(), model.rank,
        principalSingularValues.data(),
        &unusedPrincipalLeft, 1,
        &unusedPrincipalRight, 1,
        principalSuperb.data());
    if (principalInfo != 0) {
        diagnostics.explicitReferenceStatus =
            "principal_angle_svd_failed";
        return;
    }
    const double minimumSingular =
        principalSingularValues.empty()
        ? 0.0 : principalSingularValues.back();
    diagnostics.explicitMaximumPrincipalAngleRadians =
        std::acos(std::clamp(
            minimumSingular, 0.0, 1.0));
    diagnostics.explicitReferenceRan = true;
    diagnostics.explicitReferenceStatus = "completed";
#endif
}

double probeResidual(
    const std::vector<std::vector<double>>& seeds,
    LocalPortReducedSchurSolver& localSolver,
    const GlobalInterfaceCoarseModel* coarse)
{
    double maximum = 0.0;
    for (const std::vector<double>& seed : seeds) {
        std::vector<double> response;
        localSolver.localGalerkinInterfaceResponse(seed, response);
        std::vector<double> image;
        localSolver.applyFullInterface(response, image);
        if (coarse) {
            for (int mode = 0; mode < coarse->rank; ++mode) {
                const double* basis = coarse->basis.data()
                    + static_cast<std::size_t>(
                        mode * coarse->fullInterfaceDofs);
                const double* basisImage =
                    coarse->schurImages.data()
                    + static_cast<std::size_t>(
                        mode * coarse->fullInterfaceDofs);
                const double coefficient = dot(
                    basis, seed.data(), coarse->fullInterfaceDofs)
                    / dot(
                        basis, basisImage,
                        coarse->fullInterfaceDofs);
                for (int row = 0;
                     row < coarse->fullInterfaceDofs; ++row) {
                    image[static_cast<std::size_t>(row)] +=
                        coefficient * basisImage[row];
                }
            }
        }
        maximum = std::max(
            maximum, relativeDifference(image, seed));
    }
    return maximum;
}

} // namespace

std::size_t globalCoarseCurrentWorkingSetBytes()
{
    return workingSetBytes();
}

GlobalInterfaceCoarseModel buildGlobalInterfaceCoarsePrototype(
    const LocalPortModel& localModel,
    LocalPortReducedSchurSolver& localSolver,
    const ReducedDynamicSchurOperator& exactSchur,
    const std::vector<PortPatch>& physicalPatches,
    const std::vector<double>& interfaceMetricDiagonal,
    const GeneralizedTransferSourceBlocks& sources,
    const GlobalInterfaceCoarseOptions& options)
{
    const auto start = Clock::now();
    GlobalInterfaceCoarseModel model;
    model.fullInterfaceDofs = exactSchur.size();
    model.interfaceGlobalDofs = localModel.interfaceGlobalDofs;
    model.selectedInterfaceIds = options.selectedInterfaceIds;
    model.localBasisFingerprint = localBasisFingerprint(localModel);
    auto& diagnostics = model.diagnostics;
    diagnostics.requestedRank = options.requestedRank;
    diagnostics.requestedSpectralOperator = "S^{-1}M";
    diagnostics.eigenvalueMapping =
        "lambda=(phi^T S phi)/(phi^T M phi)";
    diagnostics.innerSolverRequested =
        options.inverseMode == "exact-pcg"
        ? "matrix-free-pcg" : "none";
    diagnostics.baselineWorkingSetBytes = workingSetBytes();
    diagnostics.peakWorkingSetBytes =
        diagnostics.baselineWorkingSetBytes;

    if (options.requestedRank <= 0
        || options.requestedRank > 32
        || options.maximumIterations <= 0
        || options.innerMaximumIterations <= 0
        || options.historyRank <= 0
        || options.krylovSweeps <= 0
        || !(options.ritzTolerance > 0.0)
        || !(options.innerTolerance > 0.0)
        || !(options.orthogonalityTolerance > 0.0)
        || !(options.deflationTolerance > 0.0)
        || options.selectedInterfaceIds.empty()
        || localModel.fullInterfaceDofs != exactSchur.size()
        || interfaceMetricDiagonal.size()
            != static_cast<std::size_t>(exactSchur.size())
        || sources.interfaceDofs != exactSchur.size()) {
        throw std::runtime_error(
            "[Global coarse] Invalid prototype options or dimensions.");
    }
    if (options.inverseMode == "schur-jacobi") {
        diagnostics.actualSpectralOperator =
            "D_Jacobi^{-1}M";
        diagnostics.jacobiRole =
            "direct_approximate_inverse";
        diagnostics.innerSolverActual = "not_run";
    } else if (options.inverseMode == "exact-pcg") {
        diagnostics.actualSpectralOperator =
            "S^{-1}M";
        diagnostics.jacobiRole =
            "two_level_pcg_preconditioner_only";
        diagnostics.innerSolverActual =
            "matrix-free-pcg";
    } else {
        throw std::runtime_error(
            "[Global coarse] Unsupported inverse mode: "
            + options.inverseMode);
    }
    std::vector<double> metric = interfaceMetricDiagonal;
    diagnostics.metricRawMinimum =
        std::numeric_limits<double>::infinity();
    diagnostics.metricRawMaximum = 0.0;
    for (double value : metric) {
        if (value < 0.0 || !std::isfinite(value)) {
            throw std::runtime_error(
                "[Global coarse] Interface metric must be "
                "nonnegative finite.");
        }
        diagnostics.metricRawMinimum = std::min(
            diagnostics.metricRawMinimum, value);
        diagnostics.metricRawMaximum = std::max(
            diagnostics.metricRawMaximum, value);
    }
    if (!(diagnostics.metricRawMaximum > 0.0)) {
        throw std::runtime_error(
            "[Global coarse] Interface metric is identically zero.");
    }
    // This is the same machine-scale floor used by Optimal Port. SIPG has
    // normal-derivative interface DOFs whose value trace is exactly zero, so
    // the assembled trace Gram is positive semidefinite. The floor only makes
    // the generalized-eigenproblem metric definite; it does not regularize S.
    diagnostics.metricFloor =
        128.0 * std::numeric_limits<double>::epsilon()
        * diagnostics.metricRawMaximum;
    for (double& value : metric) {
        if (value < diagnostics.metricFloor) {
            value = diagnostics.metricFloor;
            ++diagnostics.metricFlooredRows;
        }
    }

    const std::set<int> selectedIds(
        options.selectedInterfaceIds.begin(),
        options.selectedInterfaceIds.end());
    if (selectedIds.size()
        != options.selectedInterfaceIds.size()) {
        throw std::runtime_error(
            "[Global coarse] Selected interface ids must be unique.");
    }
    std::set<int> targetSet;
    std::vector<const PortPatch*> selectedPatches;
    for (const PortPatch& patch : physicalPatches) {
        if (selectedIds.count(patch.interfaceId) == 0) continue;
        selectedPatches.push_back(&patch);
        targetSet.insert(patch.target.begin(), patch.target.end());
    }
    if (selectedPatches.size() != selectedIds.size()
        || targetSet.empty()) {
        throw std::runtime_error(
            "[Global coarse] Selected interfaces are absent from the owner map.");
    }
    std::vector<int> target(targetSet.begin(), targetSet.end());

    std::vector<CompactSeed> raw;
    std::vector<int> inputColumns(
        static_cast<std::size_t>(sources.inputChannels));
    std::iota(inputColumns.begin(), inputColumns.end(), 0);
    appendSourceColumns(
        raw, sources.input, sources.interfaceDofs,
        inputColumns, target, 0);
    diagnostics.rawInputColumns = sources.inputChannels;
    std::vector<int> boundaryColumns(
        static_cast<std::size_t>(sources.boundaryChannels));
    std::iota(boundaryColumns.begin(), boundaryColumns.end(), 0);
    appendSourceColumns(
        raw, sources.boundary, sources.interfaceDofs,
        boundaryColumns, target, 1);
    diagnostics.rawBoundaryColumns = sources.boundaryChannels;
    diagnostics.rawHistoryColumns = sources.historyChannels;
    const std::vector<int> historyColumns = topHistoryColumns(
        sources, target, metric,
        options.historyRank);
    appendSourceColumns(
        raw, sources.history, sources.interfaceDofs,
        historyColumns, target, 2);
    diagnostics.compressedHistoryColumns =
        static_cast<int>(historyColumns.size());

    for (std::size_t patchIndex = 0;
         patchIndex < selectedPatches.size(); ++patchIndex) {
        CompactSeed seed;
        seed.family = 3;
        seed.sourceColumn = static_cast<int>(patchIndex);
        seed.values.assign(target.size(), 0.0);
        const std::set<int> support(
            selectedPatches[patchIndex]->target.begin(),
            selectedPatches[patchIndex]->target.end());
        for (std::size_t row = 0; row < target.size(); ++row) {
            if (support.count(target[row]) != 0) {
                seed.values[row] =
                    metric[
                        static_cast<std::size_t>(target[row])];
            }
        }
        raw.push_back(std::move(seed));
    }
    CompactSeed graphLinear;
    graphLinear.family = 3;
    graphLinear.sourceColumn =
        static_cast<int>(selectedPatches.size());
    graphLinear.values.assign(target.size(), 0.0);
    for (std::size_t patchIndex = 0;
         patchIndex < selectedPatches.size(); ++patchIndex) {
        const double coordinate = selectedPatches.size() > 1
            ? -1.0 + 2.0 * static_cast<double>(patchIndex)
                / static_cast<double>(selectedPatches.size() - 1)
            : 0.0;
        const std::set<int> support(
            selectedPatches[patchIndex]->target.begin(),
            selectedPatches[patchIndex]->target.end());
        for (std::size_t row = 0; row < target.size(); ++row) {
            if (support.count(target[row]) != 0) {
                graphLinear.values[row] =
                    coordinate
                    * metric[
                        static_cast<std::size_t>(target[row])];
            }
        }
    }
    raw.push_back(std::move(graphLinear));
    diagnostics.rawGraphColumns =
        static_cast<int>(selectedPatches.size()) + 1;

    const int candidateLimit = options.candidateDimension > 0
        ? options.candidateDimension
        : std::min(
            exactSchur.size(),
            std::max(options.requestedRank + 12,
                     3 * options.requestedRank));
    const int seedLimit = std::min(
        candidateLimit,
        std::max(options.requestedRank + 4,
                 2 * options.requestedRank));
    const std::vector<CompactSeed> compactSeeds =
        pivotedSeedQr(
            raw, target, metric,
            seedLimit, options.deflationTolerance);
    diagnostics.independentSeedColumns =
        static_cast<int>(compactSeeds.size());
    diagnostics.deflatedSeedColumns =
        static_cast<int>(raw.size() - compactSeeds.size());
    if (compactSeeds.empty()) {
        throw std::runtime_error(
            "[Global coarse] All deterministic operator seeds deflated.");
    }

    std::vector<std::vector<double>> fullSeeds;
    fullSeeds.reserve(compactSeeds.size());
    for (const CompactSeed& compact : compactSeeds) {
        std::vector<double> seed(
            static_cast<std::size_t>(exactSchur.size()), 0.0);
        for (std::size_t row = 0; row < target.size(); ++row) {
            seed[static_cast<std::size_t>(target[row])] =
                compact.values[row];
        }
        fullSeeds.push_back(std::move(seed));
    }
    diagnostics.initialProbeResidual =
        probeResidual(fullSeeds, localSolver, nullptr);

    const std::vector<double>& schurDiagonal =
        exactSchur.diagonal();
    double diagonalScale = 0.0;
    for (double value : schurDiagonal) {
        diagonalScale = std::max(
            diagonalScale, std::abs(value));
    }
    const double diagonalFloor =
        1024.0 * std::numeric_limits<double>::epsilon()
        * std::max(
            std::numeric_limits<double>::min(),
            diagonalScale);
    const auto applyPreconditioner =
        [&](const std::vector<double>& input,
            std::vector<double>& output) {
            output.resize(input.size());
            for (std::size_t row = 0;
                 row < input.size(); ++row) {
                const double diagonal =
                    schurDiagonal[row];
                if (!(std::abs(diagonal) > diagonalFloor)
                    || !std::isfinite(diagonal)) {
                    throw std::runtime_error(
                        "[Global coarse] Schur-Jacobi diagonal "
                        "contains a tiny or nonfinite entry.");
                }
                output[row] = input[row] / diagonal;
            }
        };
    const auto applySpectralInverse =
        [&](const std::vector<double>& input,
            std::vector<double>& output) {
            if (options.inverseMode == "schur-jacobi") {
                applyPreconditioner(input, output);
                ++diagnostics.preconditionerApplyCount;
                return;
            }
            const auto solveStart = Clock::now();
            ExactInnerSolveResult solve =
                solveExactSchurPcg(
                    input, localSolver, exactSchur,
                    options.innerMaximumIterations,
                    options.innerTolerance,
                    diagnostics.schurApplyCount,
                    diagnostics.preconditionerApplyCount);
            diagnostics.innerSolveSeconds +=
                secondsSince(solveStart);
            ++diagnostics.innerSolveCount;
            diagnostics.innerSolveIterations +=
                solve.iterations;
            diagnostics.maximumInnerSolveResidual =
                std::max(
                    diagnostics.maximumInnerSolveResidual,
                    solve.relativeResidual);
            diagnostics.innerSolveConverged =
                diagnostics.innerSolveConverged
                && solve.converged;
            output = std::move(solve.solution);
            std::cout
                << "[Global coarse] exact inverse solve "
                << diagnostics.innerSolveCount
                << ": iterations=" << solve.iterations
                << ", residual="
                << solve.relativeResidual << '\n';
        };

    std::vector<double> x(
        static_cast<std::size_t>(exactSchur.size()), 0.0);
    std::vector<double> y(x.size(), 0.0);
    for (int row = 0; row < exactSchur.size(); ++row) {
        x[static_cast<std::size_t>(row)] =
            std::sin(0.7548776662466927 * (row + 1));
        y[static_cast<std::size_t>(row)] =
            std::cos(0.5698402909980532 * (row + 3));
    }
    std::vector<double> sx;
    std::vector<double> sy;
    exactSchur.apply(x, sx);
    exactSchur.apply(y, sy);
    diagnostics.schurApplyCount += 2;
    diagnostics.symmetryError =
        std::abs(dot(x.data(), sy.data(), exactSchur.size())
            - dot(y.data(), sx.data(), exactSchur.size()))
        / std::max({
            std::numeric_limits<double>::epsilon(),
            std::abs(dot(x.data(), sy.data(), exactSchur.size())),
            std::abs(dot(y.data(), sx.data(), exactSchur.size()))});
    std::vector<double> projectedX;
    std::vector<double> projectedTwiceX;
    std::vector<double> projectedY;
    localSolver.projectSchurEnergyComplement(
        x, projectedX);
    localSolver.projectSchurEnergyComplement(
        projectedX, projectedTwiceX);
    localSolver.projectSchurEnergyComplement(
        y, projectedY);
    diagnostics.schurApplyCount += 3;
    diagnostics.projectorIdempotenceError =
        relativeDifference(projectedTwiceX, projectedX);
    std::vector<double> sProjectedX;
    std::vector<double> sProjectedY;
    exactSchur.apply(projectedX, sProjectedX);
    exactSchur.apply(projectedY, sProjectedY);
    diagnostics.schurApplyCount += 2;
    diagnostics.projectorSSymmetryError =
        std::abs(
            dot(x.data(), sProjectedY.data(), exactSchur.size())
            - dot(y.data(), sProjectedX.data(), exactSchur.size()))
        / std::max({
            std::numeric_limits<double>::epsilon(),
            std::abs(dot(
                x.data(), sProjectedY.data(), exactSchur.size())),
            std::abs(dot(
                y.data(), sProjectedX.data(), exactSchur.size()))});

    CandidateSpace space;
    space.rows = exactSchur.size();
    for (const std::vector<double>& seed : fullSeeds) {
        std::vector<double> localResponse;
        localSolver.localGalerkinInterfaceResponse(
            seed, localResponse);
        std::vector<double> localImage;
        exactSchur.apply(localResponse, localImage);
        ++diagnostics.schurApplyCount;
        std::vector<double> residual(seed.size(), 0.0);
        for (std::size_t row = 0; row < seed.size(); ++row) {
            residual[row] = seed[row] - localImage[row];
        }
        std::vector<double> response;
        applySpectralInverse(residual, response);
        if (!diagnostics.innerSolveConverged) break;
        appendEnergyCandidate(
            space, std::move(response), localSolver, exactSchur,
            options.deflationTolerance,
            diagnostics.schurApplyCount);
        diagnostics.peakWorkingSetBytes = std::max(
            diagnostics.peakWorkingSetBytes, workingSetBytes());
        if (static_cast<int>(
                space.vectors.size()
                / static_cast<std::size_t>(space.rows))
            >= seedLimit) {
            break;
        }
    }

    for (int sweep = 1;
         sweep < options.krylovSweeps; ++sweep) {
        const int beginRank = static_cast<int>(
            space.vectors.size()
            / static_cast<std::size_t>(space.rows));
        for (int mode = 0; mode < beginRank; ++mode) {
            if (static_cast<int>(
                    space.vectors.size()
                    / static_cast<std::size_t>(space.rows))
                >= candidateLimit) {
                break;
            }
            std::vector<double> rhs(
                static_cast<std::size_t>(space.rows), 0.0);
            const double* basis = space.vectors.data()
                + static_cast<std::size_t>(mode * space.rows);
            for (int row = 0; row < space.rows; ++row) {
                rhs[static_cast<std::size_t>(row)] =
                    metric[
                        static_cast<std::size_t>(row)]
                    * basis[row];
            }
            std::vector<double> response;
            applySpectralInverse(rhs, response);
            if (!diagnostics.innerSolveConverged) break;
            appendEnergyCandidate(
                space, std::move(response),
                localSolver, exactSchur,
                options.deflationTolerance,
                diagnostics.schurApplyCount);
        }
        diagnostics.peakWorkingSetBytes = std::max(
            diagnostics.peakWorkingSetBytes, workingSetBytes());
        if (!diagnostics.innerSolveConverged) break;
    }
    if (!diagnostics.innerSolveConverged) {
        diagnostics.candidateRank = static_cast<int>(
            space.vectors.size()
            / static_cast<std::size_t>(space.rows));
        diagnostics.basisSeconds = secondsSince(start);
        diagnostics.peakWorkingSetBytes = std::max(
            diagnostics.peakWorkingSetBytes,
            workingSetBytes());
        diagnostics.peakIncrementalMemoryBytes =
            diagnostics.peakWorkingSetBytes
                > diagnostics.baselineWorkingSetBytes
            ? diagnostics.peakWorkingSetBytes
                - diagnostics.baselineWorkingSetBytes
            : 0;
        diagnostics.status =
            "exact_inner_schur_solve_residual_gate_failed";
        return model;
    }
    if (static_cast<int>(
            space.vectors.size()
            / static_cast<std::size_t>(space.rows))
        < options.requestedRank) {
        throw std::runtime_error(
            "[Global coarse] Candidate space is smaller than requested rank.");
    }

    RitzData ritz;
    for (int iteration = 0;
         iteration < options.maximumIterations; ++iteration) {
        ritz = computeRitz(
            space, metric,
            options.requestedRank);
        diagnostics.iterations = iteration + 1;
        diagnostics.maximumRitzResidual = 0.0;
        for (double residual : ritz.residuals) {
            diagnostics.maximumRitzResidual =
                std::max(
                    diagnostics.maximumRitzResidual, residual);
        }
        std::cout
            << "[Global coarse] Ritz iteration "
            << diagnostics.iterations
            << ", candidate rank="
            << space.vectors.size()
                / static_cast<std::size_t>(space.rows)
            << ", maximum residual="
            << diagnostics.maximumRitzResidual << '\n';
        if (diagnostics.maximumRitzResidual
                <= options.ritzTolerance
            || static_cast<int>(
                space.vectors.size()
                / static_cast<std::size_t>(space.rows))
                >= candidateLimit) {
            break;
        }
        bool added = false;
        for (int mode = 0; mode < ritz.rank; ++mode) {
            if (static_cast<int>(
                    space.vectors.size()
                    / static_cast<std::size_t>(space.rows))
                >= candidateLimit) {
                break;
            }
            if (ritz.residuals[static_cast<std::size_t>(mode)]
                <= options.ritzTolerance) {
                continue;
            }
            std::vector<double> correction;
            applySpectralInverse(
                ritz.residualVectors[
                    static_cast<std::size_t>(mode)],
                correction);
            if (!diagnostics.innerSolveConverged) break;
            added = appendEnergyCandidate(
                space, std::move(correction),
                localSolver, exactSchur,
                options.deflationTolerance,
                diagnostics.schurApplyCount) || added;
        }
        diagnostics.peakWorkingSetBytes = std::max(
            diagnostics.peakWorkingSetBytes, workingSetBytes());
        if (!added) break;
    }

    finalizeRitzModes(
        ritz, metric, localSolver, exactSchur,
        diagnostics);
    const OrthogonalityData overlap =
        localCoarseOverlap(localSolver, ritz);
    diagnostics.localCoarseOrthogonalityAbsolute =
        overlap.absolute;
    diagnostics.localCoarseOrthogonalityRelative =
        overlap.relative;
    diagnostics.maximumLocalCoarseOrthogonality =
        overlap.relative;
    model.rank = ritz.rank;
    model.basis = std::move(ritz.vectors);
    model.schurImages = std::move(ritz.images);
    model.eigenvalues = std::move(ritz.values);
    model.ritzResiduals = std::move(ritz.residuals);
    model.projectedRitzResiduals =
        std::move(ritz.projectedResiduals);
    model.inverseMetricRitzResiduals =
        std::move(ritz.inverseMetricResiduals);
    diagnostics.acceptedRank = model.rank;
    diagnostics.candidateRank = static_cast<int>(
        space.vectors.size()
        / static_cast<std::size_t>(space.rows));
    model.participationRatios.resize(
        static_cast<std::size_t>(model.rank), 0.0);
    for (int mode = 0; mode < model.rank; ++mode) {
        const double* basis = model.basis.data()
            + static_cast<std::size_t>(
                mode * model.fullInterfaceDofs);
        long double total = 0.0L;
        long double squared = 0.0L;
        for (const PortPatch& patch : physicalPatches) {
            long double energy = 0.0L;
            for (int row : patch.target) {
                energy += static_cast<long double>(basis[row])
                    * metric[
                        static_cast<std::size_t>(row)]
                    * basis[row];
            }
            total += energy;
            squared += energy * energy;
        }
        model.participationRatios[
            static_cast<std::size_t>(mode)] =
            squared > 0.0L
            ? static_cast<double>(total * total / squared)
            : 0.0;
    }
    diagnostics.finalProbeResidual =
        probeResidual(fullSeeds, localSolver, &model);
    diagnostics.basisSeconds = secondsSince(start);
    if (options.explicitReference) {
        const auto referenceStart = Clock::now();
        runExplicitProjectedReference(
            localModel, localSolver, exactSchur,
            metric, model, diagnostics);
        diagnostics.explicitReferenceSeconds =
            secondsSince(referenceStart);
    }
    diagnostics.peakWorkingSetBytes = std::max(
        diagnostics.peakWorkingSetBytes, workingSetBytes());
    diagnostics.peakIncrementalMemoryBytes =
        diagnostics.peakWorkingSetBytes
            > diagnostics.baselineWorkingSetBytes
        ? diagnostics.peakWorkingSetBytes
            - diagnostics.baselineWorkingSetBytes
        : 0;
    if (!diagnostics.innerSolveConverged) {
        diagnostics.status =
            "exact_inner_schur_solve_residual_gate_failed";
    } else if (model.rank != options.requestedRank) {
        diagnostics.status = "global_coarse_rank_gate_failed";
    } else if (!(diagnostics.symmetryError < 1.0e-12)) {
        diagnostics.status =
            "global_coarse_operator_asymmetry_gate_failed";
    } else if (!(diagnostics.maximumRitzResidual
                     < options.ritzTolerance)) {
        diagnostics.status =
            "global_coarse_spectral_residual_gate_failed";
    } else if (!(
                   diagnostics.maximumLocalCoarseOrthogonality
                   < options.orthogonalityTolerance)) {
        diagnostics.status =
            "global_coarse_orthogonality_gate_failed";
    } else if (!(diagnostics.coarseSchurGramError
                     < 1.0e-10)) {
        diagnostics.status =
            "global_coarse_schur_gram_gate_failed";
    } else if (!(diagnostics.basisSeconds < 300.0)) {
        diagnostics.status = "global_coarse_time_gate_failed";
    } else if (!(diagnostics.peakIncrementalMemoryBytes
                     < UINT64_C(1024) * 1024 * 1024)) {
        diagnostics.status = "global_coarse_memory_gate_failed";
    } else {
        diagnostics.status = "passed";
    }
    return model;
}

void writeGlobalInterfaceCoarseDiagnostics(
    const GlobalInterfaceCoarseModel& model,
    const std::string& caseName,
    const std::filesystem::path& outputDirectory)
{
    std::filesystem::create_directories(outputDirectory);
    const auto& row = model.diagnostics;
    std::ofstream summary(
        outputDirectory
        / "milestone8_global_coarse_prototype.csv");
    summary
        << "case,selected_interfaces,requested_coarse_rank,"
        "coarse_rank,candidate_rank,iterations,ritz_residual,"
        "schur_residual,local_coarse_orthogonality,"
        "symmetry_error,basis_time_s,"
        "peak_incremental_memory_bytes,process_peak_memory_bytes,"
        "initial_probe_residual,final_probe_residual,"
        "schur_apply_count,preconditioner_apply_count,"
        "requested_spectral_operator,actual_spectral_operator,"
        "jacobi_role,eigenvalue_mapping,"
        "inner_solver_requested,inner_solver_actual,"
        "inner_solve_count,inner_solve_iterations,"
        "inner_solve_time_s,maximum_inner_solve_residual,"
        "inner_solve_converged,projected_ritz_residual,"
        "inverse_metric_ritz_residual,"
        "local_coarse_orthogonality_absolute,"
        "local_coarse_orthogonality_relative,"
        "coarse_schur_gram_error,coarse_mass_gram_error,"
        "projector_idempotence_error,"
        "projector_s_symmetry_error,"
        "explicit_reference_requested,"
        "explicit_reference_ran,"
        "explicit_reference_dimension,"
        "explicit_reference_seconds,"
        "explicit_reference_full_residual,"
        "explicit_reference_projected_residual,"
        "explicit_eigenvalue_relative_error,"
        "explicit_max_principal_angle_radians,"
        "explicit_reference_status,"
        "metric_raw_minimum,metric_raw_maximum,metric_floor,"
        "metric_floored_rows,"
        "raw_input_columns,raw_boundary_columns,"
        "raw_history_columns,compressed_history_columns,"
        "raw_graph_columns,independent_seed_columns,"
        "deflated_seed_columns,snapshot_used,fom_used_for_basis,"
        "pod_used,svd_used,status\n"
        << std::setprecision(17)
        << caseName << ",\"";
    for (std::size_t index = 0;
         index < model.selectedInterfaceIds.size(); ++index) {
        if (index != 0) summary << ';';
        summary << model.selectedInterfaceIds[index];
    }
    summary << "\"," << row.requestedRank
        << ',' << model.rank
        << ',' << row.candidateRank
        << ',' << row.iterations
        << ',' << row.maximumRitzResidual
        << ',' << row.maximumRitzResidual
        << ',' << row.maximumLocalCoarseOrthogonality
        << ',' << row.symmetryError
        << ',' << row.basisSeconds
        << ',' << row.peakIncrementalMemoryBytes
        << ',' << row.peakWorkingSetBytes
        << ',' << row.initialProbeResidual
        << ',' << row.finalProbeResidual
        << ',' << row.schurApplyCount
        << ',' << row.preconditionerApplyCount
        << ',' << row.requestedSpectralOperator
        << ',' << row.actualSpectralOperator
        << ',' << row.jacobiRole
        << ',' << row.eigenvalueMapping
        << ',' << row.innerSolverRequested
        << ',' << row.innerSolverActual
        << ',' << row.innerSolveCount
        << ',' << row.innerSolveIterations
        << ',' << row.innerSolveSeconds
        << ',' << row.maximumInnerSolveResidual
        << ',' << (row.innerSolveConverged ? 1 : 0)
        << ',' << row.maximumProjectedRitzResidual
        << ',' << row.maximumInverseMetricRitzResidual
        << ',' << row.localCoarseOrthogonalityAbsolute
        << ',' << row.localCoarseOrthogonalityRelative
        << ',' << row.coarseSchurGramError
        << ',' << row.coarseMassGramError
        << ',' << row.projectorIdempotenceError
        << ',' << row.projectorSSymmetryError
        << ',' << (row.explicitReferenceRequested ? 1 : 0)
        << ',' << (row.explicitReferenceRan ? 1 : 0)
        << ',' << row.explicitReferenceDimension
        << ',' << row.explicitReferenceSeconds
        << ',' << row.explicitReferenceFullResidual
        << ',' << row.explicitReferenceProjectedResidual
        << ',' << row.explicitEigenvalueRelativeError
        << ',' << row.explicitMaximumPrincipalAngleRadians
        << ',' << row.explicitReferenceStatus
        << ',' << row.metricRawMinimum
        << ',' << row.metricRawMaximum
        << ',' << row.metricFloor
        << ',' << row.metricFlooredRows
        << ',' << row.rawInputColumns
        << ',' << row.rawBoundaryColumns
        << ',' << row.rawHistoryColumns
        << ',' << row.compressedHistoryColumns
        << ',' << row.rawGraphColumns
        << ',' << row.independentSeedColumns
        << ',' << row.deflatedSeedColumns
        << ",0,0,0,0," << row.status << '\n';

    std::ofstream modes(
        outputDirectory
        / "milestone8_global_coarse_modes.csv");
    modes
        << "mode,eigenvalue,ritz_residual,"
        "projected_ritz_residual,"
        "inverse_metric_ritz_residual,"
        "participation_ratio\n"
        << std::setprecision(17);
    for (int mode = 0; mode < model.rank; ++mode) {
        modes << mode
            << ',' << model.eigenvalues[
                static_cast<std::size_t>(mode)]
            << ',' << model.ritzResiduals[
                static_cast<std::size_t>(mode)]
            << ',' << model.projectedRitzResiduals[
                static_cast<std::size_t>(mode)]
            << ',' << model.inverseMetricRitzResiduals[
                static_cast<std::size_t>(mode)]
            << ',' << model.participationRatios[
                static_cast<std::size_t>(mode)]
            << '\n';
    }
}

} // namespace mor::transient
