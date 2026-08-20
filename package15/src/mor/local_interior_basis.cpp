#include "local_interior_basis.hpp"

#include "pod_basis.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

namespace mor {

void buildCompressedInteriorBases(DeploymentResponseModel& model,
                                  int requestedRank,
                                  double discardedEnergyTolerance,
                                  double relativeSingularTolerance)
{
    const auto start = std::chrono::steady_clock::now();
    for (InteriorResponseBlock& block : model.interiors) {
        const int rows = static_cast<int>(block.globalDofs.size());
        if (rows == 0) {
            continue;
        }
#ifdef USE_MKL_PARDISO
        const int columns = model.sourceChannels;
        const int availableRank = std::min(rows, columns);
        std::vector<double> responseCopy = block.exactResponse;
        std::vector<double> singularValues(static_cast<std::size_t>(availableRank), 0.0);
        std::vector<double> leftVectors(static_cast<std::size_t>(
            rows * availableRank), 0.0);
        std::vector<double> rightVectors(static_cast<std::size_t>(
            availableRank * columns), 0.0);
        const lapack_int info = LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'S',
            static_cast<lapack_int>(rows), static_cast<lapack_int>(columns),
            responseCopy.data(), static_cast<lapack_int>(rows),
            singularValues.data(), leftVectors.data(), static_cast<lapack_int>(rows),
            rightVectors.data(), static_cast<lapack_int>(availableRank));
        if (info != 0) {
            throw std::runtime_error("[MOR Interior] Thin SVD failed for subdomain "
                + std::to_string(block.subdomain) + " with info=" + std::to_string(info));
        }

        long double totalEnergy = 0.0L;
        for (double sigma : singularValues) {
            totalEnergy += static_cast<long double>(sigma) * sigma;
        }
        block.singularValues = singularValues;
        block.retainedEnergy.resize(static_cast<std::size_t>(availableRank), 0.0);
        long double cumulativeEnergy = 0.0L;
        int numericalRank = 0;
        int energyRank = availableRank;
        bool energyRankFound = false;
        const double sigma0 = singularValues.empty() ? 0.0 : singularValues.front();
        for (int mode = 0; mode < availableRank; ++mode) {
            const double sigma = singularValues[static_cast<std::size_t>(mode)];
            if (sigma0 > 0.0 && sigma > relativeSingularTolerance * sigma0) {
                ++numericalRank;
            }
            cumulativeEnergy += static_cast<long double>(sigma) * sigma;
            const double retained = totalEnergy > 0.0L
                ? static_cast<double>(cumulativeEnergy / totalEnergy) : 1.0;
            block.retainedEnergy[static_cast<std::size_t>(mode)] = retained;
            if (!energyRankFound && 1.0 - retained <= discardedEnergyTolerance) {
                energyRank = mode + 1;
                energyRankFound = true;
            }
        }
        // A user-specified fixed rank is authoritative.  In particular, an
        // exact-rank audit must retain small but physically important response
        // directions that an energy or relative-sigma test would discard.
        block.rank = requestedRank > 0
            ? std::min(requestedRank, availableRank)
            : std::min(energyRank, std::max(1, numericalRank));
        block.localBasis.assign(leftVectors.begin(), leftVectors.begin()
            + static_cast<std::ptrdiff_t>(rows * block.rank));
        block.localCoordinateMap.assign(static_cast<std::size_t>(
            block.rank * columns), 0.0);
        for (int channel = 0; channel < columns; ++channel) {
            for (int mode = 0; mode < block.rank; ++mode) {
                block.localCoordinateMap[static_cast<std::size_t>(
                    channel * block.rank + mode)] =
                    singularValues[static_cast<std::size_t>(mode)]
                    * rightVectors[static_cast<std::size_t>(
                        channel * availableRank + mode)];
            }
        }
#else
        SnapshotDatabase snapshots;
        snapshots.rows = rows;
        snapshots.values = block.exactResponse;
        snapshots.cases.resize(static_cast<std::size_t>(model.sourceChannels));
        for (int channel = 0; channel < model.sourceChannels; ++channel) {
            snapshots.cases[static_cast<std::size_t>(channel)].index = channel;
        }
        const PodResult pod = buildGramPod(snapshots, requestedRank,
            discardedEnergyTolerance, relativeSingularTolerance);
        block.rank = pod.selectedRank;
        block.localBasis = pod.basis;
        block.singularValues = pod.singularValues;
        block.retainedEnergy = pod.retainedEnergy;
        block.localCoordinateMap.assign(static_cast<std::size_t>(
            block.rank * model.sourceChannels), 0.0);
        for (int channel = 0; channel < model.sourceChannels; ++channel) {
            for (int mode = 0; mode < block.rank; ++mode) {
                double value = 0.0;
                for (int row = 0; row < rows; ++row) {
                    value += block.localBasis[static_cast<std::size_t>(mode * rows + row)]
                        * block.exactResponse[static_cast<std::size_t>(channel * rows + row)];
                }
                block.localCoordinateMap[static_cast<std::size_t>(
                    channel * block.rank + mode)] = value;
            }
        }
#endif
    }
    model.compressionSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
}

} // namespace mor
