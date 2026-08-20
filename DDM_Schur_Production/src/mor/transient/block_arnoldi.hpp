#pragma once

// Public contracts for local Krylov basis construction and diagnostics. Dense
// basis columns are column-major so projection and BLAS calls share one layout.

#include "thermal_descriptor_system.hpp"

#include <cstddef>
#include <vector>

namespace mor::transient {

struct ArnoldiHistoryRow {
    int moment = 0;
    int inputColumns = 0;
    int addedRank = 0;
    int cumulativeRank = 0;
    int deflatedColumns = 0;
    double orthogonalityError = 0.0;
    double arnoldiResidual = 0.0;
    double solveSeconds = 0.0;
    double orthogonalizationSeconds = 0.0;
    std::size_t basisBytes = 0;
};

struct BlockArnoldiTiming {
    double symbolicAnalysisSeconds = 0.0;
    double numericalFactorizationSeconds = 0.0;
    double multiRhsSolveSeconds = 0.0;
    double orthogonalizationSeconds = 0.0;
    double totalSeconds = 0.0;
    int symbolicAnalysisCalls = 0;
    int numericalFactorizationCalls = 0;
};

struct BlockArnoldiResult {
    int rows = 0;
    int rank = 0;
    int blockSize = 0;
    int moments = 0;
    double expansionPoint = 0.0;
    double rankTolerance = 1.0e-10;
    std::vector<double> basis;
    std::vector<double> referenceTemperature;
    std::vector<ArnoldiHistoryRow> history;
    BlockArnoldiTiming timing;
};

BlockArnoldiResult buildBlockArnoldiBasis(
    const ThermalDescriptorSystem& system,
    int moments,
    double expansionPoint,
    double rankTolerance,
    double secondMomentEnergy = 1.0,
    int secondMomentMaximumColumns = 0,
    int directSolverThreads = 1);

} // namespace mor::transient
