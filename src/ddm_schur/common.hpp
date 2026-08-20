#pragma once

#include <cstddef>
#include <string>
#include <vector>

struct Mesh;
struct SparseMatrix;

namespace ddm_schur {

struct Options {
    int maxIterations = 500;
    int restart = 50;
    double relativeTolerance = 1.0e-12;
    bool coarseLinearXY = true;
    bool coarseLinearZ = true;
    bool coarseGlobalQuadraticZ = false;
    bool coarseInterfacePatches = false;
    bool coarseInterfaceLinearXY = false;
};

struct Report {
    int domains = 0;
    int totalDofs = 0;
    int interfaceDofs = 0;
    int interiorDofs = 0;
    int iterations = 0;
    int coarseDimension = 0;
    int interfacePatchCount = 0;
    int schurMatvecs = 0;
    int localSolveCalls = 0;
    int localSymbolicAnalysisCalls = 0;
    int localNumericalFactorizationCalls = 0;
    double setupSeconds = 0.0;
    double localFactorizationSeconds = 0.0;
    double localSymbolicAnalysisSeconds = 0.0;
    double localNumericalFactorizationSeconds = 0.0;
    double localSolveSeconds = 0.0;
    double coarseSolveSeconds = 0.0;
    double condensedRhsSeconds = 0.0;
    double interfaceSolveSeconds = 0.0;
    double fgmresSeconds = 0.0;
    double recoverySeconds = 0.0;
    double totalSolveSeconds = 0.0;
    double totalSeconds = 0.0;
    double interfaceRelativeResidual = 0.0;
    std::size_t memoryBytes = 0;
    std::string status = "not_run";
};

struct SolveResult {
    std::vector<double> temperature;
    Report report;
};

} // namespace ddm_schur
