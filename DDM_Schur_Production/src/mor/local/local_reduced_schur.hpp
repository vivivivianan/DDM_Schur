#pragma once

// Production reduced Schur solver.
//
// The research repository contained four interface solution strategies in
// this class: explicit Schur assembly, matrix-free FGMRES/PCG, port/core
// elimination, and an augmented sparse direct solve. Package15 measurements
// selected the last strategy, so this export keeps only that implementation.
// The solved unknown vector is
//
//     y = [T_Gamma, q_1, q_2, ..., q_N],
//
// where T_Gamma contains every physical interface temperature and q_i is the
// reduced interior coordinate vector of subdomain i. A single symmetric CSR
// matrix is assembled and factorized once for a fixed time step; every online
// step then changes only its right-hand side.

#include "local_subdomain_model.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

struct CaseConfig;
struct Mesh;
namespace ddm_schur { struct InterfacePartition; struct Options; }

namespace mor::local {

// Small dense factors are still used by construction-trace helpers in the
// Dynamic Schur workflow. Storage is row-major lower triangular; a failed
// Cholesky attempt falls back to an LDL^T factorization.
struct DenseSymmetricFactor {
    int size = 0;
    bool cholesky = true;
    std::vector<double> lower;
    std::vector<double> diagonal;
};

DenseSymmetricFactor factorDenseSymmetric(
    const std::vector<double>& matrix, int size);
void solveDenseSymmetric(
    const DenseSymmetricFactor& factor,
    std::vector<double>& rightHandSide);
void solveDenseSymmetricMultiple(
    const DenseSymmetricFactor& factor,
    std::vector<double>& rightHandSides,
    int rightHandSideCount);

class LocalReducedSchurSolver {
public:
    // The mesh, physics, and partition arguments remain in the signature so
    // the validated upstream workflow can call this clean implementation
    // without duplicating model construction. Augmented-direct itself needs
    // only the already assembled local ROM model and solver thread settings.
    LocalReducedSchurSolver(
        const Model& model,
        const ::Mesh& mesh,
        const ::CaseConfig& physics,
        const ddm_schur::InterfacePartition& partition,
        const ddm_schur::Options& options,
        int matrixFreeInterfaceThreshold,
        const std::filesystem::path& outputDirectory);
    ~LocalReducedSchurSolver();

    // Project a full global RHS into local reduced coordinates, solve the
    // coupled augmented system, and reconstruct the temperature vector.
    SolveResult solve(
        const std::vector<double>& globalRhs,
        const std::vector<double>* interfaceInitialGuess = nullptr,
        double interfaceToleranceOverride = 0.0);

    // Fast transient entry point. The caller has already formed V_i^T f_i,
    // so this avoids a repeated global-to-local projection at each step.
    SolveResult solveReducedRhs(
        const std::vector<std::vector<double>>& projectedInteriorRhs,
        const std::vector<double>& interfaceRhs,
        const std::vector<double>* interfaceInitialGuess = nullptr,
        double interfaceToleranceOverride = 0.0);

    double factorizationSeconds() const { return factorizationSeconds_; }
    double assemblySeconds() const { return assemblySeconds_; }
    double symbolicAnalysisSeconds() const;
    double numericalFactorizationSeconds() const;
    int symbolicAnalysisCalls() const;
    int numericalFactorizationCalls() const;
    std::size_t factorMemoryBytes() const;
    std::size_t interfaceNonzeros() const { return interfaceNonzeros_; }
    const std::string& interfaceSolver() const { return interfaceSolver_; }

    // The summary CSV inherited these columns from the research executable.
    // Unsupported methods are absent in this project, therefore their
    // compatibility accessors intentionally report zero/false.
    int coarseDimension() const { return 0; }
    int geometricCoarseDimension() const { return 0; }
    int operatorCoarseDimension() const { return 0; }
    bool operatorCoarseCacheHit() const { return false; }
    double operatorCoarseSetupSeconds() const { return 0.0; }
    double operatorCoarseCacheLoadSeconds() const { return 0.0; }
    double operatorCoarseCacheSaveSeconds() const { return 0.0; }
    int proxyColors() const { return 0; }
    int proxyProbingApplies() const { return 0; }
    std::size_t proxyNonzeros() const { return 0; }
    double proxySetupSeconds() const { return 0.0; }
    double proxySymbolicSeconds() const { return 0.0; }
    double proxyNumericalSeconds() const { return 0.0; }
    std::size_t proxyMemoryBytes() const { return 0; }
    bool proxyMatrixCacheHit() const { return false; }
    bool proxyFactorCacheHit() const { return false; }
    bool portCoreCacheHit() const { return false; }
    double portCoreCacheLoadSeconds() const { return 0.0; }
    double portCoreCacheSaveSeconds() const { return 0.0; }
    std::size_t portCoreCacheBytes() const { return 0; }
    double portCorePartitionSeconds() const { return 0.0; }
    double portCoreCouplingAssemblySeconds() const { return 0.0; }
    double portCoreLeafCsrSeconds() const { return 0.0; }
    double portCoreLeafFactorSeconds() const { return 0.0; }
    double portCoreEliminationSeconds() const { return 0.0; }
    double portCoreMultiRhsPortSeconds() const { return 0.0; }
    double portCoreSchurProductPortSeconds() const { return 0.0; }
    double portCoreCoreAccumulationSeconds() const { return 0.0; }
    double portCoreCoreCsrSeconds() const { return 0.0; }
    double portCoreCoreFactorSeconds() const { return 0.0; }

private:
    struct AugmentedDirectData;

    void initializeAugmentedDirect(
        const ddm_schur::Options& options,
        const std::filesystem::path& outputDirectory);
    SolveResult solveWithReducedRhs(
        std::vector<std::vector<double>> projectedInteriorRhs,
        std::vector<double> interfaceRhs);

    const Model& model_;
    std::unique_ptr<AugmentedDirectData> augmentedDirect_;
    std::string interfaceSolver_ = "augmented-pardiso";
    std::size_t interfaceNonzeros_ = 0;
    double assemblySeconds_ = 0.0;
    double factorizationSeconds_ = 0.0;
};

} // namespace mor::local
