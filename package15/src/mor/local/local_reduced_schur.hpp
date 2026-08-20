#pragma once

#include "local_subdomain_model.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

struct CaseConfig;
struct Mesh;
namespace ddm_schur { struct InterfacePartition; struct Options; }

namespace mor::local {

class PortCoreSolver;

struct DenseSymmetricFactor {
    int size = 0;
    bool cholesky = true;
    std::vector<double> lower;
    std::vector<double> diagonal;
};

DenseSymmetricFactor factorDenseSymmetric(const std::vector<double>& matrix, int size);
void solveDenseSymmetric(const DenseSymmetricFactor& factor,
                         std::vector<double>& rightHandSide);
// Solve a row-major size-by-rightHandSides dense block in place.  This is
// algebraically identical to repeated solveDenseSymmetric calls, but permits
// a level-3 BLAS triangular solve when MKL is available.
void solveDenseSymmetricMultiple(const DenseSymmetricFactor& factor,
                                 std::vector<double>& rightHandSides,
                                 int rightHandSideCount);

class LocalReducedSchurSolver {
public:
    explicit LocalReducedSchurSolver(const Model& model);
    LocalReducedSchurSolver(const Model& model,
                            const ::Mesh& mesh,
                            const ::CaseConfig& physics,
                            const ddm_schur::InterfacePartition& partition,
                            const ddm_schur::Options& options,
                            int matrixFreeInterfaceThreshold,
                            const std::filesystem::path& outputDirectory);
    ~LocalReducedSchurSolver();

    SolveResult solve(
        const std::vector<double>& globalRhs,
        const std::vector<double>* interfaceInitialGuess = nullptr,
        double interfaceToleranceOverride = 0.0);
    // Solve from already projected interior RHS blocks and a physical-order
    // interface RHS.  This is equivalent to solve(globalRhs) while avoiding
    // a redundant global-to-local basis projection in transient runs.
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
    int coarseDimension() const;
    int geometricCoarseDimension() const;
    int operatorCoarseDimension() const;
    bool operatorCoarseCacheHit() const;
    double operatorCoarseSetupSeconds() const;
    double operatorCoarseCacheLoadSeconds() const;
    double operatorCoarseCacheSaveSeconds() const;
    int proxyColors() const;
    int proxyProbingApplies() const;
    std::size_t proxyNonzeros() const;
    double proxySetupSeconds() const;
    double proxySymbolicSeconds() const;
    double proxyNumericalSeconds() const;
    std::size_t proxyMemoryBytes() const;
    bool proxyMatrixCacheHit() const;
    bool proxyFactorCacheHit() const;
    bool portCoreCacheHit() const;
    double portCoreCacheLoadSeconds() const;
    double portCoreCacheSaveSeconds() const;
    std::size_t portCoreCacheBytes() const;
    double portCorePartitionSeconds() const;
    double portCoreCouplingAssemblySeconds() const;
    double portCoreLeafCsrSeconds() const;
    double portCoreLeafFactorSeconds() const;
    double portCoreEliminationSeconds() const;
    double portCoreMultiRhsPortSeconds() const;
    double portCoreSchurProductPortSeconds() const;
    double portCoreCoreAccumulationSeconds() const;
    double portCoreCoreCsrSeconds() const;
    double portCoreCoreFactorSeconds() const;

private:
    struct SparseFactor;
    struct AugmentedDirectData;
    struct MatrixFreeData;
    void initializeExplicit();
    void initializeAugmentedDirect(
        const ddm_schur::Options& options,
        const std::filesystem::path& outputDirectory);
    void initializeMatrixFree(const ::Mesh& mesh,
                              const ::CaseConfig& physics,
                              const ddm_schur::InterfacePartition& partition,
                              const ddm_schur::Options& options,
                              const std::filesystem::path& outputDirectory);
    void initializePortCore(
        const ddm_schur::InterfacePartition& partition,
        const ddm_schur::Options& options,
        const std::filesystem::path& outputDirectory);
    void applyMatrixFree(const std::vector<double>& input,
                         std::vector<double>& output,
                         bool countMatvec);
    void applyMatrixFreePreconditioner(const std::vector<double>& residual,
                                       std::vector<double>& result);
    SolveResult solveWithReducedRhs(
        std::vector<std::vector<double>> projectedInteriorRhs,
        std::vector<double> interfaceRhs,
        const std::vector<double>* interfaceInitialGuess,
        double interfaceToleranceOverride);
    const Model& model_;
    std::vector<DenseSymmetricFactor> localFactors_;
    DenseSymmetricFactor schurFactor_;
    std::unique_ptr<SparseFactor> sparseFactor_;
    std::unique_ptr<AugmentedDirectData> augmentedDirect_;
    std::unique_ptr<MatrixFreeData> matrixFree_;
    std::unique_ptr<PortCoreSolver> portCore_;
    std::string interfaceSolver_ = "dense-llt";
    std::size_t interfaceNonzeros_ = 0;
    double assemblySeconds_ = 0.0;
    double factorizationSeconds_ = 0.0;
};

} // namespace mor::local
