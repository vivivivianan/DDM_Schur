#pragma once

// Compatibility declarations only. Matrix-free FGMRES/PCG implementations
// were removed from this export; DdmSchurSolver remains declared because the
// preserved upstream workflow contains unreachable research diagnostics.

#include "common.hpp"

#include <functional>
#include <memory>

namespace ddm_schur {

struct MatrixFreeFgmresResult {
    std::vector<double> solution;
    int iterations = 0;
    double initialRelativeResidual = 0.0;
    double relativeResidual = 0.0;
    bool converged = false;
    std::string actualSolver = "fgmres";
    bool fallbackTriggered = false;
    std::string fallbackReason;
    double operatorSeconds = 0.0;
    double preconditionerSeconds = 0.0;
    double orthogonalizationSeconds = 0.0;
    double vectorUpdateSeconds = 0.0;
    bool solutionProductAvailable = false;
};

using MatrixFreeApply = std::function<void(
    const std::vector<double>&, std::vector<double>&)>;

struct MatrixFreeFgmresWorkspace {
    std::vector<std::vector<double>> basis;
    std::vector<std::vector<double>> preconditioned;
    std::vector<std::vector<double>> hessenberg;
    std::vector<double> cosine;
    std::vector<double> sine;
    std::vector<double> rotatedRhs;
    std::vector<double> product;
    std::vector<double> residual;
    std::vector<double> work;
    std::vector<double> trialSolution;
    std::vector<double> trialProduct;
    std::vector<double> trialResidual;
};

// Stage-2 local-ROM entry point.  It deliberately uses the same restarted
// FGMRES algorithm, two-pass MGS, and true-residual acceptance gate as the
// validated Stage-1 Schur solve, while allowing a different matrix-free
// operator and preconditioner to be supplied.  When present, initialProduct
// must equal applyOperator(*initialSolution).
MatrixFreeFgmresResult solveMatrixFreeFgmres(
    const std::vector<double>& rhs,
    const Options& options,
    const MatrixFreeApply& applyOperator,
    const MatrixFreeApply& applyPreconditioner,
    const std::vector<double>* initialSolution = nullptr,
    MatrixFreeFgmresWorkspace* workspace = nullptr,
    const std::vector<double>* initialProduct = nullptr);

// Optional SPD path for the reduced full-interface Schur operator.  Every
// accepted result uses a freshly evaluated true residual.  Non-positive
// operator/preconditioner curvature is reported to the caller so it can
// fall back to the established FGMRES path without weakening the gate.
MatrixFreeFgmresResult solveMatrixFreePcg(
    const std::vector<double>& rhs,
    const Options& options,
    const MatrixFreeApply& applyOperator,
    const MatrixFreeApply& applyPreconditioner,
    const std::vector<double>* initialSolution = nullptr);

class DdmSchurSolver {
public:
    DdmSchurSolver(const Mesh& mesh,
                   const SparseMatrix& system,
                   const CaseConfig& physics,
                   Options options = {});
    ~DdmSchurSolver();

    DdmSchurSolver(DdmSchurSolver&&) noexcept;
    DdmSchurSolver& operator=(DdmSchurSolver&&) noexcept;
    DdmSchurSolver(const DdmSchurSolver&) = delete;
    DdmSchurSolver& operator=(const DdmSchurSolver&) = delete;

    SolveResult solve(const std::vector<double>& globalRhs);
    // Stage-2 MOR accessors.  These delegate to the same exact matrix-free
    // Schur operator and retained PARDISO factors used by solve().
    int interfaceDofs() const;
    const std::vector<int>& interfaceGlobalDofs() const;
    std::vector<double> condensedRhs(const std::vector<double>& globalRhs);
    void applyExactSchur(const std::vector<double>& interfaceVector,
                         std::vector<double>& result);
    std::vector<double> recover(const std::vector<double>& globalRhs,
                                const std::vector<double>& interfaceSolution);
    SolveResult solveCorrection(const std::vector<double>& globalRhs,
                                const std::vector<double>& interfaceInitialGuess);
    const Report& setupReport() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace ddm_schur
