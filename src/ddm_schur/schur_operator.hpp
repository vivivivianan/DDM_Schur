#pragma once

#include "common.hpp"

#include <memory>
#include <vector>

namespace ddm_schur {

class SchurOperator {
public:
    SchurOperator(const Mesh& mesh,
                  const SparseMatrix& system,
                  bool coarseLinearXY,
                  bool coarseLinearZ,
                  bool coarseGlobalQuadraticZ,
                  bool coarseInterfacePatches,
                  bool coarseInterfaceLinearXY);
    ~SchurOperator();

    SchurOperator(SchurOperator&&) noexcept;
    SchurOperator& operator=(SchurOperator&&) noexcept;
    SchurOperator(const SchurOperator&) = delete;
    SchurOperator& operator=(const SchurOperator&) = delete;

    int domains() const;
    int totalDofs() const;
    int interfaceDofs() const;
    int interiorDofs() const;
    int coarseDimension() const;
    int interfacePatchCount() const;
    double localFactorizationSeconds() const;
    double localSymbolicAnalysisSeconds() const;
    double localNumericalFactorizationSeconds() const;
    double localSolveSeconds() const;
    double coarseSolveSeconds() const;
    std::size_t memoryBytes() const;
    int localSolveCalls() const;
    int localSymbolicAnalysisCalls() const;
    int localNumericalFactorizationCalls() const;
    int matvecCalls() const;

    std::vector<double> condensedRhs(const std::vector<double>& globalRhs);
    void apply(const std::vector<double>& interfaceVector, std::vector<double>& result);
    void applyBlockPreconditioner(const std::vector<double>& residual, std::vector<double>& result);
    std::vector<double> recover(const std::vector<double>& globalRhs,
                                const std::vector<double>& interfaceSolution);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace ddm_schur
