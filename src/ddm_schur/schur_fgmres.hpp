#pragma once

#include "common.hpp"

#include <memory>

namespace ddm_schur {

class DdmSchurSolver {
public:
    DdmSchurSolver(const Mesh& mesh, const SparseMatrix& system, Options options = {});
    ~DdmSchurSolver();

    DdmSchurSolver(DdmSchurSolver&&) noexcept;
    DdmSchurSolver& operator=(DdmSchurSolver&&) noexcept;
    DdmSchurSolver(const DdmSchurSolver&) = delete;
    DdmSchurSolver& operator=(const DdmSchurSolver&) = delete;

    SolveResult solve(const std::vector<double>& globalRhs);
    const Report& setupReport() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace ddm_schur
