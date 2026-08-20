#include "../sipg_core.hpp"
#include "../linear_solvers.hpp"
#include "local_solver.hpp"

#include <stdexcept>

namespace ddm_schur {

struct LocalSolver::Impl {
    int size = 0;
    std::unique_ptr<SubdomainDirectSolver> solver;
};

LocalSolver::LocalSolver() : impl_(std::make_unique<Impl>()) {}

LocalSolver::LocalSolver(int size, const std::vector<Entry>& entries)
    : impl_(std::make_unique<Impl>())
{
    impl_->size = size;
    if (size == 0) {
        return;
    }
    std::vector<MatrixEntry> nativeEntries;
    nativeEntries.reserve(entries.size());
    for (const Entry& entry : entries) {
        nativeEntries.push_back({entry.row, entry.col, entry.value});
    }
    impl_->solver = std::make_unique<SubdomainDirectSolver>(size, nativeEntries);
}

LocalSolver::~LocalSolver() = default;
LocalSolver::LocalSolver(LocalSolver&&) noexcept = default;
LocalSolver& LocalSolver::operator=(LocalSolver&&) noexcept = default;

int LocalSolver::size() const { return impl_->size; }

std::size_t LocalSolver::memoryBytes() const
{
    return impl_->solver ? impl_->solver->memoryBytes() : 0;
}

double LocalSolver::symbolicAnalysisSeconds() const
{
    return impl_->solver ? impl_->solver->symbolicAnalysisSeconds() : 0.0;
}

double LocalSolver::numericalFactorizationSeconds() const
{
    return impl_->solver ? impl_->solver->numericalFactorizationSeconds() : 0.0;
}

int LocalSolver::symbolicAnalysisCalls() const
{
    return impl_->solver ? impl_->solver->symbolicAnalysisCalls() : 0;
}

int LocalSolver::numericalFactorizationCalls() const
{
    return impl_->solver ? impl_->solver->numericalFactorizationCalls() : 0;
}

void LocalSolver::solve(const std::vector<double>& rhs, std::vector<double>& solution)
{
    if (static_cast<int>(rhs.size()) != impl_->size) {
        throw std::runtime_error("[Schur] Local PARDISO right-hand side has the wrong size.");
    }
    if (impl_->size == 0) {
        solution.clear();
        return;
    }
    impl_->solver->solve(rhs, solution);
}

} // namespace ddm_schur
