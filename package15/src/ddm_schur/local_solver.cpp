#include "../sipg_core.hpp"
#include "../linear_solvers.hpp"
#include "local_solver.hpp"

#include <algorithm>
#include <stdexcept>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

namespace ddm_schur {

struct LocalSolver::Impl {
    int size = 0;
    int pardisoThreads = 1;
    std::unique_ptr<SubdomainDirectSolver> solver;
};

namespace {

class ScopedLocalMklThreads {
public:
    explicit ScopedLocalMklThreads(int threads)
    {
#ifdef USE_MKL_PARDISO
        previous_ = mkl_get_max_threads();
        mkl_set_num_threads_local(std::max(1, threads));
#else
        (void)threads;
#endif
    }
    ~ScopedLocalMklThreads()
    {
#ifdef USE_MKL_PARDISO
        mkl_set_num_threads_local(previous_);
#endif
    }
private:
#ifdef USE_MKL_PARDISO
    int previous_ = 1;
#endif
};

} // namespace

LocalSolver::LocalSolver() : impl_(std::make_unique<Impl>()) {}

LocalSolver::LocalSolver(int size, const std::vector<Entry>& entries, int pardisoThreads)
    : impl_(std::make_unique<Impl>())
{
    impl_->size = size;
    impl_->pardisoThreads = std::max(1, pardisoThreads);
    if (size == 0) {
        return;
    }
    std::vector<MatrixEntry> nativeEntries;
    nativeEntries.reserve(entries.size());
    for (const Entry& entry : entries) {
        nativeEntries.push_back({entry.row, entry.col, entry.value});
    }
    ScopedLocalMklThreads threads(impl_->pardisoThreads);
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
    ScopedLocalMklThreads threads(impl_->pardisoThreads);
    impl_->solver->solve(rhs, solution);
}

void LocalSolver::solveMultiple(const std::vector<double>& rhs,
                                int rightHandSides,
                                std::vector<double>& solution)
{
    if (rightHandSides <= 0
        || rhs.size() != static_cast<std::size_t>(impl_->size)
            * static_cast<std::size_t>(rightHandSides)) {
        throw std::runtime_error(
            "[Schur] Local PARDISO multi-RHS dimensions are invalid.");
    }
    if (impl_->size == 0) {
        solution.clear();
        return;
    }
    ScopedLocalMklThreads threads(impl_->pardisoThreads);
    impl_->solver->solveMultiple(rhs, rightHandSides, solution);
}

} // namespace ddm_schur
