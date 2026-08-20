#pragma once

#include "interface_operator.hpp"

#include <cstddef>
#include <memory>
#include <vector>

namespace ddm_schur {

class LocalSolver {
public:
    LocalSolver();
    LocalSolver(int size, const std::vector<Entry>& entries, int pardisoThreads = 1);
    ~LocalSolver();

    LocalSolver(LocalSolver&&) noexcept;
    LocalSolver& operator=(LocalSolver&&) noexcept;
    LocalSolver(const LocalSolver&) = delete;
    LocalSolver& operator=(const LocalSolver&) = delete;

    int size() const;
    std::size_t memoryBytes() const;
    double symbolicAnalysisSeconds() const;
    double numericalFactorizationSeconds() const;
    int symbolicAnalysisCalls() const;
    int numericalFactorizationCalls() const;
    void solve(const std::vector<double>& rhs, std::vector<double>& solution);
    void solveMultiple(const std::vector<double>& rhs,
                       int rightHandSides,
                       std::vector<double>& solution);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace ddm_schur
