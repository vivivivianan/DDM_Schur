#pragma once

#include "types.hpp"

#include <vector>

namespace ddm_schur { class DdmSchurSolver; }

namespace mor {

ReducedSchurModel buildReducedSchurModel(ddm_schur::DdmSchurSolver& solver,
                                         const PodResult& pod,
                                         int rank,
                                         const std::vector<double>& referenceInterface,
                                         const Fingerprints& fingerprints);

ReducedSchurModel truncateModel(const ReducedSchurModel& source, int rank);

class ReducedSchurOnlineSolver {
public:
    ReducedSchurOnlineSolver(ddm_schur::DdmSchurSolver& solver,
                             const ReducedSchurModel& model);

    SolveResult solve(const std::vector<double>& globalRhs, bool corrected);
    const char* factorizationType() const;

private:
    ddm_schur::DdmSchurSolver& solver_;
    const ReducedSchurModel& model_;
    std::vector<double> referenceSchurImage_;
    std::vector<double> denseFactor_;
    std::vector<double> diagonalFactor_;
    bool choleskyFactor_ = false;
};

} // namespace mor
