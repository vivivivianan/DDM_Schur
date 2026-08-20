#pragma once

#include <functional>
#include <string>
#include <vector>

namespace mor::transient {

struct MatrixFreeEigenOptions {
    int requestedEigenpairs = 16;
    int oversamplingVectors = 4;
    int maximumIterations = 40;
    double relativeTolerance = 1.0e-8;
    double deflationTolerance = 1.0e-12;
    // Zero disables the wall-time guard. A positive value stops at the next
    // operator boundary, without weakening residual tolerances.
    double maximumWallSeconds = 0.0;
};

struct MatrixFreeEigenResult {
    int dimension = 0;
    int rank = 0;
    int iterations = 0;
    int operatorApplies = 0;
    bool converged = false;
    std::string status = "not_run";
    std::vector<double> eigenvalues;
    std::vector<double> residuals;
    // Column-major dimension x rank.
    std::vector<double> eigenvectors;
};

using SymmetricOperatorApply =
    std::function<void(const std::vector<double>&, std::vector<double>&)>;

MatrixFreeEigenResult solveLargestSymmetricEigenpairs(
    int dimension,
    const SymmetricOperatorApply& apply,
    const MatrixFreeEigenOptions& options);

} // namespace mor::transient
