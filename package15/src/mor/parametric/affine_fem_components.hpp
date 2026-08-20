#pragma once

#include "affine_parameter.hpp"
#include "mor/source_parameterization.hpp"
#include "sipg_core.hpp"

#include <cstdint>
#include <filesystem>
#include <vector>

namespace mor::parametric {

struct AffineFemComponents {
    AffineParameter parameter;
    SparseMatrix matrixConstant;
    SparseMatrix matrixLinear;
    std::vector<SparseMatrix> matrixHarmonic;
    std::vector<double> rhsConstant;
    std::vector<double> rhsLinear;
    std::vector<std::vector<double>> rhsHarmonic;
    SourceParameterization sources;
    std::uint64_t constantMatrixHash = 0;
    std::uint64_t linearMatrixHash = 0;
    std::vector<std::uint64_t> harmonicMatrixHashes;
    std::uint64_t constantRhsHash = 0;
    std::uint64_t linearRhsHash = 0;
    std::vector<std::uint64_t> harmonicRhsHashes;
    double assemblySeconds = 0.0;
};

struct DirectParametricSystem {
    SparseMatrix matrix;
    std::vector<double> rhs;
    std::vector<double> rawSource;
    std::vector<double> heatOnlySource;
    std::vector<double> fixedAdjust;
};

struct AffineValidationRow {
    double parameterValue = 0.0;
    double matrixRelativeDifference = 0.0;
    double rhsRelativeDifference = 0.0;
    double symmetryError = 0.0;
    double minimumDiagonal = 0.0;
    std::size_t nonzeros = 0;
};

AffineFemComponents buildAffineFemComponents(
    const Mesh& mesh,
    const CaseConfig& physics,
    const SparseMatrix& referenceSystem,
    const std::vector<double>& assembledSource,
    const std::vector<double>& heatOnlySource,
    const std::vector<double>& fixedAdjust,
    const mor::Options& options);

DirectParametricSystem assembleDirectParametricSystem(
    const Mesh& mesh,
    const CaseConfig& referencePhysics,
    const AffineParameter& parameter,
    double value);

SparseMatrix composeMatrix(const AffineFemComponents& components, double value);
std::vector<double> composeRhs(const AffineFemComponents& components,
                               double value,
                               const std::vector<double>& powersW);

AffineValidationRow validateAffineOperator(
    const Mesh& mesh,
    const CaseConfig& physics,
    const AffineFemComponents& components,
    double value,
    const std::vector<double>& powersW);

void writeAffineValidation(const std::vector<AffineValidationRow>& rows,
                           const std::filesystem::path& path);

} // namespace mor::parametric
