#include "transient_reduced_model.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#include <mkl_lapacke.h>
#endif

namespace mor::transient {
namespace {

using Clock = std::chrono::steady_clock;

std::vector<double> projectSparse(const std::vector<double>& basis,
                                  int rows,
                                  int rank,
                                  const SparseMatrix& matrix)
{
    std::vector<double> reduced(static_cast<std::size_t>(rank * rank), 0.0);
    std::vector<double> column(static_cast<std::size_t>(rows), 0.0);
    std::vector<double> image;
    for (int mode = 0; mode < rank; ++mode) {
        std::copy_n(basis.begin() + static_cast<std::ptrdiff_t>(
                        static_cast<std::size_t>(mode) * rows),
                    static_cast<std::size_t>(rows), column.begin());
        image = matrix.multiply(column);
#ifdef USE_MKL_PARDISO
        cblas_dgemv(CblasColMajor, CblasTrans,
            rows, rank, 1.0, basis.data(), rows,
            image.data(), 1, 0.0,
            reduced.data() + static_cast<std::size_t>(mode) * rank, 1);
#else
        for (int test = 0; test < rank; ++test) {
            double value = 0.0;
            const double* q = basis.data() + static_cast<std::size_t>(test) * rows;
            for (int row = 0; row < rows; ++row) value += q[row] * image[row];
            reduced[static_cast<std::size_t>(test + mode * rank)] = value;
        }
#endif
    }
    return reduced;
}

std::vector<double> projectDenseInput(const std::vector<double>& basis,
                                      int rows,
                                      int rank,
                                      const std::vector<double>& input,
                                      int channels)
{
    std::vector<double> reduced(
        static_cast<std::size_t>(rank * channels), 0.0);
#ifdef USE_MKL_PARDISO
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
        rank, channels, rows, 1.0,
        basis.data(), rows, input.data(), rows,
        0.0, reduced.data(), rank);
#else
    for (int channel = 0; channel < channels; ++channel) {
        for (int mode = 0; mode < rank; ++mode) {
            double value = 0.0;
            for (int row = 0; row < rows; ++row) {
                value += basis[static_cast<std::size_t>(row + mode * rows)]
                    * input[static_cast<std::size_t>(row + channel * rows)];
            }
            reduced[static_cast<std::size_t>(mode + channel * rank)] = value;
        }
    }
#endif
    return reduced;
}

void symmetrize(std::vector<double>& matrix, int size)
{
    for (int column = 0; column < size; ++column) {
        for (int row = column + 1; row < size; ++row) {
            const std::size_t lower = static_cast<std::size_t>(row + column * size);
            const std::size_t upper = static_cast<std::size_t>(column + row * size);
            const double average = 0.5 * (matrix[lower] + matrix[upper]);
            matrix[lower] = average;
            matrix[upper] = average;
        }
    }
}

ReducedMatrixDiagnostic diagnoseReduced(const std::vector<double>& matrix,
                                        int size)
{
    ReducedMatrixDiagnostic result;
    double skew = 0.0;
    double norm = 0.0;
    for (int column = 0; column < size; ++column) {
        for (int row = 0; row < size; ++row) {
            const double value = matrix[static_cast<std::size_t>(row + column * size)];
            const double transpose = matrix[static_cast<std::size_t>(column + row * size)];
            skew += (value - transpose) * (value - transpose);
            norm += value * value;
        }
    }
    result.symmetryError = std::sqrt(skew) /
        std::max(1.0e-300, std::sqrt(norm));
#ifdef USE_MKL_PARDISO
    std::vector<double> factor = matrix;
    result.choleskySucceeded = LAPACKE_dpotrf(
        LAPACK_COL_MAJOR, 'L', size, factor.data(), size) == 0;
    std::vector<double> eigenMatrix = matrix;
    std::vector<double> eigenvalues(static_cast<std::size_t>(size), 0.0);
    const lapack_int info = LAPACKE_dsyev(
        LAPACK_COL_MAJOR, 'N', 'U', size,
        eigenMatrix.data(), size, eigenvalues.data());
    if (info != 0) {
        throw std::runtime_error("Reduced SPD diagnostic eigensolve failed.");
    }
    result.minimumEigenvalue = eigenvalues.front();
    result.maximumEigenvalue = eigenvalues.back();
#else
    result.choleskySucceeded = true;
    result.minimumEigenvalue = std::numeric_limits<double>::quiet_NaN();
    result.maximumEigenvalue = std::numeric_limits<double>::quiet_NaN();
#endif
    return result;
}

std::vector<double> projectVector(const std::vector<double>& basis,
                                  int rows,
                                  int rank,
                                  const std::vector<double>& vector)
{
    std::vector<double> result(static_cast<std::size_t>(rank), 0.0);
#ifdef USE_MKL_PARDISO
    cblas_dgemv(CblasColMajor, CblasTrans,
        rows, rank, 1.0, basis.data(), rows,
        vector.data(), 1, 0.0, result.data(), 1);
#else
    for (int mode = 0; mode < rank; ++mode) {
        for (int row = 0; row < rows; ++row) {
            result[static_cast<std::size_t>(mode)] +=
                basis[static_cast<std::size_t>(row + mode * rows)]
                * vector[static_cast<std::size_t>(row)];
        }
    }
#endif
    return result;
}

void solveSpd(const std::vector<double>& matrix,
              int size,
              std::vector<double>& rhs)
{
#ifdef USE_MKL_PARDISO
    std::vector<double> factor = matrix;
    lapack_int info = LAPACKE_dpotrf(
        LAPACK_COL_MAJOR, 'L', size, factor.data(), size);
    if (info == 0) {
        info = LAPACKE_dpotrs(
            LAPACK_COL_MAJOR, 'L', size, 1,
            factor.data(), size, rhs.data(), size);
    } else {
        std::vector<lapack_int> pivots(static_cast<std::size_t>(size));
        factor = matrix;
        info = LAPACKE_dsytrf(
            LAPACK_COL_MAJOR, 'L', size, factor.data(), size, pivots.data());
        if (info == 0) {
            info = LAPACKE_dsytrs(
                LAPACK_COL_MAJOR, 'L', size, 1,
                factor.data(), size, pivots.data(), rhs.data(), size);
        }
    }
    if (info != 0) {
        throw std::runtime_error(
            "Reduced initial-condition projection factorization failed.");
    }
#else
    (void)matrix;
    (void)size;
    (void)rhs;
    throw std::runtime_error("Reduced dense solve requires MKL in this build.");
#endif
}

} // namespace

TransientReducedModel buildTransientReducedModel(
    const ThermalDescriptorSystem& descriptor,
    BlockArnoldiResult arnoldi)
{
    const auto start = Clock::now();
    if (arnoldi.rank <= 0) {
        throw std::runtime_error("Block Arnoldi produced an empty basis.");
    }
    TransientReducedModel model;
    model.globalDofs = descriptor.dofs;
    model.rank = arnoldi.rank;
    model.blockSize = arnoldi.blockSize;
    model.moments = arnoldi.moments;
    model.sourceChannels = descriptor.sourceChannels;
    model.expansionPoint = arnoldi.expansionPoint;
    model.rankTolerance = arnoldi.rankTolerance;
    model.massType = descriptor.massType;
    model.basis = std::move(arnoldi.basis);
    model.referenceTemperature = std::move(arnoldi.referenceTemperature);
    model.arnoldiHistory = std::move(arnoldi.history);
    model.arnoldiTiming = arnoldi.timing;
    model.fingerprints = descriptor.fingerprints;
    model.nominalPowersW = descriptor.nominalPowersW;
    model.minimumPowersW = descriptor.minimumPowersW;
    model.maximumPowersW = descriptor.maximumPowersW;
    model.sourceSubdomains = descriptor.sourceSubdomains;
    model.sourceDomainEntities = descriptor.sourceDomainEntities;
    model.deploymentDofs = descriptor.deploymentDofs;

    model.reducedCapacity = projectSparse(
        model.basis, model.globalDofs, model.rank, descriptor.capacity);
    model.reducedConductivity = projectSparse(
        model.basis, model.globalDofs, model.rank, descriptor.conductivity);
    model.reducedInput = projectDenseInput(
        model.basis, model.globalDofs, model.rank,
        descriptor.input, descriptor.sourceChannels);
    symmetrize(model.reducedCapacity, model.rank);
    symmetrize(model.reducedConductivity, model.rank);

    const std::vector<double> referenceImage =
        descriptor.conductivity.multiply(model.referenceTemperature);
    std::vector<double> referenceResidual(
        static_cast<std::size_t>(model.globalDofs), 0.0);
    double residualSquared = 0.0;
    double rhsSquared = 0.0;
    for (int row = 0; row < model.globalDofs; ++row) {
        const double residual = descriptor.boundaryRhs[static_cast<std::size_t>(row)]
            - referenceImage[static_cast<std::size_t>(row)];
        referenceResidual[static_cast<std::size_t>(row)] = residual;
        residualSquared += residual * residual;
        rhsSquared += descriptor.boundaryRhs[static_cast<std::size_t>(row)]
            * descriptor.boundaryRhs[static_cast<std::size_t>(row)];
    }
    model.referenceResidual = std::sqrt(residualSquared)
        / std::max(1.0e-300, std::sqrt(rhsSquared));
    model.reducedBoundary = projectVector(
        model.basis, model.globalDofs, model.rank, referenceResidual);
    model.capacityDiagnostic = diagnoseReduced(model.reducedCapacity, model.rank);
    model.conductivityDiagnostic = diagnoseReduced(
        model.reducedConductivity, model.rank);
    if (!model.capacityDiagnostic.choleskySucceeded
        || !model.conductivityDiagnostic.choleskySucceeded
        || !(model.capacityDiagnostic.minimumEigenvalue > 0.0)
        || !(model.conductivityDiagnostic.minimumEigenvalue > 0.0)) {
        throw std::runtime_error(
            "Congruence projection produced a non-SPD reduced C or K matrix.");
    }
    for (const ArnoldiHistoryRow& row : model.arnoldiHistory) {
        model.basisOrthogonalityError = std::max(
            model.basisOrthogonalityError, row.orthogonalityError);
    }
    model.projectionSeconds = std::chrono::duration<double>(
        Clock::now() - start).count();
    return model;
}

std::vector<double> projectInitialConditionCWeighted(
    const ThermalDescriptorSystem& descriptor,
    const TransientReducedModel& model,
    const std::vector<double>& initialTemperature)
{
    if (initialTemperature.size() != static_cast<std::size_t>(model.globalDofs)) {
        throw std::runtime_error("Transient initial temperature dimension mismatch.");
    }
    std::vector<double> theta(initialTemperature.size(), 0.0);
    for (std::size_t row = 0; row < theta.size(); ++row) {
        theta[row] = initialTemperature[row] - model.referenceTemperature[row];
    }
    const std::vector<double> weighted = descriptor.capacity.multiply(theta);
    std::vector<double> coordinates = projectVector(
        model.basis, model.globalDofs, model.rank, weighted);
    solveSpd(model.reducedCapacity, model.rank, coordinates);
    return coordinates;
}

std::vector<double> reconstructTemperature(
    const TransientReducedModel& model,
    const std::vector<double>& coordinates)
{
    if (coordinates.size() != static_cast<std::size_t>(model.rank)) {
        throw std::runtime_error("Transient reduced coordinate dimension mismatch.");
    }
    std::vector<double> temperature = model.referenceTemperature;
#ifdef USE_MKL_PARDISO
    cblas_dgemv(CblasColMajor, CblasNoTrans,
        model.globalDofs, model.rank, 1.0,
        model.basis.data(), model.globalDofs,
        coordinates.data(), 1, 1.0, temperature.data(), 1);
#else
    for (int mode = 0; mode < model.rank; ++mode) {
        for (int row = 0; row < model.globalDofs; ++row) {
            temperature[static_cast<std::size_t>(row)] +=
                model.basis[static_cast<std::size_t>(row + mode * model.globalDofs)]
                * coordinates[static_cast<std::size_t>(mode)];
        }
    }
#endif
    return temperature;
}

} // namespace mor::transient
