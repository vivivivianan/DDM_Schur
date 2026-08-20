#include "transient_time_integrator.hpp"

#include "linear_solvers.hpp"

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

struct DenseFactor {
    int size = 0;
    std::vector<double> values;
#ifdef USE_MKL_PARDISO
    std::vector<lapack_int> pivots;
#endif
    bool cholesky = true;
};

DenseFactor factorDense(std::vector<double> matrix, int size)
{
    DenseFactor factor;
    factor.size = size;
#ifdef USE_MKL_PARDISO
    factor.values = matrix;
    lapack_int info = LAPACKE_dpotrf(
        LAPACK_COL_MAJOR, 'L', size, factor.values.data(), size);
    if (info == 0) return factor;
    factor.cholesky = false;
    factor.values = std::move(matrix);
    factor.pivots.resize(static_cast<std::size_t>(size));
    info = LAPACKE_dsytrf(
        LAPACK_COL_MAJOR, 'L', size, factor.values.data(), size,
        factor.pivots.data());
    if (info != 0) {
        throw std::runtime_error(
            "Reduced transient LLT and LDLT factorizations both failed.");
    }
#else
    (void)matrix;
    throw std::runtime_error("Reduced transient dense factorization requires MKL.");
#endif
    return factor;
}

void solveDense(const DenseFactor& factor, std::vector<double>& rhs)
{
#ifdef USE_MKL_PARDISO
    const lapack_int info = factor.cholesky
        ? LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', factor.size, 1,
            factor.values.data(), factor.size, rhs.data(), factor.size)
        : LAPACKE_dsytrs(LAPACK_COL_MAJOR, 'L', factor.size, 1,
            factor.values.data(), factor.size, factor.pivots.data(),
            rhs.data(), factor.size);
    if (info != 0) throw std::runtime_error("Reduced transient solve failed.");
#else
    (void)factor;
    (void)rhs;
#endif
}

void denseMatvec(const std::vector<double>& matrix,
                 int size,
                 const std::vector<double>& x,
                 std::vector<double>& y,
                 double alpha = 1.0,
                 double beta = 0.0)
{
    if (y.size() != static_cast<std::size_t>(size)) y.assign(size, 0.0);
#ifdef USE_MKL_PARDISO
    cblas_dgemv(CblasColMajor, CblasNoTrans,
        size, size, alpha, matrix.data(), size,
        x.data(), 1, beta, y.data(), 1);
#else
    for (int row = 0; row < size; ++row) {
        double value = 0.0;
        for (int column = 0; column < size; ++column) {
            value += matrix[static_cast<std::size_t>(row + column * size)]
                * x[static_cast<std::size_t>(column)];
        }
        y[static_cast<std::size_t>(row)] = alpha * value
            + beta * y[static_cast<std::size_t>(row)];
    }
#endif
}

void addReducedInput(const TransientReducedModel& model,
                     const std::vector<double>& powers,
                     double scale,
                     std::vector<double>& rhs)
{
    for (int channel = 0; channel < model.sourceChannels; ++channel) {
        const double coefficient = scale * powers[static_cast<std::size_t>(channel)];
        const double* column = model.reducedInput.data()
            + static_cast<std::size_t>(channel) * model.rank;
        for (int row = 0; row < model.rank; ++row) {
            rhs[static_cast<std::size_t>(row)] += coefficient * column[row];
        }
    }
}

void addFullInput(const ThermalDescriptorSystem& descriptor,
                  const std::vector<double>& powers,
                  double scale,
                  std::vector<double>& rhs)
{
    for (int channel = 0; channel < descriptor.sourceChannels; ++channel) {
        const double coefficient = scale * powers[static_cast<std::size_t>(channel)];
        const double* column = descriptor.input.data()
            + static_cast<std::size_t>(channel) * descriptor.dofs;
        for (int row = 0; row < descriptor.dofs; ++row) {
            rhs[static_cast<std::size_t>(row)] += coefficient * column[row];
        }
    }
}

SparseMatrix timeStepMatrix(const ThermalDescriptorSystem& descriptor,
                            double timeStep,
                            const std::string& integrator)
{
    SparseMatrix matrix(descriptor.dofs);
    matrix.appendScaledEntries(descriptor.capacity, 1.0 / timeStep);
    matrix.appendScaledEntries(descriptor.conductivity,
        integrator == "crank-nicolson" ? 0.5 : 1.0);
    matrix.finalizeCsr();
    return matrix;
}

double vectorNormSquared(const std::vector<double>& vector)
{
    double value = 0.0;
    for (double item : vector) value += item * item;
    return value;
}

double discreteReducedResidual(
    const TransientReducedModel& model,
    const PowerWaveform& waveform,
    const std::vector<double>& previous,
    const std::vector<double>& current,
    double time,
    double timeStep,
    const std::string& integrator)
{
    std::vector<double> left(static_cast<std::size_t>(model.rank), 0.0);
    std::vector<double> right(static_cast<std::size_t>(model.rank), 0.0);
    const double conductivityWeight = integrator == "crank-nicolson" ? 0.5 : 1.0;
    denseMatvec(model.reducedCapacity, model.rank, current, left,
                1.0 / timeStep, 0.0);
    denseMatvec(model.reducedConductivity, model.rank, current, left,
                conductivityWeight, 1.0);
    denseMatvec(model.reducedCapacity, model.rank, previous, right,
                1.0 / timeStep, 0.0);
    const std::vector<double> powerNext = waveform.sample(time);
    if (integrator == "crank-nicolson") {
        denseMatvec(model.reducedConductivity, model.rank, previous, right,
                    -0.5, 1.0);
        addReducedInput(model, waveform.sample(time - timeStep), 0.5, right);
        addReducedInput(model, powerNext, 0.5, right);
    } else {
        addReducedInput(model, powerNext, 1.0, right);
    }
    for (int row = 0; row < model.rank; ++row) {
        right[static_cast<std::size_t>(row)] +=
            model.reducedBoundary[static_cast<std::size_t>(row)];
        left[static_cast<std::size_t>(row)] -= right[static_cast<std::size_t>(row)];
    }
    return std::sqrt(vectorNormSquared(left)) /
        std::max(1.0e-300, std::sqrt(vectorNormSquared(right)));
}

} // namespace

ReducedTrajectory integrateReducedModel(
    const TransientReducedModel& model,
    const PowerWaveform& waveform,
    double timeStep,
    int timeSteps,
    const std::string& integrator,
    const std::vector<double>& initialCoordinates,
    bool computeMaximumTemperature)
{
    const auto totalStart = Clock::now();
    if (integrator != "backward-euler" && integrator != "crank-nicolson") {
        throw std::runtime_error(
            "Transient integrator must be backward-euler or crank-nicolson.");
    }
    if (initialCoordinates.size() != static_cast<std::size_t>(model.rank)) {
        throw std::runtime_error("Reduced initial coordinate dimension mismatch.");
    }
    ReducedTrajectory result;
    result.rank = model.rank;
    result.steps = timeSteps;
    result.timeStep = timeStep;
    result.times.resize(static_cast<std::size_t>(timeSteps + 1));
    result.states.assign(static_cast<std::size_t>(timeSteps + 1)
        * static_cast<std::size_t>(model.rank), 0.0);
    result.maximumTemperature.assign(static_cast<std::size_t>(timeSteps + 1), 0.0);
    std::copy(initialCoordinates.begin(), initialCoordinates.end(), result.states.begin());

    std::vector<double> lhs = model.reducedCapacity;
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        lhs[i] /= timeStep;
        lhs[i] += (integrator == "crank-nicolson" ? 0.5 : 1.0)
            * model.reducedConductivity[i];
    }
    const auto factorStart = Clock::now();
    const DenseFactor factor = factorDense(std::move(lhs), model.rank);
    result.timing.reducedFactorizationSeconds =
        std::chrono::duration<double>(Clock::now() - factorStart).count();
    result.timing.usedLdltFallback = !factor.cholesky;

    std::vector<double> state = initialCoordinates;
    if (computeMaximumTemperature) {
        const auto reconstructionStart = Clock::now();
        const std::vector<double> temperature = reconstructTemperature(model, state);
        result.maximumTemperature[0] = *std::max_element(
            temperature.begin(), temperature.end());
        result.timing.reconstructionSeconds += std::chrono::duration<double>(
            Clock::now() - reconstructionStart).count();
    }
    const double reconstructionBeforeStepping = result.timing.reconstructionSeconds;
    const auto steppingStart = Clock::now();
    for (int step = 1; step <= timeSteps; ++step) {
        const double time = static_cast<double>(step) * timeStep;
        const std::vector<double> powerNext = waveform.sample(time);
        std::vector<double> rhs(static_cast<std::size_t>(model.rank), 0.0);
        denseMatvec(model.reducedCapacity, model.rank, state, rhs,
                    1.0 / timeStep, 0.0);
        if (integrator == "crank-nicolson") {
            denseMatvec(model.reducedConductivity, model.rank, state, rhs,
                        -0.5, 1.0);
            const std::vector<double> powerCurrent = waveform.sample(time - timeStep);
            addReducedInput(model, powerCurrent, 0.5, rhs);
            addReducedInput(model, powerNext, 0.5, rhs);
        } else {
            addReducedInput(model, powerNext, 1.0, rhs);
        }
        for (int row = 0; row < model.rank; ++row) {
            rhs[static_cast<std::size_t>(row)] +=
                model.reducedBoundary[static_cast<std::size_t>(row)];
        }
        solveDense(factor, rhs);
        state = std::move(rhs);
        result.times[static_cast<std::size_t>(step)] = time;
        std::copy(state.begin(), state.end(), result.states.begin()
            + static_cast<std::ptrdiff_t>(static_cast<std::size_t>(step) * model.rank));
        if (computeMaximumTemperature) {
            const auto reconstructionStart = Clock::now();
            const std::vector<double> temperature = reconstructTemperature(model, state);
            result.maximumTemperature[static_cast<std::size_t>(step)] =
                *std::max_element(temperature.begin(), temperature.end());
            result.timing.reconstructionSeconds += std::chrono::duration<double>(
                Clock::now() - reconstructionStart).count();
        }
    }
    result.timing.timeSteppingSeconds = std::chrono::duration<double>(
        Clock::now() - steppingStart).count()
        - (result.timing.reconstructionSeconds - reconstructionBeforeStepping);
    result.timing.totalSeconds = std::chrono::duration<double>(
        Clock::now() - totalStart).count();
    result.timing.peakWorkingSetBytes = peakWorkingSetBytes();
    return result;
}

TransientAccuracySummary compareWithFullOrder(
    const ThermalDescriptorSystem& descriptor,
    const TransientReducedModel& model,
    const PowerWaveform& waveform,
    const ReducedTrajectory& reduced,
    double timeStep,
    int timeSteps,
    const std::string& integrator,
    const std::vector<double>& initialTemperature)
{
    if (initialTemperature.size() != static_cast<std::size_t>(descriptor.dofs)) {
        throw std::runtime_error("FOM initial temperature dimension mismatch.");
    }
    TransientAccuracySummary summary;
    summary.waveform = waveform.name;
    const SparseMatrix lhs = timeStepMatrix(descriptor, timeStep, integrator);
    const auto factorStart = Clock::now();
    SubdomainDirectSolver factor(descriptor.dofs, sparseMatrixEntries(lhs));
    summary.fomFactorizationSeconds = std::chrono::duration<double>(
        Clock::now() - factorStart).count();
    summary.fomFactorMemoryBytes = factor.memoryBytes();

    std::vector<double> theta(initialTemperature.size(), 0.0);
    for (std::size_t row = 0; row < theta.size(); ++row) {
        theta[row] = initialTemperature[row] - model.referenceTemperature[row];
    }
    double errorSquaredSum = 0.0;
    double referenceSquaredSum = 0.0;
    double maximumTemperatureSquaredError = 0.0;
    double fomPeak = -std::numeric_limits<double>::infinity();
    double romPeak = -std::numeric_limits<double>::infinity();
    double fomPeakTime = 0.0;
    double romPeakTime = 0.0;
    double fomCoreSeconds = 0.0;
    const std::vector<double> kReference =
        descriptor.conductivity.multiply(model.referenceTemperature);
    std::vector<double> boundaryOffset = descriptor.boundaryRhs;
    for (int row = 0; row < descriptor.dofs; ++row) {
        boundaryOffset[static_cast<std::size_t>(row)] -=
            kReference[static_cast<std::size_t>(row)];
    }
    for (int step = 0; step <= timeSteps; ++step) {
        const double time = static_cast<double>(step) * timeStep;
        std::vector<double> rhs;
        if (step > 0) {
            const auto coreStart = Clock::now();
            rhs = descriptor.capacity.multiply(theta);
            for (double& value : rhs) value /= timeStep;
            const std::vector<double> powerNext = waveform.sample(time);
            if (integrator == "crank-nicolson") {
                const std::vector<double> conductivityImage =
                    descriptor.conductivity.multiply(theta);
                for (int row = 0; row < descriptor.dofs; ++row) {
                    rhs[static_cast<std::size_t>(row)] -=
                        0.5 * conductivityImage[static_cast<std::size_t>(row)];
                }
                const std::vector<double> powerCurrent = waveform.sample(time - timeStep);
                addFullInput(descriptor, powerCurrent, 0.5, rhs);
                addFullInput(descriptor, powerNext, 0.5, rhs);
            } else {
                addFullInput(descriptor, powerNext, 1.0, rhs);
            }
            for (int row = 0; row < descriptor.dofs; ++row) {
                rhs[static_cast<std::size_t>(row)] +=
                    boundaryOffset[static_cast<std::size_t>(row)];
            }
            std::vector<double> next;
            factor.solve(rhs, next);
            theta = std::move(next);
            fomCoreSeconds += std::chrono::duration<double>(
                Clock::now() - coreStart).count();
        }
        std::vector<double> fom = model.referenceTemperature;
        for (std::size_t row = 0; row < fom.size(); ++row) fom[row] += theta[row];
        std::vector<double> romState(static_cast<std::size_t>(model.rank), 0.0);
        std::copy_n(reduced.states.begin() + static_cast<std::ptrdiff_t>(
                        static_cast<std::size_t>(step) * model.rank),
                    static_cast<std::size_t>(model.rank), romState.begin());
        const auto romReconstructionStart = Clock::now();
        const std::vector<double> rom = reconstructTemperature(model, romState);
        summary.romReconstructionSeconds += std::chrono::duration<double>(
            Clock::now() - romReconstructionStart).count();
        double errorSquared = 0.0;
        double referenceSquared = 0.0;
        double maximumAbsolute = 0.0;
        for (int row = 0; row < descriptor.dofs; ++row) {
            const double error = rom[static_cast<std::size_t>(row)]
                - fom[static_cast<std::size_t>(row)];
            errorSquared += error * error;
            referenceSquared += fom[static_cast<std::size_t>(row)]
                * fom[static_cast<std::size_t>(row)];
            maximumAbsolute = std::max(maximumAbsolute, std::abs(error));
        }
        const double fomMaximum = *std::max_element(fom.begin(), fom.end());
        const double romMaximum = *std::max_element(rom.begin(), rom.end());
        if (fomMaximum > fomPeak) { fomPeak = fomMaximum; fomPeakTime = time; }
        if (romMaximum > romPeak) { romPeak = romMaximum; romPeakTime = time; }
        errorSquaredSum += errorSquared;
        referenceSquaredSum += referenceSquared;
        maximumTemperatureSquaredError +=
            (romMaximum - fomMaximum) * (romMaximum - fomMaximum);
        summary.maximumAbsolute = std::max(summary.maximumAbsolute, maximumAbsolute);

        double fullResidual = 0.0;
        double energyBalance = 0.0;
        if (step > 0) {
            const std::vector<double> image = lhs.multiply(theta);
            std::vector<double> residual(image.size(), 0.0);
            for (std::size_t row = 0; row < image.size(); ++row) {
                residual[row] = image[row] - rhs[row];
            }
            fullResidual = std::sqrt(vectorNormSquared(residual)) /
                std::max(1.0e-300, std::sqrt(vectorNormSquared(rhs)));
            double balance = 0.0;
            double scale = 0.0;
            for (std::size_t row = 0; row < residual.size(); ++row) {
                balance += std::abs(residual[row]);
                scale += std::abs(rhs[row]);
            }
            energyBalance = balance / std::max(1.0e-300, scale);
        }
        double reducedResidual = 0.0;
        if (step > 0) {
            std::vector<double> previousState(static_cast<std::size_t>(model.rank), 0.0);
            std::copy_n(reduced.states.begin() + static_cast<std::ptrdiff_t>(
                            static_cast<std::size_t>(step - 1) * model.rank),
                        static_cast<std::size_t>(model.rank), previousState.begin());
            reducedResidual = discreteReducedResidual(
                model, waveform, previousState, romState,
                time, timeStep, integrator);
        }
        summary.maximumFullResidual = std::max(summary.maximumFullResidual, fullResidual);
        summary.maximumReducedResidual = std::max(
            summary.maximumReducedResidual, reducedResidual);
        summary.maximumEnergyBalanceError = std::max(
            summary.maximumEnergyBalanceError, energyBalance);
        summary.rows.push_back({waveform.name, step, time,
            std::sqrt(errorSquared) /
                std::max(1.0e-300, std::sqrt(referenceSquared)),
            maximumAbsolute, fomMaximum, romMaximum,
            fullResidual, reducedResidual, energyBalance});
        if (step == timeSteps) {
            summary.finalRelativeL2 = summary.rows.back().relativeL2;
        }
    }
    summary.fomTimeSteppingSeconds = fomCoreSeconds;
    summary.spaceTimeRelativeL2 = std::sqrt(errorSquaredSum) /
        std::max(1.0e-300, std::sqrt(referenceSquaredSum));
    summary.maximumTemperatureRmse = std::sqrt(
        maximumTemperatureSquaredError / static_cast<double>(timeSteps + 1));
    summary.peakTemperatureError = std::abs(romPeak - fomPeak);
    summary.peakTimeError = std::abs(romPeakTime - fomPeakTime);
    summary.romTimeSteppingSeconds = reduced.timing.timeSteppingSeconds;
    summary.speedup = summary.fomTimeSteppingSeconds /
        std::max(1.0e-300, reduced.timing.timeSteppingSeconds);
    return summary;
}

} // namespace mor::transient
