#include "transient_rom_diagnostics.hpp"

#include "linear_solvers.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#include <mkl_lapacke.h>
#endif

namespace mor::transient {
namespace {

void solveReduced(const std::vector<double>& matrix,
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
    }
    if (info != 0) throw std::runtime_error("Reduced DC solve failed.");
#else
    (void)matrix;
    (void)size;
    (void)rhs;
#endif
}

std::vector<double> reducedSteadyRhs(const TransientReducedModel& model,
                                     const std::vector<double>& powers)
{
    std::vector<double> rhs = model.reducedBoundary;
    for (int channel = 0; channel < model.sourceChannels; ++channel) {
        for (int row = 0; row < model.rank; ++row) {
            rhs[static_cast<std::size_t>(row)] +=
                powers[static_cast<std::size_t>(channel)]
                * model.reducedInput[static_cast<std::size_t>(
                    row + channel * model.rank)];
        }
    }
    return rhs;
}

} // namespace

std::vector<DcConsistencyRow> evaluateDcConsistency(
    const ThermalDescriptorSystem& descriptor,
    const TransientReducedModel& model)
{
    SubdomainDirectSolver factor(
        descriptor.dofs, sparseMatrixEntries(descriptor.conductivity));
    std::vector<std::pair<std::string, std::vector<double>>> cases;
    cases.push_back({"nominal", descriptor.nominalPowersW});
    std::vector<double> channel(static_cast<std::size_t>(descriptor.sourceChannels), 0.0);
    channel[0] = descriptor.nominalPowersW[0];
    cases.push_back({"single_channel", std::move(channel)});
    std::vector<DcConsistencyRow> rows;
    for (const auto& item : cases) {
        const std::vector<double> rhs = descriptorInputRhs(descriptor, item.second);
        std::vector<double> fom;
        factor.solve(rhs, fom);
        std::vector<double> coordinates = reducedSteadyRhs(model, item.second);
        solveReduced(model.reducedConductivity, model.rank, coordinates);
        const std::vector<double> rom = reconstructTemperature(model, coordinates);
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
        rows.push_back({item.first,
            std::sqrt(errorSquared) /
                std::max(1.0e-300, std::sqrt(referenceSquared)),
            maximumAbsolute,
            *std::max_element(fom.begin(), fom.end()),
            *std::max_element(rom.begin(), rom.end())});
    }
    return rows;
}

double initialProjectionOrthogonalityError(
    const ThermalDescriptorSystem& descriptor,
    const TransientReducedModel& model,
    const std::vector<double>& initialTemperature,
    const std::vector<double>& coordinates)
{
    const std::vector<double> projected = reconstructTemperature(model, coordinates);
    std::vector<double> error(static_cast<std::size_t>(model.globalDofs), 0.0);
    for (int row = 0; row < model.globalDofs; ++row) {
        error[static_cast<std::size_t>(row)] =
            initialTemperature[static_cast<std::size_t>(row)]
            - projected[static_cast<std::size_t>(row)];
    }
    const std::vector<double> weighted = descriptor.capacity.multiply(error);
    std::vector<double> projection(static_cast<std::size_t>(model.rank), 0.0);
#ifdef USE_MKL_PARDISO
    cblas_dgemv(CblasColMajor, CblasTrans,
        model.globalDofs, model.rank, 1.0,
        model.basis.data(), model.globalDofs,
        weighted.data(), 1, 0.0, projection.data(), 1);
#else
    for (int mode = 0; mode < model.rank; ++mode) {
        for (int row = 0; row < model.globalDofs; ++row) {
            projection[static_cast<std::size_t>(mode)] +=
                model.basis[static_cast<std::size_t>(row + mode * model.globalDofs)]
                * weighted[static_cast<std::size_t>(row)];
        }
    }
#endif
    double numerator = 0.0;
    double denominator = 0.0;
    for (double value : projection) numerator += value * value;
    for (double value : weighted) denominator += value * value;
    return std::sqrt(numerator) /
        std::max(1.0e-300, std::sqrt(denominator));
}

void writeArnoldiHistory(const TransientReducedModel& model,
                         const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "moment,input_columns,added_rank,cumulative_rank,deflated_columns,"
        << "orthogonality_error,arnoldi_residual,solve_seconds,orthogonalization_seconds,basis_bytes\n"
        << std::setprecision(17);
    for (const ArnoldiHistoryRow& row : model.arnoldiHistory) {
        out << row.moment << ',' << row.inputColumns << ',' << row.addedRank << ','
            << row.cumulativeRank << ',' << row.deflatedColumns << ','
            << row.orthogonalityError << ',' << row.arnoldiResidual << ','
            << row.solveSeconds << ',' << row.orthogonalizationSeconds << ','
            << row.basisBytes << '\n';
    }
}

void writeAccuracyByTime(const TransientAccuracySummary& accuracy,
                         const std::filesystem::path& path,
                         bool append)
{
    const bool header = !append || !std::filesystem::exists(path);
    std::ofstream out(path, append ? std::ios::app : std::ios::trunc);
    if (header) out << "waveform,step,time_s,relative_l2,maximum_absolute_k,"
        << "fom_maximum_k,rom_maximum_k,full_residual,reduced_residual,energy_balance_error\n";
    out << std::setprecision(17);
    for (const TransientAccuracyRow& row : accuracy.rows) {
        out << row.waveform << ',' << row.step << ',' << row.time << ','
            << row.relativeL2 << ',' << row.maximumAbsolute << ','
            << row.fomMaximumTemperature << ',' << row.romMaximumTemperature << ','
            << row.fullResidual << ',' << row.reducedResidual << ','
            << row.energyBalanceError << '\n';
    }
}

void writeAccuracyByWaveform(const std::vector<TransientAccuracySummary>& rows,
                             const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "waveform,space_time_relative_l2,maximum_absolute_k,max_temperature_rmse_k,"
        << "peak_temperature_error_k,peak_time_error_s,final_relative_l2,"
        << "maximum_full_residual,maximum_reduced_residual,maximum_energy_balance_error,"
        << "fom_factorization_seconds,fom_time_stepping_seconds,rom_time_stepping_seconds,"
        << "rom_reconstruction_seconds,speedup,fom_factor_memory_bytes\n"
        << std::setprecision(17);
    for (const TransientAccuracySummary& row : rows) {
        out << row.waveform << ',' << row.spaceTimeRelativeL2 << ','
            << row.maximumAbsolute << ',' << row.maximumTemperatureRmse << ','
            << row.peakTemperatureError << ',' << row.peakTimeError << ','
            << row.finalRelativeL2 << ',' << row.maximumFullResidual << ','
            << row.maximumReducedResidual << ',' << row.maximumEnergyBalanceError << ','
            << row.fomFactorizationSeconds << ',' << row.fomTimeSteppingSeconds << ','
            << row.romTimeSteppingSeconds << ',' << row.romReconstructionSeconds << ','
            << row.speedup << ',' << row.fomFactorMemoryBytes << '\n';
    }
}

void writeMaximumTemperatureCurves(const TransientAccuracySummary& accuracy,
                                   const std::filesystem::path& path,
                                   bool append)
{
    const bool header = !append || !std::filesystem::exists(path);
    std::ofstream out(path, append ? std::ios::app : std::ios::trunc);
    if (header) out << "waveform,step,time_s,fom_maximum_k,rom_maximum_k,error_k\n";
    out << std::setprecision(17);
    for (const TransientAccuracyRow& row : accuracy.rows) {
        out << row.waveform << ',' << row.step << ',' << row.time << ','
            << row.fomMaximumTemperature << ',' << row.romMaximumTemperature << ','
            << row.romMaximumTemperature - row.fomMaximumTemperature << '\n';
    }
}

void writeDcConsistency(const std::vector<DcConsistencyRow>& rows,
                        const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "case,relative_l2,maximum_absolute_k,fom_maximum_k,rom_maximum_k\n"
        << std::setprecision(17);
    for (const DcConsistencyRow& row : rows) {
        out << row.caseName << ',' << row.relativeL2 << ',' << row.maximumAbsolute << ','
            << row.fomMaximumTemperature << ',' << row.romMaximumTemperature << '\n';
    }
}

} // namespace mor::transient
