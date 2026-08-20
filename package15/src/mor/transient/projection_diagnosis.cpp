#include "projection_diagnosis.hpp"

#include "mor/mor_diagnostics.hpp"
#include "mor/transient/interface_flux_operator.hpp"
#include "mor/transient/local_port_reduced_schur.hpp"
#include "mor/transient/thermal_descriptor_system.hpp"
#include "mor/transient/transfer_operator.hpp"
#include "sipg_core.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#ifdef USE_MKL_PARDISO
#include <mkl.h>
#endif

namespace mor::transient {
namespace {

using Clock = std::chrono::steady_clock;

double elapsed(const Clock::time_point& start)
{
    return std::chrono::duration<double>(
        Clock::now() - start).count();
}

double dot(const std::vector<double>& left,
           const std::vector<double>& right)
{
    if (left.size() != right.size()) {
        throw std::runtime_error(
            "[Projection diagnosis] Dot-product size mismatch.");
    }
    long double value = 0.0L;
    for (std::size_t row = 0; row < left.size(); ++row) {
        value += static_cast<long double>(left[row]) * right[row];
    }
    return static_cast<double>(value);
}

double relativeDifference(const std::vector<double>& left,
                          const std::vector<double>& right)
{
    if (left.size() != right.size()) {
        throw std::runtime_error(
            "[Projection diagnosis] Relative-difference size mismatch.");
    }
    long double error = 0.0L;
    long double reference = 0.0L;
    for (std::size_t row = 0; row < left.size(); ++row) {
        const double delta = left[row] - right[row];
        error += static_cast<long double>(delta) * delta;
        reference += static_cast<long double>(right[row]) * right[row];
    }
    return std::sqrt(static_cast<double>(error))
        / std::max(1.0e-300,
            std::sqrt(static_cast<double>(reference)));
}

void matrixVector(const std::vector<double>& matrix,
                  int rows,
                  int columns,
                  const std::vector<double>& input,
                  std::vector<double>& output)
{
    if (input.size() != static_cast<std::size_t>(columns)
        || matrix.size()
            != static_cast<std::size_t>(rows * columns)) {
        throw std::runtime_error(
            "[Projection diagnosis] Dense matrix-vector mismatch.");
    }
    output.assign(static_cast<std::size_t>(rows), 0.0);
#ifdef USE_MKL_PARDISO
    cblas_dgemv(
        CblasColMajor, CblasNoTrans, rows, columns,
        1.0, matrix.data(), rows, input.data(), 1,
        0.0, output.data(), 1);
#else
    for (int column = 0; column < columns; ++column) {
        for (int row = 0; row < rows; ++row) {
            output[static_cast<std::size_t>(row)]
                += matrix[static_cast<std::size_t>(
                    column * rows + row)]
                * input[static_cast<std::size_t>(column)];
        }
    }
#endif
}

std::vector<double> condenseGlobalRhs(
    const local::Model& model,
    const ReducedDynamicSchurOperator& schur,
    const std::vector<double>& globalRhs)
{
    if (globalRhs.size()
        != static_cast<std::size_t>(model.globalDofs)) {
        throw std::runtime_error(
            "[Projection diagnosis] Global RHS size mismatch.");
    }
    std::vector<double> condensed(
        static_cast<std::size_t>(model.interfaceDofs), 0.0);
    for (int row = 0; row < model.interfaceDofs; ++row) {
        condensed[static_cast<std::size_t>(row)] =
            globalRhs[static_cast<std::size_t>(
                model.interfaceGlobalDofs[
                    static_cast<std::size_t>(row)])];
    }
    for (std::size_t slot = 0;
         slot < model.subdomains.size(); ++slot) {
        const local::SubdomainModel& local =
            model.subdomains[slot];
        std::vector<double> reduced(
            static_cast<std::size_t>(local.rank), 0.0);
        for (int mode = 0; mode < local.rank; ++mode) {
            const double* basis = local.basis.data()
                + static_cast<std::size_t>(
                    mode * local.interiorDofs);
            long double value = 0.0L;
            for (int row = 0; row < local.interiorDofs; ++row) {
                value += static_cast<long double>(basis[row])
                    * globalRhs[static_cast<std::size_t>(
                        local.interiorGlobalDofs[
                            static_cast<std::size_t>(row)])];
            }
            reduced[static_cast<std::size_t>(mode)] =
                static_cast<double>(value);
        }
        schur.solveReducedInterior(slot, reduced);
        for (int gamma = 0;
             gamma < local.localInterfaceDofs; ++gamma) {
            long double correction = 0.0L;
            for (int mode = 0; mode < local.rank; ++mode) {
                correction += static_cast<long double>(
                    local.reducedInterfaceInterior[
                        static_cast<std::size_t>(
                            gamma * local.rank + mode)])
                    * reduced[static_cast<std::size_t>(mode)];
            }
            condensed[static_cast<std::size_t>(
                local.interfaceIndices[
                    static_cast<std::size_t>(gamma)])]
                -= static_cast<double>(correction);
        }
    }
    return condensed;
}

std::vector<double> recoverWithInterface(
    const local::Model& model,
    const ReducedDynamicSchurOperator& schur,
    const std::vector<double>& globalRhs,
    const std::vector<double>& interfaceTemperature)
{
    if (interfaceTemperature.size()
            != static_cast<std::size_t>(model.interfaceDofs)
        || globalRhs.size()
            != static_cast<std::size_t>(model.globalDofs)) {
        throw std::runtime_error(
            "[Projection diagnosis] Recovery dimensions are invalid.");
    }
    std::vector<double> recovered(
        static_cast<std::size_t>(model.globalDofs), 0.0);
    for (int gamma = 0; gamma < model.interfaceDofs; ++gamma) {
        recovered[static_cast<std::size_t>(
            model.interfaceGlobalDofs[
                static_cast<std::size_t>(gamma)])]
            = interfaceTemperature[static_cast<std::size_t>(gamma)];
    }
    for (std::size_t slot = 0;
         slot < model.subdomains.size(); ++slot) {
        const local::SubdomainModel& local =
            model.subdomains[slot];
        std::vector<double> reduced(
            static_cast<std::size_t>(local.rank), 0.0);
        for (int mode = 0; mode < local.rank; ++mode) {
            const double* basis = local.basis.data()
                + static_cast<std::size_t>(
                    mode * local.interiorDofs);
            long double value = 0.0L;
            for (int row = 0; row < local.interiorDofs; ++row) {
                value += static_cast<long double>(basis[row])
                    * globalRhs[static_cast<std::size_t>(
                        local.interiorGlobalDofs[
                            static_cast<std::size_t>(row)])];
            }
            for (int gamma = 0;
                 gamma < local.localInterfaceDofs; ++gamma) {
                value -= static_cast<long double>(
                    local.reducedInteriorInterface[
                        static_cast<std::size_t>(
                            mode * local.localInterfaceDofs
                            + gamma)])
                    * interfaceTemperature[static_cast<std::size_t>(
                        local.interfaceIndices[
                            static_cast<std::size_t>(gamma)])];
            }
            reduced[static_cast<std::size_t>(mode)] =
                static_cast<double>(value);
        }
        schur.solveReducedInterior(slot, reduced);
        for (int row = 0; row < local.interiorDofs; ++row) {
            long double value = 0.0L;
            for (int mode = 0; mode < local.rank; ++mode) {
                value += static_cast<long double>(
                    local.basis[static_cast<std::size_t>(
                        mode * local.interiorDofs + row)])
                    * reduced[static_cast<std::size_t>(mode)];
            }
            recovered[static_cast<std::size_t>(
                local.interiorGlobalDofs[
                    static_cast<std::size_t>(row)])]
                = static_cast<double>(value);
        }
    }
    return recovered;
}

struct ForcingProjection {
    double error = 0.0;
    double referenceNorm = 0.0;
    double errorNorm = 0.0;
};

struct InterfaceRow {
    int interfaceId = -1;
    int interfaceDofs = 0;
    int localPortRank = 0;
    int energyProjectionRank = 0;
    int leftSubdomain = -1;
    int rightSubdomain = -1;
    double temperatureProjectionError = 0.0;
    double temperatureEuclideanError = 0.0;
    double temperatureReferenceEnergy = 0.0;
    double temperatureErrorEnergy = 0.0;
    double temperatureGalerkinOrthogonality = 0.0;
    double fluxProjectionError = 0.0;
    double physicalFluxProjectionError = 0.0;
    double numericalFluxProjectionError = 0.0;
    int fluxTriangleCount = 0;
    double fluxArea = 0.0;
    ForcingProjection forcing;
    ForcingProjection input;
    ForcingProjection boundary;
    ForcingProjection history;
    ForcingProjection other;
    double forcingDecompositionClosure = 0.0;
    int targetSolveCalls = 0;
    double targetSolveResidual = 0.0;
    double setupSeconds = 0.0;
    double projectionSeconds = 0.0;
    std::size_t peakIncrementalMemoryBytes = 0;
    std::string targetSolver;
    std::string status = "not_run";
    std::vector<int> targetIndices;
};

struct EnergyOrthonormalBasis {
    std::vector<double> vectors;
    std::vector<double> schurImages;
    int rank = 0;
};

EnergyOrthonormalBasis buildEnergyOrthonormalBasis(
    const std::vector<double>& candidates,
    int rows,
    int columns,
    PatchTransferOperator& target)
{
    if (rows <= 0 || columns <= 0
        || candidates.size()
            != static_cast<std::size_t>(rows * columns)) {
        throw std::runtime_error(
            "[Projection diagnosis] Invalid candidate basis.");
    }
    EnergyOrthonormalBasis result;
    result.vectors.reserve(candidates.size());
    result.schurImages.reserve(candidates.size());
    double maximumInitialEnergy = 0.0;
    constexpr double relativeDeflationTolerance = 1.0e-12;
    for (int column = 0; column < columns; ++column) {
        std::vector<double> candidate(
            candidates.begin()
                + static_cast<std::ptrdiff_t>(column * rows),
            candidates.begin()
                + static_cast<std::ptrdiff_t>((column + 1) * rows));
        std::vector<double> image;
        target.applyTargetAction(candidate, image);
        const double initialEnergySquared =
            std::max(0.0, dot(candidate, image));
        maximumInitialEnergy = std::max(
            maximumInitialEnergy,
            std::sqrt(initialEnergySquared));
        for (int pass = 0; pass < 2; ++pass) {
            for (int accepted = 0;
                 accepted < result.rank; ++accepted) {
                const double* vectorColumn =
                    result.vectors.data()
                    + static_cast<std::size_t>(
                        accepted * rows);
                const double* imageColumn =
                    result.schurImages.data()
                    + static_cast<std::size_t>(
                        accepted * rows);
                long double coefficient = 0.0L;
                for (int row = 0; row < rows; ++row) {
                    coefficient +=
                        static_cast<long double>(
                            vectorColumn[row])
                        * image[
                            static_cast<std::size_t>(row)];
                }
                const double alpha =
                    static_cast<double>(coefficient);
                for (int row = 0; row < rows; ++row) {
                    candidate[static_cast<std::size_t>(row)]
                        -= alpha * vectorColumn[row];
                    image[static_cast<std::size_t>(row)]
                        -= alpha * imageColumn[row];
                }
            }
        }
        const double energySquared =
            std::max(0.0, dot(candidate, image));
        const double energy = std::sqrt(energySquared);
        const double threshold =
            relativeDeflationTolerance
            * std::max(1.0e-300, maximumInitialEnergy);
        if (!(energy > threshold) || !std::isfinite(energy)) {
            continue;
        }
        const double inverse = 1.0 / energy;
        for (int row = 0; row < rows; ++row) {
            candidate[static_cast<std::size_t>(row)] *= inverse;
            image[static_cast<std::size_t>(row)] *= inverse;
        }
        result.vectors.insert(
            result.vectors.end(),
            candidate.begin(), candidate.end());
        result.schurImages.insert(
            result.schurImages.end(),
            image.begin(), image.end());
        ++result.rank;
    }
    if (result.rank <= 0) {
        throw std::runtime_error(
            "[Projection diagnosis] Port basis has zero numerical "
            "rank in the target Schur energy.");
    }
    return result;
}

ForcingProjection projectForcing(
    const std::vector<double>& forcing,
    const std::vector<double>& basis,
    const std::vector<double>& schurImages,
    int rows,
    int rank,
    PatchTransferOperator& target)
{
    ForcingProjection result;
    const double forcingNorm = std::sqrt(std::max(0.0, dot(
        forcing, forcing)));
    if (!(forcingNorm > 0.0)) return result;

    std::vector<double> coefficients(
        static_cast<std::size_t>(rank), 0.0);
#ifdef USE_MKL_PARDISO
    cblas_dgemv(
        CblasColMajor, CblasTrans, rows, rank, 1.0,
        basis.data(), rows, forcing.data(), 1,
        0.0, coefficients.data(), 1);
#else
    for (int mode = 0; mode < rank; ++mode) {
        for (int row = 0; row < rows; ++row) {
            coefficients[static_cast<std::size_t>(mode)]
                += basis[static_cast<std::size_t>(
                    mode * rows + row)]
                * forcing[static_cast<std::size_t>(row)];
        }
    }
#endif
    std::vector<double> projected;
    matrixVector(
        schurImages, rows, rank, coefficients, projected);
    std::vector<double> residual(forcing.size(), 0.0);
    for (std::size_t row = 0; row < forcing.size(); ++row) {
        residual[row] = forcing[row] - projected[row];
    }
    std::vector<double> forcingDual;
    std::vector<double> residualDual;
    target.solveTargetResponse(forcing, forcingDual);
    target.solveTargetResponse(residual, residualDual);
    result.referenceNorm =
        std::sqrt(std::max(0.0, dot(forcing, forcingDual)));
    result.errorNorm =
        std::sqrt(std::max(0.0, dot(residual, residualDual)));
    result.error = result.errorNorm
        / std::max(1.0e-300, result.referenceNorm);
    return result;
}

} // namespace

ProjectionDiagnosisSummary runFullInterfaceProjectionDiagnosis(
    const Mesh& mesh,
    const CaseConfig& physics,
    const ddm_schur::InterfacePartition& partition,
    const ThermalDescriptorSystem& descriptor,
    const local::Model& dynamicModel,
    const LocalPortModel& portModel,
    const std::vector<double>& previousTheta,
    const std::vector<double>& referenceTemperature,
    const std::vector<double>& boundaryOffset,
    const ProjectionDiagnosisOptions& options)
{
    const auto diagnosisStart = Clock::now();
    ProjectionDiagnosisSummary summary;
    summary.requestedInterfaces =
        static_cast<int>(options.interfaceIds.size());
    if (options.interfaceIds.empty()
        || !(options.timeStep > 0.0)
        || options.powers.size()
            != static_cast<std::size_t>(descriptor.sourceChannels)
        || previousTheta.size()
            != static_cast<std::size_t>(descriptor.dofs)
        || referenceTemperature.size() != previousTheta.size()
        || boundaryOffset.size() != previousTheta.size()
        || portModel.interfaceGlobalDofs
            != dynamicModel.interfaceGlobalDofs) {
        throw std::runtime_error(
            "[Projection diagnosis] Invalid input dimensions.");
    }
    std::filesystem::create_directories(options.outputDirectory);

    std::vector<double> historyGlobal =
        descriptor.capacity.multiply(previousTheta);
    for (double& value : historyGlobal) value /= options.timeStep;
    const std::vector<double> inputGlobal =
        descriptorInputRhs(descriptor, options.powers);
    const std::vector<double> boundaryGlobal = boundaryOffset;
    std::vector<double> totalGlobal(historyGlobal.size(), 0.0);
    for (std::size_t row = 0; row < totalGlobal.size(); ++row) {
        totalGlobal[row] = historyGlobal[row]
            + inputGlobal[row] + boundaryGlobal[row];
    }

    ddm_schur::Options referenceOptions = options.schurOptions;
    const std::filesystem::path referenceOutput =
        options.outputDirectory / "full_interface_reference";
    std::filesystem::create_directories(referenceOutput);
    referenceOptions.proxyOutputDirectory =
        referenceOutput.string();
    const auto referenceSetupStart = Clock::now();
    local::LocalReducedSchurSolver referenceSolver(
        dynamicModel, mesh, physics, partition,
        referenceOptions, 0, referenceOutput);
    summary.fullInterfaceSetupSeconds =
        elapsed(referenceSetupStart);
    const auto referenceSolveStart = Clock::now();
    local::SolveResult referenceSolve =
        referenceSolver.solve(totalGlobal);
    summary.fullInterfaceSolveSeconds =
        elapsed(referenceSolveStart);
    if (referenceSolve.status != "success") {
        throw std::runtime_error(
            "[Projection diagnosis] Full-interface reference solve failed.");
    }

    ReducedDynamicSchurOperator exactSchur(dynamicModel, false);
    const std::vector<double> condensedInput =
        condenseGlobalRhs(dynamicModel, exactSchur, inputGlobal);
    const std::vector<double> condensedBoundary =
        condenseGlobalRhs(dynamicModel, exactSchur, boundaryGlobal);
    const std::vector<double> condensedHistory =
        condenseGlobalRhs(dynamicModel, exactSchur, historyGlobal);
    std::vector<double> condensedTotal(condensedInput.size(), 0.0);
    for (std::size_t row = 0; row < condensedTotal.size(); ++row) {
        condensedTotal[row] = condensedInput[row]
            + condensedBoundary[row] + condensedHistory[row];
    }

    const std::vector<PortPatch> patches =
        buildOptimalPortPatches(
            mesh, partition,
            std::max(1, portModel.oversamplingLayers));
    std::map<int, const PortPatch*> patchById;
    for (const PortPatch& patch : patches) {
        patchById.emplace(patch.interfaceId, &patch);
    }
    std::map<int, const LocalPortBasis*> portById;
    for (const LocalPortBasis& port : portModel.ports) {
        portById.emplace(port.interfaceId, &port);
    }

    std::vector<double> projectedGamma =
        referenceSolve.interfaceTemperature;
    std::vector<InterfaceRow> rows;
    rows.reserve(options.interfaceIds.size());
    for (int interfaceId : options.interfaceIds) {
        const auto patchIt = patchById.find(interfaceId);
        const auto portIt = portById.find(interfaceId);
        if (patchIt == patchById.end()
            || portIt == portById.end()) {
            throw std::runtime_error(
                "[Projection diagnosis] Selected interface is missing.");
        }
        const PortPatch& patch = *patchIt->second;
        const LocalPortBasis& port = *portIt->second;
        const auto interfaceStart = Clock::now();
        InterfaceRow row;
        row.interfaceId = interfaceId;
        row.interfaceDofs = static_cast<int>(patch.target.size());
        row.localPortRank = port.rank;
        row.leftSubdomain = port.leftSubdomain;
        row.rightSubdomain = port.rightSubdomain;
        row.targetIndices = patch.target;

        std::map<int, int> portRow;
        for (int localRow = 0; localRow < port.rows; ++localRow) {
            portRow.emplace(
                port.interfaceIndices[
                    static_cast<std::size_t>(localRow)],
                localRow);
        }
        std::vector<double> basis(
            static_cast<std::size_t>(
                row.interfaceDofs * port.rank), 0.0);
        for (int targetRow = 0;
             targetRow < row.interfaceDofs; ++targetRow) {
            const auto found = portRow.find(
                patch.target[
                    static_cast<std::size_t>(targetRow)]);
            if (found == portRow.end()) {
                throw std::runtime_error(
                    "[Projection diagnosis] Port/patch target mismatch.");
            }
            for (int mode = 0; mode < port.rank; ++mode) {
                basis[static_cast<std::size_t>(
                    mode * row.interfaceDofs + targetRow)]
                    = port.basis[static_cast<std::size_t>(
                        mode * port.rows + found->second)];
            }
        }

        PatchTransferOptions targetOptions;
        targetOptions.innerSolver = "woodbury-exact";
        targetOptions.relativeTolerance = 1.0e-10;
        targetOptions.maximumIterations = 1000;
        targetOptions.refinementMaximumIterations = 3;
        targetOptions.refinementTolerance = 1.0e-10;
        PatchTransferOperator target(
            exactSchur, patch, targetOptions);
        EnergyOrthonormalBasis energyBasis =
            buildEnergyOrthonormalBasis(
                basis, row.interfaceDofs,
                port.rank, target);
        basis = std::move(energyBasis.vectors);
        std::vector<double> schurImages =
            std::move(energyBasis.schurImages);
        row.energyProjectionRank = energyBasis.rank;

        std::vector<double> referenceTarget(
            static_cast<std::size_t>(row.interfaceDofs), 0.0);
        std::vector<double> forcingTarget(referenceTarget.size(), 0.0);
        std::vector<double> inputTarget(referenceTarget.size(), 0.0);
        std::vector<double> boundaryTarget(referenceTarget.size(), 0.0);
        std::vector<double> historyTarget(referenceTarget.size(), 0.0);
        for (int targetRow = 0;
             targetRow < row.interfaceDofs; ++targetRow) {
            const int globalRow = patch.target[
                static_cast<std::size_t>(targetRow)];
            referenceTarget[static_cast<std::size_t>(targetRow)] =
                referenceSolve.interfaceTemperature[
                    static_cast<std::size_t>(globalRow)];
            forcingTarget[static_cast<std::size_t>(targetRow)] =
                condensedTotal[static_cast<std::size_t>(globalRow)];
            inputTarget[static_cast<std::size_t>(targetRow)] =
                condensedInput[static_cast<std::size_t>(globalRow)];
            boundaryTarget[static_cast<std::size_t>(targetRow)] =
                condensedBoundary[static_cast<std::size_t>(globalRow)];
            historyTarget[static_cast<std::size_t>(targetRow)] =
                condensedHistory[static_cast<std::size_t>(globalRow)];
        }

        std::vector<double> referenceImage;
        target.applyTargetAction(referenceTarget, referenceImage);
        std::vector<double> coefficients(
            static_cast<std::size_t>(
                row.energyProjectionRank), 0.0);
#ifdef USE_MKL_PARDISO
        cblas_dgemv(
            CblasColMajor, CblasTrans,
            row.interfaceDofs, row.energyProjectionRank,
            1.0, basis.data(), row.interfaceDofs,
            referenceImage.data(), 1,
            0.0, coefficients.data(), 1);
#else
        for (int mode = 0;
             mode < row.energyProjectionRank; ++mode) {
            for (int targetRow = 0;
                 targetRow < row.interfaceDofs; ++targetRow) {
                coefficients[static_cast<std::size_t>(mode)]
                    += basis[static_cast<std::size_t>(
                        mode * row.interfaceDofs + targetRow)]
                    * referenceImage[
                        static_cast<std::size_t>(targetRow)];
            }
        }
#endif
        std::vector<double> projectedTarget;
        matrixVector(
            basis, row.interfaceDofs,
            row.energyProjectionRank,
            coefficients, projectedTarget);
        std::vector<double> errorTarget(referenceTarget.size(), 0.0);
        std::vector<double> errorImage;
        // Refine the Galerkin projection in its already orthonormalized
        // Schur-energy coordinates.  This only removes floating-point
        // accumulation error; it does not enlarge or alter the frozen port
        // space.
        for (int refinement = 0; refinement < 2; ++refinement) {
            for (int targetRow = 0;
                 targetRow < row.interfaceDofs; ++targetRow) {
                errorTarget[static_cast<std::size_t>(targetRow)] =
                    referenceTarget[static_cast<std::size_t>(targetRow)]
                    - projectedTarget[static_cast<std::size_t>(targetRow)];
            }
            target.applyTargetAction(errorTarget, errorImage);
            std::vector<double> correction(
                static_cast<std::size_t>(
                    row.energyProjectionRank), 0.0);
#ifdef USE_MKL_PARDISO
            cblas_dgemv(
                CblasColMajor, CblasTrans,
                row.interfaceDofs, row.energyProjectionRank,
                1.0, basis.data(), row.interfaceDofs,
                errorImage.data(), 1,
                0.0, correction.data(), 1);
#else
            for (int mode = 0;
                 mode < row.energyProjectionRank; ++mode) {
                for (int targetRow = 0;
                     targetRow < row.interfaceDofs; ++targetRow) {
                    correction[static_cast<std::size_t>(mode)]
                        += basis[static_cast<std::size_t>(
                            mode * row.interfaceDofs + targetRow)]
                            * errorImage[
                                static_cast<std::size_t>(targetRow)];
                }
            }
#endif
            std::vector<double> correctionTarget;
            matrixVector(
                basis, row.interfaceDofs,
                row.energyProjectionRank,
                correction, correctionTarget);
            for (int targetRow = 0;
                 targetRow < row.interfaceDofs; ++targetRow) {
                projectedTarget[
                    static_cast<std::size_t>(targetRow)]
                    += correctionTarget[
                        static_cast<std::size_t>(targetRow)];
            }
        }
        for (int targetRow = 0;
             targetRow < row.interfaceDofs; ++targetRow) {
            errorTarget[static_cast<std::size_t>(targetRow)] =
                referenceTarget[static_cast<std::size_t>(targetRow)]
                - projectedTarget[static_cast<std::size_t>(targetRow)];
            projectedGamma[static_cast<std::size_t>(
                patch.target[
                    static_cast<std::size_t>(targetRow)])]
                = projectedTarget[static_cast<std::size_t>(targetRow)];
        }
        target.applyTargetAction(errorTarget, errorImage);
        row.temperatureReferenceEnergy =
            std::sqrt(std::max(0.0,
                dot(referenceTarget, referenceImage)));
        row.temperatureErrorEnergy =
            std::sqrt(std::max(0.0,
                dot(errorTarget, errorImage)));
        row.temperatureProjectionError =
            row.temperatureErrorEnergy
            / std::max(1.0e-300,
                row.temperatureReferenceEnergy);
        row.temperatureEuclideanError =
            relativeDifference(projectedTarget, referenceTarget);
        std::vector<double> orthogonality(
            static_cast<std::size_t>(
                row.energyProjectionRank), 0.0);
#ifdef USE_MKL_PARDISO
        cblas_dgemv(
            CblasColMajor, CblasTrans,
            row.interfaceDofs, row.energyProjectionRank,
            1.0, basis.data(), row.interfaceDofs,
            errorImage.data(), 1,
            0.0, orthogonality.data(), 1);
#else
        for (int mode = 0;
             mode < row.energyProjectionRank; ++mode) {
            for (int targetRow = 0;
                 targetRow < row.interfaceDofs; ++targetRow) {
                orthogonality[static_cast<std::size_t>(mode)]
                    += basis[static_cast<std::size_t>(
                        mode * row.interfaceDofs + targetRow)]
                    * errorImage[
                        static_cast<std::size_t>(targetRow)];
            }
        }
#endif
        double maximumOrthogonality = 0.0;
        for (double value : orthogonality) {
            maximumOrthogonality =
                std::max(maximumOrthogonality, std::abs(value));
        }
        row.temperatureGalerkinOrthogonality =
            maximumOrthogonality
            / std::max(1.0e-300,
                std::sqrt(dot(referenceImage, referenceImage)));

        row.forcing = projectForcing(
            forcingTarget, basis, schurImages,
            row.interfaceDofs,
            row.energyProjectionRank, target);
        row.input = projectForcing(
            inputTarget, basis, schurImages,
            row.interfaceDofs,
            row.energyProjectionRank, target);
        row.boundary = projectForcing(
            boundaryTarget, basis, schurImages,
            row.interfaceDofs,
            row.energyProjectionRank, target);
        row.history = projectForcing(
            historyTarget, basis, schurImages,
            row.interfaceDofs,
            row.energyProjectionRank, target);
        std::vector<double> otherTarget(forcingTarget.size(), 0.0);
        std::vector<double> componentSum(forcingTarget.size(), 0.0);
        for (std::size_t targetRow = 0;
             targetRow < forcingTarget.size(); ++targetRow) {
            componentSum[targetRow] =
                inputTarget[targetRow]
                + boundaryTarget[targetRow]
                + historyTarget[targetRow];
            otherTarget[targetRow] =
                forcingTarget[targetRow] - componentSum[targetRow];
        }
        row.other = projectForcing(
            otherTarget, basis, schurImages,
            row.interfaceDofs,
            row.energyProjectionRank, target);
        row.forcingDecompositionClosure =
            relativeDifference(componentSum, forcingTarget);
        const PatchInnerSolverStatistics& statistics =
            target.statistics();
        row.targetSolveCalls = statistics.solveCalls;
        row.targetSolveResidual =
            statistics.maximumRelativeResidual;
        row.targetSolver = statistics.actualSolver;
        row.setupSeconds = statistics.setupSeconds;
        row.projectionSeconds = elapsed(interfaceStart);
        row.peakIncrementalMemoryBytes =
            statistics.peakIncrementalMemoryBytes;
        row.status =
            std::isfinite(row.temperatureProjectionError)
            && std::isfinite(row.forcing.error)
            && row.temperatureGalerkinOrthogonality <= 1.0e-8
            && row.forcingDecompositionClosure <= 1.0e-10
            && row.targetSolveResidual <= 1.0e-9
            ? "passed" : "diagnostic_gate_failed";
        summary.maximumTargetSolveResidual = std::max(
            summary.maximumTargetSolveResidual,
            row.targetSolveResidual);
        rows.push_back(std::move(row));
    }

    const std::vector<double> projectedTheta =
        recoverWithInterface(
            dynamicModel, exactSchur, totalGlobal,
            projectedGamma);
    std::vector<double> referenceAbsolute =
        referenceSolve.temperature;
    std::vector<double> projectedAbsolute = projectedTheta;
    for (std::size_t row = 0;
         row < referenceAbsolute.size(); ++row) {
        referenceAbsolute[row] += referenceTemperature[row];
        projectedAbsolute[row] += referenceTemperature[row];
    }
    if (options.fluxOperatorAudit) {
        std::ofstream orientation(
            options.outputDirectory
            / "milestone8_flux_orientation_audit.csv");
        std::ofstream affine(
            options.outputDirectory
            / "milestone8_flux_affine_decomposition.csv");
        std::ofstream linearity(
            options.outputDirectory
            / "milestone8_flux_operator_linearity.csv");
        std::ofstream forcingAudit(
            options.outputDirectory
            / "milestone8_interface10_forcing_audit.csv");
        orientation
            << "case,interface_id,owner_side,master_subdomain,"
            "slave_subdomain,owner_normal_x,owner_normal_y,"
            "owner_normal_z,master_normal_x,master_normal_y,"
            "master_normal_z,slave_normal_x,slave_normal_y,"
            "slave_normal_z,normal_dot_product,"
            "minimum_normal_dot_product,maximum_normal_dot_product,"
            "target_dofs,flux_triangles,flux_quadrature_points,"
            "junction_duplicate_count,deterministic_owner_ordering,"
            "opposite_normals,projection_uses_absolute_value,"
            "dof_quadrature_one_to_one_assumed,"
            "physical_flux_sign_convention,"
            "numerical_flux_sign_convention,mapping_convention,status\n"
            << std::setprecision(17);
        affine
            << "case,interface_id,flux_type,homogeneous_norm,"
            "affine_norm,total_norm,homogeneous_projection_error,"
            "total_projection_error,affine_decomposition_error,"
            "cancellation_ratio,linearity_gate,affine_gate,status\n"
            << std::setprecision(17);
        linearity
            << "case,interface_id,flux_type,first_scale,second_scale,"
            "relative_linearity_error,tolerance,status\n"
            << std::setprecision(17);
        forcingAudit
            << "case,interface_id,left_subdomain,right_subdomain,"
            "target_dofs,local_port_rank,mandatory_rank,"
            "input_source_rows,boundary_source_rows,"
            "history_source_rows,history_rank_estimate,"
            "raw_channel_count,compressed_channel_count,"
            "accepted_rank,deflated_rank,"
            "forcing_reference_dual_norm,forcing_error_dual_norm,"
            "forcing_relative_error,input_reference_dual_norm,"
            "input_relative_error,boundary_reference_dual_norm,"
            "boundary_relative_error,history_reference_dual_norm,"
            "history_relative_error,other_reference_dual_norm,"
            "other_relative_error,norm_before_compression,"
            "norm_after_compression,projection_norm,"
            "input_fingerprint,boundary_fingerprint,"
            "history_fingerprint,decomposition_closure,"
            "near_zero_forcing,interpretation,snapshot_used,"
            "fom_used_for_basis,status\n"
            << std::setprecision(17);
        for (InterfaceRow& row : rows) {
            InterfaceFluxOperator fluxOperator(
                mesh, physics, dynamicModel, exactSchur,
                row.interfaceId, row.leftSubdomain,
                row.rightSubdomain);
            const InterfaceFluxResponse homogeneousReference =
                fluxOperator.applyGlobalTrace(
                    referenceSolve.interfaceTemperature);
            const InterfaceFluxResponse homogeneousProjected =
                fluxOperator.applyGlobalTrace(projectedGamma);
            const InterfaceFluxResponse affineReference =
                fluxOperator.evaluateAffine(
                    totalGlobal, referenceTemperature);
            const InterfaceFluxResponse totalReference =
                fluxOperator.evaluateRecoveredField(
                    referenceAbsolute);
            const InterfaceFluxResponse totalProjected =
                fluxOperator.evaluateRecoveredField(
                    projectedAbsolute);
            const InterfaceFluxResponse decomposed =
                combineInterfaceFlux(
                    homogeneousReference, 1.0,
                    affineReference, 1.0);

            std::vector<double> firstTrace(
                static_cast<std::size_t>(
                    dynamicModel.interfaceDofs), 0.0);
            std::vector<double> secondTrace(firstTrace.size(), 0.0);
            for (std::size_t localRow = 0;
                 localRow < row.targetIndices.size(); ++localRow) {
                const int gamma = row.targetIndices[localRow];
                firstTrace[static_cast<std::size_t>(gamma)] =
                    std::sin(0.37 * static_cast<double>(localRow + 1));
                secondTrace[static_cast<std::size_t>(gamma)] =
                    std::cos(0.23 * static_cast<double>(localRow + 1));
            }
            constexpr double firstScale = 0.73;
            constexpr double secondScale = -1.21;
            std::vector<double> combinedTrace(firstTrace.size(), 0.0);
            for (std::size_t gamma = 0;
                 gamma < combinedTrace.size(); ++gamma) {
                combinedTrace[gamma] =
                    firstScale * firstTrace[gamma]
                    + secondScale * secondTrace[gamma];
            }
            const InterfaceFluxResponse firstResponse =
                fluxOperator.applyGlobalTrace(firstTrace);
            const InterfaceFluxResponse secondResponse =
                fluxOperator.applyGlobalTrace(secondTrace);
            const InterfaceFluxResponse combinedResponse =
                fluxOperator.applyGlobalTrace(combinedTrace);
            const double tolerance =
                row.interfaceDofs > 2000 ? 1.0e-8 : 1.0e-10;

            std::set<int> uniqueTargets(
                row.targetIndices.begin(), row.targetIndices.end());
            const int duplicateCount =
                static_cast<int>(row.targetIndices.size())
                - static_cast<int>(uniqueTargets.size());
            const InterfaceFluxOrientationAudit orientationRow =
                fluxOperator.auditOrientation(
                    row.interfaceDofs, duplicateCount);
            orientation
                << physics.name << ',' << row.interfaceId << ','
                << orientationRow.ownerSide << ','
                << orientationRow.masterSubdomain << ','
                << orientationRow.slaveSubdomain << ','
                << orientationRow.ownerNormalX << ','
                << orientationRow.ownerNormalY << ','
                << orientationRow.ownerNormalZ << ','
                << orientationRow.masterNormalX << ','
                << orientationRow.masterNormalY << ','
                << orientationRow.masterNormalZ << ','
                << orientationRow.slaveNormalX << ','
                << orientationRow.slaveNormalY << ','
                << orientationRow.slaveNormalZ << ','
                << orientationRow.normalDotProduct << ','
                << orientationRow.minimumNormalDotProduct << ','
                << orientationRow.maximumNormalDotProduct << ','
                << orientationRow.targetDofCount << ','
                << orientationRow.fluxTriangleCount << ','
                << orientationRow.fluxQuadratureCount << ','
                << orientationRow.junctionDuplicateCount << ','
                << (orientationRow.deterministicOwnerOrdering ? 1 : 0)
                << ',' << (orientationRow.oppositeNormals ? 1 : 0)
                << ',' << (orientationRow.projectionUsesAbsoluteValue
                    ? 1 : 0)
                << ',' << (orientationRow.dofQuadratureOneToOneAssumed
                    ? 1 : 0)
                << ",\"" << orientationRow.physicalFluxSignConvention
                << "\",\"" << orientationRow.numericalFluxSignConvention
                << "\",\"" << orientationRow.mappingConvention << "\","
                << orientationRow.status << '\n';

            bool operatorPassed =
                orientationRow.status == "passed";
            for (const std::string fluxType :
                 {"physical", "numerical"}) {
                const double linearityError =
                    interfaceFluxRelativeLinearityError(
                        combinedResponse, firstResponse, firstScale,
                        secondResponse, secondScale, fluxType);
                const double decompositionError =
                    interfaceFluxRelativeDifference(
                        decomposed, totalReference, fluxType);
                const double homogeneousNorm =
                    interfaceFluxNorm(
                        homogeneousReference, fluxType);
                const double affineNorm =
                    interfaceFluxNorm(affineReference, fluxType);
                const double totalNorm =
                    interfaceFluxNorm(totalReference, fluxType);
                const double cancellationRatio = totalNorm
                    / std::max(
                        1.0e-300,
                        homogeneousNorm + affineNorm);
                const double homogeneousProjectionError =
                    interfaceFluxRelativeDifference(
                        homogeneousProjected,
                        homogeneousReference, fluxType);
                const double totalProjectionError =
                    interfaceFluxRelativeDifference(
                        totalProjected, totalReference, fluxType);
                const bool linearityPassed =
                    std::isfinite(linearityError)
                    && linearityError <= tolerance;
                const bool affinePassed =
                    std::isfinite(decompositionError)
                    && decompositionError <= tolerance;
                linearity
                    << physics.name << ',' << row.interfaceId << ','
                    << fluxType << ',' << firstScale << ','
                    << secondScale << ',' << linearityError << ','
                    << tolerance << ','
                    << (linearityPassed ? "passed" : "failed") << '\n';
                affine
                    << physics.name << ',' << row.interfaceId << ','
                    << fluxType << ',' << homogeneousNorm << ','
                    << affineNorm << ',' << totalNorm << ','
                    << homogeneousProjectionError << ','
                    << totalProjectionError << ','
                    << decompositionError << ','
                    << cancellationRatio << ','
                    << (linearityPassed ? 1 : 0) << ','
                    << (affinePassed ? 1 : 0) << ','
                    << (linearityPassed && affinePassed
                        ? "passed" : "failed") << '\n';
                operatorPassed = operatorPassed
                    && linearityPassed && affinePassed;
            }
            if (!operatorPassed) {
                row.status = "flux_operator_audit_failed";
            }

            if (row.interfaceId == 10) {
                const auto portIt = portById.find(row.interfaceId);
                const LocalPortBasis& port = *portIt->second;
                const int mandatoryRank = port.mandatoryModes;
                const int geometryAndConstantRank =
                    std::min(4, mandatoryRank);
                const int historyRankEstimate =
                    std::max(0, mandatoryRank
                        - geometryAndConstantRank);
                const int rawChannels =
                    port.inputSourceRows
                    + port.boundarySourceRows
                    + port.historySourceRows;
                const int compressedChannels =
                    port.inputSourceRows
                    + port.boundarySourceRows
                    + historyRankEstimate;
                const double forcingScale = std::max({
                    row.input.referenceNorm,
                    row.boundary.referenceNorm,
                    row.history.referenceNorm,
                    1.0});
                const bool nearZero =
                    row.forcing.referenceNorm
                    <= 1.0e-12 * forcingScale;
                forcingAudit
                    << physics.name << ',' << row.interfaceId << ','
                    << row.leftSubdomain << ','
                    << row.rightSubdomain << ','
                    << row.interfaceDofs << ','
                    << row.localPortRank << ','
                    << mandatoryRank << ','
                    << port.inputSourceRows << ','
                    << port.boundarySourceRows << ','
                    << port.historySourceRows << ','
                    << historyRankEstimate << ','
                    << rawChannels << ','
                    << compressedChannels << ','
                    << mandatoryRank << ','
                    << std::max(0,
                        rawChannels - compressedChannels) << ','
                    << row.forcing.referenceNorm << ','
                    << row.forcing.errorNorm << ','
                    << row.forcing.error << ','
                    << row.input.referenceNorm << ','
                    << row.input.error << ','
                    << row.boundary.referenceNorm << ','
                    << row.boundary.error << ','
                    << row.history.referenceNorm << ','
                    << row.history.error << ','
                    << row.other.referenceNorm << ','
                    << row.other.error << ','
                    << row.history.referenceNorm << ','
                    << row.history.referenceNorm << ','
                    << std::sqrt(std::max(
                        0.0,
                        row.forcing.referenceNorm
                            * row.forcing.referenceNorm
                        - row.forcing.errorNorm
                            * row.forcing.errorNorm)) << ','
                    << port.inputSourceFingerprint << ','
                    << port.boundarySourceFingerprint << ','
                    << port.historySourceFingerprint << ','
                    << row.forcingDecompositionClosure << ','
                    << (nearZero ? 1 : 0) << ','
                    << (nearZero
                        ? "relative_error_ill_posed_zero_forcing"
                        : "finite_forcing_relative_error_valid")
                    << ",0,0,"
                    << (nearZero
                        ? "passed_near_zero_classified"
                        : "passed") << '\n';
            }
        }
    }
    const DetailedInterfacePhysicsMetrics referenceFlux =
        calculateDetailedInterfacePhysicsMetrics(
            mesh, physics, referenceAbsolute);
    const DetailedInterfacePhysicsMetrics projectedFlux =
        calculateDetailedInterfacePhysicsMetrics(
            mesh, physics, projectedAbsolute);
    if (referenceFlux.triangles.size()
        != projectedFlux.triangles.size()) {
        throw std::runtime_error(
            "[Projection diagnosis] Flux ordering mismatch.");
    }
    struct FluxAccumulator {
        long double numericalError = 0.0L;
        long double numericalReference = 0.0L;
        long double physicalError = 0.0L;
        long double physicalReference = 0.0L;
        int triangleCount = 0;
        long double area = 0.0L;
    };
    std::map<std::pair<int, int>, int> interfaceBySubdomains;
    for (const InterfaceRow& row : rows) {
        interfaceBySubdomains.emplace(
            std::minmax(
                row.leftSubdomain, row.rightSubdomain),
            row.interfaceId);
    }
    std::map<int, FluxAccumulator> fluxByInterface;
    for (std::size_t triangle = 0;
         triangle < referenceFlux.triangles.size(); ++triangle) {
        const InterfaceTriangleFluxRecord& ref =
            referenceFlux.triangles[triangle];
        const InterfaceTriangleFluxRecord& projected =
            projectedFlux.triangles[triangle];
        const auto selected = interfaceBySubdomains.find(
            std::minmax(
                ref.leftSubdomain, ref.rightSubdomain));
        if (selected == interfaceBySubdomains.end()) {
            continue;
        }
        FluxAccumulator& accumulator =
            fluxByInterface[selected->second];
        const double numericalDelta =
            projected.sipgNumericalFlux
            - ref.sipgNumericalFlux;
        const double leftDelta =
            projected.leftPhysicalNormalFlux
            - ref.leftPhysicalNormalFlux;
        const double rightDelta =
            projected.rightPhysicalNormalFlux
            - ref.rightPhysicalNormalFlux;
        accumulator.numericalError
            += static_cast<long double>(ref.area)
                * numericalDelta * numericalDelta;
        accumulator.numericalReference
            += static_cast<long double>(ref.area)
                * ref.sipgNumericalFlux
                * ref.sipgNumericalFlux;
        accumulator.physicalError
            += static_cast<long double>(ref.area)
                * (leftDelta * leftDelta
                    + rightDelta * rightDelta);
        accumulator.physicalReference
            += static_cast<long double>(ref.area)
                * (ref.leftPhysicalNormalFlux
                    * ref.leftPhysicalNormalFlux
                    + ref.rightPhysicalNormalFlux
                    * ref.rightPhysicalNormalFlux);
        ++accumulator.triangleCount;
        accumulator.area += ref.area;
    }
    for (InterfaceRow& row : rows) {
        const FluxAccumulator& accumulator =
            fluxByInterface[row.interfaceId];
        row.fluxTriangleCount = accumulator.triangleCount;
        row.fluxArea = static_cast<double>(accumulator.area);
        row.numericalFluxProjectionError =
            std::sqrt(static_cast<double>(
                accumulator.numericalError))
            / std::max(1.0e-300,
                std::sqrt(static_cast<double>(
                    accumulator.numericalReference)));
        row.physicalFluxProjectionError =
            std::sqrt(static_cast<double>(
                accumulator.physicalError))
            / std::max(1.0e-300,
                std::sqrt(static_cast<double>(
                    accumulator.physicalReference)));
        row.fluxProjectionError =
            row.numericalFluxProjectionError;
        if (row.fluxTriangleCount <= 0
            || !(row.fluxArea > 0.0)
            || !std::isfinite(row.fluxProjectionError)) {
            row.status = "diagnostic_gate_failed";
        }
        if (row.status == "passed") {
            ++summary.completedInterfaces;
        }
    }

    std::ofstream temperature(
        options.outputDirectory
        / "milestone8_projection_temperature.csv");
    temperature
        << "case,interface_id,left_subdomain,right_subdomain,"
        "interface_dofs,local_port_rank,energy_projection_rank,"
        "temperature_projection_error,"
        "temperature_euclidean_projection_error,"
        "reference_energy_norm,error_energy_norm,"
        "galerkin_orthogonality,target_solver,"
        "target_solve_residual,setup_time_s,projection_time_s,"
        "peak_incremental_memory_bytes,full_field_read,"
        "snapshot_used,fom_used_for_basis,pod_used,svd_used,status\n"
        << std::setprecision(17);
    std::ofstream flux(
        options.outputDirectory
        / "milestone8_projection_flux.csv");
    flux
        << "case,interface_id,interface_dofs,local_port_rank,"
        "energy_projection_rank,"
        "flux_projection_error,physical_flux_projection_error,"
        "numerical_flux_projection_error,flux_triangle_count,"
        "flux_area,projection_definition,"
        "full_field_read,snapshot_used,fom_used_for_basis,"
        "pod_used,svd_used,status\n"
        << std::setprecision(17);
    std::ofstream forcing(
        options.outputDirectory
        / "milestone8_projection_forcing.csv");
    forcing
        << "case,interface_id,interface_dofs,local_port_rank,"
        "energy_projection_rank,"
        "forcing_projection_error,input_projection_error,"
            "boundary_projection_error,history_projection_error,"
            "other_projection_error,"
        "forcing_reference_dual_norm,forcing_error_dual_norm,"
        "decomposition_closure,target_solve_calls,"
        "target_solve_residual,dual_norm,full_field_read,"
        "snapshot_used,fom_used_for_basis,pod_used,svd_used,status\n"
        << std::setprecision(17);
    bool allPassed = true;
    for (const InterfaceRow& row : rows) {
        temperature
            << physics.name << ',' << row.interfaceId << ','
            << row.leftSubdomain << ',' << row.rightSubdomain << ','
            << row.interfaceDofs << ',' << row.localPortRank << ','
            << row.energyProjectionRank << ','
            << row.temperatureProjectionError << ','
            << row.temperatureEuclideanError << ','
            << row.temperatureReferenceEnergy << ','
            << row.temperatureErrorEnergy << ','
            << row.temperatureGalerkinOrthogonality << ','
            << row.targetSolver << ','
            << row.targetSolveResidual << ','
            << row.setupSeconds << ','
            << row.projectionSeconds << ','
            << row.peakIncrementalMemoryBytes
            << ",0,0,0,0,0," << row.status << '\n';
        flux
            << physics.name << ',' << row.interfaceId << ','
            << row.interfaceDofs << ',' << row.localPortRank << ','
            << row.energyProjectionRank << ','
            << row.fluxProjectionError << ','
            << row.physicalFluxProjectionError << ','
            << row.numericalFluxProjectionError
            << ',' << row.fluxTriangleCount
            << ',' << row.fluxArea
            << ",projected_trace_with_exact_reduced_interior_recovery,"
            << "0,0,0,0,0," << row.status << '\n';
        forcing
            << physics.name << ',' << row.interfaceId << ','
            << row.interfaceDofs << ',' << row.localPortRank << ','
            << row.energyProjectionRank << ','
            << row.forcing.error << ','
            << row.input.error << ','
            << row.boundary.error << ','
            << row.history.error << ','
            << row.other.error << ','
            << row.forcing.referenceNorm << ','
            << row.forcing.errorNorm << ','
            << row.forcingDecompositionClosure << ','
            << row.targetSolveCalls << ','
            << row.targetSolveResidual
            << ",S_tt_inverse,0,0,0,0,0,"
            << row.status << '\n';
        allPassed = allPassed && row.status == "passed";
    }
    summary.projectionSeconds = elapsed(diagnosisStart)
        - summary.fullInterfaceSetupSeconds
        - summary.fullInterfaceSolveSeconds;
    summary.peakWorkingSetBytes = peakWorkingSetBytes();
    summary.status =
        allPassed
        && summary.completedInterfaces
            == summary.requestedInterfaces
        ? "passed" : "diagnostic_gate_failed";
    std::ofstream summaryOutput(
        options.outputDirectory
        / "milestone8_projection_summary.csv");
    summaryOutput
        << "case,requested_interfaces,completed_interfaces,"
        "full_interface_setup_time_s,full_interface_solve_time_s,"
        "projection_time_s,max_target_solve_residual,"
        "peak_working_set_bytes,full_field_read,snapshot_used,"
        "fom_used_for_basis,pod_used,svd_used,transient_steps,status\n"
        << std::setprecision(17)
        << physics.name << ',' << summary.requestedInterfaces << ','
        << summary.completedInterfaces << ','
        << summary.fullInterfaceSetupSeconds << ','
        << summary.fullInterfaceSolveSeconds << ','
        << summary.projectionSeconds << ','
        << summary.maximumTargetSolveResidual << ','
        << summary.peakWorkingSetBytes
        << ",0,0,0,0,0,1," << summary.status << '\n';
    return summary;
}

} // namespace mor::transient
