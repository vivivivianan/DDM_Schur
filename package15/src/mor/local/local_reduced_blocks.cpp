#include "local_reduced_blocks.hpp"

#include "ddm_schur/interface_operator.hpp"
#include "local_reduced_schur.hpp"
#include "mor/pod_basis.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace mor::local {
namespace {

double relativeVectorDifference(const std::vector<double>& left,
                                const std::vector<double>& right)
{
    if (left.size() != right.size()) {
        return std::numeric_limits<double>::infinity();
    }
    double difference = 0.0;
    double scale = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        difference = std::max(difference, std::abs(left[index] - right[index]));
        scale = std::max(scale, std::max(std::abs(left[index]), std::abs(right[index])));
    }
    return difference / std::max(1.0e-300, scale);
}

void clearTemplatePayload(SubdomainModel& local)
{
    local.referenceInterior.clear();
    local.referenceInterior.shrink_to_fit();
    local.basis.clear();
    local.basis.shrink_to_fit();
    local.reducedInterior.clear();
    local.reducedInterior.shrink_to_fit();
    local.reducedInteriorInterface.clear();
    local.reducedInteriorInterface.shrink_to_fit();
    local.reducedInterfaceInterior.clear();
    local.reducedInterfaceInterior.shrink_to_fit();
    local.interiorReferenceImage.clear();
    local.interiorReferenceImage.shrink_to_fit();
    local.interfaceReferenceImage.clear();
    local.interfaceReferenceImage.shrink_to_fit();
}

} // namespace

void projectLocalReducedBlocks(Model& model,
                               const ddm_schur::InterfacePartition& partition)
{
    if (model.subdomains.size() != partition.domains.size()
        || model.interfaceDofs != static_cast<int>(partition.interfaceGlobalDofs.size())) {
        throw std::runtime_error("[Local ROM] Partition/model mismatch during projection.");
    }
    const auto projectionStart = std::chrono::steady_clock::now();
    const int gammaSize = model.interfaceDofs;
    model.interfaceEntries.clear();
    model.interfaceEntries.reserve(partition.interfaceEntries.size());
    for (const ddm_schur::Entry& entry : partition.interfaceEntries) {
        model.interfaceEntries.push_back({entry.row, entry.col, entry.value});
    }

    for (std::size_t slot = 0; slot < partition.domains.size(); ++slot) {
        const auto localProjectionStart = std::chrono::steady_clock::now();
        const ddm_schur::DomainBlocks& domain = partition.domains[slot];
        SubdomainModel& local = model.subdomains[slot];
        const int rows = local.interiorDofs;
        const int rank = local.rank;
        const int localGammaSize = local.localInterfaceDofs;
        std::vector<int> gammaToLocal(static_cast<std::size_t>(gammaSize), -1);
        for (int localGamma = 0; localGamma < localGammaSize; ++localGamma) {
            gammaToLocal[static_cast<std::size_t>(
                local.interfaceIndices[static_cast<std::size_t>(localGamma)])] = localGamma;
        }
        local.reducedInterior.assign(static_cast<std::size_t>(rank * rank), 0.0);
        local.reducedInteriorInterface.assign(
            static_cast<std::size_t>(rank * localGammaSize), 0.0);
        local.reducedInterfaceInterior.assign(
            static_cast<std::size_t>(localGammaSize * rank), 0.0);
        local.interiorReferenceImage.assign(static_cast<std::size_t>(rows), 0.0);
        local.interfaceReferenceImage.assign(static_cast<std::size_t>(localGammaSize), 0.0);

        for (const ddm_schur::Entry& entry : domain.interiorEntries) {
            local.interiorReferenceImage[static_cast<std::size_t>(entry.row)] +=
                entry.value * local.referenceInterior[static_cast<std::size_t>(entry.col)];
            for (int left = 0; left < rank; ++left) {
                const double test = local.basis[static_cast<std::size_t>(left * rows + entry.row)];
                for (int right = 0; right < rank; ++right) {
                    local.reducedInterior[static_cast<std::size_t>(left * rank + right)] +=
                        test * entry.value
                        * local.basis[static_cast<std::size_t>(right * rows + entry.col)];
                }
            }
        }
        for (int row = 0; row < rows; ++row) {
            for (const auto& entry : domain.interiorInterfaceRows[static_cast<std::size_t>(row)]) {
                const int localGamma = gammaToLocal[static_cast<std::size_t>(entry.first)];
                if (localGamma < 0) {
                    throw std::runtime_error(
                        "[Local ROM] Interior-interface entry is outside the local trace.");
                }
                for (int mode = 0; mode < rank; ++mode) {
                    local.reducedInteriorInterface[static_cast<std::size_t>(
                        mode * localGammaSize + localGamma)] +=
                        local.basis[static_cast<std::size_t>(mode * rows + row)] * entry.second;
                }
            }
        }
        for (std::size_t localGamma = 0;
             localGamma < domain.interfaceInteriorRows.size(); ++localGamma) {
            for (const auto& entry : domain.interfaceInteriorRows[localGamma]) {
                local.interfaceReferenceImage[localGamma] +=
                    entry.second * local.referenceInterior[static_cast<std::size_t>(entry.first)];
                for (int mode = 0; mode < rank; ++mode) {
                    local.reducedInterfaceInterior[static_cast<std::size_t>(
                        localGamma * static_cast<std::size_t>(rank) + mode)] += entry.second
                        * local.basis[static_cast<std::size_t>(mode * rows + entry.first)];
                }
            }
        }

        double matrixScale = 0.0;
        double symmetryDifference = 0.0;
        std::vector<double> symmetricInterior = local.reducedInterior;
        for (int row = 0; row < rank; ++row) {
            for (int column = 0; column < rank; ++column) {
                const double left = local.reducedInterior[static_cast<std::size_t>(
                    row * rank + column)];
                const double right = local.reducedInterior[static_cast<std::size_t>(
                    column * rank + row)];
                matrixScale = std::max(matrixScale, std::max(std::abs(left), std::abs(right)));
                symmetryDifference = std::max(symmetryDifference, std::abs(left - right));
                symmetricInterior[static_cast<std::size_t>(row * rank + column)] =
                    0.5 * (left + right);
            }
        }
        local.reducedInteriorSymmetryError = symmetryDifference
            / std::max(1.0e-300, matrixScale);
        const std::vector<double> eigenvalues = symmetricEigenvalues(
            std::move(symmetricInterior), rank);
        if (!eigenvalues.empty()) {
            local.reducedInteriorMaximumEigenvalue = eigenvalues.front();
            local.reducedInteriorMinimumEigenvalue = eigenvalues.back();
            local.reducedInteriorConditionEstimate = std::abs(eigenvalues.front())
                / std::max(1.0e-300, std::abs(eigenvalues.back()));
        }
        if (!(local.reducedInteriorMinimumEigenvalue > 0.0)) {
            throw std::runtime_error(
                "[Local ROM] Reduced interior block is not positive definite.");
        }

        double couplingScale = 0.0;
        double couplingDifference = 0.0;
        for (int mode = 0; mode < rank; ++mode) {
            for (int localGamma = 0; localGamma < localGammaSize; ++localGamma) {
                const double interiorInterface = local.reducedInteriorInterface[
                    static_cast<std::size_t>(mode * localGammaSize + localGamma)];
                const double interfaceInterior = local.reducedInterfaceInterior[
                    static_cast<std::size_t>(localGamma * rank + mode)];
                couplingScale = std::max(couplingScale,
                    std::max(std::abs(interiorInterface), std::abs(interfaceInterior)));
                couplingDifference = std::max(couplingDifference,
                    std::abs(interiorInterface - interfaceInterior));
            }
        }
        local.couplingSymmetryError = couplingDifference
            / std::max(1.0e-300, couplingScale);

        const DenseSymmetricFactor reducedFactor = factorDenseSymmetric(
            local.reducedInterior, rank);
        double schurScale = 0.0;
        double schurDifference = 0.0;
        constexpr int exactSchurDiagnosticLimit = 2000;
        if (localGammaSize <= exactSchurDiagnosticLimit) {
            std::vector<double> solvedColumns(static_cast<std::size_t>(rank)
                * static_cast<std::size_t>(localGammaSize), 0.0);
            for (int column = 0; column < localGammaSize; ++column) {
                std::vector<double> image(static_cast<std::size_t>(rank), 0.0);
                for (int mode = 0; mode < rank; ++mode) {
                    image[static_cast<std::size_t>(mode)] =
                        local.reducedInteriorInterface[static_cast<std::size_t>(
                            mode * localGammaSize + column)];
                }
                solveDenseSymmetric(reducedFactor, image);
                for (int mode = 0; mode < rank; ++mode) {
                    solvedColumns[static_cast<std::size_t>(
                        mode * localGammaSize + column)] =
                        image[static_cast<std::size_t>(mode)];
                }
            }
            for (int row = 0; row < localGammaSize; ++row) {
                for (int column = row; column < localGammaSize; ++column) {
                    double left = 0.0;
                    double right = 0.0;
                    for (int mode = 0; mode < rank; ++mode) {
                        left += local.reducedInterfaceInterior[static_cast<std::size_t>(
                            row * rank + mode)] * solvedColumns[static_cast<std::size_t>(
                            mode * localGammaSize + column)];
                        right += local.reducedInterfaceInterior[static_cast<std::size_t>(
                            column * rank + mode)] * solvedColumns[static_cast<std::size_t>(
                            mode * localGammaSize + row)];
                    }
                    schurScale = std::max(schurScale,
                        std::max(std::abs(left), std::abs(right)));
                    schurDifference = std::max(
                        schurDifference, std::abs(left - right));
                }
            }
        } else {
            // An explicit local Gamma-by-Gamma diagnostic would be quadratic
            // and dominates RRAM26 setup.  Test symmetry through deterministic
            // bilinear forms x^T S_i y and y^T S_i x instead.  The operator
            // itself remains exact and matrix-free.
            for (int trial = 0; trial < 6; ++trial) {
                std::vector<double> leftVector(
                    static_cast<std::size_t>(localGammaSize), 0.0);
                std::vector<double> rightVector(
                    static_cast<std::size_t>(localGammaSize), 0.0);
                for (int gamma = 0; gamma < localGammaSize; ++gamma) {
                    leftVector[static_cast<std::size_t>(gamma)] =
                        std::sin(0.013 * static_cast<double>(
                            (gamma + 1) * (trial + 1)));
                    rightVector[static_cast<std::size_t>(gamma)] =
                        std::cos(0.017 * static_cast<double>(
                            (gamma + 3) * (trial + 2)));
                }
                std::vector<double> leftReduced(static_cast<std::size_t>(rank), 0.0);
                std::vector<double> rightReduced(static_cast<std::size_t>(rank), 0.0);
                for (int mode = 0; mode < rank; ++mode) {
                    for (int gamma = 0; gamma < localGammaSize; ++gamma) {
                        const double coupling = local.reducedInteriorInterface[
                            static_cast<std::size_t>(mode * localGammaSize + gamma)];
                        leftReduced[static_cast<std::size_t>(mode)] += coupling
                            * leftVector[static_cast<std::size_t>(gamma)];
                        rightReduced[static_cast<std::size_t>(mode)] += coupling
                            * rightVector[static_cast<std::size_t>(gamma)];
                    }
                }
                solveDenseSymmetric(reducedFactor, leftReduced);
                solveDenseSymmetric(reducedFactor, rightReduced);
                double left = 0.0;
                double right = 0.0;
                for (int gamma = 0; gamma < localGammaSize; ++gamma) {
                    double imageLeft = 0.0;
                    double imageRight = 0.0;
                    for (int mode = 0; mode < rank; ++mode) {
                        const double coupling = local.reducedInterfaceInterior[
                            static_cast<std::size_t>(gamma * rank + mode)];
                        imageLeft += coupling * leftReduced[static_cast<std::size_t>(mode)];
                        imageRight += coupling * rightReduced[static_cast<std::size_t>(mode)];
                    }
                    left += leftVector[static_cast<std::size_t>(gamma)] * imageRight;
                    right += rightVector[static_cast<std::size_t>(gamma)] * imageLeft;
                }
                schurScale = std::max(schurScale,
                    std::max(std::abs(left), std::abs(right)));
                schurDifference = std::max(
                    schurDifference, std::abs(left - right));
            }
        }
        local.localSchurSymmetryError = schurDifference
            / std::max(1.0e-300, schurScale);
        local.projectionSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - localProjectionStart).count();
    }

    std::map<int, std::size_t> templateRepresentatives;
    for (std::size_t slot = 0; slot < model.subdomains.size(); ++slot) {
        SubdomainModel& local = model.subdomains[slot];
        if (!local.templateReused) {
            templateRepresentatives.emplace(local.templateId, slot);
            continue;
        }
        const auto representative = templateRepresentatives.find(local.templateId);
        if (representative == templateRepresentatives.end()) {
            throw std::runtime_error(
                "[Local ROM] Reused template precedes its representative.");
        }
        const SubdomainModel& source = model.subdomains[representative->second];
        local.templateConsistencyDifference = std::max({
            relativeVectorDifference(local.referenceInterior, source.referenceInterior),
            relativeVectorDifference(local.basis, source.basis),
            relativeVectorDifference(local.reducedInterior, source.reducedInterior),
            relativeVectorDifference(local.reducedInteriorInterface,
                source.reducedInteriorInterface),
            relativeVectorDifference(local.reducedInterfaceInterior,
                source.reducedInterfaceInterior),
            relativeVectorDifference(local.interiorReferenceImage,
                source.interiorReferenceImage),
            relativeVectorDifference(local.interfaceReferenceImage,
                source.interfaceReferenceImage),
            std::abs(local.localSchurSymmetryError - source.localSchurSymmetryError)});
        if (!(local.templateConsistencyDifference <= 1.0e-10)) {
            throw std::runtime_error(
                "[Local ROM] Template reuse rejected: projected local blocks differ.");
        }
        clearTemplatePayload(local);
    }
    model.projectionSeconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - projectionStart).count();

    model.schurConstructionSeconds = 0.0;
    model.modelBytes = estimateModelBytes(model);
}

} // namespace mor::local
