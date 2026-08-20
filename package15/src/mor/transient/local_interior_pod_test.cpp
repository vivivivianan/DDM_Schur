#include "local_interior_pod_test.hpp"

#include "../../sipg_core.hpp"
#include "../../config_io.hpp"
#include "../../mesh_loader.hpp"
#include "../../fem_assembly.hpp"
#include "../../linear_solvers.hpp"
#include "../../diagnostics_io.hpp"
#include "../../ddm_schur/interface_operator.hpp"
#include "../../ddm_schur/schur_fgmres.hpp"
#include "../pod_basis.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace {

using Clock = std::chrono::steady_clock;

double elapsed(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

double matrixValue(const SparseMatrix& matrix, int row, int col)
{
    if (matrix.n == 0) {
        return 0.0;
    }
    if (matrix.csrReady) {
        const int begin = matrix.rowPtr[static_cast<std::size_t>(row)];
        const int end = matrix.rowPtr[static_cast<std::size_t>(row + 1)];
        auto it = std::lower_bound(matrix.colInd.begin() + begin,
                                   matrix.colInd.begin() + end, col);
        if (it != matrix.colInd.begin() + end && *it == col) {
            return matrix.values[static_cast<std::size_t>(it - matrix.colInd.begin())];
        }
        return 0.0;
    }
    double result = 0.0;
    for (const MatrixEntry& entry : matrix.triplets) {
        if (entry.row == row && entry.col == col) {
            result += entry.value;
        }
    }
    return result;
}

struct DenseCholesky {
    int n = 0;
    std::vector<double> lower;

    bool factor(const std::vector<double>& a, int size)
    {
        n = size;
        lower.assign(a.begin(), a.end());
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                double value = lower[static_cast<std::size_t>(i * n + j)];
                for (int k = 0; k < j; ++k) {
                    value -= lower[static_cast<std::size_t>(i * n + k)]
                        * lower[static_cast<std::size_t>(j * n + k)];
                }
                if (i == j) {
                    if (!(value > 0.0) || !std::isfinite(value)) {
                        return false;
                    }
                    lower[static_cast<std::size_t>(i * n + j)] = std::sqrt(value);
                } else {
                    lower[static_cast<std::size_t>(i * n + j)] = value
                        / lower[static_cast<std::size_t>(j * n + j)];
                }
            }
            for (int j = i + 1; j < n; ++j) {
                lower[static_cast<std::size_t>(i * n + j)] = 0.0;
            }
        }
        return true;
    }

    std::vector<double> solve(const std::vector<double>& rhs) const
    {
        if (static_cast<int>(rhs.size()) != n) {
            throw std::runtime_error("[Local Interior MOR] dense solve dimension mismatch.");
        }
        std::vector<double> result(rhs);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                result[static_cast<std::size_t>(i)] -=
                    lower[static_cast<std::size_t>(i * n + j)]
                    * result[static_cast<std::size_t>(j)];
            }
            result[static_cast<std::size_t>(i)] /=
                lower[static_cast<std::size_t>(i * n + i)];
        }
        for (int i = n - 1; i >= 0; --i) {
            for (int j = i + 1; j < n; ++j) {
                result[static_cast<std::size_t>(i)] -=
                    lower[static_cast<std::size_t>(j * n + i)]
                    * result[static_cast<std::size_t>(j)];
            }
            result[static_cast<std::size_t>(i)] /=
                lower[static_cast<std::size_t>(i * n + i)];
        }
        return result;
    }
};

struct FomTrajectory {
    int columns = 0;
    int interfaceDofs = 0;
    std::vector<std::vector<double>> interiorSnapshots;
    std::vector<double> interfaceSnapshots;
    std::vector<double> maxTemperature;
    double squaredNorm = 0.0;
    double solveSeconds = 0.0;
    double averageIterations = 0.0;
    int maximumIterations = 0;
};

struct LocalReducedModel {
    int domainId = -1;
    int nInterior = 0;
    int nInterface = 0;
    int rank = 0;
    mor::PodResult pod;
    std::vector<double> reducedMass;
    std::vector<double> reducedStiffness;
    // Projected interior/interface blocks are retained explicitly so mass
    // coupling is not silently dropped when a case has C_IΓ terms.
    std::vector<double> reducedMassInterface;
    std::vector<double> reducedStiffnessInterface;
    std::vector<double> reducedSystem;
    DenseCholesky factor;
};

ddm_schur::Options makeSchurOptions(const ProgramOptions& options,
                                    const CaseConfig& physics)
{
    ddm_schur::Options result;
    result.maxIterations = std::max(1, options.schwarzMaxIters);
    result.restart = std::max(1, options.gmresRestart);
    result.relativeTolerance = std::min(options.pcgTolerance, 1.0e-10);
    result.coarseLinearXY = options.schurLinearXYCoarse;
    result.coarseLinearZ = options.schurLinearZCoarse;
    result.coarseGlobalQuadraticZ = options.schurGlobalQuadraticZCoarse;
    result.coarseInterfacePatches = options.schurInterfacePatchCoarse;
    result.coarseInterfaceLinearXY = options.schurInterfacePatchLinearXY;
    result.coarseEnergyAdaptive = options.schurEnergyAdaptiveCoarse;
    result.energyMaxModesPerDomain = options.schurEnergyMaxModesPerDomain;
    result.energySubspaceIterations = options.schurEnergySubspaceIterations;
    result.energyEigenvalueThreshold = options.schurEnergyEigenvalueThreshold;
    result.coarseGlobalSlow = options.schurGlobalSlowCoarse;
    result.globalSlowModes = options.schurGlobalSlowModes;
    result.globalSlowSubspaceDimension = options.schurGlobalSlowSubspaceDimension;
    result.proxyDiagnostics = options.schurProxyDiagnostics;
    result.proxyEnabled = options.schurProxyEnabled;
    result.proxyDisableCoarse = options.schurProxyDisableCoarse;
    result.proxyHighConductivityThreshold = options.schurProxyHighKThreshold;
    result.proxyUseMaterialConnectivity = options.schurProxyUseMaterialConnectivity;
    result.proxyRing = options.schurProxyRing;
    result.proxyProbeColumns = options.schurProxyProbeColumns;
    result.proxyBlockSize = options.schurProxyBlockSize;
    result.proxyValidateBlockEquivalence = options.schurProxyValidateBlockEquivalence;
    result.localSolveThreads = options.schurLocalSolveThreads;
    result.localPardisoThreads = options.schurLocalPardisoThreads;
    result.proxyCacheEnabled = options.schurProxyCacheEnabled;
    result.proxyCachePath = options.schurProxyCachePath.string();
    result.proxyOutputDirectory = options.outputDir.string();
    (void)physics;
    return result;
}

void appendSnapshot(const ddm_schur::InterfacePartition& partition,
                    const std::vector<double>& temperature,
                    FomTrajectory& trajectory)
{
    for (std::size_t d = 0; d < partition.domains.size(); ++d) {
        const auto& ids = partition.domains[d].interiorGlobalDofs;
        auto& snapshot = trajectory.interiorSnapshots[d];
        snapshot.reserve(snapshot.size() + ids.size());
        for (int global : ids) {
            snapshot.push_back(temperature[static_cast<std::size_t>(global)]);
        }
    }
    for (int global : partition.interfaceGlobalDofs) {
        trajectory.interfaceSnapshots.push_back(temperature[static_cast<std::size_t>(global)]);
    }
    trajectory.columns += 1;
}

std::vector<double> storedFullTemperature(const ddm_schur::InterfacePartition& partition,
                                          const FomTrajectory& trajectory,
                                          int column)
{
    std::vector<double> result(static_cast<std::size_t>(partition.totalDofs), 0.0);
    for (int g = 0; g < static_cast<int>(partition.interfaceGlobalDofs.size()); ++g) {
        result[static_cast<std::size_t>(partition.interfaceGlobalDofs[static_cast<std::size_t>(g)])] =
            trajectory.interfaceSnapshots[static_cast<std::size_t>(column * trajectory.interfaceDofs + g)];
    }
    for (std::size_t d = 0; d < partition.domains.size(); ++d) {
        const auto& ids = partition.domains[d].interiorGlobalDofs;
        const auto& snapshot = trajectory.interiorSnapshots[d];
        for (int i = 0; i < static_cast<int>(ids.size()); ++i) {
            result[static_cast<std::size_t>(ids[static_cast<std::size_t>(i)])] =
                snapshot[static_cast<std::size_t>(column * static_cast<int>(ids.size()) + i)];
        }
    }
    return result;
}

void fillReducedMatrix(const SparseMatrix& matrix,
                       const ddm_schur::DomainBlocks& domain,
                       const mor::PodResult& pod,
                       std::vector<double>& reduced)
{
    const int n = static_cast<int>(domain.interiorGlobalDofs.size());
    const int r = pod.selectedRank;
    reduced.assign(static_cast<std::size_t>(r * r), 0.0);
    // Traverse the sparse interior pattern, never the dense n-by-n block.
    // The dynamic partition supplies exactly the nonzero local pattern.
    for (const ddm_schur::Entry& entry : domain.interiorEntries) {
        const int i = entry.row;
        const int j = entry.col;
        if (i < 0 || i >= n || j < 0 || j >= n) {
            continue;
        }
        const int gi = domain.interiorGlobalDofs[static_cast<std::size_t>(i)];
        const int gj = domain.interiorGlobalDofs[static_cast<std::size_t>(j)];
        const double a = matrixValue(matrix, gi, gj);
        if (a == 0.0) {
            continue;
        }
        for (int p = 0; p < r; ++p) {
            const double vip = pod.basis[static_cast<std::size_t>(i + p * n)];
            for (int q = 0; q < r; ++q) {
                reduced[static_cast<std::size_t>(p * r + q)] += a * vip
                    * pod.basis[static_cast<std::size_t>(j + q * n)];
            }
        }
    }
}

void fillReducedCoupling(const SparseMatrix& matrix,
                         const ddm_schur::DomainBlocks& domain,
                         const ddm_schur::InterfacePartition& partition,
                         const mor::PodResult& pod,
                         std::vector<double>& reduced)
{
    const int n = static_cast<int>(domain.interiorGlobalDofs.size());
    const int r = pod.selectedRank;
    const int gamma = static_cast<int>(partition.interfaceGlobalDofs.size());
    reduced.assign(static_cast<std::size_t>(r * gamma), 0.0);
    for (int i = 0; i < n; ++i) {
        const int row = domain.interiorGlobalDofs[static_cast<std::size_t>(i)];
        for (const auto& entry : domain.interiorInterfaceRows[static_cast<std::size_t>(i)]) {
            if (entry.first < 0 || entry.first >= gamma) {
                continue;
            }
            const int col = partition.interfaceGlobalDofs[static_cast<std::size_t>(entry.first)];
            const double value = matrixValue(matrix, row, col);
            for (int p = 0; p < r; ++p) {
                reduced[static_cast<std::size_t>(p * gamma + entry.first)] +=
                    pod.basis[static_cast<std::size_t>(i + p * n)] * value;
            }
        }
    }
}

FomTrajectory buildFomTrajectory(const Mesh& mesh,
                                 const CaseConfig& physics,
                                 const SparseMatrix& mass,
                                 const SparseMatrix& system,
                                 const std::vector<double>& source,
                                 const std::vector<double>& fixedAdjust,
                                 const ProgramOptions& options,
                                 const ddm_schur::InterfacePartition& partition,
                                 int timeSteps)
{
    FomTrajectory trajectory;
    trajectory.interfaceDofs = static_cast<int>(partition.interfaceGlobalDofs.size());
    trajectory.interiorSnapshots.resize(partition.domains.size());
    const auto start = Clock::now();
    ddm_schur::DdmSchurSolver solver(mesh, system, physics, makeSchurOptions(options, physics));
    std::vector<double> temperature = initialTemperatureVector(mesh, physics);
    appendSnapshot(partition, temperature, trajectory);
    double iterationSum = 0.0;
    for (int step = 1; step <= timeSteps; ++step) {
        std::vector<double> rhs = makeTransientRhs(mass, temperature, source, physics.timeStep);
        applyDirichletRhs(mesh, fixedAdjust, rhs);
        const ddm_schur::SolveResult result = solver.solve(rhs);
        if (result.report.status != "success" || vectorHasNonFinite(result.temperature)) {
            throw std::runtime_error("[Local Interior MOR] FOM Dynamic Schur failed at transient step "
                + std::to_string(step) + ".");
        }
        temperature = result.temperature;
        iterationSum += result.report.iterations;
        trajectory.maximumIterations = std::max(trajectory.maximumIterations, result.report.iterations);
        const auto extrema = std::minmax_element(temperature.begin(), temperature.end());
        trajectory.maxTemperature.push_back(*extrema.second);
        double norm = 0.0;
        for (double value : temperature) {
            norm += value * value;
        }
        trajectory.squaredNorm += norm;
        appendSnapshot(partition, temperature, trajectory);
    }
    trajectory.solveSeconds = elapsed(start);
    trajectory.averageIterations = timeSteps > 0
        ? iterationSum / static_cast<double>(timeSteps) : 0.0;
    return trajectory;
}

double norm2(const std::vector<double>& values)
{
    double sum = 0.0;
    for (double value : values) {
        sum += value * value;
    }
    return std::sqrt(sum);
}

} // namespace

int runTransientLocalInteriorRomTest(
    const Mesh& mesh,
    const CaseConfig& physics,
    const SparseMatrix& mass,
    const SparseMatrix& stiffness,
    const SparseMatrix& system,
    const std::vector<double>& source,
    const std::vector<double>& fixedAdjust,
    const ProgramOptions& options,
    const std::filesystem::path& outputDir)
{
    if (physics.timeSteps <= 0 || !(physics.timeStep > 0.0)) {
        throw std::runtime_error("[Local Interior MOR] transient-local-rom-test requires positive time_steps and time_step.");
    }
    std::filesystem::create_directories(outputDir);
    const auto partition = ddm_schur::buildInterfacePartition(mesh, system);
    const int timeSteps = options.transientLocalRomMaxSteps > 0
        ? std::min(physics.timeSteps, options.transientLocalRomMaxSteps)
        : physics.timeSteps;
    std::cout << "[Local Interior MOR] interface-prescribed validation: "
              << partition.interfaceGlobalDofs.size() << " interface DOFs, "
              << partition.domains.size() << " subdomains, "
              << timeSteps << " time steps.\n";
    const FomTrajectory trajectory = buildFomTrajectory(
        mesh, physics, mass, system, source, fixedAdjust, options, partition, timeSteps);
    std::ofstream podOut(outputDir / "transient_local_interior_pod.csv");
    podOut << "subdomain_id,interior_dofs,snapshot_columns,requested_rank,chosen_rank,"
              "numerical_rank,gram_time_s,eigensolve_time_s,basis_time_s,energy_captured\n";
    std::ofstream singularOut(outputDir / "transient_local_interior_singular_values.csv");
    singularOut << "subdomain_id,mode,singular_value,retained_energy\n";
    std::ofstream summary(outputDir / "transient_local_interior_rom_summary.csv");
    summary << "rank,offline_snapshot_time_s,svd_time_s,basis_construction_time_s,"
                "reduced_matrix_time_s,offline_total_time_s,online_time_s,average_time_per_step_s,"
                "fom_time_s,interface_iterations_avg,interface_iterations_max,energy_captured,"
                "relative_l2,max_temperature_error,max_pointwise_error,max_temperature_curve_difference,"
                "maximum_full_residual,online_only_speedup,total_speedup,status\n";

    const double snapshotSeconds = trajectory.solveSeconds;
    for (int requestedRank : options.transientLocalRomRanks) {
        if (requestedRank <= 0) {
            continue;
        }
        const auto rankStart = Clock::now();
        std::vector<LocalReducedModel> models;
        models.reserve(partition.domains.size());
        double svdSeconds = 0.0;
        double basisSeconds = 0.0;
        double matrixSeconds = 0.0;
        double minimumEnergy = 1.0;
        for (std::size_t d = 0; d < partition.domains.size(); ++d) {
            const auto& domain = partition.domains[d];
            LocalReducedModel model;
            model.domainId = domain.domainId;
            model.nInterior = static_cast<int>(domain.interiorGlobalDofs.size());
            model.nInterface = static_cast<int>(domain.interfaceGlobalDofs.size());
            mor::SnapshotDatabase snapshots;
            snapshots.rows = model.nInterior;
            snapshots.cases.resize(static_cast<std::size_t>(trajectory.columns));
            snapshots.values = trajectory.interiorSnapshots[d];
            for (int column = 0; column < trajectory.columns; ++column) {
                snapshots.cases[static_cast<std::size_t>(column)].index = column;
            }
            const auto pod = mor::buildGramPod(snapshots,
                                               requestedRank,
                                               1.0e-12,
                                               options.transientLocalRomSvdTolerance);
            model.pod = pod;
            model.rank = pod.selectedRank;
            svdSeconds += pod.gramSeconds + pod.eigenSeconds;
            basisSeconds += pod.basisSeconds;
            const double energy = pod.retainedEnergy.empty()
                ? 0.0 : pod.retainedEnergy[static_cast<std::size_t>(pod.selectedRank - 1)];
            minimumEnergy = std::min(minimumEnergy, energy);
            podOut << domain.domainId << ',' << model.nInterior << ',' << trajectory.columns << ','
                   << requestedRank << ',' << pod.selectedRank << ',' << pod.numericalRank << ','
                   << pod.gramSeconds << ',' << pod.eigenSeconds << ',' << pod.basisSeconds << ','
                   << energy << '\n';
            for (int mode = 0; mode < static_cast<int>(pod.singularValues.size()); ++mode) {
                singularOut << domain.domainId << ',' << (mode + 1) << ','
                            << pod.singularValues[static_cast<std::size_t>(mode)] << ','
                            << pod.retainedEnergy[static_cast<std::size_t>(mode)] << '\n';
            }
            const auto matrixStart = Clock::now();
            fillReducedMatrix(mass, domain, model.pod, model.reducedMass);
            fillReducedMatrix(stiffness, domain, model.pod, model.reducedStiffness);
            fillReducedCoupling(mass, domain, partition, model.pod, model.reducedMassInterface);
            fillReducedCoupling(stiffness, domain, partition, model.pod, model.reducedStiffnessInterface);
            const int r = model.rank;
            model.reducedSystem.assign(static_cast<std::size_t>(r * r), 0.0);
            for (int i = 0; i < r * r; ++i) {
                model.reducedSystem[static_cast<std::size_t>(i)] =
                    model.reducedMass[static_cast<std::size_t>(i)] / physics.timeStep
                    + model.reducedStiffness[static_cast<std::size_t>(i)];
            }
            // Dirichlet elimination is applied to the assembled dynamic
            // matrix after K is built.  Use that exact A_II block for the
            // solve; Cr/Kr above remain the uneliminated reduced diagnostics.
            std::vector<double> reducedDynamic;
            fillReducedMatrix(system, domain, model.pod, reducedDynamic);
            model.reducedSystem.swap(reducedDynamic);
            if (!model.factor.factor(model.reducedSystem, r)) {
                throw std::runtime_error("[Local Interior MOR] reduced local C/dt+K is not SPD for subdomain "
                    + std::to_string(domain.domainId) + ".");
            }
            matrixSeconds += elapsed(matrixStart);
            models.push_back(std::move(model));
        }

        const auto onlineStart = Clock::now();
        double errorSquared = 0.0;
        double maximumPointwise = 0.0;
        double maximumCurveDifference = 0.0;
        double maximumResidual = 0.0;
        std::vector<double> previous = storedFullTemperature(partition, trajectory, 0);
        for (int step = 1; step <= timeSteps; ++step) {
            std::vector<double> rhs = makeTransientRhs(mass, previous, source, physics.timeStep);
            applyDirichletRhs(mesh, fixedAdjust, rhs);
            std::vector<double> current(static_cast<std::size_t>(partition.totalDofs), 0.0);
            const int gammaColumn = step;
            for (int g = 0; g < trajectory.interfaceDofs; ++g) {
                current[static_cast<std::size_t>(partition.interfaceGlobalDofs[static_cast<std::size_t>(g)])] =
                    trajectory.interfaceSnapshots[static_cast<std::size_t>(gammaColumn * trajectory.interfaceDofs + g)];
            }
            for (std::size_t d = 0; d < partition.domains.size(); ++d) {
                const auto& domain = partition.domains[d];
                const auto& model = models[d];
                const int nInterior = model.nInterior;
                std::vector<double> projected(static_cast<std::size_t>(model.rank), 0.0);
                for (int i = 0; i < nInterior; ++i) {
                    double value = rhs[static_cast<std::size_t>(domain.interiorGlobalDofs[static_cast<std::size_t>(i)])];
                    for (const auto& coupling : domain.interiorInterfaceRows[static_cast<std::size_t>(i)]) {
                        value -= coupling.second * trajectory.interfaceSnapshots[
                            static_cast<std::size_t>(gammaColumn * trajectory.interfaceDofs + coupling.first)];
                    }
                    for (int p = 0; p < model.rank; ++p) {
                        projected[static_cast<std::size_t>(p)] +=
                            model.pod.basis[static_cast<std::size_t>(i + p * nInterior)] * value;
                    }
                }
                const std::vector<double> coefficients = model.factor.solve(projected);
                for (int i = 0; i < nInterior; ++i) {
                    double value = 0.0;
                    for (int p = 0; p < model.rank; ++p) {
                        value += model.pod.basis[static_cast<std::size_t>(i + p * nInterior)]
                            * coefficients[static_cast<std::size_t>(p)];
                    }
                    current[static_cast<std::size_t>(domain.interiorGlobalDofs[static_cast<std::size_t>(i)])] = value;
                }
            }
            const std::vector<double> fom = storedFullTemperature(partition, trajectory, step);
            double stepErrorSquared = 0.0;
            for (int i = 0; i < partition.totalDofs; ++i) {
                const double difference = current[static_cast<std::size_t>(i)] - fom[static_cast<std::size_t>(i)];
                stepErrorSquared += difference * difference;
                maximumPointwise = std::max(maximumPointwise, std::abs(difference));
            }
            errorSquared += stepErrorSquared;
            const auto currentExtrema = std::max_element(current.begin(), current.end());
            const double fomMaximum = *std::max_element(fom.begin(), fom.end());
            maximumCurveDifference = std::max(maximumCurveDifference,
                                              std::abs(*currentExtrema - fomMaximum));
            maximumResidual = std::max(maximumResidual,
                                       relativeResidualNorm(system, current, rhs));
            previous = std::move(current);
        }
        const double onlineSeconds = elapsed(onlineStart);
        const double relativeL2 = trajectory.squaredNorm > 0.0
            ? std::sqrt(errorSquared / trajectory.squaredNorm) : 0.0;
        const double offlineTotal = snapshotSeconds + svdSeconds + basisSeconds + matrixSeconds;
        const double total = offlineTotal + onlineSeconds;
        const double onlineSpeedup = onlineSeconds > 0.0 ? snapshotSeconds / onlineSeconds : 0.0;
        const double totalSpeedup = total > 0.0 ? snapshotSeconds / total : 0.0;
        summary << requestedRank << ',' << snapshotSeconds << ',' << svdSeconds << ','
                << basisSeconds << ',' << matrixSeconds << ',' << offlineTotal << ','
                << onlineSeconds << ',' << onlineSeconds / std::max(1, timeSteps) << ','
                << snapshotSeconds << ',' << trajectory.averageIterations << ','
                << trajectory.maximumIterations << ',' << minimumEnergy << ',' << relativeL2 << ','
                << maximumCurveDifference << ',' << maximumPointwise << ','
                << maximumCurveDifference << ',' << maximumResidual << ',' << onlineSpeedup << ','
                << totalSpeedup << ",success\n";
        std::cout << std::setprecision(10)
                  << "[Local Interior MOR] rank " << requestedRank
                  << " chosen=" << models.front().rank
                  << " offline=" << offlineTotal << " s"
                  << " online=" << onlineSeconds << " s"
                  << " relative_L2=" << relativeL2
                  << " max_error=" << maximumPointwise
                  << " total_speedup=" << totalSpeedup << "\n";
        (void)rankStart;
    }
    std::cout << "[Local Interior MOR] wrote "
              << (outputDir / "transient_local_interior_rom_summary.csv") << "\n";
    return 0;
}
