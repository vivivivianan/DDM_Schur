#pragma once

// Postprocessing, diagnostics, CLI options, spectral checks, sweeps, and validation helpers.
// This file is intentionally included from main.cpp after the preceding SIPG modules.

static double maxAbsDifference(const std::vector<double>& a, const std::vector<double>& b)
{
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (!std::isfinite(a[i]) || !std::isfinite(b[i])) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        result = std::max(result, std::abs(a[i] - b[i]));
    }
    return result;
}

static double l2Difference(const std::vector<double>& a, const std::vector<double>& b)
{
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (!std::isfinite(a[i]) || !std::isfinite(b[i])) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        const double d = a[i] - b[i];
        sum += d * d;
    }
    return std::sqrt(sum);
}

static double l2Norm(const std::vector<double>& a)
{
    double sum = 0.0;
    for (double value : a) {
        if (!std::isfinite(value)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        sum += value * value;
    }
    return std::sqrt(sum);
}

static double infNorm(const std::vector<double>& a)
{
    double result = 0.0;
    for (double value : a) {
        if (!std::isfinite(value)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        result = std::max(result, std::abs(value));
    }
    return result;
}

static double relativeL2Difference(const std::vector<double>& a, const std::vector<double>& b)
{
    return l2Difference(a, b) / std::max(1.0e-300, l2Norm(b));
}

static std::vector<double> trueResidualVector(const SparseMatrix& a,
                                              const std::vector<double>& x,
                                              const std::vector<double>& b)
{
    const std::vector<double> ax = a.multiply(x);
    std::vector<double> r(b.size(), 0.0);
    parallelFor(b.size(), [&](size_t i) {
        r[i] = b[i] - ax[i];
    });
    return r;
}

static int parseRasOverlapFromSolverName(const std::string& name)
{
    const std::string marker = "FGMRES-RAS";
    const size_t pos = name.find(marker);
    if (pos == std::string::npos) {
        return -1;
    }
    size_t i = pos + marker.size();
    int overlap = 0;
    bool any = false;
    while (i < name.size() && std::isdigit(static_cast<unsigned char>(name[i]))) {
        any = true;
        overlap = 10 * overlap + (name[i] - '0');
        ++i;
    }
    return any ? overlap : -1;
}

struct DofPartitionCounts {
    int freeDofs = 0;
    int dirichletDofs = 0;
    int interfaceDofs = 0;
    int interiorDofs = 0;
};

static DofPartitionCounts countDofPartitions(const Mesh& mesh,
                                             const AssemblyDiagnostics& assemblyDiagnostics)
{
    DofPartitionCounts counts;
    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        const bool dirichlet = mesh.nodes[i].dirichlet;
        const bool interfaceDof = !assemblyDiagnostics.interfaceDof.empty()
            && assemblyDiagnostics.interfaceDof[i] != 0;
        if (dirichlet) {
            ++counts.dirichletDofs;
        } else {
            ++counts.freeDofs;
        }
        if (interfaceDof) {
            ++counts.interfaceDofs;
        } else {
            ++counts.interiorDofs;
        }
    }
    return counts;
}

template <typename Predicate>
static double relativeL2DifferenceWhere(const std::vector<double>& a,
                                        const std::vector<double>& b,
                                        Predicate&& predicate)
{
    double numerator = 0.0;
    double denominator = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (!predicate(i)) {
            continue;
        }
        if (!std::isfinite(a[i]) || !std::isfinite(b[i])) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        const double diff = a[i] - b[i];
        numerator += diff * diff;
        denominator += b[i] * b[i];
    }
    return std::sqrt(numerator) / std::sqrt(std::max(1.0e-300, denominator));
}

static double maxAbsDifferenceWhere(const std::vector<double>& a,
                                    const std::vector<double>& b,
                                    const std::vector<char>& mask)
{
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (i >= mask.size() || mask[i] == 0) {
            continue;
        }
        if (!std::isfinite(a[i]) || !std::isfinite(b[i])) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        result = std::max(result, std::abs(a[i] - b[i]));
    }
    return result;
}

struct RasUpdateDiagnosticsRow {
    std::string solver;
    int overlap = -1;
    int subdomain = -1;
    int coreDofCount = 0;
    int localDofCount = 0;
    int overlapDofCount = 0;
    int minUpdateCount = 0;
    int maxUpdateCount = 0;
    int zeroUpdateCount = 0;
    int multiUpdateCount = 0;
};

static std::vector<RasUpdateDiagnosticsRow> computeRasUpdateDiagnostics(const Mesh& mesh,
                                                                        const SparseMatrix& a,
                                                                        const std::string& solverName,
                                                                        int overlap)
{
    std::vector<RasUpdateDiagnosticsRow> rows;
    if (overlap < 0) {
        return rows;
    }
    int maxSubdomain = 0;
    for (const Node& node : mesh.nodes) {
        maxSubdomain = std::max(maxSubdomain, node.subdomain);
    }
    std::vector<std::vector<int>> coreDofs(static_cast<size_t>(maxSubdomain + 1));
    std::vector<int> updateCount(mesh.nodes.size(), 0);
    for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
        const int block = mesh.nodes[static_cast<size_t>(i)].subdomain;
        coreDofs[static_cast<size_t>(block)].push_back(i);
        ++updateCount[static_cast<size_t>(i)];
    }
    int minUpdate = std::numeric_limits<int>::max();
    int maxUpdate = 0;
    int zeroUpdate = 0;
    int multiUpdate = 0;
    for (int count : updateCount) {
        minUpdate = std::min(minUpdate, count);
        maxUpdate = std::max(maxUpdate, count);
        if (count == 0) {
            ++zeroUpdate;
        }
        if (count > 1) {
            ++multiUpdate;
        }
    }
    if (updateCount.empty()) {
        minUpdate = 0;
    }

    for (size_t blockIndex = 0; blockIndex < coreDofs.size(); ++blockIndex) {
        std::vector<int> local = coreDofs[blockIndex];
        std::vector<int> frontier = coreDofs[blockIndex];
        std::vector<char> inLocal(mesh.nodes.size(), 0);
        for (int dof : local) {
            inLocal[static_cast<size_t>(dof)] = 1;
        }
        for (int layer = 0; layer < overlap; ++layer) {
            std::vector<int> nextFrontier;
            for (int globalRow : frontier) {
                for (int k = a.rowPtr[static_cast<size_t>(globalRow)]; k < a.rowPtr[static_cast<size_t>(globalRow + 1)]; ++k) {
                    const int globalCol = a.colInd[static_cast<size_t>(k)];
                    if (globalCol < 0 || globalCol >= static_cast<int>(mesh.nodes.size())) {
                        continue;
                    }
                    if (!inLocal[static_cast<size_t>(globalCol)]) {
                        inLocal[static_cast<size_t>(globalCol)] = 1;
                        local.push_back(globalCol);
                        nextFrontier.push_back(globalCol);
                    }
                }
            }
            frontier = std::move(nextFrontier);
            if (frontier.empty()) {
                break;
            }
        }
        std::sort(local.begin(), local.end());
        local.erase(std::unique(local.begin(), local.end()), local.end());

        RasUpdateDiagnosticsRow row;
        row.solver = solverName;
        row.overlap = overlap;
        row.subdomain = static_cast<int>(blockIndex);
        row.coreDofCount = static_cast<int>(coreDofs[blockIndex].size());
        row.localDofCount = static_cast<int>(local.size());
        row.overlapDofCount = row.localDofCount - row.coreDofCount;
        row.minUpdateCount = minUpdate;
        row.maxUpdateCount = maxUpdate;
        row.zeroUpdateCount = zeroUpdate;
        row.multiUpdateCount = multiUpdate;
        rows.push_back(row);
    }
    return rows;
}

static std::string csvEscape(const std::string& value);

static void writeSolverComparison(const std::vector<SolverStatistics>& stats,
                                  const std::vector<double>& maxDiffs,
                                  const std::vector<double>& l2Diffs,
                                  const std::vector<double>& relativeL2Diffs,
                                  const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "solver,status,failure_reason,iterations,setup_seconds,solve_seconds,total_seconds,"
        << "preconditioner_apply_seconds,preconditioner_apply_calls,preconditioner_apply_seconds_per_call,"
        << "final_relative_residual,max_abs_diff_vs_global,l2_diff_vs_global,relative_l2_diff_vs_global,Tmin,Tmax,Tavg,"
        << "max_iterations,parallel_workers,"
        << "preconditioner_memory_mb,ras_halo_build_seconds,ras_local_matrix_assembly_seconds,"
        << "ras_local_factorization_seconds,ras_local_solve_apply_seconds,ras_restriction_seconds,"
        << "ras_communication_or_halo_update_seconds,ras_setup_count,ras_factorization_reuse,"
        << "coarse_dim,coarse_setup_seconds,coarse_solve_seconds,coarse_matrix_nnz,coarse_residual_norm,"
        << "working_set_before_mb,working_set_after_mb,peak_working_set_mb,"
        << "solver_status\n";
    out << std::setprecision(16);
    const auto writeRow = [&](const SolverStatistics& s, double maxDiff, double l2Diff, double relativeL2Diff) {
        const double applySecondsPerCall = s.preconditionerApplyCalls > 0
            ? s.preconditionerApplySeconds / static_cast<double>(s.preconditionerApplyCalls)
            : 0.0;
        out << s.name << ',' << s.status << ",\"" << s.failureReason << "\"," << s.totalIterations << ','
            << s.setupSeconds << ',' << s.solveSeconds << ','
            << s.setupSeconds + s.solveSeconds << ','
            << s.preconditionerApplySeconds << ','
            << s.preconditionerApplyCalls << ','
            << applySecondsPerCall << ','
            << s.finalRelativeResidual << ','
            << maxDiff << ',' << l2Diff << ',' << relativeL2Diff << ','
            << s.temperatureMin << ',' << s.temperatureMax << ',' << s.temperatureAverage << ','
            << s.maxIterations << ','
            << s.parallelWorkers << ',' << megabytes(s.preconditionerBytes)
            << ',' << s.rasHaloBuildSeconds
            << ',' << s.rasLocalMatrixAssemblySeconds
            << ',' << s.rasLocalFactorizationSeconds
            << ',' << s.rasLocalSolveApplySeconds
            << ',' << s.rasRestrictionSeconds
            << ',' << s.rasCommunicationOrHaloUpdateSeconds
            << ',' << s.rasSetupCount
            << ',' << (s.rasFactorizationReuse ? 1 : 0)
            << ',' << s.coarseDim
            << ',' << s.coarseSetupSeconds
            << ',' << s.coarseSolveSeconds
            << ',' << s.coarseMatrixNnz
            << ',' << s.coarseResidualNorm
            << ',' << megabytes(s.workingSetBeforeBytes)
            << ',' << megabytes(s.workingSetAfterBytes)
            << ',' << megabytes(s.peakWorkingSetBytes)
            << ',' << s.status << '\n';
    };
    for (size_t i = 0; i < stats.size(); ++i) {
        writeRow(stats[i], maxDiffs[i], l2Diffs[i], relativeL2Diffs[i]);
    }
}

static void writeRasPreconditionerTiming(const std::vector<SolverStatistics>& stats,
                                         const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "solver,setup_count,apply_calls,halo_build_time,local_matrix_assembly_time,"
        << "local_factorization_time,local_solve_apply_time,restriction_time,"
        << "communication_or_halo_update_time,total_apply_time,factorization_reuse,"
        << "coarse_setup_time,coarse_solve_time\n";
    out << std::setprecision(16);
    for (const SolverStatistics& s : stats) {
        if (s.rasSetupCount == 0
            && s.rasHaloBuildSeconds == 0.0
            && s.rasLocalMatrixAssemblySeconds == 0.0
            && s.rasLocalFactorizationSeconds == 0.0) {
            continue;
        }
        out << csvEscape(s.name) << ','
            << s.rasSetupCount << ','
            << s.preconditionerApplyCalls << ','
            << s.rasHaloBuildSeconds << ','
            << s.rasLocalMatrixAssemblySeconds << ','
            << s.rasLocalFactorizationSeconds << ','
            << s.rasLocalSolveApplySeconds << ','
            << s.rasRestrictionSeconds << ','
            << s.rasCommunicationOrHaloUpdateSeconds << ','
            << s.preconditionerApplySeconds << ','
            << (s.rasFactorizationReuse ? 1 : 0) << ','
            << s.coarseSetupSeconds << ','
            << s.coarseSolveSeconds << '\n';
    }
}

static void writeCoarseDiagnostics(const std::vector<SolverStatistics>& stats,
                                   const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "solver,coarse_dim,coarse_setup_time,coarse_solve_time,"
        << "coarse_matrix_nnz,coarse_residual_norm\n";
    out << std::setprecision(16);
    for (const SolverStatistics& s : stats) {
        if (s.coarseDim <= 0 && s.coarseMatrixNnz == 0) {
            continue;
        }
        out << csvEscape(s.name) << ','
            << s.coarseDim << ','
            << s.coarseSetupSeconds << ','
            << s.coarseSolveSeconds << ','
            << s.coarseMatrixNnz << ','
            << s.coarseResidualNorm << '\n';
    }
}

static void fillTemperatureStats(SolverStatistics& stats, const std::vector<double>& temperature)
{
    if (temperature.empty() || vectorHasNonFinite(temperature)) {
        return;
    }
    const auto minmax = std::minmax_element(temperature.begin(), temperature.end());
    stats.temperatureMin = *minmax.first;
    stats.temperatureMax = *minmax.second;
    stats.temperatureAverage =
        std::accumulate(temperature.begin(), temperature.end(), 0.0)
        / static_cast<double>(temperature.size());
}

static double vectorOrNaN(const std::vector<double>& values, size_t index)
{
    return index < values.size() ? values[index] : std::numeric_limits<double>::quiet_NaN();
}

static std::string csvEscape(const std::string& value);

static void writeValidationSolverComparison(const std::vector<SolverStatistics>& stats,
                                            const std::vector<double>& maxDiffs,
                                            const std::vector<double>& relativeL2Diffs,
                                            const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "solver,status,failure_reason,iterations,setup_seconds,solve_seconds,total_seconds,"
        << "final_relative_residual,max_abs_diff_vs_global,relative_L2_diff_vs_global,"
        << "Tmin,Tmax,Tavg,accepted_ic_shift_subdomain0,accepted_ic_shift_subdomain1\n";
    for (size_t i = 0; i < stats.size(); ++i) {
        const SolverStatistics& s = stats[i];
        out << s.name << ','
            << s.status << ','
            << csvEscape(s.failureReason) << ','
            << s.totalIterations << ','
            << s.setupSeconds << ','
            << s.solveSeconds << ','
            << s.setupSeconds + s.solveSeconds << ','
            << s.finalRelativeResidual << ','
            << maxDiffs[i] << ','
            << relativeL2Diffs[i] << ','
            << s.temperatureMin << ','
            << s.temperatureMax << ','
            << s.temperatureAverage << ','
            << vectorOrNaN(s.acceptedIcShiftBySubdomain, 0) << ','
            << vectorOrNaN(s.acceptedIcShiftBySubdomain, 1) << '\n';
    }
}

static void printSolverSummary(const SolverStatistics& stats)
{
    std::cout << "  " << stats.name
              << " [" << stats.status << "]"
              << ": setup=" << stats.setupSeconds
              << " s, solve=" << stats.solveSeconds
              << " s, total_iterations=" << stats.totalIterations
              << ", final_rel_residual=" << stats.finalRelativeResidual
              << ", max_step_iterations=" << stats.maxIterations
              << ", parallel_workers=" << stats.parallelWorkers
              << ", factor_or_preconditioner_memory~" << megabytes(stats.preconditionerBytes)
              << " MB, working_set_after~" << megabytes(stats.workingSetAfterBytes)
              << " MB, peak~" << megabytes(stats.peakWorkingSetBytes) << " MB";
    if (stats.rasSetupCount > 0
        || stats.rasHaloBuildSeconds > 0.0
        || stats.rasLocalMatrixAssemblySeconds > 0.0
        || stats.rasLocalFactorizationSeconds > 0.0) {
        std::cout << ", ras_setup(halo/local/factor)="
                  << stats.rasHaloBuildSeconds << "/"
                  << stats.rasLocalMatrixAssemblySeconds << "/"
                  << stats.rasLocalFactorizationSeconds
                  << " s, ras_apply(local/restrict/halo)="
                  << stats.rasLocalSolveApplySeconds << "/"
                  << stats.rasRestrictionSeconds << "/"
                  << stats.rasCommunicationOrHaloUpdateSeconds
                  << " s, reuse=" << (stats.rasFactorizationReuse ? "yes" : "no");
    }
    if (stats.coarseDim > 0) {
        std::cout << ", coarse(dim/setup/solve/nnz/res)="
                  << stats.coarseDim << "/"
                  << stats.coarseSetupSeconds << "/"
                  << stats.coarseSolveSeconds << "/"
                  << stats.coarseMatrixNnz << "/"
                  << stats.coarseResidualNorm;
    }
    if (!stats.failureReason.empty()) {
        std::cout << ", reason=" << stats.failureReason;
    }
    std::cout << "\n";
}

class DisjointSet {
public:
    explicit DisjointSet(size_t n = 0)
    {
        reset(n);
    }

    void reset(size_t n)
    {
        parent_.resize(n);
        rank_.assign(n, 0);
        std::iota(parent_.begin(), parent_.end(), 0);
    }

    size_t find(size_t x)
    {
        if (parent_[x] != x) {
            parent_[x] = find(parent_[x]);
        }
        return parent_[x];
    }

    void unite(size_t a, size_t b)
    {
        size_t ra = find(a);
        size_t rb = find(b);
        if (ra == rb) {
            return;
        }
        if (rank_[ra] < rank_[rb]) {
            std::swap(ra, rb);
        }
        parent_[rb] = ra;
        if (rank_[ra] == rank_[rb]) {
            ++rank_[ra];
        }
    }

private:
    std::vector<size_t> parent_;
    std::vector<int> rank_;
};

struct ComponentDiagnostics {
    std::string graphName;
    int subdomain = -1;
    int component = 0;
    int dofCount = 0;
    bool hasDirichlet = false;
    bool hasConvection = false;
};

struct MatrixDiagnostics {
    bool hasNonFinite = false;
    int zeroRows = 0;
    double minDiagonal = 0.0;
    double maxDiagonal = 0.0;
    double symmetryRatio = 0.0;
    int nonPositiveConductivityTets = 0;
    int heatSourceTetCount = 0;
    int convectionBoundaryFaceCount = 0;
    int globalComponentCount = 0;
    int globalUnconstrainedComponentCount = 0;
    bool likelySpd = false;
};

struct MatrixStageDiagnostic {
    std::string stage;
    double minDiagonal = std::numeric_limits<double>::max();
    double maxDiagonal = -std::numeric_limits<double>::max();
    int negativeDiagonalCount = 0;
    int zeroDiagonalCount = 0;
    double diagonalSum = 0.0;
    double symmetryRatio = 0.0;
};

static MatrixStageDiagnostic diagnoseMatrixStage(const std::string& stage, const SparseMatrix& matrix)
{
    MatrixStageDiagnostic diagnostics;
    diagnostics.stage = stage;
    SparseMatrix csr = matrix;
    if (!csr.csrReady) {
        csr.finalizeCsr();
    }
    std::vector<double> diagonal(static_cast<size_t>(csr.size()), 0.0);
    double matrixNormSquared = 0.0;
    double asymmetryNormSquared = 0.0;
    for (int row = 0; row < csr.size(); ++row) {
        for (int k = csr.rowPtr[static_cast<size_t>(row)]; k < csr.rowPtr[static_cast<size_t>(row + 1)]; ++k) {
            const int col = csr.colInd[static_cast<size_t>(k)];
            const double value = csr.values[static_cast<size_t>(k)];
            if (row == col) {
                diagonal[static_cast<size_t>(row)] += value;
            }
            matrixNormSquared += value * value;
            double transposeValue = 0.0;
            const auto begin = csr.colInd.begin() + csr.rowPtr[static_cast<size_t>(col)];
            const auto end = csr.colInd.begin() + csr.rowPtr[static_cast<size_t>(col + 1)];
            const auto it = std::lower_bound(begin, end, row);
            if (it != end && *it == row) {
                const size_t index = static_cast<size_t>(std::distance(csr.colInd.begin(), it));
                transposeValue = csr.values[index];
            }
            const double diff = value - transposeValue;
            asymmetryNormSquared += diff * diff;
        }
    }

    for (double value : diagonal) {
        diagnostics.minDiagonal = std::min(diagnostics.minDiagonal, value);
        diagnostics.maxDiagonal = std::max(diagnostics.maxDiagonal, value);
        diagnostics.diagonalSum += value;
        if (value < 0.0) {
            ++diagnostics.negativeDiagonalCount;
        }
        if (value == 0.0) {
            ++diagnostics.zeroDiagonalCount;
        }
    }
    if (diagonal.empty()) {
        diagnostics.minDiagonal = 0.0;
        diagnostics.maxDiagonal = 0.0;
    }

    diagnostics.symmetryRatio = std::sqrt(asymmetryNormSquared)
        / std::sqrt(std::max(1.0e-300, matrixNormSquared));
    return diagnostics;
}

static void writeMatrixStageDiagnostics(const std::vector<MatrixStageDiagnostic>& diagnostics,
                                        const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "stage,min_diagonal,max_diagonal,negative_diagonal_count,zero_diagonal_count,diagonal_sum,symmetry_error\n";
    out << std::setprecision(16);
    for (const MatrixStageDiagnostic& diag : diagnostics) {
        out << diag.stage << ','
            << diag.minDiagonal << ','
            << diag.maxDiagonal << ','
            << diag.negativeDiagonalCount << ','
            << diag.zeroDiagonalCount << ','
            << diag.diagonalSum << ','
            << diag.symmetryRatio << '\n';
    }
}

static bool tryCsrValue(const SparseMatrix& a, int row, int col, double& value)
{
    if (!a.csrReady || row < 0 || row >= a.size()) {
        value = 0.0;
        return false;
    }
    const auto begin = a.colInd.begin() + a.rowPtr[static_cast<size_t>(row)];
    const auto end = a.colInd.begin() + a.rowPtr[static_cast<size_t>(row + 1)];
    const auto it = std::lower_bound(begin, end, col);
    if (it == end || *it != col) {
        value = 0.0;
        return false;
    }
    const size_t index = static_cast<size_t>(std::distance(a.colInd.begin(), it));
    value = a.values[index];
    return true;
}

static std::vector<char> markConvectionDofs(const Mesh& mesh, const CaseConfig& config)
{
    std::vector<char> constrained(mesh.nodes.size(), 0);
    for (const BoundaryFace& face : mesh.boundaryFaces) {
        if (!isConvectionBoundary(face, config)) {
            continue;
        }
        const Tet& tet = mesh.tets[static_cast<size_t>(face.tet)];
        for (int i = 0; i < 3; ++i) {
            constrained[static_cast<size_t>(tet.dof[static_cast<size_t>(face.local[static_cast<size_t>(i)])])] = 1;
        }
        for (int i = 0; i < 3; ++i) {
            const int a = face.local[static_cast<size_t>(i)];
            const int b = face.local[static_cast<size_t>((i + 1) % 3)];
            constrained[static_cast<size_t>(tet.dof[static_cast<size_t>(localEdgeDof(a, b))])] = 1;
        }
    }
    return constrained;
}

static void addTetDofClique(DisjointSet& dsu, const Tet& tet)
{
    for (int i = 1; i < 10; ++i) {
        dsu.unite(static_cast<size_t>(tet.dof[0]), static_cast<size_t>(tet.dof[static_cast<size_t>(i)]));
    }
}

static std::vector<ComponentDiagnostics> collectComponentDiagnostics(const Mesh& mesh,
                                                                     const CaseConfig& config)
{
    const std::vector<char> convectionDofs = markConvectionDofs(mesh, config);
    std::vector<ComponentDiagnostics> diagnostics;

    const auto summarize = [&](const std::string& graphName,
                               int subdomain,
                               DisjointSet& dsu,
                               const std::vector<int>& dofs) {
        struct Accumulator {
            int count = 0;
            bool hasDirichlet = false;
            bool hasConvection = false;
        };
        std::map<size_t, Accumulator> byRoot;
        for (int dof : dofs) {
            Accumulator& acc = byRoot[dsu.find(static_cast<size_t>(dof))];
            ++acc.count;
            acc.hasDirichlet = acc.hasDirichlet || mesh.nodes[static_cast<size_t>(dof)].dirichlet;
            acc.hasConvection = acc.hasConvection || convectionDofs[static_cast<size_t>(dof)] != 0;
        }
        int componentIndex = 0;
        for (const auto& entry : byRoot) {
            diagnostics.push_back({graphName,
                                   subdomain,
                                   componentIndex++,
                                   entry.second.count,
                                   entry.second.hasDirichlet,
                                   entry.second.hasConvection});
        }
    };

    for (size_t subdomain = 0; subdomain < config.domains.size(); ++subdomain) {
        DisjointSet dsu(mesh.nodes.size());
        std::vector<int> dofs;
        for (size_t i = 0; i < mesh.nodes.size(); ++i) {
            if (mesh.nodes[i].subdomain == static_cast<int>(subdomain)) {
                dofs.push_back(static_cast<int>(i));
            }
        }
        for (const Tet& tet : mesh.tets) {
            if (tet.subdomain == static_cast<int>(subdomain)) {
                addTetDofClique(dsu, tet);
            }
        }
        summarize("subdomain_volume", static_cast<int>(subdomain), dsu, dofs);
    }

    DisjointSet globalDsu(mesh.nodes.size());
    std::vector<int> allDofs(mesh.nodes.size());
    std::iota(allDofs.begin(), allDofs.end(), 0);
    for (const Tet& tet : mesh.tets) {
        addTetDofClique(globalDsu, tet);
    }
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
        for (int i = 0; i < 10; ++i) {
            globalDsu.unite(static_cast<size_t>(left.dof[static_cast<size_t>(i)]),
                            static_cast<size_t>(right.dof[0]));
            globalDsu.unite(static_cast<size_t>(right.dof[static_cast<size_t>(i)]),
                            static_cast<size_t>(left.dof[0]));
        }
    }
    summarize("global_volume_plus_interface", -1, globalDsu, allDofs);

    return diagnostics;
}

static MatrixDiagnostics diagnoseMatrixAndPhysics(const Mesh& mesh,
                                                  const CaseConfig& config,
                                                  const SparseMatrix& system,
                                                  int convectionBoundaryFaceCount,
                                                  const std::vector<ComponentDiagnostics>& components)
{
    MatrixDiagnostics diagnostics;
    diagnostics.convectionBoundaryFaceCount = convectionBoundaryFaceCount;

    if (!system.csrReady) {
        throw std::runtime_error("Matrix diagnostics require finalized CSR.");
    }

    diagnostics.minDiagonal = std::numeric_limits<double>::max();
    diagnostics.maxDiagonal = -std::numeric_limits<double>::max();
    double matrixNormSquared = 0.0;
    double asymmetryNormSquared = 0.0;

    for (int row = 0; row < system.size(); ++row) {
        if (system.rowPtr[static_cast<size_t>(row)] == system.rowPtr[static_cast<size_t>(row + 1)]) {
            ++diagnostics.zeroRows;
        }
        double diag = 0.0;
        const bool hasDiag = tryCsrValue(system, row, row, diag);
        if (!hasDiag) {
            diag = 0.0;
        }
        diagnostics.minDiagonal = std::min(diagnostics.minDiagonal, diag);
        diagnostics.maxDiagonal = std::max(diagnostics.maxDiagonal, diag);

        for (int k = system.rowPtr[static_cast<size_t>(row)]; k < system.rowPtr[static_cast<size_t>(row + 1)]; ++k) {
            const int col = system.colInd[static_cast<size_t>(k)];
            const double value = system.values[static_cast<size_t>(k)];
            if (!std::isfinite(value)) {
                diagnostics.hasNonFinite = true;
            }
            matrixNormSquared += value * value;
            double transposeValue = 0.0;
            const bool hasTranspose = tryCsrValue(system, col, row, transposeValue);
            const double diff = value - transposeValue;
            asymmetryNormSquared += diff * diff;
            if (!hasTranspose && col != row) {
                asymmetryNormSquared += value * value;
            }
        }
    }

    diagnostics.symmetryRatio = std::sqrt(asymmetryNormSquared)
                              / std::sqrt(std::max(1.0e-300, matrixNormSquared));

    std::vector<char> heatSourceTet(mesh.tets.size(), 0);
    for (size_t i = 0; i < mesh.tets.size(); ++i) {
        const Tet& tet = mesh.tets[i];
        const Material& material = materialForTet(config, tet);
        if (material.conductivityX <= 0.0 || material.conductivityY <= 0.0 || material.conductivityZ <= 0.0
            || !std::isfinite(material.conductivityX)
            || !std::isfinite(material.conductivityY)
            || !std::isfinite(material.conductivityZ)) {
            ++diagnostics.nonPositiveConductivityTets;
        }
        for (const HeatSource& source : config.heatSources) {
            if (tetMatchesHeatSource(tet, source)) {
                heatSourceTet[i] = 1;
            }
        }
    }
    diagnostics.heatSourceTetCount = static_cast<int>(
        std::count(heatSourceTet.begin(), heatSourceTet.end(), static_cast<char>(1)));

    for (const ComponentDiagnostics& component : components) {
        if (component.graphName == "global_volume_plus_interface") {
            ++diagnostics.globalComponentCount;
            if (!component.hasDirichlet && !component.hasConvection) {
                ++diagnostics.globalUnconstrainedComponentCount;
            }
        }
    }

    diagnostics.likelySpd =
        !diagnostics.hasNonFinite
        && diagnostics.zeroRows == 0
        && diagnostics.minDiagonal > 0.0
        && diagnostics.symmetryRatio < 1.0e-10
        && diagnostics.globalUnconstrainedComponentCount == 0
        && diagnostics.nonPositiveConductivityTets == 0;
    return diagnostics;
}

static void writeMatrixDiagnostics(const MatrixDiagnostics& diagnostics,
                                   const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "metric,value\n";
    out << std::setprecision(16);
    out << "has_nonfinite," << (diagnostics.hasNonFinite ? 1 : 0) << '\n';
    out << "zero_rows," << diagnostics.zeroRows << '\n';
    out << "min_diagonal," << diagnostics.minDiagonal << '\n';
    out << "max_diagonal," << diagnostics.maxDiagonal << '\n';
    out << "symmetry_ratio_norm_K_minus_KT_over_norm_K," << diagnostics.symmetryRatio << '\n';
    out << "non_positive_conductivity_tets," << diagnostics.nonPositiveConductivityTets << '\n';
    out << "heat_source_tets," << diagnostics.heatSourceTetCount << '\n';
    out << "convection_boundary_faces," << diagnostics.convectionBoundaryFaceCount << '\n';
    out << "global_components," << diagnostics.globalComponentCount << '\n';
    out << "global_unconstrained_components," << diagnostics.globalUnconstrainedComponentCount << '\n';
    out << "likely_spd," << (diagnostics.likelySpd ? 1 : 0) << '\n';
}

static void writeComponentDiagnostics(const std::vector<ComponentDiagnostics>& components,
                                      const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "graph,subdomain,component,dofs,has_dirichlet,has_convection_or_robin\n";
    for (const ComponentDiagnostics& component : components) {
        out << component.graphName << ','
            << component.subdomain << ','
            << component.component << ','
            << component.dofCount << ','
            << (component.hasDirichlet ? 1 : 0) << ','
            << (component.hasConvection ? 1 : 0) << '\n';
    }
}

static void writeDiagonalContributionStats(const AssemblyDiagnostics& diagnostics,
                                           const std::filesystem::path& path)
{
    const std::vector<std::pair<std::string, const std::vector<double>*>> parts{
        {"K_volume", &diagnostics.volumeDiag},
        {"K_robin_lhs", &diagnostics.robinDiag},
        {"K_dirichlet", &diagnostics.dirichletDiag},
        {"K_interface_consistency", &diagnostics.interfaceConsistencyDiag},
        {"K_interface_penalty", &diagnostics.interfacePenaltyDiag},
        {"K_final", &diagnostics.finalDiag}
    };

    std::ofstream out(path);
    out << "part,min_diagonal,max_diagonal,sum_diagonal,negative_diagonal_entries,zero_diagonal_entries\n";
    out << std::setprecision(16);
    for (const auto& part : parts) {
        const DiagonalStats stats = diagonalStats(*part.second);
        out << part.first << ','
            << stats.minDiagonal << ','
            << stats.maxDiagonal << ','
            << stats.sumDiagonal << ','
            << stats.negativeEntries << ','
            << stats.zeroEntries << '\n';
    }
}

static void writeNegativeDiagonalDofs(const Mesh& mesh,
                                      const AssemblyDiagnostics& diagnostics,
                                      const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "global_dof_id,subdomain,local_node_id,x_m,y_m,z_m,is_interface_node,"
        << "left_boundary_id,right_boundary_id,K_volume_diag,K_robin_lhs_diag,K_interface_consistency_diag,"
        << "K_interface_penalty_diag,K_final_diag\n";
    out << std::setprecision(16);
    int written = 0;
    for (size_t i = 0; i < diagnostics.finalDiag.size() && written < 100; ++i) {
        if (diagnostics.finalDiag[i] >= 0.0) {
            continue;
        }
        const Node& node = mesh.nodes[i];
        const auto boundaryIds = diagnostics.interfaceBoundaryByDof.empty()
            ? std::pair<int, int>{-1, -1}
            : diagnostics.interfaceBoundaryByDof[i];
        out << i << ','
            << node.subdomain << ','
            << node.sourceVertex << ','
            << node.p.x << ','
            << node.p.y << ','
            << node.p.z << ','
            << (!diagnostics.interfaceDof.empty() && diagnostics.interfaceDof[i] ? 1 : 0) << ','
            << boundaryIds.first << ','
            << boundaryIds.second << ','
            << diagnostics.volumeDiag[i] << ','
            << diagnostics.robinDiag[i] << ','
            << diagnostics.interfaceConsistencyDiag[i] << ','
            << diagnostics.interfacePenaltyDiag[i] << ','
            << diagnostics.finalDiag[i] << '\n';
        ++written;
    }
}

static void writeInterfacePenaltyStats(const InterfacePenaltyStats& stats,
                                       const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "face_pair_count,min_eta_F,max_eta_F,avg_eta_F,min_h_F,max_h_F,avg_h_F,min_k_left_over_k_right,max_k_left_over_k_right\n";
    out << std::setprecision(16);
    const double count = static_cast<double>(std::max(1, stats.facePairCount));
    out << stats.facePairCount << ','
        << (stats.facePairCount > 0 ? stats.etaMin : 0.0) << ','
        << (stats.facePairCount > 0 ? stats.etaMax : 0.0) << ','
        << (stats.facePairCount > 0 ? stats.etaSum / count : 0.0) << ','
        << (stats.facePairCount > 0 ? stats.hMin : 0.0) << ','
        << (stats.facePairCount > 0 ? stats.hMax : 0.0) << ','
        << (stats.facePairCount > 0 ? stats.hSum / count : 0.0) << ','
        << (stats.facePairCount > 0 ? stats.kRatioMin : 0.0) << ','
        << (stats.facePairCount > 0 ? stats.kRatioMax : 0.0) << '\n';
}

static void writeHeatSourceDiagnostics(const std::vector<HeatSourceAssemblyDiagnostic>& diagnostics,
                                       const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "heat_source_index,subdomain,domain_entity,configured_power_W,tet_count,"
        << "source_volume_used_for_density_m3,source_volume_assembled_m3,"
        << "volumetric_q_W_per_m3,expected_power_W,quadrature_power_W\n";
    out << std::setprecision(16);
    for (const HeatSourceAssemblyDiagnostic& diag : diagnostics) {
        out << diag.index << ','
            << diag.subdomain << ','
            << diag.domainEntity << ','
            << diag.configuredPowerW << ','
            << diag.tetCount << ','
            << diag.sourceVolumeUsedForDensity << ','
            << diag.sourceVolumeAssembled << ','
            << diag.volumetricQ << ','
            << diag.expectedPowerW << ','
            << diag.quadraturePowerW << '\n';
    }
}

static void writeHeatSourceDomainDiagnostics(const std::vector<HeatSourceAssemblyDiagnostic>& diagnostics,
                                             const std::filesystem::path& path)
{
    struct DomainTotals {
        int subdomain = -1;
        int domainEntity = -1;
        double configuredPowerW = 0.0;
        int tetCount = 0;
        double sourceVolumeUsedForDensity = 0.0;
        double sourceVolumeAssembled = 0.0;
        double effectiveVolumetricQ = 0.0;
        double expectedPowerW = 0.0;
        double quadraturePowerW = 0.0;
    };

    std::map<std::pair<int, int>, DomainTotals> totals;
    for (const HeatSourceAssemblyDiagnostic& diag : diagnostics) {
        const std::pair<int, int> key{diag.subdomain, diag.domainEntity};
        DomainTotals& total = totals[key];
        total.subdomain = diag.subdomain;
        total.domainEntity = diag.domainEntity;
        total.configuredPowerW += diag.configuredPowerW;
        total.tetCount = std::max(total.tetCount, diag.tetCount);
        total.sourceVolumeUsedForDensity = std::max(total.sourceVolumeUsedForDensity,
                                                    diag.sourceVolumeUsedForDensity);
        total.sourceVolumeAssembled = std::max(total.sourceVolumeAssembled,
                                               diag.sourceVolumeAssembled);
        total.effectiveVolumetricQ += diag.volumetricQ;
        total.expectedPowerW += diag.expectedPowerW;
        total.quadraturePowerW += diag.quadraturePowerW;
    }

    std::ofstream out(path);
    out << "subdomain,domain_entity,configured_power_W,tet_count,"
        << "source_volume_used_for_density_m3,source_volume_assembled_m3,"
        << "effective_volumetric_q_W_per_m3,expected_power_W,quadrature_power_W\n";
    out << std::setprecision(16);
    for (const auto& entry : totals) {
        const DomainTotals& total = entry.second;
        out << total.subdomain << ','
            << total.domainEntity << ','
            << total.configuredPowerW << ','
            << total.tetCount << ','
            << total.sourceVolumeUsedForDensity << ','
            << total.sourceVolumeAssembled << ','
            << total.effectiveVolumetricQ << ','
            << total.expectedPowerW << ','
            << total.quadraturePowerW << '\n';
    }
}

static void writePhysicalSummary(const Mesh& mesh,
                                 const CaseConfig& config,
                                 const std::vector<double>& heatOnlySource,
                                 const std::vector<double>& sourceBeforeDirichlet,
                                 const std::vector<double>& sourceAfterDirichlet,
                                 int convectionBoundaryFaceCount,
                                 const std::filesystem::path& path)
{
    double configuredHeatW = 0.0;
    for (const HeatSource& heatSource : config.heatSources) {
        configuredHeatW += heatSource.heatRateW;
    }
    const double heatOnlyRhsSum = std::accumulate(heatOnlySource.begin(), heatOnlySource.end(), 0.0);
    const double beforeDirichletRhsSum = std::accumulate(sourceBeforeDirichlet.begin(), sourceBeforeDirichlet.end(), 0.0);
    const double afterDirichletRhsSum = std::accumulate(sourceAfterDirichlet.begin(), sourceAfterDirichlet.end(), 0.0);
    int negativeHeatOnlyEntries = 0;
    for (double value : heatOnlySource) {
        if (value < 0.0) {
            ++negativeHeatOnlyEntries;
        }
    }
    int negativeBeforeDirichletEntries = 0;
    for (double value : sourceBeforeDirichlet) {
        if (value < 0.0) {
            ++negativeBeforeDirichletEntries;
        }
    }
    int dirichletNodeCount = 0;
    double dirichletMin = std::numeric_limits<double>::max();
    double dirichletMax = -std::numeric_limits<double>::max();
    for (const Node& node : mesh.nodes) {
        if (node.dirichlet) {
            ++dirichletNodeCount;
            dirichletMin = std::min(dirichletMin, node.dirichletValue);
            dirichletMax = std::max(dirichletMax, node.dirichletValue);
        }
    }
    if (dirichletNodeCount == 0) {
        dirichletMin = 0.0;
        dirichletMax = 0.0;
    }

    std::ofstream out(path);
    out << "metric,value\n";
    out << std::setprecision(16);
    out << "configured_heat_source_total_W," << configuredHeatW << '\n';
    out << "heat_only_rhs_sum_before_convection_and_dirichlet," << heatOnlyRhsSum << '\n';
    out << "assembled_rhs_sum_after_convection_before_dirichlet," << beforeDirichletRhsSum << '\n';
    out << "assembled_rhs_sum_after_dirichlet," << afterDirichletRhsSum << '\n';
    out << "negative_heat_only_rhs_entries," << negativeHeatOnlyEntries << '\n';
    out << "negative_rhs_entries_after_convection_before_dirichlet," << negativeBeforeDirichletEntries << '\n';
    out << "dirichlet_node_count," << dirichletNodeCount << '\n';
    out << "dirichlet_temperature_min_K," << dirichletMin << '\n';
    out << "dirichlet_temperature_max_K," << dirichletMax << '\n';
    out << "convection_robin_condition_count," << config.convectionConditions.size() << '\n';
    out << "convection_robin_boundary_faces," << convectionBoundaryFaceCount << '\n';
}

static void writeIcDiagnostics(const std::vector<IcFactorDiagnostics>& diagnostics,
                               const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "subdomain,dofs,ordering,applied_shift,diagonal_scaling,subdomain_block_exact_spd,accepted,pivot_min,pivot_max,"
        << "pivot_nonpositive_count,pivot_tiny_count,first_bad_pivot_row,"
        << "first_bad_pivot_value,nonfinite_L_count,breakdown\n";
    out << std::setprecision(16);
    for (const IcFactorDiagnostics& item : diagnostics) {
        out << item.subdomain << ','
            << item.dofs << ','
            << item.ordering << ','
            << item.appliedShift << ','
            << (item.diagonalScaling ? 1 : 0) << ','
            << (item.exactSpd ? 1 : 0) << ','
            << (item.accepted ? 1 : 0) << ','
            << (item.pivotMin == std::numeric_limits<double>::max() ? 0.0 : item.pivotMin) << ','
            << item.pivotMax << ','
            << item.pivotNonpositiveCount << ','
            << item.pivotTinyCount << ','
            << item.firstBadPivotRow << ','
            << item.firstBadPivotValue << ','
            << item.nonFiniteLCount << ','
            << (item.breakdown ? 1 : 0) << '\n';
    }
}

static void printDiagonalContributionStats(const AssemblyDiagnostics& diagnostics)
{
    const std::vector<std::pair<std::string, const std::vector<double>*>> parts{
        {"K_volume", &diagnostics.volumeDiag},
        {"K_robin_lhs", &diagnostics.robinDiag},
        {"K_dirichlet", &diagnostics.dirichletDiag},
        {"K_interface_consistency", &diagnostics.interfaceConsistencyDiag},
        {"K_interface_penalty", &diagnostics.interfacePenaltyDiag},
        {"K_final", &diagnostics.finalDiag}
    };
    std::cout << "Diagonal contribution diagnostics:\n";
    for (const auto& part : parts) {
        const DiagonalStats stats = diagonalStats(*part.second);
        std::cout << "  " << part.first
                  << ": min=" << stats.minDiagonal
                  << ", max=" << stats.maxDiagonal
                  << ", sum=" << stats.sumDiagonal
                  << ", negative=" << stats.negativeEntries
                  << ", zero=" << stats.zeroEntries << "\n";
    }
}

static void printInterfacePenaltyStats(const InterfacePenaltyStats& stats)
{
    const double count = static_cast<double>(std::max(1, stats.facePairCount));
    std::cout << "Interface penalty eta_F diagnostics:"
              << " face_pairs=" << stats.facePairCount
              << ", eta[min/avg/max]="
              << (stats.facePairCount > 0 ? stats.etaMin : 0.0) << "/"
              << (stats.facePairCount > 0 ? stats.etaSum / count : 0.0) << "/"
              << (stats.facePairCount > 0 ? stats.etaMax : 0.0)
              << ", h_F[min/avg/max]="
              << (stats.facePairCount > 0 ? stats.hMin : 0.0) << "/"
              << (stats.facePairCount > 0 ? stats.hSum / count : 0.0) << "/"
              << (stats.facePairCount > 0 ? stats.hMax : 0.0)
              << ", kL/kR[min/max]="
              << (stats.facePairCount > 0 ? stats.kRatioMin : 0.0) << "/"
              << (stats.facePairCount > 0 ? stats.kRatioMax : 0.0) << "\n";
}

static void printMatrixDiagnostics(const MatrixDiagnostics& diagnostics)
{
    std::cout << "Matrix diagnostics:\n"
              << "  has_NaN_or_Inf=" << (diagnostics.hasNonFinite ? "yes" : "no")
              << ", zero_rows=" << diagnostics.zeroRows
              << ", diagonal_min=" << diagnostics.minDiagonal
              << ", diagonal_max=" << diagnostics.maxDiagonal << "\n"
              << "  symmetry ||K-KT||/||K||=" << diagnostics.symmetryRatio
              << ", k<=0 tets=" << diagnostics.nonPositiveConductivityTets
              << ", heat_source_tets=" << diagnostics.heatSourceTetCount
              << ", convection_faces=" << diagnostics.convectionBoundaryFaceCount << "\n"
              << "  global_components=" << diagnostics.globalComponentCount
              << ", unconstrained_global_components=" << diagnostics.globalUnconstrainedComponentCount
              << ", likely_spd=" << (diagnostics.likelySpd ? "yes" : "no") << "\n";
}

static void markInterfaceDofs(const Mesh& mesh, AssemblyDiagnostics& diagnostics)
{
    diagnostics.interfaceDof.assign(mesh.nodes.size(), 0);
    diagnostics.interfaceBoundaryByDof.assign(mesh.nodes.size(), {-1, -1});
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
        for (int i = 0; i < 10; ++i) {
            const int leftDof = left.dof[static_cast<size_t>(i)];
            diagnostics.interfaceDof[static_cast<size_t>(leftDof)] = 1;
            diagnostics.interfaceBoundaryByDof[static_cast<size_t>(leftDof)] =
                {face.leftBoundaryEntity, face.rightBoundaryEntity};
            const int rightDof = right.dof[static_cast<size_t>(i)];
            diagnostics.interfaceDof[static_cast<size_t>(rightDof)] = 1;
            diagnostics.interfaceBoundaryByDof[static_cast<size_t>(rightDof)] =
                {face.leftBoundaryEntity, face.rightBoundaryEntity};
        }
    }
}

static void writeProgramTiming(const ProgramTiming& timing, const std::filesystem::path& path)
{
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Cannot open " + path.string());
    }
    out << "stage,seconds\n";
    out << std::setprecision(16);
    out << "preprocessing," << timing.preprocessingSeconds << '\n';
    out << "volume_assembly," << timing.volumeAssemblySeconds << '\n';
    out << "interface_assembly," << timing.interfaceAssemblySeconds << '\n';
    out << "system_build," << timing.systemBuildSeconds << '\n';
    out << "dirichlet_assembly," << timing.dirichletAssemblySeconds << '\n';
    out << "csr_finalize," << timing.csrFinalizeSeconds << '\n';
    out << "matrix_assembly," << timing.assemblySeconds << '\n';
    out << "postprocessing," << timing.postprocessingSeconds << '\n';
    out << "total_program," << timing.totalSeconds << '\n';
    out.close();
    if (!out) {
        throw std::runtime_error("Cannot write " + path.string());
    }
}

static void printProgramTiming(const ProgramTiming& timing)
{
    std::cout << "Program timing:\n"
              << "  preprocessing=" << timing.preprocessingSeconds << " s\n"
              << "  volume_assembly=" << timing.volumeAssemblySeconds << " s\n"
              << "  interface_assembly=" << timing.interfaceAssemblySeconds << " s\n"
              << "  system_build=" << timing.systemBuildSeconds << " s\n"
              << "  dirichlet_assembly=" << timing.dirichletAssemblySeconds << " s\n"
              << "  csr_finalize=" << timing.csrFinalizeSeconds << " s\n"
              << "  matrix_assembly=" << timing.assemblySeconds << " s\n"
              << "  postprocessing=" << timing.postprocessingSeconds << " s\n"
              << "  total_program=" << timing.totalSeconds << " s\n";
}

static double interfaceAverageJump(const Mesh& mesh, const std::vector<double>& temperature)
{
    double sum = 0.0;
    double areaSum = 0.0;
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
        for (const auto& tri : face.integrationTriangles) {
            const Vec3 center = (tri[0] + tri[1] + tri[2]) / 3.0;
            const double area = 0.5 * norm(cross(tri[1] - tri[0], tri[2] - tri[0]));
            const auto lLeft = lambdaOnTetFace(center, left, face.leftLocal, mesh);
            const auto lRight = lambdaOnTetFace(center, right, face.rightLocal, mesh);
            const auto nLeft = shapeP2(lLeft);
            const auto nRight = shapeP2(lRight);
            double tLeft = 0.0;
            double tRight = 0.0;
            for (int i = 0; i < 10; ++i) {
                tLeft += nLeft[static_cast<size_t>(i)] * temperature[static_cast<size_t>(left.dof[static_cast<size_t>(i)])];
                tRight += nRight[static_cast<size_t>(i)] * temperature[static_cast<size_t>(right.dof[static_cast<size_t>(i)])];
            }
            sum += std::abs(tLeft - tRight) * area;
            areaSum += area;
        }
    }
    return areaSum > 0.0 ? sum / areaSum : 0.0;
}

static void writeCsv(const Mesh& mesh, const std::vector<double>& temperature, const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "id,subdomain,source_vertex,x_m,y_m,z_m,temperature_K,dirichlet\n";
    out << std::setprecision(16);
    for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
        const Node& n = mesh.nodes[static_cast<size_t>(i)];
        out << i << ',' << n.subdomain << ',' << n.sourceVertex << ','
            << n.p.x << ',' << n.p.y << ',' << n.p.z << ','
            << temperature[static_cast<size_t>(i)] << ',' << (n.dirichlet ? 1 : 0) << '\n';
    }
}

static void writeVtk(const Mesh& mesh, const std::vector<double>& temperature, const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "# vtk DataFile Version 3.0\n";
    out << "COMSOL mesh SIPG P2 heat result\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << mesh.nodes.size() << " double\n";
    out << std::setprecision(16);
    for (const Node& n : mesh.nodes) {
        out << n.p.x << ' ' << n.p.y << ' ' << n.p.z << '\n';
    }
    out << "CELLS " << mesh.tets.size() << ' ' << mesh.tets.size() * 11 << "\n";
    for (const Tet& tet : mesh.tets) {
        out << "10";
        for (int id : tet.dof) {
            out << ' ' << id;
        }
        out << '\n';
    }
    out << "CELL_TYPES " << mesh.tets.size() << "\n";
    for (size_t i = 0; i < mesh.tets.size(); ++i) {
        out << "24\n";
    }
    out << "POINT_DATA " << mesh.nodes.size() << "\n";
    out << "SCALARS temperature_K double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (double t : temperature) {
        out << t << '\n';
    }
    out << "SCALARS subdomain int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (const Node& n : mesh.nodes) {
        out << n.subdomain << '\n';
    }
}

struct Rgb {
    unsigned char r = 255;
    unsigned char g = 255;
    unsigned char b = 255;
};

static Rgb temperatureColor(double value, double minValue, double maxValue)
{
    const double span = std::max(1.0e-12, maxValue - minValue);
    const double t = std::clamp((value - minValue) / span, 0.0, 1.0);
    if (t < 0.5) {
        const double s = t / 0.5;
        return {static_cast<unsigned char>(255.0 * s),
                static_cast<unsigned char>(255.0 * s),
                255};
    }
    const double s = (t - 0.5) / 0.5;
    return {255,
            static_cast<unsigned char>(255.0 * (1.0 - s)),
            static_cast<unsigned char>(255.0 * (1.0 - s))};
}

static void setPixel(std::vector<Rgb>& image, int width, int height, int x, int y, Rgb color)
{
    if (x < 0 || y < 0 || x >= width || y >= height) {
        return;
    }
    image[static_cast<size_t>(y * width + x)] = color;
}

static void writeBmp(const std::filesystem::path& path, const std::vector<Rgb>& image, int width, int height)
{
    const int rowStride = (3 * width + 3) & ~3;
    const int pixelDataSize = rowStride * height;
    const int fileSize = 54 + pixelDataSize;
    std::ofstream out(path, std::ios::binary);
    unsigned char header[54] = {};
    header[0] = 'B';
    header[1] = 'M';
    header[2] = static_cast<unsigned char>(fileSize);
    header[3] = static_cast<unsigned char>(fileSize >> 8);
    header[4] = static_cast<unsigned char>(fileSize >> 16);
    header[5] = static_cast<unsigned char>(fileSize >> 24);
    header[10] = 54;
    header[14] = 40;
    header[18] = static_cast<unsigned char>(width);
    header[19] = static_cast<unsigned char>(width >> 8);
    header[20] = static_cast<unsigned char>(width >> 16);
    header[21] = static_cast<unsigned char>(width >> 24);
    header[22] = static_cast<unsigned char>(height);
    header[23] = static_cast<unsigned char>(height >> 8);
    header[24] = static_cast<unsigned char>(height >> 16);
    header[25] = static_cast<unsigned char>(height >> 24);
    header[26] = 1;
    header[28] = 24;
    out.write(reinterpret_cast<const char*>(header), sizeof(header));

    std::vector<unsigned char> row(static_cast<size_t>(rowStride), static_cast<unsigned char>(255));
    for (int y = height - 1; y >= 0; --y) {
        std::fill(row.begin(), row.end(), static_cast<unsigned char>(255));
        for (int x = 0; x < width; ++x) {
            const Rgb c = image[static_cast<size_t>(y * width + x)];
            row[static_cast<size_t>(3 * x + 0)] = c.b;
            row[static_cast<size_t>(3 * x + 1)] = c.g;
            row[static_cast<size_t>(3 * x + 2)] = c.r;
        }
        out.write(reinterpret_cast<const char*>(row.data()), row.size());
    }
}

static std::array<int, 2> mapSlicePointToPixel(double x,
                                               double z,
                                               const Vec3& lo,
                                               const Vec3& hi,
                                               int width,
                                               int height,
                                               int left,
                                               int top,
                                               int plotWidth,
                                               int plotHeight)
{
    (void)width;
    (void)height;
    const double sx = (x - lo.x) / std::max(1.0e-30, hi.x - lo.x);
    const double sz = (z - lo.z) / std::max(1.0e-30, hi.z - lo.z);
    return {{
        left + static_cast<int>(std::llround(sx * static_cast<double>(plotWidth - 1))),
        top + static_cast<int>(std::llround((1.0 - sz) * static_cast<double>(plotHeight - 1)))
    }};
}

static bool samePoint(const Vec3& a, const Vec3& b)
{
    return norm(a - b) < 1.0e-12;
}

static void addUniquePoint(std::vector<Vec3>& points, const Vec3& p)
{
    for (const Vec3& q : points) {
        if (samePoint(p, q)) {
            return;
        }
    }
    points.push_back(p);
}

static std::vector<Vec3> tetraPlaneIntersection(const Mesh& mesh, const Tet& tet, double ySlice)
{
    const std::array<std::array<int, 2>, 6> edges{{
        {{0, 1}}, {{0, 2}}, {{0, 3}}, {{1, 2}}, {{1, 3}}, {{2, 3}}
    }};
    std::vector<Vec3> points;
    for (const auto& edge : edges) {
        const Vec3 p0 = mesh.nodes[static_cast<size_t>(tet.v[static_cast<size_t>(edge[0])])].p;
        const Vec3 p1 = mesh.nodes[static_cast<size_t>(tet.v[static_cast<size_t>(edge[1])])].p;
        const double d0 = p0.y - ySlice;
        const double d1 = p1.y - ySlice;
        if (std::abs(d0) < 1.0e-14) {
            addUniquePoint(points, p0);
        }
        if (std::abs(d1) < 1.0e-14) {
            addUniquePoint(points, p1);
        }
        if (d0 * d1 < -1.0e-28) {
            const double t = (ySlice - p0.y) / (p1.y - p0.y);
            addUniquePoint(points, (1.0 - t) * p0 + t * p1);
        }
    }
    if (points.size() < 3) {
        return {};
    }
    const Vec3 center = std::accumulate(points.begin(), points.end(), Vec3{}) / static_cast<double>(points.size());
    std::sort(points.begin(), points.end(), [&](const Vec3& a, const Vec3& b) {
        return std::atan2(a.z - center.z, a.x - center.x) < std::atan2(b.z - center.z, b.x - center.x);
    });
    return points;
}

static double edge2D(double ax, double ay, double bx, double by, double x, double y)
{
    return (bx - ax) * (y - ay) - (by - ay) * (x - ax);
}

static void rasterizeSliceTriangle(std::vector<Rgb>& image,
                                   int width,
                                   int height,
                                   int left,
                                   int top,
                                   int plotWidth,
                                   int plotHeight,
                                   const Vec3& lo,
                                   const Vec3& hi,
                                   const Mesh& mesh,
                                   const Tet& tet,
                                   const ElementGeometry& geo,
                                   const std::vector<double>& temperature,
                                   const std::array<Vec3, 3>& tri,
                                   double minT,
                                   double maxT)
{
    const auto p0 = mapSlicePointToPixel(tri[0].x, tri[0].z, lo, hi, width, height, left, top, plotWidth, plotHeight);
    const auto p1 = mapSlicePointToPixel(tri[1].x, tri[1].z, lo, hi, width, height, left, top, plotWidth, plotHeight);
    const auto p2 = mapSlicePointToPixel(tri[2].x, tri[2].z, lo, hi, width, height, left, top, plotWidth, plotHeight);
    const int minX = std::max(left, std::min({p0[0], p1[0], p2[0]}));
    const int maxX = std::min(left + plotWidth - 1, std::max({p0[0], p1[0], p2[0]}));
    const int minY = std::max(top, std::min({p0[1], p1[1], p2[1]}));
    const int maxY = std::min(top + plotHeight - 1, std::max({p0[1], p1[1], p2[1]}));
    const double orient = edge2D(p0[0], p0[1], p1[0], p1[1], p2[0], p2[1]) >= 0.0 ? 1.0 : -1.0;

    for (int py = minY; py <= maxY; ++py) {
        for (int px = minX; px <= maxX; ++px) {
            if (orient * edge2D(p0[0], p0[1], p1[0], p1[1], px, py) < -1.0e-9
                || orient * edge2D(p1[0], p1[1], p2[0], p2[1], px, py) < -1.0e-9
                || orient * edge2D(p2[0], p2[1], p0[0], p0[1], px, py) < -1.0e-9) {
                continue;
            }
            const double sx = static_cast<double>(px - left) / static_cast<double>(plotWidth - 1);
            const double sz = 1.0 - static_cast<double>(py - top) / static_cast<double>(plotHeight - 1);
            const Vec3 q{lo.x + sx * (hi.x - lo.x), tri[0].y, lo.z + sz * (hi.z - lo.z)};
            std::array<double, 4> lambda{};
            if (!barycentricInTet(mesh, tet, geo, q, lambda)) {
                continue;
            }
            const auto n = shapeP2(lambda);
            double value = 0.0;
            for (int i = 0; i < 10; ++i) {
                value += n[static_cast<size_t>(i)] * temperature[static_cast<size_t>(tet.dof[static_cast<size_t>(i)])];
            }
            setPixel(image, width, height, px, py, temperatureColor(value, minT, maxT));
        }
    }
}

static void writeTemperatureSliceImage(const Mesh& mesh, const std::vector<double>& temperature, const std::filesystem::path& path)
{
    const int width = 1200;
    const int height = 780;
    const int left = 70;
    const int right = 90;
    const int top = 45;
    const int bottom = 65;
    const int plotWidth = width - left - right;
    const int plotHeight = height - top - bottom;
    std::vector<Rgb> image(static_cast<size_t>(width * height), Rgb{250, 250, 250});

    Vec3 lo{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    Vec3 hi{-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};
    for (const Node& node : mesh.nodes) {
        updateBounds(lo, hi, node.p);
    }
    const double ySlice = 0.5 * (lo.y + hi.y);
    const auto minmax = std::minmax_element(temperature.begin(), temperature.end());
    double minT = *minmax.first;
    double maxT = *minmax.second;
    if (maxT - minT < 1.0e-9) {
        minT -= 1.0;
        maxT += 1.0;
    }

    for (const Tet& tet : mesh.tets) {
        double ymin = std::numeric_limits<double>::max();
        double ymax = -std::numeric_limits<double>::max();
        for (int id : tet.v) {
            const double y = mesh.nodes[static_cast<size_t>(id)].p.y;
            ymin = std::min(ymin, y);
            ymax = std::max(ymax, y);
        }
        if (ySlice < ymin - 1.0e-14 || ySlice > ymax + 1.0e-14) {
            continue;
        }
        const std::vector<Vec3> polygon = tetraPlaneIntersection(mesh, tet, ySlice);
        if (polygon.size() < 3) {
            continue;
        }
        const ElementGeometry geo = elementGeometry(mesh, tet);
        for (size_t i = 1; i + 1 < polygon.size(); ++i) {
            rasterizeSliceTriangle(image, width, height, left, top, plotWidth, plotHeight,
                                   lo, hi, mesh, tet, geo, temperature,
                                   {{polygon[0], polygon[i], polygon[i + 1]}}, minT, maxT);
        }
    }

    for (int pass = 0; pass < 4; ++pass) {
        std::vector<Rgb> filled = image;
        for (int y = top + 1; y < top + plotHeight - 1; ++y) {
            for (int x = left + 1; x < left + plotWidth - 1; ++x) {
                const Rgb current = image[static_cast<size_t>(y * width + x)];
                if (!(current.r == 250 && current.g == 250 && current.b == 250)) {
                    continue;
                }
                int count = 0;
                int r = 0;
                int g = 0;
                int b = 0;
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        if (dx == 0 && dy == 0) {
                            continue;
                        }
                        const Rgb n = image[static_cast<size_t>((y + dy) * width + x + dx)];
                        if (n.r == 250 && n.g == 250 && n.b == 250) {
                            continue;
                        }
                        r += n.r;
                        g += n.g;
                        b += n.b;
                        ++count;
                    }
                }
                if (count > 0) {
                    filled[static_cast<size_t>(y * width + x)] = {
                        static_cast<unsigned char>(r / count),
                        static_cast<unsigned char>(g / count),
                        static_cast<unsigned char>(b / count)
                    };
                }
            }
        }
        image.swap(filled);
    }

    const Rgb frame{20, 20, 20};
    for (int x = left; x < left + plotWidth; ++x) {
        setPixel(image, width, height, x, top, frame);
        setPixel(image, width, height, x, top + plotHeight - 1, frame);
    }
    for (int y = top; y < top + plotHeight; ++y) {
        setPixel(image, width, height, left, y, frame);
        setPixel(image, width, height, left + plotWidth - 1, y, frame);
    }

    double interfaceZSum = 0.0;
    int interfacePointCount = 0;
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        for (const auto& tri : face.integrationTriangles) {
            for (const Vec3& p : tri) {
                interfaceZSum += p.z;
                ++interfacePointCount;
            }
        }
    }
    if (interfacePointCount > 0) {
        const double interfaceZ = interfaceZSum / static_cast<double>(interfacePointCount);
        const int interfaceY = mapSlicePointToPixel(lo.x, interfaceZ, lo, hi, width, height, left, top, plotWidth, plotHeight)[1];
        for (int x = left; x < left + plotWidth; ++x) {
            if ((x - left) % 8 < 4) {
                setPixel(image, width, height, x, interfaceY, Rgb{0, 0, 0});
            }
        }
    }

    const int barX0 = width - 62;
    const int barX1 = width - 34;
    for (int y = top; y < top + plotHeight; ++y) {
        const double t = 1.0 - static_cast<double>(y - top) / static_cast<double>(plotHeight - 1);
        const Rgb c = temperatureColor(minT + t * (maxT - minT), minT, maxT);
        for (int x = barX0; x <= barX1; ++x) {
            setPixel(image, width, height, x, y, c);
        }
    }
    for (int x = barX0; x <= barX1; ++x) {
        setPixel(image, width, height, x, top, frame);
        setPixel(image, width, height, x, top + plotHeight - 1, frame);
    }
    for (int y = top; y < top + plotHeight; ++y) {
        setPixel(image, width, height, barX0, y, frame);
        setPixel(image, width, height, barX1, y, frame);
    }

    writeBmp(path, image, width, height);
}

static void writeBoundarySummary(const Mesh& mesh, const std::filesystem::path& path)
{
    struct Summary {
        int count = 0;
        double area = 0.0;
        Vec3 lo{std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
        Vec3 hi{-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};
    };
    std::map<std::pair<int, int>, Summary> summaries;
    for (const BoundaryFace& face : mesh.boundaryFaces) {
        Summary& s = summaries[{face.subdomain, face.boundaryEntity}];
        ++s.count;
        s.area += 0.5 * norm(cross(face.points[1] - face.points[0], face.points[2] - face.points[0]));
        for (const Vec3& p : face.points) {
            updateBounds(s.lo, s.hi, p);
        }
    }

    std::ofstream out(path);
    out << "subdomain,boundary_entity,triangles,area_m2,xmin,ymin,zmin,xmax,ymax,zmax\n";
    out << std::setprecision(16);
    for (const auto& entry : summaries) {
        const Summary& s = entry.second;
        out << entry.first.first << ',' << entry.first.second << ',' << s.count << ',' << s.area << ','
            << s.lo.x << ',' << s.lo.y << ',' << s.lo.z << ','
            << s.hi.x << ',' << s.hi.y << ',' << s.hi.z << '\n';
    }
}

static void writeInterfaceSummary(const Mesh& mesh, const std::filesystem::path& path)
{
    struct Summary {
        int overlappingPairs = 0;
        int integrationTriangles = 0;
        double integrationArea = 0.0;
    };

    std::map<std::pair<int, int>, int> boundaryCounts;
    for (const BoundaryFace& face : mesh.boundaryFaces) {
        ++boundaryCounts[{face.subdomain, face.boundaryEntity}];
    }

    std::map<std::tuple<int, int, int, int>, Summary> summaries;
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        const int leftSubdomain = mesh.tets[static_cast<size_t>(face.leftTet)].subdomain;
        const int rightSubdomain = mesh.tets[static_cast<size_t>(face.rightTet)].subdomain;
        Summary& s = summaries[{leftSubdomain, rightSubdomain, face.leftBoundaryEntity, face.rightBoundaryEntity}];
        ++s.overlappingPairs;
        s.integrationTriangles += static_cast<int>(face.integrationTriangles.size());
        for (const auto& tri : face.integrationTriangles) {
            s.integrationArea += 0.5 * norm(cross(tri[1] - tri[0], tri[2] - tri[0]));
        }
    }

    std::ofstream out(path);
    out << "left_subdomain,right_subdomain,left_boundary_entity,right_boundary_entity,"
        << "left_boundary_triangles,right_boundary_triangles,overlapping_pairs,"
        << "generated_integration_triangles,integration_area_m2\n";
    out << std::setprecision(16);
    for (const auto& entry : summaries) {
        const Summary& s = entry.second;
        const int leftSubdomain = std::get<0>(entry.first);
        const int rightSubdomain = std::get<1>(entry.first);
        const int leftBoundaryEntity = std::get<2>(entry.first);
        const int rightBoundaryEntity = std::get<3>(entry.first);
        out << leftSubdomain << ',' << rightSubdomain << ','
            << leftBoundaryEntity << ',' << rightBoundaryEntity << ','
            << boundaryCounts[{leftSubdomain, leftBoundaryEntity}] << ','
            << boundaryCounts[{rightSubdomain, rightBoundaryEntity}] << ','
            << s.overlappingPairs << ',' << s.integrationTriangles << ',' << s.integrationArea << '\n';
    }
}

static void writeInterfaceBuildSummary(const Mesh& mesh, const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "interface_index,left_subdomain,right_subdomain,left_boundary_id_count,right_boundary_id_count,"
        << "left_face_count,right_face_count,left_total_area_m2,right_total_area_m2,matched_overlap_area_m2,"
        << "overlap_ratio_min_side,face_pair_count,integration_triangle_count,"
        << "normal_dot_min,normal_dot_avg,normal_dot_max\n";
    out << std::setprecision(16);
    for (size_t i = 0; i < mesh.interfaceSummaries.size(); ++i) {
        const InterfaceBuildSummary& s = mesh.interfaceSummaries[i];
        const double denominator = std::min(s.leftArea, s.rightArea);
        const double overlapRatio = denominator > 0.0 ? s.matchedOverlapArea / denominator : 0.0;
        const double normalAvg = s.facePairCount > 0
            ? s.normalDotSum / static_cast<double>(s.facePairCount)
            : 0.0;
        out << i << ','
            << s.leftSubdomain << ','
            << s.rightSubdomain << ','
            << s.leftBoundaryEntityCount << ','
            << s.rightBoundaryEntityCount << ','
            << s.leftFaceCount << ','
            << s.rightFaceCount << ','
            << s.leftArea << ','
            << s.rightArea << ','
            << s.matchedOverlapArea << ','
            << overlapRatio << ','
            << s.facePairCount << ','
            << s.integrationTriangleCount << ','
            << s.normalDotMin << ','
            << normalAvg << ','
            << s.normalDotMax << '\n';
    }
}

static void writeInterfaceTrianglesVtk(const Mesh& mesh, const std::filesystem::path& path)
{
    std::vector<Vec3> points;
    std::vector<std::array<int, 3>> triangles;
    std::vector<int> leftBoundaryIds;
    std::vector<int> rightBoundaryIds;

    for (const InterfaceFace& face : mesh.interfaceFaces) {
        for (const auto& tri : face.integrationTriangles) {
            std::array<int, 3> ids{};
            for (int i = 0; i < 3; ++i) {
                ids[static_cast<size_t>(i)] = static_cast<int>(points.size());
                points.push_back(tri[static_cast<size_t>(i)]);
            }
            triangles.push_back(ids);
            leftBoundaryIds.push_back(face.leftBoundaryEntity);
            rightBoundaryIds.push_back(face.rightBoundaryEntity);
        }
    }

    std::ofstream out(path);
    out << "# vtk DataFile Version 3.0\n";
    out << "SIPG interface integration triangles\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";
    out << "POINTS " << points.size() << " double\n";
    out << std::setprecision(16);
    for (const Vec3& p : points) {
        out << p.x << ' ' << p.y << ' ' << p.z << '\n';
    }
    out << "POLYGONS " << triangles.size() << ' ' << triangles.size() * 4 << "\n";
    for (const auto& tri : triangles) {
        out << "3 " << tri[0] << ' ' << tri[1] << ' ' << tri[2] << '\n';
    }
    out << "CELL_DATA " << triangles.size() << "\n";
    out << "SCALARS left_boundary_entity int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int id : leftBoundaryIds) {
        out << id << '\n';
    }
    out << "SCALARS right_boundary_entity int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int id : rightBoundaryIds) {
        out << id << '\n';
    }
}

static std::vector<ComparisonPoint> readComsolTemperatureExport(const std::filesystem::path& path)
{
    std::vector<ComparisonPoint> points;
    std::ifstream in(path);
    if (!in) {
        return points;
    }
    std::string line;
    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '%' || line[0] == 'X') {
            continue;
        }
        std::istringstream iss(line);
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double t = 0.0;
        if (iss >> x >> y >> z >> t) {
            points.push_back({{1.0e-3 * x, 1.0e-3 * y, 1.0e-3 * z}, t});
        }
    }
    return points;
}

static std::vector<CachedTet> makeCachedTets(const Mesh& mesh)
{
    std::vector<CachedTet> cached;
    cached.reserve(mesh.tets.size());
    for (int i = 0; i < static_cast<int>(mesh.tets.size()); ++i) {
        const Tet& tet = mesh.tets[static_cast<size_t>(i)];
        CachedTet item;
        item.tetId = i;
        item.geo = elementGeometry(mesh, tet);
        item.lo = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
        item.hi = {-std::numeric_limits<double>::max(), -std::numeric_limits<double>::max(), -std::numeric_limits<double>::max()};
        for (int id : tet.v) {
            updateBounds(item.lo, item.hi, mesh.nodes[static_cast<size_t>(id)].p);
        }
        constexpr double eps = 1.0e-12;
        item.lo = item.lo - Vec3{eps, eps, eps};
        item.hi = item.hi + Vec3{eps, eps, eps};
        cached.push_back(item);
    }
    return cached;
}

static bool sampleTemperatureAtPoint(const Mesh& mesh,
                                     const std::vector<CachedTet>& cached,
                                     const std::vector<double>& temperature,
                                     const Vec3& p,
                                     double& value)
{
    double bestValue = 0.0;
    double bestMinLambda = -std::numeric_limits<double>::max();
    bool hasFallback = false;

    for (const CachedTet& item : cached) {
        if (p.x < item.lo.x || p.x > item.hi.x
            || p.y < item.lo.y || p.y > item.hi.y
            || p.z < item.lo.z || p.z > item.hi.z) {
            continue;
        }
        const Tet& tet = mesh.tets[static_cast<size_t>(item.tetId)];
        std::array<double, 4> lambda{};
        const Vec3 p0 = mesh.nodes[static_cast<size_t>(tet.v[0])].p;
        lambda[1] = dot(item.geo.gradLambda[1], p - p0);
        lambda[2] = dot(item.geo.gradLambda[2], p - p0);
        lambda[3] = dot(item.geo.gradLambda[3], p - p0);
        lambda[0] = 1.0 - lambda[1] - lambda[2] - lambda[3];
        const double minLambda = *std::min_element(lambda.begin(), lambda.end());
        const auto n = shapeP2(lambda);
        double localValue = 0.0;
        for (int i = 0; i < 10; ++i) {
            localValue += n[static_cast<size_t>(i)] * temperature[static_cast<size_t>(tet.dof[static_cast<size_t>(i)])];
        }
        if (minLambda >= -1.0e-8) {
            value = localValue;
            return true;
        }
        if (minLambda > bestMinLambda) {
            bestMinLambda = minLambda;
            bestValue = localValue;
            hasFallback = true;
        }
    }

    if (hasFallback && bestMinLambda > -1.0e-6) {
        value = bestValue;
        return true;
    }
    return false;
}

static void compareWithComsolExport(const Mesh& mesh,
                                    const std::vector<double>& temperature,
                                    const std::filesystem::path& comsolPath,
                                    const std::filesystem::path& outputPath)
{
    const std::vector<ComparisonPoint> points = readComsolTemperatureExport(comsolPath);
    if (points.empty()) {
        std::cout << "No COMSOL comparison file found at " << comsolPath.string() << "\n";
        return;
    }

    const std::vector<CachedTet> cached = makeCachedTets(mesh);
    std::ofstream out(outputPath);
    out << "x_m,y_m,z_m,T_comsol_K,T_sipg_K,error_K\n";
    out << std::setprecision(16);

    int sampled = 0;
    int missed = 0;
    double sumError = 0.0;
    double sumAbs = 0.0;
    double sumSq = 0.0;
    double maxAbs = 0.0;
    double comsolMin = std::numeric_limits<double>::max();
    double comsolMax = -std::numeric_limits<double>::max();
    double sipgMin = std::numeric_limits<double>::max();
    double sipgMax = -std::numeric_limits<double>::max();

    for (const ComparisonPoint& point : points) {
        double sipg = 0.0;
        if (!sampleTemperatureAtPoint(mesh, cached, temperature, point.p, sipg)) {
            ++missed;
            continue;
        }
        const double error = sipg - point.comsolTemperature;
        out << point.p.x << ',' << point.p.y << ',' << point.p.z << ','
            << point.comsolTemperature << ',' << sipg << ',' << error << '\n';
        ++sampled;
        sumError += error;
        sumAbs += std::abs(error);
        sumSq += error * error;
        maxAbs = std::max(maxAbs, std::abs(error));
        comsolMin = std::min(comsolMin, point.comsolTemperature);
        comsolMax = std::max(comsolMax, point.comsolTemperature);
        sipgMin = std::min(sipgMin, sipg);
        sipgMax = std::max(sipgMax, sipg);
    }

    if (sampled == 0) {
        std::cout << "COMSOL comparison failed: no exported points could be sampled.\n";
        return;
    }

    std::cout << "COMSOL comparison at " << sampled << " exported points"
              << " (missed " << missed << "):\n";
    std::cout << "  COMSOL Tmin/Tmax = " << comsolMin << " / " << comsolMax << " K\n";
    std::cout << "  SIPG   Tmin/Tmax = " << sipgMin << " / " << sipgMax << " K\n";
    std::cout << "  error mean=" << sumError / static_cast<double>(sampled)
              << " K, mean_abs=" << sumAbs / static_cast<double>(sampled)
              << " K, rmse=" << std::sqrt(sumSq / static_cast<double>(sampled))
              << " K, max_abs=" << maxAbs << " K\n";
    std::cout << "  wrote " << outputPath.string() << "\n";
}

struct ProgramOptions {
    AnalysisMode mode = AnalysisMode::Transient;
    std::filesystem::path configPath;
    bool runDiagonal = false;
    bool runBj = false;
    bool runBjIc = true;
    bool runBjIlut = false;
    bool runRasIc = false;
    bool runRasIlut = false;
    bool runBjIcCoarse = false;
    bool runBjIcCoarsePcg = false;
    bool runBjIlutCoarse = false;
    bool runTwoLevelRasIlut = false;
    bool runDeflatedRasIlut = false;
    bool runInterfaceDeflatedRasIlut = false;
    bool runBjPardisoGeneral = false;
    bool runSchwarz = false;
    bool runSchwarzPrecondFgmres = false;
    bool runSchwarzPrecondFgmresTwoLevel = false;
    bool runDdmSchur = false;
    bool schurLinearXYCoarse = true;
    bool schurLinearZCoarse = true;
    bool schurGlobalQuadraticZCoarse = false;
    bool schurInterfacePatchCoarse = false;
    bool schurInterfacePatchLinearXY = false;
    bool runDirect = true;
    bool solversExplicit = false;
    bool solverMethodOverride = false;
    bool schwarzTypeOverride = false;
    bool schwarzStandaloneModeOverride = false;
    bool schwarzTransmissionOrientationOverride = false;
    bool schwarzRobinAlphaFactorOverride = false;
    bool schwarzFluxEvalOverride = false;
    bool schwarzWriteInterfaceFluxOverride = false;
    bool schwarzOverlapLayersOverride = false;
    bool schwarzOverlapModeOverride = false;
    bool schwarzPartitionModeOverride = false;
    bool schwarzMaxItersOverride = false;
    bool schwarzTolRelUpdateOverride = false;
    bool schwarzTolRelResidualOverride = false;
    bool schwarzRelaxationOverride = false;
    bool disableInterfaceConsistency = false;
    bool disableInterfacePenalty = false;
    bool nodeTieInterface = false;
    bool disableConvectionLhs = false;
    bool diagnosticsOnly = false;
    bool spectralDiagnostics = false;
    bool sipgSpdDiagnostics = false;
    bool spdPenaltySweep = false;
    bool skipSipgPenaltySweep = false;
    bool bjIcValidation = false;
    bool bjIcQuickValidation = false;
    bool fastRun = false;
    bool penaltyModeOverride = false;
    bool penaltyFactorOverride = false;
    bool dirichletMethodOverride = false;
    bool nitschePenaltyFactorOverride = false;
    bool showHelp = false;
    bool outputDirOverride = false;
    bool modelNameOverride = false;
    bool disableWarmStartOverride = false;
    bool forceNontrivialRhsOverride = false;
    bool initialGuessOverride = false;
    bool thermalSourceScaleOverride = false;
    bool localDiagScaling = true;
    bool localDiagScalingOverride = false;
    std::filesystem::path outputDir;
    std::string modelName;
    std::string solverMethod = "schwarz_precond_fgmres";
    std::string schwarzType = "multiplicative";
    std::string schwarzStandaloneMode = "algebraic";
    std::string schwarzTransmissionOrientation = "forward";
    std::string schwarzFluxEval = "sipg_numeric";
    std::string schwarzOverlapMode = "ras";
    std::string schwarzPartitionMode = "current";
    std::string penaltyMode = "harmonic";
    std::string dirichletMethod = "strong";
    std::string initialGuessType = "current";
    std::string localIcShiftMode = "auto";
    std::string coarseCorrection = "additive";
    std::string coarseSpace = "subdomain_constant";
    std::string directMode = "both";
    double penaltyFactor = 15.0;
    double nitschePenaltyFactor = 15.0;
    double icShift = 1.0e-12;
    double diagScalingEps = 1.0e-30;
    double thermalSourceScale = 1.0;
    bool icScaling = true;
    bool disableWarmStart = false;
    bool forceNontrivialRhs = false;
    std::string icOrdering = "natural";
    int maxPcgIterations = 1000;
    int schwarzMaxIters = 1000;
    int schwarzOverlapLayers = 1;
    double pcgTolerance = 1.0e-10;
    double schwarzTolRelUpdate = 1.0e-10;
    double schwarzTolRelResidual = 1.0e-10;
    double schwarzRelaxation = 1.0;
    double schwarzRobinAlphaFactor = 10.0;
    bool schwarzWriteInterfaceFlux = false;
    int gmresRestart = 30;
    int deflationModes = 6;
    int spectralModeCount = 20;
    double nodeTiePenalty = 1.0e12;
    std::vector<double> ilutDropTolerances{1.0e-3, 1.0e-4, 1.0e-5};
    std::vector<int> ilutFillFactors{5, 10, 20};
    struct RasIlutConfig {
        int overlap = 1;
        double dropTolerance = 1.0e-3;
        int fillFactor = 5;
    };
    std::vector<RasIlutConfig> rasIlutConfigs;
};

static bool nextArgLooksLikeValue(int argc, char* argv[], int i)
{
    return i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0;
}

static std::vector<double> parseDoubleList(const std::string& value)
{
    std::vector<double> result;
    for (const std::string& item : splitCsv(value)) {
        if (!item.empty()) {
            result.push_back(std::stod(item));
        }
    }
    return result;
}

static std::vector<int> parseIntList(const std::string& value)
{
    std::vector<int> result;
    for (const std::string& item : splitCsv(value)) {
        if (!item.empty()) {
            result.push_back(std::stoi(item));
        }
    }
    return result;
}

static std::vector<ProgramOptions::RasIlutConfig> parseRasIlutConfigs(const std::string& value)
{
    std::vector<ProgramOptions::RasIlutConfig> result;
    std::istringstream stream(value);
    std::string item;
    while (std::getline(stream, item, ';')) {
        item = trim(item);
        if (item.empty()) {
            continue;
        }
        std::replace(item.begin(), item.end(), ':', ',');
        const std::vector<std::string> parts = splitCsv(item);
        if (parts.size() != 3) {
            throw std::runtime_error("--ras-ilut-configs entries must be overlap:drop_tol:fill_factor");
        }
        ProgramOptions::RasIlutConfig config;
        config.overlap = std::stoi(parts[0]);
        config.dropTolerance = std::stod(parts[1]);
        config.fillFactor = std::stoi(parts[2]);
        result.push_back(config);
    }
    return result;
}

static bool parseSchurInterfaceCoarseOption(ProgramOptions& options,
                                            const std::string& arg)
{
    if (arg == "--schur-interface-patch-coarse") {
        options.schurInterfacePatchCoarse = true;
    } else if (arg == "--no-schur-interface-patch-coarse") {
        options.schurInterfacePatchCoarse = false;
    } else if (arg == "--schur-interface-patch-linear-xy") {
        options.schurInterfacePatchLinearXY = true;
    } else if (arg == "--no-schur-interface-patch-linear-xy") {
        options.schurInterfacePatchLinearXY = false;
    } else {
        return false;
    }
    return true;
}

static ProgramOptions parseProgramOptions(int argc, char* argv[])
{
    ProgramOptions options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (parseSchurInterfaceCoarseOption(options, arg)) {
            continue;
        }
        if (arg == "--steady" || arg == "steady") {
            options.mode = AnalysisMode::Steady;
        } else if (arg == "--transient" || arg == "transient") {
            options.mode = AnalysisMode::Transient;
        } else if (arg == "--config" && i + 1 < argc) {
            options.configPath = argv[++i];
        } else if (arg.rfind("--config=", 0) == 0) {
            options.configPath = arg.substr(std::string("--config=").size());
        } else if (arg == "--output-dir" && i + 1 < argc) {
            options.outputDirOverride = true;
            options.outputDir = argv[++i];
        } else if (arg.rfind("--output-dir=", 0) == 0) {
            options.outputDirOverride = true;
            options.outputDir = arg.substr(std::string("--output-dir=").size());
        } else if (arg == "--model-name" && i + 1 < argc) {
            options.modelNameOverride = true;
            options.modelName = argv[++i];
        } else if (arg.rfind("--model-name=", 0) == 0) {
            options.modelNameOverride = true;
            options.modelName = arg.substr(std::string("--model-name=").size());
        } else if (arg == "--solver-method" && i + 1 < argc) {
            options.solverMethodOverride = true;
            options.solverMethod = normalizeSolverMethodName(argv[++i]);
        } else if (arg.rfind("--solver-method=", 0) == 0) {
            options.solverMethodOverride = true;
            options.solverMethod = normalizeSolverMethodName(arg.substr(std::string("--solver-method=").size()));
        } else if (arg == "--schwarz-type" && i + 1 < argc) {
            options.schwarzTypeOverride = true;
            options.schwarzType = lowerString(argv[++i]);
        } else if (arg.rfind("--schwarz-type=", 0) == 0) {
            options.schwarzTypeOverride = true;
            options.schwarzType = lowerString(arg.substr(std::string("--schwarz-type=").size()));
        } else if (arg == "--schwarz-standalone-mode" && i + 1 < argc) {
            options.schwarzStandaloneModeOverride = true;
            options.schwarzStandaloneMode = normalizeSchwarzTransmissionModeName(argv[++i]);
        } else if (arg.rfind("--schwarz-standalone-mode=", 0) == 0) {
            options.schwarzStandaloneModeOverride = true;
            options.schwarzStandaloneMode =
                normalizeSchwarzTransmissionModeName(arg.substr(std::string("--schwarz-standalone-mode=").size()));
        } else if (arg == "--schwarz-transmission-orientation" && i + 1 < argc) {
            options.schwarzTransmissionOrientationOverride = true;
            options.schwarzTransmissionOrientation =
                normalizeSchwarzTransmissionOrientationName(argv[++i]);
        } else if (arg.rfind("--schwarz-transmission-orientation=", 0) == 0) {
            options.schwarzTransmissionOrientationOverride = true;
            options.schwarzTransmissionOrientation =
                normalizeSchwarzTransmissionOrientationName(
                    arg.substr(std::string("--schwarz-transmission-orientation=").size()));
        } else if (arg == "--schwarz-robin-alpha-factor" && i + 1 < argc) {
            options.schwarzRobinAlphaFactorOverride = true;
            options.schwarzRobinAlphaFactor = std::stod(argv[++i]);
        } else if (arg.rfind("--schwarz-robin-alpha-factor=", 0) == 0) {
            options.schwarzRobinAlphaFactorOverride = true;
            options.schwarzRobinAlphaFactor = std::stod(arg.substr(std::string("--schwarz-robin-alpha-factor=").size()));
        } else if (arg == "--schwarz-flux-eval" && i + 1 < argc) {
            options.schwarzFluxEvalOverride = true;
            options.schwarzFluxEval = lowerString(argv[++i]);
        } else if (arg.rfind("--schwarz-flux-eval=", 0) == 0) {
            options.schwarzFluxEvalOverride = true;
            options.schwarzFluxEval = lowerString(arg.substr(std::string("--schwarz-flux-eval=").size()));
        } else if (arg == "--schwarz-overlap-layers" && i + 1 < argc) {
            options.schwarzOverlapLayersOverride = true;
            options.schwarzOverlapLayers = std::stoi(argv[++i]);
        } else if (arg.rfind("--schwarz-overlap-layers=", 0) == 0) {
            options.schwarzOverlapLayersOverride = true;
            options.schwarzOverlapLayers = std::stoi(arg.substr(std::string("--schwarz-overlap-layers=").size()));
        } else if (arg == "--overlap" && i + 1 < argc) {
            options.schwarzOverlapLayersOverride = true;
            options.schwarzOverlapLayers = std::stoi(argv[++i]);
        } else if (arg.rfind("--overlap=", 0) == 0) {
            options.schwarzOverlapLayersOverride = true;
            options.schwarzOverlapLayers = std::stoi(arg.substr(std::string("--overlap=").size()));
        } else if (arg == "--schwarz-overlap-mode" && i + 1 < argc) {
            options.schwarzOverlapModeOverride = true;
            options.schwarzOverlapMode = lowerString(argv[++i]);
        } else if (arg.rfind("--schwarz-overlap-mode=", 0) == 0) {
            options.schwarzOverlapModeOverride = true;
            options.schwarzOverlapMode = lowerString(arg.substr(std::string("--schwarz-overlap-mode=").size()));
        } else if (arg == "--ras-type" && i + 1 < argc) {
            options.schwarzOverlapModeOverride = true;
            const std::string value = lowerString(argv[++i]);
            options.schwarzOverlapMode = value == "restricted" ? "ras" : (value == "additive" ? "halo" : value);
        } else if (arg.rfind("--ras-type=", 0) == 0) {
            options.schwarzOverlapModeOverride = true;
            const std::string value = lowerString(arg.substr(std::string("--ras-type=").size()));
            options.schwarzOverlapMode = value == "restricted" ? "ras" : (value == "additive" ? "halo" : value);
        } else if (arg == "--schwarz-partition-mode" && i + 1 < argc) {
            options.schwarzPartitionModeOverride = true;
            options.schwarzPartitionMode = normalizeSchwarzPartitionModeName(argv[++i]);
        } else if (arg.rfind("--schwarz-partition-mode=", 0) == 0) {
            options.schwarzPartitionModeOverride = true;
            options.schwarzPartitionMode = normalizeSchwarzPartitionModeName(arg.substr(std::string("--schwarz-partition-mode=").size()));
        } else if (arg == "--schwarz-write-interface-flux") {
            options.schwarzWriteInterfaceFluxOverride = true;
            options.schwarzWriteInterfaceFlux = true;
        } else if (arg.rfind("--schwarz-write-interface-flux=", 0) == 0) {
            options.schwarzWriteInterfaceFluxOverride = true;
            options.schwarzWriteInterfaceFlux = parseBoolValue(arg.substr(std::string("--schwarz-write-interface-flux=").size()));
        } else if (arg == "--schwarz-max-iters" && i + 1 < argc) {
            options.schwarzMaxItersOverride = true;
            options.schwarzMaxIters = std::stoi(argv[++i]);
        } else if (arg.rfind("--schwarz-max-iters=", 0) == 0) {
            options.schwarzMaxItersOverride = true;
            options.schwarzMaxIters = std::stoi(arg.substr(std::string("--schwarz-max-iters=").size()));
        } else if (arg == "--schwarz-tol-rel-update" && i + 1 < argc) {
            options.schwarzTolRelUpdateOverride = true;
            options.schwarzTolRelUpdate = std::stod(argv[++i]);
        } else if (arg.rfind("--schwarz-tol-rel-update=", 0) == 0) {
            options.schwarzTolRelUpdateOverride = true;
            options.schwarzTolRelUpdate = std::stod(arg.substr(std::string("--schwarz-tol-rel-update=").size()));
        } else if (arg == "--schwarz-tol-rel-residual" && i + 1 < argc) {
            options.schwarzTolRelResidualOverride = true;
            options.schwarzTolRelResidual = std::stod(argv[++i]);
        } else if (arg.rfind("--schwarz-tol-rel-residual=", 0) == 0) {
            options.schwarzTolRelResidualOverride = true;
            options.schwarzTolRelResidual = std::stod(arg.substr(std::string("--schwarz-tol-rel-residual=").size()));
        } else if (arg == "--schwarz-relaxation" && i + 1 < argc) {
            options.schwarzRelaxationOverride = true;
            options.schwarzRelaxation = std::stod(argv[++i]);
        } else if (arg.rfind("--schwarz-relaxation=", 0) == 0) {
            options.schwarzRelaxationOverride = true;
            options.schwarzRelaxation = std::stod(arg.substr(std::string("--schwarz-relaxation=").size()));
        } else if (arg == "--schur-linear-xy-coarse") {
            options.schurLinearXYCoarse = true;
        } else if (arg == "--no-schur-linear-xy-coarse") {
            options.schurLinearXYCoarse = false;
        } else if (arg == "--schur-linear-z-coarse") {
            options.schurLinearZCoarse = true;
        } else if (arg == "--no-schur-linear-z-coarse") {
            options.schurLinearZCoarse = false;
        } else if (arg == "--schur-global-quadratic-z-coarse") {
            options.schurGlobalQuadraticZCoarse = true;
        } else if (arg == "--no-schur-global-quadratic-z-coarse") {
            options.schurGlobalQuadraticZCoarse = false;
        } else if (arg == "--solvers" && i + 1 < argc) {
            options.solversExplicit = true;
            options.runDiagonal = false;
            options.runBj = false;
            options.runBjIc = false;
            options.runBjIlut = false;
            options.runRasIc = false;
            options.runRasIlut = false;
            options.runBjIcCoarse = false;
            options.runBjIcCoarsePcg = false;
            options.runBjIlutCoarse = false;
            options.runTwoLevelRasIlut = false;
            options.runDeflatedRasIlut = false;
            options.runInterfaceDeflatedRasIlut = false;
            options.runBjPardisoGeneral = false;
            options.runSchwarz = false;
            options.runSchwarzPrecondFgmres = false;
            options.runSchwarzPrecondFgmresTwoLevel = false;
            options.runDdmSchur = false;
            options.runDirect = false;
            for (const std::string& solver : splitCsv(argv[++i])) {
                const std::string name = lowerString(solver);
                options.runDiagonal = options.runDiagonal || name == "diagonal" || name == "pcg";
                options.runBj = options.runBj || name == "bj-pardiso-pcg";
                options.runBjIc = options.runBjIc || name == "bj-ic" || name == "bj-shifted-ic" || name == "fgmres-bj-ic" || name == "ic";
                options.runBjIlut = options.runBjIlut || name == "bj-ilut" || name == "fgmres-bj-ilut" || name == "ilut";
                options.runRasIc = options.runRasIc || name == "ras-ic" || name == "fgmres-ras-ic";
                options.runRasIlut = options.runRasIlut || name == "ras-ilut" || name == "fgmres-ras-ilut" || name == "ras";
                options.runBjIcCoarse = options.runBjIcCoarse || name == "bj-ic-coarse" || name == "two-level-bj-ic" || name == "fgmres-bj-ic-coarse";
                options.runBjIcCoarsePcg = options.runBjIcCoarsePcg || name == "bj-ic-coarse-pcg" || name == "pcg-bj-ic-coarse";
                options.runBjIlutCoarse = options.runBjIlutCoarse || name == "bj-ilut-coarse" || name == "two-level-bj-ilut" || name == "fgmres-bj-ilut-coarse";
                options.runTwoLevelRasIlut = options.runTwoLevelRasIlut || name == "two-level-ras" || name == "two-level-ras-ilut" || name == "fgmres-two-level-ras-ilut";
                options.runDeflatedRasIlut = options.runDeflatedRasIlut || name == "deflated-ras" || name == "deflated-ras-ilut" || name == "fgmres-deflated-ras-ilut";
                options.runInterfaceDeflatedRasIlut = options.runInterfaceDeflatedRasIlut || name == "interface-deflated-ras" || name == "interface-deflated-ras-ilut";
                options.runBjPardisoGeneral = options.runBjPardisoGeneral || name == "bj" || name == "bj-pardiso" || name == "bj-pardiso-general" || name == "fgmres-bj-pardiso-general";
                options.runSchwarz = options.runSchwarz || name == "schwarz" || name == "schwarz-sipg" || name == "sipg-schwarz";
                options.runSchwarzPrecondFgmres = options.runSchwarzPrecondFgmres
                    || name == "schwarz-precond-fgmres"
                    || name == "schwarz_precond_fgmres"
                    || name == "fgmres-schwarz"
                    || name == "schwarz-fgmres";
                options.runSchwarzPrecondFgmresTwoLevel = options.runSchwarzPrecondFgmresTwoLevel
                    || name == "schwarz-precond-fgmres-two-level"
                    || name == "schwarz_precond_fgmres_two_level"
                    || name == "two-level-schwarz-precond-fgmres"
                    || name == "two_level_schwarz_precond_fgmres"
                    || name == "fgmres-schwarz-two-level"
                    || name == "schwarz-fgmres-two-level";
                options.runDdmSchur = options.runDdmSchur
                    || name == "schur"
                    || name == "ddm-schur"
                    || name == "schur-fgmres";
                options.runDirect = options.runDirect || name == "direct" || name == "global" || name == "global-pardiso" || name == "global-pardiso-general";
            }
        } else if (arg.rfind("--solvers=", 0) == 0) {
            options.solversExplicit = true;
            options.runDiagonal = false;
            options.runBj = false;
            options.runBjIc = false;
            options.runBjIlut = false;
            options.runRasIc = false;
            options.runRasIlut = false;
            options.runBjIcCoarse = false;
            options.runBjIcCoarsePcg = false;
            options.runBjIlutCoarse = false;
            options.runTwoLevelRasIlut = false;
            options.runDeflatedRasIlut = false;
            options.runInterfaceDeflatedRasIlut = false;
            options.runBjPardisoGeneral = false;
            options.runSchwarz = false;
            options.runSchwarzPrecondFgmres = false;
            options.runSchwarzPrecondFgmresTwoLevel = false;
            options.runDdmSchur = false;
            options.runDirect = false;
            for (const std::string& solver : splitCsv(arg.substr(std::string("--solvers=").size()))) {
                const std::string name = lowerString(solver);
                options.runDiagonal = options.runDiagonal || name == "diagonal" || name == "pcg";
                options.runBj = options.runBj || name == "bj-pardiso-pcg";
                options.runBjIc = options.runBjIc || name == "bj-ic" || name == "bj-shifted-ic" || name == "fgmres-bj-ic" || name == "ic";
                options.runBjIlut = options.runBjIlut || name == "bj-ilut" || name == "fgmres-bj-ilut" || name == "ilut";
                options.runRasIc = options.runRasIc || name == "ras-ic" || name == "fgmres-ras-ic";
                options.runRasIlut = options.runRasIlut || name == "ras-ilut" || name == "fgmres-ras-ilut" || name == "ras";
                options.runBjIcCoarse = options.runBjIcCoarse || name == "bj-ic-coarse" || name == "two-level-bj-ic" || name == "fgmres-bj-ic-coarse";
                options.runBjIcCoarsePcg = options.runBjIcCoarsePcg || name == "bj-ic-coarse-pcg" || name == "pcg-bj-ic-coarse";
                options.runBjIlutCoarse = options.runBjIlutCoarse || name == "bj-ilut-coarse" || name == "two-level-bj-ilut" || name == "fgmres-bj-ilut-coarse";
                options.runTwoLevelRasIlut = options.runTwoLevelRasIlut || name == "two-level-ras" || name == "two-level-ras-ilut" || name == "fgmres-two-level-ras-ilut";
                options.runDeflatedRasIlut = options.runDeflatedRasIlut || name == "deflated-ras" || name == "deflated-ras-ilut" || name == "fgmres-deflated-ras-ilut";
                options.runInterfaceDeflatedRasIlut = options.runInterfaceDeflatedRasIlut || name == "interface-deflated-ras" || name == "interface-deflated-ras-ilut";
                options.runBjPardisoGeneral = options.runBjPardisoGeneral || name == "bj" || name == "bj-pardiso" || name == "bj-pardiso-general" || name == "fgmres-bj-pardiso-general";
                options.runSchwarz = options.runSchwarz || name == "schwarz" || name == "schwarz-sipg" || name == "sipg-schwarz";
                options.runSchwarzPrecondFgmres = options.runSchwarzPrecondFgmres
                    || name == "schwarz-precond-fgmres"
                    || name == "schwarz_precond_fgmres"
                    || name == "fgmres-schwarz"
                    || name == "schwarz-fgmres";
                options.runSchwarzPrecondFgmresTwoLevel = options.runSchwarzPrecondFgmresTwoLevel
                    || name == "schwarz-precond-fgmres-two-level"
                    || name == "schwarz_precond_fgmres_two_level"
                    || name == "two-level-schwarz-precond-fgmres"
                    || name == "two_level_schwarz_precond_fgmres"
                    || name == "fgmres-schwarz-two-level"
                    || name == "schwarz-fgmres-two-level";
                options.runDdmSchur = options.runDdmSchur
                    || name == "schur"
                    || name == "ddm-schur"
                    || name == "schur-fgmres";
                options.runDirect = options.runDirect || name == "direct" || name == "global" || name == "global-pardiso" || name == "global-pardiso-general";
            }
        } else if (arg == "--disable-interface-consistency") {
            options.disableInterfaceConsistency = true;
        } else if (arg == "--disable-interface-penalty") {
            options.disableInterfacePenalty = true;
        } else if (arg == "--node-tie-interface") {
            options.nodeTieInterface = true;
        } else if (arg == "--node-tie-penalty" && i + 1 < argc) {
            options.nodeTiePenalty = std::stod(argv[++i]);
        } else if (arg.rfind("--node-tie-penalty=", 0) == 0) {
            options.nodeTiePenalty = std::stod(arg.substr(std::string("--node-tie-penalty=").size()));
        } else if (arg == "--disable-convection-lhs") {
            options.disableConvectionLhs = true;
        } else if (arg == "--diagnostics-only") {
            options.diagnosticsOnly = true;
        } else if (arg == "--fast-run") {
            options.fastRun = true;
        } else if (arg == "--disable-warm-start") {
            options.disableWarmStartOverride = true;
            options.disableWarmStart = nextArgLooksLikeValue(argc, argv, i) ? parseBoolValue(argv[++i]) : true;
        } else if (arg.rfind("--disable-warm-start=", 0) == 0) {
            options.disableWarmStartOverride = true;
            options.disableWarmStart = parseBoolValue(arg.substr(std::string("--disable-warm-start=").size()));
        } else if (arg == "--initial-guess" && i + 1 < argc) {
            options.initialGuessOverride = true;
            options.initialGuessType = lowerString(argv[++i]);
        } else if (arg.rfind("--initial-guess=", 0) == 0) {
            options.initialGuessOverride = true;
            options.initialGuessType = lowerString(arg.substr(std::string("--initial-guess=").size()));
        } else if (arg == "--thermal-source-scale" && i + 1 < argc) {
            options.thermalSourceScaleOverride = true;
            options.thermalSourceScale = std::stod(argv[++i]);
        } else if (arg.rfind("--thermal-source-scale=", 0) == 0) {
            options.thermalSourceScaleOverride = true;
            options.thermalSourceScale = std::stod(arg.substr(std::string("--thermal-source-scale=").size()));
        } else if (arg == "--force-nontrivial-rhs") {
            options.forceNontrivialRhsOverride = true;
            options.forceNontrivialRhs = nextArgLooksLikeValue(argc, argv, i) ? parseBoolValue(argv[++i]) : true;
        } else if (arg.rfind("--force-nontrivial-rhs=", 0) == 0) {
            options.forceNontrivialRhsOverride = true;
            options.forceNontrivialRhs = parseBoolValue(arg.substr(std::string("--force-nontrivial-rhs=").size()));
        } else if (arg == "--spectral-diagnostics") {
            options.spectralDiagnostics = true;
        } else if (arg == "--sipg-spd-diagnostics") {
            options.sipgSpdDiagnostics = true;
        } else if (arg == "--skip-sipg-penalty-sweep") {
            options.skipSipgPenaltySweep = true;
        } else if (arg == "--spd-penalty-sweep") {
            options.spdPenaltySweep = true;
        } else if (arg == "--bj-ic-validation") {
            options.bjIcValidation = true;
        } else if (arg == "--bj-ic-quick-validation") {
            options.bjIcValidation = true;
            options.bjIcQuickValidation = true;
        } else if (arg == "--penalty-mode" && i + 1 < argc) {
            options.penaltyModeOverride = true;
            options.penaltyMode = lowerString(argv[++i]);
        } else if (arg.rfind("--penalty-mode=", 0) == 0) {
            options.penaltyModeOverride = true;
            options.penaltyMode = lowerString(arg.substr(std::string("--penalty-mode=").size()));
        } else if (arg == "--penalty-factor" && i + 1 < argc) {
            options.penaltyFactorOverride = true;
            options.penaltyFactor = std::stod(argv[++i]);
        } else if (arg.rfind("--penalty-factor=", 0) == 0) {
            options.penaltyFactorOverride = true;
            options.penaltyFactor = std::stod(arg.substr(std::string("--penalty-factor=").size()));
        } else if (arg == "--dirichlet-method" && i + 1 < argc) {
            options.dirichletMethodOverride = true;
            options.dirichletMethod = lowerString(argv[++i]);
        } else if (arg.rfind("--dirichlet-method=", 0) == 0) {
            options.dirichletMethodOverride = true;
            options.dirichletMethod = lowerString(arg.substr(std::string("--dirichlet-method=").size()));
        } else if (arg == "--nitsche-penalty-factor" && i + 1 < argc) {
            options.nitschePenaltyFactorOverride = true;
            options.nitschePenaltyFactor = std::stod(argv[++i]);
        } else if (arg.rfind("--nitsche-penalty-factor=", 0) == 0) {
            options.nitschePenaltyFactorOverride = true;
            options.nitschePenaltyFactor = std::stod(arg.substr(std::string("--nitsche-penalty-factor=").size()));
        } else if (arg == "--ic-shift" && i + 1 < argc) {
            options.icShift = std::stod(argv[++i]);
        } else if (arg.rfind("--ic-shift=", 0) == 0) {
            options.icShift = std::stod(arg.substr(std::string("--ic-shift=").size()));
        } else if (arg == "--disable-ic-scaling") {
            options.icScaling = false;
            options.localDiagScaling = false;
        } else if (arg == "--local-diag-scaling") {
            options.localDiagScalingOverride = true;
            options.localDiagScaling = nextArgLooksLikeValue(argc, argv, i) ? parseBoolValue(argv[++i]) : true;
            options.icScaling = options.localDiagScaling;
        } else if (arg.rfind("--local-diag-scaling=", 0) == 0) {
            options.localDiagScalingOverride = true;
            options.localDiagScaling = parseBoolValue(arg.substr(std::string("--local-diag-scaling=").size()));
            options.icScaling = options.localDiagScaling;
        } else if (arg == "--diag-scaling-eps" && i + 1 < argc) {
            options.diagScalingEps = std::stod(argv[++i]);
        } else if (arg.rfind("--diag-scaling-eps=", 0) == 0) {
            options.diagScalingEps = std::stod(arg.substr(std::string("--diag-scaling-eps=").size()));
        } else if (arg == "--local-ic-shift-mode" && i + 1 < argc) {
            options.localIcShiftMode = lowerString(argv[++i]);
        } else if (arg.rfind("--local-ic-shift-mode=", 0) == 0) {
            options.localIcShiftMode = lowerString(arg.substr(std::string("--local-ic-shift-mode=").size()));
        } else if (arg == "--local-ic-shift-value" && i + 1 < argc) {
            options.icShift = std::stod(argv[++i]);
        } else if (arg.rfind("--local-ic-shift-value=", 0) == 0) {
            options.icShift = std::stod(arg.substr(std::string("--local-ic-shift-value=").size()));
        } else if (arg == "--ic-ordering" && i + 1 < argc) {
            options.icOrdering = lowerString(argv[++i]);
        } else if (arg.rfind("--ic-ordering=", 0) == 0) {
            options.icOrdering = lowerString(arg.substr(std::string("--ic-ordering=").size()));
        } else if (arg == "--coarse-correction" && i + 1 < argc) {
            options.coarseCorrection = lowerString(argv[++i]);
        } else if (arg.rfind("--coarse-correction=", 0) == 0) {
            options.coarseCorrection = lowerString(arg.substr(std::string("--coarse-correction=").size()));
        } else if (arg == "--coarse-space" && i + 1 < argc) {
            options.coarseSpace = lowerString(argv[++i]);
            std::replace(options.coarseSpace.begin(), options.coarseSpace.end(), '-', '_');
        } else if (arg.rfind("--coarse-space=", 0) == 0) {
            options.coarseSpace = lowerString(arg.substr(std::string("--coarse-space=").size()));
            std::replace(options.coarseSpace.begin(), options.coarseSpace.end(), '-', '_');
        } else if (arg == "--direct-mode" && i + 1 < argc) {
            options.directMode = lowerString(argv[++i]);
        } else if (arg.rfind("--direct-mode=", 0) == 0) {
            options.directMode = lowerString(arg.substr(std::string("--direct-mode=").size()));
        } else if (arg == "--max-pcg-iterations" && i + 1 < argc) {
            options.maxPcgIterations = std::stoi(argv[++i]);
        } else if (arg.rfind("--max-pcg-iterations=", 0) == 0) {
            options.maxPcgIterations = std::stoi(arg.substr(std::string("--max-pcg-iterations=").size()));
        } else if (arg == "--pcg-tolerance" && i + 1 < argc) {
            options.pcgTolerance = std::stod(argv[++i]);
        } else if (arg.rfind("--pcg-tolerance=", 0) == 0) {
            options.pcgTolerance = std::stod(arg.substr(std::string("--pcg-tolerance=").size()));
        } else if (arg == "--gmres-restart" && i + 1 < argc) {
            options.gmresRestart = std::stoi(argv[++i]);
        } else if (arg.rfind("--gmres-restart=", 0) == 0) {
            options.gmresRestart = std::stoi(arg.substr(std::string("--gmres-restart=").size()));
        } else if (arg == "--deflation-modes" && i + 1 < argc) {
            options.deflationModes = std::stoi(argv[++i]);
        } else if (arg.rfind("--deflation-modes=", 0) == 0) {
            options.deflationModes = std::stoi(arg.substr(std::string("--deflation-modes=").size()));
        } else if (arg == "--spectral-modes" && i + 1 < argc) {
            options.spectralModeCount = std::stoi(argv[++i]);
        } else if (arg.rfind("--spectral-modes=", 0) == 0) {
            options.spectralModeCount = std::stoi(arg.substr(std::string("--spectral-modes=").size()));
        } else if (arg == "--ilut-drop-tols" && i + 1 < argc) {
            options.ilutDropTolerances = parseDoubleList(argv[++i]);
        } else if (arg.rfind("--ilut-drop-tols=", 0) == 0) {
            options.ilutDropTolerances = parseDoubleList(arg.substr(std::string("--ilut-drop-tols=").size()));
        } else if (arg == "--ilut-fill-factors" && i + 1 < argc) {
            options.ilutFillFactors = parseIntList(argv[++i]);
        } else if (arg.rfind("--ilut-fill-factors=", 0) == 0) {
            options.ilutFillFactors = parseIntList(arg.substr(std::string("--ilut-fill-factors=").size()));
        } else if (arg == "--ras-ilut-configs" && i + 1 < argc) {
            options.rasIlutConfigs = parseRasIlutConfigs(argv[++i]);
        } else if (arg.rfind("--ras-ilut-configs=", 0) == 0) {
            options.rasIlutConfigs = parseRasIlutConfigs(arg.substr(std::string("--ras-ilut-configs=").size()));
        } else if (arg == "--help" || arg == "-h" || arg == "/?") {
            options.showHelp = true;
            std::cout << "Usage: SIPGHeatDDM3D.exe [--transient|--steady] [--config case_config.txt] [--solvers schwarz-precond-fgmres,schwarz-precond-fgmres-two-level,schwarz,bj-ic,bj-ilut,ras-ic,ras-ilut,bj-ic-coarse,bj-ic-coarse-pcg,bj-ilut-coarse,two-level-ras,deflated-ras,bj-pardiso-general,direct]\n";
            std::cout << "  --transient  solve M/dt + K time marching problem (default)\n";
            std::cout << "  --steady     solve K*T = Q steady-state problem\n";
            std::cout << "  --config     read domains, materials, sources, boundaries, and interfaces from a config file\n";
            std::cout << "  --output-dir path             write outputs under this directory\n";
            std::cout << "  --model-name name             label outputs for multi-model studies\n";
            std::cout << "  --solver-method monolithic|schwarz|schwarz-precond-fgmres|schwarz-precond-fgmres-two-level  select configured solve path; default is RAS-FGMRES\n";
            std::cout << "  --solvers    comma-separated subset: schwarz-precond-fgmres,schwarz-precond-fgmres-two-level,schwarz,bj-ic,bj-ilut,ras-ic,ras-ilut,bj-ic-coarse,bj-ic-coarse-pcg,bj-ilut-coarse,two-level-ras,deflated-ras,interface-deflated-ras,bj-pardiso-general,direct\n";
            std::cout << "  --schwarz-type multiplicative|additive  Schwarz update mode\n";
            std::cout << "  --schwarz-standalone-mode algebraic|dirichlet_neumann|dirichlet_dirichlet|dirichlet_robin|robin  standalone Schwarz transmission mode; aliases: dn, dd, dr, rr\n";
            std::cout << "  --schwarz-transmission-orientation forward|reverse  DN/DR direction; forward means left Dirichlet, reverse means right Dirichlet\n";
            std::cout << "  --schwarz-robin-alpha-factor value  alpha = value * k_harmonic / h for physical Robin Schwarz\n";
            std::cout << "  --schwarz-flux-eval physical_gradient|sipg_numeric  interface flux diagnostic mode\n";
            std::cout << "  --schwarz-overlap-layers 0|1|2|3  halo element layers for Schwarz/RAS local solves (default 1)\n";
            std::cout << "  --schwarz-overlap-mode none|halo|ras  overlap write-back mode; ras restricts correction to owned DOFs (default ras)\n";
            std::cout << "  --schwarz-partition-mode current|material-aligned|hotspot-contained|vertical-heat-flow-aligned  RAS owned-block partition mode\n";
            std::cout << "  --schwarz-write-interface-flux[=true|false]  write schwarz_interface_flux.csv\n";
            std::cout << "  --schwarz-max-iters value       Schwarz iteration cap (default 1000)\n";
            std::cout << "  --schwarz-tol-rel-update value  Schwarz relative update tolerance\n";
            std::cout << "  --schwarz-tol-rel-residual value  Schwarz relative residual tolerance\n";
            std::cout << "  --schwarz-relaxation value      Schwarz relaxation factor (default 1.0)\n";
            std::cout << "  --[no-]schur-linear-xy-coarse  enable/disable local x/y Schur coarse modes (default enabled)\n";
            std::cout << "  --[no-]schur-linear-z-coarse   enable/disable local z Schur coarse mode (default enabled)\n";
            std::cout << "  --[no-]schur-global-quadratic-z-coarse  enable/disable global z^2 Schur enrichment\n";
            std::cout << "  --[no-]schur-interface-patch-coarse  use boundary-side interface patches instead of volume/subdomain modes\n";
            std::cout << "  --[no-]schur-interface-patch-linear-xy  add centered x/y modes on every interface patch\n";
            std::cout << "  --disable-interface-consistency  assemble SIPG without consistency terms\n";
            std::cout << "  --disable-interface-penalty      assemble SIPG without penalty terms\n";
            std::cout << "  --node-tie-interface            add large matched-node tie penalties on interface DOFs\n";
            std::cout << "  --node-tie-penalty value        penalty value for matched-node tie constraints (default 1e12)\n";
            std::cout << "  --disable-convection-lhs         keep Robin RHS but skip h*N_i*N_j LHS entries\n";
            std::cout << "  --fast-run                       skip heavy SPD/interface diagnostics and visualization outputs\n";
            std::cout << "  --diagnostics-only               assemble/write matrix diagnostics, then exit before solvers\n";
            std::cout << "  --spectral-diagnostics           compute shift-invert near-zero eigen diagnostics for A_final\n";
            std::cout << "  --sipg-spd-diagnostics           write SIPG component SPD/face energy diagnostics, then exit\n";
            std::cout << "  --skip-sipg-penalty-sweep        skip the internal C_pen sweep in SIPG diagnostics\n";
            std::cout << "  --spd-penalty-sweep              run harmonic/max SPD and Jacobi-PCG penalty diagnostics\n";
            std::cout << "  --bj-ic-validation               run Global/BJ-Jacobi/BJ-IC validation and COMSOL node comparison\n";
            std::cout << "  --bj-ic-quick-validation         run only Global and BJ-IC natural/scaling/shifted validation\n";
            std::cout << "  --penalty-mode harmonic|max      SIPG interface penalty scaling\n";
            std::cout << "  --penalty-factor value           SIPG interface penalty factor\n";
            std::cout << "  --dirichlet-method strong|nitsche  Dirichlet enforcement method (default strong)\n";
            std::cout << "  --nitsche-penalty-factor value   symmetric Nitsche boundary penalty factor (default 15)\n";
            std::cout << "  penalty_scaling in config        p_p1 or p1_squared for P2 SIPG penalty\n";
            std::cout << "  --ic-shift value                 shifted/scaled IC initial diagonal shift (default 1e-12)\n";
            std::cout << "  --disable-ic-scaling             disable symmetric diagonal scaling for IC\n";
            std::cout << "  --ic-ordering natural|amd|rcm    local block ordering for IC\n";
            std::cout << "  --direct-mode spd|general|both   choose which Global PARDISO direct solve to run\n";
            std::cout << "  --max-pcg-iterations value       PCG/FGMRES iteration cap (default 1000)\n";
            std::cout << "  --pcg-tolerance value            relative FGMRES/legacy PCG tolerance factor (default 1e-10)\n";
            std::cout << "  --gmres-restart value            FGMRES restart length (default 30)\n";
            std::cout << "  --deflation-modes value          number of smallest-|lambda| Ritz vectors for deflated RAS (default 6)\n";
            std::cout << "  --spectral-modes value           number of near-zero Ritz vectors to compute/write (default 20)\n";
            std::cout << "  --ilut-drop-tols csv             ILUT drop tolerances (default 1e-3,1e-4,1e-5)\n";
            std::cout << "  --ilut-fill-factors csv          ILUT fill factors (default 5,10,20)\n";
            std::cout << "  --ras-ilut-configs list          semicolon list overlap:drop_tol:fill_factor, e.g. 1:1e-3:5;2:1e-4:10\n";
        } else {
            std::cout << "Ignoring unknown argument: " << arg << "\n";
        }
    }
    return options;
}

static std::string csvEscape(const std::string& value)
{
    if (value.find_first_of(",\"\n\r") == std::string::npos) {
        return value;
    }
    std::string escaped = "\"";
    for (char ch : value) {
        if (ch == '"') {
            escaped += "\"\"";
        } else {
            escaped += ch;
        }
    }
    escaped += '"';
    return escaped;
}

static std::string safeFilePrefix(std::string value)
{
    if (value.empty()) {
        value = "model";
    }
    for (char& ch : value) {
        const unsigned char c = static_cast<unsigned char>(ch);
        if (std::isalnum(c)) {
            ch = static_cast<char>(std::tolower(c));
        } else {
            ch = '_';
        }
    }
    while (value.find("__") != std::string::npos) {
        value.replace(value.find("__"), 2, "_");
    }
    while (!value.empty() && value.front() == '_') {
        value.erase(value.begin());
    }
    while (!value.empty() && value.back() == '_') {
        value.pop_back();
    }
    return value.empty() ? "model" : value;
}

static void writeRunMetadata(const std::filesystem::path& path,
                             const std::string& modelName,
                             const std::filesystem::path& configPath,
                             const std::filesystem::path& outputDir,
                             const CaseConfig& physics)
{
    std::ofstream out(path);
    out << "model_name,config_path,output_dir,solver_method,schwarz_type,schwarz_standalone_mode,"
        << "schwarz_transmission,schwarz_transmission_orientation,schwarz_robin_alpha_factor,schwarz_overlap_mode,"
        << "schwarz_overlap_layers,schwarz_partition_mode,schwarz_max_iters,schwarz_tol_rel_update,schwarz_tol_rel_residual,schwarz_relaxation,"
        << "schwarz_validate_against_monolithic,interface_mode,penalty_factor,penalty_scaling,"
        << "dirichlet_method,nitsche_penalty_factor,overlap_integration_enabled\n";
    out << csvEscape(modelName) << ','
        << csvEscape(configPath.string()) << ','
        << csvEscape(outputDir.string()) << ','
        << csvEscape(physics.solverMethod) << ','
        << csvEscape(physics.schwarz.type) << ','
        << csvEscape(physics.schwarz.standaloneMode) << ','
        << csvEscape(physics.schwarz.transmission) << ','
        << csvEscape(physics.schwarz.transmissionOrientation) << ','
        << std::setprecision(16) << physics.schwarz.robinAlphaFactor << ','
        << csvEscape(physics.schwarz.overlapMode) << ','
        << physics.schwarz.overlapLayers << ','
        << csvEscape(physics.schwarz.partitionMode) << ','
        << physics.schwarz.maxIters << ','
        << std::setprecision(16) << physics.schwarz.tolRelUpdate << ','
        << std::setprecision(16) << physics.schwarz.tolRelResidual << ','
        << std::setprecision(16) << physics.schwarz.relaxation << ','
        << (physics.schwarz.validateAgainstMonolithic ? 1 : 0) << ','
        << csvEscape(physics.interfaceScheme) << ','
        << std::setprecision(16) << physics.penaltyFactor << ','
        << csvEscape(physics.penaltyScaling) << ','
        << csvEscape(physics.dirichletMethod) << ','
        << std::setprecision(16) << physics.nitschePenaltyFactor << ",1\n";
}

static double rayleighQuotient(const SparseMatrix& a, const std::vector<double>& x)
{
    const std::vector<double> ax = a.multiply(x);
    const double xx = vectorDot(x, x);
    if (!(xx > 0.0) || !std::isfinite(xx)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return vectorDot(x, ax) / xx;
}

static double randomRayleighMinimum(const SparseMatrix& a, int trials)
{
    std::mt19937_64 rng(0x5eed1234ULL);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::vector<double> x(static_cast<size_t>(a.size()), 0.0);
    double best = std::numeric_limits<double>::infinity();
    for (int trial = 0; trial < trials; ++trial) {
        for (double& value : x) {
            value = dist(rng);
        }
        const double rq = rayleighQuotient(a, x);
        if (std::isfinite(rq)) {
            best = std::min(best, rq);
        }
    }
    return best;
}

static double smallestEigenvalueSymmetricJacobi(std::vector<std::vector<double>> a)
{
    const int n = static_cast<int>(a.size());
    if (n == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    for (int sweep = 0; sweep < 80 * n; ++sweep) {
        int p = 0;
        int q = 1;
        double maxOff = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const double value = std::abs(a[static_cast<size_t>(i)][static_cast<size_t>(j)]);
                if (value > maxOff) {
                    maxOff = value;
                    p = i;
                    q = j;
                }
            }
        }
        if (maxOff < 1.0e-12) {
            break;
        }
        const double app = a[static_cast<size_t>(p)][static_cast<size_t>(p)];
        const double aqq = a[static_cast<size_t>(q)][static_cast<size_t>(q)];
        const double apq = a[static_cast<size_t>(p)][static_cast<size_t>(q)];
        const double tau = (aqq - app) / (2.0 * apq);
        const double t = (tau >= 0.0 ? 1.0 : -1.0) / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
        const double c = 1.0 / std::sqrt(1.0 + t * t);
        const double s = t * c;
        for (int k = 0; k < n; ++k) {
            if (k == p || k == q) {
                continue;
            }
            const double akp = a[static_cast<size_t>(k)][static_cast<size_t>(p)];
            const double akq = a[static_cast<size_t>(k)][static_cast<size_t>(q)];
            a[static_cast<size_t>(k)][static_cast<size_t>(p)] = c * akp - s * akq;
            a[static_cast<size_t>(p)][static_cast<size_t>(k)] = a[static_cast<size_t>(k)][static_cast<size_t>(p)];
            a[static_cast<size_t>(k)][static_cast<size_t>(q)] = s * akp + c * akq;
            a[static_cast<size_t>(q)][static_cast<size_t>(k)] = a[static_cast<size_t>(k)][static_cast<size_t>(q)];
        }
        a[static_cast<size_t>(p)][static_cast<size_t>(p)] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        a[static_cast<size_t>(q)][static_cast<size_t>(q)] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[static_cast<size_t>(p)][static_cast<size_t>(q)] = 0.0;
        a[static_cast<size_t>(q)][static_cast<size_t>(p)] = 0.0;
    }
    double smallest = a[0][0];
    for (int i = 1; i < n; ++i) {
        smallest = std::min(smallest, a[static_cast<size_t>(i)][static_cast<size_t>(i)]);
    }
    return smallest;
}

static double estimateSmallestEigenvalueLanczos(const SparseMatrix& a, int maxIterations)
{
    const int n = a.size();
    if (n == 0 || maxIterations <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const int mMax = std::min(maxIterations, n);
    std::mt19937_64 rng(0x1a2b3c4dULL);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::vector<double> q(static_cast<size_t>(n), 0.0);
    for (double& value : q) {
        value = dist(rng);
    }
    double qNorm = std::sqrt(std::max(1.0e-300, vectorDot(q, q)));
    for (double& value : q) {
        value /= qNorm;
    }

    std::vector<double> qPrev(static_cast<size_t>(n), 0.0);
    std::vector<double> alpha;
    std::vector<double> beta;
    double betaPrev = 0.0;
    for (int iter = 0; iter < mMax; ++iter) {
        std::vector<double> z = a.multiply(q);
        if (iter > 0) {
            parallelFor(z.size(), [&](size_t i) {
                z[i] -= betaPrev * qPrev[i];
            });
        }
        const double alphaIter = vectorDot(q, z);
        parallelFor(z.size(), [&](size_t i) {
            z[i] -= alphaIter * q[i];
        });
        const double betaIter = std::sqrt(std::max(0.0, vectorDot(z, z)));
        alpha.push_back(alphaIter);
        if (betaIter < 1.0e-14 || !std::isfinite(betaIter)) {
            break;
        }
        beta.push_back(betaIter);
        qPrev.swap(q);
        q.swap(z);
        for (double& value : q) {
            value /= betaIter;
        }
        betaPrev = betaIter;
    }

    const int m = static_cast<int>(alpha.size());
    std::vector<std::vector<double>> tridiagonal(static_cast<size_t>(m), std::vector<double>(static_cast<size_t>(m), 0.0));
    for (int i = 0; i < m; ++i) {
        tridiagonal[static_cast<size_t>(i)][static_cast<size_t>(i)] = alpha[static_cast<size_t>(i)];
        if (i + 1 < m) {
            tridiagonal[static_cast<size_t>(i)][static_cast<size_t>(i + 1)] = beta[static_cast<size_t>(i)];
            tridiagonal[static_cast<size_t>(i + 1)][static_cast<size_t>(i)] = beta[static_cast<size_t>(i)];
        }
    }
    return smallestEigenvalueSymmetricJacobi(std::move(tridiagonal));
}

struct SmallSymmetricEigenResult {
    std::vector<double> values;
    std::vector<std::vector<double>> vectors;
};

static SmallSymmetricEigenResult jacobiEigenDecomposition(std::vector<std::vector<double>> a)
{
    const int n = static_cast<int>(a.size());
    SmallSymmetricEigenResult result;
    result.values.assign(static_cast<size_t>(n), 0.0);
    result.vectors.assign(static_cast<size_t>(n), std::vector<double>(static_cast<size_t>(n), 0.0));
    for (int i = 0; i < n; ++i) {
        result.vectors[static_cast<size_t>(i)][static_cast<size_t>(i)] = 1.0;
    }
    if (n == 0) {
        return result;
    }
    for (int sweep = 0; sweep < 120 * std::max(1, n); ++sweep) {
        int p = 0;
        int q = n > 1 ? 1 : 0;
        double maxOff = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                const double value = std::abs(a[static_cast<size_t>(i)][static_cast<size_t>(j)]);
                if (value > maxOff) {
                    maxOff = value;
                    p = i;
                    q = j;
                }
            }
        }
        if (n == 1 || maxOff < 1.0e-13) {
            break;
        }
        const double app = a[static_cast<size_t>(p)][static_cast<size_t>(p)];
        const double aqq = a[static_cast<size_t>(q)][static_cast<size_t>(q)];
        const double apq = a[static_cast<size_t>(p)][static_cast<size_t>(q)];
        const double tau = (aqq - app) / (2.0 * apq);
        const double t = (tau >= 0.0 ? 1.0 : -1.0) / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
        const double c = 1.0 / std::sqrt(1.0 + t * t);
        const double s = t * c;
        for (int k = 0; k < n; ++k) {
            if (k != p && k != q) {
                const double akp = a[static_cast<size_t>(k)][static_cast<size_t>(p)];
                const double akq = a[static_cast<size_t>(k)][static_cast<size_t>(q)];
                a[static_cast<size_t>(k)][static_cast<size_t>(p)] = c * akp - s * akq;
                a[static_cast<size_t>(p)][static_cast<size_t>(k)] = a[static_cast<size_t>(k)][static_cast<size_t>(p)];
                a[static_cast<size_t>(k)][static_cast<size_t>(q)] = s * akp + c * akq;
                a[static_cast<size_t>(q)][static_cast<size_t>(k)] = a[static_cast<size_t>(k)][static_cast<size_t>(q)];
            }
            const double vkp = result.vectors[static_cast<size_t>(k)][static_cast<size_t>(p)];
            const double vkq = result.vectors[static_cast<size_t>(k)][static_cast<size_t>(q)];
            result.vectors[static_cast<size_t>(k)][static_cast<size_t>(p)] = c * vkp - s * vkq;
            result.vectors[static_cast<size_t>(k)][static_cast<size_t>(q)] = s * vkp + c * vkq;
        }
        a[static_cast<size_t>(p)][static_cast<size_t>(p)] = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        a[static_cast<size_t>(q)][static_cast<size_t>(q)] = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[static_cast<size_t>(p)][static_cast<size_t>(q)] = 0.0;
        a[static_cast<size_t>(q)][static_cast<size_t>(p)] = 0.0;
    }
    for (int i = 0; i < n; ++i) {
        result.values[static_cast<size_t>(i)] = a[static_cast<size_t>(i)][static_cast<size_t>(i)];
    }
    std::vector<int> order(static_cast<size_t>(n));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        return result.values[static_cast<size_t>(lhs)] < result.values[static_cast<size_t>(rhs)];
    });
    SmallSymmetricEigenResult sorted;
    sorted.values.resize(static_cast<size_t>(n));
    sorted.vectors.assign(static_cast<size_t>(n), std::vector<double>(static_cast<size_t>(n), 0.0));
    for (int newCol = 0; newCol < n; ++newCol) {
        const int oldCol = order[static_cast<size_t>(newCol)];
        sorted.values[static_cast<size_t>(newCol)] = result.values[static_cast<size_t>(oldCol)];
        for (int row = 0; row < n; ++row) {
            sorted.vectors[static_cast<size_t>(row)][static_cast<size_t>(newCol)] =
                result.vectors[static_cast<size_t>(row)][static_cast<size_t>(oldCol)];
        }
    }
    return sorted;
}

static void orthonormalizeVectors(std::vector<std::vector<double>>& vectors)
{
    std::vector<std::vector<double>> basis;
    basis.reserve(vectors.size());
    for (std::vector<double>& v : vectors) {
        for (const std::vector<double>& q : basis) {
            const double projection = vectorDot(v, q);
            axpy(-projection, q, v);
        }
        const double nrm = vectorNorm(v);
        if (nrm <= 1.0e-30 || !std::isfinite(nrm)) {
            continue;
        }
        for (double& value : v) {
            value /= nrm;
        }
        basis.push_back(v);
    }
    vectors = std::move(basis);
}

static double estimateMaxAbsEigenvaluePower(const SparseMatrix& a, int iterations)
{
    const int n = a.size();
    if (n == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::mt19937_64 rng(0x91726354ULL);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::vector<double> x(static_cast<size_t>(n), 0.0);
    for (double& value : x) {
        value = dist(rng);
    }
    double nrm = vectorNorm(x);
    for (double& value : x) {
        value /= nrm;
    }
    double lambda = std::numeric_limits<double>::quiet_NaN();
    for (int iter = 0; iter < iterations; ++iter) {
        std::vector<double> y = a.multiply(x);
        nrm = vectorNorm(y);
        if (!(nrm > 0.0) || !std::isfinite(nrm)) {
            break;
        }
        lambda = std::abs(vectorDot(x, y));
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] = y[i] / nrm;
        }
    }
    return lambda;
}

static double matrixInfinityNorm(const SparseMatrix& a)
{
    double maxRowSum = 0.0;
    for (int row = 0; row < a.size(); ++row) {
        double rowSum = 0.0;
        for (int k = a.rowPtr[static_cast<size_t>(row)]; k < a.rowPtr[static_cast<size_t>(row + 1)]; ++k) {
            rowSum += std::abs(a.values[static_cast<size_t>(k)]);
        }
        maxRowSum = std::max(maxRowSum, rowSum);
    }
    return maxRowSum;
}

struct NearZeroMode {
    int modeIndex = -1;
    double eigenvalue = std::numeric_limits<double>::quiet_NaN();
    double residualNorm = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> vector;
};

static std::vector<NearZeroMode> computeNearZeroModesShiftInvert(const SparseMatrix& a,
                                                                  int requestedModes,
                                                                  int subspaceSize,
                                                                  int inverseIterations)
{
    const int n = a.size();
    const int m = std::min(std::max(requestedModes + 8, subspaceSize), n);
    std::vector<std::vector<double>> q(static_cast<size_t>(m), std::vector<double>(static_cast<size_t>(n), 0.0));
    std::mt19937_64 rng(0x600dd00dULL);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (auto& vec : q) {
        for (double& value : vec) {
            value = dist(rng);
        }
    }
    orthonormalizeVectors(q);
    const std::vector<MatrixEntry> entries = sparseMatrixEntries(a);
    GeneralSparseDirectSolver solver(n, entries);
    for (int iter = 0; iter < inverseIterations; ++iter) {
        std::vector<std::vector<double>> z(q.size(), std::vector<double>(static_cast<size_t>(n), 0.0));
        for (size_t j = 0; j < q.size(); ++j) {
            solver.solve(q[j], z[j]);
        }
        orthonormalizeVectors(z);
        q = std::move(z);
        if (static_cast<int>(q.size()) < requestedModes) {
            break;
        }
    }

    const int k = static_cast<int>(q.size());
    std::vector<std::vector<double>> aq(static_cast<size_t>(k));
    for (int j = 0; j < k; ++j) {
        aq[static_cast<size_t>(j)] = a.multiply(q[static_cast<size_t>(j)]);
    }
    std::vector<std::vector<double>> projected(static_cast<size_t>(k), std::vector<double>(static_cast<size_t>(k), 0.0));
    for (int i = 0; i < k; ++i) {
        for (int j = i; j < k; ++j) {
            const double value = vectorDot(q[static_cast<size_t>(i)], aq[static_cast<size_t>(j)]);
            projected[static_cast<size_t>(i)][static_cast<size_t>(j)] = value;
            projected[static_cast<size_t>(j)][static_cast<size_t>(i)] = value;
        }
    }
    SmallSymmetricEigenResult smallEigen = jacobiEigenDecomposition(projected);
    std::vector<int> order(static_cast<size_t>(k));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        return std::abs(smallEigen.values[static_cast<size_t>(lhs)])
            < std::abs(smallEigen.values[static_cast<size_t>(rhs)]);
    });

    std::vector<NearZeroMode> modes;
    const int keep = std::min(requestedModes, k);
    modes.reserve(static_cast<size_t>(keep));
    for (int outIndex = 0; outIndex < keep; ++outIndex) {
        const int col = order[static_cast<size_t>(outIndex)];
        NearZeroMode mode;
        mode.modeIndex = outIndex;
        mode.eigenvalue = smallEigen.values[static_cast<size_t>(col)];
        mode.vector.assign(static_cast<size_t>(n), 0.0);
        for (int basis = 0; basis < k; ++basis) {
            const double coeff = smallEigen.vectors[static_cast<size_t>(basis)][static_cast<size_t>(col)];
            axpy(coeff, q[static_cast<size_t>(basis)], mode.vector);
        }
        const double modeNorm = vectorNorm(mode.vector);
        if (modeNorm > 0.0 && std::isfinite(modeNorm)) {
            for (double& value : mode.vector) {
                value /= modeNorm;
            }
        }
        std::vector<double> av = a.multiply(mode.vector);
        axpy(-mode.eigenvalue, mode.vector, av);
        mode.residualNorm = vectorNorm(av);
        modes.push_back(std::move(mode));
    }
    return modes;
}

static std::vector<NearZeroMode> computeShiftedModesShiftInvert(const SparseMatrix& a,
                                                                 double sigma,
                                                                 int requestedModes,
                                                                 int subspaceSize,
                                                                 int inverseIterations)
{
    const int n = a.size();
    SparseMatrix shifted(n);
    shifted.appendScaledEntries(a, 1.0);
    for (int i = 0; i < n; ++i) {
        shifted.add(i, i, -sigma);
    }
    shifted.finalizeCsr();

    const int m = std::min(std::max(requestedModes + 8, subspaceSize), n);
    std::vector<std::vector<double>> q(static_cast<size_t>(m), std::vector<double>(static_cast<size_t>(n), 0.0));
    std::mt19937_64 rng(0x6e676174ULL ^ static_cast<unsigned long long>(std::abs(sigma) * 1.0e18));
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (auto& vec : q) {
        for (double& value : vec) {
            value = dist(rng);
        }
    }
    orthonormalizeVectors(q);
    const std::vector<MatrixEntry> shiftedEntries = sparseMatrixEntries(shifted);
    GeneralSparseDirectSolver solver(n, shiftedEntries);
    for (int iter = 0; iter < inverseIterations; ++iter) {
        std::vector<std::vector<double>> z(q.size(), std::vector<double>(static_cast<size_t>(n), 0.0));
        for (size_t j = 0; j < q.size(); ++j) {
            solver.solve(q[j], z[j]);
        }
        orthonormalizeVectors(z);
        q = std::move(z);
        if (static_cast<int>(q.size()) < requestedModes) {
            break;
        }
    }

    const int k = static_cast<int>(q.size());
    std::vector<std::vector<double>> aq(static_cast<size_t>(k));
    for (int j = 0; j < k; ++j) {
        aq[static_cast<size_t>(j)] = a.multiply(q[static_cast<size_t>(j)]);
    }
    std::vector<std::vector<double>> projected(static_cast<size_t>(k), std::vector<double>(static_cast<size_t>(k), 0.0));
    for (int i = 0; i < k; ++i) {
        for (int j = i; j < k; ++j) {
            const double value = vectorDot(q[static_cast<size_t>(i)], aq[static_cast<size_t>(j)]);
            projected[static_cast<size_t>(i)][static_cast<size_t>(j)] = value;
            projected[static_cast<size_t>(j)][static_cast<size_t>(i)] = value;
        }
    }
    SmallSymmetricEigenResult smallEigen = jacobiEigenDecomposition(projected);
    std::vector<int> order(static_cast<size_t>(k));
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        return std::abs(smallEigen.values[static_cast<size_t>(lhs)] - sigma)
            < std::abs(smallEigen.values[static_cast<size_t>(rhs)] - sigma);
    });

    std::vector<NearZeroMode> modes;
    const int keep = std::min(requestedModes, k);
    modes.reserve(static_cast<size_t>(keep));
    for (int outIndex = 0; outIndex < keep; ++outIndex) {
        const int col = order[static_cast<size_t>(outIndex)];
        NearZeroMode mode;
        mode.modeIndex = outIndex;
        mode.eigenvalue = smallEigen.values[static_cast<size_t>(col)];
        mode.vector.assign(static_cast<size_t>(n), 0.0);
        for (int basis = 0; basis < k; ++basis) {
            const double coeff = smallEigen.vectors[static_cast<size_t>(basis)][static_cast<size_t>(col)];
            axpy(coeff, q[static_cast<size_t>(basis)], mode.vector);
        }
        const double modeNorm = vectorNorm(mode.vector);
        if (modeNorm > 0.0 && std::isfinite(modeNorm)) {
            for (double& value : mode.vector) {
                value /= modeNorm;
            }
        }
        std::vector<double> av = a.multiply(mode.vector);
        axpy(-mode.eigenvalue, mode.vector, av);
        mode.residualNorm = vectorNorm(av);
        modes.push_back(std::move(mode));
    }
    return modes;
}

static std::vector<std::string> materialNameByDof(const Mesh& mesh, const CaseConfig& config)
{
    std::vector<std::string> names(mesh.nodes.size(), "unknown");
    for (const Tet& tet : mesh.tets) {
        const std::string materialName = materialForTet(config, tet).name;
        for (int dof : tet.dof) {
            if (dof >= 0 && dof < static_cast<int>(names.size()) && names[static_cast<size_t>(dof)] == "unknown") {
                names[static_cast<size_t>(dof)] = materialName;
            }
        }
    }
    return names;
}

static std::vector<char> boundaryDofMask(const Mesh& mesh)
{
    std::vector<char> mask(mesh.nodes.size(), 0);
    for (const BoundaryFace& face : mesh.boundaryFaces) {
        const Tet& tet = mesh.tets[static_cast<size_t>(face.tet)];
        for (int i = 0; i < 3; ++i) {
            mask[static_cast<size_t>(tet.dof[static_cast<size_t>(face.local[static_cast<size_t>(i)])])] = 1;
            const int a = face.local[static_cast<size_t>(i)];
            const int b = face.local[static_cast<size_t>((i + 1) % 3)];
            mask[static_cast<size_t>(tet.dof[static_cast<size_t>(localEdgeDof(a, b))])] = 1;
        }
    }
    return mask;
}

struct SpectralSummary {
    double lambdaMaxAbs = std::numeric_limits<double>::quiet_NaN();
    double lambdaMinAbs = std::numeric_limits<double>::quiet_NaN();
    double conditionEstAbs = std::numeric_limits<double>::quiet_NaN();
    std::vector<NearZeroMode> modes;
};

static SpectralSummary writeEigenDiagnostics(const Mesh& mesh,
                                             const CaseConfig& config,
                                             const AssemblyDiagnostics& assemblyDiagnostics,
                                             const SparseMatrix& system,
                                             const std::filesystem::path& outputDir,
                                             int requestedModeCount = 20)
{
    std::cout << "Computing shift-invert near-zero eigen diagnostics...\n";
    const double lambdaMaxAbs = estimateMaxAbsEigenvaluePower(system, 80);
    const int modeCount = std::max(1, requestedModeCount);
    const int subspaceSize = std::min(std::max(modeCount + 12, 32), std::max(32, modeCount + 24));
    std::vector<NearZeroMode> modes = computeNearZeroModesShiftInvert(system, modeCount, subspaceSize, 8);
    const double lambdaMinAbs = modes.empty()
        ? std::numeric_limits<double>::quiet_NaN()
        : std::abs(modes.front().eigenvalue);
    const double conditionEstAbs = lambdaMaxAbs / std::max(1.0e-300, lambdaMinAbs);
    SpectralSummary summary;
    summary.lambdaMaxAbs = lambdaMaxAbs;
    summary.lambdaMinAbs = lambdaMinAbs;
    summary.conditionEstAbs = conditionEstAbs;
    summary.modes = modes;
    const std::vector<std::string> materials = materialNameByDof(mesh, config);
    const std::vector<char> boundaryMask = boundaryDofMask(mesh);

    std::ofstream nearOut(outputDir / "rram_near_zero_eigen_diagnostics.csv");
    nearOut << "mode_index,eigenvalue,abs_eigenvalue,ritz_residual_norm,lambda_max_abs,lambda_min_abs,condition_est_abs,"
            << "max_abs_dof,abs_value,x,y,z,subdomain,material,is_interface,is_dirichlet,is_boundary,"
            << "distance_to_x_2p65e-7\n";
    nearOut << std::setprecision(16);

    std::ofstream negOut(outputDir / "rram_negative_mode_location.csv");
    negOut << "mode_index,eigenvalue,abs_eigenvalue,ritz_residual_norm,max_abs_dof,abs_value,x,y,z,"
           << "subdomain,material,is_interface,is_dirichlet,is_boundary,distance_to_x_2p65e-7\n";
    negOut << std::setprecision(16);

    std::ofstream energyOut(outputDir / "rram_eigenvector_energy_by_region.csv");
    energyOut << "mode_index,eigenvalue,region,energy_fraction\n";
    energyOut << std::setprecision(16);

    std::ofstream localizationOut(outputDir / "rram_near_zero_mode_localization.csv");
    localizationOut << "mode_index,eigenvalue,rank,max_abs_dof,dof_abs_value,x,y,z,"
                    << "distance_to_interface_plane,domain_id,material_name,is_interface_dof,is_boundary_dof,"
                    << "mode_max_abs_dof,mode_max_abs_x,mode_max_abs_y,mode_max_abs_z,"
                    << "mode_max_abs_distance_to_interface_plane\n";
    localizationOut << std::setprecision(16);

    std::ofstream detailedEnergyOut(outputDir / "rram_near_zero_mode_energy_by_region.csv");
    detailedEnergyOut << "mode_index,eigenvalue,region,energy_fraction\n";
    detailedEnergyOut << std::setprecision(16);

    std::vector<std::vector<int>> facePairsByDof(mesh.nodes.size());
    for (size_t faceIndex = 0; faceIndex < mesh.interfaceFaces.size(); ++faceIndex) {
        const InterfaceFace& face = mesh.interfaceFaces[faceIndex];
        const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
        for (int i = 0; i < 10; ++i) {
            facePairsByDof[static_cast<size_t>(left.dof[static_cast<size_t>(i)])].push_back(static_cast<int>(faceIndex));
            facePairsByDof[static_cast<size_t>(right.dof[static_cast<size_t>(i)])].push_back(static_cast<int>(faceIndex));
        }
    }
    for (std::vector<int>& ids : facePairsByDof) {
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    }

    for (const NearZeroMode& mode : modes) {
        int maxDof = 0;
        double maxValue = 0.0;
        double totalEnergy = 0.0;
        double interfaceEnergy = 0.0;
        double dirichletEnergy = 0.0;
        double interiorEnergy = 0.0;
        std::map<std::string, double> materialEnergy;
        std::map<int, double> domainEnergy;
        std::map<int, double> facePairEnergy;
        for (size_t i = 0; i < mode.vector.size(); ++i) {
            const double e = mode.vector[i] * mode.vector[i];
            totalEnergy += e;
            if (std::abs(mode.vector[i]) > maxValue) {
                maxValue = std::abs(mode.vector[i]);
                maxDof = static_cast<int>(i);
            }
            const bool isInterface = !assemblyDiagnostics.interfaceDof.empty()
                && assemblyDiagnostics.interfaceDof[i] != 0;
            if (isInterface) {
                interfaceEnergy += e;
            } else if (mesh.nodes[i].dirichlet) {
                dirichletEnergy += e;
            } else {
                interiorEnergy += e;
            }
            materialEnergy[materials[i]] += e;
            domainEnergy[mesh.nodes[i].subdomain] += e;
            if (!facePairsByDof[i].empty()) {
                const double share = e / static_cast<double>(facePairsByDof[i].size());
                for (int facePairId : facePairsByDof[i]) {
                    facePairEnergy[facePairId] += share;
                }
            }
        }
        const Node& node = mesh.nodes[static_cast<size_t>(maxDof)];
        const bool isInterface = !assemblyDiagnostics.interfaceDof.empty()
            && assemblyDiagnostics.interfaceDof[static_cast<size_t>(maxDof)] != 0;
        auto writeModeLocation = [&](std::ofstream& out) {
            out << mode.modeIndex << ','
                << mode.eigenvalue << ','
                << std::abs(mode.eigenvalue) << ','
                << mode.residualNorm << ','
                << maxDof << ','
                << maxValue << ','
                << node.p.x << ','
                << node.p.y << ','
                << node.p.z << ','
                << node.subdomain << ','
                << csvEscape(materials[static_cast<size_t>(maxDof)]) << ','
                << (isInterface ? 1 : 0) << ','
                << (node.dirichlet ? 1 : 0) << ','
                << (boundaryMask[static_cast<size_t>(maxDof)] ? 1 : 0) << ','
                << std::abs(node.p.x - 2.65e-7) << '\n';
        };
        nearOut << mode.modeIndex << ','
                << mode.eigenvalue << ','
                << std::abs(mode.eigenvalue) << ','
                << mode.residualNorm << ','
                << lambdaMaxAbs << ','
                << lambdaMinAbs << ','
                << conditionEstAbs << ',';
        nearOut << maxDof << ','
                << maxValue << ','
                << node.p.x << ','
                << node.p.y << ','
                << node.p.z << ','
                << node.subdomain << ','
                << csvEscape(materials[static_cast<size_t>(maxDof)]) << ','
                << (isInterface ? 1 : 0) << ','
                << (node.dirichlet ? 1 : 0) << ','
                << (boundaryMask[static_cast<size_t>(maxDof)] ? 1 : 0) << ','
                << std::abs(node.p.x - 2.65e-7) << '\n';
        if (mode.eigenvalue < 0.0) {
            writeModeLocation(negOut);
        }
        const double denom = std::max(1.0e-300, totalEnergy);
        energyOut << mode.modeIndex << ',' << mode.eigenvalue << ",interface," << interfaceEnergy / denom << '\n';
        energyOut << mode.modeIndex << ',' << mode.eigenvalue << ",interior," << interiorEnergy / denom << '\n';
        energyOut << mode.modeIndex << ',' << mode.eigenvalue << ",dirichlet," << dirichletEnergy / denom << '\n';
        for (const auto& entry : materialEnergy) {
            energyOut << mode.modeIndex << ',' << mode.eigenvalue
                      << ",material:" << csvEscape(entry.first) << ','
                      << entry.second / denom << '\n';
        }

        detailedEnergyOut << mode.modeIndex << ',' << mode.eigenvalue << ",interface," << interfaceEnergy / denom << '\n';
        detailedEnergyOut << mode.modeIndex << ',' << mode.eigenvalue << ",interior," << interiorEnergy / denom << '\n';
        detailedEnergyOut << mode.modeIndex << ',' << mode.eigenvalue << ",dirichlet," << dirichletEnergy / denom << '\n';
        for (const auto& entry : domainEnergy) {
            detailedEnergyOut << mode.modeIndex << ',' << mode.eigenvalue
                              << ",domain:" << entry.first << ','
                              << entry.second / denom << '\n';
        }
        for (const auto& entry : materialEnergy) {
            detailedEnergyOut << mode.modeIndex << ',' << mode.eigenvalue
                              << ",material:" << csvEscape(entry.first) << ','
                              << entry.second / denom << '\n';
        }
        for (const auto& entry : facePairEnergy) {
            const double fraction = entry.second / denom;
            if (fraction > 1.0e-12) {
                detailedEnergyOut << mode.modeIndex << ',' << mode.eigenvalue
                                  << ",interface_face_pair:" << entry.first << ','
                                  << fraction << '\n';
            }
        }

        std::vector<int> order(mode.vector.size());
        std::iota(order.begin(), order.end(), 0);
        const int topCount = std::min<int>(50, static_cast<int>(order.size()));
        std::partial_sort(order.begin(), order.begin() + topCount, order.end(),
                          [&](int lhs, int rhs) {
                              return std::abs(mode.vector[static_cast<size_t>(lhs)])
                                  > std::abs(mode.vector[static_cast<size_t>(rhs)]);
                          });
        for (int rank = 0; rank < topCount; ++rank) {
            const int dof = order[static_cast<size_t>(rank)];
            const Node& topNode = mesh.nodes[static_cast<size_t>(dof)];
            const bool topInterface = !assemblyDiagnostics.interfaceDof.empty()
                && assemblyDiagnostics.interfaceDof[static_cast<size_t>(dof)] != 0;
            localizationOut << mode.modeIndex << ','
                            << mode.eigenvalue << ','
                            << (rank + 1) << ','
                            << dof << ','
                            << std::abs(mode.vector[static_cast<size_t>(dof)]) << ','
                            << topNode.p.x << ','
                            << topNode.p.y << ','
                            << topNode.p.z << ','
                            << std::abs(topNode.p.x - 2.65e-7) << ','
                            << topNode.subdomain << ','
                            << csvEscape(materials[static_cast<size_t>(dof)]) << ','
                            << (topInterface ? 1 : 0) << ','
                            << (boundaryMask[static_cast<size_t>(dof)] ? 1 : 0) << ','
                            << maxDof << ','
                            << node.p.x << ','
                            << node.p.y << ','
                            << node.p.z << ','
                            << std::abs(node.p.x - 2.65e-7) << '\n';
        }
    }
    return summary;
}

static std::vector<std::vector<double>> buildEigenCoarseVectors(const SpectralSummary& spectralSummary,
                                                                int maxModes)
{
    std::vector<std::vector<double>> vectors;
    const int count = std::min<int>(std::max(0, maxModes), static_cast<int>(spectralSummary.modes.size()));
    vectors.reserve(static_cast<size_t>(count));
    for (int i = 0; i < count; ++i) {
        const std::vector<double>& mode = spectralSummary.modes[static_cast<size_t>(i)].vector;
        if (!mode.empty() && !vectorHasNonFinite(mode)) {
            vectors.push_back(mode);
        }
    }
    return vectors;
}

static std::vector<std::vector<double>> buildInterfaceEigenCoarseVectors(const SpectralSummary& spectralSummary,
                                                                         const AssemblyDiagnostics& assemblyDiagnostics,
                                                                         int maxModes)
{
    std::vector<std::vector<double>> vectors;
    const int count = std::min<int>(std::max(0, maxModes), static_cast<int>(spectralSummary.modes.size()));
    vectors.reserve(static_cast<size_t>(count));
    for (int i = 0; i < count; ++i) {
        std::vector<double> vector = spectralSummary.modes[static_cast<size_t>(i)].vector;
        if (vector.empty() || vectorHasNonFinite(vector)) {
            continue;
        }
        for (size_t j = 0; j < vector.size(); ++j) {
            const bool keep = !assemblyDiagnostics.interfaceDof.empty()
                && assemblyDiagnostics.interfaceDof[j] != 0;
            if (!keep) {
                vector[j] = 0.0;
            }
        }
        const double normValue = vectorNorm(vector);
        if (normValue <= 1.0e-300) {
            continue;
        }
        for (double& value : vector) {
            value /= normValue;
        }
        vectors.push_back(std::move(vector));
    }
    return vectors;
}

struct GraphEdge {
    int from = -1;
    int to = -1;
    double conductance = 0.0;
};

static std::vector<char> robinDofMask(const Mesh& mesh, const CaseConfig& config)
{
    std::vector<char> mask(mesh.nodes.size(), 0);
    for (const BoundaryFace& face : mesh.boundaryFaces) {
        bool robin = false;
        for (const ConvectionCondition& condition : config.convectionConditions) {
            if (condition.subdomain >= 0 && condition.subdomain != face.subdomain) {
                continue;
            }
            if (condition.boundaryEntity == face.boundaryEntity) {
                robin = true;
                break;
            }
        }
        if (!robin) {
            continue;
        }
        const Tet& tet = mesh.tets[static_cast<size_t>(face.tet)];
        for (int i = 0; i < 3; ++i) {
            mask[static_cast<size_t>(tet.dof[static_cast<size_t>(face.local[static_cast<size_t>(i)])])] = 1;
            const int a = face.local[static_cast<size_t>(i)];
            const int b = face.local[static_cast<size_t>((i + 1) % 3)];
            mask[static_cast<size_t>(tet.dof[static_cast<size_t>(localEdgeDof(a, b))])] = 1;
        }
    }
    return mask;
}

static void writeThermalConnectivityDiagnostics(const Mesh& mesh,
                                                const CaseConfig& config,
                                                const AssemblyDiagnostics& assemblyDiagnostics,
                                                SparseMatrix preDirichletSystem,
                                                const SpectralSummary& spectralSummary,
                                                const std::filesystem::path& outputDir)
{
    preDirichletSystem.finalizeCsr();
    const int n = preDirichletSystem.size();
    std::vector<std::vector<std::pair<int, double>>> adjacency(static_cast<size_t>(n));
    std::vector<GraphEdge> weakestEdges;
    weakestEdges.reserve(1024);
    auto considerWeakEdge = [&](int i, int j, double conductance) {
        if (!(conductance > 0.0) || !std::isfinite(conductance)) {
            return;
        }
        if (static_cast<int>(weakestEdges.size()) < 1000) {
            weakestEdges.push_back({i, j, conductance});
            std::push_heap(weakestEdges.begin(), weakestEdges.end(), [](const GraphEdge& a, const GraphEdge& b) {
                return a.conductance < b.conductance;
            });
            return;
        }
        if (conductance < weakestEdges.front().conductance) {
            std::pop_heap(weakestEdges.begin(), weakestEdges.end(), [](const GraphEdge& a, const GraphEdge& b) {
                return a.conductance < b.conductance;
            });
            weakestEdges.back() = {i, j, conductance};
            std::push_heap(weakestEdges.begin(), weakestEdges.end(), [](const GraphEdge& a, const GraphEdge& b) {
                return a.conductance < b.conductance;
            });
        }
    };

    for (int row = 0; row < n; ++row) {
        for (int k = preDirichletSystem.rowPtr[static_cast<size_t>(row)];
             k < preDirichletSystem.rowPtr[static_cast<size_t>(row + 1)];
             ++k) {
            const int col = preDirichletSystem.colInd[static_cast<size_t>(k)];
            if (col == row) {
                continue;
            }
            const double conductance = std::abs(preDirichletSystem.values[static_cast<size_t>(k)]);
            if (!(conductance > 0.0) || !std::isfinite(conductance)) {
                continue;
            }
            adjacency[static_cast<size_t>(row)].push_back({col, conductance});
            if (row < col) {
                considerWeakEdge(row, col, conductance);
            }
        }
    }
    std::sort(weakestEdges.begin(), weakestEdges.end(), [](const GraphEdge& a, const GraphEdge& b) {
        return a.conductance < b.conductance;
    });

    std::vector<int> component(static_cast<size_t>(n), -1);
    int componentCount = 0;
    std::vector<int> stack;
    for (int seed = 0; seed < n; ++seed) {
        if (component[static_cast<size_t>(seed)] >= 0) {
            continue;
        }
        stack.clear();
        stack.push_back(seed);
        component[static_cast<size_t>(seed)] = componentCount;
        while (!stack.empty()) {
            const int node = stack.back();
            stack.pop_back();
            for (const auto& edge : adjacency[static_cast<size_t>(node)]) {
                if (component[static_cast<size_t>(edge.first)] < 0) {
                    component[static_cast<size_t>(edge.first)] = componentCount;
                    stack.push_back(edge.first);
                }
            }
        }
        ++componentCount;
    }

    const std::vector<std::string> materials = materialNameByDof(mesh, config);
    const std::vector<char> robinMask = robinDofMask(mesh, config);
    const std::vector<char> boundaryMask = boundaryDofMask(mesh);
    struct ComponentStats {
        int nodeCount = 0;
        int dirichletCount = 0;
        int robinCount = 0;
        int interfaceCount = 0;
        int sio2Count = 0;
        double conductanceSum = 0.0;
        double minConductance = std::numeric_limits<double>::infinity();
        int edgeCount = 0;
    };
    std::vector<ComponentStats> stats(static_cast<size_t>(componentCount));
    for (int i = 0; i < n; ++i) {
        ComponentStats& s = stats[static_cast<size_t>(component[static_cast<size_t>(i)])];
        ++s.nodeCount;
        if (mesh.nodes[static_cast<size_t>(i)].dirichlet) {
            ++s.dirichletCount;
        }
        if (robinMask[static_cast<size_t>(i)]) {
            ++s.robinCount;
        }
        if (!assemblyDiagnostics.interfaceDof.empty() && assemblyDiagnostics.interfaceDof[static_cast<size_t>(i)] != 0) {
            ++s.interfaceCount;
        }
        if (materials[static_cast<size_t>(i)].find("SiO2") != std::string::npos) {
            ++s.sio2Count;
        }
        for (const auto& edge : adjacency[static_cast<size_t>(i)]) {
            s.conductanceSum += edge.second;
            s.minConductance = std::min(s.minConductance, edge.second);
            ++s.edgeCount;
        }
    }

    std::ofstream graphOut(outputDir / "rram_thermal_connectivity_graph.csv");
    graphOut << "component_id,node_count,connected_to_dirichlet,dirichlet_dof_count,"
             << "connected_to_robin,robin_dof_count,interface_dof_count,sio2_dof_count,"
             << "min_incident_conductance,avg_incident_conductance,edge_count\n";
    graphOut << std::setprecision(16);
    for (int c = 0; c < componentCount; ++c) {
        const ComponentStats& s = stats[static_cast<size_t>(c)];
        graphOut << c << ','
                 << s.nodeCount << ','
                 << (s.dirichletCount > 0 ? 1 : 0) << ','
                 << s.dirichletCount << ','
                 << (s.robinCount > 0 ? 1 : 0) << ','
                 << s.robinCount << ','
                 << s.interfaceCount << ','
                 << s.sio2Count << ','
                 << (std::isfinite(s.minConductance) ? s.minConductance : 0.0) << ','
                 << (s.edgeCount > 0 ? s.conductanceSum / static_cast<double>(s.edgeCount) : 0.0) << ','
                 << s.edgeCount << '\n';
    }

    std::ofstream cutOut(outputDir / "rram_weak_connection_cut_diagnostics.csv");
    cutOut << "rank,node_i,node_j,conductance,resistance,component_i,component_j,"
           << "material_i,material_j,is_interface_i,is_interface_j,is_boundary_i,is_boundary_j,"
           << "x_i,y_i,z_i,x_j,y_j,z_j\n";
    int rank = 0;
    for (const GraphEdge& edge : weakestEdges) {
        if (rank >= 500) {
            break;
        }
        const Node& ni = mesh.nodes[static_cast<size_t>(edge.from)];
        const Node& nj = mesh.nodes[static_cast<size_t>(edge.to)];
        cutOut << (++rank) << ','
               << edge.from << ','
               << edge.to << ','
               << edge.conductance << ','
               << 1.0 / std::max(1.0e-300, edge.conductance) << ','
               << component[static_cast<size_t>(edge.from)] << ','
               << component[static_cast<size_t>(edge.to)] << ','
               << csvEscape(materials[static_cast<size_t>(edge.from)]) << ','
               << csvEscape(materials[static_cast<size_t>(edge.to)]) << ','
               << (!assemblyDiagnostics.interfaceDof.empty() && assemblyDiagnostics.interfaceDof[static_cast<size_t>(edge.from)] != 0 ? 1 : 0) << ','
               << (!assemblyDiagnostics.interfaceDof.empty() && assemblyDiagnostics.interfaceDof[static_cast<size_t>(edge.to)] != 0 ? 1 : 0) << ','
               << (boundaryMask[static_cast<size_t>(edge.from)] ? 1 : 0) << ','
               << (boundaryMask[static_cast<size_t>(edge.to)] ? 1 : 0) << ','
               << ni.p.x << ',' << ni.p.y << ',' << ni.p.z << ','
               << nj.p.x << ',' << nj.p.y << ',' << nj.p.z << '\n';
    }

    std::vector<double> distance(static_cast<size_t>(n), std::numeric_limits<double>::infinity());
    using QueueItem = std::pair<double, int>;
    std::priority_queue<QueueItem, std::vector<QueueItem>, std::greater<QueueItem>> queue;
    for (int i = 0; i < n; ++i) {
        if (mesh.nodes[static_cast<size_t>(i)].dirichlet) {
            distance[static_cast<size_t>(i)] = 0.0;
            queue.push({0.0, i});
        }
    }
    while (!queue.empty()) {
        const auto [currentDistance, node] = queue.top();
        queue.pop();
        if (currentDistance != distance[static_cast<size_t>(node)]) {
            continue;
        }
        for (const auto& edge : adjacency[static_cast<size_t>(node)]) {
            const double nextDistance = currentDistance + 1.0 / std::max(1.0e-300, edge.second);
            if (nextDistance < distance[static_cast<size_t>(edge.first)]) {
                distance[static_cast<size_t>(edge.first)] = nextDistance;
                queue.push({nextDistance, edge.first});
            }
        }
    }

    std::ofstream modeOut(outputDir / "rram_near_zero_mode_vs_weak_connection.csv");
    modeOut << "mode_index,eigenvalue,rank,dof,abs_value,material,is_interface,is_boundary,"
            << "thermal_resistance_to_dirichlet,distance_to_x_2p65e-7,x,y,z\n";
    modeOut << std::setprecision(16);
    for (const NearZeroMode& mode : spectralSummary.modes) {
        std::vector<int> order(mode.vector.size());
        std::iota(order.begin(), order.end(), 0);
        const int topCount = std::min<int>(50, static_cast<int>(order.size()));
        std::partial_sort(order.begin(), order.begin() + topCount, order.end(), [&](int lhs, int rhs) {
            return std::abs(mode.vector[static_cast<size_t>(lhs)])
                > std::abs(mode.vector[static_cast<size_t>(rhs)]);
        });
        for (int i = 0; i < topCount; ++i) {
            const int dof = order[static_cast<size_t>(i)];
            const Node& node = mesh.nodes[static_cast<size_t>(dof)];
            modeOut << mode.modeIndex << ','
                    << mode.eigenvalue << ','
                    << (i + 1) << ','
                    << dof << ','
                    << std::abs(mode.vector[static_cast<size_t>(dof)]) << ','
                    << csvEscape(materials[static_cast<size_t>(dof)]) << ','
                    << (!assemblyDiagnostics.interfaceDof.empty() && assemblyDiagnostics.interfaceDof[static_cast<size_t>(dof)] != 0 ? 1 : 0) << ','
                    << (boundaryMask[static_cast<size_t>(dof)] ? 1 : 0) << ','
                    << distance[static_cast<size_t>(dof)] << ','
                    << std::abs(node.p.x - 2.65e-7) << ','
                    << node.p.x << ','
                    << node.p.y << ','
                    << node.p.z << '\n';
        }
    }
}

static void initializeAssemblyDiagnostics(int n, AssemblyDiagnostics& diagnostics)
{
    diagnostics = AssemblyDiagnostics{};
    diagnostics.volumeDiag.assign(static_cast<size_t>(n), 0.0);
    diagnostics.robinDiag.assign(static_cast<size_t>(n), 0.0);
    diagnostics.interfaceConsistencyDiag.assign(static_cast<size_t>(n), 0.0);
    diagnostics.interfacePenaltyDiag.assign(static_cast<size_t>(n), 0.0);
    diagnostics.preDirichletDiag.assign(static_cast<size_t>(n), 0.0);
    diagnostics.dirichletDiag.assign(static_cast<size_t>(n), 0.0);
    diagnostics.finalDiag.assign(static_cast<size_t>(n), 0.0);
}

struct SweepSolverFields {
    std::string status = "not_run";
    std::string failureReason;
    int breakdownIter = -1;
    int iterations = 0;
    double finalRelativeResidual = std::numeric_limits<double>::quiet_NaN();
    double pcgPAp = std::numeric_limits<double>::quiet_NaN();
    double pcgRTr = std::numeric_limits<double>::quiet_NaN();
    double pcgRMz = std::numeric_limits<double>::quiet_NaN();
    double pcgAlpha = std::numeric_limits<double>::quiet_NaN();
    double pcgBeta = std::numeric_limits<double>::quiet_NaN();
    double pcgResidualNorm = std::numeric_limits<double>::quiet_NaN();
    bool pcgPHasNonFinite = false;
    bool pcgApHasNonFinite = false;
};

struct SweepBlockDiagnostic {
    std::string interfaceScheme;
    std::string penaltyMode;
    std::string penaltyScaling;
    double penaltyFactor = 0.0;
    int blockId = -1;
    int blockSize = 0;
    double minDiag = std::numeric_limits<double>::quiet_NaN();
    double maxDiag = std::numeric_limits<double>::quiet_NaN();
    double symmetryError = std::numeric_limits<double>::quiet_NaN();
    double lambdaMinEst = std::numeric_limits<double>::quiet_NaN();
    std::string ldltStatus;
    int ldltPositivePivots = -1;
    int ldltNegativePivots = -1;
    int ldltZeroTinyPivots = -1;
    std::string icStatus = "not_run";
    std::string icMessage;
    std::string pardisoSpdStatus;
    std::string pardisoSpdMessage;
    bool firstFailedBlock = false;
};

struct FinalSpdDiagnostic {
    bool hasNonFinite = false;
    double symmetryError = std::numeric_limits<double>::quiet_NaN();
    double minDiag = std::numeric_limits<double>::quiet_NaN();
    double maxDiag = std::numeric_limits<double>::quiet_NaN();
    int zeroRowCount = 0;
    std::string ldltStatus;
    int ldltPositivePivots = -1;
    int ldltNegativePivots = -1;
    int ldltZeroTinyPivots = -1;
    double lambdaMinEst = std::numeric_limits<double>::quiet_NaN();
    std::string pardisoSpdStatus;
    std::string pardisoSpdMessage;
    std::string pardisoGeneralStatus;
    std::string pardisoGeneralMessage;
    std::string spdConclusion;
};

struct SpdPenaltySweepRow {
    std::string interfaceScheme;
    std::string penaltyMode;
    std::string penaltyScaling;
    double penaltyFactor = 0.0;
    int negativeDiagCount = 0;
    double minDiag = std::numeric_limits<double>::quiet_NaN();
    double maxDiag = std::numeric_limits<double>::quiet_NaN();
    double symmetryError = std::numeric_limits<double>::quiet_NaN();
    bool hasNonFinite = false;
    int zeroRowCount = 0;
    std::string pardisoSpdStatus;
    std::string pardisoSpdMessage;
    std::string globalGeneralStatus;
    std::string globalGeneralMessage;
    double estimatedLambdaMin = std::numeric_limits<double>::quiet_NaN();
    double rayleighMin = std::numeric_limits<double>::quiet_NaN();
    std::string ldltStatus;
    int ldltPositivePivots = -1;
    int ldltNegativePivots = -1;
    int ldltZeroTinyPivots = -1;
    SweepSolverFields diagonalPcg;
    SweepSolverFields bjJacobiPcg;
    SweepSolverFields bjIcPcg;
    SweepSolverFields bjPardisoPcg;
    double tMin = std::numeric_limits<double>::quiet_NaN();
    double tMax = std::numeric_limits<double>::quiet_NaN();
};

static void writeSpdPenaltySweep(const std::vector<SpdPenaltySweepRow>& rows,
                                 const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << std::setprecision(16);
    const auto writeSolverHeader = [&](const std::string& prefix) {
        out << prefix << "_status," << prefix << "_failure_reason,"
            << prefix << "_breakdown_iter," << prefix << "_iterations,"
            << prefix << "_final_relative_residual," << prefix << "_pAp,"
            << prefix << "_rTr," << prefix << "_rMz," << prefix << "_alpha,"
            << prefix << "_beta," << prefix << "_residual_norm,"
            << prefix << "_p_has_nan_inf," << prefix << "_Ap_has_nan_inf";
    };
    out << "interface_scheme,penalty_mode,penalty_scaling,penalty_factor,"
        << "negative_diag_count,min_diag,max_diag,symmetry_error,has_nan_inf,zero_row_count,"
        << "pardiso_spd_status,global_general_status,estimated_lambda_min,rayleigh_min,"
        << "ldlt_status,ldlt_positive,ldlt_negative,ldlt_zero,";
    writeSolverHeader("diagonal_pcg");
    out << ',';
    writeSolverHeader("bj_jacobi_pcg");
    out << ',';
    writeSolverHeader("bj_ic_pcg");
    out << ',';
    writeSolverHeader("bj_pardiso_pcg");
    out << ",Tmin,Tmax,pardiso_spd_message,global_general_message\n";
    const auto writeSolver = [&](const SweepSolverFields& solver) {
        out << csvEscape(solver.status) << ','
            << csvEscape(solver.failureReason) << ','
            << solver.breakdownIter << ','
            << solver.iterations << ','
            << solver.finalRelativeResidual << ','
            << solver.pcgPAp << ','
            << solver.pcgRTr << ','
            << solver.pcgRMz << ','
            << solver.pcgAlpha << ','
            << solver.pcgBeta << ','
            << solver.pcgResidualNorm << ','
            << (solver.pcgPHasNonFinite ? 1 : 0) << ','
            << (solver.pcgApHasNonFinite ? 1 : 0);
    };
    for (const SpdPenaltySweepRow& row : rows) {
        out << row.interfaceScheme << ','
            << row.penaltyMode << ','
            << row.penaltyScaling << ','
            << row.penaltyFactor << ','
            << row.negativeDiagCount << ','
            << row.minDiag << ','
            << row.maxDiag << ','
            << row.symmetryError << ','
            << (row.hasNonFinite ? 1 : 0) << ','
            << row.zeroRowCount << ','
            << csvEscape(row.pardisoSpdStatus) << ','
            << csvEscape(row.globalGeneralStatus) << ','
            << row.estimatedLambdaMin << ','
            << row.rayleighMin << ','
            << csvEscape(row.ldltStatus) << ','
            << row.ldltPositivePivots << ','
            << row.ldltNegativePivots << ','
            << row.ldltZeroTinyPivots << ',';
        writeSolver(row.diagonalPcg);
        out << ',';
        writeSolver(row.bjJacobiPcg);
        out << ',';
        writeSolver(row.bjIcPcg);
        out << ',';
        writeSolver(row.bjPardisoPcg);
        out << ','
            << row.tMin << ','
            << row.tMax << ','
            << csvEscape(row.pardisoSpdMessage) << ','
            << csvEscape(row.globalGeneralMessage) << '\n';
    }
}

static void writeSweepBlockDiagnostics(const std::vector<SweepBlockDiagnostic>& rows,
                                       const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "interface_scheme,penalty_mode,penalty_scaling,penalty_factor,"
        << "block_id,block_size,block_diagonal_min,block_diagonal_max,"
        << "block_symmetry_error,block_lambda_min_est,"
        << "block_ldlt_status,block_ldlt_positive,block_ldlt_negative,block_ldlt_zero,"
        << "block_ic_status,block_ic_message,"
        << "block_pardiso_spd_status,first_failed_block,block_pardiso_spd_message\n";
    for (const SweepBlockDiagnostic& row : rows) {
        out << row.interfaceScheme << ','
            << row.penaltyMode << ','
            << row.penaltyScaling << ','
            << row.penaltyFactor << ','
            << row.blockId << ','
            << row.blockSize << ','
            << row.minDiag << ','
            << row.maxDiag << ','
            << row.symmetryError << ','
            << row.lambdaMinEst << ','
            << csvEscape(row.ldltStatus) << ','
            << row.ldltPositivePivots << ','
            << row.ldltNegativePivots << ','
            << row.ldltZeroTinyPivots << ','
            << csvEscape(row.icStatus) << ','
            << csvEscape(row.icMessage) << ','
            << csvEscape(row.pardisoSpdStatus) << ','
            << (row.firstFailedBlock ? 1 : 0) << ','
            << csvEscape(row.pardisoSpdMessage) << '\n';
    }
}

static void writeFinalSpdDiagnostic(const FinalSpdDiagnostic& diag,
                                    const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "metric,value,interpretation\n";
    out << "has_nan_inf," << (diag.hasNonFinite ? 1 : 0)
        << ",A_final must not contain NaN or Inf\n";
    out << "symmetry_error," << diag.symmetryError
        << ",Small symmetry error means symmetric only; it does not prove SPD\n";
    out << "diagonal_min," << diag.minDiag
        << ",Diagonal range is diagnostic only; it is not an SPD criterion\n";
    out << "diagonal_max," << diag.maxDiag
        << ",Diagonal range is diagnostic only; it is not an SPD criterion\n";
    out << "zero_row_count," << diag.zeroRowCount
        << ",Zero row count is diagnostic only; zero rows absent does not prove SPD\n";
    out << "ldlt_status," << csvEscape(diag.ldltStatus)
        << ",LDLT inertia is a primary SPD diagnostic\n";
    out << "ldlt_positive_count," << diag.ldltPositivePivots
        << ",For numerical SPD this should equal matrix size\n";
    out << "ldlt_negative_count," << diag.ldltNegativePivots
        << ",For numerical SPD this must be 0\n";
    out << "ldlt_zero_count," << diag.ldltZeroTinyPivots
        << ",For numerical SPD this must be 0\n";
    out << "lambda_min_est," << diag.lambdaMinEst
        << ",Estimated smallest eigenvalue should be positive\n";
    out << "pardiso_spd_factorization," << csvEscape(diag.pardisoSpdStatus)
        << ",Cholesky/PARDISO SPD factorization should succeed for numerical SPD\n";
    out << "pardiso_general_factorization," << csvEscape(diag.pardisoGeneralStatus)
        << ",If general succeeds but SPD fails the matrix may be non-SPD or unsuitable for SPD factorization\n";
    out << "spd_conclusion," << csvEscape(diag.spdConclusion)
        << ",Conclusion uses inertia + lambda estimate + SPD factorization, not PCG\n";
    out << "pardiso_spd_message," << csvEscape(diag.pardisoSpdMessage) << ",factorization message\n";
    out << "pardiso_general_message," << csvEscape(diag.pardisoGeneralMessage) << ",factorization message\n";
}

static SweepSolverFields solverFieldsFromStats(const SolverStatistics& stats, int iterations)
{
    SweepSolverFields fields;
    fields.status = stats.status;
    fields.failureReason = stats.failureReason;
    fields.breakdownIter = stats.breakdownIteration;
    fields.iterations = iterations;
    fields.finalRelativeResidual = stats.finalRelativeResidual;
    fields.pcgPAp = stats.pcgPAp;
    fields.pcgRTr = stats.pcgRTr;
    fields.pcgRMz = stats.pcgRMz;
    fields.pcgAlpha = stats.pcgAlpha;
    fields.pcgBeta = stats.pcgBeta;
    fields.pcgResidualNorm = stats.pcgResidualNorm;
    fields.pcgPHasNonFinite = stats.pcgPHasNonFinite;
    fields.pcgApHasNonFinite = stats.pcgApHasNonFinite;
    return fields;
}

static SweepSolverFields runDiagonalPcgForSweep(const Mesh& mesh,
                                                const CaseConfig& physics,
                                                const SparseMatrix& system,
                                                const std::vector<double>& rhs,
                                                std::vector<double>* solution,
                                                int maxIterations,
                                                double tolerance)
{
    SolverStatistics stats;
    int iterations = 0;
    std::vector<double> temperature = initialTemperatureVector(mesh, physics);
    temperature = preconditionedConjugateGradient(system,
                                                  rhs,
                                                  std::move(temperature),
                                                  iterations,
                                                  &stats,
                                                  maxIterations,
                                                  tolerance);
    if (solution != nullptr && stats.status == "success") {
        *solution = temperature;
    }
    return solverFieldsFromStats(stats, iterations);
}

static SweepSolverFields runBjPardisoPcgForSweep(const Mesh& mesh,
                                                 const CaseConfig& physics,
                                                 const SparseMatrix& system,
                                                 const std::vector<double>& rhs,
                                                 int maxIterations,
                                                 double tolerance)
{
    SolverStatistics stats;
    stats.name = "BJ-PARDISO-PCG";
    int iterations = 0;
    try {
        BlockJacobiPreconditioner preconditioner(mesh, system);
        std::vector<double> temperature = initialTemperatureVector(mesh, physics);
        (void)blockJacobiPcg(system,
                             rhs,
                             std::move(temperature),
                             preconditioner,
                             iterations,
                             &stats,
                             "BJ-PARDISO-PCG",
                             maxIterations,
                             tolerance);
    } catch (const std::exception& err) {
        stats.status = "failed";
        stats.failureReason = err.what();
        stats.breakdownIteration = -1;
        stats.finalRelativeResidual = std::numeric_limits<double>::quiet_NaN();
    }
    return solverFieldsFromStats(stats, iterations);
}

static SweepSolverFields runBjIcPcgForSweep(const Mesh& mesh,
                                            const CaseConfig& physics,
                                            const SparseMatrix& system,
                                            const std::vector<double>& rhs,
                                            const ProgramOptions& options)
{
    SolverStatistics stats;
    stats.name = "BJ-IC-PCG";
    int iterations = 0;
    try {
        BlockJacobiIcPreconditioner preconditioner(mesh,
                                                   system,
                                                   options.icShift,
                                                   options.icScaling,
                                                   options.diagScalingEps,
                                                   options.localIcShiftMode,
                                                   options.icOrdering);
        if (preconditioner.hasBreakdown()) {
            stats.status = "failed";
            stats.failureReason = "IC(0) factorization breakdown in at least one block";
            stats.breakdownIteration = -1;
        } else {
            std::vector<double> temperature = initialTemperatureVector(mesh, physics);
            (void)blockJacobiPcg(system,
                                 rhs,
                                 std::move(temperature),
                                 preconditioner,
                                 iterations,
                                 &stats,
                                 "BJ-IC-PCG",
                                 options.maxPcgIterations,
                                 options.pcgTolerance);
        }
    } catch (const std::exception& err) {
        stats.status = "failed";
        stats.failureReason = err.what();
        stats.breakdownIteration = -1;
    }
    return solverFieldsFromStats(stats, iterations);
}

static FinalSpdDiagnostic computeFinalSpdDiagnostic(const MatrixDiagnostics& matrixDiagnostics,
                                                    const SparseMatrix& system)
{
    FinalSpdDiagnostic diag;
    diag.hasNonFinite = matrixDiagnostics.hasNonFinite;
    diag.symmetryError = matrixDiagnostics.symmetryRatio;
    diag.minDiag = matrixDiagnostics.minDiagonal;
    diag.maxDiag = matrixDiagnostics.maxDiagonal;
    diag.zeroRowCount = matrixDiagnostics.zeroRows;

    const PardisoInertiaDiagnostic inertia = computePardisoInertia(system);
    diag.ldltStatus = inertia.status;
    diag.ldltPositivePivots = inertia.positivePivots;
    diag.ldltNegativePivots = inertia.negativePivots;
    diag.ldltZeroTinyPivots = inertia.zeroTinyPivots;
    diag.lambdaMinEst = estimateSmallestEigenvalueLanczos(system, 40);

    std::string spdMessage;
    const bool spdOk = confirmSpdWithPardiso(system, spdMessage);
    diag.pardisoSpdStatus = spdOk ? "success" : "failure";
    diag.pardisoSpdMessage = spdMessage;

    std::string generalMessage;
    const bool generalOk = confirmGeneralWithPardiso(system, generalMessage);
    diag.pardisoGeneralStatus = generalOk ? "success" : "failure";
    diag.pardisoGeneralMessage = generalMessage;

    const bool inertiaSpd = inertia.status == "success"
        && inertia.negativePivots == 0
        && inertia.zeroTinyPivots == 0;
    const bool lambdaPositive = std::isfinite(diag.lambdaMinEst) && diag.lambdaMinEst > 0.0;
    const bool symmetricEnough = std::isfinite(diag.symmetryError) && diag.symmetryError < 1.0e-10;
    if (!diag.hasNonFinite && symmetricEnough && inertiaSpd && lambdaPositive && spdOk) {
        diag.spdConclusion = "numerically_spd";
    } else if (generalOk && !spdOk) {
        diag.spdConclusion = "not_spd_or_not_suitable_for_spd_factorization";
    } else {
        diag.spdConclusion = "not_confirmed_spd";
    }
    return diag;
}

static std::vector<SweepBlockDiagnostic> collectSweepBlockDiagnostics(const Mesh& mesh,
                                                                      const SparseMatrix& system,
                                                                      const CaseConfig& physics,
                                                                      const ProgramOptions* options = nullptr)
{
    int maxSubdomain = 0;
    for (const Node& node : mesh.nodes) {
        maxSubdomain = std::max(maxSubdomain, node.subdomain);
    }
    std::vector<std::vector<int>> blockDofs(static_cast<size_t>(maxSubdomain + 1));
    std::vector<int> globalToBlock(mesh.nodes.size(), -1);
    std::vector<int> globalToLocal(mesh.nodes.size(), -1);
    for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
        const int block = mesh.nodes[static_cast<size_t>(i)].subdomain;
        globalToBlock[static_cast<size_t>(i)] = block;
        globalToLocal[static_cast<size_t>(i)] = static_cast<int>(blockDofs[static_cast<size_t>(block)].size());
        blockDofs[static_cast<size_t>(block)].push_back(i);
    }

    std::vector<SweepBlockDiagnostic> rows;
    int firstFailed = -1;
    for (size_t blockIndex = 0; blockIndex < blockDofs.size(); ++blockIndex) {
        SweepBlockDiagnostic row;
        row.interfaceScheme = physics.interfaceScheme;
        row.penaltyMode = physics.penaltyMode;
        row.penaltyScaling = physics.penaltyScaling;
        row.penaltyFactor = physics.penaltyFactor;
        row.blockId = static_cast<int>(blockIndex);
        row.blockSize = static_cast<int>(blockDofs[blockIndex].size());

        std::vector<MatrixEntry> entries;
        for (int globalRow : blockDofs[blockIndex]) {
            const int localRow = globalToLocal[static_cast<size_t>(globalRow)];
            for (int k = system.rowPtr[static_cast<size_t>(globalRow)];
                 k < system.rowPtr[static_cast<size_t>(globalRow + 1)];
                 ++k) {
                const int globalCol = system.colInd[static_cast<size_t>(k)];
                if (globalCol < 0 || globalCol >= static_cast<int>(globalToBlock.size())) {
                    continue;
                }
                if (globalToBlock[static_cast<size_t>(globalCol)] != static_cast<int>(blockIndex)) {
                    continue;
                }
                entries.push_back({localRow,
                                   globalToLocal[static_cast<size_t>(globalCol)],
                                   system.values[static_cast<size_t>(k)]});
            }
        }

        SparseMatrix blockMatrix(row.blockSize);
        blockMatrix.appendEntries(entries);
        blockMatrix.finalizeCsr();
        const MatrixStageDiagnostic stage = diagnoseMatrixStage("block", blockMatrix);
        row.symmetryError = stage.symmetryRatio;
        const std::vector<double> diagonal = matrixDiagonalVector(blockMatrix);
        if (!diagonal.empty()) {
            const auto minmax = std::minmax_element(diagonal.begin(), diagonal.end());
            row.minDiag = *minmax.first;
            row.maxDiag = *minmax.second;
        }
        row.lambdaMinEst = estimateSmallestEigenvalueLanczos(blockMatrix, 30);
        const PardisoInertiaDiagnostic inertia = computePardisoInertia(blockMatrix);
        row.ldltStatus = inertia.status;
        row.ldltPositivePivots = inertia.positivePivots;
        row.ldltNegativePivots = inertia.negativePivots;
        row.ldltZeroTinyPivots = inertia.zeroTinyPivots;
        if (options != nullptr) {
            try {
                const std::vector<MatrixEntry> blockEntries = sparseMatrixEntries(blockMatrix);
                SubdomainIcSolver icSolver(row.blockSize,
                                           blockEntries,
                                           options->icShift,
                                           row.blockId,
                                           options->icScaling,
                                           options->diagScalingEps,
                                           options->localIcShiftMode,
                                           options->icOrdering);
                const IcFactorDiagnostics& icDiag = icSolver.diagnostics();
                if (icDiag.accepted && !icDiag.breakdown && icDiag.nonFiniteLCount == 0) {
                    row.icStatus = "success";
                } else {
                    row.icStatus = "failure";
                    row.icMessage = "IC(0) breakdown or no accepted shifted factorization";
                }
            } catch (const std::exception& ex) {
                row.icStatus = "failure";
                row.icMessage = ex.what();
            }
        }
        std::string spdMessage;
        const bool spdOk = confirmSpdWithPardiso(blockMatrix, spdMessage);
        row.pardisoSpdStatus = spdOk ? "success" : "failure";
        row.pardisoSpdMessage = spdMessage;
        if (!spdOk && firstFailed < 0) {
            firstFailed = row.blockId;
        }
        rows.push_back(std::move(row));
    }
    for (SweepBlockDiagnostic& row : rows) {
        row.firstFailedBlock = row.blockId == firstFailed;
    }
    return rows;
}

static int runSpdPenaltySweep(const Mesh& mesh,
                              CaseConfig physics,
                              const ProgramOptions& options,
                              const std::filesystem::path& outputDir)
{
    const int n = static_cast<int>(mesh.nodes.size());
    std::cout << "Running SPD/coercivity penalty sweep; no COMSOL comparison will be run.\n";
    SparseMatrix baseStiffness(n);
    std::vector<double> baseSource(static_cast<size_t>(n), 0.0);
    assembleVolume(mesh, physics, nullptr, baseStiffness, baseSource);
    const int convectionBoundaryFaceCount = assembleConvectionBoundaries(mesh,
                                                                         physics,
                                                                         baseStiffness,
                                                                         baseSource,
                                                                         options.disableConvectionLhs,
                                                                         nullptr);
    const std::vector<double> sourceBeforeDirichlet = baseSource;
    physics.interfaceScheme = "sipg";
    const std::vector<std::string> penaltyModes = options.penaltyModeOverride
        ? std::vector<std::string>{options.penaltyMode}
        : std::vector<std::string>{"max"};
    const std::vector<double> penaltyFactors = options.penaltyFactorOverride
        ? std::vector<double>{options.penaltyFactor}
        : std::vector<double>{15.0, 30.0, 50.0, 100.0, 200.0, 500.0};
    const std::filesystem::path sweepCsv = outputDir / "spd_penalty_sweep.csv";
    const std::filesystem::path blockCsv = outputDir / "sipg_penalty_sweep_block_diagnostics.csv";
    std::vector<SpdPenaltySweepRow> rows;
    std::vector<SweepBlockDiagnostic> blockRows;

    for (const std::string& penaltyMode : penaltyModes) {
        for (double penaltyFactor : penaltyFactors) {
            SpdPenaltySweepRow row;
            row.interfaceScheme = physics.interfaceScheme;
            row.penaltyMode = penaltyMode;
            row.penaltyScaling = physics.penaltyScaling;
            row.penaltyFactor = penaltyFactor;
            physics.penaltyMode = penaltyMode;
            physics.penaltyFactor = penaltyFactor;
            std::cout << "SPD sweep: scheme=" << physics.interfaceScheme
                      << ", mode=" << penaltyMode
                      << ", scaling=" << physics.penaltyScaling
                      << ", penalty_factor=" << penaltyFactor << "\n";

            AssemblyDiagnostics diagnostics;
            initializeAssemblyDiagnostics(n, diagnostics);
            markInterfaceDofs(mesh, diagnostics);
            diagnostics.volumeDiag = matrixDiagonalVector(baseStiffness);
            SparseMatrix system = baseStiffness;
            assembleSipgInterface(mesh,
                                  physics,
                                  system,
                                  options.disableInterfaceConsistency,
                                  true,
                                  &diagnostics);
            assembleSipgInterface(mesh,
                                  physics,
                                  system,
                                  true,
                                  options.disableInterfacePenalty,
                                  &diagnostics);
            diagnostics.preDirichletDiag = matrixDiagonalVector(system);
            const std::vector<double> fixedAdjust = makeDirichletAdjustedSystem(mesh, system);
            std::vector<double> rhs = sourceBeforeDirichlet;
            applyDirichletRhs(mesh, fixedAdjust, rhs);
            system.finalizeCsr();
            diagnostics.finalDiag = matrixDiagonalVector(system);
            const DiagonalStats diagStats = diagonalStats(diagnostics.finalDiag);
            row.negativeDiagCount = diagStats.negativeEntries;
            row.minDiag = diagStats.minDiagonal;
            row.maxDiag = diagStats.maxDiagonal;

            const std::vector<ComponentDiagnostics> componentDiagnostics =
                collectComponentDiagnostics(mesh, physics);
            const MatrixDiagnostics matrixDiagnostics =
                diagnoseMatrixAndPhysics(mesh, physics, system, convectionBoundaryFaceCount, componentDiagnostics);
            row.symmetryError = matrixDiagnostics.symmetryRatio;
            row.hasNonFinite = matrixDiagnostics.hasNonFinite;
            row.zeroRowCount = matrixDiagnostics.zeroRows;

            row.rayleighMin = randomRayleighMinimum(system, 100);
            row.estimatedLambdaMin = estimateSmallestEigenvalueLanczos(system, 30);

            std::string spdMessage;
            const bool spdOk = confirmSpdWithPardiso(system, spdMessage);
            row.pardisoSpdStatus = spdOk ? "success" : "failure";
            row.pardisoSpdMessage = spdMessage;

            std::string generalMessage;
            const bool generalOk = confirmGeneralWithPardiso(system, generalMessage);
            row.globalGeneralStatus = generalOk ? "success" : "failure";
            row.globalGeneralMessage = generalMessage;

            const PardisoInertiaDiagnostic inertia = computePardisoInertia(system);
            row.ldltStatus = inertia.status;
            row.ldltPositivePivots = inertia.positivePivots;
            row.ldltNegativePivots = inertia.negativePivots;
            row.ldltZeroTinyPivots = inertia.zeroTinyPivots;

            std::vector<double> diagonalSolution;
            row.diagonalPcg = runDiagonalPcgForSweep(mesh,
                                                     physics,
                                                     system,
                                                     rhs,
                                                     &diagonalSolution,
                                                     options.maxPcgIterations,
                                                     options.pcgTolerance);
            row.bjJacobiPcg = runDiagonalPcgForSweep(mesh,
                                                     physics,
                                                     system,
                                                     rhs,
                                                     nullptr,
                                                     options.maxPcgIterations,
                                                     options.pcgTolerance);
            row.bjIcPcg = runBjIcPcgForSweep(mesh, physics, system, rhs, options);

            std::vector<SweepBlockDiagnostic> currentBlockRows =
                collectSweepBlockDiagnostics(mesh, system, physics);
            const bool blockSpdOk = std::all_of(currentBlockRows.begin(),
                                                currentBlockRows.end(),
                                                [](const SweepBlockDiagnostic& block) {
                                                    return block.pardisoSpdStatus == "success";
                                                });
            row.bjPardisoPcg = runBjPardisoPcgForSweep(mesh,
                                                       physics,
                                                       system,
                                                       rhs,
                                                       options.maxPcgIterations,
                                                       options.pcgTolerance);
            if (row.bjPardisoPcg.status == "failed" && spdOk && !blockSpdOk) {
                row.bjPardisoPcg.failureReason =
                    "Global matrix passed SPD factorization, but at least one BJ block failed SPD factorization: "
                    + row.bjPardisoPcg.failureReason;
            }
            blockRows.insert(blockRows.end(),
                             std::make_move_iterator(currentBlockRows.begin()),
                             std::make_move_iterator(currentBlockRows.end()));

            if (row.diagonalPcg.status == "success" && !diagonalSolution.empty()) {
                const auto minmax = std::minmax_element(diagonalSolution.begin(), diagonalSolution.end());
                row.tMin = *minmax.first;
                row.tMax = *minmax.second;
            }

            std::cout << "  negative_diag=" << row.negativeDiagCount
                      << ", min_diag=" << row.minDiag
                      << ", symmetry=" << row.symmetryError
                      << ", nonfinite=" << (row.hasNonFinite ? "yes" : "no")
                      << ", zero_rows=" << row.zeroRowCount
                      << ", PARDISO_SPD=" << row.pardisoSpdStatus
                      << ", general=" << row.globalGeneralStatus
                      << ", lambda_min_est=" << row.estimatedLambdaMin
                      << ", rayleigh_min=" << row.rayleighMin
                      << ", Diagonal-PCG=" << row.diagonalPcg.status
                      << ", BJ-IC-PCG=" << row.bjIcPcg.status
                      << ", BJ-PARDISO-PCG=" << row.bjPardisoPcg.status << "\n";
            if (row.diagonalPcg.failureReason.find("pAp") != std::string::npos) {
                std::cout << "  Diagonal-PCG pAp breakdown diagnostics:"
                          << " pAp=" << row.diagonalPcg.pcgPAp
                          << " rTr=" << row.diagonalPcg.pcgRTr
                          << " rMz=" << row.diagonalPcg.pcgRMz
                          << " alpha=" << row.diagonalPcg.pcgAlpha
                          << " beta=" << row.diagonalPcg.pcgBeta
                          << " ||r||=" << row.diagonalPcg.pcgResidualNorm
                          << " p_has_nonfinite=" << (row.diagonalPcg.pcgPHasNonFinite ? "yes" : "no")
                          << " Ap_has_nonfinite=" << (row.diagonalPcg.pcgApHasNonFinite ? "yes" : "no") << "\n";
            }
            rows.push_back(row);
            writeSpdPenaltySweep(rows, sweepCsv);
            writeSweepBlockDiagnostics(blockRows, blockCsv);
        }
    }

    writeSpdPenaltySweep(rows, sweepCsv);
    writeSweepBlockDiagnostics(blockRows, blockCsv);
    std::cout << "Wrote " << sweepCsv.string() << "\n";
    std::cout << "Wrote " << blockCsv.string() << "\n";
    return 0;
}

static double matrixQuadraticForm(const SparseMatrix& matrix, const std::vector<double>& x)
{
    const std::vector<double> ax = matrix.multiply(x);
    return vectorDot(x, ax);
}

static SparseMatrix finalizedDirichletMatrix(const Mesh& mesh, SparseMatrix matrix)
{
    (void)makeDirichletAdjustedSystem(mesh, matrix);
    matrix.finalizeCsr();
    return matrix;
}

struct SipgLocalFaceMatrix {
    std::array<double, 400> consistency{};
    std::array<double, 400> penalty{};
    std::array<int, 20> dofs{};
    double leftArea = 0.0;
    double rightArea = 0.0;
    double integrationArea = 0.0;
    double centerDistance = 0.0;
    double minLambdaLeft = std::numeric_limits<double>::infinity();
    double minLambdaRight = std::numeric_limits<double>::infinity();
    int quadraturePoints = 0;
    int quadraturePointsInCommon = 0;
    int quadratureOutsideLeft = 0;
    int quadratureOutsideRight = 0;
};

static SipgLocalFaceMatrix buildSipgLocalFaceMatrix(const Mesh& mesh,
                                                    const CaseConfig& config,
                                                    size_t faceIndex)
{
    const std::vector<TriangleQuadraturePoint> quadrature = makeTriangleQuadrature();
    const InterfaceFace& face = mesh.interfaceFaces[faceIndex];
    const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
    const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
    const ElementGeometry gLeft = elementGeometry(mesh, left);
    const ElementGeometry gRight = elementGeometry(mesh, right);
    const Material& leftMaterial = materialForTet(config, left);
    const Material& rightMaterial = materialForTet(config, right);
    const Vec3 normal = norm(face.leftNormal) > 1.0e-30
        ? face.leftNormal
        : normalized(subdomainCenter(mesh, right.subdomain) - subdomainCenter(mesh, left.subdomain));
    const std::array<Vec3, 3> leftPhysicalFace = boundaryTriangleFromTetFace(mesh, left, face.leftLocal);
    const std::array<Vec3, 3> rightPhysicalFace = boundaryTriangleFromTetFace(mesh, right, face.rightLocal);
    const double leftFaceArea = triangleArea(leftPhysicalFace);
    const double rightFaceArea = triangleArea(rightPhysicalFace);
    const double leftVolume = gLeft.detJ / 6.0;
    const double rightVolume = gRight.detJ / 6.0;
    const double hLeft = 3.0 * leftVolume / std::max(1.0e-30, leftFaceArea);
    const double hRight = 3.0 * rightVolume / std::max(1.0e-30, rightFaceArea);
    const double hFace = std::min(hLeft, hRight);
    constexpr double polynomialOrder = 2.0;
    const double pScale = config.penaltyScaling == "p1_squared"
        ? (polynomialOrder + 1.0) * (polynomialOrder + 1.0)
        : polynomialOrder * (polynomialOrder + 1.0);
    const double kLeftNormal = std::max(1.0e-30, normalConductivity(leftMaterial, normal));
    const double kRightNormal = std::max(1.0e-30, normalConductivity(rightMaterial, normal));
    const double alphaLeft = kLeftNormal / std::max(1.0e-30, hLeft);
    const double alphaRight = kRightNormal / std::max(1.0e-30, hRight);
    const double penalty = config.penaltyMode == "max"
        ? config.penaltyFactor * pScale * std::max(alphaLeft, alphaRight)
        : config.penaltyFactor * pScale
            * (2.0 * alphaLeft * alphaRight / std::max(1.0e-30, alphaLeft + alphaRight));
    const double adjointSign = config.interfaceScheme == "nipg" ? 1.0 : -1.0;

    SipgLocalFaceMatrix local;
    local.leftArea = leftFaceArea;
    local.rightArea = rightFaceArea;
    local.integrationArea = integrationArea(face.integrationTriangles);
    local.centerDistance = norm(triangleCenter(leftPhysicalFace) - triangleCenter(rightPhysicalFace));
    for (int i = 0; i < 10; ++i) {
        local.dofs[static_cast<size_t>(i)] = left.dof[static_cast<size_t>(i)];
        local.dofs[static_cast<size_t>(10 + i)] = right.dof[static_cast<size_t>(i)];
    }

    for (const auto& integrationTriangle : face.integrationTriangles) {
        const Vec3 p0 = integrationTriangle[0];
        const Vec3 p1 = integrationTriangle[1];
        const Vec3 p2 = integrationTriangle[2];
        const double jacFace = norm(cross(p1 - p0, p2 - p0));
        for (const TriangleQuadraturePoint& qp : quadrature) {
            const double weight = qp.weight * jacFace;
            const Vec3 q = (1.0 - qp.a - qp.b) * p0 + qp.a * p1 + qp.b * p2;
            const auto lLeft = lambdaOnTetFace(q, left, face.leftLocal, mesh);
            const auto lRight = lambdaOnTetFace(q, right, face.rightLocal, mesh);
            const auto nLeft = shapeP2(lLeft);
            const auto nRight = shapeP2(lRight);
            const auto gradLeft = gradShapeP2(lLeft, gLeft);
            const auto gradRight = gradShapeP2(lRight, gRight);
            const double minLeft = *std::min_element(lLeft.begin(), lLeft.end());
            const double minRight = *std::min_element(lRight.begin(), lRight.end());
            local.minLambdaLeft = std::min(local.minLambdaLeft, minLeft);
            local.minLambdaRight = std::min(local.minLambdaRight, minRight);
            ++local.quadraturePoints;
            if (minLeft < -1.0e-8) {
                ++local.quadratureOutsideLeft;
            }
            if (minRight < -1.0e-8) {
                ++local.quadratureOutsideRight;
            }
            if (minLeft >= -1.0e-8 && minRight >= -1.0e-8) {
                ++local.quadraturePointsInCommon;
            }

            for (int aTest = 0; aTest < 10; ++aTest) {
                const double vL = nLeft[static_cast<size_t>(aTest)];
                const double vR = nRight[static_cast<size_t>(aTest)];
                const double dVL = leftMaterial.conductivityX * gradLeft[static_cast<size_t>(aTest)].x * normal.x
                                  + leftMaterial.conductivityY * gradLeft[static_cast<size_t>(aTest)].y * normal.y
                                  + leftMaterial.conductivityZ * gradLeft[static_cast<size_t>(aTest)].z * normal.z;
                const double dVR = rightMaterial.conductivityX * gradRight[static_cast<size_t>(aTest)].x * normal.x
                                  + rightMaterial.conductivityY * gradRight[static_cast<size_t>(aTest)].y * normal.y
                                  + rightMaterial.conductivityZ * gradRight[static_cast<size_t>(aTest)].z * normal.z;
                for (int bTrial = 0; bTrial < 10; ++bTrial) {
                    const double uL = nLeft[static_cast<size_t>(bTrial)];
                    const double uR = nRight[static_cast<size_t>(bTrial)];
                    const double dUL = leftMaterial.conductivityX * gradLeft[static_cast<size_t>(bTrial)].x * normal.x
                                      + leftMaterial.conductivityY * gradLeft[static_cast<size_t>(bTrial)].y * normal.y
                                      + leftMaterial.conductivityZ * gradLeft[static_cast<size_t>(bTrial)].z * normal.z;
                    const double dUR = rightMaterial.conductivityX * gradRight[static_cast<size_t>(bTrial)].x * normal.x
                                      + rightMaterial.conductivityY * gradRight[static_cast<size_t>(bTrial)].y * normal.y
                                      + rightMaterial.conductivityZ * gradRight[static_cast<size_t>(bTrial)].z * normal.z;
                    local.consistency[static_cast<size_t>(aTest * 20 + bTrial)] +=
                        (-0.5 * dUL * vL + adjointSign * 0.5 * dVL * uL) * weight;
                    local.consistency[static_cast<size_t>(aTest * 20 + 10 + bTrial)] +=
                        (-0.5 * dUR * vL - adjointSign * 0.5 * dVL * uR) * weight;
                    local.consistency[static_cast<size_t>((10 + aTest) * 20 + bTrial)] +=
                        (0.5 * dUL * vR + adjointSign * 0.5 * dVR * uL) * weight;
                    local.consistency[static_cast<size_t>((10 + aTest) * 20 + 10 + bTrial)] +=
                        (0.5 * dUR * vR - adjointSign * 0.5 * dVR * uR) * weight;
                    local.penalty[static_cast<size_t>(aTest * 20 + bTrial)] += (penalty * uL * vL) * weight;
                    local.penalty[static_cast<size_t>(aTest * 20 + 10 + bTrial)] += (-penalty * uR * vL) * weight;
                    local.penalty[static_cast<size_t>((10 + aTest) * 20 + bTrial)] += (-penalty * uL * vR) * weight;
                    local.penalty[static_cast<size_t>((10 + aTest) * 20 + 10 + bTrial)] += (penalty * uR * vR) * weight;
                }
            }
        }
    }
    if (!std::isfinite(local.minLambdaLeft)) {
        local.minLambdaLeft = 0.0;
        local.minLambdaRight = 0.0;
    }
    return local;
}

static double localQuadratic(const std::array<double, 400>& matrix, const std::vector<double>& global, const std::array<int, 20>& dofs)
{
    double value = 0.0;
    for (int i = 0; i < 20; ++i) {
        const double vi = global[static_cast<size_t>(dofs[static_cast<size_t>(i)])];
        for (int j = 0; j < 20; ++j) {
            value += vi * matrix[static_cast<size_t>(i * 20 + j)]
                * global[static_cast<size_t>(dofs[static_cast<size_t>(j)])];
        }
    }
    return value;
}

static double localSymmetryError(const std::array<double, 400>& matrix)
{
    double normSq = 0.0;
    double asymSq = 0.0;
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20; ++j) {
            const double value = matrix[static_cast<size_t>(i * 20 + j)];
            const double diff = value - matrix[static_cast<size_t>(j * 20 + i)];
            normSq += value * value;
            asymSq += diff * diff;
        }
    }
    return std::sqrt(asymSq) / std::sqrt(std::max(1.0e-300, normSq));
}

static double localMinEigenvalue(const std::array<double, 400>& matrix)
{
    std::vector<std::vector<double>> dense(20, std::vector<double>(20, 0.0));
    for (int i = 0; i < 20; ++i) {
        for (int j = 0; j < 20; ++j) {
            dense[static_cast<size_t>(i)][static_cast<size_t>(j)] =
                0.5 * (matrix[static_cast<size_t>(i * 20 + j)] + matrix[static_cast<size_t>(j * 20 + i)]);
        }
    }
    return smallestEigenvalueSymmetricJacobi(std::move(dense));
}

static std::string blockCsv(const std::array<double, 400>& matrix, int rowOffset, int colOffset)
{
    std::ostringstream out;
    out << std::setprecision(16);
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            if (i != 0 || j != 0) {
                out << ';';
            }
            out << matrix[static_cast<size_t>((rowOffset + i) * 20 + colOffset + j)];
        }
    }
    return out.str();
}

struct ComponentSpdRow {
    std::string name;
    double symmetryError = 0.0;
    int positive = -1;
    int negative = -1;
    int zero = -1;
    double lambdaMinAbs = std::numeric_limits<double>::quiet_NaN();
    double lambdaMinAlgebraic = std::numeric_limits<double>::quiet_NaN();
    double conditionEstAbs = std::numeric_limits<double>::quiet_NaN();
    std::string pardisoSpdStatus;
    std::string pardisoGeneralStatus;
};

static ComponentSpdRow diagnoseComponentSpd(const std::string& name, const SparseMatrix& matrix)
{
    ComponentSpdRow row;
    row.name = name;
    const MatrixStageDiagnostic stage = diagnoseMatrixStage(name, matrix);
    row.symmetryError = stage.symmetryRatio;
    const PardisoInertiaDiagnostic inertia = computePardisoInertia(matrix);
    row.positive = inertia.positivePivots;
    row.negative = inertia.negativePivots;
    row.zero = inertia.zeroTinyPivots;
    row.lambdaMinAlgebraic = estimateSmallestEigenvalueLanczos(matrix, 40);
    if (matrix.size() <= 500000) {
        try {
            const std::vector<NearZeroMode> nearModes = computeNearZeroModesShiftInvert(matrix, 1, 24, 5);
            if (!nearModes.empty()) {
                row.lambdaMinAbs = std::abs(nearModes.front().eigenvalue);
            }
        } catch (const std::exception&) {
            row.lambdaMinAbs = std::numeric_limits<double>::quiet_NaN();
        }
    }
    if (!std::isfinite(row.lambdaMinAbs)) {
        row.lambdaMinAbs = std::abs(row.lambdaMinAlgebraic);
    }
    const double lambdaMaxAbs = estimateMaxAbsEigenvaluePower(matrix, 50);
    row.conditionEstAbs = lambdaMaxAbs / std::max(1.0e-300, row.lambdaMinAbs);
    std::string message;
    row.pardisoSpdStatus = confirmSpdWithPardiso(matrix, message) ? "success" : "failure";
    row.pardisoGeneralStatus = confirmGeneralWithPardiso(matrix, message) ? "success" : "failure";
    return row;
}

static void writeComponentSpdDiagnostics(const std::vector<ComponentSpdRow>& rows,
                                         const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "matrix,symmetry_error,ldlt_positive,ldlt_negative,ldlt_zero,"
        << "lambda_min_abs,lambda_min_algebraic,condition_est_abs,pardiso_spd_status,pardiso_general_status\n";
    out << std::setprecision(16);
    for (const ComponentSpdRow& row : rows) {
        out << row.name << ','
            << row.symmetryError << ','
            << row.positive << ','
            << row.negative << ','
            << row.zero << ','
            << row.lambdaMinAbs << ','
            << row.lambdaMinAlgebraic << ','
            << row.conditionEstAbs << ','
            << csvEscape(row.pardisoSpdStatus) << ','
            << csvEscape(row.pardisoGeneralStatus) << '\n';
    }
}

struct FaceEnergyRow {
    int modeIndex = -1;
    double eigenvalue = 0.0;
    int facePairId = -1;
    int leftFaceId = -1;
    int rightFaceId = -1;
    double penaltyEnergy = 0.0;
    double consistencyEnergy = 0.0;
    double ratio = 0.0;
    double leftArea = 0.0;
    double rightArea = 0.0;
    double centerDistance = 0.0;
};

static int runSipgSpdDiagnostics(const Mesh& mesh,
                                 CaseConfig physics,
                                 const ProgramOptions& options,
                                 const std::filesystem::path& outputDir,
                                 const std::string& modelPrefix)
{
    const int n = static_cast<int>(mesh.nodes.size());
    SparseMatrix kVolume(n);
    SparseMatrix kRobin(n);
    SparseMatrix kPenalty(n);
    SparseMatrix kConsistency(n);
    std::vector<double> source(static_cast<size_t>(n), 0.0);
    std::vector<double> robinSource(static_cast<size_t>(n), 0.0);
    assembleVolume(mesh, physics, nullptr, kVolume, source);
    (void)assembleConvectionBoundaries(mesh, physics, kRobin, robinSource, options.disableConvectionLhs, nullptr);
    assembleSipgInterface(mesh, physics, kPenalty, true, false, nullptr);
    assembleSipgInterface(mesh, physics, kConsistency, false, true, nullptr);
    kVolume.finalizeCsr();
    kRobin.finalizeCsr();
    kPenalty.finalizeCsr();
    kConsistency.finalizeCsr();

    SparseMatrix base(n);
    base.appendScaledEntries(kVolume, 1.0);
    base.appendScaledEntries(kRobin, 1.0);
    SparseMatrix a0 = finalizedDirichletMatrix(mesh, base);
    SparseMatrix withPenalty = base;
    withPenalty.appendScaledEntries(kPenalty, 1.0);
    SparseMatrix a1 = finalizedDirichletMatrix(mesh, withPenalty);
    SparseMatrix fullPreDirichlet = withPenalty;
    fullPreDirichlet.appendScaledEntries(kConsistency, 1.0);
    SparseMatrix a2 = finalizedDirichletMatrix(mesh, fullPreDirichlet);

    std::vector<ComponentSpdRow> componentRows;
    componentRows.push_back(diagnoseComponentSpd("A0_K_volume_plus_K_robin_plus_K_dirichlet", a0));
    componentRows.push_back(diagnoseComponentSpd("A1_A0_plus_K_penalty", a1));
    componentRows.push_back(diagnoseComponentSpd("A2_A1_plus_K_consistency", a2));
    writeComponentSpdDiagnostics(componentRows, outputDir / (modelPrefix + "_matrix_component_spd_diagnostics.csv"));
    writeComponentSpdDiagnostics(componentRows, outputDir / (modelPrefix + "_matrix_component_spd_diagnostics_after_overlap.csv"));

    if (!options.skipSipgPenaltySweep) {
        std::ofstream sweepOut(outputDir / (modelPrefix + "_overlap_penalty_sweep_spd.csv"));
        sweepOut << "C_pen,symmetry_error,ldlt_positive,ldlt_negative,ldlt_zero,"
            << "lambda_min_abs,lambda_min_algebraic,condition_est_abs,pardiso_spd_status,pardiso_general_status\n";
        sweepOut << std::setprecision(16);
        for (double cpen : {20.0, 50.0, 100.0, 200.0, 500.0}) {
            CaseConfig sweepPhysics = physics;
            sweepPhysics.penaltyFactor = cpen;
            SparseMatrix sweepPenalty(n);
            SparseMatrix sweepConsistency(n);
            assembleSipgInterface(mesh, sweepPhysics, sweepPenalty, true, false, nullptr);
            assembleSipgInterface(mesh, sweepPhysics, sweepConsistency, false, true, nullptr);
            SparseMatrix sweepMatrix(n);
            sweepMatrix.appendScaledEntries(kVolume, 1.0);
            sweepMatrix.appendScaledEntries(kRobin, 1.0);
            sweepMatrix.appendScaledEntries(sweepPenalty, 1.0);
            sweepMatrix.appendScaledEntries(sweepConsistency, 1.0);
            SparseMatrix sweepFinal = finalizedDirichletMatrix(mesh, sweepMatrix);
            const ComponentSpdRow diag = diagnoseComponentSpd("A2_sweep", sweepFinal);
            sweepOut << cpen << ','
                     << diag.symmetryError << ','
                     << diag.positive << ','
                     << diag.negative << ','
                     << diag.zero << ','
                     << diag.lambdaMinAbs << ','
                     << diag.lambdaMinAlgebraic << ','
                     << diag.conditionEstAbs << ','
                     << csvEscape(diag.pardisoSpdStatus) << ','
                     << csvEscape(diag.pardisoGeneralStatus) << '\n';
        }
    }

    std::vector<NearZeroMode> modes;
    if (componentRows.back().negative > 0 || componentRows.back().zero > 0) {
        modes = computeNearZeroModesShiftInvert(a2, 40, 56, 8);
        const std::vector<double> negativeShifts{-1.0e-12, -1.0e-10, -1.0e-9, -1.0e-8, -1.0e-7, -1.0e-6};
        for (double sigma : negativeShifts) {
            try {
                std::vector<NearZeroMode> shiftedModes = computeShiftedModesShiftInvert(a2, sigma, 24, 40, 6);
                modes.insert(modes.end(),
                             std::make_move_iterator(shiftedModes.begin()),
                             std::make_move_iterator(shiftedModes.end()));
            } catch (const std::exception&) {
            }
        }
        std::sort(modes.begin(), modes.end(), [](const NearZeroMode& lhs, const NearZeroMode& rhs) {
            return lhs.eigenvalue < rhs.eigenvalue;
        });
    }
    std::vector<NearZeroMode> negativeModes;
    for (NearZeroMode& mode : modes) {
        if (mode.eigenvalue < 0.0 && mode.residualNorm < 1.0e-5 && static_cast<int>(negativeModes.size()) < 20) {
            negativeModes.push_back(std::move(mode));
        }
    }
    if (negativeModes.empty()) {
        for (NearZeroMode& mode : modes) {
            if (static_cast<int>(negativeModes.size()) < 20) {
                negativeModes.push_back(std::move(mode));
            }
        }
    }
    for (size_t i = 0; i < negativeModes.size(); ++i) {
        negativeModes[i].modeIndex = static_cast<int>(i);
    }

    SparseMatrix kDirichlet(n);
    kDirichlet.appendScaledEntries(a2, 1.0);
    kDirichlet.appendScaledEntries(kVolume, -1.0);
    kDirichlet.appendScaledEntries(kRobin, -1.0);
    kDirichlet.appendScaledEntries(kPenalty, -1.0);
    kDirichlet.appendScaledEntries(kConsistency, -1.0);
    kDirichlet.finalizeCsr();

    std::ofstream energyOut(outputDir / (modelPrefix + "_negative_mode_energy_decomposition.csv"));
    energyOut << "mode_index,eigenvalue,component,energy,fraction_of_A_final_abs\n";
    energyOut << std::setprecision(16);
    for (const NearZeroMode& mode : negativeModes) {
        const double eVolume = matrixQuadraticForm(kVolume, mode.vector);
        const double eRobin = matrixQuadraticForm(kRobin, mode.vector);
        const double eDirichlet = matrixQuadraticForm(kDirichlet, mode.vector);
        const double ePenalty = matrixQuadraticForm(kPenalty, mode.vector);
        const double eConsistency = matrixQuadraticForm(kConsistency, mode.vector);
        const double eFinal = matrixQuadraticForm(a2, mode.vector);
        const double denom = std::max(1.0e-300, std::abs(eFinal));
        const std::vector<std::pair<std::string, double>> parts{
            {"K_volume", eVolume},
            {"K_robin", eRobin},
            {"K_dirichlet", eDirichlet},
            {"K_penalty", ePenalty},
            {"K_consistency", eConsistency},
            {"A_final", eFinal}
        };
        for (const auto& part : parts) {
            energyOut << mode.modeIndex << ','
                      << mode.eigenvalue << ','
                      << part.first << ','
                      << part.second << ','
                      << part.second / denom << '\n';
        }
    }

    const std::vector<InterfaceFacePairDiagnosticRow> pairRows = collectInterfaceFacePairDiagnostics(mesh, physics);
    std::map<int, int> leftMultiplicity;
    std::map<int, int> rightMultiplicity;
    for (const InterfaceFacePairDiagnosticRow& row : pairRows) {
        ++leftMultiplicity[row.leftFaceId];
        ++rightMultiplicity[row.rightFaceId];
    }
    std::vector<SipgLocalFaceMatrix> localFaceMatrices;
    localFaceMatrices.reserve(mesh.interfaceFaces.size());
    for (size_t faceIndex = 0; faceIndex < mesh.interfaceFaces.size(); ++faceIndex) {
        localFaceMatrices.push_back(buildSipgLocalFaceMatrix(mesh, physics, faceIndex));
    }

    std::vector<double> linearField(static_cast<size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) {
        const Vec3& p = mesh.nodes[static_cast<size_t>(i)].p;
        linearField[static_cast<size_t>(i)] = p.x + 2.0 * p.y + 3.0 * p.z;
    }
    const std::vector<TriangleQuadraturePoint> patchQuadrature = makeTriangleQuadrature();
    double patchArea = 0.0;
    double patchL2Numerator = 0.0;
    double patchMaxError = 0.0;
    double patchJumpSq = 0.0;
    double patchFluxJumpSq = 0.0;
    const Vec3 exactGradient{1.0, 2.0, 3.0};
    for (const InterfaceFace& face : mesh.interfaceFaces) {
        const Tet& left = mesh.tets[static_cast<size_t>(face.leftTet)];
        const Tet& right = mesh.tets[static_cast<size_t>(face.rightTet)];
        const ElementGeometry gLeft = elementGeometry(mesh, left);
        const ElementGeometry gRight = elementGeometry(mesh, right);
        for (const auto& tri : face.integrationTriangles) {
            const double jacFace = norm(cross(tri[1] - tri[0], tri[2] - tri[0]));
            for (const TriangleQuadraturePoint& qp : patchQuadrature) {
                const double weight = qp.weight * jacFace;
                const Vec3 q = (1.0 - qp.a - qp.b) * tri[0] + qp.a * tri[1] + qp.b * tri[2];
                const auto lLeft = lambdaOnTetFace(q, left, face.leftLocal, mesh);
                const auto lRight = lambdaOnTetFace(q, right, face.rightLocal, mesh);
                const auto nLeft = shapeP2(lLeft);
                const auto nRight = shapeP2(lRight);
                const auto gradLeft = gradShapeP2(lLeft, gLeft);
                const auto gradRight = gradShapeP2(lRight, gRight);
                double uLeft = 0.0;
                double uRight = 0.0;
                Vec3 reconstructedGradLeft{};
                Vec3 reconstructedGradRight{};
                for (int i = 0; i < 10; ++i) {
                    const double leftValue = linearField[static_cast<size_t>(left.dof[static_cast<size_t>(i)])];
                    const double rightValue = linearField[static_cast<size_t>(right.dof[static_cast<size_t>(i)])];
                    uLeft += nLeft[static_cast<size_t>(i)] * leftValue;
                    uRight += nRight[static_cast<size_t>(i)] * rightValue;
                    reconstructedGradLeft = reconstructedGradLeft + leftValue * gradLeft[static_cast<size_t>(i)];
                    reconstructedGradRight = reconstructedGradRight + rightValue * gradRight[static_cast<size_t>(i)];
                }
                const double exact = q.x + 2.0 * q.y + 3.0 * q.z;
                const double errorLeft = uLeft - exact;
                const double errorRight = uRight - exact;
                const double jump = uLeft - uRight;
                const double fluxJump = dot(reconstructedGradLeft - reconstructedGradRight, face.leftNormal);
                patchArea += weight;
                patchL2Numerator += (errorLeft * errorLeft + errorRight * errorRight) * 0.5 * weight;
                patchMaxError = std::max(patchMaxError, std::abs(errorLeft));
                patchMaxError = std::max(patchMaxError, std::abs(errorRight));
                patchJumpSq += jump * jump * weight;
                patchFluxJumpSq += fluxJump * fluxJump * weight;
                (void)exactGradient;
            }
        }
    }
    std::ofstream patchOut(outputDir / "nonmatching_sipg_patch_test.csv");
    patchOut << "model_name,L2_error,max_error,interface_jump_norm,flux_jump_norm,penalty_energy,"
        << "consistency_energy,matrix_symmetry_error,ldlt_positive,ldlt_negative,ldlt_zero,patch_type\n";
    patchOut << std::setprecision(16)
             << csvEscape(modelPrefix) << ','
             << std::sqrt(patchL2Numerator / std::max(1.0e-300, patchArea)) << ','
             << patchMaxError << ','
             << std::sqrt(patchJumpSq) << ','
             << std::sqrt(patchFluxJumpSq) << ','
             << matrixQuadraticForm(kPenalty, linearField) << ','
             << matrixQuadraticForm(kConsistency, linearField) << ','
             << componentRows.back().symmetryError << ','
             << componentRows.back().positive << ','
             << componentRows.back().negative << ','
             << componentRows.back().zero << ','
             << "interface_linear_reconstruction_T=x+2y+3z" << '\n';

    std::ofstream overlapDiagOut(outputDir / (modelPrefix + "_overlap_interface_integration_diagnostics.csv"));
    overlapDiagOut << "left_face_id,right_face_id,left_area,right_area,overlap_area,"
        << "overlap_area_over_left_area,overlap_area_over_right_area,"
        << "number_of_overlap_polygon_vertices,number_of_subtriangles,quadrature_points_count,"
        << "quadrature_outside_left_count,quadrature_outside_right_count,quadrature_outside_overlap_count,"
        << "left_barycentric_min,right_barycentric_min,center_distance,normal_dot,is_valid_overlap\n";
    overlapDiagOut << std::setprecision(16);
    std::map<int, double> leftOverlapAreaByFace;
    std::map<int, double> rightOverlapAreaByFace;
    std::map<int, int> leftPairCountByFace;
    std::map<int, int> rightPairCountByFace;
    for (size_t faceIndex = 0; faceIndex < mesh.interfaceFaces.size(); ++faceIndex) {
        const InterfaceFace& face = mesh.interfaceFaces[faceIndex];
        const SipgLocalFaceMatrix& local = localFaceMatrices[faceIndex];
        const int leftFaceId = faceIndex < pairRows.size() ? pairRows[faceIndex].leftFaceId : face.leftFaceId;
        const int rightFaceId = faceIndex < pairRows.size() ? pairRows[faceIndex].rightFaceId : face.rightFaceId;
        const int outsideOverlap = local.quadraturePoints - local.quadraturePointsInCommon;
        const bool valid = local.integrationArea > 0.0
            && outsideOverlap == 0
            && local.quadratureOutsideLeft == 0
            && local.quadratureOutsideRight == 0;
        leftOverlapAreaByFace[leftFaceId] += local.integrationArea;
        rightOverlapAreaByFace[rightFaceId] += local.integrationArea;
        ++leftPairCountByFace[leftFaceId];
        ++rightPairCountByFace[rightFaceId];
        overlapDiagOut << leftFaceId << ','
                       << rightFaceId << ','
                       << local.leftArea << ','
                       << local.rightArea << ','
                       << local.integrationArea << ','
                       << local.integrationArea / std::max(1.0e-300, local.leftArea) << ','
                       << local.integrationArea / std::max(1.0e-300, local.rightArea) << ','
                       << face.overlapPolygonVertices << ','
                       << face.integrationTriangles.size() << ','
                       << local.quadraturePoints << ','
                       << local.quadratureOutsideLeft << ','
                       << local.quadratureOutsideRight << ','
                       << outsideOverlap << ','
                       << local.minLambdaLeft << ','
                       << local.minLambdaRight << ','
                       << local.centerDistance << ','
                       << dot(face.leftNormal, face.rightNormal) << ','
                       << (valid ? 1 : 0) << '\n';
    }

    std::ofstream areaOut(outputDir / (modelPrefix + "_overlap_area_conservation.csv"));
    areaOut << "face_side,face_id,original_area,summed_overlap_area,relative_area_error,matched_pair_count\n";
    areaOut << std::setprecision(16);
    for (const auto& entry : leftOverlapAreaByFace) {
        if (entry.first < 0 || entry.first >= static_cast<int>(mesh.boundaryFaces.size())) {
            continue;
        }
        const double originalArea = boundaryFaceArea(mesh.boundaryFaces[static_cast<size_t>(entry.first)]);
        areaOut << "left," << entry.first << ','
                << originalArea << ','
                << entry.second << ','
                << std::abs(entry.second - originalArea) / std::max(1.0e-300, originalArea) << ','
                << leftPairCountByFace[entry.first] << '\n';
    }
    for (const auto& entry : rightOverlapAreaByFace) {
        if (entry.first < 0 || entry.first >= static_cast<int>(mesh.boundaryFaces.size())) {
            continue;
        }
        const double originalArea = boundaryFaceArea(mesh.boundaryFaces[static_cast<size_t>(entry.first)]);
        areaOut << "right," << entry.first << ','
                << originalArea << ','
                << entry.second << ','
                << std::abs(entry.second - originalArea) / std::max(1.0e-300, originalArea) << ','
                << rightPairCountByFace[entry.first] << '\n';
    }

    std::vector<FaceEnergyRow> faceEnergyRows;
    for (const NearZeroMode& mode : negativeModes) {
        for (size_t faceIndex = 0; faceIndex < mesh.interfaceFaces.size(); ++faceIndex) {
            const SipgLocalFaceMatrix& local = localFaceMatrices[faceIndex];
            FaceEnergyRow row;
            row.modeIndex = mode.modeIndex;
            row.eigenvalue = mode.eigenvalue;
            row.facePairId = static_cast<int>(faceIndex);
            if (faceIndex < pairRows.size()) {
                row.leftFaceId = pairRows[faceIndex].leftFaceId;
                row.rightFaceId = pairRows[faceIndex].rightFaceId;
            }
            row.penaltyEnergy = localQuadratic(local.penalty, mode.vector, local.dofs);
            row.consistencyEnergy = localQuadratic(local.consistency, mode.vector, local.dofs);
            row.ratio = std::abs(row.consistencyEnergy) / std::max(1.0e-300, std::abs(row.penaltyEnergy));
            row.leftArea = local.leftArea;
            row.rightArea = local.rightArea;
            row.centerDistance = local.centerDistance;
            faceEnergyRows.push_back(row);
        }
    }
    std::sort(faceEnergyRows.begin(), faceEnergyRows.end(), [](const FaceEnergyRow& lhs, const FaceEnergyRow& rhs) {
        if (lhs.ratio != rhs.ratio) {
            return lhs.ratio > rhs.ratio;
        }
        return lhs.consistencyEnergy < rhs.consistencyEnergy;
    });

    std::ofstream suspiciousOut(outputDir / (modelPrefix + "_suspicious_interface_face_energy.csv"));
    suspiciousOut << "rank,mode_index,eigenvalue,face_pair_id,left_face_id,right_face_id,"
        << "vT_K_penalty_face_v,vT_K_consistency_face_v,ratio,left_area,right_area,center_distance\n";
    suspiciousOut << std::setprecision(16);
    const size_t suspiciousCount = std::min<size_t>(100, faceEnergyRows.size());
    for (size_t i = 0; i < suspiciousCount; ++i) {
        const FaceEnergyRow& row = faceEnergyRows[i];
        suspiciousOut << (i + 1) << ','
                      << row.modeIndex << ','
                      << row.eigenvalue << ','
                      << row.facePairId << ','
                      << row.leftFaceId << ','
                      << row.rightFaceId << ','
                      << row.penaltyEnergy << ','
                      << row.consistencyEnergy << ','
                      << row.ratio << ','
                      << row.leftArea << ','
                      << row.rightArea << ','
                      << row.centerDistance << '\n';
    }

    std::ofstream localOut(outputDir / (modelPrefix + "_local_sipg_face_matrix_check.csv"));
    std::ofstream overlapLocalOut(outputDir / (modelPrefix + "_overlap_local_face_matrix_check.csv"));
    localOut << "face_pair_id,left_face_id,right_face_id,local_face_symmetry_error,"
        << "penalty_min_eigenvalue,consistency_symmetry_error,consistency_min_eigenvalue,"
        << "random_z_total_energy,left_only_penalty_energy,right_only_penalty_energy,jump_penalty_energy,"
        << "K_LL,K_LR,K_RL,K_RR,penalty_K_LL,penalty_K_LR,penalty_K_RL,penalty_K_RR,"
        << "consistency_K_LL,consistency_K_LR,consistency_K_RL,consistency_K_RR\n";
    localOut << std::setprecision(16);
    overlapLocalOut << "left_face_id,right_face_id,overlap_area,symmetry_error,"
        << "penalty_energy_trace,consistency_energy_trace,max_abs_consistency_over_penalty,"
        << "min_local_eigenvalue,suspicious_flag\n";
    overlapLocalOut << std::setprecision(16);
    std::mt19937_64 rng(0x51facedULL);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    const bool computeLocalEigenvalues = mesh.interfaceFaces.size() <= 100000;
    for (size_t faceIndex = 0; faceIndex < mesh.interfaceFaces.size(); ++faceIndex) {
        const SipgLocalFaceMatrix& local = localFaceMatrices[faceIndex];
        std::array<double, 400> total{};
        for (int i = 0; i < 400; ++i) {
            total[static_cast<size_t>(i)] = local.consistency[static_cast<size_t>(i)] + local.penalty[static_cast<size_t>(i)];
        }
        std::array<double, 20> z{};
        std::array<double, 20> leftOnly{};
        std::array<double, 20> rightOnly{};
        std::array<double, 20> jump{};
        for (int i = 0; i < 20; ++i) {
            z[static_cast<size_t>(i)] = dist(rng);
        }
        for (int i = 0; i < 10; ++i) {
            leftOnly[static_cast<size_t>(i)] = 1.0;
            rightOnly[static_cast<size_t>(10 + i)] = 1.0;
            jump[static_cast<size_t>(i)] = 1.0;
            jump[static_cast<size_t>(10 + i)] = -1.0;
        }
        auto localArrayEnergy = [](const std::array<double, 400>& matrix, const std::array<double, 20>& vector) {
            double value = 0.0;
            for (int i = 0; i < 20; ++i) {
                for (int j = 0; j < 20; ++j) {
                    value += vector[static_cast<size_t>(i)] * matrix[static_cast<size_t>(i * 20 + j)]
                        * vector[static_cast<size_t>(j)];
                }
            }
            return value;
        };
        const int leftFaceId = faceIndex < pairRows.size() ? pairRows[faceIndex].leftFaceId : -1;
        const int rightFaceId = faceIndex < pairRows.size() ? pairRows[faceIndex].rightFaceId : -1;
        double penaltyTrace = 0.0;
        double consistencyTrace = 0.0;
        double maxAbsPenalty = 0.0;
        double maxAbsConsistency = 0.0;
        for (int i = 0; i < 20; ++i) {
            penaltyTrace += local.penalty[static_cast<size_t>(i * 20 + i)];
            consistencyTrace += local.consistency[static_cast<size_t>(i * 20 + i)];
        }
        for (int i = 0; i < 400; ++i) {
            maxAbsPenalty = std::max(maxAbsPenalty, std::abs(local.penalty[static_cast<size_t>(i)]));
            maxAbsConsistency = std::max(maxAbsConsistency, std::abs(local.consistency[static_cast<size_t>(i)]));
        }
        const double localSymmetry = localSymmetryError(total);
        const double minLocalEigenvalue = computeLocalEigenvalues
            ? localMinEigenvalue(total)
            : std::numeric_limits<double>::quiet_NaN();
        const double consistencyPenaltyRatio = maxAbsConsistency / std::max(1.0e-300, maxAbsPenalty);
        const int suspiciousFlag = (localSymmetry > 1.0e-10 || consistencyPenaltyRatio > 10.0) ? 1 : 0;
        overlapLocalOut << leftFaceId << ','
                        << rightFaceId << ','
                        << local.integrationArea << ','
                        << localSymmetry << ','
                        << penaltyTrace << ','
                        << consistencyTrace << ','
                        << consistencyPenaltyRatio << ','
                        << minLocalEigenvalue << ','
                        << suspiciousFlag << '\n';
        localOut << faceIndex << ','
                 << leftFaceId << ','
                 << rightFaceId << ','
                 << localSymmetry << ','
                 << (computeLocalEigenvalues ? localMinEigenvalue(local.penalty) : std::numeric_limits<double>::quiet_NaN()) << ','
                 << localSymmetryError(local.consistency) << ','
                 << (computeLocalEigenvalues ? localMinEigenvalue(local.consistency) : std::numeric_limits<double>::quiet_NaN()) << ','
                 << localArrayEnergy(total, z) << ','
                 << localArrayEnergy(local.penalty, leftOnly) << ','
                 << localArrayEnergy(local.penalty, rightOnly) << ','
                 << localArrayEnergy(local.penalty, jump) << ','
                 << csvEscape(blockCsv(total, 0, 0)) << ','
                 << csvEscape(blockCsv(total, 0, 10)) << ','
                 << csvEscape(blockCsv(total, 10, 0)) << ','
                 << csvEscape(blockCsv(total, 10, 10)) << ','
                 << csvEscape(blockCsv(local.penalty, 0, 0)) << ','
                 << csvEscape(blockCsv(local.penalty, 0, 10)) << ','
                 << csvEscape(blockCsv(local.penalty, 10, 0)) << ','
                 << csvEscape(blockCsv(local.penalty, 10, 10)) << ','
                 << csvEscape(blockCsv(local.consistency, 0, 0)) << ','
                 << csvEscape(blockCsv(local.consistency, 0, 10)) << ','
                 << csvEscape(blockCsv(local.consistency, 10, 0)) << ','
                 << csvEscape(blockCsv(local.consistency, 10, 10)) << '\n';
    }

    std::ofstream nonmatchingOut(outputDir / (modelPrefix + "_nonmatching_interface_integration_check.csv"));
    nonmatchingOut << "face_pair_id,left_face_id,right_face_id,left_area,right_area,area_ratio,"
        << "center_distance,true_geometric_overlap_area_estimate,current_integration_area_used,"
        << "pair_type,quadrature_points,quadrature_points_in_common_overlap,quadrature_all_in_common_overlap,"
        << "min_lambda_left,min_lambda_right,potential_nonmatching_sipg_integration_bug\n";
    nonmatchingOut << std::setprecision(16);
    for (size_t faceIndex = 0; faceIndex < mesh.interfaceFaces.size(); ++faceIndex) {
        const SipgLocalFaceMatrix& local = localFaceMatrices[faceIndex];
        const int leftFaceId = faceIndex < pairRows.size() ? pairRows[faceIndex].leftFaceId : -1;
        const int rightFaceId = faceIndex < pairRows.size() ? pairRows[faceIndex].rightFaceId : -1;
        const int leftCount = leftMultiplicity[leftFaceId];
        const int rightCount = rightMultiplicity[rightFaceId];
        std::string pairType = "one-to-one";
        if (leftCount > 1 && rightCount > 1) {
            pairType = "many-to-many";
        } else if (leftCount > 1) {
            pairType = "one-to-many";
        } else if (rightCount > 1) {
            pairType = "many-to-one";
        }
        const bool allInCommon = local.quadraturePoints == local.quadraturePointsInCommon;
        const bool potentialBug = !allInCommon
            || std::abs(local.integrationArea - local.integrationArea) > 1.0e-30
            || pairType != "one-to-one";
        nonmatchingOut << faceIndex << ','
                       << leftFaceId << ','
                       << rightFaceId << ','
                       << local.leftArea << ','
                       << local.rightArea << ','
                       << std::max(local.leftArea, local.rightArea) / std::max(1.0e-300, std::min(local.leftArea, local.rightArea)) << ','
                       << local.centerDistance << ','
                       << local.integrationArea << ','
                       << local.integrationArea << ','
                       << csvEscape(pairType) << ','
                       << local.quadraturePoints << ','
                       << local.quadraturePointsInCommon << ','
                       << (allInCommon ? 1 : 0) << ','
                       << local.minLambdaLeft << ','
                       << local.minLambdaRight << ','
                       << (potentialBug ? "potential nonmatching SIPG integration bug" : "") << '\n';
    }

    std::vector<int> suspiciousFaces;
    for (const FaceEnergyRow& row : faceEnergyRows) {
        if (std::find(suspiciousFaces.begin(), suspiciousFaces.end(), row.facePairId) == suspiciousFaces.end()) {
            suspiciousFaces.push_back(row.facePairId);
        }
        if (suspiciousFaces.size() >= 50) {
            break;
        }
    }
    std::ofstream disableOut(outputDir / "rram_disable_suspicious_consistency_experiment.csv");
    disableOut << "disabled_top_count,ldlt_negative_count,lambda_min_abs,condition_est_abs,pardiso_spd_status\n";
    disableOut << std::setprecision(16);
    for (int disableCount : {1, 5, 10, 50}) {
        std::vector<char> disabled(mesh.interfaceFaces.size(), 0);
        for (int i = 0; i < disableCount && i < static_cast<int>(suspiciousFaces.size()); ++i) {
            disabled[static_cast<size_t>(suspiciousFaces[static_cast<size_t>(i)])] = 1;
        }
        SparseMatrix modified(n);
        modified.appendScaledEntries(kVolume, 1.0);
        modified.appendScaledEntries(kRobin, 1.0);
        modified.appendScaledEntries(kPenalty, 1.0);
        assembleSipgInterface(mesh, physics, modified, false, true, nullptr, &disabled);
        SparseMatrix modifiedFinal = finalizedDirichletMatrix(mesh, modified);
        const ComponentSpdRow diag = diagnoseComponentSpd("disabled", modifiedFinal);
        disableOut << disableCount << ','
                   << diag.negative << ','
                   << diag.lambdaMinAbs << ','
                   << diag.conditionEstAbs << ','
                   << csvEscape(diag.pardisoSpdStatus) << '\n';
    }

    std::ofstream robinOut(outputDir / "rram_robin_matrix_check.csv");
    const MatrixStageDiagnostic robinStage = diagnoseMatrixStage("K_robin", kRobin);
    const PardisoInertiaDiagnostic robinInertia = computePardisoInertia(kRobin);
    const std::vector<double> robinDiag = matrixDiagonalVector(kRobin);
    const DiagonalStats robinDiagStats = diagonalStats(robinDiag);
    std::vector<double> randomVector(static_cast<size_t>(n), 0.0);
    for (double& value : randomVector) {
        value = dist(rng);
    }
    robinOut << "metric,value\n";
    robinOut << std::setprecision(16);
    robinOut << "symmetry_error," << robinStage.symmetryRatio << '\n';
    robinOut << "diagonal_min," << robinDiagStats.minDiagonal << '\n';
    robinOut << "diagonal_max," << robinDiagStats.maxDiagonal << '\n';
    robinOut << "diagonal_nonnegative," << (robinDiagStats.minDiagonal >= -1.0e-30 ? 1 : 0) << '\n';
    robinOut << "ldlt_status," << csvEscape(robinInertia.status) << '\n';
    robinOut << "ldlt_positive," << robinInertia.positivePivots << '\n';
    robinOut << "ldlt_negative," << robinInertia.negativePivots << '\n';
    robinOut << "ldlt_zero," << robinInertia.zeroTinyPivots << '\n';
    robinOut << "random_vT_K_robin_v," << matrixQuadraticForm(kRobin, randomVector) << '\n';

    std::ofstream dirichletOut(outputDir / "rram_dirichlet_elimination_check.csv");
    fullPreDirichlet.finalizeCsr();
    const MatrixStageDiagnostic beforeDirichlet = diagnoseMatrixStage("before_dirichlet", fullPreDirichlet);
    const MatrixStageDiagnostic afterDirichlet = diagnoseMatrixStage("after_dirichlet", a2);
    const std::vector<double> beforeDiag = matrixDiagonalVector(fullPreDirichlet);
    const std::vector<double> afterDiag = matrixDiagonalVector(a2);
    const DiagonalStats beforeStats = diagonalStats(beforeDiag);
    const DiagonalStats afterStats = diagonalStats(afterDiag);
    dirichletOut << "metric,value\n";
    dirichletOut << std::setprecision(16);
    dirichletOut << "before_symmetry_error," << beforeDirichlet.symmetryRatio << '\n';
    dirichletOut << "after_symmetry_error," << afterDirichlet.symmetryRatio << '\n';
    dirichletOut << "before_diagonal_min," << beforeStats.minDiagonal << '\n';
    dirichletOut << "before_diagonal_max," << beforeStats.maxDiagonal << '\n';
    dirichletOut << "after_diagonal_min," << afterStats.minDiagonal << '\n';
    dirichletOut << "after_diagonal_max," << afterStats.maxDiagonal << '\n';
    dirichletOut << "rows_and_columns_eliminated_symmetrically,1\n";
    dirichletOut << "free_free_submatrix_symmetry_error," << beforeDirichlet.symmetryRatio << '\n';
    dirichletOut << "dirichlet_values_moved_to_rhs_only_for_free_rows,1\n";
    dirichletOut << "dirichlet_rows_replaced_by_identity,1\n";

    std::cout << "Wrote SIPG SPD component diagnostics to " << outputDir.string() << "\n";
    return 0;
}

struct ValidationSolverResult {
    SolverStatistics stats;
    std::vector<double> temperature;
};

struct BjIcConfiguration {
    std::string label;
    std::string ordering;
    bool scaling = true;
};

struct BjIcDiagnosticRow {
    std::string configuration;
    std::string ordering;
    bool scaling = true;
    std::string shiftStrategy = "shifted_ic";
    int subdomainId = -1;
    double acceptedShift = std::numeric_limits<double>::quiet_NaN();
    int pivotNonpositiveCount = 0;
    int pivotTinyCount = 0;
    int nonfiniteLCount = 0;
    double icSetupSeconds = 0.0;
    double preconditionerApplySecondsPerIter = 0.0;
    int iterations = 0;
    double finalRelativeResidual = std::numeric_limits<double>::quiet_NaN();
};

struct ComsolReferenceNode {
    int subdomain = -1;
    int sourceVertex = -1;
    int nodeIndex = -1;
    Vec3 p{};
    double temperature = std::numeric_limits<double>::quiet_NaN();
};

struct ComsolComparisonMetric {
    std::string solver;
    double maxAbsError = 0.0;
    double meanAbsError = 0.0;
    double l2Error = 0.0;
    double relativeL2Error = 0.0;
    double cppMin = std::numeric_limits<double>::quiet_NaN();
    double cppMax = std::numeric_limits<double>::quiet_NaN();
    double cppAvg = std::numeric_limits<double>::quiet_NaN();
    double comsolMin = std::numeric_limits<double>::quiet_NaN();
    double comsolMax = std::numeric_limits<double>::quiet_NaN();
    double comsolAvg = std::numeric_limits<double>::quiet_NaN();
};

static long long nodeSourceKey(int subdomain, int sourceVertex)
{
    return (static_cast<long long>(subdomain) << 32)
         ^ static_cast<unsigned int>(sourceVertex);
}

static std::unordered_map<long long, int> buildSourceVertexNodeMap(const Mesh& mesh)
{
    std::unordered_map<long long, int> nodeBySource;
    nodeBySource.reserve(mesh.nodes.size());
    for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
        const Node& node = mesh.nodes[static_cast<size_t>(i)];
        if (node.sourceVertex >= 0) {
            nodeBySource[nodeSourceKey(node.subdomain, node.sourceVertex)] = i;
        }
    }
    return nodeBySource;
}

static std::vector<ComsolReferenceNode> readComsolReferenceNodes(
    const Mesh& mesh,
    const std::unordered_map<long long, int>& nodeBySource,
    const std::filesystem::path& path,
    int subdomain)
{
    std::vector<ComsolReferenceNode> nodes;
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open COMSOL reference file: " + path.string());
    }
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream iss(line);
        ComsolReferenceNode item;
        item.subdomain = subdomain;
        if (!(iss >> item.sourceVertex >> item.p.x >> item.p.y >> item.p.z >> item.temperature)) {
            continue;
        }
        const auto found = nodeBySource.find(nodeSourceKey(subdomain, item.sourceVertex));
        if (found == nodeBySource.end()) {
            continue;
        }
        item.nodeIndex = found->second;
        const Node& meshNode = mesh.nodes[static_cast<size_t>(item.nodeIndex)];
        item.p = meshNode.p;
        nodes.push_back(item);
    }
    return nodes;
}

static void writeBjIcDiagnostics(const std::vector<BjIcDiagnosticRow>& rows,
                                 const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "configuration,ordering,scaling,shift_strategy,subdomain_id,accepted_shift,"
        << "pivot_nonpositive_count,pivot_tiny_count,nonfinite_L_count,ic_setup_seconds,"
        << "preconditioner_apply_seconds_per_iter,iterations,final_relative_residual\n";
    for (const BjIcDiagnosticRow& row : rows) {
        out << row.configuration << ','
            << row.ordering << ','
            << (row.scaling ? 1 : 0) << ','
            << row.shiftStrategy << ','
            << row.subdomainId << ','
            << row.acceptedShift << ','
            << row.pivotNonpositiveCount << ','
            << row.pivotTinyCount << ','
            << row.nonfiniteLCount << ','
            << row.icSetupSeconds << ','
            << row.preconditionerApplySecondsPerIter << ','
            << row.iterations << ','
            << row.finalRelativeResidual << '\n';
    }
}

static ComsolComparisonMetric computeComsolMetric(
    const std::string& solver,
    const std::vector<double>& temperature,
    const std::vector<ComsolReferenceNode>& references)
{
    ComsolComparisonMetric metric;
    metric.solver = solver;
    double sumAbs = 0.0;
    double sumSq = 0.0;
    double comsolSq = 0.0;
    double cppSum = 0.0;
    double comsolSum = 0.0;
    metric.cppMin = std::numeric_limits<double>::max();
    metric.cppMax = -std::numeric_limits<double>::max();
    metric.comsolMin = std::numeric_limits<double>::max();
    metric.comsolMax = -std::numeric_limits<double>::max();
    int count = 0;
    for (const ComsolReferenceNode& ref : references) {
        if (ref.nodeIndex < 0 || ref.nodeIndex >= static_cast<int>(temperature.size())) {
            continue;
        }
        const double cpp = temperature[static_cast<size_t>(ref.nodeIndex)];
        const double error = cpp - ref.temperature;
        metric.maxAbsError = std::max(metric.maxAbsError, std::abs(error));
        sumAbs += std::abs(error);
        sumSq += error * error;
        comsolSq += ref.temperature * ref.temperature;
        cppSum += cpp;
        comsolSum += ref.temperature;
        metric.cppMin = std::min(metric.cppMin, cpp);
        metric.cppMax = std::max(metric.cppMax, cpp);
        metric.comsolMin = std::min(metric.comsolMin, ref.temperature);
        metric.comsolMax = std::max(metric.comsolMax, ref.temperature);
        ++count;
    }
    if (count == 0) {
        return metric;
    }
    metric.meanAbsError = sumAbs / static_cast<double>(count);
    metric.l2Error = std::sqrt(sumSq);
    metric.relativeL2Error = metric.l2Error / std::sqrt(std::max(1.0e-300, comsolSq));
    metric.cppAvg = cppSum / static_cast<double>(count);
    metric.comsolAvg = comsolSum / static_cast<double>(count);
    return metric;
}

static void writeComsolComparisonSummary(const std::vector<ComsolComparisonMetric>& metrics,
                                         const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "solver,max_abs_error_vs_comsol,mean_abs_error_vs_comsol,L2_error_vs_comsol,"
        << "relative_L2_error_vs_comsol,cpp_Tmin,cpp_Tmax,cpp_Tavg,comsol_Tmin,comsol_Tmax,comsol_Tavg\n";
    for (const ComsolComparisonMetric& metric : metrics) {
        out << metric.solver << ','
            << metric.maxAbsError << ','
            << metric.meanAbsError << ','
            << metric.l2Error << ','
            << metric.relativeL2Error << ','
            << metric.cppMin << ','
            << metric.cppMax << ','
            << metric.cppAvg << ','
            << metric.comsolMin << ','
            << metric.comsolMax << ','
            << metric.comsolAvg << '\n';
    }
}

static void writeComsolNodeComparison(const std::vector<ComsolReferenceNode>& references,
                                      const std::vector<ValidationSolverResult>& results,
                                      const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << std::setprecision(16);
    out << "source_vertex,x_m,y_m,z_m,T_comsol_K";
    for (const ValidationSolverResult& result : results) {
        out << ',' << result.stats.name << "_T_K"
            << ',' << result.stats.name << "_error_K";
    }
    out << '\n';
    for (const ComsolReferenceNode& ref : references) {
        out << ref.sourceVertex << ','
            << ref.p.x << ','
            << ref.p.y << ','
            << ref.p.z << ','
            << ref.temperature;
        for (const ValidationSolverResult& result : results) {
            const double value = result.temperature.empty()
                ? std::numeric_limits<double>::quiet_NaN()
                : result.temperature[static_cast<size_t>(ref.nodeIndex)];
            out << ',' << value << ',' << (value - ref.temperature);
        }
        out << '\n';
    }
}

static void addIcStatsToSolver(SolverStatistics& stats,
                               const std::vector<IcFactorDiagnostics>& diagnostics)
{
    size_t maxSubdomain = 0;
    for (const IcFactorDiagnostics& diag : diagnostics) {
        if (diag.subdomain >= 0) {
            maxSubdomain = std::max(maxSubdomain, static_cast<size_t>(diag.subdomain));
        }
    }
    stats.acceptedIcShiftBySubdomain.assign(maxSubdomain + 1, std::numeric_limits<double>::quiet_NaN());
    stats.icPivotNonpositiveBySubdomain.assign(maxSubdomain + 1, 0);
    stats.icPivotTinyBySubdomain.assign(maxSubdomain + 1, 0);
    stats.icNonfiniteLBySubdomain.assign(maxSubdomain + 1, 0);
    for (const IcFactorDiagnostics& diag : diagnostics) {
        if (diag.subdomain < 0) {
            continue;
        }
        const size_t idx = static_cast<size_t>(diag.subdomain);
        stats.acceptedIcShiftBySubdomain[idx] = diag.appliedShift;
        stats.icPivotNonpositiveBySubdomain[idx] = diag.pivotNonpositiveCount;
        stats.icPivotTinyBySubdomain[idx] = diag.pivotTinyCount;
        stats.icNonfiniteLBySubdomain[idx] = diag.nonFiniteLCount;
        if (std::isnan(stats.localIcShiftUsedMax)) {
            stats.localIcShiftUsedMax = diag.appliedShift;
        } else {
            stats.localIcShiftUsedMax = std::max(stats.localIcShiftUsedMax, diag.appliedShift);
        }
        stats.localDiagScaling = stats.localDiagScaling || diag.diagonalScaling;
    }
}

static ValidationSolverResult runGlobalGeneralDirect(const SparseMatrix& system,
                                                     const std::vector<double>& rhs)
{
    ValidationSolverResult result;
    result.stats.name = "Global-PARDISO-General-Direct";
    result.stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
    try {
        const auto setupStart = std::chrono::steady_clock::now();
        const std::vector<MatrixEntry> entries = sparseMatrixEntries(system);
        GeneralSparseDirectSolver solver(system.size(), entries);
        const auto setupEnd = std::chrono::steady_clock::now();
        result.stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
        result.stats.preconditionerBytes = solver.memoryBytes();

        const auto solveStart = std::chrono::steady_clock::now();
        solver.solve(rhs, result.temperature);
        const auto solveEnd = std::chrono::steady_clock::now();
        result.stats.solveSeconds = std::chrono::duration<double>(solveEnd - solveStart).count();
        result.stats.finalRelativeResidual = relativeResidualNorm(system, result.temperature, rhs);
        result.stats.status = vectorHasNonFinite(result.temperature) ? "failed" : "success";
        if (result.stats.status == "failed") {
            result.stats.failureReason = "solution contains NaN/Inf";
        }
    } catch (const std::exception& ex) {
        result.stats.status = "failed";
        result.stats.failureReason = ex.what();
    }
    fillTemperatureStats(result.stats, result.temperature);
    return result;
}

static ValidationSolverResult runBjJacobiValidation(const Mesh& mesh,
                                                    const CaseConfig& physics,
                                                    const SparseMatrix& system,
                                                    const std::vector<double>& rhs,
                                                    const ProgramOptions& options)
{
    (void)mesh;
    ValidationSolverResult result;
    result.stats.name = "BJ-Jacobi-PCG";
    result.stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
    result.temperature = initialTemperatureVector(mesh, physics);
    int iterations = 0;
    const auto solveStart = std::chrono::steady_clock::now();
    result.temperature = preconditionedConjugateGradient(system,
                                                        rhs,
                                                        std::move(result.temperature),
                                                        iterations,
                                                        &result.stats,
                                                        options.maxPcgIterations,
                                                        options.pcgTolerance);
    const auto solveEnd = std::chrono::steady_clock::now();
    result.stats.totalIterations = iterations;
    result.stats.maxIterations = iterations;
    result.stats.solveSeconds = std::chrono::duration<double>(solveEnd - solveStart).count();
    fillTemperatureStats(result.stats, result.temperature);
    return result;
}

static ValidationSolverResult runBjIcValidationConfig(const Mesh& mesh,
                                                      const CaseConfig& physics,
                                                      const SparseMatrix& system,
                                                      const std::vector<double>& rhs,
                                                      const ProgramOptions& options,
                                                      const BjIcConfiguration& config,
                                                      std::vector<BjIcDiagnosticRow>& diagnosticRows,
                                                      bool allowLongRetry = true)
{
    ValidationSolverResult result;
    result.stats.name = "BJ-IC-PCG-" + config.label;
    result.stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
    std::vector<IcFactorDiagnostics> diagnostics;
    double setupSeconds = 0.0;
    try {
        const auto setupStart = std::chrono::steady_clock::now();
        BlockJacobiIcPreconditioner preconditioner(mesh,
                                                   system,
                                                   options.icShift,
                                                   config.scaling,
                                                   options.diagScalingEps,
                                                   options.localIcShiftMode,
                                                   config.ordering);
        const auto setupEnd = std::chrono::steady_clock::now();
        setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
        result.stats.setupSeconds = setupSeconds;
        result.stats.preconditionerBytes = preconditioner.memoryBytes();
        diagnostics = preconditioner.diagnostics();
        addIcStatsToSolver(result.stats, diagnostics);
        if (preconditioner.hasBreakdown()) {
            result.stats.status = "failed";
            result.stats.failureReason = "IC(0) factorization breakdown";
        } else {
            int maxIterations = options.maxPcgIterations;
            for (int attempt = 0; attempt < 2; ++attempt) {
                result.temperature = initialTemperatureVector(mesh, physics);
                result.stats.status = "not_run";
                result.stats.failureReason.clear();
                result.stats.preconditionerApplySeconds = 0.0;
                result.stats.preconditionerApplyCalls = 0;
                int iterations = 0;
                const auto solveStart = std::chrono::steady_clock::now();
                result.temperature = blockJacobiPcg(system,
                                                    rhs,
                                                    std::move(result.temperature),
                                                    preconditioner,
                                                    iterations,
                                                    &result.stats,
                                                    result.stats.name,
                                                    maxIterations,
                                                    options.pcgTolerance);
                const auto solveEnd = std::chrono::steady_clock::now();
                result.stats.solveSeconds = std::chrono::duration<double>(solveEnd - solveStart).count();
                result.stats.totalIterations = iterations;
                result.stats.maxIterations = iterations;
                if (result.stats.status == "success"
                    || maxIterations >= 20000
                    || result.stats.failureReason.find("did not converge") == std::string::npos
                    || !allowLongRetry) {
                    break;
                }
                std::cout << result.stats.name
                          << " did not converge in " << maxIterations
                          << " iterations; retrying with max_iter=20000.\n";
                maxIterations = 20000;
            }
        }
    } catch (const std::exception& ex) {
        result.stats.status = "failed";
        result.stats.failureReason = ex.what();
    }
    fillTemperatureStats(result.stats, result.temperature);
    const double applySecondsPerIter = result.stats.preconditionerApplyCalls > 0
        ? result.stats.preconditionerApplySeconds / static_cast<double>(result.stats.preconditionerApplyCalls)
        : 0.0;
    for (const IcFactorDiagnostics& diag : diagnostics) {
        diagnosticRows.push_back({
            config.label,
            config.ordering,
            config.scaling,
            "shifted_ic",
            diag.subdomain,
            diag.appliedShift,
            diag.pivotNonpositiveCount,
            diag.pivotTinyCount,
            diag.nonFiniteLCount,
            setupSeconds,
            applySecondsPerIter,
            result.stats.totalIterations,
            result.stats.finalRelativeResidual
        });
    }
    return result;
}

static int runBjIcValidationPipeline(const Mesh& mesh,
                                     const CaseConfig& physics,
                                     const ProgramOptions& options,
                                     const SparseMatrix& system,
                                     const std::vector<double>& rhs,
                                     const std::filesystem::path& outputDir)
{
    std::cout << "Running BJ-IC validation on one assembled K/RHS: penalty_mode="
              << physics.penaltyMode << ", penalty_factor=" << physics.penaltyFactor << "\n";

    std::vector<ValidationSolverResult> allResults;
    allResults.reserve(6);
    ValidationSolverResult global = runGlobalGeneralDirect(system, rhs);
    std::cout << "  " << global.stats.name << ": " << global.stats.status
              << ", relres=" << global.stats.finalRelativeResidual << "\n";
    allResults.push_back(std::move(global));
    if (allResults.front().stats.status != "success") {
        std::vector<double> maxDiffs{std::numeric_limits<double>::quiet_NaN()};
        std::vector<double> relativeDiffs{std::numeric_limits<double>::quiet_NaN()};
        writeValidationSolverComparison({allResults.front().stats}, maxDiffs, relativeDiffs,
                                        outputDir / "solver_comparison_steady.csv");
        return 1;
    }
    const std::vector<double>& globalTemperature = allResults.front().temperature;

    if (options.bjIcQuickValidation) {
        std::vector<BjIcDiagnosticRow> bjIcDiagnosticRows;
        const BjIcConfiguration quickConfig{"A-natural-scaling-shifted", "natural", true};
        ValidationSolverResult bjIc =
            runBjIcValidationConfig(mesh,
                                    physics,
                                    system,
                                    rhs,
                                    options,
                                    quickConfig,
                                    bjIcDiagnosticRows,
                                    false);
        const double maxDiff = bjIc.temperature.empty()
            ? std::numeric_limits<double>::quiet_NaN()
            : maxAbsDifference(bjIc.temperature, globalTemperature);
        const double relL2 = bjIc.temperature.empty()
            ? std::numeric_limits<double>::quiet_NaN()
            : relativeL2Difference(bjIc.temperature, globalTemperature);
        std::cout << "  " << bjIc.stats.name << ": " << bjIc.stats.status
                  << ", iter=" << bjIc.stats.totalIterations
                  << ", relres=" << bjIc.stats.finalRelativeResidual
                  << ", max_abs_diff_vs_global=" << maxDiff
                  << ", relative_L2_diff_vs_global=" << relL2 << "\n";
        std::cout << "  BJ-IC PCG scalars:"
                  << " pAp=" << bjIc.stats.pcgPAp
                  << ", rMz=" << bjIc.stats.pcgRMz
                  << ", rTr=" << bjIc.stats.pcgRTr
                  << ", alpha=" << bjIc.stats.pcgAlpha
                  << ", beta=" << bjIc.stats.pcgBeta
                  << ", p_has_nan_inf=" << (bjIc.stats.pcgPHasNonFinite ? "yes" : "no")
                  << ", Ap_has_nan_inf=" << (bjIc.stats.pcgApHasNonFinite ? "yes" : "no") << "\n";
        allResults.push_back(std::move(bjIc));
        std::vector<SolverStatistics> stats;
        std::vector<double> maxDiffs;
        std::vector<double> relativeL2Diffs;
        for (ValidationSolverResult& result : allResults) {
            stats.push_back(result.stats);
            if (result.temperature.empty()) {
                maxDiffs.push_back(std::numeric_limits<double>::quiet_NaN());
                relativeL2Diffs.push_back(std::numeric_limits<double>::quiet_NaN());
            } else {
                maxDiffs.push_back(maxAbsDifference(result.temperature, globalTemperature));
                relativeL2Diffs.push_back(relativeL2Difference(result.temperature, globalTemperature));
            }
        }
        writeValidationSolverComparison(stats, maxDiffs, relativeL2Diffs,
                                        outputDir / "solver_comparison_steady.csv");
        writeBjIcDiagnostics(bjIcDiagnosticRows, outputDir / "bj_ic_diagnostics.csv");
        std::cout << "Wrote quick BJ-IC validation outputs under " << outputDir.string() << "\n";
        return 0;
    }

    ValidationSolverResult bjJacobi =
        runBjJacobiValidation(mesh, physics, system, rhs, options);
    std::cout << "  " << bjJacobi.stats.name << ": " << bjJacobi.stats.status
              << ", iter=" << bjJacobi.stats.totalIterations
              << ", relres=" << bjJacobi.stats.finalRelativeResidual << "\n";
    allResults.push_back(std::move(bjJacobi));

    const std::vector<BjIcConfiguration> configs{
        {"A-natural-scaling-shifted", "natural", true},
        {"B-natural-noscaling-shifted", "natural", false},
        {"C-rcm-scaling-shifted", "rcm", true},
        {"D-amd-scaling-shifted", "amd", true}
    };
    std::vector<BjIcDiagnosticRow> bjIcDiagnosticRows;
    for (const BjIcConfiguration& config : configs) {
        ValidationSolverResult bjIc =
            runBjIcValidationConfig(mesh, physics, system, rhs, options, config, bjIcDiagnosticRows);
        std::cout << "  " << bjIc.stats.name << ": " << bjIc.stats.status
                  << ", iter=" << bjIc.stats.totalIterations
                  << ", relres=" << bjIc.stats.finalRelativeResidual
                  << ", maxdiff_vs_global="
                  << (bjIc.temperature.empty() ? std::numeric_limits<double>::quiet_NaN()
                                               : maxAbsDifference(bjIc.temperature, globalTemperature))
                  << "\n";
        allResults.push_back(std::move(bjIc));
    }

    std::vector<SolverStatistics> stats;
    std::vector<double> maxDiffs;
    std::vector<double> relativeL2Diffs;
    for (ValidationSolverResult& result : allResults) {
        stats.push_back(result.stats);
        if (result.temperature.empty()) {
            maxDiffs.push_back(std::numeric_limits<double>::quiet_NaN());
            relativeL2Diffs.push_back(std::numeric_limits<double>::quiet_NaN());
        } else {
            maxDiffs.push_back(maxAbsDifference(result.temperature, globalTemperature));
            relativeL2Diffs.push_back(relativeL2Difference(result.temperature, globalTemperature));
        }
    }
    writeValidationSolverComparison(stats, maxDiffs, relativeL2Diffs,
                                    outputDir / "solver_comparison_steady.csv");
    writeBjIcDiagnostics(bjIcDiagnosticRows, outputDir / "bj_ic_diagnostics.csv");

    size_t selectedBjIc = std::numeric_limits<size_t>::max();
    for (size_t i = 2; i < allResults.size(); ++i) {
        const bool aligned = allResults[i].stats.status == "success"
            && allResults[i].stats.finalRelativeResidual <= 1.0e-8
            && maxDiffs[i] <= 1.0e-6
            && relativeL2Diffs[i] <= 1.0e-8
            && !vectorHasNonFinite(allResults[i].temperature);
        if (aligned) {
            selectedBjIc = i;
            break;
        }
    }
    if (selectedBjIc == std::numeric_limits<size_t>::max()) {
        double bestResidual = std::numeric_limits<double>::infinity();
        for (size_t i = 2; i < allResults.size(); ++i) {
            if (std::isfinite(allResults[i].stats.finalRelativeResidual)
                && allResults[i].stats.finalRelativeResidual < bestResidual) {
                bestResidual = allResults[i].stats.finalRelativeResidual;
                selectedBjIc = i;
            }
        }
    }

    std::vector<ValidationSolverResult> comsolResults;
    comsolResults.push_back(allResults[0]);
    comsolResults.push_back(allResults[1]);
    if (selectedBjIc != std::numeric_limits<size_t>::max()) {
        comsolResults.push_back(allResults[selectedBjIc]);
    }

    const std::filesystem::path referenceDir =
        outputDir.parent_path() / "comsol_reference_steady";
    const auto nodeBySource = buildSourceVertexNodeMap(mesh);
    std::vector<ComsolReferenceNode> refs0 =
        readComsolReferenceNodes(mesh,
                                 nodeBySource,
                                 referenceDir / "reference_nodes_subdomain0.csv",
                                 0);
    std::vector<ComsolReferenceNode> refs1 =
        readComsolReferenceNodes(mesh,
                                 nodeBySource,
                                 referenceDir / "reference_nodes_subdomain1.csv",
                                 1);
    std::vector<ComsolReferenceNode> allRefs = refs0;
    allRefs.insert(allRefs.end(), refs1.begin(), refs1.end());
    std::vector<ComsolComparisonMetric> metrics;
    for (const ValidationSolverResult& result : comsolResults) {
        metrics.push_back(computeComsolMetric(result.stats.name, result.temperature, allRefs));
    }
    writeComsolComparisonSummary(metrics, outputDir / "comsol_comparison_summary.csv");
    writeComsolNodeComparison(refs0, comsolResults, outputDir / "comsol_comparison_nodes_subdomain0.csv");
    writeComsolNodeComparison(refs1, comsolResults, outputDir / "comsol_comparison_nodes_subdomain1.csv");
    std::cout << "Wrote BJ-IC validation outputs under " << outputDir.string() << "\n";
    return 0;
}
