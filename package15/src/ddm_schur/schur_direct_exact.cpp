#include "schur_direct_exact.hpp"

#include "interface_operator.hpp"
#include "local_solver.hpp"
#include "../sipg_core.hpp"
#include "../linear_solvers.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_map>

#ifdef _WIN32
#include <windows.h>
#endif

namespace ddm_schur {
namespace {

using Clock = std::chrono::steady_clock;

std::size_t physicalMemoryBytes()
{
#ifdef _WIN32
    MEMORYSTATUSEX status{};
    status.dwLength = sizeof(status);
    if (GlobalMemoryStatusEx(&status) != 0) {
        return static_cast<std::size_t>(status.ullTotalPhys);
    }
#endif
    return std::size_t(0);
}

std::uint64_t checkedAdd(std::uint64_t a, std::uint64_t b, const char* what)
{
    if (b > std::numeric_limits<std::uint64_t>::max() - a) {
        throw std::runtime_error(std::string("Exact Schur ") + what + " overflow.");
    }
    return a + b;
}

std::uint64_t checkedClique(std::size_t n)
{
    const std::uint64_t u = static_cast<std::uint64_t>(n);
    if (u != 0 && u > (std::numeric_limits<std::uint64_t>::max() - u) / 2) {
        throw std::runtime_error("Exact Schur clique nnz overflow.");
    }
    return u * (u + 1) / 2;
}

double relativeDifference(const std::vector<double>& a, const std::vector<double>& b,
                          double* maxAbsolute = nullptr)
{
    if (a.size() != b.size()) {
        return std::numeric_limits<double>::infinity();
    }
    double norm = 0.0;
    double diff = 0.0;
    double maxAbs = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        norm += a[i] * a[i];
        const double d = a[i] - b[i];
        diff += d * d;
        maxAbs = std::max(maxAbs, std::abs(d));
    }
    if (maxAbsolute != nullptr) {
        *maxAbsolute = maxAbs;
    }
    return std::sqrt(diff) / std::max(std::sqrt(norm), 1.0e-300);
}

} // namespace

struct ExactSchurDirectSolver::Impl {
    const Mesh& mesh;
    const SparseMatrix& system;
    ExactSchurDirectOptions options;
    InterfacePartition partition;
    ExactSchurDirectReport report;
    std::vector<std::vector<int>> patternRows;
    std::vector<int> rowPtr;
    std::vector<int> colInd;
    std::vector<double> values;
    std::vector<LocalSolver> localSolvers;
    SubdomainDirectSolver schurSolver;
    bool runnable = false;
    std::size_t localFactorBytes = 0;

    Impl(const Mesh& meshIn, const SparseMatrix& systemIn,
         const CaseConfig&, const ExactSchurDirectOptions& optionsIn)
        : mesh(meshIn), system(systemIn), options(optionsIn)
    {
        const auto totalStart = Clock::now();
        report.indexBits = static_cast<int>(sizeof(int) * 8);
#ifdef USE_MKL_PARDISO
        report.indexBits = static_cast<int>(sizeof(MKL_INT) * 8);
#endif
        try {
            partition = buildInterfacePartition(mesh, system);
            report.interfaceDofs = static_cast<int>(partition.interfaceGlobalDofs.size());
            auditDimensions();
            if (!runnable) {
                report.totalSeconds = std::chrono::duration<double>(Clock::now() - totalStart).count();
                return;
            }
            if (options.dryRun) {
                report.status = "dry_run";
                report.totalSeconds = std::chrono::duration<double>(Clock::now() - totalStart).count();
                return;
            }
            buildPattern();
            if (!runnable || options.patternOnly || options.dryRun) {
                report.totalSeconds = std::chrono::duration<double>(Clock::now() - totalStart).count();
                return;
            }
            assembleLocalFactors();
            checkInterfaceCouplingSymmetry();
            assembleNumericSchur();
            if (options.verifyOperator) {
                verifyOperator();
            }
            factorSchur();
            report.status = "ready";
            runnable = true;
        } catch (const std::exception& error) {
            report.status = "failed";
            report.abortReason = error.what();
            runnable = false;
        }
        report.peakWorkingSetBytes = std::max(report.peakWorkingSetBytes, peakWorkingSetBytes());
        report.totalSeconds = std::chrono::duration<double>(Clock::now() - totalStart).count();
    }

    void stop(const std::string& reason, const std::string& status = "aborted")
    {
        report.status = status;
        report.abortReason = reason;
        runnable = false;
    }

    std::size_t rawBytesFor(std::uint64_t nnz) const
    {
        const std::uint64_t rows = static_cast<std::uint64_t>(report.interfaceDofs) + 1;
        const std::uint64_t bytes = rows * sizeof(int)
            + nnz * (sizeof(int) + sizeof(double));
        if (bytes > std::numeric_limits<std::size_t>::max()) {
            throw std::runtime_error("Exact Schur raw CSR byte estimate overflow.");
        }
        return static_cast<std::size_t>(bytes);
    }

    void auditDimensions()
    {
        if (report.interfaceDofs <= 0) {
            stop("no interface DOFs in partition");
            return;
        }
        const std::uint64_t maxIndex =
#ifdef USE_MKL_PARDISO
            static_cast<std::uint64_t>(std::numeric_limits<MKL_INT>::max());
#else
            static_cast<std::uint64_t>(std::numeric_limits<int>::max());
#endif
        if (static_cast<std::uint64_t>(report.interfaceDofs) > maxIndex) {
            stop("interface dimension exceeds MKL_INT range", "infeasible");
            return;
        }
        std::uint64_t predicted = 0;
        std::unordered_map<std::uint64_t, char> interfacePattern;
        interfacePattern.reserve(partition.interfaceEntries.size());
        for (const Entry& e : partition.interfaceEntries) {
            const int r = std::min(e.row, e.col);
            const int c = std::max(e.row, e.col);
            interfacePattern[(static_cast<std::uint64_t>(r) << 32)
                             | static_cast<std::uint32_t>(c)] = 1;
        }
        predicted = static_cast<std::uint64_t>(interfacePattern.size());
        for (const DomainBlocks& domain : partition.domains) {
            ExactSchurDirectDomainReport row;
            row.domainId = domain.domainId;
            row.interiorDofs = domain.interiorGlobalDofs.size();
            row.interfaceDofs = domain.interfaceGlobalDofs.size();
            row.cliqueUpperNnz = checkedClique(row.interfaceDofs);
            row.aiGammaNnz = 0;
            for (const auto& entries : domain.interiorInterfaceRows) {
                row.aiGammaNnz += static_cast<std::uint64_t>(entries.size());
            }
            row.estimatedDenseLocalSchurBytes = row.cliqueUpperNnz * sizeof(double);
            report.domains.push_back(row);
            predicted = checkedAdd(predicted, row.cliqueUpperNnz, "predicted upper nnz");
            std::cout << "[ExactSchur] domain " << row.domainId
                      << " nI=" << row.interiorDofs
                      << " nGamma=" << row.interfaceDofs
                      << " clique_upper_nnz=" << row.cliqueUpperNnz
                      << " AIgamma_nnz=" << row.aiGammaNnz
                      << " estimated_dense_local_schur_bytes="
                      << row.estimatedDenseLocalSchurBytes << "\n";
        }
        report.predictedUpperNnz = predicted;
        const std::size_t rawUpperBytes = rawBytesFor(predicted);
        const std::size_t physical = physicalMemoryBytes();
        std::cout << "[ExactSchur] interface_dofs=" << report.interfaceDofs
                  << " predicted_upper_nnz=" << report.predictedUpperNnz
                  << " raw_csr_upper_estimate_bytes=" << rawUpperBytes
                  << " index_bits=" << report.indexBits << "\n";
        if (physical != 0 && static_cast<double>(rawUpperBytes)
                > options.memoryFraction * static_cast<double>(physical)) {
            stop("raw exact Schur CSR estimate exceeds memory fraction gate", "infeasible");
            return;
        }
        if (report.predictedUpperNnz > maxIndex) {
            stop("predicted Schur nnz exceeds MKL_INT range", "infeasible");
            return;
        }
        runnable = true;
    }

    void buildPattern()
    {
        const auto start = Clock::now();
        patternRows.assign(static_cast<std::size_t>(report.interfaceDofs), {});
        for (const Entry& e : partition.interfaceEntries) {
            const int r = std::min(e.row, e.col);
            const int c = std::max(e.row, e.col);
            patternRows[static_cast<std::size_t>(r)].push_back(c);
        }
        for (const DomainBlocks& domain : partition.domains) {
            std::vector<int> gamma;
            gamma.reserve(domain.interfaceGlobalDofs.size());
            for (int dof : domain.interfaceGlobalDofs) {
                gamma.push_back(partition.globalToInterface[static_cast<std::size_t>(dof)]);
            }
            for (std::size_t i = 0; i < gamma.size(); ++i) {
                for (std::size_t j = i; j < gamma.size(); ++j) {
                    const int a = gamma[i];
                    const int b = gamma[j];
                    patternRows[static_cast<std::size_t>(std::min(a, b))]
                        .push_back(std::max(a, b));
                }
            }
        }
        rowPtr.assign(static_cast<std::size_t>(report.interfaceDofs + 1), 0);
        std::uint64_t uniqueNnz = 0;
        for (std::size_t row = 0; row < patternRows.size(); ++row) {
            auto& cols = patternRows[row];
            std::sort(cols.begin(), cols.end());
            cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
            uniqueNnz = checkedAdd(uniqueNnz, cols.size(), "actual upper nnz");
            rowPtr[row + 1] = static_cast<int>(uniqueNnz);
        }
        report.actualUniqueUpperNnz = uniqueNnz;
        report.patternDensity = static_cast<double>(uniqueNnz)
            / std::max(1.0, static_cast<double>(report.interfaceDofs)
                * (static_cast<double>(report.interfaceDofs) + 1.0) / 2.0);
        colInd.reserve(static_cast<std::size_t>(uniqueNnz));
        for (const auto& cols : patternRows) {
            colInd.insert(colInd.end(), cols.begin(), cols.end());
        }
        report.rawCsrBytes = rowPtr.size() * sizeof(int)
            + colInd.size() * (sizeof(int) + sizeof(double));
        report.patternSeconds = std::chrono::duration<double>(Clock::now() - start).count();
        if (physicalMemoryBytes() != 0 && static_cast<double>(report.rawCsrBytes)
                > options.memoryFraction * static_cast<double>(physicalMemoryBytes())) {
            stop("constructed exact Schur CSR exceeds memory fraction gate", "infeasible");
        }
        std::cout << "[ExactSchur] pattern actual_unique_upper_nnz="
                  << report.actualUniqueUpperNnz << " density=" << report.patternDensity
                  << " seconds=" << report.patternSeconds << "\n";
    }

    std::size_t patternPosition(int row, int col) const
    {
        const auto& cols = patternRows[static_cast<std::size_t>(row)];
        auto it = std::lower_bound(cols.begin(), cols.end(), col);
        if (it == cols.end() || *it != col) {
            throw std::runtime_error("Exact Schur numeric entry is absent from exact pattern.");
        }
        const std::size_t offset = static_cast<std::size_t>(rowPtr[static_cast<std::size_t>(row)]);
        return offset + static_cast<std::size_t>(it - cols.begin());
    }

    void initializeInterfaceBlock()
    {
        values.assign(colInd.size(), 0.0);
        for (const Entry& e : partition.interfaceEntries) {
            if (e.row > e.col) {
                continue;
            }
            values[patternPosition(e.row, e.col)] += e.value;
        }
    }

    void assembleLocalFactors()
    {
        localSolvers.clear();
        localSolvers.reserve(partition.domains.size());
        for (const DomainBlocks& domain : partition.domains) {
            localSolvers.emplace_back(static_cast<int>(domain.interiorGlobalDofs.size()),
                                      domain.interiorEntries, 1);
            localFactorBytes += localSolvers.back().memoryBytes();
        }
    }

    void checkInterfaceCouplingSymmetry()
    {
        double squaredDifference = 0.0;
        double squaredNorm = 0.0;
        for (const DomainBlocks& domain : partition.domains) {
            std::unordered_map<std::uint64_t, double> lower;
            std::unordered_map<std::uint64_t, double> upper;
            for (std::size_t row = 0; row < domain.interiorInterfaceRows.size(); ++row) {
                for (const auto& entry : domain.interiorInterfaceRows[row]) {
                    const std::uint64_t key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(row)) << 32)
                        | static_cast<std::uint32_t>(entry.first);
                    lower[key] += entry.second;
                    squaredNorm += entry.second * entry.second;
                }
            }
            for (std::size_t gamma = 0; gamma < domain.interfaceInteriorRows.size(); ++gamma) {
                const int globalGamma = partition.globalToInterface[
                    static_cast<std::size_t>(domain.interfaceGlobalDofs[gamma])];
                for (const auto& entry : domain.interfaceInteriorRows[gamma]) {
                    const std::uint64_t key = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(entry.first)) << 32)
                        | static_cast<std::uint32_t>(globalGamma);
                    upper[key] += entry.second;
                }
            }
            for (const auto& item : lower) {
                const double other = upper[item.first];
                const double difference = item.second - other;
                squaredDifference += difference * difference;
            }
            for (const auto& item : upper) {
                if (lower.find(item.first) == lower.end()) {
                    squaredDifference += item.second * item.second;
                    squaredNorm += item.second * item.second;
                }
            }
        }
        report.interfaceCouplingSymmetryError = std::sqrt(squaredDifference)
            / std::max(std::sqrt(squaredNorm), 1.0e-300);
        if (report.interfaceCouplingSymmetryError > 1.0e-12) {
            throw std::runtime_error("A_GammaI and A_IGamma are not numerically symmetric");
        }
    }

    void assembleNumericSchur()
    {
        const auto start = Clock::now();
        initializeInterfaceBlock();
        int maxGamma = 0;
        for (const DomainBlocks& domain : partition.domains) {
            maxGamma = std::max(maxGamma, static_cast<int>(domain.interfaceGlobalDofs.size()));
        }
        report.batchSize = options.batchSize > 0
            ? options.batchSize
            : (maxGamma > 64 ? 64 : (maxGamma > 32 ? 32 : (maxGamma > 16 ? 16 : 8)));
        report.batchSize = std::max(1, report.batchSize);
        const std::size_t positionBase = 0;
        (void)positionBase;
        for (std::size_t domainIndex = 0; domainIndex < partition.domains.size(); ++domainIndex) {
            const DomainBlocks& domain = partition.domains[domainIndex];
            const int nI = static_cast<int>(domain.interiorGlobalDofs.size());
            const int nGamma = static_cast<int>(domain.interfaceGlobalDofs.size());
            if (nI == 0 || nGamma == 0) {
                continue;
            }
            const auto localStart = Clock::now();
            std::vector<int> gammaGlobal;
            gammaGlobal.reserve(domain.interfaceGlobalDofs.size());
            for (int dof : domain.interfaceGlobalDofs) {
                gammaGlobal.push_back(partition.globalToInterface[static_cast<std::size_t>(dof)]);
            }
            for (int first = 0; first < nGamma; first += report.batchSize) {
                const int count = std::min(report.batchSize, nGamma - first);
                std::vector<double> rhs(static_cast<std::size_t>(nI) * static_cast<std::size_t>(count), 0.0);
                for (int row = 0; row < nI; ++row) {
                    for (const auto& entry : domain.interiorInterfaceRows[static_cast<std::size_t>(row)]) {
                        const auto it = std::find(gammaGlobal.begin() + first,
                                                  gammaGlobal.begin() + first + count,
                                                  entry.first);
                        if (it != gammaGlobal.begin() + first + count) {
                            const int column = static_cast<int>(it - (gammaGlobal.begin() + first));
                            rhs[static_cast<std::size_t>(column) * static_cast<std::size_t>(nI)
                                + static_cast<std::size_t>(row)] += entry.second;
                        }
                    }
                }
                std::vector<double> x;
                const auto solveStart = Clock::now();
                localSolvers[domainIndex].solveMultiple(rhs, count, x);
                const double solveSeconds = std::chrono::duration<double>(Clock::now() - solveStart).count();
                ++report.multiRhsCalls;
                report.localSolveSeconds += solveSeconds;
                report.domains[domainIndex].phase33Calls += 1;
                report.domains[domainIndex].phase33Seconds += solveSeconds;
                for (int localRow = 0; localRow < nGamma; ++localRow) {
                    const int globalRow = gammaGlobal[static_cast<std::size_t>(localRow)];
                    for (int column = 0; column < count; ++column) {
                        double value = 0.0;
                        for (const auto& entry : domain.interfaceInteriorRows[static_cast<std::size_t>(localRow)]) {
                            value += entry.second * x[static_cast<std::size_t>(column) * static_cast<std::size_t>(nI)
                                + static_cast<std::size_t>(entry.first)];
                        }
                        const int globalColumn = gammaGlobal[static_cast<std::size_t>(first + column)];
                        if (globalRow > globalColumn) {
                            continue;
                        }
                        const int row = std::min(globalRow, globalColumn);
                        const int col = std::max(globalRow, globalColumn);
                        values[patternPosition(row, col)] -= value;
                    }
                }
            }
            const double localTotal = std::chrono::duration<double>(Clock::now() - localStart).count();
            (void)localTotal;
        }
        report.assemblySeconds = std::chrono::duration<double>(Clock::now() - start).count();
        for (std::size_t i = 0; i < report.domains.size(); ++i) {
            report.domains[i].phase11Seconds = localSolvers[i].symbolicAnalysisSeconds();
            report.domains[i].phase22Seconds = localSolvers[i].numericalFactorizationSeconds();
            report.domains[i].factorMemoryBytes = localSolvers[i].memoryBytes();
            report.localPhase11Seconds += report.domains[i].phase11Seconds;
            report.localPhase22Seconds += report.domains[i].phase22Seconds;
        }
        report.localFactorBytes = localFactorBytes;
    }

    void factorSchur()
    {
        std::vector<MatrixEntry> entries;
        entries.reserve(values.size());
        for (std::size_t row = 0; row < patternRows.size(); ++row) {
            const int begin = rowPtr[row];
            const int end = rowPtr[row + 1];
            for (int p = begin; p < end; ++p) {
                entries.push_back({static_cast<int>(row), colInd[static_cast<std::size_t>(p)],
                                   values[static_cast<std::size_t>(p)]});
            }
        }
        schurSolver = SubdomainDirectSolver(report.interfaceDofs, entries);
        report.factorPhase11Seconds = schurSolver.symbolicAnalysisSeconds();
        report.factorPhase22Seconds = schurSolver.numericalFactorizationSeconds();
        report.factorBytes = schurSolver.memoryBytes();
        if (physicalMemoryBytes() != 0 && static_cast<double>(report.factorBytes)
                > options.factorMemoryFraction * static_cast<double>(physicalMemoryBytes())) {
            stop("exact Schur factor memory exceeds factor fraction gate", "infeasible");
            return;
        }
        report.peakWorkingSetBytes = std::max(report.peakWorkingSetBytes, peakWorkingSetBytes());
    }

    void matrixFreeApply(const std::vector<double>& x, std::vector<double>& y)
    {
        y.assign(static_cast<std::size_t>(report.interfaceDofs), 0.0);
        for (const Entry& e : partition.interfaceEntries) {
            y[static_cast<std::size_t>(e.row)] += e.value * x[static_cast<std::size_t>(e.col)];
        }
        for (std::size_t d = 0; d < partition.domains.size(); ++d) {
            const DomainBlocks& domain = partition.domains[d];
            const int nI = static_cast<int>(domain.interiorGlobalDofs.size());
            if (nI == 0) {
                continue;
            }
            std::vector<double> rhs(static_cast<std::size_t>(nI), 0.0);
            for (int row = 0; row < nI; ++row) {
                for (const auto& entry : domain.interiorInterfaceRows[static_cast<std::size_t>(row)]) {
                    rhs[static_cast<std::size_t>(row)] += entry.second * x[static_cast<std::size_t>(entry.first)];
                }
            }
            std::vector<double> z;
            localSolvers[d].solve(rhs, z);
            for (std::size_t row = 0; row < domain.interfaceInteriorRows.size(); ++row) {
                for (const auto& entry : domain.interfaceInteriorRows[row]) {
                    const int globalRow = partition.globalToInterface[static_cast<std::size_t>(domain.interfaceGlobalDofs[row])];
                    y[static_cast<std::size_t>(globalRow)] -= entry.second * z[static_cast<std::size_t>(entry.first)];
                }
            }
        }
    }

    void explicitApply(const std::vector<double>& x, std::vector<double>& y) const
    {
        y.assign(static_cast<std::size_t>(report.interfaceDofs), 0.0);
        for (std::size_t row = 0; row < patternRows.size(); ++row) {
            for (int p = rowPtr[row]; p < rowPtr[row + 1]; ++p) {
                const int col = colInd[static_cast<std::size_t>(p)];
                const double value = values[static_cast<std::size_t>(p)];
                y[row] += value * x[static_cast<std::size_t>(col)];
                if (col != static_cast<int>(row)) {
                    y[static_cast<std::size_t>(col)] += value * x[row];
                }
            }
        }
    }

    void verifyOperator()
    {
        std::mt19937_64 rng(0x4e584143ULL);
        std::uniform_real_distribution<double> distribution(-1.0, 1.0);
        double maxRel = 0.0;
        double maxAbs = 0.0;
        const int checks = std::max(1, options.randomChecks);
        for (int check = 0; check < checks; ++check) {
            std::vector<double> x(static_cast<std::size_t>(report.interfaceDofs), 0.0);
            for (double& value : x) {
                value = distribution(rng);
            }
            std::vector<double> explicitY;
            std::vector<double> matrixFreeY;
            explicitApply(x, explicitY);
            matrixFreeApply(x, matrixFreeY);
            double abs = 0.0;
            const double rel = relativeDifference(explicitY, matrixFreeY, &abs);
            maxRel = std::max(maxRel, rel);
            maxAbs = std::max(maxAbs, abs);
        }
        report.operatorRelativeError = maxRel;
        report.operatorMaxAbsoluteError = maxAbs;
        report.schurSymmetryError = maxRel;
        if (maxRel > 1.0e-10) {
            throw std::runtime_error("exact Schur explicit/matrix-free operator check failed");
        }
        for (int check = 0; check < checks; ++check) {
            std::vector<double> x(static_cast<std::size_t>(report.interfaceDofs), 0.0);
            for (double& value : x) {
                value = distribution(rng);
            }
            std::vector<double> y;
            explicitApply(x, y);
            double quadratic = 0.0;
            for (std::size_t i = 0; i < x.size(); ++i) {
                quadratic += x[i] * y[i];
            }
            if (!std::isfinite(quadratic) || quadratic <= 0.0) {
                throw std::runtime_error("exact Schur SPD random check failed");
            }
        }
    }

    void condensedRhs(const std::vector<double>& rhs, std::vector<double>& g)
    {
        g.assign(static_cast<std::size_t>(report.interfaceDofs), 0.0);
        for (std::size_t i = 0; i < partition.interfaceGlobalDofs.size(); ++i) {
            g[i] = rhs[static_cast<std::size_t>(partition.interfaceGlobalDofs[i])];
        }
        const auto start = Clock::now();
        for (std::size_t d = 0; d < partition.domains.size(); ++d) {
            const DomainBlocks& domain = partition.domains[d];
            const int nI = static_cast<int>(domain.interiorGlobalDofs.size());
            if (nI == 0) {
                continue;
            }
            std::vector<double> localRhs(static_cast<std::size_t>(nI), 0.0);
            for (int i = 0; i < nI; ++i) {
                localRhs[static_cast<std::size_t>(i)] = rhs[static_cast<std::size_t>(domain.interiorGlobalDofs[static_cast<std::size_t>(i)])];
            }
            std::vector<double> z;
            localSolvers[d].solve(localRhs, z);
            for (std::size_t row = 0; row < domain.interfaceInteriorRows.size(); ++row) {
                const int globalRow = partition.globalToInterface[static_cast<std::size_t>(domain.interfaceGlobalDofs[row])];
                for (const auto& entry : domain.interfaceInteriorRows[row]) {
                    g[static_cast<std::size_t>(globalRow)] -= entry.second * z[static_cast<std::size_t>(entry.first)];
                }
            }
        }
        report.rhsSeconds = std::chrono::duration<double>(Clock::now() - start).count();
    }

    void recover(const std::vector<double>& rhs, const std::vector<double>& interfaceSolution,
                 std::vector<double>& temperature)
    {
        temperature.assign(static_cast<std::size_t>(system.size()), 0.0);
        for (std::size_t i = 0; i < partition.interfaceGlobalDofs.size(); ++i) {
            temperature[static_cast<std::size_t>(partition.interfaceGlobalDofs[i])] = interfaceSolution[i];
        }
        const auto start = Clock::now();
        for (std::size_t d = 0; d < partition.domains.size(); ++d) {
            const DomainBlocks& domain = partition.domains[d];
            const int nI = static_cast<int>(domain.interiorGlobalDofs.size());
            if (nI == 0) {
                continue;
            }
            std::vector<double> localRhs(static_cast<std::size_t>(nI), 0.0);
            for (int i = 0; i < nI; ++i) {
                localRhs[static_cast<std::size_t>(i)] = rhs[static_cast<std::size_t>(domain.interiorGlobalDofs[static_cast<std::size_t>(i)])];
                for (const auto& entry : domain.interiorInterfaceRows[static_cast<std::size_t>(i)]) {
                    localRhs[static_cast<std::size_t>(i)] -= entry.second * interfaceSolution[static_cast<std::size_t>(entry.first)];
                }
            }
            std::vector<double> localSolution;
            localSolvers[d].solve(localRhs, localSolution);
            for (int i = 0; i < nI; ++i) {
                temperature[static_cast<std::size_t>(domain.interiorGlobalDofs[static_cast<std::size_t>(i)])] = localSolution[static_cast<std::size_t>(i)];
            }
        }
        report.recoverySeconds = std::chrono::duration<double>(Clock::now() - start).count();
    }

    void solve(const std::vector<double>& rhs, std::vector<double>& temperature)
    {
        if (!runnable) {
            throw std::runtime_error("exact Schur direct solver is not runnable: " + report.abortReason);
        }
        if (rhs.size() != static_cast<std::size_t>(system.size())) {
            throw std::runtime_error("exact Schur RHS has the wrong global size");
        }
        const auto start = Clock::now();
        std::vector<double> g;
        condensedRhs(rhs, g);
        std::vector<double> interfaceSolution;
        const auto interfaceStart = Clock::now();
        schurSolver.solve(g, interfaceSolution);
        report.interfaceSolveSeconds = std::chrono::duration<double>(Clock::now() - interfaceStart).count();
        recover(rhs, interfaceSolution, temperature);
        report.peakWorkingSetBytes = std::max(report.peakWorkingSetBytes, peakWorkingSetBytes());
        report.totalSeconds += std::chrono::duration<double>(Clock::now() - start).count();
        report.status = "success";
    }

    void writeReport(const std::filesystem::path& path) const
    {
        std::filesystem::create_directories(path.parent_path());
        std::ofstream out(path);
        if (!out) {
            throw std::runtime_error("cannot write exact Schur report: " + path.string());
        }
        out << "field,value\n";
        out << "schur_direct_status," << report.status << "\n";
        out << "abort_reason,\"" << report.abortReason << "\"\n";
        out << "interface_dofs," << report.interfaceDofs << "\n";
        out << "predicted_upper_nnz," << report.predictedUpperNnz << "\n";
        out << "actual_unique_upper_nnz," << report.actualUniqueUpperNnz << "\n";
        out << "pattern_density," << report.patternDensity << "\n";
        out << "index_bits," << report.indexBits << "\n";
        out << "pattern_seconds," << report.patternSeconds << "\n";
        out << "assembly_seconds," << report.assemblySeconds << "\n";
        out << "rhs_seconds," << report.rhsSeconds << "\n";
        out << "factor_phase11_seconds," << report.factorPhase11Seconds << "\n";
        out << "factor_phase22_seconds," << report.factorPhase22Seconds << "\n";
        out << "local_phase11_seconds," << report.localPhase11Seconds << "\n";
        out << "local_phase22_seconds," << report.localPhase22Seconds << "\n";
        out << "interface_solve_seconds," << report.interfaceSolveSeconds << "\n";
        out << "recovery_seconds," << report.recoverySeconds << "\n";
        out << "total_seconds," << report.totalSeconds << "\n";
        out << "raw_csr_bytes," << report.rawCsrBytes << "\n";
        out << "factor_bytes," << report.factorBytes << "\n";
        out << "local_factor_bytes," << report.localFactorBytes << "\n";
        out << "peak_working_set_bytes," << report.peakWorkingSetBytes << "\n";
        out << "batch_size," << report.batchSize << "\n";
        out << "multi_rhs_calls," << report.multiRhsCalls << "\n";
        out << "local_solve_seconds," << report.localSolveSeconds << "\n";
        out << "operator_relative_error," << report.operatorRelativeError << "\n";
        out << "operator_max_absolute_error," << report.operatorMaxAbsoluteError << "\n";
        out << "interface_coupling_symmetry_error," << report.interfaceCouplingSymmetryError << "\n";
        out << "schur_symmetry_error," << report.schurSymmetryError << "\n";
        out << "rhs_relative_error," << report.rhsRelativeError << "\n";
        out << "true_residual," << report.trueResidual << "\n";
        out << "relative_l2," << report.relativeL2 << "\n";
        out << "max_temperature_difference," << report.maxTemperatureDifference << "\n";
        std::ofstream domains(path.parent_path() / "schur_direct_exact_domains.csv");
        domains << "domain_id,interior_dofs,interface_dofs,clique_upper_nnz,ai_gamma_nnz,estimated_dense_local_schur_bytes,phase11_seconds,phase22_seconds,phase33_calls,phase33_seconds,factor_memory_bytes\n";
        for (const auto& row : report.domains) {
            domains << row.domainId << ',' << row.interiorDofs << ',' << row.interfaceDofs << ','
                    << row.cliqueUpperNnz << ',' << row.aiGammaNnz << ','
                    << row.estimatedDenseLocalSchurBytes << ',' << row.phase11Seconds << ','
                    << row.phase22Seconds << ',' << row.phase33Calls << ',' << row.phase33Seconds
                    << ',' << row.factorMemoryBytes << '\n';
        }
    }
};

ExactSchurDirectSolver::ExactSchurDirectSolver(const Mesh& mesh,
                                               const SparseMatrix& system,
                                               const CaseConfig& physics,
                                               const ExactSchurDirectOptions& options)
    : impl_(std::make_unique<Impl>(mesh, system, physics, options)) {}

ExactSchurDirectSolver::~ExactSchurDirectSolver() = default;
bool ExactSchurDirectSolver::canSolve() const { return impl_->runnable; }
const ExactSchurDirectReport& ExactSchurDirectSolver::report() const { return impl_->report; }
void ExactSchurDirectSolver::solve(const std::vector<double>& rhs, std::vector<double>& temperature)
{
    impl_->solve(rhs, temperature);
}
void ExactSchurDirectSolver::recordTrueResidual(double value)
{
    impl_->report.trueResidual = value;
}
void ExactSchurDirectSolver::writeReport(const std::filesystem::path& path) const
{
    impl_->writeReport(path);
}

} // namespace ddm_schur
