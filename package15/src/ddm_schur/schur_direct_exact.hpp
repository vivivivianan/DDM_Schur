#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

struct Mesh;
struct SparseMatrix;
struct CaseConfig;

namespace ddm_schur {

struct ExactSchurDirectOptions {
    bool patternOnly = false;
    bool dryRun = false;
    int batchSize = 0;
    double memoryFraction = 0.50;
    double factorMemoryFraction = 0.70;
    bool verifyOperator = true;
    int randomChecks = 3;
    std::filesystem::path outputDirectory;
};

struct ExactSchurDirectDomainReport {
    int domainId = -1;
    std::size_t interiorDofs = 0;
    std::size_t interfaceDofs = 0;
    std::uint64_t cliqueUpperNnz = 0;
    std::uint64_t aiGammaNnz = 0;
    std::uint64_t estimatedDenseLocalSchurBytes = 0;
    double phase11Seconds = 0.0;
    double phase22Seconds = 0.0;
    int phase33Calls = 0;
    double phase33Seconds = 0.0;
    std::size_t factorMemoryBytes = 0;
};

struct ExactSchurDirectReport {
    std::string status = "not_run";
    std::string abortReason;
    int interfaceDofs = 0;
    std::uint64_t predictedUpperNnz = 0;
    std::uint64_t actualUniqueUpperNnz = 0;
    double patternDensity = 0.0;
    int indexBits = 32;
    double patternSeconds = 0.0;
    double assemblySeconds = 0.0;
    double rhsSeconds = 0.0;
    double factorPhase11Seconds = 0.0;
    double factorPhase22Seconds = 0.0;
    double localPhase11Seconds = 0.0;
    double localPhase22Seconds = 0.0;
    double interfaceSolveSeconds = 0.0;
    double recoverySeconds = 0.0;
    double totalSeconds = 0.0;
    std::size_t rawCsrBytes = 0;
    std::size_t factorBytes = 0;
    std::size_t localFactorBytes = 0;
    std::size_t peakWorkingSetBytes = 0;
    int batchSize = 0;
    int multiRhsCalls = 0;
    double localSolveSeconds = 0.0;
    double operatorRelativeError = 0.0;
    double operatorMaxAbsoluteError = 0.0;
    double interfaceCouplingSymmetryError = 0.0;
    double schurSymmetryError = 0.0;
    double rhsRelativeError = 0.0;
    double trueResidual = 0.0;
    double relativeL2 = 0.0;
    double maxTemperatureDifference = 0.0;
    std::vector<ExactSchurDirectDomainReport> domains;
};

class ExactSchurDirectSolver {
public:
    ExactSchurDirectSolver(const Mesh& mesh,
                           const SparseMatrix& system,
                           const CaseConfig& physics,
                           const ExactSchurDirectOptions& options);
    ~ExactSchurDirectSolver();

    ExactSchurDirectSolver(const ExactSchurDirectSolver&) = delete;
    ExactSchurDirectSolver& operator=(const ExactSchurDirectSolver&) = delete;

    bool canSolve() const;
    const ExactSchurDirectReport& report() const;
    void solve(const std::vector<double>& rhs, std::vector<double>& temperature);
    void recordTrueResidual(double value);
    void writeReport(const std::filesystem::path& path) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace ddm_schur
