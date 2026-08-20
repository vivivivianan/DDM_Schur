// Copyright (c) 2026.
//
// Production-only entry point for the validated Package15 Dynamic Schur path.
// The research repository exposed dozens of experimental switches. This
// executable intentionally exposes only parameters that do not change the
// selected algorithm: model/config paths, cold versus warm cache use, time
// horizon, thread budget, and optional summary-only FOM validation.

#include "sipg_core.hpp"
#include "config_io.hpp"
#include "mesh_loader.hpp"
#include "fem_assembly.hpp"
#include "linear_solvers.hpp"
#include "mor/transient/transient_workflow.hpp"
#include "mor/transient/local_dynamic_schur.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

namespace fs = std::filesystem;

// These values are not arbitrary defaults. They reproduce the fastest
// validated single-node configuration measured on Package15 large.
constexpr int constructionWorkers = 8;
constexpr int constructionTraceRank = 19;
constexpr double localRankTolerance = 1.0e-6;

struct CommandLine {
    fs::path config = "configs/package15_large.txt";
    fs::path output = "runs/package15_large_cold";
    fs::path cache = "cache/package15_large";
    std::string mode = "cold";
    int steps = 1;
    double timeStep = 0.05;
    int outerThreads = 16;
    int localMklThreads = 1;
    bool validateFomSummary = false;
    bool showHelp = false;
};

[[noreturn]] void optionError(const std::string& message)
{
    throw std::runtime_error("Command-line error: " + message);
}

std::string requireValue(int& index, int argc, char* argv[])
{
    if (index + 1 >= argc) {
        optionError(std::string("missing value after ") + argv[index]);
    }
    return argv[++index];
}

int parsePositiveInt(const std::string& text, const std::string& option)
{
    std::size_t consumed = 0;
    const int value = std::stoi(text, &consumed);
    if (consumed != text.size() || value <= 0) {
        optionError(option + " requires a positive integer");
    }
    return value;
}

double parsePositiveDouble(const std::string& text,
                           const std::string& option)
{
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != text.size() || !(value > 0.0) || !std::isfinite(value)) {
        optionError(option + " requires a finite positive number");
    }
    return value;
}

void printHelp()
{
    std::cout
        << "DynamicSchurProduction - validated production-only solver\n\n"
        << "Usage:\n"
        << "  DynamicSchurProduction [options]\n\n"
        << "Options:\n"
        << "  --config <file>       Case configuration file.\n"
        << "  --output <dir>        Empty output directory for this run.\n"
        << "  --cache <dir>         Descriptor/reference/local-ROM cache.\n"
        << "  --mode cold|warm      Build and save, or load a cached ROM.\n"
        << "  --steps <N>           Number of Backward-Euler steps.\n"
        << "  --dt <seconds>        Constant time-step size.\n"
        << "  --threads <N>         Augmented-factor thread budget.\n"
        << "  --local-mkl <N>       MKL threads per concurrent local ROM.\n"
        << "  --validate-fom        Add summary-only monolithic validation.\n"
        << "  --help                Print this help and exit.\n\n"
        << "Fixed algorithm:\n"
        << "  global-fom traces, M1 Local Block-Arnoldi, full interface,\n"
        << "  augmented-direct, no residual replacement, no full field.\n";
}

CommandLine parseCommandLine(int argc, char* argv[])
{
    CommandLine result;
    for (int index = 1; index < argc; ++index) {
        const std::string option = argv[index];
        if (option == "--help" || option == "-h") {
            result.showHelp = true;
        } else if (option == "--config") {
            result.config = requireValue(index, argc, argv);
        } else if (option == "--output") {
            result.output = requireValue(index, argc, argv);
        } else if (option == "--cache") {
            result.cache = requireValue(index, argc, argv);
        } else if (option == "--mode") {
            result.mode = requireValue(index, argc, argv);
        } else if (option == "--steps") {
            result.steps = parsePositiveInt(
                requireValue(index, argc, argv), option);
        } else if (option == "--dt") {
            result.timeStep = parsePositiveDouble(
                requireValue(index, argc, argv), option);
        } else if (option == "--threads") {
            result.outerThreads = parsePositiveInt(
                requireValue(index, argc, argv), option);
        } else if (option == "--local-mkl") {
            result.localMklThreads = parsePositiveInt(
                requireValue(index, argc, argv), option);
        } else if (option == "--validate-fom") {
            result.validateFomSummary = true;
        } else {
            optionError("unknown option '" + option + "'");
        }
    }

    if (result.mode != "cold" && result.mode != "warm") {
        optionError("--mode must be cold or warm");
    }
    return result;
}

fs::path absoluteNormalized(const fs::path& path)
{
    return fs::absolute(path).lexically_normal();
}

void requireEmptyOutputDirectory(const fs::path& output)
{
    std::error_code error;
    if (fs::exists(output, error)) {
        if (error) {
            throw std::runtime_error("Cannot inspect output directory: "
                                     + output.string());
        }
        if (!fs::is_directory(output) || !fs::is_empty(output)) {
            throw std::runtime_error(
                "Output directory must be absent or empty: " + output.string());
        }
    }
    fs::create_directories(output);
}

void requireWarmCache(const fs::path& cache)
{
    const fs::path required[] = {
        cache / "thermal_descriptor.bin",
        cache / "local_dynamic_reference.bin",
        cache / "local_dynamic_interior_model.bin"};
    for (const fs::path& file : required) {
        if (!fs::is_regular_file(file)) {
            throw std::runtime_error("Warm run is missing cache file: "
                                     + file.string());
        }
    }
}

void setConstructionWorkerLimit()
{
    // solverParallelWorkers() reads this process-local limit. Keeping it below
    // the 16-thread augmented-factor budget prevents memory-bandwidth
    // contention during independent construction and local-ROM tasks.
#ifdef _WIN32
    _putenv_s("SIPG_SOLVER_WORKERS",
              std::to_string(constructionWorkers).c_str());
#else
    setenv("SIPG_SOLVER_WORKERS",
           std::to_string(constructionWorkers).c_str(), 1);
#endif
}

mor::transient::Options productionOptions(const CommandLine& command,
                                           const fs::path& cache)
{
    mor::transient::Options options;
    options.generate = command.mode == "cold";
    options.loadPath = command.mode == "warm" ? cache : fs::path{};
    options.savePath = command.mode == "cold" ? cache : fs::path{};

    // Offline model construction: operator-derived traces and one local
    // Arnoldi block. No snapshot POD or data-trained interface basis is used.
    options.moments = 1;
    options.expansionPoint = 0.0;
    options.rankTolerance = localRankTolerance;
    options.secondMomentEnergy = 1.0;
    options.secondMomentMaximumColumns = 0;
    options.constructionTraceMode = "global-fom";
    options.interfaceExcitationRank = constructionTraceRank;
    options.reuseIdenticalSubdomains = true;

    // Online integration and validation policy. Full-order validation is
    // summary-only and never participates in basis generation or replacement.
    options.timeStep = command.timeStep;
    options.endTime = command.steps * command.timeStep;
    options.integrator = "backward-euler";
    options.waveform = "asynchronous_hotspots";
    options.outputMode = "max-temperature";
    options.initialMode = "ambient";
    options.compareFom = command.validateFomSummary;
    options.compareFomSummaryOnly = command.validateFomSummary;
    options.nativeReducedHistory = true;
    options.fullResidualTolerance = 1.0e-4;
    options.fullResidualFallback = false;

    // Fastest measured interface backend. The full physical interface remains
    // explicit; PARDISO eliminates it together with all local ROM coordinates.
    options.interfaceKrylov = "augmented-direct";
    options.matrixFreeInterfaceThreshold = 0;
    options.interfaceInitialGuess = "previous";
    options.interfaceTolerance = 1.0e-10;
    options.adaptiveInterfaceTolerance = 1.0e-9;
    options.localSolveThreads = command.outerThreads;
    options.localPardisoThreads = command.localMklThreads;
    options.portReduction = false;
    options.proxyEnabled = false;
    options.proxyCacheEnabled = false;
    options.portCoreCacheEnabled = false;
    return options;
}

void writeProductionMetadata(const fs::path& output,
                             const CommandLine& command,
                             const fs::path& config,
                             const fs::path& cache)
{
    std::ofstream stream(output / "production_run.csv");
    stream << "field,value\n"
           << "source_commit," << DDM_SCHUR_SOURCE_COMMIT << "\n"
           << "config," << config.string() << "\n"
           << "cache," << cache.string() << "\n"
           << "mode," << command.mode << "\n"
           << "steps," << command.steps << "\n"
           << "dt_seconds," << std::setprecision(17) << command.timeStep << "\n"
           << "construction_workers," << constructionWorkers << "\n"
           << "augmented_factor_threads," << command.outerThreads << "\n"
           << "local_mkl_threads," << command.localMklThreads << "\n"
           << "fom_summary_validation," << (command.validateFomSummary ? 1 : 0)
           << "\n";
}

} // namespace

int main(int argc, char* argv[])
{
    try {
        std::cout << std::unitbuf;
        std::cerr << std::unitbuf;
        const CommandLine command = parseCommandLine(argc, argv);
        if (command.showHelp) {
            printHelp();
            return 0;
        }

        const fs::path configPath = absoluteNormalized(command.config);
        const fs::path outputPath = absoluteNormalized(command.output);
        const fs::path cachePath = absoluteNormalized(command.cache);
        if (!fs::is_regular_file(configPath)) {
            throw std::runtime_error("Configuration file not found: "
                                     + configPath.string());
        }
        requireEmptyOutputDirectory(outputPath);
        if (command.mode == "warm") {
            requireWarmCache(cachePath);
        } else {
            fs::create_directories(cachePath);
        }

        setConstructionWorkerLimit();
        const CaseConfig physics = readCaseConfig(configPath);
        Mesh mesh = loadComsolDomainDecompositionMesh(physics);

        const mor::transient::Options options =
            productionOptions(command, cachePath);
        writeProductionMetadata(
            outputPath, command, configPath, cachePath);

        std::cout << "Dynamic Schur production path\n"
                  << "  mode=" << command.mode
                  << ", dofs=" << mesh.nodes.size()
                  << ", steps=" << command.steps
                  << ", dt=" << command.timeStep << " s\n"
                  << "  algorithm=global-fom/M1/full-interface/augmented-direct\n"
                  << "  construction_workers=" << constructionWorkers
                  << ", factor_threads=" << command.outerThreads
                  << ", local_mkl_threads=" << command.localMklThreads << "\n";

        mor::transient::runLocalBlockArnoldiDynamicSchurWorkflow(
            mesh, physics, options, outputPath);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
