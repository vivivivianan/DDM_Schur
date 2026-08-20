#pragma once

// Top-level application workflow. main.cpp intentionally stays as a small executable entry point.

static int runProgram(int argc, char* argv[])
{
    try {
        std::cout << std::unitbuf;
        std::cerr << std::unitbuf;
        const auto programStart = std::chrono::steady_clock::now();
        ProgramTiming timing;
        ProgramOptions options = parseProgramOptions(argc, argv);
        if (options.showHelp) {
            return 0;
        }
        if (options.morTransientDeploymentOnly) {
            if (options.morTransientLoadPath.empty()) {
                throw std::runtime_error(
                    "--mor-transient-deployment-only requires --mor-transient-load <directory>.");
            }
            if (options.morTransientMethod != "block-arnoldi") {
                throw std::runtime_error(
                    "Stage 2C.1 supports --mor-transient-method block-arnoldi.");
            }
            const std::filesystem::path deploymentOutput = options.outputDirOverride
                ? std::filesystem::absolute(options.outputDir).lexically_normal()
                : std::filesystem::absolute(options.morTransientLoadPath)
                    .lexically_normal() / "deployment_output";
            mor::transient::Options transientOptions;
            transientOptions.loadPath = options.morTransientLoadPath;
            transientOptions.timeStep = options.morTransientDt;
            transientOptions.endTime = options.morTransientEndTime;
            transientOptions.integrator = options.morTransientIntegrator;
            transientOptions.inputPath = options.morTransientInputPath;
            transientOptions.waveform = options.morTransientWaveform;
            transientOptions.outputMode = options.morTransientOutput;
            transientOptions.initialTemperature = options.morTransientInitialTemperature;
            transientOptions.deploymentRhsCount = options.morDeploymentRhsCount;
            transientOptions.seed = options.morSeed;
            mor::transient::runTransientDeploymentOnly(
                std::filesystem::absolute(options.morTransientLoadPath).lexically_normal(),
                transientOptions, deploymentOutput);
            return 0;
        }
        if (options.morParametricDeploymentOnly) {
            if (options.morParametricLoadPath.empty()) {
                throw std::runtime_error(
                    "--mor-parametric-deployment-only requires --mor-parametric-load <directory>.");
            }
            const std::filesystem::path deploymentOutput = options.outputDirOverride
                ? std::filesystem::absolute(options.outputDir).lexically_normal()
                : std::filesystem::absolute(options.morParametricLoadPath)
                    .lexically_normal() / "deployment_output";
            mor::Options deploymentOptions;
            deploymentOptions.parameterValue = options.morParameterValue;
            deploymentOptions.onlinePowersW = options.morOnlinePowersW;
            deploymentOptions.allowExtrapolation = options.morAllowExtrapolation;
            deploymentOptions.interfaceRank = options.morInterfaceRank;
            deploymentOptions.localRank = options.morLocalRank;
            deploymentOptions.deploymentRhsCount = options.morDeploymentRhsCount;
            deploymentOptions.seed = options.morSeed;
            deploymentOptions.saveModelPath = options.morParametricSavePath;
            mor::parametric::runParametricDeploymentOnly(
                std::filesystem::absolute(options.morParametricLoadPath).lexically_normal(),
                deploymentOutput,
                deploymentOptions);
            return 0;
        }
        if (options.morDeploymentOnly) {
            if (options.morLoadModelPath.empty()) {
                throw std::runtime_error(
                    "--mor-deployment-only requires --mor-load-model <directory>.");
            }
            const std::filesystem::path deploymentOutput = options.outputDirOverride
                ? std::filesystem::absolute(options.outputDir).lexically_normal()
                : std::filesystem::absolute(options.morLoadModelPath)
                    .lexically_normal() / "deployment_output";
            mor::Options deploymentOptions;
            deploymentOptions.interiorMode = options.morInteriorMode;
            deploymentOptions.interiorRank = options.morInteriorRank;
            deploymentOptions.interiorEnergyTolerance = options.morInteriorEnergyTolerance;
            deploymentOptions.interiorSingularValueTolerance =
                options.morInteriorSingularValueTolerance;
            deploymentOptions.interiorRankSweep = options.morInteriorRankSweep;
            deploymentOptions.saveModelPath = options.morSaveModelPath;
            deploymentOptions.storagePrecision = options.morStoragePrecision;
            deploymentOptions.reportIoTime = options.morReportIoTime;
            deploymentOptions.compareInteriorModes = options.morCompareInteriorModes;
            deploymentOptions.deploymentRhsCount = options.morDeploymentRhsCount;
            deploymentOptions.seed = options.morSeed;
            mor::runDeploymentOnly(
                std::filesystem::absolute(options.morLoadModelPath).lexically_normal(),
                deploymentOutput,
                deploymentOptions);
            return 0;
        }
        const AnalysisMode analysisMode = options.mode;
        const bool steady = analysisMode == AnalysisMode::Steady;
        const std::string modeName = analysisModeName(analysisMode);
        const std::filesystem::path projectRoot = options.configPath.empty()
            ? findProjectRoot()
            : std::filesystem::absolute(options.configPath).parent_path();
        const std::filesystem::path configPath = options.configPath.empty()
            ? projectRoot / "case_config.txt"
            : std::filesystem::absolute(options.configPath).lexically_normal();
        CaseConfig physics = readCaseConfig(configPath);
        if (options.solverMethodOverride) {
            if (!isKnownSolverMethod(options.solverMethod)) {
                throw std::runtime_error("--solver-method must be monolithic, schwarz, schwarz-precond-fgmres, schwarz-precond-fgmres-two-level, or schur-direct-exact.");
            }
            physics.solverMethod = options.solverMethod;
            physics.schwarz.enabled = options.solverMethod == "schwarz"
                || options.solverMethod == "schwarz_precond_fgmres"
                || options.solverMethod == "schwarz_precond_fgmres_two_level";
        }
        if (options.schwarzTypeOverride) {
            if (options.schwarzType != "multiplicative" && options.schwarzType != "additive") {
                throw std::runtime_error("--schwarz-type must be multiplicative or additive.");
            }
            physics.schwarz.type = options.schwarzType;
        }
        if (options.schwarzStandaloneModeOverride) {
            options.schwarzStandaloneMode =
                normalizeSchwarzTransmissionModeName(options.schwarzStandaloneMode);
            if (!isKnownSchwarzStandaloneMode(options.schwarzStandaloneMode)) {
                throw std::runtime_error(
                    "--schwarz-standalone-mode must be algebraic, dirichlet_neumann, dirichlet_dirichlet, dirichlet_robin, or robin.");
            }
            physics.schwarz.standaloneMode = options.schwarzStandaloneMode;
            physics.schwarz.transmission =
                schwarzTransmissionForStandaloneMode(options.schwarzStandaloneMode);
        }
        if (options.schwarzTransmissionOrientationOverride) {
            options.schwarzTransmissionOrientation =
                normalizeSchwarzTransmissionOrientationName(options.schwarzTransmissionOrientation);
            if (!isKnownSchwarzTransmissionOrientation(options.schwarzTransmissionOrientation)) {
                throw std::runtime_error("--schwarz-transmission-orientation must be forward or reverse.");
            }
            physics.schwarz.transmissionOrientation = options.schwarzTransmissionOrientation;
        }
        if (options.schwarzRobinAlphaFactorOverride) {
            physics.schwarz.robinAlphaFactor = options.schwarzRobinAlphaFactor;
        }
        if (options.schwarzFluxEvalOverride) {
            if (options.schwarzFluxEval != "physical_gradient" && options.schwarzFluxEval != "sipg_numeric") {
                throw std::runtime_error("--schwarz-flux-eval must be physical_gradient or sipg_numeric.");
            }
            physics.schwarz.fluxEval = options.schwarzFluxEval;
        }
        if (options.schwarzOverlapLayersOverride) {
            if (options.schwarzOverlapLayers < 0 || options.schwarzOverlapLayers > 3) {
                throw std::runtime_error("--schwarz-overlap-layers must be 0, 1, 2, or 3.");
            }
            physics.schwarz.overlapLayers = options.schwarzOverlapLayers;
        }
        if (options.schwarzOverlapModeOverride) {
            if (options.schwarzOverlapMode != "none"
                && options.schwarzOverlapMode != "halo"
                && options.schwarzOverlapMode != "ras") {
                throw std::runtime_error("--schwarz-overlap-mode must be none, halo, or ras.");
            }
            physics.schwarz.overlapMode = options.schwarzOverlapMode;
        }
        if (options.schwarzPartitionModeOverride) {
            if (!isKnownSchwarzPartitionMode(options.schwarzPartitionMode)) {
                throw std::runtime_error(
                    "--schwarz-partition-mode must be current, material-aligned, hotspot-contained, or vertical-heat-flow-aligned.");
            }
            physics.schwarz.partitionMode = options.schwarzPartitionMode;
        }
        if (physics.schwarz.overlapLayers == 0) {
            physics.schwarz.overlapMode = "none";
        } else if (physics.schwarz.overlapMode == "none") {
            physics.schwarz.overlapMode = "ras";
        }
        if (options.schwarzWriteInterfaceFluxOverride) {
            physics.schwarz.writeInterfaceFlux = options.schwarzWriteInterfaceFlux;
        }
        if (options.schwarzMaxItersOverride) {
            physics.schwarz.maxIters = options.schwarzMaxIters;
        }
        if (options.schwarzTolRelUpdateOverride) {
            physics.schwarz.tolRelUpdate = options.schwarzTolRelUpdate;
        }
        if (options.schwarzTolRelResidualOverride) {
            physics.schwarz.tolRelResidual = options.schwarzTolRelResidual;
        }
        if (options.schwarzRelaxationOverride) {
            physics.schwarz.relaxation = options.schwarzRelaxation;
        }
        if (physics.schwarz.enabled
            && physics.solverMethod != "schwarz_precond_fgmres"
            && physics.solverMethod != "schwarz_precond_fgmres_two_level") {
            physics.solverMethod = "schwarz";
        }
        if (physics.solverMethod == "schwarz" && !options.solversExplicit) {
            options.runSchwarz = true;
            options.runBjIc = false;
            options.runDirect = physics.schwarz.validateAgainstMonolithic;
        }
        if (physics.solverMethod == "schwarz_precond_fgmres" && !options.solversExplicit) {
            options.runSchwarzPrecondFgmres = true;
            options.runSchwarz = false;
            options.runBjIc = false;
            options.runDirect = physics.schwarz.validateAgainstMonolithic;
        }
        if (physics.solverMethod == "schwarz_precond_fgmres_two_level" && !options.solversExplicit) {
            options.runSchwarzPrecondFgmresTwoLevel = true;
            options.runSchwarzPrecondFgmres = false;
            options.runSchwarz = false;
            options.runBjIc = false;
            options.runDirect = physics.schwarz.validateAgainstMonolithic;
        }
        if (physics.solverMethod == "schur_direct_exact" && !options.solversExplicit) {
            options.runSchurDirectExact = true;
            options.runBjIc = false;
            options.runDirect = false;
        }
        if (options.runSchwarz) {
            if (!options.runSchwarzPrecondFgmres && !options.runSchwarzPrecondFgmresTwoLevel) {
                physics.solverMethod = "schwarz";
            }
            physics.schwarz.enabled = true;
        }
        if (options.runSchwarzPrecondFgmres) {
            physics.solverMethod = "schwarz_precond_fgmres";
            physics.schwarz.enabled = true;
        }
        if (options.runSchwarzPrecondFgmresTwoLevel) {
            physics.solverMethod = "schwarz_precond_fgmres_two_level";
            physics.schwarz.enabled = true;
        }
        if (options.penaltyModeOverride) {
            if (options.penaltyMode != "harmonic" && options.penaltyMode != "max") {
                throw std::runtime_error("--penalty-mode must be harmonic or max.");
            }
            physics.penaltyMode = options.penaltyMode;
        }
        if (options.penaltyFactorOverride) {
            physics.penaltyFactor = options.penaltyFactor;
        }
        if (options.dirichletMethodOverride) {
            if (options.dirichletMethod != "strong" && options.dirichletMethod != "nitsche") {
                throw std::runtime_error("--dirichlet-method must be strong or nitsche.");
            }
            physics.dirichletMethod = options.dirichletMethod;
        }
        if (options.nitschePenaltyFactorOverride) {
            physics.nitschePenaltyFactor = options.nitschePenaltyFactor;
        }
        if (options.initialGuessOverride) {
            if (options.initialGuessType != "current"
                && options.initialGuessType != "zero"
                && options.initialGuessType != "constant"
                && options.initialGuessType != "random-small"
                && options.initialGuessType != "random_small") {
                throw std::runtime_error("--initial-guess must be zero, current, constant, or random-small.");
            }
            if (options.initialGuessType == "random_small") {
                options.initialGuessType = "random-small";
            }
            physics.initialGuessType = options.initialGuessType;
        }
        if (options.disableWarmStartOverride) {
            physics.disableWarmStart = options.disableWarmStart;
        }
        if (options.forceNontrivialRhsOverride) {
            physics.forceNontrivialRhs = options.forceNontrivialRhs;
        }
        if (options.thermalSourceScaleOverride) {
            physics.thermalSourceScale = options.thermalSourceScale;
        }
        if (physics.interfaceScheme == "nipg") {
            throw std::runtime_error(
                "NIPG produces a nonsymmetric system and is not supported by the current "
                "SPD Schur/PARDISO and balanced two-level preconditioner. "
                "Use interface_scheme=sipg.");
        }
        if (options.localIcShiftMode != "none"
            && options.localIcShiftMode != "auto"
            && options.localIcShiftMode != "fixed") {
            throw std::runtime_error("--local-ic-shift-mode must be none, auto, or fixed.");
        }
        if (options.coarseCorrection != "additive" && options.coarseCorrection != "multiplicative") {
            throw std::runtime_error("--coarse-correction must be additive or multiplicative.");
        }
        if (options.coarseSpace != "subdomain_constant"
            && options.coarseSpace != "subdomain_material_constant"
            && options.coarseSpace != "interface_constant") {
            throw std::runtime_error("--coarse-space must be subdomain-constant, subdomain-material-constant, or interface-constant.");
        }
        const std::filesystem::path outputDir = options.outputDirOverride
            ? std::filesystem::absolute(options.outputDir).lexically_normal()
            : (physics.outputDir.empty() ? projectRoot / "output" : physics.outputDir);
        std::filesystem::create_directories(outputDir);
        const std::string modelName = options.modelNameOverride ? options.modelName : physics.name;
        const std::string modelPrefix = safeFilePrefix(modelName);
        writeRunMetadata(outputDir / "run_metadata.csv", modelName, configPath, outputDir, physics);
        const bool collectDetailedDiagnostics = !options.fastRun;

        Mesh mesh = loadComsolDomainDecompositionMesh(physics);
        const int n = static_cast<int>(mesh.nodes.size());
        AssemblyDiagnostics assemblyDiagnostics;
        if (collectDetailedDiagnostics) {
            initializeAssemblyDiagnostics(n, assemblyDiagnostics);
        }
        markInterfaceDofs(mesh, assemblyDiagnostics);

        if (options.morTransientGenerate
            || !options.morTransientLoadPath.empty()) {
            if (options.morTransientMethod != "block-arnoldi"
                && options.morTransientMethod != "local-block-arnoldi"
                && options.morTransientMethod != "local-port-block-arnoldi"
                && options.morTransientMethod != "local-interior-arnoldi-reduced-schur") {
                throw std::runtime_error(
                    "Transient MOR supports block-arnoldi, local-block-arnoldi, "
                    "local-port-block-arnoldi, or local-interior-arnoldi-reduced-schur.");
            }
            if (options.portBasisMethod != "full-interface"
                && options.portBasisMethod != "port-pod"
                && options.portBasisMethod != "steklov-schur"
                && options.portBasisMethod != "mandatory-only"
                && options.portBasisMethod != "residual-krylov"
                && options.portBasisMethod != "optimal-transfer"
                && options.portBasisMethod != "randomized-transfer"
                && options.portBasisMethod != "hybrid-randomized") {
                throw std::runtime_error(
                    "--port-basis-method is unsupported.");
            }
            if (options.optimalPortRankMode != "fixed"
                && options.optimalPortRankMode != "match-m7"
                && options.optimalPortRankMode != "eigenvalue-tolerance") {
                throw std::runtime_error(
                    "--optimal-port-rank-mode must be fixed, match-m7, or eigenvalue-tolerance.");
            }
            if (options.optimalPortInnerProduct != "trace-mass"
                && options.optimalPortInnerProduct != "penalty-weighted-mass") {
                throw std::runtime_error(
                    "--optimal-port-inner-product must be trace-mass or penalty-weighted-mass.");
            }
            if (options.optimalPortInnerSolver != "auto"
                && options.optimalPortInnerSolver != "direct"
                && options.optimalPortInnerSolver != "pcg"
                && options.optimalPortInnerSolver != "fgmres"
                && options.optimalPortInnerSolver != "iterative-schur"
                && options.optimalPortInnerSolver != "assembled-dense"
                && options.optimalPortInnerSolver != "woodbury-exact") {
                throw std::runtime_error(
                    "--optimal-port-inner-solver must be "
                    "iterative-schur, assembled-dense, woodbury-exact, "
                    "or a legacy compatibility value (auto/direct/pcg/fgmres).");
            }
            if (options.optimalPortInnerRefinementMaximumIterations < 0
                || options.optimalPortInnerRefinementMaximumIterations > 3
                || !(options.optimalPortInnerRefinementTolerance > 0.0)) {
                throw std::runtime_error(
                    "Optimal-port outer target refinement requires "
                    "0..3 iterations and a positive tolerance.");
            }
            if (options.optimalPortAblation != "mandatory-transfer"
                && options.optimalPortAblation != "mandatory-only"
                && options.optimalPortAblation != "transfer-only"
                && options.optimalPortAblation != "constant-only"
                && options.optimalPortAblation != "geometry-particular"
                && options.optimalPortAblation != "constant-geometry"
                && options.optimalPortAblation != "input-particular-only"
                && options.optimalPortAblation != "trace-transfer-only"
                && options.optimalPortAblation != "generalized-transfer-only"
                && options.optimalPortAblation != "constant-geometry-trace"
                && options.optimalPortAblation != "constant-geometry-generalized"
                && options.optimalPortAblation != "original-mandatory-trace") {
                throw std::runtime_error(
                    "--optimal-port-ablation is unsupported.");
            }
            if (options.optimalPortSourceMode != "trace-only"
                && options.optimalPortSourceMode != "trace-plus-input"
                && options.optimalPortSourceMode != "generalized-dynamic") {
                throw std::runtime_error(
                    "--optimal-port-source-mode must be trace-only, trace-plus-input, or generalized-dynamic.");
            }
            const int optimalPortStopModes =
                (options.optimalPortTopologyAudit ? 1 : 0)
                + (options.optimalPortBasisPilot ? 1 : 0)
                + (options.optimalPortTargetSolverComparison ? 1 : 0)
                + (options.optimalPortWoodburyPilot ? 1 : 0)
                + (options.optimalPortRefinementValidation ? 1 : 0)
                + (options.optimalPortRepresentativeInterfacePilot ? 1 : 0)
                + (options.optimalPortMaximumInterfaceRefinementPilot
                    ? 1 : 0)
                + (options.optimalPortAllInterfaceBasis ? 1 : 0);
            if (optimalPortStopModes > 1) {
                throw std::runtime_error(
                    "Optimal-port audit/comparison/pilot stop modes are mutually exclusive.");
            }
            if ((options.optimalPortRepresentativeInterfacePilot
                    || options.optimalPortMaximumInterfaceRefinementPilot
                    || options.optimalPortAllInterfaceBasis)
                && options.optimalPortTopologyAuditCsv.empty()) {
                throw std::runtime_error(
                    "RRAM26 basis scalability modes require "
                    "--optimal-port-topology-audit-csv.");
            }
            const bool residualKrylov =
                options.portBasisMethod == "residual-krylov"
                || options.portBasisMethod == "mandatory-only"
                || options.portBasisMethod == "hybrid-randomized";
            const std::set<int> residualRanks = {0, 1, 2, 4, 8};
            if (residualKrylov
                && (residualRanks.count(
                        options.residualKrylovMaximumRank) == 0
                    || options.residualKrylovMaximumSweeps < 1
                    || options.residualKrylovMaximumSweeps > 3
                    || !(options.residualKrylovTolerance > 0.0)
                    || options.residualKrylovBlockSize <= 0
                    || (options.residualKrylovProbeMode
                            != "operator-geometry"
                        && options.residualKrylovProbeMode
                            != "particular-only")
                    || options.residualKrylovInnerSolver
                        != "woodbury-exact")) {
                throw std::runtime_error(
                    "Residual-Krylov options are invalid.");
            }
            if ((options.residualKrylovRepresentativePilot
                    || options.residualKrylovAllInterfaceBasis)
                && (!residualKrylov
                    || options.optimalPortTopologyAuditCsv.empty())) {
                throw std::runtime_error(
                    "Residual-Krylov scalability modes require the "
                    "matching basis method and topology-audit CSV.");
            }
            if (options.residualKrylovRepresentativePilot
                && options.residualKrylovAllInterfaceBasis) {
                throw std::runtime_error(
                    "Residual-Krylov stop modes are mutually exclusive.");
            }
            const bool randomizedTransfer =
                options.portBasisMethod == "randomized-transfer"
                || options.portBasisMethod == "hybrid-randomized";
            const std::set<int> randomizedRanks = {8, 16, 32};
            if (randomizedTransfer
                && (randomizedRanks.count(options.randomizedPortRank) == 0
                    || options.randomizedPortOversampling < 0
                    || (options.randomizedPortPowerIterations != 0
                        && options.randomizedPortPowerIterations != 1))) {
                throw std::runtime_error(
                    "Randomized-transfer rank/power/oversampling options are invalid.");
            }
            if (options.randomizedPortCompareOptimal
                && options.portBasisMethod
                    != "randomized-transfer") {
                throw std::runtime_error(
                    "Randomized-versus-Optimal comparison requires "
                    "--port-basis-method randomized-transfer.");
            }
            if (options.randomizedPortRepresentativePilot
                && !randomizedTransfer) {
                throw std::runtime_error(
                    "Randomized-transfer comparison/pilot flags require "
                    "randomized-transfer or hybrid-randomized.");
            }
            if (options.randomizedPortRepresentativePilot
                && options.optimalPortTopologyAuditCsv.empty()) {
                throw std::runtime_error(
                    "Randomized-transfer representative pilot requires "
                    "--optimal-port-topology-audit-csv.");
            }
            if (options.historyCompressionMethod != "none"
                && options.historyCompressionMethod != "deterministic-rrqr") {
                throw std::runtime_error(
                    "--history-compression-method must be none or "
                    "deterministic-rrqr.");
            }
            const std::set<int> historyCompressionRanks =
                {16, 32, 64, 128, 256};
            if (options.historyCompressionMethod != "none"
                && historyCompressionRanks.count(
                    options.historyCompressionRank) == 0) {
                throw std::runtime_error(
                    "Enabled history compression requires "
                    "--history-compression-rank 16, 32, 64, 128, or 256.");
            }
            if (options.historyCompressionMethod != "none"
                && !residualKrylov) {
                throw std::runtime_error(
                    "History compression applies only to mandatory-only, "
                    "residual-krylov, or hybrid-randomized port bases.");
            }
            if (!(options.historyCompressionTolerance > 0.0)) {
                throw std::runtime_error(
                    "--history-compression-tolerance must be positive.");
            }
            if (options.historyCompressionMaximumInterfacePilot
                && (options.portBasisMethod != "hybrid-randomized"
                    || options.historyCompressionMethod
                        != "deterministic-rrqr"
                    || options.optimalPortTopologyAuditCsv.empty())) {
                throw std::runtime_error(
                    "History-compression maximum-interface pilot requires "
                    "hybrid-randomized, deterministic-rrqr, "
                    "and --optimal-port-topology-audit-csv.");
            }
            if (options.historyCompressionMaximumInterfacePilot
                && (optimalPortStopModes > 0
                    || options.residualKrylovRepresentativePilot
                    || options.residualKrylovAllInterfaceBasis
                    || options.randomizedPortRepresentativePilot
                    || options.adaptivePortLocalPilot)) {
                throw std::runtime_error(
                    "History-compression maximum-interface pilot conflicts "
                    "with another audit/pilot/all-interface stop mode.");
            }
            if (options.milestone8ProductionBasisOnly
                && (options.portBasisMethod != "hybrid-randomized"
                    || options.historyCompressionMethod
                        != "deterministic-rrqr"
                    || options.historyCompressionRank != 64
                    || !options.morTransientGenerate
                    || !options.morTransientLoadPath.empty())) {
                throw std::runtime_error(
                    "M8 production basis-only mode requires a fresh "
                    "hybrid-randomized build with deterministic history "
                    "compression rank 64.");
            }
            if (options.milestone8AdaptiveProduction
                && (options.portBasisMethod != "hybrid-randomized"
                    || options.historyCompressionMethod
                        != "deterministic-rrqr"
                    || options.historyCompressionRank != 64
                    || options.randomizedPortRank != 16
                    || options.residualKrylovMaximumRank != 4)) {
                throw std::runtime_error(
                    "M8.9 adaptive production requires the M8.7 default "
                    "CLI ranks 64/16/4; the frozen M8.8 policy applies "
                    "the five per-interface overrides internally.");
            }
            if (options.milestone8ProductionBasisOnly
                && (optimalPortStopModes > 0
                    || options.residualKrylovRepresentativePilot
                    || options.residualKrylovAllInterfaceBasis
                    || options.randomizedPortRepresentativePilot
                    || options.historyCompressionMaximumInterfacePilot
                    || options.adaptivePortLocalPilot
                    || options.morDeploymentRhsCount != 1)) {
                throw std::runtime_error(
                    "M8 production basis-only mode conflicts with another "
                    "pilot/benchmark mode.");
            }
            if (options.adaptivePortLocalPilot
                && (options.portBasisMethod != "hybrid-randomized"
                    || options.historyCompressionMethod
                        != "deterministic-rrqr"
                    || options.adaptivePortInterfaceIds.empty()
                    || options.optimalPortTopologyAuditCsv.empty()
                    || !options.morTransientGenerate
                    || !options.morTransientLoadPath.empty())) {
                throw std::runtime_error(
                    "Adaptive-port local pilot requires a fresh "
                    "hybrid-randomized build, deterministic history "
                    "compression, a nonnegative interface id, and the "
                    "topology-audit CSV.");
            }
            if (options.adaptivePortLocalPilot) {
                std::set<int> adaptiveInterfaceIds;
                for (const int interfaceId :
                     options.adaptivePortInterfaceIds) {
                    if (interfaceId < 0
                        || !adaptiveInterfaceIds.insert(
                            interfaceId).second) {
                        throw std::runtime_error(
                            "Adaptive-port interface ids must be "
                            "nonnegative and unique.");
                    }
                }
            }
            if (options.adaptivePortLocalPilot
                && (optimalPortStopModes > 0
                    || options.residualKrylovRepresentativePilot
                    || options.residualKrylovAllInterfaceBasis
                    || options.randomizedPortRepresentativePilot
                    || options.historyCompressionMaximumInterfacePilot
                    || options.milestone8ProductionBasisOnly
                    || options.morDeploymentRhsCount != 1)) {
                throw std::runtime_error(
                    "Adaptive-port local pilot conflicts with another "
                    "pilot/benchmark mode.");
            }
            if (options.globalInterfaceCoarsePrototype) {
                const std::set<int> supportedRanks{4, 8, 16, 32};
                std::set<int> selected;
                for (int interfaceId :
                     options.globalInterfaceCoarseInterfaceIds) {
                    if (interfaceId < 0
                        || !selected.insert(interfaceId).second) {
                        throw std::runtime_error(
                            "Global-coarse interface ids must be "
                            "nonnegative and unique.");
                    }
                }
                if (options.portBasisMethod != "hybrid-randomized"
                    || options.historyCompressionMethod
                        != "deterministic-rrqr"
                    || options.historyCompressionRank != 64
                    || options.morTransientLoadPath.empty()
                    || selected.empty()
                    || supportedRanks.count(
                        options.globalInterfaceCoarseRank) == 0
                    || (options.globalInterfaceCoarseInverseMode
                            != "schur-jacobi"
                        && options.globalInterfaceCoarseInverseMode
                            != "exact-pcg")
                    || options.globalInterfaceCoarseMaximumIterations <= 0
                    || options
                        .globalInterfaceCoarseInnerMaximumIterations <= 0
                    || options.globalInterfaceCoarseKrylovSweeps <= 0
                    || !(options.globalInterfaceCoarseTolerance > 0.0)
                    || !(options
                        .globalInterfaceCoarseInnerTolerance > 0.0)
                    || (options.globalInterfaceCoarseCandidateDimension > 0
                        && options.globalInterfaceCoarseCandidateDimension
                            < options.globalInterfaceCoarseRank)
                    || options.milestone8ProductionBasisOnly
                    || options.adaptivePortLocalPilot
                    || optimalPortStopModes > 0
                    || options.morDeploymentRhsCount != 1) {
                    throw std::runtime_error(
                        "M8.10.1 global coarse prototype requires the "
                        "loaded hybrid model, deterministic "
                        "history rank 64, rank 4/8/16/32, unique selected "
                        "interfaces, and a single validation RHS.");
                }
            }
            if (options.globalRandomizedSchur) {
                const std::set<int> twoCubeRanks{10, 20, 50};
                const std::set<int> tenCubeRanks{25, 50, 100};
                const std::set<int> rramRanks{50, 100, 200};
                const bool supportedRank =
                    twoCubeRanks.count(options.globalRandomizedRank) != 0
                    || tenCubeRanks.count(options.globalRandomizedRank) != 0
                    || rramRanks.count(options.globalRandomizedRank) != 0;
                if (options.globalInterfaceCoarsePrototype
                    || options.portBasisMethod != "hybrid-randomized"
                    || options.historyCompressionMethod
                        != "deterministic-rrqr"
                    || options.historyCompressionRank != 64
                    || options.morTransientLoadPath.empty()
                    || !supportedRank
                    || options.globalRandomizedInnerMaximumIterations <= 0
                    || !(options.globalRandomizedInnerTolerance > 0.0)
                    || (options.globalRandomizedComposition
                            != "global-only"
                        && options.globalRandomizedComposition
                            != "augment-local")
                    || options.milestone8ProductionBasisOnly
                    || options.adaptivePortLocalPilot
                    || optimalPortStopModes > 0
                    || options.morDeploymentRhsCount != 1) {
                    throw std::runtime_error(
                        "M8.11 global randomized Schur requires a loaded "
                        "hybrid model, deterministic history rank 64, "
                        "an approved rank, one RHS, and composition "
                        "global-only or augment-local.");
                }
            }
            if (options.projectionDiagnosis) {
                std::set<int> uniqueProjectionInterfaces(
                    options.projectionInterfaceIds.begin(),
                    options.projectionInterfaceIds.end());
                if (options.globalInterfaceCoarsePrototype
                    || options.globalRandomizedSchur
                    || options.portBasisMethod != "hybrid-randomized"
                    || options.historyCompressionMethod
                        != "deterministic-rrqr"
                    || options.historyCompressionRank != 64
                    || (options.morTransientLoadPath.empty()
                        && !options.fluxOperatorAudit)
                    || options.projectionInterfaceIds.empty()
                    || uniqueProjectionInterfaces.size()
                        != options.projectionInterfaceIds.size()
                    || options.milestone8ProductionBasisOnly
                    || options.adaptivePortLocalPilot
                    || optimalPortStopModes > 0
                    || options.morDeploymentRhsCount != 1
                    || !(options.morTransientDt > 0.0)
                    || std::abs(
                        options.morTransientEndTime
                        - options.morTransientDt)
                        > 64.0 * std::numeric_limits<double>::epsilon()
                            * std::max({
                                1.0,
                                std::abs(options.morTransientEndTime),
                                std::abs(options.morTransientDt)})) {
                    throw std::runtime_error(
                        "M8.12 projection diagnosis requires a loaded "
                        "hybrid model, deterministic history rank 64, "
                        "unique selected interfaces, and exactly one "
                        "diagnostic time step.");
                }
            }
            if (options.fluxAwarePort
                && (options.fluxAwareFluxType != "numerical"
                    && options.fluxAwareFluxType != "physical"
                    && options.fluxAwareFluxType != "both")) {
                throw std::runtime_error(
                    "Flux-aware flux type must be numerical, "
                    "physical, or both.");
            }
            mor::transient::Options transientOptions;
            transientOptions.generate = options.morTransientGenerate;
            transientOptions.loadPath = options.morTransientLoadPath;
            transientOptions.savePath = options.morTransientSavePath;
            transientOptions.moments = options.morArnoldiMoments;
            transientOptions.expansionPoint = options.morArnoldiExpansionPoint;
            transientOptions.rankTolerance = options.morArnoldiRankTolerance;
            transientOptions.basisOrthogonalization =
                options.morBasisOrthogonalization;
            transientOptions.projection = options.morProjection;
            transientOptions.fomOnly = options.morTransientFomOnly;
            transientOptions.secondMomentEnergy =
                options.morArnoldiSecondMomentEnergy;
            transientOptions.secondMomentMaximumColumns =
                options.morArnoldiSecondMomentMaximumColumns;
            if (options.morBasis != "standard"
                && options.morBasis != "shift-invert") {
                throw std::runtime_error(
                    "--mor-basis must be standard or shift-invert.");
            }
            if (options.morBasisOrthogonalization != "euclidean"
                && options.morBasisOrthogonalization != "mass"
                && options.morBasisOrthogonalization != "k-energy") {
                throw std::runtime_error(
                    "--mor-basis-orthogonalization must be euclidean, mass, or k-energy.");
            }
            if (options.morProjection != "galerkin"
                && options.morProjection != "petrov") {
                throw std::runtime_error(
                    "--mor-projection must be galerkin or petrov.");
            }
            if (options.morBasis == "shift-invert"
                && (!(options.morBasisShift > 0.0)
                    || !std::isfinite(options.morBasisShift))) {
                throw std::runtime_error(
                    "--mor-basis shift-invert requires --mor-basis-shift > 0.");
            }
            transientOptions.basisType = options.morBasis;
            transientOptions.basisShift = options.morBasisShift;
            transientOptions.boundaryAwareBasis = options.morBoundaryAwareBasis;
            if (options.morResidualTopN <= 0
                || options.morResidualEnrichmentRounds <= 0
                || options.morResidualEnrichmentRounds > 2
                || options.morResidualPilotSteps <= 0) {
                throw std::runtime_error(
                    "Residual enrichment requires positive top-N/pilot steps and one or two rounds.");
            }
            transientOptions.residualEnrichment = options.morResidualEnrichment;
            transientOptions.residualTopN = options.morResidualTopN;
            transientOptions.residualEnrichmentRounds =
                options.morResidualEnrichmentRounds;
            transientOptions.residualPilotSteps = options.morResidualPilotSteps;
            transientOptions.massType = options.morTransientMass;
            transientOptions.timeStep = options.morTransientDt;
            transientOptions.endTime = options.morTransientEndTime;
            transientOptions.integrator = options.morTransientIntegrator;
            transientOptions.inputPath = options.morTransientInputPath;
            transientOptions.waveform = options.morTransientWaveform;
            transientOptions.outputMode = options.morTransientOutput;
            transientOptions.compareFom = options.morTransientCompareFom;
            transientOptions.compareFomSummaryOnly =
                options.morTransientCompareFomSummaryOnly;
            transientOptions.nativeReducedHistory =
                options.morNativeReducedHistory;
            if (options.morInterfaceInitialGuess != "zero"
                && options.morInterfaceInitialGuess != "previous"
                && options.morInterfaceInitialGuess != "extrapolated") {
                throw std::runtime_error(
                    "MOR interface initial guess must be zero, previous, or extrapolated.");
            }
            transientOptions.interfaceInitialGuess =
                options.morInterfaceInitialGuess;
            if (options.morInterfaceKrylov != "fgmres"
                && options.morInterfaceKrylov != "pcg"
                && options.morInterfaceKrylov != "port-core"
                && options.morInterfaceKrylov != "augmented-direct") {
                throw std::runtime_error(
                    "MOR interface solver must be fgmres, pcg, port-core, "
                    "or augmented-direct.");
            }
            transientOptions.interfaceKrylov = options.morInterfaceKrylov;
            if (options.morPortCoreCacheEnabled
                && options.morPortCoreCachePath.empty()) {
                throw std::runtime_error(
                    "MOR port/core cache requires a nonempty path.");
            }
            transientOptions.portCoreCacheEnabled =
                options.morPortCoreCacheEnabled;
            transientOptions.portCoreCachePath =
                options.morPortCoreCachePath;
            transientOptions.sourceAlignedInterfaceValidation =
                options.morSourceAlignedInterfaceValidation;
            transientOptions.fullResidualTolerance =
                options.morFullResidualTolerance;
            transientOptions.fullResidualFallback =
                options.morFullResidualFallback;
            transientOptions.initialTemperature = options.morTransientInitialTemperature;
            transientOptions.initialMode = options.morTransientInitialMode;
            transientOptions.reuseIdenticalSubdomains =
                options.morTransientReuseIdenticalSubdomains;
            if (options.morConstructionTraceMode != "global-fom"
                && options.morConstructionTraceMode != "operator-coarse") {
                throw std::runtime_error(
                    "MOR construction traces must be global-fom or operator-coarse.");
            }
            transientOptions.constructionTraceMode =
                options.morConstructionTraceMode;
            transientOptions.localInteriorArnoldiReducedSchurValidation =
                options.morTransientMethod == "local-interior-arnoldi-reduced-schur";
            transientOptions.skipReducedSchurValidation =
                options.morSkipReducedSchurValidation;
            transientOptions.interfaceExcitationRank = options.morInterfaceRank;
            transientOptions.matrixFreeInterfaceThreshold =
                options.localMorMatrixFreeInterfaceThreshold;
            transientOptions.interfaceMaxIterations = options.maxPcgIterations;
            transientOptions.interfaceRestart = options.gmresRestart;
            transientOptions.interfaceTolerance = options.pcgTolerance;
            if (options.morAdaptiveInterfaceTolerance < 0.0
                || (options.morAdaptiveInterfaceTolerance > 0.0
                    && options.morAdaptiveInterfaceTolerance
                        < options.pcgTolerance)) {
                throw std::runtime_error(
                    "Adaptive interface tolerance must be zero or no stricter "
                    "than --pcg-tolerance.");
            }
            transientOptions.adaptiveInterfaceTolerance =
                options.morAdaptiveInterfaceTolerance;
            transientOptions.coarseLinearXY = options.schurLinearXYCoarse;
            transientOptions.coarseLinearZ = options.schurLinearZCoarse;
            transientOptions.proxyEnabled = options.schurProxyEnabled;
            transientOptions.proxyDisableCoarse = options.schurProxyDisableCoarse;
            transientOptions.proxyHighConductivityThreshold =
                options.schurProxyHighKThreshold;
            transientOptions.proxyUseMaterialConnectivity =
                options.schurProxyUseMaterialConnectivity;
            transientOptions.proxyRing = options.schurProxyRing;
            transientOptions.proxyProbeColumns = options.schurProxyProbeColumns;
            transientOptions.proxyBlockSize = options.schurProxyBlockSize;
            transientOptions.proxyValidateBlockEquivalence =
                options.schurProxyValidateBlockEquivalence;
            transientOptions.localSolveThreads =
                options.schurLocalSolveThreads;
            transientOptions.localPardisoThreads =
                options.schurLocalPardisoThreads;
            transientOptions.proxyCacheEnabled =
                options.schurProxyCacheEnabled;
            transientOptions.proxyCachePath =
                options.schurProxyCachePath;
            transientOptions.interfaceOperatorCoarseRank =
                options.schurInterfaceOperatorCoarseRank;
            transientOptions.interfaceOperatorCoarseSweeps =
                options.schurInterfaceOperatorCoarseSweeps;
            transientOptions.interfaceOperatorCoarsePredictor =
                options.schurInterfaceOperatorCoarsePredictor;
            transientOptions.interfaceOperatorCoarseCachePath =
                options.schurInterfaceOperatorCoarseCachePath;
            transientOptions.portReduction =
                options.morTransientMethod == "local-port-block-arnoldi"
                && options.portBasisMethod != "full-interface";
            transientOptions.portBasisMethod = options.portBasisMethod;
            transientOptions.optimalPortRank = options.optimalPortRank;
            transientOptions.optimalPortRankMode =
                options.optimalPortRankMode;
            transientOptions.optimalPortRankFile =
                options.optimalPortRankFile;
            transientOptions.optimalPortEigenvalueTolerance =
                options.optimalPortEigenvalueTolerance;
            transientOptions.optimalPortMinimumRank =
                options.optimalPortMinimumRank;
            transientOptions.optimalPortMaximumRank =
                options.optimalPortMaximumRank;
            transientOptions.optimalPortInnerProduct =
                options.optimalPortInnerProduct;
            transientOptions.optimalPortOversamplingLayers =
                options.optimalPortOversamplingLayers;
            transientOptions.optimalPortInnerSolver =
                options.optimalPortInnerSolver;
            transientOptions.optimalPortInnerTolerance =
                options.optimalPortInnerTolerance;
            transientOptions.optimalPortInnerMaximumIterations =
                options.optimalPortInnerMaximumIterations;
            transientOptions.optimalPortInnerRefinementMaximumIterations =
                options.optimalPortInnerRefinementMaximumIterations;
            transientOptions.optimalPortInnerRefinementTolerance =
                options.optimalPortInnerRefinementTolerance;
            transientOptions.optimalPortEigenTolerance =
                options.optimalPortEigenTolerance;
            transientOptions.optimalPortEigenMaximumIterations =
                options.optimalPortEigenMaximumIterations;
            transientOptions.optimalPortAblation =
                options.optimalPortAblation;
            transientOptions.optimalPortSourceMode =
                options.optimalPortSourceMode;
            transientOptions.optimalPortTopologyAudit =
                options.optimalPortTopologyAudit;
            transientOptions.optimalPortBasisPilot =
                options.optimalPortBasisPilot;
            transientOptions.optimalPortTargetSolverComparison =
                options.optimalPortTargetSolverComparison;
            transientOptions.optimalPortWoodburyPilot =
                options.optimalPortWoodburyPilot;
            transientOptions.optimalPortRefinementValidation =
                options.optimalPortRefinementValidation;
            transientOptions.optimalPortRepresentativeInterfacePilot =
                options.optimalPortRepresentativeInterfacePilot;
            transientOptions.optimalPortMaximumInterfaceRefinementPilot =
                options.optimalPortMaximumInterfaceRefinementPilot;
            transientOptions.optimalPortAllInterfaceBasis =
                options.optimalPortAllInterfaceBasis;
            transientOptions.optimalPortTopologyAuditCsv =
                options.optimalPortTopologyAuditCsv;
            transientOptions.residualKrylovMaximumRank =
                options.residualKrylovMaximumRank;
            transientOptions.residualKrylovMaximumSweeps =
                options.residualKrylovMaximumSweeps;
            transientOptions.residualKrylovTolerance =
                options.residualKrylovTolerance;
            transientOptions.residualKrylovBlockSize =
                options.residualKrylovBlockSize;
            transientOptions.residualKrylovProbeMode =
                options.residualKrylovProbeMode;
            transientOptions.residualKrylovInnerSolver =
                options.residualKrylovInnerSolver;
            transientOptions.residualKrylovRepresentativePilot =
                options.residualKrylovRepresentativePilot;
            transientOptions.residualKrylovAllInterfaceBasis =
                options.residualKrylovAllInterfaceBasis;
            transientOptions.randomizedPortRank =
                options.randomizedPortRank;
            transientOptions.randomizedPortOversampling =
                options.randomizedPortOversampling;
            transientOptions.randomizedPortPowerIterations =
                options.randomizedPortPowerIterations;
            transientOptions.randomizedPortSeed =
                options.randomizedPortSeed;
            transientOptions.randomizedPortCompareOptimal =
                options.randomizedPortCompareOptimal;
            transientOptions.randomizedPortRepresentativePilot =
                options.randomizedPortRepresentativePilot;
            transientOptions.historyCompressionMethod =
                options.historyCompressionMethod;
            transientOptions.historyCompressionRank =
                options.historyCompressionRank;
            transientOptions.historyCompressionTolerance =
                options.historyCompressionTolerance;
            transientOptions.historyCompressionMaximumInterfacePilot =
                options.historyCompressionMaximumInterfacePilot;
            transientOptions.milestone8ProductionBasisOnly =
                options.milestone8ProductionBasisOnly;
            transientOptions.milestone8AdaptiveProduction =
                options.milestone8AdaptiveProduction;
            transientOptions.adaptivePortLocalPilot =
                options.adaptivePortLocalPilot;
            transientOptions.adaptivePortInterfaceIds =
                options.adaptivePortInterfaceIds;
            transientOptions.globalInterfaceCoarsePrototype =
                options.globalInterfaceCoarsePrototype;
            transientOptions.globalInterfaceCoarseInverseMode =
                options.globalInterfaceCoarseInverseMode;
            transientOptions.globalInterfaceCoarseExplicitReference =
                options.globalInterfaceCoarseExplicitReference;
            transientOptions.globalInterfaceCoarseRank =
                options.globalInterfaceCoarseRank;
            transientOptions.globalInterfaceCoarseCandidateDimension =
                options.globalInterfaceCoarseCandidateDimension;
            transientOptions.globalInterfaceCoarseMaximumIterations =
                options.globalInterfaceCoarseMaximumIterations;
            transientOptions
                .globalInterfaceCoarseInnerMaximumIterations =
                options.globalInterfaceCoarseInnerMaximumIterations;
            transientOptions.globalInterfaceCoarseKrylovSweeps =
                options.globalInterfaceCoarseKrylovSweeps;
            transientOptions.globalInterfaceCoarseTolerance =
                options.globalInterfaceCoarseTolerance;
            transientOptions.globalInterfaceCoarseInnerTolerance =
                options.globalInterfaceCoarseInnerTolerance;
            transientOptions.globalInterfaceCoarseInterfaceIds =
                options.globalInterfaceCoarseInterfaceIds;
            transientOptions.globalRandomizedSchur =
                options.globalRandomizedSchur;
            transientOptions.globalRandomizedRank =
                options.globalRandomizedRank;
            transientOptions.globalRandomizedSeed =
                options.globalRandomizedSeed;
            transientOptions.globalRandomizedInnerMaximumIterations =
                options.globalRandomizedInnerMaximumIterations;
            transientOptions.globalRandomizedInnerTolerance =
                options.globalRandomizedInnerTolerance;
            transientOptions.globalRandomizedComposition =
                options.globalRandomizedComposition;
            transientOptions.projectionDiagnosis =
                options.projectionDiagnosis;
            transientOptions.fluxOperatorAudit =
                options.fluxOperatorAudit;
            transientOptions.fluxAwarePort =
                options.fluxAwarePort;
            transientOptions.fluxAwareFluxType =
                options.fluxAwareFluxType;
            transientOptions.projectionInterfaceIds =
                options.projectionInterfaceIds;
            transientOptions.localPortRank = options.localPortRank;
            transientOptions.localPortEnergyTolerance =
                options.localPortEnergyTolerance;
            transientOptions.localPortRankFile = options.localPortRankFile;
            transientOptions.localPortTemperatureWeight =
                options.localPortTemperatureWeight;
            transientOptions.localPortFluxWeight = options.localPortFluxWeight;
            transientOptions.localPortResidualWeight =
                options.localPortResidualWeight;
            transientOptions.localPortEnrichmentRounds =
                options.localPortEnrichmentRounds;
            transientOptions.localPortCorrected = options.localPortCorrected;
            transientOptions.deploymentRhsCount = options.morDeploymentRhsCount;
            transientOptions.seed = options.morSeed;
            if (options.morTransientMethod == "local-block-arnoldi"
                || options.morTransientMethod == "local-port-block-arnoldi"
                || options.morTransientMethod == "local-interior-arnoldi-reduced-schur") {
                mor::transient::runLocalBlockArnoldiDynamicSchurWorkflow(
                    mesh, physics, transientOptions, outputDir);
            } else {
                mor::transient::runTransientBlockArnoldiWorkflow(
                    mesh, physics, transientOptions, outputDir);
            }
            return 0;
        }

        size_t integrationTriangleCount = 0;
        for (const InterfaceFace& face : mesh.interfaceFaces) {
            integrationTriangleCount += face.integrationTriangles.size();
        }
        int dirichletNodeCount = 0;
        for (const Node& node : mesh.nodes) {
            if (node.dirichlet) {
                ++dirichletNodeCount;
            }
        }

        std::cout << "COMSOL mphtxt mesh + P2 tetra FEM + "
                  << physics.interfaceScheme << " domain decomposition\n";
        std::cout << "Project root: " << projectRoot.string() << "\n";
        std::cout << "Case config: " << configPath.string() << "\n";
        std::cout << "Case name: " << physics.name << "\n";
        std::cout << "Model name: " << modelName << "\n";
        std::cout << "Output dir: " << outputDir.string() << "\n";
        std::cout << "Analysis mode: " << modeName
                  << (steady ? " (K*T = Q)" : " (M/dt + K time marching)") << "\n";
        std::cout << "Solver method: " << physics.solverMethod << "\n";
        if (physics.schwarz.enabled || options.runSchwarz) {
            std::cout << "Schwarz options: type=" << physics.schwarz.type
                      << ", standalone_mode=" << physics.schwarz.standaloneMode
                      << ", transmission=" << physics.schwarz.transmission
                      << ", transmission_orientation=" << physics.schwarz.transmissionOrientation
                      << ", robin_alpha_factor=" << physics.schwarz.robinAlphaFactor
                      << ", flux_eval=" << physics.schwarz.fluxEval
                      << ", overlap_layers=" << physics.schwarz.overlapLayers
                      << ", overlap_mode=" << physics.schwarz.overlapMode
                      << ", partition_mode=" << physics.schwarz.partitionMode
                      << ", max_iters=" << physics.schwarz.maxIters
                      << ", tol_rel_update=" << physics.schwarz.tolRelUpdate
                      << ", tol_rel_residual=" << physics.schwarz.tolRelResidual
                      << ", relaxation=" << physics.schwarz.relaxation
                      << ", write_interface_flux="
                      << (physics.schwarz.writeInterfaceFlux ? "true" : "false")
                      << ", validate_against_monolithic="
                      << (physics.schwarz.validateAgainstMonolithic ? "true" : "false")
                      << "\n";
        }
    std::cout << "Interface consistency: " << (options.disableInterfaceConsistency ? "disabled" : "enabled")
                  << ", interface penalty: " << (options.disableInterfacePenalty ? "disabled" : "enabled") << "\n";
        std::cout << "Robin/convection LHS: " << (options.disableConvectionLhs ? "disabled" : "enabled")
                  << ", diagnostics_only=" << (options.diagnosticsOnly ? "yes" : "no") << "\n";
        std::cout << "Dirichlet enforcement: " << physics.dirichletMethod;
        if (physics.dirichletMethod == "nitsche") {
            std::cout << ", nitsche_penalty_factor=" << physics.nitschePenaltyFactor;
        }
        std::cout << "\n";
        std::cout << "Fast run: " << (options.fastRun ? "enabled" : "disabled")
                  << ", direct_mode=" << options.directMode << "\n";
        std::cout << "IC options: ordering=" << options.icOrdering
                  << ", scaling=" << (options.icScaling ? "enabled" : "disabled")
                  << ", initial_shift=" << options.icShift << "\n";
        std::cout << "Solver parallel workers: " << solverParallelWorkers()
                  << " (BJ preconditioner also solves subdomain blocks concurrently)\n";
        std::cout << "Domains: " << physics.domains.size() << "\n";
        for (size_t i = 0; i < physics.domains.size(); ++i) {
            const Material& mat = physics.domains[i].material;
            std::cout << "  domain " << i << ": " << physics.domains[i].meshPath.string()
                      << ", default_material=" << mat.name
                      << ", k=(" << mat.conductivityX << ","
                      << mat.conductivityY << "," << mat.conductivityZ << ")"
                      << " W/(m K), rho=" << mat.density
                      << " kg/m^3, cp=" << mat.heatCapacity << " J/(kg K)"
                      << ", offset_m=(" << physics.domains[i].translationMeters.x
                      << "," << physics.domains[i].translationMeters.y
                      << "," << physics.domains[i].translationMeters.z << ")\n";
            std::cout << "    COMSOL domain material mappings: "
                      << physics.domains[i].materialsByDomainEntity.size() << "\n";
            for (const auto& mapped : physics.domains[i].materialsByDomainEntity) {
                const Material& domainMat = mapped.second;
                if (domainMat.name != "Die") {
                    continue;
                }
                std::cout << "    COMSOL domain_entity=" << mapped.first
                          << ", material=" << domainMat.name
                          << ", k=(" << domainMat.conductivityX << ","
                          << domainMat.conductivityY << "," << domainMat.conductivityZ << ")"
                          << " W/(m K), rho=" << domainMat.density
                          << " kg/m^3, cp=" << domainMat.heatCapacity << " J/(kg K)\n";
            }
        }
        std::cout << "COMSOL coordinates are multiplied by coordinate_scale="
                  << physics.coordinateScale << ".\n";
        std::cout << "Nodes(P2): " << mesh.nodes.size()
                  << ", tets: " << mesh.tets.size()
                  << ", boundary triangles: " << mesh.boundaryFaces.size() << "\n";
        std::cout << physics.interfaceScheme << " interface face pairs: " << mesh.interfaceFaces.size()
                  << ", interface integration triangles: " << integrationTriangleCount << "\n";
        std::cout << physics.interfaceScheme << " penalty mode: " << physics.penaltyMode
                  << ", penalty_factor=" << physics.penaltyFactor
                  << ", penalty_scaling=" << physics.penaltyScaling << "\n";
        for (size_t i = 0; i < mesh.interfaceSummaries.size(); ++i) {
            const InterfaceBuildSummary& s = mesh.interfaceSummaries[i];
            const double overlapRatio = std::min(s.leftArea, s.rightArea) > 0.0
                ? s.matchedOverlapArea / std::min(s.leftArea, s.rightArea)
                : 0.0;
            const double normalAvg = s.facePairCount > 0
                ? s.normalDotSum / static_cast<double>(s.facePairCount)
                : 0.0;
            std::cout << "  interface " << i
                      << ": left_area=" << s.leftArea
                      << " m^2, right_area=" << s.rightArea
                      << " m^2, matched_overlap=" << s.matchedOverlapArea
                      << " m^2, overlap_ratio=" << overlapRatio
                      << ", face_pairs=" << s.facePairCount
                      << ", integration_triangles=" << s.integrationTriangleCount
                      << ", normal_dot[min/avg/max]="
                      << s.normalDotMin << "/" << normalAvg << "/" << s.normalDotMax << "\n";
            if (overlapRatio < 0.99) {
                throw std::runtime_error("Interface overlap ratio is below 0.99; stopping before solve.");
            }
        }
        std::cout << "Dirichlet conditions: " << physics.dirichletConditions.size()
                  << "; marked nodes=" << dirichletNodeCount << "\n";
        for (const BoundaryCondition& condition : physics.dirichletConditions) {
            std::cout << "  subdomain=" << (condition.subdomain < 0 ? std::string("*") : std::to_string(condition.subdomain))
                      << ", boundary_entity=" << condition.boundaryEntity
                      << ", T=" << condition.temperature << " K\n";
        }
        std::cout << "Convection/Robin conditions: " << physics.convectionConditions.size() << "\n";
        double convectionMinH = std::numeric_limits<double>::infinity();
        double convectionMaxH = -std::numeric_limits<double>::infinity();
        double convectionMinAmbient = std::numeric_limits<double>::infinity();
        double convectionMaxAmbient = -std::numeric_limits<double>::infinity();
        for (const ConvectionCondition& condition : physics.convectionConditions) {
            convectionMinH = std::min(convectionMinH, condition.coefficient);
            convectionMaxH = std::max(convectionMaxH, condition.coefficient);
            convectionMinAmbient = std::min(convectionMinAmbient, condition.ambientTemperature);
            convectionMaxAmbient = std::max(convectionMaxAmbient, condition.ambientTemperature);
        }
        if (!physics.convectionConditions.empty()) {
            std::cout << "  h range=" << convectionMinH << ".." << convectionMaxH
                      << " W/(m^2 K), Tamb range=" << convectionMinAmbient
                      << ".." << convectionMaxAmbient << " K\n";
        }
        const size_t convectionPrintLimit = std::min<size_t>(physics.convectionConditions.size(), 10);
        for (size_t i = 0; i < convectionPrintLimit; ++i) {
            const ConvectionCondition& condition = physics.convectionConditions[i];
            std::cout << "  subdomain=" << (condition.subdomain < 0 ? std::string("*") : std::to_string(condition.subdomain))
                      << ", boundary_entity=" << condition.boundaryEntity
                      << ", h=" << condition.coefficient
                      << " W/(m^2 K), Tamb=" << condition.ambientTemperature << " K\n";
        }
        if (physics.convectionConditions.size() > convectionPrintLimit) {
            std::cout << "  ... " << (physics.convectionConditions.size() - convectionPrintLimit)
                      << " additional convection/Robin conditions omitted from console log\n";
        }
        std::cout << "Heat sources: " << physics.heatSources.size() << "\n";
        double configuredHeatTotalW = 0.0;
        for (const HeatSource& heatSource : physics.heatSources) {
            configuredHeatTotalW += heatSource.heatRateW;
            std::cout << "  subdomain=" << (heatSource.subdomain < 0 ? std::string("*") : std::to_string(heatSource.subdomain))
                      << ", domain_entity=" << heatSource.domainEntity
                      << ", total_Q=" << heatSource.heatRateW << " W\n";
        }
        std::cout << "Configured heat source total: " << configuredHeatTotalW << " W\n";
        if (!(configuredHeatTotalW > 0.0)) {
            throw std::runtime_error("Configured heat source total is zero; stopping before solve.");
        }
        int preAssemblyNonPositiveConductivityTets = 0;
        for (const Tet& tet : mesh.tets) {
            const Material& material = materialForTet(physics, tet);
            if (material.conductivityX <= 0.0 || material.conductivityY <= 0.0 || material.conductivityZ <= 0.0
                || !std::isfinite(material.conductivityX)
                || !std::isfinite(material.conductivityY)
                || !std::isfinite(material.conductivityZ)) {
                ++preAssemblyNonPositiveConductivityTets;
            }
        }
        std::cout << "Pre-assembly material diagnostic: k<=0/nonfinite tets="
                  << preAssemblyNonPositiveConductivityTets << "\n";
        if (preAssemblyNonPositiveConductivityTets > 0) {
            throw std::runtime_error("Detected k<=0/nonfinite material tets before assembly; stopping before solve.");
        }
        if (options.spdPenaltySweep) {
            if (!steady) {
                throw std::runtime_error("--spd-penalty-sweep is implemented for steady K*T=Q diagnostics.");
            }
            return runSpdPenaltySweep(mesh, physics, options, outputDir);
        }
        if (options.sipgSpdDiagnostics) {
            if (!steady) {
                throw std::runtime_error("--sipg-spd-diagnostics is implemented for steady K*T=Q diagnostics.");
            }
            return runSipgSpdDiagnostics(mesh, physics, options, outputDir, modelPrefix);
        }
        if (!steady) {
            std::cout << "Time: range(" << physics.startTime << "," << physics.timeStep
                      << "," << physics.timeStep * physics.timeSteps << ") s\n";
        }
        timing.preprocessingSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - programStart).count();

        const auto assemblyStart = std::chrono::steady_clock::now();
        SparseMatrix mass(n);
        SparseMatrix stiffness(n);
        std::vector<double> source(static_cast<size_t>(n), 0.0);
        std::vector<MatrixStageDiagnostic> matrixStageDiagnostics;
        std::vector<HeatSourceAssemblyDiagnostic> heatSourceDiagnostics;
        if (collectDetailedDiagnostics) {
            heatSourceDiagnostics = collectHeatSourceDiagnostics(mesh, physics);
        }
        const auto volumeAssemblyStart = std::chrono::steady_clock::now();
        assembleVolume(mesh, physics, steady ? nullptr : &mass, stiffness, source);
        if (physics.thermalSourceScale != 1.0) {
            for (double& value : source) {
                value *= physics.thermalSourceScale;
            }
            std::cout << "Thermal source scale applied: " << physics.thermalSourceScale << "\n";
        }
        if (collectDetailedDiagnostics) {
            assemblyDiagnostics.volumeDiag = matrixDiagonalVector(stiffness);
            matrixStageDiagnostics.push_back(diagnoseMatrixStage("bulk_stiffness", stiffness));
        }
        const std::vector<double> heatOnlySource = source;
        const int convectionBoundaryFaceCount = assembleConvectionBoundaries(mesh,
                                                                             physics,
                                                                             stiffness,
                                                                             source,
                                                                             options.disableConvectionLhs,
                                                                             collectDetailedDiagnostics ? &assemblyDiagnostics : nullptr);
        if (collectDetailedDiagnostics) {
            matrixStageDiagnostics.push_back(diagnoseMatrixStage("bulk_plus_robin_lhs", stiffness));
        }
        const std::vector<double> sourceBeforeDirichlet = source;
        SparseMatrix physicalTransmissionStiffness = stiffness;
        timing.volumeAssemblySeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - volumeAssemblyStart).count();

        const auto interfaceAssemblyStart = std::chrono::steady_clock::now();
        assembleSipgInterface(mesh,
                              physics,
                              stiffness,
                              options.disableInterfaceConsistency,
                              true,
                              collectDetailedDiagnostics ? &assemblyDiagnostics : nullptr);
        if (collectDetailedDiagnostics) {
            matrixStageDiagnostics.push_back(diagnoseMatrixStage("bulk_robin_plus_sipg_consistency", stiffness));
        }
        assembleSipgInterface(mesh,
                              physics,
                              stiffness,
                              true,
                              options.disableInterfacePenalty,
                              collectDetailedDiagnostics ? &assemblyDiagnostics : nullptr);
        if (collectDetailedDiagnostics) {
            matrixStageDiagnostics.push_back(diagnoseMatrixStage("bulk_robin_plus_full_sipg", stiffness));
        }
        if (options.nodeTieInterface) {
            const int tiedPairs = assembleMatchedNodeTieInterface(mesh, stiffness, options.nodeTiePenalty, outputDir);
            std::cout << "Matched-node tie interface constraints: pairs=" << tiedPairs
                      << ", penalty=" << options.nodeTiePenalty << "\n";
            if (collectDetailedDiagnostics) {
                matrixStageDiagnostics.push_back(diagnoseMatrixStage("bulk_robin_sipg_plus_node_tie", stiffness));
            }
        }
        timing.interfaceAssemblySeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - interfaceAssemblyStart).count();

        SparseMatrix system(n);
        SparseMatrix physicalTransmissionSystem(n);
        const auto systemBuildStart = std::chrono::steady_clock::now();
        if (steady) {
            system = std::move(stiffness);
            physicalTransmissionSystem = std::move(physicalTransmissionStiffness);
        }
        if (!steady) {
            buildSystemMatrix(mass, stiffness, physics.timeStep, system);
            buildSystemMatrix(mass, physicalTransmissionStiffness, physics.timeStep, physicalTransmissionSystem);
            // The independent transient-local-rom-test needs the original K
            // blocks to form Cr/dt+Kr.  Existing transient workflows retain
            // their historical memory-saving move semantics.
            if (!options.transientLocalRomTest) {
                stiffness = SparseMatrix(n);
            }
            physicalTransmissionStiffness = SparseMatrix(n);
        }
        if (collectDetailedDiagnostics) {
            assemblyDiagnostics.preDirichletDiag = matrixDiagonalVector(system);
        }
        timing.systemBuildSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - systemBuildStart).count();
        SparseMatrix preDirichletSystem;
        if (options.spectralDiagnostics || options.runDeflatedRasIlut || options.runInterfaceDeflatedRasIlut) {
            preDirichletSystem = system;
            preDirichletSystem.finalizeCsr();
        }

        const auto dirichletAssemblyStart = std::chrono::steady_clock::now();
        std::vector<double> fixedAdjust;
        std::vector<double> sourceAfterDirichlet = source;
        if (physics.dirichletMethod == "nitsche") {
            const int nitscheFaces = assembleNitscheDirichletBoundaries(
                mesh,
                physics,
                system,
                source,
                collectDetailedDiagnostics ? &assemblyDiagnostics : nullptr);
            std::vector<double> robinBaseSource = source;
            assembleNitscheDirichletBoundaries(mesh, physics, physicalTransmissionSystem, robinBaseSource, nullptr);
            sourceAfterDirichlet = source;
            std::cout << "Nitsche Dirichlet boundary faces assembled: " << nitscheFaces
                      << ", penalty_factor=" << physics.nitschePenaltyFactor << "\n";
        } else {
            fixedAdjust = makeDirichletAdjustedSystem(mesh, system);
            makeDirichletAdjustedSystem(mesh, physicalTransmissionSystem);
            sourceAfterDirichlet = sourceBeforeDirichlet;
            applyDirichletRhs(mesh, fixedAdjust, sourceAfterDirichlet);
        }
        if (collectDetailedDiagnostics) {
            matrixStageDiagnostics.push_back(diagnoseMatrixStage(
                physics.dirichletMethod == "nitsche" ? "after_nitsche_dirichlet" : "after_symmetric_dirichlet",
                system));
        }
        if (steady && physics.forceNontrivialRhs) {
            const double trialRel = relativeResidualNorm(system,
                                                         initialTemperatureVector(mesh, physics),
                                                         sourceAfterDirichlet);
            if (!(trialRel > 1.0e-3)) {
                int targetDof = -1;
                for (int i = 0; i < static_cast<int>(mesh.nodes.size()); ++i) {
                    if (!mesh.nodes[static_cast<size_t>(i)].dirichlet) {
                        targetDof = i;
                        break;
                    }
                }
                if (targetDof >= 0) {
                    const double rhsNormBefore = vectorNorm(sourceAfterDirichlet);
                    const double delta = 1.0e-2 * std::max(1.0, rhsNormBefore);
                    source[static_cast<size_t>(targetDof)] += delta;
                    sourceAfterDirichlet = source;
                    if (!fixedAdjust.empty()) {
                        applyDirichletRhs(mesh, fixedAdjust, sourceAfterDirichlet);
                    }
                    std::cout << "force-nontrivial-rhs added point load " << delta
                              << " to dof " << targetDof
                              << " because trial initial relative residual was " << trialRel << "\n";
                } else {
                    std::cout << "WARNING: force-nontrivial-rhs found no non-Dirichlet DOF to perturb.\n";
                }
            }
        }
        timing.dirichletAssemblySeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - dirichletAssemblyStart).count();

        const auto csrFinalizeStart = std::chrono::steady_clock::now();
        if (!steady) {
            mass.finalizeCsr();
        }
        if (!steady && options.transientLocalRomTest && !stiffness.csrReady) {
            stiffness.finalizeCsr();
        }
        system.finalizeCsr();
        physicalTransmissionSystem.finalizeCsr();
        if (collectDetailedDiagnostics) {
            assemblyDiagnostics.finalDiag = matrixDiagonalVector(system);
            for (size_t i = 0; i < assemblyDiagnostics.finalDiag.size(); ++i) {
                assemblyDiagnostics.dirichletDiag[i] =
                    assemblyDiagnostics.finalDiag[i] - assemblyDiagnostics.preDirichletDiag[i];
            }
        }
        timing.csrFinalizeSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - csrFinalizeStart).count();
        timing.assemblySeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - assemblyStart).count();

        std::vector<ComponentDiagnostics> componentDiagnostics;
        MatrixDiagnostics matrixDiagnostics;
        if (collectDetailedDiagnostics) {
            componentDiagnostics = collectComponentDiagnostics(mesh, physics);
            matrixDiagnostics =
                diagnoseMatrixAndPhysics(mesh, physics, system, convectionBoundaryFaceCount, componentDiagnostics);
            printMatrixDiagnostics(matrixDiagnostics);
            printDiagonalContributionStats(assemblyDiagnostics);
            printInterfacePenaltyStats(assemblyDiagnostics.interfacePenaltyStats);
            writeMatrixDiagnostics(matrixDiagnostics, outputDir / "matrix_diagnostics.csv");
            writeMatrixStageDiagnostics(matrixStageDiagnostics, outputDir / "matrix_stage_diagnostics.csv");
            writeComponentDiagnostics(componentDiagnostics, outputDir / "component_diagnostics.csv");
            writeInterfaceBuildSummary(mesh, outputDir / "interface_build_summary.csv");
            writeDiagonalContributionStats(assemblyDiagnostics, outputDir / "diagonal_contribution_stats.csv");
            writeNegativeDiagonalDofs(mesh, assemblyDiagnostics, outputDir / "negative_diagonal_dofs.csv");
            writeInterfacePenaltyStats(assemblyDiagnostics.interfacePenaltyStats, outputDir / "interface_penalty_stats.csv");
            writeHeatSourceDiagnostics(heatSourceDiagnostics, outputDir / "heat_source_diagnostics.csv");
            writeHeatSourceDomainDiagnostics(heatSourceDiagnostics, outputDir / "heat_source_domain_diagnostics.csv");
            writePhysicalSummary(mesh,
                                 physics,
                                 heatOnlySource,
                                 sourceBeforeDirichlet,
                                 sourceAfterDirichlet,
                                 convectionBoundaryFaceCount,
                                 outputDir / "physical_summary.csv");
            const FinalSpdDiagnostic finalSpdDiagnostic =
                computeFinalSpdDiagnostic(matrixDiagnostics, system);
            const std::vector<SweepBlockDiagnostic> finalBlockSpdDiagnostics =
                collectSweepBlockDiagnostics(mesh, system, physics, &options);
            writeFinalSpdDiagnostic(finalSpdDiagnostic, outputDir / "final_matrix_spd_diagnostics.csv");
            writeSweepBlockDiagnostics(finalBlockSpdDiagnostics, outputDir / "final_bj_block_spd_diagnostics.csv");
            const bool confirmedSpd = finalSpdDiagnostic.spdConclusion == "numerically_spd";
            const bool allBjBlocksSpd = std::all_of(finalBlockSpdDiagnostics.begin(),
                                                    finalBlockSpdDiagnostics.end(),
                                                    [](const SweepBlockDiagnostic& block) {
                                                        return block.ldltStatus == "success"
                                                            && block.ldltNegativePivots == 0
                                                            && block.ldltZeroTinyPivots == 0
                                                            && block.pardisoSpdStatus == "success"
                                                            && std::isfinite(block.lambdaMinEst)
                                                            && block.lambdaMinEst > 0.0;
                                                    });
            std::string spdCase = "C";
            if (!confirmedSpd) {
                spdCase = "A";
            } else if (!allBjBlocksSpd) {
                spdCase = "B";
            }
            std::cout << "Final A_final SPD diagnostics:\n"
                      << "  NaN/Inf=" << (finalSpdDiagnostic.hasNonFinite ? "yes" : "no")
                      << ", symmetry_error=" << finalSpdDiagnostic.symmetryError
                      << ", diagonal_min/max=" << finalSpdDiagnostic.minDiag
                      << "/" << finalSpdDiagnostic.maxDiag
                      << ", zero_rows=" << finalSpdDiagnostic.zeroRowCount << "\n"
                      << "  LDLT inertia positive/negative/zero="
                      << finalSpdDiagnostic.ldltPositivePivots << "/"
                      << finalSpdDiagnostic.ldltNegativePivots << "/"
                      << finalSpdDiagnostic.ldltZeroTinyPivots
                      << ", lambda_min_est=" << finalSpdDiagnostic.lambdaMinEst
                      << ", PARDISO_SPD=" << finalSpdDiagnostic.pardisoSpdStatus
                      << ", PARDISO_general=" << finalSpdDiagnostic.pardisoGeneralStatus << "\n"
                      << "  SPD conclusion=" << finalSpdDiagnostic.spdConclusion
                      << ", classified_case=" << spdCase << "\n";
            if (finalSpdDiagnostic.pardisoGeneralStatus == "success"
                && finalSpdDiagnostic.pardisoSpdStatus != "success") {
                std::cout << "  PARDISO general succeeded but SPD mode failed: A_final may be non-SPD or numerically unsuitable for SPD factorization.\n";
            }
        } else {
            std::cout << "Fast run: skipped detailed matrix SPD checks, BJ block SPD checks, and full diagnostic CSV generation.\n";
        }
        SpectralSummary spectralSummary;
        if (options.spectralDiagnostics || options.runDeflatedRasIlut || options.runInterfaceDeflatedRasIlut) {
            const int requestedSpectralModes = std::max(options.spectralModeCount, options.deflationModes);
            spectralSummary = writeEigenDiagnostics(mesh, physics, assemblyDiagnostics, system, outputDir, requestedSpectralModes);
        }
        if (collectDetailedDiagnostics) {
            const std::vector<InterfaceFacePairDiagnosticRow> interfacePairDiagnostics =
                collectInterfaceFacePairDiagnostics(mesh, physics);
            writeInterfaceFacePairDiagnostics(interfacePairDiagnostics,
                                              outputDir / "rram_interface_face_pair_diagnostics.csv");
            writeNonmatchingInterfaceProjectionDiagnostics(interfacePairDiagnostics,
                                                           outputDir / "rram_nonmatching_interface_projection_diagnostics.csv");
            writeMissingInterfaceDiagnostics(mesh,
                                             physics,
                                             interfacePairDiagnostics,
                                             outputDir / "rram_missing_interface_diagnostics.csv");
        } else {
            std::cout << "Fast run: skipped interface face-pair and missing-interface diagnostics.\n";
        }
        if (options.spectralDiagnostics || options.runDeflatedRasIlut || options.runInterfaceDeflatedRasIlut) {
            writeThermalConnectivityDiagnostics(mesh,
                                                physics,
                                                assemblyDiagnostics,
                                                preDirichletSystem,
                                                spectralSummary,
                                                outputDir);
        }
        if (options.diagnosticsOnly) {
            std::vector<InterfaceEnergyDiagnostic> interfaceEnergyDiagnostics;
            InterfaceScaleDiagnostic interfaceScaleDiagnostics;
            std::vector<InterfaceLocalDiagnostic> localInterfaceDiagnostics;
            std::vector<InterfacePointDiagnostic> barycentricInterfaceDiagnostics;
            collectInterfaceDiagnostics(mesh,
                                        physics,
                                        interfaceEnergyDiagnostics,
                                        interfaceScaleDiagnostics,
                                        localInterfaceDiagnostics,
                                        barycentricInterfaceDiagnostics);
            printInterfaceContinuousFieldDiagnostics(interfaceEnergyDiagnostics);
            writeInterfaceContinuousFieldDiagnostics(interfaceEnergyDiagnostics,
                                                    outputDir / "interface_continuous_field_diagnostics.csv");
            writeLocalInterfacePairDiagnostics(localInterfaceDiagnostics,
                                               outputDir / "local_interface_pair_diagnostics.csv");
            writeBarycentricDiagnostics(barycentricInterfaceDiagnostics,
                                        outputDir / "barycentric_diagnostics.csv");
            SchwarzOptions rramDiagnosticOptions = physics.schwarz;
            rramDiagnosticOptions.standaloneMode = "algebraic";
            rramDiagnosticOptions.transmission = "none";
            SchwarzBlockSolver rramDiagnosticSolver(mesh,
                                                    system,
                                                    rramDiagnosticOptions,
                                                    &physics,
                                                    nullptr,
                                                    false);
            writeSchwarzBlockMatrixCheck(rramDiagnosticSolver,
                                         outputDir / "schwarz_block_matrix_check.csv");
            writeRramDiagnosticsBundle(mesh,
                                       physics,
                                       rramDiagnosticSolver,
                                       {},
                                       outputDir);
            std::vector<ValidationSolverResult> noSolverResults;
            writeRramSolverDiagnosticsSummary(noSolverResults,
                                              mesh,
                                              physics,
                                              rramDiagnosticSolver,
                                              outputDir / "rram_solver_diagnostics_summary.csv");
            writeProgramTiming(timing, outputDir / "program_timing.csv");
            std::cout << "Diagnostics-only mode: wrote matrix diagnostics and skipped solvers.\n";
            return 0;
        }
        if (collectDetailedDiagnostics
            && (matrixDiagnostics.hasNonFinite
            || matrixDiagnostics.zeroRows > 0
            || matrixDiagnostics.nonPositiveConductivityTets > 0
            || matrixDiagnostics.heatSourceTetCount == 0)) {
            throw std::runtime_error("Matrix/physics diagnostics failed stop criteria; see output diagnostics CSV files.");
        }
        if (options.bjIcValidation) {
            return runBjIcValidationPipeline(mesh,
                                             physics,
                                             options,
                                             system,
                                             sourceAfterDirichlet,
                                             outputDir);
        }
        if (options.transientLocalRomTest) {
            return runTransientLocalInteriorRomTest(mesh,
                                                    physics,
                                                    mass,
                                                    stiffness,
                                                    system,
                                                    sourceAfterDirichlet,
                                                    fixedAdjust,
                                                    options,
                                                    outputDir);
        }
        if (options.runBjIlut || options.runRasIlut || options.runTwoLevelRasIlut
            || options.runDeflatedRasIlut || options.runInterfaceDeflatedRasIlut
            || options.runBjPardisoGeneral || options.runSchwarzPrecondFgmres
            || options.runRasIc || options.runBjIcCoarse || options.runBjIlutCoarse
            || options.runSchwarzPrecondFgmresTwoLevel || options.runDdmSchur
            || options.runReducedSchur
            || (options.runLocalRom && options.localMorMode == "corrected")) {
            std::cout << "FGMRES-family solvers are selected for this run.\n";
        } else {
            std::cout << "Selected solvers avoid FGMRES-family methods.\n";
        }

        struct SolverRunResult {
            SolverStatistics stats;
            std::vector<double> temperature;
            std::vector<double> finalRhs;
        };
        std::vector<SolverRunResult> results;
        std::vector<IcFactorDiagnostics> icDiagnostics;
        std::vector<KrylovIterationRow> allKrylovRows;

        if (options.runDiagonal) {
            SolverStatistics stats;
            stats.preconditionerBytes = static_cast<size_t>(n) * sizeof(double);
            std::vector<double> temperature = runAnalysisSolver(
                steady,
                "Diagonal-PCG",
                mesh, physics, mass, system, source, fixedAdjust, stats,
                [&](const SparseMatrix& a, const std::vector<double>& b, std::vector<double> x, int& iterations) {
                    return preconditionedConjugateGradient(a,
                                                           b,
                                                           std::move(x),
                                                           iterations,
                                                           &stats,
                                                           options.maxPcgIterations,
                                                           options.pcgTolerance);
                });
            results.push_back({stats, std::move(temperature)});
        }

        if (options.runP2P1AuxPcg) {
            SolverStatistics stats;
            stats.name = "P2-P1-Auxiliary-PCG";
            stats.solverMethod = "p2_p1_auxiliary_pcg";
            stats.preconditioner = "element_additive_schwarz_plus_inherited_p1";
            stats.localSolver = "element_cholesky_reused+p1_pardiso_spd_reused";
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = peakWorkingSetBytes();
            P2P1AuxiliaryReport auxiliary;
            std::vector<double> temperature;
            std::vector<double> finalRhs;
            try {
                const bool largeCaseGates = n > 100000;
                auto preconditioner = std::make_unique<P2P1AuxiliaryPreconditioner>(
                    mesh, system, auxiliary, largeCaseGates);
                stats.setupSeconds = auxiliary.setupSeconds;
                stats.coarseEnabled = true;
                stats.coarseDim = auxiliary.p1Dofs;
                stats.coarseMatrixNnz = auxiliary.coarseNnz;
                stats.coarseSetupSeconds = auxiliary.coarseAssemblySeconds
                    + auxiliary.coarseSymbolicSeconds + auxiliary.coarseNumericalSeconds;
                stats.preconditionerBytes = preconditioner->memoryBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes,
                                                     peakWorkingSetBytes());
                const int maximumIterations = largeCaseGates
                    ? std::min(60, options.maxPcgIterations)
                    : options.maxPcgIterations;
                temperature = runAnalysisSolver(
                    steady,
                    "P2-P1-Auxiliary-PCG",
                    mesh, physics, mass, system, source, fixedAdjust, stats,
                    [&](const SparseMatrix& a,
                        const std::vector<double>& b,
                        std::vector<double> x,
                        int& iterations) {
                        return blockJacobiPcg(
                            a, b, std::move(x), *preconditioner, iterations, &stats,
                            "P2-P1-Auxiliary-PCG", maximumIterations,
                            options.pcgTolerance, {}, false);
                    },
                    &finalRhs);
                auxiliary.pcgIterations = stats.totalIterations;
                auxiliary.solveSeconds = stats.solveSeconds;
                auxiliary.totalSeconds = stats.setupSeconds + stats.solveSeconds;
                auxiliary.trueResidual = stats.trueRelativeResidual;
                auxiliary.peakMemoryBytes = std::max(auxiliary.peakMemoryBytes,
                                                     stats.peakWorkingSetBytes);
                auxiliary.status = stats.status;
                if (largeCaseGates && stats.totalIterations >= 60
                    && stats.trueRelativeResidual > options.pcgTolerance) {
                    auxiliary.status = "stopped";
                    auxiliary.stopReason = "rram5_pcg_iteration_gate_exceeded";
                    stats.status = "failed";
                    stats.failureReason = auxiliary.stopReason;
                }
            } catch (const std::exception& error) {
                stats.status = "failed";
                stats.failureReason = error.what();
                stats.setupSeconds = auxiliary.setupSeconds;
                auxiliary.totalSeconds = auxiliary.setupSeconds;
                auxiliary.peakMemoryBytes = std::max(auxiliary.peakMemoryBytes,
                                                     peakWorkingSetBytes());
                if (auxiliary.stopReason.empty()) {
                    auxiliary.stopReason = error.what();
                }
                if (auxiliary.status == "not_run") {
                    auxiliary.status = "failed";
                }
                temperature = initialTemperatureVector(mesh, physics);
            }
            writeP2P1AuxiliaryReport(auxiliary,
                outputDir / "p2_p1_auxiliary_summary.csv");
            results.push_back({stats, std::move(temperature), std::move(finalRhs)});
        }

        if (options.runSchurDirectExact) {
            SolverStatistics stats;
            stats.name = "Exact-Schur-Direct";
            stats.solverMethod = "schur_direct_exact";
            stats.preconditioner = "none_explicit_exact_schur";
            stats.localSolver = "PARDISO_local_and_interface_direct";
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = peakWorkingSetBytes();
            std::vector<double> temperature;
            std::vector<double> finalRhs;
            try {
                ddm_schur::ExactSchurDirectOptions exactOptions;
                exactOptions.patternOnly = options.schurDirectPatternOnly;
                exactOptions.dryRun = options.schurDirectDryRun;
                exactOptions.batchSize = options.schurDirectBatchSize;
                exactOptions.memoryFraction = options.schurDirectMemoryFraction;
                exactOptions.factorMemoryFraction = options.schurDirectFactorMemoryFraction;
                exactOptions.verifyOperator = options.schurDirectVerifyOperator;
                exactOptions.randomChecks = options.schurDirectRandomChecks;
                exactOptions.outputDirectory = outputDir;
                ddm_schur::ExactSchurDirectSolver solver(mesh, system, physics, exactOptions);
                solver.writeReport(outputDir / "schur_direct_exact_summary.csv");
                if (options.schurDirectPatternOnly || options.schurDirectDryRun) {
                    std::cout << "Exact Schur feasibility mode stopped before numerical solve.\n";
                    return 0;
                }
                const auto& exactReport = solver.report();
                stats.setupSeconds = exactReport.patternSeconds + exactReport.assemblySeconds
                    + exactReport.localPhase11Seconds + exactReport.localPhase22Seconds
                    + exactReport.factorPhase11Seconds + exactReport.factorPhase22Seconds;
                stats.preconditionerBytes = exactReport.rawCsrBytes + exactReport.factorBytes
                    + exactReport.localFactorBytes;
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes,
                                                     exactReport.peakWorkingSetBytes);
                stats.totalIterations = 1;
                stats.maxIterations = 1;
                if (!solver.canSolve()) {
                    stats.status = "failed";
                    stats.failureReason = exactReport.abortReason.empty()
                        ? "exact Schur feasibility gate stopped before solve"
                        : exactReport.abortReason;
                    temperature = initialTemperatureVector(mesh, physics);
                } else {
                    temperature = runAnalysisSolver(
                        steady,
                        stats.name,
                        mesh, physics, mass, system, source, fixedAdjust, stats,
                        [&](const SparseMatrix&, const std::vector<double>& b,
                            std::vector<double>, int& iterations) {
                            solver.solve(b, temperature);
                            iterations = 1;
                            stats.status = "success";
                            return temperature;
                        },
                        &finalRhs);
                    solver.recordTrueResidual(stats.trueRelativeResidual);
                    solver.writeReport(outputDir / "schur_direct_exact_summary.csv");
                }
                stats.workingSetAfterBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes,
                                                     peakWorkingSetBytes());
            } catch (const std::exception& error) {
                stats.status = "failed";
                stats.failureReason = error.what();
                temperature = initialTemperatureVector(mesh, physics);
                stats.workingSetAfterBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes,
                                                     peakWorkingSetBytes());
            }
            results.push_back({stats, std::move(temperature), std::move(finalRhs)});
        }

        if (options.runBj) {
            SolverStatistics stats;
            stats.name = "BJ-PARDISO-PCG";
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = peakWorkingSetBytes();
            std::vector<double> temperature;
            try {
                const auto setupStart = std::chrono::steady_clock::now();
                BlockJacobiPreconditioner preconditioner(mesh, system);
                const auto setupEnd = std::chrono::steady_clock::now();
                stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                stats.preconditionerBytes = preconditioner.memoryBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                temperature = runAnalysisSolver(
                    steady,
                    "BJ-PARDISO-PCG",
                    mesh, physics, mass, system, source, fixedAdjust, stats,
                    [&](const SparseMatrix& a, const std::vector<double>& b, std::vector<double> x, int& iterations) {
                        return blockJacobiPcg(a,
                                              b,
                                              std::move(x),
                                              preconditioner,
                                              iterations,
                                              &stats,
                                              "BJ-PARDISO-PCG",
                                              options.maxPcgIterations,
                                              options.pcgTolerance);
                    });
            } catch (const std::exception& err) {
                stats.status = "failed";
                stats.failureReason = err.what();
                stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                stats.workingSetAfterBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                temperature = initialTemperatureVector(mesh, physics);
            }
            results.push_back({stats, std::move(temperature)});
        }

        if (options.runBjIc) {
            SolverStatistics stats;
            stats.name = "BJ-IC-PCG";
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = peakWorkingSetBytes();
            std::vector<double> temperature;
            try {
                const auto setupStart = std::chrono::steady_clock::now();
                BlockJacobiIcPreconditioner preconditioner(mesh,
                                                           system,
                                                           options.icShift,
                                                           options.icScaling,
                                                           options.diagScalingEps,
                                                           options.localIcShiftMode,
                                                           options.icOrdering);
                const auto setupEnd = std::chrono::steady_clock::now();
                stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                stats.preconditionerBytes = preconditioner.memoryBytes();
                icDiagnostics = preconditioner.trialDiagnostics();
                stats.localDiagScaling = options.icScaling;
                stats.diagScalingEps = options.diagScalingEps;
                stats.localIcShiftMode = options.localIcShiftMode;
                addIcStatsToSolver(stats, preconditioner.diagnostics());
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                if (preconditioner.hasBreakdown()) {
                    stats.status = "failed";
                    stats.failureReason = "shifted IC factorization breakdown in at least one block";
                    stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                    stats.workingSetAfterBytes = currentWorkingSetBytes();
                    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                    temperature = initialTemperatureVector(mesh, physics);
                } else {
                    const std::filesystem::path residualHistoryPath =
                        outputDir / "bj_ic_pcg_residual_history.csv";
                    temperature = runAnalysisSolver(
                        steady,
                        "BJ-IC-PCG",
                        mesh, physics, mass, system, source, fixedAdjust, stats,
                        [&](const SparseMatrix& a, const std::vector<double>& b, std::vector<double> x, int& iterations) {
                            return blockJacobiPcg(a,
                                                  b,
                                                  std::move(x),
                                                  preconditioner,
                                                  iterations,
                                                  &stats,
                                                  "BJ-IC-PCG",
                                                  options.maxPcgIterations,
                                                  options.pcgTolerance,
                                                  residualHistoryPath,
                                                  !options.fastRun);
                        });
                }
            } catch (const std::exception& err) {
                stats.status = "failed";
                stats.failureReason = err.what();
                stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                stats.workingSetAfterBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                temperature = initialTemperatureVector(mesh, physics);
            }
            results.push_back({stats, std::move(temperature)});
        }

        const std::vector<double> dropTolerances =
            options.ilutDropTolerances.empty() ? std::vector<double>{1.0e-3, 1.0e-4, 1.0e-5} : options.ilutDropTolerances;
        const std::vector<int> fillFactors =
            options.ilutFillFactors.empty() ? std::vector<int>{5, 10, 20} : options.ilutFillFactors;

        if (options.runRasIc) {
            std::vector<ProgramOptions::RasIlutConfig> rasConfigs = options.rasIlutConfigs;
            if (rasConfigs.empty()) {
                rasConfigs.push_back({1, 0.0, 0});
            }
            for (const ProgramOptions::RasIlutConfig& rasConfig : rasConfigs) {
                SolverStatistics stats;
                stats.name = "FGMRES-RAS" + std::to_string(rasConfig.overlap) + "-IC";
                stats.workingSetBeforeBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = peakWorkingSetBytes();
                std::vector<double> temperature;
                try {
                    const auto setupStart = std::chrono::steady_clock::now();
                    RasIcPreconditioner preconditioner(mesh,
                                                       system,
                                                       rasConfig.overlap,
                                                       options.icShift,
                                                       options.icScaling,
                                                       options.diagScalingEps,
                                                       options.localIcShiftMode,
                                                       options.icOrdering);
                    const auto setupEnd = std::chrono::steady_clock::now();
                    stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                    stats.preconditionerBytes = preconditioner.memoryBytes();
                    stats.localDiagScaling = options.icScaling;
                    stats.diagScalingEps = options.diagScalingEps;
                    stats.localIcShiftMode = options.localIcShiftMode;
                    addIcStatsToSolver(stats, preconditioner.diagnostics());
                    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                    if (preconditioner.hasBreakdown()) {
                        stats.status = "failed";
                        stats.failureReason = "shifted IC factorization breakdown in at least one RAS block";
                        stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                        stats.workingSetAfterBytes = currentWorkingSetBytes();
                        stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                        temperature = initialTemperatureVector(mesh, physics);
                    } else {
                        const std::string solverName = stats.name;
                        temperature = runAnalysisSolver(
                            steady,
                            solverName,
                            mesh, physics, mass, system, source, fixedAdjust, stats,
                            [&](const SparseMatrix& a, const std::vector<double>& b, std::vector<double> x, int& iterations) {
                                return flexibleGmres(a,
                                                     b,
                                                     std::move(x),
                                                     preconditioner,
                                                     iterations,
                                                     &stats,
                                                     solverName,
                                                     options.maxPcgIterations,
                                                     options.pcgTolerance,
                                                     options.gmresRestart);
                            });
                    }
                } catch (const std::exception& err) {
                    stats.status = "failed";
                    stats.failureReason = err.what();
                    stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                    stats.workingSetAfterBytes = currentWorkingSetBytes();
                    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                    temperature = initialTemperatureVector(mesh, physics);
                }
                results.push_back({stats, std::move(temperature)});
            }
        }

        if (options.runBjIcCoarse) {
            SolverStatistics stats;
            stats.name = "FGMRES-BJ-IC-CoarseConst";
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = peakWorkingSetBytes();
            std::vector<double> temperature;
            try {
                const auto setupStart = std::chrono::steady_clock::now();
                BlockJacobiIcPreconditioner local(mesh,
                                                  system,
                                                  options.icShift,
                                                  options.icScaling,
                                                  options.diagScalingEps,
                                                  options.localIcShiftMode,
                                                  options.icOrdering);
                stats.localDiagScaling = options.icScaling;
                stats.diagScalingEps = options.diagScalingEps;
                stats.localIcShiftMode = options.localIcShiftMode;
                addIcStatsToSolver(stats, local.diagnostics());
                if (local.hasBreakdown()) {
                    throw std::runtime_error("shifted IC factorization breakdown in at least one block");
                }
                SubdomainConstantCoarseCorrectedPreconditioner<BlockJacobiIcPreconditioner> preconditioner(
                    std::move(local),
                    system,
                    buildSubdomainConstantCoarseSpace(mesh));
                const auto setupEnd = std::chrono::steady_clock::now();
                stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                stats.preconditionerBytes = preconditioner.memoryBytes();
                stats.coarseDim = preconditioner.coarseDim();
                stats.coarseSetupSeconds = preconditioner.coarseSetupSeconds();
                stats.coarseMatrixNnz = preconditioner.coarseMatrixNnz();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                const std::string solverName = stats.name;
                temperature = runAnalysisSolver(
                    steady,
                    solverName,
                    mesh, physics, mass, system, source, fixedAdjust, stats,
                    [&](const SparseMatrix& a, const std::vector<double>& b, std::vector<double> x, int& iterations) {
                        return flexibleGmres(a,
                                             b,
                                             std::move(x),
                                             preconditioner,
                                             iterations,
                                             &stats,
                                             solverName,
                                             options.maxPcgIterations,
                                             options.pcgTolerance,
                                             options.gmresRestart);
                    });
                stats.coarseSolveSeconds = preconditioner.coarseSolveSeconds();
                stats.coarseResidualNorm = preconditioner.coarseResidualNorm();
            } catch (const std::exception& err) {
                stats.status = "failed";
                stats.failureReason = err.what();
                stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                stats.workingSetAfterBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                temperature = initialTemperatureVector(mesh, physics);
            }
            results.push_back({stats, std::move(temperature)});
        }

        if (options.runBjIcCoarsePcg) {
            SolverStatistics stats;
            stats.name = "BJ-IC-Coarse-PCG";
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = peakWorkingSetBytes();
            std::vector<double> temperature;
            try {
                const auto setupStart = std::chrono::steady_clock::now();
                BlockJacobiIcPreconditioner local(mesh,
                                                  system,
                                                  options.icShift,
                                                  options.icScaling,
                                                  options.diagScalingEps,
                                                  options.localIcShiftMode,
                                                  options.icOrdering);
                stats.localDiagScaling = options.icScaling;
                stats.diagScalingEps = options.diagScalingEps;
                stats.localIcShiftMode = options.localIcShiftMode;
                addIcStatsToSolver(stats, local.diagnostics());
                if (local.hasBreakdown()) {
                    throw std::runtime_error("shifted IC factorization breakdown in at least one block");
                }
                SubdomainConstantCoarseCorrectedPreconditioner<BlockJacobiIcPreconditioner> preconditioner(
                    std::move(local),
                    system,
                    buildSubdomainConstantCoarseSpace(mesh));
                const auto setupEnd = std::chrono::steady_clock::now();
                stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                stats.preconditionerBytes = preconditioner.memoryBytes();
                stats.coarseDim = preconditioner.coarseDim();
                stats.coarseSetupSeconds = preconditioner.coarseSetupSeconds();
                stats.coarseMatrixNnz = preconditioner.coarseMatrixNnz();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                const std::filesystem::path residualHistoryPath =
                    outputDir / "bj_ic_coarse_pcg_residual_history.csv";
                temperature = runAnalysisSolver(
                    steady,
                    "BJ-IC-Coarse-PCG",
                    mesh, physics, mass, system, source, fixedAdjust, stats,
                    [&](const SparseMatrix& a, const std::vector<double>& b, std::vector<double> x, int& iterations) {
                        return blockJacobiPcg(a,
                                              b,
                                              std::move(x),
                                              preconditioner,
                                              iterations,
                                              &stats,
                                              "BJ-IC-Coarse-PCG",
                                              options.maxPcgIterations,
                                              options.pcgTolerance,
                                              residualHistoryPath,
                                              !options.fastRun);
                    });
                stats.coarseSolveSeconds = preconditioner.coarseSolveSeconds();
                stats.coarseResidualNorm = preconditioner.coarseResidualNorm();
            } catch (const std::exception& err) {
                stats.status = "failed";
                stats.failureReason = err.what();
                stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                stats.workingSetAfterBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                temperature = initialTemperatureVector(mesh, physics);
            }
            results.push_back({stats, std::move(temperature)});
        }

        if (options.runBjIlut) {
            for (double dropTolerance : dropTolerances) {
                for (int fillFactor : fillFactors) {
                    SolverStatistics stats;
                    stats.name = "FGMRES-BJ-ILUT-d" + std::to_string(dropTolerance)
                        + "-f" + std::to_string(fillFactor);
                    stats.workingSetBeforeBytes = currentWorkingSetBytes();
                    stats.peakWorkingSetBytes = peakWorkingSetBytes();
                    std::vector<double> temperature;
                    try {
                        const auto setupStart = std::chrono::steady_clock::now();
                        BlockJacobiIlutPreconditioner preconditioner(mesh,
                                                                    system,
                                                                    dropTolerance,
                                                                    fillFactor,
                                                                    options.localDiagScaling,
                                                                    options.diagScalingEps);
                        const auto setupEnd = std::chrono::steady_clock::now();
                        stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                        stats.preconditionerBytes = preconditioner.memoryBytes();
                        stats.localDiagScaling = options.localDiagScaling;
                        stats.diagScalingEps = options.diagScalingEps;
                        stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                        const std::string solverName = stats.name;
                        temperature = runAnalysisSolver(
                            steady,
                            solverName,
                            mesh, physics, mass, system, source, fixedAdjust, stats,
                            [&](const SparseMatrix& a, const std::vector<double>& b, std::vector<double> x, int& iterations) {
                                return flexibleGmres(a,
                                                     b,
                                                     std::move(x),
                                                     preconditioner,
                                                     iterations,
                                                     &stats,
                                                     solverName,
                                                     options.maxPcgIterations,
                                                     options.pcgTolerance,
                                                     options.gmresRestart);
                            });
                    } catch (const std::exception& err) {
                        stats.status = "failed";
                        stats.failureReason = err.what();
                        stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                        stats.workingSetAfterBytes = currentWorkingSetBytes();
                        stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                        temperature = initialTemperatureVector(mesh, physics);
                    }
                    results.push_back({stats, std::move(temperature)});
                }
            }
        }

        if (options.runBjIlutCoarse) {
            const std::vector<std::vector<double>> coarseVectors =
                buildSubdomainConstantCoarseVectors(mesh);
            for (double dropTolerance : dropTolerances) {
                for (int fillFactor : fillFactors) {
                    SolverStatistics stats;
                    stats.name = "FGMRES-BJ-ILUT-CoarseConst-d" + std::to_string(dropTolerance)
                        + "-f" + std::to_string(fillFactor);
                    stats.workingSetBeforeBytes = currentWorkingSetBytes();
                    stats.peakWorkingSetBytes = peakWorkingSetBytes();
                    std::vector<double> temperature;
                    try {
                        const auto setupStart = std::chrono::steady_clock::now();
                        BlockJacobiIlutPreconditioner local(mesh,
                                                           system,
                                                           dropTolerance,
                                                           fillFactor,
                                                           options.localDiagScaling,
                                                           options.diagScalingEps);
                        CoarseCorrectedPreconditioner<BlockJacobiIlutPreconditioner> preconditioner(
                            std::move(local),
                            system,
                            coarseVectors);
                        const auto setupEnd = std::chrono::steady_clock::now();
                        stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                        stats.preconditionerBytes = preconditioner.memoryBytes();
                        stats.localDiagScaling = options.localDiagScaling;
                        stats.diagScalingEps = options.diagScalingEps;
                        stats.coarseDim = preconditioner.coarseDim();
                        stats.coarseSetupSeconds = preconditioner.coarseSetupSeconds();
                        stats.coarseMatrixNnz = preconditioner.coarseMatrixNnz();
                        stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                        const std::string solverName = stats.name;
                        temperature = runAnalysisSolver(
                            steady,
                            solverName,
                            mesh, physics, mass, system, source, fixedAdjust, stats,
                            [&](const SparseMatrix& a, const std::vector<double>& b, std::vector<double> x, int& iterations) {
                                return flexibleGmres(a,
                                                     b,
                                                     std::move(x),
                                                     preconditioner,
                                                     iterations,
                                                     &stats,
                                                     solverName,
                                                     options.maxPcgIterations,
                                                     options.pcgTolerance,
                                                     options.gmresRestart);
                            });
                        stats.coarseSolveSeconds = preconditioner.coarseSolveSeconds();
                        stats.coarseResidualNorm = preconditioner.coarseResidualNorm();
                    } catch (const std::exception& err) {
                        stats.status = "failed";
                        stats.failureReason = err.what();
                        stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                        stats.workingSetAfterBytes = currentWorkingSetBytes();
                        stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                        temperature = initialTemperatureVector(mesh, physics);
                    }
                    results.push_back({stats, std::move(temperature)});
                }
            }
        }

        if (options.runRasIlut) {
            std::vector<ProgramOptions::RasIlutConfig> rasConfigs = options.rasIlutConfigs;
            if (rasConfigs.empty()) {
                for (double dropTolerance : dropTolerances) {
                    for (int fillFactor : fillFactors) {
                        ProgramOptions::RasIlutConfig config;
                        config.overlap = 1;
                        config.dropTolerance = dropTolerance;
                        config.fillFactor = fillFactor;
                        rasConfigs.push_back(config);
                    }
                }
            }
            for (const ProgramOptions::RasIlutConfig& rasConfig : rasConfigs) {
                    SolverStatistics stats;
                    stats.name = "FGMRES-RAS" + std::to_string(rasConfig.overlap)
                        + "-ILUT-d" + std::to_string(rasConfig.dropTolerance)
                        + "-f" + std::to_string(rasConfig.fillFactor);
                    stats.workingSetBeforeBytes = currentWorkingSetBytes();
                    stats.peakWorkingSetBytes = peakWorkingSetBytes();
                    std::vector<double> temperature;
                    try {
                        const auto setupStart = std::chrono::steady_clock::now();
                        RasIlutPreconditioner preconditioner(mesh,
                                                            system,
                                                            rasConfig.overlap,
                                                            rasConfig.dropTolerance,
                                                            rasConfig.fillFactor,
                                                            options.localDiagScaling,
                                                            options.diagScalingEps);
                        const auto setupEnd = std::chrono::steady_clock::now();
                        stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                        stats.preconditionerBytes = preconditioner.memoryBytes();
                        stats.localDiagScaling = options.localDiagScaling;
                        stats.diagScalingEps = options.diagScalingEps;
                        stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                        const std::string solverName = stats.name;
                        temperature = runAnalysisSolver(
                            steady,
                            solverName,
                            mesh, physics, mass, system, source, fixedAdjust, stats,
                            [&](const SparseMatrix& a, const std::vector<double>& b, std::vector<double> x, int& iterations) {
                                return flexibleGmres(a,
                                                     b,
                                                     std::move(x),
                                                     preconditioner,
                                                     iterations,
                                                     &stats,
                                                     solverName,
                                                     options.maxPcgIterations,
                                                     options.pcgTolerance,
                                                     options.gmresRestart);
                            });
                    } catch (const std::exception& err) {
                        stats.status = "failed";
                        stats.failureReason = err.what();
                        stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                        stats.workingSetAfterBytes = currentWorkingSetBytes();
                        stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                        temperature = initialTemperatureVector(mesh, physics);
                    }
                    results.push_back({stats, std::move(temperature)});
            }
        }

        if (options.runTwoLevelRasIlut) {
            std::vector<ProgramOptions::RasIlutConfig> twoLevelConfigs = options.rasIlutConfigs;
            if (twoLevelConfigs.empty()) {
                twoLevelConfigs.push_back({1, 1.0e-3, 5});
                twoLevelConfigs.push_back({2, 1.0e-4, 10});
            }
            const std::vector<std::vector<double>> coarseVectors =
                buildSubdomainConstantCoarseVectors(mesh);
            for (const ProgramOptions::RasIlutConfig& rasConfig : twoLevelConfigs) {
                SolverStatistics stats;
                stats.name = "FGMRES-RAS" + std::to_string(rasConfig.overlap)
                    + "-ILUT-CoarseConst-d" + std::to_string(rasConfig.dropTolerance)
                    + "-f" + std::to_string(rasConfig.fillFactor);
                stats.workingSetBeforeBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = peakWorkingSetBytes();
                std::vector<double> temperature;
                try {
                    const auto setupStart = std::chrono::steady_clock::now();
                    CoarseCorrectedRasIlutPreconditioner preconditioner(mesh,
                                                                        system,
                                                                        rasConfig.overlap,
                                                                        rasConfig.dropTolerance,
                                                                        rasConfig.fillFactor,
                                                                        coarseVectors,
                                                                        options.localDiagScaling,
                                                                        options.diagScalingEps);
                    const auto setupEnd = std::chrono::steady_clock::now();
                    stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                    stats.preconditionerBytes = preconditioner.memoryBytes();
                    stats.localDiagScaling = options.localDiagScaling;
                    stats.diagScalingEps = options.diagScalingEps;
                    stats.coarseDim = preconditioner.coarseDim();
                    stats.coarseSetupSeconds = preconditioner.coarseSetupSeconds();
                    stats.coarseMatrixNnz = preconditioner.coarseMatrixNnz();
                    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                    const std::string solverName = stats.name;
                    temperature = runAnalysisSolver(
                        steady,
                        solverName,
                        mesh, physics, mass, system, source, fixedAdjust, stats,
                        [&](const SparseMatrix& a, const std::vector<double>& b, std::vector<double> x, int& iterations) {
                            return flexibleGmres(a,
                                                 b,
                                                 std::move(x),
                                                 preconditioner,
                                                 iterations,
                                                 &stats,
                                                 solverName,
                                                 options.maxPcgIterations,
                                                 options.pcgTolerance,
                                                 options.gmresRestart);
                        });
                    stats.coarseSolveSeconds = preconditioner.coarseSolveSeconds();
                    stats.coarseResidualNorm = preconditioner.coarseResidualNorm();
                } catch (const std::exception& err) {
                    stats.status = "failed";
                    stats.failureReason = err.what();
                    stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                    stats.workingSetAfterBytes = currentWorkingSetBytes();
                    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                    temperature = initialTemperatureVector(mesh, physics);
                }
                results.push_back({stats, std::move(temperature)});
            }
        }

        if (options.runDeflatedRasIlut) {
            std::vector<ProgramOptions::RasIlutConfig> deflatedConfigs = options.rasIlutConfigs;
            if (deflatedConfigs.empty()) {
                deflatedConfigs.push_back({1, 1.0e-3, 5});
            }
            const std::vector<std::vector<double>> eigenVectors =
                buildEigenCoarseVectors(spectralSummary, options.deflationModes);
            for (const ProgramOptions::RasIlutConfig& rasConfig : deflatedConfigs) {
                SolverStatistics stats;
                stats.name = "FGMRES-RAS" + std::to_string(rasConfig.overlap)
                    + "-ILUT-DeflatedEig" + std::to_string(static_cast<int>(eigenVectors.size()))
                    + "-d" + std::to_string(rasConfig.dropTolerance)
                    + "-f" + std::to_string(rasConfig.fillFactor);
                stats.workingSetBeforeBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = peakWorkingSetBytes();
                std::vector<double> temperature;
                try {
                    if (eigenVectors.empty()) {
                        throw std::runtime_error("deflated RAS requested but no near-zero eigenvectors are available");
                    }
                    const auto setupStart = std::chrono::steady_clock::now();
                    CoarseCorrectedRasIlutPreconditioner preconditioner(mesh,
                                                                        system,
                                                                        rasConfig.overlap,
                                                                        rasConfig.dropTolerance,
                                                                        rasConfig.fillFactor,
                                                                        eigenVectors,
                                                                        options.localDiagScaling,
                                                                        options.diagScalingEps);
                    const auto setupEnd = std::chrono::steady_clock::now();
                    stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                    stats.preconditionerBytes = preconditioner.memoryBytes();
                    stats.localDiagScaling = options.localDiagScaling;
                    stats.diagScalingEps = options.diagScalingEps;
                    stats.coarseDim = preconditioner.coarseDim();
                    stats.coarseSetupSeconds = preconditioner.coarseSetupSeconds();
                    stats.coarseMatrixNnz = preconditioner.coarseMatrixNnz();
                    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                    const std::string solverName = stats.name;
                    temperature = runAnalysisSolver(
                        steady,
                        solverName,
                        mesh, physics, mass, system, source, fixedAdjust, stats,
                        [&](const SparseMatrix& a, const std::vector<double>& b, std::vector<double> x, int& iterations) {
                            return flexibleGmres(a,
                                                 b,
                                                 std::move(x),
                                                 preconditioner,
                                                 iterations,
                                                 &stats,
                                                 solverName,
                                                 options.maxPcgIterations,
                                                 options.pcgTolerance,
                                                 options.gmresRestart);
                        });
                    stats.coarseSolveSeconds = preconditioner.coarseSolveSeconds();
                    stats.coarseResidualNorm = preconditioner.coarseResidualNorm();
                } catch (const std::exception& err) {
                    stats.status = "failed";
                    stats.failureReason = err.what();
                    stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                    stats.workingSetAfterBytes = currentWorkingSetBytes();
                    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                    temperature = initialTemperatureVector(mesh, physics);
                }
                results.push_back({stats, std::move(temperature)});
            }
        }

        if (options.runInterfaceDeflatedRasIlut) {
            std::vector<ProgramOptions::RasIlutConfig> deflatedConfigs = options.rasIlutConfigs;
            if (deflatedConfigs.empty()) {
                deflatedConfigs.push_back({1, 1.0e-3, 5});
                deflatedConfigs.push_back({2, 1.0e-4, 10});
            }
            const std::vector<std::vector<double>> eigenVectors =
                buildInterfaceEigenCoarseVectors(spectralSummary, assemblyDiagnostics, options.deflationModes);
            for (const ProgramOptions::RasIlutConfig& rasConfig : deflatedConfigs) {
                SolverStatistics stats;
                stats.name = "FGMRES-RAS" + std::to_string(rasConfig.overlap)
                    + "-ILUT-InterfaceEig" + std::to_string(static_cast<int>(eigenVectors.size()))
                    + "-d" + std::to_string(rasConfig.dropTolerance)
                    + "-f" + std::to_string(rasConfig.fillFactor);
                stats.workingSetBeforeBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = peakWorkingSetBytes();
                std::vector<double> temperature;
                try {
                    if (eigenVectors.empty()) {
                        throw std::runtime_error("interface-deflated RAS requested but no interface eigenvectors are available");
                    }
                    const auto setupStart = std::chrono::steady_clock::now();
                    CoarseCorrectedRasIlutPreconditioner preconditioner(mesh,
                                                                        system,
                                                                        rasConfig.overlap,
                                                                        rasConfig.dropTolerance,
                                                                        rasConfig.fillFactor,
                                                                        eigenVectors,
                                                                        options.localDiagScaling,
                                                                        options.diagScalingEps);
                    const auto setupEnd = std::chrono::steady_clock::now();
                    stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                    stats.preconditionerBytes = preconditioner.memoryBytes();
                    stats.localDiagScaling = options.localDiagScaling;
                    stats.diagScalingEps = options.diagScalingEps;
                    stats.coarseDim = preconditioner.coarseDim();
                    stats.coarseSetupSeconds = preconditioner.coarseSetupSeconds();
                    stats.coarseMatrixNnz = preconditioner.coarseMatrixNnz();
                    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                    const std::string solverName = stats.name;
                    temperature = runAnalysisSolver(
                        steady,
                        solverName,
                        mesh, physics, mass, system, source, fixedAdjust, stats,
                        [&](const SparseMatrix& a, const std::vector<double>& b, std::vector<double> x, int& iterations) {
                            return flexibleGmres(a,
                                                 b,
                                                 std::move(x),
                                                 preconditioner,
                                                 iterations,
                                                 &stats,
                                                 solverName,
                                                 options.maxPcgIterations,
                                                 options.pcgTolerance,
                                                 options.gmresRestart);
                        });
                    stats.coarseSolveSeconds = preconditioner.coarseSolveSeconds();
                    stats.coarseResidualNorm = preconditioner.coarseResidualNorm();
                } catch (const std::exception& err) {
                    stats.status = "failed";
                    stats.failureReason = err.what();
                    stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                    stats.workingSetAfterBytes = currentWorkingSetBytes();
                    stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                    temperature = initialTemperatureVector(mesh, physics);
                }
                results.push_back({stats, std::move(temperature)});
            }
        }

        if (options.runBjPardisoGeneral) {
            SolverStatistics stats;
            stats.name = "FGMRES-BJ-PARDISO-General";
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = peakWorkingSetBytes();
            std::vector<double> temperature;
            try {
                const auto setupStart = std::chrono::steady_clock::now();
                BlockJacobiGeneralPardisoPreconditioner preconditioner(mesh, system);
                const auto setupEnd = std::chrono::steady_clock::now();
                stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                stats.preconditionerBytes = preconditioner.memoryBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                temperature = runAnalysisSolver(
                    steady,
                    "FGMRES-BJ-PARDISO-General",
                    mesh, physics, mass, system, source, fixedAdjust, stats,
                    [&](const SparseMatrix& a, const std::vector<double>& b, std::vector<double> x, int& iterations) {
                        return flexibleGmres(a,
                                             b,
                                             std::move(x),
                                             preconditioner,
                                             iterations,
                                             &stats,
                                             "FGMRES-BJ-PARDISO-General",
                                             options.maxPcgIterations,
                                             options.pcgTolerance,
                                             options.gmresRestart);
                    });
            } catch (const std::exception& err) {
                stats.status = "failed";
                stats.failureReason = err.what();
                stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                stats.workingSetAfterBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                temperature = initialTemperatureVector(mesh, physics);
            }
            results.push_back({stats, std::move(temperature)});
        }

        if (options.runSchwarzPrecondFgmres) {
            SolverStatistics stats;
            stats.name = "Schwarz-Precond-FGMRES";
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = peakWorkingSetBytes();
            std::vector<double> temperature;
            std::vector<double> finalRhs;
            try {
                const auto setupStart = std::chrono::steady_clock::now();
                SchwarzOptions preconditionerOptions = physics.schwarz;
                preconditionerOptions.standaloneMode = "algebraic";
                preconditionerOptions.transmission = "none";
                SchwarzBlockSolver preconditioner(mesh, system, preconditionerOptions, &physics);
                const auto setupEnd = std::chrono::steady_clock::now();
                stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                stats.preconditionerBytes = preconditioner.memoryBytes();
                preconditioner.copyTimingTo(stats);
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                writeSchwarzBlockMatrixCheck(preconditioner,
                                             outputDir / "schwarz_block_matrix_check.csv");

                SchwarzPreconditionedKrylovRunResult krylovResult =
                    runSchwarzPreconditionedFgmresAnalysis(steady,
                                                           stats.name,
                                                           mesh,
                                                           physics,
                                                           mass,
                                                           system,
                                                           source,
                                                           fixedAdjust,
                                                           preconditioner,
                                                           stats,
                                                           options.maxPcgIterations,
                                                           options.pcgTolerance,
                                                           options.gmresRestart);
                temperature = std::move(krylovResult.temperature);
                finalRhs = std::move(krylovResult.finalRhs);
                preconditioner.copyTimingTo(stats);
                allKrylovRows.insert(allKrylovRows.end(),
                                     krylovResult.iterationRows.begin(),
                                     krylovResult.iterationRows.end());
                writeKrylovIterationLog(krylovResult.iterationRows,
                                        outputDir / "krylov_iterations_schwarz_precond_fgmres.csv");
            } catch (const std::exception& err) {
                stats.status = "failed";
                stats.failureReason = err.what();
                stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                stats.workingSetAfterBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                temperature = initialTemperatureVector(mesh, physics);
            }
            results.push_back({stats, std::move(temperature), std::move(finalRhs)});
        }

        if (options.runSchwarzPrecondFgmresTwoLevel) {
            SolverStatistics stats;
            stats.name = "Schwarz-Precond-FGMRES-Two-Level";
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = peakWorkingSetBytes();
            std::vector<double> temperature;
            std::vector<double> finalRhs;
            try {
                const auto setupStart = std::chrono::steady_clock::now();
                SchwarzOptions preconditionerOptions = physics.schwarz;
                preconditionerOptions.standaloneMode = "algebraic";
                preconditionerOptions.transmission = "none";
                SchwarzBlockSolver localPreconditioner(mesh, system, preconditionerOptions, &physics);
                TwoLevelSchwarzPreconditioner preconditioner(
                    localPreconditioner,
                    mesh,
                    system,
                    buildRequestedCoarseSpace(mesh, physics, options.coarseSpace),
                    options.coarseCorrection);
                const auto setupEnd = std::chrono::steady_clock::now();
                stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                stats.preconditionerBytes = preconditioner.memoryBytes();
                stats.coarseSpace = options.coarseSpace;
                stats.coarseCorrection = options.coarseCorrection;
                preconditioner.copyTimingTo(stats);
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                writeSchwarzBlockMatrixCheck(localPreconditioner,
                                             outputDir / "schwarz_block_matrix_check_two_level.csv");

                SchwarzPreconditionedKrylovRunResult krylovResult =
                    runSchwarzPreconditionedFgmresAnalysis(steady,
                                                           stats.name,
                                                           mesh,
                                                           physics,
                                                           mass,
                                                           system,
                                                           source,
                                                           fixedAdjust,
                                                           preconditioner,
                                                           stats,
                                                           options.maxPcgIterations,
                                                           options.pcgTolerance,
                                                           options.gmresRestart);
                temperature = std::move(krylovResult.temperature);
                finalRhs = std::move(krylovResult.finalRhs);
                preconditioner.copyTimingTo(stats);
                allKrylovRows.insert(allKrylovRows.end(),
                                     krylovResult.iterationRows.begin(),
                                     krylovResult.iterationRows.end());
                writeKrylovIterationLog(krylovResult.iterationRows,
                                        outputDir / "krylov_iterations_schwarz_precond_fgmres_two_level.csv");
            } catch (const std::exception& err) {
                stats.status = "failed";
                stats.failureReason = err.what();
                stats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                stats.workingSetAfterBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                temperature = initialTemperatureVector(mesh, physics);
            }
            results.push_back({stats, std::move(temperature), std::move(finalRhs)});
        }

        if (options.runSchwarz) {
            SolverStatistics stats;
            const std::string normalizedStandaloneMode =
                normalizeSchwarzTransmissionModeName(physics.schwarz.standaloneMode);
            stats.name = isPhysicalSchwarzStandaloneMode(normalizedStandaloneMode)
                ? "Schwarz-SIPG-physical-" + normalizedStandaloneMode + "-" + physics.schwarz.type
                : "Schwarz-SIPG-algebraic-" + physics.schwarz.type;
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = peakWorkingSetBytes();
            std::vector<double> temperature;
            try {
                const auto setupStart = std::chrono::steady_clock::now();
                SchwarzBlockSolver schwarzSolver(mesh,
                                                 system,
                                                 physics.schwarz,
                                                 &physics,
                                                 &physicalTransmissionSystem);
                const auto setupEnd = std::chrono::steady_clock::now();
                const double blockAssemblySeconds =
                    std::chrono::duration<double>(setupEnd - setupStart).count();
                stats.setupSeconds = blockAssemblySeconds;
                stats.preconditionerBytes = schwarzSolver.memoryBytes();
                schwarzSolver.copyTimingTo(stats);
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                writeSchwarzBlockMatrixCheck(schwarzSolver,
                                             outputDir / "schwarz_block_matrix_check.csv");

                SchwarzRunResult schwarzResult = runSchwarzAnalysis(
                    steady,
                    stats.name,
                    mesh,
                    physics,
                    mass,
                    system,
                    source,
                    fixedAdjust,
                    schwarzSolver,
                    blockAssemblySeconds,
                    stats);
                temperature = std::move(schwarzResult.temperature);
                writeSchwarzIterationStats(schwarzResult.stepStats,
                                           outputDir / "schwarz_iterations.csv");
                writeSchwarzSubdomainSolveStats(schwarzResult.subdomainStats,
                                                outputDir / "schwarz_subdomain_solve_times.csv");
                if (physics.schwarz.writeInterfaceFlux) {
                    writeSchwarzInterfaceFluxRows(schwarzResult.interfaceFluxRows,
                                                  outputDir / "schwarz_interface_flux.csv");
                }
            } catch (const std::exception& err) {
                stats.status = "failed";
                stats.failureReason = err.what();
                stats.parallelWorkers = 1;
                stats.workingSetAfterBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                temperature = initialTemperatureVector(mesh, physics);
            }
            results.push_back({stats, std::move(temperature)});
        }

        if (options.runLocalRom) {
            SolverStatistics stats;
            stats.name = options.localMorMode == "corrected"
                ? "Local-POD-Schur-ROM-Corrected" : "Local-POD-Schur-ROM-Pure";
            stats.solverMethod = "independent_local_interior_rom_full_interface_schur";
            stats.preconditioner = options.localMorMode == "corrected"
                ? "stage1_exact_schur_residual_gated_correction"
                : "not_used_by_pure_local_rom";
            stats.localSolver = "dense_reduced_subdomain_solve";
            stats.coarseSpace = "none_full_order_physical_interface";
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            std::vector<double> temperature;
            std::vector<double> finalRhs = source;
            try {
                if (!steady) {
                    throw std::runtime_error(
                        "Local-ROM-Schur Milestones 1-3 support only the fixed steady matrix.");
                }
                ddm_schur::Options schurOptions;
                schurOptions.maxIterations = options.maxPcgIterations;
                schurOptions.restart = options.gmresRestart;
                schurOptions.relativeTolerance = options.pcgTolerance;
                schurOptions.coarseLinearXY = options.schurLinearXYCoarse;
                schurOptions.coarseLinearZ = options.schurLinearZCoarse;
                schurOptions.coarseGlobalQuadraticZ = options.schurGlobalQuadraticZCoarse;
                schurOptions.coarseInterfacePatches = options.schurInterfacePatchCoarse;
                schurOptions.coarseInterfaceLinearXY = options.schurInterfacePatchLinearXY;
                schurOptions.coarseEnergyAdaptive = options.schurEnergyAdaptiveCoarse;
                schurOptions.energyMaxModesPerDomain = options.schurEnergyMaxModesPerDomain;
                schurOptions.energySubspaceIterations = options.schurEnergySubspaceIterations;
                schurOptions.energyEigenvalueThreshold = options.schurEnergyEigenvalueThreshold;
                schurOptions.coarseGlobalSlow = options.schurGlobalSlowCoarse;
                schurOptions.globalSlowModes = options.schurGlobalSlowModes;
                schurOptions.globalSlowSubspaceDimension = options.schurGlobalSlowSubspaceDimension;
                schurOptions.proxyDiagnostics = options.schurProxyDiagnostics;
                schurOptions.proxyEnabled = options.schurProxyEnabled;
                schurOptions.proxyDisableCoarse = options.schurProxyDisableCoarse;
                schurOptions.proxyHighConductivityThreshold = options.schurProxyHighKThreshold;
                schurOptions.proxyUseMaterialConnectivity = options.schurProxyUseMaterialConnectivity;
                schurOptions.proxyRing = options.schurProxyRing;
                schurOptions.proxyProbeColumns = options.schurProxyProbeColumns;
                schurOptions.proxyBlockSize = options.schurProxyBlockSize;
                schurOptions.proxyValidateBlockEquivalence =
                    options.schurProxyValidateBlockEquivalence;
                schurOptions.localSolveThreads = options.schurLocalSolveThreads;
                schurOptions.localPardisoThreads = options.schurLocalPardisoThreads;
                schurOptions.proxyCacheEnabled = options.schurProxyCacheEnabled;
                schurOptions.proxyCachePath = options.schurProxyCachePath.string();
                schurOptions.proxyOutputDirectory = outputDir.string();

                mor::local::Options localOptions;
                localOptions.generate = options.localMorGenerate;
                localOptions.loadPath = options.localMorLoadPath;
                localOptions.savePath = options.localMorSavePath;
                localOptions.method = options.localMorMethod;
                localOptions.mode = options.localMorMode;
                localOptions.interfaceMode = options.localInterfaceMode;
                localOptions.rank = options.localMorRank;
                localOptions.rankPerSubdomain = options.localMorRankPerSubdomain;
                localOptions.rankFile = options.localMorRankFile;
                localOptions.energyTolerance = options.localMorEnergyTolerance;
                localOptions.singularValueTolerance = options.localMorSingularValueTolerance;
                localOptions.trainingCases = options.localMorTrainingCases;
                localOptions.validationCases = options.localMorValidationCases;
                localOptions.testCases = options.localMorTestCases;
                localOptions.seed = options.morSeed;
                localOptions.compareFom = options.localMorCompareFom;
                localOptions.reuseIdenticalSubdomains =
                    options.localMorReuseIdenticalSubdomains;
                localOptions.matrixFreeInterfaceThreshold =
                    options.localMorMatrixFreeInterfaceThreshold;

                const mor::local::WorkflowResult local =
                    mor::local::runLocalRomSchurWorkflow(
                        mesh, system, physics, source, heatOnlySource, fixedAdjust,
                        schurOptions, localOptions, outputDir);
                temperature = local.nominalTemperature;
                stats.coarseEnabled = local.coarseDimension > 0;
                stats.coarseDim = local.coarseDimension;
                stats.coarseSolveSeconds = local.coarseSolveSeconds;
                stats.totalIterations = local.correctionIterations;
                stats.maxIterations = local.correctionIterations;
                stats.setupSeconds = local.setupSeconds;
                stats.solveSeconds = local.onlineSeconds;
                stats.finalRelativeResidual = local.interfaceResidual;
                stats.trueRelativeResidual = local.globalResidual;
                stats.status = local.status;
                if (stats.status != "success") {
                    stats.failureReason = "Local-ROM-Schur accuracy or corrected residual gate failed";
                }
            } catch (const std::exception& err) {
                stats.status = "failed";
                stats.failureReason = err.what();
                temperature = initialTemperatureVector(mesh, physics);
            }
            applyDirichletRhs(mesh, fixedAdjust, finalRhs);
            stats.workingSetAfterBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes,
                                                peakWorkingSetBytes());
            results.push_back({stats, std::move(temperature), std::move(finalRhs)});
        }

        if (options.runReducedSchur) {
            SolverStatistics stats;
            stats.name = options.morMode == "corrected"
                ? "Reduced-Schur-Corrected" : "Reduced-Schur-Pure";
            stats.solverMethod = "exact_interface_reduced_schur";
            stats.preconditioner = options.morMode == "corrected"
                || options.morSnapshotSolver == "schur"
                ? "stage1_ring1_proxy_xyz_for_snapshots_or_correction"
                : "not_applied_in_pure_reduced_schur";
            stats.localSolver = "PARDISO_reused";
            stats.coarseSpace = "source_response_gram_pod";
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            std::vector<double> temperature;
            std::vector<double> finalRhs = source;
            try {
                if (!steady) {
                    throw std::runtime_error(
                        "reduced-schur Stage 2A supports only a fixed steady matrix.");
                }
                ddm_schur::Options schurOptions;
                schurOptions.maxIterations = options.maxPcgIterations;
                schurOptions.restart = options.gmresRestart;
                schurOptions.relativeTolerance = options.pcgTolerance;
                schurOptions.coarseLinearXY = options.schurLinearXYCoarse;
                schurOptions.coarseLinearZ = options.schurLinearZCoarse;
                schurOptions.coarseGlobalQuadraticZ = options.schurGlobalQuadraticZCoarse;
                schurOptions.coarseInterfacePatches = options.schurInterfacePatchCoarse;
                schurOptions.coarseInterfaceLinearXY = options.schurInterfacePatchLinearXY;
                schurOptions.coarseEnergyAdaptive = options.schurEnergyAdaptiveCoarse;
                schurOptions.energyMaxModesPerDomain = options.schurEnergyMaxModesPerDomain;
                schurOptions.energySubspaceIterations = options.schurEnergySubspaceIterations;
                schurOptions.energyEigenvalueThreshold = options.schurEnergyEigenvalueThreshold;
                schurOptions.coarseGlobalSlow = options.schurGlobalSlowCoarse;
                schurOptions.globalSlowModes = options.schurGlobalSlowModes;
                schurOptions.globalSlowSubspaceDimension = options.schurGlobalSlowSubspaceDimension;
                schurOptions.proxyDiagnostics = options.schurProxyDiagnostics;
                schurOptions.proxyEnabled = options.schurProxyEnabled;
                schurOptions.proxyDisableCoarse = options.schurProxyDisableCoarse;
                schurOptions.proxyHighConductivityThreshold = options.schurProxyHighKThreshold;
                schurOptions.proxyUseMaterialConnectivity = options.schurProxyUseMaterialConnectivity;
                schurOptions.proxyRing = options.schurProxyRing;
                schurOptions.proxyProbeColumns = options.schurProxyProbeColumns;
                schurOptions.proxyBlockSize = options.schurProxyBlockSize;
                schurOptions.proxyValidateBlockEquivalence = options.schurProxyValidateBlockEquivalence;
                schurOptions.localSolveThreads = options.schurLocalSolveThreads;
                schurOptions.localPardisoThreads = options.schurLocalPardisoThreads;
                schurOptions.proxyCacheEnabled = options.schurProxyCacheEnabled;
                schurOptions.proxyCachePath = options.schurProxyCachePath.string();
                schurOptions.proxyOutputDirectory = outputDir.string();

                mor::Options morOptions;
                morOptions.generateModel = options.morGenerateModel;
                morOptions.loadModelPath = options.morLoadModelPath;
                morOptions.saveModelPath = options.morSaveModelPath;
                if (morOptions.generateModel && morOptions.saveModelPath.empty()) {
                    morOptions.saveModelPath = outputDir / "mor_model";
                }
                morOptions.rank = options.morRank;
                morOptions.energyTolerance = options.morEnergyTolerance;
                morOptions.singularValueTolerance = options.morSingularValueTolerance;
                morOptions.trainingCases = options.morTrainingCases;
                morOptions.validationCases = options.morValidationCases;
                morOptions.testCases = options.morTestCases;
                morOptions.seed = options.morSeed;
                morOptions.mode = options.morMode;
                morOptions.snapshotSolver = options.morSnapshotSolver;
                morOptions.exactInteriorRecovery = options.morExactInteriorRecovery;
                morOptions.compareFom = options.morCompareFom;
                morOptions.rankSweep = options.morRankSweep;
                morOptions.interiorMode = options.morInteriorMode;
                morOptions.interiorRank = options.morInteriorRank;
                morOptions.interiorEnergyTolerance = options.morInteriorEnergyTolerance;
                morOptions.interiorSingularValueTolerance =
                    options.morInteriorSingularValueTolerance;
                morOptions.interiorRankSweep = options.morInteriorRankSweep;
                morOptions.precomputePowerResponse = options.morPrecomputePowerResponse;
                morOptions.storagePrecision = options.morStoragePrecision;
                morOptions.reportIoTime = options.morReportIoTime;
                morOptions.compareInteriorModes = options.morCompareInteriorModes;
                morOptions.deploymentRhsCount = options.morDeploymentRhsCount;
                morOptions.parametricGenerate = options.morParametricGenerate;
                morOptions.parametricLoadPath = options.morParametricLoadPath;
                morOptions.parametricSavePath = options.morParametricSavePath;
                morOptions.matrixParameter = options.morMatrixParameter;
                morOptions.parameterSubdomain = options.morParameterSubdomain;
                morOptions.parameterRegionId = options.morParameterRegionId;
                morOptions.parameterMinimum = options.morParameterMinimum;
                morOptions.parameterMaximum = options.morParameterMaximum;
                morOptions.parameterReference = options.morParameterReference;
                morOptions.parameterValue = options.morParameterValue;
                morOptions.parameterTrainingCount = options.morParameterTrainingCount;
                morOptions.parameterValidationCount = options.morParameterValidationCount;
                morOptions.parameterTestCount = options.morParameterTestCount;
                morOptions.interfaceRank = options.morInterfaceRank;
                morOptions.localRank = options.morLocalRank;
                morOptions.parametricMode = options.morParametricMode;
                morOptions.allowExtrapolation = options.morAllowExtrapolation;
                morOptions.parametricAffineValidationOnly =
                    options.morParametricAffineValidationOnly;
                morOptions.onlinePowersW = options.morOnlinePowersW;
                const bool parametric = morOptions.parametricGenerate
                    || !morOptions.parametricLoadPath.empty();
                const mor::WorkflowResult reduced = parametric
                    ? mor::parametric::runParametricReducedSchurWorkflow(
                        mesh, system, physics, source, heatOnlySource, fixedAdjust,
                        schurOptions, morOptions, outputDir)
                    : mor::runReducedSchurWorkflow(
                        mesh, system, physics, source, heatOnlySource, fixedAdjust,
                        schurOptions, morOptions, outputDir);
                temperature = reduced.nominalTemperature;
                stats.coarseEnabled = true;
                stats.coarseDim = reduced.selectedRank;
                stats.totalIterations = reduced.correctionIterations;
                stats.maxIterations = reduced.correctionIterations;
                stats.setupSeconds = reduced.setupSeconds;
                stats.solveSeconds = reduced.nominalOnlineSeconds;
                stats.finalRelativeResidual = reduced.nominalInterfaceResidual;
                stats.trueRelativeResidual = reduced.nominalGlobalResidual;
                stats.status = reduced.status;
                if (stats.status != "success") {
                    stats.failureReason = "Reduced-Schur corrected solve reached the iteration limit";
                }
            } catch (const std::exception& err) {
                stats.status = "failed";
                stats.failureReason = err.what();
                temperature = initialTemperatureVector(mesh, physics);
            }
            applyDirichletRhs(mesh, fixedAdjust, finalRhs);
            stats.workingSetAfterBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes,
                                                peakWorkingSetBytes());
            results.push_back({stats, std::move(temperature), std::move(finalRhs)});
        }

        if (options.runDdmSchur) {
            SolverStatistics stats;
            stats.name = "DDM-Schur-FGMRES";
            stats.solverMethod = "matrix_free_schur";
            const std::string localCoarseModes = options.schurLinearXYCoarse
                ? (options.schurLinearZCoarse ? "constant_plus_linear_xyz" : "constant_plus_linear_xy")
                : (options.schurLinearZCoarse ? "constant_plus_linear_z" : "constant");
            stats.preconditioner = options.schurProxyEnabled
                ? std::string("balanced_two_level_schur_proxy_ring")
                    + std::to_string(options.schurProxyRing)
                    + (options.schurProxyDisableCoarse ? "_only" : "_plus_")
                    + (options.schurProxyDisableCoarse ? "" : localCoarseModes)
                    + (options.schurGlobalSlowCoarse ? "_plus_global_slow" : "")
                : (options.schurGlobalSlowCoarse
                ? "balanced_two_level_schur_subdomain_" + localCoarseModes
                    + "_plus_global_slow"
                : (options.schurEnergyAdaptiveCoarse
                ? "balanced_two_level_schur_subdomain_" + localCoarseModes
                    + "_plus_energy_adaptive"
                : (options.schurInterfacePatchCoarse
                    ? std::string("balanced_two_level_schur_interface_patch_constant")
                        + (options.schurInterfacePatchLinearXY ? "_plus_linear_xy" : "")
                    : "balanced_two_level_schur_subdomain_" + localCoarseModes
                        + (options.schurGlobalQuadraticZCoarse ? "_plus_global_quadratic_z" : ""))));
            stats.localSolver = "PARDISO";
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = peakWorkingSetBytes();
            std::vector<double> temperature;
            std::vector<double> finalRhs;
            ddm_schur::Report lastReport;
            try {
                ddm_schur::Options schurOptions;
                schurOptions.maxIterations = options.maxPcgIterations;
                schurOptions.restart = options.gmresRestart;
                schurOptions.relativeTolerance = options.pcgTolerance;
                schurOptions.coarseLinearXY = options.schurLinearXYCoarse;
                schurOptions.coarseLinearZ = options.schurLinearZCoarse;
                schurOptions.coarseGlobalQuadraticZ = options.schurGlobalQuadraticZCoarse;
                schurOptions.coarseInterfacePatches = options.schurInterfacePatchCoarse;
                schurOptions.coarseInterfaceLinearXY = options.schurInterfacePatchLinearXY;
                schurOptions.coarseEnergyAdaptive = options.schurEnergyAdaptiveCoarse;
                schurOptions.energyMaxModesPerDomain = options.schurEnergyMaxModesPerDomain;
                schurOptions.energySubspaceIterations = options.schurEnergySubspaceIterations;
                schurOptions.energyEigenvalueThreshold = options.schurEnergyEigenvalueThreshold;
                schurOptions.coarseGlobalSlow = options.schurGlobalSlowCoarse;
                schurOptions.globalSlowModes = options.schurGlobalSlowModes;
                schurOptions.globalSlowSubspaceDimension = options.schurGlobalSlowSubspaceDimension;
                schurOptions.proxyDiagnostics = options.schurProxyDiagnostics;
                schurOptions.proxyEnabled = options.schurProxyEnabled;
                schurOptions.proxyDisableCoarse = options.schurProxyDisableCoarse;
                schurOptions.proxyHighConductivityThreshold = options.schurProxyHighKThreshold;
                schurOptions.proxyUseMaterialConnectivity = options.schurProxyUseMaterialConnectivity;
                schurOptions.proxyRing = options.schurProxyRing;
                schurOptions.proxyProbeColumns = options.schurProxyProbeColumns;
                schurOptions.proxyBlockSize = options.schurProxyBlockSize;
                schurOptions.proxyValidateBlockEquivalence =
                    options.schurProxyValidateBlockEquivalence;
                schurOptions.localSolveThreads = options.schurLocalSolveThreads;
                schurOptions.localPardisoThreads = options.schurLocalPardisoThreads;
                schurOptions.proxyCacheEnabled = options.schurProxyCacheEnabled;
                schurOptions.proxyCachePath = options.schurProxyCachePath.string();
                schurOptions.proxyOutputDirectory = outputDir.string();
                ddm_schur::DdmSchurSolver solver(mesh, system, physics, schurOptions);
                stats.setupSeconds = solver.setupReport().setupSeconds;
                stats.preconditionerBytes = solver.setupReport().memoryBytes;
                stats.coarseEnabled = true;
                stats.coarseDim = solver.setupReport().coarseDimension;
                if (options.schurMultiRhsCount < 0) {
                    throw std::runtime_error("--schur-multi-rhs-count must be nonnegative.");
                }
                if (options.schurMultiRhsCount > 0 && !steady) {
                    throw std::runtime_error(
                        "--schur-multi-rhs-count currently supports only a fixed steady matrix.");
                }

                if (options.schurMultiRhsCount == 0) {
                    temperature = runAnalysisSolver(
                        steady,
                        stats.name,
                        mesh, physics, mass, system, source, fixedAdjust, stats,
                        [&](const SparseMatrix&, const std::vector<double>& b, std::vector<double>, int& iterations) {
                            ddm_schur::SolveResult schurResult = solver.solve(b);
                            lastReport = schurResult.report;
                            iterations = schurResult.report.iterations;
                            stats.finalRelativeResidual = schurResult.report.interfaceRelativeResidual;
                            stats.rasLocalSolveApplySeconds = schurResult.report.localSolveSeconds;
                            stats.coarseSolveSeconds = schurResult.report.coarseSolveSeconds;
                            stats.status = schurResult.report.status == "success" ? "success" : "failed";
                            if (stats.status != "success") {
                                stats.failureReason = "Schur FGMRES reached the iteration limit";
                            }
                            return std::move(schurResult.temperature);
                        },
                        &finalRhs);
                } else {
                    struct MultiRhsRow {
                        int index = 0;
                        double multiplierMin = 1.0;
                        double multiplierMax = 1.0;
                        double rhsGenerationSeconds = 0.0;
                        double globalRelativeResidual = 0.0;
                        double temperatureMin = 0.0;
                        double temperatureMax = 0.0;
                        double temperatureAverage = 0.0;
                        ddm_schur::Report report;
                    };
                    std::vector<MultiRhsRow> multiRows;
                    multiRows.reserve(static_cast<std::size_t>(options.schurMultiRhsCount));
                    std::vector<double> nonHeatSource(source.size(), 0.0);
                    for (std::size_t i = 0; i < source.size(); ++i) {
                        nonHeatSource[i] = source[i] - heatOnlySource[i];
                    }

                    const auto onlineBatchStart = std::chrono::steady_clock::now();
                    for (int rhsIndex = 0; rhsIndex < options.schurMultiRhsCount; ++rhsIndex) {
                        const auto rhsStart = std::chrono::steady_clock::now();
                        std::vector<double> rhsSource;
                        std::vector<double> multipliers(physics.heatSources.size(), 1.0);
                        if (rhsIndex == 0) {
                            // Preserve the original production RHS exactly as the
                            // regression member of the ensemble.
                            rhsSource = source;
                        } else {
                            for (std::size_t sourceIndex = 0;
                                 sourceIndex < multipliers.size(); ++sourceIndex) {
                                const double phase = 0.7548776662466927
                                    * static_cast<double>((rhsIndex + 1)
                                        * static_cast<int>(sourceIndex + 1))
                                    + 0.371 * static_cast<double>(rhsIndex);
                                multipliers[sourceIndex] = 1.0 + 0.25 * std::sin(phase);
                            }
                            std::vector<double> newHeatSource;
                            assembleHeatSourceVector(mesh, physics, multipliers, newHeatSource);
                            if (physics.thermalSourceScale != 1.0) {
                                for (double& value : newHeatSource) {
                                    value *= physics.thermalSourceScale;
                                }
                            }
                            rhsSource = nonHeatSource;
                            for (std::size_t i = 0; i < rhsSource.size(); ++i) {
                                rhsSource[i] += newHeatSource[i];
                            }
                        }
                        std::vector<double> rhs = rhsSource;
                        applyDirichletRhs(mesh, fixedAdjust, rhs);
                        const double rhsGenerationSeconds = std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - rhsStart).count();

                        ddm_schur::SolveResult schurResult = solver.solve(rhs);
                        const std::vector<double> trueResidual =
                            trueResidualVector(system, schurResult.temperature, rhs);
                        const double globalRelativeResidual = l2Norm(trueResidual)
                            / std::max(1.0e-300, l2Norm(rhs));
                        const auto minmax = std::minmax_element(
                            schurResult.temperature.begin(), schurResult.temperature.end());

                        MultiRhsRow row;
                        row.index = rhsIndex;
                        row.multiplierMin = multipliers.empty()
                            ? 1.0
                            : *std::min_element(multipliers.begin(), multipliers.end());
                        row.multiplierMax = multipliers.empty()
                            ? 1.0
                            : *std::max_element(multipliers.begin(), multipliers.end());
                        row.rhsGenerationSeconds = rhsGenerationSeconds;
                        row.globalRelativeResidual = globalRelativeResidual;
                        row.temperatureMin = *minmax.first;
                        row.temperatureMax = *minmax.second;
                        row.temperatureAverage = std::accumulate(
                            schurResult.temperature.begin(),
                            schurResult.temperature.end(), 0.0)
                            / static_cast<double>(schurResult.temperature.size());
                        row.report = schurResult.report;
                        multiRows.push_back(row);

                        if (rhsIndex == 0) {
                            temperature = std::move(schurResult.temperature);
                            finalRhs = rhs;
                            lastReport = row.report;
                            stats.totalIterations = row.report.iterations;
                            stats.maxIterations = row.report.iterations;
                            stats.finalRelativeResidual = row.report.interfaceRelativeResidual;
                            stats.trueRelativeResidual = globalRelativeResidual;
                            stats.rasLocalSolveApplySeconds = row.report.localSolveSeconds;
                            stats.coarseSolveSeconds = row.report.coarseSolveSeconds;
                            stats.status = row.report.status == "success" ? "success" : "failed";
                            if (stats.status != "success") {
                                stats.failureReason = "Schur FGMRES reached the iteration limit";
                            }
                        }

                        std::cout << std::setprecision(12)
                                  << "[Schur multi-RHS] " << (rhsIndex + 1)
                                  << '/' << options.schurMultiRhsCount
                                  << " iterations=" << row.report.iterations
                                  << " interface_residual=" << row.report.interfaceRelativeResidual
                                  << " global_residual=" << row.globalRelativeResidual
                                  << " online=" << row.report.totalSolveSeconds << " s\n";
                    }
                    const double onlineBatchWallSeconds = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - onlineBatchStart).count();

                    std::ofstream runsOut(outputDir / "schur_multi_rhs_runs.csv");
                    runsOut << "rhs_index,source_multiplier_min,source_multiplier_max,iterations,"
                            << "interface_relative_residual,global_relative_residual,temperature_min,"
                            << "temperature_max,temperature_average,rhs_generation_seconds,"
                            << "condensed_rhs_seconds,fgmres_seconds,recovery_seconds,online_seconds,"
                            << "local_solve_calls,proxy_solve_calls,status\n";
                    runsOut << std::setprecision(16);
                    for (const MultiRhsRow& row : multiRows) {
                        runsOut << row.index << ','
                                << row.multiplierMin << ','
                                << row.multiplierMax << ','
                                << row.report.iterations << ','
                                << row.report.interfaceRelativeResidual << ','
                                << row.globalRelativeResidual << ','
                                << row.temperatureMin << ','
                                << row.temperatureMax << ','
                                << row.temperatureAverage << ','
                                << row.rhsGenerationSeconds << ','
                                << row.report.condensedRhsSeconds << ','
                                << row.report.fgmresSeconds << ','
                                << row.report.recoverySeconds << ','
                                << row.report.totalSolveSeconds << ','
                                << row.report.localSolveCalls << ','
                                << row.report.proxySolveCalls << ','
                                << row.report.status << '\n';
                    }

                    std::vector<int> aggregateCounts;
                    for (const int requested : {10, 50, 100}) {
                        if (requested <= options.schurMultiRhsCount) {
                            aggregateCounts.push_back(requested);
                        }
                    }
                    if (aggregateCounts.empty()
                        || aggregateCounts.back() != options.schurMultiRhsCount) {
                        aggregateCounts.push_back(options.schurMultiRhsCount);
                    }
                    std::ofstream aggregateOut(outputDir / "schur_multi_rhs_summary.csv");
                    aggregateOut << "rhs_count,setup_once_seconds,average_rhs_generation_seconds,"
                                 << "average_condensed_rhs_seconds,average_fgmres_seconds,"
                                 << "average_recovery_seconds,average_online_seconds,total_online_seconds,"
                                 << "batch_wall_seconds,average_iterations,max_interface_residual,"
                                 << "max_global_residual,local_symbolic_calls,local_numerical_calls,"
                                 << "proxy_symbolic_calls,proxy_numerical_calls,interface_graph_builds,"
                                 << "coloring_builds,proxy_pattern_builds,proxy_factorizations\n";
                    aggregateOut << std::setprecision(16);
                    for (const int count : aggregateCounts) {
                        double rhsSeconds = 0.0;
                        double condensedSeconds = 0.0;
                        double fgmresSeconds = 0.0;
                        double recoverySeconds = 0.0;
                        double onlineSeconds = 0.0;
                        double iterationSum = 0.0;
                        double maxInterfaceResidual = 0.0;
                        double maxGlobalResidual = 0.0;
                        for (int i = 0; i < count; ++i) {
                            const MultiRhsRow& row = multiRows[static_cast<std::size_t>(i)];
                            rhsSeconds += row.rhsGenerationSeconds;
                            condensedSeconds += row.report.condensedRhsSeconds;
                            fgmresSeconds += row.report.fgmresSeconds;
                            recoverySeconds += row.report.recoverySeconds;
                            onlineSeconds += row.report.totalSolveSeconds;
                            iterationSum += static_cast<double>(row.report.iterations);
                            maxInterfaceResidual = std::max(
                                maxInterfaceResidual, row.report.interfaceRelativeResidual);
                            maxGlobalResidual = std::max(
                                maxGlobalResidual, row.globalRelativeResidual);
                        }
                        const double inverseCount = 1.0 / static_cast<double>(count);
                        aggregateOut << count << ','
                                     << solver.setupReport().setupSeconds << ','
                                     << rhsSeconds * inverseCount << ','
                                     << condensedSeconds * inverseCount << ','
                                     << fgmresSeconds * inverseCount << ','
                                     << recoverySeconds * inverseCount << ','
                                     << onlineSeconds * inverseCount << ','
                                     << onlineSeconds << ','
                                     << (count == options.schurMultiRhsCount
                                         ? onlineBatchWallSeconds
                                         : std::numeric_limits<double>::quiet_NaN()) << ','
                                     << iterationSum * inverseCount << ','
                                     << maxInterfaceResidual << ','
                                     << maxGlobalResidual << ','
                                     << solver.setupReport().localSymbolicAnalysisCalls << ','
                                     << solver.setupReport().localNumericalFactorizationCalls << ','
                                     << solver.setupReport().proxySymbolicCalls << ','
                                     << solver.setupReport().proxyNumericalCalls << ','
                                     << 1 << ',' << 1 << ',' << 1 << ',' << 1 << '\n';
                    }
                    stats.solveSeconds = lastReport.totalSolveSeconds;
                    stats.workingSetAfterBytes = currentWorkingSetBytes();
                    stats.peakWorkingSetBytes = std::max(
                        stats.peakWorkingSetBytes, peakWorkingSetBytes());
                }

                // Global harmonic modes are prepared lazily from the true
                // condensed RHS.  Refresh the generic solver statistics after
                // solve so they report the final coarse dimension and split
                // deferred setup from the actual solve time.
                stats.setupSeconds = lastReport.setupSeconds;
                stats.solveSeconds = lastReport.totalSolveSeconds;
                stats.preconditionerBytes = lastReport.memoryBytes;
                stats.coarseDim = lastReport.coarseDimension;

                std::ofstream summary(outputDir / "schur_solver_summary.csv");
                summary << "domains,total_dofs,interface_dofs,interior_dofs,interface_patch_count,coarse_dimension,energy_candidate_modes,global_slow_candidate_dimension,iterations,schur_matvecs,"
                        << "local_solve_calls,local_symbolic_analysis_calls,local_numerical_factorization_calls,"
                        << "setup_seconds,energy_eigen_setup_seconds,energy_selected_eigenvalue_min,energy_selected_eigenvalue_max,"
                        << "global_slow_setup_seconds,global_slow_estimated_lambda_max,global_slow_selected_ritz_min,global_slow_selected_ritz_max,"
                        << "proxy_graph_edges,proxy_diagnostic_columns,proxy_diagnostic_schur_applies,proxy_diagnostics_seconds,"
                        << "proxy_ring3_coverage,proxy_ring3_operator_error,proxy_ring3_estimated_nnz,proxy_ring3_memory_estimate_bytes,proxy_recommended,"
                        << "proxy_nnz,proxy_density,proxy_colors,proxy_probing_schur_applies,proxy_probing_block_size,proxy_probing_block_calls,proxy_validation_schur_applies,proxy_block_maximum_difference,proxy_block_relative_difference,proxy_value_hash,proxy_symbolic_calls,proxy_numerical_calls,proxy_solve_calls,"
                        << "proxy_setup_seconds,proxy_symbolic_seconds,proxy_numerical_seconds,proxy_solve_seconds,proxy_matrix_cache_hit,proxy_factor_cache_hit,proxy_symmetry_error,"
                        << "proxy_minimum_test_rayleigh,proxy_diagonal_shift,proxy_diagonal_compensation,proxy_memory_bytes,"
                        << "local_factorization_seconds,local_symbolic_analysis_seconds,"
                        << "local_numerical_factorization_seconds,local_solve_seconds,schur_apply_seconds,coarse_solve_seconds,"
                        << "condensed_rhs_seconds,interface_solve_seconds,fgmres_seconds,recovery_seconds,"
                        << "total_solve_seconds,total_seconds,relative_residual,memory_bytes,status\n";
                summary << std::setprecision(16)
                        << lastReport.domains << ','
                        << lastReport.totalDofs << ','
                        << lastReport.interfaceDofs << ','
                        << lastReport.interiorDofs << ','
                        << lastReport.interfacePatchCount << ','
                        << lastReport.coarseDimension << ','
                        << lastReport.energyCandidateModes << ','
                        << lastReport.globalSlowCandidateDimension << ','
                        << lastReport.iterations << ','
                        << lastReport.schurMatvecs << ','
                        << lastReport.localSolveCalls << ','
                        << lastReport.localSymbolicAnalysisCalls << ','
                        << lastReport.localNumericalFactorizationCalls << ','
                        << lastReport.setupSeconds << ','
                        << lastReport.energyEigenSetupSeconds << ','
                        << lastReport.energySelectedEigenvalueMin << ','
                        << lastReport.energySelectedEigenvalueMax << ','
                        << lastReport.globalSlowSetupSeconds << ','
                        << lastReport.globalSlowEstimatedLambdaMax << ','
                        << lastReport.globalSlowSelectedRitzMin << ','
                        << lastReport.globalSlowSelectedRitzMax << ','
                        << lastReport.proxyGraphEdges << ','
                        << lastReport.proxyProbeColumns << ','
                        << lastReport.proxySchurApplies << ','
                        << lastReport.proxyDiagnosticsSeconds << ','
                        << lastReport.proxyRing3Coverage << ','
                        << lastReport.proxyRing3OperatorError << ','
                        << lastReport.proxyRing3EstimatedNnz << ','
                        << lastReport.proxyRing3MemoryEstimateBytes << ','
                        << (lastReport.proxyRecommended ? 1 : 0) << ','
                        << lastReport.proxyNnz << ','
                        << lastReport.proxyDensity << ','
                        << lastReport.proxyColors << ','
                        << lastReport.proxyProbingSchurApplies << ','
                        << lastReport.proxyProbingBlockSize << ','
                        << lastReport.proxyProbingBlockCalls << ','
                        << lastReport.proxyValidationSchurApplies << ','
                        << lastReport.proxyBlockMaximumDifference << ','
                        << lastReport.proxyBlockRelativeDifference << ','
                        << lastReport.proxyValueHash << ','
                        << lastReport.proxySymbolicCalls << ','
                        << lastReport.proxyNumericalCalls << ','
                        << lastReport.proxySolveCalls << ','
                        << lastReport.proxySetupSeconds << ','
                        << lastReport.proxySymbolicSeconds << ','
                        << lastReport.proxyNumericalSeconds << ','
                        << lastReport.proxySolveSeconds << ','
                        << (lastReport.proxyMatrixCacheHit ? 1 : 0) << ','
                        << (lastReport.proxyFactorCacheHit ? 1 : 0) << ','
                        << lastReport.proxySymmetryError << ','
                        << lastReport.proxyMinimumTestRayleigh << ','
                        << lastReport.proxyDiagonalShift << ','
                        << lastReport.proxyDiagonalCompensation << ','
                        << lastReport.proxyMemoryBytes << ','
                        << lastReport.localFactorizationSeconds << ','
                        << lastReport.localSymbolicAnalysisSeconds << ','
                        << lastReport.localNumericalFactorizationSeconds << ','
                        << lastReport.localSolveSeconds << ','
                        << lastReport.schurApplySeconds << ','
                        << lastReport.coarseSolveSeconds << ','
                        << lastReport.condensedRhsSeconds << ','
                        << lastReport.interfaceSolveSeconds << ','
                        << lastReport.fgmresSeconds << ','
                        << lastReport.recoverySeconds << ','
                        << lastReport.totalSolveSeconds << ','
                        << lastReport.totalSeconds << ','
                        << lastReport.interfaceRelativeResidual << ','
                        << lastReport.memoryBytes << ','
                        << lastReport.status << '\n';

                std::ofstream subdomainOut(outputDir / "schur_subdomain_performance.csv");
                subdomainOut << "subdomain_id,interior_dofs,interface_dofs,phase11_time,"
                             << "phase22_time,phase33_calls,phase33_total_time,"
                             << "phase33_average_time,factor_memory_bytes\n";
                subdomainOut << std::setprecision(16);
                double phase33Total = 0.0;
                int phase33Calls = 0;
                int slowestSubdomain = -1;
                double slowestPhase33 = -1.0;
                std::size_t factorMemoryTotal = 0;
                for (const ddm_schur::SubdomainPerformance& item :
                     lastReport.subdomainPerformance) {
                    const double average = item.phase33Calls > 0
                        ? item.phase33Seconds / static_cast<double>(item.phase33Calls)
                        : 0.0;
                    subdomainOut << item.subdomainId << ',' << item.interiorDofs << ','
                                 << item.interfaceDofs << ',' << item.phase11Seconds << ','
                                 << item.phase22Seconds << ',' << item.phase33Calls << ','
                                 << item.phase33Seconds << ',' << average << ','
                                 << item.factorMemoryBytes << '\n';
                    phase33Total += item.phase33Seconds;
                    phase33Calls += item.phase33Calls;
                    factorMemoryTotal += item.factorMemoryBytes;
                    if (item.phase33Seconds > slowestPhase33) {
                        slowestPhase33 = item.phase33Seconds;
                        slowestSubdomain = item.subdomainId;
                    }
                }
                std::ofstream subdomainSummaryOut(
                    outputDir / "schur_subdomain_performance_summary.csv");
                subdomainSummaryOut
                    << "subdomain_count,phase33_calls,phase33_total_time,"
                    << "phase33_average_time,slowest_subdomain,slowest_phase33_time,"
                    << "factor_memory_bytes\n"
                    << std::setprecision(16)
                    << lastReport.subdomainPerformance.size() << ',' << phase33Calls << ','
                    << phase33Total << ','
                    << (phase33Calls > 0
                        ? phase33Total / static_cast<double>(phase33Calls) : 0.0) << ','
                    << slowestSubdomain << ',' << slowestPhase33 << ','
                    << factorMemoryTotal << '\n';
            } catch (const std::exception& err) {
                stats.status = "failed";
                stats.failureReason = err.what();
                stats.workingSetAfterBytes = currentWorkingSetBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                temperature = initialTemperatureVector(mesh, physics);
            }
            results.push_back({stats, std::move(temperature), std::move(finalRhs)});
        }

        if (options.runDirect) {
            const bool runSpdDirect = options.directMode == "spd" || options.directMode == "both";
            const bool runGeneralDirect = options.directMode == "general" || options.directMode == "both";
            if (!runSpdDirect && !runGeneralDirect) {
                throw std::runtime_error("--direct-mode must be spd, general, or both.");
            }
            if (runSpdDirect) {
            SolverStatistics spdStats;
            spdStats.name = "Global-PARDISO-SPD-Direct";
            spdStats.workingSetBeforeBytes = currentWorkingSetBytes();
            spdStats.peakWorkingSetBytes = peakWorkingSetBytes();
            std::vector<double> spdTemperature;
            std::vector<double> spdFinalRhs;
            try {
                const auto setupStart = std::chrono::steady_clock::now();
                const std::vector<MatrixEntry> globalEntries = sparseMatrixEntries(system);
                SubdomainDirectSolver spdSolver(n, globalEntries);
                const auto setupEnd = std::chrono::steady_clock::now();
                spdStats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                spdStats.preconditionerBytes = spdSolver.memoryBytes();
                spdStats.peakWorkingSetBytes = std::max(spdStats.peakWorkingSetBytes, peakWorkingSetBytes());
                spdTemperature = runAnalysisSolver(
                    steady,
                    "Global-PARDISO-SPD-Direct",
                    mesh, physics, mass, system, source, fixedAdjust, spdStats,
                    [&](const SparseMatrix&, const std::vector<double>& b, std::vector<double>, int& iterations) {
                        iterations = 0;
                        std::vector<double> x;
                        spdSolver.solve(b, x);
                        spdStats.status = "success";
                        return x;
                    },
                    &spdFinalRhs);
            } catch (const std::exception& err) {
                spdStats.status = "failed";
                spdStats.failureReason = err.what();
                spdStats.parallelWorkers = static_cast<int>(solverParallelWorkers());
                spdStats.workingSetAfterBytes = currentWorkingSetBytes();
                spdStats.peakWorkingSetBytes = std::max(spdStats.peakWorkingSetBytes, peakWorkingSetBytes());
                spdTemperature = initialTemperatureVector(mesh, physics);
            }
            results.push_back({spdStats, std::move(spdTemperature), std::move(spdFinalRhs)});
            }

            if (runGeneralDirect) {
            SolverStatistics stats;
            stats.name = "Global-PARDISO-General-Direct";
            stats.workingSetBeforeBytes = currentWorkingSetBytes();
            stats.peakWorkingSetBytes = peakWorkingSetBytes();
            std::vector<double> temperature;
            std::vector<double> finalRhs;
            {
                const auto setupStart = std::chrono::steady_clock::now();
                const std::vector<MatrixEntry> globalEntries = sparseMatrixEntries(system);
                GeneralSparseDirectSolver solver(n, globalEntries);
                const auto setupEnd = std::chrono::steady_clock::now();
                stats.setupSeconds = std::chrono::duration<double>(setupEnd - setupStart).count();
                stats.preconditionerBytes = solver.memoryBytes();
                stats.peakWorkingSetBytes = std::max(stats.peakWorkingSetBytes, peakWorkingSetBytes());
                temperature = runAnalysisSolver(
                    steady,
                    "Global-PARDISO-General-Direct",
                    mesh, physics, mass, system, source, fixedAdjust, stats,
                    [&](const SparseMatrix&, const std::vector<double>& b, std::vector<double>, int& iterations) {
                        iterations = 0;
                        std::vector<double> x;
                        solver.solve(b, x);
                        stats.status = "success";
                        return x;
                    },
                    &finalRhs);
            }
            results.push_back({stats, std::move(temperature), std::move(finalRhs)});
            }
        }

        if (results.empty()) {
            throw std::runtime_error("No solvers were selected. Use --solvers local-rom,reduced-schur,schur,schur-direct-exact,direct or another supported solver.");
        }

        size_t referenceIndex = results.size() - 1;
        for (size_t i = 0; i < results.size(); ++i) {
            if (results[i].stats.name.find("Global-PARDISO") != std::string::npos) {
                referenceIndex = i;
            }
        }
        const std::vector<double>& referenceTemperature = results[referenceIndex].temperature;
        if (steady) {
            std::vector<double> steadyRhs = source;
            applyDirichletRhs(mesh, fixedAdjust, steadyRhs);
            const DofPartitionCounts dofCounts = countDofPartitions(mesh, assemblyDiagnostics);
            const double rhsL2 = std::max(1.0e-300, l2Norm(steadyRhs));
            const double rhsInf = std::max(1.0, infNorm(steadyRhs));
            const std::vector<double> directTrueResidual =
                trueResidualVector(system, referenceTemperature, steadyRhs);
            const double directTrueRelativeResidual = l2Norm(directTrueResidual) / rhsL2;
            const double directTrueInfResidual = infNorm(directTrueResidual) / rhsInf;
            const double matrixInfNormValue = collectDetailedDiagnostics
                ? matrixInfinityNorm(system)
                : std::numeric_limits<double>::quiet_NaN();
            if (collectDetailedDiagnostics && !std::isfinite(spectralSummary.lambdaMaxAbs)) {
                spectralSummary.lambdaMaxAbs = estimateMaxAbsEigenvaluePower(system, 50);
            }

            std::ofstream residualOut(outputDir / "rram_true_residual_check.csv");
            residualOut << "solver,status,reported_final_relative_residual,true_relative_residual,true_inf_residual,"
                        << "direct_true_relative_residual,direct_true_inf_residual,vector_size_direct,vector_size_iter,"
                        << "free_dof_count,dirichlet_dof_count,interface_dof_count,interior_dof_count,"
                        << "ras_overlap,ras_total_overlap_dof_count,ras_total_core_dof_count\n";
            residualOut << std::setprecision(16);

            std::ofstream errorOut(outputDir / "rram_solution_error_by_dof_type.csv");
            errorOut << "solver,status,all_dofs_relative_l2_error,free_dofs_relative_l2_error,"
                     << "dirichlet_dofs_max_error,interface_dofs_relative_l2_error,"
                     << "interior_core_dofs_relative_l2_error,direct_dirichlet_max_deviation,"
                     << "iter_dirichlet_max_deviation,free_dof_count,dirichlet_dof_count,"
                     << "interface_dof_count,interior_core_dof_count\n";
            errorOut << std::setprecision(16);

            std::ofstream updateOut(outputDir / "rram_ras_update_count_diagnostics.csv");
            updateOut << "solver,overlap,subdomain,core_dof_count,local_dof_count,overlap_dof_count,"
                      << "min_update_count,max_update_count,zero_update_count,multi_update_count\n";
            updateOut << std::setprecision(16);

            std::ofstream backwardOut(outputDir / "rram_backward_error_diagnostics.csv");
            backwardOut << "solver,status,absolute_residual_l2,true_relative_residual,true_inf_residual,"
                        << "normwise_backward_error,inf_norm_backward_error,lambda_max_abs,lambda_min_abs,"
                        << "condition_est_abs,estimated_forward_error_bound\n";
            backwardOut << std::setprecision(16);
            std::ofstream qualityOut(outputDir / "rram_iterative_convergence_quality.csv");
            qualityOut << "solver,status,true_relative_residual,absolute_residual_l2,normwise_backward_error,"
                       << "inf_norm_backward_error,deflated_residual,coarse_space_residual,"
                       << "condition_est_abs,estimated_forward_error_bound,relative_l2_error_vs_direct\n";
            qualityOut << std::setprecision(16);

            std::vector<RasUpdateDiagnosticsRow> allRasRows;
            for (const SolverRunResult& result : results) {
                const int rasOverlap = parseRasOverlapFromSolverName(result.stats.name);
                int rasTotalOverlapDofs = -1;
                int rasTotalCoreDofs = -1;
                if (rasOverlap >= 0) {
                    std::vector<RasUpdateDiagnosticsRow> rows =
                        computeRasUpdateDiagnostics(mesh, system, result.stats.name, rasOverlap);
                    rasTotalOverlapDofs = 0;
                    rasTotalCoreDofs = 0;
                    for (const RasUpdateDiagnosticsRow& row : rows) {
                        rasTotalOverlapDofs += row.overlapDofCount;
                        rasTotalCoreDofs += row.coreDofCount;
                        allRasRows.push_back(row);
                        updateOut << csvEscape(row.solver) << ','
                                  << row.overlap << ','
                                  << row.subdomain << ','
                                  << row.coreDofCount << ','
                                  << row.localDofCount << ','
                                  << row.overlapDofCount << ','
                                  << row.minUpdateCount << ','
                                  << row.maxUpdateCount << ','
                                  << row.zeroUpdateCount << ','
                                  << row.multiUpdateCount << '\n';
                    }
                }

                const std::vector<double> trueResidual =
                    trueResidualVector(system, result.temperature, steadyRhs);
                const double trueResidualL2 = l2Norm(trueResidual);
                const double trueResidualInf = infNorm(trueResidual);
                const double trueRelativeResidual = trueResidualL2 / rhsL2;
                const double normwiseBackwardError = collectDetailedDiagnostics
                    ? trueResidualL2 / (std::max(0.0, spectralSummary.lambdaMaxAbs) * l2Norm(result.temperature) + rhsL2)
                    : std::numeric_limits<double>::quiet_NaN();
                const double infBackwardError = collectDetailedDiagnostics
                    ? trueResidualInf / (matrixInfNormValue * infNorm(result.temperature) + rhsInf)
                    : std::numeric_limits<double>::quiet_NaN();
                residualOut << csvEscape(result.stats.name) << ','
                            << result.stats.status << ','
                            << result.stats.finalRelativeResidual << ','
                            << trueRelativeResidual << ','
                            << trueResidualInf / rhsInf << ','
                            << directTrueRelativeResidual << ','
                            << directTrueInfResidual << ','
                            << referenceTemperature.size() << ','
                            << result.temperature.size() << ','
                            << dofCounts.freeDofs << ','
                            << dofCounts.dirichletDofs << ','
                            << dofCounts.interfaceDofs << ','
                            << dofCounts.interiorDofs << ','
                            << rasOverlap << ','
                            << rasTotalOverlapDofs << ','
                            << rasTotalCoreDofs << '\n';
                backwardOut << csvEscape(result.stats.name) << ','
                            << result.stats.status << ','
                            << trueResidualL2 << ','
                            << trueRelativeResidual << ','
                            << trueResidualInf / rhsInf << ','
                            << normwiseBackwardError << ','
                            << infBackwardError << ','
                            << spectralSummary.lambdaMaxAbs << ','
                            << spectralSummary.lambdaMinAbs << ','
                            << spectralSummary.conditionEstAbs << ','
                            << spectralSummary.conditionEstAbs * trueRelativeResidual << '\n';
                const double relativeErrorVsDirect = relativeL2Difference(result.temperature, referenceTemperature);
                qualityOut << csvEscape(result.stats.name) << ','
                           << result.stats.status << ','
                           << trueRelativeResidual << ','
                           << trueResidualL2 << ','
                           << normwiseBackwardError << ','
                           << infBackwardError << ','
                           << std::numeric_limits<double>::quiet_NaN() << ','
                           << std::numeric_limits<double>::quiet_NaN() << ','
                           << spectralSummary.conditionEstAbs << ','
                           << spectralSummary.conditionEstAbs * trueRelativeResidual << ','
                           << relativeErrorVsDirect << '\n';

                std::vector<char> dirichletMask(mesh.nodes.size(), 0);
                for (size_t i = 0; i < mesh.nodes.size(); ++i) {
                    dirichletMask[i] = mesh.nodes[i].dirichlet ? 1 : 0;
                }
                const double directDirichletDeviation = maxAbsDifferenceWhere(
                    referenceTemperature,
                    [&]() {
                        std::vector<double> exact(mesh.nodes.size(), 0.0);
                        for (size_t i = 0; i < mesh.nodes.size(); ++i) {
                            exact[i] = mesh.nodes[i].dirichletValue;
                        }
                        return exact;
                    }(),
                    dirichletMask);
                const double iterDirichletDeviation = maxAbsDifferenceWhere(
                    result.temperature,
                    [&]() {
                        std::vector<double> exact(mesh.nodes.size(), 0.0);
                        for (size_t i = 0; i < mesh.nodes.size(); ++i) {
                            exact[i] = mesh.nodes[i].dirichletValue;
                        }
                        return exact;
                    }(),
                    dirichletMask);
                const int interiorCoreCount = static_cast<int>(std::count_if(mesh.nodes.begin(),
                                                                            mesh.nodes.end(),
                                                                            [&](const Node& node) {
                                                                                const size_t idx = static_cast<size_t>(&node - mesh.nodes.data());
                                                                                const bool interfaceDof = !assemblyDiagnostics.interfaceDof.empty()
                                                                                    && assemblyDiagnostics.interfaceDof[idx] != 0;
                                                                                return !node.dirichlet && !interfaceDof;
                                                                            }));
                errorOut << csvEscape(result.stats.name) << ','
                         << result.stats.status << ','
                         << relativeL2Difference(result.temperature, referenceTemperature) << ','
                         << relativeL2DifferenceWhere(result.temperature, referenceTemperature, [&](size_t i) {
                                return !mesh.nodes[i].dirichlet;
                            }) << ','
                         << maxAbsDifferenceWhere(result.temperature, referenceTemperature, dirichletMask) << ','
                         << relativeL2DifferenceWhere(result.temperature, referenceTemperature, [&](size_t i) {
                                return !assemblyDiagnostics.interfaceDof.empty()
                                    && assemblyDiagnostics.interfaceDof[i] != 0;
                            }) << ','
                         << relativeL2DifferenceWhere(result.temperature, referenceTemperature, [&](size_t i) {
                                const bool interfaceDof = !assemblyDiagnostics.interfaceDof.empty()
                                    && assemblyDiagnostics.interfaceDof[i] != 0;
                                return !mesh.nodes[i].dirichlet && !interfaceDof;
                            }) << ','
                         << directDirichletDeviation << ','
                         << iterDirichletDeviation << ','
                         << dofCounts.freeDofs << ','
                         << dofCounts.dirichletDofs << ','
                         << dofCounts.interfaceDofs << ','
                         << interiorCoreCount << '\n';
            }

            std::ofstream mappingOut(outputDir / "rram_mapping_dirichlet_scaling_check.md");
            mappingOut << "# RRAM Mapping, Dirichlet, Scaling, and RAS Check\n\n";
            mappingOut << "- Matrix/vector space: full global DOF vector, same `A_final` and same steady RHS used for direct and iterative solvers.\n";
            mappingOut << "- Direct vector size: " << referenceTemperature.size() << "\n";
            mappingOut << "- Free DOFs: " << dofCounts.freeDofs << "\n";
            mappingOut << "- Dirichlet DOFs: " << dofCounts.dirichletDofs << "\n";
            mappingOut << "- Interface DOFs: " << dofCounts.interfaceDofs << "\n";
            mappingOut << "- Interior/non-interface DOFs: " << dofCounts.interiorDofs << "\n";
            mappingOut << "- Direct true relative residual: " << directTrueRelativeResidual << "\n";
            mappingOut << "- Direct true infinity residual: " << directTrueInfResidual << "\n\n";
            mappingOut << "## Dirichlet Handling\n\n";
            mappingOut << "Dirichlet rows are retained in the full system and RHS values are set by `applyDirichletRhs`. "
                       << "Both direct and iterative solutions are full vectors, so no eliminated/free-only vector is mixed into the comparison.\n\n";
            mappingOut << "## Scaling\n\n";
            mappingOut << "No global diagonal, row-column, or system equilibration scaling is applied in the tested FGMRES/BJ-ILUT/RAS-ILUT path. "
                       << "Residual diagnostics in `rram_true_residual_check.csv` are unscaled true residuals computed from the original `A_final` and RHS. "
                       << "Local IC scaling, when IC is selected, is internal to that preconditioner and is not used for the RAS/BJ ILUT tests here.\n\n";
            mappingOut << "## RAS Write-Back\n\n";
            mappingOut << "RAS local solves include overlap DOFs, but correction write-back uses only core DOFs. "
                       << "The update-count diagnostics report the intended full-domain partition: every global DOF should have update count 1, "
                       << "zero-update count 0, and multi-update count 0.\n\n";
            for (const RasUpdateDiagnosticsRow& row : allRasRows) {
                mappingOut << "- " << row.solver
                           << ", subdomain " << row.subdomain
                           << ": core=" << row.coreDofCount
                           << ", local=" << row.localDofCount
                           << ", overlap_extra=" << row.overlapDofCount
                           << ", update[min/max/zero/multi]="
                           << row.minUpdateCount << '/'
                           << row.maxUpdateCount << '/'
                           << row.zeroUpdateCount << '/'
                           << row.multiUpdateCount << "\n";
            }
        }
        std::vector<SolverStatistics> comparisonStats;
        std::vector<double> maxDiffs;
        std::vector<double> l2Diffs;
        std::vector<double> relativeL2Diffs;
        const double referenceL2Norm = l2Norm(referenceTemperature);
        comparisonStats.reserve(results.size());
        maxDiffs.reserve(results.size());
        l2Diffs.reserve(results.size());
        relativeL2Diffs.reserve(results.size());
        for (SolverRunResult& result : results) {
            fillTemperatureStats(result.stats, result.temperature);
            comparisonStats.push_back(result.stats);
            if (result.temperature.empty()) {
                maxDiffs.push_back(std::numeric_limits<double>::quiet_NaN());
                l2Diffs.push_back(std::numeric_limits<double>::quiet_NaN());
                relativeL2Diffs.push_back(std::numeric_limits<double>::quiet_NaN());
            } else {
                maxDiffs.push_back(maxAbsDifference(result.temperature, referenceTemperature));
                const double l2Diff = l2Difference(result.temperature, referenceTemperature);
                l2Diffs.push_back(l2Diff);
                relativeL2Diffs.push_back(l2Diff / std::max(1.0e-300, referenceL2Norm));
            }
        }

        std::cout << "Solver comparison:\n";
        for (const SolverRunResult& result : results) {
            printSolverSummary(result.stats);
        }
        std::cout << "  final solution difference vs " << results[referenceIndex].stats.name << ":\n";
        for (size_t i = 0; i < results.size(); ++i) {
            std::cout << "    " << results[i].stats.name << ": max_abs=" << maxDiffs[i]
                      << " K, l2=" << l2Diffs[i]
                      << " K, relative_l2=" << relativeL2Diffs[i] << "\n";
        }
        writeMonolithicVsSchwarzComparison(results,
                                           physics,
                                           steady,
                                           system,
                                           outputDir / "comparison_monolithic_vs_schwarz.csv");
        writeOverlapDiagnostics(results,
                                maxDiffs,
                                relativeL2Diffs,
                                physics,
                                outputDir / "overlap_diagnostics.csv");
        writeRasPreconditionerTiming(comparisonStats,
                                     outputDir / "ras_preconditioner_timing.csv");
        writeCoarseDiagnostics(comparisonStats,
                               outputDir / "coarse_diagnostics.csv");
        writeLargeScaleSolverSummary(results,
                                     maxDiffs,
                                     relativeL2Diffs,
                                     mesh,
                                     physics,
                                     steady,
                                     outputDir / "large_scale_solver_summary.csv");
        SchwarzOptions rramDiagnosticOptions = physics.schwarz;
        rramDiagnosticOptions.standaloneMode = "algebraic";
        rramDiagnosticOptions.transmission = "none";
        SchwarzBlockSolver rramDiagnosticSolver(mesh,
                                                system,
                                                rramDiagnosticOptions,
                                                &physics,
                                                nullptr,
                                                false);
        writeRramDiagnosticsBundle(mesh,
                                   physics,
                                   rramDiagnosticSolver,
                                   referenceTemperature,
                                   outputDir);
        writeRramSolverDiagnosticsSummary(results,
                                          mesh,
                                          physics,
                                          rramDiagnosticSolver,
                                          outputDir / "rram_solver_diagnostics_summary.csv");
        if (!allKrylovRows.empty()) {
            writeKrylovIterationLog(allKrylovRows,
                                    outputDir / "krylov_iterations.csv");
        }

        const auto postStart = std::chrono::steady_clock::now();
        if (options.fastRun) {
            for (const SolverRunResult& result : results) {
                std::string fileStem = lowerString(result.stats.name);
                std::replace(fileStem.begin(), fileStem.end(), '-', '_');
                writeCsv(mesh, result.temperature, outputDir / ("temperature_" + fileStem + "_nodes.csv"));
            }
            writeSolverComparison(comparisonStats, maxDiffs, l2Diffs, relativeL2Diffs, outputDir / "solver_comparison.csv");
            writeIcDiagnostics(icDiagnostics, outputDir / "ic_factor_diagnostics.csv");
            writeInterfaceSummary(mesh, outputDir / "interface_summary.csv");
        } else {
            writeCsv(mesh, referenceTemperature, outputDir / "temperature_nodes.csv");
            writeVtk(mesh, referenceTemperature, outputDir / "temperature.vtk");
            for (const SolverRunResult& result : results) {
                std::string fileStem = lowerString(result.stats.name);
                std::replace(fileStem.begin(), fileStem.end(), '-', '_');
                writeCsv(mesh, result.temperature, outputDir / ("temperature_" + fileStem + "_nodes.csv"));
                writeCsv(mesh, result.temperature, outputDir / ("temperature_" + modeName + "_" + fileStem + "_nodes.csv"));
            }
            writeBoundarySummary(mesh, outputDir / "boundary_summary.csv");
            writeInterfaceSummary(mesh, outputDir / "interface_summary.csv");
            writeInterfaceBuildSummary(mesh, outputDir / "interface_build_summary.csv");
            writeMatrixDiagnostics(matrixDiagnostics, outputDir / "matrix_diagnostics.csv");
            writeMatrixStageDiagnostics(matrixStageDiagnostics, outputDir / "matrix_stage_diagnostics.csv");
            writeComponentDiagnostics(componentDiagnostics, outputDir / "component_diagnostics.csv");
            writeDiagonalContributionStats(assemblyDiagnostics, outputDir / "diagonal_contribution_stats.csv");
            writeNegativeDiagonalDofs(mesh, assemblyDiagnostics, outputDir / "negative_diagonal_dofs.csv");
            writeInterfacePenaltyStats(assemblyDiagnostics.interfacePenaltyStats, outputDir / "interface_penalty_stats.csv");
            writeHeatSourceDiagnostics(heatSourceDiagnostics, outputDir / "heat_source_diagnostics.csv");
            writeHeatSourceDomainDiagnostics(heatSourceDiagnostics, outputDir / "heat_source_domain_diagnostics.csv");
            writePhysicalSummary(mesh,
                                 physics,
                                 heatOnlySource,
                                 sourceBeforeDirichlet,
                                 sourceAfterDirichlet,
                                 convectionBoundaryFaceCount,
                                 outputDir / "physical_summary.csv");
            writeIcDiagnostics(icDiagnostics, outputDir / "ic_factor_diagnostics.csv");
            writeInterfaceTrianglesVtk(mesh, outputDir / "interface_integration_triangles.vtk");
            writeTemperatureSliceImage(mesh, referenceTemperature, outputDir / "temperature_slice_y_mid.bmp");
            writeSolverComparison(comparisonStats, maxDiffs, l2Diffs, relativeL2Diffs, outputDir / "solver_comparison.csv");
            writeSolverComparison(comparisonStats, maxDiffs, l2Diffs, relativeL2Diffs, outputDir / ("solver_comparison_" + modeName + ".csv"));
            if (options.runTwoLevelRasIlut || options.runDeflatedRasIlut) {
                writeSolverComparison(comparisonStats, maxDiffs, l2Diffs, relativeL2Diffs,
                                      outputDir / "rram_two_level_ras_comparison.csv");
            }
            if (options.runDeflatedRasIlut || options.runInterfaceDeflatedRasIlut) {
                writeSolverComparison(comparisonStats, maxDiffs, l2Diffs, relativeL2Diffs,
                                      outputDir / "rram_deflation_sweep.csv");
            }
        }
        if (!steady && !physics.comsolComparisonPath.empty()) {
            compareWithComsolExport(mesh, referenceTemperature,
                                    physics.comsolComparisonPath,
                                    outputDir / "comparison_comsol_t10.csv");
        } else {
            std::cout << "Skipped COMSOL comparison"
                      << (steady ? " for steady-state mode." : " because no comsol_comparison file is configured.")
                      << "\n";
        }
        timing.postprocessingSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - postStart).count();
        timing.totalSeconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - programStart).count();
        writeProgramTiming(timing, outputDir / "program_timing.csv");
        printProgramTiming(timing);
        std::cout << "Wrote solver outputs, summaries, images, and timing CSVs under "
                  << outputDir.string();
        if (!steady && !physics.comsolComparisonPath.empty()) {
            std::cout << " including comparison_comsol_t10.csv";
        }
        std::cout << "\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }
}
