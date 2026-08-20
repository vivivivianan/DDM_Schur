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
                throw std::runtime_error("--solver-method must be monolithic, schwarz, schwarz-precond-fgmres, or schwarz-precond-fgmres-two-level.");
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
        } else {
            physics.penaltyMode = "harmonic";
        }
        if (options.penaltyFactorOverride) {
            physics.penaltyFactor = options.penaltyFactor;
        } else {
            physics.penaltyFactor = 15.0;
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
            // A zero source is physically valid when temperature/flux/Robin
            // boundaries drive the problem. Matrix constraints are checked
            // independently below, so source magnitude is not a solvability
            // criterion.
            std::cout << "Configured heat source total is zero; continuing with boundary-driven problem.\n";
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
        const int heatFluxBoundaryFaceCount = assembleHeatFluxBoundaries(mesh, physics, source);
        std::cout << "Matched inward heat-flux boundary faces: " << heatFluxBoundaryFaceCount << "\n";
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
            const double transientTheta = physics.timeIntegrator == "crank_nicolson" ? 0.5 : 1.0;
            buildSystemMatrix(mass, stiffness, physics.timeStep, transientTheta, system);
            buildSystemMatrix(mass, physicalTransmissionStiffness, physics.timeStep, transientTheta, physicalTransmissionSystem);
            stiffness = SparseMatrix(n);
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
        if (options.runBjIlut || options.runRasIlut || options.runTwoLevelRasIlut
            || options.runDeflatedRasIlut || options.runInterfaceDeflatedRasIlut
            || options.runBjPardisoGeneral || options.runSchwarzPrecondFgmres
            || options.runRasIc || options.runBjIcCoarse || options.runBjIlutCoarse
            || options.runSchwarzPrecondFgmresTwoLevel || options.runDdmSchur) {
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

        if (options.runDdmSchur) {
            SolverStatistics stats;
            stats.name = "DDM-Schur-FGMRES";
            stats.solverMethod = "matrix_free_schur";
            const std::string localCoarseModes = options.schurLinearXYCoarse
                ? (options.schurLinearZCoarse ? "constant_plus_linear_xyz" : "constant_plus_linear_xy")
                : (options.schurLinearZCoarse ? "constant_plus_linear_z" : "constant");
            stats.preconditioner = options.schurInterfacePatchCoarse
                ? std::string("balanced_two_level_schur_interface_patch_constant")
                    + (options.schurInterfacePatchLinearXY ? "_plus_linear_xy" : "")
                : "balanced_two_level_schur_subdomain_" + localCoarseModes
                    + (options.schurGlobalQuadraticZCoarse ? "_plus_global_quadratic_z" : "");
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
                ddm_schur::DdmSchurSolver solver(mesh, system, schurOptions);
                stats.setupSeconds = solver.setupReport().setupSeconds;
                stats.preconditionerBytes = solver.setupReport().memoryBytes;
                stats.coarseEnabled = true;
                stats.coarseDim = solver.setupReport().coarseDimension;
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

                std::ofstream summary(outputDir / "schur_solver_summary.csv");
                summary << "domains,total_dofs,interface_dofs,interior_dofs,interface_patch_count,coarse_dimension,iterations,schur_matvecs,"
                        << "local_solve_calls,local_symbolic_analysis_calls,local_numerical_factorization_calls,"
                        << "setup_seconds,local_factorization_seconds,local_symbolic_analysis_seconds,"
                        << "local_numerical_factorization_seconds,local_solve_seconds,coarse_solve_seconds,"
                        << "condensed_rhs_seconds,interface_solve_seconds,fgmres_seconds,recovery_seconds,"
                        << "total_solve_seconds,total_seconds,relative_residual,memory_bytes,status\n";
                summary << std::setprecision(16)
                        << lastReport.domains << ','
                        << lastReport.totalDofs << ','
                        << lastReport.interfaceDofs << ','
                        << lastReport.interiorDofs << ','
                        << lastReport.interfacePatchCount << ','
                        << lastReport.coarseDimension << ','
                        << lastReport.iterations << ','
                        << lastReport.schurMatvecs << ','
                        << lastReport.localSolveCalls << ','
                        << lastReport.localSymbolicAnalysisCalls << ','
                        << lastReport.localNumericalFactorizationCalls << ','
                        << lastReport.setupSeconds << ','
                        << lastReport.localFactorizationSeconds << ','
                        << lastReport.localSymbolicAnalysisSeconds << ','
                        << lastReport.localNumericalFactorizationSeconds << ','
                        << lastReport.localSolveSeconds << ','
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
            throw std::runtime_error("No solvers were selected. Use --solvers schur,direct or another supported solver.");
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
