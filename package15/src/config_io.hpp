#pragma once

#include <cstdlib>

// Case configuration parsing and material lookup helpers.
// This file is intentionally included from main.cpp after the preceding SIPG modules.

static std::string trim(const std::string& s)
{
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

static std::string numericPart(const std::string& s)
{
    const auto hash = s.find('#');
    return trim(hash == std::string::npos ? s : s.substr(0, hash));
}

static std::vector<std::string> readLines(const std::filesystem::path& path)
{
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open " + path.string());
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)) {
        lines.push_back(line);
    }
    return lines;
}

static int firstInt(const std::string& s)
{
    std::istringstream iss(numericPart(s));
    int value = 0;
    iss >> value;
    return value;
}

static std::vector<int> parseInts(const std::string& s)
{
    std::istringstream iss(numericPart(s));
    std::vector<int> values;
    int value = 0;
    while (iss >> value) {
        values.push_back(value);
    }
    return values;
}

static std::vector<double> parseDoubles(const std::string& s)
{
    std::istringstream iss(numericPart(s));
    std::vector<double> values;
    double value = 0.0;
    while (iss >> value) {
        values.push_back(value);
    }
    return values;
}

static std::vector<std::string> splitCsv(const std::string& value)
{
    std::vector<std::string> items;
    std::string item;
    std::istringstream iss(value);
    while (std::getline(iss, item, ',')) {
        items.push_back(trim(item));
    }
    return items;
}

static std::string lowerString(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

static bool parseBoolValue(const std::string& value)
{
    const std::string v = lowerString(trim(value));
    return v == "1" || v == "true" || v == "yes" || v == "on";
}

static int parseSelectorInt(const std::string& value)
{
    const std::string v = lowerString(trim(value));
    if (v == "*" || v == "all" || v == "any") {
        return -1;
    }
    return std::stoi(v);
}

static std::vector<int> parseEntityList(const std::string& value)
{
    std::vector<int> entities;
    std::string item;
    std::istringstream semicolonStream(value);
    while (std::getline(semicolonStream, item, ';')) {
        const std::string trimmed = trim(item);
        if (!trimmed.empty()) {
            entities.push_back(std::stoi(trimmed));
        }
    }
    if (entities.empty()) {
        throw std::runtime_error("Entity list cannot be empty.");
    }
    return entities;
}

static bool containsEntity(const std::vector<int>& entities, int entity)
{
    return std::find(entities.begin(), entities.end(), entity) != entities.end();
}

static std::string expandCaseEnvironment(const std::string& input)
{
    std::string output;
    output.reserve(input.size());
    size_t position = 0;
    while (position < input.size()) {
        const size_t begin = input.find("${", position);
        if (begin == std::string::npos) {
            output.append(input, position, std::string::npos);
            break;
        }
        output.append(input, position, begin - position);
        const size_t end = input.find('}', begin + 2);
        if (end == std::string::npos) {
            throw std::runtime_error(
                "Unterminated environment variable in config path: " + input);
        }
        const std::string name = input.substr(begin + 2, end - begin - 2);
        std::string value;
#ifdef _WIN32
        char* buffer = nullptr;
        size_t size = 0;
        if (_dupenv_s(&buffer, &size, name.c_str()) == 0 && buffer != nullptr) {
            value = buffer;
        }
        std::free(buffer);
#else
        const char* buffer = std::getenv(name.c_str());
        if (buffer != nullptr) {
            value = buffer;
        }
#endif
        if (name.empty() || value.empty()) {
            throw std::runtime_error(
                "Required environment variable '" + name + "' is not set.");
        }
        output += value;
        position = end + 1;
    }
    return output;
}

static std::filesystem::path resolveCasePath(const std::filesystem::path& baseDir,
                                             const std::string& value)
{
    std::filesystem::path p(expandCaseEnvironment(trim(value)));
    if (p.is_relative()) {
        p = baseDir / p;
    }
    return p.lexically_normal();
}

static CaseConfig readCaseConfig(const std::filesystem::path& path)
{
    const std::vector<std::string> lines = readLines(path);
    const std::filesystem::path baseDir = path.parent_path();
    CaseConfig config;

    for (const std::string& rawLine : lines) {
        const std::string line = trim(numericPart(rawLine));
        if (line.empty()) {
            continue;
        }
        const auto eq = line.find('=');
        if (eq == std::string::npos) {
            continue;
        }

        const std::string key = lowerString(trim(line.substr(0, eq)));
        const std::string value = trim(line.substr(eq + 1));
        if (key == "name") {
            config.name = value;
        } else if (key == "method" || key == "solver_method") {
            config.solverMethod = normalizeSolverMethodName(value);
            if (!isKnownSolverMethod(config.solverMethod)) {
                throw std::runtime_error("solver method must be monolithic, schwarz, schwarz_precond_fgmres, or schwarz_precond_fgmres_two_level.");
            }
            if (config.solverMethod == "schwarz"
                || config.solverMethod == "schwarz_precond_fgmres"
                || config.solverMethod == "schwarz_precond_fgmres_two_level") {
                config.schwarz.enabled = true;
            }
        } else if (key == "time_integrator") {
            config.timeIntegrator = lowerString(value);
            if (config.timeIntegrator != "backward_euler"
                && config.timeIntegrator != "crank_nicolson") {
                throw std::runtime_error(
                    "time_integrator must be backward_euler or crank_nicolson.");
            }
        } else if (key == "schwarz_enabled" || key == "enabled") {
            config.schwarz.enabled = parseBoolValue(value);
            if (config.schwarz.enabled && config.solverMethod == "monolithic") {
                config.solverMethod = "schwarz";
            }
        } else if (key == "schwarz_type" || key == "type") {
            config.schwarz.type = lowerString(value);
            if (config.schwarz.type != "multiplicative" && config.schwarz.type != "additive") {
                throw std::runtime_error("Schwarz type must be multiplicative or additive.");
            }
        } else if (key == "schwarz_standalone_mode" || key == "standalone_mode") {
            config.schwarz.standaloneMode = normalizeSchwarzTransmissionModeName(value);
            if (!isKnownSchwarzStandaloneMode(config.schwarz.standaloneMode)) {
                throw std::runtime_error(
                    "Schwarz standalone_mode must be algebraic, dirichlet_neumann, dirichlet_dirichlet, dirichlet_robin, or robin.");
            }
            config.schwarz.transmission =
                schwarzTransmissionForStandaloneMode(config.schwarz.standaloneMode);
        } else if (key == "schwarz_transmission" || key == "transmission") {
            config.schwarz.transmission = normalizeSchwarzTransmissionModeName(value);
            if (config.schwarz.transmission == "algebraic") {
                config.schwarz.transmission = "none";
            }
            if (!isKnownSchwarzTransmissionMode(config.schwarz.transmission)) {
                throw std::runtime_error(
                    "Schwarz transmission must be none, dirichlet_neumann, dirichlet_dirichlet, dirichlet_robin, or robin.");
            }
            if (config.schwarz.transmission != "none") {
                config.schwarz.standaloneMode = config.schwarz.transmission;
            }
        } else if (key == "schwarz_transmission_orientation" || key == "transmission_orientation") {
            config.schwarz.transmissionOrientation =
                normalizeSchwarzTransmissionOrientationName(value);
            if (!isKnownSchwarzTransmissionOrientation(config.schwarz.transmissionOrientation)) {
                throw std::runtime_error("Schwarz transmission_orientation must be forward or reverse.");
            }
        } else if (key == "schwarz_robin_alpha_factor" || key == "robin_alpha_factor") {
            config.schwarz.robinAlphaFactor = std::stod(value);
        } else if (key == "schwarz_flux_eval" || key == "flux_eval") {
            config.schwarz.fluxEval = lowerString(value);
            if (config.schwarz.fluxEval != "physical_gradient"
                && config.schwarz.fluxEval != "sipg_numeric") {
                throw std::runtime_error("Schwarz flux_eval must be physical_gradient or sipg_numeric.");
            }
        } else if (key == "schwarz_overlap_mode" || key == "overlap_mode") {
            config.schwarz.overlapMode = lowerString(value);
            if (config.schwarz.overlapMode != "none"
                && config.schwarz.overlapMode != "halo"
                && config.schwarz.overlapMode != "ras") {
                throw std::runtime_error("Schwarz overlap_mode must be none, halo, or ras.");
            }
        } else if (key == "schwarz_overlap_layers" || key == "overlap_layers") {
            config.schwarz.overlapLayers = std::stoi(value);
            if (config.schwarz.overlapLayers < 0 || config.schwarz.overlapLayers > 3) {
                throw std::runtime_error("Schwarz overlap_layers must be 0, 1, 2, or 3.");
            }
        } else if (key == "schwarz_partition_mode" || key == "partition_mode") {
            config.schwarz.partitionMode = normalizeSchwarzPartitionModeName(value);
            if (!isKnownSchwarzPartitionMode(config.schwarz.partitionMode)) {
                throw std::runtime_error(
                    "Schwarz partition_mode must be current, material_aligned, hotspot_contained, or vertical_heat_flow_aligned.");
            }
        } else if (key == "schwarz_write_interface_flux" || key == "write_interface_flux") {
            config.schwarz.writeInterfaceFlux = parseBoolValue(value);
        } else if (key == "schwarz_max_iters" || key == "max_iters") {
            config.schwarz.maxIters = std::stoi(value);
        } else if (key == "schwarz_tol_rel_update" || key == "tol_rel_update") {
            config.schwarz.tolRelUpdate = std::stod(value);
        } else if (key == "schwarz_tol_rel_residual" || key == "tol_rel_residual") {
            config.schwarz.tolRelResidual = std::stod(value);
        } else if (key == "schwarz_relaxation" || key == "relaxation") {
            config.schwarz.relaxation = std::stod(value);
        } else if (key == "schwarz_initial_guess" || key == "initial_guess") {
            config.schwarz.initialGuess = lowerString(value);
            if (config.schwarz.initialGuess != "previous_time_step") {
                throw std::runtime_error("Schwarz initial_guess currently supports previous_time_step only.");
            }
        } else if (key == "schwarz_check_interface_jump" || key == "check_interface_jump") {
            config.schwarz.checkInterfaceJump = parseBoolValue(value);
        } else if (key == "schwarz_output_iteration_log" || key == "output_iteration_log") {
            config.schwarz.outputIterationLog = parseBoolValue(value);
        } else if (key == "schwarz_validate_against_monolithic" || key == "validate_against_monolithic") {
            config.schwarz.validateAgainstMonolithic = parseBoolValue(value);
        } else if (key == "coordinate_scale") {
            config.coordinateScale = std::stod(value);
        } else if (key == "initial_temperature") {
            config.initialTemperature = std::stod(value);
        } else if (key == "start_time") {
            config.startTime = std::stod(value);
        } else if (key == "time_step") {
            config.timeStep = std::stod(value);
        } else if (key == "time_steps") {
            config.timeSteps = std::stoi(value);
        } else if (key == "penalty_factor") {
            config.penaltyFactor = std::stod(value);
        } else if (key == "dirichlet_method") {
            config.dirichletMethod = lowerString(value);
            if (config.dirichletMethod != "strong" && config.dirichletMethod != "nitsche") {
                throw std::runtime_error("dirichlet_method must be strong or nitsche.");
            }
        } else if (key == "nitsche_penalty_factor") {
            config.nitschePenaltyFactor = std::stod(value);
        } else if (key == "penalty_mode") {
            config.penaltyMode = lowerString(value);
            if (config.penaltyMode != "harmonic" && config.penaltyMode != "max") {
                throw std::runtime_error("penalty_mode must be harmonic or max.");
            }
        } else if (key == "interface_scheme") {
            config.interfaceScheme = lowerString(value);
            if (config.interfaceScheme != "sipg" && config.interfaceScheme != "nipg") {
                throw std::runtime_error("interface_scheme must be sipg or nipg.");
            }
        } else if (key == "penalty_scaling") {
            config.penaltyScaling = lowerString(value);
            if (config.penaltyScaling == "p(p+1)" || config.penaltyScaling == "p_p1" || config.penaltyScaling == "p-p1") {
                config.penaltyScaling = "p_p1";
            } else if (config.penaltyScaling == "(p+1)^2" || config.penaltyScaling == "p1_squared" || config.penaltyScaling == "pplus1_squared") {
                config.penaltyScaling = "p1_squared";
            } else {
                throw std::runtime_error("penalty_scaling must be p_p1 or p1_squared.");
            }
        } else if (key == "auto_interfaces") {
            config.autoInterfaces = parseBoolValue(value);
        } else if (key == "output_dir") {
            config.outputDir = resolveCasePath(baseDir, value);
        } else if (key == "comsol_comparison") {
            config.comsolComparisonPath = resolveCasePath(baseDir, value);
        } else if (key == "domain") {
            const std::vector<std::string> items = splitCsv(value);
            if (items.empty()) {
                throw std::runtime_error("domain entry requires at least a mesh path.");
            }
            DomainConfig domain;
            domain.meshPath = resolveCasePath(baseDir, items[0]);
            if (items.size() >= 5) {
                domain.material.name = items[1];
                setIsotropicConductivity(domain.material, std::stod(items[2]));
                domain.material.density = std::stod(items[3]);
                domain.material.heatCapacity = std::stod(items[4]);
            }
            config.domains.push_back(domain);
        }
        else if (key == "domain_offset_m") {
            const std::vector<std::string> items = splitCsv(value);
            if (items.size() < 4) {
                throw std::runtime_error("domain_offset_m entry format: subdomain,dx_m,dy_m,dz_m.");
            }
            const int subdomain = std::stoi(items[0]);
            if (subdomain < 0 || subdomain >= static_cast<int>(config.domains.size())) {
                throw std::runtime_error("domain_offset_m references a subdomain that has not been defined.");
            }
            config.domains[static_cast<std::size_t>(subdomain)].translationMeters =
                {std::stod(items[1]), std::stod(items[2]), std::stod(items[3])};
        }
        else if (key == "domain_material") {
            const std::vector<std::string> items = splitCsv(value);
            if (items.size() < 8) {
                throw std::runtime_error(
                    "domain_material entry format: subdomain,domain_entity,material_name,kx,ky,kz,rho,cp. "
                    "subdomain can be a concrete index or * / all / any.");
            }

            const int subdomainSelector = parseSelectorInt(items[0]);
            const int domainEntity = std::stoi(items[1]);

            Material material;
            material.name = items[2];
            material.conductivityX = std::stod(items[3]);
            material.conductivityY = std::stod(items[4]);
            material.conductivityZ = std::stod(items[5]);
            material.conductivity =
                (material.conductivityX + material.conductivityY + material.conductivityZ) / 3.0;
            material.density = std::stod(items[6]);
            material.heatCapacity = std::stod(items[7]);

            const auto assignMaterialToSubdomain = [&](int subdomain) {
                config.domains[static_cast<size_t>(subdomain)]
                    .materialsByDomainEntity[domainEntity] = material;
            };

            if (subdomainSelector < 0) {
                if (config.domains.empty()) {
                    throw std::runtime_error(
                        "domain_material uses * but no domain has been defined yet. "
                        "Put all domain = ... lines before domain_material = ... lines.");
                }

                for (int subdomain = 0; subdomain < static_cast<int>(config.domains.size()); ++subdomain) {
                    assignMaterialToSubdomain(subdomain);
                }
            }
            else {
                if (subdomainSelector >= static_cast<int>(config.domains.size())) {
                    throw std::runtime_error(
                        "domain_material references a subdomain index that has not been defined.");
                }

                assignMaterialToSubdomain(subdomainSelector);
            }
        } else if (key == "dirichlet") {
            const std::vector<std::string> items = splitCsv(value);
            if (items.size() < 3) {
                throw std::runtime_error("dirichlet entry format: subdomain,boundary_entity,temperature.");
            }
            config.dirichletConditions.push_back(
                {parseSelectorInt(items[0]), std::stoi(items[1]), std::stod(items[2])});
        } else if (key == "convection" || key == "robin") {
            const std::vector<std::string> items = splitCsv(value);
            if (items.size() < 4) {
                throw std::runtime_error("convection entry format: subdomain,boundary_entity,h_W_m2K,ambient_temperature.");
            }
            config.convectionConditions.push_back(
                {parseSelectorInt(items[0]), std::stoi(items[1]), std::stod(items[2]), std::stod(items[3])});
        } else if (key == "heat_source") {
            const std::vector<std::string> items = splitCsv(value);
            if (items.size() < 3) {
                throw std::runtime_error("heat_source entry format: subdomain,domain_entity,total_W.");
            }
            config.heatSources.push_back(
                {parseSelectorInt(items[0]), std::stoi(items[1]), std::stod(items[2])});
        } else if (key == "interface") {
            const std::vector<std::string> items = splitCsv(value);
            if (items.size() < 3) {
                throw std::runtime_error(
                    "interface entry format: left_subdomain,right_subdomain,boundary_entity "
                    "or left_subdomain,right_subdomain,left_boundary_entities,right_boundary_entities.");
            }
            InterfaceConfig interfaceConfig;
            interfaceConfig.leftSubdomain = std::stoi(items[0]);
            interfaceConfig.rightSubdomain = std::stoi(items[1]);
            interfaceConfig.leftBoundaryEntities = parseEntityList(items[2]);
            interfaceConfig.rightBoundaryEntities = items.size() >= 4
                ? parseEntityList(items[3])
                : interfaceConfig.leftBoundaryEntities;
            config.interfaces.push_back(std::move(interfaceConfig));
        } else if (key == "interface_face_pair") {
            const std::vector<std::string> items = splitCsv(value);
            if (items.size() < 2) {
                throw std::runtime_error("interface_face_pair entry format: left_boundary_face_id,right_boundary_face_id.");
            }
            config.explicitInterfaceFacePairs.push_back({std::stoi(items[0]), std::stoi(items[1])});
        }
    }

    if (config.domains.empty()) {
        throw std::runtime_error("Case config has no domains. Add at least one domain = meshfile.mphtxt,... line.");
    }
    return config;
}

