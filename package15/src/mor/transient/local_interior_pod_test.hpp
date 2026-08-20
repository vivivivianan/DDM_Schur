#pragma once

#include <filesystem>
#include <vector>

struct CaseConfig;
struct Mesh;
struct ProgramOptions;
struct SparseMatrix;

int runTransientLocalInteriorRomTest(
    const Mesh& mesh,
    const CaseConfig& physics,
    const SparseMatrix& mass,
    const SparseMatrix& stiffness,
    const SparseMatrix& system,
    const std::vector<double>& source,
    const std::vector<double>& fixedAdjust,
    const ProgramOptions& options,
    const std::filesystem::path& outputDir);
