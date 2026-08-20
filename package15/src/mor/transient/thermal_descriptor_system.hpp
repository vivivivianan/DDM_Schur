#pragma once

#include "sipg_core.hpp"
#include "mor/types.hpp"

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace mor::transient {

struct TransientFingerprints {
    std::uint64_t mesh = 0;
    std::uint64_t capacity = 0;
    std::uint64_t conductivity = 0;
    std::uint64_t input = 0;
    std::uint64_t boundary = 0;
    std::uint64_t sources = 0;
};

struct MatrixDiagnostic {
    std::string name;
    std::string units;
    int rows = 0;
    std::size_t nonzeros = 0;
    double symmetryError = 0.0;
    double minimumDiagonal = 0.0;
    double maximumDiagonal = 0.0;
    int zeroDiagonalCount = 0;
    bool positiveOnFreeDofs = false;
};

struct ThermalDescriptorSystem {
    int dofs = 0;
    int sourceChannels = 0;
    std::string massType = "consistent";
    SparseMatrix capacity;
    SparseMatrix conductivity;
    // Lumped positive trace Gram diagonals assembled on the exact SIPG
    // nonmatching interface quadrature.
    std::vector<double> interfaceTraceMassDiagonal;
    std::vector<double> interfacePenaltyMassDiagonal;
    // Column-major B: input[channel * dofs + row] is the load per watt.
    std::vector<double> input;
    std::vector<double> boundaryRhs;
    std::vector<double> nominalPowersW;
    std::vector<double> minimumPowersW;
    std::vector<double> maximumPowersW;
    std::vector<int> sourceSubdomains;
    std::vector<int> sourceDomainEntities;
    std::vector<DeploymentDof> deploymentDofs;
    TransientFingerprints fingerprints;
    MatrixDiagnostic capacityDiagnostic;
    MatrixDiagnostic conductivityDiagnostic;
    double assemblySeconds = 0.0;
};

ThermalDescriptorSystem assembleThermalDescriptorSystem(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::string& massType = "consistent",
    bool collectInterfaceGram = false);

// Fingerprint every mesh/configuration input consumed by descriptor assembly.
// This is deliberately independent of the assembled matrix fingerprints so a
// persistent descriptor cache can be validated before paying assembly cost.
std::uint64_t thermalDescriptorInputFingerprint(
    const Mesh& mesh,
    const CaseConfig& physics,
    const std::string& massType,
    bool collectInterfaceGram);

void saveThermalDescriptorCache(
    const std::filesystem::path& path,
    std::uint64_t inputFingerprint,
    const ThermalDescriptorSystem& system);

ThermalDescriptorSystem loadThermalDescriptorCache(
    const std::filesystem::path& path,
    std::uint64_t expectedInputFingerprint,
    int expectedDofs);

std::vector<double> descriptorInputRhs(
    const ThermalDescriptorSystem& system,
    const std::vector<double>& powersW);

void writeTransientMatrixDiagnostics(
    const ThermalDescriptorSystem& system,
    const std::filesystem::path& path);

bool sameFingerprints(const TransientFingerprints& left,
                      const TransientFingerprints& right);

} // namespace mor::transient
