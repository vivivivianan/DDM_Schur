#pragma once

#include "mor/types.hpp"

#include <cstdint>
#include <string>
#include <vector>

struct CaseConfig;
struct Mesh;

namespace mor::parametric {

struct AffineParameter {
    std::string name;
    std::string units;
    int subdomain = -1;
    int regionId = -1;
    double minimum = 0.0;
    double maximum = 0.0;
    double reference = 0.0;
    double ambientTemperature = 0.0;
    std::vector<double> harmonicNeighborConductivities;
    bool touchesInterface = false;
    int selectedTetCount = 0;
    int selectedBoundaryFaceCount = 0;
    std::uint64_t definitionHash = 0;
};

AffineParameter resolveAffineParameter(const Mesh& mesh,
                                       const CaseConfig& physics,
                                       const mor::Options& options);

CaseConfig physicsAtParameter(const CaseConfig& physics,
                              const AffineParameter& parameter,
                              double value);

bool tetSelected(const AffineParameter& parameter,
                 int subdomain,
                 int domainEntity);

bool boundarySelected(const AffineParameter& parameter,
                      int subdomain,
                      int boundaryEntity);

double harmonicTheta(const AffineParameter& parameter,
                     double value,
                     std::size_t group);

} // namespace mor::parametric
