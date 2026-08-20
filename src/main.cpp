#include "sipg_core.hpp"
#include "config_io.hpp"
#include "mesh_loader.hpp"
#include "fem_assembly.hpp"
#include "linear_solvers.hpp"
#include "diagnostics_io.hpp"
#include "schwarz_solver.hpp"
#include "ddm_schur/schur_fgmres.hpp"
#include "application.hpp"

int main(int argc, char* argv[])
{
    return runProgram(argc, argv);
}
