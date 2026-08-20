#include "sipg_core.hpp"
#include "config_io.hpp"
#include "mesh_loader.hpp"
#include "fem_assembly.hpp"
#include "linear_solvers.hpp"
#include "auxiliary_space_solver.hpp"
#include "diagnostics_io.hpp"
#include "schwarz_solver.hpp"
#include "ddm_schur/schur_fgmres.hpp"
#include "ddm_schur/schur_direct_exact.hpp"
#include "mor/reduced_schur_workflow.hpp"
#include "mor/deployment_response_model.hpp"
#include "mor/parametric/parametric_workflow.hpp"
#include "mor/transient/transient_workflow.hpp"
#include "mor/transient/local_dynamic_schur.hpp"
#include "mor/transient/local_interior_pod_test.hpp"
#include "mor/local/local_rom_solver.hpp"
#include "application.hpp"

int main(int argc, char* argv[])
{
    return runProgram(argc, argv);
}
