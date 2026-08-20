cmake_minimum_required(VERSION 3.15)

project(SIPGHeatDDM3D LANGUAGES CXX)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(SIPGHeatDDM3D
    src/application.hpp
    src/auxiliary_space_solver.hpp
    src/config_io.hpp
    src/diagnostics_io.hpp
    src/fem_assembly.hpp
    src/linear_solvers.hpp
    src/main.cpp
    src/mesh_loader.hpp
    src/schwarz_solver.hpp
    src/sipg_core.hpp
    src/ddm_schur/common.hpp
    src/ddm_schur/interface_operator.cpp
    src/ddm_schur/interface_operator.hpp
    src/ddm_schur/local_solver.cpp
    src/ddm_schur/local_solver.hpp
    src/ddm_schur/schur_fgmres.cpp
    src/ddm_schur/schur_fgmres.hpp
    src/ddm_schur/schur_operator.cpp
    src/ddm_schur/schur_operator.hpp
    src/ddm_schur/schur_proxy.cpp
    src/ddm_schur/schur_proxy.hpp
    src/ddm_schur/schur_direct_exact.cpp
    src/ddm_schur/schur_direct_exact.hpp
    src/mor/model_io.cpp
    src/mor/model_io.hpp
    src/mor/deployment_response_model.cpp
    src/mor/deployment_response_model.hpp
    src/mor/interior_model_serialization.cpp
    src/mor/interior_model_serialization.hpp
    src/mor/local_interior_basis.cpp
    src/mor/local_interior_basis.hpp
    src/mor/local/local_subdomain_model.cpp
    src/mor/local/local_subdomain_model.hpp
    src/mor/local/local_snapshot_basis.cpp
    src/mor/local/local_snapshot_basis.hpp
    src/mor/local/local_reduced_blocks.cpp
    src/mor/local/local_reduced_blocks.hpp
    src/mor/local/local_reduced_schur.cpp
    src/mor/local/local_reduced_schur.hpp
    src/mor/local/port_core_solver.cpp
    src/mor/local/port_core_solver.hpp
    src/mor/local/local_rom_solver.cpp
    src/mor/local/local_rom_solver.hpp
    src/mor/local/local_rom_serialization.cpp
    src/mor/local/local_rom_serialization.hpp
    src/mor/mor_diagnostics.cpp
    src/mor/mor_diagnostics.hpp
    src/mor/pod_basis.cpp
    src/mor/pod_basis.hpp
    src/mor/reduced_schur_model.cpp
    src/mor/reduced_schur_model.hpp
    src/mor/reduced_schur_workflow.cpp
    src/mor/reduced_schur_workflow.hpp
    src/mor/source_parameterization.cpp
    src/mor/source_parameterization.hpp
    src/mor/types.hpp
    src/mor/parametric/affine_parameter.cpp
    src/mor/parametric/affine_parameter.hpp
    src/mor/parametric/affine_fem_components.cpp
    src/mor/parametric/affine_fem_components.hpp
    src/mor/parametric/parametric_reduced_schur.cpp
    src/mor/parametric/parametric_reduced_schur.hpp
    src/mor/parametric/parametric_workflow.cpp
    src/mor/parametric/parametric_workflow.hpp
    src/mor/transient/thermal_descriptor_system.cpp
    src/mor/transient/thermal_descriptor_system.hpp
    src/mor/transient/block_arnoldi.cpp
    src/mor/transient/block_arnoldi.hpp
    src/mor/transient/transient_reduced_model.cpp
    src/mor/transient/transient_reduced_model.hpp
    src/mor/transient/transient_time_integrator.cpp
    src/mor/transient/transient_time_integrator.hpp
    src/mor/transient/transient_rom_serialization.cpp
    src/mor/transient/transient_rom_serialization.hpp
    src/mor/transient/transient_rom_diagnostics.cpp
    src/mor/transient/transient_rom_diagnostics.hpp
    src/mor/transient/transient_input_waveform.cpp
    src/mor/transient/transient_input_waveform.hpp
    src/mor/transient/transient_workflow.cpp
    src/mor/transient/transient_workflow.hpp
    src/mor/transient/local_dynamic_schur.cpp
    src/mor/transient/local_dynamic_schur.hpp
    src/mor/transient/local_interior_pod_test.cpp
    src/mor/transient/local_interior_pod_test.hpp
    src/mor/transient/local_port_reduced_schur.cpp
    src/mor/transient/local_port_reduced_schur.hpp
    src/mor/transient/global_interface_coarse.cpp
    src/mor/transient/global_interface_coarse.hpp
    src/mor/transient/global_randomized_schur.cpp
    src/mor/transient/global_randomized_schur.hpp
    src/mor/transient/projection_diagnosis.cpp
    src/mor/transient/projection_diagnosis.hpp
    src/mor/transient/interface_flux_operator.cpp
    src/mor/transient/interface_flux_operator.hpp
    src/mor/transient/optimal_port_space.cpp
    src/mor/transient/optimal_port_space.hpp
    src/mor/transient/transfer_operator.cpp
    src/mor/transient/transfer_operator.hpp
    src/mor/transient/port_eigensolver.cpp
    src/mor/transient/port_eigensolver.hpp
    src/mor/transient/port/IPortBasisBuilder.hpp
    src/mor/transient/port/randomized_transfer_port.cpp
    src/mor/transient/port/randomized_transfer_port.hpp
    src/mor/transient/port/port_basis_factory.cpp
    src/mor/transient/port/port_basis_factory.hpp
)
target_include_directories(SIPGHeatDDM3D PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
execute_process(
    COMMAND git rev-parse --short=12 HEAD
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE DDM_SCHUR_SOURCE_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
if (NOT DDM_SCHUR_SOURCE_COMMIT)
    set(DDM_SCHUR_SOURCE_COMMIT "unknown")
endif()
target_compile_definitions(SIPGHeatDDM3D PRIVATE
    DDM_SCHUR_SOURCE_COMMIT="${DDM_SCHUR_SOURCE_COMMIT}")

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(SIPGHeatDDM3D PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP enabled for matrix/vector kernels")
else()
    message(WARNING "OpenMP was not found. Matrix/vector kernels will run with serial fallback loops.")
endif()

find_path(MKL_INCLUDE_DIR
    NAMES mkl.h
    PATHS
        "$ENV{MKLROOT}/include"
)

find_library(MKL_RT_LIBRARY
    NAMES mkl_rt
    PATHS
        "$ENV{MKLROOT}/lib"
)

find_file(MKL_RT_DLL
    NAMES mkl_rt.2.dll mkl_rt.dll
    PATHS
        "$ENV{MKLROOT}/bin"
)

find_file(MKL_CORE_DLL
    NAMES mkl_core.2.dll mkl_core.dll
    PATHS
        "$ENV{MKLROOT}/bin"
)

find_file(MKL_THREAD_DLL
    NAMES mkl_intel_thread.2.dll mkl_intel_thread.dll
    PATHS
        "$ENV{MKLROOT}/bin"
)

find_file(MKL_AVX2_DLL
    NAMES mkl_avx2.2.dll mkl_avx2.dll
    PATHS
        "$ENV{MKLROOT}/bin"
)

find_file(MKL_DEF_DLL
    NAMES mkl_def.2.dll mkl_def.dll
    PATHS
        "$ENV{MKLROOT}/bin"
)

find_file(INTEL_OPENMP_DLL
    NAMES libiomp5md.dll
    PATHS
        "$ENV{CMPLR_ROOT}/bin"
)

if (MKL_INCLUDE_DIR AND MKL_RT_LIBRARY)
    target_compile_definitions(SIPGHeatDDM3D PRIVATE USE_MKL_PARDISO)
    target_include_directories(SIPGHeatDDM3D PRIVATE "${MKL_INCLUDE_DIR}")
    target_link_libraries(SIPGHeatDDM3D PRIVATE "${MKL_RT_LIBRARY}")
    foreach(required_dll IN ITEMS
        "${MKL_RT_DLL}" "${MKL_CORE_DLL}" "${MKL_THREAD_DLL}"
        "${MKL_AVX2_DLL}" "${MKL_DEF_DLL}" "${INTEL_OPENMP_DLL}")
        if (required_dll)
            add_custom_command(TARGET SIPGHeatDDM3D POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                        "${required_dll}"
                        "$<TARGET_FILE_DIR:SIPGHeatDDM3D>"
            )
        endif()
    endforeach()
    message(STATUS "MKL PARDISO enabled: ${MKL_RT_LIBRARY}")
else()
    message(WARNING "MKL was not found. BJ solver will use the built-in dense Cholesky fallback instead of PARDISO.")
endif()

if (MSVC)
    target_compile_options(SIPGHeatDDM3D PRIVATE /W4 /permissive- /utf-8)
    target_link_libraries(SIPGHeatDDM3D PRIVATE Psapi)
else()
    target_compile_options(SIPGHeatDDM3D PRIVATE -Wall -Wextra -pedantic)
endif()

if (BUILD_TESTING AND WIN32)
    find_program(POWERSHELL_EXECUTABLE NAMES pwsh powershell REQUIRED)
    set(STAGE2B1_TEST_ROOT "${CMAKE_BINARY_DIR}/stage2b1_ctest")
    set(STAGE2B1_TEST_DRIVER "${CMAKE_SOURCE_DIR}/tests/run_stage2b1_test.ps1")

    add_test(NAME stage2b1_affine_decomposition_two_cube
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${STAGE2B1_TEST_DRIVER}" -Mode Affine
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${STAGE2B1_TEST_ROOT}")
    add_test(NAME stage2b1_parametric_training_exactness
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${STAGE2B1_TEST_DRIVER}" -Mode Generate
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${STAGE2B1_TEST_ROOT}")
    set_tests_properties(stage2b1_parametric_training_exactness PROPERTIES
        FIXTURES_SETUP stage2b1_model)

    add_test(NAME stage2b1_serialization_reload_consistency
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${STAGE2B1_TEST_DRIVER}" -Mode Reload
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${STAGE2B1_TEST_ROOT}")
    add_test(NAME stage2b1_out_of_range_rejection
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${STAGE2B1_TEST_DRIVER}" -Mode OutOfRange
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${STAGE2B1_TEST_ROOT}")
    add_test(NAME stage2b1_fingerprint_rejection
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${STAGE2B1_TEST_DRIVER}" -Mode Fingerprint
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${STAGE2B1_TEST_ROOT}")
    add_test(NAME stage2b1_corrected_residual_gate
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${STAGE2B1_TEST_DRIVER}" -Mode Corrected
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${STAGE2B1_TEST_ROOT}")
    set_tests_properties(
        stage2b1_serialization_reload_consistency
        stage2b1_out_of_range_rejection
        stage2b1_fingerprint_rejection
        stage2b1_corrected_residual_gate
        PROPERTIES FIXTURES_REQUIRED stage2b1_model)

    add_test(NAME stage2a1_fixed_matrix_rom_regression
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${STAGE2B1_TEST_DRIVER}" -Mode Stage2A1
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${STAGE2B1_TEST_ROOT}")

    set(STAGE2C1_TEST_ROOT "${CMAKE_BINARY_DIR}/stage2c1_ctest")
    set(STAGE2C1_TEST_DRIVER "${CMAKE_SOURCE_DIR}/tests/run_stage2c1_test.ps1")
    add_test(NAME stage2c1_descriptor_matrix_symmetry
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${STAGE2C1_TEST_DRIVER}" -Mode Matrix
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${STAGE2C1_TEST_ROOT}")
    set_tests_properties(stage2c1_descriptor_matrix_symmetry PROPERTIES
        FIXTURES_SETUP stage2c1_two_model)

    foreach(stage2c1_test IN ITEMS
        Dc Step Serialization InitialCondition Fingerprint ReducedSpd)
        string(TOLOWER "${stage2c1_test}" stage2c1_test_lower)
        add_test(NAME "stage2c1_${stage2c1_test_lower}"
            COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
                -File "${STAGE2C1_TEST_DRIVER}" -Mode "${stage2c1_test}"
                -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${STAGE2C1_TEST_ROOT}")
        set_tests_properties("stage2c1_${stage2c1_test_lower}" PROPERTIES
            FIXTURES_REQUIRED stage2c1_two_model)
    endforeach()

    add_test(NAME stage2c1_ten_cube_transient_regression
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${STAGE2C1_TEST_DRIVER}" -Mode TenCube
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${STAGE2C1_TEST_ROOT}")

    set(LOCAL_BLOCK_ARNOLDI_TEST_ROOT "${CMAKE_BINARY_DIR}/local_block_arnoldi_ctest")
    set(LOCAL_BLOCK_ARNOLDI_TEST_DRIVER
        "${CMAKE_SOURCE_DIR}/tests/run_local_block_arnoldi_test.ps1")
    foreach(local_dynamic_test IN ITEMS Generate InitialConditions Waveforms)
        string(TOLOWER "${local_dynamic_test}" local_dynamic_test_lower)
        add_test(NAME "local_dynamic_schur_${local_dynamic_test_lower}"
            COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
                -File "${LOCAL_BLOCK_ARNOLDI_TEST_DRIVER}" -Mode "${local_dynamic_test}"
                -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
                -Root "${LOCAL_BLOCK_ARNOLDI_TEST_ROOT}")
        set_tests_properties("local_dynamic_schur_${local_dynamic_test_lower}"
            PROPERTIES RUN_SERIAL TRUE TIMEOUT 180)
    endforeach()

    set(LOCAL_BLOCK_ARNOLDI_TEN_TEST_ROOT
        "${CMAKE_BINARY_DIR}/local_block_arnoldi_ten_ctest")
    set(LOCAL_BLOCK_ARNOLDI_TEN_TEST_DRIVER
        "${CMAKE_SOURCE_DIR}/tests/run_local_block_arnoldi_ten_test.ps1")
    add_test(NAME local_dynamic_schur_ten_cube
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${LOCAL_BLOCK_ARNOLDI_TEN_TEST_DRIVER}" -Mode Generate
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${LOCAL_BLOCK_ARNOLDI_TEN_TEST_ROOT}")
    set_tests_properties(local_dynamic_schur_ten_cube PROPERTIES
        FIXTURES_SETUP local_dynamic_ten_model RUN_SERIAL TRUE TIMEOUT 180)
    add_test(NAME local_dynamic_schur_ten_template_reuse
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${LOCAL_BLOCK_ARNOLDI_TEN_TEST_DRIVER}" -Mode Reuse
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${LOCAL_BLOCK_ARNOLDI_TEN_TEST_ROOT}")
    set_tests_properties(local_dynamic_schur_ten_template_reuse PROPERTIES
        FIXTURES_REQUIRED local_dynamic_ten_model RUN_SERIAL TRUE TIMEOUT 180)

    add_test(NAME local_dynamic_schur_matrixfree
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${CMAKE_SOURCE_DIR}/tests/run_local_dynamic_schur_matrixfree_test.ps1"
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${CMAKE_BINARY_DIR}/local_dynamic_matrixfree_ctest")
    set_tests_properties(local_dynamic_schur_matrixfree PROPERTIES
        RUN_SERIAL TRUE TIMEOUT 180)

    add_test(NAME local_dynamic_schur_production
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${CMAKE_SOURCE_DIR}/tests/run_dynamic_schur_production_test.ps1"
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${CMAKE_BINARY_DIR}/local_dynamic_production_ctest")
    set_tests_properties(local_dynamic_schur_production PROPERTIES
        RUN_SERIAL TRUE TIMEOUT 180)

    add_test(NAME local_dynamic_schur_parallel_cache
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${CMAKE_SOURCE_DIR}/tests/run_dynamic_schur_parallel_cache_test.ps1"
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${CMAKE_BINARY_DIR}/local_dynamic_parallel_cache_ctest")
    set_tests_properties(local_dynamic_schur_parallel_cache PROPERTIES
        RUN_SERIAL TRUE TIMEOUT 180)

    add_test(NAME local_dynamic_operator_coarse_construction
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${CMAKE_SOURCE_DIR}/tests/run_operator_coarse_construction_test.ps1"
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${CMAKE_BINARY_DIR}/operator_coarse_construction_ctest")
    set_tests_properties(local_dynamic_operator_coarse_construction PROPERTIES
        RUN_SERIAL TRUE TIMEOUT 180)

    add_test(NAME package15_fastest_dynamic_schur
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${CMAKE_SOURCE_DIR}/tests/run_package15_fastest_path_test.ps1"
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${CMAKE_BINARY_DIR}/package15_fastest_ctest")
    set_tests_properties(package15_fastest_dynamic_schur PROPERTIES
        RUN_SERIAL TRUE TIMEOUT 180)

    foreach(package15_profile IN ITEMS smoke medium)
        add_test(NAME "package15_geometry_${package15_profile}"
            COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
                -File "${CMAKE_SOURCE_DIR}/tests/run_package15_geometry_test.ps1"
                -Profile "${package15_profile}"
                -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
                -Root "${CMAKE_BINARY_DIR}/package15_geometry_${package15_profile}_ctest")
        set_tests_properties("package15_geometry_${package15_profile}" PROPERTIES
            RUN_SERIAL TRUE TIMEOUT 180)
    endforeach()

    set(LOCAL_PORT_TEST_ROOT "${CMAKE_BINARY_DIR}/local_port_ctest")
    set(LOCAL_PORT_TEST_DRIVER
        "${CMAKE_SOURCE_DIR}/tests/run_local_port_dynamic_schur_test.ps1")
    add_test(NAME local_port_two_cube_basis
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${LOCAL_PORT_TEST_DRIVER}" -Mode Generate
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${LOCAL_PORT_TEST_ROOT}")
    set_tests_properties(local_port_two_cube_basis PROPERTIES
        FIXTURES_SETUP local_port_model RUN_SERIAL TRUE TIMEOUT 180)
    foreach(local_port_check IN ITEMS
        Projection Assembly Enrichment Flux Corrected Fingerprint)
        string(TOLOWER "${local_port_check}" local_port_check_lower)
        add_test(NAME "local_port_${local_port_check_lower}"
            COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
                -File "${LOCAL_PORT_TEST_DRIVER}" -Mode "${local_port_check}"
                -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${LOCAL_PORT_TEST_ROOT}")
        set_tests_properties("local_port_${local_port_check_lower}" PROPERTIES
            FIXTURES_REQUIRED local_port_model RUN_SERIAL TRUE TIMEOUT 180)
    endforeach()
    add_test(NAME local_port_serialization_reload
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${LOCAL_PORT_TEST_DRIVER}" -Mode Serialization
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${LOCAL_PORT_TEST_ROOT}")
    set_tests_properties(local_port_serialization_reload PROPERTIES
        RUN_SERIAL TRUE TIMEOUT 180)
    add_test(NAME local_port_ten_cube_regression
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${LOCAL_PORT_TEST_DRIVER}" -Mode TenCube
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${LOCAL_PORT_TEST_ROOT}")
    set_tests_properties(local_port_ten_cube_regression PROPERTIES
        RUN_SERIAL TRUE TIMEOUT 240)
    add_test(NAME local_port_template_reuse
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${LOCAL_PORT_TEST_DRIVER}" -Mode TemplateReuse
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${LOCAL_PORT_TEST_ROOT}")
    set_tests_properties(local_port_template_reuse PROPERTIES
        RUN_SERIAL TRUE TIMEOUT 240)

    set(LOCAL_PORT_RRAM_CONFIG
        "D:/CPP/TEST_CHATGPT/configs/rram26_from_sim_parameter_1_bshift_minus1.txt")
    if (EXISTS "${LOCAL_PORT_RRAM_CONFIG}")
        add_test(NAME local_port_rram26_integration
            COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
                -File "${LOCAL_PORT_TEST_DRIVER}" -Mode Rram
                -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${LOCAL_PORT_TEST_ROOT}"
                -Config "${LOCAL_PORT_RRAM_CONFIG}")
        set_tests_properties(local_port_rram26_integration PROPERTIES
            RUN_SERIAL TRUE TIMEOUT 900 LABELS "large;rram26")

        set(RRAM26_DYNAMIC_PRODUCTION_TEST_ROOT
            "${CMAKE_BINARY_DIR}/rram26_dynamic_production_ctest")
        add_test(NAME local_dynamic_rram26_port_production
            COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
                -File "${CMAKE_SOURCE_DIR}/tests/run_rram26_port_dynamic_schur_production_test.ps1"
                -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
                -Root "${RRAM26_DYNAMIC_PRODUCTION_TEST_ROOT}"
                -Config "${LOCAL_PORT_RRAM_CONFIG}")
        set_tests_properties(local_dynamic_rram26_port_production PROPERTIES
            RUN_SERIAL TRUE TIMEOUT 900 LABELS "large;rram26;production")
    endif()
    set(LOCAL_PORT_CHIPLET_CONFIG
        "D:/CPP/TEST_CHATGPT/chiplet_model/case_chiplet_config_horizontal.txt")
    if (EXISTS "${LOCAL_PORT_CHIPLET_CONFIG}")
        add_test(NAME local_port_chiplet_flux_regression
            COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
                -File "${LOCAL_PORT_TEST_DRIVER}" -Mode Chiplet
                -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${LOCAL_PORT_TEST_ROOT}"
                -Config "${LOCAL_PORT_CHIPLET_CONFIG}")
        set_tests_properties(local_port_chiplet_flux_regression PROPERTIES
            RUN_SERIAL TRUE TIMEOUT 900 LABELS "large;chiplet")

        set(CHIPLET_DYNAMIC_PRODUCTION_TEST_ROOT
            "${CMAKE_BINARY_DIR}/chiplet_dynamic_production_ctest")
        add_test(NAME local_dynamic_chiplet_production
            COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
                -File "${CMAKE_SOURCE_DIR}/tests/run_chiplet_dynamic_schur_production_test.ps1"
                -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
                -Root "${CHIPLET_DYNAMIC_PRODUCTION_TEST_ROOT}"
                -Config "${LOCAL_PORT_CHIPLET_CONFIG}")
        set_tests_properties(local_dynamic_chiplet_production PROPERTIES
            RUN_SERIAL TRUE TIMEOUT 900 LABELS "large;chiplet;production")
    endif()

    set(OPTIMAL_PORT_TEST_ROOT
        "${CMAKE_BINARY_DIR}/optimal_port_ctest")
    set(OPTIMAL_PORT_TEST_DRIVER
        "${CMAKE_SOURCE_DIR}/tests/run_optimal_port_test.ps1")
    add_test(NAME optimal_port_two_cube
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile
            -ExecutionPolicy Bypass
            -File "${OPTIMAL_PORT_TEST_DRIVER}" -Mode TwoCube
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${OPTIMAL_PORT_TEST_ROOT}")
    set_tests_properties(optimal_port_two_cube PROPERTIES
        FIXTURES_SETUP optimal_port_model TIMEOUT 120)
    add_test(NAME optimal_port_serialization
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile
            -ExecutionPolicy Bypass
            -File "${OPTIMAL_PORT_TEST_DRIVER}" -Mode Serialization
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${OPTIMAL_PORT_TEST_ROOT}")
    set_tests_properties(optimal_port_serialization PROPERTIES
        FIXTURES_REQUIRED optimal_port_model TIMEOUT 120)
    foreach(optimal_port_mode IN ITEMS Provenance Steklov Pcg Fgmres Generalized TargetSolvers Refinement ResidualKrylov Hybrid HistoryCompression GlobalRandomized ProjectionDiagnosis ProductionBasisStop)
        string(TOLOWER "${optimal_port_mode}" optimal_port_mode_lower)
        add_test(NAME "optimal_port_${optimal_port_mode_lower}"
            COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile
                -ExecutionPolicy Bypass
                -File "${OPTIMAL_PORT_TEST_DRIVER}"
                -Mode "${optimal_port_mode}"
                -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
                -Root "${OPTIMAL_PORT_TEST_ROOT}")
        set_tests_properties(
            "optimal_port_${optimal_port_mode_lower}"
            PROPERTIES TIMEOUT 120)
        if(optimal_port_mode STREQUAL "Provenance")
            set_tests_properties(
                "optimal_port_${optimal_port_mode_lower}"
                PROPERTIES FIXTURES_REQUIRED optimal_port_model)
        endif()
        if(optimal_port_mode STREQUAL "Refinement")
            set_tests_properties(
                "optimal_port_${optimal_port_mode_lower}"
                PROPERTIES RUN_SERIAL TRUE TIMEOUT 300)
        endif()
    endforeach()
    add_test(NAME optimal_port_ten_transfer
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile
            -ExecutionPolicy Bypass
            -File "${OPTIMAL_PORT_TEST_DRIVER}" -Mode TenTransfer
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${OPTIMAL_PORT_TEST_ROOT}")
    set_tests_properties(optimal_port_ten_transfer PROPERTIES
        RUN_SERIAL TRUE TIMEOUT 180 LABELS "ten-cube")
    add_test(NAME optimal_port_randomized_two_cube
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile
            -ExecutionPolicy Bypass
            -File "${OPTIMAL_PORT_TEST_DRIVER}" -Mode Randomized
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${OPTIMAL_PORT_TEST_ROOT}")
    set_tests_properties(optimal_port_randomized_two_cube PROPERTIES
        FIXTURES_SETUP randomized_port_model TIMEOUT 180)
    add_test(NAME optimal_port_randomized_serialization
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile
            -ExecutionPolicy Bypass
            -File "${OPTIMAL_PORT_TEST_DRIVER}"
            -Mode RandomizedSerialization
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${OPTIMAL_PORT_TEST_ROOT}")
    add_test(NAME optimal_port_randomized_reproducibility
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile
            -ExecutionPolicy Bypass
            -File "${OPTIMAL_PORT_TEST_DRIVER}"
            -Mode RandomizedReproducibility
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${OPTIMAL_PORT_TEST_ROOT}")
    set_tests_properties(
        optimal_port_randomized_serialization
        optimal_port_randomized_reproducibility
        PROPERTIES FIXTURES_REQUIRED randomized_port_model TIMEOUT 120)

    set(LOCAL_ROM_TEST_ROOT "${CMAKE_BINARY_DIR}/local_rom_ctest")
    set(LOCAL_ROM_TEST_DRIVER "${CMAKE_SOURCE_DIR}/tests/run_local_rom_test.ps1")
    add_test(NAME local_rom_two_cube_exactness
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${LOCAL_ROM_TEST_DRIVER}" -Mode Generate
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${LOCAL_ROM_TEST_ROOT}")
    set_tests_properties(local_rom_two_cube_exactness PROPERTIES
        FIXTURES_SETUP local_rom_model)
    foreach(local_rom_test IN ITEMS Assembly Mapping Flux Serialization Fingerprint Corrected MatrixFree)
        string(TOLOWER "${local_rom_test}" local_rom_test_lower)
        add_test(NAME "local_rom_${local_rom_test_lower}"
            COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
                -File "${LOCAL_ROM_TEST_DRIVER}" -Mode "${local_rom_test}"
                -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${LOCAL_ROM_TEST_ROOT}")
        set_tests_properties("local_rom_${local_rom_test_lower}" PROPERTIES
            FIXTURES_REQUIRED local_rom_model)
    endforeach()

    set(LOCAL_ROM_TEN_TEST_ROOT "${CMAKE_BINARY_DIR}/local_rom_ten_ctest")
    set(LOCAL_ROM_TEN_TEST_DRIVER "${CMAKE_SOURCE_DIR}/tests/run_local_rom_ten_test.ps1")
    add_test(NAME local_rom_ten_generate
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${LOCAL_ROM_TEN_TEST_DRIVER}" -Mode Generate
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${LOCAL_ROM_TEN_TEST_ROOT}")
    set_tests_properties(local_rom_ten_generate PROPERTIES
        FIXTURES_SETUP local_rom_ten_model RUN_SERIAL TRUE TIMEOUT 300)
    foreach(local_rom_ten_test IN ITEMS
        Structure LocalBases ReducedBlocks InterfaceAssembly Mapping Pure Serialization Fingerprint Corrected)
        string(TOLOWER "${local_rom_ten_test}" local_rom_ten_test_lower)
        add_test(NAME "local_rom_ten_${local_rom_ten_test_lower}"
            COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
                -File "${LOCAL_ROM_TEN_TEST_DRIVER}" -Mode "${local_rom_ten_test}"
                -Exe "$<TARGET_FILE:SIPGHeatDDM3D>" -Root "${LOCAL_ROM_TEN_TEST_ROOT}")
        set_tests_properties("local_rom_ten_${local_rom_ten_test_lower}" PROPERTIES
            FIXTURES_REQUIRED local_rom_ten_model RUN_SERIAL TRUE TIMEOUT 600)
    endforeach()

    set(NUMERICAL_CORRECTNESS_TEST_ROOT
        "${CMAKE_BINARY_DIR}/numerical_correctness_ctest")
    set(NUMERICAL_CORRECTNESS_TEST_DRIVER
        "${CMAKE_SOURCE_DIR}/tests/run_numerical_correctness_test.ps1")
    foreach(correctness_mode IN ITEMS
        HeterogeneousSipg PenaltyScaling NipgReject FullResidualGate)
        string(TOLOWER "${correctness_mode}" correctness_mode_lower)
        add_test(NAME "numerical_correctness_${correctness_mode_lower}"
            COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile
                -ExecutionPolicy Bypass
                -File "${NUMERICAL_CORRECTNESS_TEST_DRIVER}"
                -Mode "${correctness_mode}"
                -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
                -Root "${NUMERICAL_CORRECTNESS_TEST_ROOT}")
        set_tests_properties(
            "numerical_correctness_${correctness_mode_lower}"
            PROPERTIES RUN_SERIAL TRUE TIMEOUT 180)
    endforeach()

    add_test(NAME schur_parallel_local_and_proxy_cache
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${CMAKE_SOURCE_DIR}/tests/run_schur_parallel_cache_test.ps1"
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${CMAKE_BINARY_DIR}/schur_parallel_cache_ctest")
    set_tests_properties(schur_parallel_local_and_proxy_cache PROPERTIES
        RUN_SERIAL TRUE TIMEOUT 180)

    add_test(NAME p2_p1_auxiliary_two_cube
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${CMAKE_SOURCE_DIR}/tests/run_p2_p1_auxiliary_test.ps1"
            -Exe "$<TARGET_FILE:SIPGHeatDDM3D>"
            -Root "${CMAKE_BINARY_DIR}/p2_p1_auxiliary_ctest")
    set_tests_properties(p2_p1_auxiliary_two_cube PROPERTIES
        RUN_SERIAL TRUE TIMEOUT 180)
endif()
