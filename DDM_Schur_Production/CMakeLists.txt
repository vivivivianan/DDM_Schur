cmake_minimum_required(VERSION 3.20)

# This project intentionally builds one executable and one algorithm. Research
# experiments live in the source repository and are not exposed here.
project(DynamicSchurProduction LANGUAGES CXX)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(DynamicSchurProduction
    src/main.cpp
    src/internal/removed_research_methods.cpp
    src/ddm_schur/interface_operator.cpp
    src/mor/model_io.cpp
    src/mor/mor_diagnostics.cpp
    src/mor/source_parameterization.cpp
    src/mor/local/local_subdomain_model.cpp
    src/mor/local/local_reduced_schur.cpp
    src/mor/transient/thermal_descriptor_system.cpp
    src/mor/transient/block_arnoldi.cpp
    src/mor/transient/transient_input_waveform.cpp
    src/mor/transient/local_dynamic_schur.cpp)

target_include_directories(DynamicSchurProduction PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Cache files record a source identifier so models from incompatible binaries
# cannot be reused silently. This clean export has its own stable identifier.
target_compile_definitions(DynamicSchurProduction PRIVATE
    DDM_SCHUR_SOURCE_COMMIT="production-clean-v1")

find_package(OpenMP REQUIRED)
target_link_libraries(DynamicSchurProduction PRIVATE OpenMP::OpenMP_CXX)

# The validated fast path relies on MKL PARDISO. A dense fallback would change
# both scalability and memory behavior, so configuration fails if MKL is absent.
set(ONEAPI_MKL_HINTS
    "$ENV{MKLROOT}"
    "C:/Program Files (x86)/Intel/oneAPI/mkl/latest")
find_path(MKL_INCLUDE_DIR NAMES mkl.h
    PATH_SUFFIXES include PATHS ${ONEAPI_MKL_HINTS} REQUIRED)
find_library(MKL_RT_LIBRARY NAMES mkl_rt
    PATH_SUFFIXES lib PATHS ${ONEAPI_MKL_HINTS} REQUIRED)
target_compile_definitions(DynamicSchurProduction PRIVATE USE_MKL_PARDISO)
target_include_directories(DynamicSchurProduction PRIVATE "${MKL_INCLUDE_DIR}")
target_link_libraries(DynamicSchurProduction PRIVATE "${MKL_RT_LIBRARY}")

if(MSVC)
    target_compile_options(DynamicSchurProduction PRIVATE /W4 /permissive- /utf-8)
    # Whole-program optimization discards dormant compatibility branches and
    # keeps the shipped binary focused on the fixed production entry point.
    target_compile_options(DynamicSchurProduction PRIVATE
        $<$<CONFIG:Release>:/GL>)
    target_link_options(DynamicSchurProduction PRIVATE
        $<$<CONFIG:Release>:/LTCG>)
    target_link_libraries(DynamicSchurProduction PRIVATE Psapi)
else()
    target_compile_options(DynamicSchurProduction PRIVATE -Wall -Wextra -pedantic)
endif()

# Copy the runtime libraries beside the executable for a self-contained
# Release directory on Windows. Missing optional architecture DLLs are skipped.
if(WIN32)
    set(ONEAPI_RUNTIME_HINTS
        "$ENV{MKLROOT}/bin"
        "$ENV{CMPLR_ROOT}/bin"
        "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/bin"
        "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin")
    set(runtime_index 0)
    foreach(runtime_name IN ITEMS
            mkl_rt.2.dll mkl_core.2.dll mkl_intel_thread.2.dll
            mkl_avx2.2.dll mkl_def.2.dll libiomp5md.dll)
        math(EXPR runtime_index "${runtime_index} + 1")
        find_file(runtime_path_${runtime_index} NAMES ${runtime_name}
            PATHS ${ONEAPI_RUNTIME_HINTS})
        if(runtime_path_${runtime_index})
            add_custom_command(TARGET DynamicSchurProduction POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${runtime_path_${runtime_index}}"
                    "$<TARGET_FILE_DIR:DynamicSchurProduction>")
        endif()
    endforeach()
endif()

if(BUILD_TESTING AND WIN32)
    find_program(POWERSHELL_EXECUTABLE NAMES pwsh powershell REQUIRED)
    add_test(NAME production_smoke
        COMMAND "${POWERSHELL_EXECUTABLE}" -NoProfile -ExecutionPolicy Bypass
            -File "${CMAKE_CURRENT_SOURCE_DIR}/tests/run_production_smoke.ps1"
            -Exe "$<TARGET_FILE:DynamicSchurProduction>"
            -Root "${CMAKE_BINARY_DIR}/production_smoke"
            -BinaryRoot "${CMAKE_BINARY_DIR}")
    set_tests_properties(production_smoke PROPERTIES RUN_SERIAL TRUE TIMEOUT 180)
endif()
