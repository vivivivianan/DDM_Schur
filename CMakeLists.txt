cmake_minimum_required(VERSION 3.15)

project(SIPGHeatDDM3D LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(SIPGHeatDDM3D
    src/application.hpp
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
)

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
