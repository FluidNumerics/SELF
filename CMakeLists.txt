#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
# 
#  Maintainers : support@fluidnumerics.com
#  Official Repository : https://github.com/FluidNumerics/self/
# 
#  Copyright © 2024 Fluid Numerics LLC
# 
#  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
#  1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 
#  2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the distribution.
# 
#  3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
#     this software without specific prior written permission.
# 
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
#  THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
#  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
cmake_minimum_required(VERSION 3.21)
cmake_policy(VERSION 3.21...3.27)

# Defaults to "none", which disables OpenMP builds
# Override with -DOMP_TARGET
# - Set to "none" to disable OpenMP
# - Set to "multicore" for multicore cpu
# - GPU offloading via openmp not supported (not by our own choice)
option(SELF_ENABLE_OPENMP "Option to enable OpenMP multicore target. Currently only tested with nvfortran (Default Off)"  OFF)
option(SELF_ENABLE_TESTING "Option to enable build of tests. (Default On)"  ON)
option(SELF_ENABLE_GPU "Option to enable GPU backend. Requires either CUDA or HIP. (Default Off)"  OFF)
option(SELF_ENABLE_APU "Option to enable APU backend. Requires either CUDA or HIP. (Default Off)"  OFF)
option(SELF_ENABLE_MPI "Option to enable MPI for distributed memory parallelism (domain decomposition). (Default Off)"  OFF)
option(SELF_ENABLE_DOUBLE_PRECISION "Option to enable double precision for floating point arithmetic. (Default On)"  ON)

if(SELF_ENABLE_OPENMP)
    set(OMP_TARGET "multicore" CACHE STRING "Target architecture for openmp.")
endif()


if(NOT DEFINED ROCM_PATH)
    if(NOT DEFINED ENV{ROCM_PATH})
        set(ROCM_PATH "/opt/rocm/" CACHE PATH "Path to which ROCm has been installed")
    else()
        set(ROCM_PATH $ENV{ROCM_PATH} CACHE PATH "Path to which ROCm has been installed")
    endif()
endif()

#SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIS_D}/cmake/modules/")
include("${CMAKE_CURRENT_LIST_DIR}/cmake/modules/OpenMPTargetOffloadFlags.cmake")
# C Language is needed in order to verify Fortran compiler is C-interoperable
# CXX language is needed to properly find "hip" package
project(self VERSION 1.0.0
	DESCRIPTION "Spectral Element Library in Fortran"
	LANGUAGES Fortran C CXX)

# Fortran compiler requirements
INCLUDE(FortranCInterface)
FortranCInterface_VERIFY()
IF(NOT FortranCInterface_VERIFIED_C)
	MESSAGE(FATAL_ERROR "Fortran compiler must support C Interface")
ENDIF(NOT FortranCInterface_VERIFIED_C)
	
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Default Fortran compiler flags
# Fortran compiler flags
if( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
    set( CMAKE_Fortran_FLAGS "${CMAKE_FORTRAN_FLAGS} -cpp -ffree-line-length-512" )
    set( CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -C -Wall -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow" )
    set( CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_DEBUG} --coverage")
    set( CMAKE_Fortran_FLAGS_PROFILE "-pg -O3")
    set( CMAKE_Fortran_FLAGS_RELEASE "-O3" )
    set( CMAKE_C_FLAGS_COVERAGE "-g -O0 --coverage")

elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" )
    set( CMAKE_Fortran_FLAGS "-fpp" )
    set( CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -debug all -check all -traceback" )
    set( CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_DEBUG} -prof-gen=srcpos")
    set( CMAKE_Fortran_FLAGS_PROFILE "-O3")
    set( CMAKE_Fortran_FLAGS_RELEASE "-O3" )
    set( CMAKE_C_FLAGS_COVERAGE "-g -O0")

elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM" )
    set( CMAKE_Fortran_FLAGS "-fpp" )
    set( CMAKE_Fortran_FLAGS_DEBUG "-debug all -check all,nouninit -traceback" )
    set( CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_DEBUG}") # Code coverage not available with ifx
    set( CMAKE_Fortran_FLAGS_PROFILE "-O3")
    set( CMAKE_Fortran_FLAGS_RELEASE "-O3" )
    set( CMAKE_C_FLAGS_COVERAGE "-g -O0")

elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Flang" )
    set( CMAKE_Fortran_FLAGS "-Mpreprocess" )
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -O0 -g -Wall -pedantic -save-temps" )
    set( CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_DEBUG}") # Code coverage not available with ifx
    set( CMAKE_Fortran_FLAGS_PROFILE "-O3")
    set( CMAKE_Fortran_FLAGS_RELEASE "-O3" )
    set( CMAKE_C_FLAGS_COVERAGE "-g -O0")

elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NVHPC" )
    set( CMAKE_Fortran_FLAGS "-Mpreprocess" )
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -O0 -Wall -C -g -traceback" )
    set( CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_DEBUG}") # Code coverage not available with ifx
    set( CMAKE_Fortran_FLAGS_PROFILE "-O3 -Mstack_arrays -Minfo")
    set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -Mstack_arrays -Minfo" )
    set( CMAKE_C_FLAGS_COVERAGE "-g -O0")
endif()

# Check for openmp support if offload target is provided
if( SELF_ENABLE_OPENMP )
    find_package(OpenMP COMPONENTS Fortran REQUIRED)
    if( ${OpenMP_Fortran_FOUND} )
        if( "${OMP_TARGET}" STREQUAL "multicore" )
            get_offload_flags(RESULT OFFLOAD_FLAGS)
        else()
            get_offload_flags(ARCH ${OMP_TARGET} RESULT OFFLOAD_FLAGS)
        endif()
        set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OFFLOAD_FLAGS}" )
        set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OFFLOAD_FLAGS}" )
        set( CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_COVERAGE} ${OFFLOAD_FLAGS}")
        set( CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_PROFILE} ${OFFLOAD_FLAGS}")
        set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OFFLOAD_FLAGS}" )
        if( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NVHPC" )
            set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Minfo=mp" )
            set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Minfo=mp" )
            set( CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_COVERAGE} -Minfo=mp")
            set( CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_PROFILE} -Minfo=mp")
            set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -Minfo=mp" )
        endif()
    else()
        message( WARNING "OpenMP target set to ${OMP_TARGET} but OpenMP_Fortran not found. GPU Offloading disabled.")
    endif()
endif()


# MPI
find_package(MPI COMPONENTS Fortran REQUIRED)

# HDF5 : See https://cmake.org/cmake/help/latest/module/FindHDF5.html
find_package(HDF5 REQUIRED Fortran)

# # JSON-Fortran
# find_library(JSONFORTRAN_LIBRARIES NAMES jsonfortran REQUIRED)
# find_path(JSONFORTRAN_INCLUDE_DIRS json_module.mod)


# FEQ-Parse
find_library(FEQPARSE_LIBRARIES NAMES feqparse REQUIRED)
find_path(FEQPARSE_INCLUDE_DIRS feqparse.mod)

if(SELF_ENABLE_DOUBLE_PRECISION)
    message("Enabling Double Precision")
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DDOUBLE_PRECISION" )
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -DDOUBLE_PRECISION" )
    set( CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_COVERAGE} -DDOUBLE_PRECISION")
    set( CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_PROFILE} -DDOUBLE_PRECISION")
    set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -DDOUBLE_PRECISION" )

    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDOUBLE_PRECISION" )
    set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDOUBLE_PRECISION" )
    set( CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_COVERAGE} -DDOUBLE_PRECISION")
    set( CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_PROFILE} -DDOUBLE_PRECISION")
    set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DDOUBLE_PRECISION" )
endif()

if(SELF_ENABLE_GPU)
    message("Enabling GPU Support")
    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DENABLE_GPU" )
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -DENABLE_GPU" )
    set( CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_COVERAGE} -DENABLE_GPU")
    set( CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_PROFILE} -DENABLE_GPU")
    set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -DENABLE_GPU" )

    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_GPU" )
    set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DENABLE_GPU" )
    set( CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_COVERAGE} -DENABLE_GPU")
    set( CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_PROFILE} -DENABLE_GPU")
    set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DENABLE_GPU" )

    # HIP (Optional)
    # IF HIP is found, we prefer HIP. IF HIP is not found, then we fall back to look for CUDA
    # only. If neither HIP nor CUDA are found, Cmake will error out.
    find_package(hip)
    if(hip_FOUND)
        find_package(hipblas REQUIRED)
        message("HIP found. Enabling HIP language.")
        enable_language(HIP)
        set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DHAVE_HIP" )
        set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -DHAVE_HIP" )
        set( CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_COVERAGE} -DHAVE_HIP")
        set( CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_PROFILE} -DHAVE_HIP")
        set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -DHAVE_HIP" )

        set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_HIP" )
        set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DHAVE_HIP" )
        set( CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_COVERAGE} -DHAVE_HIP")
        set( CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_PROFILE} -DHAVE_HIP")
        set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DHAVE_HIP" )

        set( BACKEND_LIBRARIES hip::device roc::hipblas)
    else()
        # CUDA (Optional)
        find_package(cuda)
        if(cuda_FOUND)
            message("CUDA found. Enabling CUDA language.")
            enable_language(CUDA)
            set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DHAVE_CUDA" )
            set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -DHAVE_CUDA" )
            set( CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_COVERAGE} -DHAVE_CUDA")
            set( CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_PROFILE} -DHAVE_CUDA")
            set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -DHAVE_CUDA" )
            
            set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_CUDA" )
            set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DHAVE_CUDA" )
            set( CMAKE_CXX_FLAGS_COVERAGE "${CMAKE_CXX_FLAGS_COVERAGE} -DHAVE_CUDA")
            set( CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS_PROFILE} -DHAVE_CUDA")
            set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DHAVE_CUDA" )
            
            # TO DO - need cuda libraries and hipblas libraries
        else()
            message( FATAL_ERROR "Enabling GPU support requires either HIP or CUDA." )
        endif()
    endif()
endif()
# get_cmake_property(_variableNames VARIABLES)
# list (SORT _variableNames)
# foreach (_variableName ${_variableNames})
#     message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()
# ----- end dependencies ----- #

# Libraries
add_subdirectory(${CMAKE_SOURCE_DIR}/src)


if(SELF_ENABLE_TESTING)
    enable_testing()
    add_subdirectory(${CMAKE_SOURCE_DIR}/test)
endif()
