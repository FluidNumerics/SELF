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

option(SELF_ENABLE_MULTITHREADING "Option to enable CPU multithreading for `do concurrent` loop blocks."  OFF)
option(SELF_ENABLE_TESTING "Option to enable build of tests. (Default On)"  ON)
option(SELF_ENABLE_GPU "Option to enable GPU backend. Requires either CUDA or HIP. (Default Off)"  OFF)
option(SELF_ENABLE_APU "Option to enable APU backend. Requires either CUDA or HIP. (Default Off)"  OFF)
option(SELF_ENABLE_DOUBLE_PRECISION "Option to enable double precision for floating point arithmetic. (Default On)"  ON)

set(SELF_MPIEXEC_NUMPROCS "2" CACHE STRING "The number of MPI ranks to use to launch MPI tests. Only used when launching test programs via ctest.")
set(SELF_MPIEXEC_OPTIONS "" CACHE STRING "Any additional options, such as binding options, to use for MPI tests.Only used when launching test programs via ctest. Defaults to nothing")

if(SELF_ENABLE_MULTITHREADING)
    set(SELF_MULITHREADING_NTHREADS "4" CACHE STRING "Number of threads to use for `do concurrent` loop blocks. This option is only used with GNU compilers. Other compilers use OMP_NUM_THREADS environment variable at runtime.")
endif()

if(NOT DEFINED ROCM_PATH)
    if(NOT DEFINED ENV{ROCM_PATH})
        set(ROCM_PATH "/opt/rocm/" CACHE PATH "Path to which ROCm has been installed")
    else()
        set(ROCM_PATH $ENV{ROCM_PATH} CACHE PATH "Path to which ROCm has been installed")
    endif()
endif()

# C Language is needed in order to verify Fortran compiler is C-interoperable
# CXX language is needed to properly find "hip" package
project(self VERSION 1.0.0
	DESCRIPTION "Spectral Element Library in Fortran"
	LANGUAGES Fortran C CXX)

# Fortran compiler requirements
include(FortranCInterface)
FortranCInterface_VERIFY()
if(NOT FortranCInterface_VERIFIED_C)
	message(FATAL_ERROR "Fortran compiler must support C Interface")
endif(NOT FortranCInterface_VERIFIED_C)
	
if(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
endif(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Default Fortran compiler flags
# Fortran compiler flags
if( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
    set( CMAKE_Fortran_FLAGS "${CMAKE_FORTRAN_FLAGS} -cpp -ffree-line-length-512" )
    set( CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -C -Wall -fbounds-check -fbacktrace" )
    #set( CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -C -Wall -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow" )
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
if( SELF_ENABLE_MULTITHREADING )
    if( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
        set( OFFLOAD_FLAGS "-ftree-parallelize-loops=${SELF_MULTITHREADING_NTHREADS} -fopt-info-loop" )
    elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" )
        set( OFFLOAD_FLAGS "-parallel -qopt-report -qopt-report-phase=par" )
    elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "IntelLLVM" )
        set( OFFLOAD_FLAGS "-qopenmp -qopt-report" )
    elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Flang" )
        set( OFFLOAD_FLAGS "-fopenmp" )
    elseif( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NVHPC" )
        set( OFFLOAD_FLAGS "-stdpar=multicore -Minfo=stdpar,accel" )
    endif()

    set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OFFLOAD_FLAGS}" )
    set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OFFLOAD_FLAGS}" )
    set( CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_COVERAGE} ${OFFLOAD_FLAGS}")
    set( CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_PROFILE} ${OFFLOAD_FLAGS}")
    set( CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OFFLOAD_FLAGS}" )

endif()

# MPI
find_package(MPI COMPONENTS Fortran C REQUIRED)

# HDF5 : See https://cmake.org/cmake/help/latest/module/FindHDF5.html
find_package(HDF5 REQUIRED Fortran)

# # JSON-Fortran
# find_library(JSONFORTRAN_LIBRARIES NAMES jsonfortran REQUIRED)
# find_path(JSONFORTRAN_INCLUDE_DIRS json_module.mod)


# FEQ-Parse
find_library(FEQPARSE_LIBRARIES NAMES feqparse REQUIRED)
find_path(FEQPARSE_INCLUDE_DIRS feqparse.mod)

if(SELF_ENABLE_DOUBLE_PRECISION)
    message("-- SELF Build System : Enabling Double Precision")
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
    message("-- SELF Build System : Enabling GPU Support")
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

    # Check MPI for GPU awareness
    # Add SELF's cmake module directory to the search path
    set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")

    # Include the custom module
    include(CheckMPIFeatures)
    CheckMPIFeatures()

    # HIP (Optional)
    # IF HIP is found, we prefer HIP. IF HIP is not found, then we fall back to look for CUDA
    # only. If neither HIP nor CUDA are found, Cmake will error out.
    find_package(hip)
    if(hip_FOUND)
        if(MPI_HAS_QUERY_HIP_SUPPORT)
            find_package(hipblas REQUIRED)
            #message("-- HIP found. Enabling HIP language.")
            #enable_language(HIP)
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
            message( FATAL_ERROR "MPI installation is not GPU-aware" )
        endif()
    else()
        # CUDA (Optional)
        find_package(cuda)
        if(cuda_FOUND)
            if(MPI_HAS_QUERY_CUDA_SUPPORT)
                #message("-- CUDA found. Enabling CUDA language.")
                #enable_language(CUDA)
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
                message( FATAL_ERROR "MPI installation is not GPU-aware" )
            endif()
        else()
            message( FATAL_ERROR "Enabling GPU support requires either HIP or CUDA." )
        endif()
    endif()
endif()

# ----- end dependencies ----- #

# Libraries
add_subdirectory(${CMAKE_SOURCE_DIR}/src)


if(SELF_ENABLE_TESTING)
    enable_testing()
    add_subdirectory(${CMAKE_SOURCE_DIR}/test)
endif()
