#cmake_minimum_required(VERSION 3.5)
#cmake_policy(VERSION 3.5...3.27)
#
# Invoke with
#  FC=gfortran \
#  AMDDeviceLibs_DIR=/opt/rocm/lib/cmake/AMDDeviceLibs \
#  amd_comgr_DIR=/opt/rocm/lib/cmake/amd_comgr/ \
#  cmake ../
#
cmake_minimum_required(VERSION 3.21) # HIP language support requires 3.21 or gretater
cmake_policy(VERSION 3.21...3.27)

if(NOT DEFINED ROCM_PATH)
    if(NOT DEFINED ENV{ROCM_PATH})
        set(ROCM_PATH "/opt/rocm/" CACHE PATH "Path to which ROCm has been installed")
    else()
        set(ROCM_PATH $ENV{ROCM_PATH} CACHE PATH "Path to which ROCm has been installed")
    endif()
endif()
set( CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -fPIE")

# C Language is needed in order to verify Fortran compiler is C-interoperable
# CXX language is needed to properly find "hip" package
project(self VERSION 1.0.0
	DESCRIPTION "Spectral Element Library in Fortran"
	LANGUAGES Fortran HIP C CXX)

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

endif()

# ------ dependencies ------ #
# HIP
#find_package(hipfort REQUIRED COMPONENTS hip hipblas)
find_package(hip REQUIRED)
find_package(hipblas REQUIRED)
find_library(HIPFORT_LIBRARIES NAMES libhipfort-amdgcn.a HINTS ENV HIPFORT_ROOT REQUIRED)
find_path(HIPFORT_INCLUDE_DIRS hipfort.mod HINTS ENV HIPFORT_ROOT PATH_SUFFIXES include/hipfort/amdgcn/ REQUIRED)


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


# ----- end dependencies ----- #

# Libraries
add_subdirectory(${CMAKE_SOURCE_DIR}/src)

enable_testing()

# Tests
add_subdirectory(${CMAKE_SOURCE_DIR}/test)

# Examples (additional tests; standalone models)
#add_subdirectory(${CMAKE_SOURCE_DIR}/example)
