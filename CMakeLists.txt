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
find_program(CMAKE_CXX_COMPILER "hipcc" PATHS ${ROCM_PATH})


project(self VERSION 1.0.0
	DESCRIPTION "Spectral Element Library in Fortran"
	LANGUAGES Fortran HIP C CXX)


# set( CMAKE_CXX_COMPILER "${ROCM_PATH}/hip/bin/hipcc" )
# if(NOT DEFINED CUDA_PATH)
#     if(NOT DEFINED ENV{CUDA_PATH})
#         set(CUDA_PATH "/usr/local/cuda" CACHE PATH "Path to which CUDA has been installed")
#     else()
#         set(CUDA_PATH $ENV{CUDA_PATH} CACHE PATH "Path to which CUDA has been installed")
#     endif()
# endif()

# conditional on GPU target
# if nvidia
# "-L$CUDA_PATH/targets/x86_64-linux/lib -lcudart"
# if amd
# -L$ROCM_PATH/lib -lamdhip64 -Wl,-rpath=$ROCM_PATH/lib 
#set( HIP_Fortran_LFLAGS "-L${ROCM_PATH}/lib -lamdhip64 -Wl,-rpath=${ROCM_PATH}/lib")

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
    set( CMAKE_Fortran_FLAGS "-cpp -ffree-line-length-512" )
    set( CMAKE_Fortran_FLAGS_DEBUG "-pg -g -O0 -C -Wall -fbounds-check -fbacktrace --coverage -ffpe-trap=invalid,zero,overflow" )
    set( CMAKE_Fortran_FLAGS_RELEASE "-O3" )
endif()

# ------ dependencies ------ #
# HIP
find_package(hip REQUIRED)

# MPI
find_package(MPI COMPONENTS Fortran REQUIRED)
# add_definitions(${MPI_Fortran_COMPILE_DEFINITIONS})
# include_directories(${MPI_Fortran_INCLUDE_PATH})
# link_directories(${MPI_Fortran_LIBRARIES})

# HDF5 : See https://cmake.org/cmake/help/latest/module/FindHDF5.html
find_package(HDF5 REQUIRED Fortran)
# add_definitions(${HDF5_Fortran_DEFINITIONS})
# include_directories(${HDF5_INCLUDE_DIRS})

# JSON-Fortran
find_library(JSONFORTRAN_LIBRARIES NAMES jsonfortran REQUIRED)
find_path(JSONFORTRAN_INCLUDE_DIRS json_module.mod)


# FEQ-Parse
find_library(FEQPARSE_LIBRARIES NAMES feqparse REQUIRED)
find_path(FEQPARSE_INCLUDE_DIRS feqparse.mod)


# ----- end dependencies ----- #

# Libraries
add_subdirectory(${CMAKE_SOURCE_DIR}/src)

# Driver programs
add_subdirectory(${CMAKE_SOURCE_DIR}/bin)
