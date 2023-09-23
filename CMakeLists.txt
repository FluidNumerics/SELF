#cmake_minimum_required(VERSION 3.5)
#cmake_policy(VERSION 3.5...3.27)
#
# Invoke with
#  FC=gfortran \
#  AMDDeviceLibs_DIR=/opt/rocm/lib/cmake/AMDDeviceLibs \
#  amd_comgr_DIR=/opt/rocm/lib/cmake/amd_comgr/ \
#  cmake ../
#
cmake_minimum_required(VERSION 3.21) # HIP language support requires 3.21
cmake_policy(VERSION 3.21.3...3.27)

project(self VERSION 1.0.0
	DESCRIPTION "Spectral Element Library in Fortran"
	LANGUAGES Fortran HIP C)

# Fortran compiler requirements
INCLUDE(FortranCInterface)
FortranCInterface_VERIFY()
IF(NOT FortranCInterface_VERIFIED_C)
	MESSAGE(FATAL_ERROR "Fortran compiler must support C Interface")
ENDIF(NOT FortranCInterface_VERIFIED_C)
	
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# ------ dependencies ------ #
# MPI
find_package(MPI COMPONENTS Fortran REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_DEFINITIONS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})

# HDF5 : See https://cmake.org/cmake/help/latest/module/FindHDF5.html
find_package(HDF5 REQUIRED Fortran)
add_definitions(${HDF5_Fortran_DEFINITIONS})
include_directories(${HDF5_INCLUDE_DIRS})

# JSON-Fortran
#find_package(jsonfortran REQUIRED)
find_library(libjsonfortran NAMES jsonfortran REQUIRED)
find_path(incjsonfortran json_module.mod)

# FEQ-Parse
#find_package(feqparse REQUIRED)
find_library(libfeqparse NAMES feqparse REQUIRED)
find_path(incfeqparse feqparse.mod)

# ----- end dependencies ----- #

add_subdirectory(${CMAKE_SOURCE_DIR}/src)
add_subdirectory(${CMAKE_SOURCE_DIR}/bin)



# add_executable(self bin/SELF_Main.f90)
# target_link_libraries(self PRIVATE ${libjsonfortran} ${libfeqparse})
# target_include_directories(self PRIVATE ${incjsonfortran} ${incfeqparse})