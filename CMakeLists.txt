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
	LANGUAGES Fortran HIP)

# ------ dependencies ------ #
# MPI
find_package(MPI COMPONENTS Fortran REQUIRED)

# HDF5 : See https://cmake.org/cmake/help/latest/module/FindHDF5.html
find_package(HDF5 COMPONENTS Fortran
	HINTS HDF5_PREFER_PARALLEL REQUIRED)

# JSON-Fortran
find_package(json-fortran REQUIRED)

# FEQ-Parse
find_package(feq-parse REQUIRED)

# ----- end dependencies ----- #

#add_executable(square square.hip)
#add_library(self STATIC src/ src/hip src/models src/models/hip)
#target_link_libraries(self PRIVATE hip::device)
