# =========================================================================
# SELF CMake configuration
# =========================================================================
# NOTE: use minimum CMake version required by tools/libs (Paraview, HDF5, FFTW)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.2)
PROJECT(SELF)
ENABLE_LANGUAGE(Fortran)

# =========================================================================
# Build type
# =========================================================================
# make sure that the default is a RELEASE
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release Profile Coverage."
      FORCE)
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release Profile Coverage)
ENDIF (NOT CMAKE_BUILD_TYPE)

IF (CMAKE_BUILD_TYPE MATCHES "Debug")
  ADD_DEFINITIONS("-DDEBUG")
ENDIF()

# =========================================================================
# COMPILER FLAGS
# =========================================================================

# FFLAGS depend on the compiler
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # set Flags
  SET (CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -DGNU")
  SET (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -march=native  -finline-functions")
  SET (CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS} -pg -O3 -march=native  -finline-functions")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -g -O0 -ggdb3 -fbounds-check -finit-real=nan -fbacktrace  -Wall")
  SET (CMAKE_Fortran_FLAGS_COVERAGE   "${CMAKE_Fortran_FLAGS} -g -O0 -fbounds-check -finit-real=nan -fbacktrace -coverage -Wall")
  # add flags only for compiling not linking!
  SET (SELF_COMPILE_FLAGS "")

ELSE()
  MESSAGE(SEND_ERROR "Unknown compiler")
ENDIF()

## ------- ROCm / GPU Related Settings ------- ##
# Set the GPU to compile for
set(GPU_TARGETS "gfx900" CACHE STRING "GPU targets to compile for")

# Search for rocm in common locations
list(APPEND CMAKE_PREFIX_PATH /opt/rocm/hip /opt/rocm)

# Find hip
find_package(hip REQUIRED)
## ------- ROCm / GPU Related Settings ------- ##

## --------- HDF5 Settings ----------- #
## --------- HDF5 Settings ----------- #

## --------- JSON-Fortran Settings ----------- #
find_package ( jsonfortran-${CMAKE_Fortran_COMPILER_ID} 8.3.0 REQUIRED )
## --------- JSON-Fortran Settings ----------- #

# Create library
add_library(self ...)
# Link with HIP
target_link_libraries(self hip::device)
