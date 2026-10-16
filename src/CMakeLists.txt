
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

file(GLOB SELF_FSRC "${CMAKE_CURRENT_SOURCE_DIR}/*.f*")
if(SELF_ENABLE_GPU)
    file(GLOB SELF_BACKEND_FSRC "${CMAKE_CURRENT_SOURCE_DIR}/gpu/*.f*")
    file(GLOB SELF_BACKEND_CPPSRC "${CMAKE_CURRENT_SOURCE_DIR}/gpu/*.cpp*")

    if(SELF_ACTIVE_GPU_BACKEND STREQUAL "HIP")
      set_source_files_properties(${SELF_BACKEND_CPPSRC} PROPERTIES LANGUAGE HIP)
    else()
      set_source_files_properties(${SELF_BACKEND_CPPSRC} PROPERTIES LANGUAGE CUDA)
    endif()
elseif(SELF_ENABLE_APU)
    file(GLOB SELF_BACKEND_FSRC "${CMAKE_CURRENT_SOURCE_DIR}/apu/*.f*")
    # The APU backend reuses the hand-written device kernels from the gpu tree.
    set(SELF_BACKEND_CPPSRC "${CMAKE_CURRENT_SOURCE_DIR}/gpu/SELF_MatrixMultiply.cpp")

    if(SELF_ACTIVE_GPU_BACKEND STREQUAL "HIP")
      set_source_files_properties(${SELF_BACKEND_CPPSRC} PROPERTIES LANGUAGE HIP)
    else()
      set_source_files_properties(${SELF_BACKEND_CPPSRC} PROPERTIES LANGUAGE CUDA)
    endif()
else()
    file(GLOB SELF_BACKEND_FSRC "${CMAKE_CURRENT_SOURCE_DIR}/cpu/*.f*")
endif()

file(GLOB SELF_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Enable pre-processing for source code
set_source_files_properties(
  ${SELF_FSRC}
  PROPERTIES Fortran_PREPROCESS ON
)

set_source_files_properties(
  ${SELF_BACKEND_FSRC}
  PROPERTIES Fortran_PREPROCESS ON
)

set_source_files_properties(
  ${SELF_MODEL_FSRC}
  PROPERTIES Fortran_PREPROCESS ON
)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

add_library(self SHARED ${SELF_FSRC} ${SELF_BACKEND_CPPSRC} ${SELF_BACKEND_FSRC})
#set_target_properties(self PROPERTIES OUTPUT_NAME "self") 
# The absolute library and include paths that SELF's own build has always used
# are wrapped in BUILD_INTERFACE, so the build of libself and of the in-tree tests
# is unchanged. The installed package replaces them with imported targets that
# SELFConfig.cmake re-resolves on the consumer's machine (see cmake/SELFConfig.cmake.in).
target_link_libraries(self PUBLIC
                        "$<BUILD_INTERFACE:${FEQPARSE_LIBRARIES}>"
                        "$<INSTALL_INTERFACE:SELF::feqparse>"
                        HDF5::HDF5
                        "$<BUILD_INTERFACE:${MPI_Fortran_LIBRARIES}>"
                        "$<INSTALL_INTERFACE:MPI::MPI_Fortran>"
                        ${BACKEND_LIBRARIES})

target_include_directories(self PUBLIC
                        "$<BUILD_INTERFACE:${FEQPARSE_INCLUDE_DIRS}>"
                        "$<BUILD_INTERFACE:${HDF5_INCLUDE_DIRS}>"
                        "$<BUILD_INTERFACE:${MPI_Fortran_INCLUDE_DIRS}>")

# Consumers find SELF's .mod files here, in the build tree or the install tree.
target_include_directories(self INTERFACE
                        $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
                        $<INSTALL_INTERFACE:include>)

# SELF itself is compiled with these switches through the global CMAKE_Fortran_FLAGS
# set in the top level CMakeLists.txt. The installed target hands the same switches
# to consumers, so that code which tests them sees the configuration libself was
# built with. They are INSTALL_INTERFACE only, so the in-tree tests and examples,
# which already receive them from the global flags, are compiled exactly as before.
set(SELF_INTERFACE_DEFINITIONS "")
if(SELF_ENABLE_DOUBLE_PRECISION)
    list(APPEND SELF_INTERFACE_DEFINITIONS DOUBLE_PRECISION)
endif()
if(SELF_ENABLE_MULTITHREADING)
    list(APPEND SELF_INTERFACE_DEFINITIONS MULTITHREADING)
endif()
if(SELF_ENABLE_GPU)
    list(APPEND SELF_INTERFACE_DEFINITIONS ENABLE_GPU)
elseif(SELF_ENABLE_APU)
    list(APPEND SELF_INTERFACE_DEFINITIONS ENABLE_APU)
endif()
if(SELF_ACTIVE_GPU_BACKEND STREQUAL "HIP")
    list(APPEND SELF_INTERFACE_DEFINITIONS HAVE_HIP)
elseif(SELF_ACTIVE_GPU_BACKEND STREQUAL "CUDA")
    list(APPEND SELF_INTERFACE_DEFINITIONS HAVE_CUDA)
endif()
if(SELF_INTERFACE_DEFINITIONS)
    target_compile_definitions(self INTERFACE "$<INSTALL_INTERFACE:${SELF_INTERFACE_DEFINITIONS}>")
endif()

target_compile_options(self PUBLIC -fPIC)

set_target_properties(self PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(self PROPERTIES PUBLIC_HEADER ${SELF_HEADERS})

install(TARGETS self
        EXPORT SELFTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

# CMake package: find_package(SELF CONFIG) then link SELF::self.
include(CMakePackageConfigHelpers)
set(SELF_INSTALL_CMAKEDIR lib/cmake/SELF)

install(EXPORT SELFTargets
        NAMESPACE SELF::
        DESTINATION ${SELF_INSTALL_CMAKEDIR})

# Recorded in SELFConfig.cmake as the first place to look for FEQ-Parse.
get_filename_component(SELF_FEQPARSE_LIBRARY_DIR "${FEQPARSE_LIBRARIES}" DIRECTORY)

configure_package_config_file(${CMAKE_SOURCE_DIR}/cmake/SELFConfig.cmake.in
                              ${CMAKE_CURRENT_BINARY_DIR}/SELFConfig.cmake
                              INSTALL_DESTINATION ${SELF_INSTALL_CMAKEDIR})

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/SELFConfigVersion.cmake
                                 VERSION ${PROJECT_VERSION}
                                 COMPATIBILITY SameMajorVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SELFConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/SELFConfigVersion.cmake
        DESTINATION ${SELF_INSTALL_CMAKEDIR})