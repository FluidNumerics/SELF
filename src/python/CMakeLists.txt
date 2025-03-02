
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

file(GLOB SELF_PYINT_FSRC "${CMAKE_CURRENT_SOURCE_DIR}/*.f*")

# Enable pre-processing for source code
set_source_files_properties(
  ${SELF_PYINT_FSRC}
  PROPERTIES Fortran_PREPROCESS ON
)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

add_library(self_interface SHARED ${SELF_PYINT_FSRC})
#set_target_properties(self PROPERTIES OUTPUT_NAME "self") 
target_link_libraries(self_interface PUBLIC
                        self
                        ${FEQPARSE_LIBRARIES}
                        HDF5::HDF5
                        ${MPI_Fortran_LIBRARIES}
                        ${BACKEND_LIBRARIES}
                        ${JSONFORTRAN_LIBRARIES})

target_include_directories(self_interface PUBLIC
                        ${FEQPARSE_INCLUDE_DIRS}
                        ${HDF5_INCLUDE_DIRS}
                        ${MPI_Fortran_INCLUDE_DIRS}
                        ${JSONFORTRAN_INCLUDE_DIRS})

target_compile_options(self_interface PUBLIC -fPIC)    

set_target_properties(self_interface PROPERTIES LINKER_LANGUAGE Fortran)

install(TARGETS self_interface
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

# if(JSONFORTRAN_LIBRARIES)
# add_library(self-python SHARED ${SELF_FSRC} ${SELF_BACKEND_CPPSRC} ${SELF_BACKEND_FSRC})
# target_link_libraries(self PUBLIC
#                         ${FEQPARSE_LIBRARIES}
#                         HDF5::HDF5
#                         ${MPI_Fortran_LIBRARIES}
#                         ${BACKEND_LIBRARIES})

# target_include_directories(self PUBLIC
#                         ${FEQPARSE_INCLUDE_DIRS}
#                         ${HDF5_INCLUDE_DIRS}
#                         ${MPI_Fortran_INCLUDE_DIRS})