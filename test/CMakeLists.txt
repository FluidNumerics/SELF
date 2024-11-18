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

CMAKE_MINIMUM_REQUIRED(VERSION 3.21)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
set(make_env SELF_HOME=${CMAKE_BINARY_DIR} SELF_COMPILER=${SELF_COMPILER})

function (add_fortran_tests)

    foreach (TFILE ${ARGN})
        string (REGEX REPLACE ".f90$" "" TEST_NAME ${TFILE})
        add_executable (${TEST_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${TFILE})
        target_link_libraries(${TEST_NAME}
                    ${CMAKE_BINARY_DIR}/src/libself-static.a
                    ${FEQPARSE_LIBRARIES}
                    ${HDF5_LIBRARIES}
                    ${MPI_Fortran_LIBRARIES}
                    ${BACKEND_LIBRARIES})
                           
        target_include_directories(${TEST_NAME} PRIVATE 
                        ${CMAKE_BINARY_DIR}/include
                        ${FEQPARSE_INCLUDE_DIRS}
                        ${HDF5_INCLUDE_DIRS}
                        ${MPI_Fortran_INCLUDE_DIRS})

        add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
        install(TARGETS ${TEST_NAME} DESTINATION test)
    endforeach ()

endfunction ()

function (add_mpi_fortran_tests)

    foreach (TFILE ${ARGN})
        string (REGEX REPLACE ".f90$" "" TEST_NAME ${TFILE})
        add_executable (${TEST_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${TFILE})
        target_link_libraries(${TEST_NAME}
                    ${CMAKE_BINARY_DIR}/src/libself-static.a
                    ${FEQPARSE_LIBRARIES}
                    ${HDF5_LIBRARIES}
                    ${MPI_Fortran_LIBRARIES}
                    ${BACKEND_LIBRARIES})
                           
        target_include_directories(${TEST_NAME} PRIVATE 
                        ${CMAKE_BINARY_DIR}/include
                        ${FEQPARSE_INCLUDE_DIRS}
                        ${HDF5_INCLUDE_DIRS}
                        ${MPI_Fortran_INCLUDE_DIRS})

        add_test(NAME ${TEST_NAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SELF_MPIEXEC_NUMPROCS} ${SELF_MPIEXEC_OPTIONS} ${TEST_NAME})

        install(TARGETS ${TEST_NAME} DESTINATION test)
    endforeach ()

endfunction ()

add_fortran_tests (
    "scalargridinterp_1d_constant.f90"
    "scalarderivative_1d_constant.f90"
    "scalarboundaryinterp_1d_constant.f90"
    "scalargridinterp_2d_constant.f90"
    "scalargradient_2d_constant.f90"
    "scalarboundaryinterp_2d_constant.f90"
    "vectordivergence_2d_constant.f90"
    "vectordivergence_2d_linear.f90"
    "vectorboundaryinterp_2d_constant.f90"
    "scalargridinterp_3d_constant.f90"
    "scalargradient_3d_constant.f90"
    "scalarboundaryinterp_3d_constant.f90"
    "vectordivergence_3d_constant.f90"
    "vectorgridinterp_2d_constant.f90"
    "vectorgridinterp_3d_constant.f90"
    "vectorboundaryinterp_3d_constant.f90"
    "mesh2d_setup.f90"
    "mesh2d_uniformstructured.f90"
    "mesh3d_setup.f90"
    "mesh3d_uniformstructured.f90"
    "mappedscalarderivative_1d_constant.f90"
    "mappedscalarbrderivative_1d_constant.f90"
    "mappedscalardgderivative_1d_constant.f90"
    "mappedscalargradient_2d_constant.f90"
    "mappedscalargradient_2d_linear.f90"
    "mappedscalarbrgradient_2d_constant.f90"
    "mappedscalarbrgradient_2d_constant_structuredmesh.f90"
    "mappedscalarbrgradient_2d_linear.f90"
    "mappedscalarbrgradient_2d_linear_structuredmesh.f90"
    "mappedvectordivergence_2d_constant.f90"
    "mappedvectordivergence_2d_gausslobatto_constant.f90"
    "mappedvectordivergence_2d_linear_structuredmesh.f90"
    "mappedvectordgdivergence_2d_constant.f90"
    "mappedvectordgdivergence_2d_linear.f90"
    "mappedvectordgdivergence_2d_linear_structuredmesh.f90"
    "mappedscalargradient_3d_constant.f90"
    "mappedscalargradient_3d_linear.f90"
    "mappedscalarbrgradient_3d_constant.f90"
    "mappedscalarbrgradient_3d_constant_structuredmesh.f90"
    "mappedscalarbrgradient_3d_linear.f90"
    "mappedscalarbrgradient_3d_linear_structuredmesh.f90"
    "mappedvectordivergence_3d_constant.f90"
    "mappedvectordivergence_3d_linear.f90"
    "mappedvectordivergence_3d_gausslobatto_constant.f90"
    "mappedvectordgdivergence_3d_constant.f90"
    "mappedvectordgdivergence_3d_linear.f90"
    "mappedvectordgdivergence_3d_linear_structuredmesh.f90"
    "mappedvectordgdivergence_3d_linear_sideexchange.f90"
    "nulldgmodel1d.f90"
    "nulldgmodel2d.f90"
    "nulldgmodel2d_prescribed.f90"
    "nulldgmodel3d.f90"
    "nulldgmodel3d_prescribed.f90"
    "advection_diffusion_1d_euler.f90"
    "advection_diffusion_1d_euler_pickup.f90"
    "advection_diffusion_1d_rk2.f90"
    "advection_diffusion_1d_rk3.f90"
    "advection_diffusion_1d_rk4.f90"
    "burgers1d_constant.f90"
    "burgers1d_nonormalflow.f90"
    "burgers1d_prescribed.f90"
    "burgers1d_radiation.f90"
    "advection_diffusion_2d_euler.f90"
    "advection_diffusion_2d_rk2.f90"
    "advection_diffusion_2d_rk3.f90"
    "advection_diffusion_2d_rk3_pickup.f90"
    "advection_diffusion_2d_rk4.f90"
    "advection_diffusion_3d_euler.f90"
    "advection_diffusion_3d_rk2.f90"
    "advection_diffusion_3d_rk3.f90"
    "advection_diffusion_3d_rk4.f90"
    "linear_shallow_water_2d_constant.f90"
    "linear_shallow_water_2d_nonormalflow.f90"
    "linear_shallow_water_2d_radiation.f90"
    )

add_mpi_fortran_tests( "mappedvectordgdivergence_2d_linear_mpi.f90"
                       "mappedvectordgdivergence_2d_linear_sideexchange_mpi.f90"
                       "mappedvectordgdivergence_2d_linear_structuredmesh_mpi.f90"
                       "mappedscalarbrgradient_2d_linear_mpi.f90"
                       "mappedvectordgdivergence_3d_linear_mpi.f90"
                       "mappedvectordgdivergence_3d_linear_sideexchange_mpi.f90"
                       "mappedvectordgdivergence_3d_linear_structuredmesh_mpi.f90"
                       "mappedscalarbrgradient_3d_linear_mpi.f90"
                       "advection_diffusion_2d_rk3_mpi.f90"
                       "advection_diffusion_2d_rk3_pickup_mpi.f90"
                       "advection_diffusion_3d_rk3_mpi.f90"
                       "advection_diffusion_3d_rk3_pickup_mpi.f90" )