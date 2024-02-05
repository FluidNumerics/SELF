CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
set(make_env SELF_HOME=${CMAKE_BINARY_DIR} SELF_COMPILER=${SELF_COMPILER})


function (add_fortran_tests)

    foreach (TFILE ${ARGN})
        string (REGEX REPLACE ".f90$" "" TEST_NAME ${TFILE})
        add_executable (${TEST_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${TFILE})
        target_link_libraries(${TEST_NAME}
                    ${CMAKE_BINARY_DIR}/src/libself-static.a
                    ${JSONFORTRAN_LIBRARIES}
                    ${FEQPARSE_LIBRARIES}
                    ${HDF5_LIBRARIES}
                    ${MPI_Fortran_LIBRARIES}
                    hipfort::hip
                    hipfort::hipblas
                    hip::device)
                           
        target_include_directories(${TEST_NAME} PRIVATE 
                        ${CMAKE_BINARY_DIR}/include
                        ${HDF5_INCLUDE_DIRS}
                        ${MPI_Fortran_INCLUDE_DIRS}
                        ${HIPFORT_INCLUDE_DIRS})

        add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
        install(TARGETS ${TEST_NAME} DESTINATION example)
    endforeach ()

endfunction ()


add_fortran_tests (
    "scalargridinterp_1d_cpu_constant.f90"
    "scalargridinterp_1d_gpu_constant.f90"
    "scalarderivative_1d_cpu_constant.f90"
    "scalarderivative_1d_gpu_constant.f90"
    "scalarboundaryinterp_1d_cpu_constant.f90"
    "scalarboundaryinterp_1d_gpu_constant.f90"
    "scalargridinterp_2d_cpu_constant.f90"
    "scalargridinterp_2d_gpu_constant.f90"
    "scalargradient_2d_cpu_constant.f90"
    "scalargradient_2d_gpu_constant.f90"
    )


#add_fortran_test_executable (
#    testsuite


#    "scalarboundaryinterp_2d_cpu_constant.f90"
#    "scalarboundaryinterp_2d_gpu_constant.f90"
#    "scalarboundaryinterp_3d_cpu_constant.f90"
#    "scalarboundaryinterp_3d_gpu_constant.f90"
#    "vectorboundaryinterp_2d_cpu_constant.f90"
#    "vectorboundaryinterp_2d_gpu_constant.f90"
#    "vectorboundaryinterp_3d_cpu_constant.f90"
#    "vectorboundaryinterp_3d_gpu_constant.f90"


#    "scalargridinterp_3d_cpu_constant.f90"
#    "scalargridinterp_3d_gpu_constant.f90"
#    "vectorgridinterp_2d_cpu_constant.f90"
#    "vectorgridinterp_2d_gpu_constant.f90"
#    "vectorgridinterp_3d_cpu_constant.f90"
#    "vectorgridinterp_3d_gpu_constant.f90"

#    "scalargradient_3d_cpu_constant.f90"
#    "scalargradient_3d_gpu_constant.f90"
#    "vectordivergence_2d_cpu_constant.f90"
#    "vectordivergence_2d_gpu_constant.f90"
#    "vectordivergence_3d_cpu_constant.f90"
#    "vectordivergence_3d_gpu_constant.f90"
#    "mesh2d_setup.f90"
#    "mesh3d_setup.f90"
#    "mappedscalarderivative_1d_cpu_constant.f90"
#    "mappedscalarderivative_1d_gpu_constant.f90"
#    "mappedscalarbrderivative_1d_cpu_constant.f90"
#    "mappedscalarbrderivative_1d_gpu_constant.f90"
#    "mappedscalardgderivative_1d_cpu_constant.f90"
#    "mappedscalardgderivative_1d_gpu_constant.f90"
#    "mappedvectordivergence_2d_cpu_constant.f90"
#    "mappedvectordivergence_2d_gpu_constant.f90"
#    "mappedvectordivergence_2d_cpu_gausslobatto_constant.f90"
#    "mappedscalargradient_2d_cpu_constant.f90"
#    "mappedscalargradient_2d_gpu_constant.f90"
#    "mappedscalarbrgradient_2d_cpu_constant.f90"
#    "mappedscalarbrgradient_2d_gpu_constant.f90"
#    "mappedvectordivergence_3d_cpu_constant.f90"
#    "mappedvectordivergence_3d_gpu_constant.f90"
#    "mappedvectordivergence_3d_cpu_gausslobatto_constant.f90"
#    "mappedscalargradient_3d_cpu_constant.f90"
#    "mappedscalargradient_3d_gpu_constant.f90"
#    "mappedscalarbrgradient_3d_cpu_constant.f90"
#    "mappedscalarbrgradient_3d_gpu_constant.f90"
#    )
#