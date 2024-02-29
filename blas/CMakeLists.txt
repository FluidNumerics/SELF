CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
set(make_env SELF_HOME=${CMAKE_BINARY_DIR} SELF_COMPILER=${SELF_COMPILER})


function (add_fortran_examples)

    foreach (TFILE ${ARGN})
        string (REGEX REPLACE ".f90$" "" EXAMPLE_NAME ${TFILE})
        add_executable (${EXAMPLE_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/${TFILE})
        target_link_libraries(${EXAMPLE_NAME}
                    ${CMAKE_BINARY_DIR}/src/libself-static.a
                    ${JSONFORTRAN_LIBRARIES}
                    ${FEQPARSE_LIBRARIES}
                    ${HDF5_LIBRARIES}
                    ${MPI_Fortran_LIBRARIES}
                    hip::device
                    roc::hipblas)
                           
        target_include_directories(${EXAMPLE_NAME} PRIVATE 
                        ${CMAKE_BINARY_DIR}/include
                        ${HDF5_INCLUDE_DIRS}/shared
                        ${MPI_Fortran_INCLUDE_DIRS})

        add_test(NAME ${EXAMPLE_NAME} COMMAND ${EXAMPLE_NAME})
        install(TARGETS ${EXAMPLE_NAME} DESTINATION example)
    endforeach ()

endfunction ()

add_fortran_examples (
    # Level 1
    # "axpy_test.f90"
    # Level 2
    #"test_gemv.f90"
    #"test_gemvstridedbatched.f90"
    ##### gemvXStridedBatched comparisons
    "gemvstridedbatched_op_n_32.f90"
    "gemvstridedbatched_op_t_32.f90"
    "gemvstridedbatched_op_n_64.f90"
    "gemvstridedbatched_op_t_64.f90"
    # Level 3
    #"test_gemm.f90"
    # misc
    # "example_hipSaxpy.f90"
    # "example_hipDaxpy.f90"
)