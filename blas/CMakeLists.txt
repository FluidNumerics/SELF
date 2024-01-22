file(GLOB SELF_FSRC "${CMAKE_CURRENT_SOURCE_DIR}/*.f*")

# Enable pre-processing for source code
set_source_files_properties(
  ${CMAKE_CURRENT_SOURCE_DIR}/scalargridinterp_1d_gpu_constant.f90
  PROPERTIES Fortran_PREPROCESS ON
)


add_executable (blas_test ${CMAKE_CURRENT_SOURCE_DIR}/scalargridinterp_1d_gpu_constant.f90)

#target_compile_options(self PRIVATE ${HIP_Fortran_LFLAGS})

target_link_libraries(blas_test 
                    ${CMAKE_BINARY_DIR}/src/libself-static.a
                    ${JSONFORTRAN_LIBRARIES}
                    ${FEQPARSE_LIBRARIES}
                    ${HDF5_LIBRARIES}
                    ${MPI_Fortran_LIBRARIES}
                    hip::device)

                           
target_include_directories(blas_test PRIVATE 
                        ${CMAKE_BINARY_DIR}/include
                        ${HDF5_INCLUDE_DIRS}/shared
                        ${MPI_Fortran_INCLUDE_DIRS})

install(TARGETS blas_test DESTINATION bin)
#target_compile_options(self PRIVATE)