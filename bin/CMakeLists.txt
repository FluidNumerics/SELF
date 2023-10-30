file(GLOB SELF_FSRC "${CMAKE_CURRENT_SOURCE_DIR}/*.f*")

# Enable pre-processing for source code
set_source_files_properties(
  ${CMAKE_CURRENT_SOURCE_DIR}/SELF_Main.f90
  PROPERTIES Fortran_PREPROCESS ON
)


add_executable (self ${CMAKE_CURRENT_SOURCE_DIR}/SELF_Main.f90)

#target_compile_options(self PRIVATE ${HIP_Fortran_LFLAGS})

target_link_libraries(self 
                    ${CMAKE_BINARY_DIR}/src/libself-static.a
                    ${JSONFORTRAN_LIBRARIES}
                    ${FEQPARSE_LIBRARIES}
                    ${HDF5_LIBRARIES}
                    ${MPI_Fortran_LIBRARIES}
                    hip::device)

                           
target_include_directories(self PRIVATE 
                        ${CMAKE_BINARY_DIR}/include
                        ${HDF5_INCLUDE_DIRS}/shared
                        ${MPI_Fortran_INCLUDE_DIRS})

install(TARGETS self DESTINATION bin)
#target_compile_options(self PRIVATE)
