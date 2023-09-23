file(GLOB SELF_FSRC "${CMAKE_CURRENT_SOURCE_DIR}/*.f*")

# Enable pre-processing for source code
set_source_files_properties(
  ${CMAKE_CURRENT_SOURCE_DIR}/SELF_Main.f90
  PROPERTIES Fortran_PREPROCESS ON
)

add_executable (self ${CMAKE_CURRENT_SOURCE_DIR}/SELF_Main.f90)

target_link_libraries(self ${libjsonfortran} 
                           ${libfeqparse} 
                           ${HDF5_LIBRARIES} 
                           ${MPI_Fortran_LIBRARIES}
                           ${CMAKE_BINARY_DIR}/src/libself-static.a)

target_include_directories(self PRIVATE ${incjsonfortran}
                                ${incfeqparse} 
                                ${HDF5_INCLUDE_DIRS} 
                                ${MPI_Fortran_INCLUDE_DIRS}
                                ${CMAKE_BINARY_DIR}/src)