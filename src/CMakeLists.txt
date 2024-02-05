# file(GLOB hipfort_hipblas_sources hipfort_hipblas_*.f90)
# add_library(hipfort_hipblas STATIC ${hipfort_hipblas_sources})

file(GLOB SELF_FSRC "${CMAKE_CURRENT_SOURCE_DIR}/*.f*")
file(GLOB SELF_HIPSRC "${CMAKE_CURRENT_SOURCE_DIR}/*.hip")

# Enable pre-processing for source code
set_source_files_properties(
  ${SELF_FSRC}
  PROPERTIES Fortran_PREPROCESS ON
)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

add_library(self SHARED ${SELF_FSRC} ${SELF_HIPSRC})
add_library(self-static STATIC ${SELF_FSRC} ${SELF_HIPSRC})

target_link_libraries(self PRIVATE 
                        ${FEQPARSE_LIBRARIES}
                        ${HDF5_LIBRARIES}
                        ${MPI_Fortran_LIBRARIES}
                        hipfort::hip)

target_include_directories(self PRIVATE 
                        ${FEQPARSE_INCLUDE_DIRS}
                        ${HDF5_INCLUDE_DIRS}
                        ${MPI_Fortran_INCLUDE_DIRS})
                        
target_link_libraries(self-static PRIVATE 
                        ${FEQPARSE_LIBRARIES}
                        ${HDF5_LIBRARIES}
                        ${MPI_Fortran_LIBRARIES}
                        hipfort::hip)

target_include_directories(self-static PRIVATE 
                        ${FEQPARSE_INCLUDE_DIRS}
                        ${HDF5_INCLUDE_DIRS}
                        ${MPI_Fortran_INCLUDE_DIRS})
                        
target_compile_options(self-static PRIVATE -fPIC)

install(TARGETS self-static
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

set_target_properties(self-static PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(self PROPERTIES LINKER_LANGUAGE Fortran)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)