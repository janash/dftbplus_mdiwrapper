# Locate MPI

find_package(MPI)
if(MPI_Fortran_FOUND)
   include_directories(${MPI_Fortran_INCLUDE_PATH})
else()
   add_subdirectory(STUBS_MPI)
endif()



# Link to MDI

set( MDI_LOCATION ${CMAKE_BINARY_DIR}/lib/mdi/MDI_Library/ )
link_directories( ${MDI_LOCATION} )
include_directories(${MDI_LOCATION})

link_directories( /repo/build/dftbplus/build/install/lib )
include_directories( /repo/build/dftbplus/build/install/include/dftbplus/modfiles )




# Compile the engine as an exectable

add_executable(dftbplus_mdiwrapper
               dftbplus_mdiwrapper.F90
               dftbplus_engine.f90)


# Compile the engine as a shared library

add_library(dftbplus_mdiwrapper_lib SHARED
            dftbplus_mdiwrapper.f90
            dftbplus_engine.f90)
target_link_libraries(dftbplus_mdiwrapper_lib mdi)
set_target_properties(dftbplus_mdiwrapper_lib PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    OUTPUT_NAME dftbplus_mdiwrapper)



# Ensure that MPI is properly linked

if(NOT MPI_Fortran_FOUND)
   target_include_directories(dftbplus_mdiwrapper PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/STUBS_MPI/)
   target_include_directories(dftbplus_mdiwrapper_lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/STUBS_MPI/)
   target_link_libraries(dftbplus_mdiwrapper mpi)
   target_link_libraries(dftbplus_mdiwrapper_lib mpi)
else()
   target_include_directories(dftbplus_mdiwrapper PUBLIC
      ${MPI_Fortran_INCLUDE_PATH} )
   target_include_directories(dftbplus_mdiwrapper_lib PUBLIC
      ${MPI_Fortran_INCLUDE_PATH} )
   target_link_libraries(dftbplus_mdiwrapper ${MPI_Fortran_LIBRARIES})
   target_link_libraries(dftbplus_mdiwrapper_lib ${MPI_Fortran_LIBRARIES})
endif()
if(MPI_Fortran_COMPILE_FLAGS)
   set_target_properties(dftbplus_mdiwrapper PROPERTIES
      COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
   set_target_properties(dftbplus_mdiwrapper_lib PROPERTIES
      COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}")
endif()
if(MPI_Fortran_LINK_FLAGS)
   set_target_properties(dftbplus_mdiwrapper PROPERTIES
      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
   set_target_properties(dftbplus_mdiwrapper_lib PROPERTIES
      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
endif()



target_link_libraries(dftbplus_mdiwrapper mdi
                      dftbplus
                      ${MPI_LIBRARIES})
target_link_libraries(dftbplus_mdiwrapper_lib mdi
                      dftbplus
                      ${MPI_LIBRARIES})
