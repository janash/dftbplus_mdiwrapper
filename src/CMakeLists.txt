# Locate MPI

find_package(MPI)
if(MPI_FOUND)
   include_directories(${MPI_INCLUDE_PATH})
else()
   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/STUBS_MPI/mpi.h ${CMAKE_CURRENT_BINARY_DIR}/STUBS_MPI/mpi.h COPYONLY)
endif()



# Link to MDI

set( MDI_LOCATION ${CMAKE_BINARY_DIR}/lib/mdi/MDI_Library/ )
link_directories( ${MDI_LOCATION} )
include_directories(${MDI_LOCATION})

link_directories( /repo/build/dftbplus/build/install/lib )
include_directories( /repo/build/dftbplus/build/install/include )




# Compile the driver

add_executable(dftbplus_mdiwrapper
               dftbplus_mdiwrapper.c)
target_link_libraries(dftbplus_mdiwrapper mdi
                      dftbplus
                      ${MPI_LIBRARIES})



# Ensure that MPI is properly linked

if(NOT MPI_FOUND)
   target_include_directories(dftbplus_mdiwrapper PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/STUBS_MPI/)
endif()
if(MPI_COMPILE_FLAGS)
   set_target_properties(dftbplus_mdiwrapper PROPERTIES
      COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
   set_target_properties(dftbplus_mdiwrapper PROPERTIES
      LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
