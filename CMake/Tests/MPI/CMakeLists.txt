cmake_minimum_required( VERSION 2.8 )

## PROJECT NAME and supported LANGUAGES

project ( MPITest CXX )

################################################
#   MPI stuff                                  #
################################################

find_package( MPI REQUIRED )

include_directories( ${MPI_INCLUDE_PATH} )

add_executable( MPITest cpi.cpp )

target_link_libraries( MPITest ${MPI_LIBRARIES} )


