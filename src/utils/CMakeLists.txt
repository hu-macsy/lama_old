
### set include directories ###
if ( MPI_FOUND )
    include_directories ( ${MPI_INCLUDE_PATH} )
endif ( MPI_FOUND )

include_directories ( .. ../lama ${Boost_INCLUDE_DIRS} )

### set link directories ###
if ( Boost_INCLUDE_DIR )
    if ( Boost_LIBRARY_DIRS )
        set ( LINK_PATH ${Boost_LIBRARY_DIRS} )
    endif ( Boost_LIBRARY_DIRS )
endif ( Boost_INCLUDE_DIR )

link_directories ( ${LINK_PATH} )

include_directories ( ${CUDA_INCLUDE_DIRS} )

### add executable ###

# vector converter
add_executable ( MatrixMultTest MatrixMultTest.cpp )


target_link_libraries ( MatrixMultTest ama )