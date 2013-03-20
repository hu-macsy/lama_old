cmake_minimum_required (VERSION 2.8 )

# Attention: project( GramSchmidt CXX ) does not work

project( GramSchmidt )

# make sure that the default is a RELEASE

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

include( FindLAPACK )

if ( LAPACK_FOUND )
    message( STATUS "LAPACK found, LAPACK_LIBRARIES=${LAPACK_LIBRARIES}" )
    add_executable( QR QR.cpp )
    target_link_libraries( QR ${LAPACK_LIBRARIES} )
else( LAPACK_FOUND )
    message( STATUS "LAPACK not found" )
endif( LAPACK_FOUND )

