cmake_minimum_required (VERSION 2.8 )

project( GramSchmidt )

# make sure that the default is a RELEASE

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

include( FindBLAS )

if ( BLAS_FOUND )
    message( STATUS "BLAS found, BLAS_LIBRARIES=${BLAS_LIBRARIES}" )
    add_executable( QR QR.cpp )
    target_link_libraries( QR ${BLAS_LIBRARIES} )
else( BLAS_FOUND )
    message( STATUS "BLAS not found" )
endif( BLAS_FOUND )

