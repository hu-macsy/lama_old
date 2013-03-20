cmake_minimum_required ( VERSION 2.8 )

project( PapiExample )

# make sure that the default is a RELEASE

if ( NOT CMAKE_BUILD_TYPE )
   set( CMAKE_BUILD_TYPE RELEASE CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE )
endif ( NOT CMAKE_BUILD_TYPE )

# Find PAPI and set flags

set ( PAPI /home/brandes/local/papi-5.0.1 )

find_library( PAPI_LIBRARY
              NAMES papi
              HINTS ${PAPI}/lib )

set ( PAPI_LIBRARIES ${PAPI_LIBRARY} )

include_directories( ${PAPI}/include )

add_executable( example example.cpp )

target_link_libraries ( example ${PAPI_LIBRARIES} )

