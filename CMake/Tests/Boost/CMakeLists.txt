cmake_minimum_required (VERSION 2.8 )

project( BoostExample )

# make sure that the default is a RELEASE

set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../Modules" )

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set ( Boost_COMPONENTS thread )

set ( Boost_DEBUG TRUE )
set ( Boost_DETAILED_FAILURE_MSG TRUE )

find_package ( Boost COMPONENTS ${Boost_COMPONENTS} )
# find_package ( Boost COMPONENTS ${Boost_COMPONENTS} QUIET ) 

if ( NOT Boost_THREAD_LIBRARY )
    message ( FATAL_ERROR "Boost thread library not found" )
endif ( NOT Boost_THREAD_LIBRARY )

# Alternative on WIN32
# set ( Boost_USE_STATIC_LIBS ON )
# set ( Boost_USE_MULTITHREADED ON )
# set ( Boost_USE_STATIC_RUNTIME OFF )

if ( Boost_FOUND )
    message( STATUS "Boost found, BLAS_LIBRARIES=${BLAS_LIBRARIES}" )
    add_executable( example example.cpp )
    target_link_libraries( example ${Boost_THREAD_LIBRARY} )

# Alternative on WIN32
#   set ( LAMAPP_LINK_DIRECTORIES ${LAMAPP_LINK_DIRECTORIES} ${Boost_LIBRARY_DIRS} )
#   link_directories( ${LAMAPP_LINK_DIRECTORIES} )

else( Boost_FOUND )
    message( STATUS "Boost not found" )
endif( Boost_FOUND )

