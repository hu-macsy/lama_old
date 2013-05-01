cmake_minimum_required (VERSION 2.8 )

project( BoostExample )

# make sure that the default is a RELEASE

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set ( Boost_COMPONENTS thread )

set ( Boost_DEBUG TRUE )
set ( Boost_DETAILED_FAILURE_MSG TRUE )

# Alternative on WIN32
# set ( Boost_USE_STATIC_LIBS ON )
# set ( Boost_USE_MULTITHREADED ON )
# set ( Boost_USE_STATIC_RUNTIME OFF )

# Note: set (environment) variable BOOST_ROOT to find installed Boost

find_package ( Boost COMPONENTS ${Boost_COMPONENTS} )

if ( NOT Boost_THREAD_LIBRARY )
    message ( FATAL_ERROR "Boost thread library not found" )
endif ( NOT Boost_THREAD_LIBRARY )

if ( Boost_FOUND )
    message( STATUS "Boost found, Boost_THREAD_LIBRARY=${Boost_THREAD_LIBRARY}" )
    add_executable( example example.cpp )
    target_link_libraries( example ${Boost_THREAD_LIBRARY} )

# Alternative on WIN32
#   link_directories( ${Boost_LIBRARY_DIRS} )

else( Boost_FOUND )
    message( STATUS "Boost not found" )
endif( Boost_FOUND )

