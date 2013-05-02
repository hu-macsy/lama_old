################################################
#   CMake example for simple C++ program       #
################################################

cmake_minimum_required( VERSION 2.8 )

## PROJECT NAME and supported LANGUAGES

project ( LibTest CXX )

################################################
#   Simple C++ program                         #
################################################

# default is building shared libraries, might be changed

set (BUILD_SHARED_LIBS ON CACHE BOOL

    "Build shared libraries instead of static ones.")

if ( BUILD_SHARED_LIBS )

    set ( LIBRARY_MODE SHARED )

    add_definitions( -DDLL_LINKING )

else ()

    set ( LIBRARY_MODE STATIC )

endif ()

add_library( rectangle ${LIBRARY_MODE} Rectangle.cpp Rectangle.hpp )

add_executable( example example.cpp )

target_link_libraries( example rectangle )

