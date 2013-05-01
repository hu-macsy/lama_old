################################################
#   CMake example for simple C++ program       #
################################################

cmake_minimum_required( VERSION 2.8 )

## PROJECT NAME and supported LANGUAGES

project ( LibTest CXX )

################################################
#   Simple C++ program                         #
################################################

set ( LIBRARY_MODE DYNAMIC )

add_library( rectangle ${LIBRARY_MODE} rectangle.cpp rectangle.hpp )
add_executable( example example.cpp )
target_link_libraries( example rectangle )

