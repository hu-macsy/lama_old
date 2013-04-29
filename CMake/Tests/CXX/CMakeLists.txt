################################################
#   CMake example for simple C++ program       #
################################################

cmake_minimum_required( VERSION 2.8 )

## PROJECT NAME and supported LANGUAGES

project ( CXXTest CXX )

################################################
#   Simple C++ program                         #
################################################

add_executable( example example.cpp )


