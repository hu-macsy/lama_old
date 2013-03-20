cmake_minimum_required (VERSION 2.8 )

project( OpenCLDevices )

# make sure that the default is a RELEASE

if ( NOT CMAKE_BUILD_TYPE )
   set( CMAKE_BUILD_TYPE RELEASE CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE )
endif ( NOT CMAKE_BUILD_TYPE )

# Find OpenCL and set flags

set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../../Modules" )

include( FindCUDA )
include( FindOpenCL )

if ( OPENCL_FOUND )
    message ( STATUS "OpenCL found at ${OPENCL_INCLUDE_DIRS}" )
else ( OPENCL_FOUND )
    message ( STATUS "OpenCL not found" )
endif ( OPENCL_FOUND )

include_directories( ${OPENCL_INCLUDE_DIRS} )

add_executable( OpenCLDevices main.cpp )

target_link_libraries( OpenCLDevices ${OPENCL_LIBRARY} )
