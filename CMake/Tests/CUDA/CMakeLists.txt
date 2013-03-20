cmake_minimum_required (VERSION 2.8 )

project( CUDADevices )

# make sure that the default is a RELEASE

if ( NOT CMAKE_BUILD_TYPE )
   set( CMAKE_BUILD_TYPE RELEASE CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE )
endif ( NOT CMAKE_BUILD_TYPE )

# Find CUDA and set flags

include( FindCUDA )

cuda_add_executable( CUDADevices main.cu )

