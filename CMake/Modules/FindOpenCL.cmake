# - Try to find OPENCL
#
# OPENCL_INCLUDE_PATH and OPENCL_LIBRARY_PATH can be user defined in cmake call 
#
# Once done this will define
#  OPENCL_FOUND            - System has OPENCL
#  OPENCL_INCLUDE_DIRS     - The OPENCL include directories
#  OPENCL_LIBRARIES        - The libraries needed to use OPENCL

set( OPENCL_LIBRARY_SEARCH_PATHS /usr/local/ati/stream-sdk-v2.3-lnx64/lib/x86_64/ ${OPENCL_ROOT}/lib "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v4.0\\lib\\x64" "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v3.2\\lib\\x64" )
set( OPENCL_INCLUDE_SEARCH_PATHS /usr/local/ati/stream-sdk-v2.3-lnx64/include/ ${OPENCL_ROOT}/include )


# Name on CUDA,...    
find_library(OPENCL_LIBRARY NAMES OpenCL CL PATHS ${OPENCL_LIBRARY_SEARCH_PATHS}) 

find_path(OPENCL_INCLUDE_DIR CL/cl.h PATHS ${OPENCL_INCLUDE_SEARCH_PATHS} HINTS ${CUDA_INCLUDE_DIRS})
if(EXISTS ${OPENCL_LIBRARY})
	set(OPENCL_LIBRARIES ${OPENCL_LIBRARY})
else()
	message(STATUS "WARNING OpenCL library can not be found. Please define OPENCL_LIBRARY_PATH.")
endif()

if(EXISTS ${OPENCL_INCLUDE_DIR})
	set(OPENCL_INCLUDE_DIRS ${OPENCL_INCLUDE_DIR})
else()
	message(STATUS "WARNING OPENCL_INCLUDE_DIR=${OPENCL_INCLUDE_DIR} does not exist. Please define OPENCL_INCLUDE_PATH.")
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIET and REQUIRED arguments and set OPENCL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(OPENCL  DEFAULT_MSG
                                  OPENCL_LIBRARIES 
                                  OPENCL_INCLUDE_DIRS)

mark_as_advanced(OPENCL_INCLUDE_DIRS OPENCL_LIBRARIES OPENCL_LIBRARY OPENCL_INCLUDE_DIR )