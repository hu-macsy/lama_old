###
 # @file Variables.cmake
 #
 # @license
 # Copyright (c) 2009-2013
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in
 # all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 # @endlicense
 #
 # @brief Important CMake variable definitions
 # @author Jan Ecker
 # @date 16.04.2013
 # @since 1.0.0
###

# Check if verbose mode for CMAKE is selected
if ( DEFINED LAMA_CMAKE_VERBOSE AND LAMA_CMAKE_VERBOSE )
    set ( LAMA_FIND_PACKAGE_FLAGS )
else ()
    set ( LAMA_FIND_PACKAGE_FLAGS QUIET )
endif( DEFINED LAMA_CMAKE_VERBOSE AND LAMA_CMAKE_VERBOSE )

set ( CMAKE_SYSTEM_LIBRARY_PATH )
set ( LAMA_ROOT_DIR "${CMAKE_SOURCE_DIR}/.." )

# CMAKE configuration variable that guarantees adding rpath for installed
# libraries; very useful so that installed library can be used without 
# complex settings of LD_LIBRARY_PATH

set ( CMAKE_SKIP_BUILD_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set ( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

get_property ( FIND_LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS )
message ( STATUS "FindLib64: " ${FIND_LIB64} )

# Makefile outputs
set ( CMAKE_VERBOSE_MAKEFILE OFF )

# Define default library type as SHARED
if ( NOT DEFINED BUILD_SHARED_LIBS )
    set ( BUILD_SHARED_LIBS TRUE )
endif( NOT DEFINED BUILD_SHARED_LIBS )

## BUILDTYPE

# Choose Default CMAKE_BUILD_TYPE
if ( NOT CMAKE_BUILD_TYPE )
  # Can be: (RelWithDebInfo)
  set ( CMAKE_BUILD_TYPE Release CACHE STRING 
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
endif ( NOT CMAKE_BUILD_TYPE )

message ( STATUS "Build type is set to " ${CMAKE_BUILD_TYPE} )

## set cache entries
set ( LAMA_ADDITIONAL_LINK_LIBRARIES ${LAMA_ADDITIONAL_LINK_LIBRARIES} CACHE STRING "Additional libraries for linking, separated by ;" )
set ( LAMA_ADDITIONAL_LINK_FLAGS ${LAMA_ADDITIONAL_LINK_FLAGS} CACHE STRING "Additional link flags, separated by ;" )
mark_as_advanced ( LAMA_ADDITIONAL_LINK_LIBRARIES LAMA_ADDITIONAL_LINK_FLAGS )
