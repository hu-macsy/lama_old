 ###
 # @file FindVampireTrace.cmake
 #
 # @license
 # Copyright (c) 2013
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
 # @brief Find VampirTrace
 # @author
 # @date 25.04.2013
###

#
# - Find VampirTrace
#
# This module finds the VampirTrace include directory and library
#
# It sets the following variables:
#  VAMPIRTRACE_FOUND       - Set to false, or undefined, if VampirTrace isn't found.
#  VAMPIRTRACE_INCLUDE_DIR - The VampirTrace include directory.
#  VAMPIRTRACE_LIBRARY     - The VampirTrace library to link against.
	
set ( _vampirtrace_DEBUG false )

# If VAMPIRTRACE_ROOT or VT_ROOT is defined in the environment, use it.

if ( NOT VAMPIRTRACE_ROOT AND NOT $ENV{VAMPIRTRACE_ROOT} STREQUAL "" )
  set( VAMPIRTRACE_ROOT $ENV{VAMPIRTRACE_ROOT} )
endif ( NOT VAMPIRTRACE_ROOT AND NOT $ENV{VAMPIRTRACE_ROOT} STREQUAL "" )

if ( NOT VAMPIRTRACE_ROOT AND NOT $ENV{VT_ROOT} STREQUAL "" )
  set( VAMPIRTRACE_ROOT $ENV{VT_ROOT} )
endif ( NOT VAMPIRTRACE_ROOT AND NOT $ENV{VT_ROOT} STREQUAL "" )

# If VAMPIRTRACE_INCLUDEDIR or VT_INC is defined in the environment, use it.
if ( NOT VAMPIRTRACE_INCLUDEDIR AND NOT $ENV{VAMPIRTRACE_INCLUDEDIR} STREQUAL "" )
  set( VAMPIRTRACE_INCLUDEDIR $ENV{VAMPIRTRACE_INCLUDEDIR} )
endif ( NOT VAMPIRTRACE_INCLUDEDIR AND NOT $ENV{VAMPIRTRACE_INCLUDEDIR} STREQUAL "" )

if ( NOT VAMPIRTRACE_INCLUDEDIR AND NOT $ENV{VT_INC} STREQUAL "" )
  set( VAMPIRTRACE_INCLUDEDIR $ENV{VT_INC} )
endif ( NOT VAMPIRTRACE_INCLUDEDIR AND NOT $ENV{VT_INC} STREQUAL "" )

# If VAMPIRTRACE_LIBRARYDIR was defined in the environment, use it.
if( NOT VAMPIRTRACE_LIBRARYDIR AND NOT $ENV{VAMPIRTRACE_LIBRARYDIR} STREQUAL "" )
  set( VAMPIRTRACE_LIBRARYDIR $ENV{VAMPIRTRACE_LIBRARYDIR} )
endif( NOT VAMPIRTRACE_LIBRARYDIR AND NOT $ENV{VAMPIRTRACE_LIBRARYDIR} STREQUAL "" )

if( NOT VAMPIRTRACE_LIBRARYDIR AND NOT $ENV{VT_LIB} STREQUAL "" )
  set( VAMPIRTRACE_LIBRARYDIR $ENV{VT_LIB} )
endif( NOT VAMPIRTRACE_LIBRARYDIR AND NOT $ENV{VT_LIB} STREQUAL "" )

if( VAMPIRTRACE_ROOT )
    if ( _vampirtrace_DEBUG)
       message(STATUS "VAMPIRTRACE_ROOT = ${VAMPIRTRACE_ROOT}")
    endif ( _vampirtrace_DEBUG)

    set(_vampirtrace_INCLUDE_SEARCH_DIRS ${VAMPIRTRACE_ROOT}/include )
    set(_vampirtrace_LIBRARY_SEARCH_DIRS ${VAMPIRTRACE_ROOT}/lib )
else( VAMPIRTRACE_ROOT )
    if ( _vampirtrace_DEBUG)
       message(STATUS "VAMPIRTRACE_ROOT not set")
    endif ( _vampirtrace_DEBUG)
endif( VAMPIRTRACE_ROOT )

if( VAMPIRTRACE_INCLUDEDIR )
  file(TO_CMAKE_PATH ${VAMPIRTRACE_INCLUDEDIR} VAMPIRTRACE_INCLUDEDIR)
  SET(_vampirtrace_INCLUDE_SEARCH_DIRS
    ${VAMPIRTRACE_INCLUDEDIR} )
endif( VAMPIRTRACE_INCLUDEDIR )

if( VAMPIRTRACE_LIBRARYDIR )
  file(TO_CMAKE_PATH ${VAMPIRTRACE_LIBRARYDIR} VAMPIRTRACE_LIBRARYDIR)
  SET(_vampirtrace_LIBRARY_SEARCH_DIRS
    ${VAMPIRTRACE_LIBRARYDIR} )
endif( VAMPIRTRACE_LIBRARYDIR )

# now find VAMPIRTRACE_INCLUDE_DIR

if ( _vampirtrace_DEBUG)
  message(STATUS "search include dirs for vampirtrace = ${_vampirtrace_INCLUDE_SEARCH_DIRS}")
endif ( _vampirtrace_DEBUG)

find_path(VAMPIRTRACE_INCLUDE_DIR 
          NAMES vampirtrace/vt_user.h
          HINTS ${_vampirtrace_INCLUDE_SEARCH_DIRS})

if ( _vampirtrace_DEBUG)
  message(STATUS "include dir for vampirtrace = ${VAMPIRTRACE_INCLUDE_DIR}")
endif ( _vampirtrace_DEBUG)

# now find VAMPIRTRACE_LIBRARY

if ( _vampirtrace_DEBUG)
  message(STATUS "search library dirs for vampirtrace = ${_vampirtrace_LIBRARY_SEARCH_DIRS}")
endif ( _vampirtrace_DEBUG)

find_library(VAMPIRTRACE_OTF_LIBRARY 
             NAMES otf
             HINTS ${_vampirtrace_LIBRARY_SEARCH_DIRS})

find_library(VAMPIRTRACE_VT_LIBRARY 
             NAMES vt-hyb vt.hyb
             HINTS ${_vampirtrace_LIBRARY_SEARCH_DIRS})

if (VAMPIRTRACE_VT_LIBRARY AND VAMPIRTRACE_OTF_LIBRARY)
    set( VAMPIRTRACE_LIBRARIES ${VAMPIRTRACE_VT_LIBRARY} ${VAMPIRTRACE_OTF_LIBRARY} )
endif()

if ( _vampirtrace_DEBUG)
  message(STATUS "otf library for vampirtrace = ${VAMPIRTRACE_OTF_LIBRARY}")
  message(STATUS "vt library for vampirtrace = ${VAMPIRTRACE_VT_LIBRARY}")
endif ( _vampirtrace_DEBUG)

if (VAMPIRTRACE_INCLUDE_DIR AND VAMPIRTRACE_OTF_LIBRARY AND VAMPIRTRACE_VT_LIBRARY)
   SET(VAMPIRTRACE_FOUND TRUE)
endif ()

if (VAMPIRTRACE_FOUND)
   # show which VampirTrace was found only if not quiet
   if (NOT VampirTrace_FIND_QUIETLY)
      message(STATUS "Found VampirTrace: ${VAMPIRTRACE_OTF_LIBRARY} ${VAMPIRTRACE_VT_LIBRARY}")
   endif (NOT VampirTrace_FIND_QUIETLY)
else (VAMPIRTRACE_FOUND)
   # fatal error if VampirTrace is required but not found
   if (VampirTrace_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find VampirTrace")
   endif (VampirTrace_FIND_REQUIRED)
endif (VAMPIRTRACE_FOUND)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VAMPIRTRACE  DEFAULT_MSG  VAMPIRTRACE_OTF_LIBRARY VAMPIRTRACE_VT_LIBRARY VAMPIRTRACE_INCLUDE_DIR)

mark_as_advanced(VAMPIRTRACE_INCLUDE_DIR VAMPIRTRACE_OTF_LIBRARY VAMPIRTRACE_VT_LIBRARY)
