 ###
 # @file FindACML.cmake
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
 # @brief Find ACML
 # @author
 # @date 25.04.2013
###

# - Try to find ACML
#
#  ACML_INCLUDE_DIR and ACML_LIBRARY_PATH can be user defined in cmake call 
#
#  Once done this will define
#  ACML_FOUND - System has ACML
#  ACML_INCLUDE_DIRS - The ACML include directories
#  ACML_LIBRARIES - The libraries needed to use ACML

# FILE GLOB_RECURSE do not follow symlinks
cmake_policy ( SET CMP0009 NEW )

### Search for include path for all required acml header files

# If ACML_ROOT was defined use it.

if ( ACML_ROOT AND NOT ACML_INCLUDE_DIR )
    set ( ACML_INCLUDE_DIR ${ACML_ROOT}/include )
else ( ACML_ROOT AND NOT ACML_INCLUDE_DIR )

    if ( NOT ACML_INCLUDE_DIR AND NOT $ENV{ACML_ROOT} STREQUAL "" )
        set ( ACML_INCLUDE_DIR $ENV{ACML_ROOT}/include )
    elseif ( NOT ACML_INCLUDE_DIR AND NOT $ENV{ACMLROOT} STREQUAL "" )
	    set  ( ACML_INCLUDE_DIR $ENV{ACMLROOT}/include )
	#else( ) do nothing
    endif ( NOT ACML_INCLUDE_DIR AND NOT $ENV{ACML_ROOT} STREQUAL "" )

endif ( ACML_ROOT AND NOT ACML_INCLUDE_DIR )

if ( NOT DEFINED ACML_INCLUDE_DIR )
    # Search for includes by using common paths
    set ( ACML_HINTS ${ACML_ROOT} /usr/lib/ /usr/lib32/ /usr/lib64/ )
    
    foreach ( ACML_HINT ${ACML_HINTS} )
    
        if ( EXISTS ${ACML_HINT} )
            file ( GLOB_RECURSE ACML_SEARCH_DIR ${ACML_HINT}*/acml.h )
            
            if ( DEFINED ACML_SEARCH_DIR )
                list ( LENGTH ACML_SEARCH_DIR ACML_SEARCH_DIR_COUNT )
                
                # Found more than one include dir
                if ( ${ACML_SEARCH_DIR_COUNT} GREATER 1 ) 
                    list ( GET ACML_SEARCH_DIR 0 ACML_SEARCH_DIR )
                    message ( STATUS "Found multiple ACML installations. Will choose: " ${ACML_SEARCH_DIR} )
                endif( ${ACML_SEARCH_DIR_COUNT} GREATER 1 )
                
                get_filename_component ( ACML_INCLUDE_DIR "${ACML_SEARCH_DIR}" PATH CACHE )
                break()
            endif( DEFINED ACML_SEARCH_DIR )
            
        endif( EXISTS ${ACML_HINT} )
        
    endforeach ( ACML_HINT ${ACML_HINTS} )
endif( NOT DEFINED ACML_INCLUDE_DIR )

if ( NOT EXISTS ${ACML_INCLUDE_DIR} )
    message ( STATUS "WARNING No ACML include directory found. Please define ACML_INCLUDE_DIR or ACML_ROOT." )
    unset ( ACML_INCLUDE_DIR )
else ( NOT EXISTS ${ACML_INCLUDE_DIR} )
    ### Search for available library directories

    if ( NOT DEFINED ACML_LIBRARY_PATH )
        get_filename_component ( ACML_ROOT_PATH "${ACML_INCLUDE_DIR}" PATH CACHE )
        
        if ( EXISTS ${ACML_ROOT_PATH}/lib )
            set ( ACML_LIBRARY_PATH ${ACML_ROOT_PATH}/lib )
        else( EXISTS ${ACML_ROOT_PATH}/lib )
            message ( STATUS "WARNING ACML library path not found. Please define ACML_LIBRARY_PATH." )
        endif( EXISTS ${ACML_ROOT_PATH}/lib )
    endif( NOT DEFINED ACML_LIBRARY_PATH )
    
    ### Search for libraries libacml libacml_mp_dll
    
    if ( EXISTS ${ACML_LIBRARY_PATH} )
    
        if ( WIN32 )
            # search for libacml_mp_dll lib 
	    	message ( STATUS "Using ACML library path " ${ACML_LIBRARY_PATH} )
            find_library ( ACML_MP_DLL_LIBRARY libacml_mp_dll PATHS ${ACML_LIBRARY_PATH} )
        
            if( NOT EXISTS ${ACML_MP_DLL_LIBRARY} )
                message ( STATUS "WARNING ACML library acml_mp_dll not found with ACML_LIBRARY_PATH=${ACML_LIBRARY_PATH}." )
            endif( NOT EXISTS ${ACML_MP_DLL_LIBRARY} )
        
            set ( ACML_LIBRARIES ${ACML_MP_DLL_LIBRARY} )
        else ( WIN32 )
            #search for libacml lib
            message ( STATUS "Using ACML library path " ${ACML_LIBRARY_PATH} )
            find_library ( ACML_LIBRARY acml PATHS ${ACML_LIBRARY_PATH} )
        
            if( NOT EXISTS ${ACML_LIBRARY} )
                message ( STATUS "WARNING ACML library acml not found with ACML_LIBRARY_PATH=${ACML_LIBRARY_PATH}." )
            endif( NOT EXISTS ${ACML_LIBRARY} )
        
            set ( ACML_LIBRARIES ${ACML_LIBRARY} )
        endif ( WIN32 )
        
    else( EXISTS ${ACML_LIBRARY_PATH} )
        message ( STATUS "WARNING ACML libraries not found. ACML_LIBRARY_PATH=${ACML_LIBRARY_PATH} directory does not exist." )
    endif( EXISTS ${ACML_LIBRARY_PATH} )
    
    set ( ACML_INCLUDE_DIRS ${ACML_INCLUDE_DIR} )
endif( NOT EXISTS ${ACML_INCLUDE_DIR} )

include( FindPackageHandleStandardArgs )
# handle the QUIET and REQUIRED arguments and set MKL_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args ( ACML DEFAULT_MSG ACML_LIBRARIES ACML_INCLUDE_DIRS )

mark_as_advanced( ACML_INCLUDE_DIRS ACML_LIBRARIES)
