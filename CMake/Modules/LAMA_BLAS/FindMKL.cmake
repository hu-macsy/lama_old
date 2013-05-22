 ###
 # @file FindMKL.cmake
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
 # @brief Find MKL
 # @author
 # @date 25.04.2013
 # $Id$
###

# - Try to find MKL
#
#  MKL_INCLUDE_DIR and MKL_LIBRARY_PATH can be user defined in cmake call 
#
#  Once done this will define
#  MKL_FOUND - System has MKL
#  MKL_INCLUDE_DIRS - The MKL include directories
#  MKL_LIBRARIES - The libraries needed to use MKL

# FILE GLOB_RECURSE do not follow symlinks
cmake_policy ( SET CMP0009 NEW )

### Search for include path for all required mkl header files

# If MKL_ROOT was defined in the environment, use it.
if ( MKL_ROOT AND NOT MKL_INCLUDE_DIR )
    set ( MKL_INCLUDE_DIR ${MKL_ROOT}/include )
else ( MKL_ROOT AND NOT MKL_INCLUDE_DIR )

    if ( NOT MKL_INCLUDE_DIR AND NOT $ENV{MKL_ROOT} STREQUAL "" )
        set ( MKL_INCLUDE_DIR $ENV{MKL_ROOT}/include )
    elseif ( NOT MKL_INCLUDE_DIR AND NOT $ENV{MKLROOT} STREQUAL "" )
	    set  ( MKL_INCLUDE_DIR $ENV{MKLROOT}/include )
	#else( ) do nothing
    endif ( NOT MKL_INCLUDE_DIR AND NOT $ENV{MKL_ROOT} STREQUAL "" )
    
endif ( MKL_ROOT AND NOT MKL_INCLUDE_DIR )

if ( NOT DEFINED MKL_INCLUDE_DIR )
    # Search for includes by using common paths
    set ( MKL_HINTS ${MKL_ROOT_PATH} /opt/intel/ /usr/lib/ /usr/lib32/ /usr/lib64/ )
    
    foreach ( MKL_HINT ${MKL_HINTS} )
    
        if ( EXISTS ${MKL_HINT} )
            file ( GLOB_RECURSE MKL_SEARCH_DIR ${MKL_HINT}*/mkl.h )
            
            if ( DEFINED MKL_SEARCH_DIR )
                list ( LENGTH MKL_SEARCH_DIR MKL_SEARCH_DIR_COUNT )
                
                # Found more than one include dir
                if ( ${MKL_SEARCH_DIR_COUNT} GREATER 1 ) 
                    list ( GET MKL_SEARCH_DIR 0 MKL_SEARCH_DIR )
                    message ( STATUS "Found multiple MKL installations. Will choose: " ${MKL_SEARCH_DIR} )
                endif( ${MKL_SEARCH_DIR_COUNT} GREATER 1 )
                
                get_filename_component ( MKL_INCLUDE_DIR "${MKL_SEARCH_DIR}" PATH )
                break()
            endif( DEFINED MKL_SEARCH_DIR )
            
        endif( EXISTS ${MKL_HINT} )
        
    endforeach ( MKL_HINT ${MKL_HINTS} )
endif( NOT DEFINED MKL_INCLUDE_DIR )

if ( NOT EXISTS ${MKL_INCLUDE_DIR} )
    message ( STATUS "WARNING No MKL include directory found. Please define MKL_INCLUDE_DIR or MKL_ROOT." )
    unset ( MKL_INCLUDE_DIR )
else ( NOT EXISTS ${MKL_INCLUDE_DIR} )
    ### Search for available library directories

    if ( NOT DEFINED MKL_LIBRARY_PATH )
        get_filename_component ( MKL_ROOT_PATH "${MKL_INCLUDE_DIR}" PATH )
        
        if ( LAMA_CMAKE_VERBOSE )
            message ( STATUS "MKL_ROOT_PATH = ${MKL_ROOT_PATH}" )
        endif ( LAMA_CMAKE_VERBOSE )
        
        if ( NOT DEFINED MKL_Is64 )
            get_cmake_property ( MKL_Is64 FIND_LIBRARY_USE_LIB64_PATHS )
        endif( NOT DEFINED MKL_Is64 )
        
        # search lib suffix
        if ( MKL_Is64 OR CMAKE_CL_64 ) #64 Bit libraries
            #define possible subdirectories
            set ( MKL_LIBRARY_PATH_SUFFIXES em64t intel64 )
        else( MKL_Is64 OR CMAKE_CL_64 ) #32 Bit libraries
            #define possible subdirectories
            set ( MKL_LIBRARY_PATH_SUFFIXES 32 ia32 )
        endif( MKL_Is64 OR CMAKE_CL_64 )
        
        if ( EXISTS ${MKL_ROOT_PATH}/lib )
            set( MKL_LIBRARY_PATH ${MKL_ROOT_PATH}/lib )
        else( EXISTS ${MKL_ROOT_PATH}/lib )
            message ( STATUS "WARNING MKL library path not found. Please define MKL_LIBRARY_PATH." )
        endif( EXISTS ${MKL_ROOT_PATH}/lib )
    endif( NOT DEFINED MKL_LIBRARY_PATH )

    if ( LAMA_CMAKE_VERBOSE )
        message( STATUS "MKL_LIBRARY_PATH = ${MKL_LIBRARY_PATH}, MKL_LIBRARY_PATH_SUFFIXES = ${MKL_LIBRARY_PATH_SUFFIXES}" )
    endif ( LAMA_CMAKE_VERBOSE )
    
    ### Search for libraries mkl_gnu_thread, mkl_intel_thread and mkl_core
    
    if ( EXISTS ${MKL_LIBRARY_PATH} )
        # search for mkl_intel_lp64 lib 
        find_library ( MKL_LIBRARY_LP64 mkl_intel_lp64 PATHS ${MKL_LIBRARY_PATH} PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )

        if ( NOT EXISTS ${MKL_LIBRARY_LP64} )
            message ( STATUS "WARNING MKL library mkl_intel_lp64 not found with MKL_LIBRARY_PATH=${MKL_LIBRARY_PATH}." )
        endif( NOT EXISTS ${MKL_LIBRARY_LP64} )
        set(MKL_LIBRARIES  ${MKL_LIBRARY_LP64} )
        
        if ( NOT EXISTS ${MKL_LIBRARY_BLACS} )
            message ( STATUS "WARNING MKL library mkl_blacs not found with MKL_LIBRARY_PATH=${MKL_LIBRARY_PATH}." )
        endif( NOT EXISTS ${MKL_LIBRARY_BLACS} )
        
        # search for mkl_core lib 
        find_library ( MKL_LIBRARY_CORE mkl_core PATHS ${MKL_LIBRARY_PATH} PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )
        if ( NOT EXISTS ${MKL_LIBRARY_CORE} )
            message ( STATUS "WARNING MKL library mkl_core not found with MKL_LIBRARY_PATH=${MKL_LIBRARY_PATH}." )
        endif ( NOT EXISTS ${MKL_LIBRARY_CORE} )
        
        # search for mkl_thread lib
            
        # search for gnu compiler libs
        if ( CMAKE_COMPILER_IS_GNUCC )
            find_library ( MKL_LIBRARY_GNU mkl_gnu_thread PATHS ${MKL_LIBRARY_PATH} PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )
            if ( NOT EXISTS ${MKL_LIBRARY_GNU} )
                message ( STATUS "WARNING MKL library mkl_gnu_thread not found with MKL_LIBRARY_PATH=${MKL_LIBRARY_PATH}." )
            endif( NOT EXISTS ${MKL_LIBRARY_GNU} )
            list ( APPEND MKL_LIBRARIES ${MKL_LIBRARY_GNU} )
         endif( CMAKE_COMPILER_IS_GNUCC )
        
         # search for intel compiler libs
         if ( CMAKE_C_COMPILER_ID MATCHES Intel )
            find_library ( MKL_LIBRARY_INTEL mkl_intel_thread PATHS ${MKL_LIBRARY_PATH} PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )
            if ( NOT EXISTS ${MKL_LIBRARY_INTEL} )
                message ( STATUS "WARNING MKL library mkl_intel_thread not found with MKL_LIBRARY_PATH=${MKL_LIBRARY_PATH}." )
            endif()
            list ( APPEND MKL_LIBRARIES ${MKL_LIBRARY_INTEL} )
         endif( CMAKE_C_COMPILER_ID MATCHES Intel )
        
         if ( ${CMAKE_GENERATOR} MATCHES "Visual Studio" )
            find_library ( MKL_LIBRARY_INTEL mkl_intel_thread PATHS ${MKL_LIBRARY_PATH} PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )
            if ( NOT EXISTS ${MKL_LIBRARY_INTEL} )
                message ( STATUS "WARNING MKL library mkl_intel_thread not found with MKL_LIBRARY_PATH=${MKL_LIBRARY_PATH}." )
            endif ( NOT EXISTS ${MKL_LIBRARY_INTEL} )
            list ( APPEND MKL_LIBRARIES ${MKL_LIBRARY_INTEL} )
         endif( ${CMAKE_GENERATOR} MATCHES "Visual Studio" )

        # search for gfortran. Required by older MKL versions
        if ( NOT WIN32 )
            find_library ( GFORTRAN_LIBRARY gfortran HINTS ${GFORTRAN_LIBRARY_PATH} )
            mark_as_advanced( FORCE GFORTRAN_LIBRARY )
        
            if ( EXISTS ${GFORTRAN_LIBRARY} )
                list ( APPEND MKL_LIBRARIES ${GFORTRAN_LIBRARY} )
            else ( EXISTS ${GFORTRAN_LIBRARY} )
                message ( STATUS "WARNING Library gfortran not found. Required by some older MKL libraries. Please define GFORTRAN_LIBRARY_PATH." )
            endif ( EXISTS ${GFORTRAN_LIBRARY} )
        endif ( NOT WIN32 )
        
        # conclude libs
        list ( APPEND MKL_LIBRARIES ${MKL_LIBRARY_CORE} )
        
        if ( WIN32 )
            find_library ( INTEL_IOMP_LIBRARY libiomp5md PATHS ${MKL_LIBRARY_PATH} ${MKL_ROOT_PATH}/../compiler/lib PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )
            if ( NOT EXISTS ${INTEL_IOMP_LIBRARY} )
                message ( STATUS "WARNING Intel OMP library iomp5md not found with INTEL_COMPILER_PATH=${MKL_ROOT_PATH}/../compiler/lib." )
            endif ( NOT EXISTS ${INTEL_IOMP_LIBRARY} )
            list ( APPEND MKL_LIBRARIES ${INTEL_IOMP_LIBRARY} )
        endif ( WIN32 )
        
        if ( LAMA_CMAKE_VERBOSE )
            message ( STATUS "Found MKL Libraries: ${MKL_LIBRARIES}." )
            message ( STATUS "Found MKL PLibraries: ${MKL_PLIBRARIES}." )
        endif ( LAMA_CMAKE_VERBOSE )
    else ( EXISTS ${MKL_LIBRARY_PATH} )
        message ( STATUS "WARNING MKL libraries not found. MKL_LIBRARY_PATH=${MKL_LIBRARY_PATH} directory does not exist." )
    endif ( EXISTS ${MKL_LIBRARY_PATH} )
    
    set ( MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR} )
endif ( NOT EXISTS ${MKL_INCLUDE_DIR} )

include( FindPackageHandleStandardArgs )
# handle the QUIET and REQUIRED arguments and set MKL_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args ( MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIRS )

mark_as_advanced( MKL_INCLUDE_DIRS 
                  MKL_LIBRARIES
                  MKL_ROOT_PATH
                  MKL_LIBRARY_BLACS
                  MKL_LIBRARY_CORE 
                  MKL_LIBRARY_LP64
                  MKL_LIBRARY_GNU 
                  MKL_LIBRARY_INTEL )
                  
# if MKL is found, include directories    
if ( MKL_FOUND )
    include_directories ( ${MKL_INCLUDE_DIRS} )
endif ( MKL_FOUND )
