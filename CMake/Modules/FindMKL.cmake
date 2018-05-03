###
 # @file FindMKL.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief Find MKL
 # @author 
 # @date 25.04.2013
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
if   ( MKL_ROOT AND NOT MKL_INCLUDE_DIR )
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
                    message ( "Found multiple MKL installations. Will choose: " ${MKL_SEARCH_DIR} )
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
        
        if ( SCAI_CMAKE_VERBOSE )
            message ( STATUS "MKL_ROOT_PATH = ${MKL_ROOT_PATH}" )
        endif ( SCAI_CMAKE_VERBOSE )
        
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

    if ( SCAI_CMAKE_VERBOSE )
        message( STATUS "MKL_LIBRARY_PATH = ${MKL_LIBRARY_PATH}, MKL_LIBRARY_PATH_SUFFIXES = ${MKL_LIBRARY_PATH_SUFFIXES}" )
    endif ( SCAI_CMAKE_VERBOSE )
    
    ### Search for libraries mkl_gnu_thread, mkl_intel_thread and mkl_core
    
    if ( EXISTS ${MKL_LIBRARY_PATH} )

        # search for mkl_intel_lp64 lib 

        find_library ( MKL_LIBRARY_LP64 mkl_intel_lp64 PATHS ${MKL_LIBRARY_PATH} PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )

        if ( NOT EXISTS ${MKL_LIBRARY_LP64} )
            message ( STATUS "WARNING MKL library mkl_intel_lp64 not found with MKL_LIBRARY_PATH=${MKL_LIBRARY_PATH}." )
        endif( NOT EXISTS ${MKL_LIBRARY_LP64} )

        set ( MKL_LIBRARIES  ${MKL_LIBRARY_LP64} )
        
        # search for mkl_core lib 

        find_library ( MKL_LIBRARY_CORE mkl_core PATHS ${MKL_LIBRARY_PATH} PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )

        if ( NOT EXISTS ${MKL_LIBRARY_CORE} )
            message ( STATUS "WARNING MKL library mkl_core not found with MKL_LIBRARY_PATH=${MKL_LIBRARY_PATH}." )
        endif ( NOT EXISTS ${MKL_LIBRARY_CORE} )
        
        # search for mkl_thread lib
            
        # search for gnu compiler libs
        if ( CMAKE_CXX_COMPILER_ID MATCHES GNU )
            if ( OPENMP_FOUND AND USE_OPENMP )
                find_library ( MKL_LIBRARY_GNU mkl_gnu_thread PATHS ${MKL_LIBRARY_PATH} PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )
            else ( OPENMP_FOUND AND USE_OPENMP )
                find_library ( MKL_LIBRARY_GNU mkl_sequential PATHS ${MKL_LIBRARY_PATH} PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )
            endif ( OPENMP_FOUND AND USE_OPENMP )
            if ( NOT EXISTS ${MKL_LIBRARY_GNU} )
                message ( STATUS "WARNING MKL library mkl_gnu_thread not found with MKL_LIBRARY_PATH=${MKL_LIBRARY_PATH}." )
            endif( NOT EXISTS ${MKL_LIBRARY_GNU} )
            list ( APPEND MKL_LIBRARIES ${MKL_LIBRARY_GNU} )
         endif( CMAKE_CXX_COMPILER_ID MATCHES GNU )
        
         # search for intel compiler libs ( for icc and clang (llvm) )
         if ( CMAKE_CXX_COMPILER_ID MATCHES Intel OR CMAKE_CXX_COMPILER_ID MATCHES Clang )
            if ( OPENMP_FOUND AND USE_OPENMP )
                find_library ( MKL_LIBRARY_INTEL mkl_intel_thread PATHS ${MKL_LIBRARY_PATH} PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )
            else ( OPENMP_FOUND AND USE_OPENMP )
                find_library ( MKL_LIBRARY_INTEL mkl_sequential PATHS ${MKL_LIBRARY_PATH} PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )
            endif ( OPENMP_FOUND AND USE_OPENMP )
            if ( NOT EXISTS ${MKL_LIBRARY_INTEL} )
                message ( STATUS "WARNING MKL library mkl_intel_thread not found with MKL_LIBRARY_PATH=${MKL_LIBRARY_PATH}." )
            endif()
            list ( APPEND MKL_LIBRARIES ${MKL_LIBRARY_INTEL} )
         endif( CMAKE_CXX_COMPILER_ID MATCHES Intel OR CMAKE_CXX_COMPILER_ID MATCHES Clang )
        
         if ( ${CMAKE_GENERATOR} MATCHES "Visual Studio" )
            find_library ( MKL_LIBRARY_INTEL mkl_intel_thread PATHS ${MKL_LIBRARY_PATH} PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )
            if ( NOT EXISTS ${MKL_LIBRARY_INTEL} )
                message ( STATUS "WARNING MKL library mkl_intel_thread not found with MKL_LIBRARY_PATH=${MKL_LIBRARY_PATH}." )
            endif ( NOT EXISTS ${MKL_LIBRARY_INTEL} )
            list ( APPEND MKL_LIBRARIES ${MKL_LIBRARY_INTEL} )
         endif( ${CMAKE_GENERATOR} MATCHES "Visual Studio" )

        # conclude libs
        list ( APPEND MKL_LIBRARIES ${MKL_LIBRARY_CORE} )
        
        if ( WIN32 )
            find_library ( INTEL_IOMP_LIBRARY libiomp5md PATHS ${MKL_LIBRARY_PATH} ${MKL_ROOT_PATH}/../compiler/lib PATH_SUFFIXES ${MKL_LIBRARY_PATH_SUFFIXES} )
            if ( NOT EXISTS ${INTEL_IOMP_LIBRARY} )
                message ( STATUS "WARNING Intel OMP library iomp5md not found with INTEL_COMPILER_PATH=${MKL_ROOT_PATH}/../compiler/lib." )
            endif ( NOT EXISTS ${INTEL_IOMP_LIBRARY} )
            list ( APPEND MKL_LIBRARIES ${INTEL_IOMP_LIBRARY} )
        endif ( WIN32 )
        
        if ( SCAI_CMAKE_VERBOSE )
            message ( STATUS "Found MKL Libraries: ${MKL_LIBRARIES}." )
            message ( STATUS "Found MKL PLibraries: ${MKL_PLIBRARIES}." )
        endif ( SCAI_CMAKE_VERBOSE )
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
