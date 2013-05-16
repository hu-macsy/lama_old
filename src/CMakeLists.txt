###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2011
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
 # @brief Main specification file for CMake
 # @author Thomas Brandes, Jiri Kraus
 # @date 29.03.2012
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Print LAMA Installation header
message ( STATUS "" )
message ( STATUS "LAMA Configuration:" )
message ( STATUS "===================" )
message ( STATUS "" )

### Set PROJECT NAME / supported LANGUAGES and Module Path ###
project ( LAMA CXX C )
set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../CMake/Modules" )
set ( LAMA_VERSION "1.0.0" )
set ( LAMA_VERSION_MAJOR 1 )
set ( LAMA_VERSION_MINOR 0 )
set ( LAMA_VERSION_PATCH 0 )
set ( LAMA_VERSION_COUNT 3 )



### Variables / CMAKE Functions ###
include ( Variables )
include ( SearchPackages )
include ( Functions )
include ( LAMAdefines )
include ( CompilerVersion )

### PACKAGES ###
include ( PackageMPI )
include ( PackageBoost )
include ( PackageCodeCoverage )

# Find required packages
foreach ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} ${LAMA_FIND_PACKAGE_FLAGS} REQUIRED )
endforeach ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )

# Find optional packages
foreach ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} ${LAMA_FIND_PACKAGE_FLAGS} )
endforeach ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )

# ALLOW to switch off CUDA explicitly
setAndCheckCache ( "CUDA" )

### Include Compiler Flags, after all packages are found ###
include ( CompilerFlags )

## OPENMP
if ( NOT LAMA_OMP_SCHEDULE )
    set ( LAMA_OMP_SCHEDULE "static" )
endif ( NOT LAMA_OMP_SCHEDULE )

# TODO: SUMMARY
#message ( STATUS "OpenMP schedule set to " ${LAMA_OMP_SCHEDULE} )
add_definitions ( -DLAMA_OMP_SCHEDULE=${LAMA_OMP_SCHEDULE} )

if ( OPENMP_FOUND )
    set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif( OPENMP_FOUND )



# LAMA irrelevant entries will be marked as advanced ( Remove them from default cmake GUI )
mark_as_advanced ( CUDA_BUILD_CUBIN CUDA_BUILD_EMULATION CUDA_SDK_ROOT_DIR CUDA_VERBOSE_BUILD )


## Project System-Include
include_directories ( . )

if ( Boost_INCLUDE_DIR )
    include_directories ( ${Boost_INCLUDE_DIR} )
endif ( Boost_INCLUDE_DIR )


### add subdirectories ###
add_subdirectory ( lama )
add_subdirectory ( logging )
if ( FOUND_BOOST_TEST AND LAMA_BUILD_TEST )
    add_subdirectory ( test )
endif ( FOUND_BOOST_TEST AND LAMA_BUILD_TEST )

include ( PackageDoxygen )
include ( PackageSphinx )


### install ###
configure_file( "${CMAKE_SOURCE_DIR}/lama/lama.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" )
# TODO: SUMMARY
#message ( STATUS "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" )
install ( FILES "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" DESTINATION include )


### DISTCLEAN TARGET ###
add_custom_target ( distclean )
file ( GLOB_RECURSE BUILD_GLOB_RES ${CMAKE_BINARY_DIR}/* )
add_custom_command (
    TARGET distclean
    DEPENDS clean
    COMMAND ${CMAKE_COMMAND} -E remove ${BUILD_GLOB_RES}
)

### Summary ###
message ( STATUS "" )
#message ( STATUS "        _    _")
#message ( STATUS "       ( \\__//)")
#message ( STATUS "       .'     )")
#message ( STATUS "    __/b d  .  )")
#message ( STATUS "   (_Y_`,     .)")
#message ( STATUS "    `--'-,-'  )")
#message ( STATUS "         (.  )")
#message ( STATUS "         (   )")
#message ( STATUS "        (   )")
#message ( STATUS "       ( . )         .---.")
#message ( STATUS "      (    )        (     )")
#message ( STATUS "      (   . )      (  .    )")
#message ( STATUS "      (      )    (      .  ),")
#message ( STATUS "      ( .     `\"'`  .       `)\\")
#message ( STATUS "       (      .              .)\\")
#message ( STATUS "       ((  .      .   (   .   )\\\\")
#message ( STATUS "       ((       .    (        ) \\\\")
#message ( STATUS "        ((     )     _( .   . )  \\\\")
#message ( STATUS "        ( ( .   )\"'\"`(.(     )   ( ;")
#message ( STATUS "        ( (    )      ( ( . )     \\'")
#message ( STATUS "         |~(  )        |~(  )")
#message ( STATUS "         | ||~|        | ||~|")
#message ( STATUS "    jgs  | || |        | || |  ")  
#message ( STATUS "        _| || |       _| || |")
#message ( STATUS "       /___(| |      /___(| |")
#message ( STATUS "          /___(         /___(")
message ( STATUS "" )
message ( STATUS "Summary of LAMA Configuration:" )
message ( STATUS "==============================" )
message ( STATUS "" )

lama_status_message ( HEADLINE "REQUIRED PACKAGES:" )
## C++ Compiler
if ( CMAKE_CXX_COMPILER )
    lama_status_message ( INFO "[FOUND]" "    C++ Compiler - ${CMAKE_CXX_COMPILER_ID}"
                               "${${CMAKE_CXX_COMPILER_ID}CXX_COMPILER_VERSION}" )
else ( CMAKE_CXX_COMPILER )
    lama_status_message ( ERROR "[NOT FOUND]" "C++ Compiler" )
endif ( CMAKE_CXX_COMPILER )

# C Compiler
if ( CMAKE_C_COMPILER )
    lama_status_message ( INFO "[FOUND]" "    C Compiler - ${CMAKE_C_COMPILER_ID}" 
                               "${${CMAKE_CXX_COMPILER_ID}CC_COMPILER_VERSION}" )
else ( CMAKE_C_COMPILER )
    lama_status_message ( ERROR "[NOT FOUND]" "C Compiler" )
endif ( CMAKE_C_COMPILER )

# Boost 
if ( Boost_INCLUDE_DIR )
    lama_status_message ( INFO "[FOUND]" "    Boost - Version"
                               "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}"
                               "at ${Boost_INCLUDE_DIR}" )
else ( Boost_INCLUDE_DIR )
    lama_status_message ( ERROR "[NOT FOUND]" "Boost" )
endif ( Boost_INCLUDE_DIR )

# BLAS
if ( LAMA_BLAS_LIBRARIES )
    lama_status_message ( INFO "[FOUND]" "    BLAS (${LAMA_BLAS_NAME}) with libraries: ${LAMA_BLAS_LIBRARIES}" )
    if ( LAPACK_FOUND ) 
        lama_status_message ( INFO "[FOUND]" "    LAPACK" )
    else ( LAPACK_FOUND )
        lama_status_message ( ERROR "[NOT FOUND]" "LAPACK" )
    endif ( LAPACK_FOUND )
else ( LAMA_BLAS_LIBRARIES )
    lama_status_message ( ERROR "[NOT FOUND]" "BLAS" )
endif ( LAMA_BLAS_LIBRARIES )

## OpenMP
if( OPENMP_FOUND )
    lama_status_message ( INFO "[FOUND]" "    OpenMP" )
else ( OPENMP_FOUND )
    lama_status_message ( ERROR "[NOT FOUND]" "OpenMP" )
endif ( OPENMP_FOUND )



message ( STATUS "" )
lama_status_message ( HEADLINE "OPTIONAL PACKAGES:" )

## MPI
if ( MPI_FOUND )
    if ( LAMA_USE_MPI )
        lama_status_message ( INFO "[FOUND]" "    MPI at ${MPI_INCLUDE_PATH}" )
    else ( LAMA_USE_MPI )
        lama_status_message ( WARNING "[DISABLED]" " MPI at ${MPI_INCLUDE_PATH}" )
    endif ( LAMA_USE_MPI )
else ( MPI_FOUND )
    lama_status_message ( WARNING "[NOT FOUND]" "MPI not found" )
endif ( MPI_FOUND )

## CUDA
if ( CUDA_FOUND )
    if ( LAMA_USE_CUDA )
        lama_status_message ( INFO "[FOUND]" "    CUDA ${CUDA_VERSION} at ${CUDA_INCLUDE_DIRS}" )
    else ( LAMA_USE_CUDA )
        lama_status_message ( WARNING "[DISABLED]" " CUDA ${CUDA_VERSION} at ${CUDA_INCLUDE_DIRS}" )
    endif ( LAMA_USE_CUDA )
else ( CUDA_FOUND )
    lama_status_message ( WARNING "[NOT FOUND]" "CUDA" )
endif ( CUDA_FOUND )

## DOXYGEN
if ( DOXYGEN_FOUND )
    lama_status_message ( INFO "[FOUND]" "    DOXYGEN - 'make doc' to build documentation" )
else ( DOXYGEN_FOUND )
    lama_status_message ( WARNING "[NOT FOUND]" "DOXYGEN" )
endif ( DOXYGEN_FOUND )

## Boost Unit Test Framework
if ( Boost_UNIT_TEST_FRAMEWORK_FOUND )
    lama_status_message ( INFO "[FOUND]" "    Boost Unit Test Framework" )
else ( Boost_UNIT_TEST_FRAMEWORK_FOUND )
    lama_status_message ( WARNING "[NOT FOUND]" "Boost Unit Test Framework" )
endif ( Boost_UNIT_TEST_FRAMEWORK_FOUND )

## Boost regex
if ( Boost_REGEX_FOUND )
    lama_status_message ( INFO "[FOUND]" "    Boost Regex" )
else ( Boost_REGEX_FOUND )
    lama_status_message ( WARNING "[NOT FOUND]" "Boost Regex" )
endif ( Boost_REGEX_FOUND )



message ( STATUS "" )


lama_status_message ( HEADLINE "INFO:" )
message ( STATUS "LAMA Version : ${LAMA_VERSION}" )
message ( STATUS "Build Type   : ${CMAKE_BUILD_TYPE}" )
message ( STATUS "LOG Level    : ${LAMA_LOG_LEVEL}" )
message ( STATUS "ASSERT Level : ${LAMA_ASSERT_LEVEL}" )
message ( STATUS "TRACE Level  : ${LAMA_TRACE_LEVEL}" )
message ( STATUS "" )