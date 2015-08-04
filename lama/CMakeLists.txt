###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2015
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
 # @since 1.0.0
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

set ( LAMA_VERSION_NAME "Alpamayo" )
set ( LAMA_COPYRIGHT_YEAR 2015 )
set ( LAMA_VERSION_MAJOR 1 )
set ( LAMA_VERSION_MINOR 0 )
set ( LAMA_VERSION_PATCH 1 )
set ( LAMA_VERSION_COUNT 3 )
set ( LAMA_VERSION ${LAMA_VERSION_MAJOR}.${LAMA_VERSION_MINOR}.${LAMA_VERSION_PATCH} )

### Variables / CMAKE Functions ###
include ( Variables ) #LAMA.CMakeLists
include ( SearchPackages/LAMA )
include ( Functions/allFunctions )
include ( LAMAdefines ) #LAMA.CMakeLists
include ( CompilerVersion ) #LAMA.CMakeLists

### PACKAGES ###
include ( Package/OpenMP )
include ( Package/CUDA )
include ( Package/MIC )
include ( Package/Boost ) #LAMA.CMakeLists
include ( Package/MPI )
include ( Package/CodeCoverage ) #LAMA.CMakeLists

# Find required packages
foreach    ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} ${LAMA_FIND_PACKAGE_FLAGS} REQUIRED )
endforeach ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )

# Find optional packages
foreach    ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} ${LAMA_FIND_PACKAGE_FLAGS} )
endforeach ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )

# LAMA irrelevant entries will be marked as advanced ( Remove them from default cmake GUI )
mark_as_advanced ( CUDA_BUILD_CUBIN CUDA_BUILD_EMULATION CUDA_SDK_ROOT_DIR CUDA_VERBOSE_BUILD )

### Include Compiler Flags, after all packages are found ###
# TB: removed, now already done by other packages: include ( CompilerFlags )
# TB: only flags concerning NVCC are used ( not required by other projects)
include ( SetCPPFlags ) #LAMA.CMakeLists

## MKL FOUND
if     ( LAMA_BLAS_NAME STREQUAL "MKL" )
    add_definitions( -DLAMA_MKL_BLAS )
elseif ( INTERNALBLAS_FOUND )
    add_definitions( -DLAMA_INTERNAL_BLAS )
else   ( )
    add_definitions( -DADD_ )
endif  ( )

## Project System-Include
include_directories ( . )

if    ( Boost_INCLUDE_DIR )
    include_directories ( ${Boost_INCLUDE_DIR} )
endif ( Boost_INCLUDE_DIR )

### add subdirectories ###
add_subdirectory ( lama )

if ( FOUND_BOOST_TEST AND LAMA_BUILD_TEST )
    add_subdirectory ( test )
endif ( FOUND_BOOST_TEST AND LAMA_BUILD_TEST )

### Examples ###
# copy examples in build directory + install it
file ( COPY ${CMAKE_SOURCE_DIR}/examples DESTINATION ${CMAKE_BINARY_DIR} )
install ( DIRECTORY ${CMAKE_SOURCE_DIR}/examples DESTINATION share )

### Doxygen and Sphinx Configuration ###
include ( Package/Doxygen )
include ( Package/Sphinx )

## install library.hpp
configure_file ( "${CMAKE_SOURCE_DIR}/lama/lama.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" )
install ( FILES "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" DESTINATION include )

## add custum target distClean
include ( DistClean )

include ( LAMASummary )
