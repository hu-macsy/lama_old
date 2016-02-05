###
 # @file scai/solver/CMakeLists.txt
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
 # @author Lauretta Schubert
 # @date 25.01.2016
 # @since 2.0.0
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Set PROJECT NAME / supported LANGUAGES
project ( SOLVER CXX )

### Set version of solver
# TOOD: where is the best place to define this?
set ( SCAI_SOLVER_VERSION_NAME "Alpamayo" )
set ( SCAI_SOLVER_COPYRIGHT_YEAR 2015 )
set ( SCAI_SOLVER_VERSION_MAJOR 1 )
set ( SCAI_SOLVER_VERSION_MINOR 0 )
set ( SCAI_SOLVER_VERSION_PATCH 1 )
set ( SCAI_SOLVER_VERSION_COUNT 3 )
set ( SCAI_SOLVER_VERSION ${SCAI_SOLVER_VERSION_MAJOR}.${SCAI_SOLVER_VERSION_MINOR}.${SCAI_SOLVER_VERSION_PATCH} )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

include ( VersionDefinition )

### Include Compiler Flags
if    ( NOT SCAI_COMPLETE_BUILD )
	include ( SetBuildFlags )
	include ( SCAIAssert )
	include ( SetCPPFlags )
	include ( Package/Sphinx )
else  ( NOT SCAI_COMPLETE_BUILD )
	include ( Functions/CheckValueAtCompleteBuild )
	CheckValueAtCompleteBuild ( solver )
endif ( NOT SCAI_COMPLETE_BUILD )

include ( CheckC++11 )
if    ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )
    include( Package/Boost )
endif ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )

### include all search packages
include ( InternalDependencies )
include ( ExternalDependencies )

foreach    ( module ${SCAI_SOLVER_EXTERNAL_DEPS} )
    include( Package/${module} )
endforeach ( module ${SCAI_SOLVER_EXTERNAL_DEPS} )

foreach    ( PACKAGE_TO_FIND ${SCAI_SOLVER_INTERNAL_DEPS} )
    find_package ( ${PACKAGE_TO_FIND} ${SCAI_FIND_PACKAGE_FLAGS} REQUIRED )
endforeach ( PACKAGE_TO_FIND $${SCAI_SOLVER_INTERNAL_DEPS} )

## finalize build flags
include ( ConcludeFlags )

## MKL FOUND
if     ( SCAI_BLAS_NAME STREQUAL "MKL" )
    add_definitions( -DSCAI_MKL_BLAS )
elseif ( INTERNALBLAS_FOUND )
    add_definitions( -DSCAI_INTERNAL_BLAS )
else   ( )
    add_definitions( -DADD_ )
endif  ( )

### add definitions ###

add_definitions( ${ADDITIONAL_WARNING_FLAGS} )

add_definitions ( -D${SCAI_LOGGING_FLAG} )
add_definitions ( -D${SCAI_TRACING_FLAG} )

add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )

if ( SCAI_OMP_SCHEDULE_FLAG )
    add_definitions ( -D${SCAI_OMP_SCHEDULE_FLAG} )
endif ( SCAI_OMP_SCHEDULE_FLAG )

if ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### Project System-Include
include_directories ( ../.. )

foreach    ( module ${SCAI_SOLVER_INTERNAL_DEPS} ${SCAI_SOLVER_EXTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    include_directories( ${${upper_module}_INCLUDE_DIR} )
endforeach ( module ${SCAI_SOLVER_INTERNAL_DEPS} ${SCAI_SOLVER_EXTERNAL_DEPS} )

include ( SCAIAssert )

### build solver library

include ( Functions/scaiProject )

### Add classes and header files

scai_project(

    CLASSES

        AMGSetup
        BiCG
        BiCGstab
        CG
    #   CGNE
        CGNR
        CGS
        DefaultJacobi
        GMRES
        InverseSolver
        IterativeSolver
        MINRES
        OmegaSolver
        Richardson
        SOR
        SimpleAMG
        SingleGridSetup
        SolutionProxy
        Solver
        SpecializedJacobi
        TFQMR
        TrivialPreconditioner

    HEADERS

        shared_lib
    )

### add subdirectories ###
add_subdirectory ( criteria )
add_subdirectory ( logger )

### for external amg ###
file ( GLOB AMG_SUBDIRPRESENT "amg" )

if ( AMG_SUBDIRPRESENT )
    add_subdirectory ( ${AMG_SUBDIRPRESENT} )
    message ( STATUS "amg subdir provides following solver libraries: ${SOLVER_SUB_LIBRARIES}" )
    #set ( SOLVER_SUB_LIBRARIES ${SOLVER_SUB_LIBRARIES} PARENT_SCOPE )
endif ( AMG_SUBDIRPRESENT )

## build doc
add_subdirectory( doc )

## build test
if ( FOUND_BOOST_TEST AND BUILD_TEST )
    add_subdirectory ( test )
endif ( FOUND_BOOST_TEST AND BUILD_TEST )

### add library ###

add_library ( scai_solver ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )

# Might be some additional solver libraries have been defined

message ( STATUS "SOLVER uses following solver libraries: ${SOLVER_SUB_LIBRARIES}" )
if    ( SOLVER_SUB_LIBRARIES )
    message ( STATUS "linked to solver: ${SOLVER_SUB_LIBRARIES}" )
    target_link_libraries( scai_solver ${SOLVER_SUB_LIBRARIES} )
endif ( SOLVER_SUB_LIBRARIES )

## add link libraries
list ( REVERSE SCAI_SOLVER_INTERNAL_DEPS )
foreach    ( module ${SCAI_SOLVER_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
	set ( SOLVER_LINK_LIBRARIES ${SOLVER_LINK_LIBRARIES} ${${upper_module}_LIBRARY} )
endforeach ( module ${SCAI_SOLVER_INTERNAL_DEPS} )
target_link_libraries ( scai_solver ${SCAI_START_LINK_LIBRARIES} ${SOLVER_LINK_LIBRARIES} ${SCAI_END_LINK_LIBRARIES} )

foreach    ( module ${SCAI_SOLVER_EXTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    target_link_libraries ( scai_solver ${SCAI_${upper_module}_LIBRARIES} )
endforeach ( module ${SCAI_SOLVER_EXTERNAL_DEPS} )

## install all Include files

include( Functions/relativeInstall )
relative_install( FILES ${CXX_HEADERS} DESTINATION include/scai/solver )

## install library.hpp
install ( FILES ../solver.hpp DESTINATION include/scai )

## build examples
add_subdirectory( examples )

# generate file make.inc + install it
#configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/examples/make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/examples/make.inc" )
#install ( FILES "${CMAKE_CURRENT_BINARY_DIR}/examples/make.inc" DESTINATION share/examples )

# copy examples in build directory + install it
#file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
#install ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION share )

## add custum target distClean
include ( DistClean )

### install
install ( TARGETS scai_solver DESTINATION lib )

if    ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
	include ( Summaries/solver )
endif ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
