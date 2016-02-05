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
project ( scai_solver CXX )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

### includes ###

## important Variables for project
include ( SetProjectVariables)
include ( VersionDefinition )
include ( InternalDependencies )
include ( ExternalDependencies )

## needed functions
include ( Functions/scaiFunctions )
include ( Functions/scaiProject )
include ( Functions/relativeInstall )

### Include compile dependent things
if    ( NOT SCAI_COMPLETE_BUILD )
    include ( SetBuildFlags )
    include ( CompilerVersion )
    include ( SCAIAssert )
    include ( SetCPPFlags )
    include ( Package/Sphinx )
else  ( NOT SCAI_COMPLETE_BUILD )
    include ( Functions/CheckValueAtCompleteBuild )
    CheckValueAtCompleteBuild ( ${PROJECT_SURNAME} )
endif ( NOT SCAI_COMPLETE_BUILD )

## C++11 or Boost
include ( CheckC++11 )
if    ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )
    include( Package/Boost )
endif ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )

## finalize build flags
include ( ConcludeFlags )

### add definitions ###

add_definitions( ${ADDITIONAL_WARNING_FLAGS} )
add_definitions ( -D${SCAI_LOGGING_FLAG} )
add_definitions ( -D${SCAI_TRACING_FLAG} )
add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )

## MKL FOUND
if     ( SCAI_BLAS_NAME STREQUAL "MKL" )
    add_definitions( -DSCAI_MKL_BLAS )
elseif ( INTERNALBLAS_FOUND )
    add_definitions( -DSCAI_INTERNAL_BLAS )
else   ( )
    add_definitions( -DADD_ )
endif  ( )

if ( SCAI_OMP_SCHEDULE_FLAG )
    add_definitions ( -D${SCAI_OMP_SCHEDULE_FLAG} )
endif ( SCAI_OMP_SCHEDULE_FLAG )

if ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### Project System-Include
include_directories ( ../.. )
addInternalAndExternalIncludes()

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp

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
        QMR
        Richardson
        SOR
        SimpleAMG
        SingleGridSetup
        SolutionProxy
        Solver
        SpecializedJacobi
        TFQMR
        TrivialPreconditioner

    HEADERS                  # .hpp only

        shared_lib
    )

## add subdirectories
add_subdirectory ( criteria )
add_subdirectory ( doc )
add_subdirectory ( examples )
add_subdirectory ( logger )

if ( FOUND_BOOST_TEST AND BUILD_TEST )
    add_subdirectory ( test )
endif ( FOUND_BOOST_TEST AND BUILD_TEST )

### for external amg ###
file ( GLOB AMG_SUBDIRPRESENT "amg" )

if ( AMG_SUBDIRPRESENT )
    add_subdirectory ( ${AMG_SUBDIRPRESENT} )
    message ( STATUS "amg subdir provides following solver libraries: ${SOLVER_SUB_LIBRARIES}" )
    #set ( SOLVER_SUB_LIBRARIES ${SOLVER_SUB_LIBRARIES} PARENT_SCOPE )
endif ( AMG_SUBDIRPRESENT )

### add library ###
add_library ( scai_solver ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )
set_target_properties ( ${PROJECT_NAME} PROPERTIES VERSION ${${UPPER_PROJECT_NAME}_VERSION} )

## add link libraries
addInternalAndExternalLinkLibraries()
# Might be some additional solver libraries have been defined
#message ( STATUS "SOLVER uses following solver libraries: ${SOLVER_SUB_LIBRARIES}" )
if    ( SOLVER_SUB_LIBRARIES )
#    message ( STATUS "linked to solver: ${SOLVER_SUB_LIBRARIES}" )
    target_link_libraries( scai_solver ${SOLVER_SUB_LIBRARIES} )
endif ( SOLVER_SUB_LIBRARIES )

### install ###

## install lib
install ( TARGETS scai_solver DESTINATION lib )

## install headers
relative_install( FILES ${CXX_HEADERS} DESTINATION include/scai/solver )

## install library.hpp
install ( FILES ../solver.hpp DESTINATION include/scai )

# generate file make.inc + install it
#configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/examples/make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/examples/make.inc" )
#install ( FILES "${CMAKE_CURRENT_BINARY_DIR}/examples/make.inc" DESTINATION share/examples )

# copy examples in build directory + install it
#file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
#install ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION share )

## add custum target distClean
include ( DistClean )

### summary ###
if    ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
	include ( Summaries/solver )
endif ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
