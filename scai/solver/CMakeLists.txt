###
 # @file scai/solver/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2017
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
 # @brief Main specification file for CMake
 # @author Lauretta Schubert
 # @date 25.01.2016
###

## Define name of the SCAI module project

set ( MODULE_NAME solver )

set ( INTERNAL_DEPS common logging tracing tasking hmemo kregistry blaskernel utilskernel sparsekernel dmemo lama ) 

## Define all library names with the (global) prefix SCAI_LIBRARY_PREFIX

set ( PROJECT_NAME "${SCAI_LIBRARY_PREFIX}${MODULE_NAME}" )

## finalize build flags
include ( ConcludeFlags )

### add definitions ###

include ( Settings/logging )
include ( Settings/tracing )

add_definitions ( ${ADDITIONAL_WARNING_FLAGS} )
add_definitions ( -DSCAI_LOGGING_LEVEL_${SCAI_LOGGING_LEVEL} )
add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )
add_definitions ( -DSCAI_TRACING_${SCAI_TRACING} )

if ( SCAI_OMP_SCHEDULE_FLAG )
    add_definitions ( -D${SCAI_OMP_SCHEDULE_FLAG} )
endif ( SCAI_OMP_SCHEDULE_FLAG )

if ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### Project System-Include ( also configured files )

set( CONFIGURE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include )
include_directories ( ../.. ${CONFIGURE_INCLUDE_DIR} )

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp

        AMGSetup
        BiCG
        BiCGstab
        CG
        CGNE
        CGNR
        CGS
        DecompositionSolver
        GMRES
        InverseSolver
        IterativeSolver
        Jacobi
        Kaczmarz
        MINRES
        OmegaSolver
        QMR
        Richardson
        SimpleAMG
        SingleGridSetup
        SolutionProxy
        Solver
        TFQMR
        TrivialPreconditioner

    )

## add subdirectories
add_subdirectory ( criteria )

if    ( BUILD_EXAMPLES )
add_subdirectory ( examples )
endif ( BUILD_EXAMPLES )

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

### set link directories ###
if ( WIN32 )
    link_directories ( ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

### add library ###

## Define library via macro using MODULE_NAME, INTERNAL_DEPS, EXTERNAL_DEPS
## and global variables SCAI_LIBRARY_PREFIX, SCAI_VERSION, ....

# Might be some additional solver libraries have been defined

scai_library ( OTHER_LINK ${SOLVER_SUB_LIBRARIES} )

### install ###

## install headers

relative_install ( FILES ${CXX_HEADERS} DESTINATION "include/scai/${MODULE_NAME}" )

## install library.hpp
install ( FILES ../solver.hpp DESTINATION include/scai )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 

include ( CustomCommands/SphinxDoc )

