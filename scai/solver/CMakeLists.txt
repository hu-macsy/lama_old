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

scai_module (

    MODULE_NAME   solver 

    INTERNAL_DEPS common logging tracing tasking hmemo kregistry 
                  blaskernel utilskernel sparsekernel dmemo lama 
) 

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

scai_subdirectories ( criteria examples logger test )

### for external amg ###

file ( GLOB AMG_SUBDIRPRESENT "amg" )

if ( AMG_SUBDIRPRESENT )
    add_subdirectory ( ${AMG_SUBDIRPRESENT} )
    message ( STATUS "amg subdir provides following solver libraries: ${SOLVER_SUB_LIBRARIES}" )
    #set ( SOLVER_SUB_LIBRARIES ${SOLVER_SUB_LIBRARIES} PARENT_SCOPE )
endif ( AMG_SUBDIRPRESENT )

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

