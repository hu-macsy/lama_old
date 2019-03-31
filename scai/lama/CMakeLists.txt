###
 # @file scai/lama/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief CMake configuration for module project lama ( linear algebra math )
 # @author Thomas Brandes
 # @date 04.07.2017
###

## Include required scai macros

include ( scai_macro/scai_module )
include ( scai_macro/scai_project )
include ( scai_macro/scai_subdirectories )
include ( scai_macro/scai_library )
include ( scai_function/relative_install )

## Define new SCAI module project with its internal/external dependencies

scai_module(

    MODULE_NAME   lama

    INTERNAL_DEPS common logging tracing tasking hmemo kregistry 
                  blaskernel utilskernel sparsekernel dmemo 

    EXTERNAL_DEPS OpenMP ZLIB PNG
)

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp

        DenseVector
        GridVector
        GridSection
        SparseVector

        VectorKind
        _Vector
        Vector
        VectorAssembly

        GridVector
        GridSection

    HEADERS                  # .hpp only

        Scalar
        GridReadAccess
        GridWriteAccess
        fft
        freeFunction
    )

### add subdirectories

scai_subdirectories( expression io matrix matutils mepr norm storage TEST EXAMPLES )

## Define library via macro using MODULE_NAME, INTERNAL_DEPS, EXTERNAL_DEPS
## and global variables SCAI_LIBRARY_PREFIX, SCAI_VERSION, ....

scai_library ( PREFIX  ${SCAI_LIBRARY_PREFIX}
               TYPE    ${SCAI_LIBRARY_TYPE}
               VERSION ${SCAI_VERSION}
               ${CXX_SOURCES} )

### install ###

## install headers

relative_install ( FILES ${CXX_HEADERS} DESTINATION "include/scai/${MODULE_NAME}" )

## install library.hpp

configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/lama.hpp.in" "${CMAKE_BINARY_DIR}/include/scai/lama.hpp" )
install ( FILES "${CMAKE_BINARY_DIR}/include/scai/lama.hpp" DESTINATION include/scai )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 

include ( CustomCommands/SphinxDoc )

