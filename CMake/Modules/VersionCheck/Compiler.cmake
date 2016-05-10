###
 # @file Compiler.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
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
 # @endlicense
 #
 # @brief Version variable defintions for the used compilers
 # @author Jan Ecker
 # @date 25.04.2013
###

### GNU compiler

## C Compiler
if    ( CMAKE_COMPILER_IS_GNUCC )
    execute_process ( COMMAND ${CMAKE_C_COMPILER} --version OUTPUT_VARIABLE _compiler_output )
    string ( REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" GNUCC_COMPILER_VERSION ${_compiler_output} )
endif ( CMAKE_COMPILER_IS_GNUCC )

## CXX Compiler
if ( CMAKE_COMPILER_IS_GNUCXX )
    execute_process ( COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE _compiler_output )
    string ( REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" GNUCXX_COMPILER_VERSION ${_compiler_output} )
endif ( CMAKE_COMPILER_IS_GNUCXX )

### Intel compiler

## C Compiler
if    ( CMAKE_CC_COMPILER_ID MATCHES Intel )
    execute_process ( COMMAND ${CMAKE_C_COMPILER} --version OUTPUT_VARIABLE _compiler_output )
    string ( REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" IntelCC_COMPILER_VERSION ${_compiler_output} )
endif ( CMAKE_CC_COMPILER_ID MATCHES Intel )

## CXX Compiler
if    ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
    execute_process ( COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE _compiler_output )
    string ( REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" IntelCXX_COMPILER_VERSION ${_compiler_output} )
endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
