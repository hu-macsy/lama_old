###
 # @file Compiler.cmake
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
 # @brief Version variable defintions for the used compilers
 # @author Jan Ecker
 # @date 25.04.2013
###

## C Compiler
#needs not to be enabled, because the language C is not enabled --> not CMAKE_C_COMPILER
#execute_process ( COMMAND ${CMAKE_C_COMPILER} --version OUTPUT_VARIABLE _cc_compiler_output )
#string ( REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" CC_COMPILER_VERSION ${_cc_compiler_output} )

## CXX Compiler

if ( DEFINED CMAKE_CXX_COMPILER_VERSION )
    # later CMAKE versions deliver the version of the compiler in a variable
    set ( CXX_COMPILER_VERSION "${CMAKE_CXX_COMPILER_VERSION}" )
else ( DEFINED CMAKE_CXX_COMPILER_VERSION )
    # Query the version of the compiler by a Linux command
    execute_process ( COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE _cxx_compiler_output ERROR_VARIABLE _cxx_compiler_output )
    string ( REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" CXX_COMPILER_VERSION ${_cxx_compiler_output} )
    if    ( "${CXX_COMPILER_VERSION}" STREQUAL "" )
        string ( REGEX MATCH "([0-9]+\\.[0-9])" CXX_COMPILER_VERSION ${_cxx_compiler_output} )
        if    ( "${CXX_COMPILER_VERSION}" STREQUAL "" )
    	    string ( REGEX MATCH "([0-9])" CXX_COMPILER_VERSION ${_cxx_compiler_output} )
	    endif ( "${CXX_COMPILER_VERSION}" STREQUAL "" ) 
    endif ( "${CXX_COMPILER_VERSION}" STREQUAL "" ) 
endif ( DEFINED CMAKE_CXX_COMPILER_VERSION )
