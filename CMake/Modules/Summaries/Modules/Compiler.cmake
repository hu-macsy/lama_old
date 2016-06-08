###
 # @file CMake/Modules/Summaries/Modules/Compiler.cmake
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
 # @endlicense
 #
 # @brief Summary concerning the compiler support.
 # @author Lauretta Schubert
 # @date 11.04.2016
###

heading ( "Compiler:" )

if    ( CMAKE_CXX_COMPILER AND ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR ) )
    set( REQUIRED_FOUND TRUE )
else  ( CMAKE_CXX_COMPILER AND ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR ) )
    set( REQUIRED_FOUND FALSE )
    message ( FATAL_ERROR "Compiler Configuration incomplete" )
endif ( CMAKE_CXX_COMPILER AND ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR ) )

heading2 ( "Configuration" "REQUIRED_FOUND" )
    found_message ( "C++ Compiler" "CMAKE_CXX_COMPILER" "REQUIRED" "${CMAKE_CXX_COMPILER_ID} Version ${CXX_COMPILER_VERSION}" )
    found_message ( "with C++11 support" "CXX_SUPPORTS_C11" "REQUIRED" "" )

if    ( NOT CXX_SUPPORTS_C11 )
    emptyline()
    message ( STATUS "Either compiler supporting C++11 or Boost needed." )
    found_message ( "Boost" "SCAI_BOOST_INCLUDE_DIR" "REQUIRED" "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )
endif ( NOT CXX_SUPPORTS_C11 )