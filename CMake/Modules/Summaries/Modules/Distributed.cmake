###
 # @file CMake/Modules/Summaries/Modules/Distributed.cmake
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
 # @brief Summary concerning the distributed support (MPI and GPI).
 # @author Lauretta Schubert
 # @date 11.04.2016
###

set ( REQUIRED_FOUND FALSE )
if    ( MPI_ENABLED OR GPI_ENABLED )
  set ( REQUIRED_FOUND TRUE )
endif ( MPI_ENABLED OR GPI_ENABLED )

heading3 ( "Distributed" "REQUIRED_FOUND" )
    found_message ( "MPI" "MPI_FOUND" "OPTIONAL" "Version ${MPI_VERSION} at ${SCAI_MPI_INCLUDE_DIR}" )
    found_message ( "GPI" "GPI_FOUND" "OPTIONAL" "with:" )
    if    ( GPI_FOUND )
    	message ( STATUS "                                 GPI2 Version ${GPI2_VERSION} at ${GPI2_INCLUDE_DIR}" )
    	message ( STATUS "                                 IBVERBS at ${IBVERBS_INCLUDE_DIR}" )
    # no IBVERBS_VERSION
    #foreach    ( _LIB GPI2 IBVERBS )
    #    message ( STATUS "                                 ${_LIB} Version ${${_LIB}_VERSION} at ${${_LIB}_INCLUDE_DIR}" )
    #endforeach ( _LIB GPI2 IBVERBS )
    endif ( GPI_FOUND )