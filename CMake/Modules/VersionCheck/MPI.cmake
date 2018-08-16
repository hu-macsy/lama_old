###
 # @file MPI.cmake
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
 # @author Lauretta Schubert
 # @date 07.04.2016
###

if    ( MPI_FOUND )
    execute_process ( COMMAND ${MPIEXEC} --version OUTPUT_VARIABLE _mpi_output ERROR_VARIABLE _mpi_error)
    set ( _output "${_mpi_output}${_mpi_error}" ) # some version write output to error stream
    if    ( CMAKE_COMPILER_IS_GNUCXX )
    	string ( REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" MPI_VERSION ${_output} )
    	if    ( "${MPI_VERSION}" STREQUAL "" )
	    	string ( REGEX MATCH "([0-9]+\\.[0-9])" MPI_VERSION ${_output} )
	    	if    ( "${MPI_VERSION}" STREQUAL "" )
	    		string ( REGEX MATCH "([0-9])" MPI_VERSION ${_output} )
    		endif ( "${MPI_VERSION}" STREQUAL "" )
    	endif ( "${MPI_VERSION}" STREQUAL "" )
    endif ( CMAKE_COMPILER_IS_GNUCXX )
endif ( MPI_FOUND )
