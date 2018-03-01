###
 # @file Package/MPI.cmake
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
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief findPackage and configuration of MPI
 # @author Jan Ecker
 # @date 25.04.2013
###

include ( scai_macro/scai_pragma_once )
include ( scai_macro/scai_build_variable )
include ( scai_macro/scai_summary )

### MPI_FOUND            - if MPI is found
### USE_MPI              - if MPI is enabled
### SCAI_MPI_INCLUDE_DIR - MPI include directory
### SCAI_MPI_LIBRARIES   - all needed MPI libraries
### MPI_ENABLED          - if MPI_FOUND AND USE_MPI

scai_pragma_once()

# Look for MPI first to allow SCAI_BLAS to take the correct blacs implementation
# based on the found mpi
if    ( WIN32 AND NOT ( MPI_C_INCLUDE_PATH OR MPI_CXX_INCLUDE_PATH OR MPI_C_LIBRARIES OR MPI_CXX_LIBRARIES ) )
    if    ( MPI_ROOT )
        set ( MPI_C_INCLUDE_PATH "${MPI_ROOT}/Inc" )
        set ( MPI_CXX_INCLUDE_PATH "${MPI_ROOT}/Inc" )
        set ( LAMA_MPI_LIB_DIR "${MPI_ROOT}/Lib" )
    else  ( MPI_ROOT )
        set ( MPI_C_INCLUDE_PATH "C:/Program Files/Microsoft HPC Pack 2008 R2/Inc" )
        set ( MPI_CXX_INCLUDE_PATH "C:/Program Files/Microsoft HPC Pack 2008 R2/Inc" )
        set ( LAMA_MPI_LIB_DIR "C:/Program Files/Microsoft HPC Pack 2008 R2/Lib" )
    endif ( MPI_ROOT )
    
    if    ( CMAKE_CL_64 )
        set ( LAMA_MPI_LIB_DIR "${LAMA_MPI_LIB_DIR}/amd64" )
    else  ( CMAKE_CL_64 )
        set ( LAMA_MPI_LIB_DIR "${LAMA_MPI_LIB_DIR}/i386" )
    endif ( CMAKE_CL_64 )
    
    set ( MPI_C_LIBRARIES "${LAMA_MPI_LIB_DIR}/msmpi.lib" )
    set ( MPI_CXX_LIBRARIES "${LAMA_MPI_LIB_DIR}/msmpi.lib" )
else  ( WIN32 AND NOT ( MPI_C_INCLUDE_PATH OR MPI_CXX_INCLUDE_PATH OR MPI_C_LIBRARIES OR MPI_CXX_LIBRARIES ) )
  if    ( MPI_ROOT )
      set ( MPI_COMPILER ${MPI_ROOT}/bin/mpicxx )
  endif ( MPI_ROOT )
endif ( WIN32 AND NOT ( MPI_C_INCLUDE_PATH OR MPI_CXX_INCLUDE_PATH OR MPI_C_LIBRARIES OR MPI_CXX_LIBRARIES ) )

##############################################################################
#  MPI Stuff
##############################################################################

find_package ( MPI ${SCAI_FIND_PACKAGE_FLAGS} )

# find_package ( MPI ) does not find MPIEXEC when not in PATH
# do it ourself
if    ( NOT ${MPIEXEC} )
    find_program( MPIEXEC NAMES ${_MPI_EXEC_NAMES} HINTS ${MPI_ROOT}/bin )
endif ( NOT ${MPIEXEC} )

mark_as_advanced ( MPI_EXTRA_LIBRARY MPI_LIBRARY MPI_CXX_LIBRARIES MPIEXEC_POSTFLAGS MPIEXEC_PREFLAGS )

### ALLOW to switch off MPI explicitly ###

scai_build_variable ( NAME      USE_MPI
                      BOOL 
                      DEFAULT   ${MPI_FOUND}
                      DOCSTRING "use of MPI" )

if ( MPI_FOUND AND USE_MPI )

    set ( SCAI_MPI_INCLUDE_DIR ${MPI_INCLUDE_PATH} CACHE PATH "MPI include directory" )

    # some older versions of cmake have not set MPI_CXX_LIBRARIES

    if ( DEFINED MPI_CXX_LIBRARIES )
       set ( SCAI_MPI_LIBRARIES ${MPI_CXX_LIBRARIES} CACHE PATH "MPI libraries" )
    else ()
       set ( SCAI_MPI_LIBRARIES ${MPI_LIBRARIES} CACHE PATH "MPI libraries" )
    endif ()

endif ( MPI_FOUND AND USE_MPI )

include ( VersionCheck/MPI )

if    ( USE_MPI AND NOT MPI_FOUND )
    message( FATAL_ERROR "MPI shoud be used but not found" )
endif ( USE_MPI AND NOT MPI_FOUND )

mark_as_advanced ( SCAI_MPI_LIBRARIES )
mark_as_advanced ( SCAI_MPI_INCLUDE_DIR )
 
scai_summary_external ( NAME      MPI 
                        FOUND     ${MPI_FOUND} 
                        ENABLED   ${USE_MPI}
                        VERSION   ${MPI_VERSION} 
                        INCLUDE   ${SCAI_MPI_INCLUDE_DIR} 
                        LIBRARIES ${SCAI_MPI_LIBRARIES} )
