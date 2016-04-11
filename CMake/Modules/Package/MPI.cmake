###
 # @file PackageMPI.cmake
 #
 # @license
 # Copyright (c) 2009-2013
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
 # @brief findPackage and configuration of MPI
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

### MPI_FOUND            - if MPI is found
### USE_MPI              - if MPI is enabled
### SCAI_MPI_INCLUDE_DIR - MPI include directory
### SCAI_MPI_LIBRARIES   - all needed MPI libraries
### MPI_ENABLED          - if MPI_FOUND AND USE_MPI

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

mark_as_advanced ( MPI_EXTRA_LIBRARY MPI_LIBRARY )

### ALLOW to switch off MPI explicitly ###
include ( Functions/setAndCheckCache )
setAndCheckCache ( MPI )

set ( SCAI_MPI_INCLUDE_DIR ${MPI_INCLUDE_PATH} )
set ( SCAI_MPI_LIBRARIES ${MPI_LIBRARIES} )

include ( VersionCheck/MPI )

set ( MPI_ENABLED FALSE )
if    ( USE_MPI AND MPI_FOUND )
    set ( MPI_ENABLED TRUE )
endif ( USE_MPI AND MPI_FOUND )

if    ( USE_MPI AND NOT MPI_FOUND )
    message( FATAL_ERROR "Build of LAMA MPI enabled, but configuration is incomplete!")
endif ( USE_MPI AND NOT MPI_FOUND )
 