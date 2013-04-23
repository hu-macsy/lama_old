# Look for MPI first to allow LAMA_BLAS to take the correct blacs implementation
# based on the found mpi
if ( WIN32 AND NOT ( MPI_C_INCLUDE_PATH OR MPI_CXX_INCLUDE_PATH OR MPI_C_LIBRARIES OR MPI_CXX_LIBRARIES ) )
    if ( MPI_ROOT )
        set ( MPI_C_INCLUDE_PATH "${MPI_ROOT}/Inc" )
        set ( MPI_CXX_INCLUDE_PATH "${MPI_ROOT}/Inc" )
        set ( LAMA_MPI_LIB_DIR "${MPI_ROOT}/Lib" )
    else ( MPI_ROOT )
        set ( MPI_C_INCLUDE_PATH "C:/Program Files/Microsoft HPC Pack 2008 R2/Inc" )
        set ( MPI_CXX_INCLUDE_PATH "C:/Program Files/Microsoft HPC Pack 2008 R2/Inc" )
        set ( LAMA_MPI_LIB_DIR "C:/Program Files/Microsoft HPC Pack 2008 R2/Lib" )
    endif ( MPI_ROOT )
    if ( CMAKE_CL_64 )
        set ( LAMA_MPI_LIB_DIR "${LAMA_MPI_LIB_DIR}/amd64" )
    else ( CMAKE_CL_64 )
        set ( LAMA_MPI_LIB_DIR "${LAMA_MPI_LIB_DIR}/i386" )
    endif ( CMAKE_CL_64 )
    set ( MPI_C_LIBRARIES "${LAMA_MPI_LIB_DIR}/msmpi.lib" )
    set ( MPI_CXX_LIBRARIES "${LAMA_MPI_LIB_DIR}/msmpi.lib" )
else ( WIN32 AND NOT ( MPI_C_INCLUDE_PATH OR MPI_CXX_INCLUDE_PATH OR MPI_C_LIBRARIES OR MPI_CXX_LIBRARIES ) )
  if ( MPI_ROOT )
      set ( MPI_COMPILER ${MPI_ROOT}/bin/mpicxx )
  endif ( MPI_ROOT )
endif ( WIN32 AND NOT ( MPI_C_INCLUDE_PATH OR MPI_CXX_INCLUDE_PATH OR MPI_C_LIBRARIES OR MPI_CXX_LIBRARIES ) )

##############################################################################
#  MPI Stuff
##############################################################################

find_package ( MPI ${LAMA_FIND_PACKAGE_FLAGS} )

### ALLOW to switch off MPI explicitly ###

if ( NOT LAMA_USE_MPI )
    set ( DEFAULT_USE_MPI ${MPI_FOUND} )
endif ( NOT LAMA_USE_MPI )

set ( LAMA_USE_MPI ${DEFAULT_USE_MPI} CACHE BOOL "Enable / Disable use of MPI" )

if ( MPI_FOUND AND LAMA_USE_MPI)
    include_directories ( ${MPI_INCLUDE_PATH} )
endif ( MPI_FOUND AND LAMA_USE_MPI)