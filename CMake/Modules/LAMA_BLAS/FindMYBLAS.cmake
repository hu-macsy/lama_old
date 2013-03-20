# - Try to find MYBLAS
#
# User should define MYBLAS_LIBRARY_PATH and GFORTRAN_LIBRARY_PATH
#
# Once done this will define 
#  MYBLAS_FOUND - System has LAMA_BLAS
#  MYBLAS_LIBRARIES - The libraries needed to use LAMA_BLAS

## Search MYBLAS library

set ( MYBLAS_LIBRARY_PATH "/usr/lib/" )

find_library ( MYBLAS_LIBRARY blas HINTS ${MYBLAS_LIBRARY_PATH} )
if ( EXISTS ${MYBLAS_LIBRARY} )
    set ( MYBLAS_LIBRARIES ${MYBLAS_LIBRARY} )
    
    ## Search for gfortran. Required by MYBLAS

    find_library ( GFORTRAN_LIBRARY gfortran HINTS ${GFORTRAN_LIBRARY_PATH} )
    
    if ( EXISTS ${GFORTRAN_LIBRARY} )
        list ( APPEND MYBLAS_LIBRARIES ${GFORTRAN_LIBRARY} )
    else( EXISTS ${GFORTRAN_LIBRARY} )
        message ( STATUS "WARNING Library gfortran not found. Required by MY BLAS. Please define GFORTRAN_LIBRARY_PATH." )
    endif( EXISTS ${GFORTRAN_LIBRARY} )
    
    ## Search for pthreads. Required by MYBLAS
    
    set( CMAKE_THREAD_PREFER_PTHREAD TRUE )
    
    find_package ( Threads QUIET )
    
    if ( Threads_FOUND )
        if ( CMAKE_USE_PTHREADS_INIT )
            #TODO is this correct???
            list ( APPEND MYBLAS_LIBRARIES ${CMAKE_THREAD_LIBS_INIT} )
            #add_definitions( ${CMAKE_THREAD_LIBS_INIT} ) 
        else( CMAKE_USE_PTHREADS_INIT )
            message ( STATUS "WARNING No pthreads found. Required by MYBLAS." )    
        endif( CMAKE_USE_PTHREADS_INIT )
    else( Threads_FOUND )
        message ( STATUS "WARNING No threads found. Required by MYBLAS." )
    endif( Threads_FOUND )
elseif ( NOT DEFINED MYBLAS_LIBRARY_PATH )
    message ( STATUS "WARNING MYBLAS not found. Please define MYBLAS_LIBRARY_PATH." )
else ( EXISTS ${MYBLAS_LIBRARY} )
    message ( STATUS "WARNING MYBLAS not found. MYBLAS_LIBRARY_PATH=${MYBLAS_LIBRARY_PATH} directory does not exist." )
endif ( EXISTS ${MYBLAS_LIBRARY} )

## Module footer

include ( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set MYBLAS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args ( MYBLAS DEFAULT_MSG MYBLAS_LIBRARIES )

if ( MYBLAS_FOUND )
   add_definitions ( -DLAMA_FORTRAN_BLAS_STYLE_UNDERSCORE )
endif( MYBLAS_FOUND )

mark_as_advanced( MYBLAS_LIBRARIES MYBLAS_LIBRARY )

find_package( LAPACK )
