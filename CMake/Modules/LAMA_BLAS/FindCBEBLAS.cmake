# - Try to find CBE BLAS Implementation 
#
# User can define CBEBLAS_LIBRARY_PATH
#
# Once done this will define 
# CBEBLAS_FOUND - System has LAMA_BLAS
# CBEBLAS_LIBRARIES - The libraries needed to use LAMA_BLAS

## Search CBEBLAS library

find_library ( CBEBLAS_LIBRARY blas HINTS ${CBEBLAS_LIBRARY_PATH} )
if ( EXISTS ${CBEBLAS_LIBRARY} )
    set ( CBEBLAS_LIBRARIES ${CBEBLAS_LIBRARY} )
elseif ( NOT DEFINED CBEBLAS_LIBRARY_PATH )
    message ( STATUS "WARNING CBEBLAS not found. Please define CBEBLAS_LIBRARY_PATH." )
else ()
    message ( STATUS "WARNING CBEBLAS not found. CBEBLAS_LIBRARY_PATH=${CBEBLAS_LIBRARY_PATH} directory does not exist." )
endif ()

## Module footer

include ( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set CBEBLAS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args ( CBEBLAS DEFAULT_MSG CBEBLAS_LIBRARIES )

if( CBEBLAS_FOUND )
   add_definitions ( -DLAMA_FORTRAN_BLAS_STYLE_UNDERSCORE )
endif( CBEBLAS_FOUND )

mark_as_advanced( CBEBLAS_LIBRARIES CBEBLAS_LIBRARY )