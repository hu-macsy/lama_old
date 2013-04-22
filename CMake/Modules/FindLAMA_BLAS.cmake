# - Try to find LAMA_BLAS
#   Once done this will define
#   LAMA_BLAS_FOUND - System has LAMA_BLAS
#   LAMA_BLAS_LIBRARIES - The libraries needed to use LAMA_BLAS

# try to find one of this blas libraries in this order: MKL, ACML, GOTOBLAS, FortranBLAS
include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindMKL.cmake )
setAndCheckCache( "MKL" )
	
if ( NOT MKL_FOUND )
	include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindACML.cmake )
    setAndCheckCache( "ACML" )
	
    if ( NOT ACML_FOUND )
	    include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindGOTOBLAS.cmake )
        setAndCheckCache( "GOTOBLAS" )
	
        if ( NOT GOTOBLAS_FOUND )
	        include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindFortranBLAS.cmake )
	        setAndCheckCache( "BLAS" )
	        
        endif ( NOT GOTOBLAS_FOUND )
    endif ( NOT ACML_FOUND )
endif ( NOT MKL_FOUND )
	
# at this point one blas library should be found
if ( MKL_FOUND )
	set ( LAMA_BLAS_LIBRARIES ${MKL_LIBRARIES} )
	set ( LAMA_PBLAS_LIBRARIES ${MKL_PLIBRARIES} )
		
elseif ( ACML_FOUND )
	set ( LAMA_BLAS_LIBRARIES ${ACML_LIBRARIES} )
		
elseif ( GOTOBLAS_FOUND ) 
   	set ( LAMA_BLAS_LIBRARIES ${GOTOBLAS_LIBRARIES} )
   	
#BLAS AND LAPACK found for FortranBLAS
elseif ( BLAS_FOUND )
    
    if ( LAPACK_LIBRARIES )
        set ( LAMA_BLAS_LIBRARIES ${LAPACK_LIBRARIES} )
    elseif ( LAPACK_lapack_LIBRARY AND BLAS_blas_LIBRARY )
        list ( APPEND LAMA_BLAS_LIBRARIES ${LAPACK_lapack_LIBRARY} )
    endif ()
    
endif ()    

message ( STATUS "LAMA_BLAS_LIBRARIES ${LAMA_BLAS_LIBRARIES}" )

include ( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set LAMA_BLAS_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args ( LAMA_BLAS DEFAULT_MSG LAMA_BLAS_LIBRARIES)

mark_as_advanced ( LAMA_BLAS_LIBRARIES )
mark_as_advanced ( LAMA_PBLAS_LIBRARIES )
