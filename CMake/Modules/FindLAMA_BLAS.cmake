# - Try to find LAMA_BLAS
#   Once done this will define
#   LAMA_BLAS_FOUND - System has LAMA_BLAS
#   LAMA_BLAS_LIBRARIES - The libraries needed to use LAMA_BLAS

if ( ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc64" AND CBEBLAS_FOUND )
	# cell processor

	include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindCBEBLAS.cmake )
	if ( CBEBLAS_FOUND )
		set ( LAMA_BLAS_LIBRARIES ${CBEBLAS_LIBRARIES} )
	endif ( CBEBLAS_FOUND )
	
else ( ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc64" AND CBEBLAS_FOUND )
	# no cell processor

	# try to find one of this blas libraries in this order: MKL, ACML, GOTOBLAS, FortranBLAS
	include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindMKL.cmake )
	
	if ( NOT MKL_FOUND )
			include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindACML.cmake )
	endif ( NOT MKL_FOUND )
	
	if ( NOT MKL_FOUND AND NOT ACML_FOUND )
		include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindGOTOBLAS.cmake )
	endif ( NOT MKL_FOUND AND NOT ACML_FOUND )
	
	if ( NOT MKL_FOUND AND NOT ACML_FOUND AND NOT GOTOBLAS_FOUND )
		include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindFortranBLAS.cmake )
	endif ( NOT MKL_FOUND AND NOT ACML_FOUND AND NOT GOTOBLAS_FOUND )
	
	# at this point one blas library should be found
	
	if ( MKL_FOUND )
	   message ( STATUS "BLAS library found: MKL" )
		# include is not required because lama uses its own header files 
		set ( LAMA_BLAS_LIBRARIES ${MKL_LIBRARIES} )
		set ( LAMA_PBLAS_LIBRARIES ${MKL_PLIBRARIES} )
		
	elseif ( ACML_FOUND )
		set ( LAMA_BLAS_LIBRARIES ${ACML_LIBRARIES} )
		
	elseif ( GOTOBLAS_FOUND ) 
    	set ( LAMA_BLAS_LIBRARIES ${GOTOBLAS_LIBRARIES} )
    	
	#BLAS AND LAPACK found for FortranBLAS
    elseif ( BLAS_FOUND )
    
        if ( LAPACK_LIBRARIES )
            list ( APPEND LAMA_BLAS_LIBRARIES ${LAPACK_LIBRARIES} )
            list ( APPEND LAMA_BLAS_LIBRARIES ${BLAS_LIBRARIES} )
            
        elseif ( LAPACK_lapack_LIBRARY AND BLAS_blas_LIBRARY )
            list ( APPEND LAMA_BLAS_LIBRARIES ${LAPACK_lapack_LIBRARY} )
            list ( APPEND LAMA_BLAS_LIBRARIES ${BLAS_blas_LIBRARY} )   
        endif ()
    endif ()    
endif( ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc64" AND CBEBLAS_FOUND ) # whether cell processor or not

include ( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set LAMA_BLAS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args ( LAMA_BLAS DEFAULT_MSG LAMA_BLAS_LIBRARIES)

mark_as_advanced ( LAMA_BLAS_LIBRARIES )
mark_as_advanced ( LAMA_PBLAS_LIBRARIES )
