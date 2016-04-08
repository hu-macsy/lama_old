# LAMA (core)
heading ( "External Libraries:" )

set ( REQUIRED_FOUND FALSE )
if    ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR AND SCAI_BLAS_FOUND )
    set ( REQUIRED_FOUND TRUE )
    if ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
        set( REQUIRED_FOUND FALSE )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
endif ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR AND SCAI_BLAS_FOUND )

heading2 ( "Required core" "REQUIRED_FOUND" )

    # pthreads
    found_message ( "pThreads" "SCAI_THREAD_LIBRARIES" "REQUIRED" "Version ${SCAI_THREAD_VERSION}" )
    # boost
    found_message ( "Boost" "SCAI_BOOST_INCLUDE_DIR" "REQUIRED" "Version ${BOOST_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )

    # BLAS (Lapack)
    found_message ( "BLAS" "SCAI_BLAS_FOUND" "REQUIRED" "${SCAI_BLAS_NAME} Version ${BLAS_VERSION} with:" )
    foreach    ( _B_LIB ${SCAI_SCAI_BLAS_LIBRARIES} )
        message ( STATUS "                                 ${_B_LIB}" )
    endforeach ( _B_LIB ${SCAI_SCAI_BLAS_LIBRARIES} )
    if    ( SCAI_BLAS_NAME MATCHES "BLAS" )
        found_message ( "Lapack" "LAPACK_FOUND" "REQUIRED" "" )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" )