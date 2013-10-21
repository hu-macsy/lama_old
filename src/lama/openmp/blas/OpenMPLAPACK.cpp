/**
 * @file OpenMPLAPACK.cpp
 *
 * @license
 * Copyright (c) 2009-2013
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief OpenMPLAPACK.cpp
 * @author lschubert
 * @date 02.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/openmp/OpenMPLAPACK.hpp>

// others
#include <lama/openmp/OpenMPBLAS1.hpp>

#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/BLASInterface.hpp>
#include <lama/tracing.hpp>
#include <boost/scoped_array.hpp>

// macros
#include <lama/macros/unused.hpp>

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/



#ifdef __cplusplus
} /*extern "C"*/
#endif /*__cplusplus*/

namespace lama
{

/* ------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( OpenMPLAPACK::logger, "OpenMP.LAPACK" )

/* ------------------------------------------------------------------------- */
/*      getrf<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename T>
IndexType OpenMPLAPACK::getrf(
    const enum CBLAS_ORDER UNUSED(order),
    const int m,
    const int n,
    T* const UNUSED(A),
    const int UNUSED(lda),
    int* const UNUSED(ipiv) )
{
    LAMA_REGION( "OpenMP.LAPACK.getrf<float>" )

    LAMA_LOG_INFO( logger, "getrf<float> for A of size " << m << " x " << n )

    int info = 0;
    LAMA_THROWEXCEPTION( "OpenMPLAPACK:getrf is not implemented yet")
//#ifdef F77_INT
//    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_info = info;
//#else
//#define F77_M m
//#define F77_N n
//#define F77_lda lda
//#define F77_info info
//#endif

//    if ( order == CblasColMajor )
//    {
//        F77_sgetrf( &F77_M, &F77_N, A, &F77_lda, ipiv, &F77_info );
//    }
//    else if ( m == n && n == lda )
//    {
//        for ( int i = 0; i < m; ++i )
//        {
//            for ( int j = i + 1; j < n; ++j )
//            {
//                std::swap( A[i * n + j], A[j * m + i] );
//            }
//        }
//
//        F77_sgetrf( &F77_M, &F77_N, A, &F77_lda, ipiv, &F77_info );
//
//        for ( int i = 0; i < m; ++i )
//        {
//            for ( int j = i + 1; j < n; ++j )
//            {
//                std::swap( A[i * n + j], A[j * m + i] );
//            }
//        }
//    }
//    else
//    {
//        LAMA_THROWEXCEPTION( "row major only supported for square matrices" );
//    }
//
//    for ( int i = 0; i < m; ++i )
//    {
//        --ipiv[i];   // Fortran numbering from 1 to n ->  0 to n-1
//    }
//
//    if ( info < 0 )
//    {
//        LAMA_THROWEXCEPTION( "illegal argument " << ( -info ) )
//    }
//    else if ( info > 0 )
//    {
//        LAMA_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
//    }

    return info;
}


/* ------------------------------------------------------------------------- */
/*      getinv<float>                                                        */
/* ------------------------------------------------------------------------- */

template<typename T>
void OpenMPLAPACK::getinv( const IndexType n, T* UNUSED(a), const IndexType UNUSED(lda) )
{
    LAMA_REGION( "OpenMP.LAPACK.getinv<float>" )

//    int info = 0;
    LAMA_THROWEXCEPTION( "OpenMPLAPACK:getinv is not implemented yet")
    // scoped array, will also be freed in case of exception

    boost::scoped_array<IndexType> ipiv( new IndexType[n] );
//
//#ifdef F77_INT
//    F77_INT F77_N = n, F77_lda = lda, F77_info = info;
//#else
//#define F77_N n
//#define F77_lda lda
//#define F77_info info
//#endif
//
//    LAMA_LOG_INFO( logger, "getinv<float> for " << n << " x " << n << " matrix, uses Fortran interface" )
//
//    F77_sgetrf( &F77_N, &F77_N, a, &F77_lda, ipiv.get(), &F77_info );
//
//    if ( F77_info )
//    {
//        LAMA_THROWEXCEPTION( "LAPACK sgetrf failed, info = " << F77_info )
//    }
//
//    boost::scoped_array<float> work( new float[n] );
//
//    F77_sgetri( &F77_N, a, &F77_lda, ipiv.get(), work.get(), &F77_N, &F77_info );
//
//    if ( F77_info )
//    {
//        LAMA_THROWEXCEPTION( "LAPACK sgetri failed, info = " << F77_info )
//    }
}


/* ------------------------------------------------------------------------- */
/*      getri<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename T>
int OpenMPLAPACK::getri(
    const enum CBLAS_ORDER UNUSED(order),
    const int n,
    T* const UNUSED(a),
    const int UNUSED(lda),
    int* const UNUSED(ipiv) )
{
    LAMA_REGION( "OpenMP.LAPACK.getri<float>" )

    LAMA_LOG_INFO( logger, "getri<float> for A of size " << n << " x " << n )

    int info = 0;
    LAMA_THROWEXCEPTION( "OpenMPLAPACK:getri is not implemented yet")
    // translate C indexes into  Fortran Indexes for ipiv

//    for ( int i = 0; i < n; ++i )
//    {
//        ++ipiv[i];
//    }
//
//    // transpose if not column major order
//
//    if ( order != CblasColMajor )
//    {
//        LAMA_ASSERT_EQUAL_ERROR( lda, n )
//
//        for ( int i = 0; i < n; ++i )
//        {
//            // swap row and column
//
//            for ( int j = i + 1; j < n; ++j )
//            {
//                std::swap( a[ i * n + j], a[ j * n + i ] );
//            }
//        }
//    }
//
//#ifdef F77_INT
//    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_info = info;
//#else
//#define F77_M m
//#define F77_N n
//#define F77_lda lda
//#define F77_info info
//#endif
//
//    boost::scoped_array<float> work( new float[n] );
//
//    F77_sgetri( &F77_N, a, &F77_lda, ipiv, work.get(), &F77_N, &F77_info );
//
//    if ( order != CblasColMajor )
//    {
//        // transpose back
//
//        for ( int i = 0; i < n; ++i )
//        {
//            for ( int j = i + 1; j < n; ++j )
//            {
//                std::swap( a[ i * n + j], a[ j * n + i ] );
//            }
//        }
//    }
//
//    if ( info < 0 )
//    {
//        LAMA_THROWEXCEPTION( "illegal argument " << ( -info ) )
//    }
//    else if ( info > 0 )
//    {
//        LAMA_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
//    }

    return info;
}


/* --------------------------------------------------------------------------- */

template<typename T>
int OpenMPLAPACK::tptrs(
    const enum CBLAS_ORDER UNUSED(order),
    const enum CBLAS_UPLO UNUSED(uplo),
    const enum CBLAS_TRANSPOSE UNUSED(trans),
    const enum CBLAS_DIAG UNUSED(diag),
    const int UNUSED(n),
    const int UNUSED(nrhs),
    const T* UNUSED(AP),
    T* UNUSED(B),
    const int UNUSED(ldb) )
{
    LAMA_REGION( "OpenMP.LAPACK.tptrs<float>" )

    int info = 0;
    LAMA_THROWEXCEPTION( "OpenMPLAPACK:tptrs is not implemented yet")
//    char UL = BLASHelper::lapack_uplo( uplo );
//    char TA = BLASHelper::lapack_transpose( trans );
//    char DI = BLASHelper::lapack_diag( diag );
//
//#ifdef F77_CHAR
//    F77_CHAR F77_UL, F77_TA, F77_DI;
//#else
//#   define F77_UL &UL
//#   define F77_TA &TA
//#   define F77_DI &DI
//#endif
//#ifdef F77_INT
//    F77_INT F77_n = n, F77_nrhs = nrhs, F77_ldb = ldb;
//#else
//#   define F77_n n
//#   define F77_nrhs nrhs
//#   define F77_ldb ldb
//#endif
//
//#ifdef F77_CHAR
//    F77_UL = C2F_CHAR( &UL );
//    F77_TA = C2F_CHAR( &TA );
//    F77_DI = C2F_CHAR( &DI );
//#endif
//
//    LAMA_LOG_INFO( logger, "tptrs<float>, n = " << n << ", nrhs = " << nrhs
//                   << ", order = " << order << ", UL = " << UL
//                   << ", TA = " << TA << ", DI = " << DI );
//
//    if ( order == CblasColMajor )
//    {
//        F77_stptrs( F77_UL, F77_TA, F77_DI, &F77_n, &F77_nrhs, AP, B, &F77_ldb, &info );
//    }
//    else if ( order == CblasRowMajor )
//    {
//        // TODO: transpose matrix.
//        LAMA_THROWEXCEPTION( "row major order not supported for tptrs" );
//    }

    return info;
}

/* --------------------------------------------------------------------------- */

template<typename T>
void OpenMPLAPACK::laswp(
    const enum CBLAS_ORDER UNUSED(order),
    const int UNUSED(N),
    T* UNUSED(A),
    const int UNUSED(LDA),
    const int UNUSED(K1),
    const int UNUSED(K2),
    const int* UNUSED(ipiv),
    const int UNUSED(INCX),
    SyncToken* UNUSED(syncToken) )
{
    LAMA_REGION( "OpenMP.LAPACK.laswp<float>" )

    //int i = K1;
    LAMA_THROWEXCEPTION( "OpenMPLAPACK:laswp is not implemented yet")
    if ( order == CblasRowMajor )
    {
        for ( i = K1; i < K2; ++i )
        {
            if ( ipiv[i * INCX] == i )
            {
                continue;
            }

            OpenMPBLAS1::swap<float>( N, &A[ipiv[i * INCX] * LDA], INCX, &A[i * LDA], INCX, syncToken );
        }
    }
    else if ( order == CblasColMajor )
    {
//#ifdef F77_INT
//        F77_INT F77_N = N, F77_LDA = LDA, F77_K1 = K1, F77_K2 = K2, F77_INCX = INCX;
//        F77_INT* F77_IPIV = ipiv;
//#else
//#   ifdef F77_N
//#       undef F77_N
//#   endif
//#   define F77_N    N
//#   define F77_LDA  LDA
//#   define F77_K1   K1
//#   define F77_K2   K2
//#   define F77_INCX INCX
//#   define F77_IPIV ipiv
//#endif
//        F77_slaswp( &F77_N, A, &F77_LDA, &F77_K1, &F77_K2, F77_IPIV, &F77_INCX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( 0, 1, "cblas_slaswp", "Illegal order setting, %d\n", order );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPLAPACK::setInterface( BLASInterface& BLAS )
{
    // Note: macro takes advantage of same name for routines and type definitions 
    //       ( e.g. routine CUDABLAS1::sum<T> is set for BLAS::BLAS1::sum variable

//    LAMA_INTERFACE_REGISTER_T( BLAS, getrf, float )
//    LAMA_INTERFACE_REGISTER_T( BLAS, getrf, double )
//
//    LAMA_INTERFACE_REGISTER_T( BLAS, getri, float )
//    LAMA_INTERFACE_REGISTER_T( BLAS, getri, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, getinv, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, getinv, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, tptrs, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, tptrs, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, laswp, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, laswp, double )

    // other routines are not used by LAMA yet
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the LAPACK routines                               */
/* --------------------------------------------------------------------------- */

bool OpenMPLAPACK::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::Host );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPLAPACK::initialized = registerInterface();

} /* namespace lama */
