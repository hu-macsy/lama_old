/**
 * @file OpenMPLAPACKe.cpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief OpenMPLAPACKe.cpp
 * @author lschubert
 * @date 02.07.2012
 * $Id$
 */

// hpp
#include <lama/openmp/OpenMPLAPACKe.hpp>

// others
#include <lama/openmp/OpenMPBLAS1.hpp>

#include <lama/BLASInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/openmp/BLASHelper.hpp>
#include <boost/scoped_array.hpp>

#include <mkl_lapacke.h>

namespace lama
{

/* ------------------------------------------------------------------------- */

static int lapack_order( const enum CBLAS_ORDER order )
{
    int matrix_order = LAPACK_COL_MAJOR;

    if ( order == CblasColMajor )
    {
        matrix_order = LAPACK_COL_MAJOR;
    }
    else if ( order == CblasRowMajor )
    {
        matrix_order = LAPACK_ROW_MAJOR;
    }
    else
    {
        LAMA_THROWEXCEPTION( "Illegal matrix order: " << order );
    }

    return matrix_order;
}

/* ------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( OpenMPLAPACKe::logger, "OpenMP.LAPACKe" )

/* ------------------------------------------------------------------------- */
/*      getrf<float>                                                         */
/* ------------------------------------------------------------------------- */

template<>
IndexType OpenMPLAPACKe::getrf(
    const enum CBLAS_ORDER order,
    const int m,
    const int n,
    float* const A,
    const int lda,
    int* const ipiv )
{
    LAMA_LOG_INFO( logger, "getrf<float> for A of size " << m << " x " << n )

    int info = LAPACKE_sgetrf( lapack_order( order), m, n, A, lda, ipiv );

    if ( info < 0 )
    {
        LAMA_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if ( info > 0 )
    {
        LAMA_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
    }

    return info;
}

/* ------------------------------------------------------------------------- */
/*      getrf<double>                                                        */
/* ------------------------------------------------------------------------- */

template<>
IndexType OpenMPLAPACKe::getrf(
    const enum CBLAS_ORDER order,
    const int m,
    const int n,
    double* const A,
    const int lda,
    int* const ipiv )
{
    LAMA_LOG_INFO( logger, "getrf<double> for A of size " << m << " x " << n )

    int matrix_order = lapack_order( order );

    int info = LAPACKE_dgetrf( matrix_order, m, n, A, lda, ipiv );

    if ( info < 0 )
    {
        LAMA_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if ( info > 0 )
    {
        LAMA_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
    }

    return info;
}

/* ------------------------------------------------------------------------- */
/*      getinv<float>                                                        */
/* ------------------------------------------------------------------------- */

template<>
void OpenMPLAPACKe::getinv( const IndexType n, float* a, const IndexType lda )
{
    int info = 0;

    // scoped array, will also be freed in case of exception

    boost::scoped_array<IndexType> ipiv( new IndexType[n] );

    LAMA_LOG_INFO( logger, "getinv<float> for " << n << " x " << n << " matrix, uses MKL" )

    info = LAPACKE_sgetrf( LAPACK_COL_MAJOR, n, n, a, lda, ipiv.get() );

    // return error if factorization did not work

    if ( info )
    {
        LAMA_THROWEXCEPTION( "MKL sgetrf failed, info = " << info )
    }

    info = LAPACKE_sgetri( LAPACK_COL_MAJOR, n, a, lda, ipiv.get() );

    if ( info )
    {
        LAMA_THROWEXCEPTION( "MKL sgetri failed, info = " << info )
    }
}

/* ------------------------------------------------------------------------- */
/*      getinv<double>                                                       */
/* ------------------------------------------------------------------------- */

template<>
void OpenMPLAPACKe::getinv( const IndexType n, double* a, const IndexType lda )
{
    int info = 0;

    boost::scoped_array<IndexType> ipiv( new IndexType[n] );

    LAMA_LOG_INFO( logger, "getinv<double> for " << n << " x " << n << " matrix, uses MKL" )

    int matrix_order = LAPACK_COL_MAJOR;

    info = LAPACKE_dgetrf( matrix_order, n, n, a, lda, ipiv.get() );

    // return error if factorization did not work

    if ( info )
    {
        LAMA_THROWEXCEPTION( "MKL dgetrf failed, info = " << info )
    }

    info = LAPACKE_dgetri( matrix_order, n, a, lda, ipiv.get() );

    if ( info )
    {
        LAMA_THROWEXCEPTION( "MKL dgetri failed, info = " << info )
    }
}

/* ------------------------------------------------------------------------- */
/*      getri<float>                                                         */
/* ------------------------------------------------------------------------- */

template<>
int OpenMPLAPACKe::getri(
    const enum CBLAS_ORDER order,
    const int n,
    float* const a,
    const int lda,
    int* const ipiv )
{
    LAMA_LOG_INFO( logger, "getri<float> for A of size " << n << " x " << n )

    int matrix_order = lapack_order( order );

    int info = LAPACKE_sgetri( matrix_order, n, a, lda, ipiv );

    if ( info < 0 )
    {
        LAMA_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if ( info > 0 )
    {
        LAMA_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
    }

    return info;
}

/* ------------------------------------------------------------------------- */
/*      getri<double>                                                        */
/* ------------------------------------------------------------------------- */

template<>
int OpenMPLAPACKe::getri(
    const enum CBLAS_ORDER order,
    const int n,
    double* const a,
    const int lda,
    int* const ipiv )
{
    LAMA_LOG_INFO( logger, "getri<double> for A of size " << n << " x " << n )

    int matrix_order = lapack_order( order );

    int info = LAPACKE_dgetri( matrix_order, n, a, lda, ipiv );

    if ( info < 0 )
    {
        LAMA_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if ( info > 0 )
    {
        LAMA_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
    }

    return info;
}

/* ------------------------------------------------------------------------- */
/*      tptrs<float>                                                         */
/* ------------------------------------------------------------------------- */

template<>
int OpenMPLAPACKe::tptrs(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const int n,
    const int nrhs,
    const float* AP,
    float* B,
    const int ldb )
{
    char UL = BLASHelper::lapack_uplo( uplo );
    char TA = BLASHelper::lapack_transpose( trans );
    char DI = BLASHelper::lapack_diag( diag );

    int matrix_order = lapack_order( order );

    LAMA_LOG_INFO( logger, "tptrs<float>, n = " << n << ", nrhs = " << nrhs 
                   << ", order = " << matrix_order << ", UL = " << UL 
                   << ", TA = " << TA << ", DI = " << DI );

    LAMA_ASSERT_ERROR( ldb >= std::max( 1, n ), "ldb = " << ldb << " out of range" );

    int info = LAPACKE_stptrs( matrix_order, UL, TA, DI, n, nrhs, AP, B, ldb );

    return info;
}

/* ------------------------------------------------------------------------- */
/*      tptrs<double>                                                        */
/* ------------------------------------------------------------------------- */

template<>
int OpenMPLAPACKe::tptrs(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const int n,
    const int nrhs,
    const double* AP,
    double* B,
    const int ldb )
{
    char UL = BLASHelper::lapack_uplo( uplo );
    char TA = BLASHelper::lapack_transpose( trans );
    char DI = BLASHelper::lapack_diag( diag );

    int matrix_order = lapack_order( order );

    LAMA_LOG_INFO( logger, "tptrs<double>, n = " << n << ", nrhs = " << nrhs 
                   << ", order = " << matrix_order << ", UL = " << UL 
                   << ", TA = " << TA << ", DI = " << DI );

    int info = LAPACKE_dtptrs( matrix_order, UL, TA, DI, n, nrhs, AP, B, ldb );

    return info;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPLAPACKe::setInterface( BLASInterface& BLAS )
{
    // Note: macro takes advantage of same name for routines and type definitions 
    //       ( e.g. routine CUDABLAS1::sum<T> is set for BLAS::BLAS1::sum variable

    LAMA_INTERFACE_REGISTER_T( BLAS, getrf, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, getrf, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, getri, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, getri, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, getinv, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, getinv, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, tptrs, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, tptrs, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the LAPACK routines                               */
/* --------------------------------------------------------------------------- */

bool OpenMPLAPACKe::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::Host );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPLAPACKe::initialized = registerInterface();

} /* namespace lama */
