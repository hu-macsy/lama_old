/**
 * @file LAPACK_LAPACK.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief LAPACK_LAPACK.cpp
 * @author lschubert
 * @date 02.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/blaskernel/external/LAPACK_LAPACK.hpp>

// local library
#include <scai/blaskernel/external/BLAS_BLAS1.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/external/LAPACKWrapper.hpp>
#include <scai/blaskernel/cblas.hpp>

#include <scai/hmemo/Context.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/macros/unused.hpp>

#include <boost/preprocessor.hpp>

namespace scai
{

using common::unique_ptr;

namespace blaskernel
{

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( LAPACK_LAPACK::logger, "LAPACK.LAPACK" )

/* ------------------------------------------------------------------------- */
/*      getrf<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType LAPACK_LAPACK::getrf(
    const CBLAS_ORDER order,
    const IndexType m,
    const IndexType n,
    ValueType* const A,
    const IndexType lda,
    int* const ipiv )
{
    SCAI_REGION( "LAPACK.LAPACK.getrf<float>" )

    SCAI_LOG_INFO( logger, "getrf<float> for A of size " << m << " x " << n )

    int info = 0;

#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_info = info;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_info info
#endif

    if( order == CblasColMajor )
    {
        F77_info = LAPACKWrapper::getrf( F77_M, F77_N, A, F77_lda, ipiv );
    }
    else if( m == n && n == lda )
    {
        for( int i = 0; i < m; ++i )
        {
            for( int j = i + 1; j < n; ++j )
            {
                std::swap( A[i * n + j], A[j * m + i] );
            }
        }

        F77_info = LAPACKWrapper::getrf( F77_M, F77_N, A, F77_lda, ipiv );

        for( int i = 0; i < m; ++i )
        {
            for( int j = i + 1; j < n; ++j )
            {
                std::swap( A[i * n + j], A[j * m + i] );
            }
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "row major only supported for square matrices" );
    }

    for( int i = 0; i < m; ++i )
    {
        --ipiv[i]; // Fortran numbering from 1 to n ->  0 to n-1
    }

    if( info < 0 )
    {
        COMMON_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if( info > 0 )
    {
        COMMON_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
    }

    return info;
}

/* ------------------------------------------------------------------------- */
/*      getinv<float>                                                        */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACK_LAPACK::getinv( const IndexType n, ValueType* a, const IndexType lda )
{
    SCAI_REGION( "LAPACK.LAPACK.getinv<float>" )

    int info = 0;

    // unique_ptr, delete by destructor, also done in case of exception

    common::scoped_array<IndexType> ipiv( new IndexType[n] );

#ifdef F77_INT
    F77_INT F77_N = n, F77_lda = lda, F77_info = info;
#else
#define F77_N n
#define F77_lda lda
#define F77_info info
#endif

    SCAI_LOG_INFO( logger, "getinv<float> for " << n << " x " << n << " matrix, uses Fortran interface" )

    F77_info = LAPACKWrapper::getrf( F77_N, F77_N, a, F77_lda, ipiv.get() );

    if( F77_info )
    {
        COMMON_THROWEXCEPTION( "LAPACK sgetrf failed, info = " << F77_info )
    }

    common::scoped_array<ValueType> work( new ValueType[n] );

    F77_info = LAPACKWrapper::getri( F77_N, a, F77_lda, ipiv.get(), work.get(), F77_N );

    if( F77_info )
    {
        COMMON_THROWEXCEPTION( "LAPACK sgetri failed, info = " << F77_info )
    }
}

/* ------------------------------------------------------------------------- */
/*      getri<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
int LAPACK_LAPACK::getri( const CBLAS_ORDER order, const int n, ValueType* const a, const int lda, int* const ipiv )
{
    SCAI_REGION( "LAPACK.LAPACK.getri<float>" )

    SCAI_LOG_INFO( logger, "getri<float> for A of size " << n << " x " << n )

    int info = 0;

    // translate C indexes into  Fortran Indexes for ipiv

    for( int i = 0; i < n; ++i )
    {
        ++ipiv[i];
    }

    // transpose if not column major order

    if( order != CblasColMajor )
    {
        SCAI_ASSERT_EQUAL_ERROR( lda, n )

        for( int i = 0; i < n; ++i )
        {
            // swap row and column

            for( int j = i + 1; j < n; ++j )
            {
                std::swap( a[i * n + j], a[j * n + i] );
            }
        }
    }

#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_info = info;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_info info
#endif

    common::scoped_array<ValueType> work( new ValueType[n] );

    F77_info = LAPACKWrapper::getri( F77_N, a, F77_lda, ipiv, work.get(), F77_N );

    if( order != CblasColMajor )
    {
        // transpose back

        for( int i = 0; i < n; ++i )
        {
            for( int j = i + 1; j < n; ++j )
            {
                std::swap( a[i * n + j], a[j * n + i] );
            }
        }
    }

    if( info < 0 )
    {
        COMMON_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if( info > 0 )
    {
        COMMON_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
    }

    return info;
}

template<typename ValueType>
int LAPACK_LAPACK::tptrs(
    const CBLAS_ORDER order,
    const CBLAS_UPLO uplo,
    const CBLAS_TRANSPOSE trans,
    const CBLAS_DIAG diag,
    const IndexType n,
    const IndexType nrhs,
    const ValueType* AP,
    ValueType* B,
    const IndexType ldb )
{
    SCAI_REGION( "LAPACK.LAPACK.tptrs<float>" )

    int info = 0;

    char UL = BLASHelper::lapack_uplo( uplo );
    char TA = BLASHelper::lapack_transpose( trans );
    char DI = BLASHelper::lapack_diag( diag );

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

    SCAI_LOG_INFO( logger,
                   "tptrs<float>, n = " << n << ", nrhs = " << nrhs << ", order = " << order << ", UL = " << UL << ", TA = " << TA << ", DI = " << DI );

    if( order == CblasColMajor )
    {
        LAPACKWrapper::tptrs( UL, TA, DI, n, nrhs, AP, B, ldb );
    }
    else if( order == CblasRowMajor )
    {
        // TODO: transpose matrix.
        COMMON_THROWEXCEPTION( "row major order not supported for tptrs" );
    }

    return info;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACK_LAPACK::laswp(
    const CBLAS_ORDER order,
    const IndexType N,
    ValueType* A,
    const IndexType LDA,
    const IndexType K1,
    const IndexType K2,
    const IndexType* ipiv,
    const IndexType INCX )
{
    SCAI_REGION( "LAPACK.LAPACK.laswp<float>" )

    int i = K1;

    if( order == CblasRowMajor )
    {
        for( i = K1; i < K2; ++i )
        {
            if( ipiv[i * INCX] == i )
            {
                continue;
            }

            BLAS_BLAS1::swap<ValueType>( N, &A[ipiv[i * INCX] * LDA], INCX, &A[i * LDA], INCX );
        }
    }
    else if( order == CblasColMajor )
    {
#ifdef F77_INT
        F77_INT F77_N = N, F77_LDA = LDA, F77_K1 = K1, F77_K2 = K2, F77_INCX = INCX;
        F77_INT* F77_IPIV = ipiv;
#else
#   ifdef F77_N
#       undef F77_N
#   endif
#   define F77_N    N
#   define F77_LDA  LDA
#   define F77_K1   K1
#   define F77_K2   K2
#   define F77_INCX INCX
#   define F77_IPIV ipiv
#endif
        LAPACKWrapper::laswp( F77_N, A, F77_LDA, F77_K1, F77_K2, F77_IPIV, F77_INCX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( 0, 1, "cblas_slaswp", "Illegal order setting, %d\n", order );
    }
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the LAPACK routines                               */
/* --------------------------------------------------------------------------- */

void LAPACK_LAPACK::registerKernels( bool deleteFlag )
{
    using kregistry::KernelRegistry;
    using common::context::Host;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_REPLACE;   // priority over OpenMPBLAS

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

#define LAMA_LAPACK_REGISTER(z, I, _)                                                               \
    KernelRegistry::set<BLASKernelTrait::getrf<ARITHMETIC_HOST_TYPE_##I> >( getrf, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::getri<ARITHMETIC_HOST_TYPE_##I> >( getri, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::getinv<ARITHMETIC_HOST_TYPE_##I> >( getinv, Host, flag );  \
    KernelRegistry::set<BLASKernelTrait::tptrs<ARITHMETIC_HOST_TYPE_##I> >( tptrs, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::laswp<ARITHMETIC_HOST_TYPE_##I> >( laswp, Host, flag );

    BOOST_PP_REPEAT( ARITHMETIC_HOST_EXT_TYPE_CNT, LAMA_LAPACK_REGISTER, _ )

#undef LAMA_LAPACK_REGISTER
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

LAPACK_LAPACK::LAPACK_LAPACK()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

LAPACK_LAPACK::~LAPACK_LAPACK()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

LAPACK_LAPACK LAPACK_LAPACK::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
