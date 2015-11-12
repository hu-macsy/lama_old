/**
 * @file LAPACKe_LAPACK.cpp
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
 * @brief LAPACKe_LAPACK.cpp
 * @author Lauretta Schubert
 * @date 02.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/blaskernel/external/LAPACKe_LAPACK.hpp>

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/external/BLASHelper.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/Assert.hpp>

// external
#include <mkl_lapacke.h>

// boost
#include <boost/preprocessor.hpp>
#include "../../blaskernel/openmp/OpenMPBLAS1.hpp"

namespace scai
{

using common::scoped_array;

namespace blaskernel
{

/* ---------------------------------------------------------------------------------------*/

/** mklCast converts pointers to LAMA complex numbers to
 *  ML pointers for complex numbers. This is safe as both
 *  are internally represented in the same way.
 */

static inline MKL_Complex8* mklCast( ComplexFloat* x )
{
    return reinterpret_cast<MKL_Complex8*>( x );
}

static inline const MKL_Complex8* mklCast( const ComplexFloat* x )
{
    return reinterpret_cast<const MKL_Complex8*>( x );
}

static inline MKL_Complex16* mklCast( ComplexDouble* x )
{
    return reinterpret_cast<MKL_Complex16*>( x );
}

static inline const MKL_Complex16* mklCast( const ComplexDouble* x )
{
    return reinterpret_cast<const MKL_Complex16*>( x );
}

/* ------------------------------------------------------------------------- */

static int lapack_order( const CBLAS_ORDER order )
{
    int matrix_order = LAPACK_COL_MAJOR;

    if( order == CblasColMajor )
    {
        matrix_order = LAPACK_COL_MAJOR;
    }
    else if( order == CblasRowMajor )
    {
        matrix_order = LAPACK_ROW_MAJOR;
    }
    else
    {
        COMMON_THROWEXCEPTION( "Illegal matrix order: " << order );
    }

    return matrix_order;
}

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( LAPACKe_LAPACK::logger, "LAPACKe.LAPACK" )

/* ------------------------------------------------------------------------- */
/*      getrf<float>                                                         */
/* ------------------------------------------------------------------------- */

template<>
IndexType LAPACKe_LAPACK::getrf(
    const CBLAS_ORDER order,
    const int m,
    const int n,
    float* const A,
    const int lda,
    int* const ipiv )
{
    SCAI_LOG_INFO( logger, "getrf<float> for A of size " << m << " x " << n )

    int info = LAPACKE_sgetrf( lapack_order( order ), m, n, A, lda, ipiv );

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
/*      getrf<double>                                                        */
/* ------------------------------------------------------------------------- */

template<>
IndexType LAPACKe_LAPACK::getrf(
    const CBLAS_ORDER order,
    const int m,
    const int n,
    double* const A,
    const int lda,
    int* const ipiv )
{
    SCAI_LOG_INFO( logger, "getrf<double> for A of size " << m << " x " << n )

    int matrix_order = lapack_order( order );

    int info = LAPACKE_dgetrf( matrix_order, m, n, A, lda, ipiv );

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
/*      getrf<ComplexFloat>                                                  */
/* ------------------------------------------------------------------------- */

template<>
IndexType LAPACKe_LAPACK::getrf(
    const CBLAS_ORDER order,
    const int m,
    const int n,
    ComplexFloat* const a,
    const int lda,
    int* const ipiv )
{
    SCAI_LOG_INFO( logger, "getrf<ComplexFloat> for a of size " << m << " x " << n )

    int info = LAPACKE_cgetrf( lapack_order( order ), m, n, mklCast( a ), lda, ipiv );

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
/*      getrf<ComplexDouble>                                                 */
/* ------------------------------------------------------------------------- */

template<>
IndexType LAPACKe_LAPACK::getrf(
    const CBLAS_ORDER order,
    const int m,
    const int n,
    ComplexDouble* const a,
    const int lda,
    int* const ipiv )
{
    SCAI_LOG_INFO( logger, "getrf<ComplexFloat> for a of size " << m << " x " << n )

    int info = LAPACKE_zgetrf( lapack_order( order ), m, n, mklCast( a ), lda, ipiv );

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
/*      getrf<LongDouble>                                                    */
/* ------------------------------------------------------------------------- */

template<>
IndexType LAPACKe_LAPACK::getrf( const CBLAS_ORDER, const int, const int, LongDouble* const, const int, int* const )
{
    // LongDouble not supported by LAPACK, should give a serious runtime error
    COMMON_THROWEXCEPTION( "getrf<LongDouble> unsupported" )
}

/* ------------------------------------------------------------------------- */
/*      getrf<ComplexLongDouble>                                                    */
/* ------------------------------------------------------------------------- */

template<>
IndexType LAPACKe_LAPACK::getrf( const CBLAS_ORDER, const int, const int, ComplexLongDouble* const, const int, int* const )
{
    // ComplexLongDouble not supported by LAPACK, should give a serious runtime error
    COMMON_THROWEXCEPTION( "getrf<ComplexLongDouble> unsupported" )
}

/* ------------------------------------------------------------------------- */
/*      getinv<float>                                                        */
/* ------------------------------------------------------------------------- */

template<>
void LAPACKe_LAPACK::getinv( const IndexType n, float* a, const IndexType lda )
{
    int info = 0;

    // scoped_array, will also be freed in case of exception

    scoped_array<IndexType> ipiv( new IndexType[n] );

    SCAI_LOG_INFO( logger, "getinv<float> for " << n << " x " << n << " matrix, uses MKL" )

    info = LAPACKE_sgetrf( LAPACK_COL_MAJOR, n, n, a, lda, ipiv.get() );

    // return error if factorization did not work

    if( info )
    {
        COMMON_THROWEXCEPTION( "MKL sgetrf failed, info = " << info )
    }

    info = LAPACKE_sgetri( LAPACK_COL_MAJOR, n, a, lda, ipiv.get() );

    if( info )
    {
        COMMON_THROWEXCEPTION( "MKL sgetri failed, info = " << info )
    }
}

/* ------------------------------------------------------------------------- */
/*      getinv<double>                                                       */
/* ------------------------------------------------------------------------- */

template<>
void LAPACKe_LAPACK::getinv( const IndexType n, double* a, const IndexType lda )
{
    int info = 0;

    scoped_array<IndexType> ipiv( new IndexType[n] );

    SCAI_LOG_INFO( logger, "getinv<double> for " << n << " x " << n << " matrix, uses MKL" )

    int matrix_order = LAPACK_COL_MAJOR;

    info = LAPACKE_dgetrf( matrix_order, n, n, a, lda, ipiv.get() );

    // return error if factorization did not work

    if( info )
    {
        COMMON_THROWEXCEPTION( "MKL dgetrf failed, info = " << info )
    }

    info = LAPACKE_dgetri( matrix_order, n, a, lda, ipiv.get() );

    if( info )
    {
        COMMON_THROWEXCEPTION( "MKL dgetri failed, info = " << info )
    }
}

/* ------------------------------------------------------------------------- */
/*      getinv<ComplexFloat>                                                        */
/* ------------------------------------------------------------------------- */

template<>
void LAPACKe_LAPACK::getinv( const IndexType n, ComplexFloat* a, const IndexType lda )
{
    int info = 0;

    // scoped_array: delete by destructor, works also for exceptions

    scoped_array<IndexType> ipiv( new IndexType[n] );

    SCAI_LOG_INFO( logger, "getinv<ComplexFloat> for " << n << " x " << n << " matrix, uses MKL" )

    info = LAPACKE_cgetrf( LAPACK_COL_MAJOR, n, n, mklCast( a ), lda, ipiv.get() );

    // return error if factorization did not work

    if( info )
    {
        COMMON_THROWEXCEPTION( "MKL sgetrf failed, info = " << info )
    }

    info = LAPACKE_cgetri( LAPACK_COL_MAJOR, n, mklCast( a ), lda, ipiv.get() );

    if( info )
    {
        COMMON_THROWEXCEPTION( "MKL sgetri failed, info = " << info )
    }
}

/* ------------------------------------------------------------------------- */
/*      getinv<ComplexDouble>                                                */
/* ------------------------------------------------------------------------- */

template<>
void LAPACKe_LAPACK::getinv( const IndexType n, ComplexDouble* a, const IndexType lda )
{
    int info = 0;

    scoped_array<IndexType> ipiv( new IndexType[n] );

    SCAI_LOG_INFO( logger, "getinv<ComplexDouble> for " << n << " x " << n << " matrix, uses MKL" )

    info = LAPACKE_zgetrf( LAPACK_COL_MAJOR, n, n, mklCast( a ), lda, ipiv.get() );

    // return error if factorization did not work

    if( info )
    {
        COMMON_THROWEXCEPTION( "MKL sgetrf failed, info = " << info )
    }

    info = LAPACKE_zgetri( LAPACK_COL_MAJOR, n, mklCast( a ), lda, ipiv.get() );

    if( info )
    {
        COMMON_THROWEXCEPTION( "MKL sgetri failed, info = " << info )
    }
}

/* ------------------------------------------------------------------------- */
/*      getinv<LongDouble>                                                   */
/* ------------------------------------------------------------------------- */

template<>
void LAPACKe_LAPACK::getinv( const IndexType, LongDouble*, const IndexType )
{
    // LongDouble not supported by LAPACK, should give a serious runtime error
    COMMON_THROWEXCEPTION( "getinv<LongDouble> unsupported" )
}

/* ------------------------------------------------------------------------- */
/*      getinv<ComplexLongDouble>                                                   */
/* ------------------------------------------------------------------------- */

template<>
void LAPACKe_LAPACK::getinv( const IndexType, ComplexLongDouble*, const IndexType )
{
    // ComplexLongDouble not supported by LAPACK, should give a serious runtime error
    COMMON_THROWEXCEPTION( "getinv<ComplexLongDouble> unsupported" )
}

/* ------------------------------------------------------------------------- */
/*      getri<float>                                                         */
/* ------------------------------------------------------------------------- */

template<>
int LAPACKe_LAPACK::getri( const CBLAS_ORDER order, const int n, float* const a, const int lda, int* const ipiv )
{
    SCAI_LOG_INFO( logger, "getri<float> for A of size " << n << " x " << n )

    int matrix_order = lapack_order( order );

    int info = LAPACKE_sgetri( matrix_order, n, a, lda, ipiv );

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
/*      getri<double>                                                        */
/* ------------------------------------------------------------------------- */

template<>
int LAPACKe_LAPACK::getri( const CBLAS_ORDER order, const int n, double* const a, const int lda, int* const ipiv )
{
    SCAI_LOG_INFO( logger, "getri<double> for A of size " << n << " x " << n )

    int matrix_order = lapack_order( order );

    int info = LAPACKE_dgetri( matrix_order, n, a, lda, ipiv );

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
/*      getri<ComplexFloat>                                                         */
/* ------------------------------------------------------------------------- */

template<>
int LAPACKe_LAPACK::getri( const CBLAS_ORDER order, const int n, ComplexFloat* const a, const int lda, int* const ipiv )
{
    SCAI_LOG_INFO( logger, "getri<ComplexFloat> for A of size " << n << " x " << n )

    int matrix_order = lapack_order( order );

    int info = LAPACKE_cgetri( matrix_order, n, mklCast( a ), lda, ipiv );

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
/*      getri<ComplexDouble>                                                 */
/* ------------------------------------------------------------------------- */

template<>
int LAPACKe_LAPACK::getri(
    const CBLAS_ORDER order,
    const int n,
    ComplexDouble* const a,
    const int lda,
    int* const ipiv )
{
    SCAI_LOG_INFO( logger, "getri<ComplexDouble> for A of size " << n << " x " << n )

    int matrix_order = lapack_order( order );

    int info = LAPACKE_zgetri( matrix_order, n, mklCast( a ), lda, ipiv );

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
/*      getri<LongDouble>                                                    */
/* ------------------------------------------------------------------------- */

template<>
IndexType LAPACKe_LAPACK::getri( const CBLAS_ORDER, const int, LongDouble* const, const int, int* const )
{
    // LongDouble not supported by LAPACK, should give a serious runtime error
    COMMON_THROWEXCEPTION( "getri<LongDouble> unsupported" )
}

/* ------------------------------------------------------------------------- */
/*      getri<ComplexLongDouble>                                                    */
/* ------------------------------------------------------------------------- */

template<>
IndexType LAPACKe_LAPACK::getri( const CBLAS_ORDER, const int, ComplexLongDouble* const, const int, int* const )
{
    // ComplexLongDouble not supported by LAPACK, should give a serious runtime error
    COMMON_THROWEXCEPTION( "getri<ComplexLongDouble> unsupported" )
}

/* ------------------------------------------------------------------------- */
/*      tptrs<float>                                                         */
/* ------------------------------------------------------------------------- */

template<>
int LAPACKe_LAPACK::tptrs(
    const CBLAS_ORDER order,
    const CBLAS_UPLO uplo,
    const CBLAS_TRANSPOSE trans,
    const CBLAS_DIAG diag,
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

    SCAI_LOG_INFO( logger,
                   "tptrs<float>, n = " << n << ", nrhs = " << nrhs << ", order = " << matrix_order << ", UL = " << UL << ", TA = " << TA << ", DI = " << DI );

    SCAI_ASSERT_ERROR( ldb >= std::max( 1, n ), "ldb = " << ldb << " out of range" );

    int info = LAPACKE_stptrs( matrix_order, UL, TA, DI, n, nrhs, AP, B, ldb );

    return info;
}

/* ------------------------------------------------------------------------- */
/*      tptrs<double>                                                        */
/* ------------------------------------------------------------------------- */

template<>
int LAPACKe_LAPACK::tptrs(
    const CBLAS_ORDER order,
    const CBLAS_UPLO uplo,
    const CBLAS_TRANSPOSE trans,
    const CBLAS_DIAG diag,
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

    SCAI_LOG_INFO( logger,
                   "tptrs<double>, n = " << n << ", nrhs = " << nrhs << ", order = " << matrix_order << ", UL = " << UL << ", TA = " << TA << ", DI = " << DI );

    int info = LAPACKE_dtptrs( matrix_order, UL, TA, DI, n, nrhs, AP, B, ldb );

    return info;
}

/* ------------------------------------------------------------------------- */
/*      tptrs<ComplexFloat>                                                  */
/* ------------------------------------------------------------------------- */

template<>
int LAPACKe_LAPACK::tptrs(
    const CBLAS_ORDER order,
    const CBLAS_UPLO uplo,
    const CBLAS_TRANSPOSE trans,
    const CBLAS_DIAG diag,
    const int n,
    const int nrhs,
    const ComplexFloat* AP,
    ComplexFloat* B,
    const int ldb )
{
    char UL = BLASHelper::lapack_uplo( uplo );
    char TA = BLASHelper::lapack_transpose( trans );
    char DI = BLASHelper::lapack_diag( diag );

    int matrix_order = lapack_order( order );

    SCAI_LOG_INFO( logger,
                   "tptrs<ComplexFloat>, n = " << n << ", nrhs = " << nrhs << ", order = " << matrix_order << ", UL = " << UL << ", TA = " << TA << ", DI = " << DI );

    SCAI_ASSERT_ERROR( ldb >= std::max( 1, n ), "ldb = " << ldb << " out of range" );

    int info = LAPACKE_ctptrs( matrix_order, UL, TA, DI, n, nrhs, mklCast( AP ), mklCast( B ), ldb );

    return info;
}

/* ------------------------------------------------------------------------- */
/*      tptrs<ComplexDouble>                                                 */
/* ------------------------------------------------------------------------- */

template<>
int LAPACKe_LAPACK::tptrs(
    const CBLAS_ORDER order,
    const CBLAS_UPLO uplo,
    const CBLAS_TRANSPOSE trans,
    const CBLAS_DIAG diag,
    const int n,
    const int nrhs,
    const ComplexDouble* AP,
    ComplexDouble* B,
    const int ldb )
{
    char UL = BLASHelper::lapack_uplo( uplo );
    char TA = BLASHelper::lapack_transpose( trans );
    char DI = BLASHelper::lapack_diag( diag );

    int matrix_order = lapack_order( order );

    SCAI_LOG_INFO( logger,
                   "tptrs<ComplexDouble>, n = " << n << ", nrhs = " << nrhs << ", order = " << matrix_order << ", UL = " << UL << ", TA = " << TA << ", DI = " << DI );

    SCAI_ASSERT_ERROR( ldb >= std::max( 1, n ), "ldb = " << ldb << " out of range" );

    int info = LAPACKE_ztptrs( matrix_order, UL, TA, DI, n, nrhs, mklCast( AP ), mklCast( B ), ldb );

    return info;
}

/* ------------------------------------------------------------------------- */
/*      tptrs<LongDouble>                                                 */
/* ------------------------------------------------------------------------- */

template<>
int LAPACKe_LAPACK::tptrs(
    const CBLAS_ORDER,
    const CBLAS_UPLO,
    const CBLAS_TRANSPOSE,
    const CBLAS_DIAG,
    const int,
    const int,
    const LongDouble*,
    LongDouble*,
    const int )
{
    COMMON_THROWEXCEPTION( "tptrs<LongDouble> unsupported" )
}

/* ------------------------------------------------------------------------- */
/*      tptrs<ComplexLongDouble>                                                 */
/* ------------------------------------------------------------------------- */

template<>
int LAPACKe_LAPACK::tptrs(
    const CBLAS_ORDER,
    const CBLAS_UPLO,
    const CBLAS_TRANSPOSE,
    const CBLAS_DIAG,
    const int,
    const int,
    const ComplexLongDouble*,
    ComplexLongDouble*,
    const int )
{
    COMMON_THROWEXCEPTION( "tptrs<ComplexLongDouble> unsupported" )
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void LAPACKe_LAPACK::registerKernels( bool deleteFlag )
{
    using kregistry::KernelRegistry;
    using common::context::Host;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_REPLACE;   // priority over OpenMPBLAS

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

#define LAMA_LAPACK_REGISTER(z, I, _)                                      \
    KernelRegistry::set<BLASKernelTrait::getrf<ARITHMETIC_HOST_TYPE_##I> >( getrf, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::getri<ARITHMETIC_HOST_TYPE_##I> >( getri, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::getinv<ARITHMETIC_HOST_TYPE_##I> >( getinv, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::tptrs<ARITHMETIC_HOST_TYPE_##I> >( tptrs, Host, flag );    \

    BOOST_PP_REPEAT( ARITHMETIC_HOST_EXT_TYPE_CNT, LAMA_LAPACK_REGISTER, _ )

#undef LAMA_LAPACK_REGISTER
}

LAPACKe_LAPACK::LAPACKe_LAPACK()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

LAPACKe_LAPACK::~LAPACKe_LAPACK()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

LAPACKe_LAPACK LAPACKe_LAPACK::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
