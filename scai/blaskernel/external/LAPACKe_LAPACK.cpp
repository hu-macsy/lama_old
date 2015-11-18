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
#include <scai/blaskernel/external/LAPACKeWrapper.hpp>
#include <scai/blaskernel/openmp/OpenMPBLAS1.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/ScalarType.hpp>

// external
#include <mkl_lapacke.h>

// boost
#include <boost/preprocessor.hpp>

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

template<typename ValueType>
IndexType LAPACKe_LAPACK::getrf(
    const CBLAS_ORDER order,
    const IndexType m,
    const IndexType n,
    ValueType* const A,
    const IndexType lda,
    IndexType* const ipiv )
{
    SCAI_LOG_INFO( logger, "getrf<float> for A of size " << m << " x " << n )

	if( common::getScalarType<IndexType>() != common::getScalarType<LAPACKeWrapper::LAPACKIndexType>() )
	{
		// ToDo: convert ipiv array
		COMMON_THROWEXCEPTION( "indextype mismatch" );
	}

    int info = LAPACKeWrapper::getrf<ValueType>( lapack_order( order ), static_cast<LAPACKeWrapper::LAPACKIndexType>(m), static_cast<LAPACKeWrapper::LAPACKIndexType>(n), A, static_cast<LAPACKeWrapper::LAPACKIndexType>(lda), ipiv );

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
void LAPACKe_LAPACK::getinv( const IndexType n, ValueType* a, const IndexType lda )
{
    int info = 0;

    // scoped_array, will also be freed in case of exception

    scoped_array<LAPACKeWrapper::LAPACKIndexType> ipiv( new LAPACKeWrapper::LAPACKIndexType[n] );

    SCAI_LOG_INFO( logger, "getinv<float> for " << n << " x " << n << " matrix, uses MKL" )

    info = LAPACKeWrapper::getrf( LAPACK_COL_MAJOR, static_cast<LAPACKeWrapper::LAPACKIndexType>(n), static_cast<LAPACKeWrapper::LAPACKIndexType>(n), a, static_cast<LAPACKeWrapper::LAPACKIndexType>(lda), ipiv.get() );

    // return error if factorization did not work

    if( info )
    {
        COMMON_THROWEXCEPTION( "MKL sgetrf failed, info = " << info )
    }

    info = LAPACKeWrapper::getri( LAPACK_COL_MAJOR, static_cast<LAPACKeWrapper::LAPACKIndexType>(n), a, static_cast<LAPACKeWrapper::LAPACKIndexType>(lda), ipiv.get() );

    if( info )
    {
        COMMON_THROWEXCEPTION( "MKL sgetri failed, info = " << info )
    }
}

/* ------------------------------------------------------------------------- */
/*      getri<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
int LAPACKe_LAPACK::getri( const CBLAS_ORDER order, const IndexType n, ValueType* const a, const IndexType lda, IndexType* const ipiv )
{
    SCAI_LOG_INFO( logger, "getri<float> for A of size " << n << " x " << n )

	if( common::getScalarType<IndexType>() != common::getScalarType<LAPACKeWrapper::LAPACKIndexType >() )
	{
		// ToDo: convert ipiv array
		COMMON_THROWEXCEPTION( "indextype mismatch" );
	}

    LAPACKeWrapper::LAPACKOrder matrix_order = lapack_order( order );

    int info = LAPACKeWrapper::getri( matrix_order, static_cast<LAPACKeWrapper::LAPACKIndexType>(n), a, static_cast<LAPACKeWrapper::LAPACKIndexType>(lda), ipiv );

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
/*      tptrs<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
int LAPACKe_LAPACK::tptrs(
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
    LAPACKeWrapper::LAPACKFlag UL = BLASHelper::lapack_uplo( uplo );
    LAPACKeWrapper::LAPACKFlag TA = BLASHelper::lapack_transpose( trans );
    LAPACKeWrapper::LAPACKFlag DI = BLASHelper::lapack_diag( diag );

    LAPACKeWrapper::LAPACKOrder matrix_order = lapack_order( order );

	if( common::getScalarType<IndexType>() != common::getScalarType<LAPACKeWrapper::LAPACKIndexType>() )
	{
		// ToDo: convert ipiv array
		COMMON_THROWEXCEPTION( "indextype mismatch" );
	}

    SCAI_LOG_INFO( logger,
                   "tptrs<float>, n = " << n << ", nrhs = " << nrhs << ", order = " << matrix_order << ", UL = " << UL << ", TA = " << TA << ", DI = " << DI );

    SCAI_ASSERT_ERROR( ldb >= std::max( 1, n ), "ldb = " << ldb << " out of range" );

    int info = LAPACKeWrapper::tptrs( matrix_order, UL, TA, DI, static_cast<LAPACKeWrapper::LAPACKIndexType>(n), static_cast<LAPACKeWrapper::LAPACKIndexType>(nrhs), AP, B, static_cast<LAPACKeWrapper::LAPACKIndexType>(ldb) );

    return info;
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

#define LAMA_LAPACKE_REGISTER(z, I, _)                                      \
    KernelRegistry::set<BLASKernelTrait::getrf<ARITHMETIC_HOST_TYPE_##I> >( getrf, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::getri<ARITHMETIC_HOST_TYPE_##I> >( getri, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::getinv<ARITHMETIC_HOST_TYPE_##I> >( getinv, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::tptrs<ARITHMETIC_HOST_TYPE_##I> >( tptrs, Host, flag );    \

    BOOST_PP_REPEAT( ARITHMETIC_HOST_EXT_TYPE_CNT, LAMA_LAPACKE_REGISTER, _ )

#undef LAMA_LAPACKE_REGISTER
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
