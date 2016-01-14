/*
 * BLASWrapper.hpp
 *
 *  Created on: 24.08.2015
 *      Author: eschricker
 */

/**
 * @file BLASWrapper.hpp
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
 * @brief Wrapper for BLAS functions
 * @author Eric Schricker
 * @date 24.08.2015
 * @since 2.0.0
 */

#pragma once

#include <scai/common/exception/NotSupportedValueTypeException.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/blaskernel/cblas.hpp>

namespace scai {

namespace blaskernel {

#define BLAS_PTR_DEC( returnType, name, ... ) 			\
	private: 											\
		typedef returnType (*name##_f)( __VA_ARGS__ );	\
		static name##_f name##_ptr;


#define BLAS_PTR_CRE( name )	\
	template<typename ValueType> \
	typename BLASWrapper<ValueType>::name##_f BLASWrapper<ValueType>::name##_ptr;

#define BLAS_PTR_SET( name, ValueType, pointer ) 									\
	name##_ptr = &pointer;



class COMMON_DLL_IMPORTEXPORT _BLASWrapper {
public:


#ifdef F77_INT
	typedef F77_INT BLASIndexType;
#else
	typedef int BLASIndexType;
#endif

#ifdef  F77_CHAR
    typedef F77_CHAR BLASTrans
#else
	typedef char BLASTrans;
#endif

};

extern "C"
{
void dasum_( const _BLASWrapper::BLASIndexType *, const double*, const _BLASWrapper::BLASIndexType *, double* );
}

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT BLASWrapper: public _BLASWrapper {
public:
	typedef _BLASWrapper::BLASIndexType BLASIndexType;
	typedef _BLASWrapper::BLASTrans BLASTrans;

	static void scal( const BLASIndexType n, const ValueType alpha, ValueType* x, const BLASIndexType incX )
	{
		SCAI_ASSERT_UNEQUAL( scal_ptr, 0, "scal_ptr not set for " << common::TypeTraits<ValueType>::id() )

		scal_ptr( &n, &alpha, x, &incX );
	}

	static ValueType nrm2(const BLASIndexType n, const ValueType *x,
				const BLASIndexType incX) {
		SCAI_ASSERT_UNEQUAL( nrm2_ptr, 0, "nrm2_ptr not set for " << common::TypeTraits<ValueType>::id() )

		ValueType nrm2_val;

		nrm2_ptr( &n, x, &incX, &nrm2_val );

		return nrm2_val;
	}

	static ValueType asum(const BLASIndexType n, const ValueType *x,
			BLASIndexType incX) {
		SCAI_ASSERT_UNEQUAL( asum_ptr, 0, "asum_ptr not set for " << common::TypeTraits<ValueType>::id() )

		std::cout << "asum ptr --> " << asum_ptr << std::endl;
		ValueType asum_val = -18;

		asum_ptr(&n, x, &incX, &asum_val);

		std::cout << "BLASWrapper::asum<" << common::TypeTraits<ValueType>::id() << "> --> " << asum_val << std::endl;

		return asum_val;
	}

	static BLASIndexType iamax(const BLASIndexType n, const ValueType *x,
			const BLASIndexType incX) {
		SCAI_ASSERT_UNEQUAL( iamax_ptr, 0, "iamax_ptr not set for " << common::TypeTraits<ValueType>::id() )

		BLASIndexType iamax_val = 0;

		iamax_ptr(&n, x, &incX, &iamax_val);

		return iamax_val;
	}

	static void swap(const BLASIndexType n, ValueType *x,
			const BLASIndexType incX, ValueType *y, const BLASIndexType incY) {
		SCAI_ASSERT_UNEQUAL( swap_ptr, 0, "swap_ptr not set for " << common::TypeTraits<ValueType>::id() )

		swap_ptr(&n, x, &incX, y, &incY);
	}

	static void copy(const BLASIndexType n, const ValueType *x,
			const BLASIndexType incX, ValueType *y, const BLASIndexType incY) {
		SCAI_ASSERT_UNEQUAL( copy_ptr, 0, "copy_ptr not set for " << common::TypeTraits<ValueType>::id() )

		copy_ptr(&n, x,	&incX, y, &incY);
	}

	static void axpy(const BLASIndexType n, const ValueType alpha,
			const ValueType *x, const BLASIndexType incX, ValueType *y, const BLASIndexType incY) {
		SCAI_ASSERT_UNEQUAL( axpy_ptr, 0, "axpy_ptr not set for " << common::TypeTraits<ValueType>::id() )

		axpy_ptr(&n, &alpha, x,	&incX, y, &incY);
	}

	static ValueType dot(const BLASIndexType n, const ValueType *x,
			const BLASIndexType incX, const ValueType *y, const BLASIndexType incY) {
		SCAI_ASSERT_UNEQUAL( dot_ptr, 0, "dot_ptr not set for " << common::TypeTraits<ValueType>::id() )

		ValueType dot_val = 0;

		dot_ptr(&n, x, &incX, y, &incY, &dot_val);

		return dot_val;
	}

	static void gemv( const BLASTrans transA, const BLASIndexType m, const BLASIndexType n,
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* x,
			const BLASIndexType incX, const ValueType beta, ValueType* y, const BLASIndexType incY) {
		SCAI_ASSERT_UNEQUAL( gemv_ptr, 0, "gemv_ptr not set for " << common::TypeTraits<ValueType>::id() )

		gemv_ptr( &transA, &m, &n, &alpha, A, &lda, x, &incX, &beta, y, &incY );

	}

	static void gemm( const BLASTrans transA, const BLASTrans transB,
			const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* B,
			const BLASIndexType ldb, const ValueType beta, ValueType* C, const BLASIndexType ldc) {
		SCAI_ASSERT_UNEQUAL( gemm_ptr, 0, "gemm_ptr not set for " << common::TypeTraits<ValueType>::id() )

		gemm_ptr( &transA, &transB, &m, &n, &k,	&alpha, A, &lda, B, &ldb, &beta, C, &ldc );
	}
private:

	BLAS_PTR_DEC( void, scal, const BLASIndexType* n, const ValueType* alpha, ValueType* x, const BLASIndexType* incX );
	BLAS_PTR_DEC( void, nrm2, const BLASIndexType* n, const ValueType* x, const BLASIndexType* incX, ValueType* nrm2 );
	BLAS_PTR_DEC( void, asum, const BLASIndexType* n, const ValueType *x, const BLASIndexType* incX, ValueType* asum );
	BLAS_PTR_DEC( void, iamax, const BLASIndexType* n, const ValueType *x,	const BLASIndexType* incX, BLASIndexType* iamax );
	BLAS_PTR_DEC( void, swap, const BLASIndexType* n, ValueType *x,	const BLASIndexType* incX, ValueType *y, const BLASIndexType* incY );
	BLAS_PTR_DEC( void, copy, const BLASIndexType* n, const ValueType *x, const BLASIndexType* incX, ValueType *y, const BLASIndexType* incY );
	BLAS_PTR_DEC( void, axpy, const BLASIndexType* n, const ValueType *alpha, const ValueType *x, const BLASIndexType* incX, ValueType *y, const BLASIndexType* incY );
	BLAS_PTR_DEC( void, dot, const BLASIndexType* n, const ValueType *x, const BLASIndexType* incX, const ValueType *y, const BLASIndexType* incY, ValueType* dot );
	BLAS_PTR_DEC( void, gemv, const BLASTrans* transA, const BLASIndexType* m, const BLASIndexType* n, const ValueType* alpha, const ValueType* A, const BLASIndexType* lda, const ValueType* x, const BLASIndexType* incX, const ValueType* beta, ValueType* y, const BLASIndexType* incY );
	BLAS_PTR_DEC( void, gemm, const BLASTrans* transA, const BLASTrans* transB,	const BLASIndexType* m, const BLASIndexType* n, const BLASIndexType* k, const ValueType* alpha, const ValueType* A, const BLASIndexType* lda, const ValueType* B, const BLASIndexType* ldb, const ValueType* beta, ValueType* C, const BLASIndexType* ldc );

	static bool init;
	static bool initialize();
};

} /* end namespace blaskernel */

} /* end namespace scai */

