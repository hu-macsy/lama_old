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
 * Permission is hereby granted, free of BLASTransposege, to any person obtaining a copy
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

// local library
#include <scai/blaskernel/mic/cblas.hpp>
#include <scai/blaskernel/mic/MICBLASTrait.hpp>


// external
#include <mkl.h>

namespace scai {

namespace blaskernel {


#define MICBLASWRAPPER_DEF( ValueType, MICBLASValueType, prefix1, prefix2, prefix3, DOT ) 								\
template<>																												\
class COMMON_DLL_IMPORTEXPORT BLASWrapper<ValueType>																	\
{																														\
public:																													\
	typedef BLASTrait::BLASIndexType BLASIndexType;																		\
	typedef BLASTrait::BLASTrans BLASTrans;																				\
																														\
	static void scal( const BLASIndexType n, const ValueType alpha, ValueType* x, const BLASIndexType incX )			\
	{																													\
		MIC_BLAS_NAME( scal, prefix1 )( &n, reinterpret_cast<MICBLASValueType*>( &alpha ), reinterpret_cast<MICBLASValueType*>( x ), &incX );															\
	}																													\
																														\
	static ValueType nrm2(const BLASIndexType n, const ValueType *x, const BLASIndexType incX) 							\
	{																													\
		return MIC_BLAS_NAME( nrm2, prefix2 )( &n, reinterpret_cast<MICBLASValueType*>( x ), &incX ); 															\
	}																													\
																														\
	static ValueType asum(const BLASIndexType n, const ValueType *x,BLASIndexType incX)									\
	{																													\
		return MIC_BLAS_NAME( asum, prefix2 )( &n, reinterpret_cast<MICBLASValueType*>( x ), &incX );															\
	}																													\
																														\
	static BLASIndexType iamax(const BLASIndexType n, const ValueType *x, const BLASIndexType incX) 					\
	{																													\
		return MIC_BLAS_NAME( amax, prefix3 )( &n, reinterpret_cast<MICBLASValueType*>( x ), &incX );															\
	}																													\
																														\
	static void swap(const BLASIndexType n, ValueType *x, const BLASIndexType incX, ValueType *y, 						\
			const BLASIndexType incY)																					\
	{																													\
		MIC_BLAS_NAME( swap, prefix1 )( &n, reinterpret_cast<MICBLASValueType*>( x ), &incX, reinterpret_cast<MICBLASValueType*>( y ), &incY );															\
	}																													\
																														\
	static void copy(const BLASIndexType n, const ValueType *x, const BLASIndexType incX, ValueType *y,					\
			const BLASIndexType incY)																					\
	{																													\
		MIC_BLAS_NAME( copy, prefix1 )( &n, reinterpret_cast<MICBLASValueType*>( x ),	&incX, reinterpret_cast<MICBLASValueType*>( y ), &incY);														\
	}																													\
																														\
	static void axpy(const BLASIndexType n, const ValueType alpha,														\
			const ValueType *x, const BLASIndexType incX, ValueType *y, const BLASIndexType incY) 						\
	{																													\
		MIC_BLAS_NAME( axpy, prefix1 )( &n, reinterpret_cast<MICBLASValueType*>( &alpha ), reinterpret_cast<MICBLASValueType*>( x ), &incX, reinterpret_cast<MICBLASValueType*>( y ), &incY ); 												\
	}																													\
																														\
	static ValueType dot(const BLASIndexType n, const ValueType *x,														\
			const BLASIndexType incX, const ValueType *y, const BLASIndexType incY)										\
	{																													\
		return MIC_BLAS_NAME( DOT, prefix1 )(&n, reinterpret_cast<MICBLASValueType*>( x ), &incX, reinterpret_cast<MICBLASValueType*>( y ), &incY);													\
	}																													\
																														\
	static void gemv( const BLASTrans transA, const BLASIndexType m, const BLASIndexType n,								\
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* x,						\
			const BLASIndexType incX, const ValueType beta, ValueType* y, const BLASIndexType incY) 					\
	{																													\
		MIC_BLAS_NAME( gemv, prefix1 )( &transA, &m, &n, reinterpret_cast<MICBLASValueType*>( &alpha ), reinterpret_cast<MICBLASValueType*>( A ), &lda, reinterpret_cast<MICBLASValueType*>( x ), &incX, reinterpret_cast<MICBLASValueType*>( &beta ), reinterpret_cast<MICBLASValueType*>( y ), &incY );					\
	}																													\
																														\
	static void gemm( const BLASTrans transA, const BLASTrans transB,													\
			const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,										\
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* B,						\
			const BLASIndexType ldb, const ValueType beta, ValueType* C, const BLASIndexType ldc) 						\
	{																													\
		MIC_BLAS_NAME( gemm, prefix1 )( &transA, &transB, &m, &n, &k, reinterpret_cast<MICBLASValueType*>( &alpha ), reinterpret_cast<MICBLASValueType*>( A ), &lda, reinterpret_cast<MICBLASValueType*>( B ), &ldb, reinterpret_cast<MICBLASValueType*>( &beta ), reinterpret_cast<MICBLASValueType*>( C ), &ldc );		\
	}																													\
};

MICBLASWRAPPER_DEF( float, float, s, s, is, dot );
MICBLASWRAPPER_DEF( double, float, d, d, id, dot );

#ifdef SCAI_COMPLEX_SUPPORTED
MICBLASWRAPPER_DEF( ComplexFloat, MKL_Complex8, c, sc, ic, dotc );
MICBLASWRAPPER_DEF( ComplexDouble, MKL_Complex16, z, dz, iz, dotc );
#endif

#undef MICBLASWRAPPER_DEF

} /* end namespace blaskernel */

} /* end namespace scai */

