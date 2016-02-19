/**
 * @file MICBLASWrapper.hpp
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
 * @date 25.01.2016
 * @since 2.0.0
 */

#pragma once

// local library
#include <scai/blaskernel/cblas.hpp>
#include <scai/blaskernel/mic/MICBLASTrait.hpp>

// scai internal libraries
#include <scai/common/mic/MICCallable.hpp>

// external
#include <mkl_blas.h>

namespace scai {

namespace blaskernel {

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT MICBLASWrapper;

#define MICBLASWRAPPER_DEF( ValueType, MICBLASValueType, prefix1, prefix2, prefix3, DOT ) 								\
template<>																												\
class COMMON_DLL_IMPORTEXPORT MICBLASWrapper<ValueType>																	\
{																														\
public:																													\
	typedef MICBLASTrait::BLASIndexType BLASIndexType;																	\
	typedef MICBLASTrait::BLASTrans BLASTrans;																			\
																														\
	static MIC_CALLABLE_MEMBER void gemv( const BLASTrans transA, const BLASIndexType m, const BLASIndexType n,			\
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* x,						\
			const BLASIndexType incX, const ValueType beta, ValueType* y, const BLASIndexType incY) 					\
	{																													\
		MIC_BLAS_NAME( gemv, prefix1 )( &transA, &m, &n, reinterpret_cast<const MICBLASValueType*>( &alpha ),			\
				reinterpret_cast<const MICBLASValueType*>( A ), &lda, reinterpret_cast<const MICBLASValueType*>( x ), 	\
				&incX, reinterpret_cast<const MICBLASValueType*>( &beta ), reinterpret_cast<MICBLASValueType*>( y ), 	\
				&incY );																								\
	}																													\
																														\
	static MIC_CALLABLE_MEMBER void gemm( const BLASTrans transA, const BLASTrans transB,								\
			const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,										\
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* B,						\
			const BLASIndexType ldb, const ValueType beta, ValueType* C, const BLASIndexType ldc) 						\
	{																													\
		MIC_BLAS_NAME( gemm, prefix1 )( &transA, &transB, &m, &n, &k, 													\
				reinterpret_cast<const MICBLASValueType*>( &alpha ), reinterpret_cast<const MICBLASValueType*>( A ), 	\
				&lda, reinterpret_cast<const MICBLASValueType*>( B ), &ldb, 											\
				reinterpret_cast<const MICBLASValueType*>( &beta ), reinterpret_cast<MICBLASValueType*>( C ), &ldc );	\
	}																													\
};

MICBLASWRAPPER_DEF( float, float, s, s, is, dot );
MICBLASWRAPPER_DEF( double, double, d, d, id, dot );

#ifdef SCAI_COMPLEX_SUPPORTED
MICBLASWRAPPER_DEF( ComplexFloat, MKL_Complex8, c, sc, ic, dotc );
MICBLASWRAPPER_DEF( ComplexDouble, MKL_Complex16, z, dz, iz, dotc );
#endif

#undef MICBLASWRAPPER_DEF

} /* end namespace blaskernel */

} /* end namespace scai */

