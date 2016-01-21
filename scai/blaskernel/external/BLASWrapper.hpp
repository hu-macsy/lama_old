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
 * @date 10.01.2016
 * @since 2.0.0
 */

#pragma once

// local library
#include <scai/blaskernel/external/BLASDefinitions.hpp>
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/common/macros/unused.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>

namespace scai {

namespace blaskernel {

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT BLASWrapper;

#define BLASWRAPPER_DEF( ValueType, prefix1, prefix2, prefix3, DOT ) 													\
template<>																												\
class COMMON_DLL_IMPORTEXPORT BLASWrapper<ValueType>																	\
{																														\
public:																													\
	typedef BLASDefinitions::BLASIndexType BLASIndexType;																\
	typedef BLASDefinitions::BLASTrans BLASTrans;																		\
																														\
	static void scal( const BLASIndexType n, const ValueType alpha, ValueType* x, const BLASIndexType incX )			\
	{																													\
		FORTRAN_BLAS_NAME( scal, prefix1 )( &n, &alpha, x, &incX );														\
	}																													\
																														\
	static ValueType nrm2(const BLASIndexType n, const ValueType *x, const BLASIndexType incX) 							\
	{																													\
		return FORTRAN_BLAS_NAME( nrm2, prefix2 )( &n, x, &incX ); 														\
	}																													\
																														\
	static ValueType asum(const BLASIndexType n, const ValueType *x,BLASIndexType incX)									\
	{																													\
		return FORTRAN_BLAS_NAME( asum, prefix2 )(&n, x, &incX);														\
	}																													\
																														\
	static BLASIndexType iamax(const BLASIndexType n, const ValueType *x, const BLASIndexType incX) 					\
	{																													\
		return FORTRAN_BLAS_NAME( amax, prefix3 )(&n, x, &incX);														\
	}																													\
																														\
	static void swap(const BLASIndexType n, ValueType *x, const BLASIndexType incX, ValueType *y, 						\
			const BLASIndexType incY)																					\
	{																													\
		FORTRAN_BLAS_NAME( swap, prefix1 )(&n, x, &incX, y, &incY);														\
	}																													\
																														\
	static void copy(const BLASIndexType n, const ValueType *x, const BLASIndexType incX, ValueType *y,					\
			const BLASIndexType incY)																					\
	{																													\
		FORTRAN_BLAS_NAME( copy, prefix1 )(&n, x,	&incX, y, &incY);													\
	}																													\
																														\
	static void axpy(const BLASIndexType n, const ValueType alpha,														\
			const ValueType *x, const BLASIndexType incX, ValueType *y, const BLASIndexType incY) 						\
	{																													\
		FORTRAN_BLAS_NAME( axpy, prefix1 )(&n, &alpha, x,	&incX, y, &incY); 											\
	}																													\
																														\
	static ValueType dot(const BLASIndexType n, const ValueType *x,														\
			const BLASIndexType incX, const ValueType *y, const BLASIndexType incY)										\
	{																													\
		return FORTRAN_BLAS_NAME( DOT, prefix1 )(&n, x, &incX, y, &incY);												\
	}																													\
																														\
	static void gemv( const BLASTrans transA, const BLASIndexType m, const BLASIndexType n,								\
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* x,						\
			const BLASIndexType incX, const ValueType beta, ValueType* y, const BLASIndexType incY) 					\
	{																													\
		FORTRAN_BLAS_NAME( gemv, prefix1 )( &transA, &m, &n, &alpha, A, &lda, x, &incX, &beta, y, &incY );				\
	}																													\
																														\
	static void gemm( const BLASTrans transA, const BLASTrans transB,													\
			const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,										\
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* B,						\
			const BLASIndexType ldb, const ValueType beta, ValueType* C, const BLASIndexType ldc) 						\
	{																													\
		FORTRAN_BLAS_NAME( gemm, prefix1 )( &transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );	\
	}																													\
};

BLASWRAPPER_DEF( float, s, s, is, dot );
BLASWRAPPER_DEF( double, d, d, id, dot );

#ifdef SCAI_COMPLEX_SUPPORTED
BLASWRAPPER_DEF( ComplexFloat, c, sc, ic, dotu );
BLASWRAPPER_DEF( ComplexDouble, z, dz, iz, dotu );
#endif

#undef BLASWRAPPER_DEF

} /* end namespace blaskernel */

} /* end namespace scai */

