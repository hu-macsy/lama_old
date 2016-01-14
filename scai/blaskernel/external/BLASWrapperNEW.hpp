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

#define FORTRAN_BLAS_NAME( name, prefix ) prefix##name##_

#define FORTRAN_BLAS_DEF( name, prefix, retType, definition ) 			\
        retType FORTRAN_BLAS_NAME( name, prefix )( definition );

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
#define CALL_DEF_SWAP( type ) const _BLASWrapper::BLASIndexType *, type*, const _BLASWrapper::BLASIndexType *, type*, const _BLASWrapper::BLASIndexType *

FORTRAN_BLAS_DEF( swap, s, void, CALL_DEF_SWAP( float ) )
FORTRAN_BLAS_DEF( swap, d, void, CALL_DEF_SWAP( double ) )
FORTRAN_BLAS_DEF( swap, c, void, CALL_DEF_SWAP( ComplexFloat ) )
FORTRAN_BLAS_DEF( swap, z, void, CALL_DEF_SWAP( ComplexDouble ) )

#undef CALL_DEF_SWAP
#define CALL_DEF_COPY( type ) const _BLASWrapper::BLASIndexType *, const type*, const _BLASWrapper::BLASIndexType *, type*, const _BLASWrapper::BLASIndexType *

FORTRAN_BLAS_DEF( copy, s, void, CALL_DEF_COPY( float ) );
FORTRAN_BLAS_DEF( copy, d, void, CALL_DEF_COPY( double ) );
FORTRAN_BLAS_DEF( copy, c, void, CALL_DEF_COPY( ComplexFloat ) );
FORTRAN_BLAS_DEF( copy, z, void, CALL_DEF_COPY( ComplexDouble ) );

#undef CALL_DEF_COPY
#define CALL_DEF_AXPY( type ) const _BLASWrapper::BLASIndexType *, const type*, const type*, const _BLASWrapper::BLASIndexType *, type*, const _BLASWrapper::BLASIndexType *

FORTRAN_BLAS_DEF( axpy, s, void, CALL_DEF_AXPY( float ) );
FORTRAN_BLAS_DEF( axpy, d, void, CALL_DEF_AXPY( double ) );
FORTRAN_BLAS_DEF( axpy, c, void, CALL_DEF_AXPY( ComplexFloat ) );
FORTRAN_BLAS_DEF( axpy, z, void, CALL_DEF_AXPY( ComplexDouble ) );

#undef CALL_DEF_AXPY
#define CALL_DEF_DOT( type ) const _BLASWrapper::BLASIndexType *, const type*, const _BLASWrapper::BLASIndexType *, const type*, const _BLASWrapper::BLASIndexType *

FORTRAN_BLAS_DEF( dot, s, float, CALL_DEF_DOT( float ) );
FORTRAN_BLAS_DEF( dot, d, double, CALL_DEF_DOT( double ) );
FORTRAN_BLAS_DEF( dotu, c, ComplexFloat, CALL_DEF_DOT( ComplexFloat ) );
FORTRAN_BLAS_DEF( dotu, z, ComplexDouble, CALL_DEF_DOT( ComplexDouble ) );

#undef CALL_DEF_DOT
#define CALL_DEF_SCAL( type ) const _BLASWrapper::BLASIndexType *, const type*, type*, const _BLASWrapper::BLASIndexType *

FORTRAN_BLAS_DEF( scal, s, void, CALL_DEF_SCAL( float ) );
FORTRAN_BLAS_DEF( scal, d, void, CALL_DEF_SCAL( double ) );
FORTRAN_BLAS_DEF( scal, c, void, CALL_DEF_SCAL( ComplexFloat ) );
FORTRAN_BLAS_DEF( scal, z, void, CALL_DEF_SCAL( ComplexDouble ) );

#undef CALL_DEF_SCAL
#define CALL_DEF_NRM2( type ) const _BLASWrapper::BLASIndexType *, const type*, const _BLASWrapper::BLASIndexType *

FORTRAN_BLAS_DEF( nrm2, s, float, CALL_DEF_NRM2( float ) );
FORTRAN_BLAS_DEF( nrm2, d, double, CALL_DEF_NRM2( double ) );
FORTRAN_BLAS_DEF( nrm2, sc, float, CALL_DEF_NRM2( ComplexFloat ) );
FORTRAN_BLAS_DEF( nrm2, dz, double, CALL_DEF_NRM2( ComplexDouble ) );

#undef CALL_DEF_NRM2
#define CALL_DEF_ASUM( type ) const _BLASWrapper::BLASIndexType *, const type*, const _BLASWrapper::BLASIndexType *

FORTRAN_BLAS_DEF( asum, s, float, CALL_DEF_ASUM( float ) );
FORTRAN_BLAS_DEF( asum, d, double, CALL_DEF_ASUM( double ) );
FORTRAN_BLAS_DEF( asum, sc, float, CALL_DEF_ASUM( ComplexFloat ) );
FORTRAN_BLAS_DEF( asum, dz, double, CALL_DEF_ASUM( ComplexDouble ) );

#undef CALL_DEF_ASUM
#define CALL_DEF_AMAX( type ) const _BLASWrapper::BLASIndexType *, const type*, const _BLASWrapper::BLASIndexType *

FORTRAN_BLAS_DEF( amax, is, _BLASWrapper::BLASIndexType, CALL_DEF_AMAX( float ) );
FORTRAN_BLAS_DEF( amax, id, _BLASWrapper::BLASIndexType, CALL_DEF_AMAX( double ) );
FORTRAN_BLAS_DEF( amax, ic, _BLASWrapper::BLASIndexType, CALL_DEF_AMAX( ComplexFloat ) );
FORTRAN_BLAS_DEF( amax, iz, _BLASWrapper::BLASIndexType, CALL_DEF_AMAX( ComplexDouble ) );

#undef CALL_DEF_AMAX
#define CALL_DEF_GEMV( type ) const _BLASWrapper::BLASTrans*, const _BLASWrapper::BLASIndexType *, const _BLASWrapper::BLASIndexType *, const type*, const type*, const _BLASWrapper::BLASIndexType *, const type*, const _BLASWrapper::BLASIndexType *, const type*, type*, const _BLASWrapper::BLASIndexType *

FORTRAN_BLAS_DEF( gemv, s, void, CALL_DEF_GEMV( float ) );
FORTRAN_BLAS_DEF( gemv, d, void, CALL_DEF_GEMV( double ) );
FORTRAN_BLAS_DEF( gemv, c, void, CALL_DEF_GEMV( ComplexFloat ) );
FORTRAN_BLAS_DEF( gemv, z, void, CALL_DEF_GEMV( ComplexDouble ) );

#undef CALL_DEF_GEMM
#define CALL_DEF_GEMM( type ) const _BLASWrapper::BLASTrans*, const _BLASWrapper::BLASTrans*, const _BLASWrapper::BLASIndexType *, const _BLASWrapper::BLASIndexType *, const _BLASWrapper::BLASIndexType *, const type*, const type*, const _BLASWrapper::BLASIndexType *, const type*, const _BLASWrapper::BLASIndexType *, const type*, type*, const _BLASWrapper::BLASIndexType *

FORTRAN_BLAS_DEF( gemm, s, void, CALL_DEF_GEMM( float ) );
FORTRAN_BLAS_DEF( gemm, d, void, CALL_DEF_GEMM( double ) );
FORTRAN_BLAS_DEF( gemm, c, void, CALL_DEF_GEMM( ComplexFloat ) );
FORTRAN_BLAS_DEF( gemm, z, void, CALL_DEF_GEMM( ComplexDouble ) );

#undef CALL_DEF_GEMM

} /* end extern "C" */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT BLASWrapper: public _BLASWrapper
{
public:
	typedef _BLASWrapper::BLASIndexType BLASIndexType;
	typedef _BLASWrapper::BLASTrans BLASTrans;
};

#define BLASWRAPPER_DEF( ValueType, prefix1, prefix2, prefix3, DOT ) 											\
template<>																										\
class COMMON_DLL_IMPORTEXPORT BLASWrapper<ValueType>															\
{																												\
public:																											\
	typedef _BLASWrapper::BLASIndexType BLASIndexType;															\
	typedef _BLASWrapper::BLASTrans BLASTrans;																	\
																												\
	static void scal( const BLASIndexType n, const ValueType alpha, ValueType* x, const BLASIndexType incX )	\
	{																											\
		FORTRAN_BLAS_NAME( scal, prefix1 )( &n, &alpha, x, &incX );												\
	}																											\
																												\
	static ValueType nrm2(const BLASIndexType n, const ValueType *x,											\
				const BLASIndexType incX) {																		\
		return FORTRAN_BLAS_NAME( nrm2, prefix2 )( &n, x, &incX ); 											\
	}																											\
																												\
	static ValueType asum(const BLASIndexType n, const ValueType *x,BLASIndexType incX)							\
	{																											\
		return FORTRAN_BLAS_NAME( asum, prefix2 )(&n, x, &incX);											\
	}																											\
\
	static BLASIndexType iamax(const BLASIndexType n, const ValueType *x,\
			const BLASIndexType incX) {\
		return FORTRAN_BLAS_NAME( amax, prefix3 )(&n, x, &incX);\
	}\
\
	static void swap(const BLASIndexType n, ValueType *x,\
			const BLASIndexType incX, ValueType *y, const BLASIndexType incY) {\
		FORTRAN_BLAS_NAME( swap, prefix1 )(&n, x, &incX, y, &incY);\
	}\
\
	static void copy(const BLASIndexType n, const ValueType *x,\
			const BLASIndexType incX, ValueType *y, const BLASIndexType incY) {\
		FORTRAN_BLAS_NAME( copy, prefix1 )(&n, x,	&incX, y, &incY);\
	}\
\
	static void axpy(const BLASIndexType n, const ValueType alpha,\
			const ValueType *x, const BLASIndexType incX, ValueType *y, const BLASIndexType incY) {\
		FORTRAN_BLAS_NAME( axpy, prefix1 )(&n, &alpha, x,	&incX, y, &incY); \
	}\
\
	static ValueType dot(const BLASIndexType n, const ValueType *x,\
			const BLASIndexType incX, const ValueType *y, const BLASIndexType incY) {\
		return FORTRAN_BLAS_NAME( DOT, prefix1 )(&n, x, &incX, y, &incY);\
	}\
\
	static void gemv( const BLASTrans transA, const BLASIndexType m, const BLASIndexType n,\
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* x,\
			const BLASIndexType incX, const ValueType beta, ValueType* y, const BLASIndexType incY) {\
		FORTRAN_BLAS_NAME( gemv, prefix1 )( &transA, &m, &n, &alpha, A, &lda, x, &incX, &beta, y, &incY );\
	}\
\
	static void gemm( const BLASTrans transA, const BLASTrans transB,\
			const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,\
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* B,\
			const BLASIndexType ldb, const ValueType beta, ValueType* C, const BLASIndexType ldc) {\
		FORTRAN_BLAS_NAME( gemm, prefix1 )( &transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );\
	}\
};

BLASWRAPPER_DEF( float, s, s, is, dot );
BLASWRAPPER_DEF( double, d, d, id, dot );
BLASWRAPPER_DEF( ComplexFloat, c, sc, ic, dotu );
BLASWRAPPER_DEF( ComplexDouble, z, dz, iz, dotu );

} /* end namespace blaskernel */

} /* end namespace scai */

