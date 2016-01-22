/*
 * CUBLASWrapper.hpp
 *
 *  Created on: 24.08.2015
 *      Author: eschricker
 */

#pragma once

// internal scai libraries
#include <scai/blaskernel/cuda/CUBLASDefinitions.hpp>

#include <scai/common/cuda/CUDAError.hpp>

// CUDA
#include <cuComplex.h>
#include <cublas_v2.h>

namespace scai {

extern cublasHandle_t CUDAContext_cublasHandle;

namespace blaskernel {

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT CUBLASWrapper;

#define CUBLASWRAPPER_DEF( ValueType, CUBLASValueType, CUBLASAbsValueType, prefix1, prefix2, prefix3, DOT ) 													\
template<>																												\
class COMMON_DLL_IMPORTEXPORT CUBLASWrapper<ValueType>																	\
{																														\
public:																													\
	typedef CUBLASDefinitions::BLASIndexType BLASIndexType;																\
	typedef CUBLASDefinitions::BLASTrans BLASTrans;																		\
																														\
	static void scal( const BLASIndexType n, const ValueType alpha, ValueType* x, const BLASIndexType incX )			\
	{																													\
		SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( scal, prefix1 )( CUDAContext_cublasHandle, n,								\
				reinterpret_cast<const CUBLASValueType*>( &alpha ), reinterpret_cast<CUBLASValueType*>( x ), incX ),	\
				"CUBLASWrapper::scal<" #ValueType ">()" );																\
	}																													\
																														\
	static ValueType nrm2(const BLASIndexType n, const ValueType *x, const BLASIndexType incX) 							\
	{																													\
		ValueType nrm2 = 0;																								\
		SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( nrm2, prefix2 )(CUDAContext_cublasHandle,									\
				n, reinterpret_cast<const CUBLASValueType*>(x),	incX, reinterpret_cast<CUBLASAbsValueType*>(&nrm2)),	\
				"CUBLASWrapper::nrm2<" #ValueType ">");																	\
		return nrm2;																									\
	}																													\
																														\
	static ValueType asum(const BLASIndexType n, const ValueType *x,BLASIndexType incX)									\
	{																													\
		ValueType asum = 0;																								\
		SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( asum, prefix2 )(CUDAContext_cublasHandle, n, 								\
				reinterpret_cast<const CUBLASValueType*>(x), incX, reinterpret_cast<CUBLASAbsValueType*>( &asum )),		\
				"CUBLASWrapper::asum<" #ValueType ">");																	\
		return asum;																									\
	}																													\
																														\
	static BLASIndexType iamax(const BLASIndexType n, const ValueType *x, const BLASIndexType incX) 					\
	{																													\
		BLASIndexType iamax;																							\
		SCAI_CUBLAS_CALL(CUBLAS_BLAS_NAME( amax, prefix3 )(CUDAContext_cublasHandle,									\
				n, reinterpret_cast<const CUBLASValueType*>(x),	incX, &iamax),											\
				"CUBLASWrapper::iamax<" #ValueType ">");																\
		return iamax;																									\
	}																													\
																														\
	static void swap(const BLASIndexType n, ValueType *x, const BLASIndexType incX, ValueType *y, 						\
			const BLASIndexType incY)																					\
	{																													\
		SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( swap, prefix1 )(CUDAContext_cublasHandle, n,								\
				reinterpret_cast<CUBLASValueType*>(x), incX, reinterpret_cast<CUBLASValueType*>(y), incY ),				\
				"CUBLASWrapper::swap<" #ValueType ">");																	\
	}																													\
																														\
	static void copy(const BLASIndexType n, const ValueType *x, const BLASIndexType incX, ValueType *y,					\
			const BLASIndexType incY)																					\
	{																													\
		SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( copy, prefix1 )(CUDAContext_cublasHandle, n,								\
				reinterpret_cast<const CUBLASValueType*>(x), incX, reinterpret_cast<CUBLASValueType*>(y), incY ),		\
				"CUBLASWrapper::copy<" #ValueType ">");																	\
	}																													\
																														\
	static void axpy(const BLASIndexType n, const ValueType alpha,														\
			const ValueType *x, const BLASIndexType incX, ValueType *y, const BLASIndexType incY) 						\
	{																													\
		SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( axpy, prefix1 )(CUDAContext_cublasHandle, n,								\
				reinterpret_cast<const CUBLASValueType*>(&alpha), reinterpret_cast<const CUBLASValueType*>(x),			\
				incX, reinterpret_cast<CUBLASValueType*>(y), incY),														\
				"CUBLASWrapper::axpy<" #ValueType ">");																	\
	}																													\
																														\
	static ValueType dot(const BLASIndexType n, const ValueType *x,														\
			const BLASIndexType incX, const ValueType *y, const BLASIndexType incY)										\
	{																													\
		ValueType dot;																									\
		SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( DOT, prefix1 )(CUDAContext_cublasHandle, n,									\
				reinterpret_cast<const CUBLASValueType*>(x), incX, reinterpret_cast<const CUBLASValueType*>(y), incY,	\
				reinterpret_cast<CUBLASValueType*>(&dot)), 																\
				"CUBLASWrapper::dot<" #ValueType ">");																	\
		return dot;																										\
	}																													\
																														\
	static void gemv( const BLASTrans transA, const BLASIndexType m, const BLASIndexType n,								\
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* x,						\
			const BLASIndexType incX, const ValueType beta, ValueType* y, const BLASIndexType incY) 					\
	{																													\
		SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( gemv, prefix1 )(CUDAContext_cublasHandle, transA,							\
				m, n, reinterpret_cast<const CUBLASValueType*>(&alpha), reinterpret_cast<const CUBLASValueType*>(A),	\
				lda, reinterpret_cast<const CUBLASValueType*>(x), incX, 												\
				reinterpret_cast<const CUBLASValueType*>(&beta), reinterpret_cast<CUBLASValueType*>(y), incY),			\
				"CUBLASWrapper::gemv<" #ValueType ">");																	\
	}																													\
																														\
	static void gemm( const BLASTrans transA, const BLASTrans transB,													\
			const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,										\
			const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* B,						\
			const BLASIndexType ldb, const ValueType beta, ValueType* C, const BLASIndexType ldc) 						\
	{																													\
		SCAI_CUBLAS_CALL( CUBLAS_BLAS_NAME( gemm, prefix1 )(CUDAContext_cublasHandle, transA, transB,					\
				m, n, k, reinterpret_cast<const CUBLASValueType*>(&alpha), 												\
				reinterpret_cast<const CUBLASValueType*>(A), lda, reinterpret_cast<const CUBLASValueType*>(B), 			\
				ldb, reinterpret_cast<const CUBLASValueType*>(&beta),													\
				reinterpret_cast<CUBLASValueType*>(C), ldc),															\
				"CUBLASWrapper::gemm<" #ValueType ">");																	\
	}																													\
};

CUBLASWRAPPER_DEF( float, float, float, S, S, Is, dot )
CUBLASWRAPPER_DEF( double, double, double, D, D, Id, dot )

#ifdef SCAI_COMPLEX_SUPPORTED
CUBLASWRAPPER_DEF( ComplexFloat, cuFloatComplex, float, C, Sc, Ic, dotu )
CUBLASWRAPPER_DEF( ComplexDouble, cuDoubleComplex, double, Z, Dz, Iz, dotu )
#endif

#undef CUBLASWRAPPER_DEF

} /* end namespace blaskernel */

} /* end namespace scai */

