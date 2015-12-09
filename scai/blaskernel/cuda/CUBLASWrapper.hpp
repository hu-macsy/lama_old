/*
 * CUBLASWrapper.hpp
 *
 *  Created on: 24.08.2015
 *      Author: eschricker
 */

#pragma once

// internal scai libraries
#include <scai/common/exception/NotSupportedValueTypeException.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/SCAITypes.hpp>

// CUDA
#include <cublas_v2.h>
#include <cuComplex.h>

namespace scai {

extern cublasHandle_t CUDAContext_cublasHandle;

namespace blaskernel {

//cublasHandle_t CUDAContext_cublasHandle;

class COMMON_DLL_IMPORTEXPORT CUBLASWrapper {
public:
	typedef int BLASIndexType;

	// ------------- BLAS1 -------------
	template<typename ValueType>
	static void scal(const BLASIndexType UNUSED(n), const ValueType UNUSED(alpha),
			ValueType *UNUSED(x_d), const BLASIndexType UNUSED(incX)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "scal")
	}

	template<typename ValueType>
	static ValueType nrm2(const BLASIndexType UNUSED(n),
			const ValueType *UNUSED(x_d), const BLASIndexType UNUSED(incX)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "nrm2")
	}

	template<typename ValueType>
	static ValueType asum(const BLASIndexType UNUSED(n),
			const ValueType *UNUSED(x_d), BLASIndexType UNUSED(incX)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "asum")
	}

	template<typename ValueType>
	static BLASIndexType iamax(const BLASIndexType UNUSED(n),
			const ValueType *UNUSED(x_d), const BLASIndexType UNUSED(incX)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "iamax")
	}

	template<typename ValueType>
	static void swap(const BLASIndexType UNUSED(n), ValueType *UNUSED(x_d),
			const BLASIndexType UNUSED(incX), ValueType *UNUSED(y_d),
			const BLASIndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "swap")
	}

	template<typename ValueType>
	static void copy(const BLASIndexType UNUSED(n), const ValueType *UNUSED(x_d),
			const BLASIndexType UNUSED(incX), ValueType *UNUSED(y_d),
			const BLASIndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "copy")
	}

	template<typename ValueType>
	static void axpy(const BLASIndexType UNUSED(n), const ValueType UNUSED(alpha),
			const ValueType *UNUSED(x_d), const BLASIndexType UNUSED(incX),
			ValueType *UNUSED(y_d), const BLASIndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "axpy")
	}

	template<typename ValueType>
	static ValueType dot(const BLASIndexType UNUSED(n),
			const ValueType *UNUSED(x_d), const BLASIndexType UNUSED(incX),
			const ValueType *UNUSED(y_d), const BLASIndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "dot")
	}

	// ------------- BLAS2 -------------
	template<typename ValueType>
	static void gemv(const cublasOperation_t UNUSED(transA_char),
			const BLASIndexType UNUSED(m), const BLASIndexType UNUSED(n),
			const ValueType UNUSED(alpha), const ValueType* UNUSED(A),
			const BLASIndexType UNUSED(lda), const ValueType* UNUSED(x),
			const BLASIndexType UNUSED(incX), const ValueType UNUSED(beta),
			ValueType* UNUSED(y), const BLASIndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "gemv")
	}

	// ------------- BLAS3 -------------
	template<typename ValueType>
	static void gemm(const cublasOperation_t UNUSED(transA_char),
			const cublasOperation_t UNUSED(transB_char),
			const BLASIndexType UNUSED(m), const BLASIndexType UNUSED(n),
			const BLASIndexType UNUSED(k), const ValueType UNUSED(alpha),
			const ValueType* UNUSED(A), const BLASIndexType UNUSED(lda),
			const ValueType* UNUSED(B), const BLASIndexType UNUSED(ldb),
			const ValueType UNUSED(beta), ValueType* UNUSED(C),
			const BLASIndexType UNUSED(ldc)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "gemm")
	}

private:
	/**
	 * @brief convert pointer to ComplexFloat to pointer cuComplex
	 */
	static inline cuComplex* cublasCast(ComplexFloat* x) {
		return reinterpret_cast<cuComplex*>(x);
	}

	/**
	 * @brief convert pointer to ComplexDouble to pointer cuDoubleComplex
	 */
	static inline cuDoubleComplex* cublasCast(ComplexDouble* x) {
		return reinterpret_cast<cuDoubleComplex*>(x);
	}

	/**
	 * @brief convert const pointer to ComplexFloat to const pointer cuComplex
	 */
	static inline const cuComplex* cublasCast(const ComplexFloat* x) {
		return reinterpret_cast<const cuComplex*>(x);
	}

	/**
	 * @brief convert const pointer to ComplexDouble to const pointer cuDoubleComplex
	 */
	static inline const cuDoubleComplex* cublasCast(const ComplexDouble* x) {
		return reinterpret_cast<const cuDoubleComplex*>(x);
	}

	/**
	 * @brief convert value ComplexFloat to value cuComplex
	 */
	static inline cuComplex cublasCast(ComplexFloat x) {
		return *cublasCast(&x);
	}

	/**
	 * @brief convert value ComplexDouble to value cuDoubleComplex
	 */
	static inline cuDoubleComplex cublasCast(ComplexDouble x) {
		return *cublasCast(&x);
	}

};

// -------------- scal --------------
template<>
inline void CUBLASWrapper::scal<float>(const BLASIndexType n, const float alpha,
		float *x_d, const BLASIndexType incX) {
	SCAI_CUBLAS_CALL(
			cublasSscal(CUDAContext_cublasHandle, n,
					&alpha, x_d, incX),
			"cublasWrapperScale<float>");
}

template<>
inline void CUBLASWrapper::scal<double>(const BLASIndexType n, const double alpha,
		double *x_d, const BLASIndexType incX) {
	SCAI_CUBLAS_CALL(
			cublasDscal(CUDAContext_cublasHandle, n,
					&alpha, x_d, incX),
			"cublasWrapperScale<double>");
}

template<>
inline void CUBLASWrapper::scal<ComplexFloat>(const BLASIndexType n,
		const ComplexFloat alpha, ComplexFloat *x_d, const BLASIndexType incX) {
	SCAI_CUBLAS_CALL(
			cublasCscal(CUDAContext_cublasHandle, n,
					cublasCast(&alpha), cublasCast(x_d),
					incX),
			"cublasWrapperScale<ComplexFloat>");
}

template<>
inline void CUBLASWrapper::scal<ComplexDouble>(const BLASIndexType n,
		const ComplexDouble alpha, ComplexDouble *x_d, const BLASIndexType incX) {
	SCAI_CUBLAS_CALL(
			cublasZscal(CUDAContext_cublasHandle, n,
					cublasCast(&alpha), cublasCast(x_d),
					incX),
			"cublasWrapperScale<ComplexDouble>");
}

// -------------- nrm2 --------------
template<>
inline float CUBLASWrapper::nrm2<float>(const BLASIndexType n, const float *x_d,
		const BLASIndexType incX) {
	float nrm2;
	SCAI_CUBLAS_CALL(
			cublasSnrm2(CUDAContext_cublasHandle, n,
					x_d, incX, &nrm2),
			"cublasWrapperNrm2<float>");
	return nrm2;
}

template<>
inline double CUBLASWrapper::nrm2<double>(const BLASIndexType n, const double *x_d,
		const BLASIndexType incX) {
	double nrm2;
	SCAI_CUBLAS_CALL(
			cublasDnrm2(CUDAContext_cublasHandle, n,
					x_d, incX, &nrm2),
			"cublasWrapperNrm2<double>");
	return nrm2;
}

template<>
inline ComplexFloat CUBLASWrapper::nrm2<ComplexFloat>(const BLASIndexType n,
		const ComplexFloat *x_d, const BLASIndexType incX) {
	// CUBLAS returns only float result so we convert it back to Complex
	float nrm2;
	SCAI_CUBLAS_CALL(
			cublasScnrm2(CUDAContext_cublasHandle,
					n, cublasCast(x_d),
					incX, &nrm2),
			"cublasWrapperNrm2<ComplexFloat>");
	return ComplexFloat(nrm2, 0.0f);
}

template<>
inline ComplexDouble CUBLASWrapper::nrm2<ComplexDouble>(const BLASIndexType n,
		const ComplexDouble *x_d, const BLASIndexType incX) {
	// CUBLAS returns only double result so we convert it back to Complex
	double nrm2;
	SCAI_CUBLAS_CALL(
			cublasDznrm2(CUDAContext_cublasHandle,
					n, cublasCast(x_d),
					incX, &nrm2),
			"cublasWrapperNrm2<ComplexDouble>");
	return ComplexDouble(nrm2, 0.0);
}

// -------------- asum --------------
template<>
inline float CUBLASWrapper::asum<float>(const BLASIndexType n, const float *x_d,
		BLASIndexType incX) {
	float asum;
	SCAI_CUBLAS_CALL(
			cublasSasum(CUDAContext_cublasHandle, n,
					x_d, incX, &asum),
			"cublasWrapperAsum<float>");
	return asum;
}

template<>
inline double CUBLASWrapper::asum<double>(const BLASIndexType n, const double *x_d,
		BLASIndexType incX) {
	double asum;
	SCAI_CUBLAS_CALL(
			cublasDasum(CUDAContext_cublasHandle, n,
					x_d, incX, &asum),
			"cublasWrapperAsum<double>");
	return asum;
}

template<>
inline ComplexFloat CUBLASWrapper::asum<ComplexFloat>(const BLASIndexType n,
		const ComplexFloat *x_d, BLASIndexType incX) {
	// CUBLAS returns only float result so we convert it back to Complex
	float asum;
	SCAI_CUBLAS_CALL(
			cublasScasum(CUDAContext_cublasHandle,
					n, cublasCast(x_d),
					incX, &asum),
			"cublasWrapperAsum<ComplexFloat>");
	return ComplexFloat(asum, 0.0f);
}

template<>
inline ComplexDouble CUBLASWrapper::asum<ComplexDouble>(const BLASIndexType n,
		const ComplexDouble *x_d, BLASIndexType incX) {
	// CUBLAS returns only double result so we convert it back to Complex
	double asum;
	SCAI_CUBLAS_CALL(
			cublasDzasum(CUDAContext_cublasHandle,
					n, cublasCast(x_d),
					incX, &asum),
			"cublasWrapperAsum<ComplexDouble>");
	return ComplexDouble(asum, 0.0);
}

// -------------- iamax --------------
template<>
inline CUBLASWrapper::BLASIndexType CUBLASWrapper::iamax<float>(const BLASIndexType n,
		const float *x_d, const BLASIndexType incX) {
	BLASIndexType iamax;
	SCAI_CUBLAS_CALL(
			cublasIsamax(CUDAContext_cublasHandle,
					n, x_d,
					incX, &iamax),
			"cublasWrapperIamax<float>");
	return iamax;
}

template<>
inline CUBLASWrapper::BLASIndexType CUBLASWrapper::iamax<double>(const BLASIndexType n,
		const double *x_d, const BLASIndexType incX) {
	BLASIndexType iamax;
	SCAI_CUBLAS_CALL(
			cublasIdamax(CUDAContext_cublasHandle,
					n, x_d,
					incX, &iamax),
			"cublasWrapperIamax<double>");
	return iamax;
}

template<>
inline CUBLASWrapper::BLASIndexType CUBLASWrapper::iamax<ComplexFloat>(const BLASIndexType n,
		const ComplexFloat *x_d, const BLASIndexType incX) {
	BLASIndexType iamax;
	SCAI_CUBLAS_CALL(
			cublasIcamax(CUDAContext_cublasHandle,
					n, cublasCast(x_d),
					incX, &iamax),
			"cublasWrapperIamax<ComplexFloat>");
	return iamax;
}

template<>
inline CUBLASWrapper::BLASIndexType CUBLASWrapper::iamax<ComplexDouble>(const BLASIndexType n,
		const ComplexDouble *x_d, const BLASIndexType incX) {
	BLASIndexType iamax;
	SCAI_CUBLAS_CALL(
			cublasIzamax(CUDAContext_cublasHandle,
					n, cublasCast(x_d),
					incX, &iamax),
			"cublasWrapperIamax<ComplexDouble>");
	return iamax;
}

// -------------- swap --------------
template<>
inline void CUBLASWrapper::swap<float>(const BLASIndexType n, float *x_d,
		const BLASIndexType incX, float *y_d, const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasSswap(CUDAContext_cublasHandle, n,
					x_d, incX, y_d,
					incY),
			"cublasWrapperSwap<float>");
}

template<>
inline void CUBLASWrapper::swap<double>(const BLASIndexType n, double *x_d,
		const BLASIndexType incX, double *y_d, const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasDswap(CUDAContext_cublasHandle, n,
					x_d, incX, y_d,
					incY),
			"cublasWrapperSwap<double>");
}

template<>
inline void CUBLASWrapper::swap<ComplexFloat>(const BLASIndexType n,
		ComplexFloat *x_d, const BLASIndexType incX, ComplexFloat *y_d,
		const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasCswap(CUDAContext_cublasHandle, n,
					cublasCast(x_d), incX,
					cublasCast(y_d), incY),
			"cublasWrapperSwap<ComplexFloat>");
}

template<>
inline void CUBLASWrapper::swap<ComplexDouble>(const BLASIndexType n,
		ComplexDouble *x_d, const BLASIndexType incX, ComplexDouble *y_d,
		const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasZswap(CUDAContext_cublasHandle, n,
					cublasCast(x_d), incX,
					cublasCast(y_d), incY),
			"cublasWrapperSwap<ComplexDouble>");
}

// -------------- copy --------------
template<>
inline void CUBLASWrapper::copy<float>(const BLASIndexType n, const float *x_d,
		const BLASIndexType incX, float *y_d, const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasScopy(CUDAContext_cublasHandle, n,
					x_d, incX, y_d,
					incY),
			"cublasWrapperCopy<float>");
}

template<>
inline void CUBLASWrapper::copy<double>(const BLASIndexType n, const double *x_d,
		const BLASIndexType incX, double *y_d, const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasDcopy(CUDAContext_cublasHandle, n,
					x_d, incX, y_d,
					incY),
			"cublasWrapperCopy<double>");
}

template<>
inline void CUBLASWrapper::copy<ComplexFloat>(const BLASIndexType n,
		const ComplexFloat *x_d, const BLASIndexType incX, ComplexFloat *y_d,
		const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasCcopy(CUDAContext_cublasHandle, n,
					cublasCast(x_d), incX,
					cublasCast(y_d), incY),
			"cublasWrapperCopy<ComplexFloat>");
}

template<>
inline void CUBLASWrapper::copy<ComplexDouble>(const BLASIndexType n,
		const ComplexDouble *x_d, const BLASIndexType incX, ComplexDouble *y_d,
		const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasZcopy(CUDAContext_cublasHandle, n,
					cublasCast(x_d), incX,
					cublasCast(y_d), incY),
			"cublasWrapperCopy<ComplexDouble>");
}

// -------------- axpy --------------
template<>
inline void CUBLASWrapper::axpy<float>(const BLASIndexType n, const float alpha,
		const float *x_d, const BLASIndexType incX, float *y_d,
		const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasSaxpy(CUDAContext_cublasHandle, n,
					&alpha, x_d, incX, y_d,
					incY),
			"cublasWrapperAxpy<float>");
}

template<>
inline void CUBLASWrapper::axpy<double>(const BLASIndexType n, const double alpha,
		const double *x_d, const BLASIndexType incX, double *y_d,
		const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasDaxpy(CUDAContext_cublasHandle, n,
					&alpha, x_d, incX, y_d,
					incY),
			"cublasWrapperAxpy<double>");
}

template<>
inline void CUBLASWrapper::axpy<ComplexFloat>(const BLASIndexType n,
		const ComplexFloat alpha, const ComplexFloat *x_d, const BLASIndexType incX,
		ComplexFloat *y_d, const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasCaxpy(CUDAContext_cublasHandle, n,
					cublasCast(&alpha), cublasCast(x_d),
					incX, cublasCast(y_d),
					incY),
			"cublasWrapperAxpy<ComplexFloat>");
}

template<>
inline void CUBLASWrapper::axpy<ComplexDouble>(const BLASIndexType n,
		const ComplexDouble alpha, const ComplexDouble *x_d,
		const BLASIndexType incX, ComplexDouble *y_d, const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasZaxpy(CUDAContext_cublasHandle, n,
					cublasCast(&alpha), cublasCast(x_d),
					incX, cublasCast(y_d),
					incY),
			"cublasWrapperAxpy<ComplexDouble>");
}

// -------------- dot --------------
template<>
inline float CUBLASWrapper::dot<float>(const BLASIndexType n, const float *x_d,
		const BLASIndexType incX, const float *y_d, const BLASIndexType incY) {
	float dot;
	SCAI_CUBLAS_CALL(
			cublasSdot(CUDAContext_cublasHandle, n,
					x_d, incX, y_d,
					incY, &dot),
			"cublasWrapperDot<float>");
	return dot;
}

template<>
inline double CUBLASWrapper::dot<double>(const BLASIndexType n, const double *x_d,
		const BLASIndexType incX, const double *y_d, const BLASIndexType incY) {
	double dot;
	SCAI_CUBLAS_CALL(
			cublasDdot(CUDAContext_cublasHandle, n,
					x_d, incX, y_d,
					incY, &dot),
			"cublasWrapperDot<double>");
	return dot;
}

template<>
inline ComplexFloat CUBLASWrapper::dot<ComplexFloat>(const BLASIndexType n,
		const ComplexFloat *x_d, const BLASIndexType incX, const ComplexFloat *y_d,
		const BLASIndexType incY) {
	ComplexFloat dot;
	SCAI_CUBLAS_CALL(
			cublasCdotc(CUDAContext_cublasHandle, n,
					cublasCast(x_d), incX,
					cublasCast(y_d), incY,
					cublasCast(&dot)), "cublasWrapperDot<ComplexFloat>");
	return dot;
}

template<>
inline ComplexDouble CUBLASWrapper::dot<ComplexDouble>(const BLASIndexType n,
		const ComplexDouble *x_d, const BLASIndexType incX,
		const ComplexDouble *y_d, const BLASIndexType incY) {
	ComplexDouble dot;
	SCAI_CUBLAS_CALL(
			cublasZdotc(CUDAContext_cublasHandle, n,
					cublasCast(x_d), incX,
					cublasCast(y_d), incY,
					cublasCast(&dot)), "cublasWrapperDot<ComplexDouble>");
	return dot;
}

// -------------- gemv --------------
template<>
inline void CUBLASWrapper::gemv<float>(const cublasOperation_t trans,
		const BLASIndexType m, const BLASIndexType n, const float alpha, const float* A,
		const BLASIndexType lda, const float* x, const BLASIndexType incX,
		const float beta, float* y, const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasSgemv(CUDAContext_cublasHandle, trans,
					m,
					n, &alpha, A,
					lda, x,
					incX, &beta, y,
					incY),
			"cublasWrapperGemv<float>");
}

template<>
inline void CUBLASWrapper::gemv<double>(const cublasOperation_t trans,
		const BLASIndexType m, const BLASIndexType n, const double alpha,
		const double* A, const BLASIndexType lda, const double* x,
		const BLASIndexType incX, const double beta, double* y,
		const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasDgemv(CUDAContext_cublasHandle, trans,
					m,
					n, &alpha, A,
					lda, x,
					incX, &beta, y,
					incY),
			"cublasWrapperGemv<double>");
}

template<>
inline void CUBLASWrapper::gemv<ComplexFloat>(const cublasOperation_t trans,
		const BLASIndexType m, const BLASIndexType n, const ComplexFloat alpha,
		const ComplexFloat* A, const BLASIndexType lda, const ComplexFloat* x,
		const BLASIndexType incX, const ComplexFloat beta, ComplexFloat* y,
		const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasCgemv(CUDAContext_cublasHandle, trans,
					m,
					n, cublasCast(&alpha),
					cublasCast(A), lda,
					cublasCast(x), incX,
					cublasCast(&beta), cublasCast(y),
					incY),
			"cublasWrapperGemv<ComplexFloat>");
}

template<>
inline void CUBLASWrapper::gemv<ComplexDouble>(const cublasOperation_t trans,
		const BLASIndexType m, const BLASIndexType n, const ComplexDouble alpha,
		const ComplexDouble* A, const BLASIndexType lda, const ComplexDouble* x,
		const BLASIndexType incX, const ComplexDouble beta, ComplexDouble* y,
		const BLASIndexType incY) {
	SCAI_CUBLAS_CALL(
			cublasZgemv(CUDAContext_cublasHandle, trans,
					m,
					n, cublasCast(&alpha),
					cublasCast(A), lda,
					cublasCast(x), incX,
					cublasCast(&beta), cublasCast(y),
					incY),
			"cublasWrapperGemv<ComplexDouble>");
}

// -------------- gemm --------------
template<>
inline void CUBLASWrapper::gemm<float>(const cublasOperation_t transA_char,
		const cublasOperation_t transB_char, const BLASIndexType m,
		const BLASIndexType n, const BLASIndexType k, const float alpha, const float* A,
		const BLASIndexType lda, const float* B, const BLASIndexType ldb,
		const float beta, float* C, const BLASIndexType ldc) {
	SCAI_CUBLAS_CALL(
			cublasSgemm(CUDAContext_cublasHandle, transA_char, transB_char,
					m,
					n,
					k, &alpha, A,
					lda, B,
					ldb, &beta, C,
					ldc),
			"cublasWrapperGemm<float>");
}

template<>
inline void CUBLASWrapper::gemm<double>(const cublasOperation_t transA_char,
		const cublasOperation_t transB_char, const BLASIndexType m,
		const BLASIndexType n, const BLASIndexType k, const double alpha,
		const double* A, const BLASIndexType lda, const double* B,
		const BLASIndexType ldb, const double beta, double* C,
		const BLASIndexType ldc) {
	SCAI_CUBLAS_CALL(
			cublasDgemm(CUDAContext_cublasHandle, transA_char, transB_char,
					m,
					n,
					k, &alpha, A,
					lda, B,
					ldb, &beta, C,
					ldc),
			"cublasWrapperGemm<dobule>");
}

template<>
inline void CUBLASWrapper::gemm<ComplexFloat>(
		const cublasOperation_t transA_char,
		const cublasOperation_t transB_char, const BLASIndexType m,
		const BLASIndexType n, const BLASIndexType k, const ComplexFloat alpha,
		const ComplexFloat* A, const BLASIndexType lda, const ComplexFloat* B,
		const BLASIndexType ldb, const ComplexFloat beta, ComplexFloat* C,
		const BLASIndexType ldc) {
	SCAI_CUBLAS_CALL(
			cublasCgemm(CUDAContext_cublasHandle, transA_char, transB_char,
					m,
					n,
					k, cublasCast(&alpha),
					cublasCast(A), lda,
					cublasCast(B), ldb,
					cublasCast(&beta), cublasCast(C),
					ldc),
			"cublasWrapperGemm<ComplexFloat>");
}

template<>
inline void CUBLASWrapper::gemm<ComplexDouble>(
		const cublasOperation_t transA_char,
		const cublasOperation_t transB_char, const BLASIndexType m,
		const BLASIndexType n, const BLASIndexType k, const ComplexDouble alpha,
		const ComplexDouble* A, const BLASIndexType lda, const ComplexDouble* B,
		const BLASIndexType ldb, const ComplexDouble beta, ComplexDouble* C,
		const BLASIndexType ldc) {
	SCAI_CUBLAS_CALL(
			cublasZgemm(CUDAContext_cublasHandle, transA_char, transB_char,
					m,
					n,
					k, cublasCast(&alpha),
					cublasCast(A), lda,
					cublasCast(B), ldb,
					cublasCast(&beta), cublasCast(C),
					ldc),
			"cublasWrapperGemm<double>");
}

} /* end namespace blaskernel */

} /* end namespace scai */

