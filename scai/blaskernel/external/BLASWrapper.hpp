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
#include <scai/common/SCAITypes.hpp>

namespace scai {

namespace blaskernel {

class COMMON_DLL_IMPORTEXPORT BLASWrapper {
public:
	typedef int BLASIndexType;

	// ------------- BLAS1 -------------
	template<typename ValueType>
	static void scal(const IndexType UNUSED(n), const ValueType UNUSED(alpha),
			ValueType *UNUSED(x), const IndexType UNUSED(incX)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "scal");
	}

	template<typename ValueType>
	static ValueType nrm2(const IndexType UNUSED(n), const ValueType *UNUSED(x),
			const IndexType UNUSED(incX)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "nrm2");
	}

	template<typename ValueType>
	static ValueType asum(const IndexType UNUSED(n), const ValueType *UNUSED(x),
			IndexType UNUSED(incX)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "asum");
	}

	template<typename ValueType>
	static IndexType iamax(const IndexType UNUSED(n),
			const ValueType *UNUSED(x), const IndexType UNUSED(incX)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "iamax");
	}

	template<typename ValueType>
	static void swap(const IndexType UNUSED(n), ValueType *UNUSED(x),
			const IndexType UNUSED(incX), ValueType *UNUSED(y),
			const IndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "swap");
	}

	template<typename ValueType>
	static void copy(const IndexType UNUSED(n), const ValueType *UNUSED(x),
			const IndexType UNUSED(incX), ValueType *UNUSED(y),
			const IndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "copy");
	}

	template<typename ValueType>
	static void axpy(const IndexType UNUSED(n), const ValueType UNUSED(alpha),
			const ValueType *UNUSED(x), const IndexType UNUSED(incX),
			ValueType *UNUSED(y), const IndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "axpy");
	}

	template<typename ValueType>
	static ValueType dot(const IndexType UNUSED(n), const ValueType *UNUSED(x),
			const IndexType UNUSED(incX), const ValueType *UNUSED(y),
			const IndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "dot");
	}

	// ------------- BLAS2 -------------
	template<typename ValueType>
	static void gemv(const CBLAS_ORDER order,
			const CBLAS_TRANSPOSE UNUSED(transA_char),
			const IndexType UNUSED(m), const IndexType UNUSED(n),
			const ValueType UNUSED(alpha), const ValueType* UNUSED(A),
			const IndexType UNUSED(lda), const ValueType* UNUSED(x),
			const IndexType UNUSED(incX), const ValueType UNUSED(beta),
			ValueType* UNUSED(y), const IndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "gemv");
	}

	// ------------- BLAS3 -------------
	template<typename ValueType>
	static void gemm(const CBLAS_ORDER UNUSED(order),
			const CBLAS_TRANSPOSE UNUSED(transA_char),
			const CBLAS_TRANSPOSE UNUSED(transB_char),
			const IndexType UNUSED(m), const IndexType UNUSED(n),
			const IndexType UNUSED(k), const ValueType UNUSED(alpha),
			const ValueType* UNUSED(A), const IndexType UNUSED(lda),
			const ValueType* UNUSED(B), const IndexType UNUSED(ldb),
			const ValueType UNUSED(beta), ValueType* UNUSED(C),
			const IndexType UNUSED(ldc)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "gemm");
	}

private:

};

// -------------- scal --------------
template<>
inline void BLASWrapper::scal<float>(const IndexType n, const float alpha,
		float *x, const IndexType incX) {
	cblas_sscal(static_cast<BLASIndexType>(n), alpha, x,
			static_cast<BLASIndexType>(incX));
}

template<>
inline void BLASWrapper::scal<double>(const IndexType n, const double alpha,
		double *x, const IndexType incX) {
	cblas_dscal(static_cast<BLASIndexType>(n), alpha, x,
			static_cast<BLASIndexType>(incX));
}

template<>
inline void BLASWrapper::scal<ComplexFloat>(const IndexType n,
		const ComplexFloat alpha, ComplexFloat *x, const IndexType incX) {
	// Attention: alpha is here passed by a pointer
	cblas_cscal(static_cast<BLASIndexType>(n), &alpha, x,
			static_cast<BLASIndexType>(incX));
}

template<>
inline void BLASWrapper::scal<ComplexDouble>(const IndexType n,
		const ComplexDouble alpha, ComplexDouble *x, const IndexType incX) {
	// Attention: alpha is here passed by a pointer
	cblas_zscal(static_cast<BLASIndexType>(n), &alpha, x,
			static_cast<BLASIndexType>(incX));
}

// -------------- nrm2 --------------
template<>
inline float BLASWrapper::nrm2<float>(const IndexType n, const float *x,
		const IndexType incX) {
	return cblas_snrm2(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX));
}

template<>
inline double BLASWrapper::nrm2<double>(const IndexType n, const double *x,
		const IndexType incX) {
	return cblas_dnrm2(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX));
}

template<>
inline ComplexFloat BLASWrapper::nrm2<ComplexFloat>(const IndexType n,
		const ComplexFloat *x, const IndexType incX) {
	float res = cblas_scnrm2(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX));
	return ComplexFloat(res);
}

template<>
inline ComplexDouble BLASWrapper::nrm2<ComplexDouble>(const IndexType n,
		const ComplexDouble *x, const IndexType incX) {
	double res = cblas_dznrm2(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX));
	return ComplexDouble(res);
}

// -------------- asum --------------
template<>
inline float BLASWrapper::asum<float>(const IndexType n, const float *x,
		IndexType incX) {
	return cblas_sasum(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX));
}

template<>
inline double BLASWrapper::asum<double>(const IndexType n, const double *x,
		IndexType incX) {
	return cblas_dasum(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX));
}

template<>
inline ComplexFloat BLASWrapper::asum<ComplexFloat>(const IndexType n,
		const ComplexFloat *x, IndexType incX) {
	float res = cblas_scasum(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX));
	return ComplexFloat(res);
}

template<>
inline ComplexDouble BLASWrapper::asum<ComplexDouble>(const IndexType n,
		const ComplexDouble *x, IndexType incX) {
	double res = cblas_dzasum(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX));
	return ComplexDouble(res);
}

// -------------- iamax --------------
template<>
inline IndexType BLASWrapper::iamax<float>(const IndexType n, const float *x,
		const IndexType incX) {
	return cblas_isamax(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX));
}

template<>
inline IndexType BLASWrapper::iamax<double>(const IndexType n, const double *x,
		const IndexType incX) {
	return cblas_idamax(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX));
}

template<>
inline IndexType BLASWrapper::iamax<ComplexFloat>(const IndexType n,
		const ComplexFloat *x, const IndexType incX) {
	return cblas_icamax(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX));
}

template<>
inline IndexType BLASWrapper::iamax<ComplexDouble>(const IndexType n,
		const ComplexDouble *x, const IndexType incX) {
	return cblas_izamax(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX));
}

// -------------- swap --------------
template<>
inline void BLASWrapper::swap<float>(const IndexType n, float *x,
		const IndexType incX, float *y, const IndexType incY) {
	cblas_sswap(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline void BLASWrapper::swap<double>(const IndexType n, double *x,
		const IndexType incX, double *y, const IndexType incY) {
	cblas_dswap(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline void BLASWrapper::swap<ComplexFloat>(const IndexType n, ComplexFloat *x,
		const IndexType incX, ComplexFloat *y, const IndexType incY) {
	cblas_cswap(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline void BLASWrapper::swap<ComplexDouble>(const IndexType n,
		ComplexDouble *x, const IndexType incX, ComplexDouble *y,
		const IndexType incY) {
	cblas_zswap(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

// -------------- copy --------------
template<>
inline void BLASWrapper::copy<float>(const IndexType n, const float *x,
		const IndexType incX, float *y, const IndexType incY) {
	cblas_scopy(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline void BLASWrapper::copy<double>(const IndexType n, const double *x,
		const IndexType incX, double *y, const IndexType incY) {
	cblas_dcopy(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline void BLASWrapper::copy<ComplexFloat>(const IndexType n,
		const ComplexFloat *x, const IndexType incX, ComplexFloat *y,
		const IndexType incY) {
	cblas_ccopy(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline void BLASWrapper::copy<ComplexDouble>(const IndexType n,
		const ComplexDouble *x, const IndexType incX, ComplexDouble *y,
		const IndexType incY) {
	cblas_zcopy(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

// -------------- axpy --------------
template<>
inline void BLASWrapper::axpy<float>(const IndexType n, const float alpha,
		const float *x, const IndexType incX, float *y, const IndexType incY) {
	cblas_saxpy(static_cast<BLASIndexType>(n), alpha, x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline void BLASWrapper::axpy<double>(const IndexType n, const double alpha,
		const double *x, const IndexType incX, double *y,
		const IndexType incY) {
	cblas_daxpy(static_cast<BLASIndexType>(n), alpha, x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline void BLASWrapper::axpy<ComplexFloat>(const IndexType n,
		const ComplexFloat alpha, const ComplexFloat *x, const IndexType incX,
		ComplexFloat *y, const IndexType incY) {
	// Attention: alpha is here passed by a pointer
	cblas_caxpy(static_cast<BLASIndexType>(n), &alpha, x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline void BLASWrapper::axpy<ComplexDouble>(const IndexType n,
		const ComplexDouble alpha, const ComplexDouble *x, const IndexType incX,
		ComplexDouble *y, const IndexType incY) {
	// Attention: alpha is here passed by a pointer
	cblas_zaxpy(static_cast<BLASIndexType>(n), &alpha, x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

// -------------- dot --------------
template<>
inline float BLASWrapper::dot<float>(const IndexType n, const float *x,
		const IndexType incX, const float *y, const IndexType incY) {
	return cblas_sdot(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline double BLASWrapper::dot<double>(const IndexType n, const double *x,
		const IndexType incX, const double *y, const IndexType incY) {
	return cblas_ddot(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline ComplexFloat BLASWrapper::dot<ComplexFloat>(const IndexType n,
		const ComplexFloat *x, const IndexType incX, const ComplexFloat *y,
		const IndexType incY) {
	ComplexFloat dotu;
	cblas_cdotu_sub(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY), &dotu);
	return dotu;
}

template<>
inline ComplexDouble BLASWrapper::dot<ComplexDouble>(const IndexType n,
		const ComplexDouble *x, const IndexType incX, const ComplexDouble *y,
		const IndexType incY) {
	ComplexDouble dotu;
	cblas_zdotu_sub(static_cast<BLASIndexType>(n), x,
			static_cast<BLASIndexType>(incX), y,
			static_cast<BLASIndexType>(incY), &dotu);
	return dotu;
}

// -------------- gemv --------------
template<>
inline void BLASWrapper::gemv<float>(const CBLAS_ORDER order,
		const CBLAS_TRANSPOSE transA, const IndexType m, const IndexType n,
		const float alpha, const float* A, const IndexType lda, const float* x,
		const IndexType incX, const float beta, float* y,
		const IndexType incY) {
	cblas_sgemv(order, transA, m, static_cast<BLASIndexType>(n), alpha, A, lda,
			x, static_cast<BLASIndexType>(incX), beta, y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline void BLASWrapper::gemv<double>(const CBLAS_ORDER order,
		const CBLAS_TRANSPOSE transA, const IndexType m, const IndexType n,
		const double alpha, const double* A, const IndexType lda,
		const double* x, const IndexType incX, const double beta, double* y,
		const IndexType incY) {
	cblas_dgemv(order, transA, m, static_cast<BLASIndexType>(n), alpha, A, lda,
			x, static_cast<BLASIndexType>(incX), beta, y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline void BLASWrapper::gemv<ComplexFloat>(const CBLAS_ORDER order,
		const CBLAS_TRANSPOSE transA, const IndexType m, const IndexType n,
		const ComplexFloat alpha, const ComplexFloat* A, const IndexType lda,
		const ComplexFloat* x, const IndexType incX, const ComplexFloat beta,
		ComplexFloat* y, const IndexType incY) {
	// Attention: alpha, beta must be passed here as a pointer
	cblas_cgemv(order, transA, m, static_cast<BLASIndexType>(n), &alpha, A, lda,
			x, static_cast<BLASIndexType>(incX), &beta, y,
			static_cast<BLASIndexType>(incY));
}

template<>
inline void BLASWrapper::gemv<ComplexDouble>(const CBLAS_ORDER order,
		const CBLAS_TRANSPOSE transA, const IndexType m, const IndexType n,
		const ComplexDouble alpha, const ComplexDouble* A, const IndexType lda,
		const ComplexDouble* x, const IndexType incX, const ComplexDouble beta,
		ComplexDouble* y, const IndexType incY) {
	// Attention: alpha, beta must be passed here as a pointer
	cblas_zgemv(order, transA, m, static_cast<BLASIndexType>(n), &alpha, A, lda,
			x, static_cast<BLASIndexType>(incX), &beta, y,
			static_cast<BLASIndexType>(incY));
}

// -------------- gemm --------------
template<>
inline void BLASWrapper::gemm<float>(const CBLAS_ORDER order,
		const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB,
		const IndexType m, const IndexType n, const IndexType k,
		const float alpha, const float* A, const IndexType lda, const float* B,
		const IndexType ldb, const float beta, float* C, const IndexType ldc) {
	cblas_sgemm(order, transA, transB, static_cast<BLASIndexType>(m), static_cast<BLASIndexType>(n), static_cast<BLASIndexType>(k),
			alpha, A, static_cast<BLASIndexType>(lda), B, static_cast<BLASIndexType>(ldb), beta, C, static_cast<BLASIndexType>(ldc));
}

template<>
inline void BLASWrapper::gemm<double>(const CBLAS_ORDER order,
		const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB,
		const IndexType m, const IndexType n, const IndexType k,
		const double alpha, const double* A, const IndexType lda,
		const double* B, const IndexType ldb, const double beta, double* C,
		const IndexType ldc) {
	cblas_dgemm(order, transA, transB, static_cast<BLASIndexType>(m), static_cast<BLASIndexType>(n), static_cast<BLASIndexType>(k),
			alpha, A, static_cast<BLASIndexType>(lda), B, static_cast<BLASIndexType>(ldb), beta, C, static_cast<BLASIndexType>(ldc));
}

template<>
inline void BLASWrapper::gemm<ComplexFloat>(const CBLAS_ORDER order,
		const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB,
		const IndexType m, const IndexType n, const IndexType k,
		const ComplexFloat alpha, const ComplexFloat* A, const IndexType lda,
		const ComplexFloat* B, const IndexType ldb, const ComplexFloat beta,
		ComplexFloat* C, const IndexType ldc) {
	// Attention: alpha and beta are passed by a pointer
	cblas_cgemm(order, transA, transB, static_cast<BLASIndexType>(m), static_cast<BLASIndexType>(n), static_cast<BLASIndexType>(k),
			&alpha, A, static_cast<BLASIndexType>(lda), B, static_cast<BLASIndexType>(ldb), &beta, C, static_cast<BLASIndexType>(ldc));
}

template<>
inline void BLASWrapper::gemm<ComplexDouble>(const CBLAS_ORDER order,
		const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB,
		const IndexType m, const IndexType n, const IndexType k,
		const ComplexDouble alpha, const ComplexDouble* A, const IndexType lda,
		const ComplexDouble* B, const IndexType ldb, const ComplexDouble beta,
		ComplexDouble* C, const IndexType ldc) {
	// Attention: alpha and beta are passed by a pointer
	cblas_zgemm(order, transA, transB, static_cast<BLASIndexType>(m), static_cast<BLASIndexType>(n), static_cast<BLASIndexType>(k),
			&alpha, A, static_cast<BLASIndexType>(lda), B, static_cast<BLASIndexType>(ldb), &beta, C, static_cast<BLASIndexType>(ldc));
}

} /* end namespace blaskernel */

} /* end namespace scai */

