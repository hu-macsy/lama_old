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

class COMMON_DLL_IMPORTEXPORT BLASWrapperOLD {
public:
	typedef int BLASIndexType;

	// ------------- BLAS1 -------------
	template<typename ValueType>
	static void scal(const BLASIndexType UNUSED(n), const ValueType UNUSED(alpha),
			ValueType *UNUSED(x), const BLASIndexType UNUSED(incX)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "scal");
	}

	template<typename ValueType>
	static ValueType nrm2(const BLASIndexType UNUSED(n), const ValueType *UNUSED(x),
			const BLASIndexType UNUSED(incX)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "nrm2");
	}

	template<typename ValueType>
	static ValueType asum(const BLASIndexType UNUSED(n), const ValueType *UNUSED(x),
			BLASIndexType UNUSED(incX)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "asum");
	}

	template<typename ValueType>
	static BLASIndexType iamax(const BLASIndexType UNUSED(n),
			const ValueType *UNUSED(x), const BLASIndexType UNUSED(incX)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "iamax");
	}

	template<typename ValueType>
	static void swap(const BLASIndexType UNUSED(n), ValueType *UNUSED(x),
			const BLASIndexType UNUSED(incX), ValueType *UNUSED(y),
			const BLASIndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "swap");
	}

	template<typename ValueType>
	static void copy(const BLASIndexType UNUSED(n), const ValueType *UNUSED(x),
			const BLASIndexType UNUSED(incX), ValueType *UNUSED(y),
			const BLASIndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "copy");
	}

	template<typename ValueType>
	static void axpy(const BLASIndexType UNUSED(n), const ValueType UNUSED(alpha),
			const ValueType *UNUSED(x), const BLASIndexType UNUSED(incX),
			ValueType *UNUSED(y), const BLASIndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "axpy");
	}

	template<typename ValueType>
	static ValueType dot(const BLASIndexType UNUSED(n), const ValueType *UNUSED(x),
			const BLASIndexType UNUSED(incX), const ValueType *UNUSED(y),
			const BLASIndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "dot");
	}

	// ------------- BLAS2 -------------
	template<typename ValueType>
	static void gemv(const CBLAS_ORDER order,
			const CBLAS_TRANSPOSE UNUSED(transA_char),
			const BLASIndexType UNUSED(m), const BLASIndexType UNUSED(n),
			const ValueType UNUSED(alpha), const ValueType* UNUSED(A),
			const BLASIndexType UNUSED(lda), const ValueType* UNUSED(x),
			const BLASIndexType UNUSED(incX), const ValueType UNUSED(beta),
			ValueType* UNUSED(y), const BLASIndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "gemv");
	}

	// ------------- BLAS3 -------------
	template<typename ValueType>
	static void gemm(const CBLAS_ORDER UNUSED(order),
			const CBLAS_TRANSPOSE UNUSED(transA_char),
			const CBLAS_TRANSPOSE UNUSED(transB_char),
			const BLASIndexType UNUSED(m), const BLASIndexType UNUSED(n),
			const BLASIndexType UNUSED(k), const ValueType UNUSED(alpha),
			const ValueType* UNUSED(A), const BLASIndexType UNUSED(lda),
			const ValueType* UNUSED(B), const BLASIndexType UNUSED(ldb),
			const ValueType UNUSED(beta), ValueType* UNUSED(C),
			const BLASIndexType UNUSED(ldc)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "gemm");
	}

private:

};
//
//// -------------- scal --------------
//template<>
//inline void BLASWrapper::scal<float>(const BLASIndexType n, const float alpha,
//		float *x, const BLASIndexType incX) {
//	cblas_sscal(n, alpha, x,
//			incX);
//}
//
//template<>
//inline void BLASWrapper::scal<double>(const BLASIndexType n, const double alpha,
//		double *x, const BLASIndexType incX) {
//	cblas_dscal(n, alpha, x,
//			incX);
//}
//
//#ifdef SCAI_COMPLEX_SUPPORTED
//
//template<>
//inline void BLASWrapper::scal<ComplexFloat>(const BLASIndexType n,
//		const ComplexFloat alpha, ComplexFloat *x, const BLASIndexType incX) {
//	// Attention: alpha is here passed by a pointer
//	cblas_cscal(n, &alpha, x,
//			incX);
//}
//
//template<>
//inline void BLASWrapper::scal<ComplexDouble>(const BLASIndexType n,
//		const ComplexDouble alpha, ComplexDouble *x, const BLASIndexType incX) {
//	// Attention: alpha is here passed by a pointer
//	cblas_zscal(n, &alpha, x,
//			incX);
//}
//
//#endif
//
//// -------------- nrm2 --------------
//template<>
//inline float BLASWrapper::nrm2<float>(const BLASIndexType n, const float *x,
//		const BLASIndexType incX) {
//	return cblas_snrm2(n, x,
//			incX);
//}
//
//template<>
//inline double BLASWrapper::nrm2<double>(const BLASIndexType n, const double *x,
//		const BLASIndexType incX) {
//	return cblas_dnrm2(n, x,
//			incX);
//}
//
//#ifdef SCAI_COMPLEX_SUPPORTED
//
//template<>
//inline ComplexFloat BLASWrapper::nrm2<ComplexFloat>(const BLASIndexType n,
//		const ComplexFloat *x, const BLASIndexType incX) {
//	float res = cblas_scnrm2(n, x,
//			incX);
//	return ComplexFloat(res);
//}
//
//template<>
//inline ComplexDouble BLASWrapper::nrm2<ComplexDouble>(const BLASIndexType n,
//		const ComplexDouble *x, const BLASIndexType incX) {
//	double res = cblas_dznrm2(n, x,
//			incX);
//	return ComplexDouble(res);
//}
//
//#endif
//
//// -------------- asum --------------
//template<>
//inline float BLASWrapper::asum<float>(const BLASIndexType n, const float *x,
//		BLASIndexType incX) {
//	return cblas_sasum(n, x,
//			incX);
//}
//
//template<>
//inline double BLASWrapper::asum<double>(const BLASIndexType n, const double *x,
//		BLASIndexType incX) {
//	return cblas_dasum(n, x,
//			incX);
//}
//
//#ifdef SCAI_COMPLEX_SUPPORTED
//
//template<>
//inline ComplexFloat BLASWrapper::asum<ComplexFloat>(const BLASIndexType n,
//		const ComplexFloat *x, BLASIndexType incX) {
//	float res = cblas_scasum(n, x,
//			incX);
//	return ComplexFloat(res);
//}
//
//template<>
//inline ComplexDouble BLASWrapper::asum<ComplexDouble>(const BLASIndexType n,
//		const ComplexDouble *x, BLASIndexType incX) {
//	double res = cblas_dzasum(n, x,
//			incX);
//	return ComplexDouble(res);
//}
//
//#endif
//
//// -------------- iamax --------------
//template<>
//inline BLASWrapper::BLASIndexType BLASWrapper::iamax<float>(const BLASIndexType n, const float *x,
//		const BLASIndexType incX) {
//	return cblas_isamax(n, x,
//			incX);
//}
//
//template<>
//inline BLASWrapper::BLASIndexType BLASWrapper::iamax<double>(const BLASIndexType n, const double *x,
//		const BLASIndexType incX) {
//	return cblas_idamax(n, x,
//			incX);
//}
//
//#ifdef SCAI_COMPLEX_SUPPORTED
//
//template<>
//inline BLASWrapper::BLASIndexType BLASWrapper::iamax<ComplexFloat>(const BLASIndexType n,
//		const ComplexFloat *x, const BLASIndexType incX) {
//	return cblas_icamax(n, x,
//			incX);
//}
//
//template<>
//inline BLASWrapper::BLASIndexType BLASWrapper::iamax<ComplexDouble>(const BLASIndexType n,
//		const ComplexDouble *x, const BLASIndexType incX) {
//	return cblas_izamax(n, x,
//			incX);
//}
//
//#endif
//
//// -------------- swap --------------
//template<>
//inline void BLASWrapper::swap<float>(const BLASIndexType n, float *x,
//		const BLASIndexType incX, float *y, const BLASIndexType incY) {
//	cblas_sswap(n, x,
//			incX, y,
//			incY);
//}
//
//template<>
//inline void BLASWrapper::swap<double>(const BLASIndexType n, double *x,
//		const BLASIndexType incX, double *y, const BLASIndexType incY) {
//	cblas_dswap(n, x,
//			incX, y,
//			incY);
//}
//
//#ifdef SCAI_COMPLEX_SUPPORTED
//
//template<>
//inline void BLASWrapper::swap<ComplexFloat>(const BLASIndexType n, ComplexFloat *x,
//		const BLASIndexType incX, ComplexFloat *y, const BLASIndexType incY) {
//	cblas_cswap(n, x,
//			incX, y,
//			incY);
//}
//
//template<>
//inline void BLASWrapper::swap<ComplexDouble>(const BLASIndexType n,
//		ComplexDouble *x, const BLASIndexType incX, ComplexDouble *y,
//		const BLASIndexType incY) {
//	cblas_zswap(n, x,
//			incX, y,
//			incY);
//}
//
//#endif
//
//// -------------- copy --------------
//template<>
//inline void BLASWrapper::copy<float>(const BLASIndexType n, const float *x,
//		const BLASIndexType incX, float *y, const BLASIndexType incY) {
//	cblas_scopy(n, x,
//			incX, y,
//			incY);
//}
//
//template<>
//inline void BLASWrapper::copy<double>(const BLASIndexType n, const double *x,
//		const BLASIndexType incX, double *y, const BLASIndexType incY) {
//	cblas_dcopy(n, x,
//			incX, y,
//			incY);
//}
//
//#ifdef SCAI_COMPLEX_SUPPORTED
//
//template<>
//inline void BLASWrapper::copy<ComplexFloat>(const BLASIndexType n,
//		const ComplexFloat *x, const BLASIndexType incX, ComplexFloat *y,
//		const BLASIndexType incY) {
//	cblas_ccopy(n, x,
//			incX, y,
//			incY);
//}
//
//template<>
//inline void BLASWrapper::copy<ComplexDouble>(const BLASIndexType n,
//		const ComplexDouble *x, const BLASIndexType incX, ComplexDouble *y,
//		const BLASIndexType incY) {
//	cblas_zcopy(n, x,
//			incX, y,
//			incY);
//}
//
//#endif
//
//// -------------- axpy --------------
//template<>
//inline void BLASWrapper::axpy<float>(const BLASIndexType n, const float alpha,
//		const float *x, const BLASIndexType incX, float *y, const BLASIndexType incY) {
//	cblas_saxpy(n, alpha, x,
//			incX, y,
//			incY);
//}
//
//template<>
//inline void BLASWrapper::axpy<double>(const BLASIndexType n, const double alpha,
//		const double *x, const BLASIndexType incX, double *y,
//		const BLASIndexType incY) {
//	cblas_daxpy(n, alpha, x,
//			incX, y,
//			incY);
//}
//
//#ifdef SCAI_COMPLEX_SUPPORTED
//
//template<>
//inline void BLASWrapper::axpy<ComplexFloat>(const BLASIndexType n,
//		const ComplexFloat alpha, const ComplexFloat *x, const BLASIndexType incX,
//		ComplexFloat *y, const BLASIndexType incY) {
//	// Attention: alpha is here passed by a pointer
//	cblas_caxpy(n, &alpha, x,
//			incX, y,
//			incY);
//}
//
//template<>
//inline void BLASWrapper::axpy<ComplexDouble>(const BLASIndexType n,
//		const ComplexDouble alpha, const ComplexDouble *x, const BLASIndexType incX,
//		ComplexDouble *y, const BLASIndexType incY) {
//	// Attention: alpha is here passed by a pointer
//	cblas_zaxpy(n, &alpha, x,
//			incX, y,
//			incY);
//}
//
//#endif
//
//// -------------- dot --------------
//template<>
//inline float BLASWrapper::dot<float>(const BLASIndexType n, const float *x,
//		const BLASIndexType incX, const float *y, const BLASIndexType incY) {
//	return cblas_sdot(n, x,
//			incX, y,
//			incY);
//}
//
//template<>
//inline double BLASWrapper::dot<double>(const BLASIndexType n, const double *x,
//		const BLASIndexType incX, const double *y, const BLASIndexType incY) {
//	return cblas_ddot(n, x,
//			incX, y,
//			incY);
//}
//
//#ifdef SCAI_COMPLEX_SUPPORTED
//
//template<>
//inline ComplexFloat BLASWrapper::dot<ComplexFloat>(const BLASIndexType n,
//		const ComplexFloat *x, const BLASIndexType incX, const ComplexFloat *y,
//		const BLASIndexType incY) {
//	ComplexFloat dotu;
//	cblas_cdotc_sub(n, x,
//			incX, y,
//			incY, &dotu);
//	return dotu;
//}
//
//template<>
//inline ComplexDouble BLASWrapper::dot<ComplexDouble>(const BLASIndexType n,
//		const ComplexDouble *x, const BLASIndexType incX, const ComplexDouble *y,
//		const BLASIndexType incY) {
//	ComplexDouble dotu;
//	cblas_zdotc_sub(n, x,
//			incX, y,
//			incY, &dotu);
//	return dotu;
//}
//
//#endif
//
//// -------------- gemv --------------
//template<>
//inline void BLASWrapper::gemv<float>(const CBLAS_ORDER order,
//		const CBLAS_TRANSPOSE transA, const BLASIndexType m, const BLASIndexType n,
//		const float alpha, const float* A, const BLASIndexType lda, const float* x,
//		const BLASIndexType incX, const float beta, float* y,
//		const BLASIndexType incY) {
//	cblas_sgemv(order, transA, m, n, alpha, A, lda,
//			x, incX, beta, y,
//			incY);
//}
//
//template<>
//inline void BLASWrapper::gemv<double>(const CBLAS_ORDER order,
//		const CBLAS_TRANSPOSE transA, const BLASIndexType m, const BLASIndexType n,
//		const double alpha, const double* A, const BLASIndexType lda,
//		const double* x, const BLASIndexType incX, const double beta, double* y,
//		const BLASIndexType incY) {
//	cblas_dgemv(order, transA, m, n, alpha, A, lda,
//			x, incX, beta, y,
//			incY);
//}
//
//#ifdef SCAI_COMPLEX_SUPPORTED
//
//template<>
//inline void BLASWrapper::gemv<ComplexFloat>(const CBLAS_ORDER order,
//		const CBLAS_TRANSPOSE transA, const BLASIndexType m, const BLASIndexType n,
//		const ComplexFloat alpha, const ComplexFloat* A, const BLASIndexType lda,
//		const ComplexFloat* x, const BLASIndexType incX, const ComplexFloat beta,
//		ComplexFloat* y, const BLASIndexType incY) {
//	// Attention: alpha, beta must be passed here as a pointer
//	cblas_cgemv(order, transA, m, n, &alpha, A, lda,
//			x, incX, &beta, y,
//			incY);
//}
//
//template<>
//inline void BLASWrapper::gemv<ComplexDouble>(const CBLAS_ORDER order,
//		const CBLAS_TRANSPOSE transA, const BLASIndexType m, const BLASIndexType n,
//		const ComplexDouble alpha, const ComplexDouble* A, const BLASIndexType lda,
//		const ComplexDouble* x, const BLASIndexType incX, const ComplexDouble beta,
//		ComplexDouble* y, const BLASIndexType incY) {
//	// Attention: alpha, beta must be passed here as a pointer
//	cblas_zgemv(order, transA, m, n, &alpha, A, lda,
//			x, incX, &beta, y,
//			incY);
//}
//
//#endif
//
//// -------------- gemm --------------
//template<>
//inline void BLASWrapper::gemm<float>(const CBLAS_ORDER order,
//		const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB,
//		const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,
//		const float alpha, const float* A, const BLASIndexType lda, const float* B,
//		const BLASIndexType ldb, const float beta, float* C, const BLASIndexType ldc) {
//	cblas_sgemm(order, transA, transB, m, n, k,
//			alpha, A, lda, B, ldb, beta, C, ldc);
//}
//
//template<>
//inline void BLASWrapper::gemm<double>(const CBLAS_ORDER order,
//		const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB,
//		const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,
//		const double alpha, const double* A, const BLASIndexType lda,
//		const double* B, const BLASIndexType ldb, const double beta, double* C,
//		const BLASIndexType ldc) {
//	cblas_dgemm(order, transA, transB, m, n, k,
//			alpha, A, lda, B, ldb, beta, C, ldc);
//}
//
//#ifdef SCAI_COMPLEX_SUPPORTED
//
//template<>
//inline void BLASWrapper::gemm<ComplexFloat>(const CBLAS_ORDER order,
//		const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB,
//		const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,
//		const ComplexFloat alpha, const ComplexFloat* A, const BLASIndexType lda,
//		const ComplexFloat* B, const BLASIndexType ldb, const ComplexFloat beta,
//		ComplexFloat* C, const BLASIndexType ldc) {
//	// Attention: alpha and beta are passed by a pointer
//	cblas_cgemm(order, transA, transB, m, n, k,
//			&alpha, A, lda, B, ldb, &beta, C, ldc);
//}
//
//template<>
//inline void BLASWrapper::gemm<ComplexDouble>(const CBLAS_ORDER order,
//		const CBLAS_TRANSPOSE transA, const CBLAS_TRANSPOSE transB,
//		const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,
//		const ComplexDouble alpha, const ComplexDouble* A, const BLASIndexType lda,
//		const ComplexDouble* B, const BLASIndexType ldb, const ComplexDouble beta,
//		ComplexDouble* C, const BLASIndexType ldc) {
//	// Attention: alpha and beta are passed by a pointer
//	cblas_zgemm(order, transA, transB, m, n, k,
//			&alpha, A, lda, B, ldb, &beta, C, ldc);
//}
//
//#endif

} /* end namespace blaskernel */

} /* end namespace scai */

