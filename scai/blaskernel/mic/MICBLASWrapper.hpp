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

#include <scai/blaskernel/MKLUtils.hpp>

#include <scai/common/exception/NotSupportedValueTypeException.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/SCAITypes.hpp>

#include <mkl.h>

namespace scai {

namespace blaskernel {

class COMMON_DLL_IMPORTEXPORT MICBLASWrapper {
public:
	typedef int BLASIndexType;
	typedef char BLASTranspose;

	// ------------- BLAS2 -------------
	template<typename ValueType>
	__declspec( target(mic) )
	static void gemv(const BLASTranspose UNUSED(transA_BLASTranspose),
			const BLASIndexType UNUSED(m), const BLASIndexType UNUSED(n),
			const ValueType UNUSED(alpha), const ValueType* UNUSED(A),
			const BLASIndexType UNUSED(lda), const ValueType* UNUSED(x),
			const BLASIndexType UNUSED(incX), const ValueType UNUSED(beta),
			ValueType* UNUSED(y), const BLASIndexType UNUSED(incY)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "gemv");
	}

	// ------------- BLAS3 -------------
	template<typename ValueType>
	__declspec( target(mic) ) 
	static void gemm(const BLASTranspose UNUSED(transA_BLASTranspose),
			const BLASTranspose UNUSED(transB_BLASTranspose),
			const BLASIndexType UNUSED(m), const BLASIndexType UNUSED(n),
			const BLASIndexType UNUSED(k), const ValueType UNUSED(alpha),
			const ValueType* UNUSED(A), const BLASIndexType UNUSED(lda),
			const ValueType* UNUSED(B), const BLASIndexType UNUSED(ldb),
			const ValueType UNUSED(beta), ValueType* UNUSED(C),
			const BLASIndexType UNUSED(ldc)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "gemm");
	}
};

// -------------- gemv --------------
template<>
inline void MICBLASWrapper::gemv<float>(const BLASTranspose transA, const BLASIndexType m, const BLASIndexType n,
		const float alpha, const float* A, const BLASIndexType lda, const float* x,
		const BLASIndexType incX, const float beta, float* y,
		const BLASIndexType incY) {
	sgemv( &transA, &m, &n, &alpha, A, &lda,
			x, &incX, &beta, y,
			&incY);
}

template<>
inline void MICBLASWrapper::gemv<double>(const BLASTranspose transA, const BLASIndexType m, const BLASIndexType n,
		const double alpha, const double* A, const BLASIndexType lda,
		const double* x, const BLASIndexType incX, const double beta, double* y,
		const BLASIndexType incY) {
	dgemv( &transA, &m, &n, &alpha, A, &lda,
			x, &incX, &beta, y,
			&incY);
}

template<>
inline void MICBLASWrapper::gemv<ComplexFloat>(const BLASTranspose transA, const BLASIndexType m, const BLASIndexType n,
		const ComplexFloat alpha, const ComplexFloat* A, const BLASIndexType lda,
		const ComplexFloat* x, const BLASIndexType incX, const ComplexFloat beta,
		ComplexFloat* y, const BLASIndexType incY) {
	cgemv( &transA, &m, &n, MKLUtils::cast(&alpha), MKLUtils::cast(A), &lda,
			MKLUtils::cast(x), &incX, MKLUtils::cast(&beta), MKLUtils::cast(y),
			&incY);
}

template<>
inline void MICBLASWrapper::gemv<ComplexDouble>(const BLASTranspose transA, const BLASIndexType m, const BLASIndexType n,
		const ComplexDouble alpha, const ComplexDouble* A, const BLASIndexType lda,
		const ComplexDouble* x, const BLASIndexType incX, const ComplexDouble beta,
		ComplexDouble* y, const BLASIndexType incY) {
	zgemv( &transA, &m, &n, MKLUtils::cast(&alpha), MKLUtils::cast(A), &lda,
			MKLUtils::cast(x), &incX, MKLUtils::cast(&beta), MKLUtils::cast(y),
			&incY);
}

// -------------- gemm --------------
template<>
inline void MICBLASWrapper::gemm<float>(const BLASTranspose transA, const BLASTranspose transB,
		const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,
		const float alpha, const float* A, const BLASIndexType lda, const float* B,
		const BLASIndexType ldb, const float beta, float* C, const BLASIndexType ldc) {
	sgemm(&transA, &transB, &m, &n, &k,
			&alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

template<>
inline void MICBLASWrapper::gemm<double>(const BLASTranspose transA, const BLASTranspose transB,
		const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,
		const double alpha, const double* A, const BLASIndexType lda,
		const double* B, const BLASIndexType ldb, const double beta, double* C,
		const BLASIndexType ldc) {
	dgemm(&transA, &transB, &m, &n, &k,
			&alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

template<>
inline void MICBLASWrapper::gemm<ComplexFloat>(const BLASTranspose transA, const BLASTranspose transB,
		const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,
		const ComplexFloat alpha, const ComplexFloat* A, const BLASIndexType lda,
		const ComplexFloat* B, const BLASIndexType ldb, const ComplexFloat beta,
		ComplexFloat* C, const BLASIndexType ldc) {
	cgemm(&transA, &transB, &m, &n, &k,
			MKLUtils::cast(&alpha), MKLUtils::cast(A), &lda, MKLUtils::cast(B), &ldb, MKLUtils::cast(&beta), MKLUtils::cast(C), &ldc);
}

template<>
inline void MICBLASWrapper::gemm<ComplexDouble>(const BLASTranspose transA, const BLASTranspose transB,
		const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,
		const ComplexDouble alpha, const ComplexDouble* A, const BLASIndexType lda,
		const ComplexDouble* B, const BLASIndexType ldb, const ComplexDouble beta,
		ComplexDouble* C, const BLASIndexType ldc) {
	zgemm(&transA, &transB, &m, &n, &k,
			MKLUtils::cast(&alpha), MKLUtils::cast(A), &lda, MKLUtils::cast(B), &ldb, MKLUtils::cast(&beta), MKLUtils::cast(C), &ldc);
}

} /* end namespace blaskernel */

} /* end namespace scai */

