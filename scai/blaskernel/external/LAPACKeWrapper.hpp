/**
 * @file LAPACKeWrapper.hpp
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
 * @brief Wrapper for LAPACKe functions
 * @author Eric Schricker
 * @date 12.11.2015
 * @since 2.0.0
 */

#pragma once

#include <scai/blaskernel/MKLUtils.hpp>

#include <scai/common/exception/NotSupportedValueTypeException.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/SCAITypes.hpp>

#include <mkl_lapack.h>
#include <mkl_lapacke.h>

namespace scai {

namespace blaskernel {

class COMMON_DLL_IMPORTEXPORT LAPACKeWrapper {
public:
	typedef lapack_int LAPACKIndexType;
	typedef char LAPACKFlag;
	typedef int LAPACKOrder;

	template<typename ValueType>
	static LAPACKIndexType getrf(const LAPACKOrder UNUSED(matrix_order),
			const LAPACKIndexType UNUSED(m), const LAPACKIndexType UNUSED(n),
			ValueType* const UNUSED(a), const LAPACKIndexType UNUSED(lda),
			LAPACKIndexType* const UNUSED(ipiv)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "getrf");
	}

	template<typename ValueType>
	static LAPACKIndexType getri(const LAPACKOrder UNUSED( matrix_order ),
			const LAPACKIndexType UNUSED(n), ValueType* const UNUSED(A),
			const LAPACKIndexType UNUSED(lda),
			LAPACKIndexType* const UNUSED(ipiv)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "getri");
	}

	template<typename ValueType>
	static LAPACKIndexType tptrs(const LAPACKOrder UNUSED( matrix_order ),
			const LAPACKFlag UNUSED(uplo), const LAPACKFlag UNUSED(trans),
			const LAPACKFlag UNUSED(diag), const LAPACKIndexType UNUSED(n),
			const LAPACKIndexType UNUSED(nrhs), const ValueType* UNUSED(AP),
			ValueType* UNUSED(B), const LAPACKIndexType UNUSED(ldb)) {
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "tptrs");
	}
};

template<>
LAPACKeWrapper::LAPACKIndexType LAPACKeWrapper::getrf<float>(
		const LAPACKOrder matrix_order, const LAPACKIndexType m,
		const LAPACKIndexType n, float* const a, const LAPACKIndexType lda,
		LAPACKIndexType* const ipiv) {
	return LAPACKE_sgetrf(matrix_order, m, n, a, lda, ipiv);
}

template<>
LAPACKeWrapper::LAPACKIndexType LAPACKeWrapper::getrf<double>(
		const LAPACKOrder matrix_order, const LAPACKIndexType m,
		const LAPACKIndexType n, double* const a, const LAPACKIndexType lda,
		LAPACKIndexType* const ipiv) {
	return LAPACKE_dgetrf(matrix_order, m, n, a, lda, ipiv);
}

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
LAPACKeWrapper::LAPACKIndexType LAPACKeWrapper::getrf<ComplexFloat>(
		const LAPACKOrder matrix_order, const LAPACKIndexType m,
		const LAPACKIndexType n, ComplexFloat* const a,
		const LAPACKIndexType lda, LAPACKIndexType* const ipiv) {
	return LAPACKE_cgetrf(matrix_order, m, n,
			MKLUtils::cast(a), lda, ipiv);
}

template<>
LAPACKeWrapper::LAPACKIndexType LAPACKeWrapper::getrf<ComplexDouble>(
		const LAPACKOrder matrix_order, const LAPACKIndexType m,
		const LAPACKIndexType n, ComplexDouble* const a,
		const LAPACKIndexType lda, LAPACKIndexType* const ipiv) {
	return LAPACKE_zgetrf(matrix_order, m, n,
			MKLUtils::cast(a), lda, ipiv);
}

#endif

template<>
LAPACKeWrapper::LAPACKIndexType LAPACKeWrapper::getri<float>(
		const LAPACKOrder matrix_order, const LAPACKIndexType n, float* const A,
		const LAPACKIndexType lda, LAPACKIndexType* const ipiv) {
	return LAPACKE_sgetri(matrix_order, n, A, lda, ipiv);
}

template<>
LAPACKeWrapper::LAPACKIndexType LAPACKeWrapper::getri<double>(
		const LAPACKOrder matrix_order, const LAPACKIndexType n,
		double* const A, const LAPACKIndexType lda,
		LAPACKIndexType* const ipiv) {
	return LAPACKE_dgetri(matrix_order, n, A, lda, ipiv);
}

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
LAPACKeWrapper::LAPACKIndexType LAPACKeWrapper::getri<ComplexFloat>(
		const LAPACKOrder matrix_order, const LAPACKIndexType n,
		ComplexFloat* const A, const LAPACKIndexType lda,
		LAPACKIndexType* const ipiv) {
	return LAPACKE_cgetri(matrix_order, n,
			MKLUtils::cast(A), lda, ipiv);
}

template<>
LAPACKeWrapper::LAPACKIndexType LAPACKeWrapper::getri<ComplexDouble>(
		const LAPACKOrder matrix_order, const LAPACKIndexType n,
		ComplexDouble* const A, const LAPACKIndexType lda,
		LAPACKIndexType* const ipiv) {
	return LAPACKE_zgetri(matrix_order, n,
			MKLUtils::cast(A), lda, ipiv);
}

#endif

template<>
LAPACKeWrapper::LAPACKIndexType LAPACKeWrapper::tptrs<float>(
		const LAPACKOrder matrix_order, const LAPACKFlag uplo,
		const LAPACKFlag trans, const LAPACKFlag diag, const LAPACKIndexType n,
		const LAPACKIndexType nrhs, const float* AP, float* B,
		const LAPACKIndexType ldb) {
	return LAPACKE_stptrs(matrix_order, uplo, trans, diag, n, nrhs, AP, B, ldb);
}

template<>
LAPACKeWrapper::LAPACKIndexType LAPACKeWrapper::tptrs<double>(
		const LAPACKOrder matrix_order, const LAPACKFlag uplo,
		const LAPACKFlag trans, const LAPACKFlag diag, const LAPACKIndexType n,
		const LAPACKIndexType nrhs, const double* AP, double* B,
		const LAPACKIndexType ldb) {
	return LAPACKE_dtptrs(matrix_order, uplo, trans, diag, n, nrhs, AP, B, ldb);
}

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
LAPACKeWrapper::LAPACKIndexType LAPACKeWrapper::tptrs<ComplexFloat>(
		const LAPACKOrder matrix_order, const LAPACKFlag uplo,
		const LAPACKFlag trans, const LAPACKFlag diag, const LAPACKIndexType n,
		const LAPACKIndexType nrhs, const ComplexFloat* AP, ComplexFloat* B,
		const LAPACKIndexType ldb) {
	return LAPACKE_ctptrs(matrix_order, uplo, trans, diag, n, nrhs,
			MKLUtils::cast(AP),
			MKLUtils::cast(B), ldb);
}

template<>
LAPACKeWrapper::LAPACKIndexType LAPACKeWrapper::tptrs<ComplexDouble>(
		const LAPACKOrder matrix_order, const LAPACKFlag uplo,
		const LAPACKFlag trans, const LAPACKFlag diag, const LAPACKIndexType n,
		const LAPACKIndexType nrhs, const ComplexDouble* AP, ComplexDouble* B,
		const LAPACKIndexType ldb) {
	return LAPACKE_ztptrs(matrix_order, uplo, trans, diag, n, nrhs,
			MKLUtils::cast(AP),
			MKLUtils::cast(B), ldb);
}

#endif

} /* end namespace blaskernel */

} /* end namespace scai */

