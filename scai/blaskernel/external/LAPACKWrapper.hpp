/**
 * @file LAPACKWrapper.hpp
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
 * @brief Wrapper for LAPACK functions
 * @author Eric Schricker
 * @date 12.11.2015
 * @since 2.0.0
 */

#pragma once

#include <scai/common/exception/NotSupportedValueTypeException.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai {

namespace blaskernel {

class COMMON_DLL_IMPORTEXPORT LAPACKWrapperOLD {
public:
#ifdef F77_INT
	typedef F77_INT LAPACKIndexType;
#else
	typedef int LAPACKIndexType;
#endif

#ifdef F77_CHAR
	typedef F77_CHAR LAPACKFlag;
#else
	typedef char LAPACKFlag;
#endif

	template<typename ValueType>
	static LAPACKIndexType getrf(
			const LAPACKIndexType UNUSED(m),
			const LAPACKIndexType UNUSED(n), ValueType* UNUSED(a),
			const LAPACKIndexType UNUSED(lda),
			LAPACKIndexType* UNUSED(ipivot))
	{
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "getrf");
	}

	template<typename ValueType>
	static LAPACKIndexType getri(
			const LAPACKIndexType UNUSED(n), ValueType* UNUSED(a),
			const LAPACKIndexType UNUSED(lda),
			LAPACKIndexType* UNUSED(ipivot), ValueType* UNUSED(work),
			const LAPACKIndexType UNUSED(ldwork))
	{
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "getri");
	}

	template<typename ValueType>
	static LAPACKIndexType tptrs(LAPACKFlag UNUSED(uplo),
			LAPACKFlag transa, LAPACKFlag UNUSED(diag),
			const LAPACKIndexType UNUSED(n),
			const LAPACKIndexType UNUSED(nrhs), const ValueType* UNUSED(ap),
			ValueType* UNUSED(b), const LAPACKIndexType UNUSED(ldb))
	{
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "tptrs");
	}

	template<typename ValueType>
	static void laswp(
			const LAPACKIndexType n,
			ValueType* a,
			const LAPACKIndexType lda,
			const LAPACKIndexType k1,
			const LAPACKIndexType k2,
			const LAPACKIndexType* ipiv,
			const LAPACKIndexType incx)
	{
		SCAI_THROWEXCEPTION(common::NotSupportedValueTypeException, "laswp");
	}

};

} /* end namespace blaskernel */

} /* end namespace scai */

//fallback if nothing is set in cmake
#if !defined(LAMA_FORTRAN_BLAS_STYLE_UNDERSCORE)
#if !defined(LAMA_FORTRAN_BLAS_STYLE_UPCASE)
#if !defined(LAMA_FORTRAN_BLAS_STYLE_LOWCASE)
#define LAMA_FORTRAN_BLAS_STYLE_LOWCASE
#endif
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif /*__cplusplus*/

#define F77_sgetrf sgetrf_
#define F77_dgetrf dgetrf_
#define F77_cgetrf cgetrf_
#define F77_zgetrf zgetrf_
#define F77_sgetri sgetri_
#define F77_dgetri dgetri_
#define F77_cgetri cgetri_
#define F77_zgetri zgetri_
#define F77_strtrs strtrs_
#define F77_dtrtrs dtrtrs_
#define F77_ctrtrs ctrtrs_
#define F77_ztrtrs ztrtrs_
#define F77_stptrs stptrs_
#define F77_dtptrs dtptrs_
#define F77_ctptrs ctptrs_
#define F77_ztptrs ztptrs_
#define F77_slaswp slaswp_
#define F77_dlaswp dlaswp_
#define F77_claswp claswp_
#define F77_zlaswp zlaswp_

void F77_sgetrf(const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* m,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n, float* a,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* lda,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ipivot,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* info);
void F77_dgetrf(const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* m,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n, double* a,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* lda,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ipivot,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* info);
void F77_cgetrf(const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* m,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n, ComplexFloat* a,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* lda,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ipivot,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* info);
void F77_zgetrf(const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* m,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n, ComplexDouble* a,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* lda,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ipivot,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* info);
void F77_sgetri(const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n, float* a,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* lda,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ipivot, float* work,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ldwork,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* info);
void F77_dgetri(const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n, double* a,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* lda,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ipivot, double* work,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ldwork,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* info);
void F77_cgetri(const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n, ComplexFloat* a,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* lda,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ipivot, ComplexFloat* work,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ldwork,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* info);
void F77_zgetri(const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n, ComplexDouble* a,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* lda,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ipivot, ComplexDouble* work,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ldwork,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* info);
void F77_stptrs(char* uplo, char* transa, char* diag,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* nrhs, const float* ap, float* b,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ldb,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* info);
void F77_dtptrs(char* uplo, char* transa, char* diag,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* nrhs, const double* ap, double* b,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ldb,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* info);
void F77_ctptrs(char* uplo, char* transa, char* diag,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* nrhs, const ComplexFloat* ap,
		ComplexFloat* b, const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ldb,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* info);
void F77_ztptrs(char* uplo, char* transa, char* diag,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* nrhs, const ComplexDouble* ap,
		ComplexDouble* b, const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ldb,
		scai::blaskernel::LAPACKWrapper::LAPACKIndexType* info);
scai::blaskernel::LAPACKWrapper::LAPACKIndexType F77_slaswp(
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n, float* a,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* lda,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* k1,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* k2,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ipiv,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* incx);
scai::blaskernel::LAPACKWrapper::LAPACKIndexType F77_dlaswp(
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n, double* a,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* lda,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* k1,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* k2,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ipiv,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* incx);
scai::blaskernel::LAPACKWrapper::LAPACKIndexType F77_claswp(
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n, ComplexFloat* a,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* lda,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* k1,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* k2,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ipiv,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* incx);
scai::blaskernel::LAPACKWrapper::LAPACKIndexType F77_zlaswp(
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* n, ComplexDouble* a,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* lda,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* k1,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* k2,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* ipiv,
		const scai::blaskernel::LAPACKWrapper::LAPACKIndexType* incx);

#ifdef __cplusplus
} /*extern "C"*/
#endif /*__cplusplus*/

namespace scai {

namespace blaskernel {

template<>
LAPACKWrapper::LAPACKIndexType LAPACKWrapper::getrf<float>(
		const LAPACKIndexType m,
		const LAPACKIndexType n, float* a,
		const LAPACKIndexType lda,
		LAPACKIndexType* ipivot) {
	LAPACKIndexType info;

	F77_sgetrf(&m, &n, a, &lda, ipivot, &info);

	return info;
}

template<>
LAPACKWrapper::LAPACKIndexType LAPACKWrapper::getrf<double>(
		const LAPACKIndexType m,
		const LAPACKIndexType n, double* a,
		const LAPACKIndexType lda,
		LAPACKIndexType* ipivot) {
	LAPACKIndexType info;

	F77_dgetrf(&m, &n, a, &lda, ipivot, &info);

	return info;
}

template<>
LAPACKWrapper::LAPACKIndexType LAPACKWrapper::getrf<ComplexFloat>(
		const LAPACKIndexType m,
		const LAPACKIndexType n, ComplexFloat* a,
		const LAPACKIndexType lda,
		LAPACKIndexType* ipivot) {
	LAPACKIndexType info;

	F77_cgetrf(&m, &n, a, &lda, ipivot, &info);

	return info;
}

template<>
LAPACKWrapper::LAPACKIndexType LAPACKWrapper::getrf<ComplexDouble>(
		const LAPACKIndexType m,
		const LAPACKIndexType n, ComplexDouble* a,
		const LAPACKIndexType lda,
		LAPACKIndexType* ipivot) {
	LAPACKIndexType info;

	F77_zgetrf(&m, &n, a, &lda, ipivot, &info);

	return info;
}

template<>
LAPACKWrapper::LAPACKIndexType LAPACKWrapper::getri<float>(
		const LAPACKIndexType n, float* a,
		const LAPACKIndexType lda,
		LAPACKIndexType* ipivot, float* work,
		const LAPACKIndexType ldwork) {
	LAPACKIndexType info;

	F77_sgetri(&n, a, &lda, ipivot, work, &ldwork, &info);

	return info;
}

template<>
LAPACKWrapper::LAPACKIndexType LAPACKWrapper::getri<double>(
		const LAPACKIndexType n, double* a,
		const LAPACKIndexType lda,
		LAPACKIndexType* ipivot, double* work,
		const LAPACKIndexType ldwork) {
	LAPACKIndexType info;

	F77_dgetri(&n, a, &lda, ipivot, work, &ldwork, &info);

	return info;
}

template<>
LAPACKWrapper::LAPACKIndexType LAPACKWrapper::getri<ComplexFloat>(
		const LAPACKIndexType n, ComplexFloat* a,
		const LAPACKIndexType lda,
		LAPACKIndexType* ipivot, ComplexFloat* work,
		const LAPACKIndexType ldwork) {
	LAPACKIndexType info;

	F77_cgetri(&n, a, &lda, ipivot, work, &ldwork, &info);

	return info;
}

template<>
LAPACKWrapper::LAPACKIndexType LAPACKWrapper::getri<ComplexDouble>(
		const LAPACKIndexType n, ComplexDouble* a,
		const LAPACKIndexType lda,
		LAPACKIndexType* ipivot, ComplexDouble* work,
		const LAPACKIndexType ldwork) {
	LAPACKIndexType info;

	F77_zgetri(&n, a, &lda, ipivot, work, &ldwork, &info);

	return info;
}

template<>
LAPACKWrapper::LAPACKIndexType LAPACKWrapper::tptrs<float>(
		LAPACKFlag uplo, LAPACKFlag transa,
		LAPACKFlag diag,
		const LAPACKIndexType n,
		const LAPACKIndexType nrhs, const float* ap, float* b,
		const LAPACKIndexType ldb)
{
	LAPACKIndexType info;

	F77_stptrs(&uplo, &transa, &diag, &n, &nrhs, ap, b, &ldb, &info );

	return info;
}

template<>
LAPACKWrapper::LAPACKIndexType LAPACKWrapper::tptrs<double>(
		LAPACKFlag uplo, LAPACKFlag transa,
		LAPACKFlag diag,
		const LAPACKIndexType n,
		const LAPACKIndexType nrhs, const double* ap, double* b,
		const LAPACKIndexType ldb) {
	LAPACKIndexType info;

	F77_dtptrs(&uplo, &transa, &diag, &n, &nrhs, ap, b, &ldb, &info );

	return info;
}

template<>
LAPACKWrapper::LAPACKIndexType LAPACKWrapper::tptrs<ComplexFloat>(
		LAPACKFlag uplo, LAPACKFlag transa,
		LAPACKFlag diag,
		const LAPACKIndexType n,
		const LAPACKIndexType nrhs, const ComplexFloat* ap, ComplexFloat* b,
		const LAPACKIndexType ldb) {
	LAPACKIndexType info;

	F77_ctptrs(&uplo, &transa, &diag, &n, &nrhs, ap, b, &ldb, &info );

	return info;
}

template<>
LAPACKWrapper::LAPACKIndexType LAPACKWrapper::tptrs<ComplexDouble>(
		LAPACKFlag uplo, LAPACKFlag transa,
		LAPACKFlag diag,
		const LAPACKIndexType n,
		const LAPACKIndexType nrhs, const ComplexDouble* ap, ComplexDouble* b,
		const LAPACKIndexType ldb) {
	LAPACKIndexType info;

	F77_ztptrs(&uplo, &transa, &diag, &n, &nrhs, ap, b, &ldb, &info );

	return info;
}

template<>
void LAPACKWrapper::laswp<float>(
		const LAPACKIndexType n,
		float* a,
		const LAPACKIndexType lda,
		const LAPACKIndexType k1,
		const LAPACKIndexType k2,
		const LAPACKIndexType* ipiv,
		const LAPACKIndexType incx)
{
	F77_slaswp( &n, a, &lda, &k1, &k2, ipiv, &incx);
}

template<>
void LAPACKWrapper::laswp<double>(
		const LAPACKIndexType n,
		double* a,
		const LAPACKIndexType lda,
		const LAPACKIndexType k1,
		const LAPACKIndexType k2,
		const LAPACKIndexType* ipiv,
		const LAPACKIndexType incx)
{
	F77_dlaswp( &n, a, &lda, &k1, &k2, ipiv, &incx);
}

template<>
void LAPACKWrapper::laswp<ComplexFloat>(
		const LAPACKIndexType n,
		ComplexFloat* a,
		const LAPACKIndexType lda,
		const LAPACKIndexType k1,
		const LAPACKIndexType k2,
		const LAPACKIndexType* ipiv,
		const LAPACKIndexType incx)
{
	F77_claswp( &n, a, &lda, &k1, &k2, ipiv, &incx);
}

template<>
void LAPACKWrapper::laswp<ComplexDouble>(
		const LAPACKIndexType n,
		ComplexDouble* a,
		const LAPACKIndexType lda,
		const LAPACKIndexType k1,
		const LAPACKIndexType k2,
		const LAPACKIndexType* ipiv,
		const LAPACKIndexType incx)
{
	F77_zlaswp( &n, a, &lda, &k1, &k2, ipiv, &incx);
}

} /* end namespace blaskernel */

} /* end namespace scai */

