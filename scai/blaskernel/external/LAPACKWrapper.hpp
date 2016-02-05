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
 * @date 14.01.2016
 * @since 2.0.0
 */

#pragma once

// local library
#include <scai/blaskernel/external/LAPACKTrait.hpp>

// internal scai libraries
#include <scai/common/macros/unused.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai {

namespace blaskernel {

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT LAPACKWrapper;

#define LAPACKWRAPPER_DEF( ValueType, prefix ) 															\
template<>																								\
class COMMON_DLL_IMPORTEXPORT LAPACKWrapper<ValueType>													\
{																										\
public:																									\
	typedef LAPACKTrait::LAPACKIndexType LAPACKIndexType;												\
	typedef LAPACKTrait::LAPACKFlag LAPACKFlag;															\
																										\
	static LAPACKIndexType getrf(																		\
			const LAPACKIndexType m,																	\
			const LAPACKIndexType n, ValueType* a,														\
			const LAPACKIndexType lda,																	\
			LAPACKIndexType* ipivot)																	\
	{																									\
		LAPACKIndexType info;																			\
		FORTRAN_LAPACK_NAME( getrf, prefix )(&m, &n, a, &lda, ipivot, &info);							\
		return info;																					\
	}																									\
																										\
	static LAPACKIndexType getri(																		\
			const LAPACKIndexType n, ValueType* a,														\
			const LAPACKIndexType lda,																	\
			LAPACKIndexType* ipivot, ValueType* work,													\
			const LAPACKIndexType ldwork)																\
	{																									\
		LAPACKIndexType info;																			\
		FORTRAN_LAPACK_NAME( getri, prefix )(&n, a, &lda, ipivot, work, &ldwork, &info);				\
		return info;																					\
	}																									\
																										\
	static LAPACKIndexType tptrs(LAPACKFlag uplo,														\
			LAPACKFlag transa, LAPACKFlag diag,															\
			const LAPACKIndexType n,																	\
			const LAPACKIndexType nrhs, const ValueType* ap,											\
			ValueType* b, const LAPACKIndexType ldb)													\
	{																									\
		LAPACKIndexType info;																			\
		FORTRAN_LAPACK_NAME( tptrs, prefix )(&uplo, &transa, &diag, &n, &nrhs, ap, b, &ldb, &info );	\
		return info;																					\
	}																									\
																										\
	static void laswp(																					\
			const LAPACKIndexType n,																	\
			ValueType* a,																				\
			const LAPACKIndexType lda,																	\
			const LAPACKIndexType k1,																	\
			const LAPACKIndexType k2,																	\
			const LAPACKIndexType* ipiv,																\
			const LAPACKIndexType incx)																	\
	{																									\
		FORTRAN_LAPACK_NAME( laswp, prefix )( &n, a, &lda, &k1, &k2, ipiv, &incx);						\
	}																									\
};


LAPACKWRAPPER_DEF( float, s )
LAPACKWRAPPER_DEF( double, d )

#ifdef SCAI_COMPLEX_SUPPORTED
LAPACKWRAPPER_DEF( ComplexFloat, c )
LAPACKWRAPPER_DEF( ComplexDouble, z )
#endif

#undef LAPACKWRAPPER_DEF

} /* end namespace blaskernel */

} /* end namespace scai */

