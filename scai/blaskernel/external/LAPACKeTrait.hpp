/**
 * @file LAPACKeTrait.hpp
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
 * @brief Typedefs and macros for LAPACKeWrapper
 * @author Eric Schricker
 * @date 14.01.2016
 * @since 2.0.0
 */

#pragma once

// internal scai libraries
#include <scai/common/config.hpp>

// external
#include <mkl_lapacke.h>

// macros
#define FORTRAN_LAPACKE_NAME( name, prefix ) LAPACKE_##prefix##name

namespace scai {

namespace blaskernel {

class COMMON_DLL_IMPORTEXPORT LAPACKeTrait
{
public:
	typedef lapack_int LAPACKIndexType;
	typedef char LAPACKFlag;
	typedef int LAPACKOrder;

	static inline LAPACKFlag enum2char( const CBLAS_UPLO uplo )
	{
		switch( uplo )
		{
			case CblasUpper:
				return 'U';
			case CblasLower:
				return 'L';
			default:
				COMMON_THROWEXCEPTION( "Illegal uplo: " << uplo );
		}
	}

	static inline LAPACKFlag enum2char( const CBLAS_TRANSPOSE trans )
	{
		switch( trans )
		{
			case CblasNoTrans:
				return 'N';
			case CblasTrans:
				return 'T';
			case CblasConjTrans:
				return 'C';
			default:
				COMMON_THROWEXCEPTION( "Illegal trans: " << trans );
		}
	}

	static inline LAPACKFlag enum2char( const CBLAS_DIAG diag )
	{
		switch( diag )
		{
			case CblasNonUnit:
				return 'N';
			case CblasUnit:
				return 'U';
			default:
				COMMON_THROWEXCEPTION( "Illegal diag: " << diag );
		}
	}

	static inline LAPACKOrder enum2order(const CBLAS_ORDER order) {
		switch( order )
		{
			case CblasColMajor:
				return LAPACK_COL_MAJOR;
			case CblasRowMajor:
				return LAPACK_ROW_MAJOR;
			default:
				COMMON_THROWEXCEPTION("illegal matrix order " << order )
		}
	}
};

} /* end namespace blaskernel */

} /* end namespace scai */

