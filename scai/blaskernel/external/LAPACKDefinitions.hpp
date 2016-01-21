/**
 * @file LAPACKDefinitions.hpp
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
 * @brief Definitions of LAPACK functions
 * @author Eric Schricker
 * @date 14.01.2016
 * @since 2.0.0
 */

#pragma once

// internal scai libraries
#include <scai/common/config.hpp>

//macros
#define FORTRAN_LAPACK_NAME( name, prefix ) prefix##name##_

#define FORTRAN_LAPACK_DEF( name, prefix, retType, definition ) 			\
        retType FORTRAN_LAPACK_NAME( name, prefix )( definition );

namespace scai {

namespace blaskernel {

class COMMON_DLL_IMPORTEXPORT LAPACKDefinitions
{
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
};

} /* end namespace blaskernel */

} /* end namespace scai */

extern "C"
{
#define CALL_DEF_GETRI( ValueType ) const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* n, ValueType* a, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* lda, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* ipivot, ValueType* work, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* ldwork, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* info
#define CALL_DEF_GETRF( ValueType ) const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* m,	const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* n, ValueType* a,const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* lda, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* ipivot, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* info
#define CALL_DEF_TPTRS( ValueType ) const scai::blaskernel::LAPACKDefinitions::LAPACKFlag* uplo, const scai::blaskernel::LAPACKDefinitions::LAPACKFlag* transa, const scai::blaskernel::LAPACKDefinitions::LAPACKFlag* diag, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* n, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* nrhs, const ValueType* ap, ValueType* b, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* ldb, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* info
#define CALL_DEF_LASWP( ValueType ) const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* n, ValueType* a, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* lda, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* k1, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* k2, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* ipiv, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* incx

// float
FORTRAN_LAPACK_DEF( getri, s, void, CALL_DEF_GETRI( float ) )
FORTRAN_LAPACK_DEF( getrf, s, void, CALL_DEF_GETRF( float ) )
FORTRAN_LAPACK_DEF( tptrs, s, void, CALL_DEF_TPTRS( float ) )
FORTRAN_LAPACK_DEF( laswp, s, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType, CALL_DEF_LASWP( float ))
// double
FORTRAN_LAPACK_DEF( getrf, d, void, CALL_DEF_GETRF( double ) )
FORTRAN_LAPACK_DEF( getri, d, void, CALL_DEF_GETRI( double ) )
FORTRAN_LAPACK_DEF( tptrs, d, void, CALL_DEF_TPTRS( double ) )
FORTRAN_LAPACK_DEF( laswp, d, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType, CALL_DEF_LASWP( double ))

#ifdef SCAI_COMPLEX_SUPPORTED
// ComplexFloat
FORTRAN_LAPACK_DEF( getrf, c, void, CALL_DEF_GETRF( ComplexFloat ) )
FORTRAN_LAPACK_DEF( getri, c, void, CALL_DEF_GETRI( ComplexFloat ) )
FORTRAN_LAPACK_DEF( tptrs, c, void, CALL_DEF_TPTRS( ComplexFloat ) )
FORTRAN_LAPACK_DEF( laswp, c, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType, CALL_DEF_LASWP( ComplexFloat ))
// ComplexDouble
FORTRAN_LAPACK_DEF( getri, z, void, CALL_DEF_GETRI( ComplexDouble ) )
FORTRAN_LAPACK_DEF( getrf, z, void, CALL_DEF_GETRF( ComplexDouble ) )
FORTRAN_LAPACK_DEF( tptrs, z, void, CALL_DEF_TPTRS( ComplexDouble ) )
FORTRAN_LAPACK_DEF( laswp, z, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType, CALL_DEF_LASWP( ComplexDouble ))
#endif

#undef CALL_DEF_GETRI
#undef CALL_DEF_GETRF
#undef CALL_DEF_TPTRS
#undef CALL_DEF_LASWP

#undef FORTRAN_LAPACK_DEF

} /* end extern "C" */
