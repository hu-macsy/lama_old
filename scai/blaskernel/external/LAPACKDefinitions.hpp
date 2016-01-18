/*
 * LAPACKDefinitions.hpp
 *
 *  Created on: Jan 14, 2016
 *      Author: eschricker
 */

#pragma once

#include <scai/common/config.hpp>

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

FORTRAN_LAPACK_DEF( getri, s, void, CALL_DEF_GETRI( float ) )
FORTRAN_LAPACK_DEF( getri, d, void, CALL_DEF_GETRI( double ) )
FORTRAN_LAPACK_DEF( getri, c, void, CALL_DEF_GETRI( ComplexFloat ) )
FORTRAN_LAPACK_DEF( getri, z, void, CALL_DEF_GETRI( ComplexDouble ) )

#undef CALL_DEF_GETRI
#define CALL_DEF_GETRF( ValueType ) const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* m,	const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* n, ValueType* a,const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* lda, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* ipivot, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* info

FORTRAN_LAPACK_DEF( getrf, s, void, CALL_DEF_GETRF( float ) )
FORTRAN_LAPACK_DEF( getrf, d, void, CALL_DEF_GETRF( double ) )
FORTRAN_LAPACK_DEF( getrf, c, void, CALL_DEF_GETRF( ComplexFloat ) )
FORTRAN_LAPACK_DEF( getrf, z, void, CALL_DEF_GETRF( ComplexDouble ) )

#undef CALL_DEF_GETRF
#define CALL_DEF_TPTRS( ValueType ) const scai::blaskernel::LAPACKDefinitions::LAPACKFlag* uplo, const scai::blaskernel::LAPACKDefinitions::LAPACKFlag* transa, const scai::blaskernel::LAPACKDefinitions::LAPACKFlag* diag, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* n, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* nrhs, const ValueType* ap, ValueType* b, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* ldb, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* info

FORTRAN_LAPACK_DEF( tptrs, s, void, CALL_DEF_TPTRS( float ) )
FORTRAN_LAPACK_DEF( tptrs, d, void, CALL_DEF_TPTRS( double ) )
FORTRAN_LAPACK_DEF( tptrs, c, void, CALL_DEF_TPTRS( ComplexFloat ) )
FORTRAN_LAPACK_DEF( tptrs, z, void, CALL_DEF_TPTRS( ComplexDouble ) )

#undef CALL_DEF_TPTRS
#define CALL_DEF_LASWP( ValueType ) const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* n, ValueType* a, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* lda, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* k1, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* k2, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* ipiv, const scai::blaskernel::LAPACKDefinitions::LAPACKIndexType* incx

FORTRAN_LAPACK_DEF( laswp, s, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType, CALL_DEF_LASWP( float ))
FORTRAN_LAPACK_DEF( laswp, d, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType, CALL_DEF_LASWP( double ))
FORTRAN_LAPACK_DEF( laswp, c, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType, CALL_DEF_LASWP( ComplexFloat ))
FORTRAN_LAPACK_DEF( laswp, z, scai::blaskernel::LAPACKDefinitions::LAPACKIndexType, CALL_DEF_LASWP( ComplexDouble ))

} /* end extern "C" */
