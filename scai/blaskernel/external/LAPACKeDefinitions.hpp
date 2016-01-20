/*
 * LAPACKDefinitions.hpp
 *
 *  Created on: Jan 14, 2016
 *      Author: eschricker
 */

#pragma once


#include <scai/common/config.hpp>

#include <mkl_lapacke.h>

#define FORTRAN_LAPACK_NAME( name, prefix ) prefix##name##_

#define FORTRAN_LAPACK_DEF( name, prefix, retType, definition ) 			\
        retType FORTRAN_LAPACK_NAME( name, prefix )( definition );

namespace scai {

namespace blaskernel {

class COMMON_DLL_IMPORTEXPORT LAPACKeDefinitions
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

