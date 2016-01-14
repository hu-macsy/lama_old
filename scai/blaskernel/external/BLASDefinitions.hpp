/*
 * scai::blaskernel::BLASDefinitions.hpp
 *
 *  Created on: Jan 14, 2016
 *      Author: eschricker
 */

#pragma once

#define FORTRAN_BLAS_NAME( name, prefix ) prefix##name##_

#define FORTRAN_BLAS_DEF( name, prefix, retType, definition ) 			\
        retType FORTRAN_BLAS_NAME( name, prefix )( definition );

namespace scai {

namespace blaskernel {

class COMMON_DLL_IMPORTEXPORT BLASDefinitions
{
public:
#ifdef F77_INT
	typedef F77_INT BLASIndexType;
#else
	typedef int BLASIndexType;
#endif

#ifdef  F77_CHAR
    typedef F77_CHAR BLASTrans
#else
	typedef char BLASTrans;
#endif
};

} /* end namespace blaskernel */

} /* end namespace scai */

extern "C"
{
#define CALL_DEF_SWAP( type ) const scai::blaskernel::BLASDefinitions::BLASIndexType *, type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *, type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *

FORTRAN_BLAS_DEF( swap, s, void, CALL_DEF_SWAP( float ) )
FORTRAN_BLAS_DEF( swap, d, void, CALL_DEF_SWAP( double ) )
FORTRAN_BLAS_DEF( swap, c, void, CALL_DEF_SWAP( ComplexFloat ) )
FORTRAN_BLAS_DEF( swap, z, void, CALL_DEF_SWAP( ComplexDouble ) )

#undef CALL_DEF_SWAP
#define CALL_DEF_COPY( type ) const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *, type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *

FORTRAN_BLAS_DEF( copy, s, void, CALL_DEF_COPY( float ) );
FORTRAN_BLAS_DEF( copy, d, void, CALL_DEF_COPY( double ) );
FORTRAN_BLAS_DEF( copy, c, void, CALL_DEF_COPY( ComplexFloat ) );
FORTRAN_BLAS_DEF( copy, z, void, CALL_DEF_COPY( ComplexDouble ) );

#undef CALL_DEF_COPY
#define CALL_DEF_AXPY( type ) const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, const type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *, type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *

FORTRAN_BLAS_DEF( axpy, s, void, CALL_DEF_AXPY( float ) );
FORTRAN_BLAS_DEF( axpy, d, void, CALL_DEF_AXPY( double ) );
FORTRAN_BLAS_DEF( axpy, c, void, CALL_DEF_AXPY( ComplexFloat ) );
FORTRAN_BLAS_DEF( axpy, z, void, CALL_DEF_AXPY( ComplexDouble ) );

#undef CALL_DEF_AXPY
#define CALL_DEF_DOT( type ) const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *

FORTRAN_BLAS_DEF( dot, s, float, CALL_DEF_DOT( float ) );
FORTRAN_BLAS_DEF( dot, d, double, CALL_DEF_DOT( double ) );
FORTRAN_BLAS_DEF( dotu, c, ComplexFloat, CALL_DEF_DOT( ComplexFloat ) );
FORTRAN_BLAS_DEF( dotu, z, ComplexDouble, CALL_DEF_DOT( ComplexDouble ) );

#undef CALL_DEF_DOT
#define CALL_DEF_SCAL( type ) const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *

FORTRAN_BLAS_DEF( scal, s, void, CALL_DEF_SCAL( float ) );
FORTRAN_BLAS_DEF( scal, d, void, CALL_DEF_SCAL( double ) );
FORTRAN_BLAS_DEF( scal, c, void, CALL_DEF_SCAL( ComplexFloat ) );
FORTRAN_BLAS_DEF( scal, z, void, CALL_DEF_SCAL( ComplexDouble ) );

#undef CALL_DEF_SCAL
#define CALL_DEF_NRM2( type ) const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *

FORTRAN_BLAS_DEF( nrm2, s, float, CALL_DEF_NRM2( float ) );
FORTRAN_BLAS_DEF( nrm2, d, double, CALL_DEF_NRM2( double ) );
FORTRAN_BLAS_DEF( nrm2, sc, float, CALL_DEF_NRM2( ComplexFloat ) );
FORTRAN_BLAS_DEF( nrm2, dz, double, CALL_DEF_NRM2( ComplexDouble ) );

#undef CALL_DEF_NRM2
#define CALL_DEF_ASUM( type ) const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *

FORTRAN_BLAS_DEF( asum, s, float, CALL_DEF_ASUM( float ) );
FORTRAN_BLAS_DEF( asum, d, double, CALL_DEF_ASUM( double ) );
FORTRAN_BLAS_DEF( asum, sc, float, CALL_DEF_ASUM( ComplexFloat ) );
FORTRAN_BLAS_DEF( asum, dz, double, CALL_DEF_ASUM( ComplexDouble ) );

#undef CALL_DEF_ASUM
#define CALL_DEF_AMAX( type ) const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *

FORTRAN_BLAS_DEF( amax, is, scai::blaskernel::BLASDefinitions::BLASIndexType, CALL_DEF_AMAX( float ) );
FORTRAN_BLAS_DEF( amax, id, scai::blaskernel::BLASDefinitions::BLASIndexType, CALL_DEF_AMAX( double ) );
FORTRAN_BLAS_DEF( amax, ic, scai::blaskernel::BLASDefinitions::BLASIndexType, CALL_DEF_AMAX( ComplexFloat ) );
FORTRAN_BLAS_DEF( amax, iz, scai::blaskernel::BLASDefinitions::BLASIndexType, CALL_DEF_AMAX( ComplexDouble ) );

#undef CALL_DEF_AMAX
#define CALL_DEF_GEMV( type ) const scai::blaskernel::BLASDefinitions::BLASTrans*, const scai::blaskernel::BLASDefinitions::BLASIndexType *, const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, const type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *

FORTRAN_BLAS_DEF( gemv, s, void, CALL_DEF_GEMV( float ) );
FORTRAN_BLAS_DEF( gemv, d, void, CALL_DEF_GEMV( double ) );
FORTRAN_BLAS_DEF( gemv, c, void, CALL_DEF_GEMV( ComplexFloat ) );
FORTRAN_BLAS_DEF( gemv, z, void, CALL_DEF_GEMV( ComplexDouble ) );

#undef CALL_DEF_GEMM
#define CALL_DEF_GEMM( type ) const scai::blaskernel::BLASDefinitions::BLASTrans*, const scai::blaskernel::BLASDefinitions::BLASTrans*, const scai::blaskernel::BLASDefinitions::BLASIndexType *, const scai::blaskernel::BLASDefinitions::BLASIndexType *, const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, const type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *, const type*, type*, const scai::blaskernel::BLASDefinitions::BLASIndexType *

FORTRAN_BLAS_DEF( gemm, s, void, CALL_DEF_GEMM( float ) );
FORTRAN_BLAS_DEF( gemm, d, void, CALL_DEF_GEMM( double ) );
FORTRAN_BLAS_DEF( gemm, c, void, CALL_DEF_GEMM( ComplexFloat ) );
FORTRAN_BLAS_DEF( gemm, z, void, CALL_DEF_GEMM( ComplexDouble ) );

#undef CALL_DEF_GEMM

} /* end extern "C" */
