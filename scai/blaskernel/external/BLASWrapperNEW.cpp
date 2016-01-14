#include <scai/blaskernel/external/BLASWrapperNEW.hpp>

typedef scai::blaskernel::_BLASWrapper::BLASIndexType BLASIndexType;
typedef scai::blaskernel::_BLASWrapper::BLASTrans BLASTrans;

extern "C" {

#define FORTRAN_BLAS_NAME(prefix, name ) prefix##name##_
#define FORTRAN_BLAS_DEF( name, prefix, definition ) 			\
        void FORTRAN_BLAS_NAME( prefix, name )( definition );

#define CALL_DEF_SWAP( type ) const BLASIndexType *, type*, const BLASIndexType *, type*, const BLASIndexType *

FORTRAN_BLAS_DEF( swap, s, CALL_DEF_SWAP( float ) )
FORTRAN_BLAS_DEF( swap, d, CALL_DEF_SWAP( double ) )
FORTRAN_BLAS_DEF( swap, c, CALL_DEF_SWAP( ComplexFloat ) )
FORTRAN_BLAS_DEF( swap, z, CALL_DEF_SWAP( ComplexDouble ) )

#undef CALL_DEF_SWAP
#define CALL_DEF_COPY( type ) const BLASIndexType *, const type*, const BLASIndexType *, type*, const BLASIndexType *

FORTRAN_BLAS_DEF( copy, s, CALL_DEF_COPY( float ) );
FORTRAN_BLAS_DEF( copy, d, CALL_DEF_COPY( double ) );
FORTRAN_BLAS_DEF( copy, c, CALL_DEF_COPY( ComplexFloat ) );
FORTRAN_BLAS_DEF( copy, z, CALL_DEF_COPY( ComplexDouble ) );

#undef CALL_DEF_COPY
#define CALL_DEF_AXPY( type ) const BLASIndexType *, const type*, const type*, const BLASIndexType *, type*, const BLASIndexType *

FORTRAN_BLAS_DEF( axpy, s, CALL_DEF_AXPY( float ) );
FORTRAN_BLAS_DEF( axpy, d, CALL_DEF_AXPY( double ) );
FORTRAN_BLAS_DEF( axpy, c, CALL_DEF_AXPY( ComplexFloat ) );
FORTRAN_BLAS_DEF( axpy, z, CALL_DEF_AXPY( ComplexDouble ) );

#undef CALL_DEF_AXPY
#define CALL_DEF_DOT( type ) const BLASIndexType *, const type*, const BLASIndexType *, const type*, const BLASIndexType *, type*

//FORTRAN_BLAS_DEF( dot_sub, s, CALL_DEF( float ) );
//FORTRAN_BLAS_DEF( dot_sub, d, CALL_DEF( double ) );
//FORTRAN_BLAS_DEF( dotu_sub, c, CALL_DEF( ComplexFloat ) );
//FORTRAN_BLAS_DEF( dotu_sub, z, CALL_DEF( ComplexDouble ) );
FORTRAN_BLAS_DEF( dot, s, CALL_DEF_DOT( float ) );
FORTRAN_BLAS_DEF( dot, d, CALL_DEF_DOT( double ) );
FORTRAN_BLAS_DEF( dotu, c, CALL_DEF_DOT( ComplexFloat ) );
FORTRAN_BLAS_DEF( dotu, z, CALL_DEF_DOT( ComplexDouble ) );

#undef CALL_DEF_DOT
#define CALL_DEF_SCAL( type ) const BLASIndexType *, const type*, type*, const BLASIndexType *

FORTRAN_BLAS_DEF( scal, s, CALL_DEF_SCAL( float ) );
FORTRAN_BLAS_DEF( scal, d, CALL_DEF_SCAL( double ) );
FORTRAN_BLAS_DEF( scal, c, CALL_DEF_SCAL( ComplexFloat ) );
FORTRAN_BLAS_DEF( scal, z, CALL_DEF_SCAL( ComplexDouble ) );

#undef CALL_DEF_SCAL
#define CALL_DEF_NRM2( type ) const BLASIndexType *, const type*, const BLASIndexType *, type*

//FORTRAN_BLAS_DEF( nrm2_sub, s, CALL_DEF( float ) );
//FORTRAN_BLAS_DEF( nrm2_sub, d, CALL_DEF( double ) );
//FORTRAN_BLAS_DEF( nrm2_sub, sc, CALL_DEF( ComplexFloat ) );
//FORTRAN_BLAS_DEF( nrm2_sub, dz, CALL_DEF( ComplexDouble ) );
FORTRAN_BLAS_DEF( nrm2, s, CALL_DEF_NRM2( float ) );
FORTRAN_BLAS_DEF( nrm2, d, CALL_DEF_NRM2( double ) );
FORTRAN_BLAS_DEF( nrm2, sc, CALL_DEF_NRM2( ComplexFloat ) );
FORTRAN_BLAS_DEF( nrm2, dz, CALL_DEF_NRM2( ComplexDouble ) );

#undef CALL_DEF_NRM2
#define CALL_DEF_ASUM( type ) const BLASIndexType *, const type*, const BLASIndexType *, type*

//FORTRAN_BLAS_DEF( asum_sub, s, CALL_DEF( float ) );
//FORTRAN_BLAS_DEF( asum_sub, d, CALL_DEF( double ) );
//FORTRAN_BLAS_DEF( asum_sub, sc, CALL_DEF( ComplexFloat ) );
//FORTRAN_BLAS_DEF( asum_sub, dz, CALL_DEF( ComplexDouble ) );
FORTRAN_BLAS_DEF( asum, s, CALL_DEF_ASUM( float ) );
FORTRAN_BLAS_DEF( asum, d, CALL_DEF_ASUM( double ) );
FORTRAN_BLAS_DEF( asum, sc, CALL_DEF_ASUM( ComplexFloat ) );
FORTRAN_BLAS_DEF( asum, dz, CALL_DEF_ASUM( ComplexDouble ) );

#undef CALL_DEF_ASUM
#define CALL_DEF_AMAX( type ) const BLASIndexType *, const type*, const BLASIndexType *, BLASIndexType *

//FORTRAN_BLAS_DEF( amax_sub, is, CALL_DEF( float ) );
//FORTRAN_BLAS_DEF( amax_sub, id, CALL_DEF( double ) );
//FORTRAN_BLAS_DEF( amax_sub, ic, CALL_DEF( ComplexFloat ) );
//FORTRAN_BLAS_DEF( amax_sub, iz, CALL_DEF( ComplexDouble ) );
FORTRAN_BLAS_DEF( amax, is, CALL_DEF_AMAX( float ) );
FORTRAN_BLAS_DEF( amax, id, CALL_DEF_AMAX( double ) );
FORTRAN_BLAS_DEF( amax, ic, CALL_DEF_AMAX( ComplexFloat ) );
FORTRAN_BLAS_DEF( amax, iz, CALL_DEF_AMAX( ComplexDouble ) );

#undef CALL_DEF_AMAX
#define CALL_DEF_GEMV( type ) const BLASTrans*, const BLASIndexType *, const BLASIndexType *, const type*, const type*, const BLASIndexType *, const type*, const BLASIndexType *, const type*, type*, const BLASIndexType *

FORTRAN_BLAS_DEF( gemv, s, CALL_DEF_GEMV( float ) );
FORTRAN_BLAS_DEF( gemv, d, CALL_DEF_GEMV( double ) );
FORTRAN_BLAS_DEF( gemv, c, CALL_DEF_GEMV( ComplexFloat ) );
FORTRAN_BLAS_DEF( gemv, z, CALL_DEF_GEMV( ComplexDouble ) );

#undef CALL_DEF_GEMM
#define CALL_DEF_GEMM( type ) const BLASTrans*, const BLASTrans*, const BLASIndexType *, const BLASIndexType *, const BLASIndexType *, const type*, const type*, const BLASIndexType *, const type*, const BLASIndexType *, const type*, type*, const BLASIndexType *

FORTRAN_BLAS_DEF( gemm, s, CALL_DEF_GEMM( float ) );
FORTRAN_BLAS_DEF( gemm, d, CALL_DEF_GEMM( double ) );
FORTRAN_BLAS_DEF( gemm, c, CALL_DEF_GEMM( ComplexFloat ) );
FORTRAN_BLAS_DEF( gemm, z, CALL_DEF_GEMM( ComplexDouble ) );

#undef CALL_DEF_GEMM

} /* end extern "C" */

namespace scai
{

namespace blaskernel
{

template<>
bool BLASWrapper<float>::initialize()
{
	std::cout << "BLASWrapper<float> is the universe" << std::endl;

	BLAS_PTR_SET( scal,		float, FORTRAN_BLAS_NAME( s, scal ) )
	BLAS_PTR_SET( nrm2,		float, FORTRAN_BLAS_NAME( s, nrm2 ) )
	BLAS_PTR_SET( asum,		float, FORTRAN_BLAS_NAME( s, asum ) )
	BLAS_PTR_SET( iamax,	float, FORTRAN_BLAS_NAME( is, amax ) )
	BLAS_PTR_SET( swap,		float, FORTRAN_BLAS_NAME( s, swap ) )
	BLAS_PTR_SET( copy,		float, FORTRAN_BLAS_NAME( s, copy ) )
	BLAS_PTR_SET( axpy,		float, FORTRAN_BLAS_NAME( s, axpy ) )
	BLAS_PTR_SET( dot,		float, FORTRAN_BLAS_NAME( s, dot ) )
	BLAS_PTR_SET( gemv,		float, FORTRAN_BLAS_NAME( s, gemv ) )
	BLAS_PTR_SET( gemm,		float, FORTRAN_BLAS_NAME( s, gemm ) )

	return true;
}

template<>
bool BLASWrapper<double>::initialize()
{
	std::cout << "BLASWrapper<double> is the universe" << std::endl;

	BLAS_PTR_SET( scal,		double,	FORTRAN_BLAS_NAME( d, scal ) )
	BLAS_PTR_SET( nrm2,		double,	FORTRAN_BLAS_NAME( d, nrm2 ) )
	BLAS_PTR_SET( asum,		double,	FORTRAN_BLAS_NAME( d, asum ) )
	BLAS_PTR_SET( iamax,	double,	FORTRAN_BLAS_NAME( id, amax ) )
	BLAS_PTR_SET( swap,		double,	FORTRAN_BLAS_NAME( d, swap ) )
	BLAS_PTR_SET( copy,		double,	FORTRAN_BLAS_NAME( d, copy ) )
	BLAS_PTR_SET( axpy,		double,	FORTRAN_BLAS_NAME( d, axpy ) )
	BLAS_PTR_SET( dot,		double,	FORTRAN_BLAS_NAME( d, dot ) )
	BLAS_PTR_SET( gemv,		double,	FORTRAN_BLAS_NAME( d, gemv ) )
	BLAS_PTR_SET( gemm,		double,	FORTRAN_BLAS_NAME( d, gemm ) )

	return true;
}

template<>
bool BLASWrapper<ComplexFloat>::initialize()
{
	BLAS_PTR_SET( scal,		ComplexFloat, FORTRAN_BLAS_NAME( c, scal ) )
	BLAS_PTR_SET( nrm2,		ComplexFloat, FORTRAN_BLAS_NAME( sc, nrm2 ) )
	BLAS_PTR_SET( asum,		ComplexFloat, FORTRAN_BLAS_NAME( sc, asum ) )
	BLAS_PTR_SET( iamax,	ComplexFloat, FORTRAN_BLAS_NAME( ic, amax ) )
	BLAS_PTR_SET( swap,		ComplexFloat, FORTRAN_BLAS_NAME( c, swap ) )
	BLAS_PTR_SET( copy,		ComplexFloat, FORTRAN_BLAS_NAME( c, copy ) )
	BLAS_PTR_SET( axpy,		ComplexFloat, FORTRAN_BLAS_NAME( c, axpy ) )
	BLAS_PTR_SET( dot,		ComplexFloat, FORTRAN_BLAS_NAME( c, dotu ) )
	BLAS_PTR_SET( gemv,		ComplexFloat, FORTRAN_BLAS_NAME( c, gemv ) )
	BLAS_PTR_SET( gemm,		ComplexFloat, FORTRAN_BLAS_NAME( c, gemm ) )

	return true;
}

template<>
bool BLASWrapper<ComplexDouble>::initialize()
{
	BLAS_PTR_SET( scal,		ComplexDouble,	FORTRAN_BLAS_NAME( z, scal ) )
	BLAS_PTR_SET( nrm2,		ComplexDouble,	FORTRAN_BLAS_NAME( dz, nrm2 ) )
	BLAS_PTR_SET( asum,		ComplexDouble,	FORTRAN_BLAS_NAME( dz, asum ) )
	BLAS_PTR_SET( iamax,	ComplexDouble,	FORTRAN_BLAS_NAME( iz, amax ) )
	BLAS_PTR_SET( swap,		ComplexDouble,	FORTRAN_BLAS_NAME( z, swap ) )
	BLAS_PTR_SET( copy,		ComplexDouble,	FORTRAN_BLAS_NAME( z, copy ) )
	BLAS_PTR_SET( axpy,		ComplexDouble,	FORTRAN_BLAS_NAME( z, axpy ) )
	BLAS_PTR_SET( dot,		ComplexDouble,	FORTRAN_BLAS_NAME( z, dotu ) )
	BLAS_PTR_SET( gemv,		ComplexDouble,	FORTRAN_BLAS_NAME( z, gemv ) )
	BLAS_PTR_SET( gemm,		ComplexDouble,	FORTRAN_BLAS_NAME( z, gemm ) )

	return true;
}

template<typename ValueType>
bool BLASWrapper<ValueType>::init = initialize();

BLAS_PTR_CRE( scal )
BLAS_PTR_CRE( nrm2 )
BLAS_PTR_CRE( asum )
BLAS_PTR_CRE( iamax )
BLAS_PTR_CRE( swap )
BLAS_PTR_CRE( copy )
BLAS_PTR_CRE( axpy )
BLAS_PTR_CRE( dot )
BLAS_PTR_CRE( gemv )
BLAS_PTR_CRE( gemm )

} /* end namespace blaskernel */

} /* end namespace scai */
