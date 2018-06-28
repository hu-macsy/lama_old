/**
 * @file BLASTrait.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Trait for BLAS interface
 * @author Eric Schricker
 * @date 10.01.2016
 */

#pragma once

#include <cstdio>
#include <scai/common/MatrixOp.hpp>

// macros
#define FORTRAN_BLAS_NAME( name, prefix ) prefix##name##_

#define FORTRAN_BLAS_DEF( name, prefix, retType, definition )           \
    retType FORTRAN_BLAS_NAME( name, prefix )( definition )

namespace scai
{

namespace blaskernel
{

class COMMON_DLL_IMPORTEXPORT BLASTrait
{
public:
#ifdef F77_INT
    typedef F77_INT BLASIndexType;
#else
    typedef int BLASIndexType;
#endif

#ifdef F77_CHAR
    typedef F77_CHAR BLASTrans;
#else
    typedef char BLASTrans;
#endif

    static inline BLASTrans castTrans( const common::MatrixOp op );
};

BLASTrait::BLASTrans BLASTrait::castTrans( const common::MatrixOp op )
{
    BLASTrait::BLASTrans castOp = ' ';

    if ( op == common::MatrixOp::NORMAL )
    {
        castOp = 'N';
    }
    else if ( op == common::MatrixOp::TRANSPOSE )
    {
        castOp = 'T';
    }
    else if ( op == common::MatrixOp::CONJ_TRANSPOSE )
    {
        castOp = 'C';
    }

    return castOp;
}

} /* end namespace blaskernel */

} /* end namespace scai */

extern "C"
{
#define CALL_DEF_SWAP( type ) const scai::blaskernel::BLASTrait::BLASIndexType *, type*, const scai::blaskernel::BLASTrait::BLASIndexType *, type*, const scai::blaskernel::BLASTrait::BLASIndexType *
#define CALL_DEF_COPY( type ) const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, const scai::blaskernel::BLASTrait::BLASIndexType *, type*, const scai::blaskernel::BLASTrait::BLASIndexType *
#define CALL_DEF_AXPY( type ) const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, const type*, const scai::blaskernel::BLASTrait::BLASIndexType *, type*, const scai::blaskernel::BLASTrait::BLASIndexType *
#define CALL_DEF_DOT( type ) const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, const scai::blaskernel::BLASTrait::BLASIndexType *
#define CALL_DEF_SCAL( type ) const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, type*, const scai::blaskernel::BLASTrait::BLASIndexType *
#define CALL_DEF_NRM2( type ) const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, const scai::blaskernel::BLASTrait::BLASIndexType *
#define CALL_DEF_ASUM( type ) const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, const scai::blaskernel::BLASTrait::BLASIndexType *
#define CALL_DEF_AMAX( type ) const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, const scai::blaskernel::BLASTrait::BLASIndexType *
#define CALL_DEF_GEMV( type ) const scai::blaskernel::BLASTrait::BLASTrans*, const scai::blaskernel::BLASTrait::BLASIndexType *, const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, const type*, const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, type*, const scai::blaskernel::BLASTrait::BLASIndexType *
#define CALL_DEF_GEMM( type ) const scai::blaskernel::BLASTrait::BLASTrans*, const scai::blaskernel::BLASTrait::BLASTrans*, const scai::blaskernel::BLASTrait::BLASIndexType *, const scai::blaskernel::BLASTrait::BLASIndexType *, const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, const type*, const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, const scai::blaskernel::BLASTrait::BLASIndexType *, const type*, type*, const scai::blaskernel::BLASTrait::BLASIndexType *

//float
    FORTRAN_BLAS_DEF( swap, s, void, CALL_DEF_SWAP( float ) );
    FORTRAN_BLAS_DEF( copy, s, void, CALL_DEF_COPY( float ) );
    FORTRAN_BLAS_DEF( axpy, s, void, CALL_DEF_AXPY( float ) );
    FORTRAN_BLAS_DEF( dot, s, float, CALL_DEF_DOT( float ) );
    FORTRAN_BLAS_DEF( scal, s, void, CALL_DEF_SCAL( float ) );
    FORTRAN_BLAS_DEF( nrm2, s, float, CALL_DEF_NRM2( float ) );
    FORTRAN_BLAS_DEF( asum, s, float, CALL_DEF_ASUM( float ) );
    FORTRAN_BLAS_DEF( amax, is, scai::blaskernel::BLASTrait::BLASIndexType, CALL_DEF_AMAX( float ) );
    FORTRAN_BLAS_DEF( gemv, s, void, CALL_DEF_GEMV( float ) );
    FORTRAN_BLAS_DEF( gemm, s, void, CALL_DEF_GEMM( float ) );
// double
    FORTRAN_BLAS_DEF( swap, d, void, CALL_DEF_SWAP( double ) );
    FORTRAN_BLAS_DEF( copy, d, void, CALL_DEF_COPY( double ) );
    FORTRAN_BLAS_DEF( axpy, d, void, CALL_DEF_AXPY( double ) );
    FORTRAN_BLAS_DEF( dot, d, double, CALL_DEF_DOT( double ) );
    FORTRAN_BLAS_DEF( scal, d, void, CALL_DEF_SCAL( double ) );
    FORTRAN_BLAS_DEF( nrm2, d, double, CALL_DEF_NRM2( double ) );
    FORTRAN_BLAS_DEF( asum, d, double, CALL_DEF_ASUM( double ) );
    FORTRAN_BLAS_DEF( amax, id, scai::blaskernel::BLASTrait::BLASIndexType, CALL_DEF_AMAX( double ) );
    FORTRAN_BLAS_DEF( gemm, d, void, CALL_DEF_GEMM( double ) );
    FORTRAN_BLAS_DEF( gemv, d, void, CALL_DEF_GEMV( double ) );

#ifdef SCAI_COMPLEX_SUPPORTED
    FORTRAN_BLAS_DEF( swap, c, void, CALL_DEF_SWAP( scai::ComplexFloat ) );
    FORTRAN_BLAS_DEF( copy, c, void, CALL_DEF_COPY( scai::ComplexFloat ) );
    FORTRAN_BLAS_DEF( axpy, c, void, CALL_DEF_AXPY( scai::ComplexFloat ) );

// should not be used
    FORTRAN_BLAS_DEF( dotc, c, float, CALL_DEF_DOT( scai::ComplexFloat ) );

    FORTRAN_BLAS_DEF( scal, c, void, CALL_DEF_SCAL( scai::ComplexFloat ) );
    FORTRAN_BLAS_DEF( nrm2, sc, float, CALL_DEF_NRM2( scai::ComplexFloat ) );
    FORTRAN_BLAS_DEF( asum, sc, float, CALL_DEF_ASUM( scai::ComplexFloat ) );
    FORTRAN_BLAS_DEF( amax, ic, scai::blaskernel::BLASTrait::BLASIndexType, CALL_DEF_AMAX( scai::ComplexFloat ) );
    FORTRAN_BLAS_DEF( gemv, c, void, CALL_DEF_GEMV( scai::ComplexFloat ) );
    FORTRAN_BLAS_DEF( gemm, c, void, CALL_DEF_GEMM( scai::ComplexFloat ) );
// ComplexDouble
    FORTRAN_BLAS_DEF( swap, z, void, CALL_DEF_SWAP( scai::ComplexDouble ) );
    FORTRAN_BLAS_DEF( copy, z, void, CALL_DEF_COPY( scai::ComplexDouble ) );
    FORTRAN_BLAS_DEF( axpy, z, void, CALL_DEF_AXPY( scai::ComplexDouble ) );

// should not be used
    FORTRAN_BLAS_DEF( dotc, z, double, CALL_DEF_DOT( scai::ComplexDouble ) );

    FORTRAN_BLAS_DEF( scal, z, void, CALL_DEF_SCAL( scai::ComplexDouble ) );
    FORTRAN_BLAS_DEF( nrm2, dz, double, CALL_DEF_NRM2( scai::ComplexDouble ) );
    FORTRAN_BLAS_DEF( asum, dz, double, CALL_DEF_ASUM( scai::ComplexDouble ) );
    FORTRAN_BLAS_DEF( amax, iz, scai::blaskernel::BLASTrait::BLASIndexType, CALL_DEF_AMAX( scai::ComplexDouble ) );
    FORTRAN_BLAS_DEF( gemv, z, void, CALL_DEF_GEMV( scai::ComplexDouble ) );
    FORTRAN_BLAS_DEF( gemm, z, void, CALL_DEF_GEMM( scai::ComplexDouble ) );
#endif

#undef CALL_DEF_SWAP
#undef CALL_DEF_COPY
#undef CALL_DEF_AXPY
#undef CALL_DEF_DOT
#undef CALL_DEF_SCAL
#undef CALL_DEF_NRM2
#undef CALL_DEF_ASUM
#undef CALL_DEF_AMAX
#undef CALL_DEF_GEMM
#undef CALL_DEF_GEMM

#undef FORTRAN_BLAS_DEF
} /* end extern "C" */
