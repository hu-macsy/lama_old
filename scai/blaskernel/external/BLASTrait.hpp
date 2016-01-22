/**
 * @file BLASTrait.cpp
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
 * @brief Trait for BLAS interface
 * @author Eric Schricker
 * @date 10.01.2016
 * @since 2.0.0
 */

#pragma once

// macros
#define FORTRAN_BLAS_NAME( name, prefix ) prefix##name##_

#define FORTRAN_BLAS_DEF( name, prefix, retType, definition ) 			\
        retType FORTRAN_BLAS_NAME( name, prefix )( definition );

namespace scai {

namespace blaskernel {

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
};

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
FORTRAN_BLAS_DEF( swap, s, void, CALL_DEF_SWAP( float ) )
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
FORTRAN_BLAS_DEF( swap, d, void, CALL_DEF_SWAP( double ) )
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
// ComplexFloat
FORTRAN_BLAS_DEF( swap, c, void, CALL_DEF_SWAP( ComplexFloat ) )
FORTRAN_BLAS_DEF( copy, c, void, CALL_DEF_COPY( ComplexFloat ) );
FORTRAN_BLAS_DEF( axpy, c, void, CALL_DEF_AXPY( ComplexFloat ) );
FORTRAN_BLAS_DEF( dotu, c, ComplexFloat, CALL_DEF_DOT( ComplexFloat ) );
FORTRAN_BLAS_DEF( scal, c, void, CALL_DEF_SCAL( ComplexFloat ) );
FORTRAN_BLAS_DEF( nrm2, sc, float, CALL_DEF_NRM2( ComplexFloat ) );
FORTRAN_BLAS_DEF( asum, sc, float, CALL_DEF_ASUM( ComplexFloat ) );
FORTRAN_BLAS_DEF( amax, ic, scai::blaskernel::BLASTrait::BLASIndexType, CALL_DEF_AMAX( ComplexFloat ) );
FORTRAN_BLAS_DEF( gemv, c, void, CALL_DEF_GEMV( ComplexFloat ) );
FORTRAN_BLAS_DEF( gemm, c, void, CALL_DEF_GEMM( ComplexFloat ) );
// ComplexDouble
FORTRAN_BLAS_DEF( swap, z, void, CALL_DEF_SWAP( ComplexDouble ) )
FORTRAN_BLAS_DEF( copy, z, void, CALL_DEF_COPY( ComplexDouble ) );
FORTRAN_BLAS_DEF( axpy, z, void, CALL_DEF_AXPY( ComplexDouble ) );
FORTRAN_BLAS_DEF( dotu, z, ComplexDouble, CALL_DEF_DOT( ComplexDouble ) );
FORTRAN_BLAS_DEF( scal, z, void, CALL_DEF_SCAL( ComplexDouble ) );
FORTRAN_BLAS_DEF( nrm2, dz, double, CALL_DEF_NRM2( ComplexDouble ) );
FORTRAN_BLAS_DEF( asum, dz, double, CALL_DEF_ASUM( ComplexDouble ) );
FORTRAN_BLAS_DEF( amax, iz, scai::blaskernel::BLASTrait::BLASIndexType, CALL_DEF_AMAX( ComplexDouble ) );
FORTRAN_BLAS_DEF( gemv, z, void, CALL_DEF_GEMV( ComplexDouble ) );
FORTRAN_BLAS_DEF( gemm, z, void, CALL_DEF_GEMM( ComplexDouble ) );
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
