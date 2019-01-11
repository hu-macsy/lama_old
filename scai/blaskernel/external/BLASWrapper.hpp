/**
 * @file BLASWrapper.hpp
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
 * @brief Wrapper for BLAS functions
 * @author Eric Schricker
 * @date 10.01.2016
 */

#pragma once

// local library
#include <scai/blaskernel/external/BLASTrait.hpp>
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/common/macros/unused.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>

namespace scai
{

namespace blaskernel
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT BLASWrapper;

#define BLASWRAPPER_DEF( ValueType, prefix1, prefix2, prefix3, DOT )                                                    \
    template<>                                                                                                              \
    class COMMON_DLL_IMPORTEXPORT BLASWrapper<ValueType>                                                                    \
    {                                                                                                                       \
    public:                                                                                                                 \
        typedef BLASTrait::BLASIndexType BLASIndexType;                                                                     \
        typedef BLASTrait::BLASTrans BLASTrans;                                                                             \
        \
        static void scal( const BLASIndexType n, const ValueType alpha, ValueType* x, const BLASIndexType incX )            \
        {                                                                                                                   \
            FORTRAN_BLAS_NAME( scal, prefix1 )( &n, &alpha, x, &incX );                                                     \
        }                                                                                                                   \
        \
        static ValueType nrm2(const BLASIndexType n, const ValueType *x, const BLASIndexType incX)                          \
        {                                                                                                                   \
            return FORTRAN_BLAS_NAME( nrm2, prefix2 )( &n, x, &incX );                                                      \
        }                                                                                                                   \
        \
        static ValueType asum(const BLASIndexType n, const ValueType *x,BLASIndexType incX)                                 \
        {                                                                                                                   \
            return FORTRAN_BLAS_NAME( asum, prefix2 )(&n, x, &incX);                                                        \
        }                                                                                                                   \
        \
        static BLASIndexType iamax(const BLASIndexType n, const ValueType *x, const BLASIndexType incX)                     \
        {                                                                                                                   \
            return FORTRAN_BLAS_NAME( amax, prefix3 )(&n, x, &incX);                                                        \
        }                                                                                                                   \
        \
        static void swap(const BLASIndexType n, ValueType *x, const BLASIndexType incX, ValueType *y,                       \
                         const BLASIndexType incY)                                                                                   \
        {                                                                                                                   \
            FORTRAN_BLAS_NAME( swap, prefix1 )(&n, x, &incX, y, &incY);                                                     \
        }                                                                                                                   \
        \
        static void copy(const BLASIndexType n, const ValueType *x, const BLASIndexType incX, ValueType *y,                 \
                         const BLASIndexType incY)                                                                                   \
        {                                                                                                                   \
            FORTRAN_BLAS_NAME( copy, prefix1 )(&n, x,   &incX, y, &incY);                                                   \
        }                                                                                                                   \
        \
        static void axpy(const BLASIndexType n, const ValueType alpha,                                                      \
                         const ValueType *x, const BLASIndexType incX, ValueType *y, const BLASIndexType incY)                       \
        {                                                                                                                   \
            FORTRAN_BLAS_NAME( axpy, prefix1 )(&n, &alpha, x,   &incX, y, &incY);                                           \
        }                                                                                                                   \
        \
        static ValueType dot(const BLASIndexType n, const ValueType *x,                                                     \
                             const BLASIndexType incX, const ValueType *y, const BLASIndexType incY)                                     \
        {                                                                                                                   \
            return FORTRAN_BLAS_NAME( DOT, prefix1 )(&n, x, &incX, y, &incY);                                               \
        }                                                                                                                   \
        \
        static void gemv( const BLASTrans transA, const BLASIndexType m, const BLASIndexType n,                             \
                          const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* x,                     \
                          const BLASIndexType incX, const ValueType beta, ValueType* y, const BLASIndexType incY)                     \
        {                                                                                                                   \
            FORTRAN_BLAS_NAME( gemv, prefix1 )( &transA, &m, &n, &alpha, A, &lda, x, &incX, &beta, y, &incY );              \
        }                                                                                                                   \
        \
        static void gemm( const BLASTrans transA, const BLASTrans transB,                                                   \
                          const BLASIndexType m, const BLASIndexType n, const BLASIndexType k,                                        \
                          const ValueType alpha, const ValueType* A, const BLASIndexType lda, const ValueType* B,                     \
                          const BLASIndexType ldb, const ValueType beta, ValueType* C, const BLASIndexType ldc)                       \
        {                                                                                                                   \
            FORTRAN_BLAS_NAME( gemm, prefix1 )( &transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );   \
        }                                                                                                                   \
    };

BLASWRAPPER_DEF( float, s, s, is, dot );
BLASWRAPPER_DEF( double, d, d, id, dot );

#ifdef SCAI_COMPLEX_SUPPORTED
BLASWRAPPER_DEF( ComplexFloat, c, sc, ic, dotc );
BLASWRAPPER_DEF( ComplexDouble, z, dz, iz, dotc );
#endif

#undef BLASWRAPPER_DEF

} /* end namespace blaskernel */

} /* end namespace scai */

