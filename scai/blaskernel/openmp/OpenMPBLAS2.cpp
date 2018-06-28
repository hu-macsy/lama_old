/**
 * @file OpenMPBLAS2.cpp
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
 * @brief Own implementation of BLAS2 kernels for Host using OpenMP parallelization.
 * @author Eric Schricker
 * @date 09.10.2013
 */

// hpp
#include <scai/blaskernel/openmp/OpenMPBLAS2.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

// internal scai libraries
#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/tracing.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

#include <functional>

namespace scai
{

using tasking::TaskSyncToken;
using common::TypeTraits;

namespace blaskernel
{

SCAI_LOG_DEF_LOGGER( OpenMPBLAS2::logger, "OpenMP.BLAS2" )

/** gemv */

template<typename ValueType>
void OpenMPBLAS2::gemv(
    const common::MatrixOp opA, 
    const IndexType M,
    const IndexType N,
    const ValueType alpha,
    const ValueType* A,
    const IndexType lda,
    const ValueType* X,
    const IndexType incX,
    const ValueType beta,
    ValueType* Y,
    const IndexType incY )
{
    SCAI_REGION( "OpenMP.BLAS2.gemv" )

    SCAI_LOG_INFO( logger,
                   "gemv<" << common::TypeTraits<ValueType>::id() << ">: M = " << M << ", N = " << N
                   << ", LDA = " << lda << ", incX = " << incX << ", incY = " << incY
                   << ", alpha = " << alpha << ", beta = " << beta )

    if ( M == 0 )
    {
        return; // empty X, Y, A
    }

    // N == 0: empty A, but deal with X, Y, we can handle this here
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    if ( opA == common::MatrixOp::NORMAL )
    {
        //'T'
        ValueType Z;

        if ( incX == 1 && incY == 1 )
        {
            #pragma omp parallel for private(Z) 

            for ( IndexType i = 0; i < M; i++ )
            {
                Z = static_cast<ValueType>( 0.0 );

                for ( IndexType j = 0; j < N; j++ )
                {
                    Z += A[lda * i + j] * X[j];
                }

                Y[i] = Z * alpha + Y[i] * beta;
            }
        }
        else
        {
            //incX != 1 || incY != 1
            #pragma omp parallel for private(Z) 
            for ( IndexType i = 0; i < M; i++ )
            {
                Z = static_cast<ValueType>( 0.0 );

                for ( IndexType j = 0; j < N; j++ )
                {
                    Z += A[lda * i + j] * X[j * incX];
                }

                Y[i * incY] = Z * alpha + Y[i * incY] * beta;
            }
        }
    }
    else if ( opA == common::MatrixOp::TRANSPOSE  )
    {
        //'N'
        ValueType Z;

        if ( incX == 1 && incY == 1 )
        {
            #pragma omp parallel for private(Z) 

            for ( IndexType i = 0; i < N; i++ )
            {
                Z = static_cast<ValueType>( 0.0 );

                for ( IndexType j = 0; j < M; j++ )
                {
                    Z += A[lda * j + i] * X[j];
                }

                Y[i] = Z * alpha + Y[i] * beta;
            }
        }
        else
        {
            //incX != 1 || incY != 1
            #pragma omp parallel for private(Z) 
            for ( IndexType i = 0; i < N; i++ )
            {
                Z = static_cast<ValueType>( 0.0 );

                for ( IndexType j = 0; j < M; j++ )
                {
                    Z += A[lda * j + i] * X[j * incX];
                }

                Y[i * incY] = Z * alpha + Y[i * incY] * beta;
            }
        }
    }
    else if ( opA == common::MatrixOp::CONJ_TRANSPOSE )
    {
        #pragma omp parallel for 

        for ( IndexType i = 0; i < N; i++ )
        {
            ValueType Z = 0;

            for ( IndexType j = 0; j < M; j++ )
            {
                Z += common::Math::conj( A[lda * j + i] ) * X[j * incX];
            }

            Y[i * incY] = Z * alpha + Y[i * incY] * beta;
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "illegal opA setting " << opA )
    }
}

/* ---- geam :  general add (dense) matrix  ---------------------------------- */

template<typename ValueType>
static inline void applyInPlace( ValueType& v, const ValueType alpha, bool conj )
{
    v *= alpha;
    if ( conj )
    {
        v = common::Math::conj( v );
    }
}

template<typename ValueType>
static void transposeInPlace( ValueType* C, const IndexType m, const IndexType n, const ValueType alpha, bool conj )
{
    if ( m == n )
    {
        // transpose square matrix

        for ( IndexType i = 0; i < m; ++i )
        {
            applyInPlace( C[ i * n + i ], alpha, conj );


            for ( IndexType j = i + 1; j < n; ++j )
            {
                applyInPlace( C[ i * n + j ], alpha, conj );
                applyInPlace( C[ j * n + i ], alpha, conj );

                std::swap( C[i + n * j], C[j + n * i] );
            }
        }
    }
    else
    {
        for ( IndexType start = 0; start < n * m; ++start )
        {
            IndexType next = start;
            IndexType i = 0;           // counts number of entries in permutation chain

            do
            {
                ++i;
                // get new position of element
                next = ( next % m ) * n + next / m;
            }
            while ( next > start );

            if ( next < start )
            {
                // start was already in a previous permutation chain
                continue;
            }

            if ( i == 1 )
            {
                // permutation chain of one single element
                applyInPlace( C[start], alpha, conj );
            }
            else
            {
                const ValueType tmp = C[start];   // save first element of chain

                next = start;

                do
                {
                    i = ( next % m ) * n + next / m;
                    C[next] = ( i == start ) ? tmp :  C[i];
                    applyInPlace( C[next], alpha, conj );
                    next = i;
                }
                while ( next > start );
            }
        }
    }
}

/* ---- scale :  scale (dense) matrix  ---------------------------------- */

template<typename ValueType>
void OpenMPBLAS2::scale(
    ValueType* C,
    const IndexType ldc,
    const IndexType m,
    const IndexType n,
    const ValueType alpha,
    const ValueType* A,
    const IndexType lda,
    const common::MatrixOp opA )
{
    if ( alpha == 0 )
    {
        #pragma omp parallel for
        for ( IndexType i = 0; i < m; ++i )
        {
            for ( IndexType j = 0; j < n; ++j )
            {
                C[ i * ldc + j ] = alpha;
            }
        }
    }
    else if ( opA == common::MatrixOp::NORMAL )
    {
        #pragma omp parallel for 
        for ( IndexType i = 0; i < m; ++i )
        {   
            for ( IndexType j = 0; j < n; ++j )
            {   
                C[ i * ldc + j ] = alpha * A[ i * lda + j ];
            }
        }
    }
    else if ( opA == common::MatrixOp::TRANSPOSE )
    {
        if ( C == A )
        {
            SCAI_ASSERT_EQ_ERROR( ldc, n, "illegal transpose" )
            SCAI_ASSERT_EQ_ERROR( lda, m, "illegal transpose" )
            // in-place transpose 
            transposeInPlace( C, m, n, alpha, false );
        }
        else 
        {
            #pragma omp parallel for 
            for ( IndexType i = 0; i < m; ++i )
            {   
                for ( IndexType j = 0; j < n; ++j )
                {   
                    C[ i * ldc + j ] = alpha * A[ j * lda + i ];
                }
            }
        }
    }
    else if ( opA == common::MatrixOp::CONJ )
    {
        #pragma omp parallel for 
        for ( IndexType i = 0; i < m; ++i )
        {
            for ( IndexType j = 0; j < n; ++j )
            {   
                C[ i * ldc + j ] = alpha * common::Math::conj( A[ i * lda + j ] );
            }
        }
    }
    else if ( opA == common::MatrixOp::CONJ_TRANSPOSE )
    {
        if ( C == A )
        {
            SCAI_ASSERT_EQ_ERROR( ldc, n, "illegal transpose for aliased A and C" )
            SCAI_ASSERT_EQ_ERROR( lda, m, "illegal transpose for aliased A and C" )
            transposeInPlace( C, m, n, alpha, true );
        }
        else
        {
            #pragma omp parallel for 
            for ( IndexType i = 0; i < m; ++i )
            {
                for ( IndexType j = 0; j < n; ++j )
                {
                    C[ i * ldc + j ] = alpha * common::Math::conj( A[ j * lda + i ] );
                }
            }
        }
    }
}

/* ---- geam :  general add (dense) matrix  ---------------------------------- */

template<typename ValueType>
void OpenMPBLAS2::geam(
    ValueType* C,
    const IndexType ldc,
    const IndexType m,
    const IndexType n,
    const ValueType alpha,
    const ValueType* A,
    const IndexType lda,
    const common::MatrixOp opA,
    const ValueType beta,
    const ValueType* B,
    const IndexType ldb,
    const common::MatrixOp opB )
{
    SCAI_LOG_INFO( logger, "geam: C ( " << m << " x " << n << " ) = " 
                             << alpha << " * A ( op = " << opA << " ) + " 
                             << beta << " * B ( op = " << opB << " )" )

    if ( beta == 0 )
    {
        // C = alpha * opA( A )

        scale( C, ldc, m, n, alpha, A, lda, opA );
        return;
    }

    if ( alpha == 0 )
    {
        // C = beta * opB( A )

        scale( C, ldc, m, n, beta, B, ldb, opB );
        return;
    }

    if ( opA == common::MatrixOp::NORMAL )
    {
        if ( opB == common::MatrixOp::NORMAL )
        {
            // C = alpha * A + beta * B

            for ( IndexType i = 0; i < m; ++i )
            {
                for ( IndexType j = 0; j < n; ++j )
                {
                    C[ i * ldc + j ] = alpha * A[ i * lda + j ] + beta * B[ i * ldb + j ];
                }
            }
        }
        else
        {
            // C = alpha * A + beta * transpose( B ) -> C = beta * transpose( B ), C = alpha * A + 1 * C

            scale( C, ldc, m, n, beta, B, ldb, opB );
            SCAI_ASSERT_NE_ERROR( C, A, "alias problem" )
            geam( C, ldc, m, n, alpha, A, lda, opA, ValueType( 1 ), C, ldc, common::MatrixOp::NORMAL );
        }
    }
    else if ( opB == common::MatrixOp::NORMAL )
    {
        SCAI_ASSERT_NE_ERROR( C, B, "alias problem" )
        scale( C, ldc, m, n, alpha, A, lda, opB );
        geam( C, ldc, m, n, ValueType( 1 ), C, ldc, common::MatrixOp::NORMAL, beta, B, ldb, opB );
    }
    else 
    {
        COMMON_THROWEXCEPTION( "general add matrix unsupported if not at least one op is NORMAL" )
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPBLAS2::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register BLAS2 routines for OpenMP in Kernel Registry" )
    KernelRegistry::set<BLASKernelTrait::gemv<ValueType> >( OpenMPBLAS2::gemv, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::geam<ValueType> >( OpenMPBLAS2::geam, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPBLAS2::OpenMPBLAS2()
{
    SCAI_LOG_INFO( logger, "register BLAS2 routines for OpenMP in Kernel Registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ADD );
}

OpenMPBLAS2::~OpenMPBLAS2()
{
    SCAI_LOG_INFO( logger, "unregister BLAS2 routines for OpenMP in Kernel Registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPBLAS2 OpenMPBLAS2::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
