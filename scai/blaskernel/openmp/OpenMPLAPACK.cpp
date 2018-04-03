/**
 * @file OpenMPLAPACK.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Own implementation of used LAPACK routines for Host using OpenMP parallelization.
 * @author Eric Schricker
 * @date 19.11.2013
 */

// hpp
#include <scai/blaskernel/openmp/OpenMPLAPACK.hpp>

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/openmp/OpenMPBLAS1.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

// std
#include <cmath>
#include <memory>

namespace scai
{

namespace blaskernel
{

using common::TypeTraits;

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( OpenMPLAPACK::logger, "OpenMP.LAPACK" )

/* ------------------------------------------------------------------------- */
/*      getrf<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPLAPACK::getrf(
    const IndexType m,
    const IndexType n,
    ValueType* const a,
    const IndexType lda,
    IndexType* const ipiv )
{
    SCAI_REGION( "OpenMP.LAPACK.getrf<ValueType>" )
    SCAI_LOG_INFO( logger, "getrf<" << TypeTraits<ValueType>::id() << "> for A of size " << m << " x " << n )
    IndexType index = 0;

    SCAI_ASSERT_EQ_ERROR( m, n, "getrf only for square matrices" )
    SCAI_ASSERT_GE_ERROR( lda, m, "lda too small" )

    for ( IndexType i = 0; i < m; i++ )
    {
        ipiv[i] = i;
    }

    // a( i, j ) -> a ( i * lda + j )

    for ( IndexType i = 0; i < m; i++ )
    {
        // pivoting over ( i:n, i )

        index = i;

        for ( IndexType k = i + 1; k < m; k++ )
        {
            if ( common::Math::abs( a[lda * ipiv[k] + i] ) > common::Math::abs( a[lda * ipiv[index] + i] ) )
            {
                index = k;
            }
        }

        SCAI_LOG_TRACE( logger, "Step i = " << i << ", pivot is row " << index << ", val = " << a[lda * ipiv[index] + i ] )

        std::swap( ipiv[index], ipiv[i] );

        for ( IndexType k = i + 1; k < m; k++ )
        {
            // a( k, i ) = a( k, i ) / a( i, i )

            a[lda * ipiv[k] + i] /= a[lda * ipiv[i] + i];

            for ( IndexType j = i + 1; j < n; j++ )
            {
                // a( k, j ) = a( k , j ) - a( k, i ) * a( i, j )

                a[lda * ipiv[k] + j] -= a[lda * ipiv[k] + i] * a[lda * ipiv[i] + j];

                SCAI_LOG_TRACE( logger, "a[" << k << ", " << j << "] = " << a[lda * ipiv[k] + j ] )
            }
        }
    }
}

/* ------------------------------------------------------------------------- */
/*      getinv<float>                                                        */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPLAPACK::getinv( const IndexType n, ValueType* a, const IndexType lda )
{
    SCAI_REGION( "OpenMP.LAPACK.getinv<ValueType>" )
    SCAI_LOG_INFO( logger,
                   "getinv<" << TypeTraits<ValueType>::id() << "> for " << n << " x " << n << " matrix, uses openmp" )
    // temporary array for pivot indexes needed, deleted by destructor
    std::unique_ptr<IndexType[]> ipiv( new IndexType[n] );
    getrf( n, n, a, lda, ipiv.get() );
    getri( n, a, lda, ipiv.get() );
}

/* ------------------------------------------------------------------------- */
/*      getri<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPLAPACK::getri( const IndexType n, ValueType* const A, const IndexType lda, IndexType* const ipiv )
{
    SCAI_REGION( "OpenMP.LAPACK.getri<ValueType>" )

    SCAI_LOG_INFO( logger, "getri<" << TypeTraits<ValueType>::id() << "> for A of size " << n << " x " << n << ", lda = " << lda )

    const ValueType ZERO = 0;

    SCAI_ASSERT_GE_ERROR( lda, n, "lda too small" )

    std::unique_ptr<ValueType[]> A_inv( new ValueType[n * n] );

    for ( IndexType i = 0; i < n * n; i++ )
    {
        A_inv[i] = ZERO;
    }

    for ( IndexType i = 0; i < n; i++ )
    {
        //Emulation of an identity matrix
        for ( IndexType h = 0; h < n; h++ )
        {
            A_inv[n * ipiv[h] + i] = ZERO;
        }

        A_inv[n * i + i] = 1;

        //Forwared Elimination to solve Z, which is stored in A_inv
        for ( IndexType j = 0; j < n; j++ )
        {
            for ( IndexType k = 0; k < j; k++ )
            {
                A_inv[n * ipiv[j] + i] -= A[lda * ipiv[j] + k] * A_inv[n * ipiv[k] + i];
            }
        }

        //Back Substitution to solve x, which is stored in A_inv
        for ( IndexType j = n; j-- > 0; )
        {
            for ( IndexType k = j + 1; k < n; k++ )
            {
                A_inv[n * ipiv[j] + i] -= A[lda * ipiv[j] + k] * A_inv[n * ipiv[k] + i];
            }

            A_inv[n * ipiv[j] + i] /= A[lda * ipiv[j] + j];
        }
    }

    // copy A_inv back to A

    for ( IndexType i = 0; i < n; i++ )
    {
        for ( IndexType j = 0; j < n; j++ )
        {
            A[lda * j + i] = A_inv[n * ipiv[j] + i];
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPLAPACK::tptrs(
    const CBLAS_UPLO uplo,
    const common::MatrixOp op,
    const CBLAS_DIAG diag,
    const IndexType n,
    const IndexType SCAI_UNUSED( nrhs ),
    const ValueType* AP,
    ValueType* B,
    const IndexType SCAI_UNUSED( ldb ) )
{
    SCAI_REGION( "OpenMP.LAPACK.tptrs<ValueType>" )

    bool isColOrder = true;

    if ( isColOrder )
    {
        if ( op == common::MatrixOp::NORMAL )
        {
            if ( uplo == CblasUpper )
            {
                // Back substituion
                if ( diag == CblasUnit )
                {
                    COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs is not implemented for CblasUnit yet" )
                }
                else if ( diag == CblasNonUnit )
                {
                    for ( IndexType i = n; i-- > 0; )
                    {
                        for ( IndexType j = n - 1; j > i; j-- )
                        {
                            B[i] -= AP[( j * j + j ) / 2 + i] * B[j];
                        }

                        B[i] /= AP[( i * i + i ) / 2 + i];
                    }
                }
            }
            else if ( uplo == CblasLower )
            {
                //Forward Elimination
                if ( diag == CblasUnit )
                {
                    COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs is not implemented for CblasUnit yet" )
                }
                else if ( diag == CblasNonUnit )
                {
                    for ( IndexType i = 0; i < n; i++ )
                    {
                        for ( IndexType j = 0; j < i; j++ )
                        {
                            B[i] -= AP[i + ( ( ( 2 * n - j - 1 ) * j ) / 2 )] * B[j];
                        }

                        B[i] /= AP[i + ( ( ( 2 * n - i - 1 ) * i ) / 2 )];
                    }
                }
            }
        }
        else if ( op == common::MatrixOp::TRANSPOSE )
        {
            // Todo: implement
            COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs is not implemented for CblasTrans yet" )
            /*
            if( uplo == CblasUpper )
            {
                if( diag == CblasUnit )
                {

                }
                else if( diag == CblasNonUnit )
                {

                }
            }
            else if( uplo == CblasLower )
            {
                if( diag == CblasUnit )
                {

                }
                else if( diag == CblasNonUnit )
                {

                }
            }*/
        }
        else 
        {
            COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs - colmajor - op = " << op << " not implemented" );
        }
    }
    else 
    {
        // row-major order

        if ( op == common::MatrixOp::NORMAL )
        {
            if ( uplo == CblasUpper )
            {
                if ( diag == CblasUnit )
                {
                    COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs is not implemented for CblasUnit yet" )
                }
                else if ( diag == CblasNonUnit )
                {
                    //Back substition
                    for ( IndexType i = n; i-- > 0; )
                    {
                        for ( IndexType j = n - 1; j > i; j-- )
                        {
                            B[i] -= AP[( j + ( ( ( 2 * n - i - 1 ) * i ) / 2 ) )] * B[j];
                        }

                        B[i] /= AP[( i + ( ( ( 2 * n - i - 1 ) * i ) / 2 ) )];
                    }
                }
            }
            else if ( uplo == CblasLower )
            {
                if ( diag == CblasUnit )
                {
                    COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs is not implemented for CblasUnit yet" )
                }
                else if ( diag == CblasNonUnit )
                {
                    //Forward elimination
                    for ( IndexType i = 0; i < n; i++ )
                    {
                        for ( IndexType j = 0; j < i; j++ )
                        {
                            B[i] -= AP[( ( i * i + i ) / 2 ) + j] * B[j];
                        }

                        B[i] /= AP[( ( i * i + i ) / 2 ) + i];
                    }
                }
            }
        }
        else 
        {
            // Todo: implement
            COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs - rowmajor - op = " << op << " not implemented" );
        }
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPLAPACK::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register LAPACK OpenMP-routines for Host at kernel registry [" << flag << "]" )
    KernelRegistry::set<BLASKernelTrait::getrf<ValueType> >( OpenMPLAPACK::getrf, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::getri<ValueType> >( OpenMPLAPACK::getri, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::getinv<ValueType> >( OpenMPLAPACK::getinv, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::tptrs<ValueType> >( OpenMPLAPACK::tptrs, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPLAPACK::OpenMPLAPACK()
{
    SCAI_LOG_INFO( logger, "register LAPACK OpenMP-routines for Host at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ADD );
}

OpenMPLAPACK::~OpenMPLAPACK()
{
    SCAI_LOG_INFO( logger, "unregister LAPACK OpenMP-routines for Host at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPLAPACK OpenMPLAPACK::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
