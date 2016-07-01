/**
 * @file OpenMPLAPACK.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief OpenMPLAPACK.cpp
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

#include <scai/common/unique_ptr.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

// std
#include <cmath>

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
IndexType OpenMPLAPACK::getrf(
    const CBLAS_ORDER order,
    const IndexType m,
    const IndexType n,
    ValueType* const a,
    const IndexType lda,
    IndexType* const ipiv )
{
    SCAI_REGION( "OpenMP.LAPACK.getrf<ValueType>" )
    SCAI_LOG_INFO( logger, "getrf<" << TypeTraits<ValueType>::id() << "> for A of size " << m << " x " << n )
    IndexType info = 0;
    IndexType index = 0;

    if ( m != n || n != lda )
    {
        COMMON_THROWEXCEPTION( "OpenMPLAPACK:getrf only supported for square matrices" );
    }

    for ( IndexType i = 0; i < lda; i++ )
    {
        ipiv[i] = i;
    }

    if ( order == CblasColMajor )
    {
        for ( IndexType i = 0; i < lda; i++ )
        {
            //pivoting
            index = i;

            for ( IndexType j = i; j < lda; j++ )
            {
                if ( common::Math::abs( a[lda * ipiv[i] + j] ) > common::Math::abs( a[lda * ipiv[i] + index] ) )
                {
                    index = j;
                }
            }

            IndexType temp = ipiv[index];
            ipiv[index] = ipiv[i];
            ipiv[i] = temp;
        }

        for ( IndexType i = 0; i < lda; i++ )
        {
            for ( IndexType j = i + 1; j < lda; j++ )
            {
                a[lda * ipiv[i] + j] /= a[lda * ipiv[i] + i];

                for ( IndexType k = i + 1; k < lda; k++ )
                {
                    a[lda * ipiv[k] + j] -= a[lda * ipiv[i] + j] * a[lda * ipiv[k] + i];
                }
            }
        }
    }
    else if ( order == CblasRowMajor )
    {
        for ( IndexType i = 0; i < lda; i++ )
        {
            //pivoting
            index = i;

            for ( IndexType j = i; j < lda; j++ )
            {
                if ( common::Math::abs( a[lda * ipiv[j] + i] ) > common::Math::abs( a[lda * ipiv[index] + i] ) )
                {
                    index = j;
                }
            }

            IndexType temp = ipiv[index];
            ipiv[index] = ipiv[i];
            ipiv[i] = temp;
        }

        for ( IndexType i = 0; i < lda; i++ )
        {
            for ( IndexType j = i + 1; j < lda; j++ )
            {
                a[lda * ipiv[j] + i] /= a[lda * ipiv[i] + i];

                for ( IndexType k = i + 1; k < lda; k++ )
                {
                    a[lda * ipiv[j] + k] -= a[lda * ipiv[j] + i] * a[lda * ipiv[i] + k];
                }
            }
        }
    }

    return info;
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
    common::scoped_array<IndexType> ipiv( new IndexType[n] );
    getrf( CblasRowMajor, n, n, a, lda, ipiv.get() );
    getri( CblasRowMajor, n, a, lda, ipiv.get() );
}

/* ------------------------------------------------------------------------- */
/*      getri<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType OpenMPLAPACK::getri( const CBLAS_ORDER order, const IndexType n, ValueType* const A, const IndexType lda, IndexType* const ipiv )
{
    SCAI_REGION( "OpenMP.LAPACK.getri<ValueType>" )
    SCAI_LOG_INFO( logger, "getri<" << TypeTraits<ValueType>::id() << "> for A of size " << n << " x " << n )
    IndexType info = 0;
    ValueType* A_inv = new ValueType[n * n];

    for ( IndexType i = 0; i < n * n; i++ )
    {
        A_inv[i] = static_cast<ValueType>( 0.0 );
    }

    if ( order == CblasRowMajor )
    {
        for ( IndexType i = 0; i < n; i++ )
        {
            //Emulation of an identity matrix
            for ( IndexType h = 0; h < n; h++ )
            {
                A_inv[n * ipiv[h] + i] = 0;
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
            for ( IndexType j = n - 1; j >= 0; j-- )
            {
                for ( IndexType k = j + 1; k < n; k++ )
                {
                    A_inv[n * ipiv[j] + i] -= A[lda * ipiv[j] + k] * A_inv[n * ipiv[k] + i];
                }

                A_inv[n * ipiv[j] + i] /= A[n * ipiv[j] + j];
            }
        }

        for ( IndexType i = 0; i < n; i++ )
        {
            //Emulation of an identity matrix
            for ( IndexType j = 0; j < n; j++ )
            {
                A[n * j + i] = A_inv[n * ipiv[j] + i];
            }
        }
    }
    else if ( order == CblasColMajor )
    {
        for ( IndexType i = 0; i < n; i++ )
        {
            //Columns of Matrix A^-1
            for ( IndexType h = 0; h < n; h++ )
            {
                A_inv[n * ipiv[i] + h] = static_cast<ValueType>( 0.0 );
            }

            A_inv[n * i + i] = 1;

            //Forwared Elimination to solve Z, which is stored in A_inv
            for ( IndexType j = 0; j < n; j++ )
            {
                for ( IndexType k = 0; k < j; k++ )
                {
                    A_inv[n * ipiv[i] + j] -= A[lda * ipiv[k] + j] * A_inv[n * ipiv[i] + k];
                }
            }

            //Back Substitution to solve x, which is stored in A_inv
            for ( IndexType j = n - 1; j >= 0; j-- )
            {
                for ( IndexType k = j + 1; k < n; k++ )
                {
                    A_inv[n * ipiv[i] + j] -= A[lda * ipiv[k] + j] * A_inv[n * ipiv[i] + k];
                }

                A_inv[n * ipiv[i] + j] /= A[n * ipiv[j] + j];
            }
        }

        for ( IndexType i = 0; i < n; i++ )
        {
            //Copy from A_inv to A
            for ( IndexType j = 0; j < n; j++ )
            {
                A[n * i + j] = A_inv[n * ipiv[i] + j];
            }
        }
    }

    delete[] A_inv;

    if ( info < 0 )
    {
        COMMON_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if ( info > 0 )
    {
        COMMON_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
    }

    return info;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType OpenMPLAPACK::tptrs(
    const CBLAS_ORDER order,
    const CBLAS_UPLO uplo,
    const CBLAS_TRANSPOSE trans,
    const CBLAS_DIAG diag,
    const IndexType n,
    const IndexType SCAI_UNUSED( nrhs ),
    const ValueType* AP,
    ValueType* B,
    const IndexType SCAI_UNUSED( ldb ) )
{
    SCAI_REGION( "OpenMP.LAPACK.tptrs<ValueType>" )
    IndexType info = 0;

    if ( order == CblasColMajor )
    {
        if ( trans == CblasNoTrans )
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
                    for ( IndexType i = n - 1; i >= 0; i-- )
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
        else if ( trans == CblasTrans )
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
        else if ( trans == CblasConjTrans )
        {
            COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs - colmajor - conjtrans not implemented" );
        }
    }
    else if ( order == CblasRowMajor )
    {
        if ( trans == CblasNoTrans )
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
                    for ( IndexType i = n - 1; i >= 0; i-- )
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
        else if ( trans == CblasTrans )
        {
            // Todo: implement
            COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs - colmajor - conjtrans not implemented" );
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
        else if ( trans == CblasConjTrans )
        {
            COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs - rowmajor - conjtrans not implemented" );
        }
    }

    return info;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPLAPACK::laswp(
    const CBLAS_ORDER order,
    const IndexType N,
    ValueType* A,
    const IndexType LDA,
    const IndexType K1,
    const IndexType K2,
    const IndexType* ipiv,
    const IndexType INCX )
{
    SCAI_REGION( "OpenMP.LAPACK.laswp" )
    IndexType i = K1;

    if ( order == CblasRowMajor )
    {
        for ( i = K1; i < K2; ++i )
        {
            if ( ipiv[i * INCX] == i )
            {
                continue;
            }

            OpenMPBLAS1::swap<ValueType>( N, &A[ipiv[i * INCX] * LDA], INCX, &A[i * LDA], INCX );
        }
    }
    else if ( order == CblasColMajor )
    {
        for ( i = K1; i < K2; ++i )
        {
            if ( ipiv[i * INCX] == i )
            {
                continue;
            }

            OpenMPBLAS1::swap<ValueType>( N, &A[ipiv[i * INCX] * LDA], INCX, &A[i * LDA], INCX );
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "illegal order setting " << order )
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPLAPACK::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::Host;
    SCAI_LOG_DEBUG( logger, "register LAPACK OpenMP-routines for Host at kernel registry [" << flag << "]" )
    KernelRegistry::set<BLASKernelTrait::getrf<ValueType> >( OpenMPLAPACK::getrf, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::getri<ValueType> >( OpenMPLAPACK::getri, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::getinv<ValueType> >( OpenMPLAPACK::getinv, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::tptrs<ValueType> >( OpenMPLAPACK::tptrs, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::laswp<ValueType> >( OpenMPLAPACK::laswp, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPLAPACK::OpenMPLAPACK()
{
    SCAI_LOG_INFO( logger, "register LAPACK OpenMP-routines for Host at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_HOST_LIST>::call(
        kregistry::KernelRegistry::KERNEL_ADD );
}

OpenMPLAPACK::~OpenMPLAPACK()
{
    SCAI_LOG_INFO( logger, "unregister LAPACK OpenMP-routines for Host at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_HOST_LIST>::call(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPLAPACK OpenMPLAPACK::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
