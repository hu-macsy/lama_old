/**
 * @file OpenMPLAPACK.cpp
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
 * @brief OpenMPLAPACK.cpp
 * @author Eric Schricker
 * @date 19.11.2013
 * @since 1.0.0
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
#include <scai/common/mepr/Container.hpp>

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
    const int m,
    const int n,
    ValueType* const a,
    const int lda,
    int* const ipiv )
{
    SCAI_REGION( "OpenMP.LAPACK.getrf<ValueType>" )

    SCAI_LOG_INFO( logger, "getrf<" << TypeTraits<ValueType>::id()<< "> for A of size " << m << " x " << n )

    int info = 0;
    int index = 0;

    if( m != n || n != lda )
    {
        COMMON_THROWEXCEPTION( "OpenMPLAPACK:getrf only supported for square matrices" );
    }

    for( int i = 0; i < lda; i++ )
    {
        ipiv[i] = i;
    }

    if( order == CblasColMajor )
    {
        for( int i = 0; i < lda; i++ )
        {
            //pivoting
            index = i;

            for( int j = i; j < lda; j++ )
            {
                if( common::Math::abs( a[lda * ipiv[i] + j] ) > common::Math::abs( a[lda * ipiv[i] + index] ) )
                {
                    index = j;
                }
            }

            int temp = ipiv[index];
            ipiv[index] = ipiv[i];
            ipiv[i] = temp;

        }

        for( int i = 0; i < lda; i++ )
        {
            for( int j = i + 1; j < lda; j++ )
            {
                a[lda * ipiv[i] + j] /= a[lda * ipiv[i] + i];

                for( int k = i + 1; k < lda; k++ )
                {
                    a[lda * ipiv[k] + j] -= a[lda * ipiv[i] + j] * a[lda * ipiv[k] + i];
                }
            }

        }

    }
    else if( order == CblasRowMajor )
    {
        for( int i = 0; i < lda; i++ )
        {
            //pivoting
            index = i;

            for( int j = i; j < lda; j++ )
            {
                if( common::Math::abs( a[lda * ipiv[j] + i] ) > common::Math::abs( a[lda * ipiv[index] + i] ) )
                {
                    index = j;
                }
            }

            int temp = ipiv[index];
            ipiv[index] = ipiv[i];
            ipiv[i] = temp;

        }

        for( int i = 0; i < lda; i++ )
        {
            for( int j = i + 1; j < lda; j++ )
            {
                a[lda * ipiv[j] + i] /= a[lda * ipiv[i] + i];

                for( int k = i + 1; k < lda; k++ )
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
                   "getinv<" << TypeTraits<ValueType>::id()<< "> for " << n << " x " << n << " matrix, uses openmp" )

    // temporary array for pivot indexes needed, deleted by destructor

    common::scoped_array<int> ipiv( new int[n] );

    getrf( CblasRowMajor, n, n, a, lda, ipiv.get() );
    getri( CblasRowMajor, n, a, lda, ipiv.get() );
}

/* ------------------------------------------------------------------------- */
/*      getri<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
int OpenMPLAPACK::getri( const CBLAS_ORDER order, const int n, ValueType* const A, const int lda, int* const ipiv )
{
    SCAI_REGION( "OpenMP.LAPACK.getri<ValueType>" )

    SCAI_LOG_INFO( logger, "getri<" << TypeTraits<ValueType>::id()<< "> for A of size " << n << " x " << n )
    int info = 0;
    ValueType* A_inv = new ValueType[n * n];

    for( int i = 0; i < n * n; i++ )
    {
        A_inv[i] = static_cast<ValueType>(0.0);
    }

    if( order == CblasRowMajor )
    {

        for( int i = 0; i < n; i++ )
        {
            //Emulation of an identity matrix
            for( int h = 0; h < n; h++ )
            {
                A_inv[n * ipiv[h] + i] = 0;
            }

            A_inv[n * i + i] = 1;

            //Forwared Elimination to solve Z, which is stored in A_inv
            for( int j = 0; j < n; j++ )
            {
                for( int k = 0; k < j; k++ )
                {
                    A_inv[n * ipiv[j] + i] -= A[lda * ipiv[j] + k] * A_inv[n * ipiv[k] + i];
                }
            }

            //Back Substitution to solve x, which is stored in A_inv
            for( int j = n - 1; j >= 0; j-- )
            {
                for( int k = j + 1; k < n; k++ )
                {
                    A_inv[n * ipiv[j] + i] -= A[lda * ipiv[j] + k] * A_inv[n * ipiv[k] + i];
                }

                A_inv[n * ipiv[j] + i] /= A[n * ipiv[j] + j];
            }
        }

        for( int i = 0; i < n; i++ )
        {
            //Emulation of an identity matrix
            for( int j = 0; j < n; j++ )
            {
                A[n * j + i] = A_inv[n * ipiv[j] + i];
            }
        }
    }
    else if( order == CblasColMajor )
    {
        for( int i = 0; i < n; i++ )
        {
            //Columns of Matrix A^-1
            for( int h = 0; h < n; h++ )
            {
                A_inv[n * ipiv[i] + h] = static_cast<ValueType>(0.0);
            }

            A_inv[n * i + i] = 1;

            //Forwared Elimination to solve Z, which is stored in A_inv
            for( int j = 0; j < n; j++ )
            {
                for( int k = 0; k < j; k++ )
                {
                    A_inv[n * ipiv[i] + j] -= A[lda * ipiv[k] + j] * A_inv[n * ipiv[i] + k];
                }
            }

            //Back Substitution to solve x, which is stored in A_inv
            for( int j = n - 1; j >= 0; j-- )
            {
                for( int k = j + 1; k < n; k++ )
                {
                    A_inv[n * ipiv[i] + j] -= A[lda * ipiv[k] + j] * A_inv[n * ipiv[i] + k];
                }

                A_inv[n * ipiv[i] + j] /= A[n * ipiv[j] + j];
            }
        }

        for( int i = 0; i < n; i++ )
        {
            //Copy from A_inv to A
            for( int j = 0; j < n; j++ )
            {
                A[n * i + j] = A_inv[n * ipiv[i] + j];
            }
        }
    }

    delete[] A_inv;

    if( info < 0 )
    {
        COMMON_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if( info > 0 )
    {
        COMMON_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
    }

    return info;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
int OpenMPLAPACK::tptrs(
    const CBLAS_ORDER order,
    const CBLAS_UPLO uplo,
    const CBLAS_TRANSPOSE trans,
    const CBLAS_DIAG diag,
    const int n,
    const int UNUSED(nrhs),
    const ValueType* AP,
    ValueType* B,
    const int UNUSED(ldb) )
{
    SCAI_REGION( "OpenMP.LAPACK.tptrs<ValueType>" )

    int info = 0;

    if( order == CblasColMajor )
    {
        if( trans == CblasNoTrans )
        {
            if( uplo == CblasUpper )
            {
                // Back substituion
                if( diag == CblasUnit )
                {
                    COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs is not implemented for CblasUnit yet" )
                }
                else if( diag == CblasNonUnit )
                {

                    for( int i = n - 1; i >= 0; i-- )
                    {
                        for( int j = n - 1; j > i; j-- )
                        {
                            B[i] -= AP[( j * j + j ) / 2 + i] * B[j];
                        }

                        B[i] /= AP[( i * i + i ) / 2 + i];
                    }
                }
            }
            else if( uplo == CblasLower )
            {
                //Forward Elimination
                if( diag == CblasUnit )
                {
                    COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs is not implemented for CblasUnit yet" )
                }
                else if( diag == CblasNonUnit )
                {
                    for( int i = 0; i < n; i++ )
                    {
                        for( int j = 0; j < i; j++ )
                        {
                            B[i] -= AP[i + ( ( ( 2 * n - j - 1 ) * j ) / 2 )] * B[j];
                        }

                        B[i] /= AP[i + ( ( ( 2 * n - i - 1 ) * i ) / 2 )];
                    }
                }
            }
        }
        else if( trans == CblasTrans )
        {
            COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs is not implemented for CblasTrans yet" )

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
            }
        }
        else if( trans == CblasConjTrans )
        {
            COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs - colmajor - conjtrans not implemented" );
        }
    }
    else if( order == CblasRowMajor )
    {
        if( trans == CblasNoTrans )
        {
            if( uplo == CblasUpper )
            {
                if( diag == CblasUnit )
                {
                    COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs is not implemented for CblasUnit yet" )
                }
                else if( diag == CblasNonUnit )
                {
                    //Back substition
                    for( int i = n - 1; i >= 0; i-- )
                    {
                        for( int j = n - 1; j > i; j-- )
                        {
                            B[i] -= AP[( j + ( ( ( 2 * n - i - 1 ) * i ) / 2 ) )] * B[j];
                        }

                        B[i] /= AP[( i + ( ( ( 2 * n - i - 1 ) * i ) / 2 ) )];
                    }
                }
            }
            else if( uplo == CblasLower )
            {
                if( diag == CblasUnit )
                {
                    COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs is not implemented for CblasUnit yet" )
                }
                else if( diag == CblasNonUnit )
                {
                    //Forward elimination
                    for( int i = 0; i < n; i++ )
                    {
                        for( int j = 0; j < i; j++ )
                        {
                            B[i] -= AP[( ( i * i + i ) / 2 ) + j] * B[j];
                        }

                        B[i] /= AP[( ( i * i + i ) / 2 ) + i];
                    }
                }
            }
        }
        else if( trans == CblasTrans )
        {
            COMMON_THROWEXCEPTION( "OpenMPLAPACK:tptrs - colmajor - conjtrans not implemented" );

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
            }
        }
        else if( trans == CblasConjTrans )
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
    const int N,
    ValueType* A,
    const int LDA,
    const int K1,
    const int K2,
    const int* ipiv,
    const int INCX )
{
    SCAI_REGION( "OpenMP.LAPACK.laswp" )

    int i = K1;

    if ( order == CblasRowMajor )
    {
        for( i = K1; i < K2; ++i )
        {
            if( ipiv[i * INCX] == i )
            {
                continue;
            }

            OpenMPBLAS1::swap<ValueType>( N, &A[ipiv[i * INCX] * LDA], INCX, &A[i * LDA], INCX );
        }
    }
    else if ( order == CblasColMajor )
    {
        for( i = K1; i < K2; ++i )
        {
            if( ipiv[i * INCX] == i )
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
    using common::context::Host;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register LAPACK OpenMP-routines for Host at kernel registry [" << flag << "]" )

    KernelRegistry::set<BLASKernelTrait::getrf<ValueType> >( OpenMPLAPACK::getrf, Host, flag );
    KernelRegistry::set<BLASKernelTrait::getri<ValueType> >( OpenMPLAPACK::getri, Host, flag );
    KernelRegistry::set<BLASKernelTrait::getinv<ValueType> >( OpenMPLAPACK::getinv, Host, flag );
    KernelRegistry::set<BLASKernelTrait::tptrs<ValueType> >( OpenMPLAPACK::tptrs, Host, flag );
    KernelRegistry::set<BLASKernelTrait::laswp<ValueType> >( OpenMPLAPACK::laswp, Host, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPLAPACK::OpenMPLAPACK()
{
    typedef common::mepr::ContainerV<RegistratorV, ARITHMETIC_HOST> ValueTypes;

    kregistry::instantiate( kregistry::KernelRegistry::KERNEL_ADD, ValueTypes() );
}

OpenMPLAPACK::~OpenMPLAPACK()
{
    typedef common::mepr::ContainerV<RegistratorV, ARITHMETIC_HOST> ValueTypes;

    kregistry::instantiate( kregistry::KernelRegistry::KERNEL_ERASE, ValueTypes() );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPLAPACK OpenMPLAPACK::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
