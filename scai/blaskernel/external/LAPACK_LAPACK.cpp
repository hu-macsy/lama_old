/**
 * @file LAPACK_LAPACK.cpp
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
 * @brief Implementation of LAPACK kernels for LAMA by using (Fortran) LAPACK.
 * @author Lauretta Schubert
 * @date 02.07.2012
 */

// hpp
#include <scai/blaskernel/external/LAPACK_LAPACK.hpp>

// local library
#include <scai/blaskernel/external/BLAS_BLAS1.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/external/LAPACKTrait.hpp>
#include <scai/blaskernel/external/LAPACKWrapper.hpp>
#include <scai/blaskernel/openmp/OpenMPLAPACK.hpp>
#include <scai/blaskernel/cblas.hpp>

// scai libraries
#include <scai/hmemo/Context.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>

#include <memory>

namespace scai
{

using std::unique_ptr;

namespace blaskernel
{

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( LAPACK_LAPACK::logger, "LAPACK.LAPACK" )

/* ------------------------------------------------------------------------- */
/*      getrf<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACK_LAPACK::getrf( const IndexType m,
                           const IndexType n, ValueType* const A, const IndexType lda,
                           IndexType* const ipiv )
{
    SCAI_REGION( "LAPACK.LAPACK.getrf<float>" )
    SCAI_LOG_INFO( logger, "getrf<float> for A of size " << m << " x " << n )
    typedef LAPACKTrait::LAPACKIndexType LAPACKIndexType;

    if ( common::TypeTraits<IndexType>::stype
            != common::TypeTraits<LAPACKIndexType>::stype )
    {
        // ToDo: convert ipiv array
        COMMON_THROWEXCEPTION( "indextype mismatch" );
    }

    LAPACKIndexType info = 0;

    if ( m == n )
    {
        for ( IndexType i = 0; i < m; ++i )
        {
            for ( IndexType j = i + 1; j < n; ++j )
            {
                std::swap( A[i * lda + j], A[j * lda + i] );
            }
        }

        info = LAPACKWrapper<ValueType>::getrf( static_cast<LAPACKIndexType>( m ),
                                                static_cast<LAPACKIndexType>( n ), A,
                                                static_cast<LAPACKIndexType>( lda ), ipiv );

        for ( IndexType i = 0; i < m; ++i )
        {
            for ( IndexType j = i + 1; j < n; ++j )
            {
                std::swap( A[i * lda + j], A[j * lda + i] );
            }
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "row major only supported for square matrices, current is " << m << " x " << n );
    }

    for ( IndexType i = 0; i < m; ++i )
    {
        --ipiv[i]; // Fortran numbering from 1 to n ->  0 to n-1
    }

    if ( info < 0 )
    {
        COMMON_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if ( info > 0 )
    {
        COMMON_THROWEXCEPTION(
            "value(" << info << "," << info << ")" << " is exactly zero" )
    }
}

/* ------------------------------------------------------------------------- */
/*      getinv<float>                                                        */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACK_LAPACK::getinv( const IndexType n, ValueType* a,
                            const IndexType lda )
{
    SCAI_REGION( "LAPACK.LAPACK.getinv<float>" )
    typedef LAPACKTrait::LAPACKIndexType LAPACKIndexType;
    LAPACKIndexType info = 0;
    // unique_ptr, delete by destructor, also done in case of exception
    std::unique_ptr<IndexType[]> ipiv( new IndexType[n] );
    SCAI_LOG_INFO( logger,
                   "getinv<float> for " << n << " x " << n << " matrix, uses Fortran interface" )
    info = LAPACKWrapper<ValueType>::getrf( static_cast<LAPACKIndexType>( n ),
                                            static_cast<LAPACKIndexType>( n ), a,
                                            static_cast<LAPACKIndexType>( lda ), ipiv.get() );

    if ( info )
    {
        COMMON_THROWEXCEPTION( "LAPACK sgetrf failed, info = " << info )
    }

    std::unique_ptr<ValueType[]> work( new ValueType[n] );
    info = LAPACKWrapper<ValueType>::getri( static_cast<LAPACKIndexType>( n ), a,
                                            static_cast<LAPACKIndexType>( lda ), ipiv.get(), work.get(),
                                            static_cast<LAPACKIndexType>( n ) );

    if ( info )
    {
        COMMON_THROWEXCEPTION( "LAPACK sgetri failed, info = " << info )
    }
}

/* ------------------------------------------------------------------------- */
/*      getri<float>                                                         */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACK_LAPACK::getri( const IndexType n,
                           ValueType* const a, const IndexType lda, IndexType* const ipiv )
{
    SCAI_REGION( "LAPACK.LAPACK.getri<float>" )
    SCAI_LOG_INFO( logger, "getri<float> for A of size " << n << " x " << n )
    typedef LAPACKTrait::LAPACKIndexType LAPACKIndexType;

    if ( common::TypeTraits<IndexType>::stype
            != common::TypeTraits<LAPACKIndexType>::stype )
    {
        // ToDo: convert ipiv array
        COMMON_THROWEXCEPTION( "indextype mismatch" );
    }

    LAPACKIndexType info = 0;

    // translate C indexes into  Fortran Indexes for ipiv

    for ( IndexType i = 0; i < n; ++i )
    {
        ++ipiv[i];
    }

    for ( IndexType i = 0; i < n; ++i )
    {
        // swap row and column
        for ( IndexType j = i + 1; j < n; ++j )
        {
            std::swap( a[i * lda + j], a[j * lda + i] );
        }
    }

    std::unique_ptr<ValueType[]> work( new ValueType[n] );
    info = LAPACKWrapper<ValueType>::getri( static_cast<LAPACKIndexType>( n ), a,
                                            static_cast<LAPACKIndexType>( lda ), ipiv, work.get(),
                                            static_cast<LAPACKIndexType>( n ) );

    // transpose back
    for ( IndexType i = 0; i < n; ++i )
    {
        for ( IndexType j = i + 1; j < n; ++j )
        {
            std::swap( a[i * lda + j], a[j * lda + i] );
        }
    }

    if ( info < 0 )
    {
        COMMON_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if ( info > 0 )
    {
        COMMON_THROWEXCEPTION(
            "value(" << info << "," << info << ")" << " is exactly zero" )
    }
}

template<typename ValueType>
void LAPACK_LAPACK::tptrs( const CBLAS_UPLO uplo,
                           const common::MatrixOp op, const CBLAS_DIAG diag, const IndexType n,
                           const IndexType nrhs, const ValueType* AP, ValueType* B,
                           const IndexType ldb )
{
    SCAI_REGION( "LAPACK.LAPACK.tptrs" )

    typedef LAPACKTrait::LAPACKIndexType LAPACKIndexType;
    typedef LAPACKTrait::LAPACKFlag LAPACKFlag;

    LAPACKFlag UL = LAPACKTrait::enum2char( uplo );
    LAPACKFlag TA = LAPACKTrait::enum2char( op );
    LAPACKFlag DI = LAPACKTrait::enum2char( diag );

    SCAI_LOG_INFO( logger,
                   "tptrs<" << common::TypeTraits<ValueType>::id() << ">, n = " << n
                   << ", nrhs = " << nrhs << ", UL = " << UL << ", TA = " << TA << ", DI = " << DI );

    bool colOrder = true;

    if ( colOrder )
    {
        LAPACKWrapper<ValueType>::tptrs( UL, TA, DI,
                                         static_cast<LAPACKIndexType>( n ),
                                         static_cast<LAPACKIndexType>( nrhs ), AP, B,
                                         static_cast<LAPACKIndexType>( ldb ) );
    }
    else
    {
        // row-major order not supported by LAPACK, call own implementation
        OpenMPLAPACK::tptrs( uplo, op, diag, n, nrhs, AP, B, ldb );
    }
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the LAPACK routines                               */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACK_LAPACK::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::Host;
    KernelRegistry::set<BLASKernelTrait::getrf<ValueType> >( LAPACK_LAPACK::getrf, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::getri<ValueType> >( LAPACK_LAPACK::getri, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::getinv<ValueType> >( LAPACK_LAPACK::getinv, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::tptrs<ValueType> >( LAPACK_LAPACK::tptrs, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

LAPACK_LAPACK::LAPACK_LAPACK()
{
    bool useLAPACK = false;
    common::Settings::getEnvironment( useLAPACK, "SCAI_USE_LAPACK" );

    if ( !useLAPACK )
    {
        return;
    }

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_EXT_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_REPLACE );
}

LAPACK_LAPACK::~LAPACK_LAPACK()
{
    bool useLAPACK = false;
    common::Settings::getEnvironment( useLAPACK, "SCAI_USE_LAPACK" );

    if ( !useLAPACK )
    {
        return;
    }

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_EXT_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

LAPACK_LAPACK LAPACK_LAPACK::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
