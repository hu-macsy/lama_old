/**
 * @file LAPACKe_LAPACK.cpp
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
 * @brief Implementation of LAPACK kernel routines for Host by wrapping to LAPACKe library.
 * @author Lauretta Schubert
 * @date 02.07.2012
 */

// hpp
#include <scai/blaskernel/external/LAPACKe_LAPACK.hpp>

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/cblas.hpp>
#include <scai/blaskernel/external/LAPACKeWrapper.hpp>
#include <scai/blaskernel/external/LAPACKeTrait.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Settings.hpp>
#include <scai/tracing.hpp>
#include <scai/common/Math.hpp>

// external
#include <mkl_lapacke.h>

#include <memory>

using std::unique_ptr;

namespace scai
{

using common::TypeTraits;

namespace blaskernel
{

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( LAPACKe_LAPACK::logger, "LAPACKe.LAPACK" )

/* ------------------------------------------------------------------------- */
/*      getrf                                                                */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACKe_LAPACK::getrf( const IndexType m, const IndexType n, 
                            ValueType* const A, const IndexType lda,
                            IndexType* const ipiv )
{
    SCAI_REGION( "LAPACKe.getrf" )

    SCAI_LOG_INFO( logger, "getrf<" << common::TypeTraits<ValueType>::id() << "> for A of size " << m << " x " << n )
    typedef LAPACKeTrait::LAPACKIndexType LAPACKIndexType;

    if ( TypeTraits<IndexType>::stype != TypeTraits<LAPACKIndexType>::stype )
    {
        // ToDo: convert ipiv array

        COMMON_THROWEXCEPTION( "indextype mismatch, LAMA uses " << TypeTraits<IndexType>::id()
                               << ", LAPACK uses " << TypeTraits<LAPACKIndexType>::id()  );
    }

    LAPACKIndexType* la_ipiv = ( LAPACKIndexType* ) ipiv;

    int info = LAPACKeWrapper<ValueType>::getrf( LAPACK_ROW_MAJOR, 
               static_cast<LAPACKIndexType>( m ),
               static_cast<LAPACKIndexType>( n ), A,
               static_cast<LAPACKIndexType>( lda ), la_ipiv );

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
/*      getinv                                                               */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACKe_LAPACK::getinv( const IndexType n, ValueType* a,
                             const IndexType lda )
{
    SCAI_REGION( "LAPACKe.getinv" )

    typedef LAPACKeTrait::LAPACKIndexType LAPACKIndexType;

    unique_ptr<LAPACKIndexType[]> ipiv( new LAPACKIndexType[n] );  // freed by destructor

    SCAI_LOG_INFO( logger,
                   "getinv<float> for " << n << " x " << n << " matrix, uses MKL" )

    LAPACKIndexType info = 0;

    {
        SCAI_REGION( "LAPACKe.getrf" )

        info = LAPACKeWrapper<ValueType>::getrf( LAPACK_COL_MAJOR,
                static_cast<LAPACKIndexType>( n ),
                static_cast<LAPACKIndexType>( n ),
                a,
                static_cast<LAPACKIndexType>( lda ),
                ipiv.get() );
    }

    // throw exception if factorization did not work

    if ( info )
    {
        COMMON_THROWEXCEPTION( "LAPACKe getrf failed, info = " << info )
    }

    {
        SCAI_REGION( "LAPACKe.getri" )

        info = LAPACKeWrapper<ValueType>::getri( LAPACK_COL_MAJOR,
                static_cast<LAPACKIndexType>( n ),
                a,
                static_cast<LAPACKIndexType>( lda ),
                ipiv.get() );
    }

    if ( info )
    {
        COMMON_THROWEXCEPTION( "LAPACKe getri failed, info = " << info )
    }
}

/* ------------------------------------------------------------------------- */
/*      getri                                                                */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACKe_LAPACK::getri( const IndexType n,
                            ValueType* const a, const IndexType lda, IndexType* const ipiv )
{
    SCAI_REGION( "LAPACKe.getri" )

    typedef LAPACKeTrait::LAPACKIndexType LAPACKIndexType;
    SCAI_LOG_INFO( logger, "getri<float> for A of size " << n << " x " << n )

    if ( TypeTraits<IndexType>::stype != TypeTraits<LAPACKIndexType>::stype )
    {
        // ToDo: convert ipiv array
        COMMON_THROWEXCEPTION( "indextype mismatch" );
    }

    LAPACKeTrait::LAPACKOrder matrix_order = LAPACK_ROW_MAJOR;

    LAPACKIndexType info = LAPACKeWrapper<ValueType>::getri( matrix_order,
                           static_cast<LAPACKIndexType>( n ), a,
                           static_cast<LAPACKIndexType>( lda ), reinterpret_cast<LAPACKIndexType*>( ipiv ) );

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
/*      tptrs                                                                */
/* ------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACKe_LAPACK::tptrs( const CBLAS_UPLO uplo,
                            const common::MatrixOp op, const CBLAS_DIAG diag, const IndexType n,
                            const IndexType nrhs, const ValueType* AP, ValueType* B,
                            const IndexType ldb )
{
    SCAI_REGION( "LAPACKe.tptrs" )

    typedef LAPACKeTrait::LAPACKIndexType LAPACKIndexType;
    LAPACKeTrait::LAPACKFlag UL = LAPACKeTrait::enum2char( uplo );
    LAPACKeTrait::LAPACKFlag TA = LAPACKeTrait::enum2char( op );
    LAPACKeTrait::LAPACKFlag DI = LAPACKeTrait::enum2char( diag );
    LAPACKeTrait::LAPACKOrder matrix_order = LAPACK_ROW_MAJOR;

    if ( TypeTraits<IndexType>::stype != TypeTraits<LAPACKIndexType>::stype )
    {
        // ToDo: convert ipiv array
        COMMON_THROWEXCEPTION( "indextype mismatch" );
    }

    SCAI_LOG_INFO( logger, "tptrs<" << TypeTraits<ValueType>::id() << ">"
                   << ", n = " << n << ", nrhs = " << nrhs << ", order = " << matrix_order
                   << ", UL = " << UL << ", TA = " << TA << ", DI = " << DI );

    IndexType one = 1;

    SCAI_ASSERT_ERROR( ldb >= common::Math::max( one, nrhs ), "ldb = " << ldb << " out of range" );

    int info = LAPACKeWrapper<ValueType>::tptrs( matrix_order, UL, TA, DI,
               static_cast<LAPACKIndexType>( n ),
               static_cast<LAPACKIndexType>( nrhs ), AP, B,
               static_cast<LAPACKIndexType>( ldb ) );

    if ( info )
    {
        COMMON_THROWEXCEPTION( "MKL tptrs failed, info = " << info )
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void LAPACKe_LAPACK::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_INFO( logger, "register lapack wrapper routines for Host at kernel registry" )
    KernelRegistry::set<BLASKernelTrait::getrf<ValueType> >( LAPACKe_LAPACK::getrf, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::getri<ValueType> >( LAPACKe_LAPACK::getri, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::getinv<ValueType> >( LAPACKe_LAPACK::getinv, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::tptrs<ValueType> >( LAPACKe_LAPACK::tptrs, ctx, flag );
}

LAPACKe_LAPACK::LAPACKe_LAPACK()
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

LAPACKe_LAPACK::~LAPACKe_LAPACK()
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

LAPACKe_LAPACK LAPACKe_LAPACK::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
