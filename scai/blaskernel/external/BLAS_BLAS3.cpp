/**
 * @file BLAS_BLAS3.cpp
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
 * @brief BLAS_BLAS3.cpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 */

// hpp
#include <scai/blaskernel/external/BLAS_BLAS3.hpp>

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/external/BLASWrapper.hpp>
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>

namespace scai
{

namespace blaskernel
{

using common::TypeTraits;
using tasking::TaskSyncToken;

SCAI_LOG_DEF_LOGGER( BLAS_BLAS3::logger, "BLAS.BLAS3" )

template<typename ValueType>
void BLAS_BLAS3::gemm(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE transA,
    const CBLAS_TRANSPOSE transB,
    const IndexType m,
    const IndexType n,
    const IndexType k,
    const ValueType alpha,
    const ValueType* A,
    const IndexType lda,
    const ValueType* B,
    const IndexType ldb,
    const ValueType beta,
    ValueType* C,
    const IndexType ldc )
{
    SCAI_REGION( "BLAS.BLAS3.gemm" )
    SCAI_LOG_INFO( logger,
                   "gemm<" << TypeTraits<ValueType>::id() << ">: " << "m = " << m << ", n = " << n << ", k = " << k
                   << ", lda = " << lda << ", ldb = " << ldb << ", ldc = " << ldc << ", alpha = " << alpha << ", beta = " << beta )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution not supported yet" )
    }

    BLASTrait::BLASTrans ta = '-', tb = '-';

    switch ( transA )
    {
        case CblasNoTrans:
            ta = 'N';
            break;

        case CblasTrans:
            ta = 'T';
            break;

        case CblasConjTrans:
            ta = 'C';
            break;
    }

    switch ( transB )
    {
        case CblasNoTrans:
            tb = 'N';
            break;

        case CblasTrans:
            tb = 'T';
            break;

        case CblasConjTrans:
            tb = 'C';
            break;
    }

    if ( order == CblasColMajor )
    {
        BLASWrapper<ValueType>::gemm( ta, tb, static_cast<BLASTrait::BLASIndexType>( m ), static_cast<BLASTrait::BLASIndexType>( n ), static_cast<BLASTrait::BLASIndexType>( k ), alpha, A, static_cast<BLASTrait::BLASIndexType>( lda ), B, static_cast<BLASTrait::BLASIndexType>( ldb ), beta, C, static_cast<BLASTrait::BLASIndexType>( ldc ) );
    }
    else if ( order == CblasRowMajor )
    {
        BLASWrapper<ValueType>::gemm( tb, ta, static_cast<BLASTrait::BLASIndexType>( n ), static_cast<BLASTrait::BLASIndexType>( m ), static_cast<BLASTrait::BLASIndexType>( k ), alpha, B, static_cast<BLASTrait::BLASIndexType>( ldb ), A, static_cast<BLASTrait::BLASIndexType>( lda ), beta, C, static_cast<BLASTrait::BLASIndexType>( ldc ) );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void BLAS_BLAS3::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::Host;
    bool useBLAS = false;
    int level = 0;
    useBLAS = common::Settings::getEnvironment( level, "SCAI_USE_BLAS" );

    if ( !useBLAS || ( level <= 0 ) )
    {
        SCAI_LOG_INFO( logger, "BLAS3 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS not set or 0)" )
        return;
    }
    else if ( level > 3 )
    {
        // only level 2 or level 3 wrappers might be used
        SCAI_LOG_INFO( logger,
                       "BLAS3 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS = " << level << ")" )
        return;
    }

    SCAI_LOG_DEBUG( logger, "register[" << flag << "] BLAS3 wrapper routines for Host at kernel registry: " <<
                            "T = " << common::TypeTraits<ValueType>::id() )

    KernelRegistry::set<BLASKernelTrait::gemm<ValueType> >( BLAS_BLAS3::gemm, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

BLAS_BLAS3::BLAS_BLAS3()
{
    SCAI_LOG_INFO( logger, "register BLAS3 wrapper routines for Host at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_EXT_HOST_LIST>::call(
        kregistry::KernelRegistry::KERNEL_REPLACE );
}

BLAS_BLAS3::~BLAS_BLAS3()
{
    SCAI_LOG_INFO( logger, "unregister BLAS3 wrapper routines for Host at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_EXT_HOST_LIST>::call(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

BLAS_BLAS3 BLAS_BLAS3::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
