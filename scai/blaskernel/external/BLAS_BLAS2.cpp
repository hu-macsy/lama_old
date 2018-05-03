/**
 * @file BLAS_BLAS2.cpp
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
 * @brief BLAS_BLAS2.cpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 */

// hpp
#include <scai/blaskernel/external/BLAS_BLAS2.hpp>

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/external/BLASWrapper.hpp>
#include <scai/blaskernel/openmp/OpenMPBLAS2.hpp>
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>

namespace scai
{

using namespace tasking;    // for running kernels asynchronously

namespace blaskernel
{

using common::TypeTraits;

SCAI_LOG_DEF_LOGGER( BLAS_BLAS2::logger, "BLAS.BLAS2" )

/* ---------------------------------------------------------------------------------------*/
/*    gemv                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void BLAS_BLAS2::gemv(
    const common::MatrixOp opA,
    const IndexType m,
    const IndexType n,
    const ValueType alpha,
    const ValueType* a,
    const IndexType lda,
    const ValueType* x,
    const IndexType incX,
    const ValueType beta,
    ValueType* y,
    const IndexType incY )
{
    SCAI_REGION( "BLAS.BLAS2.gemv" )
    SCAI_LOG_INFO( logger,
                   "gemv<" << TypeTraits<ValueType>::id() << ">: " << "m = " << m << ", n = " << n << ", LDA = " << lda << ", incX = " << incX << ", incY = " << incY << ", alpha = " << alpha << ", beta = " << beta )

    if ( m == 0 )
    {
        return; // empty X, Y, A
    }

    if ( common::isComplex( TypeTraits<ValueType>::stype ) && opA == common::MatrixOp::CONJ_TRANSPOSE )
    {
        // this case is not handled by Fortran BLAS, so we call own implementation

        OpenMPBLAS2::gemv( opA, m, n, alpha, a, lda, x, incX, beta, y, incY );
        return;
    }

    // N == 0: empty A, but deal with X, Y, we can handle this here

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution not supported yet." )
    }

    // BLAS routines assume data in column-major format, so we work on transposded data

    BLASTrait::BLASTrans b_opA = BLASTrait::castTrans( common::combine( opA, common::MatrixOp::TRANSPOSE ) );

    auto b_n = static_cast<BLASTrait::BLASIndexType>( n );
    auto b_m = static_cast<BLASTrait::BLASIndexType>( m );

    auto b_lda = static_cast<BLASTrait::BLASIndexType>( lda );
    auto b_incX = static_cast<BLASTrait::BLASIndexType>( incX );
    auto b_incY = static_cast<BLASTrait::BLASIndexType>( incY );

    BLASWrapper<ValueType>::gemv( b_opA, b_n, b_m,
                                  alpha, a, b_lda,
                                  x, b_incX,
                                  beta, y, b_incY );

    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void BLAS_BLAS2::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::Host;
    bool useBLAS = false;
    int level = 0;
    useBLAS = common::Settings::getEnvironment( level, "SCAI_USE_BLAS" );

    if ( !useBLAS || ( level <= 0 ) )
    {
        SCAI_LOG_INFO( logger, "BLAS2 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS not set or 0)" )
        return;
    }
    else if ( level > 2 )
    {
        // only level 3 wrappers might be used

        SCAI_LOG_INFO( logger,
                       "BLAS2 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS = " << level << ")" )
        return;
    }

    SCAI_LOG_DEBUG( logger, "register[" << flag << "] BLAS2 wrapper routines for Host at kernel registry: " <<
                    "T = " << common::TypeTraits<ValueType>::id() )

    KernelRegistry::set<BLASKernelTrait::gemv<ValueType> >( BLAS_BLAS2::gemv, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

BLAS_BLAS2::BLAS_BLAS2()
{
    SCAI_LOG_INFO( logger, "register BLAS2 wrapper routines for Host at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_EXT_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_REPLACE );
}

BLAS_BLAS2::~BLAS_BLAS2()
{
    SCAI_LOG_INFO( logger, "unregister BLAS2 wrapper routines for Host at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_EXT_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

BLAS_BLAS2 BLAS_BLAS2::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
