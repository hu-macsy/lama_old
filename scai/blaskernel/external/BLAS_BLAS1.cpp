/**
 * @file BLAS_BLAS1.cpp
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
 * @brief Implementation of LAMA BLAS1 routines using wrappers to BLAS library
 * @author Lauretta Schubert
 * @date 05.07.2012
 */

// hpp
#include <scai/blaskernel/external/BLAS_BLAS1.hpp>

// local library
#include <scai/blaskernel/external/BLASWrapper.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>

#include <functional>

namespace scai
{

using namespace common;
using namespace tasking;

namespace blaskernel
{

SCAI_LOG_DEF_LOGGER( BLAS_BLAS1::logger, "BLAS.BLAS1" )

/* ---------------------------------------------------------------------------------------*/
/*    scal                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void BLAS_BLAS1::scal(
    const IndexType n,
    const ValueType alpha,
    ValueType* x,
    const IndexType incX )
{
    if ( incX <= 0 )
    {
        return;
    }

    SCAI_REGION( "BLAS.BLAS1.scal" )
    SCAI_LOG_DEBUG( logger,
                    "scal<" << TypeTraits<ValueType>::id() << " n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    BLASWrapper<ValueType>::scal( static_cast<BLASTrait::BLASIndexType>( n ), alpha, x,
                                  static_cast<BLASTrait::BLASIndexType>( incX ) );
}

/* ---------------------------------------------------------------------------------------*/
/*    nrm2                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType BLAS_BLAS1::nrm2( const IndexType n, const ValueType* x, const IndexType incX )
{
    if ( incX <= 0 )
    {
        return static_cast<ValueType>( 0.0 );
    }

    SCAI_REGION( "BLAS.BLAS1.nrm2" )
    SCAI_LOG_DEBUG( logger,
                    "nrm2<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", x = " << x << ", incX = " << incX )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    return BLASWrapper<ValueType>::nrm2( static_cast<BLASTrait::BLASIndexType>( n ), x, static_cast<BLASTrait::BLASIndexType>( incX ) );
}

/* ---------------------------------------------------------------------------------------*/
/*    asum                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType BLAS_BLAS1::asum( const IndexType n, const ValueType* x, const IndexType incX )
{
    SCAI_REGION( "BLAS.BLAS1.asum" )
    SCAI_LOG_DEBUG( logger,
                    "asum<" << TypeTraits<ValueType>::id() << ">,  n = " << n << ", x = " << x << ", incX = " << incX )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution of BLAS1:asum not supported, do it synchronously." )
    }

    if ( incX <= 0 )
    {
        return static_cast<ValueType>( 0.0 );
    }

    return BLASWrapper<ValueType>::asum( static_cast<BLASTrait::BLASIndexType>( n ), x, static_cast<BLASTrait::BLASIndexType>( incX ) );
}

/* ---------------------------------------------------------------------------------------*/
/*    iamax                                                                               */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
IndexType BLAS_BLAS1::iamax( const IndexType n, const ValueType* x, const IndexType incX )
{
    SCAI_REGION( "BLAS.BLAS1.iamax" )
    SCAI_LOG_INFO( logger,
                   "iamax<" << TypeTraits<ValueType>::id() << ">, " << "n = " << n << ", x = " << x << ", incX = " << incX )

    if ( ( n <= 0 ) || ( incX <= 0 ) )
    {
        return 0;
    }

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    return BLASWrapper<ValueType>::iamax( static_cast<BLASTrait::BLASIndexType>( n ), x, static_cast<BLASTrait::BLASIndexType>( incX ) ) - 1;
}

/* ---------------------------------------------------------------------------------------*/
/*    swap                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void BLAS_BLAS1::swap(
    const IndexType n,
    ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY )
{
    SCAI_REGION( "BLAS.BLAS1.swap" )
    SCAI_LOG_DEBUG( logger,
                    "iamax<" << TypeTraits<ValueType>::id() << ">" << ", n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    BLASWrapper<ValueType>::swap( static_cast<BLASTrait::BLASIndexType>( n ), x, static_cast<BLASTrait::BLASIndexType>( incX ), y, static_cast<BLASTrait::BLASIndexType>( incY ) );
}

/* ---------------------------------------------------------------------------------------*/
/*    copy                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void BLAS_BLAS1::copy(
    const IndexType n,
    const ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY )
{
    SCAI_REGION( "BLAS.BLAS1.copy" )
    SCAI_LOG_DEBUG( logger,
                    "copy<" << TypeTraits<ValueType>::id() << ">, " << "n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    // wrapper: uses overloading and converts arguments if required
    BLASWrapper<ValueType>::copy( static_cast<BLASTrait::BLASIndexType>( n ), x, static_cast<BLASTrait::BLASIndexType>( incX ), y, static_cast<BLASTrait::BLASIndexType>( incY ) );
}

/* ---------------------------------------------------------------------------------------*/
/*    axpy                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void BLAS_BLAS1::axpy(
    const IndexType n,
    const ValueType alpha,
    const ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY )
{
    SCAI_REGION( "BLAS.BLAS1.axpy" )
    SCAI_LOG_DEBUG( logger,
                    "axpy<" << TypeTraits<ValueType>::id() << ">, " << "n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        // just set call to this routine as task to be executed; new thread will have syncToken == NULL
        syncToken->run( std::bind( BLASWrapper<ValueType>::axpy, n, alpha, x, incX, y, incY ) );
    }
    else
    {
        BLASWrapper<ValueType>::axpy( static_cast<BLASTrait::BLASIndexType>( n ), alpha, x, static_cast<BLASTrait::BLASIndexType>( incX ), y, static_cast<BLASTrait::BLASIndexType>( incY ) );
    }
}

/* ---------------------------------------------------------------------------------------*/
/*    dot                                                                                 */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType BLAS_BLAS1::dot(
    const IndexType n,
    const ValueType* x,
    const IndexType incX,
    const ValueType* y,
    const IndexType incY )
{
    SCAI_REGION( "BLAS.BLAS1.dot" )
    SCAI_LOG_DEBUG( logger,
                    "dot<" << TypeTraits<ValueType>::id() << ">, " << "n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return static_cast<ValueType>( 0.0 );
    }

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    ValueType res = BLASWrapper<ValueType>::dot( static_cast<BLASTrait::BLASIndexType>( n ), x,
                    static_cast<BLASTrait::BLASIndexType>( incX ), y,
                    static_cast<BLASTrait::BLASIndexType>( incY ) );
    return res;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void BLAS_BLAS1::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::Host;
    bool useBLAS = false;
    int level = 0;
    useBLAS = common::Settings::getEnvironment( level, "SCAI_USE_BLAS" );

    if ( !useBLAS || ( level <= 0 ) )
    {
        SCAI_LOG_INFO( logger, "BLAS1 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS not set or 0)" )
        return;
    }
    else if ( level > 1 )
    {
        // only level 2 or level 3 wrappers might be used
        SCAI_LOG_INFO( logger,
                       "BLAS1 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS = " << level << ")" )
        return;
    }

    SCAI_LOG_DEBUG( logger, "register[" << flag << "] BLAS1 wrapper routines for Host at kernel registry: " <<
                    "T = " << common::TypeTraits<ValueType>::id() )

    KernelRegistry::set<BLASKernelTrait::scal<ValueType> >( BLAS_BLAS1::scal, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::nrm2<ValueType> >( BLAS_BLAS1::nrm2, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::asum<ValueType> >( BLAS_BLAS1::asum, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::iamax<ValueType> >( BLAS_BLAS1::iamax, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::swap<ValueType> >( BLAS_BLAS1::swap, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::copy<ValueType> >( BLAS_BLAS1::copy, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::axpy<ValueType> >( BLAS_BLAS1::axpy, ctx, flag );

    if ( common::isComplex( common::TypeTraits<ValueType>::stype ) )
    {
        /*
         * don't register dot-product due to different fortran ABSs the handling of returning
         * ComplexTypes can't be determined
         *
         * see: http://scicomp.stackexchange.com/questions/5380/intel-mkl-difference-between-mkl-intel-lp64-and-mkl-gf-lp64
         */
    }
    else
    {
        KernelRegistry::set<BLASKernelTrait::dot<ValueType> >( BLAS_BLAS1::dot, ctx, flag );
    }
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

BLAS_BLAS1::BLAS_BLAS1()
{
    SCAI_LOG_INFO( logger, "register BLAS1 wrapper routines for Host at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_EXT_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_REPLACE );
}

BLAS_BLAS1::~BLAS_BLAS1()
{
    SCAI_LOG_INFO( logger, "unregister BLAS1 wrapper routines for Host at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_EXT_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

BLAS_BLAS1 BLAS_BLAS1::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
