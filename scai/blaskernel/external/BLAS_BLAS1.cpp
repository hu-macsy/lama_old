/**
 * @file BLAS_BLAS1.cpp
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
 * @brief Implementation of LAMA BLAS1 routines using wrappers to BLAS library
 * @author Lauretta Schubert
 * @date 05.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/blaskernel/external/BLAS_BLAS1.hpp>

// local library
#include <scai/blaskernel/external/BLASHelper.hpp>
#include <scai/blaskernel/external/BLASWrapper.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/cblas.hpp>

//Intel MPI need mpi.h to be included before stdio.h so this header comes first

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Settings.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/bind.hpp>

#include <scai/kregistry/KernelRegistry.hpp>

// boost
#include <boost/preprocessor.hpp>

// std
#include <cmath>

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
    if( incX <= 0 )
    {
        return;
    }

    SCAI_REGION( "BLAS.BLAS1.scal" )

    SCAI_LOG_DEBUG( logger,
                    "scal<" << getScalarType<ValueType>() << " n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

	BLASWrapper::scal(static_cast<BLASWrapper::IndexType>(n), alpha, x,
			static_cast<BLASWrapper::IndexType>( incX ) );
}

/* ---------------------------------------------------------------------------------------*/
/*    nrm2                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType BLAS_BLAS1::nrm2( const IndexType n, const ValueType* x, const IndexType incX )
{
    if( incX <= 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    SCAI_REGION( "BLAS.BLAS1.nrm2" )

    SCAI_LOG_DEBUG( logger,
                    "nrm2<" << getScalarType<ValueType>() << ">, n = " << n << ", x = " << x << ", incX = " << incX )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    return BLASWrapper::nrm2( static_cast<BLASWrapper::IndexType>( n ), x, static_cast<BLASWrapper::IndexType>( incX ));
}

/* ---------------------------------------------------------------------------------------*/
/*    asum                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType BLAS_BLAS1::asum( const IndexType n, const ValueType* x, const IndexType incX )
{
    SCAI_REGION( "BLAS.BLAS1.asum" )

    SCAI_LOG_DEBUG( logger,
                    "asum<" << getScalarType<ValueType>() << ">,  n = " << n << ", x = " << x << ", incX = " << incX )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution of BLAS1:asum not supported, do it synchronously." )
    }

    if( incX <= 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    return BLASWrapper::asum( static_cast<BLASWrapper::IndexType>( n ), x, static_cast<BLASWrapper::IndexType>( incX ));
}

/* ---------------------------------------------------------------------------------------*/
/*    iamax                                                                               */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
IndexType BLAS_BLAS1::iamax( const IndexType n, const ValueType* x, const IndexType incX )
{
    SCAI_REGION( "BLAS.BLAS1.iamax" )

    SCAI_LOG_INFO( logger,
                   "iamax<" << getScalarType<ValueType>() << ">, " << "n = " << n << ", x = " << x << ", incX = " << incX )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    return BLASWrapper::iamax( static_cast<BLASWrapper::IndexType>( n ), x, static_cast<BLASWrapper::IndexType>( incX ));
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
                    "iamax<" << getScalarType<ValueType>() << ">" << ", n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    BLASWrapper::swap( static_cast<BLASWrapper::IndexType>( n ), x, static_cast<BLASWrapper::IndexType>( incX ), y, static_cast<BLASWrapper::IndexType>( incY ));
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
                    "copy<" << getScalarType<ValueType>() << ">, " << "n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    // wrapper: uses overloading and converts arguments if required

    BLASWrapper::copy( static_cast<BLASWrapper::IndexType>( n ), x, static_cast<BLASWrapper::IndexType>( incX ), y, static_cast<BLASWrapper::IndexType>( incY ));
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
                    "axpy<" << getScalarType<ValueType>() << ">, " << "n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        // just set call to this routine as task to be executed; new thread will have syncToken == NULL
    	BLASWrapper::IndexType n_i = static_cast<BLASWrapper::IndexType>( n );
    	BLASWrapper::IndexType incX_i = static_cast<BLASWrapper::IndexType>( incX );
    	BLASWrapper::IndexType incY_i = static_cast<BLASWrapper::IndexType>( incY );

        syncToken->run( common::bind( BLASWrapper::axpy<ValueType>, n_i, alpha, x, incX_i, y, incY_i ) );
    }
    else
    {
        BLASWrapper::axpy( static_cast<BLASWrapper::IndexType>( n ), alpha, x, static_cast<BLASWrapper::IndexType>( incX ), y, static_cast<BLASWrapper::IndexType>( incY ));
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
                    "dot<" << getScalarType<ValueType>() << ">, " << "n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return static_cast<ValueType>(0.0);
    }

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    return BLASWrapper::dot( static_cast<BLASWrapper::IndexType>( n ), x, static_cast<BLASWrapper::IndexType>( incX ), y, static_cast<BLASWrapper::IndexType>( incY ) );
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void BLAS_BLAS1::registerKernels( bool deleteFlag )
{
    bool useBLAS = false;
    int level = 0;

    useBLAS = common::Settings::getEnvironment( level, "SCAI_USE_BLAS" );

    if( !useBLAS || ( level <= 0 ) )
    {
        SCAI_LOG_INFO( logger, "BLAS1 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS not set or 0)" )
        return;
    }
    else if( level > 1 )
    {
        // only level 2 or level 3 wrappers might be used
        SCAI_LOG_INFO( logger,
                       "BLAS1 wrapper routines for Host Interface are disabled (SCAI_USE_BLAS = " << level << ")" )
        return;
    }

    SCAI_LOG_INFO( logger, "register BLAS1 wrapper routines for Host at kernel registry" )

    // REGISTER1: give these routines priority in case of overriding

    using kregistry::KernelRegistry;
    using common::context::Host;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_REPLACE;   // priority over OpenMPBLAS

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

#define LAMA_BLAS1_REGISTER(z, I, _)                                                             \
    KernelRegistry::set<BLASKernelTrait::scal<ARITHMETIC_HOST_TYPE_##I> >( scal, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::nrm2<ARITHMETIC_HOST_TYPE_##I> >( nrm2, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::asum<ARITHMETIC_HOST_TYPE_##I> >( asum, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::iamax<ARITHMETIC_HOST_TYPE_##I> >( iamax, Host, flag );  \
    KernelRegistry::set<BLASKernelTrait::swap<ARITHMETIC_HOST_TYPE_##I> >( swap, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::copy<ARITHMETIC_HOST_TYPE_##I> >( copy, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::axpy<ARITHMETIC_HOST_TYPE_##I> >( axpy, Host, flag );    \
    KernelRegistry::set<BLASKernelTrait::dot<ARITHMETIC_HOST_TYPE_##I> >( dot, Host, flag );      \

    BOOST_PP_REPEAT( ARITHMETIC_HOST_EXT_TYPE_CNT, LAMA_BLAS1_REGISTER, _ )

#undef LAMA_BLAS1_REGISTER
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

BLAS_BLAS1::BLAS_BLAS1()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

BLAS_BLAS1::~BLAS_BLAS1()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

BLAS_BLAS1 BLAS_BLAS1::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
