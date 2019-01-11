/**
 * @file OpenMPBLAS1.cpp
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
 * @brief Own implementation of BLAS1 routines for Host using OpenMP parallelization.
 * @author Eric Schricker
 * @date 02.10.2013
 */

// hpp
#include <scai/blaskernel/openmp/OpenMPBLAS1.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

// internal scai libraries

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/OpenMP.hpp>

#include <scai/tracing.hpp>

// std
#include <cmath>
#include <complex>
#include <functional>

namespace scai
{

using tasking::TaskSyncToken;
using common::TypeTraits;

namespace blaskernel
{

SCAI_LOG_DEF_LOGGER( OpenMPBLAS1::logger, "OpenMP.BLAS1" )

/** scal */

template<typename ValueType>
void OpenMPBLAS1::scal(
    const IndexType n,
    const ValueType alpha,
    ValueType* x,
    const IndexType incX )
{
    SCAI_REGION( "OpenMP.BLAS1.scal" )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        // run this routine with same arguments by another thread, running thread will have syncToken == NULL
        syncToken->run( std::bind( scal<ValueType>, n, alpha, x, incX ) );
        return;
    }

    if ( incX <= 0 )
    {
        return;
    }

    SCAI_LOG_DEBUG( logger,
                    "scal<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    if ( incX == 1 )
    {
        #pragma omp parallel for

        for ( IndexType i = 0; i < n; i++ )
        {
            x[i] = x[i] * alpha;
        }
    }
    else
    {
        //incX != 1
        #pragma omp parallel for 
        for ( IndexType i = 0; i < n; i++ )
        {
            x[i * incX] = x[i * incX] * alpha;
        }
    }
}

/** nrm2 */

template<typename ValueType>
ValueType OpenMPBLAS1::nrm2( const IndexType n, const ValueType* x, const IndexType incX )
{
    SCAI_REGION( "OpenMP.BLAS1.nrm2" )
    SCAI_LOG_INFO( logger,
                   "nrm2<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", x = " << x << ", incX = " << incX )
    ValueType sumOfSquares = 0;  // same as ... = ValueType( 0 )

    if ( incX <= 0 )
    {
        // considered as vector x is empty
        return sumOfSquares;
    }

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution not supported yet." )
    }

    // OpenMP reduction clause cannot be used as it doesn't support complex numbers
    #pragma omp parallel shared( sumOfSquares )
    {
        ValueType tSumOfSquares = 0; // private for each thread

        if ( incX == 1 )
        {
            #pragma omp for 

            for ( IndexType i = 0; i < n; i++ )
            {
                tSumOfSquares += x[i] * common::Math::conj( x[i] );
            }
        }
        else
        {
            #pragma omp for 

            for ( IndexType i = 0; i < n; i++ )
            {
                tSumOfSquares += x[i * incX] * common::Math::conj( x[i * incX] );
            }
        }

        atomicAdd( sumOfSquares, tSumOfSquares );
    }
    return common::Math::sqrt( sumOfSquares );
}

template<>
IndexType OpenMPBLAS1::nrm2( const IndexType, const IndexType*, const IndexType )
{
    COMMON_THROWEXCEPTION( "BLAS1.nrm2 unsupported for " << TypeTraits<IndexType>::stype );
}

/** asum (l1 norm) */

template<typename ValueType>
ValueType OpenMPBLAS1::asum( const IndexType n, const ValueType* x, const IndexType incX )
{
    SCAI_REGION( "OpenMP.BLAS1.asum" )
    SCAI_LOG_INFO( logger,
                   "asum<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", x = " << x << ", incX = " << incX )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution not supported yet" )
    }

    ValueType result = 0;

    if ( incX <= 0 )
    {
        return result;
    }

    // OpenMP reduction clause cannot be used as it doesn't support complex numbers

    #pragma omp parallel shared( result )
    {
        ValueType tResult = 0; // private for each thread

        if ( incX == 1 )
        {
            #pragma omp for 

            for ( IndexType i = 0; i < n; i++ )
            {
                tResult += common::Math::abs( common::Math::real( x[i] ) ) + common::Math::abs( common::Math::imag( x[i] ) );
            }
        }
        else
        {
            #pragma omp for 

            for ( IndexType i = 0; i < n; i++ )
            {
                tResult += common::Math::abs( common::Math::real( x[i * incX] ) ) + common::Math::abs( common::Math::imag( x[i * incX] ) );
            }
        }

        atomicAdd( result, tResult );
    }
    return result;
}

template<>
IndexType OpenMPBLAS1::asum( const IndexType, const IndexType*, const IndexType )
{
    COMMON_THROWEXCEPTION( "BLAS1.asum unsupported for IndexType = " << TypeTraits<IndexType>::stype );
}

/** iamax */

template<typename ValueType>
IndexType OpenMPBLAS1::iamax( const IndexType n, const ValueType* x, const IndexType incX )
{
    SCAI_REGION( "OpenMP.BLAS1.iamax" )
    SCAI_LOG_INFO( logger,
                   "iamax<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", x = " << x << ", incX = " << incX )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution not supported here" )
    }

    if ( n <= 0 || incX <= 0 )
    {
        return 0;
    }

    IndexType maxPos = 0;
    #pragma omp parallel shared( maxPos )
    {
        IndexType tMaxPos = 0;

        if ( incX == 1 )
        {
            #pragma omp for 

            for ( IndexType i = 0; i < n; i++ )
            {
                if ( common::Math::abs( x[i] ) > common::Math::abs( x[tMaxPos] ) )
                {
                    tMaxPos = i;
                }
            }
        }
        else
        {
            #pragma omp for

            for ( IndexType i = 0; i < n; i++ )
            {
                if ( common::Math::abs( x[i * incX] ) > common::Math::abs( x[tMaxPos * incX] ) )
                {
                    tMaxPos = i;
                }
            }
        }

        #pragma omp critical
        {
            if ( ( common::Math::abs( x[tMaxPos] ) > common::Math::abs( x[maxPos] ) )
                    || ( ( common::Math::abs( x[tMaxPos] ) == common::Math::abs( x[maxPos] ) ) && tMaxPos < maxPos ) )
            {
                maxPos = tMaxPos;
            }
        }
    }
    return maxPos;
}

/** swap */

template<typename ValueType>
void OpenMPBLAS1::swap(
    const IndexType n,
    ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY )
{
    SCAI_REGION( "OpenMP.BLAS1.swap" )
    SCAI_LOG_INFO( logger,
                   "swap<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution not supported here" )
    }

    ValueType* temp = 0;
    temp = new ValueType[n];

    if ( incX == 1 && incY == 1 )
    {
        #pragma omp parallel
        {
            #pragma omp for 

            for ( IndexType i = 0; i < n; i++ )
            {
                temp[i] = x[i];
            }

            #pragma omp for 

            for ( IndexType i = 0; i < n; i++ )
            {
                x[i] = y[i];
            }

            #pragma omp for 

            for ( IndexType i = 0; i < n; i++ )
            {
                y[i] = temp[i];
            }
        }
    }
    else
    {
        //incX != 1 || incY != 1
        #pragma omp parallel
        {
            #pragma omp for 

            for ( IndexType i = 0; i < n; i++ )
            {
                temp[i] = x[i * incX];
            }

            #pragma omp for 

            for ( IndexType i = 0; i < n; i++ )
            {
                x[i * incX] = y[i * incY];
            }

            #pragma omp for 

            for ( IndexType i = 0; i < n; i++ )
            {
                y[i * incY] = temp[i];
            }
        }
    }

    delete[] temp;
}

/** copy */

template<typename ValueType>
void OpenMPBLAS1::copy(
    const IndexType n,
    const ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY )
{
    SCAI_REGION( "OpenMP.BLAS1.copy" )
    SCAI_LOG_INFO( logger,
                   "copy<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    if ( incX == 1 && incY == 1 )
    {
        #pragma omp parallel for 

        for ( IndexType i = 0; i < n; i++ )
        {
            y[i] = x[i];
        }
    }
    else
    {
        //incX != 1 || incY != 1
        #pragma omp parallel for 
        for ( IndexType i = 0; i < n; i++ )
        {
            y[i * incY] = x[i * incX];
        }
    }
}

/** axpy */

template<typename ValueType>
void OpenMPBLAS1::axpy(
    const IndexType n,
    const ValueType alpha,
    const ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY )
{
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution starts now" )
        syncToken->run( std::bind( axpy<ValueType>,  n, alpha, x, incX, y, incY ) );
        return;
    }

    SCAI_REGION( "OpenMP.BLAS1.axpy" )
    SCAI_LOG_INFO( logger,
                   "axpy<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", alpha = " << alpha
                   << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    if ( incX == 1 && incY == 1 )
    {
        #pragma omp parallel for 

        for ( IndexType i = 0; i < n; i++ )
        {
            y[i] = alpha * x[i] + y[i];
        }
    }
    else
    {
        //incX != 1 || incY != 1
        #pragma omp parallel for 
        for ( IndexType i = 0; i < n; i++ )
        {
            y[i * incY] = alpha * x[i * incX] + y[i * incY];
        }
    }
}
/** dot */

template<typename ValueType>
ValueType OpenMPBLAS1::dot(
    const IndexType n,
    const ValueType* x,
    const IndexType incX,
    const ValueType* y,
    const IndexType incY )
{
    SCAI_REGION( "OpenMP.BLAS1.dot" )

    SCAI_LOG_INFO( logger,
                   "dot<" << TypeTraits<ValueType>::id() << ">, n = " << n
                   << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return ValueType( 0 );
    }

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution not supported here" )
    }

    ValueType result = 0;
    #pragma omp parallel shared( result )
    {
        ValueType tResult = 0;

        if ( incX == 1 && incY == 1 )
        {
            #pragma omp for 

            for ( IndexType i = 0; i < n; i++ )
            {
                tResult += common::Math::conj( x[i] ) * y[i] ;
            }
        }
        else
        {
            // incX != 1 || incY != 1
            #pragma omp for 
            for ( IndexType i = 0; i < n; i++ )
            {
                tResult += common::Math::conj( x[i * incX] ) * y[i * incY] ;
            }
        }

        atomicAdd( result, tResult );
    }
    return result;
}

template<>
IndexType OpenMPBLAS1::dot(
    const IndexType,
    const IndexType*,
    const IndexType,
    const IndexType*,
    const IndexType )
{
    COMMON_THROWEXCEPTION( "BLAS1.dot unsupported for IndexType = " << TypeTraits<IndexType>::stype );
}

/** sum */

template<typename ValueType>
void OpenMPBLAS1::sum(
    const IndexType n,
    ValueType alpha,
    const ValueType* x,
    ValueType beta,
    const ValueType* y,
    ValueType* z )
{
    SCAI_REGION( "OpenMP.BLAS1.sum" )

    SCAI_LOG_INFO( logger,
                   "sum<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y << ", z = " << z )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution not supported here" )
    }

    #pragma omp parallel for 

    for ( IndexType i = 0; i < n; i++ )
    {
        z[i] = alpha * x[i] + beta * y[i];
    }

    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPBLAS1::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::Host;

    SCAI_LOG_DEBUG( logger, "register[" << flag << "] BLAS1 OpenMP-routines for Host at kernel registry: " <<
                    "T = " << common::TypeTraits<ValueType>::id() )

    KernelRegistry::set<BLASKernelTrait::scal<ValueType> >( OpenMPBLAS1::scal, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::nrm2<ValueType> >( OpenMPBLAS1::nrm2, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::asum<ValueType> >( OpenMPBLAS1::asum, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::iamax<ValueType> >( OpenMPBLAS1::iamax, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::swap<ValueType> >( OpenMPBLAS1::swap, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::copy<ValueType> >( OpenMPBLAS1::copy, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::axpy<ValueType> >( OpenMPBLAS1::axpy, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::dot<ValueType> >( OpenMPBLAS1::dot, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::sum<ValueType> >( OpenMPBLAS1::sum, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPBLAS1::OpenMPBLAS1()
{
    SCAI_LOG_INFO( logger, "register BLAS1 OpenMP-routines for Host at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ADD );
}

OpenMPBLAS1::~OpenMPBLAS1()
{
    SCAI_LOG_INFO( logger, "unregister BLAS1 OpenMP-routines for Host at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPBLAS1 OpenMPBLAS1::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
