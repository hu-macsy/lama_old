/**
 * @file OpenMPBLAS1.cpp
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
 * @brief Implementation of used BLAS1 routines in C++ wtih OpenMP parallelization.
 * @author Eric Schricker
 * @date 02.10.2013
 * @since 1.1.0
 */

// hpp
#include <scai/lama/openmp/OpenMPBLAS1.hpp>

// local library
#include <scai/lama/openmp/OpenMP.hpp>
#include <scai/lama/BLASInterface.hpp>
#include <scai/lama/LAMAInterfaceRegistry.hpp>

// internal scai libraries
#include <scai/common/macros/unused.hpp>

#include <scai/tracing.hpp>

// boost
#include <boost/preprocessor.hpp>

// std
#include <cmath>

namespace scai
{

namespace lama
{

/** Template routines use abs and sqrt, once from cmath and also as routines
 *  for complex numbers. To make this possible, all routines should be in
 *  this namespace made availale.
 */

using std::abs;
// used for float, double
using std::sqrt;
// used for float, double

using tasking::SyncToken;
using common::getScalarType;

SCAI_LOG_DEF_LOGGER( OpenMPBLAS1::logger, "OpenMP.BLAS1" )

/** scal */

template<typename ValueType>
void OpenMPBLAS1::scal(
    const IndexType n,
    const ValueType alpha,
    ValueType* x,
    const IndexType incX,
    SyncToken* syncToken )
{
    if( incX <= 0 )
    {
        return;
    }

    SCAI_REGION( "OpenMP.BLAS1.scal" )

    SCAI_LOG_DEBUG( logger,
                    "scal<" << getScalarType<ValueType>()<< ">, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    if( incX == 1 )
    {
        #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

        for( int i = 0; i < n; i++ )
        {
            x[i] = x[i] * alpha;
        }
    }
    else
    {
        //incX != 1
        #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )
        for( int i = 0; i < n; i++ )
        {
            x[i * incX] = x[i * incX] * alpha;
        }
    }
}

/** nrm2 */

template<typename ValueType>
ValueType OpenMPBLAS1::nrm2( const IndexType n, const ValueType* x, const IndexType incX, SyncToken* syncToken )
{
    SCAI_REGION( "OpenMP.BLAS1.nrm2" )

    SCAI_LOG_DEBUG( logger,
                    "nrm2<" << getScalarType<ValueType>()<< ">, n = " << n << ", x = " << x << ", incX = " << incX )

    if( incX <= 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    ValueType sumOfSquares = static_cast<ValueType>(0.0);

// OpenMP reduction clause cannot be used as it doesn't support complex numbers

    #pragma omp parallel shared( sumOfSquares )
    {
        ValueType tSumOfSquares = 0; // private for each thread

        if( incX == 1 )
        {
            #pragma omp for  schedule( SCAI_OMP_SCHEDULE )

            for( int i = 0; i < n; i++ )
            {
                tSumOfSquares += x[i] * x[i];
            }
        }
        else
        {
            #pragma omp for  schedule( SCAI_OMP_SCHEDULE )

            for( int i = 0; i < n; i++ )
            {
                tSumOfSquares += x[i * incX] * x[i * incX];
            }
        }

        atomicAdd( sumOfSquares, tSumOfSquares );
    }
    return sqrt( sumOfSquares );
}

/** asum (l1 norm) */

template<typename ValueType>
ValueType OpenMPBLAS1::asum( const IndexType n, const ValueType* x, const IndexType incX, SyncToken* syncToken )
{
    SCAI_REGION( "OpenMP.BLAS1.asum" )

    SCAI_LOG_DEBUG( logger,
                    "asum<" << getScalarType<ValueType>()<< ">, n = " << n << ", x = " << x << ", incX = " << incX )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    ValueType result = 0;

    if( incX <= 0 )
    {
        return result;
    }

// OpenMP reduction clause cannot be used as it doesn't support complex numbers

    #pragma omp parallel shared( result )
    {
        ValueType tResult = 0; // private for each thread

        if( incX == 1 )
        {
            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( int i = 0; i < n; i++ )
            {
                tResult += abs( x[i] );
            }
        }
        else
        {
            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( int i = 0; i < n; i++ )
            {
                tResult += abs( x[i * incX] );
            }
        }

        atomicAdd( result, tResult );
    }

    return result;
}

/** iamax */

template<typename ValueType>
IndexType OpenMPBLAS1::iamax( const IndexType n, const ValueType* x, const IndexType incX, SyncToken* syncToken )
{
    SCAI_REGION( "OpenMP.BLAS1.iamax" )

    SCAI_LOG_INFO( logger,
                   "iamax<" << getScalarType<ValueType>()<< ">, n = " << n << ", x = " << x << ", incX = " << incX )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    if( n <= 0 || incX <= 0 )
    {
        return 0;
    }

    IndexType maxPos = 0;

    #pragma omp parallel shared( maxPos )
    {
        IndexType tMaxPos = 0;

        if( incX == 1 )
        {
            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( int i = 0; i < n; i++ )
            {
                if( abs( x[i] ) > abs( x[tMaxPos] ) )
                {
                    tMaxPos = i;
                }
            }
        }
        else
        {
            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( int i = 0; i < n; i++ )
            {
                if( abs( x[i * incX] ) > abs( x[tMaxPos * incX] ) )
                {
                    tMaxPos = i;
                }
            }
        }

        #pragma omp critical
        {
            if( ( abs( x[tMaxPos] ) > abs( x[maxPos] ) )
                    || ( ( abs( x[tMaxPos] ) == abs( x[maxPos] ) ) && tMaxPos < maxPos ) )
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
    const IndexType incY,
    SyncToken* syncToken )
{
    SCAI_REGION( "OpenMP.BLAS1.swap" )

    SCAI_LOG_DEBUG( logger,
                    "iamax<" << getScalarType<ValueType>()<< ">, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    ValueType *temp = 0;
    temp = new ValueType[n];

    if( incX == 1 && incY == 1 )
    {
        #pragma omp parallel
        {
            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( int i = 0; i < n; i++ )
            {
                temp[i] = x[i];
            }

            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( int i = 0; i < n; i++ )
            {
                x[i] = y[i];
            }

            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( int i = 0; i < n; i++ )
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
            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( int i = 0; i < n; i++ )
            {
                temp[i] = x[i * incX];
            }

            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( int i = 0; i < n; i++ )
            {
                x[i * incX] = y[i * incY];
            }

            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( int i = 0; i < n; i++ )
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
    const IndexType incY,
    SyncToken* UNUSED(syncToken) )
{
    SCAI_REGION( "OpenMP.BLAS1.copy" )

    SCAI_LOG_DEBUG( logger,
                    "copy<" << getScalarType<ValueType>() << ">, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    if( incX == 1 && incY == 1 )
    {
        #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

        for( int i = 0; i < n; i++ )
        {
            y[i] = x[i];
        }
    }
    else
    {
        //incX != 1 || incY != 1
        #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )
        for( int i = 0; i < n; i++ )
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
    const IndexType incY,
    SyncToken* syncToken )
{
    SCAI_REGION( "OpenMP.BLAS1.axpy" )

    SCAI_LOG_DEBUG( logger,
                    "axpy<" << getScalarType<ValueType>() << ">, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    if( incX == 1 && incY == 1 )
    {
        #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

        for( int i = 0; i < n; i++ )
        {
            y[i] = alpha * x[i] + y[i];
        }
    }
    else
    {
        //incX != 1 || incY != 1
        #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )
        for( int i = 0; i < n; i++ )
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
    const IndexType incY,
    SyncToken* syncToken )
{
    SCAI_REGION( "OpenMP.BLAS1.sdot" )

    SCAI_LOG_DEBUG( logger,
                    "dot<" << getScalarType<ValueType>() << ">, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return static_cast<ValueType>(0.0);
    }

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    ValueType result = 0;

    #pragma omp parallel shared( result )
    {
        ValueType tResult = 0;

        if( incX == 1 && incY == 1 )
        {
            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( int i = 0; i < n; i++ )
            {
                tResult += x[i] * y[i];
            }
        }
        else
        {
            // incX != 1 || incY != 1
            #pragma omp for schedule( SCAI_OMP_SCHEDULE )
            for( int i = 0; i < n; i++ )
            {
                tResult += x[i * incX] * y[i * incY];
            }
        }

        atomicAdd( result, tResult );
    }
    return result;
}

/** sum */

template<typename ValueType>
void OpenMPBLAS1::sum(
    const IndexType n,
    ValueType alpha,
    const ValueType* x,
    ValueType beta,
    const ValueType* y,
    ValueType* z,
    SyncToken* syncToken )
{
    SCAI_REGION( "OpenMP.BLAS1.dot" )

    SCAI_LOG_DEBUG( logger,
                    "sum<" << getScalarType<ValueType>() << ">, n = " << n << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y << ", z = " << z )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

    for( int i = 0; i < n; i++ )
    {
        z[i] = alpha * x[i] + beta * y[i];
    }

    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPBLAS1::setInterface( BLASInterface& BLAS )
{
    SCAI_LOG_INFO( logger, "set BLAS1 routines for OpenMP in Interface" )

// Note: macro takes advantage of same name for routines and type definitions
//       ( e.g. routine CUDABLAS1::sum<ValueType> is set for BLAS::BLAS1::sum variable

#define LAMA_BLAS1_REGISTER(z, I, _)                                        \
    LAMA_INTERFACE_REGISTER_T( BLAS, scal, ARITHMETIC_HOST_TYPE_##I )       \
    LAMA_INTERFACE_REGISTER_T( BLAS, nrm2, ARITHMETIC_HOST_TYPE_##I )       \
    LAMA_INTERFACE_REGISTER_T( BLAS, asum, ARITHMETIC_HOST_TYPE_##I )       \
    LAMA_INTERFACE_REGISTER_T( BLAS, iamax, ARITHMETIC_HOST_TYPE_##I )      \
    LAMA_INTERFACE_REGISTER_T( BLAS, swap, ARITHMETIC_HOST_TYPE_##I )       \
    LAMA_INTERFACE_REGISTER_T( BLAS, copy, ARITHMETIC_HOST_TYPE_##I )       \
    LAMA_INTERFACE_REGISTER_T( BLAS, axpy, ARITHMETIC_HOST_TYPE_##I )       \
    LAMA_INTERFACE_REGISTER_T( BLAS, dot, ARITHMETIC_HOST_TYPE_##I )        \
    LAMA_INTERFACE_REGISTER_T( BLAS, sum, ARITHMETIC_HOST_TYPE_##I )        \

    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_BLAS1_REGISTER, _ )

#undef LAMA_BLAS1_REGISTER

}

/* --------------------------------------------------------------------------- */
/*    Static registration of the BLAS1 routines                                */
/* --------------------------------------------------------------------------- */

bool OpenMPBLAS1::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( hmemo::context::Host );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPBLAS1::initialized = registerInterface();

} /* end namespace lama */

} /* end namespace scai */