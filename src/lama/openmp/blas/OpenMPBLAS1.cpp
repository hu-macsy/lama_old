/**
 * @file OpenMPBLAS1.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief OpenMPBLAS1.cpp
 * @author Eric Schricker
 * @date 05.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/openmp/OpenMPBLAS1.hpp>

// others
#include <lama/openmp/BLASHelper.hpp>
#include <lama/BLASInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

// macros
#include <lama/macros/unused.hpp>

// tracing with LAMA_REGION
#include <lama/tracing.hpp>

#include <cmath>

namespace lama
{

LAMA_LOG_DEF_LOGGER( OpenMPBLAS1::logger, "OpenMP.BLAS1" )

/** scal */

template<typename T>
void OpenMPBLAS1::scal( const IndexType n, const T alpha, T* x, const IndexType incX, SyncToken* syncToken )
{

    if( incX <= 0 )
    {
        return;
    }

LAMA_REGION( "OpenMP.BLAS1.sscal" )

                            LAMA_LOG_DEBUG( logger, "scal<float>, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#pragma omp parallel for
    for( int i = 0; i < n; i++ )
    {
        x[i * incX] = x[i * incX] * alpha;
    }
}

/** nrm2 */

template<typename T>
T OpenMPBLAS1::nrm2( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "nrm2<float>, n = " << n << ", x = " << x << ", incX = " << incX )

    if( incX <= 0 )
    {
        return 0.0;
    }

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    T sumOfSquares = 0.0;
#pragma omp parallel for reduction(+:sumOfSquares)
    for( int i = 0; i < n; i++ )
    {
        sumOfSquares += ( x[i * incX] * x[i * incX] );
    }

    T result = ::sqrt( sumOfSquares );
    return result;
}

/** asum (l1 norm) */

template<typename T>
T OpenMPBLAS1::asum( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "asum<float>, n = " << n << ", x = " << x << ", incX = " << incX )

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    T result = 0;

    if( incX > 0 )
    {
#pragma omp parallel for reduction(+:result)
        for( int i = 0; i < n; i++ )
        {
            result = result + ::abs( x[i * incX] );
        }
    }
    return result;

}

/** iamax */

template<typename T>
IndexType OpenMPBLAS1::iamax( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "iamax<float>, n = " << n << ", x = " << x << ", incX = " << incX )

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    if( n <= 0 || incX <= 0 )
    {
        return 0;
    }

    IndexType max_pos = 0;

    //OpenMP 3.1 added min/max reduction operators (GCC >= 4.7)
#pragma omp parallel shared(max_pos)
    for( int i = 0; i < n; i++ )
    {
#pragma omp critical (OpenMPBLAS1_iamax)
        if( ::abs( x[i * incX] ) > ::abs( x[max_pos * incX] ) )
        {
            max_pos = i;
        }
    }

    return max_pos;
}

/** viamax */

template<typename T>
T OpenMPBLAS1::viamax( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "viamax<float>, n = " << n << ", x = " << x << ", incX = " << incX )

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    IndexType max_val = iamax(n,x,incX, syncToken);

    return x[max_val * incX];
}

/** swap */

template<typename T>
void OpenMPBLAS1::swap(
    const IndexType n,
    T* x,
    const IndexType incX,
    T* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "iamax<float>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    T *temp = 0;
    temp = new T[n];

#pragma omp parallel
    {
#pragma omp for
        for( int i = 0; i < n; i++ )
        {
            temp[i] = x[i * incX];
        }
#pragma omp for
        for( int i = 0; i < n; i++ )
        {
            x[i * incX] = y[i * incY];
        }
#pragma omp for
        for( int i = 0; i < n; i++ )
        {
            y[i * incY] = temp[i];
        }
    }
    delete[] temp;
}

/** copy */

template<typename T>
void OpenMPBLAS1::copy(
    const IndexType n,
    const T* x,
    const IndexType incX,
    T* y,
    const IndexType incY,
    SyncToken* UNUSED(syncToken) )
{
    LAMA_LOG_DEBUG( logger,
                    "copy<float>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

#pragma omp parallel for
    for( int i = 0; i < n; i++ )
    {
        y[i * incY] = x[i * incX];
    }
}

/** axpy */

template<typename T>
void OpenMPBLAS1::axpy(
    const IndexType n,
    const T alpha,
    const T* x,
    const IndexType incX,
    T* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "axpy<float>, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#pragma omp parallel for
    for( int i = 0; i < n; i++ )
    {
        y[i * incY] = alpha * x[i * incX] + y[i * incY];
    }
}
/** dot */

template<typename T>
T OpenMPBLAS1::dot(
    const IndexType n,
    const T* x,
    const IndexType incX,
    const T* y,
    const IndexType incY,
    SyncToken* syncToken )
{
LAMA_REGION( "OpenMP.BLAS1.sdot" )

	                            LAMA_LOG_DEBUG( logger,
                    "dot<float>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return 0.0;
    }

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    T result = 0;

#pragma omp parallel for reduction(+: result)
    for( int i = 0; i < n; i++ )
    {
        result += ( x[i * incX] * y[i * incY] );
    }

    return result;
}

/** sum */

template<typename T>
void OpenMPBLAS1::sum( const IndexType n, T alpha, const T* x, T beta, const T* y, T* z, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "sum<T>, n = " << n << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y << ", z = " << z )

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#pragma omp parallel for
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
    LAMA_LOG_INFO( logger, "set BLAS1 routines for OpenMP in Interface" )

    // Note: macro takes advantage of same name for routines and type definitions 
    //       ( e.g. routine CUDABLAS1::sum<T> is set for BLAS::BLAS1::sum variable

    LAMA_INTERFACE_REGISTER_T( BLAS, scal, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, scal, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, nrm2, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, nrm2, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, asum, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, asum, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, iamax, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, iamax, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, viamax, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, viamax, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, swap, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, swap, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, copy, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, copy, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, axpy, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, axpy, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, dot, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, dot, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, sum, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, sum, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the BLAS1 routines                                */
/* --------------------------------------------------------------------------- */

bool OpenMPBLAS1::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::Host );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPBLAS1::initialized = registerInterface();

} /** namespace lama */
