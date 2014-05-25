/**
 * @file blas/OpenMPBLAS1.cpp
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
 * @date 02.10.2013
 * @since 1.1.0
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

#include <boost/preprocessor.hpp>

//#include <cmath>

namespace lama
{

LAMA_LOG_DEF_LOGGER( OpenMPBLAS1::logger, "OpenMP.BLAS1" )

/** scal */

template<typename T>
void OpenMPBLAS1::scal( const IndexType n, const T alpha, T* x, const IndexType incX, SyncToken* syncToken )
{

    if ( incX <= 0 )
    {
        return;
    }

LAMA_REGION( "OpenMP.BLAS1.sscal" )

LAMA_LOG_DEBUG( logger, "scal<" << Scalar::getType<T>()<< ">, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    if ( incX == 1 )
    {
#pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
        for ( int i = 0; i < n; i++ )
        {
            x[i] = x[i] * alpha;
        }
    }
    else
    { //incX != 1
#pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
        for ( int i = 0; i < n; i++ )
        {
            x[i * incX] = x[i * incX] * alpha;
        }
    }
}

/** nrm2 */

template<typename T>
T OpenMPBLAS1::nrm2( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "nrm2<" << Scalar::getType<T>()<< ">, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( incX <= 0 )
    {
        return 0.0;
    }

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    T sumOfSquares = 0.0;
//no use of reduction, because it doesn't work with complex numbers
    if ( incX == 1 )
    {
#pragma omp parallel for reduction(+:sumOfSquares) schedule( LAMA_OMP_SCHEDULE )
        for ( int i = 0; i < n; i++ )
        {
            sumOfSquares += ( x[i] * x[i] );
        }
    }
    else
    { //incX != 1
#pragma omp parallel for reduction(+:sumOfSquares) schedule( LAMA_OMP_SCHEDULE )
        for ( int i = 0; i < n; i++ )
        {
            sumOfSquares += ( x[i * incX] * x[i * incX] );
        }
    }
    return ::sqrt( sumOfSquares );
}

/** asum (l1 norm) */

template<typename T>
T OpenMPBLAS1::asum( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "asum<" << Scalar::getType<T>()<< ">, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    T result = 0;

    if ( incX <= 0 )
    {
        return result;
    }

    if( incX == 1 )
    {
//no use of reduction, because it doesn't work with complex
#pragma omp parallel for reduction(+:result) schedule( LAMA_OMP_SCHEDULE )
        for ( int i = 0; i < n; i++ )
        {
            result = result + ::abs( x[i] );
        }
    }
    else
    { //incX != 1
#pragma omp parallel for reduction(+:result) schedule( LAMA_OMP_SCHEDULE )
        for ( int i = 0; i < n; i++ )
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
    LAMA_LOG_INFO( logger, "iamax<" << Scalar::getType<T>()<< ">, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    if ( n <= 0 || incX <= 0 )
    {
        return 0;
    }

    IndexType max_pos = 0;

    if ( incX == 1 )
    {
#pragma omp parallel shared(max_pos)
        {
            IndexType priv_max_pos = 0;
#pragma omp for schedule( LAMA_OMP_SCHEDULE )
            for ( int i = 0; i < n; i++ )
            {
                if( ::abs( x[i] ) > ::abs( x[priv_max_pos] ) )
                {
                    priv_max_pos = i;
                }
            }
#pragma omp critical
            {
                if ( ( ::abs( x[priv_max_pos] ) > ::abs( x[max_pos] ) )
                                || ( ( ::abs( x[priv_max_pos] ) == ::abs( x[max_pos] ) ) && priv_max_pos < max_pos ) )
                {
                    max_pos = priv_max_pos;
                }
            }
        }
    }
    else
    { //incX != 1
#pragma omp parallel shared(max_pos)
        {
            IndexType priv_max_pos = 0;
#pragma omp for schedule( LAMA_OMP_SCHEDULE )
            for ( int i = 0; i < n; i++ )
            {
                if ( ::abs( x[i * incX] ) > ::abs( x[priv_max_pos * incX] ) )
                {
                    priv_max_pos = i;
                }
            }
#pragma omp critical
            {
                if ( ( ::abs( x[priv_max_pos * incX] ) > ::abs( x[max_pos * incX] ) )
                                || ( ( ::abs( x[priv_max_pos * incX] ) == ::abs( x[max_pos * incX] ) )
                                                && priv_max_pos < max_pos ) )
                {
                    max_pos = priv_max_pos;
                }
            }
        }
    }

    return max_pos;
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
                    "iamax<" << Scalar::getType<T>()<< ">, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    T *temp = 0;
    temp = new T[n];

    if ( incX == 1 && incY == 1 )
    {
#pragma omp parallel
        {
#pragma omp for schedule( LAMA_OMP_SCHEDULE )
            for ( int i = 0; i < n; i++ )
            {
                temp[i] = x[i];
            }
#pragma omp for schedule( LAMA_OMP_SCHEDULE )
            for ( int i = 0; i < n; i++ )
            {
                x[i] = y[i];
            }
#pragma omp for schedule( LAMA_OMP_SCHEDULE )
            for ( int i = 0; i < n; i++ )
            {
                y[i] = temp[i];
            }
        }
    }
    else
    { //incX != 1 || incY != 1
#pragma omp parallel
        {
#pragma omp for schedule( LAMA_OMP_SCHEDULE )
            for ( int i = 0; i < n; i++ )
            {
                temp[i] = x[i * incX];
            }
#pragma omp for schedule( LAMA_OMP_SCHEDULE )
            for ( int i = 0; i < n; i++ )
            {
                x[i * incX] = y[i * incY];
            }
#pragma omp for schedule( LAMA_OMP_SCHEDULE )
            for ( int i = 0; i < n; i++ )
            {
                y[i * incY] = temp[i];
            }
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
                    "copy<" << Scalar::getType<T>() << ">, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    if ( incX == 1 && incY == 1 )
    {
#pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
        for ( int i = 0; i < n; i++ )
        {
            y[i] = x[i];
        }
    }
    else
    { //incX != 1 || incY != 1
#pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
        for ( int i = 0; i < n; i++ )
        {
            y[i * incY] = x[i * incX];
        }
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
                    "axpy<" << Scalar::getType<T>() << ">, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return;
    }

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    if ( incX == 1 && incY == 1 )
    {
#pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
        for ( int i = 0; i < n; i++ )
        {
            y[i] = alpha * x[i] + y[i];
        }
    }
    else
    { //incX != 1 || incY != 1
#pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
        for ( int i = 0; i < n; i++ )
        {
            y[i * incY] = alpha * x[i * incX] + y[i * incY];
        }
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
                    "dot<" << Scalar::getType<T>() << ">, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( ( incX <= 0 ) || ( incY <= 0 ) )
    {
        return 0.0;
    }

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    T result = 0;
//no use of reduction, because it doesn't work with complex
    if ( incX == 1 && incY == 1 )
    {
#pragma omp parallel for reduction(+: result) schedule( LAMA_OMP_SCHEDULE )
        for ( int i = 0; i < n; i++ )
        {
            result += ( x[i * incX] * y[i * incY] );
        }
    }
    else
    { //incX != 1 || incY != 1
#pragma omp parallel for reduction(+: result) schedule( LAMA_OMP_SCHEDULE )
        for ( int i = 0; i < n; i++ )
        {
            result += ( x[i * incX] * y[i * incY] );
        }
    }
    return result;
}

/** sum */

template<typename T>
void OpenMPBLAS1::sum( const IndexType n, T alpha, const T* x, T beta, const T* y, T* z, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "sum<" << Scalar::getType<T>() << ">, n = " << n << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y << ", z = " << z )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
    for ( int i = 0; i < n; i++ )
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

#define LAMA_BLAS1_REGISTER(z, I, _)                                            \
    LAMA_INTERFACE_REGISTER_T( BLAS, scal, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( BLAS, nrm2, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( BLAS, asum, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( BLAS, iamax, ARITHMETIC_TYPE##I )                \
    LAMA_INTERFACE_REGISTER_T( BLAS, swap, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( BLAS, copy, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( BLAS, axpy, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( BLAS, dot, ARITHMETIC_TYPE##I )                  \
    LAMA_INTERFACE_REGISTER_T( BLAS, sum, ARITHMETIC_TYPE##I )                  \

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_BLAS1_REGISTER, _ )

#undef LAMA_BLAS1_REGISTER

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
