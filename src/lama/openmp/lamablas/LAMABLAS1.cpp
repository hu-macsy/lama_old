/**
 * @file LAMABLAS1.cpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief LAMABLAS.cpp
 * @author lschubert
 * @date 09.08.2012
 * $Id$
 */

#include <lama/openmp/lamablas/LAMABLAS1.hpp>

namespace lama
{

#include <lama/macros/unused.hpp>
#include <lama/openmp/BLASHelper.hpp>

#include <math.h>

namespace lama
{

LAMA_LOG_DEF_LOGGER( LAMABLAS1::logger, "LAMABLAS1" );

/** scal */

template<typename T>
void LAMABLAS1::scal( const IndexType n, const T alpha, T* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "scal<T>, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    #pragma omp parallel for
    for( int i = 0; i < n; i++ )
    {
        x[i * incX] = alpha * x[i * incX];
    }
}

// instantiation
template void LAMABLAS1::scal<float>(
    const IndexType n,
    const float alpha,
    const float* x,
    const IndexType incX,
    SyncToken* syncToken );
template void LAMABLAS1::scal<double>(
    const IndexType n,
    const double alpha,
    const double* x,
    const IndexType incX,
    SyncToken* syncToken );

/** nrm2 */

template<typename T>
T LAMABLAS1::nrm2( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "nrm2<T>, n = " << n << ", x = " << x << ", incX = " << incX );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    T tmp = 0;
    #pragma omp parallel for reduction(+ : tmp)
    for( int i = 0; i < n; i++ )
    {
        tmp = std::abs( x[i * incX] ) * std::abs( x[i * incX] );
    }
    return sqrt( tmp );
}

// instantiation
template void LAMABLAS1::nrm2<float>( const IndexType n, const float* x, const IndexType incX, SyncToken* syncToken );
template void LAMABLAS1::nrm2<double>( const IndexType n, const double* x, const IndexType incX, SyncToken* syncToken );

/** asum (l1 norm) */

template<typename T>
T LAMABLAS1::asum( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "asum<T>, n = " << n << ", x = " << x << ", incX = " << incX );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    T tmp = 0;
    #pragma omp parallel for reduction(+ : tmp)
    for( int i = 0; i < n; i++ )
    {
        tmp = std::abs( x[i * incX] );
    }
    return sqrt( tmp );
}

// instantiation
template void LAMABLAS1::asum<float>( const IndexType n, const float* x, const IndexType incX, SyncToken* syncToken );
template void LAMABLAS1::asum<double>( const IndexType n, const double* x, const IndexType incX, SyncToken* syncToken );

/** iamax */

template<typename T>
IndexType LAMABLAS1::iamax( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "iamax<T>, n = " << n << ", x = " << x << ", incX = " << incX );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    int iamax = 0;
    #pragma omp parallel for
    for( int i = 0; i < n; i++ )
    {
        #pragma omp critical
        if( x[i * incX] > iamax )
        {
            iamax = x[i * incX];
        }
    }
    return iamax;
}

// instantiation
template void LAMABLAS1::iamax<float>( const IndexType n, const float* x, const IndexType incX, SyncToken* syncToken );
template void LAMABLAS1::iamax<double>(
    const IndexType n,
    const double* x,
    const IndexType incX,
    SyncToken* syncToken );

/** viamax */

template<typename T>
T LAMABLAS1::viamax( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "viamax<T>, n = " << n << ", x = " << x << ", incX = " << incX );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    int iamax = 0;
    #pragma omp parallel for reduction( max : iamax)
    for( int i = 0; i < n; i++ )
    {
        if( x[i * incX] > iamax )
        {
            iamax = x[i * incX];
        }
    }
    return iamax;
}

// instantiation
template void LAMABLAS1::viamax<float>( const IndexType n, const float* x, const IndexType incX, SyncToken* syncToken );
template void LAMABLAS1::viamax<double>(
    const IndexType n,
    const double* x,
    const IndexType incX,
    SyncToken* syncToken );

/** swap */

template<typename T>
void LAMABLAS1::swap( const IndexType n, T* x, const IndexType incX, T* y, const IndexType incY, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "iamax<T>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    #pragma omp parallel for
    for( int i = 0; i < n; i++ )
    {
        T tmp = x[i * incX];
        x[i * incX] = y[i * incY];
        y[i * incY] = tmp;
    }
}

// instantiation
template void LAMABLAS1::swap<float>(
    const IndexType n,
    float* x,
    const IndexType incX,
    float* y,
    const IndexType incY,
    SyncToken* syncToken )
template void LAMABLAS1::swap<double>(
    const IndexType n,
    double* x,
    const IndexType incX,
    double* y,
    const IndexType incY,
    SyncToken* syncToken )
/** copy */

template<typename T>
void LAMABLAS1::copy(
    const IndexType n,
    const T* x,
    const IndexType incX,
    T* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "copy<T>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    #pragma omp parallel for
    for( int i = 0; i < n; i++ )
    {
        x[i * incX] = y[i * inc];
    }
}

// instantiation
template void LAMABLAS1::copy<float>(
    const IndexType n,
    const float* x,
    const IndexType incX,
    float* y,
    const IndexType incY,
    SyncToken* syncToken );
template void LAMABLAS1::copy<double>(
    const IndexType n,
    const double* x,
    const IndexType incX,
    double* y,
    const IndexType incY,
    SyncToken* syncToken );

/** axpy */

template<typename T>
void LAMABLAS1::axpy(
    const IndexType n,
    const T alpha,
    const T* x,
    const IndexType incX,
    T* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "axpy<T>, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    #pragma omp parallel for
    for( int i = 0; i < n; i++ )
    {
        x[i] = alpha * x[i * incX] + y[i * incY];
    }
}

// instantiation
template void LAMABLAS1::axpy<float>(
    const IndexType n,
    const float alpha,
    const float* x,
    const IndexType incX,
    float* y,
    const IndexType incY,
    SyncToken* syncToken );
template void LAMABLAS1::axpy<double>(
    const IndexType n,
    const double alpha,
    const double* x,
    const IndexType incX,
    double* y,
    const IndexType incY,
    SyncToken* syncToken );

/** dot */

template<typename T>
T LAMABLAS1::dot(
    const IndexType n,
    const T* x,
    const IndexType incX,
    const T* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "dot<T>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    T dot;
    #pragma omp parallel for reduction ( +: dot )
    for( int i = 0; i < n; i++ )
    {
        dot += x[i * incX] * y[i * incY];
    }
    return dot;
}

// instantiation
template void LAMABLAS1::dot<float>(
    const IndexType n,
    const float* x,
    const IndexType incX,
    const float* y,
    const IndexType incY,
    SyncToken* syncToken );
template void LAMABLAS1::dot<double>(
    const IndexType n,
    const double* x,
    const IndexType incX,
    const double* y,
    const IndexType incY,
    SyncToken* syncToken );

/** sum */

template<typename T>
void LAMABLAS1::sum( const IndexType n, T alpha, const T* x, T beta, const T* y, T* z, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "sum<T>, n = " << n << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y << ", z = " << z );

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

// instantiation
template void LAMABLAS1::sum<float>(
    const IndexType n,
    float alpha,
    const float* x,
    float beta,
    const float* y,
    float* z,
    SyncToken* syncToken );
template void LAMABLAS1::sum<double>(
    const IndexType n,
    double alpha,
    const double* x,
    double beta,
    const double* y,
    double* z,
    SyncToken* syncToken );

/** rot */

template<typename T>
void LAMABLAS1::rot(
    const IndexType n,
    T* x,
    const IndexType incX,
    T* y,
    const IndexType incY,
    const T c,
    const T s,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "rot<T>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY << ", c = " << c << ", s = " << s );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    T stemp;

    if( n < 0 )
    {
        return
    }

        if( incX == 1 && incY == 1 )
    {
        for( int i = 0; i < n; i++ )
        {
            stemp = c * x[i] + s * y[i];
            y[i] = c * y[i] - s * x[i];
            x[i] = stemp;
        }
    }
    else
    {
        int ix = 1;
        int iy = 1;

        //TODO: check if we need +1 (assumed from fortran)
        if( incX < 0 )
        {
            ix = ( -n + 1 ) * incx + 1;
        }
        if( incY < 0 )
        {
            iy = ( -n + 1 ) * incy + 1;
        }

        for( int i = 0; i < n; i++ )
        {
            stemp = c * x[ix] + s * y[iy];
            y( iy ) = c * y[iy] - s * x[ix];
            x( ix ) = stemp;
            ix = ix + incX;
            iy = iy + incY;
        }
    }

    return;
}

// instantiation
template void LAMABLAS1::rot<float>(
    const IndexType n,
    float* x,
    const IndexType incX,
    float* y,
    const IndexType incY,
    const float c,
    const float s,
    SyncToken* syncToken );
template void LAMABLAS1::rot<double>(
    const IndexType n,
    double* x,
    const IndexType incX,
    double* y,
    const IndexType incY,
    const double c,
    const double s,
    SyncToken* syncToken );

/** rotm */

template<typename T>
void LAMABLAS1::rotm(
    const IndexType n,
    T* x,
    const IndexType incX,
    T* y,
    const IndexType incY,
    const T* P,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "rotm<T>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY << ", P = " << P );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    return;
}

// instantiation
template void LAMABLAS1::rotm<float>(
    const IndexType n,
    float* x,
    const IndexType incX,
    float* y,
    const IndexType incY,
    const float* P,
    SyncToken* syncToken )
template void LAMABLAS1::rotm<double>(
    const IndexType n,
    double* x,
    const IndexType incX,
    double* y,
    const IndexType incY,
    const double* P,
    SyncToken* syncToken )

/** rotg */

template<typename T>
void LAMABLAS1::rotg( T* a, T* b, T* c, T* s, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "rotg<T>, a = " << a << ", b = " << b << ", c = " << c << ", s = " << s );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

//TODO
    return;
}

// instantiation
template void LAMABLAS1::rotg<float>( float* a, float* b, float* c, float* s, SyncToken* syncToken )
template void LAMABLAS1::rotg<double>( double* a, double* b, double* c, double* s, SyncToken* syncToken )

/** rotmg */

template<typename T>
void LAMABLAS1::rotmg( T* d1, T* d2, T* b1, const T b2, T* p, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "rotmg<float>, d1 = " << d1 << ", d2 = " << d2 << ", b1 = " << b1 << ", b2 = " << b2 << ", p = " << p );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

//TODO:
    return;
}

// instantiation
template void LAMABLAS1::rotmg<float>( float* d1, float* d2, float* b1, const float b2, float* p, SyncToken* syncToken )
template void LAMABLAS1::rotmg<double>(
    double* d1,
    double* d2,
    double* b1,
    const double b2,
    double* p,
    SyncToken* syncToken )

/* ass */

template<typename T>
void LAMABLAS1::ass( const IndexType n, const T value, T* x, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "ass<T>, n = " << n << ", value = " << value << ", x = " << x );

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    #pragma omp parallel for
    for( int i = 0; i < n; i++ )
    {
        x[i] = value;
    }
    return;
}

// instantiation
template void LAMABLAS1::ass<float>( const IndexType n, const float value, float* x, SyncToken* syncToken );
template void LAMABLAS1::ass<double>( const IndexType n, const double value, double* x, SyncToken* syncToken );

} /* namespace lama */
