/**
 * @file OpenMPBLAS1.cpp
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
 * @brief OpenMPBLAS1.cpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 * $Id$
 */

// hpp
#include <lama/openmp/OpenMPBLAS1.hpp>

// others
#include <lama/openmp/BLASHelper.hpp>
#include <lama/BLASInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

// macros
#include <lama/macros/unused.hpp>

#include <cmath>

namespace lama
{

LAMA_LOG_DEF_LOGGER( OpenMPBLAS1::logger, "OpenMP.BLAS1" )

/** scal */

template<>
void OpenMPBLAS1::scal( const IndexType n, const float alpha, float* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "scal<float>, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX;
#else
#define F77_N n
#define F77_incX incX
#endif
    F77_sscal( &F77_N, &alpha, x, &F77_incX );
}

template<>
void OpenMPBLAS1::scal( const IndexType n, const double alpha, double* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "scal<double>, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX;
#else
#define F77_N n
#define F77_incX incX
#endif
    F77_dscal( &F77_N, &alpha, x, &F77_incX );
}

/** nrm2 */

template<>
float OpenMPBLAS1::nrm2( const IndexType n, const float* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "nrm2<float>, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    float nrm2;
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX;
#else
#define F77_N n
#define F77_incX incX
#endif
    nrm2 = F77_snrm2( &F77_N, const_cast<float*>( x ), &F77_incX );
    return nrm2;
}

template<>
double OpenMPBLAS1::nrm2( const IndexType n, const double* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "nrm2<double>, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    double nrm2;
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX;
#else
#define F77_N n
#define F77_incX incX
#endif
    nrm2 = F77_dnrm2( &F77_N, const_cast<double*>( x ), &F77_incX );
    return nrm2;
}

/** asum (l1 norm) */

template<>
float OpenMPBLAS1::asum( const IndexType n, const float* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "asum<float>, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    float asum;
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#endif
    asum = F77_sasum( &F77_N, (float*) x, &F77_incX );
    return asum;
}

template<>
double OpenMPBLAS1::asum( const IndexType n, const double* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "asum<double>, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    double asum;
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#endif
    asum = F77_dasum( &F77_N, (double*) x, &F77_incX );
    return asum;
}

/** iamax */

template<>
IndexType OpenMPBLAS1::iamax( const IndexType n, const float* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "iamax<float>, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    int iamax;
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX;
#else
#define F77_N n
#define F77_incX incX
#endif
    iamax = F77_isamax( &F77_N, const_cast<float *>( x ), &F77_incX );
    return iamax ? iamax - 1 : 0;
}

template<>
IndexType OpenMPBLAS1::iamax( const IndexType n, const double* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "iamax<double>, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    int iamax;
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX;
#else
#define F77_N n
#define F77_incX incX
#endif
    iamax = F77_idamax( &F77_N, const_cast<double *>( x ), &F77_incX );
    return iamax ? iamax - 1 : 0;
}

/** viamax */

template<>
float OpenMPBLAS1::viamax( const IndexType n, const float* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "viamax<float>, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    int iamax;
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX;
#else
#define F77_N n
#define F77_incX incX
#endif
    iamax = F77_isamax( &F77_N, const_cast<float *>( x ), &F77_incX );
    return iamax ? static_cast<float>( ::fabs( x[( iamax - 1 ) * incX] ) ) : nIndex;
}

template<>
double OpenMPBLAS1::viamax( const IndexType n, const double* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "viamax<double>, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    int iamax;
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX;
#else
#define F77_N n
#define F77_incX incX
#endif
    iamax = F77_idamax( &F77_N, const_cast<double *>( x ), &F77_incX );
    return iamax ? ::fabs( x[( iamax - 1 ) * incX] ) : nIndex;
}

/** swap */

template<>
void OpenMPBLAS1::swap(
    const IndexType n,
    float* x,
    const IndexType incX,
    float* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "iamax<float>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#define F77_incY incY
#endif
    F77_sswap( &F77_N, x, &F77_incX, y, &F77_incY );
}

template<>
void OpenMPBLAS1::swap(
    const IndexType n,
    double* x,
    const IndexType incX,
    double* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "iamax<double>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#define F77_incY incY
#endif
    F77_dswap( &F77_N, x, &F77_incX, y, &F77_incY );
}

/** copy */

template<>
void OpenMPBLAS1::copy(
    const IndexType n,
    const float* x,
    const IndexType incX,
    float* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "copy<float>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#define F77_incY incY
#endif
    F77_scopy( &F77_N, x, &F77_incX, y, &F77_incY );
}

template<>
void OpenMPBLAS1::copy(
    const IndexType n,
    const double* x,
    const IndexType incX,
    double* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "copy<double>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#define F77_incY incY
#endif
    F77_dcopy( &F77_N, x, &F77_incX, y, &F77_incY );
}

/** axpy */

template<>
void OpenMPBLAS1::axpy(
    const IndexType n,
    const float alpha,
    const float* x,
    const IndexType incX,
    float* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "axpy<float>, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#define F77_incY incY
#endif
    F77_saxpy( &F77_N, &alpha, x, &F77_incX, y, &F77_incY );
}

template<>
void OpenMPBLAS1::axpy(
    const IndexType n,
    const double alpha,
    const double* x,
    const IndexType incX,
    double* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "axpy<double>, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#define F77_incY incY
#endif
    F77_daxpy( &F77_N, &alpha, x, &F77_incX, y, &F77_incY );
}

/** dot */

template<>
float OpenMPBLAS1::dot(
    const IndexType n,
    const float* x,
    const IndexType incX,
    const float* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "dot<float>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    float dot;
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#define F77_incY incY
#endif
    dot = F77_sdot( &F77_N, const_cast<float *>( x ), &F77_incX, const_cast<float *>( y ), &F77_incY );
    return dot;
}

template<>
double OpenMPBLAS1::dot(
    const IndexType n,
    const double* x,
    const IndexType incX,
    const double* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "dot<double>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    double dot;
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#define F77_incY incY
#endif
    dot = F77_ddot( &F77_N, const_cast<double *>( x ), &F77_incX, const_cast<double *>( y ), &F77_incY );
    return dot;
}

/** sum */

template<typename T>
void OpenMPBLAS1::sum( const IndexType n, T alpha, const T* x, T beta, const T* y, T* z, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "sum<T>, n = " << n << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y << ", z = " << z )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    #pragma omp parallel for
    for ( int i = 0; i < n; i++ )
    {
        z[i] = alpha * x[i] + beta * y[i];
    }
    return;
}

/** rot */

template<>
void OpenMPBLAS1::rot(
    const IndexType n,
    float* x,
    const IndexType incX,
    float* y,
    const IndexType incY,
    const float c,
    const float s,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "rot<float>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY << ", c = " << c << ", s = " << s )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#define F77_incY incY
#endif
    F77_srot( &F77_N, x, &F77_incX, y, &F77_incY, &c, &s );
    return;
}

template<>
void OpenMPBLAS1::rot(
    const IndexType n,
    double* x,
    const IndexType incX,
    double* y,
    const IndexType incY,
    const double c,
    const double s,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "rot<double>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY << ", c = " << c << ", s = " << s )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#define F77_incY incY
#endif
    F77_drot( &F77_N, x, &F77_incX, y, &F77_incY, &c, &s );
    return;
}

/** rotm */

template<>
void OpenMPBLAS1::rotm(
    const IndexType n,
    float* x,
    const IndexType incX,
    float* y,
    const IndexType incY,
    const float* P,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "rotm<float>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY << ", P = " << P )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#define F77_incY incY
#endif
    F77_srotm( &F77_N, x, &F77_incX, y, &F77_incY, P );
    return;
}

template<>
void OpenMPBLAS1::rotm(
    const IndexType n,
    double* x,
    const IndexType incX,
    double* y,
    const IndexType incY,
    const double* P,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "rotm<double>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY << ", P = " << P )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#define F77_incY incY
#endif
    F77_drotm( &F77_N, x, &F77_incX, y, &F77_incY, P );
    //      drotm( &n, x, &incX, y, &incY, P);
    return;
}

/** rotg */

template<>
void OpenMPBLAS1::rotg( float* a, float* b, float* c, float* s, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "rotg<float>, a = " << a << ", b = " << b << ", c = " << c << ", s = " << s )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    F77_srotg( a, b, c, s );
    return;
}

template<>
void OpenMPBLAS1::rotg( double* a, double* b, double* c, double* s, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "rotg<double>, a = " << a << ", b = " << b << ", c = " << c << ", s = " << s )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    F77_drotg( a, b, c, s );
    return;
}

/** rotmg */

template<>
void OpenMPBLAS1::rotmg( float* d1, float* d2, float* b1, const float b2, float* p, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "rotmg<float>, d1 = " << d1 << ", d2 = " << d2 << ", b1 = " << b1 << ", b2 = " << b2 << ", p = " << p )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    F77_srotmg( d1, d2, b1, &b2, p );
}

template<>
void OpenMPBLAS1::rotmg( double* d1, double* d2, double* b1, const double b2, double* p, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "rotmg<double>, d1 = " << d1 << ", d2 = " << d2 << ", b1 = " << b1 << ", b2 = " << b2 << ", p = " << p )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    F77_drotmg( d1, d2, b1, &b2, p );
    //      drotm( &n, x, &incX, y, &incY, P);
    return;
}

/* ass */

template<typename T>
void OpenMPBLAS1::ass( const IndexType n, const T value, T* x, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "ass<float>, n = " << n << ", value = " << value << ", x = " << x )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    #pragma omp parallel for
    for ( int i = 0; i < n; i++ )
    {
        x[i] = value;
    }
    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPBLAS1::setInterface( BLASInterface& BLAS )
{
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

    // LAMA_INTERFACE_REGISTER_T( BLAS, swap, float )
    // LAMA_INTERFACE_REGISTER_T( BLAS, swap, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, copy, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, copy, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, axpy, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, axpy, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, dot, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, dot, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, sum, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, sum, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, ass, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, ass, double )
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
