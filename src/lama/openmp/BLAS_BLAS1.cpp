/**
 * @file BLAS_BLAS1.cpp
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
 * @brief Implementation of LAMA BLAS1 routines using wrappers to BLAS library
 * @author Lauretta Schubert
 * @date 05.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/openmp/BLAS_BLAS1.hpp>

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

LAMA_LOG_DEF_LOGGER( BLAS_BLAS1::logger, "BLAS_.BLAS1" )

/** scal */

template<typename T>
static inline void wrapperScal( IndexType n, T alpha, T* x_d, IndexType incX );

template<>
void wrapperScal( IndexType n, float alpha, float* x, IndexType incX )
{
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX;
#else
#define F77_N n
#define F77_incX incX
#endif
    F77_sscal( &F77_N, &alpha, x, &F77_incX );
}

template<>
void wrapperScal( IndexType n, double alpha, double* x, IndexType incX )
{
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX;
#else
#define F77_N n
#define F77_incX incX
#endif
    F77_dscal( &F77_N, &alpha, x, &F77_incX );
}

template<typename T>
void BLAS_BLAS1::scal( const IndexType n, const T alpha, T* x, const IndexType incX, SyncToken* syncToken )
{
	if ( incX <= 0 )
	{
		return;
	}

	LAMA_REGION( "BLAS_.BLAS1.scal" )

    LAMA_LOG_DEBUG( logger, "scal<" << Scalar::getType<T>() << " n = " << n 
                    << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    wrapperScal( n, alpha, x, incX );
}

/** nrm2 */

template<typename T>
static inline T wrapperNrm2( IndexType n, const T* x, IndexType incX );

template<>
float wrapperNrm2( IndexType n, const float* x, IndexType incX )
{
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
double wrapperNrm2( IndexType n, const double* x, IndexType incX )
{
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
}

template<typename T>
float BLAS_BLAS1::nrm2( const IndexType n, const float* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_REGION( "BLAS.BLAS1.nrm2" )

    LAMA_LOG_DEBUG( logger, "nrm2<" << Scalar::getType<T>() << ">, n = " << n 
                    << ", x = " << x << ", incX = " << incX )

    if ( incX <= 0 )
    {
    	return static_cast<T>( 0 );
    }

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    return wrapperNrm2( n, x, incX);
}

/** asum (l1 norm) */

template<>
float BLAS_BLAS1::asum( const IndexType n, const float* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "asum<float>, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    float asum = 0.0;
    if (incX > 0)
    {
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#endif
     asum = F77_sasum( &F77_N, (float*) x, &F77_incX );
    }
    return asum;
}

template<>
double BLAS_BLAS1::asum( const IndexType n, const double* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "asum<double>, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    double asum = 0.0;
    if (incX > 0)
    {
#ifdef F77_INT
    F77_INT F77_N = n, F77_incX = incX, F77_incY = incY;
#else
#define F77_N n
#define F77_incX incX
#endif
    	asum = F77_dasum( &F77_N, (double*) x, &F77_incX );
    }
    return asum;
}

/** iamax */

template<>
IndexType BLAS_BLAS1::iamax( const IndexType n, const float* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "iamax<float>, n = " << n << ", x = " << x << ", incX = " << incX )

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
IndexType BLAS_BLAS1::iamax( const IndexType n, const double* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "iamax<double>, n = " << n << ", x = " << x << ", incX = " << incX )

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

/** swap */

template<>
void BLAS_BLAS1::swap(
    const IndexType n,
    float* x,
    const IndexType incX,
    float* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "iamax<float>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

	if ( (incX <= 0) || (incY <= 0) )
	{
		return;
	}

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
void BLAS_BLAS1::swap(
    const IndexType n,
    double* x,
    const IndexType incX,
    double* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "iamax<double>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

	if ( (incX <= 0) || (incY <= 0) )
	{
		return;
	}

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
void BLAS_BLAS1::copy(
    const IndexType n,
    const float* x,
    const IndexType incX,
    float* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "copy<float>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

	if ( (incX <= 0) || (incY <= 0) )
	{
		return;
	}

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
void BLAS_BLAS1::copy(
    const IndexType n,
    const double* x,
    const IndexType incX,
    double* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "copy<double>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

	if ( (incX <= 0) || (incY <= 0) )
	{
		return;
	}

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
void BLAS_BLAS1::axpy(
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

	if ( (incX <= 0) || (incY <= 0) )
	{
		return;
	}

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
void BLAS_BLAS1::axpy(
    const IndexType n,
    const double alpha,
    const double* x,
    const IndexType incX,
    double* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_REGION( "BLAS_.BLAS1.daxpy" )

    LAMA_LOG_DEBUG( logger,
                    "axpy<double>, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

	if ( (incX <= 0) || (incY <= 0) )
	{
		return;
	}

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
float BLAS_BLAS1::dot(
    const IndexType n,
    const float* x,
    const IndexType incX,
    const float* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_REGION( "BLAS_.BLAS1.sdot" )

	LAMA_LOG_DEBUG( logger,
					"dot<float>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

	if ( (incX <= 0) || (incY <= 0) )
	{
		return 0.0;
	}

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
double BLAS_BLAS1::dot(
    const IndexType n,
    const double* x,
    const IndexType incX,
    const double* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_REGION( "BLAS_.BLAS1.ddot" )

	LAMA_LOG_DEBUG( logger,
					"dot<double>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

	if ( (incX <= 0) || (incY <= 0) )
	{
		return 0.0;
	}

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
void BLAS_BLAS1::sum( const IndexType n, T alpha, const T* x, T beta, const T* y, T* z, SyncToken* syncToken )
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

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void BLAS_BLAS1::setInterface( BLASInterface& BLAS )
{
    LAMA_LOG_INFO( logger, "set BLAS1 wrapper routines for Host Context in Interface" )

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

bool BLAS_BLAS1::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::Host );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool BLAS_BLAS1::initialized = registerInterface();

} /** namespace lama */
