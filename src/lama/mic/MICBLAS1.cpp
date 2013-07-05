/**
 * @file MICBLAS1.cpp
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
 * @brief MICBLAS1.cpp
 * @author Thomas Brandes    
 * @date 05.07.2013
 * @since 1.1.0
 */

// hpp
#include <lama/mic/MICBLAS1.hpp>

// others
#include <lama/BLASInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

// tracing with LAMA_REGION
#include <lama/tracing.hpp>

#include <cmath>

namespace lama
{

LAMA_LOG_DEF_LOGGER( MICBLAS1::logger, "MIC.BLAS1" )

/** scal */

template<typename ValueType>
void MICBLAS1::scal( const IndexType n, const ValueType alpha, ValueType* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "scal<" << Scalar::getType<ValueType>() << ">, n = " << n 
                            << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "asynchronous execution for MIC not supported yet." )
    }
 
    LAMA_ASSERT_EQUAL_ERROR( incX, 1 );

    const size_t xPtr = ( size_t ) x;

    #pragma offload target( mic ) in( xPtr, n, alpha, incX )
    {
        ValueType* x = ( ValueType* ) xPtr;

        #pragma omp parallel for 
        for ( IndexType i = 0; i < n; ++i )
        {
            x[i * incX] *= alpha;
        }
    }
}

/* ------------------------------------------------------------------------- */

/** copy */

template<typename ValueType>
void MICBLAS1::copy(
    const IndexType n,
    const ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "copy<" << Scalar::getType<ValueType>() << ">, n = " << n 
                     << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "asynchronous execution for MIC not supported yet." )
    }

    LAMA_REGION( "MIC.BLAS1.copy" )

    const void* xPtr = x;
    void* yPtr = y;

    #pragma offload target( mic ), in( xPtr, yPtr, n, incX, incY )
    {
        const ValueType* x = ( ValueType* ) xPtr;
        ValueType* y = ( ValueType* ) yPtr;

        #pragma omp parallel for 
        for ( IndexType i = 0; i < n; ++i )
        {
            y[i * incY] = x[i * incX];
        }
    }
}

/* ------------------------------------------------------------------------- */

/** axpy */

template<typename ValueType>
void MICBLAS1::axpy(
    const IndexType n,
    const ValueType alpha,
    const ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_REGION( "MIC.BLAS1.axpy" )

    LAMA_LOG_INFO( logger, "axpy<" << Scalar::getType<ValueType>() << ",  n = " << n 
                      << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX 
                      << ", y = " << y << ", incY = " << incY )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "asynchronous execution for MIC not supported yet." )
    }

    const void* xPtr = x;
    void* yPtr = y;

    #pragma offload target( mic ), in( xPtr, yPtr, n, alpha, incX, incY )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        ValueType* y = static_cast<ValueType*>( yPtr );

        #pragma omp parallel for 
        for ( IndexType i = 0; i < n; ++i )
        {
            y[i * incY] += alpha * x[i * incX];
        }
    }
}

/** dot */

template<typename ValueType>
ValueType MICBLAS1::dot(
    const IndexType n,
    const ValueType* x,
    const IndexType incX,
    const ValueType* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_REGION( "MIC.BLAS1.dot" )

    LAMA_LOG_INFO( logger, "dot<" << Scalar::getType<ValueType>() << ">, n = " << n );

    if ( syncToken )
    {
        LAMA_THROWEXCEPTION( "no asynchronous execution for MIC possible at this level." )
    }

    ValueType val = 0;

    const void* xPtr = x;
    const void* yPtr = y;

    #pragma offload target( mic ), out( val ), in( xPtr, yPtr, n, incX, incY )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const ValueType* y = static_cast<const ValueType*>( yPtr );

        val = 0;

        #pragma omp parallel for reduction( +:val )
        for ( IndexType i = 0; i < n; ++i )
        {
            val += x[i * incX] * y[i * incY];
        }
    }

    return val;
}

/** sum */

template<typename ValueType>
void MICBLAS1::sum( const IndexType n, ValueType alpha, const ValueType* x, 
                                       ValueType beta, const ValueType* y, 
                                       ValueType* z, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "sum<" << Scalar::getType<ValueType>() << ">, n = " << n 
                             << ", alpha = " << alpha << ", x = " << x 
                             << ", beta = " << beta << ", y = " << y << ", z = " << z )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "asynchronous execution for MIC not supported yet." )
    }

    LAMA_REGION( "MIC.BLAS1.sum" )

    const void* xPtr = x;
    const void* yPtr = y;
    void* zPtr = z;

    #pragma offload target( mic ), in( xPtr, yPtr, zPtr, alpha, beta )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const ValueType* y = static_cast<const ValueType*>( yPtr );
        ValueType* z = static_cast<ValueType*>( zPtr );

        #pragma omp parallel for
        for ( int i = 0; i < n; i++ )
        {
            z[i] = alpha * x[i] + beta * y[i];
        }
    }

    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void MICBLAS1::setInterface( BLASInterface& BLAS )
{
    LAMA_LOG_INFO( logger, "set BLAS1 routines for MIC in Interface" )

    // Note: macro takes advantage of same name for routines and type definitions 
    //       ( e.g. routine CUDABLAS1::sum<T> is set for BLAS::BLAS1::sum variable

    LAMA_INTERFACE_REGISTER_T( BLAS, scal, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, scal, double )

    /*
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
    */

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

bool MICBLAS1::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::MIC );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool MICBLAS1::initialized = registerInterface();

} /** namespace lama */
