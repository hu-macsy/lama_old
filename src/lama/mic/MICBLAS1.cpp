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
#include <omp.h>

#include <cmath>
#include <cstdlib>

namespace lama
{

LAMA_LOG_DEF_LOGGER( MICBLAS1::logger, "MIC.BLAS1" )

/* ------------------------------------------------------------------------- */

/** scal */

template<typename ValueType>
void MICBLAS1::scal( const IndexType n, const ValueType alpha, ValueType* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "scal<" << Scalar::getType<ValueType>() << ">, n = " << n 
                            << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    if ( n < 1 || incX < 1 )
    {
        return; 
    }

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "asynchronous execution for MIC not supported yet." )
    }
 
    void* xPtr = x;

    #pragma offload target( mic ) in( xPtr, n, alpha, incX )
    {
        ValueType* x = static_cast<ValueType*>( xPtr );

        #pragma omp parallel for 
        for ( IndexType i = 0; i < n; ++i )
        {
            x[i * incX] *= alpha;
        }
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICBLAS1::asum( const IndexType n, const ValueType* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger, "asum<" << Scalar::getType<ValueType>() << ">, n = " << n << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for mic possible at this level." )
    }

    ValueType asum = 0;

    if ( n < 1 || incX < 1 )
    {
        return asum;
    }

    const void* xPtr = x;

    #pragma offload target( mic ) in( xPtr, n, incX ), out( asum )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );

        asum = 0;

        #pragma omp parallel for reduction( + : asum )
        for ( int i = 0; i < n; ++i )
        {
            asum += std::abs( x[ i * incX ] );
        }
    }

    return asum;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType MICBLAS1::iamax( const IndexType n, const ValueType* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "iamax<" << Scalar::getType<ValueType>() << " >, n = " << n 
                            << ", x = " << x << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for mic possible at this level." )
    }

    IndexType maxIndex = 0;

    if ( n < 1 || incX < 1 )
    {
        return  maxIndex;
    }

    const void* xPtr = x;

    #pragma offload target( mic ) in( xPtr, n, incX ), out( maxIndex )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );

        maxIndex = 0;

        ValueType maxVal = - std::numeric_limits<ValueType>::max();

        #pragma omp parallel
        {  
            IndexType threadMaxIndex = -1;
            ValueType threadMaxVal   = - std::numeric_limits<ValueType>::max();

            #pragma omp for 
            for ( int i = 0; i < n; ++i )
            {
                const ValueType& val = x[ i * incX ];

                if ( val > threadMaxVal )
                {
                    threadMaxIndex = i;
                    threadMaxVal   = val;
                }
            }

            // ordered reduction needed to get smallest index
            
            #pragma omp for ordered
            for ( int nt = 0; nt < omp_get_num_threads(); ++nt )
            {
                #pragma omp ordered
                {
                    if ( threadMaxVal > maxVal )
                    {
                        maxVal = threadMaxVal;
                        maxIndex = threadMaxIndex;
                    }
                }
            }
        }
    }

    return maxIndex;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICBLAS1::viamax( const IndexType n, const ValueType* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "viamax<" << Scalar::getType<ValueType>() << " >, n = " << n  )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for mic possible at this level." )
    }

    ValueType maxVal = 1;  // default value for error cases

    if ( n < 1 || incX < 1 )
    {
        return  maxVal;
    }

    const void* xPtr = x;

    #pragma offload target( mic ) in( xPtr, n, incX ), inout( maxVal )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );

        #pragma omp parallel
        {  
            ValueType threadMaxVal   = - std::numeric_limits<ValueType>::max();

            #pragma omp for 
            for ( int i = 0; i < n; ++i )
            {
                const ValueType& val = x[ i * incX ];

                if ( val > threadMaxVal )
                {
                    threadMaxVal   = val;
                }
            }

            #pragma omp critical
            {
                if ( threadMaxVal > maxVal )
                {
                    maxVal = threadMaxVal;
                }
            }
        }
    }

    return maxVal;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MICBLAS1::swap(
    const IndexType n,
    ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_DEBUG( logger,
                    "iamax<long double>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    if ( n < 1 || incX < 1 || incY < 1 )
    {
        return;
    }

    void* xPtr = x;
    void* yPtr = y;

    #pragma offload target( MIC ) in( xPtr, yPtr, incX, incY, n )
    {
        ValueType* x = static_cast<ValueType*>( xPtr );
        ValueType* y = static_cast<ValueType*>( yPtr );

        #pragma omp parallel for
        for( int i = 0; i < n; ++i )
        {
            std::swap( x[ i * incX ], y[ i * incY ] );
        }
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICBLAS1::nrm2( const IndexType n, const ValueType* x, const IndexType incX, SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "nrm2<" << Scalar::getType<ValueType>() << ">( n = " << n << " )" )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for mic possible at this level." )
    }

    const void* xPtr = x;

    ValueType sum = 0;

    if ( n < 1 || incX < 1 )
    {
        return sum;
    }

    #pragma offload target( mic ) in( xPtr, n, incX ), out( sum )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );

        sum = 0;

        #pragma omp parallel for reduction( + : sum )
        for( int i = 0; i < n; ++i )
        {
            sum += x[ i * incX ] * x[ i * incX ];
        }
    }

    return std::sqrt( sum );
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

    if ( n < 1 || incX < 1 || incY < 1 )
    {
        return;
    }

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

    if ( n < 1 || incX < 1 || incY < 1 )
    {
        return;
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
            y[ i * incY ] += alpha * x[ i * incX ];
        }
    }
}

/* ------------------------------------------------------------------------- */
/*  dot                                                                      */
/* ------------------------------------------------------------------------- */

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

    if ( n < 1 || incX < 1 || incY < 1 )
    {
        return val;
    }

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
