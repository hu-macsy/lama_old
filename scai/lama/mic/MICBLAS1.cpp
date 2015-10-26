/**
 * @file MICBLAS1.cpp
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
 * @brief MICBLAS1.cpp
 * @author Thomas Brandes
 * @date 05.07.2013
 * @since 1.1.0
 */

// hpp
#include <scai/lama/mic/MICBLAS1.hpp>

// local library
#include <scai/lama/BLASKernelTrait.hpp>

// internal scai libraries
#include <scai/hmemo/mic/MICContext.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tasking/SyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/ScalarType.hpp>

// external
#include <omp.h>

// std
#include <cmath>
#include <cstdlib>

namespace scai
{

using namespace hmemo;

namespace lama
{

SCAI_LOG_DEF_LOGGER( MICBLAS1::logger, "MIC.BLAS1" )

/* ------------------------------------------------------------------------- */

/** scal */

template<typename ValueType>
void MICBLAS1::scal(
    const IndexType n,
    const ValueType alpha,
    ValueType* x,
    const IndexType incX,
    tasking::SyncToken* syncToken )
{
    SCAI_LOG_DEBUG( logger,
                    "scal<" << common::getScalarType<ValueType>() << ">, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    if( n < 1 || incX < 1 )
    {
        return;
    }

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution for MIC not supported yet." )
    }

    void* xPtr = x;

    int device = MICContext::getCurrentDevice();

#pragma offload target( MIC : device ) in( xPtr, n, alpha, incX )
    {
        ValueType* x = static_cast<ValueType*>( xPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < n; ++i )
        {
            x[i * incX] *= alpha;
        }
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICBLAS1::asum( const IndexType n, const ValueType* x, const IndexType incX, tasking::SyncToken* syncToken )
{
    SCAI_LOG_DEBUG( logger,
                    "asum<" << common::getScalarType<ValueType>() << ">, n = " << n << ", x = " << x << ", incX = " << incX )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for mic possible at this level." )
    }

    ValueType asum = static_cast<ValueType>(0.0);

    if( n < 1 || incX < 1 )
    {
        return asum;
    }

    const void* xPtr = x;

    int device = MICContext::getCurrentDevice();

#pragma offload target( MIC : device ) in( xPtr, n, incX ), out( asum )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );

        asum = static_cast<ValueType>(0.0);

        #pragma omp parallel for reduction( + : asum )

        for( int i = 0; i < n; ++i )
        {
            asum += std::abs( x[i * incX] );
        }
    }

    return asum;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType MICBLAS1::iamax( const IndexType n, const ValueType* x, const IndexType incX, tasking::SyncToken* syncToken )
{
    SCAI_LOG_INFO( logger,
                   "iamax<" << common::getScalarType<ValueType>() << " >, n = " << n << ", x = " << x << ", incX = " << incX )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for mic possible at this level." )
    }

    IndexType maxIndex = 0;

    if( n < 1 || incX < 1 )
    {
        return maxIndex;
    }

    const void* xPtr = x;

    int device = MICContext::getCurrentDevice();

#pragma offload target( MIC : device ) in( xPtr, n, incX ), out( maxIndex )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );

        maxIndex = 0;

        ValueType maxVal = -std::numeric_limits<ValueType>::max();

        #pragma omp parallel
        {
            IndexType threadMaxIndex = -1;
            ValueType threadMaxVal = -std::numeric_limits<ValueType>::max();

            #pragma omp for

            for( int i = 0; i < n; ++i )
            {
                const ValueType& val = x[i * incX];

                if( val > threadMaxVal )
                {
                    threadMaxIndex = i;
                    threadMaxVal = val;
                }
            }

            // ordered reduction needed to get smallest index

            #pragma omp for ordered

            for( int nt = 0; nt < omp_get_num_threads(); ++nt )
            {
                #pragma omp ordered
                {
                    if( threadMaxVal > maxVal )
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
void MICBLAS1::swap(
    const IndexType n,
    ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY,
    tasking::SyncToken* syncToken )
{
    SCAI_LOG_DEBUG( logger,
                    "iamax<long double>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    if( n < 1 || incX < 1 || incY < 1 )
    {
        return;
    }

    void* xPtr = x;
    void* yPtr = y;

    int device = MICContext::getCurrentDevice();

#pragma offload target( MIC : device ) in( xPtr, yPtr, incX, incY, n )
    {
        ValueType* x = static_cast<ValueType*>( xPtr );
        ValueType* y = static_cast<ValueType*>( yPtr );

        #pragma omp parallel for

        for( int i = 0; i < n; ++i )
        {
            std::swap( x[i * incX], y[i * incY] );
        }
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICBLAS1::nrm2( const IndexType n, const ValueType* x, const IndexType incX, tasking::SyncToken* syncToken )
{
    SCAI_LOG_INFO( logger, "nrm2<" << common::getScalarType<ValueType>() << ">( n = " << n << " )" )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for mic possible at this level." )
    }

    const void* xPtr = x;

    ValueType sum = static_cast<ValueType>(0.0);

    if( n < 1 || incX < 1 )
    {
        return sum;
    }

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( xPtr, n, incX ), out( sum )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );

        sum = static_cast<ValueType>(0.0);

        #pragma omp parallel for reduction( + : sum )

        for( int i = 0; i < n; ++i )
        {
            sum += x[i * incX] * x[i * incX];
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
    tasking::SyncToken* syncToken )
{
    SCAI_LOG_DEBUG( logger,
                    "copy<" << common::getScalarType<ValueType>() << ">, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution for MIC not supported yet." )
    }

    // SCAI_REGION( "MIC.BLAS1.copy" )

    if( n < 1 || incX < 1 || incY < 1 )
    {
        return;
    }

    const void* xPtr = x;
    void* yPtr = y;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( xPtr, yPtr, n, incX, incY )
    {
        const ValueType* x = (ValueType*) xPtr;
        ValueType* y = (ValueType*) yPtr;

        #pragma omp parallel for

        for( IndexType i = 0; i < n; ++i )
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
    tasking::SyncToken* syncToken )
{
    // SCAI_REGION( "MIC.BLAS1.axpy" )

    SCAI_LOG_INFO( logger,
                   "axpy<" << common::getScalarType<ValueType>() << ",  n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution for MIC not supported yet." )
    }

    if( n < 1 || incX < 1 || incY < 1 )
    {
        return;
    }

    const void* xPtr = x;
    void* yPtr = y;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( xPtr, yPtr, n, alpha, incX, incY )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        ValueType* y = static_cast<ValueType*>( yPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < n; ++i )
        {
            y[i * incY] += alpha * x[i * incX];
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
    tasking::SyncToken* syncToken )
{
    // SCAI_REGION( "MIC.BLAS1.dot" )

    SCAI_LOG_INFO( logger, "dot<" << common::getScalarType<ValueType>() << ">" 
                           << ", n = " << n << ", incX = " << incX << ", incY = " << incY );

    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "no asynchronous execution for MIC possible at this level." )
    }

    ValueType val = static_cast<ValueType>(0.0);

    if ( n < 1 || incX < 1 || incY < 1 )
    {
        return val;
    }

    int device = MICContext::getCurrentDevice();

    const void* xPtr = x;
    const void* yPtr = y;

#pragma offload target( MIC : device ), out( val ), in( xPtr, yPtr, n, incX, incY )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const ValueType* y = static_cast<const ValueType*>( yPtr );

        val = static_cast<ValueType>(0.0);

        #pragma omp parallel for reduction( +:val )

        for( IndexType i = 0; i < n; ++i )
        {
            val += x[i * incX] * y[i * incY];
        }
    }

    SCAI_LOG_INFO( logger, "dot: result = " << val )

    return val;
}

/** sum */

template<typename ValueType>
void MICBLAS1::sum(
    const IndexType n,
    ValueType alpha,
    const ValueType* x,
    ValueType beta,
    const ValueType* y,
    ValueType* z,
    tasking::SyncToken* syncToken )
{
    SCAI_LOG_DEBUG( logger,
                    "sum<" << common::getScalarType<ValueType>() << ">, n = " << n << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y << ", z = " << z )

    if( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution for MIC not supported yet." )
    }

    // SCAI_REGION( "MIC.BLAS1.sum" )

    const void* xPtr = x;
    const void* yPtr = y;
    void* zPtr = z;

    int device = MICContext::getCurrentDevice();

#pragma offload target( MIC : device ), in( xPtr, yPtr, zPtr, alpha, beta )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const ValueType* y = static_cast<const ValueType*>( yPtr );
        ValueType* z = static_cast<ValueType*>( zPtr );

        #pragma omp parallel for

        for( int i = 0; i < n; i++ )
        {
            z[i] = alpha * x[i] + beta * y[i];
        }
    }

    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void MICBLAS1::registerKernels()
{
    SCAI_LOG_INFO( logger, "register BLAS1 kernels for MIC in Kernel Registry" )

    using namespace scai::kregistry;

    // ctx will contain the context for which registration is done, here MIC

    common::ContextType ctx = common::context::MIC;

    KernelRegistry::set<BLASKernelTrait::scal<float> >( scal, ctx );
    KernelRegistry::set<BLASKernelTrait::scal<double> >( scal, ctx );

    KernelRegistry::set<BLASKernelTrait::nrm2<float> >( nrm2, ctx );
    KernelRegistry::set<BLASKernelTrait::nrm2<double> >( nrm2, ctx );

    KernelRegistry::set<BLASKernelTrait::asum<float> >( asum, ctx );
    KernelRegistry::set<BLASKernelTrait::asum<double> >( asum, ctx );

    KernelRegistry::set<BLASKernelTrait::iamax<float> >( iamax, ctx );
    KernelRegistry::set<BLASKernelTrait::iamax<double> >( iamax, ctx );

    KernelRegistry::set<BLASKernelTrait::swap<float> >( swap, ctx );
    KernelRegistry::set<BLASKernelTrait::swap<double> >( swap, ctx );

    KernelRegistry::set<BLASKernelTrait::copy<float> >( copy, ctx );
    KernelRegistry::set<BLASKernelTrait::copy<double> >( copy, ctx );

    KernelRegistry::set<BLASKernelTrait::axpy<float> >( axpy, ctx );
    KernelRegistry::set<BLASKernelTrait::axpy<double> >( axpy, ctx );

    KernelRegistry::set<BLASKernelTrait::dot<float> >( dot, ctx );
    KernelRegistry::set<BLASKernelTrait::dot<double> >( dot, ctx );

    KernelRegistry::set<BLASKernelTrait::sum<float> >( sum, ctx );
    KernelRegistry::set<BLASKernelTrait::sum<double> >( sum, ctx );
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the BLAS1 routines                                */
/* --------------------------------------------------------------------------- */

bool MICBLAS1::registerInterface()
{
    registerKernels();
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool MICBLAS1::initialized = registerInterface();

} /* end namespace lama */

} /* end namespace scai */
