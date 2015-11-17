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
#include <scai/blaskernel/mic/MICBLAS1.hpp>

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>

// internal scai libraries
#include <scai/hmemo/mic/MICContext.hpp>
#include <scai/hmemo/mic/MICSyncToken.hpp>

#include <scai/kregistry/KernelRegistry.hpp>

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
using tasking::MICSyncToken;

namespace blaskernel
{

SCAI_LOG_DEF_LOGGER( MICBLAS1::logger, "MIC.BLAS1" )

/* ------------------------------------------------------------------------- */

/** scal */

template<typename ValueType>
void MICBLAS1::scal(
    const IndexType n,
    const ValueType alpha,
    ValueType* x,
    const IndexType incX )
{
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
    }

    SCAI_LOG_DEBUG( logger,
                    "scal<" << common::getScalarType<ValueType>() << ">, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    if( n < 1 || incX < 1 )
    {
        return;
    }

    void* xPtr = x;

    int device = MICContext::getCurrentDevice();

    const ValueType* alphaPtr = &alpha;

#pragma offload target( MIC : device ) in( xPtr, n, alphaPtr[0:1], incX )
    {
        ValueType* x = static_cast<ValueType*>( xPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < n; ++i )
        {
            x[i * incX] *= *alphaPtr;
        }
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICBLAS1::asum( const IndexType n, const ValueType* x, const IndexType incX )
{
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
    }

    SCAI_LOG_DEBUG( logger,
                    "asum<" << common::getScalarType<ValueType>() << ">, n = " << n << ", x = " << x << ", incX = " << incX )

    ValueType asum = static_cast<ValueType>(0.0);

    if( n < 1 || incX < 1 )
    {
        return asum;
    }

    const void* xPtr = x;

    int device = MICContext::getCurrentDevice();

    ValueType* asumPtr = &asum;

#pragma offload target( MIC : device ) in( xPtr, n, incX ), out( asumPtr[0:1] )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );

        *asumPtr = static_cast<ValueType>(0.0);

        #pragma omp parallel
	{
		ValueType local_asum = static_cast<ValueType>( 0.0 );

		#pragma omp for
	        for( int i = 0; i < n; ++i )
	        {
	            local_asum += ::fabs( x[i * incX] );
        	}

		#pragma omp critical
		*asumPtr += local_asum;
	}
    }

    return asum;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType MICBLAS1::iamax( const IndexType n, const ValueType* x, const IndexType incX )
{
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
    }

    SCAI_LOG_INFO( logger,
                   "iamax<" << common::getScalarType<ValueType>() << " >, n = " << n << ", x = " << x << ", incX = " << incX )

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
    const IndexType incY )
{
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
    }

    SCAI_LOG_DEBUG( logger,
                    "iamax<long double>, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

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
ValueType MICBLAS1::nrm2( const IndexType n, const ValueType* x, const IndexType incX )
{
    SCAI_LOG_INFO( logger, "nrm2<" << common::getScalarType<ValueType>() << ">( n = " << n << " )" )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
    }

    const void* xPtr = x;

    ValueType sum = static_cast<ValueType>(0.0);

    if( n < 1 || incX < 1 )
    {
        return sum;
    }

    int device = MICContext::getCurrentDevice();
    ValueType* sumPtr = &sum;

#pragma offload target( mic : device ) in( xPtr, n, incX ), out( sumPtr[0:1] )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );

        *sumPtr = static_cast<ValueType>(0.0);

        #pragma omp parallel 
	{
		ValueType local_sum = static_cast<ValueType>( 0.0 );
		#pragma omp for
        	for( int i = 0; i < n; ++i )
	        {
	            local_sum += x[i * incX] * x[i * incX];
	        }

		#pragma omp critical
		*sumPtr += local_sum;
	}
    }

    return sqrt( sum );
}

/* ------------------------------------------------------------------------- */

/** copy */

template<typename ValueType>
void MICBLAS1::copy(
    const IndexType n,
    const ValueType* x,
    const IndexType incX,
    ValueType* y,
    const IndexType incY )
{
    SCAI_LOG_DEBUG( logger,
                    "copy<" << common::getScalarType<ValueType>() << ">, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
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
    const IndexType incY )
{
    // SCAI_REGION( "MIC.BLAS1.axpy" )

    SCAI_LOG_INFO( logger,
                   "axpy<" << common::getScalarType<ValueType>() << ",  n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
    }

    if( n < 1 || incX < 1 || incY < 1 )
    {
        return;
    }

    const void* xPtr = x;
    void* yPtr = y;

    int device = MICContext::getCurrentDevice();

    const ValueType* alphaPtr = &alpha;

#pragma offload target( mic : device ), in( xPtr, yPtr, n, alphaPtr[0:1], incX, incY )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        ValueType* y = static_cast<ValueType*>( yPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < n; ++i )
        {
            y[i * incY] += *alphaPtr * x[i * incX];
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
    const IndexType incY )
{
    // SCAI_REGION( "MIC.BLAS1.dot" )

    SCAI_LOG_INFO( logger, "dot<" << common::getScalarType<ValueType>() << ">" 
                           << ", n = " << n << ", incX = " << incX << ", incY = " << incY );

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
    }

    ValueType val = static_cast<ValueType>(0.0);

    if ( n < 1 || incX < 1 || incY < 1 )
    {
        return val;
    }

    int device = MICContext::getCurrentDevice();

    const void* xPtr = x;
    const void* yPtr = y;

    ValueType* valPtr = &val;

#pragma offload target( MIC : device ), out( valPtr[0:1] ), in( xPtr, yPtr, n, incX, incY )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const ValueType* y = static_cast<const ValueType*>( yPtr );

        *valPtr = static_cast<ValueType>(0.0);

        #pragma omp parallel
	{
		ValueType local_val = static_cast<ValueType>( 0.0 );
		*valPtr  = static_cast<ValueType>( 0.0 );
		
		#pragma omp for
        	for( IndexType i = 0; i < n; ++i )
        	{
	            local_val += x[i * incX] * y[i * incY];
	        }

		#pragma omp critical
		*valPtr += local_val;
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
    ValueType* z )
{
    SCAI_LOG_DEBUG( logger,
                    "sum<" << common::getScalarType<ValueType>() << ">, n = " << n << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y << ", z = " << z )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
    }

    // SCAI_REGION( "MIC.BLAS1.sum" )

    const void* xPtr = x;
    const void* yPtr = y;
    void* zPtr = z;
    const ValueType* alphaPtr = &alpha;
    const ValueType* betaPtr = &beta;

    int device = MICContext::getCurrentDevice();

#pragma offload target( MIC : device ), in( xPtr, yPtr, zPtr, alphaPtr[0:1], betaPtr[0:1] )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const ValueType* y = static_cast<const ValueType*>( yPtr );
        ValueType* z = static_cast<ValueType*>( zPtr );

        #pragma omp parallel for

        for( int i = 0; i < n; i++ )
        {
            z[i] = *alphaPtr * x[i] + *betaPtr * y[i];
        }
    }

    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void MICBLAS1::registerKernels( bool deleteFlag )
{
    SCAI_LOG_INFO( logger, "register BLAS1 kernels for MIC in Kernel Registry" )

    using kregistry::KernelRegistry;
    using common::context::MIC;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD ;   // add it or delete it

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

/*    KernelRegistry::set<BLASKernelTrait::scal<float> >( scal, MIC, flag );
    KernelRegistry::set<BLASKernelTrait::scal<double> >( scal, MIC, flag );

    KernelRegistry::set<BLASKernelTrait::scal<ComplexFloat> >( scal, MIC, flag );
    KernelRegistry::set<BLASKernelTrait::scal<ComplexDouble> >( scal, MIC, flag );

    KernelRegistry::set<BLASKernelTrait::nrm2<float> >( nrm2, MIC, flag );
    KernelRegistry::set<BLASKernelTrait::nrm2<double> >( nrm2, MIC, flag );

    KernelRegistry::set<BLASKernelTrait::asum<float> >( asum, MIC, flag );
    KernelRegistry::set<BLASKernelTrait::asum<double> >( asum, MIC, flag );

    KernelRegistry::set<BLASKernelTrait::iamax<float> >( iamax, MIC, flag );
    KernelRegistry::set<BLASKernelTrait::iamax<double> >( iamax, MIC, flag );

    KernelRegistry::set<BLASKernelTrait::swap<float> >( swap, MIC, flag );
    KernelRegistry::set<BLASKernelTrait::swap<double> >( swap, MIC, flag );

    KernelRegistry::set<BLASKernelTrait::copy<float> >( copy, MIC, flag );
    KernelRegistry::set<BLASKernelTrait::copy<double> >( copy, MIC, flag );

    KernelRegistry::set<BLASKernelTrait::axpy<float> >( axpy, MIC, flag );
    KernelRegistry::set<BLASKernelTrait::axpy<double> >( axpy, MIC, flag );

    KernelRegistry::set<BLASKernelTrait::dot<float> >( dot, MIC, flag );
    KernelRegistry::set<BLASKernelTrait::dot<double> >( dot, MIC, flag );

    KernelRegistry::set<BLASKernelTrait::sum<float> >( sum, MIC, flag );
    KernelRegistry::set<BLASKernelTrait::sum<double> >( sum, MIC, flag );
*/
    
#define LAMA_BLAS1_REGISTER(z, I, _)                                                             \
        KernelRegistry::set<BLASKernelTrait::scal<ARITHMETIC_HOST_TYPE_##I> >( scal, MIC, flag );    \
        KernelRegistry::set<BLASKernelTrait::nrm2<ARITHMETIC_HOST_TYPE_##I> >( nrm2, MIC, flag );    \
        KernelRegistry::set<BLASKernelTrait::asum<ARITHMETIC_HOST_TYPE_##I> >( asum, MIC, flag );    \
        KernelRegistry::set<BLASKernelTrait::iamax<ARITHMETIC_HOST_TYPE_##I> >( iamax, MIC, flag );  \
        KernelRegistry::set<BLASKernelTrait::swap<ARITHMETIC_HOST_TYPE_##I> >( swap, MIC, flag );    \
        KernelRegistry::set<BLASKernelTrait::copy<ARITHMETIC_HOST_TYPE_##I> >( copy, MIC, flag );    \
        KernelRegistry::set<BLASKernelTrait::axpy<ARITHMETIC_HOST_TYPE_##I> >( axpy, MIC, flag );    \
        KernelRegistry::set<BLASKernelTrait::dot<ARITHMETIC_HOST_TYPE_##I> >( dot, MIC, flag );      \
    	KernelRegistry::set<BLASKernelTrait::sum<ARITHMETIC_HOST_TYPE_##I> >( sum, MIC, flag );

        BOOST_PP_REPEAT( ARITHMETIC_HOST_EXT_TYPE_CNT, LAMA_BLAS1_REGISTER, _ )

#undef LAMA_BLAS1_REGISTER

}

/* --------------------------------------------------------------------------- */
/*    Static intialization with registration                                   */
/* --------------------------------------------------------------------------- */

MICBLAS1::RegisterGuard::RegisterGuard()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

MICBLAS1::RegisterGuard::~RegisterGuard()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

MICBLAS1::RegisterGuard MICBLAS1::guard;    // guard variable for registration

} /* end namespace blaskernel */

} /* end namespace scai */
