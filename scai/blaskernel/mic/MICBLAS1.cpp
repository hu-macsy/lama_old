/**
 * @file MICBLAS1.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief MICBLAS1.cpp
 * @author Thomas Brandes
 * @date 05.07.2013
 */

// hpp
#include <scai/blaskernel/mic/MICBLAS1.hpp>

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>

// internal scai libraries
#include <scai/hmemo/mic/MICContext.hpp>
#include <scai/tasking/mic/MICSyncToken.hpp>

#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/OpenMP.hpp>

// std
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
                    "scal<" << common::TypeTraits<ValueType>::id() << ">, n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX )

    if ( n < 1 || incX < 1 )
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

        for ( IndexType i = 0; i < n; ++i )
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
                    "asum<" << common::TypeTraits<ValueType>::id() << ">, n = " << n << ", x = " << x << ", incX = " << incX )
    ValueType asum = static_cast<ValueType>( 0.0 );

    if ( n < 1 || incX < 1 )
    {
        return asum;
    }

    const void* xPtr = x;
    int device = MICContext::getCurrentDevice();
    ValueType* asumPtr = &asum;
#pragma offload target( MIC : device ) in( xPtr, n, incX ), out( asumPtr[0:1] )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        *asumPtr = static_cast<ValueType>( 0.0 );
        #pragma omp parallel
        {
            ValueType local_asum = static_cast<ValueType>( 0.0 );
            #pragma omp for

            for ( int i = 0; i < n; ++i )
            {
                local_asum += common::Math::abs( common::Math::real( x[i * incX] ) ) + common::Math::abs( common::Math::imag( x[i * incX] ) );
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
                   "iamax<" << common::TypeTraits<ValueType>::id() << " >, n = " << n << ", x = " << x << ", incX = " << incX )
    IndexType maxIndex = 0;

    if ( n < 1 || incX < 1 )
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

            for ( int i = 0; i < n; ++i )
            {
                const ValueType& val = x[i * incX];

                if ( val > threadMaxVal )
                {
                    threadMaxIndex = i;
                    threadMaxVal = val;
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

    if ( n < 1 || incX < 1 || incY < 1 )
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

        for ( int i = 0; i < n; ++i )
        {
            std::swap( x[i * incX], y[i * incY] );
        }
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICBLAS1::nrm2( const IndexType n, const ValueType* x, const IndexType incX )
{
    SCAI_LOG_INFO( logger, "nrm2<" << common::TypeTraits<ValueType>::id() << ">( n = " << n << " )" )
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
    }

    const void* xPtr = x;
    ValueType sum = static_cast<ValueType>( 0.0 );

    if ( n < 1 || incX < 1 )
    {
        return sum;
    }

    int device = MICContext::getCurrentDevice();
    ValueType* sumPtr = &sum;
#pragma offload target( mic : device ) in( xPtr, n, incX ), out( sumPtr[0:1] )
    {
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        *sumPtr = static_cast<ValueType>( 0.0 );
        #pragma omp parallel
        {
            ValueType local_sum = static_cast<ValueType>( 0.0 );
            #pragma omp for

            for ( int i = 0; i < n; ++i )
            {
                local_sum += x[i * incX] * common::Math::conj( x[i * incX] );
            }

            #pragma omp critical
            *sumPtr += local_sum;
        }
    }
    return common::Math::sqrt( sum );
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
                    "copy<" << common::TypeTraits<ValueType>::id() << ">, n = " << n << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
    }

    // SCAI_REGION( "MIC.BLAS1.copy" )

    if ( n < 1 || incX < 1 || incY < 1 )
    {
        return;
    }

    const void* xPtr = x;
    void* yPtr = y;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ), in( xPtr, yPtr, n, incX, incY )
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
    const IndexType incY )
{
    // SCAI_REGION( "MIC.BLAS1.axpy" )
    SCAI_LOG_INFO( logger,
                   "axpy<" << common::TypeTraits<ValueType>::id() << ",  n = " << n << ", alpha = " << alpha << ", x = " << x << ", incX = " << incX << ", y = " << y << ", incY = " << incY )
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
    }

    if ( n < 1 || incX < 1 || incY < 1 )
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

        for ( IndexType i = 0; i < n; ++i )
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
    SCAI_LOG_INFO( logger, "dot<" << common::TypeTraits<ValueType>::id() << ">"
                   << ", n = " << n << ", incX = " << incX << ", incY = " << incY );
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for MIC not supported yet." )
    }

    ValueType val = static_cast<ValueType>( 0.0 );

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
        *valPtr = static_cast<ValueType>( 0.0 );
        #pragma omp parallel
        {
            ValueType local_val = static_cast<ValueType>( 0.0 );
            *valPtr  = static_cast<ValueType>( 0.0 );
            #pragma omp for

            for ( IndexType i = 0; i < n; ++i )
            {
                local_val += common::Math::conj( x[i * incX] ) * y[i * incY];
            }

            #pragma omp critical
            *valPtr += local_val;
        }
    }

    SCAI_LOG_DEBUG( logger, "dot: result = " << val )

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
    SCAI_LOG_INFO( logger,
                    "sum<" << common::TypeTraits<ValueType>::id() << ">, z[" << n << "] = " << alpha << " * x + " << beta << " * y "
                    ", x = " << x << ", y = " << y << ", z = " << z )
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
        const ValueType* x = reinterpret_cast<const ValueType*>( xPtr );
        const ValueType* y = reinterpret_cast<const ValueType*>( yPtr );
        ValueType* z = reinterpret_cast<ValueType*>( zPtr );

        #pragma omp parallel for
        for ( int i = 0; i < n; i++ )
        {
            z[i] = *alphaPtr * x[i] + *betaPtr * y[i];
        }
    }
    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICBLAS1::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_DEBUG( logger, "register[" << flag << "], BLAS1<" << common::TypeTraits<ValueType>::id() << ">" )

    KernelRegistry::set<BLASKernelTrait::scal<ValueType> >( MICBLAS1::scal, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::nrm2<ValueType> >( MICBLAS1::nrm2, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::asum<ValueType> >( MICBLAS1::asum, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::iamax<ValueType> >( MICBLAS1::iamax, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::swap<ValueType> >( MICBLAS1::swap, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::copy<ValueType> >( MICBLAS1::copy, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::axpy<ValueType> >( MICBLAS1::axpy, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::dot<ValueType> >( MICBLAS1::dot, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::sum<ValueType> >( MICBLAS1::sum, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static intialization with registration                                   */
/* --------------------------------------------------------------------------- */

MICBLAS1::RegisterGuard::RegisterGuard()
{
    SCAI_LOG_INFO( logger, "register BLAS1 OpenMP-routines for MIC at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_MIC_LIST>::call(
        kregistry::KernelRegistry::KERNEL_ADD );
}

MICBLAS1::RegisterGuard::~RegisterGuard()
{
    SCAI_LOG_INFO( logger, "unregister BLAS1 OpenMP-routines for MIC at kernel registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_MIC_LIST>::call(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

MICBLAS1::RegisterGuard MICBLAS1::guard;    // guard variable for registration

} /* end namespace blaskernel */

} /* end namespace scai */
