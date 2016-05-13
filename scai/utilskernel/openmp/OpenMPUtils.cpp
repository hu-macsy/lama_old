/**
 * @file OpenMPUtils.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Implementation of CSR utilities with OpenMP for the Inteface
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// hpp
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>

// local library
#include <scai/utilskernel/UtilKernelTrait.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tracing.hpp>

#include <scai/common/Constants.hpp>
#include <scai/common/Complex.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/OpenMP.hpp>

namespace scai
{

using common::TypeTraits;

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( OpenMPUtils::logger, "OpenMP.Utils" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::conj( ValueType mValues[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.conj" )

    if ( n > 0 && common::isComplex( TypeTraits<ValueType>::stype ) )
    {
        SCAI_LOG_INFO( logger, "conj, #n = " << n )

        #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType i = 0; i < n; i++ )
        {
            mValues[i] = common::Math::conj( mValues[i] );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType,typename OtherValueType>
void OpenMPUtils::setScale(
    ValueType outValues[],
    const ValueType value,
    const OtherValueType inValues[],
    const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.setScale" )

    SCAI_LOG_INFO( logger, "setScale, #n = " << n << ", value = " << value )

    // alias of outValues == inValues is no problem

    if ( value == common::constants::ZERO )
    {
        // Important : inValues might be undefined
        setVal( outValues, n, value, common::reduction:: COPY );
        return;
    }

    if ( value == common::constants::ONE )
    {
        set( outValues, inValues, n, common::reduction:: COPY );
        return;
    }

    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

    for( IndexType i = 0; i < n; i++ )
    {
        outValues[i] = static_cast<ValueType>( inValues[i] ) * value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::reduceSum( const ValueType array[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.reduceSum" )

    ValueType val( 0 );

    #pragma omp parallel shared( val )
    {
        ValueType threadVal( 0 );

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType i = 0; i < n; ++i )
        {
            threadVal += array[i];
        }

        atomicAdd( val, threadVal );
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::reduceMaxVal( const ValueType array[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.reduceMaxVal" )

    ValueType val( TypeTraits<ValueType>::getMin() );

    #pragma omp parallel
    {
        ValueType threadVal( TypeTraits<ValueType>::getMin() );

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for ( IndexType i = 0; i < n; ++i )
        {
            if ( array[i] > threadVal )
            {
                threadVal = array[i];
            }
        }

        #pragma omp critical
        {
            if ( threadVal > val )
            {
                val = threadVal;
            }
        }
    }
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::reduceMinVal( const ValueType array[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.reduceMinVal" )

    ValueType val( TypeTraits<ValueType>::getMax() );

    #pragma omp parallel
    {
        ValueType threadVal( TypeTraits<ValueType>::getMax() );

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for ( IndexType i = 0; i < n; ++i )
        {
            if ( array[i] < threadVal )
            {
                threadVal = array[i];
            }
        }

        #pragma omp critical
        {
            if( threadVal < val )
            {
                val = threadVal;
            }
        }
    }
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::reduceAbsMaxVal( const ValueType array[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.reduceAbsMaxVal" )

    ValueType val( 0 );

    #pragma omp parallel
    {
        ValueType threadVal( 0 );

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for ( IndexType i = 0; i < n; ++i )
        {
            ValueType elem = common::Math::abs( array[i] );

            if ( elem > threadVal )
            {
                threadVal = elem;
            }
        }

        #pragma omp critical
        {
            if( threadVal > val )
            {
                val = threadVal;
            }
        }
    }
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::reduce( const ValueType array[], const IndexType n, const common::reduction::ReductionOp op )
{
    SCAI_LOG_INFO ( logger, "reduce # array<" << TypeTraits<ValueType>::id() << ">[" << n << "], op = " << op )

    switch ( op )
    {
        case common::reduction::ADD :
            return reduceSum( array, n );
        case common::reduction::MAX :
            return reduceMaxVal( array, n );
        case common::reduction::MIN :
            return reduceMinVal( array, n );
        case common::reduction::ABS_MAX :
            return reduceAbsMaxVal( array, n );
        default:
            COMMON_THROWEXCEPTION( "Unsupported reduce op " << op )
    }

    return ValueType( 0 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::setVal( ValueType array[], const IndexType n, const ValueType val, const common::reduction::ReductionOp op )
{
    SCAI_REGION( "OpenMP.Utils.setVal" )

    SCAI_LOG_INFO( logger, "setVal<" << TypeTraits<ValueType>::id() << ">: " << "array[" << n << "] = "
                            << val << ", op = " << op )

    ValueType value = static_cast<ValueType>( val );
    switch ( op )
    {
        case common::reduction::COPY :
        {
            #pragma omp parallel
            {
                #pragma omp for schedule(SCAI_OMP_SCHEDULE)
                for ( IndexType i = 0; i < n; i++ )
                {
                    array[i] = value;
                }
            }
            break;
        }
        case common::reduction::ADD :
        {
            if ( val == common::constants::ZERO ) 
            {
                return;
            }

            #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)
            for ( IndexType i = 0; i < n; i++ )
            {
                array[i] += value;
            }
            break;
        }
        case common::reduction::SUB :
        {
            if ( val == common::constants::ZERO ) 
            {
                return;
            }

            #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)
            for ( IndexType i = 0; i < n; i++ )
            {
                array[i] -= value;
            }
            break;
        }
        case common::reduction::MULT :
        {
            // scale all values of the array 

            if ( val == common::constants::ONE ) 
            {
                // skip it
            }
            else if ( val == common::constants::ZERO )
            {
                setVal( array, n, ValueType( 0 ), common::reduction::COPY );
            }
            else
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)
    
                for ( IndexType i = 0; i < n; i++ )
                {
                    array[i] *= value;
                }
            }
            break;
        }
        case common::reduction::DIVIDE :
        {
            // scale all values of the array 

            if ( val == common::constants::ONE ) 
            {
                // skip it
            }
            else if ( val == common::constants::ZERO )
            {
                COMMON_THROWEXCEPTION( "DIVIDE by ZEROR, val = " << val )
            }
            else
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)
    
                for ( IndexType i = 0; i < n; i++ )
                {
                    array[i] /= value;
                }
            }
            break;
        }
        default:
            COMMON_THROWEXCEPTION( "Unsupported reduction op : " << op )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::setOrder( ValueType array[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.setOrder" )

    SCAI_LOG_DEBUG( logger,
                    "setOrder<" << TypeTraits<ValueType>::id() << ">: " << "array[" << n << "] = 0, 1, 2, ..., " << ( n - 1 ) )

    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

    for ( IndexType i = 0; i < n; ++i )
    {
        array[i] = static_cast<ValueType>( i );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::getValue( const ValueType* array, const IndexType i )
{
    SCAI_LOG_DEBUG( logger, "getValue<" << TypeTraits<ValueType>::id() << ">: i = " << i )

    return array[i];
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.absMaxDiffVal" )

    SCAI_LOG_DEBUG( logger, "absMaxDiffVal<" << TypeTraits<ValueType>::id() << ">: " << "array[" << n << "]" )

    ValueType val = static_cast<ValueType>(0.0);

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>(0.0);

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for ( IndexType i = 0; i < n; ++i )
        {
            ValueType elem = common::Math::abs( array1[i] - array2[i] );

            if( elem > threadVal )
            {
                threadVal = elem;
            }
        }

        #pragma omp critical
        {
            SCAI_LOG_TRACE( logger, "max val of thread  = " << threadVal << ", global was " << val )

            if( threadVal > val )
            {
                val = threadVal;
            }
        }
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool OpenMPUtils::isSorted( const ValueType array[], const IndexType n, bool ascending )
{
    SCAI_REGION( "OpenMP.Utils.isSorted" )

    SCAI_LOG_INFO( logger,
                   "isSorted<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", ascending = " << ascending )

    bool sorted = true; //!< will be set to false at violations

    if( ascending )
    {
        for ( IndexType i = 1; i < n; i++ )
        {
            if( array[i - 1] > array[i] )
            {
                sorted = false;
                break;
            }
        }
    }
    else
    {
        for ( IndexType i = 1; i < n; i++ )
        {
            if( array[i - 1] < array[i] )
            {
                sorted = false;
                break;
            }
        }
    }

    return sorted;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1,typename ValueType2>
void OpenMPUtils::set( ValueType1 out[], const ValueType2 in[], const IndexType n, const common::reduction::ReductionOp op )
{
    SCAI_REGION( "OpenMP.Utils.set" )

    SCAI_LOG_DEBUG( logger,
                    "set: out<" << TypeTraits<ValueType1>::id() << "[" << n << "]"
                    << ", op = " << op << "  in<" << TypeTraits<ValueType2>::id() << ">[" << n << "]" )

    switch ( op ) 
    {
        case common::reduction::COPY :
        {
            if( in != reinterpret_cast<ValueType2*> ( out ) )
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = static_cast<ValueType1>( in[i] );
                }
            }
            break;
        }
        case common::reduction::ADD :
        {
            #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] += static_cast<ValueType1>( in[i] );
            }
            break;
        }
        case common::reduction::SUB :
        {
            #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] -= static_cast<ValueType1>( in[i] );
            }
            break;
        }
        case common::reduction::DIVIDE :
        {
            #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] /= static_cast<ValueType1>( in[i] );
            }
            break;
        }
        case common::reduction::MULT :
        {
            #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] *= static_cast<ValueType1>( in[i] );
            }
            break;
        }
        default:
        {
            COMMON_THROWEXCEPTION( "unsupported reduction op in set: " << op )
        }
    }
}

/* --------------------------------------------------------------------------- */

bool OpenMPUtils::validIndexes( const IndexType array[], const IndexType n, const IndexType size )
{
    SCAI_REGION( "OpenMP.Utils.validIndexes" )

    SCAI_LOG_DEBUG( logger, "validIndexes: array[" << n << "], size " << size )

    bool validFlag = true;

    // #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE) reduction( & : validFlag )

    for ( IndexType i = 0; i < n; i++ )
    {
        SCAI_LOG_INFO( logger, "validIndexes, array[ " << i << " ] = " << array[i] )

        if( size <= array[i] || 0 > array[i] )
        {
            // exception only in debug mode

            /*
             COMMON_THROWEXCEPTION( "array[" << i << "] = " << array[i]
             << " is illegal index, size = " << size )
             */

            validFlag = false;
        }
    }

    SCAI_LOG_INFO( logger, "validFlag = " << validFlag )

    return validFlag;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1,typename ValueType2>
void OpenMPUtils::setGather( ValueType1 out[], const ValueType2 in[], const IndexType indexes[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.setGather" )

    SCAI_LOG_DEBUG( logger,
                    "setGather: out<" << TypeTraits<ValueType1>::id() << ">[" << n << "]"
                     << " = in<" << TypeTraits<ValueType2>::id() << ">[ indexes[" << n << "] ]" )

    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

    for( IndexType i = 0; i < n; i++ )
    {
        out[i] = static_cast<ValueType1>( in[indexes[i]] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::scatterVal( ValueType out[], const IndexType indexes[], const ValueType value, const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.scatterVal" )

    SCAI_LOG_DEBUG( logger,
                    "scatterVal: out<" << TypeTraits<ValueType>::id() << ">"
                     << "[ indexes[" << n << "] ]" << " = " << value )

    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

    for ( IndexType i = 0; i < n; i++ )
    {
        out[indexes[i]] = value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1,typename ValueType2>
void OpenMPUtils::setScatter( ValueType1 out[], const IndexType indexes[], const ValueType2 in[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.setScatter" )

    SCAI_LOG_DEBUG( logger,
                    "setScatter: out<" << TypeTraits<ValueType1>::id() << ">"
                     << "[ indexes[" << n << "] ]" << " = in<" << TypeTraits<ValueType2>::id() << ">[" << n << "]" )

    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

    for( IndexType i = 0; i < n; i++ )
    {
        out[indexes[i]] = static_cast<ValueType1>( in[i] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::invert( ValueType array[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.invert" )

    SCAI_LOG_INFO( logger, "invert array[ " << n << " ]" )

    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

    for( IndexType i = 0; i < n; ++i )
    {
        array[i] = static_cast<ValueType>(1.0) / array[i];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::scanSerial( ValueType array[], const IndexType n )
{
    SCAI_LOG_DEBUG( logger, "scanSerial: " << n << " entries" )

    // In this case we do it just serial, probably faster

    ValueType runningSum = 0;

    for( IndexType i = 0; i < n; i++ )
    {
        ValueType tmp = runningSum;
        runningSum += array[i];
        SCAI_LOG_TRACE( logger, "scan, row = " << i << ", size = " << array[i] << ", offset = " << runningSum )
        array[i] = tmp;
    }

    return runningSum;;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::scanParallel( PartitionId numThreads, ValueType array[], const IndexType n )
{
    // std::cout << "Scan with " << numThreads << " in parallel" << std::endl;

    // For more threads, we do it in parallel
    // Attention: MUST USE schedule(static)

    common::scoped_array<ValueType> threadCounter( new ValueType[numThreads] );

    SCAI_LOG_DEBUG( logger, "scanParallel: " << n << " entries for " << numThreads << " threads" )

    #pragma omp parallel
    {
        ValueType myCounter = 0;

        #pragma omp for schedule(static)

        for( IndexType i = 0; i < n; i++ )
        {
            myCounter += array[i];
        }

        threadCounter[omp_get_thread_num()] = myCounter;
    }

    ValueType runningSum = scanSerial( threadCounter.get(), numThreads );

    // Each thread sets now its offsets

    #pragma omp parallel
    {
        ValueType myRunningSum = threadCounter[omp_get_thread_num()];

        #pragma omp for schedule(static)

        for( IndexType i = 0; i < n; i++ )
        {
            ValueType tmp = myRunningSum;
            myRunningSum += array[i];
            array[i] = tmp;
        }
    }

    return runningSum;;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::scan( ValueType array[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.scan" )

    SCAI_LOG_INFO( logger, "scan array[ " << n << " ]" )

    int numThreads = 1; // will be set to available threads in parallel region

    #pragma omp parallel
    #pragma omp master
    {
        numThreads = omp_get_num_threads();
    }

    SCAI_LOG_INFO( logger, "scan " << n << " entries, #threads = " << numThreads )

    static int minThreads = 3;  

    ValueType total;

    if ( numThreads < minThreads )
    {
        total = scanSerial( array, n );
    }
    else
    {
        total = scanParallel( numThreads, array, n );
    }

    array[n] = total;

    return total;
}

/* --------------------------------------------------------------------------- */

static void* ptr = NULL;

template<typename ValueType>
void OpenMPUtils::sort( ValueType array[], IndexType perm[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.sort" )

    for ( int i = 0; i < n; ++i )
    {
        perm[i] = i;
    }

    // sort using a custom function object

    struct compare
    {
        static bool f( IndexType a, IndexType b )
        {   
            ValueType* arr = reinterpret_cast<ValueType*>( ptr );
            return arr[a] < arr[b];
        }   
    };

    ptr = array;

    std::sort( perm, perm + n, compare::f );

    common::scoped_array<ValueType> tmp( new ValueType[n] );

    for ( int i = 0; i < n; ++i )
    {
        tmp[i] = array[i];
    }

    for ( int i = 0; i < n; ++i )
    {
        array[i] = tmp[perm[i]];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType OpenMPUtils::countNonZeros( const ValueType denseArray[], const IndexType n, const ValueType eps )
{
    IndexType nonZeros = 0;

    for ( IndexType i = 0; i < n; ++i )
    {
        if ( common::Math::abs( denseArray[i] ) > eps )
        {
            nonZeros++;
        }
    }

    return nonZeros;
}


/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType OpenMPUtils::compress(
        ValueType sparseArray[],
        IndexType sparseIndexes[],
        const ValueType denseArray[],
        const IndexType n,
        const ValueType eps )
{
    IndexType nonZeros = 0;

    for ( IndexType i = 0; i < n; ++i )
    {
        if ( common::Math::abs( denseArray[i] ) > eps )
        {
            sparseArray[nonZeros] = denseArray[i];
            sparseIndexes[nonZeros] = i;
            nonZeros++;
        }
    }

    return nonZeros;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPUtils::Registrator::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::Host;

    SCAI_LOG_INFO( logger, "register UtilsKernel OpenMP-routines for Host at kernel registry [" << flag << "]" )

    // we keep the registrations for IndexType as we do not need conversions

    kregistry::KernelRegistry::set<UtilKernelTrait::validIndexes>( validIndexes, ctx, flag );
}

template<typename ValueType>
void OpenMPUtils::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::Host;

    SCAI_LOG_INFO( logger, "register UtilsKernel OpenMP-routines for Host at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << "]" )

    // we keep the registrations for IndexType as we do not need conversions

    KernelRegistry::set<UtilKernelTrait::conj<ValueType> >( conj, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::reduce<ValueType> >( reduce, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setOrder<ValueType> >( setOrder, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::getValue<ValueType> >( getValue, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setVal<ValueType> >( setVal, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::absMaxDiffVal<ValueType> >( absMaxDiffVal, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::isSorted<ValueType> >( isSorted, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::invert<ValueType> >( invert, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::scan<ValueType> >( scan, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::sort<ValueType> >( sort, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::countNonZeros<ValueType> >( countNonZeros, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::compress<ValueType> >( compress, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void OpenMPUtils::RegistratorVO<ValueType, OtherValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::Host;

    SCAI_LOG_INFO( logger, "register UtilsKernel OpenMP-routines for Host at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )

    KernelRegistry::set<UtilKernelTrait::setScale<ValueType, OtherValueType> >( setScale, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setGather<ValueType, OtherValueType> >( setGather, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setScatter<ValueType, OtherValueType> >( setScatter, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::set<ValueType, OtherValueType> >( set, ctx, flag );

}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPUtils::OpenMPUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;

    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_ARRAY_HOST_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARITHMETIC_ARRAY_HOST_LIST, SCAI_ARITHMETIC_ARRAY_HOST_LIST>::call( flag );

}

OpenMPUtils::~OpenMPUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_ARRAY_HOST_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARITHMETIC_ARRAY_HOST_LIST, SCAI_ARITHMETIC_ARRAY_HOST_LIST>::call( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPUtils OpenMPUtils::guard;   

} /* end namespace utilskernel */

} /* end namespace scai */
