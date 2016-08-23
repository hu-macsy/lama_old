/**
 * @file MICUtils.cpp
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
 * @brief Implementation of CSR utilities with MIC for the Inteface
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// hpp
#include <scai/utilskernel/mic/MICUtils.hpp>

// local project
#include <scai/utilskernel/UtilKernelTrait.hpp>

// other SCAI projects
#include <scai/tracing.hpp>

#include <scai/hmemo/mic/MICContext.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/OpenMP.hpp>

namespace scai
{

using namespace hmemo;

using  common::TypeTraits;

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( MICUtils::logger, "MIC.Utils" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::conj( ValueType array[], const IndexType n )
{
    SCAI_LOG_DEBUG( logger,
                    "conj<" << common::getScalarType<ValueType>() << ">: " << "n =" << n )
    void* arrayPtr = array;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device), in ( arrayPtr, n )
    {
        ValueType* array = static_cast<ValueType*>( arrayPtr );
        #pragma omp parallel for

        for ( IndexType i = 0; i < n; ++i )
        {
            array[i] = common::Math::conj( array[i] );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::exp( ValueType array[], const IndexType n )
{
    SCAI_LOG_DEBUG( logger,
                    "exp<" << common::getScalarType<ValueType>() << ">: " << "n =" << n )
    void* arrayPtr = array;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device), in ( arrayPtr, n )
    {
        ValueType* array = static_cast<ValueType*>( arrayPtr );
        #pragma omp parallel for

        for ( IndexType i = 0; i < n; ++i )
        {
            array[i] = common::Math::exp( array[i] );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::vectorScale( ValueType result[], const ValueType x[], const ValueType y[], const IndexType n )
{
    SCAI_LOG_DEBUG( logger,
                    "vectorScale<" << common::getScalarType<ValueType>() << ">: " << "n =" << n )
    void* rPtr = result;
    const void* xPtr = x;
    const void* yPtr = y;
    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device), in ( rPtr, xPtr, yPtr, n )
    {
        ValueType* result = static_cast<ValueType*>( rPtr );
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const ValueType* y = static_cast<const ValueType*>( yPtr );
        #pragma omp parallel for

        for ( IndexType i = 0; i < n; ++i )
        {
            result[i] = x[i] * y[i];
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, typename OtherValueType>
void MICUtils::setScale(
    ValueType outValues[],
    const ValueType value,
    const OtherValueType inValues[],
    const IndexType n )
{
    SCAI_LOG_INFO( logger, "setScale, #n = " << n << ", value = " << value )

    // alias of outValues == inValues is no problem

    if ( value == scai::common::constants::ZERO )
    {
        // Important : inValues might be undefined
        setVal( outValues, n, value, reduction::COPY );
        return;
    }

    if ( value == scai::common::constants::ONE )
    {
        set( outValues, inValues, n, reduction::COPY );
        return;
    }

    void* outPtr = outValues;
    const void* inPtr = inValues;
    const ValueType* valPtr = &value;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ) in( outPtr, inPtr, valPtr[0:1], n )
    {
        ValueType* outValues = static_cast<ValueType*>( outPtr );
        const OtherValueType* inValues = static_cast<const OtherValueType*>( inPtr );
	const ValueType& valRef = *valPtr;

        #pragma omp parallel for

        for ( IndexType i = 0; i < n; i++ )
        {
            outValues[i] = static_cast<ValueType>( inValues[i] ) * valRef;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICUtils::reduceSum( const ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "sum # array = " << array << ", n = " << n )
    ValueType val( 0 );
    const void* arrayPtr = array;
    ValueType* valPtr = &val;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ) in( arrayPtr, n ), out( valPtr[0:1] )
    {
        *valPtr = ValueType( 0 );
        const ValueType* array = reinterpret_cast<const ValueType*>( arrayPtr );
        #pragma omp parallel
        {
            ValueType threadVal = 0;
            #pragma omp for

            for ( IndexType i = 0; i < n; ++i )
            {
                threadVal += array[i];
            }

            #pragma omp critical
            {
                *valPtr += threadVal;
            }
        }
    }
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICUtils::reduce( const ValueType array[], const IndexType n, const reduction::ReductionOp op )
{
    SCAI_LOG_INFO( logger, "reduce # array = " << array << ", n = " << n << ", op = " << op )

    switch ( op )
    {
        case reduction::ADD :
            return reduceSum( array, n );

        case reduction::MAX :
            return reduceMaxVal( array, n );

        case reduction::MIN :
            return reduceMinVal( array, n );

        case reduction::ABS_MAX :
            return reduceAbsMaxVal( array, n );

        default:
            COMMON_THROWEXCEPTION( "Unsupported reduce op " << op )
    }

    return ValueType( 0 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::setVal( ValueType array[], const IndexType n, const ValueType val, const reduction::ReductionOp op )
{
    SCAI_LOG_DEBUG( logger, "setVal<" << common::getScalarType<ValueType>() << ">: " << "array[" << n << "] = " << val )
    int device = MICContext::getCurrentDevice();
    void* arrayPtr = array;
    const ValueType* valPtr = &val;

    switch ( op )
    {
        case reduction::COPY :
#pragma offload target( mic : device ), in( arrayPtr, n, valPtr[0:1] )
            {
                ValueType* array = static_cast<ValueType*>( arrayPtr );
                ValueType value = static_cast<ValueType>( *valPtr );
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    array[i] = value;
                }
            }
            break;

        case reduction::ADD :
        {
            if ( val == common::constants::ZERO )
            {
                // skip it
            }
            else
            {
#pragma offload target( mic : device ), in( arrayPtr, n, valPtr[0:1] )
                {
                    ValueType* array = static_cast<ValueType*>( arrayPtr );
                    ValueType value = static_cast<ValueType>( *valPtr );

                    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                    for ( IndexType i = 0; i < n; i++ )
                    {
                        array[i] += value;
                    }
                }
            }

            break;
        }

        case reduction::SUB :
        {
            if ( val == common::constants::ZERO )
            {
                // skip it
            }
            else
            {
#pragma offload target( mic : device ), in( arrayPtr, n, valPtr[0:1] )
                {
                    ValueType* array = static_cast<ValueType*>( arrayPtr );
                    ValueType value = static_cast<ValueType>( *valPtr );

                    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                    for ( IndexType i = 0; i < n; i++ )
                    {
                        array[i] -= value;
                    }
                }
            }

            break;
        }

        case reduction::MULT :
        {
            // scale all values of the array
            if ( val == common::constants::ONE )
            {
                // skip it
            }
            else if ( val == common::constants::ZERO )
            {
                setVal( array, n, ValueType( 0 ), reduction::COPY );
            }
            else
            {
#pragma offload target( mic : device ), in( arrayPtr, n, valPtr[0:1], op )
                {
                    ValueType* array = static_cast<ValueType*>( arrayPtr );
                    ValueType value = static_cast<ValueType>( *valPtr );

                    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                    for ( IndexType i = 0; i < n; i++ )
                    {
                        array[i] *= value;
                    }
                }
            }

            break;
        }

        case reduction::DIVIDE :
        {
            // scale all values of the array
            if ( val == common::constants::ONE )
            {
                // skip it
            }
            else if ( val == common::constants::ZERO )
            {
                setVal( array, n, ValueType( 0 ), reduction::COPY );
            }
            else
            {
#pragma offload target( mic : device ), in( arrayPtr, n, valPtr[0:1], op )
                {
                    ValueType* array = static_cast<ValueType*>( arrayPtr );
                    ValueType value = static_cast<ValueType>( *valPtr );

                    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                    for ( IndexType i = 0; i < n; i++ )
                    {
                        array[i] /= value;
                    }
                }
            }

            break;
        }

        default:
            COMMON_THROWEXCEPTION( "Unsupported reduction op : " << op )
    }

    SCAI_LOG_DEBUG( logger, "Ready:: setVal<" << common::getScalarType<ValueType>() << ">: " << "array[" << n << "] = " << val )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::setOrder( ValueType array[], const IndexType n )
{
    SCAI_LOG_DEBUG( logger,
                    "setOrder<" << common::getScalarType<ValueType>() << ">: " << "array[" << n << "] = 0, 1, 2, ..., " << ( n - 1 ) )
    void* arrayPtr = array;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device), in ( arrayPtr, n )
    {
        ValueType* array = static_cast<ValueType*>( arrayPtr );
        #pragma omp parallel for

        for ( IndexType i = 0; i < n; ++i )
        {
            array[i] = static_cast<ValueType>( i );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::setSequence( ValueType array[], const ValueType startValue, const ValueType inc, const IndexType n )
{
    SCAI_LOG_DEBUG( logger, "setSequence<" << TypeTraits<ValueType>::id() << ">: " << "array[" << n << "] = " << startValue
                            << "..., " << ( startValue + (n-1) * inc ) )
    
    void* arrayPtr = array;
    const ValueType* sPtr = &startValue;
    const ValueType* iPtr = &inc;
    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device), in ( arrayPtr, n, sPtr[0:1], iPtr[0:1] )
    {
        ValueType* array = static_cast<ValueType*>( arrayPtr );
	const ValueType& startValue = *sPtr;
	const ValueType& inc = *iPtr;
        #pragma omp parallel for
        for ( IndexType i = 0; i < n; ++i )
        {
            array[i] = static_cast<ValueType>( startValue + i * inc );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICUtils::getValue( const ValueType* array, const IndexType i )
{
    SCAI_LOG_DEBUG( logger, "getValue<" << common::getScalarType<ValueType>() << ">: i = " << i )
    ValueType val = static_cast<ValueType>( 0.0 );
    const void* arrayPtr = array;
    ValueType* valPtr = &val;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ), in( arrayPtr, i ), out( valPtr[0:1] )
    {
        const ValueType* array = static_cast<const ValueType*>( arrayPtr );
        *valPtr = array[i];
    }
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICUtils::reduceMaxVal( const ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "maxval<" << common::getScalarType<ValueType>() << ">: " << "array[" << n << "]" )
    const ValueType zero( common::TypeTraits<ValueType>::getMin() );
    ValueType val = zero;

    if ( n < 1 )
    {
        return val;
    }

    const void* arrayPtr = array;
    ValueType* valPtr = &val;
    const ValueType* zeroPtr = &zero;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ) in( arrayPtr, n, zeroPtr[0:1] ), out( valPtr[0:1] )
    {
        *valPtr = *zeroPtr;
        const ValueType* array = static_cast<const ValueType*>( arrayPtr );
        #pragma omp parallel
        {
            ValueType threadVal( *zeroPtr );
            #pragma omp for

            for ( IndexType i = 0; i < n; ++i )
            {
                if ( array[i] > threadVal )
                {
                    threadVal = array[i];
                }
            }

            #pragma omp critical
            {
                if ( threadVal > ( *valPtr ) )
                {
                    *valPtr = threadVal;
                }
            }
        }
    }
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICUtils::reduceMinVal( const ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "maxval<" << common::getScalarType<ValueType>() << ">: " << "array[" << n << "]" )
    ValueType zero( common::TypeTraits<ValueType>::getMax() );
    ValueType val = zero;

    if ( n < 1 )
    {
        return val;
    }

    const void* arrayPtr = array;
    ValueType* zeroPtr = &zero;
    ValueType* valPtr = &val;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ) in( arrayPtr, n, zeroPtr[0:1] ), out( valPtr[0:1] )
    {
        *valPtr = *zeroPtr;
        const ValueType* array = static_cast<const ValueType*>( arrayPtr );
        #pragma omp parallel
        {
            ValueType threadVal( *zeroPtr );
            #pragma omp for

            for ( IndexType i = 0; i < n; ++i )
            {
                if ( array[i] < threadVal )
                {
                    threadVal = array[i];
                }
            }

            #pragma omp critical
            {
                if ( threadVal < *valPtr )
                {
                    *valPtr = threadVal;
                }
            }
        }
    }
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICUtils::reduceAbsMaxVal( const ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "absMaxVal<" << common::getScalarType<ValueType>() << ">: " << "array[" << n << "]" )
    SCAI_LOG_INFO( logger, "array = " << array )
    ValueType val( 0 );
    // array is already on MIC device
    const void* ptr = array;
    ValueType* valPtr = &val;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ) in( n, ptr ), inout( valPtr[0:1] )
    {
        const ValueType* array = reinterpret_cast<const ValueType*>( ptr );
        #pragma omp parallel
        {
            ValueType threadVal( 0 );
            #pragma omp for

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
                if ( threadVal > *valPtr )
                {
                    *valPtr = threadVal;
                }
            }
        }
    }
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICUtils::absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n )
{
    SCAI_LOG_DEBUG( logger, "absMaxDiffVal<" << common::getScalarType<ValueType>() << ">: " << "array[" << n << "]" )
    ValueType val( 0 );
    const void* ptr1 = array1;
    const void* ptr2 = array2;
    ValueType* valPtr = &val;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ) in( n, ptr1, ptr2 ), inout( valPtr[0:1] )
    {
        const ValueType* arrayPtr1 = reinterpret_cast<const ValueType*>( ptr1 );
        const ValueType* arrayPtr2 = reinterpret_cast<const ValueType*>( ptr2 );
        #pragma omp parallel
        {
            ValueType threadVal( 0 );
            #pragma omp for

            for ( IndexType i = 0; i < n; ++i )
            {
                ValueType elem = common::Math::abs( arrayPtr1[i] - arrayPtr2[i] );

                if ( elem > threadVal )
                {
                    threadVal = elem;
                }
            }

            #pragma omp critical
            {
                if ( threadVal > *valPtr )
                {
                    *valPtr = threadVal;
                }
            }
        }
    }
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool MICUtils::isSorted( const ValueType array[], const IndexType n, bool ascending )
{
    SCAI_LOG_INFO( logger,
                   "isSorted<" << common::getScalarType<ValueType>() << ">, n = " << n << ", ascending = " << ascending )
    bool sorted = true; //!< will be set to false at violations
    const void* arrayPtr = array;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ) in( arrayPtr, n, ascending ), out( sorted )
    {
        const ValueType* array = static_cast<const ValueType*>( arrayPtr );
        sorted = true;

        if ( ascending )
        {
            #pragma omp parallel for

            for ( IndexType i = 1; i < n; i++ )
            {
                if ( sorted && ( array[i - 1] > array[i] ) )
                {
                    sorted = false;
                }
            }
        }
        else
        {
            #pragma omp parallel for

            for ( IndexType i = 1; i < n; i++ )
            {
                if ( sorted && ( array[i - 1] < array[i] ) )
                {
                    sorted = false;
                }
            }
        }
    }
    return sorted;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
void MICUtils::set( ValueType1 out[], const ValueType2 in[], const IndexType n, const reduction::ReductionOp op )
{
    SCAI_LOG_DEBUG( logger,
                    "set: out<" << TypeTraits<ValueType1>::id() << "[" << n << "]"
                    << ", op = " << op << "  in<" << TypeTraits<ValueType2>::id() << ">[" << n << "]" )
    void* outPtr = out;
    const void* inPtr = in;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ), in( outPtr, inPtr, n, op )
    {
        ValueType1* out = static_cast<ValueType1*>( outPtr );
        const ValueType2* in = static_cast<const ValueType2*>( inPtr );

        switch ( op )
        {
            case reduction::COPY :
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = static_cast<ValueType1>( in[i] );
                }

                break;
            }

            case reduction::ADD :
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] += static_cast<ValueType1>( in[i] );
                }

                break;
            }

            case reduction::SUB :
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] -= static_cast<ValueType1>( in[i] );
                }

                break;
            }

            case reduction::MULT :
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] *= static_cast<ValueType1>( in[i] );
                }

                break;
            }

            case reduction::DIVIDE :
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] /= static_cast<ValueType1>( in[i] );
                }

                break;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

bool MICUtils::validIndexes( const IndexType array[], const IndexType n, const IndexType size )
{
    SCAI_LOG_DEBUG( logger, "validIndexes: array[" << n << "], size " << size )
    bool validFlag = true;
    const void* arrayPtr = array;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ), in ( n, arrayPtr, size ), inout( validFlag )
    {
        const IndexType* array = static_cast<const IndexType*>( arrayPtr );
        #pragma omp parallel for reduction( & : validFlag )

        for ( IndexType i = 0; i < n; i++ )
        {
            if ( size <= array[i] || 0 > array[i] )
            {
                validFlag = false;
            }
        }
    }
    return validFlag;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
void MICUtils::setGather( ValueType1 out[], const ValueType2 in[], const IndexType indexes[], const IndexType n )
{
    SCAI_LOG_DEBUG( logger,
                    "setGather: out<" << common::getScalarType<ValueType1>() << ">[" << n << "]" << " = in<" << common::getScalarType<ValueType2>() << ">[ indexes[" << n << "] ]" )
    void* outPtr = out;
    const void* inPtr = in;
    const void* indexesPtr = indexes;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ) in( inPtr, outPtr, indexesPtr, n )
    {
        ValueType1* out = static_cast<ValueType1*>( outPtr );
        const ValueType2* in = static_cast<const ValueType2*>( inPtr );
        const IndexType* indexes = static_cast<const IndexType*>( indexesPtr );
        #pragma omp parallel for

        for ( IndexType i = 0; i < n; i++ )
        {
            out[i] = static_cast<ValueType1>( in[indexes[i]] );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
void MICUtils::setScatter( ValueType1 out[], const IndexType indexes[], const ValueType2 in[], const reduction::ReductionOp op, const IndexType n )
{
    SCAI_LOG_DEBUG( logger,
                    "setScatter: out<" << common::getScalarType<ValueType1>() << ">" << "[ indexes[" << n << "] ]" << " = in<" << common::getScalarType<ValueType2>() << ">[" << n << "]" )
    void* outPtr = out;
    const void* indexesPtr = indexes;
    const void* inPtr = in;
    int device = MICContext::getCurrentDevice();

    #pragma offload target( mic : device ) in( outPtr, indexesPtr, inPtr, n, op )
    {
        ValueType1* out = static_cast<ValueType1*>( outPtr );
        const ValueType2* in = static_cast<const ValueType2*>( inPtr );
        const IndexType* indexes = static_cast<const IndexType*>( indexesPtr );

        if ( op == reduction::COPY )
        {
            #pragma omp parallel for

            for ( IndexType i = 0; i < n; i++ )
            {
                out[indexes[i]] = static_cast<ValueType1>( in[i] );
            }
        }
        else if ( op == reduction::ADD )
        {
            // No parallelization at this time

            for ( IndexType i = 0; i < n; i++ )
            {
                out[indexes[i]] += static_cast<ValueType1>( in[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::scatterVal( ValueType out[], const IndexType indexes[], const ValueType value, const IndexType n )
{
    SCAI_LOG_DEBUG( logger,
                    "scatterVal: out<" << TypeTraits<ValueType>::id() << ">"
                    << "[ indexes[" << n << "] ]" << " = " << value )
    void* outPtr = out;
    const void* indexesPtr = indexes;
    const ValueType* vPtr = &value;
    int device = MICContext::getCurrentDevice();

    #pragma offload target( mic : device ) in( outPtr, indexesPtr, n, vPtr[0:1] )
    {
        ValueType* out = static_cast<ValueType*>( outPtr );
        const IndexType* indexes = static_cast<const IndexType*>( indexesPtr );
	const ValueType& value = *vPtr; 
        #pragma omp parallel for

        for ( IndexType i = 0; i < n; i++ )
        {
            out[indexes[i]] = value;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::invert( ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "invert array[ " << n << " ]" )
    void* array_ptr = array;
    int device = MICContext::getCurrentDevice();
#pragma offload target( MIC : device ), in ( n, array_ptr )
    {
        ValueType* array = static_cast<ValueType*>( array_ptr );
        #pragma omp parallel for

        for ( IndexType i = 0; i < n; ++i )
        {
            array[i] = static_cast<ValueType>( 1.0 ) / array[i];
        }
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void MICUtils::Registrator::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::MIC;
    SCAI_LOG_INFO( logger, "register[flag=" << flag << "] UtilsKernel OpenMP-routines for MIC at kernel registry" )
    // we keep the registrations for IndexType as we do not need conversions
    KernelRegistry::set<UtilKernelTrait::validIndexes>( validIndexes, ctx, flag );
}

template<typename ValueType>
void MICUtils::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::MIC;
    SCAI_LOG_DEBUG( logger, "register UtilsKernel OpenMP-routines for MIC at kernel registry [" << flag
                            << " --> " << common::getScalarType<ValueType>() << "]" )
    // we keep the registrations for IndexType as we do not need conversions
    KernelRegistry::set<UtilKernelTrait::reduce<ValueType> >( reduce, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setOrder<ValueType> >( setOrder, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setSequence<ValueType> >( setSequence, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::getValue<ValueType> >( getValue, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::absMaxDiffVal<ValueType> >( absMaxDiffVal, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::isSorted<ValueType> >( isSorted, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setVal<ValueType> >( setVal, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::invert<ValueType> >( invert, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::exp<ValueType> >( exp, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::conj<ValueType> >( conj, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::vectorScale<ValueType> >( vectorScale, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::scatterVal<ValueType> >( scatterVal, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void MICUtils::RegistratorVO<ValueType, OtherValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::MIC;
    SCAI_LOG_DEBUG( logger, "register UtilsKernel OpenMP-routines for MIC at kernel registry [" << flag
                             << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )
    // we keep the registrations for IndexType as we do not need conversions
    KernelRegistry::set<UtilKernelTrait::setScale<ValueType, OtherValueType> >( setScale, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setGather<ValueType, OtherValueType> >( setGather, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setScatter<ValueType, OtherValueType> >( setScatter, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::set<ValueType, OtherValueType> >( set, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICUtils::RegisterGuard::RegisterGuard()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_ARRAY_MIC_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARITHMETIC_ARRAY_MIC_LIST, SCAI_ARITHMETIC_ARRAY_MIC_LIST>::call( flag );
}

MICUtils::RegisterGuard::~RegisterGuard()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_ARRAY_MIC_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARITHMETIC_ARRAY_MIC_LIST, SCAI_ARITHMETIC_ARRAY_MIC_LIST>::call( flag );
}

MICUtils::RegisterGuard MICUtils::guard;    // guard variable for registration

} /* end namespace utilskernel */

} /* end namespace scai */
