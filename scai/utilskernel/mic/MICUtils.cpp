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
ValueType MICUtils::reduce( const ValueType array[], const IndexType n, const binary::BinaryOp op )
{
    SCAI_LOG_INFO( logger, "reduce # array = " << array << ", n = " << n << ", op = " << op )

    switch ( op )
    {
        case binary::ADD :
            return reduceSum( array, n );

        case binary::MAX :
            return reduceMaxVal( array, n );

        case binary::MIN :
            return reduceMinVal( array, n );

        case binary::ABS_MAX :
            return reduceAbsMaxVal( array, n );

        default:
            COMMON_THROWEXCEPTION( "Unsupported reduce op " << op )
    }

    return ValueType( 0 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::setVal( ValueType array[], const IndexType n, const ValueType val, const binary::BinaryOp op )
{
    SCAI_LOG_DEBUG( logger, "setVal<" << common::getScalarType<ValueType>() << ">: " << "array[" << n << "] = " << val )
    int device = MICContext::getCurrentDevice();
    void* arrayPtr = array;
    const ValueType* valPtr = &val;

    switch ( op )
    {
        case binary::COPY :
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

        case binary::ADD :
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

        case binary::SUB :
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

        case binary::MULT :
        {
            // scale all values of the array
            if ( val == common::constants::ONE )
            {
                // skip it
            }
            else if ( val == common::constants::ZERO )
            {
                setVal( array, n, ValueType( 0 ), binary::COPY );
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

        case binary::DIVIDE :
        {
            // scale all values of the array
            if ( val == common::constants::ONE )
            {
                // skip it
            }
            else if ( val == common::constants::ZERO )
            {
                setVal( array, n, ValueType( 0 ), binary::COPY );
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
            COMMON_THROWEXCEPTION( "Unsupported binary op : " << op )
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

void MICUtils::setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "compute inverse perm, n = " << n )

    void* inversePermPtr = inversePerm;
    const void* permPtr = perm;
    // Parallel execution is safe as perm does not contain a value twice
    int device = MICContext::getCurrentDevice();

    #pragma offload target( mic : device ) in( permPtr, inversePermPtr, n )
    {
        IndexType* inversePerm = reinterpret_cast<IndexType*>( inversePermPtr );
        const IndexType* perm = reinterpret_cast<const IndexType*>( permPtr );

        #pragma omp parallel for
        for ( IndexType ii = 0; ii < n; ii++ )
        {
            IndexType i = perm[ii];
            inversePerm[i] = ii;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
void MICUtils::set( ValueType1 out[], const ValueType2 in[], const IndexType n, const binary::BinaryOp op )
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
            case binary::COPY :
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = static_cast<ValueType1>( in[i] );
                }

                break;
            }

            case binary::ADD :
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] += static_cast<ValueType1>( in[i] );
                }

                break;
            }

            case binary::SUB :
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] -= static_cast<ValueType1>( in[i] );
                }

                break;
            }

            case binary::MULT :
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] *= static_cast<ValueType1>( in[i] );
                }

                break;
            }

            case binary::DIVIDE :
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
void MICUtils::setGather( 
    ValueType1 out[], 
    const ValueType2 in[], 
    const IndexType indexes[], 
    const utilskernel::binary::BinaryOp,
    const IndexType n )
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
void MICUtils::setScatter( ValueType1 out[], const IndexType indexes[], const ValueType2 in[], const binary::BinaryOp op, const IndexType n )
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

        if ( op == binary::COPY )
        {
            #pragma omp parallel for

            for ( IndexType i = 0; i < n; i++ )
            {
                out[indexes[i]] = static_cast<ValueType1>( in[i] );
            }
        }
        else if ( op == binary::ADD )
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
void MICUtils::unaryOp( ValueType out[], const ValueType in[], const IndexType n, const unary::UnaryOp op )
{
    SCAI_LOG_INFO( logger, "unaryOp<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", op = " << op )

    if ( n <= 0 )
    {
        return;
    }

    bool done;  //!< indicates whether op is supported and has been done

    void* out_ptr = out;
    const void* in_ptr = in;

    int device = MICContext::getCurrentDevice();

#pragma offload target( MIC : device ), in ( n, out_ptr, in_ptr, op )
    {
        ValueType* out = static_cast<ValueType*>( out_ptr );
        const ValueType* in  = static_cast<const ValueType*>( in_ptr );

        done = true;

        switch ( op )
        {
            case unary::CONJ :
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = common::Math::conj( in[i] );
                }
                break;
            }

            case unary::MINUS :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = -in[i];
                }
                

                break;
            }

            case unary::ABS :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )
   
                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = common::Math::abs( in[i] );
                }
    
                break;
            }

            case unary::EXP :
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = common::Math::exp( in[i] );
                }
                break;
            }

            case unary::SQRT :
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = common::Math::sqrt( in[i] );
                }
                break;
            }

            case unary::SIN :
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = common::Math::sin( in[i] );
                }
                break;
            }

            case unary::COS :
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = common::Math::cos( in[i] );
                }
                break;
            }

            case unary::TAN :
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = common::Math::tan( in[i] );
                }
                break;
            }

            case unary::ATAN :
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = common::Math::atan( in[i] );
                }
                break;
            }

            case unary::LOG :
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = common::Math::log( in[i] );
                }
                break;
            }

            case unary::FLOOR :
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = common::Math::floor( in[i] );
                }
                break;
            }

            case unary::CEIL :
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = common::Math::ceil( in[i] );
                }
                break;
            }

            default:
            {
                done = false;
            }
        }
    }

    if ( !done )
    {
        COMMON_THROWEXCEPTION( "out[] = op( in[] ), op = " << op << ", unsupported" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::binaryOp( ValueType out[], const ValueType in1[], const ValueType in2[], const IndexType n, const binary::BinaryOp op )
{
    SCAI_LOG_INFO( logger, "applyBiaryOp<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", op = " << op )

    if ( n <= 0 )
    {
        return;
    }

    bool done;   //!< indicates whether op is supported and has been done

    void* out_ptr = out;
    const void* in1_ptr = in1;
    const void* in2_ptr = in2;

    int device = MICContext::getCurrentDevice();

#pragma offload target( MIC : device ), in( n, out_ptr, in1_ptr, in2_ptr, op ), out( done )
    {
        ValueType* out = static_cast<ValueType*>( out_ptr );
        const ValueType* in1  = static_cast<const ValueType*>( in1_ptr );
        const ValueType* in2  = static_cast<const ValueType*>( in2_ptr );

        done = true;

        switch ( op )
        {
            case binary::ADD:
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = in1[i] + in2[i];
                }
                break;
            }

            case binary::SUB:
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = in1[i] - in2[i];
                }
                break;
            }

            case binary::MULT:
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = in1[i] * in2[i];
                }
                break;
            }

            case binary::DIVIDE:
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = in1[i] / in2[i];
                }
                break;
            }

            case binary::POW:
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = common::Math::pow( in1[i], in2[i] );
                }
                break;
            }

            case binary::COPY_SIGN:
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = common::Math::copysign( in1[i], in2[i] );
                }
                break;
            }

            case binary::MIN :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )
    
                for ( IndexType i = 0; i < n; i++ )
                {   
                    out[i] = common::Math::min( in1[i], in2[i] );
                }

                break;
            }

            case binary::MAX :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = common::Math::max( in1[i], in2[i] );
                }

                break;
            }

            case binary::ABS_MAX :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )
    
                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = common::Math::max( common::Math::abs( in1[i] ),
                                                common::Math::abs( in2[i] ) );
                }

                break;
            }

            default:
            {
                done = false;
            }
        }
    }

    if ( !done )
    {
        COMMON_THROWEXCEPTION( "out[] = in1[] <op> in2[], op = " << op << ", unsupported" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::binaryOpScalar1( 
    ValueType out[], 
    const ValueType value, 
    const ValueType in[], 
    const IndexType n, 
    const binary::BinaryOp op )
{
    SCAI_LOG_INFO( logger, "binaryOpScalar1<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", op = " << op )

    if ( n <= 0 )
    {
        return;
    }

    void* out_ptr = out;
    const void* in_ptr = in;
    const ValueType* vPtr = &value;

    bool done;   //!< indicates whether op is supported and has been done

    int device = MICContext::getCurrentDevice();

#pragma offload target( MIC : device ), in( n, out_ptr, in_ptr, vPtr[0:1], op ), out( done )
    {
        ValueType* out = static_cast<ValueType*>( out_ptr );
        const ValueType* in = static_cast<const ValueType*>( in_ptr );
        const ValueType& value = *vPtr; 

        done = true;

        switch ( op )
        {
            case binary::ADD :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = value + in[i];
                }

                break;
            }

            case binary::SUB :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )
    
                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = value - in[i];
                }
    
                break;
            }

            case binary::MULT:
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = value * in[i];
                }
                break;
            }

            case binary::DIVIDE:
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = value / in[i];
                }
                break;
            }

            case binary::POW:
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = common::Math::pow( value, in[i] );
                }
                break;
            }

            case binary::MIN :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = common::Math::min( value, in[i] );
                }
    
                break;
            }

            case binary::MAX :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )
    
                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = common::Math::max( value, in[i] );
                }

                break;
            }

            default:
            {
                done = false;
            }
        }
    }

    if ( !done )
    {
        COMMON_THROWEXCEPTION( "out[] = value <op> in[], op = " << op << ", unsupported" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::binaryOpScalar2( 
    ValueType out[], 
    const ValueType in[], 
    const ValueType value, 
    const IndexType n, 
    const binary::BinaryOp op )
{
    SCAI_LOG_INFO( logger, "binaryOpScalar2<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", op = " << op )

    if ( n <= 0 )
    {
        return;
    }

    bool done;

    void*            out_ptr = out;
    const void*      in_ptr  = in;
    const ValueType* vPtr    = &value;

    int device = MICContext::getCurrentDevice();

#pragma offload target( MIC : device ), in( n, out_ptr, in_ptr, vPtr[0:1], op ), out( done )
    {
        ValueType*       out   = static_cast<ValueType*>( out_ptr );
        const ValueType* in    = static_cast<const ValueType*>( in_ptr );
        const ValueType& value = *vPtr; 

        done = true;

        switch ( op )
        {
            case binary::ADD :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = in[i] + value;
                }

                break;
            }

            case binary::SUB :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = in[i] - value;
                }
    
                break;
            }

            case binary::MULT:
            {
                #pragma omp parallel for

                for ( IndexType i = 0; i < n; ++i )
                {
                    out[i] = in[i] * value;
                }
                break;
            }

            case binary::DIVIDE :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

                for ( IndexType i = 0; i < n; i++ )
                {   
                    out[i] = in[i] / value;
                }

                break;
            }

            case binary::POW :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )
                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = common::Math::pow( in[i], value );
                }

                break;
            }

            case binary::MIN :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )
    
                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = common::Math::min( in[i], value );
                }

                break;
            }

            case binary::MAX :
            {
                #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )
    
                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = common::Math::max( in[i], value );
                }

                break;
            }

            default:
            {
                done = false;
            }
        }
    }

    if ( !done )
    {
        COMMON_THROWEXCEPTION( "out[] = in[] <op> value, op = " << op << ", unsupported" )
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void MICUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::MIC;
    SCAI_LOG_INFO( logger, "register[flag=" << flag << "] UtilsKernel OpenMP-routines for MIC at kernel registry" )
    // we keep the registrations for IndexType as we do not need conversions
    KernelRegistry::set<UtilKernelTrait::validIndexes>( validIndexes, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setInversePerm>( setInversePerm, ctx, flag );
}

template<typename ValueType>
void MICUtils::RegArrayKernels<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
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
    KernelRegistry::set<UtilKernelTrait::scatterVal<ValueType> >( scatterVal, ctx, flag );
}

template<typename ValueType>
void MICUtils::RegNumericKernels<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::MIC;
    SCAI_LOG_DEBUG( logger, "register UtilsKernel OpenMP-routines for MIC at kernel registry [" << flag
                            << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<UtilKernelTrait::unaryOp<ValueType> >( unaryOp, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::binaryOp<ValueType> >( binaryOp, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::binaryOpScalar1<ValueType> >( binaryOpScalar1, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::binaryOpScalar2<ValueType> >( binaryOpScalar2, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void MICUtils::RegistratorVO<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::MIC;
    SCAI_LOG_DEBUG( logger, "register UtilsKernel OpenMP-routines for MIC at kernel registry [" << flag
                             << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )
    // we keep the registrations for IndexType as we do not need conversions
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
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegArrayKernels, SCAI_ARRAY_TYPES_MIC_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegNumericKernels, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARRAY_TYPES_MIC_LIST, SCAI_ARRAY_TYPES_MIC_LIST>::registerKernels( flag );
}

MICUtils::RegisterGuard::~RegisterGuard()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegArrayKernels, SCAI_ARRAY_TYPES_MIC_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegNumericKernels, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARRAY_TYPES_MIC_LIST, SCAI_ARRAY_TYPES_MIC_LIST>::registerKernels( flag );
}

MICUtils::RegisterGuard MICUtils::guard;    // guard variable for registration

} /* end namespace utilskernel */

} /* end namespace scai */
