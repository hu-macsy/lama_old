/**
 * @file MICUtils.cpp
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
 * @brief Implementation of CSR utilities with MIC for the Inteface
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.1.0
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
        setVal( outValues, n, value, common::reduction::COPY );
        return;
    }

    if ( value == scai::common::constants::ONE )
    {
        set( outValues, inValues, n, common::reduction::COPY );
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

        #pragma omp parallel for

        for ( IndexType i = 0; i < n; i++ )
        {
            outValues[i] = static_cast<ValueType>( inValues[i] ) * (*valPtr);
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

		ValueType threadVal = 0;
        #pragma omp parallel private(threadVal)
        {
			#pragma omp parallel for
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
ValueType MICUtils::reduce( const ValueType array[], const IndexType n, const common::reduction::ReductionOp op )
{
    SCAI_LOG_INFO( logger, "reduce # array = " << array << ", n = " << n << ", op = " << op )

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
void MICUtils::setVal( ValueType array[], const IndexType n, const ValueType val, const common::reduction::ReductionOp op )
{
    SCAI_LOG_DEBUG( logger, "setVal<" << common::getScalarType<ValueType>() << ">: " << "array[" << n << "] = " << val )

    int device = MICContext::getCurrentDevice();

    void* arrayPtr = array;
    const ValueType* valPtr = &val;

    switch ( op )
    {
        case common::reduction::COPY :

            #pragma offload target( mic : device ), in( arrayPtr, n, valPtr[0:1] )
            {
                ValueType* array = static_cast<ValueType*>( arrayPtr );

                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    array[i] = *valPtr;
                }
            }

            break;

        case common::reduction::ADD :
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

                    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                    for ( IndexType i = 0; i < n; i++ )
                    {
                        array[i] += *valPtr;
                    }
                }
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
                #pragma offload target( mic : device ), in( arrayPtr, n, valPtr[0:1], op )
                {
                    ValueType* array = static_cast<ValueType*>( arrayPtr );

                    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                    for ( IndexType i = 0; i < n; i++ )
                    {
                        array[i] *= *valPtr;
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

    const ValueType zero( - common::TypeTraits<ValueType>::getMax() );

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
                if ( threadVal > (*valPtr) )
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
    ValueType *valPtr = &val;

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

    #pragma omp parallel
    {
        ValueType threadVal( 0 );

        #pragma omp for

        for ( IndexType i = 0; i < n; ++i )
        {
            ValueType elem = common::Math::abs( array1[i] - array2[i] );

            if ( elem > threadVal )
            {
                threadVal = elem;
            }
        }

        #pragma omp critical
        {
            SCAI_LOG_TRACE( logger, "max val of thread  = " << threadVal << ", global was " << val )

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
void MICUtils::set( ValueType1 out[], const ValueType2 in[], const IndexType n, const common::reduction::ReductionOp op )
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
            case common::reduction::COPY :
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = static_cast<ValueType1>( in[i] );
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
            case common::reduction::MULT :
            {
                #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] *= static_cast<ValueType1>( in[i] );
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
void MICUtils::setScatter( ValueType1 out[], const IndexType indexes[], const ValueType2 in[], const IndexType n )
{
    SCAI_LOG_DEBUG( logger,
                    "setScatter: out<" << common::getScalarType<ValueType1>() << ">" << "[ indexes[" << n << "] ]" << " = in<" << common::getScalarType<ValueType2>() << ">[" << n << "]" )

    void* outPtr = out;
    const void* indexesPtr = indexes;
    const void* inPtr = in;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( outPtr, indexesPtr, inPtr, n )
    {
        ValueType1* out = static_cast<ValueType1*>( outPtr );
        const ValueType2* in = static_cast<const ValueType2*>( inPtr );
        const IndexType* indexes = static_cast<const IndexType*>( indexesPtr );

        #pragma omp parallel for

        for ( IndexType i = 0; i < n; i++ )
        {
            out[indexes[i]] = static_cast<ValueType1>( in[i] );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::invert( ValueType array[], const IndexType n )
{
    // SCAI_REGION( "MIC.invert" )

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

void MICUtils::registerKernels( bool deleteFlag )
{
    SCAI_LOG_INFO( logger, "register Utils kernels for MIC in Kernel Registry" )

    using kregistry::KernelRegistry;
    using common::context::MIC;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD ;   // add it or delete it

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

    // Instantations for IndexType, not done by ARITHMETIC_TYPE macrods

    KernelRegistry::set<UtilKernelTrait::validIndexes>( validIndexes, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::reduce<IndexType> >( reduce, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::setVal<IndexType> >( setVal, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::setOrder<IndexType> >( setOrder, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::getValue<IndexType> >( getValue, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::isSorted<IndexType> >( isSorted, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::setScatter<int, int> >( setScatter, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::setGather<int, int> >( setGather, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::set<int, int> >( set, MIC, flag );


#define LAMA_UTILS2_REGISTER(z, J, TYPE )                                                                        \
    KernelRegistry::set<UtilKernelTrait::setScale<TYPE, ARITHMETIC_MIC_TYPE_##J> >( setScale, MIC, flag );     \
    KernelRegistry::set<UtilKernelTrait::setGather<TYPE, ARITHMETIC_MIC_TYPE_##J> >( setGather, MIC, flag );   \
    KernelRegistry::set<UtilKernelTrait::setScatter<TYPE, ARITHMETIC_MIC_TYPE_##J> >( setScatter, MIC, flag ); \
    KernelRegistry::set<UtilKernelTrait::set<TYPE, ARITHMETIC_MIC_TYPE_##J> >( set, MIC, flag );               \

#define LAMA_UTILS_REGISTER(z, I, _)                                                                             \
    KernelRegistry::set<UtilKernelTrait::reduce<ARITHMETIC_MIC_TYPE_##I> >( reduce, MIC, flag );               \
    KernelRegistry::set<UtilKernelTrait::setVal<ARITHMETIC_MIC_TYPE_##I> >( setVal, MIC, flag );               \
    KernelRegistry::set<UtilKernelTrait::setOrder<ARITHMETIC_MIC_TYPE_##I> >( setOrder, MIC, flag );           \
    KernelRegistry::set<UtilKernelTrait::getValue<ARITHMETIC_MIC_TYPE_##I> >( getValue, MIC, flag );           \
    KernelRegistry::set<UtilKernelTrait::absMaxDiffVal<ARITHMETIC_MIC_TYPE_##I> >( absMaxDiffVal, MIC, flag ); \
    KernelRegistry::set<UtilKernelTrait::isSorted<ARITHMETIC_MIC_TYPE_##I> >( isSorted, MIC, flag );           \
    KernelRegistry::set<UtilKernelTrait::invert<ARITHMETIC_MIC_TYPE_##I> >( invert, MIC, flag );               \
    BOOST_PP_REPEAT( ARITHMETIC_MIC_TYPE_CNT,                                                                   \
                     LAMA_UTILS2_REGISTER,                                                                       \
                     ARITHMETIC_MIC_TYPE_##I )

    BOOST_PP_REPEAT( ARITHMETIC_MIC_TYPE_CNT, LAMA_UTILS_REGISTER, _ )

#undef LAMA_UTILS_REGISTER
#undef LAMA_UTILS2_REGISTER
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICUtils::RegisterGuard::RegisterGuard()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

MICUtils::RegisterGuard::~RegisterGuard()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

MICUtils::RegisterGuard MICUtils::guard;    // guard variable for registration

} /* end namespace utilskernel */

} /* end namespace scai */