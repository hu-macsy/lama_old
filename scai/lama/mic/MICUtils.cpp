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
#include <scai/lama/mic/MICUtils.hpp>

// local project
#include <scai/lama/UtilKernelTrait.hpp>

// other SCAI projects
#include <scai/tracing.hpp>

#include <scai/hmemo/mic/MICContext.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/Constants.hpp>

namespace scai
{

using namespace hmemo;

namespace lama
{

SCAI_LOG_DEF_LOGGER( MICUtils::logger, "MIC.Utils" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::scale( ValueType array[], const ValueType value, const IndexType n )
{
    SCAI_LOG_INFO( logger, "scale, #n = " << n << ", value = " << value )

    if( value == scai::common::constants::ONE )
    {
        return;
    }

    if( value == scai::common::constants::ZERO )
    {
        setVal( array, n, static_cast<ValueType>(0.0) );
    }
    else
    {
        void* arrayPtr = array;

        int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( arrayPtr, value, n )
        {
            ValueType* array = static_cast<ValueType*>( arrayPtr );

            #pragma omp parallel for

            for( IndexType i = 0; i < n; i++ )
            {
                array[i] *= value;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType,typename OtherValueType>
void MICUtils::setScale(
    ValueType outValues[],
    const ValueType value,
    const OtherValueType inValues[],
    const IndexType n )
{
    SCAI_LOG_INFO( logger, "setScale, #n = " << n << ", value = " << value )

    // alias of outValues == inValues is no problem

    if( value == scai::common::constants::ZERO )
    {
        // Important : inValues might be undefined
        setVal( outValues, n, value );
        return;
    }

    if( value == scai::common::constants::ONE )
    {
        set( outValues, inValues, n );
        return;
    }

    void* outPtr = outValues;
    const void* inPtr = inValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( outPtr, inPtr, value, n )
    {
        ValueType* outValues = static_cast<ValueType*>( outPtr );
        const OtherValueType* inValues = static_cast<const OtherValueType*>( inPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < n; i++ )
        {
            outValues[i] = static_cast<ValueType>( inValues[i] ) * value;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICUtils::sum( const ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "sum # array = " << array << ", n = " << n )

    ValueType val = static_cast<ValueType>(0.0);

    const void* arrayPtr = array;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( arrayPtr, n ), out( val )
    {
        val = static_cast<ValueType>(0.0);

        const ValueType* array = static_cast<const ValueType*>( arrayPtr );

        #pragma omp parallel for reduction( +:val )

        for( IndexType i = 0; i < n; ++i )
        {
            val += array[i];
        }
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICUtils::setVal( ValueType array[], const IndexType n, const ValueType val, const common::reduction::ReductionOp op )
{
    SCAI_LOG_DEBUG( logger, "setVal<" << common::getScalarType<ValueType>() << ">: " << "array[" << n << "] = " << val )

    int device = MICContext::getCurrentDevice();

    void* arrayPtr = array;

    switch ( op )
    {
        case common::reduction::COPY :
        {
			#pragma offload target( mic : device ), in( arrayPtr, n, val )
			{
				ValueType* array = static_cast<ValueType*>( arrayPtr );

				#pragma omp parallel for

				for( IndexType i = 0; i < n; ++i )
				{
					array[i] = val;
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

			#pragma offload target( mic : device ), in( arrayPtr, n, val )
			{
				ValueType* array = static_cast<ValueType*>( arrayPtr );

				#pragma omp parallel for

				for( IndexType i = 0; i < n; ++i )
				{
					array[i] += val;
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
				#pragma offload target( mic : device ), in( arrayPtr, n, val )
				{
					ValueType* array = static_cast<ValueType*>( arrayPtr );

					#pragma omp parallel for

					for( IndexType i = 0; i < n; ++i )
					{
						array[i] *= val;
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

        for( IndexType i = 0; i < n; ++i )
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

    ValueType val = static_cast<ValueType>(0.0);

    const void* arrayPtr = array;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( arrayPtr, i ), out( val )
    {
        const ValueType* array = static_cast<const ValueType*>( arrayPtr );
        val = array[i];
    }
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICUtils::maxval( const ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "maxval<" << common::getScalarType<ValueType>() << ">: " << "array[" << n << "]" )

    ValueType val = static_cast<ValueType>(0.0);

    if( n < 1 )
    {
        return val;
    }

    const void* arrayPtr = array;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( arrayPtr, n ), out( val )
    {
        val = static_cast<ValueType>(0.0);

        const ValueType* array = static_cast<const ValueType*>( arrayPtr );

        #pragma omp parallel
        {
            ValueType threadVal = static_cast<ValueType>(0.0);

            #pragma omp for

            for( IndexType i = 0; i < n; ++i )
            {
                if( array[i] > threadVal )
                {
                    threadVal = array[i];
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
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICUtils::absMaxVal( const ValueType array[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "absMaxVal<" << common::getScalarType<ValueType>() << ">: " << "array[" << n << "]" )
    SCAI_LOG_INFO( logger, "array = " << array )

    ValueType val = static_cast<ValueType>(0.0);

    // array is already on MIC device

    const void* ptr = array;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( n, ptr ), inout( val )
    {
        const ValueType* array = static_cast<const ValueType*>( ptr );

        #pragma omp parallel
        {
            ValueType threadVal = static_cast<ValueType>(0.0);

            #pragma omp for

            for( IndexType i = 0; i < n; ++i )
            {
                ValueType elem = std::abs( array[i] );

                if( elem > threadVal )
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
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICUtils::absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n )
{
    SCAI_LOG_DEBUG( logger, "absMaxDiffVal<" << common::getScalarType<ValueType>() << ">: " << "array[" << n << "]" )

    ValueType val = static_cast<ValueType>(0.0);

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>(0.0);

        #pragma omp for

        for( IndexType i = 0; i < n; ++i )
        {
            ValueType elem = std::abs( array1[i] - array2[i] );

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

        if( ascending )
        {
            #pragma omp parallel for

            for( IndexType i = 1; i < n; i++ )
            {
                if( sorted && ( array[i - 1] > array[i] ) )
                {
                    sorted = false;
                }
            }
        }
        else
        {
            #pragma omp parallel for

            for( IndexType i = 1; i < n; i++ )
            {
                if( sorted && ( array[i - 1] < array[i] ) )
                {
                    sorted = false;
                }
            }
        }
    }

    return sorted;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1,typename ValueType2>
void MICUtils::set( ValueType1 out[], const ValueType2 in[], const IndexType n )
{
    SCAI_LOG_INFO( logger,
                   "set: out<" << common::getScalarType<ValueType1>() << "[" << n << "]" << " = in<" << common::getScalarType<ValueType2>() << ">[" << n << "]" )

    void* outPtr = out;
    const void* inPtr = in;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( outPtr, inPtr, n )
    {
        ValueType1* out = static_cast<ValueType1*>( outPtr );
        const ValueType2* in = static_cast<const ValueType2*>( inPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < n; i++ )
        {
            out[i] = static_cast<ValueType1>( in[i] );
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

        for( IndexType i = 0; i < n; i++ )
        {
            if( size <= array[i] || 0 > array[i] )
            {
                validFlag = false;
            }
        }
    }

    return validFlag;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1,typename ValueType2>
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

        for( IndexType i = 0; i < n; i++ )
        {
            out[i] = static_cast<ValueType1>( in[indexes[i]] );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1,typename ValueType2>
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

        for( IndexType i = 0; i < n; i++ )
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

        for( IndexType i = 0; i < n; ++i )
        {
            array[i] = static_cast<ValueType>(1.0) / array[i];
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


    KernelRegistry::set<UtilKernelTrait::setScale<float, float> >( setScale, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::setScale<double, float> >( setScale, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::setScale<float, double> >( setScale, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::setScale<double, double> >( setScale, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::sum<IndexType> >( sum, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::sum<float> >( sum, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::sum<double> >( sum, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::setVal<IndexType> >( setVal, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::setVal<float> >( setVal, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::setVal<double> >( setVal, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::setOrder<IndexType> >( setOrder, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::getValue<IndexType> >( getValue, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::getValue<float> >( getValue, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::getValue<double> >( getValue, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::maxval<IndexType> >( maxval, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::maxval<float> >( maxval, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::maxval<double> >( maxval, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::absMaxVal<float> >( absMaxVal, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::absMaxVal<double> >( absMaxVal, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::absMaxDiffVal<float> >( absMaxDiffVal, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::absMaxDiffVal<double> >( absMaxDiffVal, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::isSorted<IndexType> >( isSorted, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::isSorted<float> >( isSorted, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::isSorted<double> >( isSorted, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::setScatter<int, int> >( setScatter, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::setScatter<float, float> >( setScatter, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::setScatter<double, float> >( setScatter, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::setScatter<float, double> >( setScatter, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::setScatter<double, double> >( setScatter, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::setGather<int, int> >( setGather, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::setGather<float, float> >( setGather, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::setGather<double, float> >( setGather, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::setGather<float, double> >( setGather, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::setGather<double, double> >( setGather, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::set<int, int> >( set, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::set<float, float> >( set, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::set<double, float> >( set, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::set<float, double> >( set, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::set<double, double> >( set, MIC, flag );

    KernelRegistry::set<UtilKernelTrait::invert<float> >( invert, MIC, flag );
    KernelRegistry::set<UtilKernelTrait::invert<double> >( invert, MIC, flag );
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

} /* end namespace lama */

} /* end namespace scai */
