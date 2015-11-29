/**
 * @file OpenMPUtils.cpp
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
 * @brief Implementation of CSR utilities with OpenMP for the Inteface
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/lama/openmp/OpenMPUtils.hpp>

// local library
#include <scai/lama/UtilKernelTrait.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tracing.hpp>

#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/OpenMP.hpp>

// boost
#include <boost/preprocessor.hpp>

using std::abs;

namespace scai
{

using common::TypeTraits;

namespace lama
{

SCAI_LOG_DEF_LOGGER( OpenMPUtils::logger, "OpenMP.Utils" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::scale( ValueType mValues[], const ValueType value, const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.scale" )

    SCAI_LOG_INFO( logger, "scale, #n = " << n << ", value = " << value )

    if( value == scai::common::constants::ONE )
    {
        return;
    }

    if( value == scai::common::constants::ZERO )
    {
        #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType i = 0; i < n; i++ )
        {
            mValues[i] = static_cast<ValueType>(0.0);
        }
    }
    else
    {
        #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType i = 0; i < n; i++ )
        {
            mValues[i] *= value;
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

    if ( value == scai::common::constants::ZERO )
    {
        // Important : inValues might be undefined
        setVal( outValues, n, value );
        return;
    }

    if ( value == scai::common::constants::ONE )
    {
        set( outValues, inValues, n );
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
ValueType OpenMPUtils::sum( const ValueType array[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.sum" )

    SCAI_LOG_INFO( logger, "sum # array = " << array << ", n = " << n )
    ValueType val = static_cast<ValueType>(0.0);

    #pragma omp parallel shared( val )
    {
        ValueType tVal = static_cast<ValueType>(0.0);

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType i = 0; i < n; ++i )
        {
            tVal += array[i];
        }

        atomicAdd( val, tVal );
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::setVal( ValueType array[], const IndexType n, const ValueType val )
{
    SCAI_REGION( "OpenMP.Utils.setVal" )

    SCAI_LOG_DEBUG( logger, "setVal<" << TypeTraits<ValueType>::id() << ">: " << "array[" << n << "] = " << val )

    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

    for( IndexType i = 0; i < n; ++i )
    {
        array[i] = val;
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

    for( IndexType i = 0; i < n; ++i )
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
ValueType OpenMPUtils::maxval( const ValueType array[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.maxVal" )

    SCAI_LOG_DEBUG( logger, "maxval<" << TypeTraits<ValueType>::id() << ">: " << "array[" << n << "]" )

    ValueType val = static_cast<ValueType>(0.0);

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>(0.0);

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType i = 0; i < n; ++i )
        {
            if( array[i] > threadVal )
            {
                threadVal = array[i];
            }
        }

        #pragma omp critical
        {
            SCAI_LOG_TRACE( logger, "max val of thread = " << threadVal << ", global was " << val )

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
ValueType OpenMPUtils::absMaxVal( const ValueType array[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.absMaxVal" )

    SCAI_LOG_DEBUG( logger, "absMaxVal<" << TypeTraits<ValueType>::id() << ">: " << "array[" << n << "]" )

    ValueType val = static_cast<ValueType>(0.0);

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>(0.0);

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType i = 0; i < n; ++i )
        {
            ValueType elem = abs( array[i] );

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
ValueType OpenMPUtils::absMaxDiffVal( const ValueType array1[], const ValueType array2[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.absMaxDiffVal" )

    SCAI_LOG_DEBUG( logger, "absMaxDiffVal<" << TypeTraits<ValueType>::id() << ">: " << "array[" << n << "]" )

    ValueType val = static_cast<ValueType>(0.0);

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>(0.0);

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType i = 0; i < n; ++i )
        {
            ValueType elem = abs( array1[i] - array2[i] );

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
        for( IndexType i = 1; i < n; i++ )
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
        for( IndexType i = 1; i < n; i++ )
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
void OpenMPUtils::set( ValueType1 out[], const ValueType2 in[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.set" )

    SCAI_LOG_DEBUG( logger,
                    "set: out<" << TypeTraits<ValueType1>::id() << "[" << n << "]" << " = in<" << TypeTraits<ValueType2>::id() << ">[" << n << "]" )

    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

    for( IndexType i = 0; i < n; i++ )
    {
        out[i] = static_cast<ValueType1>( in[i] );
    }
}

/* --------------------------------------------------------------------------- */

bool OpenMPUtils::validIndexes( const IndexType array[], const IndexType n, const IndexType size )
{
    SCAI_REGION( "OpenMP.Utils.validIndexes" )

    SCAI_LOG_DEBUG( logger, "validIndexes: array[" << n << "], size " << size )

    bool validFlag = true;

    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE) reduction( & : validFlag )

    for( IndexType i = 0; i < n; i++ )
    {
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
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPUtils::registerKernels( bool deleteFlag )
{
    using namespace scai::kregistry;
    using common::context::Host;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD ;   // add it or delete it
 
    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

    // Instantations for IndexType, not done by ARITHMETIC_TYPE macrods

    KernelRegistry::set<UtilKernelTrait::validIndexes>( validIndexes, Host, flag );

    // we keep the registrations for IndexType as we do not need conversions

    KernelRegistry::set<UtilKernelTrait::sum<IndexType> >( sum, Host, flag );

    KernelRegistry::set<UtilKernelTrait::setVal<IndexType> >( setVal, Host, flag );
    KernelRegistry::set<UtilKernelTrait::setOrder<IndexType> >( setOrder, Host, flag );
    KernelRegistry::set<UtilKernelTrait::getValue<IndexType> >( getValue, Host, flag );

    KernelRegistry::set<UtilKernelTrait::maxval<IndexType> >( maxval, Host, flag );
    KernelRegistry::set<UtilKernelTrait::isSorted<IndexType> >( isSorted, Host, flag );

    KernelRegistry::set<UtilKernelTrait::setScatter<IndexType, IndexType> >( setScatter, Host, flag );
    KernelRegistry::set<UtilKernelTrait::setGather<IndexType, IndexType> >( setGather, Host, flag );
    KernelRegistry::set<UtilKernelTrait::set<IndexType, IndexType> >( set, Host, flag );

#define LAMA_UTILS2_REGISTER(z, J, TYPE )                                                                        \
    KernelRegistry::set<UtilKernelTrait::setScale<TYPE, ARITHMETIC_HOST_TYPE_##J> >( setScale, Host, flag );     \
    KernelRegistry::set<UtilKernelTrait::setGather<TYPE, ARITHMETIC_HOST_TYPE_##J> >( setGather, Host, flag );   \
    KernelRegistry::set<UtilKernelTrait::setScatter<TYPE, ARITHMETIC_HOST_TYPE_##J> >( setScatter, Host, flag ); \
    KernelRegistry::set<UtilKernelTrait::set<TYPE, ARITHMETIC_HOST_TYPE_##J> >( set, Host, flag );               \

#define LAMA_UTILS_REGISTER(z, I, _)                                                                             \
    KernelRegistry::set<UtilKernelTrait::scale<ARITHMETIC_HOST_TYPE_##I> >( scale, Host, flag );                 \
    KernelRegistry::set<UtilKernelTrait::sum<ARITHMETIC_HOST_TYPE_##I> >( sum, Host, flag );                     \
    KernelRegistry::set<UtilKernelTrait::setVal<ARITHMETIC_HOST_TYPE_##I> >( setVal, Host, flag );               \
    KernelRegistry::set<UtilKernelTrait::setOrder<ARITHMETIC_HOST_TYPE_##I> >( setOrder, Host, flag );           \
    KernelRegistry::set<UtilKernelTrait::getValue<ARITHMETIC_HOST_TYPE_##I> >( getValue, Host, flag );           \
    KernelRegistry::set<UtilKernelTrait::maxval<ARITHMETIC_HOST_TYPE_##I> >( maxval, Host, flag );               \
    KernelRegistry::set<UtilKernelTrait::absMaxVal<ARITHMETIC_HOST_TYPE_##I> >( absMaxVal, Host, flag );         \
    KernelRegistry::set<UtilKernelTrait::absMaxDiffVal<ARITHMETIC_HOST_TYPE_##I> >( absMaxDiffVal, Host, flag ); \
    KernelRegistry::set<UtilKernelTrait::isSorted<ARITHMETIC_HOST_TYPE_##I> >( isSorted, Host, flag );           \
    KernelRegistry::set<UtilKernelTrait::invert<ARITHMETIC_HOST_TYPE_##I> >( invert, Host, flag );               \
    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT,                                                                   \
                     LAMA_UTILS2_REGISTER,                                                                       \
                     ARITHMETIC_HOST_TYPE_##I )

    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_UTILS_REGISTER, _ )

#undef LAMA_UTILS_REGISTER
#undef LAMA_UTILS2_REGISTER

}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPUtils::OpenMPUtils()
{
    bool deleteFlag = false;  
    registerKernels( deleteFlag );
}

OpenMPUtils::~OpenMPUtils()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPUtils OpenMPUtils::guard;   

} /* end namespace lama */

} /* end namespace scai */
