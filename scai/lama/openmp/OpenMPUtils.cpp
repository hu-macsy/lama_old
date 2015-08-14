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
#include <scai/lama/openmp/OpenMP.hpp>

// others
#include <scai/lama/LAMAInterfaceRegistry.hpp>
#include <scai/tracing.hpp>

// boost
#include <boost/preprocessor.hpp>

namespace lama
{

using std::abs;
using common::getScalarType;

LAMA_LOG_DEF_LOGGER( OpenMPUtils::logger, "OpenMP.Utils" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::scale( ValueType mValues[], const ValueType value, const IndexType n )
{
    LAMA_REGION( "OpenMP.Utils.scale" )

    LAMA_LOG_INFO( logger, "scale, #n = " << n << ", value = " << value )

    if( value == static_cast<ValueType>( 1 ) )
    {
        return;
    }

    if( value == 0 )
    {
        #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )

        for( IndexType i = 0; i < n; i++ )
        {
            mValues[i] = 0;
        }
    }
    else
    {
        #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )

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
    LAMA_REGION( "OpenMP.Utils.setScale" )

    LAMA_LOG_INFO( logger, "setScale, #n = " << n << ", value = " << value )

    // alias of outValues == inValues is no problem

    if( value == static_cast<ValueType>( 0 ) )
    {
        // Important : inValues might be undefined
        setVal( outValues, n, value );
        return;
    }

    if( value == static_cast<ValueType>( 1 ) )
    {
        set( outValues, inValues, n );
        return;
    }

    #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )

    for( IndexType i = 0; i < n; i++ )
    {
        outValues[i] = static_cast<ValueType>( inValues[i] ) * value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::sum( const ValueType array[], const IndexType n )
{
    LAMA_REGION( "OpenMP.Utils.sum" )

    LAMA_LOG_INFO( logger, "sum # array = " << array << ", n = " << n )
    ValueType val = static_cast<ValueType>( 0 );

    #pragma omp parallel shared( val )
    {
        ValueType tVal = static_cast<ValueType>( 0 );

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )

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
    LAMA_REGION( "OpenMP.Utils.setVal" )

    LAMA_LOG_DEBUG( logger, "setVal<" << getScalarType<ValueType>() << ">: " << "array[" << n << "] = " << val )

    #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )

    for( IndexType i = 0; i < n; ++i )
    {
        array[i] = val;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::setOrder( ValueType array[], const IndexType n )
{
    LAMA_REGION( "OpenMP.Utils.setOrder" )

    LAMA_LOG_DEBUG( logger,
                    "setOrder<" << getScalarType<ValueType>() << ">: " << "array[" << n << "] = 0, 1, 2, ..., " << ( n - 1 ) )

    #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )

    for( IndexType i = 0; i < n; ++i )
    {
        array[i] = static_cast<ValueType>( i );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::getValue( const ValueType* array, const IndexType i )
{
    LAMA_LOG_DEBUG( logger, "getValue<" << getScalarType<ValueType>() << ">: i = " << i )

    return array[i];
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::maxval( const ValueType array[], const IndexType n )
{
    LAMA_REGION( "OpenMP.Utils.maxVal" )

    LAMA_LOG_DEBUG( logger, "maxval<" << getScalarType<ValueType>() << ">: " << "array[" << n << "]" )

    ValueType val = static_cast<ValueType>( 0.0 );

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>( 0.0 );

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )

        for( IndexType i = 0; i < n; ++i )
        {
            if( array[i] > threadVal )
            {
                threadVal = array[i];
            }
        }

        #pragma omp critical
        {
            LAMA_LOG_TRACE( logger, "max val of thread = " << threadVal << ", global was " << val )

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
    LAMA_REGION( "OpenMP.Utils.absMaxVal" )

    LAMA_LOG_DEBUG( logger, "absMaxVal<" << getScalarType<ValueType>() << ">: " << "array[" << n << "]" )

    ValueType val = static_cast<ValueType>( 0.0 );

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>( 0.0 );

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )

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
            LAMA_LOG_TRACE( logger, "max val of thread  = " << threadVal << ", global was " << val )

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
    LAMA_REGION( "OpenMP.Utils.absMaxDiffVal" )

    LAMA_LOG_DEBUG( logger, "absMaxDiffVal<" << getScalarType<ValueType>() << ">: " << "array[" << n << "]" )

    ValueType val = static_cast<ValueType>( 0.0 );

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>( 0.0 );

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )

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
            LAMA_LOG_TRACE( logger, "max val of thread  = " << threadVal << ", global was " << val )

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
    LAMA_REGION( "OpenMP.Utils.isSorted" )

    LAMA_LOG_INFO( logger,
                   "isSorted<" << getScalarType<ValueType>() << ">, n = " << n << ", ascending = " << ascending )

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
    LAMA_REGION( "OpenMP.Utils.set" )

    LAMA_LOG_DEBUG( logger,
                    "set: out<" << getScalarType<ValueType1>() << "[" << n << "]" << " = in<" << getScalarType<ValueType2>() << ">[" << n << "]" )

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

    for( IndexType i = 0; i < n; i++ )
    {
        out[i] = static_cast<ValueType1>( in[i] );
    }
}

/* --------------------------------------------------------------------------- */

bool OpenMPUtils::validIndexes( const IndexType array[], const IndexType n, const IndexType size )
{
    LAMA_REGION( "OpenMP.Utils.validIndexes" )

    LAMA_LOG_DEBUG( logger, "validIndexes: array[" << n << "], size " << size )

    bool validFlag = true;

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE) reduction( & : validFlag )

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
    LAMA_REGION( "OpenMP.Utils.setGather" )

    LAMA_LOG_DEBUG( logger,
                    "setGather: out<" << getScalarType<ValueType1>() << ">[" << n << "]" << " = in<" << getScalarType<ValueType2>() << ">[ indexes[" << n << "] ]" )

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

    for( IndexType i = 0; i < n; i++ )
    {
        out[i] = static_cast<ValueType1>( in[indexes[i]] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1,typename ValueType2>
void OpenMPUtils::setScatter( ValueType1 out[], const IndexType indexes[], const ValueType2 in[], const IndexType n )
{
    LAMA_REGION( "OpenMP.Utils.setScatter" )

    LAMA_LOG_DEBUG( logger,
                    "setScatter: out<" << getScalarType<ValueType1>() << ">" << "[ indexes[" << n << "] ]" << " = in<" << getScalarType<ValueType2>() << ">[" << n << "]" )

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

    for( IndexType i = 0; i < n; i++ )
    {
        out[indexes[i]] = static_cast<ValueType1>( in[i] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::invert( ValueType array[], const IndexType n )
{
    LAMA_REGION( "OpenMP.Utils.invert" )

    LAMA_LOG_INFO( logger, "invert array[ " << n << " ]" )

    ValueType one = static_cast<ValueType>( 1.0 );

    #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )

    for( IndexType i = 0; i < n; ++i )
    {
        array[i] = one / array[i];
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPUtils::setInterface( UtilsInterface& Utils )
{
    // Instantations for IndexType, not done by ARITHMETIC_TYPE macrods

    LAMA_INTERFACE_REGISTER( Utils, validIndexes )

    // we keep the registrations for IndexType as we do not need conversions

    LAMA_INTERFACE_REGISTER_T( Utils, sum, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, setVal, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, setOrder, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, getValue, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, maxval, IndexType )
    LAMA_INTERFACE_REGISTER_T( Utils, isSorted, IndexType )

    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, IndexType, IndexType )
    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, IndexType, IndexType )
    LAMA_INTERFACE_REGISTER_TT( Utils, set, IndexType, IndexType )

#define LAMA_UTILS2_REGISTER(z, J, TYPE )                                       \
    LAMA_INTERFACE_REGISTER_TT( Utils, setScale, TYPE, ARITHMETIC_TYPE##J )     \
    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, TYPE, ARITHMETIC_TYPE##J )    \
    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, TYPE, ARITHMETIC_TYPE##J )   \
    LAMA_INTERFACE_REGISTER_TT( Utils, set, TYPE, ARITHMETIC_TYPE##J )          \

#define LAMA_UTILS_REGISTER(z, I, _)                                                 \
    LAMA_INTERFACE_REGISTER_T( Utils, scale, ARITHMETIC_TYPE##I )                    \
    LAMA_INTERFACE_REGISTER_T( Utils, sum, ARITHMETIC_TYPE##I )                      \
    LAMA_INTERFACE_REGISTER_T( Utils, setVal, ARITHMETIC_TYPE##I )                   \
    LAMA_INTERFACE_REGISTER_T( Utils, getValue, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( Utils, maxval, ARITHMETIC_TYPE##I )                   \
    LAMA_INTERFACE_REGISTER_T( Utils, absMaxVal, ARITHMETIC_TYPE##I )                \
    LAMA_INTERFACE_REGISTER_T( Utils, absMaxDiffVal, ARITHMETIC_TYPE##I )            \
    LAMA_INTERFACE_REGISTER_T( Utils, isSorted, ARITHMETIC_TYPE##I )                 \
    LAMA_INTERFACE_REGISTER_T( Utils, invert, ARITHMETIC_TYPE##I )                   \
    BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_UTILS2_REGISTER, ARITHMETIC_TYPE##I )

    BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_UTILS_REGISTER, _ )

#undef LAMA_UTILS_REGISTER
#undef LAMA_UTILS2_REGISTER

}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool OpenMPUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( memory::context::Host );
    setInterface( interface.Utils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPUtils::initialized = registerInterface();

} // namespace lama
