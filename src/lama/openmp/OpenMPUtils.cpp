/**
 * @file OpenMPUtils.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * $Id$
 */

// hpp
#include <lama/openmp/OpenMPUtils.hpp>

// others
#include <lama/LAMAInterface.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

#include <omp.h>

namespace lama
{

LAMA_LOG_DEF_LOGGER( OpenMPUtils::logger, "OpenMP.Utils" );

template<typename ValueType,typename OtherValueType>
void OpenMPUtils::scale( ValueType mValues[], const IndexType n, const OtherValueType value )
{
    LAMA_LOG_INFO( logger, "scale, #n = " << n << ", value = " << value );

    #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
    for( IndexType i = 0; i < n; i++ )
    {
        mValues[i] *= static_cast<ValueType>( value );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::sum( const ValueType array[], const IndexType n )
{
    LAMA_LOG_INFO( logger, "sum # array = " << array << ", n = " << n );
    ValueType val = static_cast<ValueType>( 0.0 );

    #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE ) reduction( +:val )
    for( IndexType i = 0; i < n; ++i )
    {
        val += array[i];
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::setVal( ValueType array[], const IndexType n, const ValueType val )
{
    LAMA_LOG_DEBUG( logger, "setVal<" << typeid(ValueType).name() << ">: " << "array[" << n << "] = " << val );

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
    LAMA_LOG_DEBUG( logger,
                    "setOrder<" << typeid(ValueType).name() << ">: " << "array[" << n << "] = 0, 1, 2, ..., " << ( n - 1 ) );

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
    LAMA_LOG_DEBUG( logger, "getValue<" << typeid(ValueType).name() << ">: i = " << i );

    return array[i];
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::maxval( const ValueType array[], const IndexType n )
{
    LAMA_LOG_DEBUG( logger, "maxval<" << typeid(ValueType).name() << ">: " << "array[" << n << "]" );

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
            LAMA_LOG_TRACE( logger, "max val of thread = " << threadVal << ", global was " << val );

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
    LAMA_LOG_DEBUG( logger, "absMaxVal<" << typeid(ValueType).name() << ">: " << "array[" << n << "]" );

    ValueType val = static_cast<ValueType>( 0.0 );

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>( 0.0 );

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
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
            LAMA_LOG_TRACE( logger, "max val of thread  = " << threadVal << ", global was " << val );

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
    LAMA_LOG_DEBUG( logger, "absMaxDiffVal<" << typeid(ValueType).name() << ">: " << "array[" << n << "]" );

    ValueType val = static_cast<ValueType>( 0.0 );

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>( 0.0 );

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
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
            LAMA_LOG_TRACE( logger, "max val of thread  = " << threadVal << ", global was " << val );

            if( threadVal > val )
            {
                val = threadVal;
            }
        }
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1,typename ValueType2>
void OpenMPUtils::set( ValueType1 out[], const ValueType2 in[], const IndexType n )
{
    LAMA_LOG_DEBUG( logger,
                    "set: out<" << typeid(ValueType1).name() << "[" << n << "]" << " = in<" << typeid(ValueType2).name() << ">[" << n << "]" );

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
    for( IndexType i = 0; i < n; i++ )
    {
        out[i] = static_cast<ValueType1>( in[i] );
    }
}

/* --------------------------------------------------------------------------- */

bool OpenMPUtils::validIndexes( const IndexType array[], const IndexType n, const IndexType size )
{
    LAMA_LOG_DEBUG( logger, "validIndexes: array[" << n << "], size " << size );

    bool validFlag = true;

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE) reduction( & : validFlag )
    for( IndexType i = 0; i < n; i++ )
    {
        if( size <= array[i] || 0 > array[i] )
        {
            // exception only in debug mode

            /*
             LAMA_THROWEXCEPTION( "array[" << i << "] = " << array[i]
             << " is illegal index, size = " << size );
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
    LAMA_LOG_DEBUG( logger,
                    "setGather: out<" << typeid(ValueType1).name() << ">[" << n << "]" << " = in<" << typeid(ValueType2).name() << ">[ indexes[" << n << "] ]" );

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
    LAMA_LOG_DEBUG( logger,
                    "setScatterr: out<" << typeid(ValueType1).name() << ">" << "[ indexes[" << n << "] ]" << " = in<" << typeid(ValueType2).name() << ">[" << n << "]" );

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
    ValueType one = static_cast<ValueType>( 1.0 );

    #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
    for( IndexType i = 0; i < n; ++i )
    {
        array[i] = one / array[i];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::scaleVal( ValueType array[], const IndexType n, const ValueType val )
{
    #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
    for( IndexType i = 0; i < n; ++i )
    {
        array[i] *= val;
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPUtils::setInterface( UtilsInterface& Utils )
{
    LAMA_INTERFACE_REGISTER( Utils, validIndexes );

    LAMA_INTERFACE_REGISTER_TT( Utils, scale, float, float );
    LAMA_INTERFACE_REGISTER_TT( Utils, scale, double, float );
    LAMA_INTERFACE_REGISTER_TT( Utils, scale, float, double );
    LAMA_INTERFACE_REGISTER_TT( Utils, scale, double, double );

    LAMA_INTERFACE_REGISTER_T( Utils, sum, IndexType );
    LAMA_INTERFACE_REGISTER_T( Utils, sum, float );
    LAMA_INTERFACE_REGISTER_T( Utils, sum, double );

    LAMA_INTERFACE_REGISTER_T( Utils, setVal, IndexType );
    LAMA_INTERFACE_REGISTER_T( Utils, setVal, float );
    LAMA_INTERFACE_REGISTER_T( Utils, setVal, double );

    LAMA_INTERFACE_REGISTER_T( Utils, setOrder, IndexType );

    LAMA_INTERFACE_REGISTER_T( Utils, getValue, IndexType );
    LAMA_INTERFACE_REGISTER_T( Utils, getValue, float );
    LAMA_INTERFACE_REGISTER_T( Utils, getValue, double );

    LAMA_INTERFACE_REGISTER_T( Utils, maxval, IndexType );
    LAMA_INTERFACE_REGISTER_T( Utils, maxval, float );
    LAMA_INTERFACE_REGISTER_T( Utils, maxval, double );

    LAMA_INTERFACE_REGISTER_T( Utils, absMaxVal, float );
    LAMA_INTERFACE_REGISTER_T( Utils, absMaxVal, double );

    LAMA_INTERFACE_REGISTER_T( Utils, absMaxDiffVal, float );
    LAMA_INTERFACE_REGISTER_T( Utils, absMaxDiffVal, double );

    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, int, int );

    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, float, float );
    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, double, float );
    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, float, double );
    LAMA_INTERFACE_REGISTER_TT( Utils, setScatter, double, double );

    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, int, int );

    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, float, float );
    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, float, double );
    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, double, float );
    LAMA_INTERFACE_REGISTER_TT( Utils, setGather, double, double );

    LAMA_INTERFACE_REGISTER_TT( Utils, set, int, int );

    LAMA_INTERFACE_REGISTER_TT( Utils, set, float, float );
    LAMA_INTERFACE_REGISTER_TT( Utils, set, float, double );
    LAMA_INTERFACE_REGISTER_TT( Utils, set, double, float );
    LAMA_INTERFACE_REGISTER_TT( Utils, set, double, double );

    LAMA_INTERFACE_REGISTER_T( Utils, invert, float );
    LAMA_INTERFACE_REGISTER_T( Utils, invert, double );
}

} // namespace lama
