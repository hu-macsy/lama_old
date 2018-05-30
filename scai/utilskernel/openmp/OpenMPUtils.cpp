/**
 * @file OpenMPUtils.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Implementation of CSR utilities with OpenMP for the Inteface
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// hpp
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>

// local library
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/SparseKernelTrait.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tracing.hpp>

#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/Utils.hpp>
#include <scai/common/OpenMP.hpp>

#include <algorithm>
#include <memory>

//#include <parallel/sort.h>

namespace scai
{

using common::BinaryOp;
using common::CompareOp;
using common::UnaryOp;
using common::TypeTraits;

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( OpenMPUtils::logger, "OpenMP.Utils" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::reduceSum( const ValueType array[], const IndexType n, const ValueType zero )
{
    SCAI_REGION( "OpenMP.Utils.reduceSum" )

    ValueType val = zero;

    #pragma omp parallel shared( val )
    {
        ValueType threadVal = zero;

        #pragma omp for 

        for ( IndexType i = 0; i < n; ++i )
        {
            threadVal += array[i];
        }

        atomicAdd( val, threadVal );
    }
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::reduceMaxVal( const ValueType array[], const IndexType n, const ValueType zero )
{
    SCAI_REGION( "OpenMP.Utils.reduceMaxVal" )

    ValueType val = zero;

    #pragma omp parallel
    {
        ValueType threadVal = zero;

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
            if ( threadVal > val )
            {
                val = threadVal;
            }
        }
    }
    return val;
}

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
ComplexFloat OpenMPUtils::reduceMaxVal( const ComplexFloat[], const IndexType, const ComplexFloat )
{
    COMMON_THROWEXCEPTION( "minval not supported for complex arrays." )
}

template<>
ComplexDouble OpenMPUtils::reduceMaxVal( const ComplexDouble[], const IndexType, const ComplexDouble )
{
    COMMON_THROWEXCEPTION( "minval not supported for complex arrays." )
}

template<>
ComplexLongDouble OpenMPUtils::reduceMaxVal( const ComplexLongDouble[], const IndexType, const ComplexLongDouble )
{
    COMMON_THROWEXCEPTION( "minval not supported for complex arrays." )
}

#endif 

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::reduceMinVal( const ValueType array[], const IndexType n, const ValueType zero )
{
    SCAI_REGION( "OpenMP.Utils.reduceMinVal" )

    ValueType val = zero;

    #pragma omp parallel
    {
        ValueType threadVal = zero;

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
            if ( threadVal < val )
            {
                val = threadVal;
            }
        }
    }
    return val;
}

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
ComplexFloat OpenMPUtils::reduceMinVal( const ComplexFloat[], const IndexType, const ComplexFloat )
{
    COMMON_THROWEXCEPTION( "minval not supported for complex arrays." )
}

template<>
ComplexDouble OpenMPUtils::reduceMinVal( const ComplexDouble[], const IndexType, const ComplexDouble )
{
    COMMON_THROWEXCEPTION( "minval not supported for complex arrays." )
}

template<>
ComplexLongDouble OpenMPUtils::reduceMinVal( const ComplexLongDouble[], const IndexType, const ComplexLongDouble )
{
    COMMON_THROWEXCEPTION( "minval not supported for complex arrays." )
}

#endif 

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::reduceBinOp(
    const ValueType array[],
    const IndexType n,
    const ValueType zero,
    const BinaryOp op )
{
    SCAI_REGION( "OpenMP.Utils.reduceBinOp" )

    ValueType val = zero;

    #pragma omp parallel
    {
        ValueType threadVal = zero;

        #pragma omp for 

        for ( IndexType i = 0; i < n; ++i )
        {
            threadVal = applyBinary( threadVal, op, array[i] );
        }

        #pragma omp critical
        {
            val = applyBinary( val, op, threadVal );
        }
    }
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::reduceAbsMaxVal( const ValueType array[], const IndexType n, const ValueType zero )
{
    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    SCAI_REGION( "OpenMP.Utils.reduceAbsMaxVal" )

    RealType val = zero;

    #pragma omp parallel
    {
        RealType threadVal = zero;

        #pragma omp for 

        for ( IndexType i = 0; i < n; ++i )
        {
            RealType elem = common::Math::abs( array[i] );

            if ( elem > threadVal )
            {
                threadVal = elem;
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
ValueType OpenMPUtils::reduce(
    const ValueType array[],
    const IndexType n,
    const ValueType zero,
    const BinaryOp op )
{
    SCAI_LOG_INFO ( logger, "reduce # array<" << TypeTraits<ValueType>::id() << ">[" << n << "], op = " << op )

    // Own reduction methods for specific binary op are more efficient than a general solution.
    // Performance gain is a factor of 1.3 - 1.4

    switch ( op )
    {
        case BinaryOp::ADD :
            return reduceSum( array, n, zero );

        case BinaryOp::MIN :
            return reduceMinVal( array, n, zero );

        case BinaryOp::MAX :
            return reduceMaxVal( array, n, zero );

        case BinaryOp::ABS_MAX :
            return reduceAbsMaxVal( array, n, zero );

        default:
            return reduceBinOp( array, n, zero, op );
    }

    return ValueType( 0 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::setVal( ValueType array[], const IndexType n, const ValueType val, const BinaryOp op )
{
    SCAI_REGION( "OpenMP.Utils.setVal" )
    SCAI_LOG_INFO( logger, "setVal<" << TypeTraits<ValueType>::id() << ">: " << "array[" << n << "] = "
                   << val << ", op = " << op )
    ValueType value = static_cast<ValueType>( val );

    switch ( op )
    {
        case BinaryOp::COPY :
        {
            #pragma omp parallel
            {
                #pragma omp for 

                for ( IndexType i = 0; i < n; i++ )
                {
                    array[i] = value;
                }
            }
            break;
        }

        case BinaryOp::ADD :
        {
            if ( val == common::Constants::ZERO )
            {
                return;
            }

            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                array[i] += value;
            }

            break;
        }

        case BinaryOp::SUB :
        {
            if ( val == common::Constants::ZERO )
            {
                return;
            }

            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                array[i] -= value;
            }

            break;
        }

        case BinaryOp::MULT :
        {
            // scale all values of the array
            if ( val == common::Constants::ONE )
            {
                // skip it
            }
            else if ( val == common::Constants::ZERO )
            {
                setVal( array, n, ValueType( 0 ), BinaryOp::COPY );
            }
            else
            {
                #pragma omp parallel for 

                for ( IndexType i = 0; i < n; i++ )
                {
                    array[i] *= value;
                }
            }

            break;
        }

        case BinaryOp::DIVIDE :
        {
            // scale all values of the array
            if ( val == common::Constants::ONE )
            {
                // skip it
            }
            else if ( val == common::Constants::ZERO )
            {
                COMMON_THROWEXCEPTION( "DIVIDE by ZERO, val = " << val )
            }
            else
            {
                #pragma omp parallel for 

                for ( IndexType i = 0; i < n; i++ )
                {
                    array[i] /= value;
                }
            }

            break;
        }

        default:
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                array[i] = applyBinary( array[i], op, value );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::scaleVectorAddScalar( ValueType array1[], const ValueType array2[], const IndexType n,
                                        const ValueType alpha, const ValueType beta )
{
    SCAI_REGION( "OpenMP.Utils.scaleVectorAddScalar" )
    SCAI_LOG_DEBUG( logger,
                    "scaleVectorAddScalar<" << TypeTraits<ValueType>::id() << ">: " << "array1[" << n << "] = alpha( " << alpha
                    << ") * array2 + beta (" << beta << ")" )
    #pragma omp parallel for 

    for ( IndexType i = 0; i < n; ++i )
    {
        array1[i] = alpha * array2[i] + beta;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::setOrder( ValueType array[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.setOrder" )
    SCAI_LOG_DEBUG( logger,
                    "setOrder<" << TypeTraits<ValueType>::id() << ">: " << "array[" << n << "] = 0, 1, 2, ..., " << ( n - 1 ) )
    #pragma omp parallel for 

    for ( IndexType i = 0; i < n; ++i )
    {
        array[i] = static_cast<ValueType>( i );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::setSequence( ValueType array[], const ValueType startValue, const ValueType inc, const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.setSequence" )
    SCAI_LOG_DEBUG( logger, "setSequence<" << TypeTraits<ValueType>::id() << ">: " << "array[" << n << "] = " << startValue
                    << "..., " << ( startValue + static_cast<ValueType>( n - 1 ) * inc ) )

    #pragma omp parallel for 

    for ( IndexType i = 0; i < n; ++i )
    {
        // use static cast to avoid any mixed binary operation
        array[i] = startValue + static_cast<ValueType>( i ) * inc;
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
    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    SCAI_REGION( "OpenMP.Utils.absMaxDiffVal" )
    SCAI_LOG_DEBUG( logger, "absMaxDiffVal<" << TypeTraits<ValueType>::id() << ">: " << "array[" << n << "]" )

    RealType val = static_cast<RealType>( 0 );

    #pragma omp parallel
    {
        RealType threadVal = static_cast<RealType>( 0 );

        #pragma omp for 

        for ( IndexType i = 0; i < n; ++i )
        {
            RealType elem = common::Math::abs( array1[i] - array2[i] );

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

template <typename ValueType>
ValueType OpenMPUtils::reduce2(
    const ValueType array1[],
    const ValueType array2[],
    const IndexType n,
    const BinaryOp binOp,
    const ValueType zero,
    const BinaryOp redOp )
{
    if ( binOp == BinaryOp::SUB && redOp == BinaryOp::ABS_MAX )
    {
        return absMaxDiffVal( array1, array2, n );
    }

    SCAI_REGION( "OpenMP.Utils.reduce2" )

    SCAI_LOG_DEBUG( logger, "reduce2<" << TypeTraits<ValueType>::id() << ">: " << "arr1,2[" << n << "]" )

    ValueType val = zero;

    #pragma omp parallel
    {
        ValueType threadVal = zero;

        #pragma omp for 

        for ( IndexType i = 0; i < n; ++i )
        {
            ValueType elem = applyBinary( array1[i], binOp, array2[i] );
            threadVal = applyBinary( threadVal, redOp, elem );
        }

        #pragma omp critical
        {
            val = applyBinary( val, redOp, threadVal );
        }
    }
    return val;
}

/* --------------------------------------------------------------------------- */

template <typename ValueType>
bool OpenMPUtils::allCompare(
    const ValueType array1[],
    const ValueType array2[],
    const IndexType n,
    const CompareOp op )
{
    bool val = true;
 
    #pragma omp parallel
    {
        bool threadVal = true;

        #pragma omp for 

        for ( IndexType i = 0; i < n; ++i )
        {
            bool elem = compare( array1[i], op, array2[i] );
            threadVal = threadVal && elem;
        }

        #pragma omp critical
        {
            val = val && threadVal;
        }
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template <typename ValueType>
bool OpenMPUtils::allCompareScalar(
    const ValueType array[],
    const ValueType scalar,
    const IndexType n,
    const CompareOp op )
{
    bool val = true;

    #pragma omp parallel
    {
        bool threadVal = true;

        #pragma omp for 

        for ( IndexType i = 0; i < n; ++i )
        {
            bool elem = compare( array[i], op, scalar );
            threadVal = threadVal && elem;
        }

        #pragma omp critical
        {
            val = val && threadVal;
        }
    }

    return val;
}

/* ------------------------------------------------------------------------------------------------------------------ */

void OpenMPUtils::setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.inversePerm" )

    SCAI_LOG_INFO( logger, "compute inverse perm, n = " << n )

    // Parallel execution is safe as perm does not contain a value twice

    #pragma omp parallel for 

    for ( IndexType ii = 0; ii < n; ii++ )
    {
        IndexType i = perm[ii];
        SCAI_ASSERT_VALID_INDEX_DEBUG( i, n, "permutation value out of range, perm[" << ii << "]" )
        inversePerm[i] = ii;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool OpenMPUtils::isSorted( const ValueType array[], const IndexType n, const CompareOp op )
{
    SCAI_REGION( "OpenMP.Utils.isSorted" )
    SCAI_LOG_INFO( logger,
                   "isSorted<" << TypeTraits<ValueType>::id() << ">, n = " << n << ", op = " << op )
    bool sorted = true; //!< will be set to false at violations

    for ( IndexType i = 1; i < n; i++ )
    {
        if ( !compare( array[i - 1], op, array[i] ) )
        {
            sorted = false;
            break;
        }
    }

    return sorted;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
void OpenMPUtils::set( ValueType1 out[], const ValueType2 in[], const IndexType n, const BinaryOp op )
{
    SCAI_REGION( "OpenMP.Utils.set" )

    SCAI_LOG_DEBUG( logger,
                    "set: out<" << TypeTraits<ValueType1>::id() << "[" << n << "]"
                    << ", op = " << op << "  in<" << TypeTraits<ValueType2>::id() << ">[" << n << "]" )

    switch ( op )
    {
        case BinaryOp::COPY :
        {
            if ( in != reinterpret_cast<ValueType2*> ( out ) )
            {
                #pragma omp parallel for 

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = static_cast<ValueType1>( in[i] );
                }
            }

            break;
        }

        case BinaryOp::ADD :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] += static_cast<ValueType1>( in[i] );
            }

            break;
        }

        case BinaryOp::SUB :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] -= static_cast<ValueType1>( in[i] );
            }

            break;
        }

        case BinaryOp::DIVIDE :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] /= static_cast<ValueType1>( in[i] );
            }

            break;
        }

        case BinaryOp::MULT :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] *= static_cast<ValueType1>( in[i] );
            }

            break;
        }

        default:
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] = applyBinary( out[i], op, static_cast<ValueType1>( in[i] ) );
            }

            break;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
void OpenMPUtils::setSection(
    ValueType1 out[],
    const IndexType inc1,
    const ValueType2 in[],
    const IndexType inc2,
    const IndexType n,
    const BinaryOp op )
{
    SCAI_REGION( "OpenMP.Utils.setSection" )

    SCAI_LOG_DEBUG( logger,
                    "setSection: out<" << TypeTraits<ValueType1>::id() << "[" << n << "]"
                    << ", op = " << op << "  in<" << TypeTraits<ValueType2>::id() << ">[" << n << "]" )

    switch ( op )
    {
        case BinaryOp::COPY :
        {
            if ( in != reinterpret_cast<ValueType2*> ( out ) )
            {
                #pragma omp parallel for 

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i * inc1] = static_cast<ValueType1>( in[i * inc2] );
                }
            }

            break;
        }

        case BinaryOp::ADD :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i * inc1] += static_cast<ValueType1>( in[i * inc2] );
            }

            break;
        }

        case BinaryOp::SUB :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i * inc1] -= static_cast<ValueType1>( in[i * inc2] );
            }

            break;
        }

        case BinaryOp::DIVIDE :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i * inc1] /= static_cast<ValueType1>( in[i * inc2] );
            }

            break;
        }

        case BinaryOp::MULT :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i * inc1] *= static_cast<ValueType1>( in[i * inc2] );
            }

            break;
        }

        default:
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i * inc1] = applyBinary( out[i * inc1], op, static_cast<ValueType1>( in[i * inc2] ) );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::fillSection(
    ValueType out[],
    const IndexType inc,
    const ValueType val,
    const IndexType n,
    const BinaryOp op )
{
    SCAI_REGION( "OpenMP.Utils.fillSection" )

    SCAI_LOG_DEBUG( logger,
                    "fillSection: out<" << TypeTraits<ValueType>::id() << "[" << n << "]"
                    << ", op = " << op << ", val = " << val )

    switch ( op )
    {
        case BinaryOp::COPY :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i * inc] = val;
            }

            break;
        }

        case BinaryOp::ADD :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i * inc] += val;
            }

            break;
        }

        case BinaryOp::SUB :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i * inc] -= val;
            }

            break;
        }

        case BinaryOp::DIVIDE :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i * inc] /= val;
            }

            break;
        }

        case BinaryOp::MULT :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i * inc] *= val;
            }

            break;
        }

        default:
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i * inc] = applyBinary( out[i * inc], op, val );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
static void floorOp( ValueType out[], const ValueType in[], const IndexType n )
{
    #pragma omp parallel for 

    for ( IndexType i = 0; i < n; i++ )
    {
        out[i] = common::Math::floor( in[i] );
    }
}

// template specialization for IndexType as floor not available on IndexType

template<>
void floorOp( IndexType out[], const IndexType in[], const IndexType n )
{
    #pragma omp parallel for 

    for ( IndexType i = 0; i < n; i++ )
    {
        out[i] = in[i];
    }
}

template<typename ValueType>
void OpenMPUtils::unaryOp( ValueType out[], const ValueType in[], const IndexType n, const UnaryOp op )
{
    SCAI_REGION( "OpenMP.Utils.unaryOp" )

    SCAI_LOG_DEBUG( logger, "unaryOp<" << TypeTraits<ValueType>::id() << ", op = " << op << ">, n = " << n )

    if ( n <= 0 )
    {
        return;
    }

    switch ( op )
    {
        case UnaryOp::MINUS :
        {
            {
                #pragma omp parallel for 

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = -in[i];
                }
            }

            break;
        }

        case UnaryOp::ABS :
        {
            {
                #pragma omp parallel for 

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = common::Math::abs( in[i] );
                }
            }

            break;
        }

        case UnaryOp::FLOOR :
        {
            floorOp( out, in, n );

            break;
        }

        default:
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] = applyUnary( op, in[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::binaryOpScalar(
    ValueType out[],
    const ValueType in[],
    const ValueType value,
    const IndexType n,
    const BinaryOp op,
    bool swapScalar )
{
    SCAI_REGION( "OpenMP.Utils.binOpScalar" )

    SCAI_LOG_DEBUG( logger, "binaryOpScalar1<" << TypeTraits<ValueType>::id() << ", op = " << op << ">, n = " << n )

    if ( n <= 0 )
    {
        return;
    }

    switch ( op )
    {
        case BinaryOp::ADD :
        {
            // ignore swapScalar, does not matter here

            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] = in[i] + value;
            }

            break;
        }

        case BinaryOp::SUB :
        {
            if ( swapScalar )
            {
                #pragma omp parallel for 
                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = value - in[i];
                }
            }
            else
            {
                #pragma omp parallel for 
                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = in[i] - value;
                }
            }
            break;
        }

        case BinaryOp::MULT :
        {
            if ( value == common::Constants::ZERO )
            {
                // Important : in might be undefined

                setVal( out, n, value, BinaryOp:: COPY );
            }
            else if ( value == common::Constants::ONE )
            {
                set( out, in, n, BinaryOp:: COPY );
            }
            else
            {
                #pragma omp parallel for 

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = in[i] * value;
                }
            }

            break;
        }

        case BinaryOp::DIVIDE :
        {
            if ( swapScalar )
            {
                #pragma omp parallel for 

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = value / in[i];
                }
            }
            else
            {
                #pragma omp parallel for 

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = in[i] / value;
                }
            }
            break;
        }
        default:
        {
            if ( swapScalar )
            {
                #pragma omp parallel for 

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = applyBinary( value, op, in[i] );
                }
            }
            else
            {
                #pragma omp parallel for 

                for ( IndexType i = 0; i < n; i++ )
                {
                    out[i] = applyBinary( in[i], op, value );
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::binaryOp( ValueType out[], const ValueType in1[], const ValueType in2[], const IndexType n, const BinaryOp op )
{
    SCAI_REGION( "OpenMP.Utils.binOp" )

    SCAI_LOG_DEBUG( logger, "binaryOp<" << TypeTraits<ValueType>::id() << ", op = " << op << ">, n = " << n )

    if ( n <= 0 )
    {
        return;
    }

    // using parallel loops with a specific operator might give a slight performance gain
    // but did not show any relevant gain when full compile time optimization is enabled.
    // But make sure: applyBinary must be inlined

    switch ( op )
    {
        case BinaryOp::ADD :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] = in1[i] + in2[i];
            }

            break;
        }

        case BinaryOp::MULT :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] = in1[i] * in2[i];
            }

            break;
        }

        default:
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] = applyBinary( in1[i], op, in2[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

bool OpenMPUtils::validIndexes( const IndexType array[], const IndexType n, const IndexType size )
{
    SCAI_REGION( "OpenMP.Utils.validIndexes" )
    SCAI_LOG_DEBUG( logger, "validIndexes: array[" << n << "], size " << size )

    IndexType invalid = 0;

    // #pragma omp parallel for 

    #pragma omp parallel
    {
        IndexType tInvalid = 0;  // each thread counts invalid indexes for its part

        #pragma omp parallel for 

        for ( IndexType i = 0; i < n; i++ )
        {
            if ( ! common::Utils::validIndex( array[i], size ) )
            {
                tInvalid++;
            }
        }

        atomicAdd( invalid, tInvalid );
    }

    SCAI_LOG_INFO( logger, "#invalid indexes = " << invalid << " of " << n )

    return invalid == 0;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
void OpenMPUtils::setGather(
    ValueType1 out[],
    const ValueType2 in[],
    const IndexType indexes[],
    const BinaryOp op,
    const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.setGather" )

    SCAI_LOG_DEBUG( logger,
                    "setGather: out<" << TypeTraits<ValueType1>::id() << ">[" << n << "] "
                    << op << " = in<" << TypeTraits<ValueType2>::id() << ">[ indexes[" << n << "] ]" )

    switch ( op )
    {
        case BinaryOp::COPY :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] = static_cast<ValueType1>( in[indexes[i]] );
            }

            break;
        }

        case BinaryOp::ADD :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] += static_cast<ValueType1>( in[indexes[i]] );
            }

            break;
        }

        case BinaryOp::SUB :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] -= static_cast<ValueType1>( in[indexes[i]] );
            }

            break;
        }

        case BinaryOp::MULT :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] *= static_cast<ValueType1>( in[indexes[i]] );
            }

            break;
        }

        case BinaryOp::DIVIDE :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] /= static_cast<ValueType1>( in[indexes[i]] );
            }

            break;
        }

        default:
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[i] = applyBinary( out[i], op, static_cast<ValueType1>( in[indexes[i]] ) );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

static IndexType _binarySearch( const IndexType indexes[], const IndexType n, const IndexType pos )
{
    IndexType first = 0;
    IndexType last  = n;

    while ( first < last )
    {
        IndexType middle = first + ( last - first ) / 2;

        if ( indexes[middle] == pos )
        {
            return middle;
        }
        else if ( indexes[middle] > pos )
        {
            last = middle;
        }
        else
        {
            first = middle + 1;
        }
    }

    return invalidIndex;
}

/* --------------------------------------------------------------------------- */

void OpenMPUtils::binarySearch( IndexType outPos[], 
                                const IndexType indexes[], const IndexType m, 
                                const IndexType inPos[], const IndexType n )
{
    #pragma omp parallel for 

    for ( IndexType i = 0; i < n; i++ )
    {
        outPos[i] = _binarySearch( indexes, m, inPos[i] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
void OpenMPUtils::setGatherSparse(
    ValueType1 target[],
    const ValueType2 sourceZeroValue,
    const ValueType2 sourceNonZeroValues[],
    const IndexType sourceNonZeroIndexes[],
    const IndexType sourceNNZ,
    const IndexType indexes[],
    const BinaryOp op,
    const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.setGatherSparse" )

    SCAI_LOG_INFO( logger,
                   "setGatherSparse: target<" << TypeTraits<ValueType1>::id() << ">[" << n << "] "
                    << op << " = sourceSparse<" << TypeTraits<ValueType2>::id() << ">[ indexes[" << n << "] ]" )

    #pragma omp parallel for 

    for ( IndexType i = 0; i < n; i++ )
    {
        // Binary search possible as source index array must be sorted

        IndexType k = _binarySearch( sourceNonZeroIndexes, sourceNNZ, indexes[i] );

        SCAI_LOG_TRACE( logger, "binarySearch( " << indexes[i] << " ) -> " << k )

        // default value if value not availabe in sparse input array

        ValueType2 sourceVal = k == invalidIndex ? sourceZeroValue : sourceNonZeroValues[k];

        target[i] = applyBinary( target[i], op, static_cast<ValueType1>( sourceVal ) );
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
    #pragma omp parallel for 

    for ( IndexType i = 0; i < n; i++ )
    {
        out[indexes[i]] = value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
void OpenMPUtils::setScatter(
    ValueType1 out[],
    const IndexType indexes[],
    const bool unique,
    const ValueType2 in[],
    const BinaryOp op,
    const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.setScatter" )

    SCAI_LOG_INFO( logger, "setScatter: out<" << TypeTraits<ValueType1>::id() << ">"
                   << "[ indexes[" << n << "], unique = " << unique << " ]"
                   << op << " = in<" << TypeTraits<ValueType2>::id() << ">[" << n << "]" )

    SCAI_LOG_DEBUG( logger, "addresses: out = " << out << ", indexes = " << indexes << ", in = " << in )

    if ( op == BinaryOp::COPY )
    {
        #pragma omp parallel for 

        for ( IndexType i = 0; i < n; i++ )
        {
            out[indexes[i]] = static_cast<ValueType1>( in[i] );
        }
    }
    else if ( op == BinaryOp::ADD )
    {
        #pragma omp parallel for 

        for ( IndexType i = 0; i < n; i++ )
        {
            atomicAdd( out[indexes[i]], static_cast<ValueType1>( in[i] ) );
        }
    }
    else if ( op == BinaryOp::SUB )
    {
        #pragma omp parallel for 

        for ( IndexType i = 0; i < n; i++ )
        {
            atomicAdd( out[indexes[i]], - static_cast<ValueType1>( in[i] ) );
        }
    }
    else
    {
        if ( unique )
        {
            // no double indexes, we can do it parallel

            #pragma omp parallel for 

            for ( IndexType i = 0; i < n; i++ )
            {
                out[indexes[i]] = applyBinary( out[indexes[i]], op, static_cast<ValueType1>( in[i] ) );
            }
        }
        else
        {
            // there might be double indexes, just do it serially

            for ( IndexType i = 0; i < n; i++ )
            {
                out[indexes[i]] = applyBinary( out[indexes[i]], op, static_cast<ValueType1>( in[i] ) );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::scanSerial( ValueType array[], const IndexType n, ValueType first, bool exclusive )
{
    SCAI_LOG_DEBUG( logger, "scanSerial: " << n << " entries, first = " << first << ", exclusive = " << exclusive )

    ValueType runningSum = first;

    if ( exclusive )
    {
        for ( IndexType i = 0; i < n; i++ )
        {
            ValueType tmp = runningSum;
            runningSum += array[i];
            array[i] = tmp;
        }
    }
    else
    {
        for ( IndexType i = 0; i < n; i++ )
        {
            runningSum += array[i];
            array[i] = runningSum;
        }
    }

    return runningSum;;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::scanParallel( PartitionId numThreads, ValueType array[], const IndexType n, const ValueType zero, const bool exclusive )
{
    std::unique_ptr<ValueType[]> threadValues( new ValueType[numThreads] );

    SCAI_LOG_DEBUG( logger, "scanParallel: " << n << " entries for " << numThreads << " threads" )

    ValueType runningSum;

    #pragma omp parallel
    {
        IndexType lb, ub;

        omp_get_my_range( lb, ub, n );

        ValueType myLocalSum = 0;

        for ( IndexType i = lb; i < ub; i++ )
        {
            myLocalSum += array[i];
        }

        SCAI_LOG_TRACE( logger, "scanParallel: local sum on " << lb << " - " << ub << " is " << myLocalSum )

        threadValues[omp_get_thread_num()] = myLocalSum;

        // Important: barrier before and after serial scan are mandatory

        #pragma omp barrier

        #pragma omp master
        {
            bool threadExclusive = true;
            runningSum = scanSerial( threadValues.get(), numThreads, zero, threadExclusive );
            SCAI_LOG_TRACE( logger, "runningSum for all = " << runningSum )
        }

        #pragma omp barrier

        // Each thread has now its start value

        ValueType myRunningSum = threadValues[omp_get_thread_num()];

        if ( lb < ub )
        {
            // local part on this thread has at least one element
            scanSerial( array + lb, ub - lb, myRunningSum, exclusive );
        }
    }

    return runningSum;;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::scan( ValueType array[], const IndexType n, ValueType zero, const bool exclusive, const bool append )
{
    SCAI_REGION( "OpenMP.Utils.scan" )

    int numThreads = 1; // will be set to available threads in parallel region

    #pragma omp parallel
    #pragma omp master
    {
        numThreads = omp_get_num_threads();
    }

    SCAI_LOG_DEBUG( logger, "scan " << n << " entries, #threads = " << numThreads << ", append = " << append )

    static int minThreads = 3;
    ValueType total;

    if ( numThreads < minThreads )
    {
        total = scanSerial( array, n, zero, exclusive );
    }
    else
    {
        total = scanParallel( numThreads, array, n, zero, exclusive );
    }

    if ( append )
    {
        SCAI_LOG_TRACE( logger, "append total = " << total << " at pos = " << n )

        array[n] = total;
    }

    return total;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPUtils::unscan( ValueType array[], const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.unscan" )
    SCAI_LOG_DEBUG( logger, "unscan array[ " << n << " ]" )

    if ( n == 0 )
    {
        return 0;
    }

    ValueType first = array[0];

    // due to anti-dependencies we do it serial, otherwise temporary required

    for ( IndexType i = 0; i < n; i++ )
    {
        array[i] = array[i + 1] - array[i];
    }

    return first;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::sortValues(
    ValueType outValues[],
    const ValueType inValues[],
    const IndexType n,
    const bool ascending )
{
    if ( outValues != inValues )
    {
        set( outValues, inValues, n, BinaryOp::COPY );
        sortValues( outValues, outValues, n, ascending );
        return;
    }

    SCAI_REGION( "OpenMP.Utils.sortValues" )

    // sort using a custom function object

    struct compare
    {
        static bool greater ( const ValueType a, const ValueType b )
        {
            return a > b;
        }
        static bool less ( const ValueType a, const ValueType b )
        {
            return a < b;
        }
    };

    if ( ascending )
    {
        std::sort( outValues, outValues + n, compare::less );
    }
    else
    {
        std::sort( outValues, outValues + n, compare::greater );
    }
}

static const void* ptr = NULL;

template<typename ValueType>
void OpenMPUtils::sort(
    IndexType perm[],
    ValueType outValues[],
    const ValueType inValues[],
    const IndexType n,
    const bool ascending )
{
    if ( perm == NULL )
    {
        sortValues( outValues, inValues, n, ascending );
        return;
    }

    SCAI_REGION( "OpenMP.Utils.sort" )

    for ( IndexType i = 0; i < n; ++i )
    {
        perm[i] = i;
    }

    // sort using a custom function object
    struct compare
    {
        static bool isAscending ( const IndexType a, const IndexType b )
        {
            const ValueType* arr = reinterpret_cast<const ValueType*>( ptr );
            return arr[a] < arr[b];
        }
        static bool isDescending ( const IndexType a, const IndexType b )
        {
            const ValueType* arr = reinterpret_cast<const ValueType*>( ptr );
            return arr[a] > arr[b];
        }
    };

    // comparison function needs a global array

    ptr = inValues;

    if ( ascending )
    {
        std::stable_sort( perm, perm + n, compare::isAscending );
    }
    else
    {
        std::stable_sort( perm, perm + n, compare::isDescending );
    }

    // now use perm to compute outValues if required

    if ( outValues == inValues )
    {
        std::unique_ptr<ValueType[]> tmp( new ValueType[n] );
        set( tmp.get(), inValues, n, BinaryOp::COPY );
        setGather( outValues, tmp.get(), perm, BinaryOp::COPY, n );
    }
    else if ( outValues != NULL )
    {
        setGather( outValues, inValues, perm, BinaryOp::COPY, n );
    }
}

#ifdef SCAI_COMPLEX_SUPPORTED

template<>
void OpenMPUtils::sort( IndexType[], ComplexFloat[], const ComplexFloat[], const IndexType, const bool )
{
    COMMON_THROWEXCEPTION( "sort unsupported for complex values" )
}

template<>
void OpenMPUtils::sort( IndexType[], ComplexDouble[], const ComplexDouble[], const IndexType, const bool )
{
    COMMON_THROWEXCEPTION( "sort unsupported for complex values" )
}

template<>
void OpenMPUtils::sort( IndexType[], ComplexLongDouble[], const ComplexLongDouble[], const IndexType, const bool )
{
    COMMON_THROWEXCEPTION( "sort unsupported for complex values" )
}

#endif

/* --------------------------------------------------------------------------- */

template<typename KeyType, typename ValueType>
void OpenMPUtils::qsort( KeyType keys[], ValueType values[], IndexType left, IndexType right, bool ascending )
{
    SCAI_LOG_INFO( logger, "qsort( " << left << " : " << right << " )" )

    if ( left >= right )
    {
        return;
    }

    IndexType i = left;
    IndexType j = right;
    IndexType m = ( left + right ) / 2;

    KeyType pivot = keys[ m ];

    while ( i <= j )
    {
        if ( ascending ) 
        {
            while ( keys[i] < pivot )
            {
                ++i;
            }
            while ( keys[j] > pivot )
            {
                --j;
            }
        }
        else
        {
            while ( keys[i] > pivot )
            {
                ++i;
            }
            while ( keys[j] < pivot )
            {
                --j;
            }
        }

        if ( i < j )
        {
            SCAI_LOG_TRACE( logger, "swap i = " << i << ", j = " << j << ", " << keys[i] << ", " << keys[j] )
            std::swap( keys[i], keys[j] );
            std::swap( values[i], values[j] );
            i++;
            j--;
        } 
        else if ( i == j )
        {
            // no swap, just to finish the loop with the right partition index
            i++;
        }
    }

    SCAI_LOG_TRACE( logger, "partition index i = " << i )

    if ( left < i - 1 )
    {
        qsort( keys, values, left, i - 1, ascending );
    }

    if ( i < right )
    {
        qsort( keys, values, i, right, ascending );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPUtils::sortInPlace(
    IndexType indexes[],
    ValueType values[],
    const IndexType n,
    const bool ascending )
{
    SCAI_REGION( "OpenMP.Utils.sortInPlace" )

    if ( n >= 2 )
    {
        qsort( indexes, values, 0, n - 1, ascending );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType OpenMPUtils::countNonZeros( const ValueType denseArray[], const IndexType n, const ValueType zero, ValueType eps )
{
    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    SCAI_REGION( "OpenMP.Utils.countNZ" )

    IndexType nonZeros = 0;

    RealType absEps = eps;

    #pragma omp parallel for reduction( +:nonZeros )

    for ( IndexType i = 0; i < n; ++i )
    {
        RealType absVal = common::Math::abs( denseArray[i] - zero );

        if ( absVal > absEps )
        {
            nonZeros++;
        }
    }

    SCAI_LOG_INFO( logger, "countNonZeros<" << TypeTraits<ValueType>::id() << "> of array[" << n << "] -> " << nonZeros )

    return nonZeros;
}

/* --------------------------------------------------------------------------- */

template<typename TargetType, typename SourceType>
IndexType OpenMPUtils::compress(
    TargetType sparseArray[],
    IndexType sparseIndexes[],
    const SourceType denseArray[],
    const IndexType n,
    const SourceType zero,
    const SourceType eps )
{
    typedef typename common::TypeTraits<SourceType>::RealType RealType;

    SCAI_REGION( "OpenMP.Utils.compress" )

    IndexType nonZeros = 0;

    // use of parallel for + atomicInc might be possible but would give an arbitrary order of sparse indexes

    RealType absEps = eps;

    for ( IndexType i = 0; i < n; ++i )
    {
        RealType sourceVal = common::Math::abs( denseArray[i] - zero );

        if ( sourceVal > absEps )
        {
            IndexType k = nonZeros++;  // parallel: atomicInc( nonZeros );

            if ( sparseArray )
            {
                sparseArray[k] = static_cast<TargetType>( denseArray[i] );
            }

            if ( sparseIndexes )
            {
                sparseIndexes[k] = i;
            }
        }
    }

    return nonZeros;
}

/* --------------------------------------------------------------------------- */

template<typename BucketType>
void OpenMPUtils::countBuckets(
    IndexType bucketSizes[],
    const BucketType nBuckets,
    const BucketType bucketMap[],
    const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.countBuckets" )

    SCAI_LOG_INFO( logger, "countBuckets, #elems = " << n << ", #buckets = " << nBuckets )

    // initialize size array for each bucket

    for ( BucketType i = 0; i < nBuckets; i++ )
    {
        bucketSizes[i] = 0;
    }

    // increment size of a bucket for each value mapped to this bucket
    // OpenMP parallelization requires atomicInc to be safe

    #pragma omp parallel for

    for ( IndexType k = 0; k < n; k++ )
    {
        BucketType iBucket = bucketMap[k];

        // No error message here, just count

        if ( scai::common::Utils::validIndex( iBucket, nBuckets ) )
        {
            atomicInc( bucketSizes[iBucket] );
        }
    }

    // Note: out of range values if bucketSizes.sum() != n
}

/* --------------------------------------------------------------------------- */

template<typename BucketType>
void OpenMPUtils::sortInBuckets( IndexType sortedIndexes[],
                                 IndexType offsets[],           // used as tmp, remains unchanged
                                 const BucketType nBuckets,
                                 const BucketType bucketMap[],
                                 const IndexType n )
{
    SCAI_REGION( "OpenMP.Utils.sortInBuckets" )

    SCAI_LOG_INFO( logger, "sortInBuckets, #elems = " << n << ", #buckets = " << nBuckets )

    // OpenMP parallelization not helpful as it destroys stable sorting

    for ( IndexType k = 0; k < n; k++ )
    {
        BucketType iBucket = bucketMap[k];

        if ( !scai::common::Utils::validIndex( iBucket, nBuckets ) )
        {
            continue;
        }

        IndexType& offset = offsets[iBucket];
  
        SCAI_ASSERT_LT_DEBUG( offset, n, "out of range offset" )
        sortedIndexes[offset] = k;
        offset++;
        SCAI_LOG_TRACE( logger, k << " is in bucket " << iBucket << " of " << nBuckets << ", offset = " << offsets[iBucket] )
    }

    // set back the old offsets

    for ( IndexType i = nBuckets; i > 0; --i )
    {
        offsets[i] = offsets[i - 1];
    }

    offsets[0] = 0;
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPUtils::countAddSparse(
    const IndexType indexes1[],
    const IndexType n1,
    const IndexType indexes2[],
    const IndexType n2 )
{
    SCAI_LOG_DEBUG( logger, "countAddSparse( n1 = " << n1 << ", n2 = " << n2 << ")" )

    IndexType n = 0;

    IndexType i1 = 0;
    IndexType i2 = 0;

    // merge via the sorted indexes

    while ( i1 < n1 && i2 < n2 )
    {
        if ( indexes1[i1] == indexes2[i2] )
        {
            ++n;
            ++i1;
            ++i2;
        }
        else if ( indexes1[i1] < indexes2[i2] )
        {
            ++n;
            ++i1;
        }
        else
        {
            ++n;
            ++i2;
        }
    }

    // add remaining indexes from array1 or array2

    n += ( n1 - i1 ) + ( n2 - i2 );

    return n;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType OpenMPUtils::addSparse(
    IndexType indexes[],
    ValueType values[],
    const IndexType indexes1[],
    const ValueType values1[],
    const ValueType zero1,
    const IndexType n1,
    const ValueType alpha,
    const IndexType indexes2[],
    const ValueType values2[],
    const ValueType zero2,
    const IndexType n2,
    const ValueType beta )
{
    SCAI_REGION( "OpenMP.Utils.addSparse" )

    SCAI_LOG_DEBUG( logger, "addSparse: "<< alpha << "* Sp( n1 = " << n1 << ", zero1 = " << zero1 << "), "
                             << beta << " * Sp( n2 = " << n2 << ", zero2 = " << zero2 )

    IndexType n = 0;

    IndexType i1 = 0;
    IndexType i2 = 0;

    // merge via the sorted indexes

    while ( i1 < n1 && i2 < n2 )
    {
        if ( indexes1[i1] == indexes2[i2] )
        {
            // entry at same position

            indexes[n] = indexes1[i1];
            values[n]  = alpha * values1[i1] + beta * values2[i2];
            ++n;
            ++i1;
            ++i2;
        }
        else if ( indexes1[i1] < indexes2[i2] )
        {
            // entry only in array1

            indexes[n] = indexes1[i1];
            values[n]  = alpha * values1[i1] + beta * zero2;
            ++n;
            ++i1;
        }
        else
        {
            // entry only in array2

            indexes[n] = indexes2[i2];
            values[n]  = beta * values2[i2] + alpha * zero1;
            ++n;
            ++i2;
        }
    }

    while ( i1 < n1 )
    {
        indexes[n] = indexes1[i1];
        values[n]  = alpha * values1[i1] + beta * zero2;
        ++n;
        ++i1;
    }

    while ( i2 < n2 )
    {
        indexes[n] = indexes2[i2];
        values[n]  = beta * values2[i2] + alpha * zero1;
        ++n;
        ++i2;
    }

    return n;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType OpenMPUtils::joinSparse(
    IndexType indexes[],
    ValueType values[],
    const IndexType indexes1[],
    const ValueType values1[],
    const IndexType n1,
    const IndexType indexes2[],
    const ValueType values2[],
    const IndexType n2 )
{
    SCAI_REGION( "OpenMP.Utils.joinSparse" )

    SCAI_LOG_DEBUG( logger, "joinSparse: Sp( n1 = " << n1 << " ), "
                             << " with Sp( n2 = " << n2 << " )" )

    IndexType n = 0;

    IndexType i1 = 0;
    IndexType i2 = 0;

    // join merge via the sorted indexes

    while ( i1 < n1 && i2 < n2 )
    {
        if ( indexes1[i1] == indexes2[i2] )
        {
            // entry at same position, take the new one

            indexes[n] = indexes1[i1];
            values[n]  = values2[i2];
            ++n;
            ++i1;
            ++i2;
        }
        else if ( indexes1[i1] < indexes2[i2] )
        {
            // entry only in array1

            indexes[n] = indexes1[i1];
            values[n]  = values1[i1];
            ++n;
            ++i1;
        }
        else
        {
            // entry only in array2

            indexes[n] = indexes2[i2];
            values[n]  = values2[i2];
            ++n;
            ++i2;
        }
    }

    while ( i1 < n1 )
    {
        indexes[n] = indexes1[i1];
        values[n]  = values1[i1];
        ++n;
        ++i1;
    }

    while ( i2 < n2 )
    {
        indexes[n] = indexes2[i2];
        values[n]  = values2[i2];
        ++n;
        ++i2;
    }

    return n;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType OpenMPUtils::binopSparse(
    IndexType indexes[],
    ValueType values[],
    const IndexType indexes1[],
    const ValueType values1[],
    const ValueType zero1,
    const IndexType n1,
    const IndexType indexes2[],
    const ValueType values2[],
    const ValueType zero2,
    const IndexType n2,
    const BinaryOp op )
{
    if ( op == common::BinaryOp::COPY )
    {
        return joinSparse( indexes, values, indexes1, values1, n1, indexes2, values2, n2 );
    }

    SCAI_REGION( "OpenMP.Utils.binopSparse" )

    SCAI_LOG_DEBUG( logger, "binopSparse: Sp( n1 = " << n1 << ", zero1 = " << zero1 << "), "
                             << op << " Sp( n2 = " << n2 << ", zero2 = " << zero2 )

    IndexType n = 0;

    IndexType i1 = 0;
    IndexType i2 = 0;

    // merge via the sorted indexes

    while ( i1 < n1 && i2 < n2 )
    {
        if ( indexes1[i1] == indexes2[i2] )
        {
            // entry at same position

            indexes[n] = indexes1[i1];
            values[n]  = applyBinary( values1[i1], op, values2[i2] );
            ++n;
            ++i1;
            ++i2;
        }
        else if ( indexes1[i1] < indexes2[i2] )
        {
            // entry only in array1

            indexes[n] = indexes1[i1];
            values[n]  = applyBinary( values1[i1], op, zero2 );
            ++n;
            ++i1;
        }
        else
        {
            // entry only in array2

            indexes[n] = indexes2[i2];
            values[n]  = applyBinary( zero1, op, values2[i2] );
            ++n;
            ++i2;
        }
    }

    while ( i1 < n1 )
    {
        indexes[n] = indexes1[i1];
        values[n]  = applyBinary( values1[i1], op, zero2 );
        ++n;
        ++i1;
    }

    while ( i2 < n2 )
    {
        indexes[n] = indexes2[i2];
        values[n]  = applyBinary( zero1, op, values2[i2] );
        ++n;
        ++i2;
    }

    return n;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType OpenMPUtils::allCompareSparse(
    bool& allFlag,
    const IndexType indexes1[],
    const ValueType values1[],
    const ValueType zero1,
    const IndexType n1,
    const IndexType indexes2[],
    const ValueType values2[],
    const ValueType zero2,
    const IndexType n2,
    const CompareOp op )
{
    SCAI_REGION( "OpenMP.Utils.binopSparse" )

    SCAI_LOG_DEBUG( logger, "binopSparse: Sp( n1 = " << n1 << ", zero1 = " << zero1 << "), "
                             << op << " Sp( n2 = " << n2 << ", zero2 = " << zero2 )

    IndexType i1 = 0;
    IndexType i2 = 0;
    IndexType n  = 0;

    allFlag = true;

    // merge via the sorted indexes

    while ( i1 < n1 && i2 < n2 )
    {
        if ( indexes1[i1] == indexes2[i2] )
        {
            // entry at same position

            allFlag = allFlag && compare( values1[i1], op, values2[i2] );
            ++i1;
            ++i2;
            ++n;
        }
        else if ( indexes1[i1] < indexes2[i2] )
        {
            // entry only in array1

            allFlag = allFlag && compare( values1[i1], op, zero2 );
            ++i1;
            ++n;
        }
        else
        {
            // entry only in array2

            allFlag = allFlag && compare( zero1, op, values2[i2] );
            ++i2;
            ++n;
        }
    }

    while ( i1 < n1 )
    {
        allFlag = allFlag && compare( values1[i1], op, zero2 );
        ++i1;
        ++n;
    }

    while ( i2 < n2 )
    {
        allFlag = allFlag && compare( zero1, op, values2[i2] );
        ++i2;
        ++n;
    }

    return n;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
static inline
void appendSparse( IndexType indexes[], ValueType values[], IndexType& n,
                   const IndexType ia, const ValueType v, const BinaryOp op )
{
    if ( n > 0 && ( indexes[n - 1] == ia ) )
    {
        values[n - 1]  = applyBinary( values[n - 1], op, v );
    }
    else
    {
        values[n]  = v;
        indexes[n] = ia;
        n++;
    }
}

template<typename ValueType>
IndexType OpenMPUtils::mergeSparse(
    IndexType indexes[],
    ValueType values[],
    const IndexType indexes1[],
    const ValueType values1[],
    const IndexType n1,
    const IndexType indexes2[],
    const ValueType values2[],
    const IndexType n2,
    const BinaryOp op )
{
    SCAI_REGION( "OpenMP.Utils.mergeSparse" )

    SCAI_LOG_DEBUG( logger, "mergeSparse: Sp( n1 = " << n1 << " ), "
                             << " with Sp( n2 = " << n2 << " )" )

    IndexType n = 0;

    IndexType i1 = 0;
    IndexType i2 = 0;

    // merge via the sorted indexes

    while ( i1 < n1 && i2 < n2 )
    {
        if ( indexes1[i1] <= indexes2[i2] )
        {
            appendSparse( indexes, values, n, indexes1[i1], values1[i1], op );
            ++i1;
        }
        else 
        {
            appendSparse( indexes, values, n, indexes2[i2], values2[i2], op );
            ++i2;
        }
    }

    while ( i1 < n1 )
    {
        appendSparse( indexes, values, n, indexes1[i1], values1[i1], op );
        ++i1;
    }

    while ( i2 < n2 )
    {
        appendSparse( indexes, values, n, indexes2[i2], values2[i2], op );
        ++i2;
    }

    return n;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPUtils::BaseKernels::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_INFO( logger, "register UtilsKernel OpenMP-routines for Host at kernel registry [" << flag << "]" )
    // we keep the registrations for IndexType as we do not need conversions
    KernelRegistry::set<UtilKernelTrait::validIndexes>( validIndexes, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::countBuckets<IndexType> >( countBuckets, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::sortInBuckets<IndexType> >( sortInBuckets, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setInversePerm>( setInversePerm, ctx, flag );
    KernelRegistry::set<SparseKernelTrait::countAddSparse>( countAddSparse, ctx, flag );
}

template<typename ValueType>
void OpenMPUtils::ArrayKernels<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register UtilsKernel OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    // we keep the registrations for IndexType as we do not need conversions
    KernelRegistry::set<UtilKernelTrait::reduce<ValueType> >( reduce, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::reduce2<ValueType> >( reduce2, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::allCompare<ValueType> >( allCompare, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::allCompareScalar<ValueType> >( allCompareScalar, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setOrder<ValueType> >( setOrder, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setSequence<ValueType> >( setSequence, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::getValue<ValueType> >( getValue, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setVal<ValueType> >( setVal, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::scaleVectorAddScalar<ValueType> >( scaleVectorAddScalar, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::isSorted<ValueType> >( isSorted, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::scan<ValueType> >( scan, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::unscan<ValueType> >( unscan, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::sort<ValueType> >( sort, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::sortInPlace<ValueType> >( sortInPlace, ctx, flag );
    KernelRegistry::set<SparseKernelTrait::countNonZeros<ValueType> >( countNonZeros, ctx, flag );
    KernelRegistry::set<SparseKernelTrait::addSparse<ValueType> >( addSparse, ctx, flag );
    KernelRegistry::set<SparseKernelTrait::binopSparse<ValueType> >( binopSparse, ctx, flag );
    KernelRegistry::set<SparseKernelTrait::allCompareSparse<ValueType> >( allCompareSparse, ctx, flag );
    KernelRegistry::set<SparseKernelTrait::mergeSparse<ValueType> >( mergeSparse, ctx, flag );

    KernelRegistry::set<UtilKernelTrait::unaryOp<ValueType> >( unaryOp, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::binaryOp<ValueType> >( binaryOp, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::binaryOpScalar<ValueType> >( binaryOpScalar, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::scatterVal<ValueType> >( scatterVal, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::fillSection<ValueType> >( fillSection, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void OpenMPUtils::BinOpKernels<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register UtilsKernel OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )
    KernelRegistry::set<UtilKernelTrait::setGather<ValueType, OtherValueType> >( setGather, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setGatherSparse<ValueType, OtherValueType> >( setGatherSparse, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setScatter<ValueType, OtherValueType> >( setScatter, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::set<ValueType, OtherValueType> >( set, ctx, flag );
    KernelRegistry::set<UtilKernelTrait::setSection<ValueType, OtherValueType> >( setSection, ctx, flag );
    KernelRegistry::set<SparseKernelTrait::compress<ValueType, OtherValueType> >( compress, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPUtils::OpenMPUtils()
{
    SCAI_LOG_INFO( logger, "register UtilsKernel OpenMP-routines for Host" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    BaseKernels::registerKernels( flag );
    kregistry::mepr::RegistratorV<ArrayKernels, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<BinOpKernels, SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels( flag );
}

OpenMPUtils::~OpenMPUtils()
{
    SCAI_LOG_INFO( logger, "unregister UtilsKernel OpenMP-routines for Host" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    BaseKernels::registerKernels( flag );
    kregistry::mepr::RegistratorV<ArrayKernels, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<BinOpKernels, SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPUtils OpenMPUtils::guard;

} /* end namespace utilskernel */

} /* end namespace scai */
