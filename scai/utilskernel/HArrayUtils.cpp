/**
 * @file HArrayUtils.cpp
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
 * @brief Implemenation of the utility routines for HArrays using corresponding kernel routines.
 * @author Thomas Brandes
 * @date 10.10.2011
 */

// hpp
#include <scai/utilskernel/HArrayUtils.hpp>

// local library
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/utilskernel/mepr/UtilsWrapper.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/ScalarType.hpp>

#include <scai/common/macros/loop.hpp>

// std
#include <iostream>

namespace scai
{

using namespace hmemo;

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( HArrayUtils::logger, "HArrayUtils" )

/* --------------------------------------------------------------------------- */

void HArrayUtils::assign( _HArray& target, const _HArray& source, const ContextPtr prefLoc )
{
    if ( &target == &source )
    {
        SCAI_LOG_INFO( logger, "assign skipped, target and source are aliased" )
    }
    else
    {
        setArray( target, source, binary::COPY, prefLoc );
    }
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::setArray(
    _HArray& target,
    const _HArray& source,
    const binary::BinaryOp op,
    const ContextPtr prefLoc )
{
    mepr::UtilsWrapperTT1<SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::setArray( target, source, op, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::setArrayImpl(
    HArray<TargetValueType>& target,
    const HArray<SourceValueType>& source,
    const binary::BinaryOp op,
    const ContextPtr prefLoc )
{
    SCAI_ASSERT_ERROR( isBinarySupported<TargetValueType>( op ),
                       op << " not supported for " << common::TypeTraits<TargetValueType>::id() )

    // set should be available on interface for each loc

    static LAMAKernel<UtilKernelTrait::set<TargetValueType, SourceValueType> > set;

    ContextPtr loc = prefLoc;

    if ( !loc )
    {
        if ( op == binary::COPY )
        {
            // if no context is given we assign where source has a valid copy available
            loc = source.getValidContext();
        }
        else
        {
            // if no context is given we assign where target has a valid copy available
            loc = target.getValidContext();
        }
    }

    set.getSupportedContext( loc );

    const IndexType n = source.size();

    SCAI_CONTEXT_ACCESS( loc )

    if ( op == binary::COPY )
    {
        ReadAccess<SourceValueType> sourceVals( source, loc );
        WriteOnlyAccess<TargetValueType> targetVals( target, loc, n );
        set[loc]( targetVals.get(), sourceVals.get(), n, op );
    }
    else
    {
        ReadAccess<SourceValueType> sourceVals( source, loc );
        WriteAccess<TargetValueType> targetVals( target, loc );
        set[loc]( targetVals.get(), sourceVals.get(), n, op );
    }
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::setArraySection(
    _HArray& target,
    const IndexType targetOffset,
    const IndexType targetStride,
    const _HArray& source,
    const IndexType sourceOffset,
    const IndexType sourceStride,
    const IndexType n,
    const binary::BinaryOp op,
    const ContextPtr prefLoc )
{
    mepr::UtilsWrapperTT1<SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::setArraySection(
        target, targetOffset, targetStride,
        source, sourceOffset, sourceStride,
        n, op, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::setArraySectionImpl(
    HArray<TargetValueType>& target,
    const IndexType targetOffset,
    const IndexType targetStride,
    const HArray<SourceValueType>& source,
    const IndexType sourceOffset,
    const IndexType sourceStride,
    const IndexType n,
    const binary::BinaryOp op,
    const ContextPtr prefLoc )
{
    SCAI_ASSERT_ERROR( isBinarySupported<TargetValueType>( op ),
                       op << " not supported for " << common::TypeTraits<TargetValueType>::id() )

    // in contrary to setArray we assume correctly allocated target

    static LAMAKernel<UtilKernelTrait::setSection<TargetValueType, SourceValueType> > setSection;

    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        // optional argument has not been set, take default

        loc = source.getValidContext();
    }

    setSection.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<SourceValueType> sourceVals( source, loc );
    WriteAccess<TargetValueType> targetVals( target, loc );

    setSection[loc]( targetVals.get() + targetOffset, targetStride,
                     sourceVals.get() + sourceOffset, sourceStride,
                     n, op );
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::gather(
    _HArray& target,
    const _HArray& source,
    const HArray<IndexType>& indexes,
    const binary::BinaryOp op,
    const ContextPtr prefLoc )
{
    // use metaprogramming to call the gather version with the correct value types for target and source
    mepr::UtilsWrapperTT1<SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::gather( target, source, indexes, op, prefLoc );
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::scatter(
    _HArray& target,
    const HArray<IndexType>& indexes,
    const _HArray& source,
    const binary::BinaryOp op,
    const ContextPtr prefLoc )
{
    // use metaprogramming to call the scatter version with the correct value types for target and source
    mepr::UtilsWrapperTT1<SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::scatter( target, indexes, source, op, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::gatherImpl(
    HArray<TargetValueType>& target,
    const HArray<SourceValueType>& source,
    const HArray<IndexType>& indexes,
    const binary::BinaryOp op,
    const ContextPtr prefLoc )
{
    SCAI_ASSERT_ERROR( isBinarySupported<TargetValueType>( op ),
                       op << " not supported for " << common::TypeTraits<TargetValueType>::id() )

    SCAI_REGION( "HArray.gather" )

    SCAI_LOG_INFO( logger, "target[ " << target.size() << " ] " << op << " = source [ indexes [ "
                   << indexes.size() << " ] : " << source.size() << " ]" )

    // choose location for the operation where source array is currently valid
    static LAMAKernel<UtilKernelTrait::setGather<TargetValueType, SourceValueType> > setGather;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = source.getValidContext();
    }

    setGather.getSupportedContext( loc );
    const IndexType n = indexes.size();

    if ( op == binary::COPY )
    {
        //  target[i] = source[ indexes[i] ]

        WriteOnlyAccess<TargetValueType> wTarget( target, loc, n );
        SCAI_CONTEXT_ACCESS( loc )
        ReadAccess<SourceValueType> rSource( source, loc );
        ReadAccess<IndexType> rIndexes( indexes, loc );
        setGather[loc] ( wTarget.get(), rSource.get(), rIndexes.get(), op, n );
    }
    else
    {
        //  target[i] op= source[ indexes[i] ]

        SCAI_ASSERT_GE_ERROR( target.size(), n, "target array too small" )

        WriteAccess<TargetValueType> wTarget( target, loc );
        SCAI_CONTEXT_ACCESS( loc )
        ReadAccess<SourceValueType> rSource( source, loc );
        ReadAccess<IndexType> rIndexes( indexes, loc );
        setGather[loc] ( wTarget.get(), rSource.get(), rIndexes.get(), op, n );
    }
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::scatterImpl(
    HArray<TargetValueType>& target,
    const HArray<IndexType>& indexes,
    const HArray<SourceValueType>& source,
    const binary::BinaryOp op,
    const ContextPtr prefLoc )
{
    SCAI_ASSERT_ERROR( isBinarySupported<TargetValueType>( op ),
                       op << " not supported for " << common::TypeTraits<TargetValueType>::id() )

    SCAI_REGION( "HArray.scatter" )

    SCAI_LOG_INFO( logger, "target[ indexes[ " << indexes.size() << " ] : " << target.size()
                   << " ] = source[ " << source.size() << " ]" )

    SCAI_ASSERT(      ( op == binary::COPY )
                      || ( op == binary::ADD )
                      || ( op == binary::SUB ) , "Unsupported reduction op " << op  )

    SCAI_ASSERT( HArrayUtils::validIndexes( indexes, target.size(), prefLoc ),
                 "illegal scatter index, target has size " << target.size() )
    static LAMAKernel<UtilKernelTrait::setScatter<TargetValueType, SourceValueType> > setScatter;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = source.getValidContext();
    }

    setScatter.getSupportedContext( loc );
    const IndexType n = indexes.size();
    WriteAccess<TargetValueType> wTarget( target, loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<SourceValueType> rSource( source, loc );
    ReadAccess<IndexType> rIndexes( indexes, loc );
    //  target[ indexes[i] ] = source[i]
    setScatter[loc] ( wTarget.get(), rIndexes.get(), rSource.get(), op, n );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::assignScalar(
    _HArray& target,
    const ValueType value,
    const binary::BinaryOp op,
    ContextPtr prefLoc )
{
    mepr::UtilsWrapperT<ValueType, SCAI_ARRAY_TYPES_HOST_LIST>::setScalar( target, value, op, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::setScalar(
    HArray<ValueType>& target,
    const ValueType value,
    const binary::BinaryOp op,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_ERROR( isBinarySupported<ValueType>( op ),
                       op << " not supported for " << common::TypeTraits<ValueType>::id() )

    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        // default location: where we have valid copy of the data
        loc = target.getValidContext();
    }

    setVal.getSupportedContext( loc );
    const IndexType n = target.size();
    SCAI_LOG_INFO( logger, target << " = " << value << ", to do at " << *loc << ", n = " << n )

    if ( op == binary::COPY )
    {
        // Note: very important is to specify the size n here as it might not have been allocated
        WriteOnlyAccess<ValueType> wTarget( target, loc, n );
        SCAI_CONTEXT_ACCESS( loc )
        setVal[loc]( wTarget.get(), n, value, op );
    }
    else
    {
        WriteAccess<ValueType> wTarget( target, loc );
        SCAI_CONTEXT_ACCESS( loc )
        setVal[loc]( wTarget.get(), n, value, op );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::setVal(
    _HArray& target,
    const IndexType index,
    const ValueType val )
{
    SCAI_ASSERT_DEBUG( index < target.size(), "index = " << index << " out of range for target = " << target );
    mepr::UtilsWrapperT< ValueType, SCAI_ARRAY_TYPES_HOST_LIST>::setValImpl( target, index, val );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::setValImpl(
    HArray<ValueType>& target,
    const IndexType index,
    const ValueType val,
    const binary::BinaryOp op )
{
    SCAI_ASSERT_ERROR( isBinarySupported<ValueType>( op ),
                       op << " not supported for " << common::TypeTraits<ValueType>::id() )

    // setting single value will directly copy to the device with the valid incarnation

    ContextPtr loc = target.getValidContext();

    if ( op == binary::COPY )
    {
        WriteAccess<ValueType> wTarget( target, loc );
        wTarget.setValue( val, index );
    }
    else
    {
        static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;
        setVal.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        WriteAccess<ValueType> wTarget( target, loc );
        setVal[loc]( wTarget.get() + index, 1, val, op );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::getVal(
    const _HArray& array,
    const IndexType index )
{
    ValueType val = mepr::UtilsWrapperT< ValueType, SCAI_ARRAY_TYPES_HOST_LIST>::getValImpl( array, index );
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::getValImpl(
    const HArray<ValueType>& array,
    const IndexType index )
{
    SCAI_ASSERT_DEBUG( index < array.size(), "index = " << index << " out of range for array = " << array );
    // get the data from a valid context, avoids any memory copy.
    ContextPtr loc = array.getValidContext();
    ReadAccess<ValueType> rArray( array, loc );
    ValueType val;
    rArray.getValue( val, index );
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::reduce(
    const HArray<ValueType>& array,
    const binary::BinaryOp redOp,
    const ContextPtr prefLoc )
{
    static LAMAKernel<UtilKernelTrait::reduce<ValueType> > reduce;
    // preferred location: where valid values of the array are available
    ContextPtr loc = array.getValidContext( prefLoc );
    reduce.getSupportedContext( loc );
    ReadAccess<ValueType> readArray( array, loc );
    SCAI_CONTEXT_ACCESS( loc )
    ValueType zero = zeroBinary<ValueType>( redOp );
    ValueType redVal = reduce[loc]( readArray.get(), readArray.size(), zero, redOp );
    return redVal;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::asum( const HArray<ValueType>& array, const ContextPtr prefLoc )
{
    const IndexType n = array.size();

    if ( n == 0 )
    {
        return ValueType( 0 );
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::asum<ValueType> > asum;
    // preferred location: where valid values of the array are available
    ContextPtr loc = array.getValidContext( prefLoc );
    asum.getSupportedContext( loc );
    ReadAccess<ValueType> readArray( array, loc );
    SCAI_CONTEXT_ACCESS( loc )
    ValueType result = asum[loc]( n, readArray.get(), 1 );
    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::nrm2( const HArray<ValueType>& array, const ContextPtr prefLoc )
{
    const IndexType n = array.size();

    if ( n == 0 )
    {
        return ValueType( 0 );
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::nrm2<ValueType> > nrm2;
    // preferred location: where valid values of the array are available
    ContextPtr loc = array.getValidContext( prefLoc );
    nrm2.getSupportedContext( loc );
    ReadAccess<ValueType> readArray( array, loc );
    SCAI_CONTEXT_ACCESS( loc )
    ValueType result = nrm2[loc]( n, readArray.get(), 1 );
    return result;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::absMaxDiffVal(
    const HArray<ValueType>& array1,
    const HArray<ValueType>& array2,
    ContextPtr prefLoc )
{
    // is just a special case of a reduction with two combined arrays

    return reduce2( array1, array2, binary::SUB, binary::ABS_MAX, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::reduce2(
    const HArray<ValueType>& array1,
    const HArray<ValueType>& array2,
    const binary::BinaryOp binOp,
    const binary::BinaryOp redOp,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_ERROR( isBinarySupported<ValueType>( binOp ),
                       binOp << " not supported for " << common::TypeTraits<ValueType>::id() )

    ValueType zero = zeroBinary<ValueType>( redOp ); // checks also for valid op

    SCAI_ASSERT_EQUAL( array1.size(), array2.size(), "array size mismatch for reduce2" )

    const IndexType n = array1.size();

    static LAMAKernel<UtilKernelTrait::reduce2<ValueType> > reduce2;

    ContextPtr loc = prefLoc;

    // Rule for default location: where array1 has valid values

    if ( loc == ContextPtr() )
    {
        loc = array1.getValidContext();
    }

    reduce2.getSupportedContext( loc );

    ReadAccess<ValueType> readArray1( array1, loc );
    ReadAccess<ValueType> readArray2( array2, loc );

    SCAI_CONTEXT_ACCESS( loc )

    ValueType redVal = reduce2[loc]( readArray1.get(),
                                     readArray2.get(),
                                     n,
                                     binOp,
                                     zero,
                                     redOp );
    return redVal;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::dotProduct(
    const HArray<ValueType>& array1,
    const HArray<ValueType>& array2,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQUAL( array1.size(), array2.size(), "array size mismatch for dotproduct" )
    static LAMAKernel<blaskernel::BLASKernelTrait::dot<ValueType> > dot;
    ContextPtr loc = prefLoc;

    // Rule for default location: where array1 has valid values

    if ( loc == ContextPtr() )
    {
        loc = array1.getValidContext();
    }

    dot.getSupportedContext( loc );
    ReadAccess<ValueType> readArray1( array1, loc );
    ReadAccess<ValueType> readArray2( array2, loc );
    SCAI_CONTEXT_ACCESS( loc )
    const ValueType res = dot[loc]( readArray1.size(), readArray1.get(), 1, readArray2.get(), 1 );
    return res;

    // Note: other solution: reduce2( array1, array2, n, binary::MULT, 0, binary::ADD )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::axpy(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    ContextPtr prefLoc )
{
    SCAI_LOG_INFO( logger, "result += " << alpha << " * x"
                   ", x = " << x << ", result = " << result )

    if ( alpha == common::constants::ZERO )
    {
        return;
    }

    SCAI_ASSERT_EQUAL( result.size(), x.size(), "size mismatch" )
    IndexType n = x.size();
    static LAMAKernel<blaskernel::BLASKernelTrait::axpy<ValueType> > axpy;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = result.getValidContext();
    }

    axpy.getSupportedContext( loc );  // axpy must be available at loc
    ReadAccess<ValueType> xAccess( x, loc );
    WriteAccess<ValueType> resultAccess( result, loc, true );
    SCAI_CONTEXT_ACCESS( loc )
    axpy[loc]( n, alpha, xAccess.get(), 1, resultAccess.get(), 1 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::arrayPlusArray(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    ContextPtr prefLoc )
{
    SCAI_LOG_INFO( logger, "result = " << alpha << " * x + " << beta << " * y" <<
                   ", x = " << x << ", y = " << y << ", result = " << result )

    // check for zero terms as we do not need read access and assert correct sizes

    if ( beta == common::constants::ZERO )
    {
        if ( alpha == common::constants::ZERO )
        {
            setScalar( result, ValueType( 0 ), binary::COPY, prefLoc );
        }
        else
        {
            binaryOpScalar1( result, alpha, x, binary::MULT, prefLoc );
        }

        return;
    }

    if ( alpha == common::constants::ZERO )
    {
        binaryOpScalar1( result, beta, y, binary::MULT, prefLoc );
        return;
    }

    if ( &y == &result && beta == common::constants::ONE )
    {
        axpy( result, alpha, x, prefLoc );
        return;
    }

    if ( &x == &result && alpha == common::constants::ONE )
    {
        axpy( result, beta, y, prefLoc );
        return;
    }

    SCAI_ASSERT_EQUAL( x.size(), y.size(), "size mismatch" )
    const IndexType n = x.size();
    static LAMAKernel<blaskernel::BLASKernelTrait::sum<ValueType> > sum;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = x.getValidContext();
    }

    sum.getSupportedContext( loc );
    // no alias checks are done here
    ReadAccess<ValueType> xAccess( x, loc );
    ReadAccess<ValueType> yAccess( y, loc );
    WriteOnlyAccess<ValueType> resultAccess( result, loc, x.size() );
    SCAI_CONTEXT_ACCESS( loc )
    sum[loc]( n, alpha, xAccess.get(), beta, yAccess.get(), resultAccess.get() );
}

template<typename ValueType>
void HArrayUtils::arrayPlusScalar(
    hmemo::HArray<ValueType>& result,
    const ValueType alpha,
    const hmemo::HArray<ValueType>& x,
    const ValueType beta,
    hmemo::ContextPtr prefLoc )
{
    // operation is result = alpha * x + b (elementwise)
    SCAI_LOG_INFO( logger, "result = " << alpha << " * x + " << beta << ", x = " << x
                   << ", result = " << result )

    // check for zero terms as we do not need read access and assert correct sizes

    if ( alpha == common::constants::ZERO ) // result = b
    {
        result.resize( x.size() );
        assignScalar( result, beta, binary::COPY, prefLoc );
        return;
    }

    if ( &x == &result && alpha == common::constants::ONE ) // result += b (elementwise)
    {
        assignScalar( result, beta, binary::ADD, prefLoc );
        return;
    }

    if ( alpha == common::constants::ONE ) // result = x + b (elementwise)
    {
        binaryOpScalar1( result, beta, x, binary::ADD, prefLoc );
        return;
    }

    if ( beta == common::constants::ZERO )
    {
        if ( alpha == common::constants::ZERO ) // result = 0
        {
            setScalar( result, ValueType( 0 ), binary::COPY, prefLoc );
        }
        else // result = alpha * x
        {
            binaryOpScalar1( result, alpha, x, binary::MULT, prefLoc );
        }

        return;
    }

    // else full operation: result = alpha * x + b (elementwise)

    const IndexType n = x.size();

    static LAMAKernel<UtilKernelTrait::scaleVectorAddScalar<ValueType> > scaleVectorAddScalar;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = x.getValidContext();
    }

    scaleVectorAddScalar.getSupportedContext( loc );

    // due to possible alias of result and x, write access must follow read

    ReadAccess<ValueType> rX( x, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, n );

    SCAI_CONTEXT_ACCESS( loc )

    scaleVectorAddScalar[loc]( wResult.get(), rX.get(), n, alpha, beta );
}

template<typename ValueType>
void HArrayUtils::arrayTimesArray(
    hmemo::HArray<ValueType>& result,
    const ValueType alpha,
    const hmemo::HArray<ValueType>& x,
    const hmemo::HArray<ValueType>& y,
    hmemo::ContextPtr prefLoc )
{
    SCAI_ASSERT_EQUAL( x.size(), y.size(), "size mismatch" )
    const IndexType n = x.size();

    if ( alpha == common::constants::ZERO )
    {
        setScalar( result, ValueType( 0 ), binary::COPY, prefLoc );
        return;
    }

    static LAMAKernel<UtilKernelTrait::binaryOp<ValueType> > binaryOp;

    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = x.getValidContext();
    }

    binaryOp.getSupportedContext( loc );
    // no alias checks are done here
    {
        ReadAccess<ValueType> xAccess( x, loc );
        ReadAccess<ValueType> yAccess( y, loc );
        WriteOnlyAccess<ValueType> resultAccess( result, loc, x.size() );
        SCAI_CONTEXT_ACCESS( loc )
        binaryOp[loc]( resultAccess.get(), xAccess.get(), yAccess.get(), n, binary::MULT );
    }

    if ( alpha != common::constants::ONE )
    {
        setScalar( result, alpha, binary::MULT, prefLoc );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::unaryOp(
    hmemo::HArray<ValueType>& result,
    const hmemo::HArray<ValueType>& x,
    const unary::UnaryOp op,
    ContextPtr prefLoc )
{
    SCAI_ASSERT( isUnarySupported<ValueType>( op ),
                 op << " not supported for " << common::TypeTraits<ValueType>::id() )

    const IndexType n = x.size();

    static LAMAKernel<UtilKernelTrait::unaryOp<ValueType> > unaryOpKernel;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = x.getValidContext();
    }

    unaryOpKernel.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    // due to possible alias of result and x, write access must follow read

    ReadAccess<ValueType> rX( x, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, n );

    unaryOpKernel[loc]( wResult.get(), rX.get(), n, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::binaryOp(
    hmemo::HArray<ValueType>& result,
    const hmemo::HArray<ValueType>& x,
    const hmemo::HArray<ValueType>& y,
    const binary::BinaryOp op,
    ContextPtr prefLoc )
{
    SCAI_ASSERT( isBinarySupported<ValueType>( op ),
                 op << " not supported for " << common::TypeTraits<ValueType>::id() )

    const IndexType n = x.size();

    SCAI_ASSERT_GE_ERROR( y.size(), n, "binaryOp: y has to small size" )

    static LAMAKernel<UtilKernelTrait::binaryOp<ValueType> > binaryOp;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = x.getValidContext();
    }

    binaryOp.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    // due to possible alias of result and x, write access must follow read

    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, n );

    binaryOp[loc]( wResult.get(), rX.get(), rY.get(), n, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::binaryOpScalar1(
    hmemo::HArray<ValueType>& result,
    const ValueType x,
    const hmemo::HArray<ValueType>& y,
    const binary::BinaryOp op,
    ContextPtr prefLoc )
{
    SCAI_ASSERT( isBinarySupported<ValueType>( op ),
                 op << " not supported for " << common::TypeTraits<ValueType>::id() )

    const IndexType n = y.size();

    static LAMAKernel<UtilKernelTrait::binaryOpScalar1<ValueType> > binaryOpScalar1;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = y.getValidContext();
    }

    binaryOpScalar1.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    // due to possible alias of result and y, write access must follow read

    ReadAccess<ValueType> rY( y, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, n );

    binaryOpScalar1[loc]( wResult.get(), x, rY.get(), n, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::binaryOpScalar2(
    hmemo::HArray<ValueType>& result,
    const hmemo::HArray<ValueType>& x,
    const ValueType y,
    const binary::BinaryOp op,
    ContextPtr prefLoc )
{
    SCAI_ASSERT( isBinarySupported<ValueType>( op ),
                 op << " not supported for " << common::TypeTraits<ValueType>::id() )

    const IndexType n = x.size();

    static LAMAKernel<UtilKernelTrait::binaryOpScalar2<ValueType> > binaryOpScalar2;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = x.getValidContext();
    }

    binaryOpScalar2.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    // due to possible alias of result and x, write access must follow read

    ReadAccess<ValueType> rX( x, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, n );

    binaryOpScalar2[loc]( wResult.get(), rX.get(), y, n, op );
}

/* --------------------------------------------------------------------------- */

bool HArrayUtils::validIndexes(
    const hmemo::HArray<IndexType>& array,
    const IndexType size,
    const hmemo::ContextPtr prefLoc )
{
    const IndexType n = array.size();
    static LAMAKernel<UtilKernelTrait::validIndexes> validIndexes;
    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    validIndexes.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> rValues( array, loc );
    bool valid = validIndexes[loc]( rValues.get(), n, size );
    return valid;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool HArrayUtils::isSorted(
    const hmemo::HArray<ValueType>& array,
    const bool ascending,
    hmemo::ContextPtr prefLoc )
{
    const IndexType n = array.size();
    static LAMAKernel<UtilKernelTrait::isSorted<ValueType> > isSorted;
    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    isSorted.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<ValueType> rValues( array, loc );
    bool sorted = isSorted[loc]( rValues.get(), n, ascending );
    return sorted;
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::setOrder( hmemo::HArray<IndexType>& array, IndexType n, hmemo::ContextPtr prefLoc )
{
    static LAMAKernel<UtilKernelTrait::setOrder<IndexType> > setOrder;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();   // should be first context
    }

    setOrder.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteOnlyAccess<IndexType> wArray( array, loc, n );
    setOrder[loc]( wArray.get(), n );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::setSequence( hmemo::HArray<ValueType>& array, ValueType startValue, ValueType inc, IndexType n, hmemo::ContextPtr prefLoc )
{
    static LAMAKernel<UtilKernelTrait::setSequence<ValueType> > setSequence;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();   // should be first context
    }

    setSequence.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteOnlyAccess<ValueType> wArray( array, loc, n );
    setSequence[loc]( wArray.get(), startValue, inc, n );
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::setRandom( hmemo::_HArray& array,
                             const IndexType n,
                             const float fillRate,
                             const hmemo::ContextPtr prefLoc )
{
    // use meta-programming to call setRandomImpl<ValueType> with the type of array

    mepr::UtilsWrapper< SCAI_ARRAY_TYPES_HOST_LIST>::setRandom( array, n, fillRate, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::setRandomImpl( hmemo::HArray<ValueType>& array,
                                 const IndexType n,
                                 const float fillRate,
                                 const hmemo::ContextPtr prefLoc )
{
    ContextPtr loc = Context::getHostPtr();   // currently only available on host
    WriteOnlyAccess<ValueType> wArray( array, loc, n );

    if ( fillRate >= 1.0f )
    {
        for ( IndexType i = 0; i < n; ++i )
        {
            common::Math::random( wArray[i] );
        }
    }
    else
    {
        for ( IndexType i = 0; i < n; ++i )
        {
            float x = static_cast<float>( rand() ) / static_cast<float>( RAND_MAX );

            if ( x < fillRate )
            {
                common::Math::random( wArray[i] );
            }
            else
            {
                wArray[i] = ValueType( 0 );
            }
        }
    }

    if ( prefLoc != ContextPtr() )
    {
        array.prefetch( prefLoc );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::scan(
    hmemo::HArray<ValueType>& array,
    hmemo::ContextPtr prefLoc )
{
    const IndexType n = array.size();
    static LAMAKernel<UtilKernelTrait::scan<ValueType> > scan;
    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    scan.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteAccess<ValueType> wValues( array, loc );
    // One additional element will be added at the end
    wValues.resize( n + 1 );
    ValueType total = scan[loc]( wValues.get(), n );
    return total;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::unscan(
    hmemo::HArray<ValueType>& array,
    hmemo::ContextPtr prefLoc )
{
    const IndexType n = array.size();

    if ( n < 1 )
    {
        return ValueType( 0 );
    }

    static LAMAKernel<UtilKernelTrait::unscan<ValueType> > unscan;
    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    unscan.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteAccess<ValueType> wValues( array, loc );
    ValueType first = unscan[loc]( wValues.get(), n );

    // One additional element will be removed from end
    wValues.resize( n - 1 );

    return first;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::sort(
    hmemo::HArray<IndexType>* perm,
    hmemo::HArray<ValueType>* outValues,
    const hmemo::HArray<ValueType>& inValues,
    const bool ascending,
    hmemo::ContextPtr prefLoc )
{
    const IndexType n = inValues.size();

    if ( n == 0 )
    {
        if ( perm )
        {
            perm->clear();
        }

        return;
    }

    static LAMAKernel<UtilKernelTrait::sort<ValueType> > sort;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = inValues.getValidContext();
    }

    sort.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    // ReadAccess for inValues is always before WriteAccess for outValues, so alias is okay

    ReadAccess<ValueType> rValues( inValues, loc );

    if ( perm == NULL )
    {
        WriteOnlyAccess<ValueType> wValues( *outValues, loc, n );
        sort[loc]( NULL, wValues.get(), rValues.get(), n, ascending );
    }
    else if ( outValues == NULL )
    {
        WriteOnlyAccess<IndexType> wPerm( *perm, loc, n );
        sort[loc]( wPerm.get(), NULL, rValues.get(), n, ascending );
    }
    else
    {
        WriteOnlyAccess<ValueType> wValues( *outValues, loc, n );
        WriteOnlyAccess<IndexType> wPerm( *perm, loc, n );
        sort[loc]( wPerm.get(), wValues.get(), rValues.get(), n, ascending );
    }
}

/* --------------------------------------------------------------------------- */

template<typename BucketType>
void HArrayUtils::bucketSort(
    hmemo::HArray<IndexType>& offsets,
    hmemo::HArray<IndexType>& perm,
    const hmemo::HArray<BucketType>& array,
    const BucketType nb,
    hmemo::ContextPtr prefLoc )
{
    const IndexType n = array.size();

    if ( n == 0 )
    {
        perm.clear();
        IndexType zero = 0;
        offsets.init( zero, nb + 1 );  // offsets = { 0, 0, ..., 0 }
        return;
    }

    static LAMAKernel<UtilKernelTrait::countBuckets<BucketType> > countBuckets;
    static LAMAKernel<UtilKernelTrait::scan<IndexType> > scan;
    static LAMAKernel<UtilKernelTrait::sortInBuckets<BucketType> > sortInBuckets;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    countBuckets.getSupportedContext( loc, scan, sortInBuckets );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<BucketType> bucketMap( array, loc );
    WriteOnlyAccess<IndexType> sizes( offsets, loc, static_cast<IndexType>( nb + 1 ) );

    countBuckets[loc]( sizes.get(), nb, bucketMap, n );
    IndexType total = scan[loc]( sizes.get(), nb );

    // Note: total can be < n if array contains values < 0 or >= nb

    WriteOnlyAccess<IndexType> wPerm( perm, loc, total );
    sortInBuckets[loc]( wPerm, sizes, nb, bucketMap, n );
}

/* --------------------------------------------------------------------------- */

template<typename BucketType>
void HArrayUtils::bucketCount(
    hmemo::HArray<IndexType>& bucketSizes,
    const hmemo::HArray<BucketType>& array,
    const BucketType nb,
    hmemo::ContextPtr prefLoc )
{
    const IndexType n = array.size();

    if ( n == 0 )
    {
        IndexType zeroVal = 0;
        bucketSizes.init( zeroVal, nb );
        return;
    }

    static LAMAKernel<UtilKernelTrait::countBuckets<BucketType> > countBuckets;

    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    countBuckets.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<BucketType> bucketMap( array, loc );

    bucketSizes.reserve( loc, nb + 1 );

    WriteOnlyAccess<IndexType> wSizes( bucketSizes, loc, nb );

    countBuckets[loc]( wSizes, nb, bucketMap, n );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::mergeSort(
    hmemo::HArray<ValueType>& values,
    const hmemo::HArray<IndexType>& offsets,
    bool ascending,
    hmemo::ContextPtr ctx )
{
    mergeSortOptional( values, NULL, offsets, ascending, ctx );
}

template<typename ValueType>
void HArrayUtils::mergeSort(
    hmemo::HArray<ValueType>& values,
    hmemo::HArray<IndexType>& perm,
    const hmemo::HArray<IndexType>& offsets,
    bool ascending,
    hmemo::ContextPtr ctx )
{
    mergeSortOptional( values, &perm, offsets, ascending, ctx );
}

template<typename ValueType>
void HArrayUtils::mergeSortOptional(
    hmemo::HArray<ValueType>& values,
    hmemo::HArray<IndexType>* perm,
    const hmemo::HArray<IndexType>& offsets,
    bool ascending,
    hmemo::ContextPtr )
{
    IndexType n  = values.size();       // number of values to sort
    IndexType nb = offsets.size() - 1;  // number of sorted subarray

    SCAI_LOG_INFO( logger, "mergeSort of values[" << n << "] with " << nb << " sorted subarrays." )

    if ( perm )
    {
        SCAI_ASSERT_EQ_ERROR( values.size(), perm->size(), "array size mismatch" )
    }

    // we need a (writeable) copy of the offset array for running offsets

    HArray<IndexType> copyOffsets( offsets );

    // Temporary array required for the sorted values, in-place not supported

    HArray<ValueType> sorted;
    HArray<IndexType> newPerm;

    {
        WriteAccess<IndexType> wOffsets( copyOffsets );
        ReadAccess<IndexType> rOffsets( offsets );
        WriteOnlyAccess<ValueType> wSorted( sorted, n );
        ReadAccess<ValueType> rValues( values );

        common::unique_ptr<WriteOnlyAccess<IndexType> > wPerm;
        common::unique_ptr<ReadAccess<IndexType> > rPerm;

        if ( perm )
        {
            wPerm.reset( new WriteOnlyAccess<IndexType>( newPerm, n ) );
            rPerm.reset( new ReadAccess<IndexType>( *perm ) );
        }

        for ( IndexType i = 0; i < n; ++i )
        {
            // find the next minimal element

            IndexType nextIndex = nIndex;

            for ( IndexType k = 0; k < nb; ++k )
            {
                if ( wOffsets[k] == rOffsets[k + 1] )
                {
                    continue;  // no more values in subarray k
                }

                if ( nextIndex == nIndex )
                {
                    nextIndex = k;
                }
                else if ( ascending && rValues[wOffsets[k]] < rValues[wOffsets[nextIndex]] )
                {
                    nextIndex = k;
                }
                else if ( !ascending && rValues[wOffsets[k]] > rValues[wOffsets[nextIndex]] )
                {
                    nextIndex = k;
                }
            }

            SCAI_ASSERT_NE_ERROR( nextIndex, nIndex, "no more elements found" )

            IndexType& pos = wOffsets[nextIndex];

            SCAI_LOG_DEBUG( logger, "Next elem in subarray = " << nextIndex
                            << ", pos = " << pos << ", val = " << rValues[pos] );

            // sort in the next minimal element

            wSorted[i] = rValues[pos];

            if ( perm )
            {
                ( *wPerm )[i] = ( *rPerm )[pos];
            }

            pos++;
        }

        // Proof: wOffsets[k] == rOffsets[k+1]

        for ( IndexType k = 0; k < nb; ++k )
        {
            SCAI_ASSERT_EQ_ERROR( rOffsets[k + 1], wOffsets[k], "serious problem during mergesort" )
        }
    }

    // set the sorted values back in the input array, also for perm if available

    values.swap( sorted );

    if ( perm )
    {
        perm->swap( newPerm );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::buildSparseIndexes(
    hmemo::HArray<IndexType>& sparseIndexes,
    const hmemo::HArray<ValueType>& denseArray,
    hmemo::ContextPtr prefLoc )
{
    const IndexType n = denseArray.size();
    static LAMAKernel<UtilKernelTrait::countNonZeros<ValueType> > countNonZeros;
    static LAMAKernel<UtilKernelTrait::compress<ValueType, ValueType> > compress;
    ContextPtr loc = prefLoc;

    // default location for conversion: where we have the dense values

    if ( loc == ContextPtr() )
    {
        loc = denseArray.getValidContext();
    }

    compress.getSupportedContext( loc, countNonZeros );
    SCAI_CONTEXT_ACCESS( loc )
    ValueType eps = common::TypeTraits<ValueType>::eps1();
    ReadAccess<ValueType> rDenseArray( denseArray, loc );
    // we count the non-zeros at first to have sizes for sparse data
    IndexType sparseN = countNonZeros[loc]( rDenseArray.get(), n, eps );
    WriteOnlyAccess<IndexType> wSparseIndexes( sparseIndexes, loc, sparseN );
    sparseN = compress[loc]( NULL, wSparseIndexes.get(), rDenseArray.get(), n, eps );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::buildSparseArray(
    hmemo::HArray<ValueType>& sparseArray,
    hmemo::HArray<IndexType>& sparseIndexes,
    const hmemo::_HArray& denseArray,
    hmemo::ContextPtr prefLoc )
{
    // call buildSparseArrayImpl<ValueType, SourceType> with SourceType = denseArray.getValueType()

    mepr::UtilsWrapperT< ValueType, SCAI_ARRAY_TYPES_HOST_LIST>::buildSparse( sparseArray, sparseIndexes, denseArray, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, typename OtherType>
void HArrayUtils::buildSparseArrayImpl(
    hmemo::HArray<ValueType>& sparseArray,
    hmemo::HArray<IndexType>& sparseIndexes,
    const hmemo::HArray<OtherType>& denseArray,
    hmemo::ContextPtr prefLoc )
{
    const IndexType n = denseArray.size();
    static LAMAKernel<UtilKernelTrait::countNonZeros<OtherType> > countNonZeros;
    static LAMAKernel<UtilKernelTrait::compress<ValueType, OtherType> > compress;
    ContextPtr loc = prefLoc;

    // default location for conversion: where we have the dense values

    if ( loc == ContextPtr() )
    {
        loc = denseArray.getValidContext();
    }

    compress.getSupportedContext( loc, countNonZeros );
    SCAI_CONTEXT_ACCESS( loc )
    ValueType eps = common::TypeTraits<ValueType>::eps1();
    ReadAccess<OtherType> rDenseArray( denseArray, loc );
    // we count the non-zeros at first to have sizes for sparse data
    IndexType sparseN = countNonZeros[loc]( rDenseArray.get(), n, eps );
    WriteOnlyAccess<ValueType> wSparseArray( sparseArray, loc, sparseN );
    WriteOnlyAccess<IndexType> wSparseIndexes( sparseIndexes, loc, sparseN );
    sparseN = compress[loc]( wSparseArray.get(), wSparseIndexes.get(), rDenseArray.get(), n, eps );
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::buildDenseArray(
    hmemo::_HArray& denseArray,
    const IndexType denseN,
    const hmemo::_HArray& sparseArray,
    const hmemo::HArray<IndexType>& sparseIndexes,
    hmemo::ContextPtr prefLoc )
{
    const IndexType sparseN = sparseArray.size();
    SCAI_ASSERT_LE( sparseN, denseN, "sparse array cannot have more entries than dense array." )
    SCAI_ASSERT( HArrayUtils::validIndexes( sparseIndexes, denseN, prefLoc ),
                 "Illegal index for sparse indexes, dense array size = " << denseN );
    // use of existent HArray utilities, even if we have two write accesses for denseArray
    denseArray.clear();
    denseArray.resize( denseN );
    HArrayUtils::assignScalar( denseArray, 0.0, binary::COPY, prefLoc );
    HArrayUtils::scatter( denseArray, sparseIndexes, sparseArray, binary::COPY, prefLoc );
}

/* --------------------------------------------------------------------------- */

IndexType HArrayUtils::findPosInSortedIndexes( const hmemo::HArray<IndexType> indexes, const IndexType pos )
{
    ReadAccess<IndexType> rIndexes( indexes );

    IndexType first = 0;
    IndexType last  = indexes.size();

    while ( first < last )
    {
        IndexType middle = first + ( last - first ) / 2;

        if ( rIndexes[middle] == pos )
        {
            return middle;
        }
        else if ( rIndexes[middle] > pos )
        {
            last = middle;
        }
        else
        {
            first = middle + 1;
        }
    }

    return nIndex;
}

/* --------------------------------------------------------------------------- */

#define HARRAUTILS_SPECIFIER_LVL2( ValueType, OtherValueType )                                                    \
    template void HArrayUtils::gatherImpl<ValueType, OtherValueType>(                                             \
            hmemo::HArray<ValueType>&,                                                                            \
            const hmemo::HArray<OtherValueType>&,                                                                 \
            const hmemo::HArray<IndexType>&,                                                                      \
            const binary::BinaryOp,                                                                               \
            const hmemo::ContextPtr );                                                                            \
    template void HArrayUtils::setArrayImpl<ValueType, OtherValueType>( hmemo::HArray<ValueType>&,                \
            const hmemo::HArray<OtherValueType>&,                                                                 \
            const binary::BinaryOp,                                                                               \
            hmemo::ContextPtr );                                                                                  \
    template void HArrayUtils::setArraySectionImpl<ValueType, OtherValueType>(                                    \
            hmemo::HArray<ValueType>&, const IndexType, const IndexType,                                          \
            const hmemo::HArray<OtherValueType>&, const IndexType, const IndexType,                               \
            const IndexType,                                                                                      \
            const binary::BinaryOp,                                                                               \
            hmemo::ContextPtr );                                                                                  \
    template void HArrayUtils::scatterImpl<ValueType, OtherValueType>(                                            \
            hmemo::HArray<ValueType>&,                                                                            \
            const hmemo::HArray<IndexType>&,                                                                      \
            const hmemo::HArray<OtherValueType>&,                                                                 \
            const binary::BinaryOp,                                                                               \
            const hmemo::ContextPtr );

/*
*/

#define HARRAYUTILS_SPECIFIER( ValueType )                                                                                   \
    template void HArrayUtils::setVal<ValueType>( hmemo::_HArray&, const IndexType, const ValueType );                       \
    template void HArrayUtils::setValImpl<ValueType>( hmemo::HArray<ValueType>&, const IndexType, const ValueType,           \
            const binary::BinaryOp );                                              \
    template ValueType HArrayUtils::getVal<ValueType>( const hmemo::_HArray&, const IndexType );                             \
    template ValueType HArrayUtils::getValImpl<ValueType>( const hmemo::HArray<ValueType>&, const IndexType );               \
    template void HArrayUtils::assignScalar<ValueType>( hmemo::_HArray&, const ValueType,                                    \
            const binary::BinaryOp, hmemo::ContextPtr);                          \
    template void HArrayUtils::setScalar<ValueType>( hmemo::HArray<ValueType>&, const ValueType,                             \
            const binary::BinaryOp, hmemo::ContextPtr);                             \
    template ValueType HArrayUtils::reduce<ValueType>( const hmemo::HArray<ValueType>&,                                      \
            const binary::BinaryOp, hmemo::ContextPtr );                                                                     \
    template ValueType HArrayUtils::reduce2<ValueType>(                                                                      \
            const hmemo::HArray<ValueType>&,                                                                                 \
            const hmemo::HArray<ValueType>&,                                                                                 \
            const binary::BinaryOp,                                                                                          \
            const binary::BinaryOp,                                                                                          \
            hmemo::ContextPtr );                                                                                             \
    template ValueType HArrayUtils::absMaxDiffVal<ValueType>( const hmemo::HArray<ValueType>&,                               \
            const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                                                            \
    template void HArrayUtils::axpy<ValueType>( hmemo::HArray<ValueType>&, const ValueType,                                  \
            const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                                                            \
    template void HArrayUtils::arrayPlusArray<ValueType>( hmemo::HArray<ValueType>&, const ValueType,                        \
            const hmemo::HArray<ValueType>&, const ValueType,                                                                \
            const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                                                            \
    template void HArrayUtils::arrayPlusScalar<ValueType>( hmemo::HArray<ValueType>&, const ValueType,                       \
            const hmemo::HArray<ValueType>&, const ValueType, hmemo::ContextPtr );                                           \
    template void HArrayUtils::arrayTimesArray<ValueType>( hmemo::HArray<ValueType>&, const ValueType,                       \
            const hmemo::HArray<ValueType>&, const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                           \
    template ValueType HArrayUtils::dotProduct<ValueType>( const hmemo::HArray<ValueType>&,                                  \
            const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                                                            \
    template ValueType HArrayUtils::asum<ValueType>( const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                   \
    template ValueType HArrayUtils::nrm2<ValueType>( const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                   \
    template bool HArrayUtils::isSorted<ValueType>( const hmemo::HArray<ValueType>&, const bool, hmemo::ContextPtr );        \
    template ValueType HArrayUtils::scan<ValueType>( hmemo::HArray<ValueType>&, hmemo::ContextPtr );                         \
    template ValueType HArrayUtils::unscan<ValueType>( hmemo::HArray<ValueType>&, hmemo::ContextPtr );                       \
    template void HArrayUtils::sort<ValueType>(                                                                              \
            hmemo::HArray<IndexType>*,                                                                                           \
            hmemo::HArray<ValueType>*,                                                                                           \
            const hmemo::HArray<ValueType>&,                                                                                     \
            const bool,                                                                                                          \
            hmemo::ContextPtr );                                                                                                 \
    template void HArrayUtils::mergeSort<ValueType>(                                                                         \
            hmemo::HArray<ValueType>&,                                                                                           \
            const hmemo::HArray<IndexType>&,                                                                                     \
            const bool,                                                                                                          \
            hmemo::ContextPtr );                                                                                                 \
    template void HArrayUtils::mergeSort<ValueType>(                                                                         \
            hmemo::HArray<ValueType>&,                                                                                           \
            hmemo::HArray<IndexType>&,                                                                                           \
            const hmemo::HArray<IndexType>&,                                                                                     \
            const bool,                                                                                                          \
            hmemo::ContextPtr );                                                                                                 \
    template void HArrayUtils::setSequence<ValueType>( hmemo::HArray<ValueType>&, ValueType, ValueType, IndexType,           \
            hmemo::ContextPtr );                                                                                             \
    template void HArrayUtils::setRandomImpl<ValueType>( hmemo::HArray<ValueType>&, IndexType, float, hmemo::ContextPtr );   \
    template void HArrayUtils::buildSparseArray<ValueType>( hmemo::HArray<ValueType>&, hmemo::HArray<IndexType>&,            \
            const hmemo::_HArray&, hmemo::ContextPtr );                                                            \
    \
    SCAI_COMMON_LOOP_LVL2( ValueType, HARRAUTILS_SPECIFIER_LVL2, SCAI_ARRAY_TYPES_HOST )

SCAI_COMMON_LOOP( HARRAYUTILS_SPECIFIER, SCAI_ARRAY_TYPES_HOST )

#undef HARRAYUTILS_SPECIFIER
#undef HARRAUTILS_SPECIFIER_LVL2

// The unary/binary operations are instantiated for all array types
// But keep in mind that some operations are not supported for IndexType
// and will therefore cause runtime errors, i.e. throw an exception

#define HARRAYUTILS_SPECIFIER( ValueType )                    \
    template void HArrayUtils::unaryOp<ValueType>(            \
            hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<ValueType>&,                      \
            const unary::UnaryOp, hmemo::ContextPtr);             \
    template void HArrayUtils::binaryOp<ValueType>(           \
            hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<ValueType>&,                     \
            const hmemo::HArray<ValueType>&,                     \
            const binary::BinaryOp,                              \
            hmemo::ContextPtr);                                  \
    template void HArrayUtils::binaryOpScalar1<ValueType>(    \
            hmemo::HArray<ValueType>&,                            \
            const ValueType,                                      \
            const hmemo::HArray<ValueType>&,                      \
            const binary::BinaryOp,                               \
            hmemo::ContextPtr);                                   \
    template void HArrayUtils::binaryOpScalar2<ValueType>(    \
            hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<ValueType>&,                      \
            const ValueType,                                      \
            const binary::BinaryOp,                               \
            hmemo::ContextPtr);

// Note: it is now safe to instantiate unary/binary kernels also for IndexType

SCAI_COMMON_LOOP( HARRAYUTILS_SPECIFIER, SCAI_ARRAY_TYPES_HOST )

#undef HARRAYUTILS_SPECIFIER

// ToDo: template instantiation of bucketSort/bucketCount for PartitionId but only if PartitionId != IndexType
// ToDo: template instantiation of bucketSort/bucketCount for PartitionId but only if PartitionId != IndexType

template void HArrayUtils::bucketSort( hmemo::HArray<IndexType>& offsets, hmemo::HArray<IndexType>& perm,
                                       const hmemo::HArray<IndexType>& array, const IndexType nb, hmemo::ContextPtr prefLoc );

template void HArrayUtils::bucketCount( hmemo::HArray<IndexType>& sizes, const hmemo::HArray<IndexType>& array,
                                        const IndexType nb, hmemo::ContextPtr prefLoc );

/* --------------------------------------------------------------------------- */

} /* end namespace utilskernel */

} /* end namespace scai */
