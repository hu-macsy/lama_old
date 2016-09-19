/**
 * @file HArrayUtils.cpp
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
 * @brief Definition of class with utility routines.
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
        assignOp( target, source, reduction::COPY, prefLoc );
    }
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::assignOp(
    _HArray& target,
    const _HArray& source,
    const reduction::ReductionOp op,
    const ContextPtr prefLoc )
{
    ContextPtr loc = prefLoc;

    if ( !loc )
    {
        if ( op == reduction::COPY )
        {
            // if no context is given we assign where source has a valid copy available
            loc = source.getValidContext();
        }
        else
        {
            // if no context is given we reduce where target has a valid copy available
            loc = target.getValidContext();
        }
    }

    mepr::UtilsWrapperTT1<SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::setArray( target, source, op, loc );
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::setArray(
    HArray<TargetValueType>& target,
    const HArray<SourceValueType>& source,
    const reduction::ReductionOp op,
    const ContextPtr prefLoc )
{
    // verify that dynamic cast operations went okay before
//    SCAI_ASSERT_ERROR( &target, "NULL target" )
//    SCAI_ASSERT_ERROR( &source, "NULL source" )
    // set should be available on interface for each loc
    static LAMAKernel<UtilKernelTrait::set<TargetValueType, SourceValueType> > set;
    ContextPtr loc = prefLoc;
    set.getSupportedContext( loc );
    const IndexType n = source.size();
    SCAI_CONTEXT_ACCESS( loc )

    if ( op == reduction::COPY )
    {
        ReadAccess<SourceValueType> sourceVals( source, loc );
        WriteOnlyAccess<TargetValueType> targetVals( target, loc, n );
        // Implemenation of set @ loc is available
        set[loc]( targetVals.get(), sourceVals.get(), n, op );
    }
    else
    {
        ReadAccess<SourceValueType> sourceVals( source, loc );
        WriteAccess<TargetValueType> targetVals( target, loc );
        // Implemenation of set @ loc is available
        set[loc]( targetVals.get(), sourceVals.get(), n, op );
    }
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::assignGather(
    _HArray& target,
    const _HArray& source,
    const HArray<IndexType>& indexes,
    const ContextPtr prefLoc )
{
    // use metaprogramming to call the gather version with the correct value types for target and source
    mepr::UtilsWrapperTT1<SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::gather( target, source, indexes, prefLoc );
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::assignScatter(
    _HArray& target,
    const HArray<IndexType>& indexes,
    const _HArray& source,
    const reduction::ReductionOp op,
    const ContextPtr prefLoc )
{
    // use metaprogramming to call the scatter version with the correct value types for target and source
    mepr::UtilsWrapperTT1<SCAI_ARRAY_TYPES_HOST_LIST, SCAI_ARRAY_TYPES_HOST_LIST>::scatter( target, indexes, source, op, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::gather(
    HArray<TargetValueType>& target,
    const HArray<SourceValueType>& source,
    const HArray<IndexType>& indexes,
    const ContextPtr prefLoc )
{
    SCAI_REGION( "HArray.gather" )
    // choose location for the operation where source array is currently valid
    static LAMAKernel<UtilKernelTrait::setGather<TargetValueType, SourceValueType> > setGather;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = source.getValidContext();
    }

    setGather.getSupportedContext( loc );
    const IndexType n = indexes.size();
    WriteOnlyAccess<TargetValueType> wTarget( target, loc, n );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<SourceValueType> rSource( source, loc );
    ReadAccess<IndexType> rIndexes( indexes, loc );
    //  target[i] = source[ indexes[i] ]
    setGather[loc] ( wTarget.get(), rSource.get(), rIndexes.get(), n );
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType, typename SourceValueType>
void HArrayUtils::scatter(
    HArray<TargetValueType>& target,
    const HArray<IndexType>& indexes,
    const HArray<SourceValueType>& source,
    const reduction::ReductionOp op,
    const ContextPtr prefLoc )
{
    SCAI_REGION( "HArray.scatter" )

    SCAI_ASSERT( ( op == reduction::COPY ) || ( op == reduction::ADD ), "Unsupported reduction op " << op  )

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
    const reduction::ReductionOp op,
    ContextPtr prefLoc )
{
    mepr::UtilsWrapperT<ValueType, SCAI_ARRAY_TYPES_HOST_LIST>::setScalar( target, value, op, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::setScalar(
    HArray<ValueType>& target,
    const ValueType value,
    const reduction::ReductionOp op,
    ContextPtr prefLoc )
{
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

    if ( op == reduction::COPY )
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
    const ValueType val )
{
    // setting single value will directly copy to the device with the valid incarnation
    ContextPtr loc = target.getValidContext();  
    WriteAccess<ValueType> wTarget( target, loc );
    wTarget.setValue( val, index );
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
void HArrayUtils::assignScaled(
    HArray<ValueType>& result,
    const ValueType beta,
    const HArray<ValueType>& y,
    ContextPtr prefLoc )
{
    const IndexType n = y.size();  // determines size of result

    // beta = 0    : saves the need of a read access for y
    // result == y : only one write access needed ( write + read not possible)

    if ( beta == common::constants::ZERO )
    {
        setScalar( result, beta, reduction::COPY, prefLoc );
    }
    else if ( &result == &y )
    {
        if ( beta == common::constants::ONE )
        {
            return;
        }

        // result := beta * result, use setScalar, op == MULT
        setScalar( result, beta, reduction::MULT, prefLoc );
    }
    else
    {
        // Note: we do not use BLAS1:axpy here to guarantee same LAMA OpenMP schedule
        //       and to support type conversions in place for multiprecision support
        static LAMAKernel<UtilKernelTrait::setScale<ValueType, ValueType> > setScale;
        ContextPtr loc = prefLoc;

        if ( loc == ContextPtr() )
        {
            loc = y.getValidContext();
        }

        setScale.getSupportedContext( loc );
        SCAI_CONTEXT_ACCESS( loc )
        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( result, loc, n );
        setScale[loc]( wResult.get(), beta, rY.get(), n );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::scale( HArray<ValueType>& array, const ValueType beta, ContextPtr prefLoc )
{
    setScalar( array, beta, reduction::MULT, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::conj( HArray<ValueType>& array, ContextPtr prefLoc )
{
    IndexType n = array.size();

    if ( ( n > 0 ) && common::isComplex( common::TypeTraits<ValueType>::stype ) )
    {
        static LAMAKernel<UtilKernelTrait::conj<ValueType> > conj;
        ContextPtr loc = prefLoc;

        if ( loc == ContextPtr() )
        {
            loc = array.getValidContext();
        }

        conj.getSupportedContext( loc );
        SCAI_CONTEXT_ACCESS( loc )
        WriteAccess<ValueType> values( array, loc );
        conj[loc]( values.get(), n );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::exp( HArray<ValueType>& array, ContextPtr prefLoc )
{
    IndexType n = array.size();

    if ( n > 0 )
    {
        static LAMAKernel<UtilKernelTrait::exp<ValueType> > exp;
        // preferred location: where valid values of the array are available
        ContextPtr loc = array.getValidContext( prefLoc);
        exp.getSupportedContext( loc );
        SCAI_CONTEXT_ACCESS( loc )
        WriteAccess<ValueType> values( array, loc );
        exp[loc]( values.get(), n );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::reduce(
    const HArray<ValueType>& array,
    const reduction::ReductionOp redOp,
    const ContextPtr prefLoc )
{
    static LAMAKernel<UtilKernelTrait::reduce<ValueType> > reduce;
    // preferred location: where valid values of the array are available
    ContextPtr loc = array.getValidContext( prefLoc );
    reduce.getSupportedContext( loc );
    ReadAccess<ValueType> readArray( array, loc );
    SCAI_CONTEXT_ACCESS( loc )
    ValueType redVal = reduce[loc]( readArray.get(), readArray.size(), redOp );
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
    SCAI_ASSERT_EQUAL( array1.size(), array2.size(), "array size mismatch for building differences" )
    static LAMAKernel<UtilKernelTrait::absMaxDiffVal<ValueType> > absMaxDiffVal;
    ContextPtr loc = prefLoc;

    // Rule for default location: where array1 has valid values

    if ( loc == ContextPtr() )
    {
        loc = array1.getValidContext();
    }

    absMaxDiffVal.getSupportedContext( loc );
    ReadAccess<ValueType> readArray1( array1, loc );
    ReadAccess<ValueType> readArray2( array2, loc );
    SCAI_CONTEXT_ACCESS( loc )
    ValueType redVal = absMaxDiffVal[loc]( readArray1.get(), readArray2.get(), readArray1.size() );
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

    if ( alpha == scai::common::constants::ZERO )
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

    if ( beta == scai::common::constants::ZERO )
    {
        assignScaled( result, alpha, x, prefLoc );
        return;
    }

    if ( alpha == scai::common::constants::ZERO )
    {
        assignScaled( result, beta, y, prefLoc );
        return;
    }

    if ( &y == &result && beta == scai::common::constants::ONE )
    {
        axpy( result, alpha, x, prefLoc );
        return;
    }

    if ( &x == &result && alpha == scai::common::constants::ONE )
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
void HArrayUtils::arrayTimesArray(
    hmemo::HArray<ValueType>& result,
    const ValueType alpha,
    const hmemo::HArray<ValueType>& x,
    const hmemo::HArray<ValueType>& y,
    hmemo::ContextPtr prefLoc )
{
    SCAI_ASSERT_EQUAL( x.size(), y.size(), "size mismatch" )
    const IndexType n = x.size();
    
    if ( alpha == scai::common::constants::ZERO )
    {
        setScalar( result, ValueType(0.0), reduction::COPY, prefLoc );
        return;
    }

    static LAMAKernel<UtilKernelTrait::vectorScale<ValueType> > vectorScale;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = x.getValidContext();
    }

    vectorScale.getSupportedContext( loc );
    // no alias checks are done here
    {
        ReadAccess<ValueType> xAccess( x, loc );
        ReadAccess<ValueType> yAccess( y, loc );
        WriteOnlyAccess<ValueType> resultAccess( result, loc, x.size() );
        SCAI_CONTEXT_ACCESS( loc )
        vectorScale[loc]( resultAccess.get(), xAccess.get(), yAccess.get(), n );
    }

    if ( alpha != scai::common::constants::ONE )
    {
        setScalar( result, alpha, reduction::MULT, prefLoc );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::invert( HArray<ValueType>& array, ContextPtr prefLoc )
{
    const IndexType size = array.size();
    static LAMAKernel<UtilKernelTrait::invert<ValueType> > invert;
    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    invert.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteAccess<ValueType> wValues( array, loc );
    invert[loc]( wValues.get(), size );
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
    hmemo::HArray<ValueType>& array,
    hmemo::HArray<IndexType>& perm,
    hmemo::ContextPtr prefLoc )
{
    const IndexType n = array.size();

    if ( n == 0 )
    {
        perm.clear();
        return;
    }

    static LAMAKernel<UtilKernelTrait::sort<ValueType> > sort;
    ContextPtr loc = prefLoc;

    // default location for check: where we have valid entries

    if ( loc == ContextPtr() )
    {
        loc = array.getValidContext();
    }

    sort.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteAccess<ValueType> wValues( array, loc );
    WriteOnlyAccess<IndexType> wPerm( perm, loc, n );
    sort[loc]( wValues.get(), wPerm.get(), n );
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
void HArrayUtils::buildSparseArray(
    hmemo::HArray<ValueType>& sparseArray,
    hmemo::HArray<IndexType>& sparseIndexes,
    const hmemo::HArray<ValueType>& denseArray,
    hmemo::ContextPtr prefLoc )
{
    const IndexType n = denseArray.size();
    static LAMAKernel<UtilKernelTrait::countNonZeros<ValueType> > countNonZeros;
    static LAMAKernel<UtilKernelTrait::compress<ValueType> > compress;
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
    WriteOnlyAccess<ValueType> wSparseArray( sparseArray, loc, sparseN );
    WriteOnlyAccess<IndexType> wSparseIndexes( sparseIndexes, loc, sparseN );
    sparseN = compress[loc]( wSparseArray.get(), wSparseIndexes.get(), rDenseArray.get(), n, eps );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::buildDenseArray(
    hmemo::HArray<ValueType>& denseArray,
    const IndexType denseN,
    const hmemo::HArray<ValueType>& sparseArray,
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
    HArrayUtils::setScalar( denseArray, ValueType( 0 ), reduction::COPY, prefLoc );
    HArrayUtils::scatter( denseArray, sparseIndexes, sparseArray, reduction::COPY, prefLoc );
}

/* --------------------------------------------------------------------------- */

#define HARRAUTILS_SPECIFIER_LVL2( ValueType, OtherValueType )                                                          \
    template void HArrayUtils::gather<ValueType, OtherValueType>( hmemo::HArray<ValueType>&,                            \
            const hmemo::HArray<OtherValueType>&,                                                                       \
            const hmemo::HArray<IndexType>&,                                                                            \
            const hmemo::ContextPtr );                                                                                  \
    template void HArrayUtils::setArray<ValueType, OtherValueType>( hmemo::HArray<ValueType>&,                          \
            const hmemo::HArray<OtherValueType>&,                                                                       \
            const reduction::ReductionOp,                                                                               \
            hmemo::ContextPtr );                                                                                        \
    template void HArrayUtils::scatter<ValueType, OtherValueType>( hmemo::HArray<ValueType>&,                           \
            const hmemo::HArray<IndexType>&,                                                                            \
            const hmemo::HArray<OtherValueType>&,                                                                       \
            const reduction::ReductionOp,                                                                               \
            const hmemo::ContextPtr );

#define HARRAYUTILS_SPECIFIER( ValueType )                                                                                        \
    template void HArrayUtils::setVal<ValueType>( hmemo::_HArray&, const IndexType, const ValueType );                            \
    template ValueType HArrayUtils::getVal<ValueType>( const hmemo::_HArray&, const IndexType );                                  \
    template ValueType HArrayUtils::getValImpl<ValueType>( const hmemo::HArray<ValueType>&, const IndexType );                    \
    template void HArrayUtils::assignScalar<ValueType>( hmemo::_HArray&, const ValueType,                                         \
                                                        const reduction::ReductionOp, hmemo::ContextPtr);                         \
    template void HArrayUtils::setScalar<ValueType>( hmemo::HArray<ValueType>&, const ValueType,                                  \
                                                     const reduction::ReductionOp, hmemo::ContextPtr);                            \
    template void HArrayUtils::assignScaled<ValueType>( hmemo::HArray<ValueType>&, const ValueType,                               \
            const hmemo::HArray<ValueType>&, hmemo::ContextPtr);                                                                  \
    template void HArrayUtils::scale<ValueType>( hmemo::HArray<ValueType>&, const ValueType, hmemo::ContextPtr );                 \
    template void HArrayUtils::conj<ValueType>( hmemo::HArray<ValueType>&, hmemo::ContextPtr);                                    \
    template void HArrayUtils::exp<ValueType>( hmemo::HArray<ValueType>&, hmemo::ContextPtr);                                     \
    template ValueType HArrayUtils::reduce<ValueType>( const hmemo::HArray<ValueType>&,                                           \
            const reduction::ReductionOp, hmemo::ContextPtr );                                                                    \
    template ValueType HArrayUtils::absMaxDiffVal<ValueType>( const hmemo::HArray<ValueType>&,                                    \
            const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                                                                 \
    template void HArrayUtils::axpy<ValueType>( hmemo::HArray<ValueType>&, const ValueType,                                       \
            const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                                                                 \
    template void HArrayUtils::arrayPlusArray<ValueType>( hmemo::HArray<ValueType>&, const ValueType,                             \
            const hmemo::HArray<ValueType>&, const ValueType,                                                                     \
            const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                                                                 \
    template void HArrayUtils::arrayTimesArray<ValueType>( hmemo::HArray<ValueType>&, const ValueType,                            \
            const hmemo::HArray<ValueType>&, const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                                \
    template ValueType HArrayUtils::dotProduct<ValueType>( const hmemo::HArray<ValueType>&,                                       \
            const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                                                                 \
    template ValueType HArrayUtils::asum<ValueType>( const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                        \
    template ValueType HArrayUtils::nrm2<ValueType>( const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                        \
    template void HArrayUtils::invert<ValueType>( hmemo::HArray<ValueType>&, hmemo::ContextPtr );                                 \
    template bool HArrayUtils::isSorted<ValueType>( const hmemo::HArray<ValueType>&, const bool, hmemo::ContextPtr );             \
    template ValueType HArrayUtils::scan<ValueType>( hmemo::HArray<ValueType>&, hmemo::ContextPtr );                              \
    template ValueType HArrayUtils::unscan<ValueType>( hmemo::HArray<ValueType>&, hmemo::ContextPtr );                            \
    template void HArrayUtils::sort<ValueType>( hmemo::HArray<ValueType>&, hmemo::HArray<IndexType>&, hmemo::ContextPtr );        \
    template void HArrayUtils::setSequence<ValueType>( hmemo::HArray<ValueType>&, ValueType, ValueType, IndexType,                \
            hmemo::ContextPtr );                                                                                                  \
    template void HArrayUtils::setRandomImpl<ValueType>( hmemo::HArray<ValueType>&, IndexType, float, hmemo::ContextPtr );        \
    template void HArrayUtils::buildSparseArray<ValueType>( hmemo::HArray<ValueType>&, hmemo::HArray<IndexType>&,                 \
            const hmemo::HArray<ValueType>&, hmemo::ContextPtr );                                                                 \
    template void HArrayUtils::buildDenseArray<ValueType>( hmemo::HArray<ValueType>&, const IndexType,                            \
            const hmemo::HArray<ValueType>&,                                                                                      \
            const hmemo::HArray<IndexType>&, hmemo::ContextPtr );                                                                 \
                                                                                                                                  \
    SCAI_COMMON_LOOP_LVL2( ValueType, HARRAUTILS_SPECIFIER_LVL2, SCAI_ARRAY_TYPES_HOST )

SCAI_COMMON_LOOP( HARRAYUTILS_SPECIFIER, SCAI_ARRAY_TYPES_HOST )

#undef HARRAYUTILS_SPECIFIER
#undef HARRAUTILS_SPECIFIER_LVL2

// ToDo: template instantiation of bucketSort/bucketCount for PartitionId but only if PartitionId != IndexType

template void HArrayUtils::bucketSort( hmemo::HArray<IndexType>& offsets, hmemo::HArray<IndexType>& perm,
                                       const hmemo::HArray<IndexType>& array, const IndexType nb, hmemo::ContextPtr prefLoc );

template void HArrayUtils::bucketCount( hmemo::HArray<IndexType>& sizes, const hmemo::HArray<IndexType>& array,
                                        const IndexType nb, hmemo::ContextPtr prefLoc );

/* --------------------------------------------------------------------------- */

} /* end namespace utilskernel */

} /* end namespace scai */
