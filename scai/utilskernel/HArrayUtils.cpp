/**
 * @file HArrayUtils.cpp
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
 * @brief Definition of class with utility routines.
 * @author Thomas Brandes
 * @date 10.10.2011
 * @since 1.0.0
 */

// hpp
#include <scai/utilskernel/HArrayUtils.hpp>

// local library
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/utilskernel/mepr/UtilsWrapper.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/ScalarType.hpp>

// std
#include <iostream>

namespace scai
{

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( HArrayUtils::logger, "HArrayUtils" )

void HArrayUtils::assign( hmemo::_HArray& target, const hmemo::_HArray& source, const hmemo::ContextPtr loc /* = ContextPtr() */)
{
    hmemo::ContextPtr validLoc = loc;

    if( !validLoc )
    {
        // if no context is given we assign where source has a valid copy available

        validLoc = source.getValidContext();
    }

    set( target, source, common::reduction::COPY, validLoc );
}

void HArrayUtils::set(
    hmemo::_HArray& target,
    const hmemo::_HArray& source,
    const common::reduction::ReductionOp op,
    const hmemo::ContextPtr loc )
{
    hmemo::ContextPtr validLoc = loc;

    if ( !validLoc )
    {
        if ( op == common::reduction::COPY )
        {
            // if no context is given we assign where source has a valid copy available

            validLoc = source.getValidContext();
        }
        else
        {
            // if no context is given we reduce where target has a valid copy available

            validLoc = target.getValidContext();
        }
    }

    mepr::UtilsWrapper<ARITHMETIC_ARRAY_HOST_LIST>::setImpl( target, source, op, validLoc );
}

template<typename ValueType1,typename ValueType2>
void HArrayUtils::setImpl(
    hmemo::HArray<ValueType1>& target,
    const hmemo::HArray<ValueType2>& source,
    const common::reduction::ReductionOp op,
    const hmemo::ContextPtr prefContext )
{
    // verify that dynamic cast operations went okay before

    SCAI_ASSERT_ERROR( &target, "NULL target" )
    SCAI_ASSERT_ERROR( &source, "NULL source" )

    // set should be available on interface for each loc

    static LAMAKernel<UtilKernelTrait::set<ValueType1, ValueType2> > set;

    hmemo::ContextPtr loc = set.getValidContext( prefContext );

    const IndexType n = source.size();

    SCAI_CONTEXT_ACCESS( loc )

    if ( op == common::reduction::COPY )
    {
        hmemo::WriteOnlyAccess<ValueType1> targetVals( target, loc, n );
        hmemo::ReadAccess<ValueType2> sourceVals( source, loc );

        // Implemenation of set @ loc is available

        set[loc]( targetVals.get(), sourceVals.get(), n, op );
    }
    else
    {
        hmemo::WriteAccess<ValueType1> targetVals( target, loc );
        hmemo::ReadAccess<ValueType2> sourceVals( source, loc );

        // Implemenation of set @ loc is available

        set[loc]( targetVals.get(), sourceVals.get(), n, op );
    }
}

void HArrayUtils::gather(
    hmemo::_HArray& target,
    const hmemo::_HArray& source,
    const hmemo::HArray<IndexType>& indexes )
{
    SCAI_REGION( "HArray.gather" )
    mepr::UtilsWrapper<ARITHMETIC_ARRAY_HOST_LIST>::gatherImpl( target, source, indexes );
}

template<typename ValueType1,typename ValueType2>
void HArrayUtils::gatherImpl(
    hmemo::HArray<ValueType1>& target,
    const hmemo::HArray<ValueType2>& source,
    const hmemo::HArray<IndexType>& indexes )
{
    SCAI_REGION( "HArray.gatherImpl" )

    // choose location for the operation where source array is currently valid

    static LAMAKernel<UtilKernelTrait::setGather<ValueType1, ValueType2> > setGather;

    hmemo::ContextPtr context = setGather.getValidContext( source.getValidContext() );

    const IndexType n = indexes.size();

    hmemo::WriteOnlyAccess<ValueType1> wTarget( target, context, n );

    SCAI_CONTEXT_ACCESS( context )

    hmemo::ReadAccess<ValueType2> rSource( source, context );
    hmemo::ReadAccess<IndexType> rIndexes( indexes, context );

    //  target[i] = source[ indexes[i] ]

    setGather[context] ( wTarget.get(), rSource.get(), rIndexes.get(), n );
}

template<typename ValueType>
void HArrayUtils::setScalar(
    hmemo::_HArray& target,
    const ValueType value,
    const common::reduction::ReductionOp op,
    hmemo::ContextPtr prefContext )
{
    mepr::UtilsWrapperT< ValueType, ARITHMETIC_ARRAY_HOST_LIST>::setScalarImpl( target, value, op, prefContext );
}

template<typename ValueType, typename OtherValueType>
void HArrayUtils::setScalarImpl(
    hmemo::HArray<ValueType>& target,
    const OtherValueType value,
    const common::reduction::ReductionOp op,
    hmemo::ContextPtr prefContext )
{
    static LAMAKernel<UtilKernelTrait::setVal<ValueType, OtherValueType> > setVal;

    hmemo::ContextPtr context = setVal.getValidContext( prefContext );

    const IndexType n = target.size();

    SCAI_LOG_INFO( logger, target << " = " << value << ", to do at " << *context << ", n = " << n )

    if ( op == common::reduction::COPY )
    {
        // Note: very important is to specify the size n here as it might not have been allocated

        hmemo::WriteOnlyAccess<ValueType> wTarget( target, context, n );
        SCAI_CONTEXT_ACCESS( context )
        setVal[context]( wTarget.get(), n, value, op );
    }
    else
    {
        hmemo::WriteAccess<ValueType> wTarget( target, context );
        SCAI_CONTEXT_ACCESS( context )
        setVal[context]( wTarget.get(), n, value, op );
    }
}

template<typename ValueType>
void HArrayUtils::setVal(
    hmemo::_HArray& target,
    const IndexType index,
    const ValueType val )
{
    SCAI_ASSERT_DEBUG( index < target.size(), "index = " << index << " out of range for target = " << target );

    mepr::UtilsWrapperT< ValueType, ARITHMETIC_ARRAY_HOST_LIST>::setValImpl( target, index, val );
}

template<typename ValueType, typename OtherValueType>
void HArrayUtils::setValImpl(
    hmemo::HArray<ValueType>& target,
    const IndexType index,
    const OtherValueType val )
{
    hmemo::ContextPtr loc = target.getValidContext();   // preferred location where to fill

    static LAMAKernel<UtilKernelTrait::setVal<ValueType, OtherValueType> > setVal;

    loc = setVal.getValidContext( loc );

    SCAI_LOG_INFO( logger, "setVal<" << common::TypeTraits<ValueType>::id() << ">[" << index
                           << "] = " << val << " @ " << *loc )

    hmemo::WriteAccess<ValueType> wTarget( target, loc );

    SCAI_CONTEXT_ACCESS( loc )

    setVal[loc]( wTarget.get() + index, 1, val, common::reduction::COPY );
}

template<typename ValueType>
ValueType HArrayUtils::getVal(
    const hmemo::_HArray& array,
    const IndexType index )
{
    ValueType val = mepr::UtilsWrapperT< ValueType, ARITHMETIC_ARRAY_HOST_LIST>::getValImpl( array, index );
    return val;
}

template<typename ValueType, typename OtherValueType>
ValueType HArrayUtils::getValImpl(
    const hmemo::HArray<OtherValueType>& array,
    const IndexType index )
{
    SCAI_ASSERT_DEBUG( index < array.size(), "index = " << index << " out of range for array = " << array );

    // get the data from a valid context, avoids any memory copy.

    hmemo::ContextPtr loc = array.getValidContext();

    static LAMAKernel<UtilKernelTrait::getValue<OtherValueType> > getValue;

    loc = getValue.getValidContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    hmemo::ReadAccess<OtherValueType> rArray( array, loc );

    OtherValueType val = getValue[loc]( rArray.get(), index );

    return static_cast<ValueType>( val );
}

template<typename ValueType>
void HArrayUtils::assignScaled(
    hmemo::HArray<ValueType>& result,
    const ValueType beta,
    const hmemo::HArray<ValueType>& y,
    hmemo::ContextPtr prefLoc )
{
    const IndexType n = y.size();  // determines size of result

    // beta = 0    : saves the need of a read access for y
    // result == y : only one write access needed ( write + read not possible)

    if( beta == common::constants::ZERO )
    {
        // result := 0

        static LAMAKernel<UtilKernelTrait::setVal<ValueType, ValueType> > setVal;

        hmemo::ContextPtr loc = setVal.getValidContext( prefLoc );

        hmemo::WriteOnlyAccess<ValueType> wResult( result, loc, n );

        SCAI_CONTEXT_ACCESS( loc )

        setVal[loc]( wResult.get(), n, ValueType( 0 ), common::reduction::COPY );
    }
    else if( &result == &y )
    {
        if( beta == common::constants::ONE )
        {
            return;
        }

        // result := beta * result, use setVal, op == MULL

        static LAMAKernel<UtilKernelTrait::setVal<ValueType, ValueType> > setVal;

        hmemo::ContextPtr loc = setVal.getValidContext( prefLoc );

        hmemo::WriteAccess<ValueType> wResult( result, loc );

        SCAI_CONTEXT_ACCESS( loc )

        setVal[loc]( wResult.get(), n, beta, common::reduction::MULT );
    }
    else
    {
        // result := beta * y
        
        // Note: we do not use BLAS1:axpy here to guarantee same LAMA OpenMP schedule
        //       and to support type conversions in place for multiprecision support
        
        static LAMAKernel<UtilKernelTrait::setScale<ValueType, ValueType> > setScale;

        hmemo::ContextPtr loc = setScale.getValidContext( prefLoc );

        SCAI_CONTEXT_ACCESS( loc )

        hmemo::ReadAccess<ValueType> rY( y, loc );
        hmemo::WriteOnlyAccess<ValueType> wResult( result, loc, n );

        setScale[loc]( wResult.get(), beta, rY.get(), n );
    }
}

template<typename ValueType>
void HArrayUtils::scale( hmemo::HArray<ValueType>& array, const ValueType beta, hmemo::ContextPtr prefLoc )
{
    setScalar( array, beta, common::reduction::MULT, prefLoc );
}

template<typename ValueType>
void HArrayUtils::conj( hmemo::HArray<ValueType>& array, hmemo::ContextPtr prefLoc )
{
    IndexType n = array.size();

    if ( ( n > 0 ) && common::scalar::isComplex( common::TypeTraits<ValueType>::stype ) )
    {
        static LAMAKernel<UtilKernelTrait::conj<ValueType> > conj;

        hmemo::ContextPtr loc = conj.getValidContext( prefLoc );

        SCAI_CONTEXT_ACCESS( loc )

        hmemo::WriteAccess<ValueType> values( array, loc );

        conj[loc]( values.get(), n );
    }
}

template<typename ValueType>
void HArrayUtils::addScalar( hmemo::HArray<ValueType>& array, const ValueType scalar, hmemo::ContextPtr prefLoc )
{
    IndexType n = array.size();

    static LAMAKernel<UtilKernelTrait::addScalar<ValueType> > addScalar;

    hmemo::ContextPtr loc = addScalar.getValidContext( prefLoc );

    SCAI_CONTEXT_ACCESS( loc )

    hmemo::WriteAccess<ValueType> values( array, loc );

    addScalar[loc]( values.get(), n, scalar );
}

template<typename ValueType>
ValueType HArrayUtils::reduce( const hmemo::HArray<ValueType>& array, const common::reduction::ReductionOp redOp )
{
    static LAMAKernel<UtilKernelTrait::reduce<ValueType> > reduce;

    // preferred location: where valid values of the array are available

    hmemo::ContextPtr loc = reduce.getValidContext( array.getValidContext() );

    hmemo::ReadAccess<ValueType> readArray( array, loc );

    SCAI_CONTEXT_ACCESS( loc )

    ValueType redVal = reduce[loc]( readArray.get(), readArray.size(), redOp );

    return redVal;
}

template<typename ValueType>
ValueType HArrayUtils::absMaxDiffVal(
    const hmemo::HArray<ValueType>& array1,
    const hmemo::HArray<ValueType>& array2 )
{
    SCAI_ASSERT_EQUAL( array1.size(), array2.size(), "array size mismatch for building differences" )

    static LAMAKernel<UtilKernelTrait::absMaxDiffVal<ValueType> > absMaxDiffVal;

    // Rule for location: where array1 has valid values

    hmemo::ContextPtr loc = absMaxDiffVal.getValidContext( array1.getValidContext() );

    hmemo::ReadAccess<ValueType> readArray1( array1, loc );
    hmemo::ReadAccess<ValueType> readArray2( array2, loc );

    SCAI_CONTEXT_ACCESS( loc )

    ValueType redVal = absMaxDiffVal[loc]( readArray1.get(), readArray2.get(), readArray1.size() );

    return redVal;
}

template<typename ValueType>
void HArrayUtils::SpecifierV<ValueType>::specify()
{
    using common::mepr::TemplateSpecifier;

    TemplateSpecifier::set( HArrayUtils::setVal<ValueType> );
    TemplateSpecifier::set( HArrayUtils::getVal<ValueType> );
    TemplateSpecifier::set( HArrayUtils::setScalar<ValueType> );
    TemplateSpecifier::set( HArrayUtils::assignScaled<ValueType> );
    TemplateSpecifier::set( HArrayUtils::scale<ValueType> );
    TemplateSpecifier::set( HArrayUtils::conj<ValueType> );
    TemplateSpecifier::set( HArrayUtils::reduce<ValueType> );
    TemplateSpecifier::set( HArrayUtils::absMaxDiffVal<ValueType> );
}

HArrayUtils::HArrayUtils()
{
    common::mepr::TemplateSpecifierV<SpecifierV, ARITHMETIC_ARRAY_HOST_LIST>::call();
}

HArrayUtils HArrayUtils::guard;

} /* end namespace utilskernel */

} /* end namespace scai */
