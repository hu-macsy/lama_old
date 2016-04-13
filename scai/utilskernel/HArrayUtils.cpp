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

#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/ScalarType.hpp>

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
    set( target, source, common::reduction::COPY, prefLoc );
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::set(
    _HArray& target,
    const _HArray& source,
    const common::reduction::ReductionOp op,
    const ContextPtr prefLoc )
{
    ContextPtr loc = prefLoc;

    if ( !loc )
    {
        if ( op == common::reduction::COPY )
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

    mepr::UtilsWrapper<ARITHMETIC_ARRAY_HOST_LIST>::setImpl( target, source, op, loc );
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType,typename SourceValueType>
void HArrayUtils::setImpl(
    HArray<TargetValueType>& target,
    const HArray<SourceValueType>& source,
    const common::reduction::ReductionOp op,
    const ContextPtr prefContext )
{
    // verify that dynamic cast operations went okay before

    SCAI_ASSERT_ERROR( &target, "NULL target" )
    SCAI_ASSERT_ERROR( &source, "NULL source" )

    // set should be available on interface for each loc

    static LAMAKernel<UtilKernelTrait::set<TargetValueType, SourceValueType> > set;

    ContextPtr loc = prefContext;

    set.getSupportedContext( loc );

    const IndexType n = source.size();

    SCAI_CONTEXT_ACCESS( loc )

    if ( op == common::reduction::COPY )
    {
        WriteOnlyAccess<TargetValueType> targetVals( target, loc, n );
        ReadAccess<SourceValueType> sourceVals( source, loc );

        // Implemenation of set @ loc is available

        set[loc]( targetVals.get(), sourceVals.get(), n, op );
    }
    else
    {
        WriteAccess<TargetValueType> targetVals( target, loc );
        ReadAccess<SourceValueType> sourceVals( source, loc );

        // Implemenation of set @ loc is available

        set[loc]( targetVals.get(), sourceVals.get(), n, op );
    }
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::gather(
    _HArray& target,
    const _HArray& source,
    const HArray<IndexType>& indexes,
    const ContextPtr prefLoc )
{
    mepr::UtilsWrapper<ARITHMETIC_ARRAY_HOST_LIST>::gatherImpl( target, source, indexes, prefLoc );
}

/* --------------------------------------------------------------------------- */

void HArrayUtils::scatter(
    _HArray& target,
    const HArray<IndexType>& indexes,
    const _HArray& source,
    const ContextPtr prefLoc )
{
    mepr::UtilsWrapper<ARITHMETIC_ARRAY_HOST_LIST>::scatterImpl( target, indexes, source, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename TargetValueType,typename SourceValueType>
void HArrayUtils::gatherImpl(
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

template<typename TargetValueType,typename SourceValueType>
void HArrayUtils::scatterImpl(
    HArray<TargetValueType>& target,
    const HArray<IndexType>& indexes,
    const HArray<SourceValueType>& source,
    const ContextPtr prefLoc )
{
    SCAI_REGION( "HArray.scatter" )

    // choose location for the operation where source array is currently valid

    static LAMAKernel<UtilKernelTrait::setScatter<TargetValueType, SourceValueType> > setScatter;

    ContextPtr loc = prefLoc;

    if ( loc == ContextPtr() )
    {
        loc = source.getValidContext();
    }

    setScatter.getSupportedContext( loc );

    const IndexType n = indexes.size();

    WriteOnlyAccess<TargetValueType> wTarget( target, loc, n );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<SourceValueType> rSource( source, loc );
    ReadAccess<IndexType> rIndexes( indexes, loc );

    //  target[ indexes[i] ] = source[i]

    setScatter[loc] ( wTarget.get(), rIndexes.get(), rSource.get(), n );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::setScalar(
    _HArray& target,
    const ValueType value,
    const common::reduction::ReductionOp op,
    ContextPtr prefLoc )
{
    mepr::UtilsWrapperT<ValueType, ARITHMETIC_ARRAY_HOST_LIST>::setScalarImpl( target, value, op, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::setScalarImpl(
    HArray<ValueType>& target,
    const ValueType value,
    const common::reduction::ReductionOp op,
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

    if ( op == common::reduction::COPY )
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

    mepr::UtilsWrapperT< ValueType, ARITHMETIC_ARRAY_HOST_LIST>::setValImpl( target, index, val );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::setValImpl(
    HArray<ValueType>& target,
    const IndexType index,
    const ValueType val )
{
    ContextPtr loc = target.getValidContext();   // preferred location where to fill

    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

    setVal.getSupportedContext( loc );

    SCAI_LOG_INFO( logger, "setVal<" << common::TypeTraits<ValueType>::id() << ">[" << index
                           << "] = " << val << " @ " << *loc )

    WriteAccess<ValueType> wTarget( target, loc );

    SCAI_CONTEXT_ACCESS( loc )

    setVal[loc]( wTarget.get() + index, 1, val, common::reduction::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType HArrayUtils::getVal(
    const _HArray& array,
    const IndexType index )
{
    ValueType val = mepr::UtilsWrapperT< ValueType, ARITHMETIC_ARRAY_HOST_LIST>::getValImpl( array, index );
    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType, typename OtherValueType>
ValueType HArrayUtils::getValImpl(
    const HArray<OtherValueType>& array,
    const IndexType index )
{
    SCAI_ASSERT_DEBUG( index < array.size(), "index = " << index << " out of range for array = " << array );

    // get the data from a valid context, avoids any memory copy.

    ContextPtr loc = array.getValidContext();

    static LAMAKernel<UtilKernelTrait::getValue<OtherValueType> > getValue;

    getValue.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<OtherValueType> rArray( array, loc );

    OtherValueType val = getValue[loc]( rArray.get(), index );

    return static_cast<ValueType>( val );
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
        setScalarImpl( result, beta, common::reduction::COPY, prefLoc );
    }
    else if( &result == &y )
    {
        if ( beta == common::constants::ONE )
        {
            return;
        }

        // result := beta * result, use setScalar, op == MULT

        setScalarImpl( result, beta, common::reduction::MULT, prefLoc );
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
    setScalarImpl( array, beta, common::reduction::MULT, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::conj( HArray<ValueType>& array, ContextPtr prefLoc )
{
    IndexType n = array.size();

    if ( ( n > 0 ) && common::scalar::isComplex( common::TypeTraits<ValueType>::stype ) )
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
ValueType HArrayUtils::reduce( 
    const HArray<ValueType>& array, 
    const common::reduction::ReductionOp redOp,
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
    ContextPtr prefLoc)
{
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
    ContextPtr prefLoc)
{
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

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void HArrayUtils::invert( HArray<ValueType>& array, ContextPtr prefContext )
{
    const IndexType size = array.size();

    static LAMAKernel<UtilKernelTrait::invert<ValueType> > invert;

    ContextPtr loc = prefContext;

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
    TemplateSpecifier::set( HArrayUtils::axpy<ValueType> );
    TemplateSpecifier::set( HArrayUtils::arrayPlusArray<ValueType> );
    TemplateSpecifier::set( HArrayUtils::dotProduct<ValueType> );
    TemplateSpecifier::set( HArrayUtils::asum<ValueType> );
    TemplateSpecifier::set( HArrayUtils::invert<ValueType> );
}

/* --------------------------------------------------------------------------- */

HArrayUtils::HArrayUtils()
{
    common::mepr::TemplateSpecifierV<SpecifierV, ARITHMETIC_ARRAY_HOST_LIST>::call();
}

/* --------------------------------------------------------------------------- */

HArrayUtils HArrayUtils::guard;

} /* end namespace utilskernel */

} /* end namespace scai */
