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
#include <scai/common/preprocessor.hpp>
#include <scai/common/ScalarType.hpp>

// std
#include <iostream>

using namespace scai::hmemo;

namespace scai
{

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( HArrayUtils::logger, "HArrayUtils" )

void HArrayUtils::assign( _HArray& target, const _HArray& source, const ContextPtr loc /* = ContextPtr() */)
{
    ContextPtr validLoc = loc;

    if( !validLoc )
    {
        // if no context is given we assign where source has a valid copy available

        validLoc = source.getValidContext();
    }

    set( target, source, common::reduction::COPY, validLoc );
}

template<typename ValueType1,typename ValueType2>
void HArrayUtils::setImpl(
    HArray<ValueType1>& target,
    const HArray<ValueType2>& source,
    const common::reduction::ReductionOp op,
    const ContextPtr prefContext )
{
    // verify that dynamic cast operations went okay before

    SCAI_ASSERT_ERROR( &target, "NULL target" )
    SCAI_ASSERT_ERROR( &source, "NULL source" )

    // set should be available on interface for each loc

    static LAMAKernel<UtilKernelTrait::set<ValueType1, ValueType2> > set;

    ContextPtr loc = set.getValidContext( prefContext );

    const IndexType n = source.size();

    SCAI_CONTEXT_ACCESS( loc )

    if ( op == common::reduction::COPY )
    {
        WriteOnlyAccess<ValueType1> targetVals( target, loc, n );
        ReadAccess<ValueType2> sourceVals( source, loc );

        // Implemenation of set @ loc is available

        set[loc]( targetVals.get(), sourceVals.get(), n, op );
    }
    else
    {
        WriteAccess<ValueType1> targetVals( target, loc );
        ReadAccess<ValueType2> sourceVals( source, loc );

        // Implemenation of set @ loc is available

        set[loc]( targetVals.get(), sourceVals.get(), n, op );
    }
}



void HArrayUtils::set( 
    _HArray& target, 
    const _HArray& source, 
    const common::reduction::ReductionOp op, 
    const ContextPtr loc )
{
    ContextPtr validLoc = loc;

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

void HArrayUtils::gather(
    _HArray& target,
    const _HArray& source,
    const HArray<IndexType>& indexes )
{
    SCAI_REGION( "HArray.gather" )

    mepr::UtilsWrapper<ARITHMETIC_ARRAY_HOST_LIST>::gatherImpl( target, source, indexes );
}

template<typename ValueType1,typename ValueType2>
void HArrayUtils::gatherImpl(
    HArray<ValueType1>& target,
    const HArray<ValueType2>& source,
    const HArray<IndexType>& indexes )
{
    SCAI_REGION( "HArray.gatherImpl" )

    // choose location for the operation where source array is currently valid

    static LAMAKernel<UtilKernelTrait::setGather<ValueType1, ValueType2> > setGather;

    ContextPtr context = setGather.getValidContext( source.getValidContext() );

    const IndexType n = indexes.size();

    WriteOnlyAccess<ValueType1> wTarget( target, context, n );

    SCAI_CONTEXT_ACCESS( context )

    ReadAccess<ValueType2> rSource( source, context );
    ReadAccess<IndexType> rIndexes( indexes, context );

    //  target[i] = source[ indexes[i] ]

    setGather[context] ( wTarget.get(), rSource.get(), rIndexes.get(), n );
}

template<typename ValueType>
void HArrayUtils::setScalar( _HArray& target, const ValueType value, const common::reduction::ReductionOp op, ContextPtr prefContext )
{
    mepr::UtilsWrapperT< ValueType, ARITHMETIC_ARRAY_HOST_LIST>::setScalarImpl( target, value, op, prefContext );
}

template<typename ValueType, typename OtherValueType>
void HArrayUtils::setScalarImpl( HArray<ValueType>& target, const OtherValueType value, const common::reduction::ReductionOp op, ContextPtr prefContext )
{
    static LAMAKernel<UtilKernelTrait::setVal<ValueType, OtherValueType> > setVal;

    ContextPtr context = setVal.getValidContext( prefContext );

    const IndexType n = target.size();

    SCAI_LOG_INFO( logger, target << " = " << value << ", to do at " << *context << ", n = " << n )

    if ( op == common::reduction::COPY )
    {
        // Note: very important is to specify the size n here as it might not have been allocated

        WriteOnlyAccess<ValueType> wTarget( target, context, n );
        SCAI_CONTEXT_ACCESS( context )
        setVal[context]( wTarget.get(), n, value, op );
    }
    else
    {
        WriteAccess<ValueType> wTarget( target, context );
        SCAI_CONTEXT_ACCESS( context )
        setVal[context]( wTarget.get(), n, value, op );
    }
}

template<typename ValueType>
void HArrayUtils::setVal( HArray<ValueType>& target, const IndexType index, ValueType val )
{
    SCAI_ASSERT_DEBUG( index < target.size(), "index = " << index << " out of range for target = " << target );

    ContextPtr loc = target.getValidContext();   // preferred location where to fill

    static LAMAKernel<UtilKernelTrait::setVal<ValueType, ValueType> > setVal;

    loc = setVal.getValidContext( loc );

    SCAI_LOG_INFO( logger, "setVal<" << common::TypeTraits<ValueType>::id() << ">[" << index 
                           << "] = " << val << " @ " << *loc )

    WriteAccess<ValueType> wTarget( target, loc );

    SCAI_CONTEXT_ACCESS( loc )

    setVal[loc]( wTarget.get() + index, 1, val, common::reduction::COPY );
}

template<typename ValueType>
ValueType HArrayUtils::getVal( const HArray<ValueType>& array, const IndexType index )
{
    SCAI_ASSERT_DEBUG( index < array.size(), "index = " << index << " out of range for array = " << array );

    // get the data from a valid context, avoids any memory copy.

    ContextPtr loc = array.getValidContext();

    static LAMAKernel<UtilKernelTrait::getValue<ValueType> > getValue;

    loc = getValue.getValidContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<ValueType> rArray( array, loc );

    return getValue[loc]( rArray.get(), index );
}

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

    if( beta == common::constants::ZERO )
    {
        // result := 0

        static LAMAKernel<UtilKernelTrait::setVal<ValueType, ValueType> > setVal;

        ContextPtr loc = setVal.getValidContext( prefLoc );

        WriteOnlyAccess<ValueType> wResult( result, loc, n );

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

        ContextPtr loc = setVal.getValidContext( prefLoc );

        WriteAccess<ValueType> wResult( result, loc );

        SCAI_CONTEXT_ACCESS( loc )

        setVal[loc]( wResult.get(), n, beta, common::reduction::MULT );
    }
    else
    {
        // result := beta * y
        
        // Note: we do not use BLAS1:axpy here to guarantee same LAMA OpenMP schedule
        //       and to support type conversions in place for multiprecision support
        
        static LAMAKernel<UtilKernelTrait::setScale<ValueType, ValueType> > setScale;

        ContextPtr loc = setScale.getValidContext( prefLoc );

        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( result, loc, n );

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

        ContextPtr loc = conj.getValidContext( prefLoc );

        SCAI_CONTEXT_ACCESS( loc )

        WriteAccess<ValueType> values( array, loc );

        conj[loc]( values.get(), n );
    }
}

// template instantiation for the supported data types

template void HArrayUtils::setVal( hmemo::HArray<IndexType>& , const IndexType , IndexType );
template void HArrayUtils::setScalar( hmemo::_HArray& , const IndexType , const common::reduction::ReductionOp op, const ContextPtr ctx );

template IndexType HArrayUtils::getVal( const hmemo::HArray<IndexType>& , const IndexType );

template void HArrayUtils::assignScaled(
    hmemo::HArray<IndexType>& ,
    const IndexType ,
    const HArray<IndexType>& ,
    hmemo::ContextPtr  );

/** Macro instantiates operations that have also type conversion */

/** Macro instantiates operations for supported arithmetic types */

#define HARRAY_UTILS_INSTANTIATE(z, I, _)                                       \
    template                                                                        \
    void HArrayUtils::setVal(                                                    \
            HArray<ARITHMETIC_HOST_TYPE_##I>& target,                            \
            const IndexType index,                                                  \
            ARITHMETIC_HOST_TYPE_##I val );                                         \
                                                                                    \
    template                                                                        \
    void HArrayUtils::scale(                                                         \
            HArray<ARITHMETIC_HOST_TYPE_##I>& array,                                 \
            const ARITHMETIC_HOST_TYPE_##I beta,                                     \
            ContextPtr prefLoc );                                                   \
                                                                                    \
    template                                                                        \
    void HArrayUtils::conj(                                                         \
            HArray<ARITHMETIC_HOST_TYPE_##I>& target,                            \
            ContextPtr loc );                                                       \
                                                                                    \
    template                                                                        \
    ARITHMETIC_HOST_TYPE_##I HArrayUtils::getVal(                                   \
            const HArray<ARITHMETIC_HOST_TYPE_##I>&,                               \
            const IndexType );                                                      \
                                                                                    \
    template                                                                        \
    void HArrayUtils::assignScaled(                                              \
            HArray<ARITHMETIC_HOST_TYPE_##I>& result,                            \
            const ARITHMETIC_HOST_TYPE_##I beta,                                    \
            const HArray<ARITHMETIC_HOST_TYPE_##I>& y,                           \
            ContextPtr loc );                                                       \
                                                                                    \

BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, HARRAY_UTILS_INSTANTIATE, _ )

#undef HARRAY_UTILS2_INSTANTIATE
#undef HARRAY_UTILS_INSTANTIATE

} /* end namespace utilskernel */

} /* end namespace scai */
