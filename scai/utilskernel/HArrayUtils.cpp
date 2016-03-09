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

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/preprocessor.hpp>

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

template<typename ValueType>
void HArrayUtils::setImpl1( hmemo::HArray<ValueType>& target, const hmemo::_HArray& source, const common::reduction::ReductionOp op, const hmemo::ContextPtr loc )
{
    common::scalar::ScalarType sourceType = source.getValueType();

    // Different types -> select for corresponding template routine

    switch ( sourceType )
    {
        case common::scalar::INDEX_TYPE:
            setImpl( target, dynamic_cast<const hmemo::HArray<IndexType>&>( source ), op, loc );
            break;

        // for all the other arithmetic types we use BOOST_PP_REPEAT to loop over supported types

#define HARRAY_SET( z, I, _ )                                                                             \
                                                                                                          \
        case common::TypeTraits<ARITHMETIC_HOST_TYPE_##I>::stype :                                        \
            setImpl( target, dynamic_cast<const hmemo::HArray<ARITHMETIC_HOST_TYPE_##I>& >( source ), op, loc ); \
            break;                                                                                        \

        BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, HARRAY_SET, _ )

#undef HARRAY_SET

default        :
        COMMON_THROWEXCEPTION( "unsupported source type : " )
    }
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

    switch ( target.getValueType() )
    {
        case common::scalar::INDEX_TYPE:
            setImpl1( dynamic_cast<hmemo::HArray<IndexType>&>( target ), source, op, validLoc );
            break;

        // for all the other arithmetic types we use BOOST_PP_REPEAT to loop over supported types

#define HARRAY_SET1( z, I, _ )                                                                                    \
        case common::TypeTraits<ARITHMETIC_HOST_TYPE_##I>::stype :                                                \
            setImpl1( dynamic_cast<hmemo::HArray< ARITHMETIC_HOST_TYPE_##I>& >( target ), source, op, validLoc ); \
            break;

        BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, HARRAY_SET1, _ )

#undef HARRAY_SET1

default        :
        COMMON_THROWEXCEPTION( "unsupported target type : " << target.getValueType() )
    }
}

template<typename ValueType1,typename ValueType2>
void HArrayUtils::gather(
    hmemo::HArray<ValueType1>& target,
    const hmemo::HArray<ValueType2>& source,
    const hmemo::HArray<IndexType>& indexes )
{
    SCAI_REGION( "HArray.gather" )

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
void HArrayUtils::setScalar( hmemo::HArray<ValueType>& target, const ValueType value, const common::reduction::ReductionOp op, hmemo::ContextPtr prefContext )
{
    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

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

/*
void HArrayUtils::setScalar( hmemo::_HArray& target, const Scalar& value, const common::reduction::ReductionOp op, hmemo::ContextPtr context )
{
    common::scalar::ScalarType arrayType = target.getValueType();

    switch( arrayType )
    {
        case common::scalar::INDEX_TYPE :
        {
            hmemo::HArray<IndexType>& typedTarget = dynamic_cast<HArray<IndexType>&>( target );
            setScalar( typedTarget, value.getValue<IndexType>(), op, context );
            break;
        }

        // for all supported arithmetic types generate it

#define HARRAY_SET_SCALAR( z, I, _ )                                                              \
        case common::TypeTraits<ARITHMETIC_HOST_TYPE_##I>::stype :                                \
        {                                                                                         \
            hmemo::HArray<ARITHMETIC_HOST_TYPE_##I>& typedTarget =                                \
                    dynamic_cast<hmemo::HArray<ARITHMETIC_HOST_TYPE_##I>&>( target );             \
            setScalar( typedTarget, value.getValue<ARITHMETIC_HOST_TYPE_##I>(), op, context );    \
            break;                                                                                \
        }
        
        BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, HARRAY_SET_SCALAR, _ )

#undef HARRAY_SET_SCALAR

        default :
        {
            COMMON_THROWEXCEPTION( target << ": setScalar for value type " << arrayType << " not supported" )
        }
    }
}
*/

template<typename ValueType>
void HArrayUtils::setVal( hmemo::HArray<ValueType>& target, const IndexType index, ValueType val )
{
    SCAI_ASSERT_DEBUG( index < target.size(), "index = " << index << " out of range for target = " << target );

    hmemo::ContextPtr loc = target.getValidContext();   // preferred location where to fill

    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

    loc = setVal.getValidContext( loc );

    SCAI_LOG_INFO( logger, "setVal<" << common::TypeTraits<ValueType>::id() << ">[" << index 
                           << "] = " << val << " @ " << *loc )

    hmemo::WriteAccess<ValueType> wTarget( target, loc );

    SCAI_CONTEXT_ACCESS( loc )

    setVal[loc]( wTarget.get() + index, 1, val, common::reduction::COPY );
}

template<typename ValueType>
ValueType HArrayUtils::getVal( const hmemo::HArray<ValueType>& array, const IndexType index )
{
    SCAI_ASSERT_DEBUG( index < array.size(), "index = " << index << " out of range for array = " << array );

    // get the data from a valid context, avoids any memory copy.

    hmemo::ContextPtr loc = array.getValidContext();

    static LAMAKernel<UtilKernelTrait::getValue<ValueType> > getValue;

    loc = getValue.getValidContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    hmemo::ReadAccess<ValueType> rArray( array, loc );

    return getValue[loc]( rArray.get(), index );
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

        static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

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

        static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

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

// template instantiation for the supported data types

template void HArrayUtils::setVal( hmemo::HArray<IndexType>& , const IndexType , IndexType );

template void HArrayUtils::setScalar( hmemo::HArray<IndexType>& target, const IndexType value,
                                      const common::reduction::ReductionOp op, hmemo::ContextPtr prefContext );

template IndexType HArrayUtils::getVal( const hmemo::HArray<IndexType>& , const IndexType );

template void HArrayUtils::assignScaled(
    hmemo::HArray<IndexType>& ,
    const IndexType ,
    const hmemo::HArray<IndexType>& ,
    hmemo::ContextPtr  );

template void HArrayUtils::gather(
    hmemo::HArray<IndexType>& ,
    const hmemo::HArray<IndexType>& ,
    const hmemo::HArray<IndexType>& );

/** Macro instantiates operations that have also type conversion */

#define HARRAY_UTILS2_INSTANTIATE(z, J, TYPE)                                    \
    template                                                                     \
    void HArrayUtils::gather(                                                    \
            hmemo::HArray<TYPE>& target,                                         \
            const hmemo::HArray<ARITHMETIC_HOST_TYPE_##J>& source,               \
            const hmemo::HArray<IndexType>& indexes );                           \

/** Macro instantiates operations for supported arithmetic types */

#define HARRAY_UTILS_INSTANTIATE(z, I, _)                                           \
    template                                                                        \
    void HArrayUtils::setVal(                                                       \
            hmemo::HArray<ARITHMETIC_HOST_TYPE_##I>& target,                        \
            const IndexType index,                                                  \
            ARITHMETIC_HOST_TYPE_##I val );                                         \
                                                                                    \
    template                                                                        \
    void HArrayUtils::setScalar( hmemo::HArray<ARITHMETIC_HOST_TYPE_##I>& target,   \
            const ARITHMETIC_HOST_TYPE_##I value,                                   \
            const common::reduction::ReductionOp op,                                \
            hmemo::ContextPtr prefContext );                                        \
    template                                                                        \
    void HArrayUtils::scale(                                                        \
            hmemo::HArray<ARITHMETIC_HOST_TYPE_##I>& array,                         \
            const ARITHMETIC_HOST_TYPE_##I beta,                                    \
            hmemo::ContextPtr prefLoc );                                            \
                                                                                    \
    template                                                                        \
    void HArrayUtils::conj(                                                         \
            hmemo::HArray<ARITHMETIC_HOST_TYPE_##I>& target,                        \
            hmemo::ContextPtr loc );                                                \
                                                                                    \
    template                                                                        \
    ARITHMETIC_HOST_TYPE_##I HArrayUtils::getVal(                                   \
            const hmemo::HArray<ARITHMETIC_HOST_TYPE_##I>&,                         \
            const IndexType );                                                      \
                                                                                    \
    template                                                                        \
    void HArrayUtils::assignScaled(                                                 \
            hmemo::HArray<ARITHMETIC_HOST_TYPE_##I>& result,                        \
            const ARITHMETIC_HOST_TYPE_##I beta,                                    \
            const hmemo::HArray<ARITHMETIC_HOST_TYPE_##I>& y,                       \
            hmemo::ContextPtr loc );                                                \
                                                                                    \
    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT,                                      \
                     HARRAY_UTILS2_INSTANTIATE, ARITHMETIC_HOST_TYPE_##I )

BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, HARRAY_UTILS_INSTANTIATE, _ )

#undef HARRAY_UTILS2_INSTANTIATE
#undef HARRAY_UTILS_INSTANTIATE

} /* end namespace utilskernel */

} /* end namespace scai */
