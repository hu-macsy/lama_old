/**
 * @file LAMAArrayUtils.cpp
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
#include <scai/lama/LAMAArrayUtils.hpp>

// local library
#include <scai/lama/UtilKernelTrait.hpp>
#include <scai/lama/LAMAKernel.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/Constants.hpp>

// boost
#include <boost/preprocessor.hpp>

// std
#include <iostream>

using namespace scai::hmemo;

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( LAMAArrayUtils::logger, "LAMAArrayUtils" )

template<typename ValueType1,typename ValueType2>
void LAMAArrayUtils::assignImpl(
    LAMAArray<ValueType1>& target,
    const LAMAArray<ValueType2>& source,
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

    WriteOnlyAccess<ValueType1> targetVals( target, loc, n );
    ReadAccess<ValueType2> sourceVals( source, loc );

    // Implemenation of set @ loc is available

    set[loc]( targetVals.get(), sourceVals.get(), n );
}

template<typename ValueType>
void LAMAArrayUtils::assignImpl1( LAMAArray<ValueType>& target, const ContextArray& source, const ContextPtr loc )
{
    common::scalar::ScalarType sourceType = source.getValueType();

    if( sourceType == target.getValueType() )
    {
        // dynamic cast is absolutely safe

        const LAMAArray<ValueType>& typedSource = dynamic_cast<const LAMAArray<ValueType>&>( source );

        // use assign method of LAMA array, more efficient

        target.assign( typedSource, loc );

        return;
    }

    // Different types -> select for corresponding template routine

    switch( sourceType )
    {

        case common::scalar::INDEX_TYPE:
            assignImpl( target, dynamic_cast<const LAMAArray<IndexType>&>( source ), loc );
            break;

#define LAMA_ARRAY_ASSIGN( z, I, _ )                                                                \
case SCALAR_ARITHMETIC_TYPE##I:                                                                     \
    assignImpl( target, dynamic_cast<const LAMAArray<ARITHMETIC_HOST_TYPE_##I>& >( source ), loc ); \
    break;                                                                                          \

        BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_ARRAY_ASSIGN, _ )

#undef LAMA_ARRAY_ASSIGN

default        :
        COMMON_THROWEXCEPTION( "unsupported source type : " )
    }
}

void LAMAArrayUtils::assign( ContextArray& target, const ContextArray& source, const ContextPtr loc /* = ContextPtr() */)
{
    ContextPtr validLoc = loc;

    if( !validLoc )
    {
// if no context is given we assign where source has a valid copy available

        validLoc = source.getValidContext();
    }

    switch ( target.getValueType() )
    {
        case common::scalar::INDEX_TYPE:
            assignImpl1( dynamic_cast<LAMAArray<IndexType>&>( target ), source, validLoc );
            break;

#define LAMA_ARRAY_ASSIGN1( z, I, _ )                                                                \
case SCALAR_ARITHMETIC_TYPE##I:                                                                      \
    assignImpl1( dynamic_cast<LAMAArray< ARITHMETIC_HOST_TYPE_##I>& >( target ), source, validLoc ); \
    break;

        BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_ARRAY_ASSIGN1, _ )

#undef LAMA_ARRAY_ASSIGN1

default        :
        COMMON_THROWEXCEPTION( "unsupported target type : " << target.getValueType() )
    }
}

template<typename ValueType1,typename ValueType2>
void LAMAArrayUtils::gather(
    LAMAArray<ValueType1>& target,
    const LAMAArray<ValueType2>& source,
    const LAMAArray<IndexType>& indexes )
{
    SCAI_REGION( "LAMAArray.gather" )

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
void LAMAArrayUtils::assignScalar( LAMAArray<ValueType>& target, const Scalar& value, ContextPtr prefContext )
{
    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

    ContextPtr context = setVal.getValidContext( prefContext );

    SCAI_LOG_INFO( logger, target << " = " << value << ", to do at " << *context )

    // assignment takes place at the specified context, no check here

    const IndexType n = target.size();

    const ValueType val = value.getValue<ValueType>();

    WriteOnlyAccess<ValueType> values( target, context );

    SCAI_CONTEXT_ACCESS( context )

    setVal[context]( values.get(), n, val );
}

void LAMAArrayUtils::assignScalar( hmemo::ContextArray& target, const Scalar& value, hmemo::ContextPtr context )
{
    common::scalar::ScalarType arrayType = target.getValueType();

    switch( arrayType )
    {
        case common::scalar::INDEX_TYPE:
        {
            LAMAArray<IndexType>& typedTarget = dynamic_cast<LAMAArray<IndexType>&>( target );
            assignScalar( typedTarget, value, context );
            break;
        }

            // for all supported arithmetic types generate it

#define LAMA_ARRAY_ASSIGN_SCALAR( z, I, _ )                                   \
case SCALAR_ARITHMETIC_TYPE##I:                                               \
{                                                                             \
    LAMAArray<ARITHMETIC_HOST_TYPE_##I>& typedTarget =                        \
            dynamic_cast<LAMAArray<ARITHMETIC_HOST_TYPE_##I>&>( target );     \
    assignScalar( typedTarget, value, context );                              \
    break;                                                                    \
}

BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_ARRAY_ASSIGN_SCALAR, _ )

#undef LAMA_ARRAY_ASSIGN_SCALAR

default        :
        {
            COMMON_THROWEXCEPTION( target << ": assignScalar for value type " << arrayType << " not supported" )
        }
    }
}

template<typename ValueType>
void LAMAArrayUtils::setVal( LAMAArray<ValueType>& target, const IndexType index, ValueType val )
{
    SCAI_ASSERT_DEBUG( index < target.size(), "index = " << index << " out of range for target = " << target );

    ContextPtr loc = target.getValidContext();   // preferred location where to fill

    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

    loc = setVal.getValidContext( loc );

    WriteAccess<ValueType> wTarget( target, loc );

    SCAI_CONTEXT_ACCESS( loc )

    setVal[loc]( wTarget.get() + index, 1, val );
}

template<typename ValueType>
void LAMAArrayUtils::assignScaled(
    LAMAArray<ValueType>& result,
    const ValueType beta,
    const LAMAArray<ValueType>& y,
    ContextPtr prefLoc )
{
    const IndexType n = result.size();

    SCAI_ASSERT_EQUAL_ERROR( n, y.size() );

    // beta = 0    : saves the need of a read access for y
    // result == y : only one write access needed ( write + read not possible)

    if( beta == scai::common::constants::ZERO )
    {
        // result := 0

        static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

        ContextPtr loc = setVal.getValidContext( prefLoc );

        WriteAccess<ValueType> wResult( result, loc );

        SCAI_CONTEXT_ACCESS( loc )

        setVal[loc]( wResult.get(), n, static_cast<ValueType>(0.0) );
    }
    else if( &result == &y )
    {
        if( beta == scai::common::constants::ONE )
        {
            return;
        }

        // result := beta * result, is result *= beta

        static LAMAKernel<UtilKernelTrait::scale<ValueType> > scale;

        ContextPtr loc = scale.getValidContext( prefLoc );

        WriteAccess<ValueType> wResult( result, loc );

        SCAI_CONTEXT_ACCESS( loc )

        scale[loc]( wResult.get(), beta, n );
    }
    else
    {
        // result := beta * y
        
        // Note: we do not use BLAS1:axpy here to guarantee same LAMA OpenMP schedule
        //       and to support type conversions in place for multiprecision support
        
        static LAMAKernel<UtilKernelTrait::setScale<ValueType, ValueType> > setScale;

        ContextPtr loc = setScale.getValidContext( prefLoc );

        SCAI_CONTEXT_ACCESS( loc )

        WriteAccess<ValueType> wResult( result, loc );
        ReadAccess<ValueType> rY( y, loc );

        setScale[loc]( wResult.get(), beta, rY.get(), n );
    }
}

// template instantiation for the supported data types

template void LAMAArrayUtils::setVal( hmemo::LAMAArray<IndexType>& , const IndexType , IndexType );

template void LAMAArrayUtils::assignScaled(
    hmemo::LAMAArray<IndexType>& ,
    const IndexType ,
    const LAMAArray<IndexType>& ,
    hmemo::ContextPtr  );

template void LAMAArrayUtils::gather(
    hmemo::LAMAArray<IndexType>& ,
    const hmemo::LAMAArray<IndexType>& ,
    const hmemo::LAMAArray<IndexType>& );

/** Macro instantiates operations that have also type conversion */

#define LAMA_ARRAY_UTILS2_INSTANTIATE(z, J, TYPE)                                   \
    template                                                                        \
    void LAMAArrayUtils::gather(                                                    \
            LAMAArray<TYPE>& target,                                                \
            const LAMAArray<ARITHMETIC_HOST_TYPE_##J>& source,                      \
            const LAMAArray<IndexType>& indexes );                                  \

/** Macro instantiates operations for supported arithmetic types */

#define LAMA_ARRAY_UTILS_INSTANTIATE(z, I, _)                                       \
    template                                                                        \
    void LAMAArrayUtils::setVal(                                                    \
            LAMAArray<ARITHMETIC_HOST_TYPE_##I>& target,                            \
            const IndexType index,                                                  \
            ARITHMETIC_HOST_TYPE_##I val );                                         \
                                                                                    \
    template                                                                        \
    void LAMAArrayUtils::assignScaled(                                              \
            LAMAArray<ARITHMETIC_HOST_TYPE_##I>& result,                            \
            const ARITHMETIC_HOST_TYPE_##I beta,                                    \
            const LAMAArray<ARITHMETIC_HOST_TYPE_##I>& y,                           \
            ContextPtr loc );                                                       \
                                                                                    \
    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT,                                      \
                     LAMA_ARRAY_UTILS2_INSTANTIATE, ARITHMETIC_HOST_TYPE_##I )

BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_ARRAY_UTILS_INSTANTIATE, _ )

#undef LAMA_ARRAY_UTILS2_INSTANTIATE
#undef LAMA_ARRAY_UTILS_INSTANTIATE

} /* end namespace lama */

} /* end namespace scai */
