/**
 * @file LAMAArrayUtils.cpp
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
 * @brief Definition of class with utility routines.
 * @author Thomas Brandes
 * @date 10.10.2011
 * $Id$
 */

// hpp
#include <lama/LAMAArrayUtils.hpp>

// others
#include <lama/ContextFactory.hpp>
#include <lama/ReadAccess.hpp>
#include <lama/WriteAccess.hpp>
#include <lama/ContextAccess.hpp>
#include <lama/LAMAInterface.hpp>

// tracing
#include <lama/tracing.hpp>

#include <iostream>

namespace lama
{

LAMA_LOG_DEF_LOGGER( LAMAArrayUtils::logger, "LAMAArrayUtils" )

template<typename ValueType1,typename ValueType2>
void LAMAArrayUtils::assignImpl2( LAMAArray<ValueType1>& target, 
                                  const LAMAArray<ValueType2>& source,
                                  const ContextPtr loc )
{
    // verify that dynamic cast operations went okay before

    LAMA_ASSERT_ERROR( &target, "NULL target" )
    LAMA_ASSERT_ERROR( &source, "NULL source" )

    // set should be available on interface for each loc

    LAMA_INTERFACE_FN_TT( set, loc, Utils, Copy, ValueType1, ValueType2 )

    const IndexType n = source.size();

    WriteOnlyAccess<ValueType1> targetVals( target, loc, n );
    ReadAccess<ValueType2> sourceVals( source, loc );

    LAMA_CONTEXT_ACCESS( loc )

    set( targetVals.get(), sourceVals.get(), n );
}

template<typename ValueType>
void LAMAArrayUtils::assignImpl1( LAMAArray<ValueType>& target, 
                                  const _LAMAArray& source,
                                  const ContextPtr loc )
{
    Scalar::ScalarType sourceType = source.getValueType();

    if ( sourceType == target.getValueType() )
    {
        // dynamic cast is absolutely safe

        const LAMAArray<ValueType>& typedSource = dynamic_cast<const LAMAArray<ValueType>&>( source );

        // use assign method of LAMA array, more efficient

        target.assign( typedSource, loc );

        return;
    }

    // Different types -> select for corresponding template routine

    switch ( sourceType )
    {
    case Scalar::FLOAT:
        assignImpl2( target, dynamic_cast<const LAMAArray<float>&>( source ), loc );
        break;
    case Scalar::DOUBLE:
        assignImpl2( target, dynamic_cast<const LAMAArray<double>&>( source ), loc );
        break;
    case Scalar::INDEX_TYPE:
        assignImpl2( target, dynamic_cast<const LAMAArray< IndexType>& >( source ), loc );
        break;
    default:
        LAMA_THROWEXCEPTION( "unsupported source type : " )
    }
}

void LAMAArrayUtils::assign( _LAMAArray& target, const _LAMAArray& source, const ContextPtr loc )
{
    switch ( target.getValueType() )
    {
    case Scalar::FLOAT:
        assignImpl1( dynamic_cast<LAMAArray<float>&>( target ), source, loc );
        break;
    case Scalar::DOUBLE:
        assignImpl1( dynamic_cast<LAMAArray<double>&>( target ), source, loc );
        break;
    case Scalar::INDEX_TYPE:
        assignImpl1( dynamic_cast<LAMAArray< IndexType>& >( target ), source, loc );
        break;
    default:
        LAMA_THROWEXCEPTION( "unsupported target type : " )
    }
}

void LAMAArrayUtils::assign( _LAMAArray& target, const _LAMAArray& source )
{
    // if no context is given we assign where source has a valid copy available

    ContextPtr loc = source.getValidContext( Context::Host );
    assign( target, source, loc );
}

template<typename ValueType1,typename ValueType2>
void LAMAArrayUtils::gather(
    LAMAArray<ValueType1>& target,
    const LAMAArray<ValueType2>& source,
    const LAMAArray<IndexType>& indexes )
{
    LAMA_REGION( "LAMAArray.gather" )

    // choose location for the operation where source array is currently valid

    ContextPtr context = source.getValidContext( Context::Host );

    LAMA_INTERFACE_FN_TT( setGather, context, Utils, Copy, ValueType1, ValueType2 )

    const IndexType n = indexes.size();

    WriteOnlyAccess<ValueType1> wTarget( target, context, n );

    ReadAccess<ValueType2> rSource( source, context );
    ReadAccess<IndexType> rIndexes( indexes, context );

    LAMA_CONTEXT_ACCESS( context )

    // target[i] = source[ indexes[i] ]
    setGather( wTarget.get(), rSource.get(), rIndexes.get(), n );
}

template<typename ValueType>
void LAMAArrayUtils::assign( LAMAArray<ValueType>& target, const Scalar& value, ContextPtr context )
{
    LAMA_LOG_INFO( logger, target << " = " << value )

    // assignment takes place at the given context

    LAMA_INTERFACE_FN_T( setVal, context, Utils, Setter, ValueType )

    const IndexType n = target.size();

    const ValueType val = value.getValue<ValueType>();

    WriteAccess<ValueType> values( target, context );

    LAMA_CONTEXT_ACCESS( context )

    // values[i] = val, i = 0, .., n-1

    setVal( values.get(), n, val );
}

template
void LAMAArrayUtils::gather(
    LAMAArray<int>& target,
    const LAMAArray<int>& source,
    const LAMAArray<IndexType>& indexes );

template
void LAMAArrayUtils::gather(
    LAMAArray<float>& target,
    const LAMAArray<float>& source,
    const LAMAArray<IndexType>& indexes );

template
void LAMAArrayUtils::gather(
    LAMAArray<double>& target,
    const LAMAArray<double>& source,
    const LAMAArray<IndexType>& indexes );

template
void LAMAArrayUtils::gather(
    LAMAArray<float>& target,
    const LAMAArray<double>& source,
    const LAMAArray<IndexType>& indexes );

template
void LAMAArrayUtils::gather(
    LAMAArray<double>& target,
    const LAMAArray<float>& source,
    const LAMAArray<IndexType>& indexes );

template
void LAMAArrayUtils::assign( LAMAArray<float>& target, const Scalar& value, ContextPtr context );

template
void LAMAArrayUtils::assign( LAMAArray<double>& target, const Scalar& value, ContextPtr context );

template
void LAMAArrayUtils::assign( LAMAArray<IndexType>& target, const Scalar& value, ContextPtr context );

} // namespace
