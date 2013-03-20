/**
 * @file UtilsTest.cpp
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
 * @brief Contains tests for the class CUDAUtils and OpenMPUtils
 * @author: Jan Ecker
 * @date 19.11.2012
 * $
 **/

// boost
#include <boost/test/unit_test.hpp>

// others
#include <lama/ContextAccess.hpp>
#include <lama/HostReadAccess.hpp>
#include <lama/LAMAArray.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/ReadAccess.hpp>
#include <lama/WriteAccess.hpp>

#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

/* ------------------------------------------------------------------------------------------------------------------ */

// Dummy type, needed to use the lama interface
typedef bool NoType;

/* ------------------------------------------------------------------------------------------------------------------ */

namespace lama
{
namespace UtilsTest
{

template<typename ValueType,typename OtherValueType>
void scaleTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_TT( scale, loc, Utils, Transform, ValueType, OtherValueType );

    ValueType valuesValues[] =
    { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
    const IndexType nValues = sizeof( valuesValues ) / sizeof(ValueType);
    OtherValueType expectedValues[] =
    { 0, 2, 4, 6, 8, 0, 2, 4, 6, 8, 0, 2, 4, 6, 8 };

    const OtherValueType mult = 2.0;

    LAMAArray<ValueType> values( nValues, valuesValues );
    {
        WriteAccess<ValueType> wValues( values, loc );

        LAMA_CONTEXT_ACCESS( loc );

        scale( wValues.get(), nValues, mult );
    }

    HostReadAccess<ValueType> rValues( values );

    for( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expectedValues[i], rValues[i] );
    }

}

template<typename ValueType>
void sumTest( ContextPtr loc )
{
    LAMA_INTERFACE_FN_T( sum, loc, Utils, Reductions, ValueType );

    {
        ValueType valuesValues[] =
        { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
        const IndexType nValues = sizeof( valuesValues ) / sizeof(ValueType);

        const ValueType expectedSum = 30;

        LAMAArray<ValueType> values( nValues, valuesValues );

        ReadAccess<ValueType> rValues( values, loc );

        LAMA_CONTEXT_ACCESS( loc );

        const ValueType resultSum = sum( rValues.get(), nValues );

        BOOST_CHECK_EQUAL( expectedSum, resultSum );
    }

    {
        ValueType valuesValues[] =
            { };
        const IndexType nValues = 0;

        const ValueType expectedSum = 0;

        LAMAArray<ValueType> values( nValues, valuesValues );

        ReadAccess<ValueType> rValues( values, loc );

        LAMA_CONTEXT_ACCESS( loc );

        const ValueType resultSum = sum( rValues.get(), nValues );

        BOOST_CHECK_EQUAL( expectedSum, resultSum );
    }
}

template<typename ValueType>
void setValTest( ContextPtr loc )
{

    LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, ValueType );

    {
        const IndexType n = 20;

        LAMAArray<ValueType> values;

        {
            WriteOnlyAccess<ValueType> wValues( values, loc, n );

            LAMA_CONTEXT_ACCESS( loc );

            setVal( wValues.get(), n, 7 );
        }

        HostReadAccess<ValueType> rValues( values );

        for( IndexType i = 0; i < n; i++ )
        {
            BOOST_CHECK_EQUAL( 7, rValues.get()[i] );
        }
    }

    {
        const IndexType n = 0;

        LAMAArray<ValueType> values;

        {
            WriteOnlyAccess<ValueType> wValues( values, loc, n );

            LAMA_CONTEXT_ACCESS( loc );

            setVal( wValues.get(), n, 7 );
        }
    }
}

template<typename NoType>
void setOrderTest( ContextPtr loc )
{

    LAMA_INTERFACE_FN_T( setOrder, loc, Utils, Setter, IndexType );

    {
        const IndexType n = 20;

        LAMAArray<IndexType> values;

        {
            WriteOnlyAccess<IndexType> wValues( values, loc, n );

            LAMA_CONTEXT_ACCESS( loc );

            setOrder( wValues.get(), n );
        }

        HostReadAccess<IndexType> rValues( values );

        for( IndexType i = 0; i < n; i++ )
        {
            BOOST_CHECK_EQUAL( i, rValues.get()[i] );
        }
    }

    {
        const IndexType n = 0;

        LAMAArray<IndexType> values;

        {
            WriteOnlyAccess<IndexType> wValues( values, loc, n );

            LAMA_CONTEXT_ACCESS( loc );

            setOrder( wValues.get(), n );
        }
    }

}

template<typename ValueType>
void invertTest( ContextPtr loc )
{

    LAMA_INTERFACE_FN_T( invert, loc, Utils, Math, ValueType );

    {
        // TODO: should it be possible to pass 0 elements? What should be the result?
        ValueType valuesValues[] =
        { 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4 };
        const IndexType nValues = sizeof( valuesValues ) / sizeof(ValueType);

        LAMAArray<ValueType> values( nValues, valuesValues );

        {
            WriteAccess<ValueType> wValues( values, loc );

            LAMA_CONTEXT_ACCESS( loc );

            invert( wValues.get(), nValues );
        }

        HostReadAccess<ValueType> rValues( values );

        for( IndexType i = 0; i < nValues; i++ )
        {
            BOOST_CHECK_CLOSE( 1 / valuesValues[i], rValues.get()[i], 1 );
        }
    }

    {
        const IndexType n = 0;

        LAMAArray<ValueType> values;

        {
            WriteOnlyAccess<ValueType> wValues( values, loc, n );

            LAMA_CONTEXT_ACCESS( loc );

            invert( wValues.get(), n );
        }
    }

}

// TODO: add SPMV tests

}//namespace UtilsTest
} //namespace lama

/* ------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( UtilsTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.UtilsTest" );

LAMA_AUTO_TEST_CASE_T( sumTest, UtilsTest );
LAMA_AUTO_TEST_CASE_T( setValTest, UtilsTest );
LAMA_AUTO_TEST_CASE_T( invertTest, UtilsTest );

LAMA_AUTO_TEST_CASE_TDUMMY( setOrderTest, UtilsTest );

LAMA_AUTO_TEST_CASE_TT( scaleTest, UtilsTest );

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END();
