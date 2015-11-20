/**
 * @file UtilsTest.cpp
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
 * @brief Contains tests for the class CUDAUtils and OpenMPUtils
 * @author: Jan Ecker
 * @date 19.11.2012
 * @since 1.0.0
 **/

// boost
#include <boost/test/unit_test.hpp>

// others
#include <scai/lama/LAMAKernel.hpp>
#include <scai/lama/LAMAArray.hpp>
#include <scai/lama/UtilKernelTrait.hpp>
#include <scai/hmemo.hpp>

#include <scai/common/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

/* ------------------------------------------------------------------------------------------------------------------ */

// Dummy type, needed to use the lama interface
typedef bool NoType;

/* ------------------------------------------------------------------------------------------------------------------ */

namespace scai
{

namespace lama
{

namespace UtilsTest
{

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void scaleTest( ContextPtr loc )
{
    static LAMAKernel<UtilKernelTrait::scale<ValueType> > scale;

    ValueType valuesValues[] =
    { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
    const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
    ValueType expectedValues[] =
    { 0, 2, 4, 6, 8, 0, 2, 4, 6, 8, 0, 2, 4, 6, 8 };
    const ValueType mult = 2.0;
    LAMAArray<ValueType> values( nValues, valuesValues );
    {
        WriteAccess<ValueType> wValues( values, loc );
        SCAI_CONTEXT_ACCESS( loc );
        scale[loc]( wValues.get(), mult, nValues );
    }
    ReadAccess<ValueType> rValues( values );

    for ( IndexType i = 0; i < nValues; i++ )
    {
        BOOST_CHECK_EQUAL( expectedValues[i], rValues[i] );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void sumTest( ContextPtr loc )
{
    static LAMAKernel<UtilKernelTrait::sum<ValueType> > sum;
    {
        ValueType valuesValues[] =
        { 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4 };
        const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
        const ValueType expectedSum = 30;
        LAMAArray<ValueType> values( nValues, valuesValues );
        ReadAccess<ValueType> rValues( values, loc );
        SCAI_CONTEXT_ACCESS( loc );
        const ValueType resultSum = sum[loc]( rValues.get(), nValues );
        BOOST_CHECK_EQUAL( expectedSum, resultSum );
    }
    {
        const ValueType expectedSum = 0;
        LAMAArray<ValueType> values;
        ReadAccess<ValueType> rValues( values, loc );
        SCAI_CONTEXT_ACCESS( loc );
        const ValueType resultSum = sum[loc]( rValues.get(), values.size() );
        BOOST_CHECK_EQUAL( expectedSum, resultSum );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void setValTest( ContextPtr loc )
{
    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;
    {
        const IndexType n = 20;
        LAMAArray<ValueType> values;
        {
            WriteOnlyAccess<ValueType> wValues( values, loc, 3 * n );
            SCAI_CONTEXT_ACCESS( loc );
            setVal[loc]( wValues.get(), 3 * n, 0 );
            // overwrite in the middle to check that there is no out-of-range set
            setVal[loc]( wValues.get() + n, n, 10 );
        }
        ReadAccess<ValueType> rValues( values );

        for ( IndexType i = 0; i < n; i++ )
        {
            BOOST_CHECK_EQUAL( 0, rValues.get()[i + 0 * n] );
            BOOST_CHECK_EQUAL( 10, rValues.get()[i + 1 * n] );
            BOOST_CHECK_EQUAL( 0, rValues.get()[i + 2 * n] );
        }
    }
    {
        const IndexType n = 0;
        LAMAArray<ValueType> values;
        {
            WriteOnlyAccess<ValueType> wValues( values, loc, n );
            SCAI_CONTEXT_ACCESS( loc );
            setVal[loc]( wValues.get(), n, 7 );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void isSortedTest( ContextPtr loc )
{
    static LAMAKernel<UtilKernelTrait::isSorted<ValueType> > isSorted;
    {
        ValueType values1[] =
        { 1, 2, 2, 2, 5, 8 };
        ValueType values2[] =
        { 2, 2, 1, 0 };
        ValueType values3[] =
        { 1, 0, 1 };
        const IndexType nValues1 = sizeof( values1 ) / sizeof( ValueType );
        const IndexType nValues2 = sizeof( values2 ) / sizeof( ValueType );
        const IndexType nValues3 = sizeof( values3 ) / sizeof( ValueType );
        LAMAArray<ValueType> valueArray1( nValues1, values1 );
        LAMAArray<ValueType> valueArray2( nValues2, values2 );
        LAMAArray<ValueType> valueArray3( nValues3, values3 );
        ReadAccess<ValueType> rValues1( valueArray1, loc );
        ReadAccess<ValueType> rValues2( valueArray2, loc );
        ReadAccess<ValueType> rValues3( valueArray3, loc );
        SCAI_CONTEXT_ACCESS( loc );
        // values1 are sorted, ascending = true
        BOOST_CHECK( isSorted[loc]( rValues1.get(), nValues1, true ) );
        BOOST_CHECK( ! isSorted[loc]( rValues1.get(), nValues1, false ) );
        // values2 are sorted, ascending = false
        BOOST_CHECK( isSorted[loc]( rValues2.get(), nValues2, false ) );
        BOOST_CHECK( ! isSorted[loc]( rValues2.get(), nValues2, true ) );
        BOOST_CHECK( isSorted[loc]( rValues2.get(), 0, true ) );
        // only first two values are sorted
        BOOST_CHECK( isSorted[loc]( rValues2.get(), 1, true ) );
        // only first two values are sorted
        BOOST_CHECK( isSorted[loc]( rValues2.get(), 2, true ) );
        // only first two values are sorted
        // values3 are not sorted, neither ascending nor descending
        BOOST_CHECK( ! isSorted[loc]( rValues3.get(), nValues3, false ) );
        BOOST_CHECK( ! isSorted[loc]( rValues3.get(), nValues3, true ) );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename NoType>
void setOrderTest( ContextPtr loc )
{
    static LAMAKernel<UtilKernelTrait::setOrder<IndexType> > setOrder;
    {
        const IndexType n = 20;
        LAMAArray<IndexType> values;
        {
            WriteOnlyAccess<IndexType> wValues( values, loc, n );
            SCAI_CONTEXT_ACCESS( loc );
            setOrder[loc]( wValues.get(), n );
        }
        ReadAccess<IndexType> rValues( values );

        for ( IndexType i = 0; i < n; i++ )
        {
            BOOST_CHECK_EQUAL( i, rValues.get()[i] );
        }
    }
    {
        const IndexType n = 0;
        LAMAArray<IndexType> values;
        {
            WriteOnlyAccess<IndexType> wValues( values, loc, n );
            SCAI_CONTEXT_ACCESS( loc );
            setOrder[loc]( wValues.get(), n );
        }
    }
}

template<typename ValueType>
void invertTest( ContextPtr loc )
{
    static LAMAKernel<UtilKernelTrait::invert<ValueType> > invert;
    {
        // TODO: should it be possible to pass 0 elements? What should be the result?
        ValueType valuesValues[] =
        { 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4 };
        const IndexType nValues = sizeof( valuesValues ) / sizeof( ValueType );
        LAMAArray<ValueType> values( nValues, valuesValues );
        {
            WriteAccess<ValueType> wValues( values, loc );
            SCAI_CONTEXT_ACCESS( loc );
            invert[loc]( wValues.get(), nValues );
        }
        ReadAccess<ValueType> rValues( values );

        for ( IndexType i = 0; i < nValues; i++ )
        {
            SCAI_CHECK_CLOSE( 1 / valuesValues[i], rValues.get()[i], 1 );
        }
    }
    {
        const IndexType n = 0;
        LAMAArray<ValueType> values;
        {
            WriteOnlyAccess<ValueType> wValues( values, loc, n );
            SCAI_CONTEXT_ACCESS( loc );
            invert[loc]( wValues.get(), n );
        }
    }
}

// TODO: add SPMV tests

}//namespace UtilsTest

} /* end namespace lama */

} /* end namespace scai */

/* ------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( UtilsTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.UtilsTest" )

LAMA_AUTO_TEST_CASE_CT( sumTest, UtilsTest, scai::lama )
LAMA_AUTO_TEST_CASE_CT( isSortedTest, UtilsTest, scai::lama )
LAMA_AUTO_TEST_CASE_CT( setValTest, UtilsTest, scai::lama )
LAMA_AUTO_TEST_CASE_CT( invertTest, UtilsTest, scai::lama )

LAMA_AUTO_TEST_CASE_CTDUMMY( setOrderTest, UtilsTest )

LAMA_AUTO_TEST_CASE_CT( scaleTest, UtilsTest, scai::lama )

/* ------------------------------------------------------------------------------------------------------------------ */BOOST_AUTO_TEST_SUITE_END()
