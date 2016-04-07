/**
 * @file LArrayTest.cpp
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
 * @brief Tests for the class LArray
 * @author: Thomas Brandes
 * @date 22.02.2016
 **/

#include <boost/test/unit_test.hpp>

#include <scai/utilskernel/LArray.hpp>
#include <scai/utilskernel/LArray.hpp>

#include <scai/common/test/TestMacros.hpp>

using namespace scai::utilskernel;
using namespace scai::hmemo;
using namespace scai::common;

extern ContextPtr testContext;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( LArrayTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.LArrayTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( indexTest )
{
    testContext = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "indexTest on " << *testContext )

    // the LArray allows indexed access, but attention: can be very slow

    LArray<IndexType> array( testContext );

    const IndexType N = 10;

    array.resize( N );

    for ( IndexType i = 0; i < N; ++i )
    {
        array[i] = i;
    }

    for ( IndexType i = 0; i < N; ++i )
    {
        BOOST_CHECK_EQUAL( i, array[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( assignTest )
{
    testContext = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "assignTest on " << *testContext )

    // the LArray allows indexed access, but attention: can be very slow

    LArray<IndexType> array( testContext );

    const IndexType N = 10;
    const IndexType val = 3;

    array.resize( N );

    array = val;

    for ( IndexType i = 0; i < N; ++i )
    {
        BOOST_CHECK_EQUAL( val, array[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    testContext = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "constructorTest on " << *testContext )

    // the LArray allows indexed access, but attention: can be very slow

    const IndexType N = 10;
    const IndexType val = 3;

    LArray<IndexType> array( N, val, testContext );

    BOOST_CHECK( array.isValid( testContext ) );

    for ( IndexType i = 0; i < N; ++i )
    {
        BOOST_CHECK_EQUAL( val, array[i] );
    }

    const IndexType myVals[N] = { 1, 5, 9, 4, 6, 3, 7, 8, 0, 2 };

    LArray<IndexType> array1( N, myVals, testContext );

    BOOST_CHECK( array.isValid( testContext ) );

    for ( IndexType i = 0; i < N; ++i )
    {
        BOOST_CHECK_EQUAL( myVals[i], array1[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( reductionTest, ValueType, scai_array_test_types )
{
    testContext = Context::getContextPtr();

    SCAI_LOG_INFO( logger, "reductionTest on " << *testContext )

    // the LArray allows indexed access, but attention: can be very slow

    const ValueType myVals[] = { 9, 5, 1, 4, 6, 3, 7, 8, 2, 0 };

    const IndexType N = sizeof( myVals ) / sizeof( ValueType );

    LArray<ValueType> array( N, myVals, testContext );

    BOOST_CHECK( array.isValid( testContext ) );

    BOOST_CHECK_EQUAL( 0, array.min() );
    BOOST_CHECK_EQUAL( 9, array.max() );
    BOOST_CHECK_EQUAL( 45, array.sum() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
