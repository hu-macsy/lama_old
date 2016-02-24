/**
 * @file HArrayUtilsTest.cpp
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
 * @brief Tests for the class HArrayUtils
 * @author: Thomas Brandes
 * @date 22.01.2016
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/test/TestMacros.hpp>
#include <scai/common/ReductionOp.hpp>

using namespace scai::utilskernel;
using namespace scai::hmemo;
using namespace scai::common;

typedef boost::mpl::list<IndexType, float, double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( HArrayUtilsTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.HArrayUtilsTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SetScalarTest, ValueType, test_types )
{
    const IndexType N = 10;

    const ValueType a = 1;
    const ValueType b = 2;
    const ValueType c = 3;

    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();

    HArray<ValueType> array( N );

    HArrayUtils::setScalar( array, a, reduction::COPY, ctx );  // array = a
    HArrayUtils::setScalar( array, b, reduction::ADD, ctx  );  // array += b
    HArrayUtils::setScalar( array, c, reduction::MULT, ctx  ); // array *= c

    ValueType expectedValue = ( a + b ) * c;

    {
        ReadAccess<ValueType> read( array, host );

        for ( IndexType i = 0; i < N; ++i )
        {
            BOOST_CHECK_EQUAL( expectedValue, read[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SetValueTest, ValueType, test_types )
{
    const IndexType N = 10;
    const IndexType k = 3;    // one value in range 0 .. N-1

    const ValueType a = 1;
    const ValueType b = 2;

    ContextPtr ctx  = Context::getContextPtr();
    ContextPtr host = Context::getHostPtr();

    HArray<ValueType> array( N );

    HArrayUtils::setScalar( array, a, reduction::COPY, ctx );  // array = a

    HArrayUtils::setVal( array, k, b );  // array[k] = b

    {
        ReadAccess<ValueType> read( array, host );

        for ( IndexType i = 0; i < N; ++i )
        {
            if ( i == k )
            {
                BOOST_CHECK_EQUAL( b, read[i] );
            }
            else
            {
                BOOST_CHECK_EQUAL( a, read[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
