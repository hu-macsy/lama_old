/**
 * @file LAMAArrayRefTest.cpp
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
 * @brief Contains all relevant tests for class LAMAArrayRef
 * @author: Thomas Brandes
 * @date 07.07.2015
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <memory/LAMAArrayRef.hpp>
#include <memory/HostWriteAccess.hpp>
#include <memory/HostReadAccess.hpp>

#include "MockContext.hpp"

using namespace boost;
using namespace common;
using namespace memory;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( LAMAArrayRefTest )

/* --------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( logger, "Test.memory" )

typedef boost::mpl::list<float, double> test_types;

template<typename ValueType>
static void update( ValueType* data, size_t size )
{
    for ( size_t i = 0; i < size; ++i )
    {
        data[i] = 1;
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( refTest, ValueType, test_types )
{
    const IndexType n = 10;
    ValueType myData[10] = {   1, 2, 3, 4, 5, 5, 4, 3, 2, 1};

    {
        // LAMA array keeps myData on Host

        LAMAArrayRef<ValueType> lamaArray( myData, 10 );

        ContextPtr context = Context::create( context::UserContext, 0 );
        ContextPtr host = Context::create( context::Host, -1 );

        BOOST_ASSERT( context );

        {
            WriteAccess<ValueType> write( lamaArray, context );
            update( write.get(), 10 );
        }

        {
            ReadAccess<ValueType> read( lamaArray, host );
        }
    }

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( myData[i], 1 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
