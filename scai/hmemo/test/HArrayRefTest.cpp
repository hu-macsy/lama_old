/**
 * @file HArrayRefTest.cpp
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
 * @brief Contains all relevant tests for class HArrayRef
 * @author: Thomas Brandes
 * @date 07.07.2015
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/hmemo.hpp>

using namespace boost;
using namespace scai::common;
using namespace scai::hmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( HArrayRefTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.memory" )

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

template<typename ValueType>
static void update( ValueType* data, size_t size )
{
    SCAI_LOG_DEBUG( logger, "update at data = " << data )

    for ( size_t i = 0; i < size; ++i )
    {
        data[i] = 1;
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( refTest, ValueType, test_types )
{
    ContextPtr mock = Context::create( context::UserContext, 0 );
    ContextPtr host = Context::create( context::Host, -1 );

    BOOST_CHECK( mock );
    BOOST_CHECK( host );

    const IndexType N = 10;

    ValueType myData[ N ] = { 1, 2, 3, 4, 5, 5, 4, 3, 2, 1};

    const ValueType* myConstData = myData;

    SCAI_LOG_INFO( logger, "myData = " << myData )

    {
        // hArray keeps myData on Host

        HArrayRef<ValueType> hArray( N, myData );

        SCAI_LOG_INFO( logger, "hArray = " << hArray )

        {
            // modify the data @ mock

            WriteAccess<ValueType> write( hArray, mock );
            update( write.get(), N );
        }

        {
            WriteAccess<ValueType> write( hArray, mock );
        }

        SCAI_LOG_INFO( logger, "modified at userContext: hArray = " << hArray )

        {
            // get valid data back @ host

            ReadAccess<ValueType> read( hArray, host );
        }

        {
            WriteAccess<ValueType> write( hArray, host );
            BOOST_CHECK_THROW(
            {
                write.resize( 2 * N );
            },
            Exception )
        }

        for ( IndexType i = 0; i < N; ++i )
        {
            BOOST_CHECK_EQUAL( myData[i], 1 );
        }

        {
            // this will create a LAMA array with a const reference,

            HArrayRef<ValueType> hArray( N, myConstData );

            BOOST_CHECK_EQUAL( hArray.size(), static_cast<IndexType>( N ) );
            // Write access should not be allowed
            BOOST_CHECK_THROW(
            {
                WriteAccess<ValueType> hArrayWAccess( hArray, host);
            }
            , Exception );
        }
    }
}

/* --------------------------------------------------------------------- */


/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
