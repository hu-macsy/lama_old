/**
 * @file HArrayRefTest.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Contains all relevant tests for class HArrayRef
 * @author Thomas Brandes
 * @date 07.07.2015
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/hmemo.hpp>

using namespace boost;
using namespace scai;
using namespace common;
using namespace hmemo;

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
    ContextPtr mock = Context::create( ContextType::UserContext, 0 );
    ContextPtr host = Context::create( ContextType::Host, -1 );
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
            Exception );
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
                WriteAccess<ValueType> hArrayWAccess( hArray, host );
            }
            , Exception );
        }
    }
}

/* --------------------------------------------------------------------- */


/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
