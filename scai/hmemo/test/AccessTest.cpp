/**
 * @file AccessTest.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Basic tests for LAMA arrays.
 * @author Thomas Brandes
 * @date 08.07.2015
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/assert.hpp>

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/HArrayRef.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/hmemo/WriteOnlyAccess.hpp>
#include <scai/hmemo/ReadAccess.hpp>

#include <functional>

using namespace boost;
using namespace scai;
using namespace scai::hmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( AccessTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.AccessTest" )

/* --------------------------------------------------------------------- */

template<typename ValueType>
void readTest( const ValueType values[], const IndexType N, const ValueType sum )
{
    ValueType mySum = 0;

    for ( IndexType i = 0; i < N; ++i )
    {
        mySum += values[i];
    }

    BOOST_CHECK_EQUAL( sum, mySum );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( releaseTest )
{
    ContextPtr hostContext = Context::getHostPtr();
    HArray<IndexType> hArray; // default, not allocated at all
    ReadAccess<IndexType> readTestAccess( hArray, hostContext );
    readTestAccess.release();
    WriteAccess<IndexType> writeAccess( hArray, hostContext );
    writeAccess.resize( 10 );
    IndexType* data = writeAccess.get();

    for ( IndexType i = 0; i < 10; i++ )
    {
        data[i] = 3;
    }

    writeAccess.release();
    BOOST_CHECK_THROW( { writeAccess.resize( 20 ); }, common::Exception );
    BOOST_CHECK_THROW( { writeAccess.get()[0] = static_cast<IndexType> ( 5.0 ); }, common::Exception );
    ReadAccess<IndexType> readAccess( hArray, hostContext );
    const IndexType* readData = readAccess.get();

    for ( IndexType i = 0; i < 5; i++ )
    {
        BOOST_CHECK_EQUAL( static_cast<IndexType>( 3 ), readData[i] );
    }

    readAccess.release();
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( accessTest )
{
    const IndexType n = 10;
    const double value = 1.0;
    const double value2 = 2.0;
    ContextPtr hostContext = Context::getHostPtr();
    HArray<double> hArray( n, value );
    {
        ReadAccess<double> hArrayRAccess( hArray, hostContext );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value, hArrayRAccess.get()[i] );
        }

        WriteAccess<double> tmpWriteAccess( hArray, hostContext );
    }
    {
        WriteAccess<double> hArrayWAccess( hArray, hostContext );

        for ( IndexType i = 0; i < n; ++i )
        {
            hArrayWAccess.get()[i] = value2;
        }

        hArrayWAccess.release();
        ReadAccess<double> tmpReadAccess( hArray, hostContext );
        ReadAccess<double> hArrayRAccess( hArray, hostContext );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( value2, hArrayRAccess.get()[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( aliasTest )
{
    const IndexType N = 10;
    const double value = 1.0;
    HArray<double> hArray( N, value );
    ContextPtr hostContext = Context::getHostPtr();
    {
        // read and write access at same time by same thread
        ReadAccess<double> read( hArray, hostContext );
        WriteAccess<double> write( hArray, hostContext );

        for ( IndexType i = 0; i < N; ++i )
        {
            write.get()[i] = 2.0 * read.get()[i];
        }
    }
    {
        // verify that operation was really on the same array
        ReadAccess<double> read( hArray, hostContext );
        readTest<double>( read.get(), N, N * value * 2.0 );
    }
    {
        // with a single write access resize is possilbe
        WriteAccess<double> write( hArray, hostContext );
        write.resize( 2 * N );
    }
    {
        // read and write access at same time by same thread
        WriteOnlyAccess<double> write( hArray, hostContext, 2 * N );
        BOOST_CHECK_THROW(
        {
            // read access no more possible as write only did not take care about valid data

            ReadAccess<double> read( hArray, hostContext );
        }, common::Exception );
    }
    {
        // with read and write at the same time resize throws Exception
        ReadAccess<double> read( hArray, hostContext );
        WriteAccess<double> write( hArray, hostContext );
        BOOST_CHECK_THROW(
        {
            write.resize( 3 * N );
        }, common::Exception );
    }
    {
        // read and write access at same time by same thread
        ReadAccess<double> read( hArray, hostContext );
        WriteAccess<double> write( hArray, hostContext );
        // a clear is not possible as it affects the other access
        // Note: clear is the same as resize( 0 )
        BOOST_CHECK_THROW(
        {
            write.clear();
        }, common::Exception );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( delayedReleaseTest )
{
    ContextPtr testContext  = Context::getContextPtr();
    ContextPtr hostContext  = Context::getContextPtr( common::ContextType::Host );

    // host context: we use mock context

    if ( testContext.get() == hostContext.get() )
    {
        testContext  = Context::getContextPtr( common::ContextType::UserContext, 1 );
    }

    SCAI_LOG_INFO( logger, "delayedReleaseTest, test context = " << *testContext )
    HArray<double> X( 10, 5.0 );
    // write access on test Context
    {
        WriteAccess<double> write( X, testContext );
    }
    // read access on Host, okay as write is released
    {
        ReadAccess<double> read( X, hostContext );
    }
    std::function<void()> delay;
    // write access on test context, but delay the release
    {
        WriteAccess<double> write( X, testContext );
        delay = write.releaseDelayed();
        // write access can no more be used
        BOOST_CHECK_THROW(
        {
            double* ptr = write.get();
            write.getMemory().memset( ptr, 0, 1 );
        },
        common::Exception );
    }
    // read access on Host, not possible as release on write is not done
    // this is the reason why we use on host a mock context
    BOOST_CHECK_THROW(
    {
        ReadAccess<double> read( X, hostContext );
    },
    common::Exception  );
    delay();  // now it is okay
    {
        ReadAccess<double> read( X, hostContext );
    }
    // read access on testContext, but delay the release
    {
        ReadAccess<double> read( X, testContext );
        delay = read.releaseDelayed();
    }
    // write access on Host, not possible as release on read is not done
    BOOST_CHECK_THROW(
    {
        WriteAccess<double> write( X, hostContext );
    },
    common::Exception );
    delay();  // now it is okay
    {
        WriteAccess<double> write( X, hostContext );
    }
}

BOOST_AUTO_TEST_CASE( readAccessMoveConstructorTest )
{
    HArray<int> array { 2, 3, -1 };

    ReadAccess<int> read1(array);
    const auto data = read1.get();
    const auto size = read1.size();

    ReadAccess<int> read2(std::move(read1));

    BOOST_TEST(read2.get() == data);
    BOOST_TEST(read2.size() == size);

    // Trying to access the data should trigger an assertion
    // (note: a more specific and documented exception would be preferable here)
    BOOST_CHECK_THROW( read1.get(), scai::common::AssertException);
}

BOOST_AUTO_TEST_CASE( readAccessFunctionTest )
{
    const auto context = Context::getContextPtr();
    HArray<int> array ({ 2, 3, -5, 1 }, context);

    const auto read = readAccess(array);

    BOOST_TEST( read.size() == 4 );
    BOOST_TEST( read[0] ==  2 );
    BOOST_TEST( read[1] ==  3 );
    BOOST_TEST( read[2] == -5 );
    BOOST_TEST( read[3] ==  1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
