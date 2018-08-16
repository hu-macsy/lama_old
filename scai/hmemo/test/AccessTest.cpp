/**
 * @file AccessTest.cpp
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
#include <scai/hmemo/HostReadAccess.hpp>
#include <scai/hmemo/HostWriteAccess.hpp>
#include <scai/hmemo/HostWriteOnlyAccess.hpp>

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

    ReadAccess<int> read1( array );
    const auto data = read1.get();
    const auto size = read1.size();

    ReadAccess<int> read2( std::move( read1 ) );

    BOOST_TEST( read2.get() == data );
    BOOST_TEST( read2.size() == size );

    // Trying to access the data should trigger an assertion
    // (note: a more specific and documented exception would be preferable here)
    BOOST_CHECK_THROW( read1.get(), scai::common::AssertException );
}

BOOST_AUTO_TEST_CASE( readAccessFunctionTest )
{
    const auto context = Context::getContextPtr();
    HArray<int> array ( { 2, 3, -5, 1 }, context );

    const auto read = readAccess( array );

    BOOST_TEST( read.size() == 4 );
    BOOST_TEST( read[0] ==  2 );
    BOOST_TEST( read[1] ==  3 );
    BOOST_TEST( read[2] == -5 );
    BOOST_TEST( read[3] ==  1 );
}

BOOST_AUTO_TEST_CASE( readAccessFunctionGivesAccessForCorrectContext )
{
    const auto context = Context::getContextPtr();
    HArray<int> array ( { 2, 3, -5, 1 }, context );

    const auto expectedDataPointer = [&] ()
    {
        return ReadAccess<int>( array, context ).get();
    }
    ();

    auto read = readAccess( array, context );
    BOOST_TEST( read.get() == expectedDataPointer );
    BOOST_TEST( read.getMemory().getContextPtr() == context );
}

BOOST_AUTO_TEST_CASE( writeAccessMoveConstructorTest )
{
    HArray<int> array { 2, 3, -1 };

    WriteAccess<int> write1( array );
    const auto data = write1.get();
    const auto size = write1.size();

    WriteAccess<int> write2( std::move( write1 ) );

    BOOST_TEST( write2.get() == data );
    BOOST_TEST( write2.size() == size );

    // Trying to access the data should trigger an assertion
    // (note: a more specific and documented exception would be preferable here)
    BOOST_CHECK_THROW( write1.get(), scai::common::AssertException );
}

BOOST_AUTO_TEST_CASE( writeAccessFunctionTest )
{
    const auto context = Context::getContextPtr();
    HArray<int> array ( { 2, 3, -5, 1 }, context );

    auto write = writeAccess( array );

    BOOST_TEST( write.size() == 4 );
    BOOST_TEST( write[0] ==  2 );
    BOOST_TEST( write[1] ==  3 );
    BOOST_TEST( write[2] == -5 );
    BOOST_TEST( write[3] ==  1 );
}

BOOST_AUTO_TEST_CASE( writeAccessFunctionGivesAccessForCorrectContext )
{
    const auto context = Context::getContextPtr();
    HArray<int> array ( { 2, 3, -5, 1 }, context );

    const auto expectedDataPointer = [&] ()
    {
        return WriteAccess<int>( array, context ).get();
    }
    ();

    auto write = writeAccess( array, context );
    BOOST_TEST( write.get() == expectedDataPointer );
    BOOST_TEST( write.getMemory().getContextPtr() == context );
}

BOOST_AUTO_TEST_CASE( writeOnlyAccessMoveConstructorTest )
{
    HArray<int> array { 2, 3, -1 };

    WriteOnlyAccess<int> write1( array, 4 );
    const auto data = write1.get();
    const auto size = write1.size();

    WriteOnlyAccess<int> write2( std::move( write1 ) );

    BOOST_TEST( write2.get() == data );
    BOOST_TEST( write2.size() == size );

    // Trying to access the data should trigger an assertion
    // (note: a more specific and documented exception would be preferable here)
    BOOST_CHECK_THROW( write1.get(), scai::common::AssertException );
}

BOOST_AUTO_TEST_CASE( writeOnlyAccessFunctionTest )
{
    const auto context = Context::getContextPtr();
    HArray<int> array ( { 2, 3, -5, 1 }, context );

    {
        auto write = writeOnlyAccess( array, 3 );
        BOOST_TEST( write.size() == 3 );
        write[0] =  4;
        write[1] =  6;
        write[2] = -3;
    }

    {
        HostReadAccess<int> read( array );
        BOOST_TEST( read.size() == 3 );
        BOOST_TEST( read[0] ==  4 );
        BOOST_TEST( read[1] ==  6 );
        BOOST_TEST( read[2] == -3 );
    }
}

BOOST_AUTO_TEST_CASE( writeOnlyAccessFunctionGivesAccessForCorrectContext )
{
    const auto context = Context::getContextPtr();
    HArray<int> array ( { 2, 3, -5, 1 }, context );

    // For the writeOnlyAccess, we can not check that the pointer is the same
    // like we did for ReadAccess and WriteAccess, because WriteOnly is permitted
    // to allocate different memory, so we only verify that the context is correct.
    BOOST_TEST( writeOnlyAccess( array, context, 3 ).getMemory().getContextPtr() == context );
}

BOOST_AUTO_TEST_CASE( hostReadAccessMoveConstructorTest )
{
    HArray<int> array { 2, 3, -1 };

    HostReadAccess<int> read1( array );
    const auto data = read1.get();
    const auto size = read1.size();

    HostReadAccess<int> read2( std::move( read1 ) );

    BOOST_TEST( read2.get() == data );
    BOOST_TEST( read2.size() == size );

    // Trying to access the data should trigger an assertion
    // (note: a more specific and documented exception would be preferable here)
    BOOST_CHECK_THROW( read1.get(), scai::common::AssertException );
}

BOOST_AUTO_TEST_CASE( hostReadAccessFunctionTest )
{
    const auto context = Context::getContextPtr();
    HArray<int> array ( { 2, 3, -5, 1 }, context );

    const auto read = hostReadAccess( array );

    BOOST_TEST( read.size() == 4 );
    BOOST_TEST( read[0] ==  2 );
    BOOST_TEST( read[1] ==  3 );
    BOOST_TEST( read[2] == -5 );
    BOOST_TEST( read[3] ==  1 );
}

BOOST_AUTO_TEST_CASE( hostWriteAccessMoveConstructorTest )
{
    HArray<int> array { 2, 3, -1 };

    HostWriteAccess<int> write1( array );
    const auto data = write1.get();
    const auto size = write1.size();

    HostWriteAccess<int> write2( std::move( write1 ) );

    BOOST_TEST( write2.get() == data );
    BOOST_TEST( write2.size() == size );

    // Trying to access the data should trigger an assertion
    // (note: a more specific and documented exception would be preferable here)
    BOOST_CHECK_THROW( write1.get(), scai::common::AssertException );
}

BOOST_AUTO_TEST_CASE( hostWriteAccessFunctionTest )
{
    const auto context = Context::getContextPtr();
    HArray<int> array ( { 2, 3, -5, 1 }, context );

    auto write = hostWriteAccess( array );

    BOOST_TEST( write.size() == 4 );
    BOOST_TEST( write[0] ==  2 );
    BOOST_TEST( write[1] ==  3 );
    BOOST_TEST( write[2] == -5 );
    BOOST_TEST( write[3] ==  1 );
}

BOOST_AUTO_TEST_CASE( hostWriteOnlyAccessMoveConstructorTest )
{
    HArray<int> array { 2, 3, -1 };

    HostWriteOnlyAccess<int> write1( array, 4 );
    const auto data = write1.get();
    const auto size = write1.size();

    HostWriteOnlyAccess<int> write2( std::move( write1 ) );

    BOOST_TEST( write2.get() == data );
    BOOST_TEST( write2.size() == size );

    // Trying to access the data should trigger an assertion
    // (note: a more specific and documented exception would be preferable here)
    BOOST_CHECK_THROW( write1.get(), scai::common::AssertException );
}

BOOST_AUTO_TEST_CASE( hostWriteOnlyAccessFunctionTest )
{
    const auto context = Context::getContextPtr();
    HArray<int> array ( { 2, 3, -5, 1 }, context );

    {
        auto write = hostWriteOnlyAccess( array, 3 );
        BOOST_TEST( write.size() == 3 );
        write[0] =  4;
        write[1] =  6;
        write[2] = -3;
    }

    {
        HostReadAccess<int> read( array );
        BOOST_TEST( read.size() == 3 );
        BOOST_TEST( read[0] ==  4 );
        BOOST_TEST( read[1] ==  6 );
        BOOST_TEST( read[2] == -3 );
    }

}

BOOST_AUTO_TEST_CASE( hostReadAccessSTLCompatibilityTest )
{
    const auto array = HArray<int> { 5, 3, 1, 2 };
    const std::vector<int> expectedElements { 5, 3, 1, 2 };

    // Range-based for
    {
        std::vector<int> elements;

        for ( auto element : hostReadAccess( array ) )
        {
            elements.push_back( element );
        }

        BOOST_TEST( elements == expectedElements, boost::test_tools::per_element() );
    }

    // copy, begin/end
    {
        std::vector<int> elements;
        const auto access = hostReadAccess( array );
        std::copy( access.begin(), access.end(), std::back_inserter( elements ) );

        BOOST_TEST( elements == expectedElements, boost::test_tools::per_element() );
    }

    // copy, cbegin/cend
    {
        std::vector<int> elements;
        const auto access = hostReadAccess( array );
        std::copy( access.cbegin(), access.cend(), std::back_inserter( elements ) );

        BOOST_TEST( elements == expectedElements, boost::test_tools::per_element() );
    }

    // is_heap, which requires a random access iterator (unlike copy which only requires an input iterator)
    {
        const auto access = hostReadAccess( array );

        BOOST_TEST ( std::is_heap( access.begin(), access.end() ) );
    }
}

BOOST_AUTO_TEST_CASE( hostWriteAccessSTLCompatibilityTest )
{
    const auto array = HArray<int> { 5, 3, 1, 2 };
    const std::vector<int> expectedElements { 5, 3, 1, 2 };

    // Range-based for, read-only
    {
        auto mutableArray = array;
        std::vector<int> elements;

        for ( auto element : hostWriteAccess( mutableArray ) )
        {
            elements.push_back( element );
        }

        BOOST_TEST( elements == expectedElements, boost::test_tools::per_element() );
    }

    // Range-based for, write
    {
        auto mutableArray = array;

        for ( auto& element : hostWriteAccess( mutableArray ) )
        {
            element = 2;
        }

        const std::vector<int> expectedMutated { 2, 2, 2, 2 };
        BOOST_TEST( hostWriteAccess( mutableArray ) == expectedMutated, boost::test_tools::per_element() );
    }

    // copy, begin/end
    {
        auto mutableArray = array;
        std::vector<int> elements;
        const auto access = hostWriteAccess( mutableArray );
        std::copy( access.begin(), access.end(), std::back_inserter( elements ) );

        BOOST_TEST( elements == expectedElements, boost::test_tools::per_element() );
    }

    // copy, cbegin/cend
    {
        auto mutableArray = array;
        std::vector<int> elements;
        const auto access = hostWriteAccess( mutableArray );
        std::copy( access.cbegin(), access.cend(), std::back_inserter( elements ) );

        BOOST_TEST( elements == expectedElements, boost::test_tools::per_element() );
    }

    // is_heap, which requires a random access iterator (unlike copy which only requires an input iterator)
    {
        auto mutableArray = array;
        const auto access = hostWriteAccess( mutableArray );

        BOOST_TEST ( std::is_heap( access.begin(), access.end() ) );
    }

    // copy into HostWriteAccess, begin/end
    {
        auto mutableArray = array;
        auto access = hostWriteAccess( mutableArray );
        std::copy( expectedElements.begin(), expectedElements.end(), access.begin() );

        BOOST_TEST( access == expectedElements, boost::test_tools::per_element() );
    }

    // sort, requires random access
    {
        auto mutableArray = array;
        auto access = hostWriteAccess( mutableArray );
        std::sort( access.begin(), access.end() );

        const std::vector<int> expectedSorted { 1, 2, 3, 5 };
        BOOST_TEST ( access == expectedSorted, boost::test_tools::per_element() );
    }
}

BOOST_AUTO_TEST_CASE( hostWriteOnlyAccessSTLCompatibilityTest )
{
    const auto array = HArray<int> { 5, 3, 1, 2 };
    const std::vector<int> expectedElements { 5, 3, 1, 2 };

    // Range-based for, write
    {
        auto mutableArray = array;

        for ( auto& element : hostWriteOnlyAccess( mutableArray, 3 ) )
        {
            element = 2;
        }

        const std::vector<int> expectedMutated { 2, 2, 2 };
        BOOST_TEST( hostReadAccess( mutableArray ) == expectedMutated, boost::test_tools::per_element() );
    }

    // copy into HostWriteOnlyAccess, begin/end
    {
        auto mutableArray = array;
        auto access = hostWriteOnlyAccess( mutableArray, 4 );
        std::copy( expectedElements.begin(), expectedElements.end(), access.begin() );

        BOOST_TEST( access == expectedElements, boost::test_tools::per_element() );
    }

    // copy into followed by sort, requires random access
    {
        auto mutableArray = array;
        auto access = hostWriteOnlyAccess( mutableArray, 4 );
        std::copy( expectedElements.begin(), expectedElements.end(), access.begin() );
        std::sort( access.begin(), access.end() );

        const std::vector<int> expectedSorted { 1, 2, 3, 5 };
        BOOST_TEST ( access == expectedSorted, boost::test_tools::per_element() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
