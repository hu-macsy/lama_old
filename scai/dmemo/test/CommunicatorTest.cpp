/**
 * @file CommunicatorTest.cpp
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
 * @brief Test methods for the communicator
 * @author: Thomas Brandes
 * @date 09.05.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/logging.hpp>

#include <scai/dmemo.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>

#include <scai/tasking.hpp>
#include <scai/hmemo.hpp>

#include <scai/common/unique_ptr.hpp>
#include <scai/common/exception/Exception.hpp>

using namespace scai;
using namespace hmemo;
using namespace dmemo;
using namespace tasking;

using common::Exception;
using common::unique_ptr;
using common::scoped_array;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CommunicatorTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.CommunicatorTest" )

/* --------------------------------------------------------------------- */

// ToDo: take ARITHMETIC_HOST_TYPE for the tests

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( basicTest )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    SCAI_LOG_INFO( logger, "basicTest for comm = " << *comm )

    BOOST_CHECK( comm->getRank() >= 0 );
    BOOST_CHECK( comm->getRank() < comm->getSize() );

    std::stringstream out;
 
    out << *comm;

    BOOST_CHECK( out.str().length() > 0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( computeOwnersTest )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    IndexType rank = comm->getRank();
    IndexType size = comm->getSize();

    IndexType n = 17;

    std::vector<IndexType> localIndexes;
    std::vector<IndexType> nonLocalIndexes;

    for ( IndexType i = 0; i < n; ++i )
    {
        localIndexes.push_back( rank * n + i );
    }

    GeneralDistribution dist( n * size, localIndexes, comm );

    for ( PartitionId p = 0; p < size; ++p )
    {
        if ( p == rank )
        {
            for ( IndexType i = 0; i < n; ++i )
            {
                BOOST_CHECK( dist.isLocal( p * n + i ) );
            }
        }
        else
        {
            for ( IndexType i = 0; i < n; ++i )
            {
                nonLocalIndexes.push_back( p * n + i );
            }
        }
    }

    std::vector < PartitionId > owners;
    comm->computeOwners( nonLocalIndexes, dist, owners );
    std::vector<PartitionId>::size_type currentIndex = 0;

    for ( PartitionId p = 0; p < size; ++p )
    {
        if ( p != rank )
        {
            for ( IndexType i = 0; i < n; ++i )
            {
                BOOST_CHECK_EQUAL( p, owners[currentIndex++] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( allocatePlanTest )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    IndexType rank = comm->getRank();
    IndexType size = comm->getSize();

    std::vector<IndexType> reqQuantities( size );

    for ( PartitionId p = 0; p < size; ++p )
    {
        if ( p != rank )
        {
            reqQuantities[p] = ( 2 * p + rank ) % 3;
        }
        else
        {
            reqQuantities[p] = 0;
        }
    }

    CommunicationPlan requiredPlan( reqQuantities.data(), reqQuantities.size() );

    // verify that requiredPlan is correctly set up

    IndexType offsetCheck = 0;

    for ( PartitionId p = 0; p < requiredPlan.size(); ++p )
    {
        IndexType n = requiredPlan[p].quantity;
        PartitionId partitionId = requiredPlan[p].partitionId;
        IndexType nExpected = ( 2 * partitionId + rank ) % 3;
        BOOST_CHECK_EQUAL( n, nExpected );
        BOOST_CHECK_EQUAL( requiredPlan[p].offset, offsetCheck );
        offsetCheck += n;
    }

    BOOST_CHECK_EQUAL( offsetCheck, requiredPlan.totalQuantity() );
    CommunicationPlan providesPlan;
    providesPlan.allocateTranspose( requiredPlan, *comm );
    offsetCheck = 0;

    for ( PartitionId p = 0; p < providesPlan.size(); ++p )
    {
        IndexType n = providesPlan[p].quantity;
        PartitionId partitionId = providesPlan[p].partitionId;
        IndexType nExpected = ( partitionId + 2 * rank ) % 3;
        BOOST_CHECK_EQUAL( n, nExpected );
        BOOST_CHECK_EQUAL( providesPlan[p].offset, offsetCheck );
        offsetCheck += n;
    }

    BOOST_CHECK_EQUAL( offsetCheck, providesPlan.totalQuantity() );

}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( bcastStringTest )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    std::string val = "Dummy";

    IndexType root = 0;

    if ( comm->getRank() == root )
    {
        val = "Hello";
    }

    SCAI_LOG_INFO( logger, *comm << ": val = " << val );

    comm->bcast( val, root );

    SCAI_LOG_INFO( logger, *comm << ": val = " << val );

    BOOST_CHECK_EQUAL( "Hello", val );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildHaloTest )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    IndexType rank = comm->getRank();
    IndexType size = comm->getSize();

    IndexType vectorSize = size;
    BlockDistribution distribution( vectorSize, comm );
    std::vector<IndexType> requiredIndexes;
    const PartitionId leftNeighbor = comm->getNeighbor( -1 );
    const PartitionId rightNeighbor = comm->getNeighbor( 1 );
// Each processor requires values from left and right neighbor

    if ( !distribution.isLocal( leftNeighbor ) )
    {
        requiredIndexes.push_back( leftNeighbor );
    }

    if ( rightNeighbor != leftNeighbor && !distribution.isLocal( rightNeighbor ) )
    {
        requiredIndexes.push_back( rightNeighbor );
    }

    const IndexType noReqIndexes = static_cast<IndexType>( requiredIndexes.size() );

    Halo halo;

    HaloBuilder::build( distribution, requiredIndexes, halo );

    const Halo& haloRef = halo;

    const CommunicationPlan& requiredPlan = haloRef.getRequiredPlan();

    const CommunicationPlan& providesPlan = haloRef.getProvidesPlan();

    // check for a correct provide plan

    IndexType offsetCheck = 0;

    for ( PartitionId p = 0; p < requiredPlan.size(); ++p )
    {
        IndexType n = requiredPlan[p].quantity;
        BOOST_CHECK_EQUAL( ( IndexType ) 1, n );
        PartitionId neighbor = requiredPlan[p].partitionId;
        BOOST_CHECK( neighbor == leftNeighbor || neighbor == rightNeighbor );
        BOOST_CHECK_EQUAL( requiredPlan[p].offset, offsetCheck );
        offsetCheck += n;
    }

    BOOST_CHECK_EQUAL( noReqIndexes, requiredPlan.totalQuantity() );

    offsetCheck = 0;
    PartitionId nProvides = providesPlan.size();

    for ( PartitionId p = 0; p < nProvides; ++p )
    {
        IndexType n = providesPlan[p].quantity;
        BOOST_CHECK_EQUAL( n, static_cast<IndexType>( 1 ) );
        PartitionId neighbor = providesPlan[p].partitionId;
        BOOST_CHECK( neighbor == leftNeighbor || neighbor == rightNeighbor );
        BOOST_CHECK_EQUAL( providesPlan[p].offset, offsetCheck );
        offsetCheck += n;
    }

    BOOST_CHECK_EQUAL( noReqIndexes, providesPlan.totalQuantity() );

    const ReadAccess<IndexType> providesIndexes( haloRef.getProvidesIndexes() );

    for ( PartitionId p = 0; p < providesPlan.size(); ++p )
    {
        const IndexType* indexes = providesIndexes + providesPlan[p].offset;
        IndexType expectedLocalIndex = rank;
        BOOST_CHECK_EQUAL( expectedLocalIndex, distribution.local2global( indexes[0] ) );
    }

    BOOST_CHECK_EQUAL( noReqIndexes, halo.getHaloSize() );

    IndexType nIndexes = static_cast<IndexType>( requiredIndexes.size() );

    for ( IndexType i = 0; i < nIndexes; ++i )
    {
        const IndexType haloIndex = halo.global2halo( requiredIndexes[i] );
        BOOST_CHECK( 0 <= haloIndex && haloIndex < halo.getHaloSize() );
    }

}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( updateHaloTest, ValueType, test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    IndexType rank = comm->getRank();
    IndexType size = comm->getSize();

    SCAI_LOG_INFO( logger, "updateHaloTest<" << common::getScalarType<ValueType>() << ">" );
    const IndexType factor = 4;
    const IndexType vectorSize = factor * size;
    BlockDistribution distribution( vectorSize, comm );
    std::vector<IndexType> requiredIndexes;

    for ( IndexType i = 0; i < factor; ++i )
    {
        const IndexType requiredIndex = ( ( rank + 1 ) * factor + i ) % vectorSize;

        if ( distribution.isLocal( requiredIndex ) )
        {
            continue;
        }

        requiredIndexes.push_back( requiredIndex );
    }

    SCAI_LOG_INFO( logger, "build the Halo" );
    Halo halo;
    HaloBuilder::build( distribution, requiredIndexes, halo );
    SCAI_LOG_INFO( logger, "halo is now available: " << halo );
    HArray<ValueType> localData;
    {
        WriteOnlyAccess<ValueType> localDataAccess( localData, distribution.getLocalSize() );

        for ( IndexType i = 0; i < localData.size(); ++i )
        {
            localDataAccess[i] = static_cast<ValueType>( distribution.local2global( i ) );
        }
    }
    SCAI_LOG_INFO( logger, "update halo data by communicator" );
    HArray<ValueType> haloData;
    comm->updateHalo( haloData, localData, halo );
    BOOST_CHECK_EQUAL( static_cast<IndexType>( requiredIndexes.size() ), haloData.size() );
    {
        ReadAccess<ValueType> haloDataAccess( haloData );

        for ( IndexType i = 0; i < static_cast<IndexType>( requiredIndexes.size() ); ++i )
        {
            ValueType expectedValue = static_cast<ValueType>( requiredIndexes[i] );
            BOOST_CHECK_EQUAL( expectedValue, haloDataAccess[i] );
        }
    }
    requiredIndexes.clear();

    for ( IndexType i = 0; i < vectorSize; ++i )
    {
        if ( distribution.isLocal( i ) || ( i + rank ) % 2 == 0 )
        {
            continue;
        }

        requiredIndexes.push_back( i );
    }

    HaloBuilder::build( distribution, requiredIndexes, halo );
    comm->updateHalo( haloData, localData, halo );
    BOOST_CHECK_EQUAL( static_cast<IndexType>( requiredIndexes.size() ), haloData.size() );
    {
        ReadAccess<ValueType> haloDataAccess( haloData );

        for ( IndexType i = 0; i < static_cast<IndexType>( requiredIndexes.size() ); ++i )
        {
            ValueType expectedValue = static_cast<ValueType>( requiredIndexes[i] );
            BOOST_CHECK_EQUAL( expectedValue, haloDataAccess[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( shiftTest, ValueType, test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    IndexType rank = comm->getRank();
    IndexType size = comm->getSize();

    // Idea of this shift Test:
    // - allocate on each processor an array with one element for each processor
    // - shift this array around all processors and each processor writes one value at its rank
    // - verify that each processor has written the right value
    if ( size > 1 )
    {
        const IndexType vectorSize = size;
        HArray<ValueType> sendBuffer( vectorSize, static_cast<ValueType>( rank ) );
        HArray<ValueType> recvBuffer;

        for ( PartitionId rounds = 0; rounds < size; ++rounds )
        {
            comm->shiftArray( recvBuffer, sendBuffer, 1 );
            {
                WriteAccess<ValueType> recvBufferAccess( recvBuffer );
                recvBufferAccess[rank] = static_cast<ValueType>( rank );
            }
            sendBuffer.swap( recvBuffer );
        }

        {
            ReadAccess<ValueType> recvBufferAccess( recvBuffer );

            for ( IndexType i = 0; i < size; ++i )
            {
                ValueType value = static_cast<ValueType>( i );
                BOOST_CHECK_EQUAL( value, recvBufferAccess[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( shiftAsyncTest, ValueType, test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    IndexType rank = comm->getRank();

    HArray<ValueType> sendBuffer( 2, static_cast<ValueType>( rank ) );
    HArray<ValueType> recvBuffer;
    {
        WriteAccess<ValueType> sbuffer( sendBuffer );
        sbuffer[ 0 ] = static_cast<ValueType>( rank );
        sbuffer[ 1 ] = static_cast<ValueType>( comm->getNeighbor( 1 ) );
    }
    // if we do not keep the token, it will be synchronized immeadiately
    unique_ptr<SyncToken> token1( comm->shiftAsync( recvBuffer, sendBuffer, 1 ) );
    SCAI_LOG_INFO( logger, "token for shiftAsync before wait : " << *token1 );
    token1->wait();
    SCAI_LOG_INFO( logger, "token for shiftAsync after wait : " << *token1 );
    SCAI_LOG_INFO( logger, "async shift without token, should have been synchronized here" );
    BOOST_CHECK_EQUAL( 2, recvBuffer.size() );
    {
        ReadAccess<ValueType> rbuffer( recvBuffer );
    }
    // dir = 0: should be like assignment
    delete comm->shiftAsync( recvBuffer, sendBuffer, 0 );
    SCAI_LOG_INFO( logger, "async shift dir = 0, self assing" );
    BOOST_CHECK_EQUAL( 2, recvBuffer.size() );
    {
        ReadAccess<ValueType> sbuffer( sendBuffer );
        ReadAccess<ValueType> rbuffer( recvBuffer );
        BOOST_CHECK_EQUAL( sbuffer[0], rbuffer[ 0 ] );
        BOOST_CHECK_EQUAL( sbuffer[1], rbuffer[ 1 ] );
    }
    // Buffers must be different
    SCAI_LOG_INFO( logger, "async shift : using same send and recv buffer should fail" );
    BOOST_CHECK_THROW( delete comm->shiftAsync( sendBuffer, sendBuffer, 1 ), Exception );
    // We also verify that the exception is thrown before changing the send Buffer
    BOOST_CHECK_EQUAL( 2, sendBuffer.size() );
    SCAI_LOG_INFO( logger, "async shift : try to access send / receive buffer before synchronization" );
    unique_ptr<SyncToken> token( comm->shiftAsync( recvBuffer, sendBuffer, 1 ) );
    // read access on send buffer should be possible
    {
        ReadAccess<ValueType> sbuffer( sendBuffer );
        ValueType value0 = static_cast<ValueType>( rank );
        BOOST_CHECK_EQUAL( value0, sbuffer[ 0 ] );
    }
    SCAI_LOG_INFO( logger, *token << ": test for correct locks" );

    // Note: multiple accesses on HArray are possible, the following exceptions  are no more thrown

    if ( !token->isSynchronized() )
    {
        // write access on send buffer should be locked
        // SCAI_CHECK_THROW( WriteAccess<ValueType> sbuffer( sendBuffer ), Exception );
    }

    if ( !token->isSynchronized() )
    {
        // read access on recv buffer should be locked if communication is not finished yet
        // SCAI_CHECK_THROW( ReadAccess<ValueType> recvBufferAccess( recvBuffer ), Exception );
    }

    token->wait();
    token->wait(); // second wait should not harm
    {
        WriteAccess<ValueType> sendBufferWrite( sendBuffer );
        sendBufferWrite.resize( 0 );
    }
    BOOST_CHECK_EQUAL( 2, recvBuffer.size() );
    {
        ReadAccess<ValueType> rbuffer( recvBuffer );
        ValueType value0 = static_cast<ValueType>( comm->getNeighbor( -1 ) );
        ValueType value1 = static_cast<ValueType>( rank );
        BOOST_CHECK_EQUAL( value0, rbuffer[ 0 ] );
        BOOST_CHECK_EQUAL( value1, rbuffer[ 1 ] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( bcastTest, ValueType, test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    IndexType rank = comm->getRank();
    IndexType size = comm->getSize();

    SCAI_LOG_INFO( logger, "bcastTest<" << common::getScalarType<ValueType>() << ">" )
    IndexType N = 5;
    ValueType dummyVal = 13;
    scoped_array<ValueType> vector( new ValueType[N + 1] );
    vector[N] = dummyVal;

    for ( PartitionId p = 0; p < size; p++ )
    {
        if ( p == rank )
        {
            for ( IndexType i = 0; i < N; i++ )
            {
                vector[i] = static_cast<ValueType>( i + p );
            }
        }

        // processor p bcast the new vector
        comm->bcast( vector.get(), N, p );

        // all processors must now have the correct values

        for ( IndexType i = 0; i < N; i++ )
        {
            ValueType value = static_cast<ValueType>( i + p );
            BOOST_CHECK_EQUAL( value, vector[i] );
        }

        // make sure that next value has not been overwritten
        BOOST_CHECK_EQUAL( dummyVal, vector[N] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scatterTest, ValueType, test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    IndexType rank = comm->getRank();
    IndexType size = comm->getSize();

    const PartitionId root = 0;
    IndexType n = 2;
    IndexType allN = 0; // only root will have full size

    if ( rank == root )
    {
        allN = size * n;
    }

    ValueType dummyVal = 13;
    scoped_array<ValueType> myvals( new ValueType[n + 1] );
    scoped_array<ValueType> allvals( new ValueType[allN] );
    myvals[0] = 0;
    myvals[1] = 1;
    myvals[2] = size + dummyVal;

    if ( rank == root )
    {
        // fill the send data
        for ( int i = 0; i < size; i++ )
        {
            allvals[2 * i] = static_cast<ValueType>( i );
            allvals[2 * i + 1] = dummyVal + i;
        }
    }

    comm->scatter( myvals.get(), n, root, allvals.get() );
    // also verify that myvals has not been overwritten at end
    BOOST_CHECK_EQUAL( myvals[0], static_cast<ValueType>( rank ) );
    BOOST_CHECK_EQUAL( myvals[1], dummyVal + rank );
    BOOST_CHECK_EQUAL( myvals[2], size + dummyVal );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scatterVTest, ValueType, test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    IndexType rank = comm->getRank();
    IndexType size = comm->getSize();

    const PartitionId root = 0;
    IndexType n = rank; // number of elements I receive
    IndexType allN = 0;

    if ( rank == root )
    {
        allN = size * ( size - 1 ) / 2;
    }

    ValueType dummyVal = 13;
    scoped_array<ValueType> myvals( new ValueType[n + 1] );
    scoped_array<ValueType> allvals( new ValueType[allN] );
    scoped_array<int> sizes( new int[size] );

    for ( int i = 0; i <= n; i++ )
    {
        myvals[i] = dummyVal;
    }

    if ( rank == root )
    {
        // fill the send data
        int offset = 0;

        for ( int p = 0; p < size; p++ )
        {
            for ( int i = 0; i < p; i++ )
            {
                allvals[offset] = dummyVal + i + 1;
                ++offset;
            }

            sizes[p] = p;
        }
    }

    comm->scatterV( myvals.get(), n, root, allvals.get(), sizes.get() );

    // also verify that myvals has not been overwritten at end
    for ( int i = 0; i < n; i++ )
    {
        BOOST_CHECK_EQUAL( myvals[i], dummyVal + i + 1 );
    }

    BOOST_CHECK_EQUAL( myvals[n], dummyVal );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gatherTest, ValueType, test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    IndexType rank = comm->getRank();
    IndexType size = comm->getSize();

    const PartitionId root = 0;
    IndexType allN = 0; // only root will have full size

    if ( rank == root )
    {
        allN = size;
    }

    ValueType dummyVal = -1;
    scoped_array<ValueType> myvals( new ValueType[2] );
    scoped_array<ValueType> allvals( new ValueType[allN] );

    if ( rank == root )
    {
        for ( IndexType i = 0; i < size; ++i )
        {
            allvals[i] = static_cast<ValueType>( i );
        }
    }

    myvals[0] = static_cast<ValueType>( rank );
    myvals[1] = dummyVal;
    comm->gather( allvals.get(), 1, root, myvals.get() );

    // also verify that myvals has not been overwritten at end

    if ( rank == root )
    {
        for ( IndexType i = 0; i < size; ++i )
        {
            BOOST_CHECK_EQUAL( allvals[i], static_cast<ValueType>( i ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gatherVTest, ValueType, test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    IndexType rank = comm->getRank();
    IndexType size = comm->getSize();

    const PartitionId root = 0;
    IndexType n = rank; // number of elements I send
    IndexType allN = 0;

    if ( rank == root )
    {
        allN = size * ( size - 1 ) / 2;
    }

    ValueType dummyVal = 13;
    scoped_array<ValueType> myvals( new ValueType[n + 1] );
    scoped_array<ValueType> allvals( new ValueType[allN] );
    scoped_array<int> sizes( new int[size] );

    for ( int i = 0; i < n; i++ )
    {
        myvals[i] = static_cast<ValueType>( rank );
    }

    myvals[n] = dummyVal;

    if ( rank == root )
    {
        // fill the recv data
        int offset = 0;

        for ( int p = 0; p < size; p++ )
        {
            for ( int i = 0; i < p; i++ )
            {
                allvals[offset] = -1;
                ++offset;
            }

            sizes[p] = p;
        }
    }

    comm->gatherV( allvals.get(), n, root, myvals.get(), sizes.get() );

    if ( rank == root )
    {
        int offset = 0;

        for ( int p = 0; p < size; ++p )
        {
            for ( int i = 0; i < p; ++i )
            {
                BOOST_CHECK_EQUAL( allvals[offset], static_cast<ValueType>( p ) );
                ++offset;
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( swapTest, ValueType, test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicator();

    BOOST_REQUIRE( comm );

    IndexType rank = comm->getRank();
    IndexType size = comm->getSize();

    int n = 10;
    scoped_array<ValueType> vector( new ValueType[n] );

    // initialize vector individually for each processor

    for ( IndexType i = 0; i < n; i++ )
    {
        vector[i] = static_cast<ValueType>( 2 * rank + 3 * i );
    }

    PartitionId partner = size - 1 - rank;
    comm->swap( vector.get(), n, partner );

    for ( IndexType i = 0; i < n; i++ )
    {
        ValueType value = static_cast<ValueType>( 2 * partner + 3 * i );
        BOOST_CHECK_EQUAL( value, vector[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */