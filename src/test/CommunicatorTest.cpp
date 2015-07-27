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
 * @brief Contains the implementation of the class CommunicatorTest
 * @author: Alexander Büchel, Thomas Brandes
 * @date 09.05.2012
 * @since 1.0.0
 **/

#include <test/distributed/CommunicatorTest.hpp>

#include <test/TestMacros.hpp>

#include <lama/LAMATypes.hpp>

#include <lama/Communicator.hpp>

#include <lama/distribution/HaloBuilder.hpp>
#include <lama/distribution/GeneralDistribution.hpp>
#include <lama/distribution/BlockDistribution.hpp>

#include <common/unique_ptr.hpp>

using namespace lama;
using namespace memory;
using namespace tasking;
using common::Exception;
using common::unique_ptr;
using common::scoped_array;

LAMA_LOG_DEF_LOGGER( logger, "Test.CommunicatorTest" )

/* --------------------------------------------------------------------- */

CommunicatorTest::CommunicatorTest( const char* communicatorType )
    :
    mCommunicatorType( communicatorType )
{
    comm = Communicator::get( mCommunicatorType.c_str() );
    rank = comm->getRank();
    size = comm->getSize();
}

CommunicatorTest::~CommunicatorTest()
{
    comm = CommunicatorPtr();
}

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE( CommunicatorTest, CommunicatorCtrTest )
// get a communicator and give it free (shared pointer)
comm = Communicator::get( mCommunicatorType.c_str() );
comm = CommunicatorPtr();
// get again the MPI communicator (might be 2nd call of MPI_Init)
comm = Communicator::get( mCommunicatorType.c_str() );
LAMA_COMMON_TEST_CASE_END()

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE( CommunicatorTest, computeOwnersTest )
IndexType n = 17;
std::vector<IndexType> localIndexes, nonLocalIndexes;

for ( IndexType i = 0; i < n; ++i )
{
    localIndexes.push_back( rank * n + i );
}

GeneralDistribution dist( n* size, localIndexes, comm );

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

LAMA_COMMON_TEST_CASE_END()

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE( CommunicatorTest, allocatePlanTest )

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
LAMA_COMMON_TEST_CASE_END()

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE( CommunicatorTest, bcastStringTest )
std::string val = "Dummy";

if ( comm->getRank() == 0 )
{
    val = "Hello";
}

LAMA_LOG_INFO( logger, *comm << ": val = " << val );

comm->bcast( val, 0 );

LAMA_LOG_INFO( logger, *comm << ": val = " << val );

BOOST_CHECK_EQUAL( "Hello", val );
LAMA_COMMON_TEST_CASE_END()

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE( CommunicatorTest, buildHaloTest )
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

LAMA_COMMON_TEST_CASE_END()

/* --------------------------------------------------------------------- */

template<typename ValueType>
void CommunicatorTest::updateHaloTest()
{
    LAMA_LOG_INFO( logger, "updateHaloTest<" << common::getScalarType<ValueType>() << ">" );
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

    LAMA_LOG_INFO( logger, "build the Halo" );
    Halo halo;
    HaloBuilder::build( distribution, requiredIndexes, halo );
    LAMA_LOG_INFO( logger, "halo is now available: " << halo );
    LAMAArray<ValueType> localData;
    {
        WriteOnlyAccess<ValueType> localDataAccess( localData, distribution.getLocalSize() );

        for ( IndexType i = 0; i < localData.size(); ++i )
        {
            localDataAccess[i] = static_cast<ValueType>( distribution.local2global( i ) );
        }
    }
    LAMA_LOG_INFO( logger, "update halo data by communicator" );
    LAMAArray<ValueType> haloData;
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

template<typename ValueType>
void CommunicatorTest::shiftTest()
{
    // Idea of this shift Test:
    // - allocate on each processor an array with one element for each processor
    // - shift this array around all processors and each processor writes one value at its rank
    // - verify that each processor has written the right value
    if ( size > 1 )
    {
        const IndexType vectorSize = size;
        LAMAArray<ValueType> sendBuffer( vectorSize, static_cast<ValueType>( rank ) );
        LAMAArray<ValueType> recvBuffer;

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

LAMA_COMMON_TEST_CASE_TM( CommunicatorTest, ValueType, shiftASyncTest )
{
    LAMAArray<ValueType> sendBuffer( 2, static_cast<ValueType>( rank ) );
    LAMAArray<ValueType> recvBuffer;
    {
        WriteAccess<ValueType> sbuffer( sendBuffer );
        sbuffer[ 0 ] = static_cast<ValueType>( rank );
        sbuffer[ 1 ] = static_cast<ValueType>( comm->getNeighbor( 1 ) );
    }
    // if we do not keep the token, it will be synchronized immeadiately
    unique_ptr<SyncToken> token1( comm->shiftAsync( recvBuffer, sendBuffer, 1 ) );
    LAMA_LOG_INFO( logger, "token for shiftAsync before wait : " << *token1 );
    token1->wait();
    LAMA_LOG_INFO( logger, "token for shiftAsync after wait : " << *token1 );
    LAMA_LOG_INFO( logger, "async shift without token, should have been synchronized here" );
    BOOST_CHECK_EQUAL( 2, recvBuffer.size() );
    {
        ReadAccess<ValueType> rbuffer( recvBuffer );
    }
    // dir = 0: should be like assignment
    delete comm->shiftAsync( recvBuffer, sendBuffer, 0 );
    LAMA_LOG_INFO( logger, "async shift dir = 0, self assing" );
    BOOST_CHECK_EQUAL( 2, recvBuffer.size() );
    {
        ReadAccess<ValueType> sbuffer( sendBuffer );
        ReadAccess<ValueType> rbuffer( recvBuffer );
        BOOST_CHECK_EQUAL( sbuffer[0], rbuffer[ 0 ] );
        BOOST_CHECK_EQUAL( sbuffer[1], rbuffer[ 1 ] );
    }
    // Buffers must be different
    LAMA_LOG_INFO( logger, "async shift : using same send and recv buffer should fail" );
    LAMA_CHECK_THROW( delete comm->shiftAsync( sendBuffer, sendBuffer, 1 ), Exception );
    // We also verify that the exception is thrown before changing the send Buffer
    BOOST_CHECK_EQUAL( 2, sendBuffer.size() );
    LAMA_LOG_INFO( logger, "async shift : try to access send / receive buffer before synchronization" );
    unique_ptr<SyncToken> token( comm->shiftAsync( recvBuffer, sendBuffer, 1 ) );
    // read access on send buffer should be possible
    {
        ReadAccess<ValueType> sbuffer( sendBuffer );
        ValueType value0 = static_cast<ValueType>( rank );
        BOOST_CHECK_EQUAL( value0, sbuffer[ 0 ] );
    }
    LAMA_LOG_INFO( logger, *token << ": test for correct locks" );

    if ( !token->isSynchronized() )
    {
        // write access on send buffer should be locked
        LAMA_CHECK_THROW( WriteAccess<ValueType> sbuffer( sendBuffer ), Exception );
    }

    if ( !token->isSynchronized() )
    {
        // read access on recv buffer should be locked if communication is not finished yet
        LAMA_CHECK_THROW( ReadAccess<ValueType> recvBufferAccess( recvBuffer ), Exception );
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
LAMA_COMMON_TEST_CASE_TM_END();

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TM( CommunicatorTest, ValueType, bcastTest )
{
    LAMA_LOG_INFO( logger, "bcastTest<" << common::getScalarType<ValueType>() << ">" )
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
LAMA_COMMON_TEST_CASE_TM_END();

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TM( CommunicatorTest, ValueType, scatterTest )
{
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
LAMA_COMMON_TEST_CASE_TM_END();

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TM( CommunicatorTest, ValueType, scatterVTest )
{
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
LAMA_COMMON_TEST_CASE_TM_END();

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TM( CommunicatorTest, ValueType, gatherTest )
{
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
LAMA_COMMON_TEST_CASE_TM_END();

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TM( CommunicatorTest, ValueType, gatherVTest )
{
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
LAMA_COMMON_TEST_CASE_TM_END();

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TM( CommunicatorTest, ValueType, swapTest )
{
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
LAMA_COMMON_TEST_CASE_TM_END()

/* ----------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE( CommunicatorTest, writeAtTest )
{
    LAMA_WRITEAT_PTR_TEST( comm );
}
LAMA_COMMON_TEST_CASE_END()

/* ----------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_RUNNER( CommunicatorTest )
{
// disable inlining otherwise derived classes will not find it
#define LAMA_COMM_TEST( z, I, _ )           \
    swapTest<ARRAY_TYPE##I>();              \
    gatherTest<ARRAY_TYPE##I>();            \
    gatherVTest<ARRAY_TYPE##I>();           \
    scatterTest<ARRAY_TYPE##I>();           \
    scatterVTest<ARRAY_TYPE##I>();          \
    bcastTest<ARRAY_TYPE##I>();             \
    shiftTest<ARRAY_TYPE##I>();             \
    shiftASyncTest<ARRAY_TYPE##I>();        \
    updateHaloTest<ARRAY_TYPE##I>();        \
    // instantiate methods for all supported data types
    BOOST_PP_REPEAT( ARRAY_TYPE_CNT, LAMA_COMM_TEST, _ )
#undef LAMA_COMM_TEST
    bcastStringTest();
    buildHaloTest();
    allocatePlanTest();
    computeOwnersTest();
    CommunicatorCtrTest();
}
