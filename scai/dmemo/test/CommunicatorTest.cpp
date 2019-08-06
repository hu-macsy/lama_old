/**
 * @file CommunicatorTest.cpp
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
 * @brief Test methods for the communicator
 * @author Thomas Brandes
 * @date 09.05.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/logging.hpp>

#include <scai/dmemo.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>

#include <scai/tasking.hpp>
#include <scai/hmemo.hpp>
#include <scai/utilskernel.hpp>

#include <scai/common/exception/Exception.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/test/TestMacros.hpp>

#include <memory>

using namespace scai;
using namespace hmemo;
using namespace dmemo;
using namespace tasking;
using namespace utilskernel;

using common::Exception;
using std::unique_ptr;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CommunicatorTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.CommunicatorTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( basicTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );
    SCAI_LOG_INFO( logger, "basicTest for comm = " << *comm )

    BOOST_CHECK( common::Utils::validIndex( comm->getRank(), comm->getSize() ) );
    std::stringstream out;
    out << *comm;
    BOOST_CHECK( out.str().length() > 0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( computeOwnersTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();
    IndexType n = 17;

    HArray<IndexType> localIndexes;

    IndexType first = static_cast<IndexType>( rank ) * n;
    IndexType inc   = 1;

    HArrayUtils::setSequence( localIndexes, first, inc, n );

    GeneralDistribution dist( n * size, localIndexes, false, comm );

    HArray<IndexType> nonLocalIndexes;

    IndexType pos = 0;

    {
        hmemo::WriteOnlyAccess<IndexType> wNonLocalIndexes( nonLocalIndexes, ( size - 1 ) * n );

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
                    wNonLocalIndexes[pos++] = p * n + i;
                }
            }
        }
    }

    BOOST_CHECK_EQUAL( pos, nonLocalIndexes.size() );

    HArray<PartitionId> owners;

    comm->computeOwners( owners, dist, nonLocalIndexes );

    pos = 0;

    hmemo::ReadAccess<PartitionId> rOwners( owners );

    for ( PartitionId p = 0; p < size; ++p )
    {
        if ( p != rank )
        {
            for ( IndexType i = 0; i < n; ++i )
            {
                BOOST_CHECK_EQUAL( p, rOwners[pos++] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( bcastStringTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );
    std::string val = "Dummy";
    PartitionId root = 0;

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

BOOST_AUTO_TEST_CASE_TEMPLATE( shiftTest, ValueType, scai_numeric_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();

    // Idea of this shift Test:
    // - allocate on each processor an array with one element for each processor
    // - shift this array around all processors and each processor writes one value at its rank
    // - verify that each processor has written the right value

    if ( size > 1 )
    {
        const IndexType vectorSize = size;
        HArray<ValueType> sendBuffer( vectorSize, static_cast<ValueType>( -1 ) );
        HArray<ValueType> recvBuffer( vectorSize );

        for ( PartitionId rounds = 0; rounds < size; ++rounds )
        {
            {
                WriteAccess<ValueType> sendBufferAccess( sendBuffer );
                sendBufferAccess[rank] = static_cast<ValueType>( rank );
            }
            // Note: recvBuffer will be allocated with same size as sendBuffer
            comm->shiftArray( recvBuffer, sendBuffer, 1 );
            sendBuffer.swap( recvBuffer );
        }

        {
            ReadAccess<ValueType> sendBufferAccess( sendBuffer );

            for ( IndexType i = 0; i < size; ++i )
            {
                ValueType value = static_cast<ValueType>( i );
                BOOST_CHECK_EQUAL( value, sendBufferAccess[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( joinTest, ValueType, scai_numeric_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );

    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();

    ValueType val = static_cast<ValueType>( rank );

    HArray<ValueType> localArray( IndexType( rank ) + 1, val );
    HArray<ValueType> globalArray;

    comm->joinArray( globalArray, localArray );

    // globalArray: 0 1 1 2 2 2 3 3 3 3 4 4 4 4 4 ...

    IndexType expectedSize = comm->getSize();
    expectedSize = expectedSize * ( expectedSize + 1 ) / 2;

    BOOST_REQUIRE_EQUAL( expectedSize, globalArray.size() );

    ReadAccess<ValueType> rGlobal( globalArray );

    IndexType pos = 0;

    for ( IndexType rank = 0; rank < size; ++rank )
    {
        ValueType expectedVal = static_cast<ValueType>( rank );

        for ( IndexType k = 0; k <= rank; ++k )
        {
            BOOST_CHECK_EQUAL( expectedVal, rGlobal[pos] );
            pos++;
        }
    }

    BOOST_CHECK_EQUAL( pos, expectedSize );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( shiftAsyncTest, ValueType, scai_numeric_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );
    PartitionId rank = comm->getRank();
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
    BOOST_CHECK_EQUAL( IndexType( 2 ), recvBuffer.size() );
    {
        ReadAccess<ValueType> rbuffer( recvBuffer );
    }
    // dir = 0: should be like assignment
    delete comm->shiftAsync( recvBuffer, sendBuffer, 0 );
    SCAI_LOG_INFO( logger, "async shift dir = 0, self assing" );
    BOOST_CHECK_EQUAL( IndexType( 2 ), recvBuffer.size() );
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
    BOOST_CHECK_EQUAL( IndexType( 2 ), sendBuffer.size() );
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
    BOOST_CHECK_EQUAL( IndexType( 2 ), recvBuffer.size() );
    {
        ReadAccess<ValueType> rbuffer( recvBuffer );
        ValueType value0 = static_cast<ValueType>( comm->getNeighbor( -1 ) );
        ValueType value1 = static_cast<ValueType>( rank );
        BOOST_CHECK_EQUAL( value0, rbuffer[ 0 ] );
        BOOST_CHECK_EQUAL( value1, rbuffer[ 1 ] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( bcastTest, ValueType, scai_numeric_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();
    SCAI_LOG_INFO( logger, "bcastTest<" << common::getScalarType<ValueType>() << ">" )
    IndexType N = 5;
    ValueType dummyVal = 13;
    unique_ptr<ValueType[]> vector( new ValueType[N + 1] );
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

BOOST_AUTO_TEST_CASE( bcastFailTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    PartitionId root = comm->getSize() + 1;

    IndexType dummyVal = 13;

    // Illegal root for bcast should throw an exception

    BOOST_CHECK_THROW (
    {
        comm->bcast( &dummyVal, 1, root );
    }, common::Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scatterTest, ValueType, scai_numeric_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();
    const PartitionId root = 0;
    IndexType n = 2;
    IndexType allN = 0; // only root will have full size

    if ( rank == root )
    {
        allN = size * n;
    }

    ValueType dummyVal = 13;
    unique_ptr<ValueType[]> myvals( new ValueType[n + 1] );
    unique_ptr<ValueType[]> allvals( new ValueType[allN] );
    myvals[0] = 0;
    myvals[1] = 1;
    myvals[2] = static_cast<ValueType>( size ) + dummyVal;

    if ( rank == root )
    {
        // fill the send data
        for ( IndexType i = 0; i < size; i++ )
        {
            allvals[2 * i] = static_cast<ValueType>( i );
            allvals[2 * i + 1] = dummyVal + static_cast<ValueType>( i );
        }
    }

    comm->scatter( myvals.get(), n, root, allvals.get() );
    // also verify that myvals has not been overwritten at end
    BOOST_CHECK_EQUAL( myvals[0], static_cast<ValueType>( rank ) );
    BOOST_CHECK_EQUAL( myvals[1], dummyVal + static_cast<ValueType>( rank ) );
    BOOST_CHECK_EQUAL( myvals[2], dummyVal + static_cast<ValueType>( size ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scatterVTest, ValueType, scai_numeric_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();
    const PartitionId root = 0;
    IndexType n = rank; // number of elements I receive
    IndexType allN = 0;

    if ( rank == root )
    {
        allN = size * ( size - 1 ) / 2;
    }

    ValueType dummyVal = 13;
    unique_ptr<ValueType[]> myvals( new ValueType[n + 1] );
    unique_ptr<ValueType[]> allvals( new ValueType[allN] );
    unique_ptr<IndexType[]> sizes( new IndexType[size] );

    for ( IndexType i = 0; i <= n; i++ )
    {
        myvals[i] = dummyVal;
    }

    if ( rank == root )
    {
        // fill the send data
        IndexType offset = 0;

        for ( PartitionId p = 0; p < size; p++ )
        {
            for ( PartitionId i = 0; i < p; i++ )
            {
                allvals[offset] = dummyVal + static_cast<ValueType>( i + 1 );
                ++offset;
            }

            sizes[p] = p;
        }
    }

    comm->scatterV( myvals.get(), n, root, allvals.get(), sizes.get() );

    // also verify that myvals has not been overwritten at end
    for ( IndexType i = 0; i < n; i++ )
    {
        BOOST_CHECK_EQUAL( myvals[i], dummyVal + static_cast<ValueType>( i + 1 ) );
    }

    BOOST_CHECK_EQUAL( myvals[n], dummyVal );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gatherTest, ValueType, scai_numeric_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();
    const PartitionId root = 0;
    IndexType allN = 0; // only root will have full size

    if ( rank == root )
    {
        allN = size;
    }

    ValueType dummyVal = -1;
    unique_ptr<ValueType[]> myvals( new ValueType[2] );
    unique_ptr<ValueType[]> allvals( new ValueType[allN] );

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

BOOST_AUTO_TEST_CASE_TEMPLATE( maxLocTest, ValueType, scai_array_test_types )
{
    using common::Math;
    using common::TypeTraits;

    typedef typename TypeTraits<ValueType>::RealType RealType;   // only here comparison

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    Math::srandom( 1751 + comm->getRank() * 17 );

    // test default ( kind == 0 ) and specific ( kind == 1 ) implementation

    for ( int kind = 0; kind < 2; ++kind )
    {
        // test it for each processor to be the root

        for ( PartitionId root = 0; root < comm->getSize(); ++root )
        {
            const IndexType N = 10;

            auto vals = utilskernel::randomHArray<RealType>( N, 1 );
            RealType localMax = vals[0];
            IndexType localMaxLoc = 0;

            for ( IndexType i = 0; i < N; ++i )
            {
                RealType v = vals[i];

                if ( v > localMax )
                {
                    localMax = v;
                    localMaxLoc = i;
                }
            }

            SCAI_LOG_INFO( logger, *comm << ": checkMaxLoc, local " << localMax << " @ " << localMaxLoc )

            RealType globalMax1 = comm->max( localMax );

            BOOST_CHECK( Math::abs( globalMax1 ) >= Math::abs( localMax ) );

            IndexType globalMaxLoc = localMaxLoc;
            RealType globalMax    = localMax;

            if ( kind == 0 )
            {
                comm->maxlocDefault( globalMax, globalMaxLoc, root );
            }
            else
            {
                comm->maxloc( globalMax, globalMaxLoc, root );
            }

            comm->bcast( &globalMax, 1, root );
            comm->bcast( &globalMaxLoc, 1, root );

            BOOST_CHECK_EQUAL( globalMax1, globalMax );

            SCAI_LOG_INFO( logger, *comm << ": checkMaxLoc, global " << globalMax << " @ " << globalMaxLoc )

            BOOST_CHECK( Math::abs( globalMax ) >= Math::abs( localMax ) );

            bool any = globalMaxLoc == localMaxLoc;

            any = comm->any( any );

            BOOST_CHECK( any );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( minLocTest, ValueType, scai_array_test_types )
{
    using common::Math;

    typedef typename common::TypeTraits<ValueType>::RealType RealType;   // for comparison

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    Math::srandom( 11171 + comm->getRank() * 17 );

    // test default ( kind == 0 ) and specific ( kind == 1 ) implementation
    for ( int kind = 0; kind < 2; ++kind )
    {
        // test it for each processor to be the root
        for ( PartitionId root = 0; root < comm->getSize(); ++root )
        {
            IndexType N = 5;
            auto vals = randomHArray<RealType>( N, 1 );
            RealType localMin = vals[0];
            IndexType localMinLoc = 0;

            for ( IndexType i = 0; i < N; ++i )
            {
                RealType v = vals[i];

                if ( v < localMin )
                {
                    localMin = v;
                    localMinLoc = i;
                }
            }

            SCAI_LOG_INFO( logger, *comm << ": checkMinLoc, local " << localMin << " @ " << localMinLoc )

            RealType globalMin1 = comm->min( localMin );

            bool error1 = localMin < globalMin1;

            BOOST_CHECK( !error1 );

            // BOOST_CHECK( Math::abs( globalMax1 ) >= Math::abs( localMax ) );

            IndexType globalMinLoc = localMinLoc;
            RealType globalMin = localMin;

            if ( kind == 0 )
            {
                comm->minlocDefault( globalMin, globalMinLoc, root );
            }
            else
            {
                comm->minloc( globalMin, globalMinLoc, root );
            }

            comm->bcast( &globalMin, 1, root );
            comm->bcast( &globalMinLoc, 1, root );

            BOOST_CHECK_EQUAL( globalMin1, globalMin );

            SCAI_LOG_INFO( logger, *comm << ": checkMinLoc, global " << globalMin << " @ " << globalMinLoc )

            bool error = RealType( localMin ) < RealType( globalMin );

            BOOST_CHECK( !error );

            bool any = globalMinLoc == localMinLoc;

            any = comm->any( any );

            BOOST_CHECK( any );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( minTest, ValueType, scai_array_test_types )
{
    using common::Math;

    typedef typename common::TypeTraits<ValueType>::RealType RealType;   // for comparison

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    Math::srandom( 1751 + comm->getRank() * 17 );

    IndexType N = 5;
    auto vals = randomHArray<RealType>( N, 10 );
    RealType localMin = vals[0];

    for ( IndexType i = 0; i < N; ++i )
    {
        RealType v = vals[i];
        v = Math::abs( v );

        if ( v > localMin )
        {
            localMin = v;
        }
    }

    RealType globalMin = comm->min( localMin );

    BOOST_CHECK( Math::abs( globalMin ) <= Math::abs( localMin ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sumTest, ValueType, scai_array_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    ValueType i = common::TypeTraits<ValueType>::imaginaryUnit();

    ValueType rank = static_cast<ValueType>( comm->getRank() );
    ValueType size = static_cast<ValueType>( comm->getSize() );

    ValueType val = ( rank + 1 ) + rank * i;

    ValueType expected = size * ( size + 1 ) / 2 + size * ( size - 1 ) / 2 * i;

    ValueType sum = comm->sum( val );

    SCAI_LOG_DEBUG( logger, "sumTest<" << common::TypeTraits<ValueType>::id() << ">: comm = "
                    << *comm << ", val = " << val << ", sum = " << sum << ", expected = " << expected )

    BOOST_CHECK_EQUAL( expected, sum );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sumArrayTest, ValueType, scai_array_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    PartitionId size = comm->getSize();

    common::Math::srandom( 17511 );  // same vals on each processor

    IndexType N = 5;
    auto vals = randomHArray<ValueType>( N, 1 );

    HArray<ValueType> saveVals( vals );

    comm->sumArray( vals );

    HArrayUtils::setScalar<ValueType>( saveVals, size, common::BinaryOp::MULT );

    // be carful: values are not exactly the same

    BOOST_CHECK( HArrayUtils::maxDiffNorm( vals, saveVals ) < 0.0001 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( scanTest, ValueType, scai_array_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    struct
    {
        inline ValueType operator() ( const PartitionId r )
        {
            return ValueType( r * 2 + 1 );
        }
    } f;

    PartitionId rank = comm->getRank();   // used to compute my value

    for ( int kind = 0; kind < 2; ++kind )
    {
        ValueType v = f( rank );    // ValueType v = rank * 2 + 1;

        ValueType scanV = kind == 0 ? comm->scanDefault( v ) : comm->scan( v );

        ValueType expected = 0;

        for ( IndexType i = 0; i <= rank; ++i )
        {
            expected += f( i );  // note: is inclusive scan
        }
    
        SCAI_LOG_DEBUG( logger, *comm << ": scan, kind = " << kind << ", v = " << v << ", scanV = " 
                                << scanV << ", expected = " << expected )

        BOOST_CHECK_EQUAL( expected, scanV );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gatherVTest, ValueType, scai_numeric_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();
    const PartitionId root = 0;
    IndexType n = rank; // number of elements I send
    IndexType allN = 0;

    if ( rank == root )
    {
        allN = size * ( size - 1 ) / 2;
    }

    ValueType dummyVal = 13;
    unique_ptr<ValueType[]> myvals( new ValueType[n + 1] );
    unique_ptr<ValueType[]> allvals( new ValueType[allN] );
    unique_ptr<IndexType[]> sizes( new IndexType[size] );

    for ( IndexType i = 0; i < n; i++ )
    {
        myvals[i] = static_cast<ValueType>( rank );
    }

    myvals[n] = dummyVal;

    if ( rank == root )
    {
        // fill the recv data
        IndexType offset = 0;

        for ( PartitionId p = 0; p < size; p++ )
        {
            for ( IndexType i = 0; i < p; i++ )
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
        IndexType offset = 0;

        for ( PartitionId p = 0; p < size; ++p )
        {
            for ( PartitionId i = 0; i < p; ++i )
            {
                BOOST_CHECK_EQUAL( allvals[offset], static_cast<ValueType>( p ) );
                ++offset;
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( all2allTest, ValueType, scai_array_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    struct
    {
        inline ValueType operator() ( const PartitionId source, const PartitionId target )
        {
            return ValueType( source * 10 + target + 3 );
        }
    } f;

    PartitionId rank = comm->getRank();   // used to compute my value
    PartitionId size = comm->getSize();   // used for size of buffers

    SCAI_LOG_INFO( logger, "all2allTest<" << common::TypeTraits<ValueType>::id() << ">" )

    std::vector<ValueType> sendBuffer( size );
    std::vector<ValueType> recvBuffer( size );
    std::vector<ValueType> expRecvBuffer( size );

    // fill the send buffer and expected recv buffer
 
    for ( PartitionId i = 0; i < size; ++i )
    {
        sendBuffer[i] = f( rank, i );
        expRecvBuffer[i] = f( i, rank );
    }

    comm->all2all( recvBuffer.data(), sendBuffer.data() );

    BOOST_TEST( recvBuffer == expRecvBuffer, boost::test_tools::per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( all2allvTest, ValueType, scai_array_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    // functional computes for ( source, target ) the number of values to send

    struct
    {
        inline IndexType operator() ( const PartitionId source, const PartitionId target )
        {
            return 1 + ( source % 2 ) + ( target % 3 );
        }
    } n;

    // functional computes for ( source, target, k ) the k-th value sent from source to target

    struct
    {
        inline ValueType operator() ( const PartitionId source, const PartitionId target, const IndexType k )
        {
            return ValueType( source * 100 + target + 12 - k );
        }
    } f;

    PartitionId rank = comm->getRank();   // used to compute my value
    PartitionId size = comm->getSize();   // used for size of buffers

    std::vector<ValueType> sendValues;

    std::vector<IndexType> sendSizes( size );
    std::vector<IndexType> recvSizes( size );

    // compute the send values

    for ( IndexType i = 0; i < size; ++i )
    {
        sendSizes[i] = n( rank, i );

        for ( IndexType k = 0; k < sendSizes[i]; ++k )
        {
            sendValues.push_back( f( rank, i, k ) );
        }
    }

    comm->all2all( recvSizes.data(), sendSizes.data() );

    IndexType recvSizeSum = 0;

    for ( IndexType i = 0; i < size; ++i )
    { 
        auto nv = recvSizes[i];
        BOOST_CHECK_EQUAL( nv, n( i, rank ) );
        recvSizeSum += nv;
    }

    std::vector<ValueType> recvValues( recvSizeSum, ValueType( 0 ) );

    // set up pointers for send and receive buffers

    std::vector<ValueType*> recvBuffer( size );
    std::vector<const ValueType*> sendBuffer( size );

    IndexType offsetSend = 0;
    IndexType offsetRecv= 0;

    for ( IndexType i = 0; i < size; ++i )
    {
        sendBuffer[i] = &sendValues[ offsetSend ];
        recvBuffer[i] = &recvValues[ offsetRecv];
        offsetSend += sendSizes[i];
        offsetRecv += recvSizes[i];
    }

    BOOST_REQUIRE_EQUAL( offsetRecv, recvSizeSum );
    BOOST_REQUIRE_EQUAL( offsetSend, static_cast<IndexType>( sendValues.size() ) );

    SCAI_LOG_INFO( logger, "all2allvTest<" << common::TypeTraits<ValueType>::id() << ">"
                           << ", total send = " << sendValues.size() << ", total recv = " << recvValues.size() )

    comm->all2allv( recvBuffer.data(), recvSizes.data(), sendBuffer.data(), sendSizes.data() );

    std::vector<ValueType> expRecvValues;

    for ( IndexType i = 0; i < size; ++i )
    {
        for ( IndexType k = 0; k < recvSizes[i]; ++k )
        {
            expRecvValues.push_back( f( i, rank, k ) );
        }
    }

    BOOST_TEST( recvValues == expRecvValues, boost::test_tools::per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( swapTest, ValueType, scai_numeric_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();
    IndexType n = 10;
    unique_ptr<ValueType[]> vector( new ValueType[n] );

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

BOOST_AUTO_TEST_CASE( nodeTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    PartitionId nodeSize = comm->getNodeSize();
    PartitionId nodeRank = comm->getNodeRank();

    BOOST_CHECK( nodeRank < nodeSize );

    Communicator::ThreadSafetyLevel level = comm->getThreadSafetyLevel();

    BOOST_CHECK( int( level ) > 0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    std::ostringstream outDerived;

    outDerived << *comm;

    BOOST_CHECK( outDerived.str().length() > 0 );

    // verify that a derived communicator class has overridden the
    // default implementation of the base class Distriution

    std::ostringstream outBase;
    comm->Communicator::writeAt( outBase );
    BOOST_CHECK( outDerived.str() != outBase.str() );
}


/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( procArrayTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    PartitionId np = comm->getSize();

    PartitionId procArray[3];
    PartitionId posArray[3];

    bool replace = true;

    common::Settings::putEnvironment( "SCAI_NP", "1x3x2", replace );

    Communicator::getUserProcArray( procArray );

    BOOST_CHECK_EQUAL( procArray[0], PartitionId( 1 ) );
    BOOST_CHECK_EQUAL( procArray[1], PartitionId( 3 ) );
    BOOST_CHECK_EQUAL( procArray[2], PartitionId( 2 ) );

    common::Settings::putEnvironment( "SCAI_NP", "2", replace );

    Communicator::getUserProcArray( procArray );

    BOOST_CHECK_EQUAL( procArray[0], PartitionId( 2 ) );
    BOOST_CHECK_EQUAL( procArray[1], PartitionId( 0 ) );
    BOOST_CHECK_EQUAL( procArray[2], PartitionId( 0 ) );

    common::Settings::putEnvironment( "SCAI_NP", "5_1", replace );

    Communicator::getUserProcArray( procArray );

    BOOST_CHECK_EQUAL( procArray[0], PartitionId( 5 ) );
    BOOST_CHECK_EQUAL( procArray[1], PartitionId( 1 ) );
    BOOST_CHECK_EQUAL( procArray[2], PartitionId( 0 ) );

    common::Settings::putEnvironment( "SCAI_NP", "", replace );

    Communicator::factorize2( procArray, np, 16, 1 );

    BOOST_CHECK_EQUAL( np, procArray[0] * procArray[1] );

    if ( np < 5 )
    {
        // all processors should be dedicated to 1st factor

        BOOST_CHECK_EQUAL( np, procArray[0] );
    }

    posArray[0] = posArray[1] = posArray[2] = invalidPartition;

    comm->getGrid2Rank( posArray, procArray );

    BOOST_CHECK_EQUAL( comm->getRank(), posArray[1] * procArray[0] + posArray[0] );
    BOOST_CHECK_EQUAL( invalidPartition, posArray[2] );

    Communicator::factorize3( procArray, np, 1, 16, 1 );

    BOOST_CHECK_EQUAL( np, procArray[0] * procArray[1] * procArray[2] );

    if ( np < 5 )
    {
        // all processors should be dedicated to 2nd factor

        BOOST_CHECK_EQUAL( np, procArray[1] );
    }

    posArray[0] = posArray[1] = posArray[2] = invalidPartition;
    comm->getGrid3Rank( posArray, procArray );

    BOOST_CHECK_EQUAL( comm->getRank(), posArray[2] * procArray[0] * procArray[1] + posArray[1] * procArray[0] + posArray[0] );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( randomSeedTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    comm->setSeed( 5 );
    double val1 = common::Math::random<double>( 1 );
    comm->setSeed( 11 );
    common::Math::random<double>( 100 );
    comm->setSeed( 5 );
    double val3 = common::Math::random<double>( 1 );

    BOOST_CHECK_EQUAL( val1, val3 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */
