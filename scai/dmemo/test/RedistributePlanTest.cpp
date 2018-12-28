/**
 * @file RedistributePlanTest.cpp
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
 * @brief Tests for the class RedistributePlan.
 * @author Thomas Brandes
 * @date 14.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/RedistributePlan.hpp>

#include <scai/common/test/TestMacros.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/hmemo/HostReadAccess.hpp>

#include <memory>
#include <algorithm>
#include <numeric>

using std::shared_ptr;

using namespace scai;
using namespace hmemo;
using namespace dmemo;
using namespace common;

using boost::test_tools::per_element;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( RedistributePlanTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.RedistributePlanTest" );

/* --------------------------------------------------------------------- */

/** Simple function to define a value in a distributed array. 
 *  It is used to verify correct values in the halo.
 */
template<typename ValueType>
static ValueType globalValue( IndexType globalIndex )
{
    // return static_cast<ValueType>( 2 * globalIndex + 1 );
    return static_cast<ValueType>( globalIndex );
}

template<typename ValueType>
static HArray<ValueType> distributedArray( const Distribution& dist )
{
    HArray<ValueType> localArray;  // local part of the distributed 'global' array

    // use own scope for write access to make sure that access is closed before return

    {
        IndexType localIndex = 0;   // running local index

        for ( auto& entry : hostWriteOnlyAccess( localArray, dist.getLocalSize() ) )
        {
            entry = globalValue<ValueType>( dist.local2Global( localIndex++ ) );
        }

    }  // filled the local array with 'global' values

    return localArray;    // each processor gets its local part
}

/* --------------------------------------------------------------------- */

template<typename ValueType>
static HArray<ValueType> distributedDenseMatrix( const Distribution& dist, const IndexType N )
{
    HArray<ValueType> localArray;  // local part of the distributed 'global' array

    // use own scope for write access to make sure that access is closed before return

    {
        IndexType localIndex = 0;   // running local row index
        IndexType colIndex = 0;     // running column index

        for ( auto& entry : hostWriteOnlyAccess( localArray, dist.getLocalSize() * N ) )
        {
            entry = ValueType( 2 * dist.local2Global( localIndex ) + 3 * colIndex );
            colIndex++;
            if ( colIndex == N )
            {
                colIndex = 0;
                localIndex++;
            }
        }

    }  // filled the local matrix with 'global' values

    return localArray;    // each processor gets its local part
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( redistributeTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    const IndexType N = 20;
    const IndexType CHUNK_SIZE = 2;

    // just define two arbitrary distributons

    auto sourceDist = blockDistribution( N );
    auto targetDist = cyclicDistribution( N, CHUNK_SIZE );

    // fill the 'distributed' arrays with 'same' values

    const auto sourceArray = distributedArray<ValueType>( *sourceDist );
    const auto expTargetArray = distributedArray<ValueType>( *targetDist );

    RedistributePlan plan( targetDist, sourceDist );

    HArray<ValueType> targetArray;
    plan.redistribute( targetArray, sourceArray );

    BOOST_REQUIRE_EQUAL( targetArray.size(), targetDist->getLocalSize() );

    BOOST_TEST( hostReadAccess( expTargetArray ) == hostReadAccess( targetArray ), per_element() );

    plan.reverse();

    HArray<ValueType> newSourceArray;
    plan.redistribute( newSourceArray, targetArray );
    BOOST_TEST( hostReadAccess( sourceArray ) == hostReadAccess( newSourceArray ), per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( redistributeNTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    const IndexType N = 20;
    const IndexType CHUNK_SIZE = 2;

    const IndexType M = 4;

    // just define two arbitrary distributons

    auto sourceDist = blockDistribution( N );
    auto targetDist = cyclicDistribution( N, CHUNK_SIZE );

    // fill the 'distributed' arrays with 'same' values

    const auto sourceArray = distributedDenseMatrix<ValueType>( *sourceDist, M );
    const auto expTargetArray = distributedDenseMatrix<ValueType>( *targetDist, M );

    RedistributePlan plan( targetDist, sourceDist );

    HArray<ValueType> targetArray;
    plan.redistributeN( targetArray, sourceArray, M );

    BOOST_REQUIRE_EQUAL( targetArray.size(), targetDist->getLocalSize() * M );

    BOOST_TEST( hostReadAccess( expTargetArray ) == hostReadAccess( targetArray ), per_element() );

    plan.reverse();

    HArray<ValueType> newSourceArray;
    plan.redistributeN( newSourceArray, targetArray, M );
    BOOST_TEST( hostReadAccess( sourceArray ) == hostReadAccess( newSourceArray ), per_element() );
}

/* --------------------------------------------------------------------- */

static HArray<IndexType> distributedOffsets( const Distribution& dist )
{
    const IndexType localN = dist.getLocalSize();

    HArray<IndexType> localArray;  // local part of the distributed 'global' array

    {
        auto wLocal = hostWriteOnlyAccess( localArray, localN + 1 );

        // use own scope for write access to make sure that access is closed before return

        for ( IndexType localIndex = 0; localIndex < localN; ++localIndex )
        {
            wLocal[localIndex] = dist.local2Global( localIndex ) % 5;
        }
    }

    localArray.resize( localN );

    utilskernel::HArrayUtils::scan1( localArray );

    return localArray;    // each processor gets its local part
}

template<typename ValueType>
static HArray<ValueType> distributedRaggedArray( const Distribution& dist, const HArray<IndexType> offsets )
{
    HArray<ValueType> localRaggedArray;

    const IndexType localN = dist.getLocalSize();

    SCAI_ASSERT_EQ_ERROR( offsets.size(), localN + 1, "serious mismatch" )

    {
        auto rOffsets  = hostReadAccess( offsets );
        auto wLocal    = hostWriteOnlyAccess( localRaggedArray, rOffsets[localN] );

        IndexType localIndex = 0;

        for ( IndexType i = 0; i < localN; ++i )
        {
            IndexType sizeI = rOffsets[i + 1] - rOffsets[i];
            IndexType globalI = dist.local2Global( i );

            for ( IndexType j = 0; j < sizeI; ++j )
            {
                wLocal[localIndex++] = globalI + 2 * j;
            }
        }

        SCAI_ASSERT_EQ_ERROR( localIndex, rOffsets[localN], "serious mismatch" )
    }

    return localRaggedArray;
}

BOOST_AUTO_TEST_CASE( redistributeVTest )
{
    typedef SCAI_TEST_TYPE ValueType;

    const IndexType N = 30;
    const IndexType CHUNK_SIZE = 2;

    // just define two arbitrary distributons

    auto sourceDist = blockDistribution( N );
    auto targetDist = cyclicDistribution( N, CHUNK_SIZE );

    // fill the 'distributed' arrays with 'same' values

    const auto sourceOffsets = distributedOffsets( *sourceDist );
    SCAI_LOG_DEBUG( logger, "source offsets = " << printIt( sourceOffsets ) )

    const auto targetOffsets = distributedOffsets( *targetDist );
    SCAI_LOG_DEBUG( logger, "target offsets = " << printIt( targetOffsets ) )

    const auto sourceArray = distributedRaggedArray<ValueType>( *sourceDist, sourceOffsets );
    SCAI_LOG_DEBUG( logger, "source array = " << printIt( sourceArray ) )
    const auto expTargetArray = distributedRaggedArray<ValueType>( *targetDist, targetOffsets );
    SCAI_LOG_DEBUG( logger, "exp target array = " << printIt( expTargetArray ) )

    RedistributePlan plan( targetDist, sourceDist );

    HArray<ValueType> targetArray;
    plan.redistributeV( targetArray, targetOffsets, sourceArray, sourceOffsets );

    BOOST_TEST( hostReadAccess( expTargetArray ) == hostReadAccess( targetArray ), per_element() );

    plan.reverse();

    HArray<ValueType> newSourceArray;
    plan.redistributeV( newSourceArray, sourceOffsets, targetArray, targetOffsets );
    BOOST_TEST( hostReadAccess( sourceArray ) == hostReadAccess( newSourceArray ), per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( redistributorTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    IndexType size = comm->getSize();
    IndexType rank = comm->getRank();

    const IndexType N = 2 * size;

    //              P = 0     P = 1   P = 2   P = 3
    //  dist1        0, 7     1, 6    2, 5    3, 4
    //  dist2        1, 4     2, 5    3, 6    0, 7

    HArray<IndexType> myIndexes1  = { rank, 2 * size - 1 - rank };
    HArray<IndexType> myIndexes2  = { rank + 1 == size ? 0 : rank + 1, size + rank };

    // the following array contains the new owners of myIndexes1, first element goes to left processor

    HArray<IndexType> newOwners = { rank == 0 ? size - 1 : rank - 1, size - 1 - rank };

    {
        auto rI1 = hostReadAccess( myIndexes1 );
        auto rI2 = hostReadAccess( myIndexes2 );
        auto rO = hostReadAccess( newOwners );

        SCAI_LOG_DEBUG( logger, *comm << ", indexes1 = " << rI1[0] << ", " << rI1[1]
                        << ", indexes2 = " << rI2[0] << ", " << rI2[1]
                        << ", new owners = " << rO[0] << ", " << rO[1] )
    }

    auto sourceDist = std::make_shared<GeneralDistribution>( N, myIndexes1, true, comm );
    auto targetDist = std::make_shared<GeneralDistribution>( N, myIndexes2, true, comm );

    RedistributePlan r1( targetDist, sourceDist );
    RedistributePlan r2( newOwners, sourceDist );

    // both redistributions should be the same, we prove by redistribution of some data

    HArray<IndexType> sourceData( { 100 + rank, 200 + rank } );

    HArray<IndexType> targetData1;
    HArray<IndexType> targetData2;

    r1.redistribute( targetData1, sourceData );
    r2.redistribute( targetData2, sourceData );

    {
        auto rS = hostReadAccess( sourceData );
        auto rT1 = hostReadAccess( targetData1 );
        auto rT2 = hostReadAccess( targetData2 );

        SCAI_LOG_DEBUG( logger, *comm << ", source = " << rS[0] << ", " << rS[1]
                        << ", target1 = " << rT1[0] << ", " << rT1[1]
                        << ", target2 = " << rT2[0] << ", " << rT2[1] )
    }

    BOOST_TEST( hostReadAccess( targetData1 ) == hostReadAccess( targetData2 ), per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    const IndexType N = 100;
    const IndexType CHUNK_SIZE = 2;

    // just define two arbitrary distributons

    auto sourceDist = blockDistribution( N );
    auto targetDist = cyclicDistribution( N, CHUNK_SIZE );

    RedistributePlan r( targetDist, sourceDist );
    std::ostringstream out;
    out << r ;
    BOOST_CHECK( out.str().length() >  0 );
}

BOOST_AUTO_TEST_CASE( constructorTest )
{
    using boost::test_tools::per_element;

    // Helper struct to make assigning test data for individual processors simpler
    struct ExpectedResult
    {
        std::vector<IndexType> keepSourceIndexes;
        std::vector<IndexType> keepTargetIndexes;
        std::vector<IndexType> exchangeSourceIndexes;
        std::vector<IndexType> exchangeTargetIndexes;
        std::vector<IndexType> local2global;
    };

    const auto checkRedistributePlanAgainstExpected = [] ( const RedistributePlan & redist, const ExpectedResult & expected )
    {
        BOOST_TEST( hostReadAccess( redist.getKeepSourceIndexes() ) == expected.keepSourceIndexes, per_element() );
        BOOST_TEST( hostReadAccess( redist.getKeepTargetIndexes() ) == expected.keepTargetIndexes, per_element() );
        BOOST_TEST( hostReadAccess( redist.getExchangeSourceIndexes() ) == expected.exchangeSourceIndexes, per_element() );
        BOOST_TEST( hostReadAccess( redist.getExchangeTargetIndexes() ) == expected.exchangeTargetIndexes, per_element() );

        const auto sourceDist = redist.getSourceDistributionPtr();
        const auto targetDist = redist.getTargetDistributionPtr();

        // Collect local2global from target distribution
        auto local2global = std::vector<IndexType>( targetDist->getLocalSize() );
        std::iota( local2global.begin(), local2global.end(), 0 );
        std::transform( local2global.cbegin(), local2global.cend(), local2global.begin(),
                        [targetDist] ( IndexType localIndex )
        {
            return targetDist->local2Global( localIndex );
        } );
        BOOST_TEST( local2global == expected.local2global, per_element() );
        BOOST_TEST( targetDist->getGlobalSize() == sourceDist->getGlobalSize() );
    };

    const auto comm = Communicator::getCommunicatorPtr();
    const auto rank = comm->getRank();
    const auto numPartitions = comm->getSize();
    auto expected = ExpectedResult();

    if ( numPartitions == 1 )
    {
        const auto sourceDist = DistributionPtr ( new CyclicDistribution( 4, 1, comm ) );
        const auto newOwnersOfLocal = HArray<IndexType> { 0, 0, 0, 0 };

        expected.keepSourceIndexes = { 0, 1, 2, 3 };
        expected.keepTargetIndexes = { 0, 1, 2, 3 };
        expected.exchangeSourceIndexes = { };
        expected.exchangeTargetIndexes = { };
        expected.local2global = { 0, 1, 2, 3 };

        const auto redist = RedistributePlan( newOwnersOfLocal, sourceDist );
        checkRedistributePlanAgainstExpected( redist, expected );
    }
    else if ( numPartitions == 2 )
    {
        const auto sourceDist = DistributionPtr ( new CyclicDistribution( 8, 1, comm ) );

        HArray<IndexType> newOwnersOfLocal;

        // Input data
        switch ( rank )
        {
            case 0:
                newOwnersOfLocal = { 1, 0, 1, 1 };
                break;
            case 1:
                newOwnersOfLocal = { 1, 0, 0, 1 };
                break;
        }

        // Expected data
        if ( rank == 0 )
        {
            expected.keepSourceIndexes = { 1 };
            expected.keepTargetIndexes = { 0 };
            expected.exchangeSourceIndexes = { 0, 2, 3 };
            expected.exchangeTargetIndexes = { 1, 2 };
            expected.local2global = { 2, 3, 5 };
        }
        else if ( rank == 1 )
        {
            expected.keepSourceIndexes = { 0, 3 };
            expected.keepTargetIndexes = { 1, 4 };
            expected.exchangeSourceIndexes = { 1, 2 };
            expected.exchangeTargetIndexes = { 0, 2, 3 };
            expected.local2global = { 0, 1, 4, 6, 7 };
        }

        const auto redist = RedistributePlan( newOwnersOfLocal, sourceDist );

        checkRedistributePlanAgainstExpected( redist, expected );
    }
    else if ( numPartitions == 3 )
    {
        const auto sourceDist = DistributionPtr ( new CyclicDistribution( 13, 1, comm ) );

        HArray<IndexType> newOwnersOfLocal;

        switch ( rank )
        {
            case 0:
                newOwnersOfLocal = { 2, 0, 1, 0, 2};
                break;
            case 1:
                newOwnersOfLocal = { 2, 2, 2, 0 };
                break;
            case 2:
                newOwnersOfLocal = { 2, 2, 0, 1 };
                break;
        }

        if ( rank == 0 )
        {
            expected.keepSourceIndexes = { 1, 3 };
            expected.keepTargetIndexes = { 0, 2 };
            expected.exchangeSourceIndexes = { 2, 0, 4 };
            expected.exchangeTargetIndexes = { 3, 1 };
            expected.local2global = { 3, 8, 9, 10 };
        }
        else if ( rank == 1 )
        {
            expected.keepSourceIndexes = { };
            expected.keepTargetIndexes = { };
            expected.exchangeSourceIndexes = { 3, 0, 1, 2 };
            expected.exchangeTargetIndexes = { 0, 1 };
            expected.local2global = { 6, 11 };
        }
        else if ( rank == 2 )
        {
            expected.keepSourceIndexes = { 0, 1 };
            expected.keepTargetIndexes = { 2, 4 };
            expected.exchangeSourceIndexes = { 2, 3 };
            expected.exchangeTargetIndexes = { 0, 6, 1, 3, 5 };
            expected.local2global = { 0, 1, 2, 4, 5, 7, 12 };
        }

        const auto redist = RedistributePlan( newOwnersOfLocal, sourceDist );

        checkRedistributePlanAgainstExpected( redist, expected );
    }
    else
    {
        BOOST_TEST_MESSAGE( "No test data for " << numPartitions << " partitions." );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
