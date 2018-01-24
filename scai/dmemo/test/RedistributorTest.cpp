/**
 * @file RedistributorTest.cpp
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
 * @brief Tests for the class Redistributor.
 * @author Thomas Brandes
 * @date 14.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/Redistributor.hpp>

#include <scai/common/test/TestMacros.hpp>

#include <scai/hmemo/HostReadAccess.hpp>

#include <memory>
#include <algorithm>
#include <numeric>

using std::shared_ptr;

using namespace scai;
using namespace hmemo;
using namespace dmemo;
using namespace common;

/* --------------------------------------------------------------------- */

struct RedistributorTestConfig
{
    RedistributorTestConfig()
    {
        comm = Communicator::getCommunicatorPtr();
    }

    ~RedistributorTestConfig()
    {
    }

    CommunicatorPtr comm;
};

BOOST_FIXTURE_TEST_SUITE( RedistributorTest, RedistributorTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.RedistributorTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( redistributeTest )
{
    typedef SCAI_TEST_TYPE ValueType;
    IndexType size = 10;
    IndexType chunkSize = 1;
    shared_ptr<Distribution> distBlock( new BlockDistribution( size, comm ) );
    shared_ptr<Distribution> distCyclic( new CyclicDistribution( size, chunkSize, comm ) );
    IndexType blockLocalSize = distBlock->getLocalSize();
    IndexType cyclicLocalSize = distCyclic->getLocalSize();
    HArray<ValueType> myData1( blockLocalSize );
    {
        WriteAccess<ValueType> data ( myData1 );

        for ( IndexType i = 0; i < blockLocalSize; i++ )
        {
            data[i] = static_cast<ValueType>( 100 * comm->getRank() + i );
        }
    }
    HArray<ValueType> myData2( cyclicLocalSize );
    Redistributor r1( distCyclic, distBlock );
    Redistributor r2( distBlock, distCyclic );
    SCAI_LOG_DEBUG( logger, "redistribute 1" );
    r1.redistribute( myData2, myData1 );
    SCAI_LOG_DEBUG( logger, "redistribute 2" );
    r2.redistribute( myData1, myData2 );
    {
        ReadAccess<ValueType> data ( myData1 );

        for ( IndexType i = 0; i < blockLocalSize; i++ )
        {
            ValueType expected = static_cast<ValueType>( 100 * comm->getRank() + i );
            SCAI_CHECK_CLOSE( data[i], expected, 1 );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildRowPlansTest )
{
    typedef SCAI_TEST_TYPE ValueType;
    IndexType size = 10;
    IndexType chunkSize = 1;
    shared_ptr<Distribution> distBlock( new BlockDistribution( size, comm ) );
    shared_ptr<Distribution> distCyclic( new CyclicDistribution( size, chunkSize, comm ) );
    IndexType blockLocalSize = distBlock->getLocalSize();
    IndexType cyclicLocalSize = distCyclic->getLocalSize();
    HArray<ValueType> myData1( blockLocalSize );
    {
        WriteAccess<ValueType> data ( myData1 );

        for ( IndexType i = 0; i < blockLocalSize; i++ )
        {
            data[i] = static_cast<ValueType>( 100 * comm->getRank() + i );
        }
    }
    HArray<ValueType> myData2( cyclicLocalSize );

    Redistributor r1( distCyclic, distBlock );
    Redistributor r2( distBlock, distCyclic );

    // For testing buildRowPlans we set the sizes just as 1

    HArray<IndexType> targetSizes( myData2.size(), 1 );  // myData2 is target of r1
    HArray<IndexType> sourceSizes( myData1.size(), 1 );  // myData1 is target of r1

    r1.buildRowPlans( targetSizes, sourceSizes );

    SCAI_LOG_DEBUG( logger, "redistribute 1" );
    r1.redistribute( myData2, myData1 );
    SCAI_LOG_DEBUG( logger, "redistribute 2" );
    r2.redistribute( myData1, myData2 );

    {
        ReadAccess<ValueType> data ( myData1 );

        for ( IndexType i = 0; i < blockLocalSize; i++ )
        {
            ValueType expected = static_cast<ValueType>( 100 * comm->getRank() + i );
            SCAI_CHECK_CLOSE( data[i], expected, 1 );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    IndexType size = 10;
    IndexType chunkSize = 1;
    shared_ptr<Distribution> distBlock( new BlockDistribution( size, comm ) );
    shared_ptr<Distribution> distCyclic( new CyclicDistribution( size, chunkSize, comm ) );
    Redistributor r( distCyclic, distBlock );
    std::ostringstream out;
    out << r ;
    BOOST_CHECK( out.str().length() >  0 );
}

BOOST_AUTO_TEST_CASE( redistributorConstructorFromNewLocalOwnersTest )
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

    const auto checkRedistributorAgainstExpected = [] (const Redistributor & redist, const ExpectedResult & expected)
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
                        [targetDist] (IndexType localIndex) { return targetDist->local2global(localIndex); } );
        BOOST_TEST( local2global == expected.local2global, per_element() );
        BOOST_TEST( targetDist->getGlobalSize() == sourceDist->getGlobalSize() );
    };

    const auto comm = Communicator::getCommunicatorPtr();
    const auto rank = comm->getRank();
    const auto numPartitions = comm->getSize();
    auto expected = ExpectedResult();

    if ( numPartitions == 1 )
    {
         const auto sourceDist = DistributionPtr ( new CyclicDistribution( 4, 1, comm) );
         const auto newOwnersOfLocal = HArray<IndexType> { 0, 0, 0, 0 };

         expected.keepSourceIndexes = { 0, 1, 2, 3 };
         expected.keepTargetIndexes = { 0, 1, 2, 3 };
         expected.exchangeSourceIndexes = { };
         expected.exchangeTargetIndexes = { };
         expected.local2global = { 0, 1, 2, 3 };

         const auto redist = Redistributor( newOwnersOfLocal, sourceDist );
         checkRedistributorAgainstExpected( redist, expected );
    }
    else if ( numPartitions == 2 )
    {
        const auto sourceDist = DistributionPtr ( new CyclicDistribution( 8, 1, comm ) );

        HArray<IndexType> newOwnersOfLocal;

        // Input data
        switch ( rank )
        {
            case 0: newOwnersOfLocal = { 1, 0, 1, 1 }; break;
            case 1: newOwnersOfLocal = { 1, 0, 0, 1 }; break;
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

        const auto redist = Redistributor( newOwnersOfLocal, sourceDist );
        checkRedistributorAgainstExpected( redist, expected );
    }
    else if ( numPartitions == 3 )
    {
        const auto sourceDist = DistributionPtr ( new CyclicDistribution( 13, 1, comm ) );

        HArray<IndexType> newOwnersOfLocal;

        switch ( rank )
        {
            case 0: newOwnersOfLocal = { 2, 0, 1, 0, 2}; break;
            case 1: newOwnersOfLocal = { 2, 2, 2, 0 }; break;
            case 2: newOwnersOfLocal = { 2, 2, 0, 1 }; break;
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

        const auto redist = Redistributor( newOwnersOfLocal, sourceDist );
        checkRedistributorAgainstExpected( redist, expected );
    }
    else
    {
        BOOST_TEST_MESSAGE("No test data for " << numPartitions << " partitions.");
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
