/**
 * @file GeneralDistributionTest.cpp
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
 * @brief Contains the implementation of the class GeneralDistributionTest.
 * @author: Alexander Büchel, Thomas Brandes
 * @date 30.07.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/distribution/GeneralDistribution.hpp>

#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/test/distributed/DistributionTest.hpp>

using namespace scai::lama;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct GeneralDistributionTestConfig
{
    GeneralDistributionTestConfig()
    {
        comm = Communicator::getCommunicator( scai::lama::communicator::MPI );
        rank = comm->getRank();
        size = comm->getSize();
        elemsPerPartition = 10;
        globalSize = elemsPerPartition * size;

        for ( IndexType k = 0; k < elemsPerPartition; ++k )
        {
            localIndexes.push_back( k * size + rank );
        }

        dist = DistributionPtr( new GeneralDistribution( globalSize, localIndexes, comm ) );
    }

    ~GeneralDistributionTestConfig()
    {
        comm = CommunicatorPtr();
    }

    PartitionId rank;
    PartitionId size;

    IndexType elemsPerPartition;
    IndexType globalSize;

    std::vector<IndexType> localIndexes;

    DistributionPtr dist;
};

BOOST_FIXTURE_TEST_SUITE( GeneralDistributionTest, GeneralDistributionTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.GeneralDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( commonTestCases )
{
    DistributionTest disttest( dist );

    if ( base_test_case )
    {
        SCAI_LOG_INFO( logger, "Run test method " << testcase << " in GeneralDistributionTest." );
        DISTRIBUTION_COMMONTESTCASES( disttest );
    }
    else
    {
        disttest.runTests();
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( generalSizeTest )
{
    BOOST_CHECK( dist->getLocalSize() == elemsPerPartition );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    std::vector<IndexType> localIndexes;
    DistributionPtr generaldist1( new GeneralDistribution( 1, localIndexes, comm ) );
    DistributionPtr generaldist2( generaldist1 );
    DistributionPtr generaldist3( new GeneralDistribution( 1, localIndexes, comm ) );
    DistributionPtr generaldist4( new GeneralDistribution( 3, localIndexes, comm ) );
    BOOST_CHECK( ( *generaldist1 ).isEqual( *generaldist2 ) );
    BOOST_CHECK( !( *generaldist1 ).isEqual( *generaldist3 ) );
    BOOST_CHECK( !( *generaldist1 ).isEqual( *generaldist4 ) );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
