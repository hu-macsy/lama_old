/**
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * @file MetisDistributionTest.cpp
 *
 * @license
 * Copyright (c) 2013
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
 * @brief MetisDistributionTest.cpp
 * @author Lauretta Schubert
 * @date 01.07.2013
 * since 1.1.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/MetisDistribution.hpp>
#include <scai/dmemo.hpp>

#include <scai/common/shared_ptr.hpp>

#include <vector>

using namespace scai::dmemo;

using scai::common::shared_ptr;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct MetisDistributionTestConfig
{
    MetisDistributionTestConfig()
    {
        comm = Communicator::getCommunicator();

        rank = comm->getRank();
        size = comm->getSize();

        globalSize = 17;

        DistributionPtr d ( new NoDistribution( globalSize ) );

        matrix.reset( new Distributed( d ) );

        // weights

        float weight = static_cast<float>( 1.0 / size );
        parts.reserve( size );

        for ( int i = 0; i < size - 1; ++i )
        {
            parts[i] = weight;
        }

        parts[ size - 1 ] = 1.0f - ( size - 1 ) * weight;
        dist = DistributionPtr( new MetisDistribution( comm, *matrix, parts ) );
    }

    ~MetisDistributionTestConfig()
    {
        comm = CommunicatorPtr();
    }

    PartitionId rank;
    PartitionId size;

    IndexType globalSize;

    DistributionPtr dist;

    shared_ptr<Distributed> matrix;

    std::vector<float> parts;
};

BOOST_FIXTURE_TEST_SUITE( MetisDistributionTest, MetisDistributionTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.MetisDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    Distributed distMatrix( *matrix );
    DistributionPtr generaldist1( new MetisDistribution( comm, distMatrix, parts ) );
    DistributionPtr generaldist2( generaldist1 );
    DistributionPtr generaldist3( new MetisDistribution( comm, distMatrix, parts ) );
    BOOST_CHECK(  ( *generaldist1 ).isEqual( *generaldist2 ) );
    BOOST_CHECK( !( *generaldist1 ).isEqual( *generaldist3 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
