/**
 * @file MetisDistributionTest.cpp
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
 * @brief Specific tests for MetisDistribution
 * @author Thomas Brandes
 * @date 30.07.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/MetisDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/utilskernel/LArray.hpp>

using namespace scai;
using namespace dmemo;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

BOOST_AUTO_TEST_SUITE( MetisDistributionTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.MetisDistributionTest" );

/* --------------------------------------------------------------------- */

class CSRGraph : public Distributed
{
public:

    CSRGraph( DistributionPtr  dist ) : Distributed( dist )
    {
    }

    IndexType getCSRGraphSize() const
    {
        return 12;  // number of edges
    }

    virtual void buildCSRGraph( IndexType ia[], IndexType ja[], IndexType vwgt[], const IndexType* ) const
    {
        ia[0] = 0;
        ja[0] = 1;
        ja[1] = 3;
        ia[1] = 2;
        ja[2] = 0;
        ja[3] = 2;
        ia[2] = 4;
        ja[4] = 1;
        ja[5] = 4;
        ia[3] = 6;
        ja[6] = 0;
        ja[7] = 4;
        ia[4] = 8;
        ja[8] = 3;
        ja[9] = 2;
        ja[10] = 5;
        ia[5] = 11;
        ja[11] = 4;
        ia[6]  = 12;

        vwgt[0] = 2;
        vwgt[1] = 2;
        vwgt[2] = 2;
        vwgt[3] = 2;
        vwgt[4] = 3;
        vwgt[5] = 1;
    }
};

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    // CSR random Matrix

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    if ( comm->getSize() < 2 )
    {
        SCAI_LOG_DEBUG( logger, "skipped test for MetisDistribution, only 1 processor." )
        return;
    }

    DistributionPtr repD( new NoDistribution( 6 ) );

    CSRGraph csrgraph( repD );

    float weight = 1.0f;

    DistributionPtr dist( new MetisDistribution( comm, csrgraph, weight ) );

    SCAI_LOG_DEBUG( logger, *comm << ", metis dist = " << *dist )

    BOOST_CHECK( dist->getLocalSize() < dist->getGlobalSize() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
