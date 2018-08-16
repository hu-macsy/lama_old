/**
 * @file MetisPartitioningTest.cpp
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
 * @brief Specific tests for MetisPartitioning
 * @author Thomas Brandes
 * @date 30.07.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/partitioning/MetisPartitioning.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

using namespace scai;
using namespace hmemo;
using namespace dmemo;
using namespace lama;
using namespace utilskernel;
using namespace partitioning;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

BOOST_AUTO_TEST_SUITE( MetisPartitioningTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.MetisPartitioningTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    // CSR random matrix

    IndexType ia[]  = { 0,    2,    4,    6,    8,      11, 12 };
    IndexType ja[]  = { 1, 3, 0, 2, 1, 4, 0, 4, 3, 2, 5, 4 };
    DefaultReal vals[] = { 1, 3, 0, 2, 1, 4, 0, 4, 3, 2, 5, 4 };

    IndexType numRows = sizeof( ia ) / sizeof( IndexType ) - 1;
    IndexType numColumns = numRows;
    IndexType numValues = sizeof( ja ) / sizeof( IndexType );

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    if ( comm->getSize() < 2 )
    {
        SCAI_LOG_DEBUG( logger, "skipped test for MetisDistribution, only 1 processor." )
        return;
    }

    CSRStorage<DefaultReal> csrStorage;
    csrStorage.setRawCSRData( numRows, numColumns, numValues, ia, ja, vals );
    CSRSparseMatrix<DefaultReal> csrMatrix( csrStorage );

    float weight = 1.0f;

    MetisPartitioning partitioning;

    DistributionPtr dist = partitioning.partitionIt( comm, csrMatrix, weight );

    HArray<IndexType> newLocalOwners;

    partitioning.squarePartitioning( newLocalOwners, csrMatrix, weight );

    BOOST_CHECK_EQUAL( newLocalOwners.size(), csrMatrix.getRowDistribution().getLocalSize() );

    BOOST_CHECK( HArrayUtils::validIndexes( newLocalOwners, comm->getSize() ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
