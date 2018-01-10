/**
 * @file PartitioningTest.cpp
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
 * @brief Tests that will be applied to all registered partitionings of the factory.
 * @author Thomas Brandes
 * @date 23.08.2017
 */

#include <boost/test/unit_test.hpp>

#include <scai/partitioning.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>

using namespace scai;
using namespace dmemo;
using namespace partitioning;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( PartitioningTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.PartitioningTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( partitionTest )
{

    dmemo::CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    float weight = 1.0f; // same weight for all processors

    std::vector<std::string> values;  // string is create type for the factory

    Partitioning::getCreateValues( values );

    for ( size_t i = 0; i < values.size(); ++i )
    {
        PartitioningPtr part( Partitioning::create( values[i] ) );

        common::Stencil3D<DefaultReal> stencil( 27 );
        common::Grid3D grid( 20, 20, 20 );
        StencilMatrix<DefaultReal> stencilMatrix( grid, stencil );
        CSRSparseMatrix<DefaultReal> csrMatrix( stencilMatrix );

        dmemo::DistributionPtr dist( part->partitionIt( comm, csrMatrix, weight ) );

        BOOST_REQUIRE( dist.get() );
        BOOST_CHECK_EQUAL( dist->getGlobalSize(), grid.size() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */

