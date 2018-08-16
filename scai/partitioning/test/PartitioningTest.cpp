/**
 * @file PartitioningTest.cpp
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
 * @brief Tests that will be applied to all registered partitionings of the factory.
 * @author Thomas Brandes
 * @date 23.08.2017
 */

#include <boost/test/unit_test.hpp>

#include <scai/partitioning.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>

using namespace scai;
using namespace dmemo;
using namespace partitioning;
using namespace utilskernel;
using namespace lama;

using common::Grid3D;
using common::Stencil3D;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( PartitioningTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.PartitioningTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( serialPartitioning )
{
    // This test is only for a single node where an examples matrix is partitioned for a certain number of processors

    dmemo::CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    if ( comm->getSize() > 1 )
    {
        return;     //  this is a test for partitioning on a single node
    }

    hmemo::HArray<float> weights( { 1.0f, 1.3f, 1.5f } );

    IndexType nPart = weights.size();   // number of partitions used

    std::vector<std::string> values;  // string is create type for the factory

    Partitioning::getCreateValues( values );

    for ( size_t i = 0; i < values.size(); ++i )
    {
        PartitioningPtr part( Partitioning::create( values[i] ) );

        // generate a 27-point stencil matrix on 3D grid 20 x 20 x 20 

        StencilMatrix<DefaultReal> stencilMatrix( Grid3D( 20, 20, 20 ), Stencil3D<DefaultReal>( 27 ) );

        // generate a (replicated) CSR matrix from the stencil matrix

        CSRSparseMatrix<DefaultReal> csrMatrix( stencilMatrix );

        hmemo::HArray<PartitionId> newLocalOwners;

        SCAI_LOG_DEBUG( logger, *comm << ": partitioning of matrix " << csrMatrix << " via " << *part )

        // call the virtual partitioning method

        part->squarePartitioning( newLocalOwners, csrMatrix, weights );

        BOOST_CHECK_EQUAL( newLocalOwners.size(), csrMatrix.getRowDistribution().getLocalSize() );

        BOOST_CHECK( HArrayUtils::validIndexes( newLocalOwners, nPart ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */

