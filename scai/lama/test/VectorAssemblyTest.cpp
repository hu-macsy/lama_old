/**
 * @file test/VectorAssemblyTest.cpp
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
 * @brief Test routines for VectorAssemblyAccess
 * @author Thomas Brandes
 * @date 26.09.2017
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/VectorAssemblyAccess.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/test/TestDistributions.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( VectorAssemblyTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.VectorAssemblyTest" );

/* ------------------------------------------------------------------------- */

/** For the matrix tests here it is sufficient to take only one of the possible value types. */

typedef RealType ValueType;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( vectorAssemblyTest )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    PartitionId commSize = comm->getSize();
    PartitionId commRank = comm->getRank();

    IndexType n  = 100;

    IndexType raw_ia[] =  { 0, 11, 12, 13, 34, 51, 60, 17, 18, 29, 31, 54, 73 };
    ValueType raw_val[] = { 1, -5,  1, -1,  2,  1,  3,  4,  5,  8, -3, -2, -4 };

    IndexType nnz = sizeof( raw_ia ) / sizeof( IndexType );

    SCAI_ASSERT_EQ_ERROR( nnz, sizeof( raw_val ) / sizeof( ValueType ), "" )

    dmemo::TestDistributions dists( n );

    for ( size_t j = 0; j < dists.size(); ++j )
    {

        dmemo::DistributionPtr dist = dists[j];

        SparseVector<ValueType> vector1;

        vector1.setSameValue( dist, 0 );

        // Assemble vector data by arbitrary processors, entries at same location are ignored

        {
            VectorAssemblyAccess<ValueType> assembly( vector1 );

            for ( IndexType i = 0; i < nnz; ++i )
            {
                // choose one 'arbitrary' processor that pushes the element

                PartitionId p1 = i % commSize;
                PartitionId p2 = ( i * 5 - 2 ) % commSize;

                if ( p1 == commRank )
                {
                    assembly.push( raw_ia[i], raw_val[i] );
                }

                if ( p2 == commRank )
                {
                    assembly.push( raw_ia[i], raw_val[i] );
                }
            }
        }

        SCAI_LOG_INFO( logger, *comm << ": assembled this vector with REPLACE: " << vector1 )

        // Assemble vector data by arbitrary processors, entries at same location are summed up

        {
            VectorAssemblyAccess<ValueType> assembly( vector1, common::BinaryOp::ADD );

            for ( IndexType i = 0; i < nnz; ++i )
            {
                PartitionId p1 = ( i + 2 ) % commSize;
                PartitionId p2 = ( i * 3 + 7 ) % commSize;

                if ( p1 == commRank )
                {
                    assembly.push( raw_ia[i], raw_val[i] );
                }

                if ( p2 == commRank )
                {
                    assembly.push( raw_ia[i], raw_val[i] );
                }

                // assembly called by each processor, but only inserted once by owner

                assembly.pushReplicated( raw_ia[i], raw_val[i] );
            }
        }

        SCAI_LOG_INFO( logger, *comm << ": assembled this vector with SUM: " << vector1 )

        SparseVector<ValueType> vector2;

        vector2.setSparseRawData( n, nnz, raw_ia, raw_val, 0 );

        vector2.redistribute( dist );

        vector2 *= 4;

        // both vectors must be exaclty the same

        BOOST_CHECK_EQUAL( vector1.maxDiffNorm( vector2 ), 0 );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
