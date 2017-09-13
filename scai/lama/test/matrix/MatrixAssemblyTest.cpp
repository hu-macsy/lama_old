/**
 * @file test/matrix/MatrixAssemblyTest.cpp
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
 * @brief Test routines for MatrixAssemblyAccess
 * @author Thomas Brandes
 * @date 07.09.2017
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/matrix/MatrixAssemblyAccess.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( MatrixAssemblyTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.MatrixAssemblyTest" );

/* ------------------------------------------------------------------------- */

/** For the matrix tests here it is sufficient to take only one of the possible value types. */

typedef RealType ValueType;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( simpleTest )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    PartitionId commSize = comm->getSize();
    PartitionId commRank = comm->getRank();

    IndexType numRows    = 10;
    IndexType numColumns = 15;

    IndexType raw_ia[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2 };
    IndexType raw_ja[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3 };
    ValueType raw_val[] = { 1, -1, 1, -1, 2, 1, 3, 4, 5, 8, -3, -2, -4 };

    IndexType nnz = sizeof( raw_ia ) / sizeof( IndexType );

    SCAI_ASSERT_EQ_ERROR( nnz, sizeof( raw_ja ) / sizeof( IndexType ), "" )
    SCAI_ASSERT_EQ_ERROR( nnz, sizeof( raw_val ) / sizeof( ValueType ), "" )

    dmemo::DistributionPtr rowDist( new dmemo::BlockDistribution( numRows, comm ) );
    dmemo::DistributionPtr colDist( new dmemo::NoDistribution( numColumns ) );

    CSRSparseMatrix<ValueType> matrix1( rowDist, colDist );

    // Assemble matrix data by arbitrary processors, entries at same location are ignored

    {
        MatrixAssemblyAccess<ValueType> assembly( matrix1 );

        for ( IndexType i = 0; i < nnz; ++i )
        {
            // choose two processors that push the elements

            PartitionId p1 = i % commSize;
            PartitionId p2 = ( i * 5 - 2 ) % commSize;

            if ( p1 == commRank )
            {
                assembly.push( raw_ia[i], raw_ja[i], raw_val[i] );
            }
            if ( p2 == commRank )
            {
                assembly.push( raw_ia[i], raw_ja[i], raw_val[i] );
            }
        }
    }

    SCAI_LOG_INFO( logger, *comm << ": assembled this matrix with REPLACE: " << matrix1 )

    // Assemble matrix data by arbitrary processors, entries at same location are summed up

    {
        MatrixAssemblyAccess<ValueType> assembly( matrix1, common::binary::ADD );

        for ( IndexType i = 0; i < nnz; ++i )
        {
            PartitionId p1 = ( i + 2 ) % commSize;
            PartitionId p2 = ( i * 3 + 7 ) % commSize;

            if ( p1 == commRank )
            {
                assembly.push( raw_ia[i], raw_ja[i], raw_val[i] );
            }
            if ( p2 == commRank )
            {
                assembly.push( raw_ia[i], raw_ja[i], raw_val[i] );
            }
        }
    }

    SCAI_LOG_INFO( logger, *comm << ": assembled this matrix with SUM: " << matrix1 )

    hmemo::HArrayRef<IndexType> cooIA( nnz, raw_ia );
    hmemo::HArrayRef<IndexType> cooJA( nnz, raw_ja );
    hmemo::HArrayRef<ValueType> cooValues( nnz, raw_val );

    COOStorage<ValueType> coo( numRows, numColumns, cooIA, cooJA, cooValues );

    CSRSparseMatrix<ValueType> matrix2( coo, rowDist, colDist );

    matrix2 *= 3;    // entries have been inserted once and added twice

    // verify that both matrices are same

    const CSRStorage<ValueType>& local1 = matrix1.getLocalStorage();
    const CSRStorage<ValueType>& local2 = matrix2.getLocalStorage();

    BOOST_CHECK_EQUAL( local1.maxDiffNorm( local2 ), 0 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
