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
 * @brief Test routines for MatrixAssembly
 * @author Thomas Brandes
 * @date 07.09.2017
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama.hpp>
#include <scai/lama/matrix/MatrixAssembly.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/test/TestDistributions.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( MatrixAssemblyTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.MatrixAssemblyTest" );

/* ------------------------------------------------------------------------- */

/** For the matrix tests here it is sufficient to take only one of the possible value types. */

typedef DefaultReal ValueType;

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

    hmemo::HArrayRef<IndexType> cooIA( nnz, raw_ia );
    hmemo::HArrayRef<IndexType> cooJA( nnz, raw_ja );
    hmemo::HArrayRef<ValueType> cooValues( nnz, raw_val );

    COOStorage<ValueType> coo( numRows, numColumns, cooIA, cooJA, cooValues );

    SCAI_ASSERT_EQ_ERROR( nnz, sizeof( raw_ja ) / sizeof( IndexType ), "" )
    SCAI_ASSERT_EQ_ERROR( nnz, sizeof( raw_val ) / sizeof( ValueType ), "" )

    MatrixAssembly<ValueType> assembly;

    for ( IndexType i = 0; i < nnz; ++i )
    {
        // choose two processors that push the elements

        if ( i % commSize == commRank )
        {
            assembly.push( raw_ia[i], raw_ja[i], raw_val[i] );
        }
    }

    BOOST_CHECK_EQUAL( 10, assembly.getNumRows() );
    BOOST_CHECK_EQUAL( 10, assembly.getNumColumns() );

    dmemo::TestDistributions rowDists( numRows );

    auto colDist = std::make_shared<dmemo::NoDistribution>( numColumns );

    for ( size_t j = 0; j < rowDists.size(); ++j )
    {
        dmemo::DistributionPtr rowDist = rowDists[j];

        if ( rowDist->getCommunicator() != assembly.getCommunicator() )
        {
            continue;
        }

        auto matrixAssembled = zero<CSRSparseMatrix<ValueType>>( rowDist, colDist );

        matrixAssembled.fillFromAssembly( assembly );

        auto matrixExpected = distribute<CSRSparseMatrix<ValueType>>( coo, rowDist, colDist );

        // verify that both matrices are same

        const CSRStorage<ValueType>& local1 = matrixAssembled.getLocalStorage();
        const CSRStorage<ValueType>& local2 = matrixExpected.getLocalStorage();

        BOOST_CHECK_EQUAL( local1.maxDiffNorm( local2 ), 0 );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( globalTest )
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

    hmemo::HArrayRef<IndexType> cooIA( nnz, raw_ia );
    hmemo::HArrayRef<IndexType> cooJA( nnz, raw_ja );
    hmemo::HArrayRef<ValueType> cooValues( nnz, raw_val );

    COOStorage<ValueType> coo( numRows, numColumns, cooIA, cooJA, cooValues );

    SCAI_ASSERT_EQ_ERROR( nnz, sizeof( raw_ja ) / sizeof( IndexType ), "" )
    SCAI_ASSERT_EQ_ERROR( nnz, sizeof( raw_val ) / sizeof( ValueType ), "" )

    MatrixAssembly<ValueType> assembly;

    for ( IndexType i = 0; i < nnz; ++i )
    {
        // choose two processors that push the elements

        if ( i % commSize == commRank )
        {
            assembly.push( raw_ia[i], raw_ja[i], raw_val[i] );
        }
    }

    COOStorage<ValueType> cooGlobal = assembly.buildGlobalCOO( numRows, numColumns, common::BinaryOp::COPY );
 
    BOOST_REQUIRE_EQUAL( numRows, cooGlobal.getNumRows() );
    BOOST_REQUIRE_EQUAL( numColumns, cooGlobal.getNumColumns() );

    BOOST_CHECK_EQUAL( coo.maxDiffNorm( cooGlobal ), 0 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( failTest )
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

    MatrixAssembly<ValueType> assembly( comm );

    for ( IndexType i = 0; i < nnz; ++i )
    {
        // choose two processors that push the elements

        if ( i % commSize == commRank )
        {
            assembly.push( raw_ia[i], raw_ja[i], raw_val[i] );
        }
    }

    auto dist = std::make_shared<dmemo::BlockDistribution>( numRows, comm);

    COOStorage<ValueType> localStorage;

    localStorage = assembly.buildLocalCOO( *dist, numColumns, common::BinaryOp::COPY );

    BOOST_CHECK_EQUAL( comm->sum( localStorage.getNumValues() ), nnz );

    // number of columns too small

    BOOST_CHECK_THROW( {
        localStorage = assembly.buildLocalCOO( *dist, 5, common::BinaryOp::COPY );
    }, common::Exception );

    // number of rows too small

    dist = std::make_shared<dmemo::BlockDistribution>( 5, comm);

    BOOST_CHECK_THROW( {
        localStorage = assembly.buildLocalCOO( *dist, numColumns, common::BinaryOp::COPY );
    }, common::Exception );

    // building local data only with distribution that has same communicator

    auto noDist = std::make_shared<dmemo::NoDistribution>( numRows );

    if ( noDist->getCommunicator() != assembly.getCommunicator() )
    {
        SCAI_LOG_INFO( logger, "noDist::comm = " << noDist->getCommunicator() <<
                               " != assembly::comm = " << assembly.getCommunicator() )
        BOOST_CHECK_THROW( {
            localStorage = assembly.buildLocalCOO( *noDist, numColumns, common::BinaryOp::COPY );
        }, common::Exception );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( operatorTest )
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

    hmemo::HArrayRef<IndexType> cooIA( nnz, raw_ia );
    hmemo::HArrayRef<IndexType> cooJA( nnz, raw_ja );
    hmemo::HArrayRef<ValueType> cooValues( nnz, raw_val );

    COOStorage<ValueType> coo( numRows, numColumns, cooIA, cooJA, cooValues );

    SCAI_ASSERT_EQ_ERROR( nnz, sizeof( raw_ja ) / sizeof( IndexType ), "" )
    SCAI_ASSERT_EQ_ERROR( nnz, sizeof( raw_val ) / sizeof( ValueType ), "" )

    MatrixAssembly<ValueType> assembly;

    for ( IndexType i = 0; i < nnz; ++i )
    {
        // choose two processors that push the elements

        if ( i % commSize == commRank )
        {
            assembly.push( raw_ia[i], raw_ja[i], raw_val[i] );
        }
        if ( ( i * 5 + 2 ) % commSize == commRank )
        {
            assembly.push( raw_ia[i], raw_ja[i], raw_val[i] );
        }
    }

    BOOST_CHECK_EQUAL( 10, assembly.getNumRows() );
    BOOST_CHECK_EQUAL( 10, assembly.getNumColumns() );
    // BOOST_CHECK_EQUAL( 2 * nnz, assembly.getNumValues() );

    dmemo::TestDistributions rowDists( numRows );

    auto rowDist = std::make_shared<dmemo::BlockDistribution>( numRows );
    auto colDist = std::make_shared<dmemo::NoDistribution>( numColumns );

    auto matrixAssembled = zero<CSRSparseMatrix<ValueType>>( rowDist, colDist );

    // fillFromAssembly will redistribute the assembled data and add it to the matrix

    matrixAssembled.fillFromAssembly( assembly, common::BinaryOp::ADD );

    auto matrixExpected = distribute<CSRSparseMatrix<ValueType>>( coo, rowDist, colDist );

    matrixExpected = 2 * matrixExpected;

    // verify that both matrices are same
  
    BOOST_CHECK_EQUAL( matrixExpected.maxDiffNorm( matrixAssembled ), 0 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
