/**
 * @file CollectiveFileTest.cpp
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
 * @brief Test for read/write of collective files.
 * @author Thomas Brandes
 * @date 15.01.2019
 */

#include <boost/test/unit_test.hpp>

#include <scai/dmemo.hpp>
#include <scai/dmemo/CollectiveFile.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/dmemo/test/TestMacros.hpp>
#include <scai/dmemo/test/distributedArray.hpp>

using namespace scai;
using namespace dmemo;
using namespace dmemo;

using boost::test_tools::per_element;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CollectiveFileTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.CollectiveFileTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( FailTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const char fileName[] = "non_exsitingFile.data";

    SCAI_LOG_DEBUG( logger, "open non-existing file " << fileName )

    auto colFile = comm->collectiveFile();

    BOOST_CHECK_THROW(
    {
        colFile->open( fileName, "r" );
    },
    common::Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( readWriteTest )
{
    typedef DefaultReal ValueType;

    const char fileName[] = "colFile.data";

    const IndexType N = 100;

    auto fillArray = []( IndexType k ) { return ValueType( 2 * k + 1 ); };

    hmemo::HArray<IndexType> testData( { 4711, 4912 } );

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    auto dist = blockDistribution( N, comm );

    // fill a distributed array that will be written to the collective file

    hmemo::HArray<ValueType> writeData = distributedArray<ValueType>( *dist, fillArray );

    auto colFile = comm->collectiveFile();
    colFile->open( fileName, "w" );
    colFile->writeSingle( N );
    BOOST_CHECK_EQUAL( colFile->currentPos(), sizeof( N ) );
    colFile->writeAll( writeData );
    colFile->writeSingle( testData );

    // verify that all processors have the same file position
    BOOST_CHECK_EQUAL( comm->max( colFile->currentPos() ), colFile->currentPos() );

    colFile->close();

    hmemo::HArray<ValueType> readDistData;
    hmemo::HArray<IndexType> readSingleData;

    IndexType N1;

    colFile->open( fileName, "r" );
    colFile->readSingle( N1 );
    BOOST_REQUIRE_EQUAL( N, N1 );

    colFile->readAll( readDistData, dist->getLocalSize(), dist->lb() );
    colFile->readSingle( readSingleData, testData.size() );

    colFile->close();

    BOOST_TEST( hmemo::hostReadAccess( writeData ) == hmemo::hostReadAccess( readDistData ), per_element() );
    BOOST_TEST( hmemo::hostReadAccess( testData ) == hmemo::hostReadAccess( readSingleData ), per_element() );

    if ( comm->getRank() == 0 )
    {
        std::remove( fileName );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */

