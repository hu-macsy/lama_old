/**
 * @file PartitionIOTest.cpp
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
 * @brief Test of PartitionIO methods for write/read of partitioned data
 * @author Thomas Brandes
 * @date 16.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/lama/io/PartitionIO.hpp>
#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/utilskernel.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>

#include <scai/testsupport/uniquePathComm.hpp>
#include <scai/testsupport/GlobalTempDir.hpp>

using namespace scai;
using namespace lama;
using namespace dmemo;

using scai::testsupport::uniquePathSharedAmongNodes;
using scai::testsupport::GlobalTempDir;

using hmemo::HArray;
using utilskernel::HArrayUtils;

using boost::test_tools::per_element;

/** Output files should be deleted unless for debugging it might be useful to check them. */

#define DELETE_OUTPUT_FILES

/* ------------------------------------------------------------------------- */

/** This function removes for a parallel execution all replicated distributions.
 *
 *  Read/write of replicated data should be done more carefully.
 */

static void removeReplicatedDistributions( std::vector<dmemo::DistributionPtr>& distributions )
{
    // skip this test for a replicated distribution with multiple processes
    // has a special handling

    CommunicatorPtr currentComm = Communicator::getCommunicatorPtr();

    if ( currentComm->getSize() == 1 )
    {
        return;
    }

    IndexType pos = 0;

    for ( size_t i = 0; i < distributions.size(); ++i )
    {
        if ( distributions[i]->isReplicated() )
        {
            continue;
        }

        distributions[pos++] = distributions[i];
    }

    distributions.resize( pos );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( PartitionIOTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.PartitionIOTest" );

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( DistributionSingleIO )
{
    // Purpose: write and read of different distributions into a single file

    const IndexType n = 31;

    dmemo::TestDistributions testDists( n );
    removeReplicatedDistributions( testDists );

    for ( size_t i = 0; i < testDists.size(); ++i )
    {
        DistributionPtr dist = testDists[i];

        SCAI_LOG_INFO( logger, "DistributionSingleIO: dist[" << i <<  "] = " << *dist )

        CommunicatorPtr comm = dist->getCommunicatorPtr();
        const std::string distFileName = uniquePathSharedAmongNodes(
                                            GlobalTempDir::getPath(),
                                            *comm,
                                            "TestDist"
                                         ) + ".txt";

        BOOST_TEST_MESSAGE("DistributionSingleIO: distFilename = " << distFileName);

        SCAI_LOG_INFO( logger, *comm << ": writeDistribution " << *dist )
        PartitionIO::write( *dist, distFileName );
        // Hint:      We assume a common file system for all processors
        // Attention: write should have an implicit synchronization
        BOOST_CHECK( FileIO::fileExists( distFileName ) );
        DistributionPtr newDist = PartitionIO::readDistribution( distFileName, comm );
        SCAI_LOG_INFO( logger, *comm << ": readDistribution " << *newDist )
        // should be equal
        HArray<IndexType> myIndexes1;
        HArray<IndexType> myIndexes2;
        dist->getOwnedIndexes( myIndexes1 );
        newDist->getOwnedIndexes( myIndexes2 );

        BOOST_TEST( hostReadAccess( myIndexes1 ) == hostReadAccess( myIndexes2 ), per_element() );

#ifdef DELETE_OUTPUT_FILES
        // only one processor should delete the file
        int rc = PartitionIO::removeFile( distFileName, *comm );
        BOOST_CHECK_EQUAL( 0, rc );
        BOOST_CHECK( !PartitionIO::fileExists( distFileName, *comm ) );
#endif

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( DistributionMultipleIO )
{
    // Purpose: write and read of different distributions into a partitioned file

    const IndexType n = 25;

    dmemo::TestDistributions testDists( n );
    removeReplicatedDistributions( testDists );

    for ( size_t i = 0; i < testDists.size(); ++i )
    {
        DistributionPtr dist = testDists[i];
        CommunicatorPtr comm = dist->getCommunicatorPtr();

        const std::string fileName = uniquePathSharedAmongNodes(
                                        GlobalTempDir::getPath(),
                                        *comm,
                                        "TestDist%r"
                                        ) + ".txt";

        BOOST_TEST_MESSAGE("DistributionMultipleIO: fileName = " << fileName);

        std::string pFileName = fileName;
        bool isPartitioned;
        PartitionIO::getPartitionFileName( pFileName, isPartitioned, *comm );

        if ( !isPartitioned )
        {
            return;  // otherwise would be the same as DistributionSingleIO
        }

        SCAI_LOG_INFO( logger, *comm << ": writeDistribution " << *dist << " to " << pFileName )
        PartitionIO::write( *dist, fileName );
        BOOST_CHECK( FileIO::fileExists( pFileName ) );
        DistributionPtr newDist = PartitionIO::readDistribution( fileName, comm );
        SCAI_LOG_INFO( logger, *comm << ": readDistribution " << *newDist << " from " << pFileName )
        // collect all owners on root processor and then compare
        HArray<IndexType> owners1;
        HArray<IndexType> owners2;
        dist->allOwners( owners1, 0 );
        newDist->allOwners( owners2, 0 );

        BOOST_TEST( hostReadAccess( owners1 ) == hostReadAccess( owners2 ), per_element() );

#ifdef DELETE_OUTPUT_FILES
        int rc = PartitionIO::removeFile( fileName, *comm );
        BOOST_CHECK_EQUAL( 0, rc );
        BOOST_CHECK( !PartitionIO::fileExists( fileName, *comm ) );
#endif

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( VectorSingleIO )
{
    typedef DefaultReal ValueType;   // no focus here on type

    const IndexType n = 18;

    dmemo::TestDistributions testDists( n );
    removeReplicatedDistributions( testDists );

    const std::string distFileName   = "TestDist.txt";
    const std::string vectorFileName = "TestVector.frv";

    for ( size_t i = 0; i < testDists.size(); ++i )
    {
        DistributionPtr dist = testDists[i];
        CommunicatorPtr comm = dist->getCommunicatorPtr();

        SCAI_LOG_DEBUG( logger, "Test case " << i << " of " << testDists.size() << ": VectorSingleIO with dist = " << *dist )

        DenseVector<ValueType> vector;

        vector.setRandom( dist, 1 );

        // now write the distributed vector and its distribution, each to a single file

        vector.writeToFile( vectorFileName, "", common::ScalarType::INTERNAL, FileIO::BINARY );
        PartitionIO::write( *dist, distFileName );

        DenseVector<ValueType> readVector;

        // read vector and reconstruct its old distribution

        readVector.readFromFile( vectorFileName, distFileName );

        SCAI_LOG_DEBUG( logger, "Read vector from file " << vectorFileName << ", dist = " << distFileName 
                                << ", read vector = " << readVector );

        // The local parts of the two vectors must be exactly the same

        const hmemo::HArray<ValueType>& local = vector.getLocalValues();
        const hmemo::HArray<ValueType>& readLocal = readVector.getLocalValues();

        BOOST_REQUIRE_EQUAL( local.size(), readLocal.size() );

        ValueType diff = utilskernel::HArrayUtils::maxDiffNorm( local, readLocal );

        BOOST_CHECK( diff == ValueType( 0 ) );

        // remove the files for vector data + distribution

        int rc = PartitionIO::removeFile( vectorFileName, *comm );
        BOOST_REQUIRE_EQUAL( 0, rc );
        BOOST_CHECK( !PartitionIO::fileExists( vectorFileName, *comm ) );

        rc = PartitionIO::removeFile( distFileName, *comm );
        BOOST_REQUIRE_EQUAL( 0, rc );
        BOOST_CHECK( !PartitionIO::fileExists( distFileName, *comm ) );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( VectorPartitionIO )
{
    typedef DefaultReal ValueType;   // no focus here on type

    const IndexType n = 18;

    dmemo::TestDistributions testDists( n );
    removeReplicatedDistributions( testDists );

    const std::string distFileName   = "TestDist.txt";
    const std::string vectorFileName = "TestVector%r.frv";

    for ( size_t i = 0; i < testDists.size(); ++i )
    {
        DistributionPtr dist = testDists[i];
        CommunicatorPtr comm = dist->getCommunicatorPtr();

        DenseVector<ValueType> vector;

        vector.setRandom( dist, 1 );

        bool withDist = dist->getBlockDistributionSize() == invalidIndex;

        if ( withDist )
        {
            PartitionIO::write( *dist, distFileName );
        }

        vector.writeToFile( vectorFileName, "", common::ScalarType::INTERNAL, FileIO::BINARY );

        DenseVector<ValueType> readVector;

        if ( withDist )
        {
            readVector.readFromFile( vectorFileName, distFileName );

            SCAI_LOG_INFO( logger, "Read vector ( " << vectorFileName
                           << " ) with dist ( " << distFileName << " ): " << readVector )
        }
        else
        {
            readVector.readFromFile( vectorFileName );

            SCAI_LOG_INFO( logger, "Read vector ( " << vectorFileName << " ): " << readVector )
        }

        // we replicate now the vector, proves same distribution and same values

        DistributionPtr repDist( new NoDistribution( vector.size() ) );

        vector.redistribute( repDist );
        readVector.redistribute( repDist );

        // The two vector must be exactly the same

        const hmemo::HArray<ValueType>& local = vector.getLocalValues();
        const hmemo::HArray<ValueType>& readLocal = readVector.getLocalValues();

        BOOST_TEST( hostReadAccess( local ) == hostReadAccess( readLocal ), per_element() );

        int rc = PartitionIO::removeFile( vectorFileName, *comm );
        BOOST_CHECK_EQUAL( 0, rc );
        BOOST_CHECK( !PartitionIO::fileExists( vectorFileName, *comm ) );

        if ( withDist )
        {
            SCAI_LOG_DEBUG( logger, "remove file " << distFileName )
            int rc = PartitionIO::removeFile( distFileName, *comm );
            BOOST_CHECK_EQUAL( 0, rc );
            BOOST_CHECK( !PartitionIO::fileExists( distFileName, *comm ) );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SparseVectorPartitionIO )
{
    typedef DefaultReal ValueType;   // no focus here on type

    const IndexType n = 18;

    dmemo::TestDistributions testDists( n );
    removeReplicatedDistributions( testDists );

    const std::string distFileName   = "TestDist.txt";
    const std::string vectorFileName = "TestVector%r.frv";

    for ( size_t i = 0; i < testDists.size(); ++i )
    {
        DistributionPtr dist = testDists[i];
        CommunicatorPtr comm = dist->getCommunicatorPtr();

        SparseVector<ValueType> vector;

        float fillRate = 0.2;

        vector.setSparseRandom( dist, 0, fillRate, 1 );

        // find out it we have also to write a distribution file 

        bool withDist = dist->getBlockDistributionSize() == invalidIndex;

        if ( withDist )
        {
            PartitionIO::write( *dist, distFileName );
        }

        vector.writeToFile( vectorFileName, "", common::ScalarType::INTERNAL, FileIO::BINARY );

        SparseVector<ValueType> readVector;

        if ( withDist )
        {
            readVector.readFromFile( vectorFileName, distFileName );

            SCAI_LOG_INFO( logger, "Read vector ( " << vectorFileName
                           << " ) with dist ( " << distFileName << " ): " << readVector )
        }
        else
        {
            readVector.readFromFile( vectorFileName );

            SCAI_LOG_INFO( logger, "Read vector ( " << vectorFileName << " ): " << readVector )
        }

        // we replicate now the vector, proves same distribution and same values

        DistributionPtr repDist( new NoDistribution( vector.size() ) );

        vector.redistribute( repDist );
        readVector.redistribute( repDist );

        // The two vector must be exactly the same

        const hmemo::HArray<ValueType>& local = vector.getNonZeroValues();
        const hmemo::HArray<ValueType>& readLocal = readVector.getNonZeroValues();

        BOOST_TEST( hostReadAccess( local ) == hostReadAccess( readLocal ), per_element() );

        int rc = PartitionIO::removeFile( vectorFileName, *comm );
        BOOST_CHECK_EQUAL( 0, rc );
        BOOST_CHECK( !PartitionIO::fileExists( vectorFileName, *comm ) );

        if ( withDist )
        {
            SCAI_LOG_DEBUG( logger, "remove file " << distFileName )
            int rc = PartitionIO::removeFile( distFileName, *comm );
            BOOST_CHECK_EQUAL( 0, rc );
            BOOST_CHECK( !PartitionIO::fileExists( distFileName, *comm ) );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( _MatrixSingleIO )
{
    typedef DefaultReal ValueType;   // no focus here on type

    const IndexType numRows = 15;   // # rows for global matrix
    const IndexType numCols = 15;   // # cols for global matrix

    dmemo::TestDistributions testDists( numRows );
    removeReplicatedDistributions( testDists );

    DistributionPtr colDist( new NoDistribution( numCols ) );

    const std::string distFileName   = "TestDist.txt";
    const std::string matrixFileName = "TestMatrix.frm";

    for ( size_t i = 0; i < testDists.size(); ++i )
    {
        DistributionPtr rowDist = testDists[i];
        const Communicator& comm = rowDist->getCommunicator();

        auto matrix = zero<CSRSparseMatrix<ValueType>>( rowDist, colDist );

        float fillRate = 0.2f;

        MatrixCreator::fillRandom( matrix, fillRate );

        matrix.writeToFile( matrixFileName, "", common::ScalarType::INTERNAL, common::ScalarType::INTERNAL, FileIO::BINARY );

        PartitionIO::write( *rowDist, distFileName );

        CSRSparseMatrix<ValueType> readMatrix;

        // read matrix and reconstruct its old distribution

        readMatrix.readFromFile( matrixFileName, distFileName );

        // The local parts of the two matrices must be exactly the same

        const CSRStorage<ValueType>& local = matrix.getLocalStorage();
        const CSRStorage<ValueType>& readLocal = readMatrix.getLocalStorage();

        BOOST_REQUIRE_EQUAL( local.getNumRows(), readLocal.getNumRows() );
        BOOST_REQUIRE_EQUAL( local.getNumColumns(), readLocal.getNumColumns() );

        BOOST_CHECK_EQUAL( local.maxDiffNorm( readLocal ), 0 );

        int rc = PartitionIO::removeFile( matrixFileName, comm );
        BOOST_CHECK_EQUAL( 0, rc );
        BOOST_CHECK( !PartitionIO::fileExists( matrixFileName, comm ) );

        rc = PartitionIO::removeFile( distFileName, comm );
        BOOST_CHECK_EQUAL( 0, rc );
        BOOST_CHECK( !PartitionIO::fileExists( distFileName, comm ) );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( _MatrixPartitionIO )
{
    typedef DefaultReal ValueType;   // no focus here on type

    IndexType numRows = 15;   // # rows for global matrix
    IndexType numCols = 15;   // # cols for global matrix

    dmemo::TestDistributions testDists( numRows );
    removeReplicatedDistributions( testDists );

    DistributionPtr colDist( new NoDistribution( numCols ) );

    const std::string distFileName   = "TestDist.txt";
    const std::string matrixFileName = "TestMatrix%r.frm";

    for ( size_t i = 0; i < testDists.size(); ++i )
    {
        DistributionPtr rowDist = testDists[i];
        CommunicatorPtr comm = rowDist->getCommunicatorPtr();

        auto matrix = zero<CSRSparseMatrix<ValueType>>( rowDist, colDist );

        float fillRate = 0.2f;

        MatrixCreator::fillRandom( matrix, fillRate );

        matrix.writeToFile( matrixFileName, "", common::ScalarType::INTERNAL, common::ScalarType::INTERNAL, FileIO::BINARY );

        bool withDist = rowDist->getBlockDistributionSize() == invalidIndex;

        if ( withDist )
        {
            PartitionIO::write( *rowDist, distFileName );
        }

        SCAI_LOG_INFO( logger, "written matrix " << matrix << " to partitioned file " << distFileName )

        CSRSparseMatrix<ValueType> readMatrix;

        if ( withDist )
        {
            readMatrix.readFromFile( matrixFileName, distFileName );

            SCAI_LOG_INFO( logger, "Read matrix ( " << matrixFileName
                           << " ) with dist ( " << distFileName << " ): " << readMatrix )
        }
        else
        {
            readMatrix.readFromFile( matrixFileName );

            SCAI_LOG_INFO( logger, "Read matrix ( " << matrixFileName << " ): " << readMatrix )
        }

        BOOST_REQUIRE_EQUAL( numRows, readMatrix.getNumRows() );

        // It might happen that the read matrix has less columns if last columns have only zeros

        BOOST_REQUIRE( readMatrix.getNumColumns() <= numCols );

        // we replicate now the matrix, proves same distribution and same values

        DistributionPtr repDist( new NoDistribution( numRows ) );

        matrix.redistribute( repDist, colDist );
        readMatrix.redistribute( repDist, colDist );

        // The local parts of the two matrices must be exactly the same

        const CSRStorage<ValueType>& local = matrix.getLocalStorage();
        const CSRStorage<ValueType>& readLocal = readMatrix.getLocalStorage();

        BOOST_REQUIRE_EQUAL( local.getNumRows(), readLocal.getNumRows() );
        BOOST_REQUIRE_EQUAL( local.getNumColumns(), readLocal.getNumColumns() );

        BOOST_CHECK_EQUAL( local.maxDiffNorm( readLocal ), 0 );

        int rc = PartitionIO::removeFile( matrixFileName, *comm );
        BOOST_CHECK_EQUAL( 0, rc );
        BOOST_CHECK( !PartitionIO::fileExists( matrixFileName, *comm ) );

        if ( withDist )
        {
            rc = PartitionIO::removeFile( distFileName, *comm );
            BOOST_CHECK_EQUAL( 0, rc );
            BOOST_CHECK( !PartitionIO::fileExists( distFileName, *comm ) );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */

