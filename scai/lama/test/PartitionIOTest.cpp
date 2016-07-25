/**
 * @file PartitionIOTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Test of PartitionIO methods for write/read of partitioned data
 * @author Thomas Brandes
 * @date 16.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/TestMacros.hpp>

#include <scai/lama/io/PartitionIO.hpp>
#include <scai/lama/io/FileIO.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/utilskernel/LArray.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>

using namespace scai;
using namespace lama;
using namespace dmemo;

/** Output files should be deleted unless for debugging it might be useful to check them. */

#undef DELETE_OUTPUT_FILES

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( PartitionIOTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.PartitionIOTest" );

/* ------------------------------------------------------------------------- */

class TestDistributions : public std::vector<DistributionPtr>
{
public:

    TestDistributions()
    {
        const IndexType globalSize = 31;

        CommunicatorPtr comm = Communicator::getCommunicatorPtr();
        push_back( DistributionPtr( new CyclicDistribution( globalSize, 3, comm ) ) );
        push_back( DistributionPtr( new BlockDistribution( globalSize, comm ) ) );

        utilskernel::LArray<PartitionId> owners;

        utilskernel::HArrayUtils::setRandom( owners, globalSize, 1.0 );

        {
            hmemo::WriteAccess<PartitionId> wOwners( owners );

            PartitionId nPartitions = comm->getSize();

            for ( IndexType i = 0; i < globalSize; ++i )
            {
                wOwners[i] = wOwners[i] % nPartitions;
            }
        }

        push_back( DistributionPtr( new GeneralDistribution( owners, comm ) ) );

        float weight = static_cast<float>( comm->getRank() + 1 );

        push_back( DistributionPtr( new GenBlockDistribution( globalSize, weight, comm ) ) );
    }
};

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( DistributionSingleIO )
{
    TestDistributions testDists;

    for ( size_t i = 0; i < testDists.size(); ++i )
    {
        DistributionPtr dist = testDists[i];

        SCAI_LOG_INFO( logger, "DistributionSingleIO: dist[" << i <<  "] = " << *dist )

        CommunicatorPtr comm = dist->getCommunicatorPtr();
        const std::string fileName = "TestDist.txt";
        SCAI_LOG_INFO( logger, *comm << ": writeDistribution " << *dist )
        PartitionIO::write( *dist, fileName );
        // Hint:      We assume a common file system for all processors
        // Attention: write should have an implicit synchronization
        BOOST_CHECK( FileIO::fileExists( fileName ) );
        DistributionPtr newDist = PartitionIO::readDistribution( fileName, comm );
        SCAI_LOG_INFO( logger, *comm << ": readDistribution " << *newDist )
        // should be equal
        utilskernel::LArray<IndexType> myIndexes1;
        utilskernel::LArray<IndexType> myIndexes2;
        dist->getOwnedIndexes( myIndexes1 );
        newDist->getOwnedIndexes( myIndexes2 );
        BOOST_REQUIRE_EQUAL( myIndexes1.size(), myIndexes2.size() );
        BOOST_CHECK_EQUAL( 0, myIndexes1.maxDiffNorm( myIndexes2 ) );

#ifdef DELETE_OUTPUT_FILES
        if ( comm->getRank() == 0 )
        {
            // only one processor should delete the file
            int rc = FileIO::removeFile( fileName );
            BOOST_CHECK_EQUAL( 0, rc );
            BOOST_CHECK( !FileIO::fileExists( fileName ) );
        }
#endif

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( DistributionMultipleIO )
{
    TestDistributions testDists;

    for ( size_t i = 0; i < testDists.size(); ++i )
    {
        DistributionPtr dist = testDists[i];
        CommunicatorPtr comm = dist->getCommunicatorPtr();

        const std::string fileName  = "TestDist%r.txt";
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
        utilskernel::LArray<IndexType> owners1;
        utilskernel::LArray<IndexType> owners2;
        dist->allOwners( owners1, 0 );
        newDist->allOwners( owners2, 0 );
        BOOST_REQUIRE_EQUAL( owners1.size(), owners2.size() );
        BOOST_CHECK_EQUAL( 0, owners1.maxDiffNorm( owners2 ) );

#ifdef DELETE_OUTPUT_FILES
        int rc = FileIO::removeFile( pFileName );
        BOOST_CHECK_EQUAL( 0, rc );
        BOOST_CHECK( !FileIO::fileExists( pFileName ) );
#endif

    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( VectorSingleIO )
{
    typedef RealType ValueType;   // no focus here on type

    TestDistributions testDists;

    const std::string distFileName   = "TestDist.txt";
    const std::string vectorFileName = "TestVector.frv";

    for ( size_t i = 0; i < testDists.size(); ++i )
    {
        DistributionPtr dist = testDists[i];
        CommunicatorPtr comm = dist->getCommunicatorPtr();

        // skip this test for a replicated distribution with multiple processes
        // has a special handling

        if ( dist->isReplicated() )
        {
            CommunicatorPtr currentComm = Communicator::getCommunicatorPtr();

            if ( currentComm->getSize() > 1 )
            {
                continue;
            }
        }

        DenseVector<ValueType> vector;

        float fillRate = 1;

        vector.setRandom( dist, fillRate );

        vector.writeToFile( vectorFileName, "", common::scalar::INTERNAL, FileIO::BINARY );
        PartitionIO::write( *dist, distFileName );
        
        DenseVector<ValueType> readVector;

        // read vector and reconstruct its old distribution

        readVector.readFromFile( vectorFileName, distFileName );

        // The local parts of the two vectors must be exactly the same

        const hmemo::HArray<ValueType>& local = vector.getLocalValues();
        const hmemo::HArray<ValueType>& readLocal = readVector.getLocalValues();

        BOOST_REQUIRE_EQUAL( local.size(), readLocal.size() );

        ValueType diff = utilskernel::HArrayUtils::absMaxDiffVal( local, readLocal );

        BOOST_CHECK( diff == ValueType( 0 ) );

        int rc = PartitionIO::removeFile( vectorFileName, *comm );
        BOOST_CHECK_EQUAL( 0, rc );
        BOOST_CHECK( !PartitionIO::fileExists( vectorFileName, *comm ) );

        rc = PartitionIO::removeFile( distFileName, *comm );
        BOOST_CHECK_EQUAL( 0, rc );
        BOOST_CHECK( !PartitionIO::fileExists( distFileName, *comm ) );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( VectorPartitionIO )
{
    typedef RealType ValueType;   // no focus here on type

    TestDistributions testDists;

    const std::string distFileName   = "TestDist.txt";
    const std::string vectorFileName = "TestVector%r.frv";

    for ( size_t i = 0; i < testDists.size(); ++i )
    {
        DistributionPtr dist = testDists[i];
        CommunicatorPtr comm = dist->getCommunicatorPtr();

        // skip this test for a replicated distribution with multiple processes
        // has a special handling

        if ( dist->isReplicated() )
        {
            CommunicatorPtr currentComm = Communicator::getCommunicatorPtr();

            if ( currentComm->getSize() > 1 )
            {
                continue;
            }
        }

        DenseVector<ValueType> vector;

        float fillRate = 1;

        vector.setRandom( dist, fillRate );

        // use binary IO to avoid loss of precision

        bool withDist = dist->getBlockDistributionSize() == nIndex;

        if ( withDist )
        {
            PartitionIO::write( *dist, distFileName );
        }

        vector.writeToFile( vectorFileName, "", common::scalar::INTERNAL, FileIO::BINARY );
        
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

        BOOST_REQUIRE_EQUAL( local.size(), readLocal.size() );

        ValueType diff = utilskernel::HArrayUtils::absMaxDiffVal( local, readLocal );

        BOOST_CHECK( diff == ValueType( 0 ) );

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

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */

