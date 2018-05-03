/**
 * @file lama/test/matrix/Redistribution.cpp
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
 * @brief Redistribution of storages for matrices
 * @author Thomas Brandes
 * @date 03.08.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/Halo.hpp>
#include <scai/dmemo/HaloBuilder.hpp>
#include <scai/dmemo/Redistributor.hpp>

#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/JDSStorage.hpp>
#include <scai/lama/storage/COOStorage.hpp>
#include <scai/lama/storage/DIAStorage.hpp>
#include <scai/lama/storage/ELLStorage.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/common/test/TestMacros.hpp>
#include <scai/lama/test/storage/Storages.hpp>

#include <scai/logging.hpp>

using namespace scai;
using namespace lama;
using namespace hmemo;
using namespace dmemo;

using std::unique_ptr;

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void setDenseData( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 4;
    const IndexType numColumns = 4;
    static ValueType values[] =
    { 6.0, 0.0, 0.0, 4.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 4.0, 2.0, 5.0, 0.0, 3.0 };
    // just make sure that number of entries in values matches the matrix size
    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof( values ) / sizeof ( ValueType ) ) );
    // Note: diagonal property of sparse matrices will be set due to square matrix
    storage.setRawDenseData( numRows, numColumns, values );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( Redistribution )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.P_MatrixStorageTest" );

/* ------------------------------------------------------------------------- */

typedef boost::mpl::list < CSRStorage<float>,
        ELLStorage<double>,
        JDSStorage<float>,
        COOStorage<double>,
        DenseStorage<float>,
        DIAStorage<double>
        > StorageTypes;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildHaloTest )
{
    // it should be okay just to test one ValueType
    typedef SCAI_TEST_TYPE ValueType;
    TypedStorages<ValueType> allMatrixStorages;

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& storage = *allMatrixStorages[s];
        SCAI_LOG_DEBUG( logger, "buildHaloTest, storage = " << storage << " @ " << *storage.getContextPtr() )
        CommunicatorPtr comm = Communicator::getCommunicatorPtr();
        SCAI_LOG_INFO( logger, *comm << ": buildHaloTest" );
        setDenseData( storage );
        CSRStorage<ValueType> compare;
        setDenseData( compare );
        const IndexType numRows = storage.getNumRows();
        const IndexType numColumns = storage.getNumColumns();
        DistributionPtr colDist = DistributionPtr( new BlockDistribution( numColumns, comm ) );
        // create matrix storage for local and halo part of same type as storage
        unique_ptr<MatrixStorage<ValueType> > localStorage ( storage.newMatrixStorage() );
        unique_ptr<MatrixStorage<ValueType> > haloStorage ( storage.newMatrixStorage() );
        Halo halo;
        SCAI_LOG_INFO( logger, *comm << ", split halo : " << storage )
        storage.splitHalo( *localStorage, *haloStorage, halo, *colDist, NULL );
        SCAI_LOG_DEBUG( logger, *comm << ": split done, local = " << *localStorage
                        << ", halo = " << *haloStorage << ", halo exchg = " << halo );
        BOOST_CHECK_EQUAL( localStorage->getNumRows(), storage.getNumRows() );
        BOOST_CHECK_EQUAL( haloStorage->getNumRows(), storage.getNumRows() );
        bool keepDiagonalFlag = false;  // does not matter here
        storage.joinHalo( *localStorage, *haloStorage, halo, *colDist, keepDiagonalFlag );
        SCAI_LOG_DEBUG( logger, *comm << ": join done, result = " << storage );
        BOOST_REQUIRE_EQUAL( storage.getNumRows(), numRows );
        BOOST_REQUIRE_EQUAL( storage.getNumColumns(), numColumns );

        // Now we should have the original matrix

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType j = 0; j < numColumns; ++j )
            {
                BOOST_CHECK_EQUAL( storage.getValue( i, j ),
                                   compare.getValue( i, j ) );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( replicateTest )
{
    // it should be okay just to test one ValueType
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    typedef SCAI_TEST_TYPE ValueType;
    TypedStorages<ValueType> allMatrixStorages;

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& matrixStorage = *allMatrixStorages[s];
        SCAI_LOG_INFO( logger, *comm << ": replicateTest" );
        setDenseData( matrixStorage );
        CSRStorage<ValueType> compare;
        setDenseData( compare );
        const IndexType numRows = matrixStorage.getNumRows();
        const IndexType numColumns = matrixStorage.getNumColumns();
        BlockDistribution rowDist( numRows, comm );
        SCAI_LOG_INFO( logger, matrixStorage << ": localize"
                       " for row distribution = " << rowDist );
        // Localize the matrix data according to the source distribution
        matrixStorage.localize( matrixStorage, rowDist );
        SCAI_LOG_INFO( logger, matrixStorage << ": is local storage" );
        matrixStorage.replicate( matrixStorage, rowDist );
        SCAI_LOG_INFO( logger, matrixStorage << ": is global/replicated storage" );
        BOOST_REQUIRE_EQUAL( matrixStorage.getNumRows(), numRows );
        BOOST_REQUIRE_EQUAL( matrixStorage.getNumColumns(), numColumns );

        // Now we should have the original matrix

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType j = 0; j < numColumns; ++j )
            {
                BOOST_CHECK_EQUAL( matrixStorage.getValue( i, j ),
                                   compare.getValue( i, j ) );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( redistributeTest )
{
    // it should be okay just to test one ValueType
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    typedef SCAI_TEST_TYPE ValueType;
    TypedStorages<ValueType> allMatrixStorages;

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& matrixStorage = *allMatrixStorages[s];
        SCAI_LOG_INFO( logger, *comm << ": redistributeTest" );
        setDenseData( matrixStorage );
        CSRStorage<ValueType> compare;
        setDenseData( compare );
        const IndexType numRows = matrixStorage.getNumRows();
        const IndexType numColumns = matrixStorage.getNumColumns();
        DistributionPtr rowDist1 = DistributionPtr ( new BlockDistribution ( numRows, comm ) );
        DistributionPtr rowDist2 = DistributionPtr ( new CyclicDistribution ( numRows, 1, comm ) );
        SCAI_LOG_INFO( logger, matrixStorage << ": localize"
                       " for row distribution = " << *rowDist1 );
        // Localize the matrix data according to the source distribution
        matrixStorage.localize( matrixStorage, *rowDist1 );
        SCAI_LOG_INFO( logger, matrixStorage << ": is local storage" );
        Redistributor redistributor( rowDist2, rowDist1 );
        matrixStorage.redistribute( matrixStorage, redistributor );
        SCAI_LOG_INFO( logger, matrixStorage << ": is redistributed storage" );
        matrixStorage.replicate( matrixStorage, *rowDist2 );
        SCAI_LOG_INFO( logger, matrixStorage << ": is global/replicated storage" );
        BOOST_REQUIRE_EQUAL( matrixStorage.getNumRows(), numRows );
        BOOST_REQUIRE_EQUAL( matrixStorage.getNumColumns(), numColumns );
        // Now we should have the original matrix
        ValueType maxDiff = matrixStorage.maxDiffNorm( compare );
        BOOST_CHECK_EQUAL( 0, maxDiff );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( exchangeHaloTest )
{
    // it should be okay just to test one ValueType
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    typedef SCAI_TEST_TYPE ValueType;
    TypedStorages<ValueType> allMatrixStorages;

    for ( size_t s = 0; s < allMatrixStorages.size(); ++s )
    {
        MatrixStorage<ValueType>& matrixStorage = *allMatrixStorages[s];
        SCAI_LOG_INFO( logger, *comm << ": buildHaloTest" );
        setDenseData( matrixStorage );
        CSRStorage<ValueType> compare;
        setDenseData( compare );
        const IndexType numRows = matrixStorage.getNumRows();
        const IndexType numColumns = matrixStorage.getNumColumns();
        DistributionPtr rowDist = DistributionPtr( new BlockDistribution( numColumns, comm ) );
        // create local matrix storage
        matrixStorage.localize( matrixStorage, *rowDist );
        // build a vector of required indexes
        Halo halo;
        std::vector<IndexType> requiredIndexes;// will keep ALL non-local indexes

        for ( IndexType i = 0; i < numRows; ++i )
        {
            if ( ! rowDist->isLocal( i ) )
            {
                requiredIndexes.push_back( i );
            }
        }

        SCAI_LOG_INFO( logger, *comm << ": #required rows = " << requiredIndexes.size() );

        {
            HArrayRef<IndexType> haloIndexes( requiredIndexes );  // does not copy values
            HaloBuilder::build( *rowDist, haloIndexes, halo );
        }

        unique_ptr<MatrixStorage<ValueType> > haloMatrix( matrixStorage.newMatrixStorage() );
        haloMatrix->exchangeHalo( halo, matrixStorage, *comm );
        SCAI_LOG_INFO( logger, *comm << ": halo matrix = " << *haloMatrix );
        BOOST_REQUIRE_EQUAL( haloMatrix->getNumRows(), static_cast<IndexType>( requiredIndexes.size() ) );
        BOOST_REQUIRE_EQUAL( haloMatrix->getNumColumns(), numColumns );

        // Halo matrix must contain the right values

        for ( IndexType i = 0; i < haloMatrix->getNumRows(); ++i )
        {
            IndexType globalRow = requiredIndexes[i];

            for ( IndexType j = 0; j < numColumns; ++j )
            {
                BOOST_CHECK_EQUAL( haloMatrix->getValue( i, j ),
                                   compare.getValue( globalRow, j ) );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */
