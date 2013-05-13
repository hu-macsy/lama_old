/**
 * @file P_MatrixStorageTest.cpp
 *
 * @license
 * Copyright (c) 2012
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief P_MatrixStorageTest.cpp
 * @author Thomas Brandes
 * @date 03.08.2012
 * $Id$
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/HostReadAccess.hpp>
#include <lama/distribution/Halo.hpp>
#include <lama/distribution/HaloBuilder.hpp>
#include <lama/distribution/Redistributor.hpp>
#include <lama/CommunicatorFactory.hpp>

#include <lama/storage/DenseStorage.hpp>
#include <lama/storage/CSRStorage.hpp>
#include <lama/storage/JDSStorage.hpp>
#include <lama/storage/COOStorage.hpp>
#include <lama/storage/DIAStorage.hpp>
#include <lama/storage/ELLStorage.hpp>

#include <lama/distribution/BlockDistribution.hpp>
#include <lama/distribution/CyclicDistribution.hpp>

#include <test/TestMacros.hpp>

#include <logging/logging.hpp>

using namespace boost;
using namespace lama;

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void setDenseData( MatrixStorage<ValueType>& storage )
{
    const IndexType numRows = 4;
    const IndexType numColumns = 4;

    static ValueType values[] =
    { 6.0, 0.0, 0.0, 4.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 4.0, 2.0, 5.0, 0.0, 3.0 };

    // just make sure that number of entries in values matches the matrix size

    BOOST_CHECK_EQUAL( numRows * numColumns, IndexType( sizeof(values) / sizeof (ValueType) ) );

    ValueType eps = static_cast<ValueType>( 1E-5 );

    // Note: diagonal property of sparse matrices will be set due to square matrix

    storage.setRawDenseData( numRows, numColumns, values, eps );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( P_MatrixStorageTest )
;

/* ------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( logger, "Test.P_MatrixStorageTest" );

/* ------------------------------------------------------------------------- */

typedef boost::mpl::list< CSRStorage<float>,
                          ELLStorage<double>,
                          JDSStorage<float>,
                          COOStorage<double>,
                          DenseStorage<float>,
                          DIAStorage<double>
                        > StorageTypes;

/* ------------------------------------------------------------------------- */

using boost::shared_ptr;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( buildHaloTest, StorageType, StorageTypes ) 
{
    StorageType matrixStorage;

    typedef typename StorageType::ValueType ValueType;

    CommunicatorPtr comm = CommunicatorFactory::get();

    LAMA_LOG_INFO( logger, *comm << ": buildHaloTest" );

    setDenseData( matrixStorage );

    CSRStorage<ValueType> compare;
    setDenseData( compare );

    const IndexType numRows = matrixStorage.getNumRows();
    const IndexType numColumns = matrixStorage.getNumColumns();

    DistributionPtr colDist = DistributionPtr( new BlockDistribution( numColumns, comm ) );

    // create matrix storage for local and halo part of same type as matrixStorage

    shared_ptr<MatrixStorage<ValueType> > localStorage ( matrixStorage.create() );
    shared_ptr<MatrixStorage<ValueType> > haloStorage ( matrixStorage.create() );

    Halo halo;

    LAMA_LOG_INFO( logger, *comm << ", split halo : " << matrixStorage )

    matrixStorage.splitHalo( *localStorage, *haloStorage, halo, *colDist, NULL );

    LAMA_LOG_DEBUG( logger, *comm << ": split done, local = " << *localStorage
                    << ", halo = " << *haloStorage << ", halo exchg = " << halo );

    BOOST_CHECK_EQUAL( localStorage->getNumRows(), matrixStorage.getNumRows() );
    BOOST_CHECK_EQUAL( haloStorage->getNumRows(), matrixStorage.getNumRows() );

    bool keepDiagonalFlag = false;  // does not matter here

    matrixStorage.joinHalo( *localStorage, *haloStorage, halo, *colDist, keepDiagonalFlag );

    LAMA_LOG_DEBUG( logger, *comm << ": join done, result = " << matrixStorage );

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

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( replicateTest, StorageType, StorageTypes ) {
    typedef typename StorageType::ValueType ValueType;

    StorageType matrixStorage;

    CommunicatorPtr comm = CommunicatorFactory::get();

    LAMA_LOG_INFO( logger, *comm << ": replicateTest" );

    setDenseData( matrixStorage );

    CSRStorage<ValueType> compare;
    setDenseData( compare );

    const IndexType numRows = matrixStorage.getNumRows();
    const IndexType numColumns = matrixStorage.getNumColumns();

    BlockDistribution rowDist( numRows, comm );

    LAMA_LOG_INFO( logger, matrixStorage << ": localize"
                   " for row distribution = " << rowDist );

    // Localize the matrix data according to the source distribution

    matrixStorage.localize( matrixStorage, rowDist );

    LAMA_LOG_INFO( logger, matrixStorage << ": is local storage" );

    matrixStorage.replicate( matrixStorage, rowDist );

    LAMA_LOG_INFO( logger, matrixStorage << ": is global/replicated storage" );

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

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( redistributeTest, StorageType, StorageTypes ) 
{
    typedef typename StorageType::ValueType ValueType;

    StorageType matrixStorage;

    CommunicatorPtr comm = CommunicatorFactory::get();

    LAMA_LOG_INFO( logger, *comm << ": redistributeTest" );

    setDenseData( matrixStorage );

    CSRStorage<ValueType> compare;
    setDenseData( compare );

    const IndexType numRows = matrixStorage.getNumRows();
    const IndexType numColumns = matrixStorage.getNumColumns();

    DistributionPtr rowDist1 = DistributionPtr ( new BlockDistribution ( numRows, comm ) );
    DistributionPtr rowDist2 = DistributionPtr ( new CyclicDistribution ( numRows, 1, comm ) );

    LAMA_LOG_INFO( logger, matrixStorage << ": localize"
                   " for row distribution = " << *rowDist1 );

    // Localize the matrix data according to the source distribution

    matrixStorage.localize( matrixStorage, *rowDist1 );

    LAMA_LOG_INFO( logger, matrixStorage << ": is local storage" );

    Redistributor redistributor( rowDist2, rowDist1 );

    matrixStorage.redistribute( matrixStorage, redistributor );

    LAMA_LOG_INFO( logger, matrixStorage << ": is redistributed storage" );

    matrixStorage.replicate( matrixStorage, *rowDist2 );

    LAMA_LOG_INFO( logger, matrixStorage << ": is global/replicated storage" );

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

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( exchangeHaloTest, StorageType, StorageTypes ) 
{
    // Test for method MatrixStorage::exchangeHalo

    StorageType matrixStorage;

    typedef typename StorageType::ValueType ValueType;

    CommunicatorPtr comm = CommunicatorFactory::get();

    LAMA_LOG_INFO( logger, *comm << ": buildHaloTest" );

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
        if ( ! rowDist->isLocal(i) )
        {
            requiredIndexes.push_back( i );
        }
    }

    LAMA_LOG_INFO( logger, *comm << ": #required rows = " << requiredIndexes.size() );

    HaloBuilder::build( *rowDist, requiredIndexes, halo );

    StorageType haloMatrix;
    haloMatrix.exchangeHalo( halo, matrixStorage, *comm );

    LAMA_LOG_INFO( logger, *comm << ": halo matrix = " << haloMatrix );

    BOOST_REQUIRE_EQUAL( haloMatrix.getNumRows(), static_cast<IndexType>( requiredIndexes.size() ) );
    BOOST_REQUIRE_EQUAL( haloMatrix.getNumColumns(), numColumns );

    // Halo matrix must contain the right values

    for ( IndexType i = 0; i < haloMatrix.getNumRows(); ++i )
    {
        IndexType globalRow = requiredIndexes[i];

        for ( IndexType j = 0; j < numColumns; ++j )
        {
            BOOST_CHECK_EQUAL( haloMatrix.getValue( i, j ),
                               compare.getValue( globalRow, j ) );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* ------------------------------------------------------------------------- */
