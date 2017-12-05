/**
 * @file test/matrix/SparseMatrixTest.cpp
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
 * @brief Test routines for specific methods/constructors of SparseMatrix
 * @author Thomas Brandes
 * @date 24.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/storage/Storages.hpp>
#include <scai/lama/test/storage/TestStorages.hpp>

#include <scai/lama/matrix/SparseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( SparseMatrixTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.SparseMatrixTest" );

/* ------------------------------------------------------------------------- */

/** For the matrix tests here it is sufficient to take only one of the possible value types. */

typedef RealType ValueType;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( localConstructorTest )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    TypedStorages<ValueType> storages( ctx );

    for ( size_t i = 0; i < storages.size(); ++i )
    {
        std::shared_ptr<MatrixStorage<ValueType> > storage = storages[i];

        setDenseData( *storage );  // fill it with some data

        SCAI_LOG_DEBUG( logger, *comm << ": localConstructorTest " << i << " of " << storages.size() 
                                << ", storage = " << *storage )

        dmemo::DistributionPtr rowDist( new dmemo::BlockDistribution( storage->getNumRows(), comm ) );

        SparseMatrix<ValueType> matrix( storage, rowDist );  // local constructor

        BOOST_CHECK_EQUAL( matrix.getContextPtr(), storage->getContextPtr() );

        BOOST_CHECK( matrix.isConsistent() );

        if ( rowDist->getGlobalSize() > 0 )
        {
            // be careful: all processors must call it or all not
            BOOST_CHECK( matrix.getMemoryUsage() > 0 );
        }

        BOOST_CHECK_EQUAL( matrix.getNumRows(), rowDist->getGlobalSize() );
        BOOST_CHECK_EQUAL( matrix.getNumColumns(), storage->getNumColumns() );
        BOOST_CHECK_EQUAL( matrix.getMatrixKind(), MatrixKind::SPARSE );

        // Note: storage has been modified (as shared pointer)

        BOOST_CHECK_EQUAL( storage->getNumRows(), rowDist->getLocalSize() );
    }
}

/* ------------------------------------------------------------------------- */

typedef boost::mpl::list < CSRSparseMatrix<ValueType>,
        ELLSparseMatrix<ValueType>,
        DIASparseMatrix<ValueType>,
        JDSSparseMatrix<ValueType>,
        COOSparseMatrix<ValueType>
        > SparseMatrixTypes;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SwapLocalStorageTest, _MatrixType, SparseMatrixTypes )
{
    typedef typename _MatrixType::StorageType StorageType;

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    const IndexType n = 100;  // size of the square matrix

    dmemo::DistributionPtr rowDist( new dmemo::BlockDistribution( n, comm ) );
    dmemo::DistributionPtr colDist( new dmemo::BlockDistribution( n, comm ) );

    _MatrixType matrix( rowDist, colDist );

    MatrixCreator::fillRandom( matrix, 0.2f );

    // only row block distributed

    const StorageType& haloStorage = matrix.getHaloStorage();

    StorageType storage;

    storage.allocate( matrix.getRowDistribution().getLocalSize(), matrix.getColDistribution().getLocalSize() );

    matrix.swapLocalStorage( storage );

    BOOST_CHECK( storage.getNumValues() > 0 );

    storage.scale( 2 );

    matrix.swapLocalStorage( storage );

    // halo must not have been changed

    BOOST_CHECK_EQUAL( &haloStorage, &matrix.getHaloStorage() );

    BOOST_CHECK_EQUAL( IndexType( 0 ), storage.getNumValues() );

    storage.clear(); // resizes to 0 x 0

    BOOST_CHECK_THROW(
    {
        matrix.swapLocalStorage( storage );
    }, common::Exception );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( FlatCopyTest, _MatrixType, SparseMatrixTypes )
{
    typedef typename _MatrixType::StorageType StorageType;

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    const IndexType n = 100;  // size of the square matrix

    dmemo::DistributionPtr rowDist( new dmemo::BlockDistribution( n, comm ) );
    dmemo::DistributionPtr colDist( new dmemo::BlockDistribution( n, comm ) );

    _MatrixType matrix( rowDist, colDist );

    MatrixCreator::fillRandom( matrix, 0.2f );

    // only row block distributed

    StorageType storage;
    storage.allocate( matrix.getRowDistribution().getLocalSize(), matrix.getColDistribution().getLocalSize() );

    matrix.swapLocalStorage( storage );

    BOOST_CHECK( storage.getNumValues() > 0 );

    storage.scale( 2 );

    matrix.swapLocalStorage( storage );

    BOOST_CHECK_EQUAL( IndexType( 0 ), storage.getNumValues() );

    storage.clear(); // resizes to 0 x 0

    BOOST_CHECK_THROW(
    {
        matrix.swapLocalStorage( storage );
    }, common::Exception );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( CopyConstructorTest, _MatrixType, SparseMatrixTypes )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    typedef typename _MatrixType::StorageType StorageType;

    const IndexType numRows = 4;
    const IndexType numCols = 7;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData( numRows * numCols, ValueType( 0 ) );

    common::Math::srandom( 1317 );    // makes sure that all processors generate same data

    utilskernel::HArrayUtils::setSparseRandom( denseData, fillRate, 1 );

    StorageType globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );

    dmemo::DistributionPtr rowDist( new dmemo::BlockDistribution( numRows, comm ) );
    dmemo::DistributionPtr colDist( new dmemo::BlockDistribution( numCols, comm ) );

    _MatrixType matrix1( globalStorage );
    matrix1.redistribute( rowDist, colDist );

    SCAI_LOG_INFO( logger, "Test copy constructor SparseMatrix( matrix1 ), with matrix1 = " << matrix1 )

    SparseMatrix<ValueType> matrix2( matrix1 );

    const StorageType& localStorage1 = matrix1.getLocalStorage();
    const MatrixStorage<ValueType>& localStorage2 = matrix2.getLocalStorage();

    BOOST_REQUIRE_EQUAL( localStorage1.getNumRows(), localStorage2.getNumRows() );
    BOOST_REQUIRE_EQUAL( localStorage1.getNumColumns(), localStorage2.getNumColumns() );

    BOOST_CHECK_EQUAL( 0, matrix1.maxDiffNorm( matrix2 ) );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
