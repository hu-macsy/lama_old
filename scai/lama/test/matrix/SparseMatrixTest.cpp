/**
 * @file test/matrix/SparseMatrixTest.cpp
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
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>

#include <scai/lama/storage/DenseStorage.hpp>

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

typedef DefaultReal ValueType;

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

        PartitionId np = comm->getSize();

        auto rowDist = std::make_shared<dmemo::BlockDistribution>( storage->getNumRows() * np, comm );

        SparseMatrix<ValueType> matrix( rowDist, storage );  // local constructor

        BOOST_CHECK_EQUAL( matrix.getContextPtr(), storage->getContextPtr() );

        BOOST_CHECK( matrix.isConsistent() );

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

BOOST_AUTO_TEST_CASE_TEMPLATE( CopyConstructorTest, MatrixType, SparseMatrixTypes )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    typedef typename MatrixType::StorageType StorageType;

    const IndexType numRows = 4;
    const IndexType numCols = 7;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData( numRows * numCols, ValueType( 0 ) );

    common::Math::srandom( 1317 );    // makes sure that all processors generate same data

    utilskernel::HArrayUtils::setSparseRandom( denseData, fillRate, 1 );

    // wrap the (2D) array denseData as a dense storage

    DenseStorage<ValueType> denseStorage( numRows, numCols, std::move( denseData ) );

    auto rowDist = std::make_shared<dmemo::BlockDistribution>( numRows, comm );
    auto colDist = std::make_shared<dmemo::BlockDistribution>( numCols, comm );

    auto matrix1 = distribute<MatrixType>( denseStorage, rowDist, colDist );

    SCAI_LOG_INFO( logger, "Test copy constructor SparseMatrix( matrix1 ), with matrix1 = " << matrix1 )

    SparseMatrix<ValueType> matrix2( matrix1 );

    const StorageType& localStorage1 = matrix1.getLocalStorage();
    const MatrixStorage<ValueType>& localStorage2 = matrix2.getLocalStorage();

    BOOST_REQUIRE_EQUAL( localStorage1.getNumRows(), localStorage2.getNumRows() );
    BOOST_REQUIRE_EQUAL( localStorage1.getNumColumns(), localStorage2.getNumColumns() );

    BOOST_CHECK_EQUAL( 0, matrix1.maxDiffNorm( matrix2 ) );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( matrixAddTest, MatrixType, SparseMatrixTypes )
{
    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    ValueType alpha = 2;
    ValueType beta  = -1;

    const IndexType numRows = 3;
    const IndexType numCols = 3;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData1( numRows * numCols, ValueType( 0 ) );
    hmemo::HArray<ValueType> denseData2( numRows * numCols, ValueType( 0 ) );
    hmemo::HArray<ValueType> expDenseData( numRows * numCols );

    common::Math::srandom( 1317 );    // makes sure that all processors generate same data

    utilskernel::HArrayUtils::setSparseRandom( denseData1, fillRate, 1 );
    utilskernel::HArrayUtils::setSparseRandom( denseData2, fillRate, 1 );
    utilskernel::HArrayUtils::arrayPlusArray( expDenseData, alpha, denseData1, beta, denseData2 );

    // wrap the (2D) array denseData as a dense storage

    DenseStorage<ValueType> denseStorage1( numRows, numCols, std::move( denseData1 ) );
    DenseStorage<ValueType> denseStorage2( numRows, numCols, std::move( denseData2 ) );

    DenseStorage<ValueType> expectedResult( numRows, numCols, std::move( expDenseData ) );

    auto rowDist = std::make_shared<dmemo::BlockDistribution>( numRows, comm );
    auto colDist = std::make_shared<dmemo::BlockDistribution>( numCols, comm );

    auto matrix1 = distribute<MatrixType>( denseStorage1, rowDist, colDist );
    auto matrix2 = distribute<MatrixType>( denseStorage2, rowDist, colDist );
    auto expMatrix = distribute<MatrixType>( expectedResult, rowDist, colDist );

    // parallel matrix add on row/col distributed matrix

    auto matrix  = eval<MatrixType>( alpha * matrix1 + beta * matrix2 );

    RealType<ValueType> maxDiff = matrix.maxDiffNorm( expMatrix );
    BOOST_CHECK( maxDiff < 0.0001 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
