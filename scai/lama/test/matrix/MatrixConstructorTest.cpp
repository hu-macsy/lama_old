/**
 * @file test/matrix/MatrixConstructorTest.cpp
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
 * @brief Contains constructor tests for all matrices
 * @author Thomas Brandes
 * @date 24.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/test/matrix/Matrices.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( _MatrixConstructorTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.MatrixConstructorTest" );

/* ------------------------------------------------------------------------- */

/** For the matrix tests here it is sufficient to take only one of the possible value types. */

typedef RealType ValueType;

/** Define a list of all matrix types. */

typedef boost::mpl::list < CSRSparseMatrix<ValueType>,
        DIASparseMatrix<ValueType>,
        COOSparseMatrix<ValueType>,
        JDSSparseMatrix<ValueType>,
        ELLSparseMatrix<ValueType>,
        DenseMatrix<ValueType>
        > MatrixTypes;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( defaultConstructorTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::StorageType StorageType;

    MatrixType matrix;   // default constructor

    // check zero sizes
    BOOST_CHECK_EQUAL( IndexType( 0 ), matrix.getNumRows() );
    BOOST_CHECK_EQUAL( IndexType( 0 ), matrix.getNumColumns() );

    // check correct format / type
    BOOST_CHECK_EQUAL( common::TypeTraits<ValueType>::stype, matrix.getValueType() );

    const StorageType& local = matrix.getLocalStorage();

    BOOST_CHECK_EQUAL( matrix.getFormat(), local.getFormat() );

    if ( matrix.getMatrixKind() == MatrixKind::SPARSE )
    {
        SparseMatrix<ValueType>& spMatrix = reinterpret_cast<SparseMatrix<ValueType>& >( matrix );
        spMatrix.prefetch();
    }

    BOOST_CHECK_EQUAL( IndexType( 0 ), local.getNumRows() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( sizeConstructorTest, MatrixType, MatrixTypes )
{
    const IndexType numRows = 13;
    const IndexType numCols = 17;

    MatrixType matrix( numRows, numCols );

    // check correct sizes, distributions

    BOOST_CHECK_EQUAL( numRows, matrix.getNumRows() );
    BOOST_CHECK_EQUAL( numCols, matrix.getNumColumns() );

    // local storage is full matrix, as it is replicated

    const _MatrixStorage& localStorage = matrix.getLocalStorage();

    BOOST_CHECK_EQUAL( numRows, localStorage.getNumRows() );
    BOOST_CHECK_EQUAL( numCols, localStorage.getNumColumns() );

    BOOST_CHECK( matrix.isConsistent() );

    matrix.clear();
    matrix.purge();

    BOOST_CHECK_EQUAL( IndexType( 0 ), matrix.getNumRows() );
    BOOST_CHECK_EQUAL( IndexType( 0 ), matrix.getNumColumns() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( distConstructorTest, MatrixType, MatrixTypes )
{
    const IndexType numRows = 13;
    const IndexType numCols = 17;

    scai::dmemo::TestDistributions rowDists( numRows );
    scai::dmemo::TestDistributions colDists( numCols );

    for ( size_t irow = 0; irow < rowDists.size(); ++irow )
    {
        for ( size_t icol = 0; icol < colDists.size(); ++icol )
        {
            dmemo::DistributionPtr rowDist = rowDists[irow];
            dmemo::DistributionPtr colDist = colDists[icol];

            MatrixType matrix( rowDist, colDist );

            SCAI_LOG_INFO( logger, "Matrix constructed with distributions: " << matrix )

            // check correct sizes, distributions

            BOOST_CHECK_EQUAL( numRows, matrix.getNumRows() );
            BOOST_CHECK_EQUAL( numCols, matrix.getNumColumns() );

            const _MatrixStorage& localStorage = matrix.getLocalStorage();

            BOOST_CHECK_EQUAL( rowDist->getLocalSize(), localStorage.getNumRows() );
            BOOST_CHECK_EQUAL( colDist->getLocalSize(), localStorage.getNumColumns() );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( storageConstructorTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::StorageType StorageType;

    const IndexType numRows = 4;
    const IndexType numCols = 7;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData( numRows * numCols, ValueType( 0 ) );

    common::Math::srandom( 1317 );  // makes sure that all processors generate same data

    utilskernel::HArrayUtils::setSparseRandom( denseData, fillRate, 1 );

    StorageType globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );

    dmemo::TestDistributions rowDists( numRows );
    dmemo::TestDistributions colDists( numCols );

    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( numCols ) );

    for ( size_t irow = 0; irow < rowDists.size(); ++irow )
    {
        for ( size_t icol = 0; icol < colDists.size(); ++icol )
        {
            dmemo::DistributionPtr rowDist = rowDists[irow];
            dmemo::DistributionPtr colDist = colDists[icol];

            MatrixType matrix( globalStorage );
            MatrixType matrix1( matrix, rowDist, colDist );
            MatrixType matrix2_tmp( matrix, rowDist, repColDist );
            MatrixType matrix2( matrix2_tmp.getLocalStorage(), rowDist, colDist );

            const StorageType& localStorage1 = matrix1.getLocalStorage();
            const StorageType& localStorage2 = matrix2.getLocalStorage();

            BOOST_REQUIRE_EQUAL( localStorage1.getNumRows(), localStorage2.getNumRows() );
            BOOST_REQUIRE_EQUAL( localStorage1.getNumColumns(), localStorage2.getNumColumns() );

            BOOST_CHECK_EQUAL( ValueType( 0 ), localStorage1.maxDiffNorm( localStorage2 ) );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( convertConstructorTest, MatrixType, MatrixTypes )
{
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    typedef typename MatrixType::StorageType StorageType;

    const IndexType n = 11;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData( n * n, ValueType( 0 ), ctx );

    common::Math::srandom( 1317 );  // makes sure that all processors generate same data

    utilskernel::HArrayUtils::setSparseRandom( denseData, fillRate, 1 );

    StorageType globalStorage;

    globalStorage.setDenseData( n, n, denseData );

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        _Matrices testMatrices( ctx );

        for ( size_t k = 0; k < testMatrices.size(); ++k )
        {
            _Matrix& otherMatrix = *testMatrices[k];

            otherMatrix.assign( globalStorage );

            otherMatrix.redistribute( dist, dist );

            MatrixType matrix1( otherMatrix );  // copy constructor arbitrary matrix

            MatrixType matrix2( otherMatrix, dist, dist ); // same as before but with redist

            const StorageType& localStorage1 = matrix1.getLocalStorage();
            const StorageType& localStorage2 = matrix2.getLocalStorage();

            BOOST_REQUIRE_EQUAL( localStorage1.getNumRows(), localStorage2.getNumRows() );
            BOOST_REQUIRE_EQUAL( localStorage1.getNumColumns(), localStorage2.getNumColumns() );

            BOOST_CHECK_EQUAL( ValueType( 0 ), localStorage1.maxDiffNorm( localStorage2 ) );

            // ToDo: not yet available for all types
            // MatrixType matrix3( otherMatrix, true );  // transpose
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( CopyConstructorTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::StorageType StorageType;

    const IndexType numRows = 4;
    const IndexType numCols = 7;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData( numRows * numCols, ValueType( 0 ) );
    common::Math::srandom( 1317 );  // makes sure that all processors generate same data
    utilskernel::HArrayUtils::setSparseRandom( denseData, fillRate, 1 );

    StorageType globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );

    dmemo::TestDistributions rowDists( numRows );
    dmemo::TestDistributions colDists( numCols );

    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( numCols ) );

    for ( size_t irow = 0; irow < rowDists.size(); ++irow )
    {
        for ( size_t icol = 0; icol < colDists.size(); ++icol )
        {
            dmemo::DistributionPtr rowDist = rowDists[irow];
            dmemo::DistributionPtr colDist = colDists[icol];

            MatrixType matrix1( globalStorage );
            matrix1.redistribute( rowDist, colDist );

            MatrixType matrix2( matrix1 );

            const StorageType& localStorage1 = matrix1.getLocalStorage();
            const StorageType& localStorage2 = matrix2.getLocalStorage();

            BOOST_REQUIRE_EQUAL( localStorage1.getNumRows(), localStorage2.getNumRows() );
            BOOST_REQUIRE_EQUAL( localStorage1.getNumColumns(), localStorage2.getNumColumns() );

            BOOST_CHECK_EQUAL( 0, matrix1.maxDiffNorm( matrix2 ) );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( AssignmentOpTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::StorageType StorageType;

    const IndexType numRows = 4;
    const IndexType numCols = 7;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData( numRows * numCols, ValueType( 0 ) );
    common::Math::srandom( 1317 );  // makes sure that all processors generate same data
    utilskernel::HArrayUtils::setSparseRandom( denseData, fillRate, 1 );

    StorageType globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );

    dmemo::TestDistributions rowDists( numRows );
    dmemo::TestDistributions colDists( numCols );

    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( numCols ) );

    for ( size_t irow = 0; irow < rowDists.size(); ++irow )
    {
        for ( size_t icol = 0; icol < colDists.size(); ++icol )
        {
            dmemo::DistributionPtr rowDist = rowDists[irow];
            dmemo::DistributionPtr colDist = colDists[icol];

            MatrixType matrix1( globalStorage );
            matrix1.redistribute( rowDist, colDist );

            MatrixType matrix2;

            matrix2 = matrix1;

            const StorageType& localStorage1 = matrix1.getLocalStorage();
            const StorageType& localStorage2 = matrix2.getLocalStorage();

            BOOST_REQUIRE_EQUAL( localStorage1.getNumRows(), localStorage2.getNumRows() );
            BOOST_REQUIRE_EQUAL( localStorage1.getNumColumns(), localStorage2.getNumColumns() );

            BOOST_CHECK_EQUAL( 0, matrix1.maxDiffNorm( matrix2 ) );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( transposeConstructorTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::StorageType StorageType;

    const IndexType numRows = 4;
    const IndexType numCols = 7;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData( numRows * numCols, ValueType( 0 ) );
    common::Math::srandom( 1317 );  // makes sure that all processors generate same data
    utilskernel::HArrayUtils::setSparseRandom( denseData, fillRate, 10 );

    StorageType globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );

    dmemo::TestDistributions rowDists( numRows );
    dmemo::TestDistributions colDists( numCols );

    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( numCols ) );

    for ( size_t irow = 0; irow < rowDists.size(); ++irow )
    {
        for ( size_t icol = 0; icol < colDists.size(); ++icol )
        {
            dmemo::DistributionPtr rowDist = rowDists[irow];
            dmemo::DistributionPtr colDist = colDists[icol];

            if ( rowDist->isReplicated() != colDist->isReplicated() )
            {
                continue; // transpose not supported for replicated matrices with distributed columns
            }

            MatrixType matrix1( globalStorage );
            matrix1.redistribute( rowDist, colDist );
            bool transposeFlag = true;

            SCAI_LOG_INFO( logger, "transposeConstructorTest with matrix1 = " << matrix1 )

            MatrixType matrix2( matrix1, transposeFlag );

            const StorageType& localStorage1 = matrix1.getLocalStorage();
            const StorageType& localStorage2 = matrix2.getLocalStorage();

            BOOST_REQUIRE_EQUAL( localStorage1.getNumRows(), localStorage2.getNumColumns() );
            BOOST_REQUIRE_EQUAL( localStorage1.getNumColumns(), localStorage2.getNumRows() );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( virtualConstructorTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::StorageType StorageType;

    const IndexType numRows = 4;
    const IndexType numCols = numRows;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData( numRows * numCols, ValueType( 0 ) );

    common::Math::srandom( 1317 );  // makes sure that all processors generate same data

    utilskernel::HArrayUtils::setSparseRandom( denseData, fillRate, 1 );

    StorageType globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );

    dmemo::TestDistributions dists( numRows );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        MatrixType matrix( globalStorage );
        matrix.redistribute( dist, dist );

        // virtual default constructor generates default matrix

        common::unique_ptr<MatrixType> newMatrix( matrix.newMatrix() );

        const StorageType& newLocalStorage = newMatrix->getLocalStorage();

        BOOST_CHECK_EQUAL( IndexType( 0 ), newLocalStorage.getNumRows() );
        BOOST_CHECK_EQUAL( IndexType( 0 ), newLocalStorage.getNumColumns() );

        common::unique_ptr<MatrixType > copyMatrix( matrix.copy() );

        const StorageType& localStorage = matrix.getLocalStorage();
        const StorageType& copyLocalStorage = copyMatrix->getLocalStorage();

        BOOST_CHECK_EQUAL( localStorage.getNumRows(), copyLocalStorage.getNumRows() );
        BOOST_CHECK_EQUAL( copyLocalStorage.getNumColumns(), copyLocalStorage.getNumColumns() );

        BOOST_CHECK_EQUAL( ValueType( 0 ), localStorage.maxDiffNorm( copyLocalStorage ) );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( expConstructorTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::StorageType StorageType;

    const IndexType numRows = 16;
    const IndexType numCols = numRows;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData( numRows * numCols, ValueType( 0) );

    common::Math::srandom( 1317 );   // makes sure that all processors generate same data

    utilskernel::HArrayUtils::setSparseRandom( denseData, fillRate, 1 );

    StorageType globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );
    dmemo::TestDistributions dists( numRows );
    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( numCols ) );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        MatrixType matrix1( globalStorage );

        matrix1.redistribute( dist, repColDist );     // only row distribution

        const StorageType& localStorage1 = matrix1.getLocalStorage();

        BOOST_CHECK_EQUAL( localStorage1.getNumRows(), dist->getLocalSize() );
        BOOST_CHECK_EQUAL( localStorage1.getNumColumns(), numCols );

        ValueType mult = 2;

        MatrixType matrix2( mult * matrix1 );
        MatrixType matrix3( matrix2 - matrix1 );

        // Note: matrix3 = 2 * matrix1 - matrix1 = matrix1, so matrix1 and matrix3 must be equal

        const StorageType& localStorage3 = matrix3.getLocalStorage();

        BOOST_CHECK_EQUAL( localStorage1.getNumRows(), localStorage3.getNumRows() );
        BOOST_CHECK_EQUAL( localStorage1.getNumColumns(), localStorage3.getNumColumns() );

        BOOST_CHECK_EQUAL( ValueType( 0 ), localStorage1.maxDiffNorm( localStorage3 ) );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ExpMMConstructorTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::StorageType StorageType;

    const IndexType n = 13;

    hmemo::HArray<ValueType> denseData( n * n, ValueType( 0 ) );

    common::Math::srandom( 1317 );     // makes sure that all processors generate same data

    float fillRate = 0.2;
    utilskernel::HArrayUtils::setSparseRandom( denseData, fillRate, 1 );

    StorageType globalStorage;
    globalStorage.setDenseData( n, n, denseData );
    dmemo::TestDistributions dists( n );
    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( n ) );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        MatrixType matrix1( globalStorage );

        if ( matrix1.getMatrixKind() == MatrixKind::DENSE )
        {
            return;   // DENSE supports some other distribution
        }

        MatrixType unity;
        unity.setIdentity( n );

        matrix1.redistribute( dist, repColDist );     // only row distribution

        MatrixType matrix2( matrix1 * unity );

        const StorageType& localStorage1 = matrix1.getLocalStorage();
        const StorageType& localStorage2 = matrix2.getLocalStorage();

        BOOST_CHECK_EQUAL( ValueType( 0 ), localStorage1.maxDiffNorm( localStorage2 ) );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fileConstructorTest, MatrixType, MatrixTypes )
{
    // Note: here we only test constructor MatrixType( "fileName" )
    //       as readFromFile is already tested in PartitionIO

    typedef typename MatrixType::StorageType StorageType;

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    const IndexType numRows = 16;
    const IndexType numCols = 16;

    std::string fileName = "myMatrix.psc";   // binary type, so no loss of precision

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData( numRows * numCols, ValueType( 0 ) );

    common::Math::srandom( 31991 );     // makes sure that all processors generate same data

    utilskernel::HArrayUtils::setSparseRandom( denseData, fillRate, 1 );

    StorageType globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );

    if ( comm->getRank() == 0 )
    {
        globalStorage.writeToFile( fileName );
    }

    comm->synchronize();

    dmemo::DistributionPtr rowDist( new dmemo::NoDistribution( numRows ) );
    dmemo::DistributionPtr colDist( new dmemo::NoDistribution( numCols ) );

    // comm->synchronize();

    MatrixType matrix1( fileName );

    BOOST_CHECK( matrix1.isConsistent() );
    BOOST_CHECK_EQUAL( numRows, matrix1.getNumRows() );
    BOOST_CHECK( matrix1.getNumColumns() <= numCols );

    matrix1.redistribute( rowDist, colDist );

    MatrixType matrix2( globalStorage );

    matrix2.redistribute( rowDist, colDist );

    const StorageType& localStorage1 = matrix1.getLocalStorage();
    const StorageType& localStorage2 = matrix2.getLocalStorage();

    BOOST_CHECK_EQUAL( ValueType( 0 ), localStorage1.maxDiffNorm( localStorage2 ) );

    if ( comm->getRank() == 0 )
    {
        int rc = FileIO::removeFile( fileName );
        BOOST_CHECK_EQUAL( 0, rc );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( swapTest, MatrixType, MatrixTypes )
{
    // Note: swap can only be done with matrices of same type, same format

    const IndexType n1 = 10;
    const IndexType n2 = 10;

    MatrixType matrix1;
    MatrixCreator::buildPoisson2D( matrix1, 5, n1, n1 );

    MatrixType matrix0;

    matrix0.swap( matrix1 );

    BOOST_CHECK_EQUAL( IndexType( 0 ), matrix1.getNumRows() );
    BOOST_CHECK_EQUAL( IndexType( 0 ), matrix1.getNumColumns() );

    BOOST_CHECK_EQUAL( n1 * n2, matrix0.getNumRows() );
    BOOST_CHECK_EQUAL( n1 * n2, matrix0.getNumColumns() );

    const IndexType m1 = 10;
    const IndexType m2 = 10;
    const IndexType m3 = 10;

    MatrixType matrix2;
    MatrixCreator::buildPoisson3D( matrix2, 7, m1, m2, m3 );

    SCAI_LOG_DEBUG( logger, "before swap: matrix2 = " << matrix2 << ", matrix0 = " << matrix0 )
    matrix2.swap( matrix0 );
    SCAI_LOG_DEBUG( logger, "after swap: matrix2 = " << matrix2 << ", matrix0 = " << matrix0 )

    BOOST_CHECK_EQUAL( n1 * n2, matrix2.getNumRows() );
    BOOST_CHECK_EQUAL( n1 * n2, matrix2.getNumColumns() );

    BOOST_CHECK_EQUAL( m1 * m2 * m3, matrix0.getNumRows() );
    BOOST_CHECK_EQUAL( m1 * m2 * m3, matrix0.getNumColumns() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
