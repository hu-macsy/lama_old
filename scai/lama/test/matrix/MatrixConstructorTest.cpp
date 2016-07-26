/**
 * @file test/matrix/MatrixConstructorTest.cpp
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

#include <scai/lama/test/TestMacros.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( MatrixConstructorTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.MatrixConstructorTest" );

/* ------------------------------------------------------------------------- */

/** For the matrix tests here it is sufficient to take only one of the possible value types. */

typedef RealType ValueType;

/** Define a list of all matrix types. */

typedef boost::mpl::list<CSRSparseMatrix<ValueType>,
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
    BOOST_CHECK_EQUAL( 0, matrix.getNumRows() );
    BOOST_CHECK_EQUAL( 0, matrix.getNumColumns() );
    // check correct format / type
    BOOST_CHECK_EQUAL( common::TypeTraits<ValueType>::stype, matrix.getValueType() );
    const StorageType& local = matrix.getLocalStorage();
    BOOST_CHECK_EQUAL( 0, local.getNumRows() );
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
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( distConstructorTest, MatrixType, MatrixTypes )
{
    typedef SCAI_TEST_TYPE ValueType;    // type itself does not matter

    const IndexType numRows = 13;
    const IndexType numCols = 17;

    scai::dmemo::TestDistributions rowDists( numRows );
    scai::dmemo::TestDistributions colDists( numCols );

    typedef SCAI_TEST_TYPE ValueType;    // type itself does not matter
    
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

    hmemo::HArray<ValueType> denseData;
    std::srand( 1317 );                   // makes sure that all processors generate same data
    utilskernel::HArrayUtils::setRandom( denseData, numRows * numCols, fillRate );

    StorageType globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );

    dmemo::TestDistributions rowDists( numRows );
    dmemo::TestDistributions colDists( numCols );

    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( numCols ) );

    typedef SCAI_TEST_TYPE ValueType;    // type itself does not matter

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

BOOST_AUTO_TEST_CASE_TEMPLATE( transposeConstructorTest, MatrixType, MatrixTypes )
{
    typedef typename MatrixType::StorageType StorageType;

    const IndexType numRows = 4;
    const IndexType numCols = 7;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData;
    std::srand( 1317 );                   // makes sure that all processors generate same data
    utilskernel::HArrayUtils::setRandom( denseData, numRows * numCols, fillRate );

    StorageType globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );

    dmemo::TestDistributions rowDists( numRows );
    dmemo::TestDistributions colDists( numCols );

    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( numCols ) );

    typedef SCAI_TEST_TYPE ValueType;    // type itself does not matter

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

    hmemo::HArray<ValueType> denseData;
    std::srand( 1317 );                   // makes sure that all processors generate same data
    utilskernel::HArrayUtils::setRandom( denseData, numRows * numCols, fillRate );

    StorageType globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );

    dmemo::TestDistributions dists( numRows );

    typedef SCAI_TEST_TYPE ValueType;    // type itself does not matter

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        MatrixType matrix( globalStorage );
        matrix.redistribute( dist, dist );

        // virtual default constructor generates default matrix

        common::unique_ptr<MatrixType> newMatrix( matrix.newMatrix() );
   
        const StorageType& newLocalStorage = newMatrix->getLocalStorage();

        BOOST_CHECK_EQUAL( 0, newLocalStorage.getNumRows() );
        BOOST_CHECK_EQUAL( 0, newLocalStorage.getNumColumns() );

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

    hmemo::HArray<ValueType> denseData;
    std::srand( 1317 );                   // makes sure that all processors generate same data
    utilskernel::HArrayUtils::setRandom( denseData, numRows * numCols, fillRate );

    StorageType globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );
    dmemo::TestDistributions dists( numRows );
    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( numCols ) );

    typedef SCAI_TEST_TYPE ValueType;    // type itself does not matter

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

BOOST_AUTO_TEST_SUITE_END();
