/**
 * @file test/matrix/ELLSparseMatrixTest.cpp
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
 * @brief Contains only test specific for the ELL Sparse matrix
 * @author Thomas Brandes
 * @date 24.03.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/test/TestMacros.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ELLSparseMatrixTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.ELLSparseMatrixTest" );

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( defaultConstructorTest, ValueType, scai_arithmetic_test_types )
{
    ELLSparseMatrix<ValueType> matrix;
    // check zero sizes
    BOOST_CHECK_EQUAL( 0, matrix.getNumRows() );
    BOOST_CHECK_EQUAL( 0, matrix.getNumColumns() );
    // check correct format / type
    BOOST_CHECK_EQUAL( common::TypeTraits<ValueType>::stype, matrix.getValueType() );
    BOOST_CHECK_EQUAL( Matrix::ELL, matrix.getFormat() );
    const ELLStorage<ValueType>& local = matrix.getLocalStorage();
    const ELLStorage<ValueType>& halo = matrix.getHaloStorage();
    BOOST_CHECK_EQUAL( local.getNumRows(), halo.getNumRows() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( sizeConstructorTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // type itself does not matter

    const IndexType numRows = 13;
    const IndexType numCols = 17;

    ELLSparseMatrix<ValueType> matrix( numRows, numCols );

    // check correct sizes, distributions

    BOOST_CHECK_EQUAL( numRows, matrix.getNumRows() );
    BOOST_CHECK_EQUAL( numCols, matrix.getNumColumns() );

    // local storage is full matrix, as it is replicated

    const _MatrixStorage& localStorage = matrix.getLocalStorage();

    BOOST_CHECK_EQUAL( numRows, localStorage.getNumRows() );
    BOOST_CHECK_EQUAL( numCols, localStorage.getNumColumns() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( distConstructorTest )
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
         
            ELLSparseMatrix<ValueType> matrix( rowDist, colDist );
    
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

BOOST_AUTO_TEST_CASE( storageConstructorTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // type itself does not matter

    const IndexType numRows = 4;
    const IndexType numCols = 7;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData;
    std::srand( 1317 );                   // makes sure that all processors generate same data
    utilskernel::HArrayUtils::setRandom( denseData, numRows * numCols, fillRate );

    ELLStorage<ValueType> globalStorage;
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

            ELLSparseMatrix<ValueType> matrix( globalStorage );
            ELLSparseMatrix<ValueType> matrix1( matrix, rowDist, colDist );
            ELLSparseMatrix<ValueType> matrix2_tmp( matrix, rowDist, repColDist );
            ELLSparseMatrix<ValueType> matrix2( matrix2_tmp.getLocalStorage(), rowDist, colDist );

            const ELLStorage<ValueType>& localStorage1 = matrix1.getLocalStorage();
            const ELLStorage<ValueType>& localStorage2 = matrix2.getLocalStorage();

            BOOST_REQUIRE_EQUAL( localStorage1.getNumRows(), localStorage2.getNumRows() );
            BOOST_REQUIRE_EQUAL( localStorage1.getNumColumns(), localStorage2.getNumColumns() );

            BOOST_CHECK_EQUAL( ValueType( 0 ), localStorage1.maxDiffNorm( localStorage2 ) );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( transposeConstructorTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // type itself does not matter

    const IndexType numRows = 4;
    const IndexType numCols = 7;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData;
    std::srand( 1317 );                   // makes sure that all processors generate same data
    utilskernel::HArrayUtils::setRandom( denseData, numRows * numCols, fillRate );

    ELLStorage<ValueType> globalStorage;
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

            ELLSparseMatrix<ValueType> matrix1( globalStorage );
            matrix1.redistribute( rowDist, colDist );
            bool transposeFlag = true;
 
            SCAI_LOG_INFO( logger, "transposeConstructorTest with matrix1 = " << matrix1 )

            ELLSparseMatrix<ValueType> matrix2( matrix1, transposeFlag );

            const ELLStorage<ValueType>& localStorage1 = matrix1.getLocalStorage();
            const ELLStorage<ValueType>& localStorage2 = matrix2.getLocalStorage();

            BOOST_REQUIRE_EQUAL( localStorage1.getNumRows(), localStorage2.getNumColumns() );
            BOOST_REQUIRE_EQUAL( localStorage1.getNumColumns(), localStorage2.getNumRows() );
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( virtualConstructorTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // type itself does not matter

    const IndexType numRows = 4;
    const IndexType numCols = numRows;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData;
    std::srand( 1317 );                   // makes sure that all processors generate same data
    utilskernel::HArrayUtils::setRandom( denseData, numRows * numCols, fillRate );

    ELLStorage<ValueType> globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );

    dmemo::TestDistributions dists( numRows );

    typedef SCAI_TEST_TYPE ValueType;    // type itself does not matter

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        ELLSparseMatrix<ValueType> matrix( globalStorage );
        matrix.redistribute( dist, dist );

        // virtual default constructor generates default matrix

        common::unique_ptr<ELLSparseMatrix<ValueType> > newMatrix( matrix.newMatrix() );
   
        const ELLStorage<ValueType>& newLocalStorage = newMatrix->getLocalStorage();

        BOOST_CHECK_EQUAL( 0, newLocalStorage.getNumRows() );
        BOOST_CHECK_EQUAL( 0, newLocalStorage.getNumColumns() );

        common::unique_ptr<ELLSparseMatrix<ValueType> > copyMatrix( matrix.copy() );

        const ELLStorage<ValueType>& localStorage = matrix.getLocalStorage();
        const ELLStorage<ValueType>& copyLocalStorage = copyMatrix->getLocalStorage();

        BOOST_CHECK_EQUAL( localStorage.getNumRows(), copyLocalStorage.getNumRows() );
        BOOST_CHECK_EQUAL( copyLocalStorage.getNumColumns(), copyLocalStorage.getNumColumns() );

        BOOST_CHECK_EQUAL( ValueType( 0 ), localStorage.maxDiffNorm( copyLocalStorage ) );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( expConstructorTest )
{
    typedef SCAI_TEST_TYPE ValueType;    // type itself does not matter

    const IndexType numRows = 4;
    const IndexType numCols = numRows;

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData;
    std::srand( 1317 );                   // makes sure that all processors generate same data
    utilskernel::HArrayUtils::setRandom( denseData, numRows * numCols, fillRate );

    ELLStorage<ValueType> globalStorage;
    globalStorage.setDenseData( numRows, numCols, denseData );

    dmemo::TestDistributions dists( numRows );
    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( numCols ) );

    typedef SCAI_TEST_TYPE ValueType;    // type itself does not matter

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        ELLSparseMatrix<ValueType> matrix1( globalStorage );

        matrix1.redistribute( dist, repColDist );     // only row distribution

        ValueType mult = 2;

        ELLSparseMatrix<ValueType> matrix2( mult * matrix1 );
        ELLSparseMatrix<ValueType> matrix3( matrix2 - matrix1 );

        const ELLStorage<ValueType>& localStorage1 = matrix1.getLocalStorage();
        const ELLStorage<ValueType>& localStorage3 = matrix3.getLocalStorage();

        BOOST_CHECK_EQUAL( localStorage1.getNumRows(), localStorage3.getNumRows() );
        BOOST_CHECK_EQUAL( localStorage1.getNumColumns(), localStorage3.getNumColumns() );

        BOOST_CHECK_EQUAL( ValueType( 0 ), localStorage1.maxDiffNorm( localStorage3 ) );
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
