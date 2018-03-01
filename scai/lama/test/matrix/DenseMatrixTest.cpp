/**
 * @file test/matrix/DenseMatrixTest.cpp
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
 * @brief Test routines for specific methods/constructors of DenseMatrix
 * @author Thomas Brandes
 * @date 24.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/test/TestMacros.hpp>

#include <scai/lama/test/storage/Storages.hpp>
#include <scai/lama/test/storage/TestStorages.hpp>

#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/expression/all.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DenseMatrixTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.SparseMatrixTest" );

/* ------------------------------------------------------------------------- */

/** For the matrix tests here it is sufficient to take only one of the possible value types. */

typedef DefaultReal ValueType;

/* ------------------------------------------------------------------------- */

typedef boost::mpl::list < CSRSparseMatrix<ValueType>,
        ELLSparseMatrix<ValueType>,
        JDSSparseMatrix<ValueType>,
        DIASparseMatrix<ValueType>,
        COOSparseMatrix<ValueType>,
        DenseMatrix<ValueType>
        > SparseMatrixTypes;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( matrixTimesVectorN, MatrixType, SparseMatrixTypes )
{
    // Test vector = Matrix * vector, where vector stands for multiple vectors
    // i.e. vector is a dense matrix
    // Note: not yet available for distributed matrix

    auto comm = dmemo::Communicator::getCommunicatorPtr();

    const IndexType n = 20;  // size of the square matrix

    auto matrix = zero<MatrixType>( 20, 20 );

    MatrixCreator::fillRandom( matrix, 0.2f );

    const IndexType k = 3;  // number of vectors

    DenseMatrix<ValueType> vectorK( n, k );

    MatrixCreator::fillRandom( vectorK, 1.0f );

    DenseMatrix<ValueType> vectorK1;

    vectorK1 = matrix * vectorK;  // calls matrixTimesVectorN

    BOOST_CHECK_EQUAL( n, vectorK1.getNumRows() );
    BOOST_CHECK_EQUAL( k, vectorK1.getNumColumns() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( invertTest, ValueType, scai_numeric_test_types )
{
    // inverse of a matrix not yet for distributed matrices

    const IndexType n = 20;  // size of the square matrix

    auto matrix = zero<DenseMatrix<ValueType>>( n, n );

    MatrixCreator::fillRandom( matrix, 1.0f );

    DenseMatrix<ValueType> invMatrix;

    invMatrix.invert( matrix );

    auto e1 = identity<DenseMatrix<ValueType>>( n );
    
    auto e2 = eval<DenseMatrix<ValueType>>( matrix * invMatrix );

    BOOST_CHECK( e1.maxDiffNorm( e2 ) < 1e-3 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( columnDistributionTest )
{
    // Motivation: dedicated test to check that split/join of distributed column
    //             data works which is typique for dense matrices.

    typedef DefaultReal ValueType;

    // In this method we just test that split/join of columns works fine

    const IndexType numRows = 5;
    const IndexType numColumns = 12;
 
    dmemo::DistributionPtr repRowDist( new dmemo::NoDistribution( numRows ) );
    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( numColumns ) );

    DenseStorage<ValueType> denseData( numRows, numColumns );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            denseData.setValue( i, j, ValueType( i * 100 + j ) );
        }
    }

    dmemo::TestDistributions dists( numColumns );

    for ( size_t k = 0; k < dists.size(); ++k )
    {
        dmemo::DistributionPtr dist = dists[k];

        DenseMatrix<ValueType> denseM( denseData );

        // split the column data

        denseM.redistribute( repRowDist, dist );

        // join the column data

        denseM.redistribute( repRowDist, repColDist );

        DenseStorage<ValueType> dense1 = denseM.getLocalStorage();

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType j = 0; j < numColumns; ++j )
            {
                BOOST_CHECK_EQUAL( denseData.getValue( i, j ), ValueType( i * 100 + j ) );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( dense2SparseTest )
{
    // Motivation: dedicated test to check that split/join of distributed column
    //             data works which is typique for dense matrices.

    typedef DefaultReal ValueType;

    // In this method we just test that split/join of columns works fine

    const IndexType numRows = 5;
    const IndexType numColumns = 12;
 
    dmemo::DistributionPtr repRowDist( new dmemo::NoDistribution( numRows ) );
    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( numColumns ) );

    DenseStorage<ValueType> denseData( numRows, numColumns );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            denseData.setValue( i, j, ValueType( i * 100 + j ) );
        }
    }

    dmemo::TestDistributions dists( numColumns );

    for ( size_t k = 0; k < dists.size(); ++k )
    {
        dmemo::DistributionPtr dist = dists[k];

        DenseMatrix<ValueType> denseM( denseData );

        denseM.redistribute( repRowDist, dist );

        CSRSparseMatrix<ValueType> sparseM( denseM );

        BOOST_REQUIRE_EQUAL( denseM.getColDistribution(), sparseM.getColDistribution() );
        BOOST_REQUIRE_EQUAL( denseM.getRowDistribution(), sparseM.getRowDistribution() );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType j = 0; j < numColumns; ++j )
            {
                BOOST_CHECK_EQUAL( sparseM.getValue( i, j ), ValueType( i * 100 + j ) );
            }
        }

        sparseM.redistribute( repRowDist, repColDist );

        CSRStorage<ValueType> csrData = sparseM.getLocalStorage();

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType j = 0; j < numColumns; ++j )
            {
                BOOST_CHECK_EQUAL( csrData.getValue( i, j ), ValueType( i * 100 + j ) );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
