/**
 * @file test/matrix/DenseMatrixTest.cpp
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

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DenseMatrixTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.SparseMatrixTest" );

/* ------------------------------------------------------------------------- */

/** For the matrix tests here it is sufficient to take only one of the possible value types. */

typedef RealType ValueType;

/* ------------------------------------------------------------------------- */

typedef boost::mpl::list<CSRSparseMatrix<ValueType>,
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

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    const IndexType n = 20;  // size of the square matrix

    MatrixType matrix( 20, 20 );

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

BOOST_AUTO_TEST_CASE_TEMPLATE( invertTest, ValueType, scai_arithmetic_test_types )
{
    // inverse of a matrix not yet for distributed matrices

    const IndexType n = 20;  // size of the square matrix

    DenseMatrix<ValueType> matrix( n, n );

    MatrixCreator::fillRandom( matrix, 1.0f );

    DenseMatrix<ValueType> invMatrix;

    invMatrix.invert( matrix );

    DenseMatrix<ValueType> e1;
    e1.setIdentity( n );
    DenseMatrix<ValueType> e2( matrix * invMatrix );

    BOOST_CHECK( e1.maxDiffNorm( e2 ) < 1e-3 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
