/**
 * @file test/matrix/HybridMatrixTest.cpp
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
 * @brief Test routines for the operator matrix class GramianMatrix
 * @author Thomas Brandes
 * @date 24.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/matrix/HybridMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/common/test/TestMacros.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;
using namespace dmemo;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( HybridMatrixTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.HybridMatrixTest" );

/* ------------------------------------------------------------------------- */

/** For the matrix tests here it is sufficient to take only one of the possible value types. */

typedef DefaultReal ValueType;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( gemvTest, ValueType, scai_numeric_test_types )
{
    const IndexType numRows = 10;
    const IndexType numCols = 15;

    auto input1 = zero<CSRSparseMatrix<ValueType>>( numRows, numCols );
    auto input2 = zero<CSRSparseMatrix<ValueType>>( numRows, numCols );

    common::Math::srandom( 15 );     // all processor must generate same random numbers

    MatrixCreator::fillRandom( input1, 0.1 );
    MatrixCreator::fillRandom( input2, 0.1 );

    TestDistributions rowDists( numRows );
    TestDistributions colDists( numCols );

    for ( size_t ir = 0; ir < rowDists.size(); ++ir )
    {
        for ( size_t ic = 0; ic < colDists.size(); ++ic )
        {
            DistributionPtr rowDist = rowDists[ir];
            DistributionPtr colDist = colDists[ic];

            auto csr1 = distribute<CSRSparseMatrix<ValueType>>( input1, rowDist, colDist );
            auto csr2 = distribute<CSRSparseMatrix<ValueType>>( input2, rowDist, colDist );

            // build: csr1 + csr2, hybrid1 : explcitly, hybrid2 implicitly

            auto hybrid1 = eval<CSRSparseMatrix<ValueType>>( input1 + input2 );
            hybrid1.redistribute( rowDist, colDist );

            HybridMatrix<ValueType> hybrid2( csr1, csr2 );

            auto x = denseVectorLinear( colDist, ValueType( 1 ), ValueType( 0.2 ) );
            auto y = denseVectorLinear( rowDist, ValueType( 5 ), ValueType( -0.2 ) );

            auto y1 = denseVectorEval( 2 * hybrid1 * x  + 3 * y );
            auto y2 = denseVectorEval( 2 * hybrid2 * x  + 3 * y );

            RealType<ValueType> eps = 0.0001;

            BOOST_CHECK( y1.maxDiffNorm( y2 ) < eps );

            auto x1 = denseVectorEval( 2 * transpose( hybrid1 ) * y  + 3 * x );
            auto x2 = denseVectorEval( 2 * transpose( hybrid2 ) * y  + 3 * x );

            BOOST_CHECK( x1.maxDiffNorm( x2 ) < eps );
        }
    }    
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
