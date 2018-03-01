/**
 * @file test/matrix/GramianMatrixTest.cpp
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
 * @brief Test routines for the operator matrix class GramianMatrix
 * @author Thomas Brandes
 * @date 24.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/matrix/GramianMatrix.hpp>
#include <scai/lama/matrix/MatrixWithT.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/common/test/TestMacros.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;
using namespace dmemo;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( GramianMatrixTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.SparseMatrixTest" );

/* ------------------------------------------------------------------------- */

/** For the matrix tests here it is sufficient to take only one of the possible value types. */

typedef DefaultReal ValueType;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( constructorTest, ValueType, scai_numeric_test_types )
{
    IndexType numRows = 3;
    IndexType numCols = 4;

    // build a rectangular matrix
 
    hmemo::HArray<IndexType> ia    ( { 0,     2,     4,     6 } );
    hmemo::HArray<IndexType> ja    ( { 0, 1,  1, 2,  2,  3    } );
    hmemo::HArray<ValueType> values( { 1, 2, -1, -1, 2, 1     } );

    CSRStorage<ValueType> csr( numRows, numCols, ia, ja, values );

    TestDistributions rowDists( numRows );
    TestDistributions colDists( numCols );

    auto repRowDist = std::make_shared<NoDistribution>( numRows );
    auto repColDist = std::make_shared<NoDistribution>( numCols );

    for ( size_t ir = 0; ir < rowDists.size(); ++ir )
    {
        for ( size_t ic = 0; ic < colDists.size(); ++ic )
        {
            DistributionPtr rowDist = rowDists[ir];
            DistributionPtr colDist = colDists[ic];

            auto a = CSRSparseMatrix<ValueType>( csr );
            CSRSparseMatrix<ValueType> aT;
            aT.assignTranspose( a );

            // matrix - matrix mutliplication not well supported for distributed matrices

            auto aTaExplicit = eval<CSRSparseMatrix<ValueType>>( aT * a );

            BOOST_CHECK( aTaExplicit.checkSymmetry() );

            aTaExplicit.redistribute( colDist, colDist );
            a.redistribute( rowDist, colDist );

            GramianMatrix<ValueType> aTaImplicit( a );

            BOOST_CHECK_EQUAL( aTaImplicit.getColDistribution(), aTaExplicit.getColDistribution() );
            BOOST_CHECK_EQUAL( aTaImplicit.getRowDistribution(), aTaExplicit.getRowDistribution() );


            auto x = linearDenseVector<ValueType>( a.getColDistributionPtr(), 1.0, 0.00001 );

            auto y1 = eval<DenseVector<ValueType>>( aTaExplicit * x );
            auto y2 = eval<DenseVector<ValueType>>( aTaImplicit * x );

            RealType<ValueType> eps = 0.0001;

            BOOST_CHECK( y1.maxDiffNorm( y2 ) < eps );

            y1 = transpose( aTaExplicit ) * x;
            y2 = transpose( aTaImplicit ) * x;

            BOOST_CHECK( y1.maxDiffNorm( y2 ) < eps );
        }
    }    
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( combinedTest, ValueType, scai_numeric_test_types )
{
    IndexType numRows = 3;
    IndexType numCols = 4;

    // build a rectangular matrix
 
    hmemo::HArray<IndexType> ia    ( { 0,     2,     4,     6 } );
    hmemo::HArray<IndexType> ja    ( { 0, 1,  1, 2,  2,  3    } );
    hmemo::HArray<ValueType> values( { 1, 2, -1, -1, 2, 1     } );

    CSRSparseMatrix<ValueType> a( CSRStorage<ValueType>( numRows, numCols, ia, ja, values ) );
    CSRSparseMatrix<ValueType> aT;
    aT.assignTranspose( a );
    auto aTaExplicit = eval<CSRSparseMatrix<ValueType>>( aT * a );
    
    MatrixWithT<ValueType> optA( a );
    GramianMatrix<ValueType> aTaImplicit( optA );

    BOOST_CHECK_EQUAL( aTaImplicit.getColDistribution(), aTaExplicit.getColDistribution() );
    BOOST_CHECK_EQUAL( aTaImplicit.getRowDistribution(), aTaExplicit.getRowDistribution() );

    auto x = fill<DenseVector<ValueType>>( a.getColDistributionPtr(), 0 );
    x.fillRandom( 5 );

    auto y1 = eval<DenseVector<ValueType>>( aTaExplicit * x );
    auto y2 = eval<DenseVector<ValueType>>( aTaImplicit * x );

    RealType<ValueType> eps = 0.0001;

    y1 = transpose( aTaExplicit ) * x;
    y2 = transpose( aTaImplicit ) * x;

    BOOST_CHECK( y1.maxDiffNorm( y2 ) < eps );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
