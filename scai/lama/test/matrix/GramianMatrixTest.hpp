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
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/common/test/TestMacros.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;

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
    // build a rectangular matrix
 
    hmemo::HArray<ValueType> ia ( { 0,     2,     4,     6 } );
    hmemo::HArray<ValueType> ja ( { 0, 1,  1, 2,  2,  3    } );
    hmemo::HArray<ValueType> val( { 1, 2, -1, -1, 2, 1    } );

    CSRSparseMatrix<ValueType> a( CSRStorage<ValueType>( m, n, ia, ja, values ) );
    CSRSparseMatrix<ValueType> aT;
    aT.assignTranspose( a );
    auto aTaExplicit = eval<CSRSparseMatrix<ValueType>>( aT * a );
    
    GramianMatrix<ValueType> aTaImplicit( a );

    BOOST_CHECK_EQUAL( aTaImplict.getColDistribution(), aTaExplicit.getColDistribution() );
    BOOST_CHECK_EQUAL( aTaImplict.getRowDistribution(), aTaExplicit.getRowDistribution() );

    BOOST_CHECK_EQUAL( ataImplicit.isSymmetric() );
    BOOST_CHECK_EQUAL( ataExplicit.isSymmetric() );

    auto x = fill<DenseVector<ValueType>>( a.getColDistributionPtr(), 0 );
    x.fillRandom( 5, 1.0f );

    auto y1 = eval<DenseVector<ValueType>>( aTaExplicit * x );
    auto y2 = eval<DenseVector<ValueType>>( aTaImplicit * x );

    RealType<ValueType> eps = 0.0001;

    BOOST_CHECK( y1.maxDiffNorm( y2 ) < eps );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
