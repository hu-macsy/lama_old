/**
 * @file DenseVectorTest.cpp
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
 * @brief Contains specific tests for class DenseVector (mainly constructors)
 * @author Thomas Brandes
 * @date 27.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/TestMacros.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>

using namespace scai;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DenseVectorTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.DenseVectorTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( cTorTest, ValueType, scai_arithmetic_test_types )
{
    // replicated dense vector

    IndexType n = 4;
    DenseVector<ValueType> v( n, ValueType( 1 ) );

    BOOST_CHECK_EQUAL( n, v.size() );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( v.getValue( i ), 1.0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( matExpConstructorTest )
{
    // Note: it is sufficient to consider one matrix type and one value type
 
    typedef RealType ValueType;

    IndexType nCols = 4;
    IndexType nRows = 5;

    dmemo::TestDistributions colDists( nCols );
    dmemo::TestDistributions rowDists( nRows );

    for ( size_t i = 0; i < rowDists.size(); ++i )
    {
        dmemo::DistributionPtr rowDist = rowDists[i];

        for ( size_t j = 0; j < colDists.size(); ++j )
        {
            dmemo::DistributionPtr colDist = colDists[i];

            CSRSparseMatrix<double> mat( rowDist, colDist );
            DenseVector<ValueType> x( colDist, 3 );
            SCAI_LOG_INFO( logger, "linear algebra expression: a*A*x" );
            DenseVector<ValueType> y( 2 * mat * x );
            DenseVector<ValueType> r( rowDist, 0 );
            BOOST_CHECK( r.getLocalValues().maxDiffNorm( y.getLocalValues() ) == 0 );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( vecExpConstructorTest )
{
    typedef RealType ValueType;

    IndexType n = 4;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        DenseVector<ValueType> b( dist , 3 );
        DenseVector<ValueType> a( dist , 3 );
        DenseVector<ValueType> c( a + b );
        DenseVector<ValueType> r( dist , 6 );

        // prove same distribution, same values of r and c

        BOOST_CHECK( r.getLocalValues().maxDiffNorm( c.getLocalValues() ) == 0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
