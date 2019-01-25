/**
 * @file VectorConstructorTest.cpp
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
 * @brief Contains tests for all vector constructors that are same for dense
 *        and sparse vectors.
 * @author Thomas Brandes
 * @date 27.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/matrix/Matrices.hpp>
#include <scai/common/test/TestMacros.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/SparseVector.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

using namespace scai;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( VectorConstructorTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.VectorConstructorTest" )

/* --------------------------------------------------------------------- */

typedef SCAI_TEST_TYPE ValueType;

/** Define a list of all matrix types. */

typedef boost::mpl::list < DenseVector<ValueType>, SparseVector<ValueType> > VectorTypes;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SizeValueConstructorTest, VectorType, VectorTypes )
{
    // replicated vector with initial value

    IndexType n = 4;

    auto v = fill<VectorType>( n, 1 );

    BOOST_CHECK_EQUAL( n, v.size() );
    BOOST_CHECK( v.getDistribution().isReplicated() );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( v.getValue( i ), 1.0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SVConstructorTest, VectorType, VectorTypes )
{
    IndexType n = 4;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        auto v1 = fill<VectorType>( dist , 3 );
        auto v2 = eval<VectorType>( 2 * v1 );

        BOOST_CHECK_EQUAL( v2.getDistribution(), v1.getDistribution() );

        auto result = fill<VectorType>( dist , 6 );

        // prove same distribution, same values of v2 and v3

        BOOST_CHECK_EQUAL( result.maxDiffNorm( v2 ),  0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SVSVConstructorTest, VectorType, VectorTypes )
{
    IndexType n = 4;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        auto in1 = fill<VectorType>( dist , 3 );
        auto in2 = fill<VectorType>( dist , 2 );

        auto out1 = eval<VectorType>( - ( 3 * in1 + ( -in2 ) * 6 ) );

        BOOST_CHECK_EQUAL( out1.getDistribution(), in1.getDistribution() );

        auto result1 = fill<VectorType>( dist , 3 );

        // prove same distribution, same values of v2 and v3

        BOOST_CHECK_EQUAL( result1.maxDiffNorm( out1 ),  0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SVVConstructorTest, VectorType, VectorTypes )
{
    IndexType n = 4;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        auto in1 = fill<VectorType>( dist , 3 );
        auto in2 = fill<VectorType>( dist , 2 );

        auto out1 = eval<VectorType>( in1 * in2 ); 
        auto out2 = eval<VectorType>( 3 * in1 * in2 ); 
        auto out3 = eval<VectorType>( in1 * 3 * in2 ); 
        auto out4 = eval<VectorType>( in1 * in2 * 3 ); 
        auto out5 = eval<VectorType>( in1 / 3 * in2 ); 
        auto out6 = eval<VectorType>( in1 * in2 / 3 ); 
        auto out7 = eval<VectorType>( ( -in1 ) * ( -in2 ) );
        auto out8 = eval<VectorType>( 2 * in1 * in2 / 3 ); 
        auto out9 = eval<VectorType>( 2 * ( -in1 ) * ( -in2 ) / ( -1 ) ); 

        BOOST_CHECK_EQUAL( out1.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out2.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out3.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out4.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out5.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out6.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out7.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out8.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out9.getDistribution(), in1.getDistribution() );

        auto result1 = fill<VectorType>( dist , 6 );
        auto result2 = fill<VectorType>( dist , 18 );
        auto result3 = fill<VectorType>( dist , 18 );
        auto result4 = fill<VectorType>( dist , 18 );
        auto result5 = fill<VectorType>( dist , 2 );
        auto result6 = fill<VectorType>( dist , 2 );
        auto result7 = fill<VectorType>( dist , 6 );
        auto result8 = fill<VectorType>( dist , 4 );
        auto result9 = fill<VectorType>( dist , -12 );

        // prove same distribution, same values of v2 and v3

        BOOST_CHECK_EQUAL( result1.maxDiffNorm( out1 ),  0 );
        BOOST_CHECK_EQUAL( result2.maxDiffNorm( out2 ),  0 );
        BOOST_CHECK_EQUAL( result3.maxDiffNorm( out3 ),  0 );
        BOOST_CHECK_EQUAL( result4.maxDiffNorm( out4 ),  0 );
        BOOST_CHECK_EQUAL( result5.maxDiffNorm( out5 ),  0 );
        BOOST_CHECK_EQUAL( result6.maxDiffNorm( out6 ),  0 );
        BOOST_CHECK_EQUAL( result7.maxDiffNorm( out7 ),  0 );
        BOOST_CHECK_EQUAL( result8.maxDiffNorm( out8 ),  0 );
        BOOST_CHECK_EQUAL( result9.maxDiffNorm( out9 ),  0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SV_S_ConstructorTest, VectorType, VectorTypes )
{
    IndexType n = 4;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        auto in = fill<VectorType>( dist , 4 );

        auto out1 = eval<VectorType>( in + 2 );
        auto out2 = eval<VectorType>( in - 2 );
        auto out3 = eval<VectorType>( 2 + in );
        auto out4 = eval<VectorType>( 2 - in );
        auto out5 = eval<VectorType>( 2 * in + 3); 
        auto out6 = eval<VectorType>( in * 2 - 3 );
        auto out7 = eval<VectorType>( - ( 3 + in / 2 ) );
        auto out8 = eval<VectorType>( 3  - in / 2 );
        auto out9 = eval<VectorType>( 3 * ( 2 * in + 4 ) + 5 );

        // prove same distribution

        BOOST_CHECK_EQUAL( out1.getDistribution(), in.getDistribution() );
        BOOST_CHECK_EQUAL( out2.getDistribution(), in.getDistribution() );
        BOOST_CHECK_EQUAL( out3.getDistribution(), in.getDistribution() );
        BOOST_CHECK_EQUAL( out4.getDistribution(), in.getDistribution() );
        BOOST_CHECK_EQUAL( out5.getDistribution(), in.getDistribution() );
        BOOST_CHECK_EQUAL( out6.getDistribution(), in.getDistribution() );
        BOOST_CHECK_EQUAL( out7.getDistribution(), in.getDistribution() );
        BOOST_CHECK_EQUAL( out8.getDistribution(), in.getDistribution() );
        BOOST_CHECK_EQUAL( out9.getDistribution(), in.getDistribution() );

        auto result1 = fill<VectorType>( dist , 6 );
        auto result2 = fill<VectorType>( dist , 2 );
        auto result3 = fill<VectorType>( dist , 6 );
        auto result4 = fill<VectorType>( dist , -2 );
        auto result5 = fill<VectorType>( dist , 11 );
        auto result6 = fill<VectorType>( dist , 5 );
        auto result7 = fill<VectorType>( dist , -5 );
        auto result8 = fill<VectorType>( dist , 1 );
        auto result9 = fill<VectorType>( dist , 41 );

        BOOST_CHECK_EQUAL( result1.maxDiffNorm( out1 ),  0 );
        BOOST_CHECK_EQUAL( result2.maxDiffNorm( out2 ),  0 );
        BOOST_CHECK_EQUAL( result3.maxDiffNorm( out3 ),  0 );
        BOOST_CHECK_EQUAL( result4.maxDiffNorm( out4 ),  0 );
        BOOST_CHECK_EQUAL( result5.maxDiffNorm( out5 ),  0 );
        BOOST_CHECK_EQUAL( result6.maxDiffNorm( out6 ),  0 );
        BOOST_CHECK_EQUAL( result7.maxDiffNorm( out7 ),  0 );
        BOOST_CHECK_EQUAL( result8.maxDiffNorm( out8 ),  0 );
        BOOST_CHECK_EQUAL( result9.maxDiffNorm( out9 ),  0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SMV_ConstructorTest, VectorType, VectorTypes )
{
    IndexType n = 4;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        auto in = fill<VectorType>( dist , 4 );

        auto m = identity<CSRSparseMatrix<ValueType>>( dist );
        m *= 3;

        auto out1 = eval<VectorType>( m * in );
        auto out2 = eval<VectorType>( m * in * 2 );
        auto out3 = eval<VectorType>( m * 2 * in );
        auto out4 = eval<VectorType>( 2 * m * in );
        auto out5 = eval<VectorType>( - ( m * in ) );
        auto out6 = eval<VectorType>( - ( 3 * m * in ) );

        // prove same distribution

        BOOST_CHECK_EQUAL( out1.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out2.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out3.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out4.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out5.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out6.getDistribution(), *dist );

        auto result1 = fill<VectorType>( dist , 12 );
        auto result2 = fill<VectorType>( dist , 24 );
        auto result3 = fill<VectorType>( dist , 24 );
        auto result4 = fill<VectorType>( dist , 24 );
        auto result5 = fill<VectorType>( dist , -12 );
        auto result6 = fill<VectorType>( dist , -36 );

        BOOST_CHECK_EQUAL( result1.maxDiffNorm( out1 ),  0 );
        BOOST_CHECK_EQUAL( result2.maxDiffNorm( out2 ),  0 );
        BOOST_CHECK_EQUAL( result3.maxDiffNorm( out3 ),  0 );
        BOOST_CHECK_EQUAL( result4.maxDiffNorm( out4 ),  0 );
        BOOST_CHECK_EQUAL( result5.maxDiffNorm( out5 ),  0 );
        BOOST_CHECK_EQUAL( result6.maxDiffNorm( out6 ),  0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SMV_Transpose_ConstructorTest, VectorType, VectorTypes )
{
    IndexType n = 4;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        auto in = fill<VectorType>( dist , 4 );

        auto m = identity<CSRSparseMatrix<ValueType>>( dist );
        m *= 3;

        auto out1 = eval<VectorType>( transpose( m ) * in );
        auto out2 = eval<VectorType>( transpose( m ) * in * 2 );
        auto out3 = eval<VectorType>( transpose( m ) * 2 * in );
        auto out4 = eval<VectorType>( 2 * transpose( m ) * in );
        auto out5 = eval<VectorType>( - ( transpose( m ) * in ) );
        auto out6 = eval<VectorType>( - ( 3 * transpose( m ) * in ) );

        // prove same distribution

        BOOST_CHECK_EQUAL( out1.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out2.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out3.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out4.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out5.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out6.getDistribution(), *dist );

        auto result1 = fill<VectorType>( dist , 12 );
        auto result2 = fill<VectorType>( dist , 24 );
        auto result3 = fill<VectorType>( dist , 24 );
        auto result4 = fill<VectorType>( dist , 24 );
        auto result5 = fill<VectorType>( dist , -12 );
        auto result6 = fill<VectorType>( dist , -36 );

        BOOST_CHECK_EQUAL( result1.maxDiffNorm( out1 ),  0 );
        BOOST_CHECK_EQUAL( result2.maxDiffNorm( out2 ),  0 );
        BOOST_CHECK_EQUAL( result3.maxDiffNorm( out3 ),  0 );
        BOOST_CHECK_EQUAL( result4.maxDiffNorm( out4 ),  0 );
        BOOST_CHECK_EQUAL( result5.maxDiffNorm( out5 ),  0 );
        BOOST_CHECK_EQUAL( result6.maxDiffNorm( out6 ),  0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SMV_SV_ConstructorTest, VectorType, VectorTypes )
{
    IndexType n = 4;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        auto in1 = fill<VectorType>( dist , 4 );
        auto in2 = fill<VectorType>( dist , 2 );

        auto m = identity<CSRSparseMatrix<ValueType>>( dist );
        m *= 3;

        auto out1 = eval<VectorType>( m * in1 + in2 );
        auto out2 = eval<VectorType>( in2 - m * in1 );
        auto out3 = eval<VectorType>( ( -m ) * in1 - in2 );
        auto out4 = eval<VectorType>( 2 * ( m * in1 + in2 )  );
        auto out5 = eval<VectorType>( ( m * in1 + 2 * in2 ) * 3 );
        auto out6 = eval<VectorType>( ( 2 * in2  + m * in1 ) / 2 );

        // prove same distribution

        BOOST_CHECK_EQUAL( out1.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out2.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out3.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out4.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out5.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out6.getDistribution(), *dist );

        auto result1 = fill<VectorType>( dist , 14 );
        auto result2 = fill<VectorType>( dist , -10 );
        auto result3 = fill<VectorType>( dist , -14 );
        auto result4 = fill<VectorType>( dist , 28 );
        auto result5 = fill<VectorType>( dist , 48 );
        auto result6 = fill<VectorType>( dist , 8 );

        BOOST_CHECK_EQUAL( result1.maxDiffNorm( out1 ),  0 );
        BOOST_CHECK_EQUAL( result2.maxDiffNorm( out2 ),  0 );
        BOOST_CHECK_EQUAL( result3.maxDiffNorm( out3 ),  0 );
        BOOST_CHECK_EQUAL( result4.maxDiffNorm( out4 ),  0 );
        BOOST_CHECK_EQUAL( result5.maxDiffNorm( out5 ),  0 );
        BOOST_CHECK_EQUAL( result6.maxDiffNorm( out6 ),  0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SMTV_SV_ConstructorTest, VectorType, VectorTypes )
{
    IndexType n = 4;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        auto in1 = fill<VectorType>( dist , 4 );
        auto in2 = fill<VectorType>( dist , 2 );

        auto m = identity<CSRSparseMatrix<ValueType>>( dist );
        m *= 3;

        auto out1 = eval<VectorType>( transpose( m ) * in1 + in2 );
        auto out2 = eval<VectorType>( in2 - transpose( m) * in1 );
        auto out3 = eval<VectorType>( ( - transpose( m ) ) * in1 - in2 );
        auto out4 = eval<VectorType>( 2 * ( transpose( m ) * in1 + in2 )  );
        auto out5 = eval<VectorType>( ( transpose( m ) * in1 + 2 * in2 ) * 3 );
        auto out6 = eval<VectorType>( ( 2 * in2  + transpose( m ) * in1 ) / 2 );

        // prove same distribution

        BOOST_CHECK_EQUAL( out1.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out2.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out3.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out4.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out5.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out6.getDistribution(), *dist );

        auto result1 = fill<VectorType>( dist , 14 );
        auto result2 = fill<VectorType>( dist , -10 );
        auto result3 = fill<VectorType>( dist , -14 );
        auto result4 = fill<VectorType>( dist , 28 );
        auto result5 = fill<VectorType>( dist , 48 );
        auto result6 = fill<VectorType>( dist , 8 );

        BOOST_CHECK_EQUAL( result1.maxDiffNorm( out1 ),  0 );
        BOOST_CHECK_EQUAL( result2.maxDiffNorm( out2 ),  0 );
        BOOST_CHECK_EQUAL( result3.maxDiffNorm( out3 ),  0 );
        BOOST_CHECK_EQUAL( result4.maxDiffNorm( out4 ),  0 );
        BOOST_CHECK_EQUAL( result5.maxDiffNorm( out5 ),  0 );
        BOOST_CHECK_EQUAL( result6.maxDiffNorm( out6 ),  0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
