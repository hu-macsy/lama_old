/**
 * @file VectorConstructorTest.cpp
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

#include <scai/lama/expression/all.hpp>

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

    VectorType v( n, ValueType( 1 ) );

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

        VectorType v1( dist , 3 );
        VectorType v2( 2 * v1 );

        BOOST_CHECK_EQUAL( v2.getDistribution(), v1.getDistribution() );

        VectorType result( dist , 6 );

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

        VectorType in1( dist , 3 );
        VectorType in2( dist , 2 );

        VectorType out1( - ( 3 * in1 + ( -in2 ) * 6 ) );

        BOOST_CHECK_EQUAL( out1.getDistribution(), in1.getDistribution() );

        VectorType result1( dist , 3 );

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

        VectorType in1( dist , 3 );
        VectorType in2( dist , 2 );

        VectorType out1( in1 * in2 ); 
        VectorType out2( 3 * in1 * in2 ); 
        VectorType out3( in1 * 3 * in2 ); 
        VectorType out4( in1 * in2 * 3 ); 
        VectorType out5( in1 / 3 * in2 ); 
        VectorType out6( in1 * in2 / 3 ); 
        VectorType out7( ( -in1 ) * ( -in2 ) );
        VectorType out8( 2 * in1 * in2 / 3 ); 
        VectorType out9( 2 * ( -in1 ) * ( -in2 ) / ( -1 ) ); 

        BOOST_CHECK_EQUAL( out1.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out2.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out3.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out4.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out5.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out6.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out7.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out8.getDistribution(), in1.getDistribution() );
        BOOST_CHECK_EQUAL( out9.getDistribution(), in1.getDistribution() );

        VectorType result1( dist , 6 );
        VectorType result2( dist , 18 );
        VectorType result3( dist , 18 );
        VectorType result4( dist , 18 );
        VectorType result5( dist , 2 );
        VectorType result6( dist , 2 );
        VectorType result7( dist , 6 );
        VectorType result8( dist , 4 );
        VectorType result9( dist , -12 );

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

        VectorType in( dist , 4 );

        VectorType out1( in + 2 );
        VectorType out2( in - 2 );
        VectorType out3( 2 + in );
        VectorType out4( 2 - in );
        VectorType out5( 2 * in + 3); 
        VectorType out6( in * 2 - 3 );
        VectorType out7( - ( 3 + in / 2 ) );
        VectorType out8( 3  - in / 2 );
        VectorType out9( 3 * ( 2 * in + 4 ) + 5 );

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

        VectorType result1( dist , 6 );
        VectorType result2( dist , 2 );
        VectorType result3( dist , 6 );
        VectorType result4( dist , -2 );
        VectorType result5( dist , 11 );
        VectorType result6( dist , 5 );
        VectorType result7( dist , -5 );
        VectorType result8( dist , 1 );
        VectorType result9( dist , 41 );

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

        VectorType in( dist , 4 );

        CSRSparseMatrix<ValueType> m;
        m.setIdentity( dist );
        m *= 3;

        VectorType out1( m * in );
        VectorType out2( m * in * 2 );
        VectorType out3( m * 2 * in );
        VectorType out4( 2 * m * in );
        VectorType out5( - ( m * in ) );
        VectorType out6( - ( 3 * m * in ) );

        // prove same distribution

        BOOST_CHECK_EQUAL( out1.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out2.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out3.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out4.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out5.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out6.getDistribution(), *dist );

        VectorType result1( dist , 12 );
        VectorType result2( dist , 24 );
        VectorType result3( dist , 24 );
        VectorType result4( dist , 24 );
        VectorType result5( dist , -12 );
        VectorType result6( dist , -36 );

        BOOST_CHECK_EQUAL( result1.maxDiffNorm( out1 ),  0 );
        BOOST_CHECK_EQUAL( result2.maxDiffNorm( out2 ),  0 );
        BOOST_CHECK_EQUAL( result3.maxDiffNorm( out3 ),  0 );
        BOOST_CHECK_EQUAL( result4.maxDiffNorm( out4 ),  0 );
        BOOST_CHECK_EQUAL( result5.maxDiffNorm( out5 ),  0 );
        BOOST_CHECK_EQUAL( result6.maxDiffNorm( out6 ),  0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SVM_ConstructorTest, VectorType, VectorTypes )
{
    IndexType n = 4;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        VectorType in( dist , 4 );

        CSRSparseMatrix<ValueType> m;
        m.setIdentity( dist );
        m *= 3;

        VectorType out1( in * m );
        VectorType out2( in * m * 2 );
        VectorType out3( in * 2 * m );
        VectorType out4( 2 * in * m );
        VectorType out5( - ( in * m ) );
        VectorType out6( - ( 3 * in * m ) );

        // prove same distribution

        BOOST_CHECK_EQUAL( out1.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out2.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out3.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out4.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out5.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out6.getDistribution(), *dist );

        VectorType result1( dist , 12 );
        VectorType result2( dist , 24 );
        VectorType result3( dist , 24 );
        VectorType result4( dist , 24 );
        VectorType result5( dist , -12 );
        VectorType result6( dist , -36 );

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

        VectorType in1( dist , 4 );
        VectorType in2( dist , 2 );

        CSRSparseMatrix<ValueType> m;
        m.setIdentity( dist );
        m *= 3;

        VectorType out1( m * in1 + in2 );
        VectorType out2( in2 - m * in1 );
        VectorType out3( ( -m ) * in1 - in2 );
        VectorType out4( 2 * ( m * in1 + in2 )  );
        VectorType out5( ( m * in1 + 2 * in2 ) * 3 );
        VectorType out6( ( 2 * in2  + m * in1 ) / 2 );

        // prove same distribution

        BOOST_CHECK_EQUAL( out1.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out2.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out3.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out4.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out5.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out6.getDistribution(), *dist );

        VectorType result1( dist , 14 );
        VectorType result2( dist , -10 );
        VectorType result3( dist , -14 );
        VectorType result4( dist , 28 );
        VectorType result5( dist , 48 );
        VectorType result6( dist , 8 );

        BOOST_CHECK_EQUAL( result1.maxDiffNorm( out1 ),  0 );
        BOOST_CHECK_EQUAL( result2.maxDiffNorm( out2 ),  0 );
        BOOST_CHECK_EQUAL( result3.maxDiffNorm( out3 ),  0 );
        BOOST_CHECK_EQUAL( result4.maxDiffNorm( out4 ),  0 );
        BOOST_CHECK_EQUAL( result5.maxDiffNorm( out5 ),  0 );
        BOOST_CHECK_EQUAL( result6.maxDiffNorm( out6 ),  0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( SVM_SV_ConstructorTest, VectorType, VectorTypes )
{
    IndexType n = 4;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        VectorType in1( dist , 4 );
        VectorType in2( dist , 2 );

        CSRSparseMatrix<ValueType> m;
        m.setIdentity( dist );
        m *= 3;

        VectorType out1( in1 * m + in2 );
        VectorType out2( in2 - in1 * m );
        VectorType out3( in1 * ( -m ) - in2 );
        VectorType out4( 2 * ( in1 * m + in2 )  );
        VectorType out5( ( in1 * m + 2 * in2 ) * 3 );
        VectorType out6( ( 2 * in2  + in1 * m ) / 2 );

        // prove same distribution

        BOOST_CHECK_EQUAL( out1.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out2.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out3.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out4.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out5.getDistribution(), *dist );
        BOOST_CHECK_EQUAL( out6.getDistribution(), *dist );

        VectorType result1( dist , 14 );
        VectorType result2( dist , -10 );
        VectorType result3( dist , -14 );
        VectorType result4( dist , 28 );
        VectorType result5( dist , 48 );
        VectorType result6( dist , 8 );

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
