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

BOOST_AUTO_TEST_CASE( SetGetValueTest )
{
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

    IndexType n = 5;   // let it really small, this test is very inefficient

    utilskernel::LArray<ValueType> data;

    utilskernel::HArrayUtils::setRandom( data, n, 1.0 );

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        utilskernel::LArray<ValueType> data1( n );

        dmemo::DistributionPtr dist = dists[i];
  
        DenseVector<ValueType> distV( dist );
 
        // set the value in the distributed vector

        for ( IndexType k = 0; k < n; ++k )
        {
            distV.setValue( k, ValueType( data[k] ) );
        }

        // get the value from the distributed vector

        for ( IndexType k = 0; k < n; ++k )
        {
            Scalar a ( distV.getValue( k ) );
            data1[k] = a.getValue<ValueType>();
        }

        BOOST_CHECK( data.maxDiffNorm( data1 ) == 0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetAndBuildTest )
{
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    IndexType n = 16;

    DenseVector<ValueType> repV;

    dmemo::DistributionPtr repDist( new dmemo::NoDistribution( n ) );

    // Note: all processors must have the same random numbers

    repV.setRandom( repDist, 1.0 );

    BOOST_CHECK_EQUAL( n, repV.getLocalValues().size() );

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];
  
        // distV = repV, but distributed

        DenseVector<ValueType> distV( repV, dist );
        
        DenseVector<ValueType> newV( dist );
 
        hmemo::HArray<ValueType> tmp;

        distV.buildValues( tmp );
        newV.setValues( tmp );

        // replicate newV, so we can compare with repV

        newV.redistribute( repV.getDistributionPtr() );

        BOOST_CHECK_EQUAL( n, newV.getLocalValues().size() );

        BOOST_CHECK_EQUAL( 0, newV.getLocalValues().maxDiffNorm( repV.getLocalValues() ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SequenceTest )
{
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    IndexType n = 16;

    DenseVector<ValueType> repV( n, 0, 1, ctx );  // just a sequence: 0, 1, ...

    BOOST_CHECK_EQUAL( n, repV.getLocalValues().size() );

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];
  
        // distV = repV, but distributed

        DenseVector<ValueType> distV( repV, dist );
        
        BOOST_CHECK_EQUAL( distV.min() , 0 );
        BOOST_CHECK_EQUAL( distV.max() , n - 1 );

        BOOST_CHECK_EQUAL( dist->getLocalSize(), distV.getLocalValues().size() );

        // distV1 = [ 0, ..., n-1] distributed, must be same

        DenseVector<ValueType> distV1( dist, 0, 1, ctx );

        BOOST_CHECK_EQUAL( 0, distV1.getLocalValues().maxDiffNorm( distV.getLocalValues() ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( vecAddExpConstructorTest )
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

BOOST_AUTO_TEST_CASE( vecMultExpConstructorTest )
{
    typedef RealType ValueType;

    IndexType n = 4;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        DenseVector<ValueType> b( dist , 3 );
        DenseVector<ValueType> a( dist , 3 );
        DenseVector<ValueType> c( a * b ); 
        DenseVector<ValueType> r( dist , 9 );

        // prove same distribution, same values of r and c

        BOOST_CHECK( r.getLocalValues().maxDiffNorm( c.getLocalValues() ) == 0 );
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
            SCAI_LOG_INFO( logger, "linear algebra expression: alpha * Matrix * Vector" );
            DenseVector<ValueType> y( 2 * mat * x );
            
            BOOST_CHECK_EQUAL( y.getDistribution(), *rowDist );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
