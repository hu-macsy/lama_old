/**
 * @file SparseVectorTest.cpp
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
 * @brief Contains specific tests for class SparseVector (mainly constructors)
 * @author Thomas Brandes
 * @date 17.01.2017
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/test/TestMacros.hpp>
#include <scai/lama/test/matrix/Matrices.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/lama/SparseVector.hpp>
#include <scai/lama/matrix/SparseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/lama/expression/all.hpp>

using namespace scai;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( SparseVectorTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.SparseVectorTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( cTorTest, ValueType, scai_numeric_test_types )
{
    // replicated sparse vector of a certain size only with zero elements

    ValueType zeroValues[] = { 0, 1 };

    const IndexType n = 4;

    for ( IndexType icase = 0; icase < 2; ++icase )
    {
        ValueType zero = zeroValues[icase];

        SparseVector<ValueType> v( n, zero );

        BOOST_CHECK_EQUAL( n, v.size() );

        for ( IndexType i = 0; i < n; ++i )
        {
            BOOST_CHECK_EQUAL( v.getValue( i ), zero );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( consistencyTest )
{
    typedef RealType ValueType;

    IndexType n = 10;

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    // create distributed sparse vector, initialized with 0

    SparseVector<ValueType> v( dist );

    BOOST_CHECK( v.isConsistent() );

    // usually it is not simple to make a vector inconsistent

    hmemo::HArray<IndexType>& indexes = const_cast<hmemo::HArray<IndexType>&>( v.getNonZeroIndexes() );

    // make it inconsistent on one processor only, isConsistent can deal with it

    if ( comm->getRank() == ( comm->getSize() / 2 ) )
    {
        indexes.resize( 1 );
    }

    // consistency check fails on all processors

    BOOST_CHECK( ! v.isConsistent() );

    hmemo::HArray<ValueType>& values = const_cast<hmemo::HArray<ValueType>&>( v.getNonZeroValues() );

    if ( dist->getLocalSize() > 0 )
    {
        indexes.setSameValue( 1, IndexType( 0 ) );
        values.setSameValue( 1, ValueType( 3 ) );
    }
    else
    {
        indexes.clear();
        values.clear();
    }

    BOOST_CHECK( v.isConsistent() );

    if ( dist->getLocalSize() > 0 )
    {
        indexes.setSameValue( dist->getLocalSize(), IndexType( 1 ) );
    }

    BOOST_CHECK( ! v.isConsistent() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CopyConstructorTest )
{
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

    IndexType n = 20;   // let it really small, this test is very inefficient
    IndexType bound = 1000;

    utilskernel::LArray<ValueType> data( n, ValueType( 0 ) );

    common::Math::srandom( 13151 );   // This test only works if all processors have same random numbers

    const float fillRate = 0.2;   // makes it worth to be sparse

    data.setSparseRandom( fillRate, bound ); // ~ 20% values are non-zero with value between 0 and 1

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        utilskernel::LArray<ValueType> data1( n );

        dmemo::DistributionPtr dist = dists[i];

        DenseVector<ValueType> denseV( data );

        denseV.redistribute( dist );

        SparseVector<ValueType> sparseV( denseV );

        // get each value from the distributed sparse vector

        for ( IndexType k = 0; k < n; ++k )
        {
            data1[k] = sparseV.getValue( k ).getValue<ValueType>();
        }

        BOOST_CHECK_EQUAL( 0, data.maxDiffNorm( data1 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SparseConstructorTest )
{
    typedef RealType ValueType;

    ValueType values_raw[] = { 3, 4, 3, 1, 2 };
    IndexType indexes_raw[] = { 7, 3, 11, 13, 5 };

    IndexType nnz = sizeof( values_raw ) / sizeof( ValueType );
    IndexType n   = 20;

    hmemo::HArray<ValueType> values( nnz, values_raw );
    hmemo::HArray<IndexType> indexes( nnz, indexes_raw );

    Scalar zero = 1;

    dmemo::DistributionPtr dist( new dmemo::NoDistribution( n ) );

    // The constructor copies the arrays and sorts the entries

    SparseVector<ValueType> s( dist, indexes, values, zero );

    const hmemo::HArray<IndexType>& spIndexes = s.getNonZeroIndexes();

    // indexes must be sorted:  spIndexes[0] < spIndexes[1] < ... < spIndexes[nnz-1]

    BOOST_CHECK( utilskernel::HArrayUtils::isSorted( spIndexes, common::CompareOp::LT ) );
    BOOST_CHECK_EQUAL( s.getZero(), zero );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( RedistributeTest )
{
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

    IndexType n = 20;   // let it really small, this test is very inefficient

    utilskernel::LArray<ValueType> data( n, ValueType( 0 ) );

    common::Math::srandom( 13151 );   // This test only works if all processors have same random numbers

    const float fillRate = 0.2;   // makes it worth to be sparse

    data.setSparseRandom( fillRate, 1 );

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        utilskernel::LArray<ValueType> data1( n );

        dmemo::DistributionPtr dist = dists[i];

        SparseVector<ValueType> sparseV( data );

        sparseV.redistribute( dist );

        // get each value from the distributed sparse vector

        for ( IndexType k = 0; k < n; ++k )
        {
            data1[k] = sparseV.getValue( k ).getValue<ValueType>();
        }

        BOOST_CHECK_EQUAL( 0, data.maxDiffNorm( data1 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( diffTest )
{
    // Test of different binary operations with sparse vectors
    // For comparison the same operations are computed with dense vectors
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

    IndexType n = 10;

    IndexType rawNonZeroIndexes[] = { 0, 5, 6 };
    ValueType rawNonZeroValues[] = { 5, 6, 7 };
    ValueType zero = 2;

    SparseVector<ValueType> xS1( n, 3, rawNonZeroIndexes, rawNonZeroValues, zero );
    SparseVector<ValueType> xS2( xS1 );

    BOOST_CHECK_EQUAL( xS2.getZero(), xS1.getZero() );

    xS1 -= xS2;

    BOOST_CHECK( xS1.maxNorm() < Scalar( 1e-4 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( binOpSparseTest )
{
    // Test of different binary operations with sparse vectors
    // For comparison the same operations are computed with dense vectors
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

    IndexType n = 10;

    IndexType rawNonZeroIndexes1[] = { 0, 5, 6 };
    ValueType rawNonZeroValues1[] = { 5, 6, 7 };

    IndexType rawNonZeroIndexes2[] = { 0, 4, 6 };
    ValueType rawNonZeroValues2[] = { 5, 4, 9 };

    // Note: binary operations should give sparse vector with maximal 4 elements

    for ( IndexType icase = 0; icase < 5; ++icase )
    {
        ValueType zero1 = 1.0;
        ValueType zero2 = 1.0;

        SparseVector<ValueType> xS1( n, 3, rawNonZeroIndexes1, rawNonZeroValues1, zero1 );
        SparseVector<ValueType> xS2( n, 3, rawNonZeroIndexes2, rawNonZeroValues2, zero2 );

        // use copy constructor to build equivalent dense vectors

        DenseVector<ValueType> xD1( xS1 );
        DenseVector<ValueType> xD2( xS2 );

        SCAI_LOG_DEBUG( logger, "Run test case " << icase << " for binop on sparse vectors" )

        switch ( icase )
        {
            case 0 : xS1 += xS2; 
                     BOOST_CHECK( xS1.getNonZeroIndexes().size() <= 4 );
                     BOOST_CHECK( abs( xS1.getZero() - Scalar( zero1 + zero2 ) ) < Scalar( 0.0001 ) );
                     xD1 += xD2;
                     break;
            case 1 : xS1.cwiseProduct( xS2 );
                     BOOST_CHECK( xS1.getNonZeroIndexes().size() <= 4 );
                     BOOST_CHECK( abs( xS1.getZero() - Scalar( zero1 * zero2 ) ) < Scalar( 0.0001 ) );
                     xD1 *= xD2;
                     break;
            case 2 : xS1 = 5 * xS1 - 2 * xS2; 
                     BOOST_CHECK( xS1.getNonZeroIndexes().size() <= 4 );
                     BOOST_CHECK( abs( xS1.getZero() - Scalar( 5 * zero1 - 2 * zero2 ) ) < Scalar( 0.0001 ) );
                     xD1 = 5 * xD1 - 2 * xD2;
                     break;
            case 3 : xS1.invert();   // this op is okay if zero element is not 0
                     xD1.invert();
                     break;
            case 4 : xS1 += xD2;     // add with dense vector okay, but might be entry for each elem
                     xD1 += xD2;
                     break;
            default: 
                     BOOST_CHECK( false );  // just fail if not handled
        }

        xD1 -= xS1;

        BOOST_CHECK( xD1.maxNorm() < Scalar( 1e-4 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( binOpDenseTest )
{
    // Test of different binary operations where result is dense vector
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

    IndexType n = 10;

    IndexType rawNonZeroIndexes1[] = { 0, 5, 6 };
    ValueType rawNonZeroValues1[] = { 5, 6, 7 };

    IndexType rawNonZeroIndexes2[] = { 0, 4, 6 };
    ValueType rawNonZeroValues2[] = { 5, 4, 9 };

    // Note: binary operations should give sparse vector with maximal 4 elements

    for ( IndexType icase = 0; icase < 4; ++icase )
    {
        ValueType zero1 = 0.0;
        ValueType zero2 = 1.0;

        SparseVector<ValueType> xS1( n, 3, rawNonZeroIndexes1, rawNonZeroValues1, zero1 );
        SparseVector<ValueType> xS2( n, 3, rawNonZeroIndexes2, rawNonZeroValues2, zero2 );

        DenseVector<ValueType> xD1( xS1 );
        DenseVector<ValueType> xD2( xS2 );

        // use copy constructor to build equivalent dense vectors

        DenseVector<ValueType> result1( n, ValueType( 5 ) );
        DenseVector<ValueType> result2( n, ValueType( 5 ) );

        SCAI_LOG_DEBUG( logger, "Run test case " << icase << " for binop on sparse vectors" )

        switch ( icase )
        {
            case 0 : result1 += xS1;
                     result2 += xD1;
                     break;
            case 1 : result1 += xS2;
                     result2 += xD2;
                     break;
            case 2 : result1 *= xS1;
                     result2 *= xD1;
                     break;
            case 3 : result1 *= xS2;
                     result2 *= xD2;
                     break;
            default: 
                     BOOST_CHECK( false );  // just fail if not handled
        }

        result1 -= result2;

        BOOST_CHECK( result1.maxNorm() < Scalar( 1e-4 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( reduceTest )
{
    // Test of different reduction operations on sparse vector
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

    IndexType n = 10;

    IndexType rawNonZeroIndexes1[] = { 0, 5, 6 };
    ValueType rawNonZeroValues1[] = { 5, -6, 7 };

    IndexType rawNonZeroIndexes2[] = { 0, 4, 6 };
    ValueType rawNonZeroValues2[] = { -5, 4, 9 };

    // Note: binary operations should give sparse vector with maximal 4 elements

    ValueType zero1 = 0;
    ValueType zero2 = 0;

    SparseVector<ValueType> xS1( n, 3, rawNonZeroIndexes1, rawNonZeroValues1, zero1 );
    SparseVector<ValueType> xS2( n, 3, rawNonZeroIndexes2, rawNonZeroValues2, zero2 );

    DenseVector<ValueType> xD1( xS1 );
    DenseVector<ValueType> xD2( xS2 );

    BOOST_CHECK_EQUAL( xS1.min(), xD1.min() );
    BOOST_CHECK_EQUAL( xS2.min(), xD2.min() );
    BOOST_CHECK_EQUAL( xS1.max(), xD1.max() );
    BOOST_CHECK_EQUAL( xS2.max(), xD2.max() );
    BOOST_CHECK_EQUAL( xS1.sum(), xD1.sum() );
    BOOST_CHECK_EQUAL( xS1.l1Norm(), xD1.l1Norm() );
    BOOST_CHECK_EQUAL( xS1.l2Norm(), xD1.l2Norm() );
    BOOST_CHECK_EQUAL( xS2.maxNorm(), xD2.maxNorm() );
    BOOST_CHECK_EQUAL( xS1.dotProduct( xS1 ), xD1.dotProduct(xD1 ) );
    BOOST_CHECK_EQUAL( xS1.dotProduct( xS2 ), xD1.dotProduct(xD2 ) );
    BOOST_CHECK_EQUAL( xS1.dotProduct( xD2 ), xD1.dotProduct(xD2 ) );
    BOOST_CHECK_EQUAL( xD1.dotProduct( xS2 ), xD1.dotProduct(xD2 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
