/**
 * @file DenseVectorTest.cpp
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
 * @brief Contains specific tests for class DenseVector (mainly constructors)
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
#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/fft.hpp>

#include <scai/utilskernel.hpp>

#include <scai/testsupport/uniquePath.hpp>
#include <scai/testsupport/GlobalTempDir.hpp>

using namespace scai;
using namespace lama;

using hmemo::HArray;

using scai::testsupport::uniquePath;
using scai::testsupport::GlobalTempDir;

using boost::test_tools::per_element;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DenseVectorTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.DenseVectorTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fillRepTest, ValueType, scai_numeric_test_types )
{
    // replicated dense vector

    const IndexType N   = 12;
    const ValueType VAL = 2;
   
    auto dV = denseVector<ValueType>( N, VAL );

    BOOST_CHECK_EQUAL( N, dV.size() );

    std::vector<ValueType> expected( N, VAL );

    BOOST_TEST( hmemo::hostReadAccess( dV.getLocalValues() ) == expected, per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fillDistTest, ValueType, scai_numeric_test_types )
{
    // distributed dense vector

    const IndexType N   = 15;
    const ValueType VAL = 1;
   
    auto dist = dmemo::blockDistribution( N );
    auto dV = denseVector<ValueType>( dist, VAL );

    BOOST_CHECK_EQUAL( N, dV.size() );

    std::vector<ValueType> expected( dist->getLocalSize(), VAL );

    BOOST_TEST( hmemo::hostReadAccess( dV.getLocalValues() ) == expected, per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fillRepFnTest, ValueType, scai_numeric_test_types )
{
    // fill replicated dense vector with function

    const IndexType N   = 12;

    auto lambda = []( IndexType i ) { return ValueType( i ); };

    auto dV = denseVectorFill<ValueType>( N, lambda );

    BOOST_CHECK_EQUAL( N, dV.size() );

    std::vector<ValueType> expected;

    for ( IndexType i = 0; i < N; ++i )
    {
        expected.push_back( ValueType( i ) );
    }

    BOOST_TEST( hmemo::hostReadAccess( dV.getLocalValues() ) == expected, per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fillDistFnTest, ValueType, scai_numeric_test_types )
{
    // distributed dense vector

    const IndexType N   = 15;

    auto lambda = []( IndexType i ) { return ValueType( i ); };

    auto dist = dmemo::blockDistribution( N );
    auto dV = denseVectorFill<ValueType>( dist, lambda );

    BOOST_CHECK_EQUAL( N, dV.size() );

    std::vector<ValueType> expected;

    for ( IndexType i = 0; i < dist->getLocalSize(); ++i )
    {
        expected.push_back( ValueType( dist->local2Global( i ) ) );
    }

    BOOST_TEST( hmemo::hostReadAccess( dV.getLocalValues() ) == expected, per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( consistencyTest )
{
    typedef DefaultReal ValueType;

    IndexType n = 10;

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    // create distributed dense vector

    auto v = denseVector<ValueType>( dist, 1 );

    BOOST_CHECK( v.isConsistent() );

    // usually it is not simple to make a vector inconsistent

    hmemo::HArray<ValueType>& localData = v.getLocalValues();

    // make it inconsistent on one processor only, isConsistent can deal with it

    if ( comm->getRank() == ( comm->getSize() / 2 ) )
    {
        localData.resize( 17 );
    }
 
    // consistency check fails on all processors

    BOOST_CHECK( ! v.isConsistent() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetGetValueTest )
{
    // Note: it is sufficient to consider one value type

    typedef DefaultReal ValueType;

    IndexType n = 5;   // let it really small, this test is very inefficient

    common::Math::srandom( 13151 );   // This test only works if all processors have same random numbers

    auto data = utilskernel::randomHArray<ValueType>( n, 1 );

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        HArray<ValueType> data1( n );

        dmemo::DistributionPtr dist = dists[i];

        DenseVector<ValueType> distV( dist, 0 );

        // set each value in the distributed vector

        for ( IndexType k = 0; k < n; ++k )
        {
            distV.setValue( k, ValueType( data[k] ) );
        }

        // get each value from the distributed vector

        for ( IndexType k = 0; k < n; ++k )
        {
            data1[k] = distV.getValue( k );
        }

        BOOST_TEST( hostReadAccess( data ) == hostReadAccess( data1 ), per_element() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetAndBuildTest )
{
    // Note: it is sufficient to consider one value type

    typedef DefaultReal ValueType;

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    IndexType n = 16;

    dmemo::DistributionPtr repDist( new dmemo::NoDistribution( n ) );

    DenseVector<ValueType> repV;

    // Note: all processors must have the same random numbers

    common::Math::srandom( 13141 );   // This test only works if all processors have same random numbers

    repV.setRandom( repDist, 1 );

    BOOST_CHECK_EQUAL( n, repV.getLocalValues().size() );

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        // distV = repV, but distributed

        DenseVector<ValueType> distV( repV );
        distV.redistribute( dist );

        DenseVector<ValueType> newV;
        newV.allocate( dist );

        hmemo::HArray<ValueType> tmp;

        distV.buildLocalValues( tmp );
        newV.setDenseValues( tmp );

        // replicate newV, so we can compare with repV

        newV.redistribute( repV.getDistributionPtr() );

        BOOST_CHECK_EQUAL( n, newV.getLocalValues().size() );

        BOOST_TEST( hostReadAccess( newV.getLocalValues() ) == hostReadAccess( repV.getLocalValues() ), per_element() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( RangeTest )
{
    // Note: it is sufficient to consider one value type

    typedef DefaultReal ValueType;

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    IndexType n = 16;

    DenseVector<ValueType> repV = denseVectorLinear<ValueType>( n, 0, 1, ctx );  // just a sequence: 0, 1, ...

    BOOST_CHECK_EQUAL( n, repV.getLocalValues().size() );

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        // distV = repV, but distributed

        auto distV = distribute<DenseVector<ValueType>>( repV, dist );

        BOOST_CHECK_EQUAL( distV.min() , 0 );
        BOOST_CHECK_EQUAL( distV.max() , n - 1 );

        BOOST_CHECK_EQUAL( dist->getLocalSize(), distV.getLocalValues().size() );

        // distV1 = [ 0, ..., n-1] distributed, must be same

        DenseVector<ValueType> distV1 = denseVectorLinear<ValueType>( dist, 0, 1, ctx );

        BOOST_TEST( hostReadAccess( distV1.getLocalValues() ) == hostReadAccess( distV.getLocalValues() ), per_element() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ScanTest )
{
    // Note: it is sufficient to consider one value type

    typedef DefaultReal ValueType;

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType n = 100;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        DenseVector<ValueType> distV = denseVectorLinear<ValueType>( dist, 1, 1, ctx );

        try
        {
            distV.scan();

            BOOST_CHECK_EQUAL( distV.size(), n );
            BOOST_CHECK_EQUAL( distV.getLocalValues().size(), dist->getLocalSize() );
 
            distV.replicate();

            {
                hmemo::ReadAccess<ValueType> rVector( distV.getLocalValues() );

                for ( IndexType i = 0; i < n; ++i )
                {
                    ValueType expected = static_cast<ValueType>( ( i + 1 ) * ( i + 2 ) / 2  );
                    ValueType computed = rVector[i];
                    BOOST_CHECK_EQUAL( expected, computed );
                }
            }
        }
        catch ( common::Exception& e )
        {
            SCAI_LOG_INFO( logger, "scan unsupported for this distribution: " << *dist << ", no block distribution" )
            BOOST_CHECK_EQUAL( dist->getBlockDistributionSize(), invalidIndex );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fileConstructorTest, ValueType, scai_numeric_test_types )
{
    // Note: here we only test constructor DenseVector( "fileName" )
    //       as readFromFile is already tested in PartitionIO

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType n = 10;

    const auto fileName = uniquePath(GlobalTempDir::getPath(), "myVector") + ".psc";

    float fillRate = 0.2;

    hmemo::HArray<ValueType> denseData( n, ValueType( 0 ), ctx );

    common::Math::srandom( 31991 );   // makes sure that all processors generate same data

    utilskernel::HArrayUtils::setSparseRandom( denseData, fillRate, 1 );

    if ( comm->getRank() == 0 )
    {
        FileIO::write( denseData, fileName );
    }

    comm->synchronize();

    auto vector1 = read<DenseVector<ValueType>>( fileName );

    BOOST_CHECK_EQUAL( n, vector1.size() );

    auto dist = std::make_shared<dmemo::BlockDistribution>( n, comm );

    vector1.redistribute( dist );
    vector1.prefetch( ctx );
    vector1.wait();

    DenseVector<ValueType> vector2( ctx );
    vector2.assign( denseData );
    vector2.redistribute( dist );

    // vector1 and vector2 must be equal

    const HArray<ValueType>& localValues1 = vector1.getLocalValues();
    const HArray<ValueType>& localValues2 = vector2.getLocalValues();

    BOOST_TEST( hostReadAccess( localValues1 ) == hostReadAccess( localValues2 ), per_element() );

    if ( comm->getRank() == 0 )
    {
        int rc = FileIO::removeFile( fileName );
        BOOST_CHECK_EQUAL( 0, rc );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( matExpConstructorTest )
{
    // Note: it is sufficient to consider one matrix type and one value type

    typedef DefaultReal ValueType;

    IndexType nCols = 4;
    IndexType nRows = 5;

    dmemo::TestDistributions colDists( nCols );
    dmemo::TestDistributions rowDists( nRows );

    for ( size_t i = 0; i < rowDists.size(); ++i )
    {
        auto rowDist = rowDists[i];

        for ( size_t j = 0; j < colDists.size(); ++j )
        {
            auto colDist = colDists[i];

            auto mat = zero<CSRSparseMatrix<DefaultReal>>( rowDist, colDist );

            auto x   = denseVector<ValueType>( colDist, 3 );

            SCAI_LOG_INFO( logger, "linear algebra expression: alpha * Matrix * Vector" );

            auto y = denseVectorEval( 2 * mat * x );

            BOOST_CHECK_EQUAL( y.getDistribution(), *rowDist );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( scalarExpConstructorTest )
{
    typedef DefaultReal ValueType;

    IndexType n = 5;

    dmemo::TestDistributions dists( n );

    const ValueType v1 = 3;
    const ValueType v2 = 2;
    const ValueType v3 = v1 + v2;

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        auto x = denseVector<ValueType>( dist, v1 );
        SCAI_LOG_INFO( logger, "linear algebra expression: alpha + Vector" );
        auto y = denseVectorEval( v2 + x );
        auto z = denseVectorEval( x + v2 );
        auto r = denseVector<ValueType>( dist, v3 );

        // prove same distribution, same values of r and y/z

        BOOST_TEST( hostReadAccess( r.getLocalValues() ) == hostReadAccess( y.getLocalValues() ), per_element() );
        BOOST_TEST( hostReadAccess( r.getLocalValues() ) == hostReadAccess( z.getLocalValues() ), per_element() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( swapTest )
{
    typedef DefaultReal ValueType;

    const IndexType n = 10;

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    auto dist1 = std::make_shared<dmemo::BlockDistribution>( n, comm );
    auto dist2 = std::make_shared<dmemo::CyclicDistribution>( n, 1, comm );

    auto x1 = denseVector<ValueType>( dist1, 1 );
    auto x2 = denseVector<ValueType>( dist2, 2 );

    x1.swap( x2 );

    BOOST_CHECK_EQUAL( x1.getValue( 8 ), ValueType( 2 ) );
    BOOST_CHECK_EQUAL( x2.getValue( 1 ), ValueType( 1 ) );

    BOOST_CHECK_EQUAL( x1.getLocalValues().size(), dist2->getLocalSize() );
    BOOST_CHECK_EQUAL( x2.getLocalValues().size(), dist1->getLocalSize() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( assignTest )
{
    typedef DefaultReal ValueType;

    const IndexType n = 10;

    auto comm = dmemo::Communicator::getCommunicatorPtr();
    auto dist = std::make_shared<dmemo::BlockDistribution>( n, comm );

    DenseVector<ValueType> x;
    x.allocate( dist );          // be careful, x contains undefined values

    x = 5;  // calls DenseVector::operator= ( Scalar )

    BOOST_CHECK_EQUAL( x.l1Norm(), 5 * n );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( MatrixVectorMultTest, ValueType, scai_numeric_test_types )
{
    // test  vector = scalar * matrix * vector + scalar * vector with all distributions, formats

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType nRows = 16;
    const IndexType nCols = 13;

    // generate random input data, same on all processors

    common::Math::srandom( 51413 );

    auto A = zero<DenseMatrix<ValueType>>( nRows, nCols, ctx );

    MatrixCreator::fillRandom( A, 0.1 );

    DenseVector<ValueType> x( ctx );
    DenseVector<ValueType> y( ctx );

    x.setRandom( A.getColDistributionPtr(), 1 );
    y.setRandom( A.getRowDistributionPtr(), 1 );

    auto res = denseVectorEval( 2 * A * x - y );

    // Now we do the same with all other matrices and all kind of distributions

    dmemo::TestDistributions colDists( nCols );
    dmemo::TestDistributions rowDists( nRows );

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    for ( size_t i = 0; i < rowDists.size(); ++i )
    {
        dmemo::DistributionPtr rowDist = rowDists[i];

        for ( size_t j = 0; j < colDists.size(); ++j )
        {
            dmemo::DistributionPtr colDist = colDists[j];

            Matrices<ValueType> matrices( ctx );   // operations only with same ValueType

            for ( size_t k = 0; k < matrices.size(); ++k )
            {
                Matrix<ValueType>& A1 = *matrices[k];

                A1.assignDistribute( A, rowDist, colDist );
                 
                auto x1 = distribute<DenseVector<ValueType>>( x, colDist );
                auto y1 = distribute<DenseVector<ValueType>>( y, rowDist );

                SCAI_LOG_DEBUG( logger, "matrixTimesVector with this matrix: " << A1 )

                auto res1 = denseVectorEval( 2 * A1 * x1 - y1 );

                // redistribute res1 with same dist as res for comparison

                res1.redistribute( res.getDistributionPtr() );

                BOOST_CHECK( res.maxDiffNorm( res1 ) < eps );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( VectorMatrixMultTest )
{
    typedef DefaultReal ValueType;

    // test of: vector = scalar * vector * matrix + scalar * vector with all distributions, formats

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType nRows = 2;
    const IndexType nCols = 3;

    // generate random input data, same on all processors

    common::Math::srandom( 51413 );
    const IndexType bound = 1;

    // Do vector = scalar * vector * matrix + scalar * vector with replicated storage data

    DenseMatrix<ValueType> A;
    A.setContextPtr( ctx );
    A.allocate( nRows, nCols );
    MatrixCreator::fillRandom( A, 0.1 );

    // Before call of setRandom, x and y must be allocated, but not initialized

    DenseVector<ValueType> x( ctx );
    DenseVector<ValueType> y( ctx );

    x.setRandom( A.getRowDistributionPtr(), bound );
    y.setRandom( A.getColDistributionPtr(), bound );

    auto res = denseVectorEval( transpose( A ) * 2 * x - y );

    // Now we do the same with all distributed matrices/vectors and all kind of distributions

    dmemo::TestDistributions colDists( nCols );
    dmemo::TestDistributions rowDists( nRows );

    for ( size_t i = 0; i < rowDists.size(); ++i )
    {
        SCAI_LOG_DEBUG( logger, "Row distribution [ " << i << " ] = " << *rowDists[i] )
    }

    for ( size_t j = 0; j < colDists.size(); ++j )
    {
        SCAI_LOG_DEBUG( logger, "Col distribution [ " << j << " ] = " << *colDists[j] )
    }

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    for ( size_t i = 0; i < rowDists.size(); ++i )
    {
        auto rowDist = rowDists[i];

        for ( size_t j = 0; j < colDists.size(); ++j )
        {
            auto colDist = colDists[j];

            Matrices<ValueType> matrices( ctx );  // currently restricted, only of ValueType

            for ( size_t k = 0; k < matrices.size(); ++k )
            {
                Matrix<ValueType>& A1 = *matrices[k];

                A1.setCommunicationKind( SyncKind::SYNCHRONOUS );

                A1.assignDistribute( A, rowDist, colDist );

                auto x1 = distribute<DenseVector<ValueType>>( x, rowDist );
                auto y1 = distribute<DenseVector<ValueType>>( y, colDist );

                SCAI_LOG_INFO( logger, "vectorTimesMatrix[" << i << "," << j << "," << k << "] with this matrix: " << A1 << ", y1 = " << y1 )

                auto res1 =  denseVectorEval( 2 * transpose( A1 ) * x1 - y1 );

                SCAI_LOG_INFO( logger, "res1 = 2 * transpose( A1 ) * x1 - y1: " << res1 )

                BOOST_CHECK_EQUAL( res1.size(), A1.getNumColumns() );

                // A1.vectorTimesMatrix( res1, Scalar( 2 ), x1, Scalar( -1 ), y1 );

                res1.redistribute( res.getDistributionPtr() );
                res1 -= res;

                BOOST_CHECK( res1.maxNorm() < eps );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( VectorMatrixMultTransposeTest, ValueType, scai_numeric_test_types )
{
    // test transpose( A ) * x, once with implicitliy tranposed matrix and one with explicitly transposed

    const IndexType nRows = 7;
    const IndexType nCols = 4;

    // generate random input data, same on all processors

    common::Math::srandom( 1413 );

    DenseMatrix<ValueType> A;
    A.allocate( nRows, nCols );
    MatrixCreator::fillRandom( A, 0.1 );

    DenseVector<ValueType> x;
    DenseVector<ValueType> y;

    x.setRandom( A.getRowDistributionPtr(), 1 );
    y.setRandom( A.getColDistributionPtr(), 1 );

    DenseMatrix<ValueType> At;
    At.assignTranspose( A );   // builds explicitly transposed matrix

    auto res1 = denseVectorEval( 2 * transpose( A ) * x - y );
    auto res2 = denseVectorEval( 2 * At * x - y );

    const HArray<ValueType>& v1 = res1.getLocalValues();
    const HArray<ValueType>& v2 = res2.getLocalValues();

    // results might be slightly different due to rounding errors
    // BOOST_TEST( hostReadAccess( v1 ) == hostReadAccess( v2 ), per_element() );

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    BOOST_CHECK( utilskernel::HArrayUtils::maxDiffNorm( v1, v2 ) < eps );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE ( VectorPlusScalarExpressionTest )
{
    typedef DefaultReal ValueType;

    const IndexType n = 4;
    ValueType sourceVals[] = { 3, 1, 4, 2 };
    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    DenseVector<ValueType> x( ctx );
    x.setRawData( n, sourceVals );

    ValueType alpha = 34.7;
    ValueType beta  = 5.2;

    // test expression
    DenseVector<ValueType> res1;
    res1 = alpha * x + beta;
    // test constructor with expression
    auto res2 = denseVectorEval( alpha * x + beta );
    // test alias
    x = alpha * x + beta;

    for ( IndexType i = 0; i < n; ++i )
    {
        ValueType erg = alpha * sourceVals[i] + beta;
        BOOST_CHECK_EQUAL( erg, res1.getValue( i ) );
        BOOST_CHECK_EQUAL( erg, res2.getValue( i ) );
        BOOST_CHECK_EQUAL( erg,    x.getValue( i ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE ( sortTest, ValueType, scai_array_test_types )
{
    if ( common::isComplex( common::TypeTraits<ValueType>::stype ) )
    {
        return;    // sort of complex numbers is not relevant
    }

    ValueType sortValues[]   = { 5, 9, 4, 8, 1, 2, 3 };
    ValueType sortedValues[] = { 1, 2, 3, 4, 5, 8, 9 };
    IndexType sortPerm[]     = { 4, 5, 6, 2, 0, 3, 1 };

    const IndexType n = sizeof( sortValues ) / sizeof( ValueType );

    DenseVector<ValueType> sortVector;
    hmemo::HArray<ValueType> sortArray( n, sortValues );
    sortVector.assign( sortArray );

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    dmemo::DistributionPtr blockDist( new dmemo::BlockDistribution( n, comm ) );

    sortVector.redistribute( blockDist );

    {
        bool descending = false;
        DenseVector<ValueType> tmp( sortVector );
        tmp.sort( descending );    // parallel sorting with global permutation
        BOOST_CHECK( tmp.isSorted( descending ) );
    }

    bool ascending = true;

    DenseVector<IndexType> perm;

    sortVector.sort( perm, ascending );    // parallel sorting with global permutation

    BOOST_REQUIRE_EQUAL( n, perm.size() );
    BOOST_CHECK( sortVector.isSorted( ascending ) );

    dmemo::DistributionPtr repDist( new dmemo::NoDistribution( n ) );

    sortVector.redistribute( repDist );
    perm.redistribute( repDist );

    BOOST_REQUIRE( utilskernel::HArrayUtils::isSorted( sortVector.getLocalValues(), common::CompareOp::LE ) );

    hmemo::ReadAccess<ValueType> rSorted( sortVector.getLocalValues() );
    hmemo::ReadAccess<IndexType> rPerm( perm.getLocalValues() );

    BOOST_REQUIRE_EQUAL( n, rSorted.size() );
    BOOST_REQUIRE_EQUAL( n, rPerm.size() );

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( rPerm[i], sortPerm[i] );
        BOOST_CHECK_EQUAL( rSorted[i], sortedValues[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE ( gatherTest, ValueType, scai_array_test_types )
{
    DenseVector<ValueType> source( HArray<ValueType>( { 5, 9, 4, 8, 1, 2, 3 } ) );
    DenseVector<IndexType> index( HArray<IndexType>( { 3, 4, 1, 0, 6, 2 } ) );

    HArray<ValueType> expTargetValues( { 8, 1, 9, 5, 3, 4 } );

    dmemo::TestDistributions sourceDistributions( source.size() );
    dmemo::TestDistributions indexDistributions( index.size() );

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    for ( size_t sd = 0; sd < sourceDistributions.size(); ++sd )
    {
        for ( size_t id = 0; id < indexDistributions.size(); ++id )
        {
            dmemo::DistributionPtr sourceDist = sourceDistributions[sd];
            dmemo::DistributionPtr indexDist = indexDistributions[id];

            source.redistribute( sourceDist );
            index.redistribute( indexDist );

            DenseVector<ValueType> target( indexDist, 0 );

            SCAI_LOG_INFO( logger, "gather source[index] with source = " << source << ", index = " << index )

            source.gatherFrom( target, index, common::BinaryOp::ADD );

            BOOST_CHECK_EQUAL( target.size(), index.size() );
            BOOST_CHECK_EQUAL( target.getDistribution(), index.getDistribution() );
            BOOST_CHECK_EQUAL( target.getDistribution().getLocalSize(), target.getLocalValues().size() );

            auto rTarget = hostReadAccess( target.getLocalValues() );
            auto rExpTarget = hostReadAccess( expTargetValues );

            for ( IndexType i = 0; i < target.size(); ++i )
            {
                IndexType localIndex = target.getDistribution().global2Local( i );

                if ( localIndex != invalidIndex )
                {
                    BOOST_CHECK_MESSAGE( rTarget[localIndex] == rExpTarget[i],
                                         *comm << ": targetLocal[" << localIndex << "] = " << rTarget[localIndex]
                                         << " must be equal to targetValues[" << i << "] = " << rExpTarget[i] );
                }
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE ( gatherByPlanTest )
{
    typedef DefaultReal ValueType;

    DenseVector<ValueType> source( HArray<ValueType>( { 5, 9, 4, 8, 1, 2, 3 } ) );
    DenseVector<IndexType> index( HArray<IndexType>( { 3, 4, 1, 0, 6, 2 } ) );

    HArray<ValueType> expTargetValues( { 8, 1, 9, 5, 3, 4 } );

    dmemo::TestDistributions sourceDistributions( source.size() );
    dmemo::TestDistributions indexDistributions( index.size() );

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    for ( size_t sd = 0; sd < sourceDistributions.size(); ++sd )
    {
        for ( size_t id = 0; id < indexDistributions.size(); ++id )
        {
            dmemo::DistributionPtr sourceDist = sourceDistributions[sd];
            dmemo::DistributionPtr indexDist = indexDistributions[id];

            if ( sourceDist->getCommunicator() != indexDist->getCommunicator() )
            {
                // global addressing plan can only be used for distribution on same processors
                continue;
            }

            source.redistribute( sourceDist );
            index.redistribute( indexDist );

            DenseVector<ValueType> target( indexDist, 0 );

            SCAI_LOG_INFO( logger, "gather source[index] with source = " << source << ", index = " << index )

            auto plan = source.globalAddressingPlan( index );
            source.gatherByPlan( target, plan, common::BinaryOp::ADD );

            DenseVector<ValueType> expTarget( expTargetValues );
            expTarget.redistribute( indexDist );

            BOOST_TEST( hmemo::hostReadAccess( target.getLocalValues() ) == hmemo::hostReadAccess( expTarget.getLocalValues() ), per_element() );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( scatterTest )
{
    typedef DefaultReal ValueType;

    // targetVector[ indexVector ] = sourceVector

    hmemo::HArray<ValueType> sourceArray( { 5, 9, 4, 8, 1, 2, 2, 3, 4 } );
    hmemo::HArray<IndexType> indexArray(  { 3, 4, 1, 0, 5, 2, 1, 0, 2 } );

    const auto op = common::BinaryOp::ADD;   // reduction operator for multiple entries
    const bool unique = false;

    BOOST_REQUIRE_EQUAL( sourceArray.size(), indexArray.size() );

    const IndexType m = utilskernel::HArrayUtils::max( indexArray ) + 1;

    hmemo::HArray<ValueType> expTargetArray( m, ValueType( 0 ) );
    utilskernel::HArrayUtils::scatter( expTargetArray, indexArray, unique, sourceArray, op );

    DenseVector<ValueType> source( std::move( sourceArray ) );
    DenseVector<IndexType> index( std::move( indexArray ) );
    DenseVector<ValueType> expTarget( std::move( expTargetArray ) );

    dmemo::TestDistributions targetDistributions( expTarget.size() );
    dmemo::TestDistributions indexDistributions( source.size() );

    for ( size_t sd = 0; sd < targetDistributions.size(); ++sd )
    {
        for ( size_t id = 0; id < indexDistributions.size(); ++id )
        {
            dmemo::DistributionPtr targetDist = targetDistributions[sd];
            dmemo::DistributionPtr indexDist = indexDistributions[id];

            SCAI_LOG_DEBUG( logger, "target( dist = " << *targetDist << ")[ index ] = source( dist = " << *indexDist << " )" )

            source.redistribute( indexDist );
            index.redistribute( indexDist );

            auto target = denseVector<ValueType>( targetDist, 0 );

            if ( targetDist->isReplicated() != indexDist->isReplicated() )
            {
                // either both are replicated or both distributed

                BOOST_CHECK_THROW(
                {
                    target.scatter( index, unique, source, common::BinaryOp::ADD );
                }, common::Exception );

                continue;
            }

            target.scatter( index, unique, source, common::BinaryOp::ADD );
            expTarget.redistribute( targetDist );

            BOOST_TEST( hmemo::hostReadAccess( target.getLocalValues() ) == hmemo::hostReadAccess( expTarget.getLocalValues() ), per_element() );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( scatterByPlanTest )
{
    typedef DefaultReal ValueType;

    // targetVector[ indexVector ] = sourceVector

    hmemo::HArray<ValueType> sourceArray( { 5, 9, 4, 8, 1, 2, 2, 3, 4 } );
    hmemo::HArray<IndexType> indexArray(  { 3, 4, 1, 0, 5, 2, 1, 0, 2 } );

    const auto op = common::BinaryOp::ADD;   // reduction operator for multiple entries
    const bool unique = false;

    BOOST_REQUIRE_EQUAL( sourceArray.size(), indexArray.size() );

    const IndexType m = utilskernel::HArrayUtils::max( indexArray ) + 1;

    hmemo::HArray<ValueType> expTargetArray( m, ValueType( 0 ) );
    utilskernel::HArrayUtils::scatter( expTargetArray, indexArray, unique, sourceArray, op );

    DenseVector<ValueType> source( std::move( sourceArray ) );
    DenseVector<IndexType> index( std::move( indexArray ) );
    DenseVector<ValueType> expTarget( std::move( expTargetArray ) );

    dmemo::TestDistributions targetDistributions( expTarget.size() );
    dmemo::TestDistributions indexDistributions( source.size() );

    for ( size_t sd = 0; sd < targetDistributions.size(); ++sd )
    {
        for ( size_t id = 0; id < indexDistributions.size(); ++id )
        {
            dmemo::DistributionPtr targetDist = targetDistributions[sd];
            dmemo::DistributionPtr indexDist = indexDistributions[id];

            if ( targetDist->getCommunicator() != indexDist->getCommunicator() )
            {
                continue;
            }

            SCAI_LOG_DEBUG( logger, "target( dist = " << *targetDist << ")[ index ] = source( dist = " << *indexDist << " )" )

            source.redistribute( indexDist );
            index.redistribute( indexDist );

            auto target = denseVector<ValueType>( targetDist, 0 );

            auto plan = target.globalAddressingPlan( index, unique );
            target.scatterByPlan( plan, source, common::BinaryOp::ADD );

            expTarget.redistribute( targetDist );

            BOOST_TEST( hmemo::hostReadAccess( target.getLocalValues() ) == hmemo::hostReadAccess( expTarget.getLocalValues() ), per_element() );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( moveTest )
{
    // Test of different reduction operations on sparse vector
    // Note: it is sufficient to consider one value type

    typedef DefaultReal ValueType;

    IndexType n = 10;
    ValueType rawValues[] = { 1, 5, 4, 3, 2, 0, 1, 5, 3, 4, 2 };

    hmemo::HArray<ValueType> values( n, rawValues );

    const ValueType* ptr0 = hmemo::ReadAccess<ValueType>( values ).get();

    DenseVector<ValueType> xD( std::move( values ) );

    const ValueType* ptr1 = hmemo::ReadAccess<ValueType>( xD.getLocalValues() ).get();

    BOOST_CHECK_EQUAL( ptr0, ptr1 );  // that is great, all the same data used

    ValueType s = xD.sum();

    DenseVector<ValueType> xD1( std::move( xD ) );
    
    BOOST_CHECK_EQUAL( 0, xD.sum() );
    BOOST_CHECK_EQUAL( s, xD1.sum() );

    DenseVector<ValueType> xD2;
    xD2 = std::move( xD1 );

    BOOST_CHECK_EQUAL( 0, xD1.sum() );
    BOOST_CHECK_EQUAL( s, xD2.sum() );

    const ValueType* ptr2 = hmemo::ReadAccess<ValueType>( xD2.getLocalValues() ).get();

    BOOST_CHECK_EQUAL( ptr1, ptr2 );  // that is great, all the same data used
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( resizeTest )
{
    typedef SCAI_TEST_TYPE ValueType;     // value type does not matter for resize

    using namespace hmemo;

    const IndexType n = 13;

    dmemo::TestDistributions oldDists( n );

    dmemo::TestDistributions fillDists( n + 2 );
    dmemo::TestDistributions truncDists( n - 2 );

    std::srand( 1311 );   // same random numbers on all processors

    DenseVector<ValueType> v1;

    HArray<ValueType> data1( n );
    utilskernel::HArrayUtils::setRandom( data1, 1 );

    for ( size_t j = 0; j < oldDists.size(); ++j )
    {
        for ( size_t k = 0; k < fillDists.size(); ++k )
        {
            v1.assign( data1 );
            v1.resize( fillDists[k] );
            ValueType v = v1[ n + 1 ];
            BOOST_CHECK_EQUAL( v, 0 );
            v = v1[3];
            ValueType r = data1[3];
            BOOST_CHECK_EQUAL( v, r );
        }

        for ( size_t k = 0; k < truncDists.size(); ++k )
        {
            v1.assign( data1 );
            v1.resize( truncDists[k] );
            ValueType x = v1[ n - 3 ];      // get distributed value
            ValueType y = data1[ n - 3 ];   // get replicated value
            BOOST_CHECK_EQUAL( x, y );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( purgeTest )
{
    auto hostCtx = hmemo::Context::getHostPtr();
    auto hostMem = hostCtx->getMemoryPtr();

    typedef SCAI_TEST_TYPE ValueType;     // value type does not matter for resize

    const IndexType n = 13;

    DenseVector<ValueType> v1( n, ValueType( 0 ), hostCtx );

    size_t nAllocates = hostMem->allocates();
    size_t nAllocatedBytes = hostMem->allocatedBytes();

    v1.clear();

    // clear does not free any memory, can be verified by statistics of host memory

    BOOST_CHECK_EQUAL( nAllocates, hostMem->allocates() );
    BOOST_CHECK_EQUAL( nAllocatedBytes, hostMem->allocatedBytes() );

    v1.purge();

    // purge has freed the allocated, can be verified by statistics of host memory

    BOOST_CHECK_EQUAL( nAllocates, hostMem->allocates() + 1 );
    BOOST_CHECK_EQUAL( nAllocatedBytes, hostMem->allocatedBytes() + n * sizeof( ValueType ) );
}

#ifdef SCAI_COMPLEX_SUPPORTED

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fftTest, ValueType, scai_fft_test_types )
{
    typedef common::Complex<RealType<ValueType>> FFTType;

    DenseVector<FFTType> y( HArray<FFTType>( { 0.2, 0.16, 0, 0, 0, 0, 0, 0 } ) );

    const IndexType n = y.size();

    fft( y );

    FFTType yS1 = y[0];
    FFTType yS2 = y[n - 2];

    RealType<ValueType> eps = 0.00001;

    BOOST_CHECK( common::Math::abs( yS1 - FFTType( 0.2 + 0.16 ) ) < eps );
    BOOST_CHECK( common::Math::abs( yS2 - FFTType( 0.2, 0.16 ) ) < eps );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ifftTest, ValueType, scai_fft_test_types )
{
    typedef common::Complex<RealType<ValueType>> FFTType;

    DenseVector<FFTType> y( HArray<FFTType>( { 0.2, 0.16, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } ) );

    const IndexType n = y.size();

    ifft( y );

    BOOST_CHECK_EQUAL( y.size(), n );

    FFTType yS1 = y[0];
    FFTType yS2 = y[n - 2];

    RealType<ValueType> eps = 0.00001;

    BOOST_CHECK( common::Math::abs( yS1 - FFTType( 0.2 + 0.16 ) ) < eps );
    BOOST_CHECK( common::Math::abs( yS2 - FFTType( 0.2, -0.16 ) ) < eps );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( allFFTTest, ValueType, scai_fft_test_types )
{
    typedef common::Complex<RealType<ValueType>> FFTType;

    const IndexType N = 13;   

    // generate random vector

    HArray<ValueType> randomValues( N );
    utilskernel::HArrayUtils::setRandom( randomValues, 1 );

    DenseVector<ValueType> x( randomValues );

    // fill up to N2 with is next power 2 value
 
    const IndexType N2 = 1 << common::Math::nextpow2( N );

    SCAI_LOG_INFO( logger, "allFFTTest, N = " << N << ", N2 = " << N2 )

    DenseVector<FFTType> y;
    y = cast<FFTType>( x );
    y.resize( std::make_shared<dmemo::NoDistribution>( N2 ) );

    fft( y );        // FFT forward
    ifft( y );       // FFT backward

    y.resize( x.getDistributionPtr() );

    // divide by N after fft - ifft to get the original result

    y /= ValueType( N2 );

    auto x1 = convert<DenseVector<ValueType>>( y );

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    // due to rounding errors: not possible is
    // BOOST_TEST( hostReadAccess( x.getLocalValues() ) == hostReadAccess( x1.getLocalValues() ), per_element() );

    BOOST_CHECK( utilskernel::HArrayUtils::maxDiffNorm( x.getLocalValues(), x1.getLocalValues() ) < eps );
}

#endif

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

