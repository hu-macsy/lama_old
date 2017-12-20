/**
 * @file DenseVectorTest.cpp
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

#include <scai/lama/expression/all.hpp>
#include <scai/testsupport/uniquePath.hpp>
#include <scai/testsupport/GlobalTempDir.hpp>

using namespace scai;
using namespace lama;

using scai::testsupport::uniquePath;
using scai::testsupport::GlobalTempDir;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DenseVectorTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.DenseVectorTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( cTorTest, ValueType, scai_numeric_test_types )
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

BOOST_AUTO_TEST_CASE( consistencyTest )
{
    typedef RealType ValueType;

    IndexType n = 10;

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();
    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    // create distributed dense vector

    DenseVector<ValueType> v( dist, ValueType( 1 ) );

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

    typedef RealType ValueType;

    IndexType n = 5;   // let it really small, this test is very inefficient

    utilskernel::LArray<ValueType> data( n );

    common::Math::srandom( 13151 );   // This test only works if all processors have same random numbers

    data.setRandom( 1 );

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        utilskernel::LArray<ValueType> data1( n );

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

        BOOST_CHECK_EQUAL( 0, data.maxDiffNorm( data1 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetAndBuildTest )
{
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

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

        DenseVector<ValueType> distV( repV, dist );

        DenseVector<ValueType> newV( dist );

        hmemo::HArray<ValueType> tmp;

        distV.buildLocalValues( tmp );
        newV.setDenseValues( tmp );

        // replicate newV, so we can compare with repV

        newV.redistribute( repV.getDistributionPtr() );

        BOOST_CHECK_EQUAL( n, newV.getLocalValues().size() );

        BOOST_CHECK_EQUAL( 0, newV.getLocalValues().maxDiffNorm( repV.getLocalValues() ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( RangeTest )
{
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    IndexType n = 16;

    DenseVector<ValueType> repV = linearValuesVector<ValueType>( n, 0, 1, ctx );  // just a sequence: 0, 1, ...

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

        DenseVector<ValueType> distV1 = linearValuesVector<ValueType>( dist, 0, 1, ctx );

        BOOST_CHECK_EQUAL( 0, distV1.getLocalValues().maxDiffNorm( distV.getLocalValues() ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ScanTest )
{
    // Note: it is sufficient to consider one value type

    typedef RealType ValueType;

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType n = 100;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        DenseVector<ValueType> distV = linearValuesVector<ValueType>( dist, 1, 1, ctx );

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
                    Scalar computed = rVector[i];
                    BOOST_CHECK_EQUAL( expected, computed.getValue<ValueType>() );
                }
            }
        }
        catch ( common::Exception& e )
        {
            SCAI_LOG_INFO( logger, "scan unsupported for this distribution: " << *dist << ", no block distribution" )
            BOOST_CHECK_EQUAL( dist->getBlockDistributionSize(), nIndex );
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

    DenseVector<ValueType> vector1( fileName );

    BOOST_CHECK_EQUAL( n, vector1.size() );

    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    vector1.redistribute( dist );
    vector1.prefetch( ctx );
    vector1.wait();

    DenseVector<ValueType> vector2( ctx );
    vector2.assign( denseData );
    vector2.redistribute( dist );

    // vector1 and vector2 must be equal

    const utilskernel::LArray<ValueType> localValues1 = vector1.getLocalValues();
    const utilskernel::LArray<ValueType> localValues2 = vector2.getLocalValues();

    BOOST_CHECK_EQUAL( ValueType( 0 ), localValues1.maxDiffNorm( localValues2 ) );

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

            CSRSparseMatrix<RealType> mat( rowDist, colDist );

            DenseVector<ValueType> x( colDist, 3 );
            SCAI_LOG_INFO( logger, "linear algebra expression: alpha * _Matrix * Vector" );
            DenseVector<ValueType> y( 2 * mat * x );

            BOOST_CHECK_EQUAL( y.getDistribution(), *rowDist );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( scalarExpConstructorTest )
{
    typedef RealType ValueType;

    IndexType n = 5;

    dmemo::TestDistributions dists( n );

    for ( size_t i = 0; i < dists.size(); ++i )
    {
        dmemo::DistributionPtr dist = dists[i];

        DenseVector<ValueType> x( dist, 3 );
        SCAI_LOG_INFO( logger, "linear algebra expression: alpha + Vector" );
        DenseVector<ValueType> y( ValueType( 2 ) + x );
        DenseVector<ValueType> z( x + ValueType( 2 ) );
        DenseVector<ValueType> r( dist, 5 );

        // prove same distribution, same values of r and y/z

        BOOST_CHECK( r.getLocalValues().maxDiffNorm( y.getLocalValues() ) == 0 );
        BOOST_CHECK( r.getLocalValues().maxDiffNorm( z.getLocalValues() ) == 0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( swapTest )
{
    typedef RealType ValueType;

    const IndexType n = 10;

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    dmemo::DistributionPtr dist1( new dmemo::BlockDistribution( n, comm ) );
    dmemo::DistributionPtr dist2( new dmemo::CyclicDistribution( n, 1, comm ) );

    DenseVector<ValueType> x1( dist1, 1 );
    DenseVector<ValueType> x2( dist2, 2 );

    x1.swap( x2 );

    BOOST_CHECK_EQUAL( x1.getValue( 8 ), ValueType( 2 ) );
    BOOST_CHECK_EQUAL( x2.getValue( 1 ), ValueType( 1 ) );

    BOOST_CHECK_EQUAL( x1.getLocalValues().size(), dist2->getLocalSize() );
    BOOST_CHECK_EQUAL( x2.getLocalValues().size(), dist1->getLocalSize() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( assignTest )
{
    typedef RealType ValueType;

    const IndexType n = 10;

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    dmemo::DistributionPtr dist( new dmemo::BlockDistribution( n, comm ) );

    DenseVector<ValueType> x( dist );

    x = 5;  // calls DenseVector::operator= ( Scalar )

    BOOST_CHECK_EQUAL( x.l1Norm(), 5 * n );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( _MatrixVectorMultTest, ValueType, scai_numeric_test_types )
{
    // test  vector = scalar * matrix * vector + scalar * vector with all distributions, formats

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType nRows = 16;
    const IndexType nCols = 13;

    // generate random input data, same on all processors

    common::Math::srandom( 51413 );

    DenseMatrix<ValueType> A;
    A.setContextPtr( ctx );
    A.allocate( nRows, nCols );
    MatrixCreator::fillRandom( A, 0.1 );
    DenseVector<ValueType> x( ctx );
    DenseVector<ValueType> y( ctx );
    x.setRandom( A.getColDistributionPtr(), 1 );
    y.setRandom( A.getRowDistributionPtr(), 1 );
    DenseVector<ValueType> res( ValueType( 2 ) * A * x - y );

    // Now we do the same with all other matrices and all kind of distributions

    dmemo::TestDistributions colDists( nCols );
    dmemo::TestDistributions rowDists( nRows );

    NormType<ValueType> eps = common::TypeTraits<ValueType>::small();

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

                A1.assign( A );
                A1.redistribute( rowDist, colDist );

                DenseVector<ValueType> x1( x, colDist );
                DenseVector<ValueType> y1( y, rowDist );
                DenseVector<ValueType> res1;

                SCAI_LOG_INFO( logger, "matrixTimesVector with this matrix: " << A1 )

                res1 = ValueType( 2 ) * A1 * x1 - y1;

                res1.redistribute( res.getDistributionPtr() );
                res1 -= res;

                BOOST_CHECK( res1.maxNorm() < eps );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( VectorMatrixMultTest )
{
    typedef float ValueType;

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

    DenseVector<ValueType> x( A.getRowDistributionPtr(), ctx );
    DenseVector<ValueType> y( A.getColDistributionPtr(), ctx );

    x.setRandom( A.getRowDistributionPtr(), bound );
    y.setRandom( A.getColDistributionPtr(), bound );

    DenseVector<ValueType> res( ValueType( 2 ) * x * A - y );

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

    NormType<ValueType> eps = common::TypeTraits<ValueType>::small();

    for ( size_t i = 0; i < rowDists.size(); ++i )
    {
        dmemo::DistributionPtr rowDist = rowDists[i];

        for ( size_t j = 0; j < colDists.size(); ++j )
        {
            dmemo::DistributionPtr colDist = colDists[j];

            Matrices<ValueType> matrices( ctx );  // currently restricted, only of ValueType

            for ( size_t k = 0; k < matrices.size(); ++k )
            {
                Matrix<ValueType>& A1 = *matrices[k];

                A1.setCommunicationKind( SyncKind::SYNCHRONOUS );

                A1.assign( A );
                A1.redistribute( rowDist, colDist );

                DenseVector<ValueType> x1( x, rowDist );
                DenseVector<ValueType> y1( y, colDist );
                DenseVector<ValueType> res1;

                SCAI_LOG_INFO( logger, "vectorTimesMatrix[" << i << "," << j << "," << k << "] with this matrix: " << A1 << ", y1 = " << y1 )

                res1 = ValueType( 2 ) * x1 * A1 - y1;

                SCAI_LOG_INFO( logger, "res1 = 2 * x1 * A1 - y1: " << res1 )

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

BOOST_AUTO_TEST_CASE( VectorMatrixMult1Test )
{
    // only serial

    typedef float ValueType;

    // test  vector = scalar * matrix * vector + scalar * vector with all distributions, formats

    hmemo::ContextPtr ctx = hmemo::Context::getContextPtr();

    const IndexType nRows = 7;
    const IndexType nCols = 4;

    // generate random input data, same on all processors

    common::Math::srandom( 1413 );

    DenseMatrix<ValueType> A;
    A.setContextPtr( ctx );
    A.allocate( nRows, nCols );
    MatrixCreator::fillRandom( A, 0.1 );

    DenseVector<ValueType> x( ctx );
    DenseVector<ValueType> y( ctx );

    x.setRandom( A.getRowDistributionPtr(), 1 );
    y.setRandom( A.getColDistributionPtr(), 1 );

    DenseMatrix<ValueType> At( A, true );

    DenseVector<ValueType> res1( ValueType( 2 ) * x * A - y );
    DenseVector<ValueType> res2( ValueType( 2 ) * At * x - y );

    const utilskernel::LArray<ValueType>& v1 = res1.getLocalValues();
    const utilskernel::LArray<ValueType>& v2 = res2.getLocalValues();

    BOOST_CHECK_EQUAL( 0, v1.maxDiffNorm( v2 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE ( VectorPlusScalarExpressionTest )
{
    typedef float ValueType;

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
    DenseVector<ValueType> res2( alpha * x + beta );
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

// BOOST_AUTO_TEST_CASE_TEMPLATE ( gatherTest, ValueType, scai_array_test_types )

BOOST_AUTO_TEST_CASE( gatherTest )
{
    typedef float ValueType;

    ValueType sourceValues[] = { 5, 9, 4, 8, 1, 2, 3 };
    ValueType indexValues[]  = { 3, 4, 1, 0, 6, 2 };
    ValueType targetValues[] = { 8, 1, 9, 5, 3, 4 };

    const IndexType m  = sizeof( sourceValues ) / sizeof( ValueType );
    const IndexType n  = sizeof( indexValues ) / sizeof( ValueType );
    const IndexType n1 = sizeof( targetValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n, n1 );

    DenseVector<ValueType> source;
    hmemo::HArray<ValueType> sourceArray( m, sourceValues );
    source.assign( sourceArray );

    DenseVector<IndexType> index;
    hmemo::HArray<ValueType> indexArray( n, indexValues );
    index.assign( indexArray );

    dmemo::TestDistributions sourceDistributions( m );
    dmemo::TestDistributions indexDistributions( n );

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    for ( size_t sd = 0; sd < sourceDistributions.size(); ++sd )
    {
        for ( size_t id = 0; id < indexDistributions.size(); ++id )
        {
            dmemo::DistributionPtr sourceDist = sourceDistributions[sd];
            dmemo::DistributionPtr indexDist = indexDistributions[id];

            source.redistribute( sourceDist );
            index.redistribute( indexDist );

            DenseVector<ValueType> target( indexDist );
            target = 0;

            SCAI_LOG_INFO( logger, "gather source[index] with source = " << source << ", index = " << index )

            target.gather( source, index, common::BinaryOp::ADD );

            BOOST_CHECK_EQUAL( target.size(), index.size() );
            BOOST_CHECK_EQUAL( target.getDistribution(), index.getDistribution() );
            BOOST_CHECK_EQUAL( target.getDistribution().getLocalSize(), target.getLocalValues().size() );

            hmemo::ReadAccess<ValueType> rTarget( target.getLocalValues() );

            for ( IndexType i = 0; i < n; ++i )
            {
                IndexType localIndex = target.getDistribution().global2local( i );

                if ( localIndex != nIndex )
                {
                    BOOST_CHECK_MESSAGE( rTarget[localIndex] == targetValues[i],
                                         *comm << ": targetLocal[" << localIndex << "] = " << rTarget[localIndex]
                                         << " must be equal to targetValues[" << i << "] = " << targetValues[i] );
                }
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( scatterTest )
{
    typedef float ValueType;

    ValueType sourceValues[] = { 5, 9, 4, 8, 1, 2 };
    ValueType indexValues[]  = { 3, 4, 1, 0, 5, 2 };
    ValueType targetValues[] = { 8, 4, 2, 5, 9, 1 };

    const IndexType n1 = sizeof( sourceValues ) / sizeof( ValueType );
    const IndexType n  = sizeof( indexValues ) / sizeof( ValueType );
    const IndexType m  = sizeof( targetValues ) / sizeof( ValueType );

    BOOST_REQUIRE_EQUAL( n, n1 );

    DenseVector<ValueType> source;
    hmemo::HArray<ValueType> sourceArray( m, sourceValues );
    source.assign( sourceArray );

    DenseVector<IndexType> index;
    hmemo::HArray<ValueType> indexArray( n, indexValues );
    index.assign( indexArray );

    dmemo::CommunicatorPtr comm = dmemo::Communicator::getCommunicatorPtr();

    dmemo::DistributionPtr targetDist( new dmemo::BlockDistribution( m, comm ) );
    dmemo::DistributionPtr indexDist( new dmemo::BlockDistribution( n, comm ) );

    dmemo::TestDistributions targetDistributions( m );
    dmemo::TestDistributions indexDistributions( n );

    for ( size_t sd = 0; sd < targetDistributions.size(); ++sd )
    {
        for ( size_t id = 0; id < indexDistributions.size(); ++id )
        {
            dmemo::DistributionPtr targetDist = targetDistributions[sd];
            dmemo::DistributionPtr indexDist = indexDistributions[id];

            source.redistribute( indexDist );
            index.redistribute( indexDist );

            DenseVector<ValueType> target( targetDist );
            target = 0;

            if ( targetDist->isReplicated() != indexDist->isReplicated() )
            {
                // either both are replicated or both distributed

                BOOST_CHECK_THROW(
                {
                    target.scatter( index, source, common::BinaryOp::ADD );
                }, common::Exception );

                continue;
            }

            target.scatter( index, source, common::BinaryOp::ADD );

            hmemo::ReadAccess<ValueType> rTarget( target.getLocalValues() );

            for ( IndexType i = 0; i < n; ++i )
            {
                IndexType localIndex = target.getDistribution().global2local( i );

                if ( localIndex != nIndex )
                {
                    BOOST_CHECK_MESSAGE( rTarget[localIndex] == targetValues[i],
                                         *comm << ": targetLocal[" << localIndex << "] = " << rTarget[localIndex]
                                         << " must be equal to targetValues[" << i << "] = " << targetValues[i] );
                }
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
