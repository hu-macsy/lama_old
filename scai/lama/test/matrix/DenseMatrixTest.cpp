/**
 * @file test/matrix/DenseMatrixTest.cpp
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
 * @brief Test routines for specific methods/constructors of DenseMatrix
 * @author Thomas Brandes
 * @date 24.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/test/TestMacros.hpp>

#include <scai/lama/test/storage/Storages.hpp>
#include <scai/lama/test/storage/TestStorages.hpp>

#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/fft.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;

using boost::test_tools::per_element;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DenseMatrixTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.DenseMatrixTest" );

/* ------------------------------------------------------------------------- */

/** For the matrix tests here it is sufficient to take only one of the possible value types. */

typedef DefaultReal ValueType;

/* ------------------------------------------------------------------------- */

typedef boost::mpl::list < CSRSparseMatrix<ValueType>,
        ELLSparseMatrix<ValueType>,
        JDSSparseMatrix<ValueType>,
        DIASparseMatrix<ValueType>,
        COOSparseMatrix<ValueType>,
        DenseMatrix<ValueType>
        > SparseMatrixTypes;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( matrixTimesVectorN, MatrixType, SparseMatrixTypes )
{
    return;

    // Test vector = Matrix * vector, where vector stands for multiple vectors
    // i.e. vector is a dense matrix
    // Note: not yet available for distributed matrix

    auto comm = dmemo::Communicator::getCommunicatorPtr();

    const IndexType n = 20;  // size of the square matrix

    auto matrix = zero<MatrixType>( 20, 20 );

    MatrixCreator::fillRandom( matrix, 0.2f );

    const IndexType k = 3;  // number of vectors

    DenseMatrix<ValueType> vectorK( n, k );

    MatrixCreator::fillRandom( vectorK, 1.0f );

    DenseMatrix<ValueType> vectorK1;

    vectorK1 = matrix * vectorK;  // calls matrixTimesVectorN

    BOOST_CHECK_EQUAL( n, vectorK1.getNumRows() );
    BOOST_CHECK_EQUAL( k, vectorK1.getNumColumns() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( invertTest, ValueType, scai_numeric_test_types )
{
    // inverse of a matrix not yet for distributed matrices

    const IndexType n = 20;  // size of the square matrix

    auto matrix = zero<DenseMatrix<ValueType>>( n, n );

    MatrixCreator::fillRandom( matrix, 1.0f );

    DenseMatrix<ValueType> invMatrix;

    invMatrix.invert( matrix );

    auto e1 = identity<DenseMatrix<ValueType>>( n );
    
    auto e2 = eval<DenseMatrix<ValueType>>( matrix * invMatrix );

    BOOST_CHECK( e1.maxDiffNorm( e2 ) < 1e-3 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( columnDistributionTest )
{
    // Motivation: dedicated test to check that split/join of distributed column
    //             data works which is typique for dense matrices.

    typedef DefaultReal ValueType;

    // In this method we just test that split/join of columns works fine

    const IndexType numRows = 5;
    const IndexType numColumns = 12;
 
    dmemo::DistributionPtr repRowDist( new dmemo::NoDistribution( numRows ) );
    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( numColumns ) );

    DenseStorage<ValueType> denseData( numRows, numColumns );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            denseData.setValue( i, j, ValueType( i * 100 + j ) );
        }
    }

    dmemo::TestDistributions dists( numColumns );

    for ( size_t k = 0; k < dists.size(); ++k )
    {
        dmemo::DistributionPtr dist = dists[k];

        DenseMatrix<ValueType> denseM( denseData );

        // split the column data

        denseM.redistribute( repRowDist, dist );

        // join the column data

        denseM.redistribute( repRowDist, repColDist );

        DenseStorage<ValueType> dense1 = denseM.getLocalStorage();

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType j = 0; j < numColumns; ++j )
            {
                BOOST_CHECK_EQUAL( denseData.getValue( i, j ), ValueType( i * 100 + j ) );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( dense2SparseTest )
{
    // Motivation: dedicated test to check that split/join of distributed column
    //             data works which is typique for dense matrices.

    typedef DefaultReal ValueType;

    // In this method we just test that split/join of columns works fine

    const IndexType numRows = 5;
    const IndexType numColumns = 12;
 
    dmemo::DistributionPtr repRowDist( new dmemo::NoDistribution( numRows ) );
    dmemo::DistributionPtr repColDist( new dmemo::NoDistribution( numColumns ) );

    DenseStorage<ValueType> denseData( numRows, numColumns );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            denseData.setValue( i, j, ValueType( i * 100 + j ) );
        }
    }

    dmemo::TestDistributions dists( numColumns );

    for ( size_t k = 0; k < dists.size(); ++k )
    {
        dmemo::DistributionPtr dist = dists[k];

        DenseMatrix<ValueType> denseM( denseData );

        denseM.redistribute( repRowDist, dist );

        CSRSparseMatrix<ValueType> sparseM( denseM );

        BOOST_REQUIRE_EQUAL( denseM.getColDistribution(), sparseM.getColDistribution() );
        BOOST_REQUIRE_EQUAL( denseM.getRowDistribution(), sparseM.getRowDistribution() );

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType j = 0; j < numColumns; ++j )
            {
                BOOST_CHECK_EQUAL( sparseM.getValue( i, j ), ValueType( i * 100 + j ) );
            }
        }

        sparseM.redistribute( repRowDist, repColDist );

        CSRStorage<ValueType> csrData = sparseM.getLocalStorage();

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType j = 0; j < numColumns; ++j )
            {
                BOOST_CHECK_EQUAL( csrData.getValue( i, j ), ValueType( i * 100 + j ) );
            }
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( resizeTest )
{
    typedef SCAI_TEST_TYPE ValueType;     // value type does not matter for resize

    using namespace hmemo;

    const IndexType m = 8;
    const IndexType n = 13;

    std::srand( 1311 );   // same random numbers on all processors

    DenseMatrix<ValueType> matrix1;

    HArray<ValueType> data1( m * n );
    utilskernel::HArrayUtils::setRandom( data1, 1 );
    DenseStorage<ValueType> repMatrix( m, n, data1 );

    dmemo::DistributionPtr rowDistBig = std::make_shared<dmemo::BlockDistribution>( m + 1 );
    dmemo::DistributionPtr rowDistSmall = std::make_shared<dmemo::CyclicDistribution>( m - 1, 2 );
    dmemo::DistributionPtr colDistBig = std::make_shared<dmemo::NoDistribution>( n + 1 );
    dmemo::DistributionPtr colDistSmall = std::make_shared<dmemo::NoDistribution>( n - 1 );

    matrix1.assign( repMatrix );

    matrix1.resize( rowDistBig, colDistSmall );

    ValueType v1 = matrix1.getValue( 2, 5 );
    ValueType v2 = repMatrix.getValue( 2, 5 );
    BOOST_CHECK_EQUAL( v1, v2 );

    v1 = matrix1.getValue( m, 3 );   // row m has been filled with zero
    BOOST_CHECK_EQUAL( v1, 0 );
 
    matrix1.resize( rowDistSmall, colDistBig );

    v1 = matrix1.getValue( 4, n );   // column n has been filled with zero
    BOOST_CHECK_EQUAL( v1, 0 );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( matrixMultTest )
{
    typedef SCAI_TEST_TYPE ValueType;     // value type does not matter for resize

    using namespace hmemo;

    const IndexType m = 5;
    const IndexType k = 4;
    const IndexType n = 3;

    HArray<ValueType> dataA( { 1, 2, 3, 4,
                               0, 1, 0, -2,
                               2, 1, -1, 1,
                               -3, 1, -1, 1,
                               2, 0, 1, -1 } );

    HArray<ValueType> dataB( { 2, 1, 0,
                               1, 0, -1,
                               2, 0, 3,
                               1, 1, 2  } );

    HArray<ValueType> dataC( { 2, 1, 0,
                               1, 0, -1,
                               3, 1, 2,
                               2, 0, 3,
                               1, 1, 2  } );

    DenseStorage<ValueType> storageA( m, k, dataA );
    DenseStorage<ValueType> storageB( k, n, dataB );
    DenseStorage<ValueType> storageC( m, n, dataC );

    ValueType alpha = 2;
    ValueType beta  = -1;

    // serial computation to compare later with distributed computation

    DenseStorage<ValueType> expected;
    expected.matrixTimesMatrix( alpha, storageA, storageB, beta, storageC );
   
    dmemo::TestDistributions distributions( m );

    for ( size_t id = 0; id < distributions.size(); ++id )
    {
        auto distM = distributions[m];

        SCAI_LOG_DEBUG( logger, "matrixTimesMatrixDense: dist for result, a = " << *distM );

        auto repM = std::make_shared<dmemo::NoDistribution>( m );
        auto repN = std::make_shared<dmemo::NoDistribution>( n );
        auto repK = std::make_shared<dmemo::NoDistribution>( k );

        auto matrixA = distribute<DenseMatrix<ValueType>>( storageA, distM, repK );
        auto matrixB = distribute<DenseMatrix<ValueType>>( storageB, repK, repN );
        auto matrixC = distribute<DenseMatrix<ValueType>>( storageC, distM, repN );

        SCAI_LOG_DEBUG( logger, "matrixTimesMatrixDense: a = " << matrixA << ", b = " << matrixB << ", c = " << matrixC )

        DenseMatrix<ValueType> matrixResult;

        matrixResult = alpha * matrixA * matrixB + beta * matrixC;

        SCAI_LOG_DEBUG( logger, "matrixResult: res = " << matrixResult )

        matrixResult.redistribute( repM, repN );

        BOOST_TEST( hostReadAccess( matrixResult.getLocalStorage().getValues() ) == 
                    hostReadAccess( expected.getValues() ), per_element() );
    }
}

/* ------------------------------------------------------------------------- */

#ifdef SCAI_COMPLEX_SUPPORTED

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fftTestRow, ValueType, scai_fft_test_types )
{
    using hmemo::HArray;
    using hmemo::hostReadAccess;

    typedef common::Complex<RealType<ValueType>> FFTType;

    HArray<FFTType> input( { 0.5, 1.0, 0.0, 0.0,
                             1,   2,   3,   4   } );

    const IndexType numRows = 2;
    const IndexType numColumns = 4;

    DenseMatrix<FFTType> matrix( DenseStorage<FFTType>( numRows, numColumns, input ) );

    fft( matrix, 1 );    // fft along columns, is parallel for each row

    HArray<FFTType> expResult( { 1.5, FFTType( 0.5, -1 ), -0.5, FFTType( 0.5, 1 ),
                                 10, FFTType( -2, 2 ), -2, FFTType( -2, -2 ) } );

    BOOST_TEST( hostReadAccess( matrix.getLocalStorage().getValues() ) == hostReadAccess( expResult ), boost::test_tools::per_element() );

    ifft( matrix, 1 );   // inverse fft
 
    matrix *= ValueType( 1 ) / ValueType( numColumns );

    BOOST_TEST( hostReadAccess( matrix.getLocalStorage().getValues() ) == hostReadAccess( input ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fftTestCol, ValueType, scai_fft_test_types )
{
    using hmemo::HArray;
    using hmemo::hostReadAccess;

    typedef common::Complex<RealType<ValueType>> FFTType;

    HArray<FFTType> input( { 0.5,  1, 
                             1.0,  2, 
                             0.0,  3,
                             0.0,  4 } );

    const IndexType numRows = 4;
    const IndexType numColumns = 2;

    DenseMatrix<FFTType> matrix( DenseStorage<FFTType>( numRows, numColumns, input ) );

    fft( matrix, 0 );    // fft along rows, is parallel for each column

    HArray<FFTType> expResult( { 1.5,                10,
                                 FFTType( 0.5, -1 ), FFTType( -2, 2 ),
                                 -0.5,               -2, 
                                 FFTType( 0.5, 1 ),  FFTType( -2, -2 ) } );

    BOOST_TEST( hostReadAccess( matrix.getLocalStorage().getValues() ) == hostReadAccess( expResult ), boost::test_tools::per_element() );

    ifft( matrix, 0 );   // inverse fft
 
    matrix *= ValueType( 1 ) / ValueType( numRows );

    BOOST_TEST( hostReadAccess( matrix.getLocalStorage().getValues() ) == hostReadAccess( input ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( allFFTTest, ValueType, scai_fft_test_types )
{
    typedef common::Complex<RealType<ValueType>> FFTType;

    const IndexType M = 13; 
    const IndexType N = 19;

    // generate random vector

    hmemo::HArray<ValueType> randomValues( M * N );
    utilskernel::HArrayUtils::setRandom( randomValues, 1 );

    DenseStorage<ValueType> storage( M, N, randomValues );
    DenseMatrix<ValueType> x( storage );

    // fill up to M2 x N2, dimensions that are the next power 2 values

    const IndexType M2 = 1 << common::Math::nextpow2( M );
    const IndexType N2 = 1 << common::Math::nextpow2( N );

    auto rowDist = std::make_shared<dmemo::BlockDistribution>( M2 );
    auto colDist = std::make_shared<dmemo::BlockDistribution>( N2 );

    // convert to complex matrix, fill it up and distribute it

    DenseMatrix<FFTType> y;
    y = cast<FFTType>( x );

    SCAI_LOG_INFO( logger, "resize y = " << y << " with row dist = " << *rowDist << ", col dist = " << *colDist )

    y.resize( rowDist, colDist );   // fill up and distribute

    fft( y );        // FFT forward
    ifft( y );       // FFT backward

    y *= ValueType( 1 ) / ValueType( M2 * N2 );

    // resize back to original data

    y.resize( x.getRowDistributionPtr(), x.getColDistributionPtr() );

    // divide by M * N after fft - ifft to get the original result

    auto x1 = convert<DenseMatrix<ValueType>>( y );

    RealType<ValueType> eps = common::TypeTraits<ValueType>::small();

    if ( false )
    {
        // helpful for debugging, but we have only close equality here

        const hmemo::HArray<ValueType>& denseData = x1.getLocalStorage().getValues();
        BOOST_TEST( hostReadAccess( denseData ) == hostReadAccess( randomValues ), boost::test_tools::per_element() );
    }

    BOOST_CHECK( x.maxDiffNorm( x1 ) < eps );
}

/* ------------------------------------------------------------------------- */

#endif

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
