/**
 * @file test/matrix/DenseMatrixTest.cpp
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
#include <scai/lama/expression/all.hpp>
#include <scai/lama/fft.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DenseMatrixTest )

/* ------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.SparseMatrixTest" );

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

#ifdef SCAI_COMPLEX_SUPPORTED

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fftTestRow1, ValueType, scai_fft_test_types )
{
    using hmemo::HArray;
    using hmemo::hostReadAccess;

    typedef common::Complex<RealType<ValueType>> FFTType;

    HArray<ValueType> input( { 0.5, 1.0 } );

    DenseMatrix<ValueType> x( DenseStorage<ValueType>( 1, 2, input ) );

    HArray<FFTType> result( { 1.5, FFTType( 0.5, -1 ), -0.5, FFTType( 0.5, 1 ) } );

    DenseMatrix<FFTType> res;

    const IndexType dim = 1;

    fft( res, x, dim, 4 );    // fft for each row

    SCAI_LOG_DEBUG( logger, "Result of fft : " << res )

    BOOST_CHECK_EQUAL( res.getNumRows(), 1 );
    BOOST_CHECK_EQUAL( res.getNumColumns(), 4 );

    BOOST_TEST( hostReadAccess( res.getLocalStorage().getValues() ) == hostReadAccess( result ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fftTestRow2, ValueType, scai_fft_test_types )
{
    using hmemo::HArray;
    using hmemo::hostReadAccess;

    typedef common::Complex<RealType<ValueType>> FFTType;

    HArray<ValueType> input( { 1.0, 1.5 } );

    DenseMatrix<ValueType> x( DenseStorage<ValueType>( 1, 2, input ) );

    HArray<FFTType> result( { 2.5, FFTType( 1, -1.5 ), -0.5, FFTType( 1, 1.5 )  } );

    DenseMatrix<FFTType> res;

    const IndexType dim = 1;

    fft( res, x, dim, 4 );    // fft for each row

    SCAI_LOG_ERROR( logger, "Result of fft : " << res )

    BOOST_CHECK_EQUAL( res.getNumRows(), 1 );
    BOOST_CHECK_EQUAL( res.getNumColumns(), 4 );

    BOOST_TEST( hostReadAccess( res.getLocalStorage().getValues() ) == hostReadAccess( result ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fftTestRowArray, ValueType, scai_fft_test_types )
{
    using hmemo::HArray;
    using hmemo::hostReadAccess;

    typedef common::Complex<RealType<ValueType>> FFTType;

    HArray<ValueType> input( { 0.5, 1.0, 1.0, 1.5 } );

    HArray<FFTType> result( { 1.5, FFTType( 0.5, -1 ), -0.5, FFTType( 0.5, 1 ) ,
                              2.5, FFTType( 1, -1.5 ), -0.5, FFTType( 1, 1.5 )  } );

    hmemo::HArray<FFTType> res;

    utilskernel::FFTUtils::fft_many( res, input, 2, 4, 1 );

    SCAI_LOG_ERROR( logger, "Result of fft : " << res )

    BOOST_TEST( hostReadAccess( res ) == hostReadAccess( result ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fftTestRowMatrix, ValueType, scai_fft_test_types )
{
    using hmemo::HArray;
    using hmemo::hostReadAccess;

    typedef common::Complex<RealType<ValueType>> FFTType;

    HArray<ValueType> input( { 0.5, 1.0, 1.0, 1.5 } );

    DenseMatrix<ValueType> x( DenseStorage<ValueType>( 2, 2, input ) );

    HArray<FFTType> result( { 1.5, FFTType( 0.5, -1 ), -0.5, FFTType( 0.5, 1 ) ,
                              2.5, FFTType( 1, -1.5 ), -0.5, FFTType( 1, 1.5 )  } );

    DenseMatrix<FFTType> res;

    const IndexType dim = 1;

    fft( res, x, dim, 4 );    // fft for each row

    SCAI_LOG_ERROR( logger, "Result of fft : " << res )

    BOOST_CHECK_EQUAL( res.getNumRows(), 2 );
    BOOST_CHECK_EQUAL( res.getNumColumns(), 4 );

    BOOST_TEST( hostReadAccess( res.getLocalStorage().getValues() ) == hostReadAccess( result ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( fftTestCol, ValueType, scai_fft_test_types )
{
    return;

    using hmemo::HArray;
    using hmemo::hostReadAccess;

    typedef common::Complex<RealType<ValueType>> FFTType;

    HArray<ValueType> input( { 0.5, 1.0, 1.0, 1.5 } );

    DenseMatrix<ValueType> x( DenseStorage<ValueType>( 2, 2, input ) );

    HArray<FFTType> result( { 1.5, 2.5, FFTType( 0.5, -1 ), FFTType( 1, -1.5 ) ,
                              -0.5, -0.5, FFTType( 0.5, 1 ), FFTType( 1, 1.5 )  } );

    DenseMatrix<FFTType> res;

    const IndexType dim = 0;

    fft( res, x, dim, 4 );    // fft for each column

    SCAI_LOG_ERROR( logger, "Result of fft : " << res )

    BOOST_CHECK_EQUAL( res.getNumRows(), 4 );
    BOOST_CHECK_EQUAL( res.getNumColumns(), 2 );

    BOOST_TEST( hostReadAccess( res.getLocalStorage().getValues() ) == hostReadAccess( result ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ifftTestRow, ValueType, scai_fft_test_types )
{
    return;

    using hmemo::HArray;
    using hmemo::hostReadAccess;

    typedef common::Complex<RealType<ValueType>> FFTType;

    HArray<ValueType> input( { 0.5, 1.0, 2.0, 2.5 } );

    DenseMatrix<ValueType> x( DenseStorage<ValueType>( 2, 2, input ) );

    HArray<FFTType> result( { 1.5, FFTType( 0.5, 1.0 ), -0.5, FFTType( 0.5, -1.0 ) ,
                              4.5, FFTType( 2.0, 2.5 ), -0.5, FFTType( 2.0, -2.5 )  } );

    DenseMatrix<FFTType> res;

    const IndexType dim = 1;

    ifft( res, x, dim, 4 );    // ifft for each row

    SCAI_LOG_ERROR( logger, "Result of fft : " << res )

    BOOST_CHECK_EQUAL( res.getNumRows(), 2 );
    BOOST_CHECK_EQUAL( res.getNumColumns(), 4 );

    // BOOST_TEST( hostReadAccess( res.getLocalStorage().getValues() ) == hostReadAccess( result ), boost::test_tools::per_element() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( ifftTestCol, ValueType, scai_fft_test_types )
{
    return;

    using hmemo::HArray;
    using hmemo::hostReadAccess;

    typedef common::Complex<RealType<ValueType>> FFTType;

    HArray<ValueType> input( { 0.5, 1.0, 1.0, 1.5 } );

    DenseMatrix<ValueType> x( DenseStorage<ValueType>( 2, 2, input ) );

    HArray<FFTType> result( { 2.5, 3.5, 
                              FFTType( 0.5, 2.0 ), FFTType( 1.0, 2.5 ) ,
                              -1.5, -1.5, 
                              FFTType( 0.5, -2.0 ), FFTType( 1.0, -2.5 )  } );

    DenseMatrix<FFTType> res;

    const IndexType dim = 0;

    ifft( res, x, dim, 4 );    // ifft for each column

    SCAI_LOG_ERROR( logger, "Result of fft : " << res )

    BOOST_CHECK_EQUAL( res.getNumRows(), 4 );
    BOOST_CHECK_EQUAL( res.getNumColumns(), 2 );

    BOOST_TEST( hostReadAccess( res.getLocalStorage().getValues() ) == hostReadAccess( result ), boost::test_tools::per_element() );
}

#endif

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
