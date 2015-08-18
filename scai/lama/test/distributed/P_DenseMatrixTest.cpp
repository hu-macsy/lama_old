/**
 * @file P_DenseMatrixTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Contains the implementation of the class P_DenseMatrixTest.
 * @author: Alexander Büchel, Michael Drost
 * @date 25.04.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/common/shared_ptr.hpp>
#include <scai/common/unique_ptr.hpp>

#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/distribution/BlockDistribution.hpp>
#include <scai/lama/distribution/GeneralDistribution.hpp>
#include <scai/lama/distribution/CyclicDistribution.hpp>
#include <scai/lama/distribution/NoDistribution.hpp>

#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>

#include <scai/lama/LAMAInterfaceRegistry.hpp>

#include <test/TestSparseMatrices.hpp>
#include <test/TestMacros.hpp>

using namespace scai::common;
using namespace scai::lama;
using namespace scai::hmemo;

typedef boost::mpl::list<double, float> test_types;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct P_DenseMatrixTestConfig
{
    P_DenseMatrixTestConfig()
    {
        comm = Communicator::get( "MPI" );
    }

    ~P_DenseMatrixTestConfig()
    {
        comm = CommunicatorPtr();
    }
};

BOOST_FIXTURE_TEST_SUITE( P_DenseMatrixTest, P_DenseMatrixTestConfig );

SCAI_LOG_DEF_LOGGER( logger, "Test.P_DenseMatrixTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( buildTest, ValueType, test_types )
{
    PartitionId size = comm->getSize();
    int numRows = 3 * size;
    int numCols = 5 * size;
    scoped_array<ValueType> values( new ValueType[ numRows * numCols ] );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numCols; ++j )
        {
            ValueType value = static_cast<ValueType> ( i * numCols + j + 1.0 );
            values[ i * numCols + j ] = value;
        }
    }

    DenseMatrix<ValueType> repM;
    repM.setRawDenseData( numRows, numCols, values.get() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numCols; ++j )
        {
            Scalar value = repM.getValue( i , j );
            Scalar expectedvalue = Scalar( static_cast<ValueType>( values[i * numCols + j] ) );
            SCAI_CHECK_SCALAR_SMALL( value - expectedvalue, ValueType, eps<ValueType>() );
        }
    }

    shared_ptr<Distribution> dist( new BlockDistribution( numRows, comm ) );
    shared_ptr<Distribution> distCol( new BlockDistribution( numCols, comm ) );
    DenseMatrix<ValueType> distM( repM, dist, distCol );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numCols; ++j )
        {
            Scalar value = distM.getValue( i, j );
            Scalar expectedvalue = repM.getValue( i, j );
            SCAI_CHECK_SCALAR_SMALL( value - expectedvalue , ValueType, eps<ValueType>() );
        }
    }
}

/* --------------------------------------------------------------------- */

void cyclicMultVectorTest( const IndexType chunkSize, const IndexType n )
{
    scoped_array<float> values( new float[n * n] );

    for ( IndexType i = 0; i < n; ++i )
    {
        for ( IndexType j = 0; j < n; ++j )
        {
            values[i * n + j] = 1.0f / n;
        }
    }

    DenseMatrix<float> repMatrix;
    repMatrix.setRawDenseData( n, n, values.get() );
    scai::common::shared_ptr<Distribution> dist( new CyclicDistribution( n, chunkSize, comm ) );
    DenseMatrix<float> distMatrix( repMatrix, dist, dist );
    distMatrix.setCommunicationKind( Matrix::ASYNCHRONOUS );
    DenseVector<float> x( dist, 1.0f );
    DenseVector<float> res( dist, 0.0f );
    res = distMatrix * x;

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_CLOSE( 1.0f, res.getValue( i ).getValue<float>(), 1.0 / ( 10.0 * n ) );
    }

    res = 0.0;

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( 0.0f, res.getValue( i ) );
    }

    distMatrix.setCommunicationKind( Matrix::SYNCHRONOUS );
    res = distMatrix * x;

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_CLOSE( 1.0f, res.getValue( i ).getValue<float>(), 1.0 / ( 10.0 * n ) );
    }
}

/* -------------------------------------------------------------------------- */

static Distribution* createDistribution( const IndexType n, CommunicatorPtr comm, int kind )
{
    Distribution* dist;

    if ( kind == 0 )
    {
        dist = new BlockDistribution( n, comm );
    }
    else if ( kind == 1 )
    {
        dist = new NoDistribution( n );
    }
    else if ( kind == 2 )
    {
        std::vector<IndexType> localIndexes;
        IndexType size = comm->getSize();
        IndexType rank = comm->getRank();

        for ( int k = 0; k < n; k++ )
        {
            if ( k % size == rank )
            {
                localIndexes.push_back( k );
            }
        }

        dist = new GeneralDistribution( n, localIndexes, comm );
    }
    else if ( kind == 3 )
    {
        IndexType chunkSize = comm->getSize();
        dist = new CyclicDistribution( n, chunkSize, comm );
    }
    else
    {
        COMMON_THROWEXCEPTION( "kind = " << kind << " unsupported here" )
    }

    return dist;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( maxDiffNormTest, ValueType, test_types )
{
    DenseMatrix<ValueType> n6m4SMatrix  = TestSparseMatrices::n6m4MatrixE1<ValueType>();
    DenseMatrix<ValueType> n6m4SMatrix1 = 2.5 * n6m4SMatrix;
    shared_ptr<Distribution> dist( new BlockDistribution( n6m4SMatrix1.getNumRows(), comm ) );
    shared_ptr<Distribution> distCol( new BlockDistribution( n6m4SMatrix1.getNumColumns(), comm ) );
    DenseMatrix<ValueType> distM( n6m4SMatrix1, dist, distCol );
    // check matrices with different distributions
    Scalar maxDiffNorm = n6m4SMatrix.maxDiffNorm( distM );
    BOOST_CHECK_EQUAL( maxDiffNorm.getValue<ValueType>(), 13.5 );
    maxDiffNorm = distM.maxDiffNorm( n6m4SMatrix );
    BOOST_CHECK_EQUAL( maxDiffNorm.getValue<ValueType>(), 13.5 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multVectorTest, ValueType, test_types )
{
    int numRows = 20;
    int numCols = 31;
    // definition of raw data for setup and comparison
    scoped_array<ValueType> valuesA( new ValueType[numRows * numCols] );
    scoped_array<ValueType> valuesX( new ValueType[numCols] );
    scoped_array<ValueType> valuesY( new ValueType[numRows] );

    // intialise data for the matrix

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numCols; ++j )
        {
            ValueType value = static_cast<ValueType> ( 100.0 - ::fabs( 2.0 * i - j ) );
            valuesA[i * numCols + j] = value;
        }
    }

    // initialize the vector x

    for ( IndexType j = 0; j < numCols; ++j )
    {
        valuesX[j] = static_cast<ValueType> ( 1.2 * ( j + 1 ) );
    }

// compute Y = A * X for comparison

    for ( IndexType i = 0; i < numRows; ++i )
    {
        valuesY[i] = 0.0;

        for ( IndexType j = 0; j < numCols; ++j )
        {
            valuesY[i] += static_cast<ValueType> ( valuesA[i * numCols + j] * valuesX[j] );
        }
    }

    // construct replicated matrix A and vector X to be redistributed
    DenseMatrix<ValueType> rA;
    rA.setRawDenseData( numRows, numCols, valuesA.get() );
    DenseVector<ValueType> rX( numCols, valuesX.get() );
    // try different distributions for rows and colums
    DistributionPtr rowDist;
    DistributionPtr colDist;

    for ( int i = 0; i < 4; i++ )
    {
        rowDist = DistributionPtr( createDistribution( numRows, comm, i ) );

        for ( int j = 0; j < 4; j++ )
        {
            colDist = DistributionPtr( createDistribution( numCols, comm, j ) );

            for ( int k = 0; k < 2; k++ )
            {
                // redistribute A, X, setup result vector Y
                DenseMatrix<ValueType> A( rA, rowDist, colDist );
                DenseVector<ValueType> X( rX, colDist );

                if ( k == 0 )
                {
                    A.setCommunicationKind( Matrix::SYNCHRONOUS );
                }
                else
                {
                    A.setCommunicationKind( Matrix::ASYNCHRONOUS );
                }

                SCAI_LOG_INFO( logger, "mult matrix A = " << A << " with vector X = " << X );
                DenseVector<ValueType> result ( A * X );
                BOOST_REQUIRE_EQUAL( result.size(), numRows );
                BOOST_REQUIRE_EQUAL( result.getDistribution(), *rowDist );

                // compare the result vector Y with the precomputed results valuesY

                for ( IndexType m = 0; m < numRows; ++m )
                {
                    Scalar value = Scalar( valuesY[m] );
                    Scalar expectedvalue = result.getValue( m );
                    //1e-1 is used because there are 8 of 1280 cases which fails with eps<1e-1
                    SCAI_CHECK_SCALAR_SMALL( value - expectedvalue , ValueType, 1e-1 );
                }
            }
        }
    }

    //Tests for Cyclic Dist speciliatzion in case of a square matrix
    {
        PartitionId size = comm->getSize();
        IndexType chunkSize = 7;
        IndexType n = 2 * size * chunkSize;
        cyclicMultVectorTest( chunkSize, n );
        n = 3 * size * chunkSize + chunkSize - 2;
        cyclicMultVectorTest( chunkSize, n );
        //Calculate the size so that each process gets more than one chunk and that not
        //all chunks have the same numbers of chunks and that we have at least von chunk that is not
        // square
        n = 3 * size * chunkSize + size / 2 * chunkSize + chunkSize - 2;
        cyclicMultVectorTest( chunkSize, n );
        //not all process get a chunk
        n = size / 2 * chunkSize + chunkSize - 2;
        cyclicMultVectorTest( chunkSize, n );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildSquareTest )
{
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();
    int numRows = 3 * size;
    int numCols = 3 * size;
    scoped_array<double> values( new double[numRows * numCols] );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numCols; ++j )
        {
            double value = 0.0;

            if ( j == i || j + size == i || j - size == i || j + 2 * size == i || j - 2 * size == i
                    || j + ( numRows - 1 ) == i || j - ( numRows - 1 ) == i )
            {
                value = 1000.0 * ( i + 1 ) + ( j + 1 );
            }

            values[i * numCols + j] = value;
        }
    }

    DenseMatrix<double> repM;
    repM.setRawDenseData( numRows, numCols, values.get() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numCols; ++j )
        {
            BOOST_CHECK_CLOSE( values[i * numCols + j], repM.getValue( i, j ).getValue<double>(), 1e-16 );
        }
    }

    std::vector<IndexType> localIndexes;

    // take a cyclic(1) distribution as general distribution

    for ( IndexType i = 0; i < numRows; ++i )
    {
        if ( i % size == rank )
        {
            localIndexes.push_back( i );
        }
    }

    shared_ptr<Distribution> dist( new GeneralDistribution( numRows, localIndexes, comm ) );
    DenseMatrix<double> distM( repM, dist, dist );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numCols; ++j )
        {
            BOOST_CHECK_CLOSE( repM.getValue( i, j ).getValue<double>(), distM.getValue( i, j ).getValue<double>(),
                               1e-16 );
        }
    }
}

/* --------------------------------------------------------------------- */

void cyclicDistTestImpl( const IndexType chunkSize, const IndexType n )
{
    scoped_array<float> values( new float[n * n] );

    for ( IndexType i = 0; i < n; ++i )
    {
        for ( IndexType j = 0; j < n; ++j )
        {
            values[i * n + j] = static_cast<float>( i * n + j );
        }
    }

    DenseMatrix<float> repMatrix;
    repMatrix.setRawDenseData( n, n, values.get() );
    scai::common::shared_ptr<Distribution> dist( new CyclicDistribution( n, chunkSize, comm ) );
    DenseMatrix<float> distMatrix( repMatrix, dist, dist );

    for ( IndexType i = 0; i < n; i++ )
    {
        for ( IndexType j = 0; j < n; j++ )
        {
            Scalar expectedvalue = repMatrix.getValue( i, j );
            Scalar value = distMatrix.getValue( i, j );
            Scalar diff = expectedvalue - value;
            SCAI_CHECK_SCALAR_SMALL( diff, float, 1E-8 );
        }
    }
}

BOOST_AUTO_TEST_CASE( cyclicDistTest )
{
    PartitionId size = comm->getSize();
    const IndexType chunkSize = 7;
    IndexType n = 2 * size * chunkSize;
    cyclicDistTestImpl( chunkSize, n );
    n = 3 * size * chunkSize + chunkSize - 2;
    cyclicDistTestImpl( chunkSize, n );
    //Calculate the size so that each process gets more than one chunk and that not
    //all chunks have the same numbers of chunks and that we have at least von chunk that is not
    // square
    n = 3 * size * chunkSize + size / 2 * chunkSize + chunkSize - 2;
    cyclicDistTestImpl( chunkSize, n );
    //not all process get a chunk
    n = size / 2 * chunkSize + chunkSize - 2;
    cyclicDistTestImpl( chunkSize, n );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
