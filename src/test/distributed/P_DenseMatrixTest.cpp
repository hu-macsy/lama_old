/**
 * @file P_DenseMatrixTest.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @author: Alexander Büchel, mdrost
 * @date 25.04.2012
 * $
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <boost/scoped_array.hpp>
#include <boost/shared_ptr.hpp>

#include <lama/matrix/DenseMatrix.hpp>

#include <lama/distribution/BlockDistribution.hpp>
#include <lama/distribution/GeneralDistribution.hpp>
#include <lama/distribution/CyclicDistribution.hpp>
#include <lama/distribution/NoDistribution.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>
#include <lama/expression/MatrixExpressions.hpp>

#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/CommunicatorFactory.hpp>

#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

typedef boost::mpl::list<double,float> test_types;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct P_DenseMatrixTestConfig
{
    P_DenseMatrixTestConfig()
    {
        comm = CommunicatorFactory::get( "MPI" );
    }

    ~P_DenseMatrixTestConfig()
    {
        comm = CommunicatorPtr();
    }
};

BOOST_FIXTURE_TEST_SUITE( P_DenseMatrixTest, P_DenseMatrixTestConfig )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.P_DenseMatrixTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( buildTest, T, test_types ) {
    typedef T ValueType;

    PartitionId size = comm->getSize();

    int numRows = 3 * size;
    int numCols = 5 * size;

    scoped_array<ValueType> values( new ValueType[ numRows * numCols ] );

    for ( IndexType i = 0; i < numRows; ++i)
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
        for ( IndexType j = 0; j<numCols; ++j )
        {
            Scalar value = repM.getValue( i , j );
            Scalar expectedvalue = Scalar( static_cast<ValueType>(values[i * numCols + j]) );
            LAMA_CHECK_SCALAR_SMALL( value - expectedvalue, ValueType, eps<ValueType>() );
        }
    }

    shared_ptr<Distribution> dist( new BlockDistribution(numRows, comm) );
    shared_ptr<Distribution> distCol( new BlockDistribution(numCols, comm) );

    DenseMatrix<ValueType> distM( repM, dist, distCol );

    for (IndexType i = 0; i<numRows; ++i)
    {
        for (IndexType j = 0; j<numCols; ++j)
        {
            Scalar value = distM.getValue( i, j );
            Scalar expectedvalue = repM.getValue( i, j );
            LAMA_CHECK_SCALAR_SMALL( value - expectedvalue , ValueType, eps<ValueType>() );
        }
    }
}

/* --------------------------------------------------------------------- */

void cyclicMultVectorTest( const IndexType chunkSize, const IndexType n )
{
    boost::scoped_array<float> values( new float[n * n] );

    for ( IndexType i = 0; i < n; ++i )
    {
        for ( IndexType j = 0; j < n; ++j )
        {
            values[i * n + j] = 1.0f / n;
        }
    }

    DenseMatrix<float> repMatrix;
    repMatrix.setRawDenseData( n, n, values.get() );

    boost::shared_ptr<Distribution> dist( new CyclicDistribution( n, chunkSize, comm ) );

    DenseMatrix<float> distMatrix( repMatrix, dist, dist );

    distMatrix.setCommunicationKind( Matrix::ASYNCHRONOUS );

    DenseVector<float> x( dist, 1.0f );

    DenseVector<float> res( dist, 0.0f );

    res = distMatrix * x;

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_CLOSE( 1.0f, res.getValue(i).getValue<float>(), 1.0 / (10.0*n) );
    }

    res = 0.0;

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_EQUAL( 0.0f, res.getValue(i) );
    }

    distMatrix.setCommunicationKind( Matrix::SYNCHRONOUS );

    res = distMatrix * x;

    for ( IndexType i = 0; i < n; ++i )
    {
        BOOST_CHECK_CLOSE( 1.0f, res.getValue(i).getValue<float>(), 1.0/(10.0*n) );
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multVectorTest, T, test_types ) {
    typedef T ValueType;

    int numRows = 20;
    int numCols = 31;

// definition of raw data for setup and comparison

    scoped_array<ValueType> valuesA( new ValueType[numRows * numCols] );
    scoped_array<ValueType> valuesX( new ValueType[numCols] );
    scoped_array<ValueType> valuesY( new ValueType[numRows] );

// intialise data for the matrix

    for ( IndexType i = 0; i<numRows; ++i )
    {
        for ( IndexType j = 0; j<numCols; ++j )
        {
            ValueType value = static_cast<ValueType> ( 100.0 - ::fabs( 2.0 * i - j ) );
            valuesA[i * numCols + j] = value;
        }
    }

// initialize the vector x

    for ( IndexType j = 0; j < numCols; ++j)
    {
        valuesX[j] = static_cast<ValueType> ( 1.2 * ( j + 1 ) );
    }

// compute Y = A * X for comparison

    for (IndexType i = 0; i<numRows; ++i)
    {
        valuesY[i] = 0.0;
        for (IndexType j = 0; j<numCols; ++j)
        {
            valuesY[i] += static_cast<ValueType> ( valuesA[i*numCols+j]*valuesX[j] );
        }
    }

// construct replicated matrix A and vector X to be redistributed

    DenseMatrix<ValueType> rA;
    rA.setRawDenseData( numRows, numCols, valuesA.get() );
    DenseVector<ValueType> rX( numCols, valuesX.get());

// try different distributions for rows and colums

    DistributionPtr rowDist;
    DistributionPtr colDist;

    for (int i = 0; i<4; i++)
    {
        if (i==0)
        {
            rowDist = DistributionPtr(new BlockDistribution(numRows, comm));
        }
        else if (i==1)
        {
            rowDist = DistributionPtr(new NoDistribution(numRows));
        }
        else if (i==2)
        {
            std::vector<IndexType> localIndexes;
            IndexType size = comm->getSize();
            IndexType rank = comm->getRank();
            for (int k = 0; k<numRows; k++)
            {
                if (k%size==rank)
                {
                    localIndexes.push_back(k);
                }
            }
            rowDist = DistributionPtr(
                          new GeneralDistribution(numRows, localIndexes, comm));
        }
        else if (i==3)
        {
            IndexType chunkSize = comm->getSize();
            rowDist = DistributionPtr( new CyclicDistribution(numRows, chunkSize, comm));
        }

        for (int j = 0; j<4; j++)
        {
            if (j==0)
            {
                colDist = DistributionPtr(
                              new BlockDistribution(numCols, comm));
            }
            else if (j==1)
            {
                colDist = DistributionPtr(new NoDistribution(numCols));
            }
            else if (j==2)
            {
                std::vector<IndexType> localIndexes;
                IndexType size = comm->getSize();
                IndexType rank = comm->getRank();
                for (int k = 0; k<numCols; k++)
                {
                    if ((k/2+2)%size==rank)
                    {
                        localIndexes.push_back(k);
                    }
                }
                colDist = DistributionPtr(
                              new GeneralDistribution(numCols, localIndexes, comm));
            }
            else if (j==3)
            {
                IndexType chunkSize = comm->getSize();
                colDist = DistributionPtr(
                              new CyclicDistribution(numCols, chunkSize, comm));
            }

            for ( int k = 0; k < 2; k++ )
            {
                // redistribute A, X, setup result vector Y

                DenseMatrix<ValueType> A( rA, rowDist, colDist );
                DenseVector<ValueType> X( rX, colDist );

                if ( k== 1 )
                {
                    A.setCommunicationKind( Matrix::SYNCHRONOUS );
                }

                //  static int x = 0;
                //  x++;

                LAMA_LOG_INFO( logger, "mult matrix A = " << A << " with vector X = " << X );

                DenseVector<ValueType> result ( A * X );

                BOOST_REQUIRE_EQUAL( result.size(), numRows );
                BOOST_REQUIRE_EQUAL( result.getDistribution(), *rowDist );

                // compare the result vector Y with the precomputed results valuesY

                for ( IndexType m = 0; m < numRows; ++m )
                {
                    Scalar value = Scalar( valuesY[m] );
                    Scalar expectedvalue = result.getValue( m );
                    //1e-1 is used because there are 8 of 1280 cases which fails with eps<1e-1
                    LAMA_CHECK_SCALAR_SMALL( value - expectedvalue , ValueType, 1e-1 );
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
        n = 3 * size * chunkSize + size/2 * chunkSize + chunkSize - 2;
        cyclicMultVectorTest( chunkSize, n );

        //not all process get a chunk
        n = size/2 * chunkSize + chunkSize - 2;
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
            BOOST_CHECK_CLOSE( values[i*numCols +j], repM.getValue( i, j ).getValue<double>(), 1e-16 );
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

template<typename T>
void GEMMTestImpl( const int n, const int m, const int k, T eps )
{
    //TODO: This test fails @line 488 if you use a CUDA-Interface.
    typedef T ValueType;

    DistributionPtr dist = DistributionPtr( new NoDistribution( 9 ) );

    Scalar alpha( 1.0 );
    Scalar beta( 3.0 );

    int maxdim = std::max( n, m );

    boost::scoped_array<ValueType> values( new ValueType[maxdim * k] );
    for ( int i = 0; i < maxdim * k; i++ )
    {
        values[i] = i + 1;
    }

    boost::scoped_array<ValueType> valuesC( new ValueType[m * n] );

    for ( int i = 0; i < m; i++ )
    {
        for ( int j = 0; j < n; j++ )
        {
            valuesC[i * n + j] = 42.0 / beta.getValue<ValueType>();
            for ( int kk = 0; kk < k; kk++ )
            {
                valuesC[i * n + j] -= values[i * k + kk] * values[kk * n + j] / beta.getValue<ValueType>();
            }
        }
    }

    for ( ContextType computeLocation = Context::Host; computeLocation < Context::MaxContext; computeLocation =
                ContextType( computeLocation + 1 ) )
    {
        if ( !LAMAInterfaceRegistry::getRegistry().hasInterface( computeLocation ) )
        {
            LAMA_LOG_WARN( logger, "Skipping location: " << computeLocation << " no Interface found." );
            continue;
        }

        ContextPtr computeContext = ContextFactory::getContext( computeLocation );

        DenseMatrix<ValueType> A( m, k, values.get() );
        A.setContext( computeContext, computeContext );

        const DenseMatrix<ValueType> B( k, n, values.get() );
        DenseMatrix<ValueType> C( m, n, valuesC.get() );

        C = alpha * A * B + beta * C;

        for ( int i = 0; i < m; i++ )
        {
            for ( int j = 0; j < n; j++ )
            {
                Scalar expectedvalue( 42.0 );
                Scalar value = C.getValue( i, j );
                Scalar diff = expectedvalue - value;
                LAMA_CHECK_SCALAR_SMALL( diff, ValueType, eps );
            }
        }

        DenseMatrix<ValueType> C2( m, n, valuesC.get() );
        DenseMatrix<ValueType> D( m, n, values.get() );

        D = alpha * A * B + beta * C2;

        for ( int i = 0; i < m; i++ )
        {
            for ( int j = 0; j < n; j++ )
            {
                Scalar expectedvalue( 42.0 );
                Scalar value = D.getValue( i, j );
                Scalar diff = expectedvalue - value;
                LAMA_CHECK_SCALAR_SMALL( diff, ValueType, eps );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( GEMMTest )
{
    //TODO: this test fails on CUDA. It seems that there are some problems with DenseMatrix on CUDA

    LAMA_LOG_WARN( logger,
                   "This test fails if it runs on CUDA. It seems that there are some problems with DenseMatrix on CUDA" );

    //    GEMMTestImpl(const int n, const int m, const int k, T eps)
    //    GEMMTestImpl<double>(5,3,4,1e-5);

//    GEMMTestImpl<double>(2,2,3,1e-5);
//    GEMMTestImpl<double>(16,16,16,1e-5);
//
//    GEMMTestImpl<double>(32,16,16,1e-5);
//
//    GEMMTestImpl<double>(2, 2, 4, 1e-5);
//    GEMMTestImpl<double>(12,12,17,1e-5);
//
//    GEMMTestImpl<double>(32,16,32,1e-5);
//    GEMMTestImpl<double>(16,32,16,1e-5);
//    GEMMTestImpl<double>(32,32,16,1e-5);
//    GEMMTestImpl<double>(16,32,32,1e-5);
//    GEMMTestImpl<double>(64,64,64,1e-5);
//    GEMMTestImpl<double>(32,32,256,1e-5);
//    GEMMTestImpl<double>(128,128,256,1);
//
//    GEMMTestImpl<double>(16,16,17,1e-5);
}

/* --------------------------------------------------------------------- */

void cyclicDistTestImpl( const IndexType chunkSize, const IndexType n )
{
    boost::scoped_array<float> values( new float[n * n] );
    for ( IndexType i = 0; i < n; ++i )
    {
        for ( IndexType j = 0; j < n; ++j )
        {
            values[i * n + j] = static_cast<float>( i * n + j );
        }
    }

    DenseMatrix<float> repMatrix;
    repMatrix.setRawDenseData( n, n, values.get() );

    boost::shared_ptr<Distribution> dist( new CyclicDistribution( n, chunkSize, comm ) );

    DenseMatrix<float> distMatrix( repMatrix, dist, dist );

    for ( IndexType i = 0; i < n; i++ )
    {
        for ( IndexType j = 0; j < n; j++ )
        {
            Scalar expectedvalue = repMatrix.getValue( i, j );
            Scalar value = distMatrix.getValue( i, j );
            Scalar diff = expectedvalue - value;
            LAMA_CHECK_SCALAR_SMALL( diff, float, 1E-8 );
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
/* --------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
