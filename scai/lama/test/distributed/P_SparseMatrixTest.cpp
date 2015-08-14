/**
 * @file P_SparseMatrixTest.cpp
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
 * @brief Contains the implementation of the class P_SparseMatrixTest.
 * @author: Alexander Büchel, Thomas Brandes
 * @date 26.04.2012
 * @since 1.0.0
 **/

#include <test/distributed/P_SparseMatrixTest.hpp>
#include <test/SameMatrixHelper.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>

#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/distribution/NoDistribution.hpp>
#include <scai/lama/distribution/BlockDistribution.hpp>
#include <scai/lama/distribution/GenBlockDistribution.hpp>
#include <scai/lama/distribution/CyclicDistribution.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <test/TestMacros.hpp>

using namespace lama;
using namespace memory;

/* --------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename MatrixType>, P_SparseMatrixTest<MatrixType>::logger,
                              "Test.P_SparseMatrixTest" );

/* --------------------------------------------------------------------- */

template<typename MatrixType>
P_SparseMatrixTest<MatrixType>::P_SparseMatrixTest()
{
    comm = Communicator::get( );
}

template<typename MatrixType>
P_SparseMatrixTest<MatrixType>::~P_SparseMatrixTest()
{
    comm = CommunicatorPtr();
}

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( P_SparseMatrixTest, MatrixType, cTorTest )
{
    LAMA_LOG_INFO( logger, "cTorTest" );
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType n = 4;
    CommunicatorPtr comm = Communicator::get( "MPI" );
    DistributionPtr dist( new BlockDistribution( n, comm ) );
    CSRSparseMatrix<double> matrix( dist, dist );
    BOOST_CHECK_EQUAL( matrix.getNumRows() , n );
    BOOST_CHECK_EQUAL( matrix.getNumColumns() , n );
    CSRSparseMatrix<ValueType> matrix2( n, n );
    BOOST_CHECK_EQUAL( matrix2.getNumRows() , n );
    BOOST_CHECK_EQUAL( matrix2.getNumColumns() , n );
    CSRSparseMatrix<ValueType> matrix3;
    BOOST_CHECK_EQUAL( matrix3.getNumRows() , 0 );
    BOOST_CHECK_EQUAL( matrix3.getNumColumns() , 0 );
}
LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* ------------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( P_SparseMatrixTest, MatrixType, repDistTest )
{
    const IndexType N1 = 5;
    const IndexType N2 = 5;
    LAMA_LOG_INFO( logger, "RepDistTest: Problem size = " << N1 << " x " << N2 );
    std::ostringstream inputSetName;
    inputSetName << "Poisson( 2D9P_" << N1 << "_" << N2 << " )";
    LAMA_LOG_INFO( logger, "repDistTest, input set = " << inputSetName.str() );
    CSRSparseMatrix<double> inputAhelp;
    MatrixCreator<double>::buildPoisson2D( inputAhelp, 9, N1, N2 );
    MatrixType inputA( inputAhelp );
    DistributionPtr repDist( new NoDistribution( inputA.getNumRows() ) );
    DistributionPtr blockDist( new BlockDistribution( inputA.getNumRows(), comm ) );
    LAMA_LOG_INFO( logger, "Construct replicated C from " << inputA );
    MatrixType C ( inputA, repDist, repDist );
    LAMA_LOG_INFO( logger, "Replicated C: " << C );
    testSameMatrix( inputA, C );
    LAMA_LOG_INFO( logger, "Construct distributed D from " << C );
    MatrixType D ( C, blockDist, blockDist );
    testSameMatrix( C, D );
    LAMA_LOG_INFO( logger, "Replicate columns of D " << D );
    D.redistribute( blockDist, repDist );
    testSameMatrix( C, D );
    LAMA_LOG_INFO( logger, "Replicate rows and distribte columns of D (is nonsense)" << D );
// This case is nonsense as each processor will require all values from
// other processors, but should work
    D.redistribute( repDist, blockDist );
    testSameMatrix( C, D );
}
LAMA_COMMON_TEST_CASE_TEMPLATE_END();
/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( P_SparseMatrixTest, MatrixType, replicateTest )
{
    const IndexType N1 = 5;
    const IndexType N2 = 5;
    LAMA_LOG_INFO( logger, "ReplicateTest: Problem size = " << N1 << " x " << N2 );
    std::ostringstream inputSetName;
    inputSetName << "Poisson( 2D9P_" << N1 << "_" << N2 << " )";
    LAMA_LOG_INFO( logger, "get input set for " << inputSetName.str() );
    CSRSparseMatrix<double> inputAhelp;
    MatrixCreator<double>::buildPoisson2D( inputAhelp, 9, N1, N2 );
    MatrixType inputA ( inputAhelp );
    DistributionPtr repDist( new NoDistribution( inputA.getNumRows() ) );
    MatrixType B ( inputA, repDist, repDist );
    testSameMatrix( inputA, B );
    MatrixType C ( inputA );
    testSameMatrix( inputA, C );
}
LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( P_SparseMatrixTest, MatrixType, assignTest )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType N1 = 5;
    const IndexType N2 = 5;
    LAMA_LOG_INFO( logger, "AssignTest: Problem size = " << N1 << " x " << N2 );
    std::ostringstream inputSetName;
    inputSetName << "Poisson( 2D9P_" << N1 << "_" << N2 << " )";
    LAMA_LOG_INFO( logger, "get input set for " << inputSetName.str() );
    CSRSparseMatrix<ValueType> inputAhelp;
    MatrixCreator<ValueType>::buildPoisson2D( inputAhelp, 9, N1, N2 );
    MatrixType inputA ( inputAhelp );
//const IndexType n = inputA.getNumRows();
//DistributionPtr dist = makeDistribution( n, comm, 1 );
// ELLSparseMatrix<double> ellA ( n, dist );  // constructor square matrix
    ELLSparseMatrix<ValueType> ellA;
    LAMA_LOG_INFO( logger, "zero distributed ELL Matrix = " << ellA );
    ellA = inputA;
    LAMA_LOG_INFO( logger, "copied distributed ELL Matrix = " << ellA );
    testSameMatrix( inputA, ellA );
}
LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( P_SparseMatrixTest, MatrixType, transposeTest )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType N1 = 2;
    const IndexType N2 = 3;
    common::scoped_array<ValueType> rawDense( new ValueType[ N1 * N2 ] );

// Initialize matrix with some data

    for ( IndexType i = 0; i < N1; ++i )
    {
        for ( IndexType j = 0; j < N2; ++ j )
        {
            int diff = i - j;

            if ( ( diff < -3 ) || ( diff > 3 ) )
            {
                return;
            }

            rawDense[ i * N2 + j ] = static_cast<ValueType>( 2 * i - j );
        }
    }

    LAMA_LOG_INFO( logger, "TransposeTest: Problem size = " << N1 << " x " << N2 );
    MatrixType repMatrix;
    repMatrix.setRawDenseData( N1, N2, rawDense.get() );
    LAMA_LOG_INFO( logger, "Replicated matrix = " << repMatrix );
    DistributionPtr bdist1 ( new BlockDistribution( N1, comm ) );
    DistributionPtr bdist2 ( new BlockDistribution( N2, comm ) );
    MatrixType distMatrix( repMatrix, bdist1, bdist2 );
    LAMA_LOG_INFO( logger, "Distributed matrix = " << distMatrix );
    bool transposeFlag = true;
    MatrixType transMatrix( distMatrix, transposeFlag );
    LAMA_LOG_INFO( logger, "Distributed transpose matrix = " << transMatrix );
    BOOST_CHECK_EQUAL( distMatrix.getDistribution(), transMatrix.getColDistribution() );
    BOOST_CHECK_EQUAL( distMatrix.getColDistribution(), transMatrix.getDistribution() );
    DistributionPtr rdist1 ( new NoDistribution( N1 ) );
    DistributionPtr rdist2 ( new NoDistribution( N2 ) );
    transMatrix.redistribute( rdist2, rdist1 );
    LAMA_LOG_INFO( logger, "Replicated transpose matrix = " << transMatrix );
    DenseStorage<ValueType> transStorage( transMatrix.getLocalStorage() );
    LAMA_LOG_INFO( logger, "Replicated transposed storage = " << transStorage );
    ReadAccess<ValueType> rValues( transStorage.getData() );

    for ( IndexType i = 0; i < N1; ++i )
    {
        for ( IndexType j = 0; j < N2; ++ j )
        {
            BOOST_CHECK_EQUAL( rawDense[ i * N2 + j ], rValues[ j * N1 + i ] );
        }
    }
}
LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* --------------------------------------------------------------------- */

static DistributionPtr makeDistribution( const IndexType n, CommunicatorPtr comm, int kind )
{
    if ( kind == 0 )
    {
        return DistributionPtr( new BlockDistribution( n, comm ) );
    }
    else if ( kind == 1 )
    {
        return DistributionPtr( new CyclicDistribution( n, 3, comm ) );
    }
    else if ( kind == 2 )
    {
        float weight = 1.0f;

        if ( comm->getRank() % 2 == 1 )
        {
            weight = 0.0001f;
        }

        return DistributionPtr( new GenBlockDistribution( n, weight, comm ) );
    }
    else if ( kind == 3 )
    {
        float weight = 1.0;

        if ( comm->getRank() % 2 == 0 )
        {
            weight = 0.0001;
        }

        return DistributionPtr( new GenBlockDistribution( n, weight, comm ) );
    }
    else if ( kind == 4 )
    {
        PartitionId size = comm->getSize();
        PartitionId rank = comm->getRank();
        IndexType val = 1713;
        std::vector<IndexType> localIndexes;

        for ( IndexType i = 0; i < n; i++ )
        {
            PartitionId owner = val % size;

            if ( owner == rank )
            {
                localIndexes.push_back( i );
            }

            val = ( val * 21 ) % 1913;
        }

        return DistributionPtr( new GeneralDistribution( n, localIndexes, comm ) );
    }
    else
    {
        COMMON_THROWEXCEPTION( "unsupported kind of makeDistribution" );
    }
}

LAMA_COMMON_TEST_CASE_TEMPLATE( P_SparseMatrixTest, MatrixType, createPoissonTest )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType N1 = 5;
    const IndexType N2 = 5;
    const IndexType n = N1 * N2;
    LAMA_LOG_INFO( logger, "CreatePossonTest: Problem size = " << N1 << " x " << N2 );
    std::ostringstream inputSetName;
// generate string for input set, e.g. 2D9P_3_3
    inputSetName << "Poisson( 2D9P_" << N1 << "_" << N2 << " )";
    LAMA_LOG_INFO( logger, "get input set for " << inputSetName.str() );
    CSRSparseMatrix<ValueType> inputAhelp;
    MatrixCreator<ValueType>::buildPoisson2D( inputAhelp, 9, N1, N2 );
    MatrixType inputA ( inputAhelp );
    LAMA_LOG_INFO( logger, inputSetName.str() << ": inputA = " << inputA );
    BOOST_CHECK_EQUAL( inputA.getNumRows(), n );
    BOOST_CHECK_EQUAL( inputA.getNumColumns(), n );
    LAMA_LOG_INFO( logger, "row distribution : " << inputA.getDistribution() );
    LAMA_LOG_INFO( logger, "col distribution : " << inputA.getColDistribution() );
    DistributionPtr blockDist( new BlockDistribution( n , comm ) );
// Version 1 : redistribute in constructor
    LAMA_LOG_INFO( logger, "B1 ( inputA, blockDist, blockDist )" );
    MatrixType B1 ( inputA, blockDist, blockDist );
    LAMA_LOG_INFO( logger, "B2 ( inputA )" );
    MatrixType B2 ( inputA );
// Version 2 : redistribute as method
    LAMA_LOG_INFO( logger, "B2.redistribute( blockDist, blockDist );" );
    B2.redistribute( blockDist, blockDist );
    testSameMatrix( inputA, B2 );

// Now test some other distributions

    for ( int kind = 1; kind <= 4; kind++ )
    {
        DistributionPtr dist = makeDistribution ( n, comm, kind );
        LAMA_LOG_INFO( logger, *comm << ": redistribute to dist = " << *dist );
        B2.redistribute( dist, dist );
        testSameMatrix( inputA, B2 );
        // generate a new one, pointer equality will fail, no comm necessary
        dist = makeDistribution ( n, comm, kind );
        B2.redistribute( dist, dist );
        testSameMatrix( inputA, B2 );
    }
}
LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_RUNNER_TEMPLATE( P_SparseMatrixTest )
{
    repDistTest();
    replicateTest();
    assignTest();
    transposeTest();
    createPoissonTest();
}

/* --------------------------------------------------------------------- */

template class P_SparseMatrixTest<CSRSparseMatrix<double> > ;
template class P_SparseMatrixTest<CSRSparseMatrix<float> > ;
template class P_SparseMatrixTest<ELLSparseMatrix<double> > ;
template class P_SparseMatrixTest<ELLSparseMatrix<float> > ;
template class P_SparseMatrixTest<JDSSparseMatrix<float> > ;
template class P_SparseMatrixTest<JDSSparseMatrix<double> > ;
template class P_SparseMatrixTest<DIASparseMatrix<float> > ;
template class P_SparseMatrixTest<DIASparseMatrix<double> > ;
template class P_SparseMatrixTest<COOSparseMatrix<float> > ;
template class P_SparseMatrixTest<COOSparseMatrix<double> > ;
