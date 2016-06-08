/**
 * @file P_SparseMatrixTest.cpp
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
 * @brief Contains the implementation of the class P_SparseMatrixTest.
 * @author Alexander Büchel, Thomas Brandes
 * @date 26.04.2012
 */

#include <scai/lama/test/distributed/P_SparseMatrixTest.hpp>
#include <scai/lama/test/SameMatrixHelper.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>

#include <scai/lama/storage/DenseStorage.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;
using namespace scai::dmemo;

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename MatrixType>, P_SparseMatrixTest<MatrixType>::logger,
                              "Test.P_SparseMatrixTest" );

/* --------------------------------------------------------------------- */

template<typename MatrixType>
P_SparseMatrixTest<MatrixType>::P_SparseMatrixTest()
{
    comm = Communicator::getCommunicator( );
}

template<typename MatrixType>
P_SparseMatrixTest<MatrixType>::~P_SparseMatrixTest()
{
    comm = CommunicatorPtr();
}

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( P_SparseMatrixTest, MatrixType, cTorTest )
{
    SCAI_LOG_INFO( logger, "cTorTest" );
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType n = 4;
    CommunicatorPtr comm = Communicator::getCommunicator();
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
    SCAI_LOG_INFO( logger, "RepDistTest: Problem size = " << N1 << " x " << N2 );
    std::ostringstream inputSetName;
    inputSetName << "Poisson( 2D9P_" << N1 << "_" << N2 << " )";
    SCAI_LOG_INFO( logger, "repDistTest, input set = " << inputSetName.str() );
    CSRSparseMatrix<double> inputAhelp;
    MatrixCreator<double>::buildPoisson2D( inputAhelp, 9, N1, N2 );
    MatrixType inputA( inputAhelp );
    DistributionPtr repDist( new NoDistribution( inputA.getNumRows() ) );
    DistributionPtr blockDist( new BlockDistribution( inputA.getNumRows(), comm ) );
    SCAI_LOG_INFO( logger, "Construct replicated C from " << inputA );
    MatrixType C ( inputA, repDist, repDist );
    SCAI_LOG_INFO( logger, "Replicated C: " << C );
    testSameMatrix( inputA, C );
    SCAI_LOG_INFO( logger, "Construct distributed D from " << C );
    MatrixType D ( C, blockDist, blockDist );
    testSameMatrix( C, D );
    SCAI_LOG_INFO( logger, "Replicate columns of D " << D );
    D.redistribute( blockDist, repDist );
    testSameMatrix( C, D );
    SCAI_LOG_INFO( logger, "Replicate rows and distribte columns of D (is nonsense)" << D );
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
    SCAI_LOG_INFO( logger, "ReplicateTest: Problem size = " << N1 << " x " << N2 );
    std::ostringstream inputSetName;
    inputSetName << "Poisson( 2D9P_" << N1 << "_" << N2 << " )";
    SCAI_LOG_INFO( logger, "get input set for " << inputSetName.str() );
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
    SCAI_LOG_INFO( logger, "AssignTest: Problem size = " << N1 << " x " << N2 );
    std::ostringstream inputSetName;
    inputSetName << "Poisson( 2D9P_" << N1 << "_" << N2 << " )";
    SCAI_LOG_INFO( logger, "get input set for " << inputSetName.str() );
    CSRSparseMatrix<ValueType> inputAhelp;
    MatrixCreator<ValueType>::buildPoisson2D( inputAhelp, 9, N1, N2 );
    MatrixType inputA ( inputAhelp );
//const IndexType n = inputA.getNumRows();
//DistributionPtr dist = makeDistribution( n, comm, 1 );
// ELLSparseMatrix<double> ellA ( n, dist );  // constructor square matrix
    ELLSparseMatrix<ValueType> ellA;
    SCAI_LOG_INFO( logger, "zero distributed ELL Matrix = " << ellA );
    ellA = inputA;
    SCAI_LOG_INFO( logger, "copied distributed ELL Matrix = " << ellA );
    testSameMatrix( inputA, ellA );
}
LAMA_COMMON_TEST_CASE_TEMPLATE_END();

/* --------------------------------------------------------------------- */

LAMA_COMMON_TEST_CASE_TEMPLATE( P_SparseMatrixTest, MatrixType, transposeTest )
{
    typedef typename MatrixType::MatrixValueType ValueType;
    const IndexType N1 = 2;
    const IndexType N2 = 3;
    scai::common::scoped_array<ValueType> rawDense( new ValueType[ N1 * N2 ] );

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

    SCAI_LOG_INFO( logger, "TransposeTest: Problem size = " << N1 << " x " << N2 );
    MatrixType repMatrix;
    repMatrix.setRawDenseData( N1, N2, rawDense.get() );
    SCAI_LOG_INFO( logger, "Replicated matrix = " << repMatrix );
    DistributionPtr bdist1 ( new BlockDistribution( N1, comm ) );
    DistributionPtr bdist2 ( new BlockDistribution( N2, comm ) );
    MatrixType distMatrix( repMatrix, bdist1, bdist2 );
    SCAI_LOG_INFO( logger, "Distributed matrix = " << distMatrix );
    bool transposeFlag = true;
    MatrixType transMatrix( distMatrix, transposeFlag );
    SCAI_LOG_INFO( logger, "Distributed transpose matrix = " << transMatrix );
    BOOST_CHECK_EQUAL( distMatrix.getDistribution(), transMatrix.getColDistribution() );
    BOOST_CHECK_EQUAL( distMatrix.getColDistribution(), transMatrix.getDistribution() );
    DistributionPtr rdist1 ( new NoDistribution( N1 ) );
    DistributionPtr rdist2 ( new NoDistribution( N2 ) );
    transMatrix.redistribute( rdist2, rdist1 );
    SCAI_LOG_INFO( logger, "Replicated transpose matrix = " << transMatrix );
    DenseStorage<ValueType> transStorage( transMatrix.getLocalStorage() );
    SCAI_LOG_INFO( logger, "Replicated transposed storage = " << transStorage );
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
    SCAI_LOG_INFO( logger, "CreatePossonTest: Problem size = " << N1 << " x " << N2 );
    std::ostringstream inputSetName;
// generate string for input set, e.g. 2D9P_3_3
    inputSetName << "Poisson( 2D9P_" << N1 << "_" << N2 << " )";
    SCAI_LOG_INFO( logger, "get input set for " << inputSetName.str() );
    CSRSparseMatrix<ValueType> inputAhelp;
    MatrixCreator<ValueType>::buildPoisson2D( inputAhelp, 9, N1, N2 );
    MatrixType inputA ( inputAhelp );
    SCAI_LOG_INFO( logger, inputSetName.str() << ": inputA = " << inputA );
    BOOST_CHECK_EQUAL( inputA.getNumRows(), n );
    BOOST_CHECK_EQUAL( inputA.getNumColumns(), n );
    SCAI_LOG_INFO( logger, "row distribution : " << inputA.getDistribution() );
    SCAI_LOG_INFO( logger, "col distribution : " << inputA.getColDistribution() );
    DistributionPtr blockDist( new BlockDistribution( n , comm ) );
// Version 1 : redistribute in constructor
    SCAI_LOG_INFO( logger, "B1 ( inputA, blockDist, blockDist )" );
    MatrixType B1 ( inputA, blockDist, blockDist );
    SCAI_LOG_INFO( logger, "B2 ( inputA )" );
    MatrixType B2 ( inputA );
// Version 2 : redistribute as method
    SCAI_LOG_INFO( logger, "B2.redistribute( blockDist, blockDist );" );
    B2.redistribute( blockDist, blockDist );
    testSameMatrix( inputA, B2 );

// Now test some other distributions

    for ( int kind = 1; kind <= 4; kind++ )
    {
        DistributionPtr dist = makeDistribution ( n, comm, kind );
        SCAI_LOG_INFO( logger, *comm << ": redistribute to dist = " << *dist );
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
