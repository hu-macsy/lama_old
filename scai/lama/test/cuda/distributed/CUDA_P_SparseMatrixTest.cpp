/**
 * @file CUDA_P_SparseMatrixTest.cpp
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
 * @author Alexander Büchel, Thomas Brandes
 * @brief Contains the implementation of the class CUDA_P_SparseMatrixTest.
 * @date 02.05.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/CommunicatorFactory.hpp>
#include <scai/lama/ContextFactory.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/dmemo/distribution/BlockDistribution.hpp>
#include <scai/dmemo/distribution/GenBlockDistribution.hpp>
#include <scai/dmemo/distribution/NoDistribution.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/lama/matutils/MatrixCreator.hpp>

#include <scai/lama/test/cuda/CUDAContext.hpp>

#include <scai/lama/test/TestMacros.hpp>

using namespace scai::lama;

typedef boost::mpl::list<double, float> test_types;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct CUDA_P_SparseMatrixTestConfig
{
    CUDA_P_SparseMatrixTestConfig()
    {
        comm = CommunicatorFactory::get( "MPI" );
    }

    ~CUDA_P_SparseMatrixTestConfig()
    {
        comm = CommunicatorPtr();
    }
};

BOOST_FIXTURE_TEST_SUITE( CUDA_P_SparseMatrixTest, CUDA_P_SparseMatrixTestConfig )
;

SCAI_LOG_DEF_LOGGER( logger, "Test.CUDA_P_SparseMatrixTest" );

/* --------------------------------------------------------------------- */

void testEqualMatrix( const _SparseMatrix& m1, const _SparseMatrix& m2 )
{
    // make sure that m1 and m2 are both square matrices of same size
    const IndexType n = m1.getNumRows();
    BOOST_CHECK_EQUAL( n, m1.getNumColumns() );
    BOOST_CHECK_EQUAL( n, m2.getNumRows() );
    BOOST_CHECK_EQUAL( n, m2.getNumColumns() );
    DenseVector<double> XA( m1.getColDistributionPtr(), 1.0 );
    DenseVector<double> XB( m2.getColDistributionPtr(), 1.0 );
    DenseVector<double> YA( m1.getDistributionPtr(), 0.0 );
    DenseVector<double> YB( m2.getDistributionPtr(), 0.0 );
    DenseVector<double> y1( YA );
    const DenseVector<double> x1( XA );
    DenseVector<double> y2( YB );
    const DenseVector<double> x2( XB );
    y1 = m1 * x1; // YA = mA * XA;
    y2 = m2 * x2; // YB = mB * XB;

    for ( IndexType i = 0; i < n; i++ )
    {
        SCAI_CHECK_SCALAR_CLOSE( YA.getValue( i ), YB.getValue( i ), double, 1.0e-5 );
    }
}

BOOST_AUTO_TEST_CASE( createPoissonTest )
{
    const IndexType N1 = 100;
    const IndexType N2 = 100;
    const IndexType n = N1 * N2;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    std::ostringstream inputSetName;
    // generate string for input set, e.g. 2D9P_3_3
    inputSetName << "Poisson( 2D9P_" << N1 << "_" << N2 << " )";
    SCAI_LOG_INFO( logger, "get input set for " << inputSetName.str() );
    CSRSparseMatrix<double> inputA;
    MatrixCreator<double>::buildPoisson2D( inputA, 9, N1, N2 );
    SCAI_LOG_INFO( logger, inputSetName.str() << ": inputA = " << inputA );
    BOOST_CHECK_EQUAL( inputA.getNumRows(), n );
    BOOST_CHECK_EQUAL( inputA.getNumColumns(), n );
    SCAI_LOG_INFO( logger, "row distribution : " << inputA.getDistribution() );
    SCAI_LOG_INFO( logger, "col distribution : " << inputA.getColDistribution() );
    CommunicatorPtr comm = inputA.getDistribution().getCommunicatorPtr();
    DistributionPtr blockDist( new BlockDistribution( n, comm ) );
    SCAI_LOG_INFO( logger, "B1 ( inputA, blockDist, blockDist )" );
    CSRSparseMatrix<double> B1( inputA, blockDist, blockDist );
    SCAI_LOG_INFO( logger, "B2 ( inputA )" );
    CSRSparseMatrix<double> B2( inputA );
    ContextPtr context = lama_test::CUDAContext::getContext();
    B2.setContext( context, context );
    B2.setCommunicationKind( Matrix::ASYNCHRONOUS );
    // Version 2 : redistribute as method
    SCAI_LOG_INFO( logger, "B2.redistribute( blockDist, blockDist );" );
    B2.redistribute( blockDist, blockDist );
    testEqualMatrix( inputA, B2 );
    float weight = 1.0;

    if ( comm->getRank() % 2 == 1 )
    {
        weight = 0.0001;
    }

    DistributionPtr genDist( new GenBlockDistribution( inputA.getNumRows(), weight, comm ) );
    SCAI_LOG_INFO( logger,
                   *comm << ": redistribute to " << *genDist << ", my local size = " << genDist->getLocalSize() );
    SCAI_LOG_INFO( logger, "B2.redistribute( genDist, genDist );" );
    B2.redistribute( genDist, genDist );
    testEqualMatrix( inputA, B2 );
    weight = 1.0;

    if ( comm->getRank() % 2 == 0 )
    {
        weight = 0.001;
    }

    DistributionPtr genDist1( new GenBlockDistribution( inputA.getNumRows(), weight, comm ) );
    SCAI_LOG_INFO( logger,
                   *comm << ": redistribute to " << *genDist1 << ", my local size = " << genDist1->getLocalSize() );
    SCAI_LOG_INFO( logger, "B2.redistribute( genDist1, genDist1 );" );
    B2.redistribute( genDist1, genDist1 );
    testEqualMatrix( inputA, B2 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( repDistTest )
{
    const IndexType N1 = 5;
    const IndexType N2 = 5;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    std::ostringstream inputSetName;
    inputSetName << "Poisson( 2D9P_" << N1 << "_" << N2 << " )";
    SCAI_LOG_INFO( logger, "repDistTest, input set = " << inputSetName.str() );
    CSRSparseMatrix<double> inputA;
    MatrixCreator<double>::buildPoisson2D( inputA, 9, N1, N2 );
    CommunicatorPtr comm = CommunicatorFactory::get( "MPI" );
    DistributionPtr repDist( new NoDistribution( inputA.getNumRows() ) );
    DistributionPtr blockDist( new BlockDistribution( inputA.getNumRows(), comm ) );
    SCAI_LOG_INFO( logger, "Construct replicated C from " << inputA );
    CSRSparseMatrix<double> C( inputA, repDist, repDist );
    testEqualMatrix( inputA, C );
    SCAI_LOG_INFO( logger, "Construct distributed D from " << C );
    CSRSparseMatrix<double> D( C, blockDist, blockDist );
    testEqualMatrix( C, D );
    SCAI_LOG_INFO( logger, "Replicate columns of D " << D );
    D.redistribute( blockDist, repDist );
    testEqualMatrix( C, D );
    SCAI_LOG_INFO( logger, "Replicate rows and distribte columns of D (is nonsense)" << D );
    // This case is nonsense as each processor will require all values from
    // other processors, but should work
    D.redistribute( repDist, blockDist );
    testEqualMatrix( C, D );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( replicateTest )
{
    const IndexType N1 = 5;
    const IndexType N2 = 5;
    SCAI_LOG_INFO( logger, "Problem size = " << N1 << " x " << N2 );
    std::ostringstream inputSetName;
    inputSetName << "Poisson( 2D9P_" << N1 << "_" << N2 << " )";
    SCAI_LOG_INFO( logger, "get input set for " << inputSetName.str() );
    CSRSparseMatrix<double> inputA;
    MatrixCreator<double>::buildPoisson2D( inputA, 9, N1, N2 );
    DistributionPtr repDist( new NoDistribution( inputA.getNumRows() ) );
    CSRSparseMatrix<double> B( inputA, repDist, repDist );
    testEqualMatrix( inputA, B );
    CSRSparseMatrix<double> C( inputA );
    testEqualMatrix( inputA, C );
}
/* --------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
