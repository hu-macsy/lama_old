/**
 * @file CUDA_P_SparseMatrixTest.cpp
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
 * @endlicense
 *
 * @brief Contains the implementation of the class CUDA_P_SparseMatrixTest.
 * @author Alexander Büchel, Thomas Brandes
 * @date 02.05.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/CommunicatorFactory.hpp>
#include <scai/lama/ContextFactory.hpp>
#include <scai/lama/DenseVector.hpp>

#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>

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
