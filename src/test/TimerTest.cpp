/**
 * @file TimerTest.cpp
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
 * @brief Contains the implementation of the class TimerTest.
 * @author Alexander Büchel, schubert
 * @date 02.02.2012
 * $Id$
 */

#include <boost/test/unit_test.hpp>
#include <boost/shared_ptr.hpp>

#include <lama/tracing.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>

#include <lama/distribution/BlockDistribution.hpp>

#include <lama/CSRSparseMatrix.hpp>
#include <lama/ELLSparseMatrix.hpp>
#include <lama/LAMATypes.hpp>

#include <lama/matutils/MatrixCreator.hpp>

#include <lama/ContextFactory.hpp>
#include <lama/Context.hpp>

#include <lama/Communicator.hpp>
#include <lama/CommunicationPlan.hpp>
#include <lama/CommunicatorFactory.hpp>

#include <omp.h>

using namespace boost;
using namespace lama;

/* --------------------------------------------------------------------- */

static CommunicatorPtr mComm;

struct TimerTestConfig
{
    TimerTestConfig()
    {
        mComm = CommunicatorFactory::getCommunicator( Communicator::MPI );
    }

    ~TimerTestConfig()
    {
        mComm = CommunicatorPtr();
    }
};

BOOST_FIXTURE_TEST_SUITE( TimerTest, TimerTestConfig )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.TimerTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( testTimer )
{

    IndexType myRank = mComm->getRank();
    if ( myRank == 0 )
    {
        std::cout << "this test should not print stack warning of LAMATimeTracer." << std::endl;
    }

    int noThreads = 1;
    int devNo = 0;

    DistributionPtr mDistribution;
    ContextPtr mLocalContext = ContextFactory::getContext( Context::Host );
    ContextPtr mHaloContext = ContextFactory::getContext( Context::Host );
    Matrix::SyncKind mCommunicationKind = Matrix::ASYNCHRONOUS;
    Matrix::SyncKind mComputeKind = Matrix::ASYNCHRONOUS;

    if ( myRank == 1 )
    {
        mLocalContext = ContextFactory::getContext( Context::CUDA, devNo );
        mHaloContext = ContextFactory::getContext( Context::CUDA, devNo );
    }
    omp_set_num_threads( noThreads );

    const IndexType N1 = 100;
    const IndexType N2 = 100;
    CSRSparseMatrix<double> tmp;
    MatrixCreator<double>::buildPoisson2D( tmp, 9, N1, N2 );

    ELLSparseMatrix<double> tmp2( tmp );
    IndexType numRows = tmp.getNumRows();

    for ( IndexType i = 0; i < 10; ++i )
    {
        DistributionPtr dist( new BlockDistribution( numRows, mComm ) );
        ELLSparseMatrix<double> A( tmp2, dist, dist );

        DenseVector<double> rhs( 1.0, dist );
        DenseVector<double> solution( 0.0, dist );

        A.setCommunicationKind( mCommunicationKind );
        A.setComputeKind( mComputeKind );
        A.setContext( mLocalContext, mHaloContext );

        solution = A * rhs;
        if ( myRank == 0 )
        {
            std::cout << "iteration " << i << " of " << 10 << " time: " << LAMA_TIMETRACER("SparseMV") << std::endl;
        }
    }

    CommunicatorFactory::release();
}
;
/* ---------------------------------------------------------------------- */
BOOST_AUTO_TEST_SUITE_END();
