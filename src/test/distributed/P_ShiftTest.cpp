/**
 * @file P_ShiftTest.cpp
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
 * @brief Contains the implementation of the class P_ShiftTest.
 * @author: Alexander Büchel, schubert
 * @date 12.03.2012
 * $
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/Distribution.hpp>
#include <lama/ContextFactory.hpp>

#include <lama/CommunicatorFactory.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>

#include <lama/CSRSparseMatrix.hpp>

#include <lama/distributed/loadBalancing/DistributionStrategyFactory.hpp>
#include <lama/distributed/loadBalancing/BalanceDistribution.hpp>
#include <lama/distributed/BlockDistribution.hpp>

#include <lama/matutils/MatrixCreator.hpp>

#include <omp.h>
#include <mpi.h>

#include <iomanip>
#include <ios>

using namespace lama;
using namespace boost;

const int MASTER = 0;

typedef boost::mpl::list<double,float> test_types;

/* ------------------------------------------------------------------------- */

struct P_ShiftTestConfig
{
    P_ShiftTestConfig()
    {
        mComm = CommunicatorFactory::getCommunicator( Communicator::MPI );

        //TODO: Do we need this stuff here anymore?
        // getting problem: "Poisson( 3D7P_50_50_50 )", "Poisson( 2D9P_10_10 )"
        IndexType pDim = 1; // 1d, 2d, 3d
        IndexType pPoints = 3; // 1d: 3; 2d: 5, 9; 3d: 7, 27
        IndexType pSize[] =
        { 20, 100, 100 };

        std::stringstream dimension;
        for ( IndexType i = 0; i < pDim; ++i )
        {
            dimension << "_" << pSize[i];
        }
        dim = dimension.str();

        problem << "Poisson( " << pDim << "D" << pPoints << "P" << dim << " )";

        eps = 0.05; // = 5%
    }

    ~P_ShiftTestConfig()
    {
        mComm = CommunicatorPtr();
        CommunicatorFactory::release();
    }

    CommunicatorPtr mComm;
    std::stringstream problem;
    double eps;
    std::string dim;
};

BOOST_FIXTURE_TEST_SUITE( P_ShiftTest, P_ShiftTestConfig )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.P_ShiftTest" );

/* ------------------------------------------------------------------------- */

bool checkDeviation( std::vector<double> times, double eps )
{
    IndexType size = times.size();
    double sum = 0.0;
    for ( IndexType i = 0; i < size; ++i )
    {
        sum += times[i];
    }
    double average = sum / size;

    // check balance

    for ( IndexType i = 0; i < size; ++i )
    {
        if ( fabs( times[i] - average ) / average > eps )
        {
            return false;
        }
    }
    return true;
}

template<typename MatrixType>
void printInformations( CommunicatorPtr mComm, MatrixType /*matrix*/, BalanceDistributionPtr dist, double myTotalTime )
{
    IndexType size = mComm->getSize();
    IndexType myRank = mComm->getRank();

    std::vector<double> totalTimes( size );
    mComm->gather( &totalTimes[0], 1, MASTER, &myTotalTime );

    double maxTotalTime = totalTimes[0];
    for ( IndexType i = 1; i < size; ++i )
    {
        if ( totalTimes[i] > maxTotalTime )
        {
            maxTotalTime = totalTimes[i];
        }
    }

    std::vector<float> weights( size );

    if ( myRank == MASTER )
    {
        weights = dist->getWeights();
        std::stringstream weightstr;
        std::stringstream timestr;
        for ( IndexType i = 0; i < size; ++i )
        {
            weightstr << weights[i] << " ";
            timestr << totalTimes[i] << " ";
        }

        std::cout << ::std::setprecision( 2 ) << weightstr.str() << std::endl << ::std::setprecision( 8 )
                  << ::std::fixed << timestr.str() << std::endl << "max " << maxTotalTime << std::endl;
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( shiftTest, T, test_types ) {
//TODO:unknown location(0): fatal error in "shiftTest<d>": std::exception:
//LAMA error in line 226 of file ~/CUDACSRInterface.cpp:
//lama_CSPBLAS_DCSRAGEMVPBVAsync_cuda
//Error: LAMA_STATUS_SCSRAGEMVPBV_CUDAKERNEL_FAILED

    IndexType size = mComm->getSize();
    IndexType myRank = mComm->getRank();

    /* test configuration */

    /* default: 1 Thread, CUDA device 0 */
    int noThreads = 1;
    int devNo = 0;

    DistributionPtr mDistribution;
    ContextPtr mLocalContext = ContextFactory::getContext( Context::Host);
    ContextPtr mHaloContext = ContextFactory::getContext( Context::Host);
    Matrix::SyncKind mCommunicationKind = Matrix::ASYNCHRONOUS;
    Matrix::SyncKind mComputeKind = Matrix::ASYNCHRONOUS;

//    if(size != 8)
//    {
//        std::cerr << "Start test suite for configuration 1 with 8 processes.";
//        exit(0);
//    }

    if( myRank % 2 == 0 )// even processes calculate on CUDA device
    {
        if( myRank % 4 == 0 ) // every second GPU is device Number 1
        {
            devNo = 1;
        }

        mLocalContext = ContextFactory::getContext( Context::CUDA, devNo);
        mHaloContext = ContextFactory::getContext( Context::CUDA, devNo );
    }
//    else // odd processes calculate on CUDA with 5 Threads
//    {
//        noThreads = 5;
//    }

    omp_set_num_threads( noThreads );

    /* init */
    std::stringstream numProc;
    numProc << size;

    std::string numProcesses = numProc.str();

    const IndexType N1 = 100;
    const IndexType N2 = 100;
    CSRSparseMatrix<double> inputA;
    MatrixCreator<double>::buildPoisson2D( inputA, 9, N1, N2);
//    const SBLASInputSet& inputSet = bf::InputSetRegistry<SBLASInputSet>::getRegistry().get(problem.str());
//    const CSRSparseMatrix<double>& inputA = inputSet.getA();

//    std::string prefix = "/home/lschubert/workspace/LAMA/res/testfiles";//Configuration::getInstance( ).getPath( );
//    std::string filename = prefix + "/dwt___59.mtx";
//    const sblas::CSRSparseMatrix<double> inputA(filename);

    CSRSparseMatrix<double> repM( inputA );

    if(myRank == MASTER)
    {
        std::cout << "total num Rows " << repM.getNumRows() << " total non zeros " << repM.getLocalNumValues() << std::endl;
    }

//    /* calc */
    std::string typeNames[] =
    {   "Default", "Metis"};
    DistributionStrategyType types[] =
    {   DistributionStrategy::Default, DistributionStrategy::Metis};

    double start, myExecutionTime;
    for(IndexType i = 0; i < 2; i++)
    {
        mComm->synchronize();
        if(myRank == MASTER)
        {
            std::cout << "**********" << typeNames[i] << "**********" <<std::endl;
        }
        BalanceDistributionPtr distC =
            DistributionStrategyFactory::getDistribution( types[i], mComm, repM );

        // redistribute from BlockDistribution to BalanceDistribution Default Recalc Strategy
        CSRSparseMatrix<double> distMC( repM, distC, distC );

        DenseVector<double> rhs( 1.0, distC );
        DenseVector<double> solution( 0.0, distC );
        Vector& sol = solution;

        distMC.setCommunicationKind( mCommunicationKind );
        distMC.setComputeKind( mComputeKind );
        distMC.setContext( mLocalContext, mHaloContext );

        // SpMV
        mComm->synchronize();
        start = MPI_Wtime();
        sol = distMC * rhs;
        myExecutionTime = MPI_Wtime() - start;

        printInformations( mComm, distMC, distC, myExecutionTime );

        std::vector<double> executionTimes( size );

        mComm->gather( &executionTimes[0], 1, MASTER, &myExecutionTime );
        mComm->bcast( &executionTimes[0], size, MASTER );

        std::vector<float> weights(size);

        distC->printDistributionVector( problem.str() + "_" + typeNames[i] + "_" + numProcesses + "_proc_init" );

        // check balance
        bool balanced = checkDeviation(executionTimes, eps);

        IndexType count = 1;
        while( !balanced && count < 2 )//10)
        {
            mComm->synchronize();
            if(myRank == MASTER)
            {
                std::cout << "*** it " << count << " ***" <<std::endl;
            }

            {
                LAMA_REGION("shift");
                distC = DistributionStrategyFactory::shift( types[i], mComm, distMC, distC, executionTimes );
            }

            CSRSparseMatrix<double> distMC2(repM, distC, distC);
            DenseVector<double> rhs2(1.0, distC);
            DenseVector<double> solution2(0.0, distC);
            Vector& sol2 = solution2;

            distMC2.setCommunicationKind( mCommunicationKind );
            distMC2.setComputeKind( mComputeKind );
            distMC2.setContext( mLocalContext, mHaloContext );

            mComm->synchronize();
            start = MPI_Wtime();
            sol2 = distMC2 * rhs2;
            myExecutionTime = MPI_Wtime() - start;

            printInformations(mComm, distMC2, distC, myExecutionTime);

            mComm->gather(&executionTimes[0], 1, MASTER, &myExecutionTime);
            mComm->bcast(&executionTimes[0], size, MASTER);

            // check balance
            balanced = checkDeviation(executionTimes, eps);

            std::stringstream iter;
            iter << count;

            distC->printDistributionVector( problem.str() + dim + "_" + typeNames[i] + "_" + numProcesses + "_proc_" + iter.str() + "_it" );
            count++;
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
