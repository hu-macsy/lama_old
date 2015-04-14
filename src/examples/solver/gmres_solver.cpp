/**
 * @file gmres_solver.cpp
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
 * @brief Example driver program for GMRES solver
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.0.0
 */

#include "LamaConfig.hpp"
#include "LamaTiming.hpp"

#include <lama/matrix/all.hpp>

#include <lama/DenseVector.hpp>
#include <lama/CommunicatorFactory.hpp>
#include <lama/distribution/GenBlockDistribution.hpp>
#include <lama/distribution/MetisDistribution.hpp>

#include <lama/solver/GMRES.hpp>
#include <lama/solver/TrivialPreconditioner.hpp>
#include <lama/solver/logger/CommonLogger.hpp>
#include <lama/solver/criteria/ResidualThreshold.hpp>
#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/norm/L2Norm.hpp>

using namespace std;
using namespace lama;

int main( int argc, char* argv[] )
{
    LamaConfig lamaconf;

    // Get (default) communicator, will be MPI if available

    const Communicator& comm = lamaconf.getCommunicator();

    int myRank   = comm.getRank();
    int numProcs = comm.getSize();

    const char* filename;

    if ( argc < 2 )
    {
        if ( myRank == 0 )
        {
            cout << "Usage: " << argv[0] << " <filename> [Host|CUDA] [CSR|ELL|JDS]" << endl;
        }
        exit( 1 );
    }

    filename = argv[1];

    // take the remaining arguments for configuration

    for ( int i = 2; i < argc; ++i )
    {
        lamaconf.setArg( argv[i] );
    }

    // use auto pointer so that matrix will be deleted at program exit

    auto_ptr<Matrix> matrixPtr;
    auto_ptr<Vector> rhsPtr;

    if ( lamaconf.getValueType() == Scalar::FLOAT )
    {
        matrixPtr.reset( lamaconf.createSparseMatrix<float>() );
        rhsPtr.reset( new DenseVector<float>() );
    }
    else
    {
        matrixPtr.reset( lamaconf.createSparseMatrix<double>() );
        rhsPtr.reset( new DenseVector<double>() );
    }

    Matrix& matrix = *matrixPtr;
    Vector& rhs = *rhsPtr;

    CSRSparseMatrix<double> inMatrix;

    // Each processor should print its configuration

    cout << lamaconf << endl;

    {
        LamaTiming timer( comm, "Loading data" );

        // read matrix + rhs from disk

        inMatrix.readFromFile( filename );
        rhs.readFromFile( filename );

        // only square matrices are accetpted

        LAMA_ASSERT_EQUAL( inMatrix.getNumRows(), inMatrix.getNumColumns() )
        LAMA_ASSERT_EQUAL( inMatrix.getNumRows(), rhs.size() )
    }

    // for solutin create vector with same format/type as rhs, size = numRows, init = 0.0

    auto_ptr<Vector> solutionPtr( rhs.create( rhs.getDistributionPtr() ) );
    Vector& solution = *solutionPtr;

    int numRows = inMatrix.getNumRows();

    solution = 0.0;   // intialize of a vector

    // distribute data (trivial block partitioning)

    if ( numProcs > 1 )
    {
        LamaTiming timer( comm, "Redistribution" );

        // determine a new distribution so that each processor gets part of the matrix according to its weight

        float weight = lamaconf.getWeight();

        DistributionPtr dist;

        if ( lamaconf.useMetis() )
        {
            LamaTiming timer( comm, "Metis" );
            dist.reset( new MetisDistribution( lamaconf.getCommunicatorPtr(), inMatrix, weight ) );
        }
        else
        {
            dist.reset( new GenBlockDistribution( numRows, weight, lamaconf.getCommunicatorPtr() ) );
        }

        inMatrix.redistribute( dist, dist );
        rhs.redistribute ( dist );
        solution.redistribute ( dist );

        cout << comm << ": matrix = " << inMatrix ;
    }

    {
        LamaTiming timer( comm, "Type conversion from CSR<double> to target format" );
        matrix = inMatrix;
    }

    inMatrix.clear();

    double matrixSize  = matrix.getMemoryUsage() / 1024.0 / 1024.0;

    if ( myRank == 0 )
    {
        cout << "Matrix Size = " << matrixSize << " MB" << endl;
    }

    {
        LamaTiming timer( comm, "Prefetching" );

        matrix.setCommunicationKind( lamaconf.getCommunicationKind() );
        matrix.setContext( lamaconf.getContextPtr() );
        rhs.setContext( lamaconf.getContextPtr() );
        solution.setContext( lamaconf.getContextPtr() );

        rhs.prefetch();
        matrix.prefetch();
        matrix.wait();
        rhs.wait();
    }

    // setting up solver from file "solveconfig.txt"

    std::ostringstream loggerName;

    loggerName << "<GMRES>, " << lamaconf.getCommunicator() << ": ";

    LoggerPtr logger( new CommonLogger ( loggerName.str(), lamaconf.getLogLevel(),
                   LoggerWriteBehaviour::toConsoleOnly,
                   std::auto_ptr<Timer>( new Timer() ) ) );

    GMRES mySolver( "GMResSolver", logger );
    mySolver.setKrylovDim( 30 );

    Scalar eps = 0.001;
    NormPtr norm = NormPtr( new L2Norm() );

    CriterionPtr rt( new ResidualThreshold( norm, eps, ResidualThreshold::Relative ) );

    if ( lamaconf.hasMaxIter() )
    {
        CriterionPtr it( new IterationCount( lamaconf.getMaxIter() ) );
 
        // stop if iteration count reached OR residual threshold is reached

        rt.reset( new Criterion ( it, rt, Criterion::OR ) );
    }

    mySolver.setStoppingCriterion( rt );

    // SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    // mySolver.setPreconditioner( preconditioner );

    {
        LamaTiming timer( comm, "Solver setup" );

        mySolver.initialize( matrix );
    }

    {
        LamaTiming timer( comm, "Solver solve" );

        mySolver.solve( solution, rhs );
    }

    bool writeFlag = false;

    if ( writeFlag )
    {
        LamaTiming timer( comm, "Writing solution" );

        solution.writeToFile( "CG_solution" );
    }
}

