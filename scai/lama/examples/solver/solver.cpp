/**
 * @file solver.cpp
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
 * @brief Example driver program for solvers in LAMA
 * @author Thomas Brandes
 * @date 04.12.2015
 */

#include "LamaConfig.hpp"
#include "LamaTiming.hpp"

#include <scai/lama/matrix/all.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/distribution/GenBlockDistribution.hpp>

#include <scai/lama/solver/GMRES.hpp>
#include <scai/lama/solver/TrivialPreconditioner.hpp>
#include <scai/lama/solver/logger/CommonLogger.hpp>
#include <scai/lama/solver/criteria/ResidualThreshold.hpp>
#include <scai/lama/solver/criteria/IterationCount.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/common/unique_ptr.hpp>

using namespace std;
using namespace scai;
using namespace lama;
using common::unique_ptr;

/**
 *  Main program 
 *
 *  - first arg is filename for input matrix
 *  - all other arguments are passed to the configuration lamaconf
 *  - configuration will contain all information to setup the solver for the input matrix
 */
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
            cout << "Usage: " << argv[0] << " <filename> [Host|CUDA] [CSR|ELL|JDS|DIA|COO] [CG|BiCG|...]" << endl;
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

    unique_ptr<Matrix> matrixPtr;
    unique_ptr<Vector> rhsPtr;

    if ( lamaconf.getValueType() == common::scalar::FLOAT )
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
 
        std::cout << "Matrix from file " << filename << " : " << inMatrix << std::endl;

        try
        {
            rhs.readFromFile( filename );
        }
        catch ( const std::exception& )
        {
            std::cout << "reading vector from file " << filename << " failed, take sum( Matrix, 2 ) " << std::endl;

            {
                scai::common::unique_ptr<Vector> xPtr( rhs.clone() );
                Vector& x = *xPtr;
                x.resize( inMatrix.getColDistributionPtr() );
                x = Scalar( 1 );
                rhs.resize( inMatrix.getDistributionPtr() );
                rhs = inMatrix * x;
            }
        }

        // only square matrices are accetpted

        SCAI_ASSERT_EQUAL( inMatrix.getNumRows(), inMatrix.getNumColumns(), "size mismatch" )
        SCAI_ASSERT_EQUAL( inMatrix.getNumRows(), rhs.size(), "size mismatch" )
    }

    // for solutin create vector with same format/type as rhs, size = numRows, init = 0.0

    unique_ptr<Vector> solutionPtr( rhs.clone( rhs.getDistributionPtr() ) );
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
            // dist.reset( new MetisDistribution( lamaconf.getCommunicatorPtr(), inMatrix, weight ) );
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
        matrix.setContextPtr( lamaconf.getContextPtr() );
        rhs.setContextPtr( lamaconf.getContextPtr() );
        solution.setContextPtr( lamaconf.getContextPtr() );

        rhs.prefetch();
        matrix.prefetch();
        matrix.wait();
        rhs.wait();
    }

    // setting up solver from file "solveconfig.txt"

    std::ostringstream solverName;

    solverName << "<" << lamaconf.getSolverName() << ">";

    std::ostringstream loggerName;

    loggerName << solverName.str() << ", " << lamaconf.getCommunicator() << ": ";

    LoggerPtr logger( new CommonLogger ( loggerName.str(), 
                                         lamaconf.getLogLevel(),
                                         LoggerWriteBehaviour::toConsoleOnly ) );

    common::unique_ptr<Solver> mySolver( Solver::create( lamaconf.getSolverName(), solverName.str() ) );

    IterativeSolver* itSolver = dynamic_cast<IterativeSolver*>( mySolver.get() );

    SCAI_ASSERT( itSolver, "Not an iterative solver: " << *mySolver )

    mySolver->setLogger( logger );

    Scalar eps = 0.001;
    NormPtr norm = NormPtr( new L2Norm() );

    CriterionPtr rt( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );

    if ( lamaconf.hasMaxIter() )
    {
        CriterionPtr it( new IterationCount( lamaconf.getMaxIter() ) );
 
        // stop if iteration count reached OR residual threshold is reached

        rt.reset( new Criterion ( it, rt, Criterion::OR ) );
    }

    itSolver->setStoppingCriterion( rt );

    // SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    // mySolver.setPreconditioner( preconditioner );

    {
        LamaTiming timer( comm, "Solver setup" );

        mySolver->initialize( matrix );
    }

    {
        LamaTiming timer( comm, "Solver solve" );

        mySolver->solve( solution, rhs );
    }

    bool writeFlag = false;

    if ( writeFlag )
    {
        LamaTiming timer( comm, "Writing solution" );

        solution.writeToFile( "CG_solution" );
    }
}

