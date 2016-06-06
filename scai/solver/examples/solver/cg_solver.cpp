/**
 * @file solver/examples/solver/cg_solver.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @brief Example driver program for CG solver
 * @author Thomas Brandes
 * @date 02.07.2012
 */

#include "LamaConfig.hpp"
#include "LamaTiming.hpp"

#include <scai/lama/matrix/all.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/dmemo/Communicator.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

#include <scai/common/unique_ptr.hpp>
#include <scai/common/exception/Exception.hpp>
#include <scai/tracing.hpp>
#include <scai/logging.hpp>

using namespace std;
using namespace scai::dmemo;
using namespace scai::lama;
using namespace scai::solver;

int main( int argc, char* argv[] )
{
    SCAI_REGION( "cg_solver" )

    SCAI_LOG_THREAD( "main" )

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

    // using pointer classes MatrixPtr, VectorPtr take care of delete

    MatrixPtr matrixPtr( lamaconf.getMatrix() );
    VectorPtr rhsPtr( matrixPtr->newDenseVector() );

    Matrix& matrix = *matrixPtr;
    Vector& rhs = *rhsPtr;

    // Each processor should print its configuration

    cout << lamaconf << endl;

    {
        LamaTiming timer( comm, "Loading data" );

        // read matrix + rhs from disk

        matrix.readFromFile( filename );

        try
        {
            rhs.readFromFile( filename );
        }
        catch ( const std::exception& )
        {
            std::cout << "reading vector from file " << filename << " failed, take sum( Matrix, 2 ) " << std::endl;

            {
                scai::common::unique_ptr<Vector> xPtr( rhs.newVector() );
                Vector& x = *xPtr;
                x.allocate( matrix.getColDistributionPtr() );
                x = Scalar( 1 );
                rhs = matrix * x;
            }
        }

        // only square matrices are accetpted

        SCAI_ASSERT_EQUAL( matrix.getNumRows(), matrix.getNumColumns(), "" )
        SCAI_ASSERT_EQUAL( matrix.getNumRows(), rhs.size(), "" )
    }

    // for solution create vector with same format/type as rhs, size = numRows, init = 0.0

    scai::common::unique_ptr<Vector> solutionPtr( Vector::create( rhs.getCreateValue() ) );
    Vector& solution = *solutionPtr;

    int numRows = matrix.getNumRows();

    solution.allocate( matrix.getColDistributionPtr() );
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
            dist.reset( Distribution::getDistribution( "METIS", lamaconf.getCommunicatorPtr(), matrix, weight ) );
        }
        else
        {
            dist.reset( new GenBlockDistribution( numRows, weight, lamaconf.getCommunicatorPtr() ) );
        }

        matrix.redistribute( dist, dist );
        rhs.redistribute ( dist );
        solution.redistribute ( dist );

        cout << comm << ": matrix = " << matrix ;
    }

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

    std::ostringstream loggerName;

    loggerName << "<CG>, " << lamaconf.getCommunicator() << ": ";

    LoggerPtr logger( new CommonLogger ( loggerName.str(),
                                         lamaconf.getLogLevel(),
                                         LoggerWriteBehaviour::toConsoleOnly
                                       ) );

    CG mySolver( "CGSolver", logger );

    Scalar eps = 0.00001;
    NormPtr norm = NormPtr( new L2Norm() );

    CriterionPtr rt( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );

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

