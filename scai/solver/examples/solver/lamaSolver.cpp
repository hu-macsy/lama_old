/**
 * @file lamaSolver.cpp
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
 * @brief Example driver program for solvers in LAMA
 * @author Thomas Brandes
 * @date 04.12.2015
 */

#include "LamaConfig.hpp"
#include "LamaTiming.hpp"

#include <scai/lama/matrix/all.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/lama/norm/Norm.hpp>
#include <scai/lama/StorageIO.hpp>

#include <scai/solver/GMRES.hpp>
#include <scai/solver/SimpleAMG.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

#include <scai/common/unique_ptr.hpp>

using namespace std;
using namespace scai;
using namespace dmemo;
using namespace lama;
using namespace solver;

typedef RealType ValueType;

#define HOST_PRINT( rank, msg )             \
{                                           \
    if ( rank == 0 )                        \
    {                                       \
        std::cout << msg << std::endl;      \
    }                                       \
}                                           \

/** Help routine to combine criterions. */

static void orCriterion( CriterionPtr& crit, const CriterionPtr& add )
{
    if ( crit )
    { 
        // combine current one with the new one by an OR

        crit.reset( new Criterion ( crit, add, Criterion::OR ) );
    }
    else
    { 
        crit = add;
    }
}

/**
 *  Main program
 *
 *  - first arg is filename for input matrix
 *  - all other arguments are passed to the configuration lamaconf
 *  - configuration will contain all information to setup the solver for the input matrix
 */
int main( int argc, const char* argv[] )
{
    common::Settings::parseArgs( argc, argv );

    LamaConfig lamaconf;   // must be defined after parseArgs

    const Communicator& comm = lamaconf.getCommunicator();

    int myRank   = comm.getRank();
    int numProcs = comm.getSize();

    // accept only 2 - 4 arguments

    if ( argc < 2  || argc > 4 )
    {
        if ( myRank == 0 )
        {
            CONFIG_ERROR( "Illegal number of arguments" )
        }

        return -1;
    }

    std::string matrixFilename = argv[1];
    std::string rhsFilename    = argc <= 2 ? "" : argv[2];
    std::string solFilename    = argc <= 3 ? "" : argv[3];

    // use auto pointer so that matrix will be deleted at program exit

    scai::common::unique_ptr<Matrix> matrixPtr( lamaconf.getMatrix() );
    scai::common::unique_ptr<Vector> rhsPtr( matrixPtr->newDenseVector() );

    Matrix& matrix = *matrixPtr;
    Vector& rhs = *rhsPtr;

    // input matrix will be CSR format

    scai::common::unique_ptr<Matrix> inMatrixPtr( Matrix::getMatrix( Matrix::CSR, lamaconf.getValueType() ) );
    Matrix& inMatrix = *inMatrixPtr;

    // Here each processor should print its configuration

    cout << lamaconf << endl;

    // Now read in matrix and rhs

    {
        LamaTiming timer( comm, "Loading data" );

        // read matrix + rhs from disk

        inMatrix.readFromFile( matrixFilename );

        HOST_PRINT( myRank, "Matrix from file " << matrixFilename << " : " << inMatrix )

        if ( rhsFilename.size() == 0 && _StorageIO::hasSuffix( matrixFilename, "frm" ) )
        {
            // this filename can also be taken for vector

            rhsFilename = matrixFilename.substr( 0, matrixFilename.size() - 4 ) + ".frv";
        }

        if ( ! _StorageIO::fileExists( rhsFilename ) )
        {
            HOST_PRINT( myRank, "rhs file " << rhsFilename << " does not exist, take default rhs" )
            rhsFilename = "";
        }

        if ( rhsFilename.size() )
        {
            rhs.readFromFile( rhsFilename );
            HOST_PRINT( myRank, "rhs from file " << rhsFilename << " : " << rhs )
        }
        else 
        {
            // build default rhs as rhs = A * x with x = 1

            scai::common::unique_ptr<Vector> xPtr( rhs.newVector() );
            Vector& x = *xPtr;
            x.allocate( inMatrix.getColDistributionPtr() );
            x = Scalar( 1 );
            rhs = inMatrix * x;

            HOST_PRINT( myRank, "rhs is sum( Matrix, 2) : " << rhs )
        }

        // only square matrices are accetpted
        SCAI_ASSERT_EQUAL( inMatrix.getNumRows(), inMatrix.getNumColumns(), "size mismatch" )
        SCAI_ASSERT_EQUAL( inMatrix.getNumRows(), rhs.size(), "size mismatch" )
    }

    // for solution create vector with same format/type as rhs, size = numRows, init = 0.0

    scai::common::unique_ptr<Vector> solutionPtr( rhs.newVector() );
    Vector& solution = *solutionPtr;
    int numRows = inMatrix.getNumRows();
    solution.allocate( inMatrix.getColDistributionPtr() );
    solution = 0.0;   // initialize of a vector

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

    // Now convert to th desired matrix format

    {
        LamaTiming timer( comm, "Type conversion from CSR<ValueType> to target format" );
        matrix = inMatrix;
    }

    cout << "matrix = " << matrix << endl;

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

    // setting up solver, get it from solver factory

    string solverName = "<" + lamaconf.getSolverName() + ">";

    scai::common::unique_ptr<Solver> mySolver( Solver::create( lamaconf.getSolverName(), solverName ) );

    // setting up a common logger, prints also rank of communicator

    ostringstream loggerName;

    loggerName << solverName << ", " << lamaconf.getCommunicator() << ": ";

    LoggerPtr logger( new CommonLogger ( loggerName.str(),
                                         lamaconf.getLogLevel(),
                                         LoggerWriteBehaviour::toConsoleOnly ) );

    mySolver->setLogger( logger );

    // Set up stopping criterion, take thresholds and max iterations if specified

    CriterionPtr crit;

    NormPtr norm( Norm::create( lamaconf.getNorm() ) );   // Norm from factory

    double eps = lamaconf.getAbsoluteTolerance();

    if ( eps > 0.0 )
    {
        crit.reset( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );
    }

    eps = lamaconf.getRelativeTolerance();

    if ( eps > 0.0 )
    {
        CriterionPtr rt( new ResidualThreshold( norm, eps, ResidualThreshold::Relative ) );

        orCriterion( crit, rt );
    }

    eps = lamaconf.getDivergenceTolerance();

    if ( eps > 0.0 )
    {
        CriterionPtr dt( new ResidualThreshold( norm, eps, ResidualThreshold::Divergence ) );

        orCriterion( crit, dt );
    }

    if ( lamaconf.hasMaxIter() )
    {
        CriterionPtr it( new IterationCount( lamaconf.getMaxIter() ) );

        orCriterion( crit, it );
    }

    if ( !crit )
    {
        HOST_PRINT( myRank, "No criterion set, take default" )

        eps = 0.001;

        crit.reset( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );
    }

    // Stopping criterion can ony be set for an iterative solver

    IterativeSolver* itSolver = dynamic_cast<IterativeSolver*>( mySolver.get() );

    if ( itSolver != NULL )
    {
        itSolver->setStoppingCriterion( crit );
    }
    else
    {
        cout << "Not iterative solver, stopping criterion ignored" << endl;
    }

    // Allow individual settings for GMRES solver

    GMRES* gmresSolver = dynamic_cast<GMRES*>( mySolver.get() );

    if ( gmresSolver )
    {
        // here is the possibility to set solver specific values

        int dim = 5;
 
        common::Settings::getEnvironment( dim, "SCAI_KRYLOV_DIM" );

        HOST_PRINT( myRank, "GMRES solver, krylov dim = " << dim )

        gmresSolver->setKrylovDim( dim );
    }

    // Allow individual settings for AMG solver

    SimpleAMG* amgSolver = dynamic_cast<SimpleAMG*>( mySolver.get() );

    if ( amgSolver )
    {
        amgSolver->setHostOnlyLevel( 4 );
        amgSolver->setReplicatedLevel( 5 );
        amgSolver->setMaxLevels( 25 );
        amgSolver->setMinVarsCoarseLevel( 200 );
    }

    // Initialization with timing

    {
        LamaTiming timer( comm, "Solver setup" );
        mySolver->initialize( matrix );
    }

    // Solve system with timing

    double solverTime;   // saves run-time spent in solve

    {
        LamaTiming timer( comm, "Solver solve" );
        mySolver->solve( solution, rhs );
        solverTime = timer.getTime();
    }

    if ( itSolver != NULL )
    {
        double iterTime = solverTime / itSolver->getIterationCount();
        HOST_PRINT( myRank, "Time per Iteration = " << ( iterTime * 1000.0 ) << " ms" )
    }

    // if 3rd argument for solution file is specified, write it or compare it

    if ( solFilename.size() )
    {
        if ( _StorageIO::fileExists( solFilename ) )
        {
            HOST_PRINT( myRank, "Compare solution with vector in " << solFilename )
            LamaTiming timer( comm, "Comparing solution" );
            scai::common::unique_ptr<Vector> compSolutionPtr( rhs.newVector() );
            Vector& compSolution = *compSolutionPtr;
            compSolution.readFromFile( solFilename );
            compSolution.redistribute( solution.getDistributionPtr() );
            compSolution -= solution;
            Scalar maxDiff = compSolution.maxNorm();
            HOST_PRINT( myRank, "Maximal difference between solution in " << solFilename << ": " << maxDiff )
        }
        else
        {
            HOST_PRINT( myRank, "Write solution to output file " << solFilename << ".mtx (Matrix Market)" )
            LamaTiming timer( comm, "Writing solution" );
            // replicate solution on all processors
            DistributionPtr repDist( new NoDistribution( solution.size() ) );
            solution.redistribute( repDist );
            if ( myRank == 0 )
            {
                solution.writeToFile( solFilename, File::MATRIX_MARKET );
            }
        }
    }
}