/**
 * @file lamaSolver.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
#include <scai/lama/io/FileIO.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/lama/norm/Norm.hpp>
#include <scai/lama/matrix/StencilMatrix.hpp>

#include <scai/solver/GMRES.hpp>
#include <scai/solver/SimpleAMG.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

#include <scai/tracing.hpp>

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
     
static bool isNumeric( RealType& val, const string& str )
{
    bool is = true;

    for ( std::string::const_iterator p = str.begin(); str.end() != p; ++p )
    {
        if ( isdigit( *p ) )
        {
            continue;
        }

        if ( *p == '-' )
        {
            continue;
        }

        if ( *p == '.' )
        {
            continue;
        }

        is = false;
        break;
    }

    if ( is )
    {
        istringstream valStr( str );
        valStr >> val;
        is = ! valStr.fail();
    }

    return is;
}

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
    SCAI_REGION( "Main.lamaSolver" )

    common::Settings::parseArgs( argc, argv );

    LamaConfig lamaconf;   // must be defined after parseArgs

    const Communicator& comm = lamaconf.getCommunicator();

    int myRank   = comm.getRank();
    int numProcs = comm.getSize();

    // accept only 2 - 4 arguments, --SCAI_xxx do not count here

    if ( argc < 2  || argc > 5 )
    {
        if ( myRank == 0 )
        {
            CONFIG_ERROR( "Illegal number of arguments" )
        }

        return -1;
    }

    try
    {
        std::istringstream stencilSpecification( argv[1] );
        std::string rhsFilename           = argc <= 2 ? "" : argv[2];
        std::string startSolutionFilename = argc <= 3 ? "" : argv[3];
        std::string finalSolutionFilename = argc <= 4 ? "" : argv[4];

        // use auto pointer so that matrix will be deleted at program exit

        IndexType nDims = 0;
        IndexType nPoints = 0;
        IndexType n1 = 0;

        stencilSpecification >> nDims >> nPoints >> n1;

        scai::common::unique_ptr<Matrix> matrixPtr;

        switch ( nDims )
        {
            case 1:  
            {
                common::Stencil1D<RealType> stencil( nPoints );
                common::Grid1D grid( n1 );
                matrixPtr.reset( new StencilMatrix<RealType>( grid, stencil ) );
                break;
            }
         
            case 2:  
            {
                common::Stencil2D<RealType> stencil( nPoints );
                IndexType n2 = n1;
                stencilSpecification >> n2;
                common::Grid2D grid( n1, n2 );
                matrixPtr.reset( new StencilMatrix<RealType>( grid, stencil ) );
                break;
            }
         
            case 3:  
            {
                common::Stencil3D<RealType> stencil( nPoints );
                IndexType n2 = n1;
                stencilSpecification >> n2;
                IndexType n3 = n2;
                stencilSpecification >> n3;
                common::Grid3D grid( n1, n2, n3 );
                matrixPtr.reset( new StencilMatrix<RealType>( grid, stencil ) );
                break;
            }
        }

        std::cout << "Stencil matrix = " << *matrixPtr << std::endl;

        scai::common::unique_ptr<Vector> rhsPtr( matrixPtr->newVector() );
        scai::common::unique_ptr<Vector> solutionPtr( rhsPtr->newVector() );

        Matrix& matrix   = *matrixPtr;
        Vector& rhs      = *rhsPtr;
        Vector& solution = *solutionPtr;

        HOST_PRINT( 0, "Setup matrix, rhs, initial solution" )

        {
            SCAI_REGION( "Main.loadData" )

            RealType val;

            if ( isNumeric( val, rhsFilename ) )
            {
                rhs.allocate( matrix.getRowDistributionPtr() );
                rhs = Scalar( val );

                HOST_PRINT( myRank, "Set rhs = " << val )
            }
            else if ( rhsFilename.size() )
            {
                rhs.readFromFile( rhsFilename, matrix.getRowDistributionPtr() );

                HOST_PRINT( myRank, "Read rhs from file " << rhsFilename << " : " << rhs )
            }
            else
            {
                // build default rhs as rhs = A * x with x = 1

                scai::common::unique_ptr<Vector> xPtr( rhs.newVector() );
                Vector& x = *xPtr;
                x.allocate( matrix.getColDistributionPtr() );
                x = Scalar( 1 );
                rhs = matrix * x;

                HOST_PRINT( myRank, "Set rhs = sum( Matrix, 2) : " << rhs )
            }

            if ( isNumeric( val, startSolutionFilename ) )
            {
                solution.allocate( matrix.getRowDistributionPtr() );
                solution = Scalar( val );
                HOST_PRINT( myRank, "Set initial solution = " << val )
            }
            else if ( startSolutionFilename.size() )
            {
                solution.readFromFile( startSolutionFilename, matrix.getRowDistributionPtr() );
                HOST_PRINT( myRank, "Read initial solution from file " << startSolutionFilename << ": " << solution )
            }
            else
            {
                solution.allocate( matrix.getRowDistributionPtr() );
                solution = Scalar( 0 );   // initialize of a vector
                HOST_PRINT( myRank, "Set initial solution = 0" )
            }

            // only square matrices are accetpted

            SCAI_ASSERT_EQUAL( matrix.getNumRows(), rhs.size(), "size mismatch: #rows of matrix must be equal size of rhs" )
            SCAI_ASSERT_EQUAL( matrix.getNumColumns(), solution.size(),
                               "size mismatch: #cols of matrix must be equal size of initial solution" )
        }

        // distribute data (trivial block partitioning)

        if ( numProcs > 1 )
        {
            SCAI_REGION( "Main.redistribute" )

            LamaTiming timer( comm, "Redistribution" );

            DistributionPtr oldDist = matrix.getRowDistributionPtr();

            DistributionPtr dist = matrix.getRowDistributionPtr();

            rhs.redistribute ( dist );
            solution.redistribute ( dist );

            HOST_PRINT( myRank, "matrix redistributed = " << matrix )
        }

        // Now convert to th desired matrix format

        RealType matrixSize  = matrix.getMemoryUsage() / 1024.0 / 1024.0;

        HOST_PRINT( myRank,  "Matrix Size = " << matrixSize << " MB" )

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

        RealType eps = lamaconf.getAbsoluteTolerance();

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
            SCAI_REGION( "Main.solveSetup" )
            LamaTiming timer( comm, "Solver setup" );
            mySolver->initialize( matrix );
        }

        // Solve system with timing

        double solverTime;   // saves run-time spent in solve

        {
            SCAI_REGION( "Main.solveIt" )
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

        if ( finalSolutionFilename.size() )
        {
            if ( FileIO::fileExists( finalSolutionFilename ) )
            {
                HOST_PRINT( myRank, "Compare solution with vector in " << finalSolutionFilename )
                LamaTiming timer( comm, "Comparing solution" );
                scai::common::unique_ptr<Vector> compSolutionPtr( rhs.newVector() );
                Vector& compSolution = *compSolutionPtr;
                compSolution.readFromFile( finalSolutionFilename );
                compSolution.redistribute( solution.getDistributionPtr() );
                compSolution -= solution;
                Scalar maxDiff = compSolution.maxNorm();
                HOST_PRINT( myRank, "Maximal difference between solution in " << finalSolutionFilename << ": " << maxDiff )
            }
            else
            {
                HOST_PRINT( myRank, "Write solution to output file " << finalSolutionFilename )
                LamaTiming timer( comm, "Writing solution" );
                solution.writeToFile( finalSolutionFilename );
            }
        }
    }
    catch ( common::Exception& e )
    {
        HOST_PRINT( myRank, "Terminate due to error: " << e.what() )
    }
}
