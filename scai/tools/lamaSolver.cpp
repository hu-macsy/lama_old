/**
 * @file tools/lamaSolver.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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
#include <scai/dmemo/SingleDistribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/partitioning/Partitioning.hpp>
#include <scai/lama/norm/Norm.hpp>

#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/CG.hpp>
#include <scai/solver/SimpleAMG.hpp>
#include <scai/solver/GMRES.hpp>
#include <scai/solver/Kaczmarz.hpp>
#include <scai/solver/Richardson.hpp>
#include <scai/solver/BiCGstab.hpp>
#include <scai/solver/DecompositionSolver.hpp>
#include <scai/solver/QMR.hpp>
#include <scai/solver/TFQMR.hpp>

#include <scai/tracing.hpp>

#include <memory>

using namespace std;
using namespace scai;
using namespace dmemo;
using namespace lama;
using namespace solver;

typedef DefaultReal ValueType;

#define HOST_PRINT( rank, msg )                 \
    {                                           \
        if ( rank == 0 )                        \
        {                                       \
            std::cout << msg << std::endl;      \
        }                                       \
    }                                           \
     
static bool isNumeric( double& val, const string& str )
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

template<typename ValueType>
static void orCriterion( CriterionPtr<ValueType>& crit, const CriterionPtr<ValueType>& add )
{
    if ( crit )
    {
        // combine current one with the new one by an OR

        crit.reset( new Criterion<ValueType> ( crit, add, BooleanOp::OR ) );
    }
    else
    {
        crit = add;
    }
}

template<typename ValueType>
void doPartitioning( Matrix<ValueType>& matrix, Vector<ValueType>& rhs, Vector<ValueType>& solution )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    if ( comm->getSize() < 2 )
    {
        return;
    }

    SCAI_REGION( "lamaSolver.redistribute" )

    LamaTiming timer( *comm, "Redistribution" );

    // determine a new distribution so that each processor gets part of the matrix according to its weight

    float weight = 1.0f;

    common::Settings::getEnvironment( weight, "SCAI_WEIGHT" );
  
    std::string parKind = "BLOCK";     // default partitioning strategy

    common::Settings::getEnvironment( parKind, "SCAI_PARTITIONING" );

    using namespace partitioning;

    DistributionPtr dist;

    if ( parKind == "OFF" )
    {
        // let it unchanged
 
        dist = matrix.getRowDistributionPtr();
    }
    else 
    {
        PartitioningPtr graphPartitioning( Partitioning::create( parKind ) );

        std::cout << "Partitioning, kind = " << parKind << ", weight = " << weight << ", partitioner = " << *graphPartitioning << std::endl;

        dist = graphPartitioning->partitionIt( comm, matrix, weight );
    }

    matrix.redistribute( dist, dist );

    rhs.redistribute ( dist );
    solution.redistribute ( dist );
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
        std::string matrixFilename        = argv[1];
        std::string rhsFilename           = argc <= 2 ? "" : argv[2];
        std::string startSolutionFilename = argc <= 3 ? "" : argv[3];
        std::string finalSolutionFilename = argc <= 4 ? "" : argv[4];

        // use auto pointer so that matrix will be deleted at program exit

        MatrixPtr<ValueType> matrixPtr( lamaconf.getMatrix<ValueType>() );
        Matrix<ValueType>& matrix = *matrixPtr;

        DenseVector<ValueType> rhs( matrixPtr->getRowDistributionPtr(), 0 );
        DenseVector<ValueType> solution;

        CSRSparseMatrix<ValueType> inMatrix;  // matrix read fromfile

        // Here each processor should print its configuration

        cout << lamaconf << endl;

        // Now read in matrix and rhs

        HOST_PRINT( 0, "Setup matrix, rhs, initial solution" )

        {
            SCAI_REGION( "Main.loadData" )

            LamaTiming timer( comm, "Loading data" );

            // read matrix + rhs from disk

            std::string distFilename = "";

            if ( common::Settings::getEnvironment( distFilename, "SCAI_DISTRIBUTION" ) )
            {
                if ( distFilename.size() < 2 )
                {
                    distFilename = "";
                    HOST_PRINT( 0 , "Read matrix from file " << matrixFilename << ", with mapping by cols" )
                }
                else
                {
                    HOST_PRINT( 0, "Read matrix from file " << matrixFilename << ", with mapping from file " << distFilename )
                }

                inMatrix.readFromFile( matrixFilename, distFilename );
            }
            else
            {
                HOST_PRINT( myRank, "Read matrix from file " << matrixFilename )
                inMatrix.readFromFile( matrixFilename );
            }

            SCAI_ASSERT_EQUAL( inMatrix.getNumRows(), inMatrix.getNumColumns(), "solver only with square matrices" )

            HOST_PRINT( myRank, "Matrix from file " << matrixFilename << " : " << inMatrix )

            if ( rhsFilename.size() == 0 && FileIO::hasSuffix( matrixFilename, "frm" ) )
            {
                // this filename can also be taken for vector

                rhsFilename = matrixFilename.substr( 0, matrixFilename.size() - 4 ) + ".frv";

                if ( ! FileIO::fileExists( rhsFilename ) )
                {
                    HOST_PRINT( myRank, "rhs file " << rhsFilename << " does not exist, take default rhs" )
                    rhsFilename = "";
                }
            }

            double val;

            if ( isNumeric( val, rhsFilename ) )
            {
                rhs.setSameValue( inMatrix.getRowDistributionPtr(), val );

                HOST_PRINT( myRank, "Set rhs = " << val )
            }
            else if ( rhsFilename.size() )
            {
                rhs.readFromFile( rhsFilename, inMatrix.getRowDistributionPtr() );

                HOST_PRINT( myRank, "Read rhs from file " << rhsFilename << " : " << rhs )
            }
            else
            {
                // build default rhs as rhs = A * x with x = 1

                auto x = fill<DenseVector<ValueType>>( inMatrix.getColDistributionPtr(), ValueType( 1 ) );

                rhs = inMatrix * x;

                HOST_PRINT( myRank, "Set rhs = sum( Matrix, 2) : " << rhs )
            }

            if ( isNumeric( val, startSolutionFilename ) )
            {
                solution.setSameValue( inMatrix.getRowDistributionPtr(), val );
                HOST_PRINT( myRank, "Set initial solution = " << val )
            }
            else if ( startSolutionFilename.size() )
            {
                solution.readFromFile( startSolutionFilename, inMatrix.getRowDistributionPtr() );
                HOST_PRINT( myRank, "Read initial solution from file " << startSolutionFilename << ": " << solution )
            }
            else
            {
                solution.setSameValue( inMatrix.getRowDistributionPtr(), 0 );
                HOST_PRINT( myRank, "Set initial solution = 0" )
            }

            // only square matrices are accetpted

            SCAI_ASSERT_EQUAL( inMatrix.getNumRows(), rhs.size(), "size mismatch: #rows of matrix must be equal size of rhs" )
            SCAI_ASSERT_EQUAL( inMatrix.getNumColumns(), solution.size(),
                               "size mismatch: #cols of matrix must be equal size of initial solution" )
        }

        doPartitioning( inMatrix, solution, rhs );

        // Now convert to th desired matrix format

        {
            LamaTiming timer( comm, "Type conversion from CSR<ValueType> to target format" );
            matrix = inMatrix;
        }

        cout << "matrix = " << matrix << endl;

        inMatrix.clear();
        double matrixSize  = matrix.getMemoryUsage() / 1024.0 / 1024.0;

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

        std::unique_ptr<Solver<ValueType> > mySolver( Solver<ValueType>::getSolver( lamaconf.getSolverName() ) );

        // setting up a common logger, prints also rank of communicator

        ostringstream loggerName;

        loggerName << solverName << ", " << lamaconf.getCommunicator() << ": ";

        bool suppressWriting = comm.getRank() > 0;   // only host processor logs

        LoggerPtr logger( new CommonLogger ( loggerName.str(),
                                             lamaconf.getLogLevel(),
                                             LoggerWriteBehaviour::toConsoleOnly,
                                             suppressWriting ) );

        mySolver->setLogger( logger );

        // Set up stopping criterion, take thresholds and max iterations if specified

        CriterionPtr<ValueType> crit;

        NormPtr<ValueType> norm( Norm<ValueType>::create( lamaconf.getNorm() ) );   // Norm from factory

        double eps = lamaconf.getAbsoluteTolerance();

        if ( eps > 0.0 )
        {
            crit.reset( new ResidualThreshold<ValueType>( norm, eps, ResidualCheck::Absolute ) );
        }

        eps = lamaconf.getRelativeTolerance();

        if ( eps > 0.0 )
        {
            CriterionPtr<ValueType> rt( new ResidualThreshold<ValueType>( norm, eps, ResidualCheck::Relative ) );

            orCriterion( crit, rt );
        }

        eps = lamaconf.getDivergenceTolerance();

        if ( eps > 0.0 )
        {
            CriterionPtr<ValueType> dt( new ResidualThreshold<ValueType>( norm, eps, ResidualCheck::Divergence ) );

            orCriterion( crit, dt );
        }

        if ( lamaconf.hasMaxIter() )
        {
            CriterionPtr<ValueType> it( new IterationCount<ValueType>( lamaconf.getMaxIter() ) );

            orCriterion( crit, it );
        }

        if ( !crit )
        {
            HOST_PRINT( myRank, "No criterion set, take default" )

            eps = 0.001;

            crit.reset( new ResidualThreshold<ValueType>( norm, eps, ResidualCheck::Absolute ) );
        }

        // Stopping criterion can ony be set for an iterative solver

        IterativeSolver<ValueType>* itSolver = dynamic_cast<IterativeSolver<ValueType>*>( mySolver.get() );

        if ( itSolver != NULL )
        {
            itSolver->setStoppingCriterion( crit );
        }
        else
        {
            cout << "Not iterative solver, stopping criterion ignored" << endl;
        }

        // Allow individual settings for GMRES solver

        // ToDo: enable GMRES

        /*
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

        */

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

                DenseVector<ValueType> compSolution;

                compSolution.readFromFile( finalSolutionFilename );
                compSolution.redistribute( solution.getDistributionPtr() );
                compSolution -= solution;
                RealType<ValueType> maxDiff = compSolution.maxNorm();
                HOST_PRINT( myRank, "Maximal difference between solution in " << finalSolutionFilename << ": " << maxDiff )
            }
            else
            {
                HOST_PRINT( myRank, "Write solution to output file " << finalSolutionFilename )
                LamaTiming timer( comm, "Writing solution" );
                solution.writeToFile( finalSolutionFilename );
            }
        }

        HOST_PRINT( myRank, "lamaSolver finished, solver = " << *mySolver << ", A = " << matrix )
    }
    catch ( common::Exception& e )
    {
        HOST_PRINT( myRank, "Terminate due to error: " << e.what() )
    }
}
