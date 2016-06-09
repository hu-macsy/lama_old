/**
 * @file solver/examples/solver/solver.cpp
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
#include <scai/lama/norm/L2Norm.hpp>

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

    // only one argument for the filename should remain

    if ( argc < 1 || argc > 3 )
    {
        if ( myRank == 0 )
        {
            CONFIG_ERROR( "Illegal number of arguments" )
        }

        exit( 1 );
    }

    const char* matrix_filename = argv[1];
    const char* vector_filename = argc <= 2 ? argv[1] : argv[2];

    // use auto pointer so that matrix will be deleted at program exit

    scai::common::unique_ptr<Matrix> matrixPtr( lamaconf.getMatrix() );
    scai::common::unique_ptr<Vector> rhsPtr( matrixPtr->newDenseVector() );

    Matrix& matrix = *matrixPtr;
    Vector& rhs = *rhsPtr;

    CSRSparseMatrix<ValueType> inMatrix;
    // Each processor should print its configuration
    cout << lamaconf << endl;
    {
        LamaTiming timer( comm, "Loading data" );
        // read matrix + rhs from disk
        inMatrix.readFromFile( matrix_filename );
        std::cout << "Matrix from file " << matrix_filename << " : " << inMatrix << std::endl;

        try
        {
            rhs.readFromFile( vector_filename );
        }
        catch ( const std::exception& )
        {
            std::cout << "reading vector from file " << vector_filename << " failed, take sum( Matrix, 2 ) " << std::endl;
            {
                scai::common::unique_ptr<Vector> xPtr( rhs.newVector() );
                Vector& x = *xPtr;
                x.allocate( inMatrix.getColDistributionPtr() );
                x = Scalar( 1 );
                rhs.allocate( inMatrix.getRowDistributionPtr() );
                rhs = inMatrix * x;
            }
        }

        // only square matrices are accetpted
        SCAI_ASSERT_EQUAL( inMatrix.getNumRows(), inMatrix.getNumColumns(), "size mismatch" )
        SCAI_ASSERT_EQUAL( inMatrix.getNumRows(), rhs.size(), "size mismatch" )
    }
    // for solutin create vector with same format/type as rhs, size = numRows, init = 0.0
    scai::common::unique_ptr<Vector> solutionPtr( Vector::create( rhs.getCreateValue() ) );
    Vector& solution = *solutionPtr;
    int numRows = inMatrix.getNumRows();
    solution.allocate( inMatrix.getColDistributionPtr() );
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
        LamaTiming timer( comm, "Type conversion from CSR<ValueType> to target format" );
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
    scai::common::unique_ptr<Solver> mySolver( Solver::create( lamaconf.getSolverName(), solverName.str() ) );
    IterativeSolver* itSolver = dynamic_cast<IterativeSolver*>( mySolver.get() );
    SCAI_ASSERT( itSolver, "Not an iterative solver: " << *mySolver )

    mySolver->setLogger( logger );

    CriterionPtr crit;

    NormPtr norm( new L2Norm() );

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
        std::cout << "No criterion set, take default" << std::endl;

        eps = 0.001;

        crit.reset( new ResidualThreshold( norm, eps, ResidualThreshold::Absolute ) );
    }

    itSolver->setStoppingCriterion( crit );

    GMRES* gmresSolver = dynamic_cast<GMRES*>( mySolver.get() );

    if ( gmresSolver )
    {
        // here is the possibility to set solver specific values

        int dim = 5;
 
        common::Settings::getEnvironment( dim, "SCAI_KRYLOV_DIM" );

        std::cout << "GMRES solver, krylov dim = " << dim << std::endl;

        gmresSolver->setKrylovDim( dim );
    }

    SimpleAMG* amgSolver = dynamic_cast<SimpleAMG*>( mySolver.get() );

    if ( amgSolver )
    {
        amgSolver->setHostOnlyLevel( 4 );
        amgSolver->setReplicatedLevel( 5 );
        amgSolver->setMaxLevels( 25 );
        amgSolver->setMinVarsCoarseLevel( 200 );
    }

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
