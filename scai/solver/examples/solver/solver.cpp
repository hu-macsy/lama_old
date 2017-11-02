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
#include <scai/dmemo/Partitioning.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/lama/norm/L2Norm.hpp>
#include <scai/lama/norm/L1Norm.hpp>
#include <scai/lama/norm/MaxNorm.hpp>
#include <scai/lama/StorageIO.hpp>

#include <scai/solver/GMRES.hpp>
#include <scai/solver/SimpleAMG.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

#include <memory>

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

    // accept only 2 - 4 arguments

    if ( argc < 2  || argc > 4 )
    {
        if ( myRank == 0 )
        {
            CONFIG_ERROR( "Illegal number of arguments" )
        }

        return -1;
    }

    const char* matrix_filename = argv[1];
    const char* rhs_filename    = argc <= 2 ? argv[1] : argv[2];
    const char* sol_filename    = argc <= 3 ? NULL : argv[3];

    // use auto pointer so that matrix will be deleted at program exit

    std::unique_ptr<Matrix> matrixPtr( lamaconf.getMatrix() );
    std::unique_ptr<Vector> rhsPtr( matrixPtr->newDenseVector() );

    _Matrix& matrix = *matrixPtr;
    Vector& rhs = *rhsPtr;

    // input matrix will be CSR format

<<<<<<< HEAD
    scai::common::unique_ptr<Matrix> in_MatrixPtr( _Matrix::getMatrix( _Matrix::CSR, lamaconf.getValueType() ) );
    _Matrix& inMatrix = *in_MatrixPtr;
=======
    std::unique_ptr<Matrix> inMatrixPtr( Matrix::getMatrix( Matrix::CSR, lamaconf.getValueType() ) );
    Matrix& inMatrix = *inMatrixPtr;
>>>>>>> lama_intern/feature/remove_boost

    // Each processor should print its configuration

    cout << lamaconf << endl;

    {
        LamaTiming timer( comm, "Loading data" );

        // read matrix + rhs from disk

        inMatrix.readFromFile( matrix_filename );
        cout << "Matrix from file " << matrix_filename << " : " << inMatrix << endl;

        try
        {
            rhs.readFromFile( rhs_filename );
            cout << "rhs from file " << rhs_filename << " : " << rhs << endl;
        }
        catch ( const exception& )
        {
            cout << "reading vector from file " << rhs_filename << " failed, take sum( _Matrix, 2 ) " << endl;
            {
                std::unique_ptr<Vector> xPtr( rhs.newVector() );
                Vector& x = *xPtr;
                x.allocate( inMatrix.getColDistributionPtr() );
                x = Scalar( 1 );
                rhs = inMatrix * x;
            }
        }

        // only square matrices are accetpted
        SCAI_ASSERT_EQUAL( inMatrix.getNumRows(), inMatrix.getNumColumns(), "size mismatch" )
        SCAI_ASSERT_EQUAL( inMatrix.getNumRows(), rhs.size(), "size mismatch" )
    }

    // for solution create vector with same format/type as rhs, size = numRows, init = 0.0

    std::unique_ptr<Vector> solutionPtr( Vector::create( rhs.getCreateValue() ) );
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

        CommunicatorPtr comm = lamaconf.getCommunicatorPtr();
        DistributionPtr dist;

        if ( lamaconf.useMetis() )
        {
            LamaTiming timer( comm, "Metis" );

            PartitioningPtr partitioning( Partitioning::create( "METIS" ) );

            SCAI_ASSERT_ERROR( partitioning.get(), "METIS partitioning not available" )

            dist = partitioning->partitionIt( comm, inMatrix, weight );
        }
        else
        {
            // same as using Partitioning::create( "BLOCK" )

            dist.reset( new GenBlockDistribution( numRows, weight, comm ) );
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

    // setting up solver from file "solveconfig.txt"
    ostringstream solverName;
    solverName << "<" << lamaconf.getSolverName() << ">";
    ostringstream loggerName;
    loggerName << solverName.str() << ", " << lamaconf.getCommunicator() << ": ";
    LoggerPtr logger( new CommonLogger ( loggerName.str(),
                                         lamaconf.getLogLevel(),
                                         LoggerWriteBehaviour::toConsoleOnly ) );
    std::unique_ptr<Solver> mySolver( Solver::create( lamaconf.getSolverName(), solverName.str() ) );
    IterativeSolver* itSolver = dynamic_cast<IterativeSolver*>( mySolver.get() );
    SCAI_ASSERT( itSolver, "Not an iterative solver: " << *mySolver )

    mySolver->setLogger( logger );

    CriterionPtr crit;

    NormPtr norm;

    if ( lamaconf.getNorm() == "L1" )
    {
        norm.reset( new L1Norm() );
    }
    else if ( lamaconf.getNorm() == "Max" )
    {
        norm.reset( new MaxNorm() );
    }
    else
    {
        norm.reset( new L2Norm() );
    }

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
        cout << "No criterion set, take default" << endl;

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

        cout << "GMRES solver, krylov dim = " << dim << endl;

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

    {
        LamaTiming timer( comm, "Solver setup" );
        mySolver->initialize( matrix );
    }
    {
        LamaTiming timer( comm, "Solver solve" );
        mySolver->solve( solution, rhs );
    }

    // if 3rd argument for solution file is specified, write it or compare it

    if ( sol_filename != NULL )
    {
        if ( _StorageIO::fileExists( sol_filename ) )
        {
            cout << "Compare solution with vector in " << sol_filename << endl;
            LamaTiming timer( comm, "Comparing solution" );
            std::unique_ptr<Vector> compSolutionPtr( Vector::create( rhs.getCreateValue() ) );
            Vector& compSolution = *compSolutionPtr;
            compSolution.readFromFile( sol_filename );
            compSolution.redistribute( solution.getDistributionPtr() );
            compSolution -= solution;
            Scalar maxDiff = compSolution.maxNorm();
            cout << "Maximal difference between solution in " << sol_filename << ": " << maxDiff << endl;
        }
        else
        {
            cout << "Write solution to output file " << sol_filename << ".mtx (Matrix Market)" << endl;
            LamaTiming timer( comm, "Writing solution" );
            // replicate solution on all processors
            DistributionPtr repDist( new NoDistribution( solution.size() ) );
            solution.redistribute( repDist );
            if ( myRank == 0 )
            {
                solution.writeToFile( sol_filename, File::MATRIX_MARKET );
            }
        }
    }
}
