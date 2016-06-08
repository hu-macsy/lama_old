/**
 * @file solver/examples/solver/amg_solver.cpp
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
 * @brief Example driver program for CG solver
 * @author Thomas Brandes
 * @date 02.07.2012
 */

#include "LamaConfig.hpp"

#include <scai/lama/matrix/all.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/lama/expression/all.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/SimpleAMG.hpp>
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>

using namespace scai;
using namespace dmemo;
using namespace lama;
using namespace solver;

using namespace std;

using common::Walltime;

typedef RealType ValueType;

void dummy( const LamaConfig& lamaconf )
{
    CSRSparseMatrix<ValueType> m;
    m.setIdentity( 100 );
    DenseVector<ValueType> v1( 100, 1.0  );
    m.setContextPtr( lamaconf.getContextPtr() );
    m.setCommunicationKind( lamaconf.getCommunicationKind() );
    v1.setContextPtr( lamaconf.getContextPtr() );
    DenseVector<ValueType> v2( m *  v1 );
}

/** Read in a partitioning file for the input data if available. */

vector<IndexType>* readPartitionVector( const char* filename, int commSize, int expectedSize )
{
    // read in general distribution if available for number of processors

    ostringstream partitionFileName;

    partitionFileName << filename << "." << commSize ;

    ifstream pfile ( partitionFileName.str().c_str() );

    if ( pfile.is_open() )
    {
        vector<IndexType> *pvector = new vector<IndexType>;

        pvector->reserve( expectedSize );

        IndexType elem;

        while ( pfile >> elem )
        {
            pvector->push_back (elem);
        }

        cout << "Read partitioning file " << partitionFileName.str()
             << ", actual size = " << pvector->size() << ", expected = " << expectedSize << endl;

        pfile.close();

        return pvector;
    }
    else
    {
        cout << "Could not open any partitioning file " << partitionFileName.str() << endl;

        return NULL;
    }
}

int main( int argc, char* argv[] )
{
    LamaConfig lamaconf;

    // Get (default) communicator, will be MPI if available

    const CommunicatorPtr& comm = lamaconf.getCommunicatorPtr();

    int myRank   = comm->getRank();
    int numProcs = comm->getSize();

    string filename;

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

    MatrixPtr matrixPtr( lamaconf.getMatrix() );
    VectorPtr rhsPtr( matrixPtr->newDenseVector() );

    Matrix& matrix = *matrixPtr;
    Vector& rhs = *rhsPtr;

    // Each processor should print its configuration

    cout << lamaconf << endl;

    double start = Walltime::get();   // start timing of reading

    // read matrix + rhs from disk

    matrix.readFromFile( filename + ".frm" );
    rhs.readFromFile( filename + ".frv" );

    // only square matrices are accetpted

    SCAI_ASSERT_EQUAL( matrix.getNumRows(), matrix.getNumColumns(), "size mismatch" )
    SCAI_ASSERT_EQUAL( matrix.getNumRows(), rhs.size(), "size mismatch" )

    int numRows = matrix.getNumRows();

    // for solution create vector with same format/type as rhs, size = numRows, init = 0.0

    scai::common::unique_ptr<Vector> solutionPtr( rhs.newVector() );
    Vector& solution = *solutionPtr;

    solution.allocate( rhs.getDistributionPtr() );
    solution = 0.0;   // intialize of a vector

    double stop = Walltime::get();  // stop timing for reading

    if ( myRank == 0 )
    {
        cout << "loading data took " << stop - start << " secs." << endl;
    }

    if ( numProcs > 1 )
    {
        // determine a new distribution so that each processor gets part of the matrix according to its weight

        float weight = lamaconf.getWeight();

        DistributionPtr dist;

        if ( lamaconf.useMetis() )
        {
            // maybue there is already a partition vector storred as file

            common::shared_ptr<vector<IndexType> > mapVector;

            mapVector.reset( readPartitionVector( filename.c_str(), numProcs, numRows ) );

            if ( mapVector )
            {
                dist.reset( new GeneralDistribution( *mapVector, numRows, comm ) );
            }
            else
            {
                start = Walltime::get();   // start timing of Metis

                dist.reset( Distribution::getDistributionPtr( "METIS", lamaconf.getCommunicatorPtr(), matrix, weight ) );
                stop = Walltime::get();   // stop timing of Metis

                if ( myRank == 0 )
                {
                    cout << "Metis call took " << stop - start << " secs." << endl;
                }
            }
        }
        else
        {
            dist.reset( new GenBlockDistribution( numRows, weight, lamaconf.getCommunicatorPtr() ) );
        }

        start = Walltime::get();   // start timing of redistribution

        matrix.redistribute( dist, dist );
        rhs.redistribute ( dist );
        solution.redistribute ( dist );

        cout << comm << ": matrix = " << matrix ;

        stop = Walltime::get();   // stop timing of redistribution

        if ( myRank == 0 )
        {
            cout << "redistributing data took " << stop - start << " secs." << endl;
        }
    }

    double matrixSize  = matrix.getMemoryUsage() / 1024.0 / 1024.0;

    if ( myRank == 0 )
    {
        cout << "Matrix Size = " << matrixSize << " MB" << endl;
    }

    start = Walltime::get();  // start time of data transfer

    dummy( lamaconf );

    matrix.setCommunicationKind( lamaconf.getCommunicationKind() );
    matrix.setContextPtr( lamaconf.getContextPtr() );
    rhs.setContextPtr( lamaconf.getContextPtr() );
    solution.setContextPtr( lamaconf.getContextPtr() );

    rhs.prefetch();
    matrix.prefetch();
    matrix.wait();
    rhs.wait();

    stop = Walltime::get();

    if ( myRank == 0 )
    {
        cout << "prefetching data took " << stop - start << " secs." << endl;
    }

    // setting up solver from file "solveconfig.txt"

    std::ostringstream loggerName;

    loggerName << "<CG>, " << lamaconf.getCommunicator() << ": ";

    LoggerPtr logger( new CommonLogger ( loggerName.str(), LogLevel::advancedInformation,
                                         LoggerWriteBehaviour::toConsoleOnly ) );

    CG mySolver( "CGSolver", logger );

    Scalar relEps = 1e-14;
    Scalar absEps = 1e-16;

    NormPtr relNorm = NormPtr( new L2Norm() );
    NormPtr absNorm = NormPtr( new L2Norm() );

    CriterionPtr relRT( new ResidualThreshold( relNorm, relEps, ResidualThreshold::Relative ) );
    CriterionPtr absRT( new ResidualThreshold( absNorm, absEps, ResidualThreshold::Absolute ) );

    CriterionPtr rt( new Criterion( relRT, absRT, Criterion::OR ) );

    if ( lamaconf.hasMaxIter() )
    {
        CriterionPtr it( new IterationCount( lamaconf.getMaxIter() ) );

        // stop if iteration count reached OR residual threshold is reached

        rt.reset( new Criterion ( it, rt, Criterion::OR ) );
    }

    mySolver.setStoppingCriterion( rt );

    LoggerPtr amgLogger( new CommonLogger ( loggerName.str(), LogLevel::solverInformation,
                                            LoggerWriteBehaviour::toConsoleOnly ) );

    common::shared_ptr<SimpleAMG> amgSolver( new SimpleAMG( "SimpleAMG solver", amgLogger ) );

    amgSolver->setHostOnlyLevel( 4 );
    amgSolver->setReplicatedLevel( 5 );
    amgSolver->setMaxLevels( 25 );
    amgSolver->setMinVarsCoarseLevel( 200 );

    mySolver.setPreconditioner( amgSolver );

    start = Walltime::get();

    // initialize solver

    mySolver.initialize( matrix );

    stop = Walltime::get();

    if ( myRank == 0 )
    {
        cout << "Solver setup took " << stop - start << " secs." << endl;
    }

    dummy( lamaconf );

    // run solver

    start = Walltime::get();

    mySolver.solve( solution, rhs );

    stop = Walltime::get();

    if ( myRank == 0 )
    {
        cout << "Solution phase took " << stop - start << " secs." << endl;
    }

    bool writeFlag = false;

    if ( writeFlag )
    {
        start = Walltime::get();
        solution.writeToFile( "CG_solution" );
        stop = Walltime::get();

        if ( myRank == 0 )
        {
            cout << "Writing solution: " << stop - start << " secs." << endl;
        }
    }
}

