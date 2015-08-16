/**
 * @file cg_solver.cpp
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
 * @brief Example driver program for CG solver
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.0.0
 */

#include "LamaConfig.hpp"

#include <scai/lama/matrix/all.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/distribution/NoDistribution.hpp>
#include <scai/lama/distribution/GeneralDistribution.hpp>
#include <scai/lama/distribution/GenBlockDistribution.hpp>
#include <scai/lama/expression/all.hpp>

#include <scai/lama/solver/CG.hpp>
#include <scai/lama/solver/SimpleAMG.hpp>
#include <scai/lama/solver/logger/CommonLogger.hpp>
#include <scai/lama/solver/criteria/ResidualThreshold.hpp>
#include <scai/lama/solver/criteria/IterationCount.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/unique_ptr.hpp>

using namespace std;
using namespace scai::lama;

using scai::common::Walltime;
using scai::common::unique_ptr;

void dummy( const LamaConfig& lamaconf )
{
    CSRSparseMatrix<double> m;
    m.setIdentity( 100 );
    DenseVector<double> v1( 100, 1.0  );
    m.setContext( lamaconf.getContextPtr() );
    m.setCommunicationKind( lamaconf.getCommunicationKind() );
    v1.setContext( lamaconf.getContextPtr() );
    DenseVector<double> v2( m *  v1 );
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
        vector<IndexType> *pvector = new vector<scai::lama::IndexType>;

        pvector->reserve( expectedSize );

        scai::lama::IndexType elem;

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

    // Each processor should print its configuration

    cout << lamaconf << endl;

    double start = Walltime::get();   // start timing of reading

    // read matrix + rhs from disk

    matrix.readFromFile( filename );
    rhs.readFromFile( filename );

    // only square matrices are accetpted

    SCAI_ASSERT_EQUAL( matrix.getNumRows(), matrix.getNumColumns() )
    SCAI_ASSERT_EQUAL( matrix.getNumRows(), rhs.size() )

    int numRows = matrix.getNumRows();

    // for solutin create vector with same format/type as rhs, size = numRows, init = 0.0

    unique_ptr<Vector> solutionPtr( rhs.clone( rhs.getDistributionPtr() ) );
    Vector& solution = *solutionPtr;

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

            mapVector.reset( readPartitionVector( filename, numProcs, numRows ) );

            if ( mapVector )
            {
                dist.reset( new GeneralDistribution( *mapVector, numRows, comm ) );
            }
            else
            {
                start = Walltime::get();   // start timing of Metis

                dist.reset( Distribution::getDistribution( "METIS", lamaconf.getCommunicatorPtr(), matrix, weight ) );
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
    matrix.setContext( lamaconf.getContextPtr() );
    rhs.setContext( lamaconf.getContextPtr() );
    solution.setContext( lamaconf.getContextPtr() );

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

    LoggerPtr logger( new CommonLogger ( loggerName.str(), scai::lama::LogLevel::advancedInformation,
                   LoggerWriteBehaviour::toConsoleOnly,
                    new Timer() ) );

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

    LoggerPtr amgLogger( new CommonLogger ( loggerName.str(), scai::lama::LogLevel::solverInformation,
                   LoggerWriteBehaviour::toConsoleOnly,
                   new Timer() ) );

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

