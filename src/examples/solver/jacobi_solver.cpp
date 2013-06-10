/**
 * @file jacobi_solver.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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

#include "LAMAConfig.hpp"

#include <lama/matrix/all.hpp>

#include <lama/DenseVector.hpp>
#include <lama/CommunicatorFactory.hpp>
#include <lama/distribution/BlockDistribution.hpp>
#include <lama/distribution/NoDistribution.hpp>

#include <lama/solver/CG.hpp>
#include <lama/solver/SpecializedJacobi.hpp>
#include <lama/solver/logger/CommonLogger.hpp>
#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/norm/L2Norm.hpp>

#include <lama/Walltime.hpp>

using namespace std;
using namespace lama;

/** ValueType is the type used for matrix and vector elements. */
typedef double ValueType;

typedef DenseVector<ValueType> VectorType;

SparseMatrix<ValueType>* createMatrix( const char* type )
{
    if ( strcmp( type, "CSR" ) == 0 )
    {
        return new CSRSparseMatrix<ValueType>();
    }
    else if ( strcmp( type, "ELL" ) == 0 )
    {
        return new ELLSparseMatrix<ValueType>();
    }
    else if ( strcmp( type, "JDS" ) == 0 )
    {
        return new JDSSparseMatrix<ValueType>();
    }
    else if ( strcmp( type, "DIA" ) == 0 )
    {
        return new DIASparseMatrix<ValueType>();
    }
    else
    {
        return NULL;
    }
}

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

    auto_ptr<SparseMatrix<ValueType> > matrixPtr( createMatrix( lamaconf.getFormat() ) );

    if ( myRank == 0 )
    {
        cout << "mContext = " << *lamaconf.mContext << endl;
        cout << lamaconf << endl;
    }

    double start = Walltime::get();   // start timing of reading

    SparseMatrix<ValueType>& matrix = *matrixPtr;

    // read matrix + rhs from disk

    matrix.readFromFile( filename );

    VectorType rhs( filename );

    // only square matrices are accetpted

    LAMA_ASSERT_EQUAL( matrix.getNumRows(), matrix.getNumColumns() )
    LAMA_ASSERT_EQUAL( matrix.getNumRows(), rhs.size() )

    int numRows = matrix.getNumRows();
    VectorType solution( numRows, 0.0 );

    double stop = Walltime::get();  // stop timing for reading

    if ( myRank == 0 )
    {
        cout << "loading data took " << stop - start << " secs." << endl;
    }

    start = Walltime::get();   // start timing of redistribution

    // distribute data (trivial block partitioning)

    DistributionPtr dist( new BlockDistribution( numRows, lamaconf.getCommunicatorPtr() ) );

    matrix.redistribute( dist, dist );

    rhs.redistribute ( dist );
    solution.redistribute ( dist );

    stop = Walltime::get();   // stop timing of redistribution

    if ( myRank == 0 )
    {
        cout << "redistributing data took " << stop - start << " secs." << endl;
    }

    double matrixSize  = matrix.getMemoryUsage() / 1024.0 / 1024.0;

    if ( myRank == 0 )
    {
        cout << "Matrix Size = " << matrixSize << " MB" << endl;
    }

    start = Walltime::get();  // start time of data transfer

    matrix.setCommunicationKind( lamaconf.mCommunicationKind );
    matrix.setContext( lamaconf.mContext );
    rhs.setContext( lamaconf.mContext );
    solution.setContext( lamaconf.mContext );

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

    LoggerPtr logger( new CommonLogger ( loggerName.str(), LogLevel::solverInformation,
                   LoggerWriteBehaviour::toConsoleOnly,
                   std::auto_ptr<Timer>( new Timer() ) ) );

    SpecializedJacobi mySolver( "CGSolver", logger );

    Scalar eps = 0.00001;
    NormPtr norm = NormPtr( new L2Norm() );

    CriterionPtr rt( new IterationCount( 100 ) );

    mySolver.setStoppingCriterion( rt );

    start = Walltime::get();

    // initialize solver

    mySolver.initialize( matrix );

    stop = Walltime::get();

    if ( myRank == 0 )
    {
        cout << "Solver setup took " << stop - start << " secs." << endl;
    }

    // run solver

    start = Walltime::get();

    mySolver.solve( solution, rhs );

    stop = Walltime::get();

    if ( myRank == 0 )
    {
        cout << "Solution phase took " << stop - start << " secs." << endl;
    }

    start = Walltime::get();

    // solution.writeToFile( "CG_solution" );

    stop = Walltime::get();

    if ( myRank == 0 )
    {
        cout << "Writing solution: " << stop - start << " secs." << endl;
    }
}

