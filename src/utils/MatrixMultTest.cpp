/**
 * @file MatrixMultTest.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief MatrixMultTest.cpp
 * @author Jan Ecker
 * @date 03.01.2013
 * $Id$
 */

#include <lama/Scalar.hpp>
#include <lama/Context.hpp>
#include <lama/matrix/CSRSparseMatrix.hpp>
#include <lama/expression/MatrixExpressions.hpp>
#include <lama/HostReadAccess.hpp>
#include <lama/storage/CSRStorage.hpp>
#include <lama/matutils/MatrixCreator.hpp>
#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

//using namespace lama;

template<typename ValueType>
double multiply(
    lama::CSRSparseMatrix<ValueType> *matrixA,
    lama::CSRSparseMatrix<ValueType> *matrixB,
    lama::CSRSparseMatrix<ValueType> *matrixC,
    lama::ContextType context,
    const bool prefetch,
    const bool silentFlag )
{
    lama::ContextPtr loc = lama::ContextFactory::getContext( context );

    matrixA->setContext( loc );
    matrixB->setContext( loc );

    if( prefetch )
    {
        if( !silentFlag )
        {
            std::cout << "Prefetching matrices... " << std::flush;
            matrixA->prefetch();
            std::cout << "..." << std::flush;
        }
        matrixB->prefetch();
        if( !silentFlag )
        {
            std::cout << "done!" << std::endl;
        }
    }

    double time;
    double start;

    if( !silentFlag )
    {
        std::cout << "Starting multiplication on " << context << " Context with size " << matrixA->getNumRows() << "x"
                  << matrixB->getNumColumns() << std::endl;
    }

    start = MPI_Wtime();

    *matrixC = 1.0 * *matrixA * *matrixB;

    time = MPI_Wtime() - start;

    if( !silentFlag )
    {
        std::cout << "Multiplication finished, time: " << time << "s" << std::endl << std::endl;
    }
    else
    {
        std::cout << time << std::endl;
    }

    return time;
}

template<typename ValueType>
void maxRow(
    lama::CSRSparseMatrix<ValueType> *matrix,
    const bool silentFlag)
{
    lama::ContextPtr locHost = lama::ContextFactory::getContext( lama::Context::Host );

    matrix->setContext( locHost );
    if( !silentFlag )
    {
        std::cout << "Checking for longest row..." << std::flush;
    }

    lama::CSRStorage<ValueType> localStorage;
    localStorage = matrix->getLocalStorage();

    lama::HostReadAccess<lama::IndexType> ia( localStorage.getIA() );

    lama::IndexType longestRow = 0;
    for( int i = 0; i < matrix->getNumRows(); ++i )
    {
        if ( ia[i+1]-ia[i] > longestRow ){
            longestRow = ia[i+1]-ia[i];
        }
    }
    if( !silentFlag )
    {
        std::cout << std::endl << "finished, longest row has " << longestRow << " Elements!" << std::endl;
    }
}

template<typename ValueType>
void validate(
    lama::CSRSparseMatrix<ValueType> *matrixGiven,
    lama::CSRSparseMatrix<ValueType> *matrixCorrect,
    const bool silentFlag,
    const double epsilon = 1e-100 )
{
    lama::ContextPtr locHost = lama::ContextFactory::getContext( lama::Context::Host );

    matrixGiven->setContext( locHost );
    matrixCorrect->setContext( locHost );

    int numErrors = 0;

    if( !silentFlag )
    {
        std::cout << "Running validation with epsilon = " << epsilon << std::flush;
    }

    lama::CSRStorage<ValueType> localStorage;
    localStorage = matrixCorrect->getLocalStorage();

    lama::HostReadAccess<lama::IndexType> correctIa( localStorage.getIA() );
    lama::HostReadAccess<lama::IndexType> correctJA( localStorage.getJA() );

    if( localStorage.getNumValues() != matrixGiven->getLocalStorage().getNumValues() )
    {
        std::cout << "numValues is different, is " << matrixGiven->getLocalStorage().getNumValues() << " but should "
                  << localStorage.getNumValues() << std::flush;
    }

//#pragma omp parallel for
    for( int i = 0; i < matrixCorrect->getNumRows(); i++ )
    {
        lama::IndexType start = correctIa[i];
        lama::IndexType end = correctIa[i + 1];
        for( lama::IndexType k = start; k < end; k++ )
        {
            lama::IndexType j = correctJA[k];
            lama::Scalar valCorrect = matrixCorrect->getValue( i, j );
            lama::Scalar valGiven = matrixGiven->getValue( i, j );

            if( abs( valGiven.getValue<ValueType>() ) - abs( valCorrect.getValue<ValueType>() ) > epsilon )
            {
                std::cout << "Error in Matrix on position (" << i << ", " << j << ") value is " << valGiven
                          << " but should " << valCorrect << std::endl;
                numErrors++;
            }
        }
    }
    if( numErrors > 0 )
    {
        std::cout << std::endl << "Validation finished, " << numErrors << " errors found!" << std::endl;
    }
    else
    {
        std::cout << std::endl << "Validation finished, no errors found!" << std::endl;
    }
}

int main( int argc, char **argv )
{
    bool validateFlag = false;
    bool prefetchFlag = false;
    bool silentFlag = false;
    bool benchmarkFlag = false;
    bool randomFlag = false;
    int size = 1000;
    float density = 1.0;
    std::string inputMatrix;
    int matrix = -1;
    int parameter;

    while( ( parameter = getopt( argc, argv, "vpm:i:sbrg:d:" ) ) != -1 )
    {
        switch( parameter )
        {
        case 'v':
            validateFlag = true;
            break;
        case 'p':
            prefetchFlag = true;
            break;
        case 'm':
            matrix = atoi( optarg );
            break;
        case 'i':
            inputMatrix = optarg;
            break;
        case 's':
            silentFlag = true;
            break;
        case 'b':
            benchmarkFlag = true;
            break;
        case 'r':
            randomFlag = true;
            break;
        case 'g':
            size = atoi( optarg );
            break;
        case 'd':
            density = atof( optarg );
            break;
        default:
            abort();
            break;
        }
    }

    std::string matrixFile;

    switch( matrix )
    {
    case 0:
        matrixFile = "../../../res/benchfiles/bcspwr01.mtx";
        break;
    case 1:
        matrixFile = "../../../res/benchfiles/bcsstk01.mtx";
        break;
    case 2:
        matrixFile = "../../../res/benchfiles/bcsstk18.mtx";
        break;
    case 3:
        matrixFile = "../../../res/benchfiles/e40r5000.mtx";
        break;
    case 4:
        matrixFile = "../../../res/benchfiles/impcol_b.mtx";
        break;
    case 5:
        matrixFile = "../../../res/benchfiles/mc2depi.mtx";
        break;
    case 6:
        matrixFile = "../../../res/benchfiles/s3dkt3m2.mtx";
        break;
    case 7:
        matrixFile = "../../../res/benchfiles/tols90.mtx";
        break;
    case 8:
        matrixFile = "../../../res/benchfiles/amazon-2008.mtx";
        break;
    default:
        matrixFile = inputMatrix;
        break;
    }

    if( !silentFlag )
    {
        std::cout << "matrixMultTest running" << std::endl << std::endl;
        std::cout << "Loading matrix \"" << matrixFile << "\"... " << std::flush;
    }

    lama::CSRSparseMatrix<double> matrixA;
    lama::CSRSparseMatrix<double> matrixB;
    if ( randomFlag ) {
        matrixA = lama::CSRSparseMatrix<double>( size, size );
        matrixB = lama::CSRSparseMatrix<double>( size, size );
        lama::MatrixCreator<double>::fillRandom( matrixA, density );
        lama::MatrixCreator<double>::fillRandom( matrixB, density );
    } else {
        matrixA = lama::CSRSparseMatrix<double>( inputMatrix );
        matrixB =  lama::CSRSparseMatrix<double>( matrixA );
    }

    float alloc;
    if( !silentFlag )
    {
        std::cout << "done!" << std::endl << std::endl;
        std::cout << "Matrix A info: " << std::endl;
        std::cout << "Size: \t\t" << matrixA.getNumRows() << "x" << matrixA.getNumColumns() << std::endl;
        std::cout << "NNZ: \t\t" << matrixA.getNumValues() << std::endl;
        std::cout << "NNZ/Row: \t" << matrixA.getNumValues() / matrixA.getNumRows() << std::endl;
        std::cout << "Density: \t" << (float)matrixA.getNumValues() / (float)(matrixA.getNumRows() * matrixA.getNumRows()) << std::endl;
        std::cout << "Filesize: \t" << ( ( matrixA.getMemoryUsage() / 1024.0 ) / 1024.0 ) << " MB" << std::endl
                  << std::endl;

        std::cout << "Matrix B info: " << std::endl;
        std::cout << "Size: \t\t" << matrixB.getNumRows() << "x" << matrixB.getNumColumns() << std::endl;
        std::cout << "NNZ: \t\t" << matrixB.getNumValues() << std::endl;
        std::cout << "NNZ/Row: \t" << matrixB.getNumValues() / matrixB.getNumRows() << std::endl;
        std::cout << "Density: \t" << (float)matrixB.getNumValues() / (float)(matrixB.getNumRows() * matrixB.getNumRows()) << std::endl;
        std::cout << "Filesize: \t" << ( ( matrixB.getMemoryUsage() / 1024.0 ) / 1024.0 ) << " MB" << std::endl
                  << std::endl;

        float density = (float)matrixA.getNumValues() / ( (float)matrixA.getNumRows() * (float)matrixA.getNumColumns() );
        float prop    = 1.0 - pow( ( 1.0 - density ), 2 * matrixA.getNumColumns() );
        alloc   = prop * matrixA.getNumRows() * matrixA.getNumColumns();
        std::cout << "prop: " << prop << std::endl << std::endl;
        std::cout << "Alloc: " << alloc << std::endl << std::endl;

    }

    lama::CSRSparseMatrix<double> matrixCHost;
    lama::CSRSparseMatrix<double> matrixCCuda;

    if( benchmarkFlag )
    {
        double time1 = multiply<double>( &matrixA, &matrixB, &matrixCCuda, lama::Context::CUDA, prefetchFlag,
                                         silentFlag );
        double time2 = multiply<double>( &matrixA, &matrixB, &matrixCHost, lama::Context::Host, prefetchFlag,
                                         silentFlag );

        if( !silentFlag )
        {
            std::cout << "Speedup: " << time2 / time1 << std::endl << std::endl;

        }
    }
    else
    {
        multiply<double>( &matrixA, &matrixB, &matrixCCuda, lama::Context::CUDA, prefetchFlag, silentFlag );
    }

    if( !silentFlag )
    {
        std::cout << "Done!" << std::endl << std::endl;

        std::cout << "Matrix info: " << std::endl;
        std::cout << "Size: \t\t" << matrixCCuda.getNumRows() << "x" << matrixCCuda.getNumColumns() << std::endl;
        std::cout << "NNZ: \t\t" << matrixCCuda.getNumValues() << std::endl;
        std::cout << "NNZ/Row: \t" << matrixCCuda.getNumValues() / matrixCCuda.getNumRows() << std::endl;
        std::cout << "Filesize: \t" << ( ( matrixCCuda.getMemoryUsage() / 1024.0 ) / 1024.0 ) << " MB" << std::endl
                  << std::endl;

        std::cout << "Abweichung: " << (1 - (float)matrixCCuda.getNumValues() / alloc) * 100 << "%" << std::endl;

        maxRow ( &matrixCCuda, silentFlag );
    }

    if( validateFlag )
    {
        validate<double>( &matrixCCuda, &matrixCHost, silentFlag );
    }

    return 0;
}
