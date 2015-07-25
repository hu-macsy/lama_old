/**
 * @file OpenMPCSRUtils.cpp
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
 * @brief Implementation of CSR utilities with OpenMP
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.0.0
 */

// for dll_import
#include <lama/openmp/OpenMPCSRUtils.hpp>
#include <lama/openmp/OpenMP.hpp>

// others
#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

// trace
#include <tracing/tracing.hpp>

// boost
#include <common/bind.hpp>
#include <boost/scoped_array.hpp>
#include <boost/preprocessor.hpp>

#include <lama/macros/unused.hpp>
#include <vector>

namespace lama
{

using std::abs;
// is used for abs( float ), abs( double )

using common::getScalarType;

LAMA_LOG_DEF_LOGGER( OpenMPCSRUtils::logger, "OpenMP.CSRUtils" )

/** Number of minimal threads for which parallelization is effective. */

static int minThreads = 3;

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::scanSerial( IndexType array[], const IndexType numValues )
{
    LAMA_LOG_DEBUG( logger, "scanSerial: " << numValues << " entries" )

    // In this case we do it just serial, probably faster

    IndexType runningSum = 0;

    for( IndexType i = 0; i < numValues; i++ )
    {
        IndexType tmp = runningSum;
        runningSum += array[i];
        LAMA_LOG_TRACE( logger, "scan, row = " << i << ", size = " << array[i] << ", offset = " << runningSum )
        array[i] = tmp;
    }

    return runningSum;;
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::scanParallel( PartitionId numThreads, IndexType array[], const IndexType numValues )
{
    // std::cout << "Scan with " << numThreads << " in parallel" << std::endl;

    // For more threads, we do it in parallel
    // Attention: MUST USE schedule(static)

    boost::scoped_array<IndexType> threadCounter( new IndexType[numThreads] );

    LAMA_LOG_DEBUG( logger, "scanParallel: " << numValues << " entries for " << numThreads << " threads" )

    #pragma omp parallel
    {
        IndexType myCounter = 0;

        #pragma omp for schedule(static)

        for( IndexType i = 0; i < numValues; i++ )
        {
            myCounter += array[i];
        }

        threadCounter[omp_get_thread_num()] = myCounter;
    }

    IndexType runningSum = scanSerial( threadCounter.get(), numThreads );

    // Each thread sets now its offsets

    #pragma omp parallel
    {
        IndexType myRunningSum = threadCounter[omp_get_thread_num()];

        #pragma omp for schedule(static)

        for( IndexType i = 0; i < numValues; i++ )
        {
            IndexType tmp = myRunningSum;
            myRunningSum += array[i];
            array[i] = tmp;
        }
    }

    return runningSum;;
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::scan( IndexType array[], const IndexType numValues )
{
    int numThreads = 1; // will be set to available threads in parallel region

    #pragma omp parallel
    #pragma omp master
    {
        numThreads = omp_get_num_threads();
    }

    LAMA_LOG_INFO( logger, "scan " << numValues << " entries, #threads = " << numThreads )

    if( numThreads < minThreads )
    {
        return scanSerial( array, numValues );
    }
    else
    {
        return scanParallel( numThreads, array, numValues );
    }
}

/* --------------------------------------------------------------------------- */

bool OpenMPCSRUtils::validOffsets( const IndexType array[], const IndexType n, const IndexType total )
{
    LAMA_LOG_INFO( logger, "check offset array[ " << n << "] for validity, total = " << total )

    bool validFlag = true;

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE) reduction( && : validFlag )

    for( IndexType i = 0; i < n; i++ )
    {
        if( !( array[i] <= array[i + 1] ) )
        {
            LAMA_LOG_DEBUG( logger, "illegal offsets at pos " << i << ": " << array[i] << " - " << array[i+1] )
            validFlag = false;
        }
    }

    if( array[n] != total )
    {
        LAMA_LOG_DEBUG( logger, "last entry in offset array = " << array[n] << " illegal, should be " << total )
        validFlag = false;
    }

    return validFlag;
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::sizes2offsets( IndexType array[], const IndexType numValues )
{
    IndexType totalValues = scan( array, numValues );
    array[numValues] = totalValues;

    LAMA_LOG_INFO( logger, "sizes2offsets, #values = " << numValues << ", total = " << totalValues )

    return totalValues;
}

/* --------------------------------------------------------------------------- */

void OpenMPCSRUtils::offsets2sizes( IndexType sizes[], const IndexType offsets[], const IndexType numRows )
{
    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

    for( IndexType i = 0; i < numRows; i++ )
    {
        sizes[i] = offsets[i + 1] - offsets[i];
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPCSRUtils::offsets2sizesGather(
    IndexType sizes[],
    const IndexType offsets[],
    const IndexType rowIndexes[],
    const IndexType numRows )
{
    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

    for( IndexType i = 0; i < numRows; i++ )
    {
        IndexType row = rowIndexes[i];
        sizes[i] = offsets[row + 1] - offsets[row];
    }
}

/* --------------------------------------------------------------------------- */

bool OpenMPCSRUtils::hasDiagonalProperty(
    const IndexType numDiagonals,
    const IndexType csrIA[],
    const IndexType csrJA[] )
{
    LAMA_LOG_INFO( logger, "hasDiagonalProperty, #numDiagonals = " << numDiagonals )

    bool diagonalProperty = true;

    #pragma omp parallel for reduction( && : diagonalProperty )

    for( IndexType i = 0; i < numDiagonals; ++i )
    {
        if( !diagonalProperty )
        {
            continue;
        }

        if( csrIA[i] == csrIA[i + 1] )
        {
            diagonalProperty = false;
        }
        else if( csrJA[csrIA[i]] != i )
        {
            diagonalProperty = false;
        }
    }

    LAMA_LOG_DEBUG( logger, "hasDiagonalProperty = " << diagonalProperty )

    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::sortRowElements(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const bool diagonalFlag )
{

    LAMA_LOG_INFO( logger, "sort elements in each of " << numRows << " rows, diagonal flag = " << diagonalFlag )

    #pragma omp parallel for

    for( IndexType i = 0; i < numRows; ++i )
    {
        // use bubble sort as sort algorithm

        const IndexType start = csrIA[i];
        IndexType end = csrIA[i + 1] - 1;

        LAMA_LOG_TRACE( logger, "row " << i << ": sort " << start << " - " << end )

        bool sorted = false;

        while( !sorted )
        {
            sorted = true; // will be reset if any wrong order appears

            LAMA_LOG_TRACE( logger, "sort from " << start << " - " << end )

            for( IndexType jj = start; jj < end; ++jj )
            {
                bool swapIt = false;

                LAMA_ASSERT_DEBUG(
                    csrJA[jj] != csrJA[jj + 1],
                    "row " << i << ": " << jj << ", column index " << csrJA[jj] << ", " << csrJA[ jj+1] << " appears twice" );

                // if diagonalFlag is set, column i is the smallest one

                if( diagonalFlag && ( csrJA[jj + 1] == i ) && ( csrJA[jj] != i ) )
                {
                    swapIt = true;
                }
                else if( diagonalFlag && ( csrJA[jj] == i ) )
                {
                    swapIt = false;
                }
                else
                {
                    swapIt = csrJA[jj] > csrJA[jj + 1];
                }

                if( swapIt )
                {
                    LAMA_LOG_TRACE( logger, "swap at pos " << jj << " : " << csrJA[jj] << " - " << csrJA[ jj+1 ] )
                    sorted = false;
                    std::swap( csrJA[jj], csrJA[jj + 1] );
                    std::swap( csrValues[jj], csrValues[jj + 1] );
                }
            }

            --end;
        }
    }
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::countNonEmptyRowsByOffsets( const IndexType offsets[], const IndexType numRows )
{
    IndexType counter = 0;

    #pragma omp parallel for reduction( +:counter )

    for( IndexType i = 0; i < numRows; ++i )
    {
        const IndexType nzRow = offsets[i + 1] - offsets[i];

        if( nzRow > 0 )
        {
            counter++;
        }
    }

    LAMA_LOG_INFO( logger, "#non-zero rows = " << counter << ", counted by offsets" )

    return counter;
}

/* --------------------------------------------------------------------------- */

void OpenMPCSRUtils::setNonEmptyRowsByOffsets(
    IndexType rowIndexes[],
    const IndexType numNonEmptyRows,
    const IndexType offsets[],
    const IndexType numRows )
{
    IndexType counter = 0;

    // Note: this routine is not easy to parallelize, no offsets for rowIndexes available

    for( IndexType i = 0; i < numRows; ++i )
    {
        const IndexType nzRow = offsets[i + 1] - offsets[i];

        if( nzRow > 0 )
        {
            rowIndexes[counter] = i;
            counter++;
        }
    }

    LAMA_ASSERT_EQUAL_DEBUG( counter, numNonEmptyRows )

    LAMA_LOG_INFO( logger, "#non-zero rows = " << counter << ", set by offsets" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1,typename ValueType2>
void OpenMPCSRUtils::scaleRows(
    ValueType1 csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const ValueType2 values[] )
{
    #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )

    for( IndexType i = 0; i < numRows; ++i )
    {
        ValueType1 tmp = static_cast<ValueType1>( values[i] );

        for( IndexType j = csrIA[i]; j < csrIA[i + 1]; ++j )
        {
            csrValues[j] *= tmp;
        }
    }
}

/* --------------------------------------------------------------------------- */

static inline IndexType atomicInc( IndexType& var )
{
    return __sync_fetch_and_add( &var, 1 );
}

template<typename ValueType>
void OpenMPCSRUtils::convertCSR2CSC(
    IndexType cscIA[],
    IndexType cscJA[],
    ValueType cscValues[],
    const IndexType rIA[],
    const IndexType rJA[],
    const ValueType rValues[],
    IndexType numRows,
    IndexType numColumns,
    IndexType numValues )
{
    LAMA_REGION( "OpenMP.CSRUtils.CSR2CSC" )

    LAMA_LOG_INFO( logger, "convertCSR2CSC of matrix " << numRows << " x " << numColumns )

    LAMA_ASSERT_EQUAL_DEBUG( numValues, rIA[numRows] )

    // initialization of column counters with 0

    #pragma omp parallel for

    for( IndexType i = 0; i < numColumns; ++i )
    {
        cscIA[i] = 0;
    }

    // loop over all rows of the row matrix to count columns, not yet OpenMP parallelized

    #pragma omp parallel for

    for( IndexType i = 0; i < numRows; ++i )
    {
        // loop over all none zero elements of column i

        for( IndexType jj = rIA[i]; jj < rIA[i + 1]; ++jj )
        {
            IndexType j = rJA[jj];
            #pragma omp atomic
            cscIA[j]++;
        }
    }

    sizes2offsets( cscIA, numColumns );

    LAMA_LOG_INFO( logger, "convertCSR2CSC, #num values counted = " << cscIA[ numColumns ] )

    LAMA_ASSERT_EQUAL_DEBUG( numValues, cscIA[numColumns] )

    // temporary copy neeeded of cscIA

    std::vector<IndexType> cscIA1( numColumns );

    #pragma omp parallel for

    for( IndexType i = 0; i < numColumns; ++i )
    {
        cscIA1[i] = cscIA[i];
    }

    // fill in the array cscJA and cscValues

    #pragma omp parallel for

    for( IndexType i = 0; i < numRows; ++i )
    {
        for( IndexType jj = rIA[i]; jj < rIA[i + 1]; ++jj )
        {
            IndexType j = rJA[jj];
            IndexType k = atomicInc( cscIA1[j] ); // k keeps old value
            cscJA[k] = i;
            cscValues[k] = rValues[jj];
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType UNUSED( numColumns ),
    const IndexType UNUSED( nnz ),
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "normalGEMV<" << getScalarType<ValueType>() << ", #threads = " << omp_get_max_threads() 
                    << ">, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )

    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    // ToDo: for efficiency the following cases should be considered
    // result = y, beta = 1.0 : result += alpha * A * x
    // otherwise: alpha = -1.0, 1.0, beta = 1.0, -1.0, 0.0

    // Note: for beta = 0.0, y could be uninitialized.
    //       0.0 * undefined should deliver 0.0, but valgrind cannot deal with it

    #pragma omp parallel
    {
        // Note: region will be entered by each thread

        LAMA_REGION( "OpenMP.CSR.normalGEMV" )

        #pragma omp for schedule(LAMA_OMP_SCHEDULE)

        for( IndexType i = 0; i < numRows; ++i )
        {
            ValueType temp = 0.0;

            for( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
            {
                IndexType j = csrJA[jj];
                temp += csrValues[jj] * x[j];
            }

            if( 0 == beta )
            {
                result[i] = alpha * temp;
            }
            else
            {
                result[i] = alpha * temp + beta * y[i];
            }
        }
    }

    if( LAMA_LOG_TRACE_ON( logger ) )
    {
        std::cout << "NormalGEMV: result = ";

        for( IndexType i = 0; i < numRows; ++i )
        {
            std::cout << " " << result[i];
        }

        std::cout << std::endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::normalGEVM(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "normalGEVM<" << getScalarType<ValueType>() << ", #threads = " << omp_get_max_threads() << ">, result[" << numColumns << "] = " << alpha << " * x * A + " << beta << " * y " )

    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    // ToDo: for efficiency the cases of alpha and beta = 1.0 / 0.0 should be considered

    // Note: for beta = 0.0, y could be uninitialized.
    //       0.0 * undefined should deliver 0.0, but valgrind cannot deal with it

    #pragma omp parallel
    {
        // Note: region will be entered by each thread
        LAMA_REGION( "OpenMP.CSR.normalGEVM" )
        #pragma omp for schedule(LAMA_OMP_SCHEDULE)

        for( IndexType i = 0; i < numColumns; ++i )
        {
            ValueType sum = 0.0;
            bool diag = false;

            if( i < numRows && csrIA[i] != csrIA[i + 1] && csrJA[csrIA[i]] == i )
            {
                sum += csrValues[csrIA[i]] * x[i];
                diag = true;
            }

            for( IndexType j = 0; j < numRows; ++j )
            {
                for( IndexType k = csrIA[j]; k < csrIA[j + 1]; ++k )
                {
                    if( !( diag && i == j ) && csrJA[k] == i )
                    {
                        sum += csrValues[k] * x[j];
                        break;
                    }
                }
            }

            result[i] = alpha * sum;
        }
    }

    if( 0 != beta )
    {
        #pragma omp for schedule(LAMA_OMP_SCHEDULE)

        for( IndexType i = 0; i < numColumns; ++i )
        {
            result[i] += beta * y[i];
        }
    }

    if( LAMA_LOG_TRACE_ON( logger ) )
    {
        std::cout << "NormalGEVM: result = ";

        for( IndexType i = 0; i < numColumns; ++i )
        {
            std::cout << " " << result[i];
        }

        std::cout << std::endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::sparseGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType numNonZeroRows,
    const IndexType rowIndexes[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    SyncToken* syncToken )
{
    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    #pragma omp parallel
    {
        // Note: region will be entered by each thread

        LAMA_REGION( "OpenMP.CSR.sparseGEMV" )

        #pragma omp for schedule(LAMA_OMP_SCHEDULE)

        for( IndexType ii = 0; ii < numNonZeroRows; ++ii )
        {
            ValueType temp = 0.0;
            IndexType i = rowIndexes[ii];

            for( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
            {
                IndexType j = csrJA[jj];
                temp += csrValues[jj] * x[j];
            }

            result[i] += alpha * temp;
        }
    }

    if( LAMA_LOG_TRACE_ON( logger ) )
    {
        std::cout << "sparseGEMV: result = ";

        for( IndexType ii = 0; ii < numNonZeroRows; ++ii )
        {
            IndexType i = rowIndexes[ii];
            std::cout << " " << i << ":" << result[i];
        }

        std::cout << std::endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::sparseGEVM(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType numColumns,
    const IndexType numNonZeroRows,
    const IndexType rowIndexes[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    SyncToken* syncToken )
{
    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    #pragma omp parallel
    {
        // Note: region will be entered by each thread

        LAMA_REGION( "OpenMP.CSR.normalGEVM" )

        #pragma omp for schedule(LAMA_OMP_SCHEDULE)

        for( IndexType i = 0; i < numColumns; ++i )
        {
            ValueType sum = 0.0;

            for( IndexType jj = 0; jj < numNonZeroRows; ++jj )
            {
                IndexType j = rowIndexes[jj];

                for( IndexType k = csrIA[j]; k < csrIA[j + 1]; ++k )
                {
                    if( csrJA[k] == i )
                    {
                        sum += csrValues[k] * x[i];
                        break;
                    }
                }
            }

            result[i] = alpha * sum;
        }
    }

    if( LAMA_LOG_TRACE_ON( logger ) )
    {
        std::cout << "sparseGEMV: result = ";

        for( IndexType i = 0; i < numColumns; ++i )
        {
            std::cout << " " << result[i];
        }

        std::cout << std::endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::gemm(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType m,
    const IndexType n,
    const IndexType p,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "gemm<" << getScalarType<ValueType>() << ">, " << " result " << m << " x " << n << " CSR " << m << " x " << p )

    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

    for( IndexType i = 0; i < m; ++i )
    {
        for( IndexType k = 0; k < n; ++k )
        {
            ValueType temp = 0.0;

            for( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
            {
                IndexType j = csrJA[jj];

                // LAMA_ASSERT_DEBUG( j < p , "index j = " << j << " out of range " << p )

                // csrValues[jj] stands for CSR( i, j )

                temp += csrValues[jj] * x[j * n + k]; // x(j,k)

            }

            result[i * n + k] = alpha * temp + beta * y[i * n + k];
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::jacobi(
    ValueType* solution,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega,
    const IndexType numRows,
    class SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "jacobi<" << getScalarType<ValueType>() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    const ValueType oneMinusOmega = static_cast<ValueType>( 1.0 ) - omega;

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.CSR.jacobi" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )

        for( IndexType i = 0; i < numRows; i++ )
        {
            ValueType temp = rhs[i];
            const ValueType diag = csrValues[csrIA[i]];

            for( IndexType j = csrIA[i] + 1; j < csrIA[i + 1]; j++ )
            {
                temp -= csrValues[j] * oldSolution[csrJA[j]];
            }

            // here we take advantange of a good branch precondiction

            if( omega == 1.0 )
            {
                solution[i] = temp / diag;
            }
            else if( omega == 0.5 )
            {
                solution[i] = omega * ( temp / diag + oldSolution[i] );
            }
            else
            {
                solution[i] = omega * ( temp / diag ) + oneMinusOmega * oldSolution[i];
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::jacobiHalo(
    ValueType solution[],
    const IndexType localIA[],
    const ValueType localValues[],
    const IndexType haloIA[],
    const IndexType haloJA[],
    const ValueType haloValues[],
    const IndexType haloRowIndexes[],
    const ValueType oldSolution[],
    const ValueType omega,
    const IndexType numNonEmptyRows )
{
    LAMA_LOG_INFO( logger,
                   "jacobiHalo<" << getScalarType<ValueType>() << ">" << ", #rows (not empty) = " << numNonEmptyRows << ", omega = " << omega );

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.CSR.jacabiHalo" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )

        for( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            IndexType i = ii; // default: rowIndexes == NULL stands for identity

            if( haloRowIndexes )
            {
                i = haloRowIndexes[ii];
            }

            ValueType temp = 0.0;

            const ValueType diag = localValues[localIA[i]];

            for( IndexType j = haloIA[i]; j < haloIA[i + 1]; j++ )
            {
                temp += haloValues[j] * oldSolution[haloJA[j]];
            }

            if( omega == 1.0 )
            {
                solution[i] -= temp / diag;
            }
            else
            {
                solution[i] -= omega * ( temp / diag );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::jacobiHaloWithDiag(
    ValueType solution[],
    const ValueType localDiagValues[],
    const IndexType haloIA[],
    const IndexType haloJA[],
    const ValueType haloValues[],
    const IndexType haloRowIndexes[],
    const ValueType oldSolution[],
    const ValueType omega,
    const IndexType numNonEmptyRows )
{
    LAMA_LOG_INFO( logger,
                   "jacobiHaloWithDiag<" << getScalarType<ValueType>() << ">" << ", #rows (not empty) = " << numNonEmptyRows << ", omega = " << omega );

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.CSR.jacabiHaloWithDiag" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )

        for( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            IndexType i = ii; // default: rowIndexes == NULL stands for identity

            if( haloRowIndexes )
            {
                i = haloRowIndexes[ii];
            }

            ValueType temp = 0.0;

            const ValueType diag = localDiagValues[i];

            for( IndexType j = haloIA[i]; j < haloIA[i + 1]; j++ )
            {
                temp += haloValues[j] * oldSolution[haloJA[j]];
            }

            if( omega == 1.0 )
            {
                solution[i] -= temp / diag;
            }
            else
            {
                solution[i] -= omega * ( temp / diag );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::matrixAddSizes(
    IndexType cSizes[],
    const IndexType numRows,
    const IndexType numColumns,
    bool diagonalProperty,
    const IndexType aIA[],
    const IndexType aJA[],
    const IndexType bIA[],
    const IndexType bJA[] )
{
    LAMA_LOG_INFO( logger,
                   "matrixAddSizes for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )

    LAMA_REGION( "OpenMP.CSR.matrixAddSizes" )

    // determine the number of entries in output matrix

    const IndexType NINIT = numColumns + 1; // marks unused colums
    const IndexType END = numColumns + 2; // marks end of list

    #pragma omp parallel
    {
        boost::scoped_array<IndexType> indexList( new IndexType[numColumns] );

        for( IndexType j = 0; j < numColumns; j++ )
        {
            indexList[j] = NINIT;
        }

        #pragma omp for

        for( IndexType i = 0; i < numRows; ++i )
        {
            IndexType length = 0;
            IndexType firstCol = END;

            bool diagonal = false; // will be set if diagonal element is set

            // loop over all none zero elements of row i of input matrix a

            for( IndexType jj = aIA[i]; jj < aIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix a
                // so we are at position a(i,j)

                IndexType j = aJA[jj];

                // element a(i,j) will generate an output element c(i,j)

                if( indexList[j] == NINIT )
                {
                    LAMA_LOG_TRACE( logger, "entry for [" << i << "," << j << "] by a" )

                    // Add column position j to the indexList

                    indexList[j] = firstCol;
                    firstCol = j;
                    ++length;

                    if( j == i )
                    {
                        diagonal = true;
                    }
                }
            }

            for( IndexType jj = bIA[i]; jj < bIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix b
                // so we are at position b(i,j)

                IndexType j = bJA[jj];

                // element a(i,j) will generate an output element c(i,j)

                if( indexList[j] == NINIT )
                {
                    LAMA_LOG_TRACE( logger, "entry for [" << i << "," << j << "] by b" )

                    // Add column position j to the indexList

                    indexList[j] = firstCol;
                    firstCol = j;
                    ++length;

                    if( j == i )
                    {
                        diagonal = true;
                    }
                }
            }

            if( diagonalProperty && !diagonal )
            {
                ++length; // diagaonal needed, but not filled yet
            }

            // so we have now the correct length

            LAMA_LOG_TRACE( logger, "row " << i << " will have " << length << " entries" )

            cSizes[i] = length;

            // reset indexList for next use

            while( firstCol != END )
            {
                LAMA_LOG_TRACE( logger, "entry [" << i << "," << firstCol << "] will be set" )
                IndexType nextCol = indexList[firstCol];
                indexList[firstCol] = NINIT;
                firstCol = nextCol;
            }

        } // end loop over all rows of input matrices
    }

    return sizes2offsets( cSizes, numRows );
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::matrixMultiplySizes(
    IndexType cSizes[],
    const IndexType m,
    const IndexType n,
    const IndexType /* k */,
    bool diagonalProperty,
    const IndexType aIA[],
    const IndexType aJA[],
    const IndexType bIA[],
    const IndexType bJA[] )
{
    LAMA_REGION( "OpenMP.CSR.matrixMultiplySizes" )

    LAMA_LOG_INFO( logger,
                   "matrixMutliplySizes for " << m << " x " << n << " matrix" << ", diagonalProperty = " << diagonalProperty )

    // determine the number of entries in output matrix

    const IndexType NINIT = n + 1; // marks unused colums
    const IndexType END = n + 2; // marks end of list

    IndexType newElems = 0;
    IndexType doubleElems = 0;

    #pragma omp parallel reduction( + : newElems, doubleElems )
    {
        boost::scoped_array<IndexType> indexList( new IndexType[n] );

        for( IndexType j = 0; j < n; j++ )
        {
            indexList[j] = NINIT;
        }

        #pragma omp for

        for( IndexType i = 0; i < m; ++i )
        {
            IndexType length = 0;
            IndexType firstCol = END;

            bool diagonal = false; // will be set if diagonal element is set

            // loop over all none zero elements of row i of input matrix a

            for( IndexType jj = aIA[i]; jj < aIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix a
                // so we are at position a(i,j)

                IndexType j = aJA[jj];

                // loop over all none zero elements of row j of input matrix b
                // that is the row of b that corresponds to the column of the current
                // element of a

                for( IndexType kk = bIA[j]; kk < bIA[j + 1]; ++kk )
                {
                    // k is the column of none zero element kk of row j of input matrix b
                    // so we are looking at position b(j,k)

                    IndexType k = bJA[kk];

                    LAMA_ASSERT_DEBUG( 0 <= k && k < n, "invalid k = " << k )

                    // element a(i,j) an b(j,k) will generate the output element c(i,k)

                    if( indexList[k] == NINIT )
                    {
                        LAMA_LOG_TRACE( logger, "entry for [" << i << "," << k << "]" )

                        // Add column position k to the indexList

                        indexList[k] = firstCol;
                        firstCol = k;
                        ++length;

                        if( k == i )
                        {
                            diagonal = true;
                        }

                        newElems++;
                    }
                    else
                    {
                        doubleElems++;
                    }
                }
            }

            if( diagonalProperty && !diagonal )
            {
                ++length; // diagaonal needed, but not filled yet
            }

            // so we have now the correct length

            LAMA_LOG_TRACE( logger, "row " << i << " will have " << length << " entries" )

            cSizes[i] = length;

            // reset indexList for next use

            while( firstCol != END )
            {
                LAMA_LOG_TRACE( logger, "entry [" << i << "," << firstCol << "] will be set" )
                IndexType nextCol = indexList[firstCol];
                indexList[firstCol] = NINIT;
                firstCol = nextCol;
            }

        } //end loop over all rows of input matrix

    }

    sizes2offsets( cSizes, m );

    return cSizes[m];
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::matrixAdd(
    IndexType cJA[],
    ValueType cValues[],
    const IndexType cIA[],
    const IndexType numRows,
    const IndexType numColumns,
    bool diagonalProperty,
    const ValueType alpha,
    const IndexType aIA[],
    const IndexType aJA[],
    const ValueType aValues[],
    const ValueType beta,
    const IndexType bIA[],
    const IndexType bJA[],
    const ValueType bValues[] )
{
    LAMA_REGION( "OpenMP.CSR.matrixAdd" )

    LAMA_LOG_INFO( logger,
                   "matrixAddJA for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )

    const IndexType NINIT = numColumns + 1;
    const IndexType END = numColumns + 2;

    // determine the number of entries in output matrix

    #pragma omp parallel
    {
        boost::scoped_array<IndexType> indexList( new IndexType[numColumns] );
        boost::scoped_array<ValueType> valueList( new ValueType[numColumns] );

        for( IndexType j = 0; j < numColumns; j++ )
        {
            indexList[j] = NINIT;
            valueList[j] = 0.0;
        }

        #pragma omp for

        for( IndexType i = 0; i < numRows; ++i )
        {
            IndexType length = 0;
            IndexType firstCol = END;

            for( IndexType jj = aIA[i]; jj < aIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix a
                // so we are at position a(i,j)

                IndexType j = aJA[jj];

                valueList[j] += alpha * aValues[jj];

                LAMA_LOG_TRACE( logger, "entry for [" << i << "," << j << "] by a" << ", new val = " << valueList[j] )

                // element a(i,j) will generate an output element c(i,j)

                if( indexList[j] == NINIT )
                {
                    // Add column position j to the indexList

                    indexList[j] = firstCol;

                    firstCol = j;
                    ++length;
                }
            }

            for( IndexType jj = bIA[i]; jj < bIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix b
                // so we are at position b(i,j)

                IndexType j = bJA[jj];

                valueList[j] += beta * bValues[jj];

                LAMA_LOG_TRACE( logger, "entry for [" << i << "," << j << "] by b" << ", new val = " << valueList[j] )

                // element b(i,j) will generate an output element c(i,j)

                if( indexList[j] == NINIT )
                {
                    // Add column position j to the indexList

                    indexList[j] = firstCol;

                    firstCol = j;
                    ++length;
                }
            }

            IndexType offset = cIA[i];

            if( diagonalProperty )
            {
                // first element is reserved for diagonal element
                LAMA_LOG_TRACE( logger, "entry for [" << i << "," << i << "] as diagonal" )
                cJA[offset] = i;
                cValues[offset] = 0.0;
                ++offset;
            }

            // fill in csrJA, csrValues and reset indexList, valueList for next use

            while( firstCol != END )
            {
                IndexType nextCol = indexList[firstCol];
                ValueType val = valueList[firstCol];

                indexList[firstCol] = NINIT;
                valueList[firstCol] = 0.0; // reset for next time

                if( diagonalProperty && firstCol == i )
                {
                    LAMA_LOG_TRACE( logger, "diagonal already added before" )
                    LAMA_LOG_TRACE( logger, "entry for [" << i << "," << i << "] = " << val )
                    cValues[cIA[i]] = val;
                }
                else
                {
                    LAMA_LOG_TRACE( logger, "entry for [" << i << "," << firstCol << "] = " << val )
                    cJA[offset] = firstCol;
                    cValues[offset] = val;
                    ++offset;
                }

                firstCol = nextCol;
            }

            // make sure that we have still the right offsets

            LAMA_ASSERT_EQUAL_DEBUG( offset, cIA[i + 1] )

        } //end loop over all rows of input matrix a
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPCSRUtils::matrixMultiplyJA(
    IndexType cJA[],
    const IndexType cIA[],
    const IndexType numRows,
    const IndexType numColumns,
    bool diagonalProperty,
    const IndexType aIA[],
    const IndexType aJA[],
    const IndexType bIA[],
    const IndexType bJA[] )
{
    LAMA_LOG_INFO( logger,
                   "matrixMutliplyJA for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )

    const IndexType NINIT = numColumns + 1;
    const IndexType END = numColumns + 2;

    // determine the number of entries in output matrix

    #pragma omp parallel
    {
        boost::scoped_array<IndexType> indexList( new IndexType[numColumns] );

        for( IndexType j = 0; j < numColumns; j++ )
        {
            indexList[j] = NINIT;
        }

        #pragma omp for

        for( IndexType i = 0; i < numRows; ++i )
        {
            IndexType length = 0;
            IndexType firstCol = END;

            // loop over all none zero elements of row i of input matrix a

            for( IndexType jj = aIA[i]; jj < aIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix a
                // so we are at position a(i,j)

                IndexType j = aJA[jj];

                // loop over all none zero elements of row j of input matrix b
                // that is the row of b that corresponds to the column of the current
                // element of a

                for( IndexType kk = bIA[j]; kk < bIA[j + 1]; ++kk )
                {
                    // k is the column of none zero element kk of row j of input matrix b
                    // so we are looking at position b(j,k)

                    IndexType k = bJA[kk];

                    // element a(i,j) an b(j,k) will generate the output element c(i,k)

                    if( indexList[k] == NINIT )
                    {
                        // Add column position k to the indexList

                        indexList[k] = firstCol;
                        firstCol = k;
                        ++length;
                        LAMA_LOG_TRACE( logger, "entry for [" << i << "," << k << "]" )
                    }
                }
            }

            IndexType offset = cIA[i];

            if( diagonalProperty )
            {
                // first element is reserved for diagonal element
                LAMA_LOG_TRACE( logger, "entry for [" << i << "," << i << "] as diagonal" )
                cJA[offset] = i;
                ++offset;
            }

            // fill in csrJA and reset indexList for next use

            while( firstCol != END )
            {
                IndexType nextCol = indexList[firstCol];

                indexList[firstCol] = NINIT;

                if( diagonalProperty && firstCol == i )
                {
                    LAMA_LOG_TRACE( logger, "diagonal already added before" )
                }
                else
                {
                    LAMA_LOG_TRACE( logger, "entry for [" << i << "," << firstCol << "]" )
                    cJA[offset] = firstCol;
                    ++offset;
                }

                firstCol = nextCol;
            }

            // make sure that we have still the right offsets

            LAMA_ASSERT_EQUAL_DEBUG( offset, cIA[i + 1] )

        } //end loop over all rows of input matrix a
    }
}

/* --------------------------------------------------------------------------- */
/*
 template<typename ValueType>
 void OpenMPCSRUtils::matrixMultiply (
 ValueType cValues[], const IndexType cIA[], const IndexType cJA[],
 const IndexType numRows,
 const ValueType alpha,
 const IndexType aIA[], const IndexType aJA[], const ValueType aValues[],
 const IndexType bIA[], const IndexType bJA[], const ValueType bValues[] )
 {
 // loop over all rows of output matrix c
 // TODO: Test if this loop is faster with OpenMP than the loop mentioned in
 // Sparse Matrix Multiplication Package (SMMP)
 // Randolph E. Bank Craig C. Douglas April 23, 2001

 #pragma omp parallel for
 for ( IndexType i = 0; i < numRows; ++i)
 {
 //loop over all none zero elements of row i of output matrix c
 for ( IndexType jj = cIA[i]; jj<cIA[i+1]; ++jj )
 {
 IndexType j = cJA[jj];
 cValues[jj] = 0.0;

 if( j == -1 ){
 continue;
 }

 for (IndexType kk = aIA[i]; kk < aIA[i+1]; ++kk)
 {
 ValueType b_kj = 0.0;

 IndexType k = aJA[kk];

 //TODO: See above this loop can be avoided if we choose to work
 //      with a temporarie array like it is done in SMMP.
 //      But we will have less parallelisem in this case

 for ( IndexType ll = bIA[k]; ll < bIA[k+1]; ++ll)
 {
 if ( bJA[ll] ==j )
 {
 b_kj = bValues[ll];
 break;
 }
 }

 cValues[jj] += aValues[kk] * b_kj;
 }

 cValues[jj] *= alpha;

 LAMA_LOG_TRACE( logger, "Computed output Element ( "
 << i << ", " << j << " ) = " << cValues[jj] )
 }
 }
 }
 */

template<typename ValueType>
void OpenMPCSRUtils::matrixMultiply(
    const IndexType cIA[],
    IndexType cJA[],
    ValueType cValues[],
    const IndexType m,
    const IndexType n,
    const IndexType /* k */,
    const ValueType alpha,
    bool diagonalProperty,
    const IndexType aIA[],
    const IndexType aJA[],
    const ValueType aValues[],
    const IndexType bIA[],
    const IndexType bJA[],
    const ValueType bValues[] )
{
    LAMA_REGION( "OpenMP.CSR.matrixMultiply" )

    //TODO: Rewrite this!
    const IndexType NINIT = n + 1;
    const IndexType END = n + 2;

    // determine the number of entries in output matrix
    #pragma omp parallel
    {
        boost::scoped_array<IndexType> indexList( new IndexType[n] );

        for( IndexType j = 0; j < n; j++ )
        {
            indexList[j] = NINIT;
        }

        #pragma omp for

        for( IndexType i = 0; i < m; ++i )
        {
            IndexType length = 0;
            IndexType firstCol = END;

            // loop over all none zero elements of row i of input matrix a

            for( IndexType jj = aIA[i]; jj < aIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix a
                // so we are at position a(i,j)

                IndexType j = aJA[jj];

                // loop over all none zero elements of row j of input matrix b
                // that is the row of b that corresponds to the column of the current
                // element of a

                for( IndexType kk = bIA[j]; kk < bIA[j + 1]; ++kk )
                {
                    // k is the column of none zero element kk of row j of input matrix b
                    // so we are looking at position b(j,k)

                    IndexType k = bJA[kk];

                    // element a(i,j) an b(j,k) will generate the output element c(i,k)

                    if( indexList[k] == NINIT )
                    {
                        // Add column position k to the indexList

                        indexList[k] = firstCol;
                        firstCol = k;
                        ++length;
                        LAMA_LOG_TRACE( logger, "entry for [" << i << "," << k << "]" )
                    }
                }
            }

            IndexType offset = cIA[i];

            if( diagonalProperty )
            {
                // first element is reserved for diagonal element
                LAMA_LOG_TRACE( logger, "entry for [" << i << "," << i << "] as diagonal" )
                cJA[offset] = i;
                ++offset;
            }

            // fill in csrJA and reset indexList for next use

            while( firstCol != END )
            {
                IndexType nextCol = indexList[firstCol];

                indexList[firstCol] = NINIT;

                if( diagonalProperty && firstCol == i )
                {
                    LAMA_LOG_TRACE( logger, "diagonal already added before" )
                }
                else
                {
                    LAMA_LOG_TRACE( logger, "entry for [" << i << "," << firstCol << "]" )
                    cJA[offset] = firstCol;
                    ++offset;
                }

                firstCol = nextCol;
            }

            // make sure that we have still the right offsets

            // LAMA_ASSERT_EQUAL_DEBUG( offset, cIA[i+1] )

        } //end loop over all rows of input matrix a
    }

    #pragma omp parallel for

    for( IndexType i = 0; i < m; ++i )
    {
        //loop over all none zero elements of row i of output matrix c
        for( IndexType jj = cIA[i]; jj < cIA[i + 1]; ++jj )
        {
            IndexType j = cJA[jj];
            cValues[jj] = 0.0;

            if( j == -1 )
            {
                continue;
            }

            for( IndexType kk = aIA[i]; kk < aIA[i + 1]; ++kk )
            {
                ValueType b_kj = 0.0;

                IndexType k = aJA[kk];

                //TODO: See above this loop can be avoided if we choose to work
                //      with a temporarie array like it is done in SMMP.
                //      But we will have less parallelisem in this case

                for( IndexType ll = bIA[k]; ll < bIA[k + 1]; ++ll )
                {
                    if( bJA[ll] == j )
                    {
                        b_kj = bValues[ll];
                        break;
                    }
                }

                cValues[jj] += aValues[kk] * b_kj;
            }

            cValues[jj] *= alpha;

            LAMA_LOG_TRACE( logger, "Computed output Element ( " << i << ", " << j << " ) = " << cValues[jj] )
        }
    }
}

/* --------------------------------------------------------------------------- */

static IndexType findCol( const IndexType columns[], const IndexType n, const IndexType j, IndexType lastPos )
{
    // lastPos is index from where to start

    for( IndexType i = 0; i < n; ++i )
    {
        if( columns[lastPos] == j )
        {
            // found, return lastPos, increment it for next search

            return lastPos++;
        }

        lastPos++;

        if( lastPos > n )
        {
            lastPos = 0;
        }
    }

    // not found, return size

    return n;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPCSRUtils::absMaxDiffRowUnsorted(
    const IndexType n1,
    const IndexType csrJA1[],
    const ValueType csrValues1[],
    const IndexType n2,
    const IndexType csrJA2[],
    const ValueType csrValues2[] )
{
    // No assumption about any sorting in a row

    ValueType val = static_cast<ValueType>( 0.0 );

    IndexType helpIndex = 0; // some kind of thread-safe global value for findCol

    for( IndexType i1 = 0; i1 < n1; ++i1 )
    {
        ValueType diff = csrValues1[i1];

        IndexType j = csrJA1[i1];
        IndexType i2 = findCol( csrJA2, n2, j, helpIndex );

        if( i2 < n2 )
        {
            diff -= csrValues2[i2];
        }

        diff = abs( diff );

        if( diff > val )
        {
            val = diff;
        }
    }

    // check for row elements in 2nd matrix that are not available in first one

    helpIndex = 0;

    for( IndexType i2 = 0; i2 < n2; ++i2 )
    {
        IndexType j = csrJA2[i2];
        IndexType i1 = findCol( csrJA1, n1, j, helpIndex );

        if( i1 < n1 )
        {
            continue; // already compare in first loop
        }

        ValueType diff = abs( csrValues2[i2] );

        if( diff > val )
        {
            val = diff;
        }
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPCSRUtils::absMaxDiffRowSorted(
    const IndexType n1,
    const IndexType csrJA1[],
    const ValueType csrValues1[],
    const IndexType n2,
    const IndexType csrJA2[],
    const ValueType csrValues2[] )
{
    // Note: the implementation assumes that rows are sorted according to column indexes

    ValueType val = static_cast<ValueType>( 0.0 );

    IndexType i2 = 0;
    IndexType i1 = 0;

    while( i1 < n1 || i2 < n2 )
    {
        ValueType diff = 0;

        if( i1 >= n1 ) // row of matrix1 completely traversed
        {
            diff = csrValues2[i2];
            ++i2;
        }
        else if( i2 >= n2 ) // row of matrix2 completely traversed
        {
            diff = csrValues1[i1];
            ++i1;
        }
        else
        {
            // so we have still values in both rows

            IndexType j1 = csrJA1[i1];
            IndexType j2 = csrJA2[i2];

            if( j1 == j2 )
            {
                diff = csrValues1[i1] - csrValues2[i2];
                ++i1;
                ++i2;
            }
            else if( j1 < j2 )
            {
                diff = csrValues1[i1];
                ++i1;

                if( i1 < n1 )
                {
                    LAMA_ASSERT_ERROR( csrJA1[i1 - 1] < csrJA1[i1], "unsorted col indexes at csrJA1[" << i1 << "]" )
                }
            }
            else if( j1 > j2 )
            {
                diff = csrValues2[i2];
                ++i2;

                if( i2 < n2 )
                {
                    LAMA_ASSERT_ERROR( csrJA2[i2 - 1] < csrJA2[i2], "unsorted col indexes at csrJA2[" << i2 << "]" )
                }
            }
        }

        diff = abs( diff );

        if( diff > val )
        {
            val = diff;
        }
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPCSRUtils::absMaxDiffVal(
    IndexType numRows,
    bool sortedRows,
    const IndexType csrIA1[],
    const IndexType csrJA1[],
    const ValueType csrValues1[],
    const IndexType csrIA2[],
    const IndexType csrJA2[],
    const ValueType csrValues2[] )
{
    LAMA_LOG_INFO( logger,
                   "absMaxDiffVal<" << getScalarType<ValueType>() << ">: " << "csr[" << numRows << "], sorted = " << sortedRows )

    ValueType (*absMaxDiffRow)(
        const IndexType,
        const IndexType[],
        const ValueType[],
        const IndexType,
        const IndexType[],
        const ValueType[] );

    if( sortedRows )
    {
        absMaxDiffRow = OpenMPCSRUtils::absMaxDiffRowSorted<ValueType>;
    }
    else
    {
        absMaxDiffRow = OpenMPCSRUtils::absMaxDiffRowUnsorted<ValueType>;
    }

    ValueType val = static_cast<ValueType>( 0.0 );

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>( 0.0 );

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )

        for( IndexType i = 0; i < numRows; ++i )
        {
            IndexType offs1 = csrIA1[i];
            IndexType offs2 = csrIA2[i];
            IndexType n1 = csrIA1[i + 1] - offs1;
            IndexType n2 = csrIA2[i + 1] - offs2;

            ValueType maxRow = 0.0;

            maxRow = absMaxDiffRow( n1, &csrJA1[offs1], &csrValues1[offs1], n2, &csrJA2[offs2], &csrValues2[offs2] );

            if( maxRow > threadVal )
            {
                threadVal = maxRow;
            }
        }

        #pragma omp critical
        {
            LAMA_LOG_TRACE( logger, "max val of thread  = " << threadVal << ", global was " << val )

            if( threadVal > val )
            {
                val = threadVal;
            }
        }
    }

    return val;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPCSRUtils::setInterface( CSRUtilsInterface& CSRUtils )
{
    LAMA_INTERFACE_REGISTER( CSRUtils, sizes2offsets )
    LAMA_INTERFACE_REGISTER( CSRUtils, offsets2sizes )
    LAMA_INTERFACE_REGISTER( CSRUtils, validOffsets )
    LAMA_INTERFACE_REGISTER( CSRUtils, hasDiagonalProperty )

    LAMA_INTERFACE_REGISTER( CSRUtils, matrixAddSizes )
    LAMA_INTERFACE_REGISTER( CSRUtils, matrixMultiplySizes )
    LAMA_INTERFACE_REGISTER( CSRUtils, matrixMultiplyJA )

#define LAMA_CSR_UTILS2_REGISTER(z, J, TYPE )                                        \
    LAMA_INTERFACE_REGISTER_TT( CSRUtils, scaleRows, TYPE, ARITHMETIC_TYPE##J )      \

#define LAMA_CSR_UTILS_REGISTER(z, I, _)                                             \
    LAMA_INTERFACE_REGISTER_T( CSRUtils, convertCSR2CSC, ARITHMETIC_TYPE##I )        \
    LAMA_INTERFACE_REGISTER_T( CSRUtils, sortRowElements, ARITHMETIC_TYPE##I )       \
    LAMA_INTERFACE_REGISTER_T( CSRUtils, normalGEMV, ARITHMETIC_TYPE##I )            \
    LAMA_INTERFACE_REGISTER_T( CSRUtils, sparseGEMV, ARITHMETIC_TYPE##I )            \
    LAMA_INTERFACE_REGISTER_T( CSRUtils, normalGEVM, ARITHMETIC_TYPE##I )            \
    LAMA_INTERFACE_REGISTER_T( CSRUtils, sparseGEVM, ARITHMETIC_TYPE##I )            \
    LAMA_INTERFACE_REGISTER_T( CSRUtils, gemm, ARITHMETIC_TYPE##I )                  \
    LAMA_INTERFACE_REGISTER_T( CSRUtils, matrixAdd, ARITHMETIC_TYPE##I )             \
    LAMA_INTERFACE_REGISTER_T( CSRUtils, matrixMultiply, ARITHMETIC_TYPE##I )        \
    LAMA_INTERFACE_REGISTER_T( CSRUtils, jacobi, ARITHMETIC_TYPE##I )                \
    LAMA_INTERFACE_REGISTER_T( CSRUtils, jacobiHalo, ARITHMETIC_TYPE##I )            \
    LAMA_INTERFACE_REGISTER_T( CSRUtils, jacobiHaloWithDiag, ARITHMETIC_TYPE##I )    \
    LAMA_INTERFACE_REGISTER_T( CSRUtils, absMaxDiffVal, ARITHMETIC_TYPE##I )         \
    \
    BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT,                                            \
                     LAMA_CSR_UTILS2_REGISTER,                                       \
                     ARITHMETIC_TYPE##I )                                            \

    BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_CSR_UTILS_REGISTER, _ )

#undef LAMA_CSR_UTILS_REGISTER
#undef LAMA_CSR_UTILS2_REGISTER

}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool OpenMPCSRUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( memory::context::Host );
    setInterface( interface.CSRUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPCSRUtils::initialized = registerInterface();

} // namespace lama
