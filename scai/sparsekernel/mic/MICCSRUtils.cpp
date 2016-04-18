/**
 * @file MICCSRUtils.cpp
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
 * @brief Implementation of CSR utilities with MIC
 * @author Thomas Brandes
 * @date 02.07.2013
 * @since 1.1.0
 */

// for dll_import
#include <scai/sparsekernel/mic/MICCSRUtils.hpp>

// local libray
#include <scai/sparsekernel/CSRKernelTrait.hpp>

// internal scai libraries
#include <scai/tasking/mic/MICSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/OpenMP.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/bind.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>

namespace scai
{

using tasking::SyncToken;
using tasking::MICSyncToken;

using hmemo::MICContext;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( MICCSRUtils::logger, "MIC.CSRUtils" )

/** Number of minimal threads for which parallelization is effective. */

static int minThreads = 3;

/* --------------------------------------------------------------------------- */

IndexType MICCSRUtils::scanSerial( IndexType array[], const IndexType numValues )
{
    SCAI_LOG_DEBUG( logger, "scanSerial: " << numValues << " entries" )

    // In this case we do it just serial, probably faster

    IndexType runningSum = 0;

    void* arrayPtr = array;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( arrayPtr, numValues ), out( runningSum )
    {
        IndexType* array = static_cast<IndexType*>( arrayPtr );

        runningSum = 0;

        for( IndexType i = 0; i < numValues; i++ )
        {
            IndexType tmp = runningSum;
            runningSum += array[i];
            array[i] = tmp;
        }
    }

    return runningSum;;
}

/* --------------------------------------------------------------------------- */

IndexType MICCSRUtils::scanParallel( PartitionId numThreads, IndexType array[], const IndexType numValues )
{
    // std::cout << "Scan with " << numThreads << " in parallel" << std::endl;

    SCAI_LOG_DEBUG( logger, "scanParallel: " << numValues << " entries for " << numThreads << " threads" )

    void* arrayPtr = array;

    void* threadCounterPtr; // temporary array on MIC device

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( arrayPtr, numValues, numThreads ), out( threadCounterPtr )
    {
        IndexType* threadCounter = new IndexType[numThreads];

        const IndexType* array = (IndexType*) arrayPtr;

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

        threadCounterPtr = threadCounter;
    }

    IndexType* threadCounter = static_cast<IndexType*>( threadCounterPtr );

    IndexType runningSum = scanSerial( threadCounter, numThreads );

    // Each thread sets now its offsets

#pragma offload target( mic : device ), in( threadCounterPtr, arrayPtr, numValues )
    {
        IndexType* threadCounter = static_cast<IndexType*>( threadCounterPtr );
        IndexType* array = static_cast<IndexType*>( arrayPtr );

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

        delete[] threadCounter;
    }

    return runningSum;;
}

/* --------------------------------------------------------------------------- */

IndexType MICCSRUtils::scan( IndexType array[], const IndexType numValues )
{
    int numThreads = 1; // will be set to available threads in parallel region

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), out( numThreads )
    {
        #pragma omp parallel
        #pragma omp master
        {
            numThreads = omp_get_num_threads();
        }
    }

    SCAI_LOG_INFO( logger, "scan " << numValues << " entries, #threads = " << numThreads )

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

bool MICCSRUtils::validOffsets( const IndexType array[], const IndexType n, const IndexType total )
{
    SCAI_LOG_INFO( logger, "check offset array[ " << n << "] for validity, total = " << total )

    bool validFlag = true;

    int device = MICContext::getCurrentDevice();

    const void* arrayPtr = array;

#pragma offload target( mic : device ) in( arrayPtr, n, total ) out( validFlag )
    {
        validFlag = true;

        const IndexType* array = static_cast<const IndexType*>( arrayPtr );

        #pragma omp parallel for reduction( && : validFlag )

        for( IndexType i = 0; i < n; i++ )
        {
            if( array[i] > array[i + 1] )
            {
                validFlag = false;
            }
        }

        if( array[n] != total )
        {
            validFlag = false;
        }
    }

    return validFlag;
}

/* --------------------------------------------------------------------------- */

IndexType MICCSRUtils::sizes2offsets( IndexType array[], const IndexType numValues )
{
    IndexType totalValues = scan( array, numValues );

    // now we have to set array[numValues] = totalValues

    void* arrayPtr = array;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( arrayPtr, numValues, totalValues )
    {
        IndexType* array = static_cast<IndexType*>( arrayPtr );

        array[numValues] = totalValues;
    }

    SCAI_LOG_INFO( logger, "sizes2offsets, #values = " << numValues << ", total = " << totalValues )

    return totalValues;
}

/* --------------------------------------------------------------------------- */

void MICCSRUtils::offsets2sizes( IndexType sizes[], const IndexType offsets[], const IndexType numRows )
{
    void* sizesPtr = sizes;
    const void* offsetsPtr = offsets;

    int device = MICContext::getCurrentDevice();

#pragma offload target( MIC : device ), in( sizesPtr, offsetsPtr, numRows )
    {
        IndexType* sizes = static_cast<IndexType*>( sizesPtr );
        const IndexType* offsets = static_cast<const IndexType*>( offsetsPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; i++ )
        {
            sizes[i] = offsets[i + 1] - offsets[i];
        }
    }
}

/* --------------------------------------------------------------------------- */

void MICCSRUtils::offsets2sizesGather(
    IndexType sizes[],
    const IndexType offsets[],
    const IndexType rowIndexes[],
    const IndexType numRows )
{
    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

    for( IndexType i = 0; i < numRows; i++ )
    {
        IndexType row = rowIndexes[i];
        sizes[i] = offsets[row + 1] - offsets[row];
    }
}

/* --------------------------------------------------------------------------- */

bool MICCSRUtils::hasDiagonalProperty( const IndexType numDiagonals, const IndexType csrIA[], const IndexType csrJA[] )
{
    SCAI_LOG_INFO( logger, "hasDiagonalProperty, #numDiagonals = " << numDiagonals )

    bool diagonalProperty = true;

    size_t csrIAPtr = (size_t) csrIA;
    size_t csrJAPtr = (size_t) csrJA;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( csrIAPtr, csrJAPtr, numDiagonals ), out( diagonalProperty )
    {
        const IndexType* csrIA = (IndexType*) csrIAPtr;
        const IndexType* csrJA = (IndexType*) csrJAPtr;

        diagonalProperty = true;

        #pragma omp parallel for reduction( && : diagonalProperty )

        for( IndexType i = 0; i < numDiagonals; ++i )
        {
            if( diagonalProperty )
            {
                if( csrIA[i] == csrIA[i + 1] )
                {
                    diagonalProperty = false;
                }
                else if( csrJA[csrIA[i]] != i )
                {
                    diagonalProperty = false;
                }
            }
        }
    }

    SCAI_LOG_DEBUG( logger, "hasDiagonalProperty = " << diagonalProperty )

    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICCSRUtils::sortRowElements(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const bool diagonalFlag )
{

    SCAI_LOG_INFO( logger, "sort elements in each of " << numRows << " rows, diagonal flag = " << diagonalFlag )

    #pragma omp parallel for

    for( IndexType i = 0; i < numRows; ++i )
    {
        // use bubble sort as sort algorithm

        const IndexType start = csrIA[i];
        IndexType end = csrIA[i + 1] - 1;

        bool sorted = false;

        while( !sorted )
        {
            sorted = true; // will be reset if any wrong order appears

            for( IndexType jj = start; jj < end; ++jj )
            {
                bool swapIt = false;

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

IndexType MICCSRUtils::countNonEmptyRowsByOffsets( const IndexType offsets[], const IndexType numRows )
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

    SCAI_LOG_INFO( logger, "#non-zero rows = " << counter << ", counted by offsets" )

    return counter;
}

/* --------------------------------------------------------------------------- */

void MICCSRUtils::setNonEmptyRowsByOffsets(
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

    SCAI_ASSERT_EQUAL_DEBUG( counter, numNonEmptyRows )

    SCAI_LOG_INFO( logger, "#non-zero rows = " << counter << ", set by offsets" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1,typename ValueType2>
void MICCSRUtils::scaleRows(
    ValueType1 csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const ValueType2 values[] )
{
    void* csrValuesPtr = csrValues;
    const void* csrIAPtr = csrIA;
    const void* valuesPtr = values;

    int device = MICContext::getCurrentDevice();

#pragma offload target( MIC : device ) in( numRows, csrValuesPtr, csrIAPtr, valuesPtr )
    {
        ValueType1* csrValues = static_cast<ValueType1*>( csrValuesPtr );
        const ValueType2* values = static_cast<const ValueType2*>( valuesPtr );
        const IndexType* csrIA = static_cast<const IndexType*>( csrIAPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; ++i )
        {
            ValueType1 tmp = static_cast<ValueType1>( values[i] );

            for( IndexType j = csrIA[i]; j < csrIA[i + 1]; ++j )
            {
                csrValues[j] *= tmp;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICCSRUtils::convertCSR2CSC(
    IndexType cIA[],
    IndexType cJA[],
    ValueType cValues[],
    const IndexType rIA[],
    const IndexType rJA[],
    const ValueType rValues[],
    IndexType numRows,
    IndexType numColumns,
    IndexType numValues )
{
    SCAI_LOG_INFO( logger, "convertCSR2CSC of matrix " << numRows << " x " << numColumns )

    SCAI_ASSERT_EQUAL_DEBUG( numValues, rIA[numRows] )

    // initialization of column counters with 0

    for( IndexType i = 0; i < numColumns; ++i )
    {
        cIA[i] = 0;
    }

    // loop over all rows of the row matrix to count columns, not yet MIC parallelized

    for( IndexType i = 0; i < numRows; ++i )
    {
        // loop over all none zero elements of column i

        for( IndexType jj = rIA[i]; jj < rIA[i + 1]; ++jj )
        {
            IndexType j = rJA[jj];
            SCAI_ASSERT_DEBUG( j < numColumns, "column index " << j << " out of range, #cols = " << numColumns )
            cIA[j]++;
        }
    }

    sizes2offsets( cIA, numColumns );

    SCAI_LOG_INFO( logger, "convertCSR2CSC, #num values counted = " << cIA[ numColumns ] )

    SCAI_ASSERT_EQUAL_DEBUG( numValues, cIA[numColumns] )

    // fill in the array cJA and cValues

    for( IndexType i = 0; i < numRows; ++i )
    {
        for( IndexType jj = rIA[i]; jj < rIA[i + 1]; ++jj )
        {
            IndexType j = rJA[jj];
            cJA[cIA[j]] = i;
            cValues[cIA[j]] = rValues[jj];
            cIA[j]++;
        }
    }

    SCAI_LOG_INFO( logger, "convertCSR2CSC, #num values counted = " << cIA[ numColumns ] )

    // set back the old offsets

    for( IndexType i = numColumns; i > 0; --i )
    {
        cIA[i] = cIA[i - 1];
    }

    cIA[0] = 0;

    SCAI_ASSERT_EQUAL_DEBUG( cIA[numColumns], numValues )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICCSRUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType /* numColumns */,
    const IndexType /* nnz */,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    SCAI_REGION( "MIC.CSRUtils.normalGEMV" )

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << common::getScalarType<ValueType>() << ">, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
    }

    // SCAI_REGION( "MIC.CSR.normalGEMV" )

    void* resultPtr = result;
    const void* xPtr = x;
    const void* yPtr = y;
    const void* csrIAPtr = csrIA;
    const void* csrJAPtr = csrJA;
    const void* csrValuesPtr = csrValues;
    const ValueType* alphaPtr = &alpha;
    const ValueType* betaPtr = &beta;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device  )in( resultPtr, xPtr, yPtr, \
                          csrIAPtr, csrJAPtr, csrValuesPtr, alphaPtr[0:1], betaPtr[0:1], numRows )
    {
        ValueType* result = (ValueType*) resultPtr;
        const ValueType* x = (ValueType*) xPtr;
        const ValueType* y = (ValueType*) yPtr;
        const IndexType* csrIA = (IndexType*) csrIAPtr;
        const IndexType* csrJA = (IndexType*) csrJAPtr;
        const ValueType* csrValues = (ValueType*) csrValuesPtr;

        #pragma omp parallel
        {
            #pragma omp for schedule(SCAI_OMP_SCHEDULE)

            for( IndexType i = 0; i < numRows; ++i )
            {
                ValueType temp = static_cast<ValueType>(0.0);

                for( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
                {
                    IndexType j = csrJA[jj];
                    temp += csrValues[jj] * x[j];
                }

                if( *betaPtr == static_cast<ValueType>( 0.0 ) )
                {
                    result[i] = (*alphaPtr) * temp;
                }
                else
                {
                    result[i] = (*alphaPtr) * temp + (*betaPtr) * y[i];
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICCSRUtils::sparseGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType numNonZeroRows,
    const IndexType rowIndexes[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
    }

    // SCAI_REGION( "MIC.CSR.sparseGEMV" )

    // conversion of pointer to size_t to cheat offload

    const size_t resultPtr = (size_t) result;
    const size_t xPtr = (size_t) x;
    const size_t rowIndexesPtr = (size_t) rowIndexes;
    const size_t csrIAPtr = (size_t) csrIA;
    const size_t csrJAPtr = (size_t) csrJA;
    const size_t csrValuesPtr = (size_t) csrValues;
    const ValueType* alphaPtr = &alpha;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( resultPtr, xPtr, rowIndexesPtr, csrIAPtr, csrJAPtr, csrValuesPtr, alphaPtr[0:1], numNonZeroRows )
    {
        ValueType* result = (ValueType*) resultPtr;
        const ValueType* x = (ValueType*) xPtr;
        const IndexType* rowIndexes = (IndexType*) rowIndexesPtr;
        const IndexType* csrIA = (IndexType*) csrIAPtr;
        const IndexType* csrJA = (IndexType*) csrJAPtr;
        const ValueType* csrValues = (ValueType*) csrValuesPtr;

        #pragma omp parallel
        {
            // Note: region will be entered by each thread

            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( IndexType ii = 0; ii < numNonZeroRows; ++ii )
            {
                ValueType temp = static_cast<ValueType>(0.0);
                IndexType i = rowIndexes[ii];

                for( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
                {
                    IndexType j = csrJA[jj];
                    temp += csrValues[jj] * x[j];
                }

                result[i] += (*alphaPtr) * temp;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICCSRUtils::gemm(
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
    const ValueType csrValues[] )
{
    SCAI_LOG_INFO( logger,
                   "gemm<" << common::getScalarType<ValueType>() << ">, " << " result " << m << " x " << n << " CSR " << m << " x " << p )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
    }

    void* resultPtr = result;
    const void* xPtr = x;
    const void* yPtr = y;

    const void* csrIAPtr = csrIA;
    const void* csrJAPtr = csrJA;
    const void* csrValuesPtr = csrValues;
    const ValueType* alphaPtr = &alpha;
    const ValueType* betaPtr = &beta;

    int device = MICContext::getCurrentDevice();

    // gemm is  dense = sparse * dense

#pragma offload target( MIC : device ) in( alphaPtr[0:1], betaPtr[0:1], m, n, p, xPtr, yPtr, \
                                               csrIAPtr, csrJAPtr, csrValuesPtr, resultPtr )
    {
        ValueType* result = static_cast<ValueType*>( resultPtr );

        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const ValueType* y = static_cast<const ValueType*>( yPtr );
        const IndexType* csrIA = static_cast<const IndexType*>( csrIAPtr );
        const IndexType* csrJA = static_cast<const IndexType*>( csrJAPtr );
        const ValueType* csrValues = static_cast<const ValueType*>( csrValuesPtr );

        #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

        for( IndexType i = 0; i < m; ++i )
        {
            for( IndexType k = 0; k < n; ++k )
            {
                ValueType temp = static_cast<ValueType>(0.0);

                for( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
                {
                    IndexType j = csrJA[jj];

                    // SCAI_ASSERT_DEBUG( j < p , "index j = " << j << " out of range " << p )

                    // csrValues[jj] stands for CSR( i, j )

                    temp += csrValues[jj] * x[j * n + k]; // x(j,k)

                }

                result[i * n + k] = (*alphaPtr) * temp + (*betaPtr) * y[i * n + k];
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICCSRUtils::jacobi(
    ValueType* solution,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega,
    const IndexType numRows )
{
    SCAI_LOG_INFO( logger,
                   "jacobi<" << common::getScalarType<ValueType>() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
    }

    // SCAI_REGION( "MIC.CSR.jacobi" )

    const size_t solutionPtr = (size_t) solution;
    const size_t oldSolutionPtr = (size_t) oldSolution;
    const size_t rhsPtr = (size_t) rhs;
    const size_t csrIAPtr = (size_t) csrIA;
    const size_t csrJAPtr = (size_t) csrJA;
    const size_t csrValuesPtr = (size_t) csrValues;

    const ValueType* omegaPtr = &omega;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( solutionPtr, oldSolutionPtr, rhsPtr, csrIAPtr, csrJAPtr, csrValuesPtr, omegaPtr[0:1], numRows )
    {
        ValueType* solution = (ValueType*) solutionPtr;
        const ValueType* oldSolution = (ValueType*) oldSolutionPtr;
        const ValueType* rhs = (ValueType*) rhsPtr;
        const IndexType* csrIA = (IndexType*) csrIAPtr;
        const IndexType* csrJA = (IndexType*) csrJAPtr;
        const ValueType* csrValues = (ValueType*) csrValuesPtr;

        const ValueType oneMinusOmega = static_cast<ValueType>(1.0) - (*omegaPtr);

        #pragma omp parallel
        {
            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( IndexType i = 0; i < numRows; i++ )
            {
                ValueType temp = rhs[i];
                const ValueType diag = csrValues[csrIA[i]];

                for( IndexType j = csrIA[i] + 1; j < csrIA[i + 1]; j++ )
                {
                    temp -= csrValues[j] * oldSolution[csrJA[j]];
                }

                // here we take advantange of a good branch precondiction

                if( (*omegaPtr) == static_cast<ValueType>( 1.0 ) )
                {
                    solution[i] = temp / diag;
                }
                else if( (*omegaPtr) == 0.5 )
                {
                    solution[i] = (*omegaPtr) * ( temp / diag + oldSolution[i] );
                }
                else
                {
                    solution[i] = (*omegaPtr) * ( temp / diag ) + oneMinusOmega * oldSolution[i];
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICCSRUtils::jacobiHalo(
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
    SCAI_LOG_INFO( logger,
                   "jacobiHalo<" << common::getScalarType<ValueType>() << ">" << ", #rows (not empty) = " << numNonEmptyRows << ", omega = " << omega );

    // SCAI_REGION( "MIC.CSR.jacabiHalo" )

    const size_t solutionPtr = (size_t) solution;
    const size_t oldSolutionPtr = (size_t) oldSolution;
    const size_t localValuesPtr = (size_t) localValues;
    const size_t localIAPtr = (size_t) localIA;
    const size_t haloIAPtr = (size_t) haloIA;
    const size_t haloJAPtr = (size_t) haloJA;
    const size_t haloValuesPtr = (size_t) haloValues;
    const size_t haloRowIndexesPtr = (size_t) haloRowIndexes;

    const ValueType* omegaPtr = &omega;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( solutionPtr, oldSolutionPtr, localValuesPtr, localIAPtr, \
                                            haloIAPtr, haloJAPtr, haloValuesPtr, haloRowIndexesPtr, \
                                            omegaPtr[0:1], numNonEmptyRows )
    {
        ValueType* solution = (ValueType*) solutionPtr;
        const ValueType* oldSolution = (ValueType*) oldSolutionPtr;
        const ValueType* localValues = (ValueType*) localValuesPtr;
        const IndexType* localIA = (IndexType*) localIAPtr;
        const IndexType* haloIA = (IndexType*) haloIAPtr;
        const IndexType* haloJA = (IndexType*) haloJAPtr;
        const ValueType* haloValues = (ValueType*) haloValuesPtr;
        const IndexType* haloRowIndexes = (IndexType*) haloRowIndexesPtr;

        #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            IndexType i = ii; // default: rowIndexes == NULL stands for identity

            if( haloRowIndexes )
            {
                i = haloRowIndexes[ii];
            }

            ValueType temp = static_cast<ValueType>(0.0);

            const ValueType diag = localValues[localIA[i]];

            for( IndexType j = haloIA[i]; j < haloIA[i + 1]; j++ )
            {
                temp += haloValues[j] * oldSolution[haloJA[j]];
            }

            if( ( *omegaPtr ) == static_cast<ValueType>( 1.0 ) )
            {
                solution[i] -= temp / diag;
            }
            else
            {
                solution[i] -= (*omegaPtr) * ( temp / diag );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICCSRUtils::jacobiHaloWithDiag(
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
    SCAI_LOG_INFO( logger,
                   "jacobiHaloWithDiag<" << common::getScalarType<ValueType>() << ">" << ", #rows (not empty) = " << numNonEmptyRows << ", omega = " << omega );

    // SCAI_REGION( "MIC.CSR.jacabiHaloWithDiag" )

    const size_t solutionPtr = (size_t) solution;
    const size_t oldSolutionPtr = (size_t) oldSolution;
    const size_t localDiagValuesPtr = (size_t) localDiagValues;
    const size_t haloIAPtr = (size_t) haloIA;
    const size_t haloJAPtr = (size_t) haloJA;
    const size_t haloValuesPtr = (size_t) haloValues;
    const size_t haloRowIndexesPtr = (size_t) haloRowIndexes;

    const ValueType* omegaPtr = &omega;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( solutionPtr, oldSolutionPtr, localDiagValuesPtr, \
                                            haloIAPtr, haloJAPtr, haloValuesPtr, haloRowIndexesPtr, \
                                            omegaPtr[0:1], numNonEmptyRows )
    {
        ValueType* solution = (ValueType*) solutionPtr;
        const ValueType* oldSolution = (ValueType*) oldSolutionPtr;
        const ValueType* localDiagValues = (ValueType*) localDiagValuesPtr;
        const IndexType* haloIA = (IndexType*) haloIAPtr;
        const IndexType* haloJA = (IndexType*) haloJAPtr;
        const ValueType* haloValues = (ValueType*) haloValuesPtr;
        const IndexType* haloRowIndexes = (IndexType*) haloRowIndexesPtr;

        const ValueType oneMinusOmega = static_cast<ValueType>(1.0) - (*omegaPtr);

        #pragma omp parallel
        {
            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
            {
                IndexType i = ii; // default: rowIndexes == NULL stands for identity

                if( haloRowIndexes )
                {
                    i = haloRowIndexes[ii];
                }

                ValueType temp = static_cast<ValueType>(0.0);

                const ValueType diag = localDiagValues[i];

                for( IndexType j = haloIA[i]; j < haloIA[i + 1]; j++ )
                {
                    temp += haloValues[j] * oldSolution[haloJA[j]];
                }

                if( (*omegaPtr) == static_cast<ValueType>( 1.0 ) )
                {
                    solution[i] -= temp / diag;
                }
                else
                {
                    solution[i] -= (*omegaPtr) * ( temp / diag );
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

IndexType MICCSRUtils::matrixAddSizes(
    IndexType cSizes[],
    const IndexType numRows,
    const IndexType numColumns,
    bool diagonalProperty,
    const IndexType aIA[],
    const IndexType aJA[],
    const IndexType bIA[],
    const IndexType bJA[] )
{
    // SCAI_REGION( "MIC.CSR.matrixAddSizes" )

    SCAI_LOG_INFO( logger,
                   "matrixAddSizes for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )

    // determine the number of entries in output matrix

    const IndexType NINIT = numColumns + 1; // marks unused colums
    const IndexType END = numColumns + 2; // marks end of list

    void* cSizesPtr = cSizes;

    const void* aIAPtr = aIA;
    const void* aJAPtr = aJA;
    const void* bIAPtr = bIA;
    const void* bJAPtr = bJA;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in ( numRows, numColumns, diagonalProperty, cSizesPtr, \
                                            aIAPtr, aJAPtr, bIAPtr, bJAPtr )
    {
        IndexType* cSizes = static_cast<IndexType*>( cSizesPtr );

        const IndexType* aIA = static_cast<const IndexType*>( aIAPtr );
        const IndexType* aJA = static_cast<const IndexType*>( aJAPtr );
        const IndexType* bIA = static_cast<const IndexType*>( bIAPtr );
        const IndexType* bJA = static_cast<const IndexType*>( bJAPtr );

        #pragma omp parallel
        {
            IndexType* indexList = new IndexType[numColumns];

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

                cSizes[i] = length;

                // reset indexList for next use

                while( firstCol != END )
                {
                    IndexType nextCol = indexList[firstCol];
                    indexList[firstCol] = NINIT;
                    firstCol = nextCol;
                }

            } // end loop over all rows of input matrices

            delete[] indexList;
        }
    }

    IndexType totalSize = sizes2offsets( cSizes, numRows );

    return totalSize;
}

/* --------------------------------------------------------------------------- */

IndexType MICCSRUtils::matrixMultiplySizes(
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
    // SCAI_REGION( "MIC.CSR.matrixMultiplySizes" )

    SCAI_LOG_INFO( logger,
                   "matrixMutliplySizes for " << m << " x " << n << " matrix" << ", diagonalProperty = " << diagonalProperty )

    // determine the number of entries in output matrix

    IndexType newElems = 0;
    IndexType doubleElems = 0;

    void* cSizesPtr = cSizes;

    const void* aIAPtr = aIA;
    const void* aJAPtr = aJA;
    const void* bIAPtr = bIA;
    const void* bJAPtr = bJA;

    SCAI_ASSERT_ERROR( n <= 8192, "problem too big for add" );

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in ( m, n, diagonalProperty, cSizesPtr, \
                                                aIAPtr, aJAPtr, bIAPtr, bJAPtr ), out( newElems, doubleElems )
    {
        const IndexType NINIT = n + 1; // marks unused colums
        const IndexType END = n + 2; // marks end of list

        IndexType* cSizes = static_cast<IndexType*>( cSizesPtr );

        const IndexType* aIA = static_cast<const IndexType*>( aIAPtr );
        const IndexType* aJA = static_cast<const IndexType*>( aJAPtr );
        const IndexType* bIA = static_cast<const IndexType*>( bIAPtr );
        const IndexType* bJA = static_cast<const IndexType*>( bJAPtr );

        newElems = 0;
        doubleElems = 0;

        #pragma omp parallel reduction( + : newElems, doubleElems )
        {
            IndexType* indexList = new IndexType[n];

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

                        // if ( k < 0 || k >= n ) std::cout << "k = " << k << " illegal" << std::endl;

                        // element a(i,j) an b(j,k) will generate the output element c(i,k)

                        if( indexList[k] == NINIT )
                        {
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

                // std::cout << "row will have " << length << " entries " << std::endl;

                cSizes[i] = length;

                // reset indexList for next use

                while( firstCol != END )
                {
                    IndexType nextCol = indexList[firstCol];
                    indexList[firstCol] = NINIT;
                    firstCol = nextCol;
                }

            } //end loop over all rows of input matrix

            delete[] indexList;
        }
    }

    IndexType totalSize = sizes2offsets( cSizes, m );

    return totalSize;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICCSRUtils::matrixAdd(
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
    // SCAI_REGION( "MIC.CSR.matrixAdd" )

    SCAI_LOG_INFO( logger,
                   "matrixAddJA for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )

    const IndexType NINIT = numColumns + 1;
    const IndexType END = numColumns + 2;

    void* cJAPtr = cJA;
    void* cValuesPtr = cValues;
    const void* cIAPtr = cIA;

    const void* aIAPtr = aIA;
    const void* aJAPtr = aJA;
    const void* aValuesPtr = aValues;
    const void* bIAPtr = bIA;
    const void* bJAPtr = bJA;
    const void* bValuesPtr = bValues;

    const ValueType* alphaPtr = &alpha;
    const ValueType* betaPtr = &beta;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in ( numRows, numColumns, diagonalProperty, \
                                               cJAPtr, cValuesPtr, cIAPtr, \
                                               aIAPtr, aJAPtr, aValuesPtr, bIAPtr, \
											   bJAPtr, bValuesPtr, alphaPtr[0:1], betaPtr[0:1] )
    {
        ValueType* cValues = static_cast<ValueType*>( cValuesPtr );
        IndexType* cJA = static_cast<IndexType*>( cJAPtr );
        const IndexType* cIA = static_cast<const IndexType*>( cIAPtr );

        const IndexType* aIA = static_cast<const IndexType*>( aIAPtr );
        const IndexType* aJA = static_cast<const IndexType*>( aJAPtr );
        const ValueType* aValues = static_cast<const ValueType*>( aValuesPtr );
        const IndexType* bIA = static_cast<const IndexType*>( bIAPtr );
        const IndexType* bJA = static_cast<const IndexType*>( bJAPtr );
        const ValueType* bValues = static_cast<const ValueType*>( bValuesPtr );

        #pragma omp parallel
        {
            IndexType* indexList = new IndexType[numColumns];
            ValueType* valueList = new ValueType[numColumns];

            for( IndexType j = 0; j < numColumns; j++ )
            {
                indexList[j] = NINIT;
                valueList[j] = static_cast<ValueType>(0.0);
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

                    valueList[j] += (*alphaPtr) * aValues[jj];

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

                    valueList[j] += (*betaPtr) * bValues[jj];

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
                    cJA[offset] = i;
                    cValues[offset] = static_cast<ValueType>(0.0);
                    ++offset;
                }

                // fill in csrJA, csrValues and reset indexList, valueList for next use

                while( firstCol != END )
                {
                    IndexType nextCol = indexList[firstCol];
                    ValueType val = valueList[firstCol];

                    indexList[firstCol] = NINIT;
                    valueList[firstCol] = static_cast<ValueType>(0.0); // reset for next time

                    if( diagonalProperty && firstCol == i )
                    {
                        cValues[cIA[i]] = val;
                    }
                    else
                    {
                        cJA[offset] = firstCol;
                        cValues[offset] = val;
                        ++offset;
                    }

                    firstCol = nextCol;
                }

            } //end loop over all rows of input matrix a

            delete[] indexList;
            delete[] valueList;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICCSRUtils::matrixMultiply(
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
    // SCAI_REGION( "MIC.CSR.matrixMultiply" )

    const IndexType NINIT = n + 1;
    const IndexType END = n + 2;

    void* cJAPtr = cJA;
    void* cValuesPtr = cValues;
    const void* cIAPtr = cIA;

    const void* aIAPtr = aIA;
    const void* aJAPtr = aJA;
    const void* aValuesPtr = aValues;
    const void* bIAPtr = bIA;
    const void* bJAPtr = bJA;
    const void* bValuesPtr = bValues;

    const ValueType* alphaPtr = &alpha;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in ( m, n, diagonalProperty, \
                                                cJAPtr, cValuesPtr, cIAPtr, \
                                                aIAPtr, aJAPtr, aValuesPtr, \
						bIAPtr, bJAPtr, bValuesPtr, \
						alphaPtr[0:1] )
    {
        ValueType* cValues = static_cast<ValueType*>( cValuesPtr );
        IndexType* cJA = static_cast<IndexType*>( cJAPtr );
        const IndexType* cIA = static_cast<const IndexType*>( cIAPtr );

        const IndexType* aIA = static_cast<const IndexType*>( aIAPtr );
        const IndexType* aJA = static_cast<const IndexType*>( aJAPtr );
        const ValueType* aValues = static_cast<const ValueType*>( aValuesPtr );
        const IndexType* bIA = static_cast<const IndexType*>( bIAPtr );
        const IndexType* bJA = static_cast<const IndexType*>( bJAPtr );
        const ValueType* bValues = static_cast<const ValueType*>( bValuesPtr );

        // determine the number of entries in output matrix
        #pragma omp parallel
        {
            IndexType* indexList = new IndexType[n];

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
                        }
                    }
                }

                IndexType offset = cIA[i];

                if( diagonalProperty )
                {
                    // first element is reserved for diagonal element
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
                        // diagonal already added before
                    }
                    else
                    {
                        cJA[offset] = firstCol;
                        ++offset;
                    }

                    firstCol = nextCol;
                }

                // make sure that we have still the right offsets

            } //end loop over all rows of input matrix a

            delete[] indexList;
        }

        #pragma omp parallel for

        for( IndexType i = 0; i < m; ++i )
        {
            //loop over all none zero elements of row i of output matrix c
            for( IndexType jj = cIA[i]; jj < cIA[i + 1]; ++jj )
            {
                IndexType j = cJA[jj];
                cValues[jj] = static_cast<ValueType>(0.0);

                if( j == -1 )
                {
                    continue;
                }

                for( IndexType kk = aIA[i]; kk < aIA[i + 1]; ++kk )
                {
                    ValueType b_kj = static_cast<ValueType>(0.0);

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

                cValues[jj] *= (*alphaPtr);
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

// MIC device needs nIndex value
__attribute__( ( target( mic ) ))
static IndexType targetNIndex = nIndex;

__attribute__( ( target( mic ) ))
static IndexType findCol(
    const IndexType columns[],
    const IndexType n,
    const IndexType j,
    IndexType lastPos )
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

    // not found

    return targetNIndex;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICCSRUtils::absMaxDiffRowUnsorted(
    const IndexType n1,
    const IndexType csrJA1[],
    const ValueType csrValues1[],
    const IndexType n2,
    const IndexType csrJA2[],
    const ValueType csrValues2[] )
{
    // No assumption about any sorting in a row

    ValueType val = static_cast<ValueType>(0.0);

    IndexType helpIndex = 0; // some kind of thread-safe global value for findCol

    for( IndexType i1 = 0; i1 < n1; ++i1 )
    {
        ValueType diff = csrValues1[i1];

        IndexType j = csrJA1[i1];
        IndexType i2 = findCol( csrJA2, n2, j, helpIndex );

        if( i2 != targetNIndex )
        {
            diff -= csrValues2[i2];
        }

        diff = common::Math::abs( diff );

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

        if( i1 != targetNIndex )
        {
            continue; // already compare in first loop
        }

        ValueType diff = common::Math::abs( csrValues2[i2] );

        if( diff > val )
        {
            val = diff;
        }
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICCSRUtils::absMaxDiffRowSorted(
    const IndexType n1,
    const IndexType csrJA1[],
    const ValueType csrValues1[],
    const IndexType n2,
    const IndexType csrJA2[],
    const ValueType csrValues2[] )
{
    // Note: the implementation assumes that rows are sorted according to column indexes

    ValueType val = static_cast<ValueType>(0.0);

    IndexType i2 = 0;
    IndexType i1 = 0;

    while( i1 < n1 || i2 < n2 )
    {
        ValueType diff = static_cast<ValueType>(0.0);

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
                    // SCAI_ASSERT_ERROR( csrJA1[i1-1] < csrJA1[i1], "unsorted col indexes at csrJA1[" << i1 << "]" )
                }
            }
            else if( j1 > j2 )
            {
                diff = csrValues2[i2];
                ++i2;

                if( i2 < n2 )
                {
                    // SCAI_ASSERT_ERROR( csrJA2[i2-1] < csrJA2[i2], "unsorted col indexes at csrJA2[" << i2 << "]" )
                }
            }
        }

        diff = common::Math::abs( diff );

        if( diff > val )
        {
            val = diff;
        }
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICCSRUtils::absMaxDiffVal(
    IndexType numRows,
    bool sortedRows,
    const IndexType csrIA1[],
    const IndexType csrJA1[],
    const ValueType csrValues1[],
    const IndexType csrIA2[],
    const IndexType csrJA2[],
    const ValueType csrValues2[] )
{
    SCAI_LOG_INFO( logger,
                   "absMaxDiffVal<" << common::getScalarType<ValueType>() << ">: " << "csr[" << numRows << "], sorted = " << sortedRows )

    ValueType val = 0;

    const void* csrIA1Ptr = csrIA1;
    const void* csrJA1Ptr = csrJA1;
    const void* csrValues1Ptr = csrValues1;

    const void* csrIA2Ptr = csrIA2;
    const void* csrJA2Ptr = csrJA2;
    const void* csrValues2Ptr = csrValues2;

    ValueType* valPtr = &val;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( csrIA1Ptr, csrJA1Ptr, csrValues1Ptr,  \
                                            csrIA2Ptr, csrJA2Ptr, csrValues2Ptr,  \
                                            numRows, sortedRows ), out ( valPtr[0:1] )
    {
        const IndexType* csrIA1 = static_cast<const IndexType*>( csrIA1Ptr );
        const IndexType* csrJA1 = static_cast<const IndexType*>( csrJA1Ptr );
        const ValueType* csrValues1 = static_cast<const ValueType*>( csrValues1Ptr );

        const IndexType* csrIA2 = static_cast<const IndexType*>( csrIA2Ptr );
        const IndexType* csrJA2 = static_cast<const IndexType*>( csrJA2Ptr );
        const ValueType* csrValues2 = static_cast<const ValueType*>( csrValues2Ptr );

        ValueType (*absMaxDiffRow)(
            const IndexType,
            const IndexType[],
            const ValueType[],
            const IndexType,
            const IndexType[],
            const ValueType[] );

        if( sortedRows )
        {
            absMaxDiffRow = MICCSRUtils::absMaxDiffRowSorted<ValueType>;
        }
        else
        {
            absMaxDiffRow = MICCSRUtils::absMaxDiffRowUnsorted<ValueType>;
        }

        *valPtr = 0;

        #pragma omp parallel
        {
            ValueType threadVal = 0;

            #pragma omp for schedule( SCAI_OMP_SCHEDULE )

            for( IndexType i = 0; i < numRows; ++i )
            {
                IndexType offs1 = csrIA1[i];
                IndexType offs2 = csrIA2[i];
                IndexType n1 = csrIA1[i + 1] - offs1;
                IndexType n2 = csrIA2[i + 1] - offs2;

                ValueType maxRow = static_cast<ValueType>(0.0);

                maxRow = absMaxDiffRow( n1, &csrJA1[offs1], &csrValues1[offs1], n2, &csrJA2[offs2],
                                        &csrValues2[offs2] );

                if( maxRow > threadVal )
                {
                    threadVal = maxRow;
                }
            }

            #pragma omp critical
            {
                if( threadVal > *valPtr )
                {
                    *valPtr = threadVal;
                }
            }
        }
    }

    return val;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void MICCSRUtils::Registrator::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_INFO( logger, "register CSSUtils OpenMP-routines for MIC at kernel registry [" << flag << "]" )

    KernelRegistry::set<CSRKernelTrait::sizes2offsets>( sizes2offsets, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::validOffsets>( validOffsets, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::offsets2sizes>( offsets2sizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::hasDiagonalProperty>( hasDiagonalProperty, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixAddSizes>( matrixAddSizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiplySizes>( matrixMultiplySizes, ctx, flag );
}

template<typename ValueType>
void MICCSRUtils::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_INFO( logger, "register CSSUtils OpenMP-routines for MIC at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<CSRKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::sparseGEMV<ValueType> >( sparseGEMV, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixAdd<ValueType> >( matrixAdd, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiply<ValueType> >( matrixMultiply, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::jacobiHaloWithDiag<ValueType> >( jacobiHaloWithDiag, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::absMaxDiffVal<ValueType> >( absMaxDiffVal, ctx, flag );

    /*
      KernelRegistry::set<CSRKernelTrait::convertCSR2CSC<float> >( convertCSR2CSC, ctx, flag );
      KernelRegistry::set<CSRKernelTrait::convertCSR2CSC<double> >( convertCSR2CSC, ctx, flag );

      KernelRegistry::set<CSRKernelTrait::sortRowElements<float> >( sortRowElements, ctx, flag );
      KernelRegistry::set<CSRKernelTrait::sortRowElements<double> >( sortRowElements, ctx, flag );
     */
}

template<typename ValueType, typename OtherValueType>
void MICCSRUtils::RegistratorVO<ValueType, OtherValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_INFO( logger, "register CSSUtils OpenMP-routines for MIC at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )

    KernelRegistry::set<CSRKernelTrait::scaleRows<ValueType, OtherValueType> >( scaleRows, ctx, flag );
}
/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICCSRUtils::RegisterGuard::RegisterGuard()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;

    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_MIC_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARITHMETIC_MIC_LIST, SCAI_ARITHMETIC_MIC_LIST>::call( flag );
}

MICCSRUtils::RegisterGuard::~RegisterGuard()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_MIC_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARITHMETIC_MIC_LIST, SCAI_ARITHMETIC_MIC_LIST>::call( flag );
}

MICCSRUtils::RegisterGuard MICCSRUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */
