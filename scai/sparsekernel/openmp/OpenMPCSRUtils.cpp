/**
 * @file OpenMPCSRUtils.cpp
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
 * @brief Implementation of CSR utilities with OpenMP
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// for dll_import
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>

// local library
#include <scai/sparsekernel/CSRKernelTrait.hpp>

// internal scai libraries
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>
#include <scai/blaskernel/openmp/OpenMPBLAS2.hpp>
#include <scai/blaskernel/openmp/OpenMPLAPACK.hpp>
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/bind.hpp>
#include <scai/common/function.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

// std
#include <vector>

namespace scai
{

using common::scoped_array;
using common::TypeTraits;

using tasking::TaskSyncToken;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( OpenMPCSRUtils::logger, "OpenMP.CSRUtils" )

/** Number of minimal threads for which parallelization is effective. */

static int minThreads = 3;

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::scanSerial( IndexType array[], const IndexType numValues )
{
    SCAI_LOG_DEBUG( logger, "scanSerial: " << numValues << " entries" )
    // In this case we do it just serial, probably faster
    IndexType runningSum = 0;

    for ( IndexType i = 0; i < numValues; i++ )
    {
        IndexType tmp = runningSum;
        runningSum += array[i];
        SCAI_LOG_TRACE( logger, "scan, row = " << i << ", size = " << array[i] << ", offset = " << runningSum )
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
    scoped_array<IndexType> threadCounter( new IndexType[numThreads] );
    SCAI_LOG_DEBUG( logger, "scanParallel: " << numValues << " entries for " << numThreads << " threads" )
    #pragma omp parallel
    {
        IndexType myCounter = 0;
        #pragma omp for schedule(static)

        for ( IndexType i = 0; i < numValues; i++ )
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

        for ( IndexType i = 0; i < numValues; i++ )
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
    SCAI_LOG_INFO( logger, "scan " << numValues << " entries, #threads = " << numThreads )

    if ( numThreads < minThreads )
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
    SCAI_LOG_INFO( logger, "check offset array[ " << n << "] for validity, total = " << total )
    bool validFlag = true;
    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE) reduction( && : validFlag )

    for ( IndexType i = 0; i < n; i++ )
    {
        if ( !( array[i] <= array[i + 1] ) )
        {
            SCAI_LOG_DEBUG( logger, "illegal offsets at pos " << i << ": " << array[i] << " - " << array[i + 1] )
            validFlag = false;
        }
    }

    if ( array[n] != total )
    {
        SCAI_LOG_DEBUG( logger, "last entry in offset array = " << array[n] << " illegal, should be " << total )
        validFlag = false;
    }

    return validFlag;
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::sizes2offsets( IndexType array[], const IndexType numValues )
{
    IndexType totalValues = scan( array, numValues );
    array[numValues] = totalValues;
    SCAI_LOG_INFO( logger, "sizes2offsets, #values = " << numValues << ", total = " << totalValues )
    return totalValues;
}

/* --------------------------------------------------------------------------- */

void OpenMPCSRUtils::offsets2sizes( IndexType sizes[], const IndexType offsets[], const IndexType numRows )
{
    if ( sizes == offsets )
    {
        // when using the same array we do it sequential
        for ( IndexType i = 0; i < numRows; i++ )
        {
            sizes[i] = sizes[i + 1] - sizes[i];
        }
    }
    else
    {
        #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

        for ( IndexType i = 0; i < numRows; i++ )
        {
            sizes[i] = offsets[i + 1] - offsets[i];
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPCSRUtils::offsets2sizesGather(
    IndexType sizes[],
    const IndexType offsets[],
    const IndexType rowIndexes[],
    const IndexType numRows )
{
    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

    for ( IndexType i = 0; i < numRows; i++ )
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
    SCAI_LOG_INFO( logger, "hasDiagonalProperty, #numDiagonals = " << numDiagonals )
    bool diagonalProperty = true;
    #pragma omp parallel for reduction( && : diagonalProperty )

    for ( IndexType i = 0; i < numDiagonals; ++i )
    {
        if ( !diagonalProperty )
        {
            continue;
        }

        if ( csrIA[i] == csrIA[i + 1] )
        {
            diagonalProperty = false;
        }
        else if ( csrJA[csrIA[i]] != i )
        {
            diagonalProperty = false;
        }
    }

    SCAI_LOG_DEBUG( logger, "hasDiagonalProperty = " << diagonalProperty )
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
    SCAI_REGION( "OpenMP.CSR.sortRow" )

    SCAI_LOG_INFO( logger, "sort elements in each of " << numRows << " rows, diagonal flag = " << diagonalFlag )

    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        // use bubble sort as sort algorithm
        const IndexType start = csrIA[i];
        IndexType end = csrIA[i + 1] - 1;
        SCAI_LOG_DEBUG( logger, "row " << i << ": sort " << start << " - " << end )
        bool sorted = false;

        while ( !sorted )
        {
            sorted = true; // will be reset if any wrong order appears
            SCAI_LOG_TRACE( logger, "sort from " << start << " - " << end )

            for ( IndexType jj = start; jj < end; ++jj )
            {
                bool swapIt = false;
                SCAI_ASSERT_DEBUG(
                    csrJA[jj] != csrJA[jj + 1],
                    "row " << i << ": " << jj << ", column index " << csrJA[jj] << ", " << csrJA[ jj + 1] << " appears twice" );

                // if diagonalFlag is set, column i is the smallest one

                if ( diagonalFlag && ( csrJA[jj + 1] == i ) && ( csrJA[jj] != i ) )
                {
                    swapIt = true;
                }
                else if ( diagonalFlag && ( csrJA[jj] == i ) )
                {
                    swapIt = false;
                }
                else
                {
                    swapIt = csrJA[jj] > csrJA[jj + 1];
                }

                if ( swapIt )
                {
                    SCAI_LOG_TRACE( logger, "swap at pos " << jj << " : " << csrJA[jj] << " - " << csrJA[ jj + 1 ] )
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

    for ( IndexType i = 0; i < numRows; ++i )
    {
        const IndexType nzRow = offsets[i + 1] - offsets[i];

        if ( nzRow > 0 )
        {
            counter++;
        }
    }

    SCAI_LOG_INFO( logger, "#non-zero rows = " << counter << ", counted by offsets" )
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

    for ( IndexType i = 0; i < numRows; ++i )
    {
        const IndexType nzRow = offsets[i + 1] - offsets[i];

        if ( nzRow > 0 )
        {
            rowIndexes[counter] = i;
            counter++;
        }
    }

    SCAI_ASSERT_EQUAL_DEBUG( counter, numNonEmptyRows )
    SCAI_LOG_INFO( logger, "#non-zero rows = " << counter << ", set by offsets" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::countNonZeros(
    IndexType sizes[],
    const IndexType ia[],
    const IndexType ja[],
    const ValueType values[],
    const IndexType numRows,
    const ValueType eps,
    const bool diagonalFlag )
{
    SCAI_REGION( "OpenMP.CSRUtils.countNonZeros" )
    SCAI_LOG_INFO( logger, "countNonZeros of CSR<" << TypeTraits<ValueType>::id() << ">( " << numRows
                   << "), eps = " << eps << ", diagonal = " << diagonalFlag )
    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType cnt = 0;

        for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            bool isDiagonal = diagonalFlag && ( ja[jj] == i );
            bool nonZero    = common::Math::abs( values[jj] ) > eps;
            SCAI_LOG_TRACE( logger, "i = " << i << ", j = " << ja[jj] << ", val = " << values[jj]
                            << ", isDiagonal = " << isDiagonal << ", nonZero = " << nonZero )

            if ( nonZero || isDiagonal )
            {
                ++cnt;
            }
        }

        SCAI_LOG_TRACE( logger, "sizes[" << i << "] = " << cnt )
        sizes[i] = cnt;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::compress(
    IndexType newJA[],
    ValueType newValues[],
    const IndexType newIA[],
    const IndexType ia[],
    const IndexType ja[],
    const ValueType values[],
    const IndexType numRows,
    const ValueType eps,
    const bool diagonalFlag )
{
    SCAI_REGION( "OpenMP.CSR.compress" )
    SCAI_LOG_INFO( logger, "compress of CSR<" << TypeTraits<ValueType>::id() << ">( " << numRows
                   << "), eps = " << eps << ", diagonal = " << diagonalFlag )
    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType offs = newIA[i];
        SCAI_LOG_TRACE( logger, "row i: " << ia[i] << ":" << ia[i + 1] << " -> " << newIA[i] << ":" << newIA[i + 1] )

        for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            bool isDiagonal = diagonalFlag && ( ja[jj] == i );
            bool nonZero    = common::Math::abs( values[jj] ) > eps;

            if ( nonZero || isDiagonal )
            {
                newJA[ offs ]     = ja[jj];
                newValues[ offs ] = values[jj];
                ++offs;
            }
        }

        // make sure that filling the compressed data fits to the computed offsets
        SCAI_ASSERT_EQUAL_DEBUG( offs, newIA[i + 1] )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType1, typename ValueType2>
void OpenMPCSRUtils::scaleRows(
    ValueType1 csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const ValueType2 values[] )
{
    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

    for ( IndexType i = 0; i < numRows; ++i )
    {
        ValueType1 tmp = static_cast<ValueType1>( values[i] );

        for ( IndexType j = csrIA[i]; j < csrIA[i + 1]; ++j )
        {
            csrValues[j] *= tmp;
        }
    }
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::getValuePos( const IndexType i, const IndexType j, const IndexType csrIA[], const IndexType csrJA[] )
{
    IndexType pos = nIndex;

    for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
    {
        if ( csrJA[jj] == j )
        {
            pos = jj;
            break;
        }
    }

    return pos;
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::getValuePosCol( IndexType row[], IndexType pos[],
        const IndexType j,
        const IndexType csrIA[], const IndexType numRows,
        const IndexType csrJA[], const IndexType )
{
    SCAI_REGION( "OpenMP.CSRUtils.getValuePosCol" )

    IndexType cnt  = 0;   // counts number of available row entries in column j

    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
        {
            if ( csrJA[jj] == j )
            {
                IndexType k = atomicInc( cnt );
                row[k] = i;
                pos[k] = jj;
                break;
            }
        }
    }

    return cnt;
}

/* --------------------------------------------------------------------------- */

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
    SCAI_REGION( "OpenMP.CSRUtils.CSR2CSC" )
    SCAI_LOG_INFO( logger, "convertCSR2CSC of matrix " << numRows << " x " << numColumns )
    SCAI_ASSERT_EQUAL_DEBUG( numValues, rIA[numRows] )
    // initialization of column counters with 0
    #pragma omp parallel for

    for ( IndexType i = 0; i < numColumns; ++i )
    {
        cscIA[i] = 0;
    }

    // loop over all rows of the row matrix to count columns, not yet OpenMP parallelized
    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        // loop over all none zero elements of column i
        for ( IndexType jj = rIA[i]; jj < rIA[i + 1]; ++jj )
        {
            IndexType j = rJA[jj];
            #pragma omp atomic
            cscIA[j]++;
        }
    }

    sizes2offsets( cscIA, numColumns );
    SCAI_LOG_INFO( logger, "convertCSR2CSC, #num values counted = " << cscIA[ numColumns ] )
    SCAI_ASSERT_EQUAL_DEBUG( numValues, cscIA[numColumns] )
    // temporary copy neeeded of cscIA
    std::vector<IndexType> cscIA1( numColumns );
    #pragma omp parallel for

    for ( IndexType i = 0; i < numColumns; ++i )
    {
        cscIA1[i] = cscIA[i];
    }

    // fill in the array cscJA and cscValues
    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType jj = rIA[i]; jj < rIA[i + 1]; ++jj )
        {
            IndexType j = rJA[jj];
            IndexType k = atomicInc( cscIA1[j] ); // k keeps old value
            cscJA[k] = i;
            cscValues[k] = rValues[jj];
        }
    }
}

/* --------------------------------------------------------------------------- */

// Alternative routine is needed as bind can only deal with up to 9 arguments

template<typename ValueType>
void OpenMPCSRUtils::normalGEMV_s(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads()
                   << ">, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y " )
    // ToDo: for efficiency the following cases should be considered
    // result = y, beta = 1.0 : result += alpha * A * x
    // otherwise: alpha = -1.0, 1.0, beta = 1.0, -1.0, 0.0
    // Note: for beta = 0.0, y could be uninitialized.
    //       0.0 * undefined should deliver 0.0, but valgrind cannot deal with it
    #pragma omp parallel
    {
        // Note: region will be entered by each thread
        SCAI_REGION( "OpenMP.CSR.normalGEMV" )
        #pragma omp for schedule(SCAI_OMP_SCHEDULE)

        for ( IndexType i = 0; i < numRows; ++i )
        {
            ValueType temp = static_cast<ValueType>( 0.0 );

            for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
            {
                IndexType j = csrJA[jj];
                temp += csrValues[jj] * x[j];
            }

            if ( beta == common::constants::ZERO )
            {
                result[i] = alpha * temp;
            }
            else
            {
                result[i] = alpha * temp + beta * y[i];
            }
        }
    }

    if ( SCAI_LOG_TRACE_ON( logger ) )
    {
        std::cout << "NormalGEMV: result = ";

        for ( IndexType i = 0; i < numRows; ++i )
        {
            std::cout << " " << result[i];
        }

        std::cout << std::endl;
    }
}

template<typename ValueType>
void OpenMPCSRUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType SCAI_UNUSED( numColumns ),
    const IndexType SCAI_UNUSED( nnz ),
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        syncToken->run( common::bind( normalGEMV_s<ValueType>, result, alpha, x, beta, y,
                                      numRows, csrIA, csrJA, csrValues ) );
    }
    else
    {
        normalGEMV_s( result, alpha, x, beta, y, numRows, csrIA, csrJA, csrValues );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::normalGEVM_s(
    ValueType result[],
    std::pair<ValueType, const ValueType*> ax,
    std::pair<ValueType, const ValueType*> by,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    normalGEVM( result, ax.first, ax.second, by.first, by.second, numRows, numColumns, csrIA, csrJA, csrValues );
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
    const ValueType csrValues[] )
{
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        // bind takes maximal 9 arguments, so we put (alpha, x) and (beta, y) in a struct
        SCAI_LOG_INFO( logger, "normalGEVM<" << TypeTraits<ValueType>::id() << ", launch it as an asynchronous task" )
        syncToken->run( common::bind( normalGEVM_s<ValueType>, result,
                                      std::pair<ValueType, const ValueType*>( alpha, x ),
                                      std::pair<ValueType, const ValueType*>( beta, y ),
                                      numRows, numColumns, csrIA, csrJA, csrValues ) );
        return;
    }

    SCAI_LOG_INFO( logger,
                   "normalGEVM<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads() << ">, result[" << numColumns << "] = "
                   << alpha << " * x * A + " << beta << " * y " )

    // result := alpha * x * A + beta * y -> result:= beta * y; result += alpha * x * A

    utilskernel::OpenMPUtils::binaryOpScalar1( result, beta, y, numColumns, utilskernel::binary::MULT );

    #pragma omp parallel
    {
        // Note: region will be entered by each thread

        SCAI_REGION( "OpenMP.CSR.normalGEVM" )

        #pragma omp for 

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
            {
                IndexType j = csrJA[jj];
                ValueType v = alpha * csrValues[jj] * x[i];
            
                atomicAdd( result[j], v );
            }
        }
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
    const ValueType csrValues[] )
{
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        syncToken->run( common::bind( sparseGEMV<ValueType>, result, alpha,
                                      x, numNonZeroRows, rowIndexes, csrIA, csrJA, csrValues ) );
        return;
    }

    #pragma omp parallel
    {
        // Note: region will be entered by each thread
        SCAI_REGION( "OpenMP.CSR.sparseGEMV" )
        #pragma omp for schedule(SCAI_OMP_SCHEDULE)

        for ( IndexType ii = 0; ii < numNonZeroRows; ++ii )
        {
            ValueType temp = static_cast<ValueType>( 0.0 );
            IndexType i = rowIndexes[ii];

            for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
            {
                IndexType j = csrJA[jj];
                temp += csrValues[jj] * x[j];
            }

            result[i] += alpha * temp;
        }
    }

    if ( SCAI_LOG_TRACE_ON( logger ) )
    {
        std::cout << "sparseGEMV: result = ";

        for ( IndexType ii = 0; ii < numNonZeroRows; ++ii )
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
    const ValueType csrValues[] )
{
    SCAI_LOG_INFO( logger,
                   "sparseGEVM<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads()
                   << ">, result[" << numColumns << "] += " << alpha << " * x * A" )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_ERROR( logger, "asynchronous execution not supported here" )
    }

    // While GEMV gathers x values needed for each row, GEVM scatters values in result

    #pragma omp parallel
    {
        // Note: region will be entered by each thread

        SCAI_REGION( "OpenMP.CSR.sparseGEVM" )

        #pragma omp for 

        for ( IndexType ii = 0; ii < numNonZeroRows; ++ii )
        {
            IndexType i = rowIndexes[ii];

            for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
            {
                IndexType j = csrJA[jj];
                ValueType v = alpha * csrValues[jj] * x[i];

                atomicAdd( result[j], v );
            }
        }
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
    const ValueType csrValues[] )
{
    SCAI_LOG_INFO( logger,
                   "gemm<" << TypeTraits<ValueType>::id() << ">, " << " result " << m << " x " << n << " CSR " << m << " x " << p )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_ERROR( logger, "asynchronous execution not supported here" )
    }

    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

    for ( IndexType i = 0; i < m; ++i )
    {
        for ( IndexType k = 0; k < n; ++k )
        {
            ValueType temp = static_cast<ValueType>( 0.0 );

            for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
            {
                IndexType j = csrJA[jj];
                // SCAI_ASSERT_DEBUG( j < p , "index j = " << j << " out of range " << p )
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
    const IndexType numRows )
{
    SCAI_LOG_INFO( logger,
                   "jacobi<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", omega = " << omega )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        syncToken->run( common::bind( jacobi<ValueType>, solution, csrIA, csrJA, csrValues, oldSolution, rhs, omega, numRows ) );
        return;
    }

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.CSR.jacobi" )
        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for ( IndexType i = 0; i < numRows; i++ )
        {
            ValueType temp = rhs[i];
            const ValueType diag = csrValues[csrIA[i]];

            for ( IndexType j = csrIA[i] + 1; j < csrIA[i + 1]; j++ )
            {
                temp -= csrValues[j] * oldSolution[csrJA[j]];
            }

            // here we take advantange of a good branch precondiction

            if ( omega == common::constants::ONE )
            {
                solution[i] = temp / diag;
            }
            else if ( omega == 0.5 )
            {
                solution[i] = omega * ( temp / diag + oldSolution[i] );
            }
            else
            {
                solution[i] = omega * ( temp / diag ) + ( static_cast<ValueType>( 1.0 ) - omega ) * oldSolution[i];
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
    SCAI_LOG_INFO( logger,
                   "jacobiHalo<" << TypeTraits<ValueType>::id() << ">" << ", #rows (not empty) = " << numNonEmptyRows << ", omega = " << omega );
    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.CSR.jacabiHalo" )
        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for ( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            IndexType i = ii; // default: rowIndexes == NULL stands for identity

            if ( haloRowIndexes )
            {
                i = haloRowIndexes[ii];
            }

            ValueType temp = static_cast<ValueType>( 0.0 );
            const ValueType diag = localValues[localIA[i]];

            for ( IndexType j = haloIA[i]; j < haloIA[i + 1]; j++ )
            {
                temp += haloValues[j] * oldSolution[haloJA[j]];
            }

            if ( omega == common::constants::ONE )
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
    SCAI_LOG_INFO( logger,
                   "jacobiHaloWithDiag<" << TypeTraits<ValueType>::id() << ">" << ", #rows (not empty) = " << numNonEmptyRows << ", omega = " << omega );
    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.CSR.jacabiHaloWithDiag" )
        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for ( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            IndexType i = ii; // default: rowIndexes == NULL stands for identity

            if ( haloRowIndexes )
            {
                i = haloRowIndexes[ii];
            }

            ValueType temp = static_cast<ValueType>( 0.0 );
            const ValueType diag = localDiagValues[i];

            for ( IndexType j = haloIA[i]; j < haloIA[i + 1]; j++ )
            {
                temp += haloValues[j] * oldSolution[haloJA[j]];
            }

            if ( omega == common::constants::ONE )
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
void OpenMPCSRUtils::decomposition(
    ValueType* const solution,
    const IndexType* csrIA,
    const IndexType* csrJA,
    const ValueType* csrValues,
    const ValueType* rhs,
    const IndexType numRows,
    const IndexType /*nnz*/,
    const bool /*isSymmetic*/ )
{
    // current workaround without MKL: inverse solver of dense matrix

    common::scoped_array<ValueType> denseA( new ValueType[ numRows * numRows ] );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numRows; ++j )
        {
            denseA[ i * numRows + j ] = static_cast<ValueType>( 0 );
        }
    }

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType jj = csrIA[i]; jj < csrIA[i+1]; ++jj )
        {
            SCAI_ASSERT_VALID_INDEX_DEBUG( csrJA[jj], numRows, "illegal col index, square matrix assumed" );
            denseA[ i * numRows + csrJA[jj] ] = csrValues[ jj ];
        }
    }

    SCAI_LOG_INFO( logger, "decomposition<" << common::TypeTraits<ValueType>::id() << "> via inverse" )

    // now call inverse solver of lapack 

    blaskernel::OpenMPLAPACK::getinv( numRows, denseA.get(), numRows );

    ValueType alpha = 1;
    ValueType beta  = 0;
    IndexType inc1 = 1;

    // initialize solution even if beta == 0

    for ( IndexType i = 0; i < numRows; ++i )
    {
        solution[i] = beta;  
    }

    blaskernel::OpenMPBLAS2::gemv( CblasRowMajor, CblasNoTrans, 
                                   numRows, numRows, alpha, denseA.get(), numRows, rhs, inc1, beta, solution, inc1 );

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
    SCAI_LOG_INFO( logger,
                   "matrixAddSizes for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )
    SCAI_REGION( "OpenMP.CSR.matrixAddSizes" )
    // determine the number of entries in output matrix
    const IndexType NINIT = numColumns + 1; // marks unused colums
    const IndexType END = numColumns + 2; // marks end of list
    #pragma omp parallel
    {
        scoped_array<IndexType> indexList( new IndexType[numColumns] );

        for ( IndexType j = 0; j < numColumns; j++ )
        {
            indexList[j] = NINIT;
        }

        #pragma omp for

        for ( IndexType i = 0; i < numRows; ++i )
        {
            IndexType length = 0;
            IndexType firstCol = END;
            bool diagonal = false; // will be set if diagonal element is set

            // loop over all none zero elements of row i of input matrix a

            for ( IndexType jj = aIA[i]; jj < aIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix a
                // so we are at position a(i,j)
                IndexType j = aJA[jj];

                // element a(i,j) will generate an output element c(i,j)

                if ( indexList[j] == NINIT )
                {
                    SCAI_LOG_TRACE( logger, "entry for [" << i << "," << j << "] by a" )
                    // Add column position j to the indexList
                    indexList[j] = firstCol;
                    firstCol = j;
                    ++length;

                    if ( j == i )
                    {
                        diagonal = true;
                    }
                }
            }

            for ( IndexType jj = bIA[i]; jj < bIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix b
                // so we are at position b(i,j)
                IndexType j = bJA[jj];

                // element a(i,j) will generate an output element c(i,j)

                if ( indexList[j] == NINIT )
                {
                    SCAI_LOG_TRACE( logger, "entry for [" << i << "," << j << "] by b" )
                    // Add column position j to the indexList
                    indexList[j] = firstCol;
                    firstCol = j;
                    ++length;

                    if ( j == i )
                    {
                        diagonal = true;
                    }
                }
            }

            if ( diagonalProperty && !diagonal )
            {
                ++length; // diagaonal needed, but not filled yet
            }

            // so we have now the correct length
            SCAI_LOG_TRACE( logger, "row " << i << " will have " << length << " entries" )
            cSizes[i] = length;

            // reset indexList for next use

            while ( firstCol != END )
            {
                SCAI_LOG_TRACE( logger, "entry [" << i << "," << firstCol << "] will be set" )
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
    SCAI_REGION( "OpenMP.CSR.matrixMultiplySizes" )
    SCAI_LOG_INFO( logger,
                   "matrixMutliplySizes for " << m << " x " << n << " matrix" << ", diagonalProperty = " << diagonalProperty )
    // determine the number of entries in output matrix
    const IndexType NINIT = n + 1; // marks unused colums
    const IndexType END = n + 2; // marks end of list
    IndexType newElems = 0;
    IndexType doubleElems = 0;
    #pragma omp parallel reduction( + : newElems, doubleElems )
    {
        scoped_array<IndexType> indexList( new IndexType[n] );

        for ( IndexType j = 0; j < n; j++ )
        {
            indexList[j] = NINIT;
        }

        #pragma omp for

        for ( IndexType i = 0; i < m; ++i )
        {
            IndexType length = 0;
            IndexType firstCol = END;
            bool diagonal = false; // will be set if diagonal element is set

            // loop over all none zero elements of row i of input matrix a

            for ( IndexType jj = aIA[i]; jj < aIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix a
                // so we are at position a(i,j)
                IndexType j = aJA[jj];

                // loop over all none zero elements of row j of input matrix b
                // that is the row of b that corresponds to the column of the current
                // element of a

                for ( IndexType kk = bIA[j]; kk < bIA[j + 1]; ++kk )
                {
                    // k is the column of none zero element kk of row j of input matrix b
                    // so we are looking at position b(j,k)

                    IndexType k = bJA[kk];

                    SCAI_ASSERT_VALID_INDEX_DEBUG( k, n, "invalid value" )

                    // element a(i,j) an b(j,k) will generate the output element c(i,k)

                    if ( indexList[k] == NINIT )
                    {
                        SCAI_LOG_TRACE( logger, "entry for [" << i << "," << k << "]" )
                        // Add column position k to the indexList
                        indexList[k] = firstCol;
                        firstCol = k;
                        ++length;

                        if ( k == i )
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

            if ( diagonalProperty && !diagonal )
            {
                ++length; // diagaonal needed, but not filled yet
            }

            // so we have now the correct length
            SCAI_LOG_TRACE( logger, "row " << i << " will have " << length << " entries" )
            cSizes[i] = length;

            // reset indexList for next use

            while ( firstCol != END )
            {
                SCAI_LOG_TRACE( logger, "entry [" << i << "," << firstCol << "] will be set" )
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
    SCAI_REGION( "OpenMP.CSR.matrixAdd" )
    SCAI_LOG_INFO( logger,
                   "matrixAddJA for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )
    const IndexType NINIT = numColumns + 1;
    const IndexType END = numColumns + 2;
    // determine the number of entries in output matrix
    #pragma omp parallel
    {
        scoped_array<IndexType> indexList( new IndexType[numColumns] );
        scoped_array<ValueType> valueList( new ValueType[numColumns] );

        for ( IndexType j = 0; j < numColumns; j++ )
        {
            indexList[j] = NINIT;
            valueList[j] = static_cast<ValueType>( 0.0 );
        }

        #pragma omp for

        for ( IndexType i = 0; i < numRows; ++i )
        {
            IndexType length = 0;
            IndexType firstCol = END;

            for ( IndexType jj = aIA[i]; jj < aIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix a
                // so we are at position a(i,j)
                IndexType j = aJA[jj];
                valueList[j] += alpha * aValues[jj];
                SCAI_LOG_TRACE( logger, "entry for [" << i << "," << j << "] by a" << ", new val = " << valueList[j] )

                // element a(i,j) will generate an output element c(i,j)

                if ( indexList[j] == NINIT )
                {
                    // Add column position j to the indexList
                    indexList[j] = firstCol;
                    firstCol = j;
                    ++length;
                }
            }

            for ( IndexType jj = bIA[i]; jj < bIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix b
                // so we are at position b(i,j)
                IndexType j = bJA[jj];
                valueList[j] += beta * bValues[jj];
                SCAI_LOG_TRACE( logger, "entry for [" << i << "," << j << "] by b" << ", new val = " << valueList[j] )

                // element b(i,j) will generate an output element c(i,j)

                if ( indexList[j] == NINIT )
                {
                    // Add column position j to the indexList
                    indexList[j] = firstCol;
                    firstCol = j;
                    ++length;
                }
            }

            IndexType offset = cIA[i];

            if ( diagonalProperty )
            {
                // first element is reserved for diagonal element
                SCAI_LOG_TRACE( logger, "entry for [" << i << "," << i << "] as diagonal" )
                cJA[offset] = i;
                cValues[offset] = static_cast<ValueType>( 0.0 );
                ++offset;
            }

            // fill in csrJA, csrValues and reset indexList, valueList for next use

            while ( firstCol != END )
            {
                IndexType nextCol = indexList[firstCol];
                ValueType val = valueList[firstCol];
                indexList[firstCol] = NINIT;
                valueList[firstCol] = static_cast<ValueType>( 0.0 ); // reset for next time

                if ( diagonalProperty && firstCol == i )
                {
                    SCAI_LOG_TRACE( logger, "diagonal already added before" )
                    SCAI_LOG_TRACE( logger, "entry for [" << i << "," << i << "] = " << val )
                    cValues[cIA[i]] = val;
                }
                else
                {
                    SCAI_LOG_TRACE( logger, "entry for [" << i << "," << firstCol << "] = " << val )
                    cJA[offset] = firstCol;
                    cValues[offset] = val;
                    ++offset;
                }

                firstCol = nextCol;
            }

            // make sure that we have still the right offsets
            SCAI_ASSERT_EQUAL_DEBUG( offset, cIA[i + 1] )
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
    SCAI_LOG_INFO( logger,
                   "matrixMutliplyJA for " << numRows << " x " << numColumns << " matrix" << ", diagonalProperty = " << diagonalProperty )
    const IndexType NINIT = numColumns + 1;
    const IndexType END = numColumns + 2;
    // determine the number of entries in output matrix
    #pragma omp parallel
    {
        scoped_array<IndexType> indexList( new IndexType[numColumns] );

        for ( IndexType j = 0; j < numColumns; j++ )
        {
            indexList[j] = NINIT;
        }

        #pragma omp for

        for ( IndexType i = 0; i < numRows; ++i )
        {
            IndexType length = 0;
            IndexType firstCol = END;

            // loop over all none zero elements of row i of input matrix a

            for ( IndexType jj = aIA[i]; jj < aIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix a
                // so we are at position a(i,j)
                IndexType j = aJA[jj];

                // loop over all none zero elements of row j of input matrix b
                // that is the row of b that corresponds to the column of the current
                // element of a

                for ( IndexType kk = bIA[j]; kk < bIA[j + 1]; ++kk )
                {
                    // k is the column of none zero element kk of row j of input matrix b
                    // so we are looking at position b(j,k)
                    IndexType k = bJA[kk];

                    // element a(i,j) an b(j,k) will generate the output element c(i,k)

                    if ( indexList[k] == NINIT )
                    {
                        // Add column position k to the indexList
                        indexList[k] = firstCol;
                        firstCol = k;
                        ++length;
                        SCAI_LOG_TRACE( logger, "entry for [" << i << "," << k << "]" )
                    }
                }
            }

            IndexType offset = cIA[i];

            if ( diagonalProperty )
            {
                // first element is reserved for diagonal element
                SCAI_LOG_TRACE( logger, "entry for [" << i << "," << i << "] as diagonal" )
                cJA[offset] = i;
                ++offset;
            }

            // fill in csrJA and reset indexList for next use

            while ( firstCol != END )
            {
                IndexType nextCol = indexList[firstCol];
                indexList[firstCol] = NINIT;

                if ( diagonalProperty && firstCol == i )
                {
                    SCAI_LOG_TRACE( logger, "diagonal already added before" )
                }
                else
                {
                    SCAI_LOG_TRACE( logger, "entry for [" << i << "," << firstCol << "]" )
                    cJA[offset] = firstCol;
                    ++offset;
                }

                firstCol = nextCol;
            }

            // make sure that we have still the right offsets
            SCAI_ASSERT_EQUAL_DEBUG( offset, cIA[i + 1] )
        } //end loop over all rows of input matrix a
    }
}

/* --------------------------------------------------------------------------- */

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
    //TODO: Rewrite this!
    const IndexType NINIT = n + 1;
    const IndexType END = n + 2;
    // determine the number of entries in output matrix
    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.CSR.matrixMultiply1" )

        scoped_array<IndexType> indexList( new IndexType[n] );

        for ( IndexType j = 0; j < n; j++ )
        {
            indexList[j] = NINIT;
        }

        #pragma omp for

        for ( IndexType i = 0; i < m; ++i )
        {
            IndexType length = 0;
            IndexType firstCol = END;

            // loop over all none zero elements of row i of input matrix a

            for ( IndexType jj = aIA[i]; jj < aIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix a
                // so we are at position a(i,j)
                IndexType j = aJA[jj];

                // loop over all none zero elements of row j of input matrix b
                // that is the row of b that corresponds to the column of the current
                // element of a

                for ( IndexType kk = bIA[j]; kk < bIA[j + 1]; ++kk )
                {
                    // k is the column of none zero element kk of row j of input matrix b
                    // so we are looking at position b(j,k)
                    IndexType k = bJA[kk];

                    // element a(i,j) an b(j,k) will generate the output element c(i,k)

                    if ( indexList[k] == NINIT )
                    {
                        // Add column position k to the indexList
                        indexList[k] = firstCol;
                        firstCol = k;
                        ++length;
                        SCAI_LOG_TRACE( logger, "entry for [" << i << "," << k << "]" )
                    }
                }
            }

            IndexType offset = cIA[i];

            if ( diagonalProperty )
            {
                // first element is reserved for diagonal element
                SCAI_LOG_TRACE( logger, "entry for [" << i << "," << i << "] as diagonal" )
                cJA[offset] = i;
                ++offset;
            }

            // fill in csrJA and reset indexList for next use

            while ( firstCol != END )
            {
                IndexType nextCol = indexList[firstCol];
                indexList[firstCol] = NINIT;

                if ( diagonalProperty && firstCol == i )
                {
                    SCAI_LOG_TRACE( logger, "diagonal already added before" )
                }
                else
                {
                    SCAI_LOG_TRACE( logger, "entry for [" << i << "," << firstCol << "]" )
                    cJA[offset] = firstCol;
                    ++offset;
                }

                firstCol = nextCol;
            }

            // make sure that we have still the right offsets
            // SCAI_ASSERT_EQUAL_DEBUG( offset, cIA[i+1] )
        } //end loop over all rows of input matrix a
    }

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.CSR.matrixMultiply2" )

        // temporary array for cRow, allocated by each thread

        scoped_array<ValueType> cRow( new ValueType[n] );

        for ( IndexType j = 0; j < n; j++ )
        {
            cRow[j] = static_cast<ValueType>( 0 );
        }

        #pragma omp for

        for ( IndexType i = 0; i < m; ++i )
        {
            // compute row of result matrix C, loop over all non-zero element of row i of A

            for ( IndexType jj = aIA[i]; jj < aIA[i + 1]; ++jj )
            {
                IndexType j = aJA[jj];

                for ( IndexType kk = bIA[j]; kk < bIA[j + 1]; ++kk )
                {
                    IndexType k = bJA[kk];

                    // element a(i,j) an b(j,k) will be added to  c(i,k)

                    cRow[k] += aValues[jj] * bValues[kk];
                }
            }

            // now take all values from the row as needed

            for ( IndexType jj = cIA[i]; jj < cIA[i + 1]; ++jj )
            {
                IndexType j = cJA[jj];

                if ( j == nIndex )
                {
                    continue;
                }

                cValues[jj] = alpha * cRow[j];

                cRow[j] = static_cast<ValueType>( 0 );
            }

            // just make sure that we really got all non-zero entries

            // for ( IndexType j = 0; j < n; j++ )
            // {
            //     SCAI_ASSERT_EQ_ERROR( cRow[j], ValueType( 0 ), "serious mismatch" )
            // }
        }
    }
}

/* --------------------------------------------------------------------------- */

static IndexType findCol( const IndexType columns[], const IndexType n, const IndexType j, IndexType lastPos )
{
    // lastPos is index from where to start
    for ( IndexType i = 0; i < n; ++i )
    {
        if ( columns[lastPos] == j )
        {
            // found, return lastPos, increment it for next search
            return lastPos++;
        }

        lastPos++;

        if ( lastPos > n )
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

    for ( IndexType i1 = 0; i1 < n1; ++i1 )
    {
        ValueType diff = csrValues1[i1];
        IndexType j = csrJA1[i1];
        IndexType i2 = findCol( csrJA2, n2, j, helpIndex );

        if ( i2 < n2 )
        {
            diff -= csrValues2[i2];
        }

        diff = common::Math::abs( diff );

        if ( diff > val )
        {
            val = diff;
        }
    }

    // check for row elements in 2nd matrix that are not available in first one
    helpIndex = 0;

    for ( IndexType i2 = 0; i2 < n2; ++i2 )
    {
        IndexType j = csrJA2[i2];
        IndexType i1 = findCol( csrJA1, n1, j, helpIndex );

        if ( i1 < n1 )
        {
            continue; // already compare in first loop
        }

        ValueType diff = common::Math::abs( csrValues2[i2] );

        if ( diff > val )
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

    while ( i1 < n1 || i2 < n2 )
    {
        ValueType diff = 0;

        if ( i1 >= n1 ) // row of matrix1 completely traversed
        {
            diff = csrValues2[i2];
            ++i2;
        }
        else if ( i2 >= n2 ) // row of matrix2 completely traversed
        {
            diff = csrValues1[i1];
            ++i1;
        }
        else
        {
            // so we have still values in both rows
            IndexType j1 = csrJA1[i1];
            IndexType j2 = csrJA2[i2];

            if ( j1 == j2 )
            {
                diff = csrValues1[i1] - csrValues2[i2];
                ++i1;
                ++i2;
            }
            else if ( j1 < j2 )
            {
                diff = csrValues1[i1];
                ++i1;

                if ( i1 < n1 )
                {
                    SCAI_ASSERT_ERROR( csrJA1[i1 - 1] < csrJA1[i1], "unsorted col indexes at csrJA1[" << i1 << "]" )
                }
            }
            else if ( j1 > j2 )
            {
                diff = csrValues2[i2];
                ++i2;

                if ( i2 < n2 )
                {
                    SCAI_ASSERT_ERROR( csrJA2[i2 - 1] < csrJA2[i2], "unsorted col indexes at csrJA2[" << i2 << "]" )
                }
            }
        }

        diff = common::Math::abs( diff );

        if ( diff > val )
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
    SCAI_LOG_INFO( logger,
                   "absMaxDiffVal<" << TypeTraits<ValueType>::id() << ">: " << "csr[" << numRows << "], sorted = " << sortedRows )
    ValueType ( *absMaxDiffRow )(
        const IndexType,
        const IndexType[],
        const ValueType[],
        const IndexType,
        const IndexType[],
        const ValueType[] );

    if ( sortedRows )
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
        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for ( IndexType i = 0; i < numRows; ++i )
        {
            IndexType offs1 = csrIA1[i];
            IndexType offs2 = csrIA2[i];
            IndexType n1 = csrIA1[i + 1] - offs1;
            IndexType n2 = csrIA2[i + 1] - offs2;
            ValueType maxRow = static_cast<ValueType>( 0.0 );
            maxRow = absMaxDiffRow( n1, &csrJA1[offs1], &csrValues1[offs1], n2, &csrJA2[offs2], &csrValues2[offs2] );

            if ( maxRow > threadVal )
            {
                threadVal = maxRow;
            }
        }

        #pragma omp critical
        {
            SCAI_LOG_TRACE( logger, "max val of thread  = " << threadVal << ", global was " << val )

            if ( threadVal > val )
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

void OpenMPCSRUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::context::ContextType ctx = common::context::Host;
    SCAI_LOG_DEBUG( logger, "register CSRUtils OpenMP-routines for Host at kernel registry [" << flag << "]" )
    KernelRegistry::set<CSRKernelTrait::getValuePos>( getValuePos, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::getValuePosCol>( getValuePosCol, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::sizes2offsets>( sizes2offsets, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::offsets2sizes>( offsets2sizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::validOffsets>( validOffsets, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::hasDiagonalProperty>( hasDiagonalProperty, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixAddSizes>( matrixAddSizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiplySizes>( matrixMultiplySizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiplyJA>( matrixMultiplyJA, ctx, flag );
}

template<typename ValueType>
void OpenMPCSRUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::context::ContextType ctx = common::context::Host;
    SCAI_LOG_DEBUG( logger, "register CSRUtils OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<CSRKernelTrait::convertCSR2CSC<ValueType> >( convertCSR2CSC, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::sortRowElements<ValueType> >( sortRowElements, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::sparseGEMV<ValueType> >( sparseGEMV, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::normalGEVM<ValueType> >( normalGEVM, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::sparseGEVM<ValueType> >( sparseGEVM, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::gemm<ValueType> >( gemm, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixAdd<ValueType> >( matrixAdd, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiply<ValueType> >( matrixMultiply, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::jacobiHaloWithDiag<ValueType> >( jacobiHaloWithDiag, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::absMaxDiffVal<ValueType> >( absMaxDiffVal, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::countNonZeros<ValueType> >( countNonZeros, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::compress<ValueType> >( compress, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::decomposition<ValueType> >( decomposition, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void OpenMPCSRUtils::RegistratorVO<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::context::ContextType ctx = common::context::Host;
    SCAI_LOG_DEBUG( logger, "register CSRUtils OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )
    KernelRegistry::set<CSRKernelTrait::scaleRows<ValueType, OtherValueType> >( scaleRows, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPCSRUtils::OpenMPCSRUtils()
{
    SCAI_LOG_INFO( logger, "register CSRUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_HOST_LIST, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

OpenMPCSRUtils::~OpenMPCSRUtils()
{
    SCAI_LOG_INFO( logger, "unregister CSRUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_HOST_LIST, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPCSRUtils OpenMPCSRUtils::guard;

} /* end namespace sparsekernel */

} /* end namespace scai */
