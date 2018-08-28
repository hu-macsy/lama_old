/**
 * @file OpenMPCSRUtils.cpp
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
 * @brief Implementation of CSR utilities with OpenMP
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// for dll_import
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/sparsekernel/openmp/BuildSparseIndexes.hpp>

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
#include <scai/common/macros/unused.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

// std
#include <vector>
#include <memory>
#include <functional>

using std::unique_ptr;

namespace scai
{

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
    unique_ptr<IndexType[]> threadCounter( new IndexType[numThreads] );
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
    #pragma omp parallel for reduction( && : validFlag )

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
        #pragma omp parallel for 

        for ( IndexType i = 0; i < numRows; i++ )
        {
            sizes[i] = offsets[i + 1] - offsets[i];
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPCSRUtils::gatherSizes(
    IndexType sizes[],
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType rowIndexes[],
    const IndexType nIndexes )
{
    #pragma omp parallel for 

    for ( IndexType i = 0; i < nIndexes; i++ )
    {
        IndexType row = rowIndexes[i];
        SCAI_ASSERT_VALID_INDEX_DEBUG( row, numRows, "illegal row index" );
        sizes[i] = csrIA[row + 1] - csrIA[row];
    }
}

/* --------------------------------------------------------------------------- */

bool OpenMPCSRUtils::hasDiagonalProperty(
    const IndexType numDiagonals,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const bool )
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

        bool found = false;

        // ToDo: binary search on sorted rows might be more efficient

        for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
        {
            if ( csrJA[jj] == i )
            {
                found = true;
                break;
            }
        }
   
        if ( !found )
        {
            diagonalProperty = false;
        }
    }

    SCAI_LOG_DEBUG( logger, "hasDiagonalProperty = " << diagonalProperty )

    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

bool OpenMPCSRUtils::hasSortedRows(
    const IndexType csrIA[],
    const IndexType csrJA[],
    const IndexType numRows,
    const IndexType,
    const IndexType )
{
    bool isSorted = true;

    SCAI_LOG_DEBUG( logger, "hasSortedRows: #rows = " << numRows )

    #pragma omp parallel for
    for ( IndexType i = 0; i < numRows; ++i )
    {
        if ( !isSorted )
        {
            continue;   // might be break, but not allowed for parallel for
        }

        IndexType start = csrIA[i];
        const IndexType end = csrIA[i + 1];

        SCAI_LOG_TRACE( logger, "check row " << i << " : " << start << " - " << end )

        if ( start + 1 >= end )
        {
            continue;  // less than 2 elements, is sorted
        }

        if ( ! utilskernel::OpenMPUtils::isSorted( &csrJA[start], end - start, common::CompareOp::LT ) )
        {
            isSorted = false;
        }
    }

    return isSorted;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::sortRows(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType,
    const IndexType nnz )
{
    // argument numColumns is not needed here, we make no checks here

    if ( nnz != csrIA[numRows] )
    {
        COMMON_THROWEXCEPTION( "serious error for arguments, nnz = " << nnz << ", csrIA[" << numRows << " ] = " << csrIA[ numRows ] )
    }

    SCAI_REGION( "OpenMP.CSR.sortRows" )

    SCAI_LOG_INFO( logger, "sort elements in each of " << numRows << " rows" )

    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType start = csrIA[i];

        const IndexType end = csrIA[i + 1];

        SCAI_LOG_TRACE( logger, "sort row " << i << ", start = " << start << ", end = " << end )

        if ( start + 1 < end )
        { 
            // at least two elements in row, call quicksort implementation of utilskernel

            utilskernel::OpenMPUtils::sortInPlace( &csrJA[start], &csrValues[start], end - start, true );
        }
    }
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::nonEmptyRows( IndexType rowIndexes[], const IndexType csrIA[], const IndexType numRows )
{
    IndexType counter = 0;

    if ( rowIndexes == NULL )
    {
        #pragma omp parallel for reduction( + : counter )
       
        for ( IndexType i = 0; i < numRows; ++i )
        {
            if ( csrIA[i + 1] > csrIA[i]  )
            {
                counter++;
            }
        }
    }
    else
    {
        // OpenMP parallelization not possible if indexes are needed

        for ( IndexType i = 0; i < numRows; ++i )
        {
            if ( csrIA[i + 1] > csrIA[i]  )
            {
                rowIndexes[counter++] = i;
            }
        }
    }

    SCAI_LOG_INFO( logger, "#non-empty rows = " << counter << ", counted by offsets" )

    return counter;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::countNonZeros(
    IndexType sizes[],
    const IndexType ia[],
    const IndexType[],
    const ValueType values[],
    const IndexType numRows,
    const RealType<ValueType> eps )
{
    SCAI_REGION( "OpenMP.CSRUtils.countNonZeros" )

    SCAI_LOG_INFO( logger, "countNonZeros of CSR<" << TypeTraits<ValueType>::id() << ">( " 
                            << numRows << "), eps = " << eps  )
    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType cnt = 0;

        for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            if ( common::Math::abs( values[jj] ) <= eps )
            {
                continue;  // skip this zero element
            }

            ++cnt;
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
    const RealType<ValueType> eps )
{
    SCAI_REGION( "OpenMP.CSR.compress" )

    SCAI_LOG_INFO( logger, "compress of CSR<" << TypeTraits<ValueType>::id() << ">( " 
                            << numRows << "), eps = " << eps )

    #pragma omp parallel for
    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType offs = newIA[i];

        SCAI_LOG_TRACE( logger, "row i: " << ia[i] << ":" << ia[i + 1] << " -> " << newIA[i] << ":" << newIA[i + 1] )

        for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            if ( common::Math::abs( values[jj] ) <= eps )
            {
                continue;  // skip this zero element
            }

            // take over the non-zero element

            newJA[ offs ]     = ja[jj];
            newValues[ offs ] = values[jj];
            ++offs;
        }

        // make sure that filling the compressed data fits to the computed offsets
        SCAI_ASSERT_EQUAL_ERROR( offs, newIA[i + 1] )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::getDiagonal(
    ValueType diagonal[],
    const IndexType numDiagonals,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    const bool )
{
    #pragma omp parallel for 

    for ( IndexType i = 0; i < numDiagonals; ++i )
    {
        diagonal[i] = ValueType( 0 );

        for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
        {
            if ( csrJA[jj] == i )
            {
                diagonal[i] = csrValues[jj];
                break;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool OpenMPCSRUtils::setDiagonal(
    ValueType csrValues[],
    const ValueType diagonal,
    const IndexType numDiagonals,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const bool )
{
    bool okay = true;

    #pragma omp parallel for 

    for ( IndexType i = 0; i < numDiagonals; ++i )
    {
        bool found = false;

        for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
        {   
            if ( csrJA[jj] == i )
            {
                csrValues[jj] = diagonal;
                found = true;
                break;
            }
        }

        if ( !found )
        {
            okay = false;
        }
    }

    if ( diagonal == ValueType( 0 ) )
    {
        okay = true;  // missing diagonal entries have the correct value
    }

    return okay;
}
/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool OpenMPCSRUtils::setDiagonalV(
    ValueType csrValues[],
    const ValueType diagonal[],
    const IndexType numDiagonals,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const bool )
{
    bool okay = true;

    const ValueType ZERO = 0;  

    #pragma omp parallel for 

    for ( IndexType i = 0; i < numDiagonals; ++i )
    {
        bool found = false;

        for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
        {
            if ( csrJA[jj] == i )
            {
                csrValues[jj] = diagonal[i];
                found = true;
                break;
            }
        }

        if ( !found && diagonal[i] != ZERO )
        {
            okay = false;  // so we could not set a non-zero value
        }
    }

    return okay;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::setRows(
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const ValueType values[],
    const common::BinaryOp op )
{
    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; ++i )
    {
        if ( op == common::BinaryOp::MULT )
        {
            for ( IndexType j = csrIA[i]; j < csrIA[i + 1]; ++j )
            {
                csrValues[j] *= values[i];
            }
        }
        else
        {
            for ( IndexType j = csrIA[i]; j < csrIA[i + 1]; ++j )
            {
                csrValues[j] = common::applyBinary( csrValues[j], op, values[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::setColumns(
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const IndexType numRows,
    const ValueType values[],
    const common::BinaryOp op )
{
    SCAI_REGION( "OpenMP.CSR.setColumns" )

    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; ++i )
    {
        if ( op == common::BinaryOp::MULT )
        {
            for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
            {
                csrValues[jj] *= values[csrJA[jj]];
            }
        }
        else
        {
            for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
            {
                csrValues[jj] = common::applyBinary( csrValues[jj], op, values[csrJA[jj]] );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::getValuePos( const IndexType i, const IndexType j, const IndexType csrIA[], const IndexType csrJA[] )
{
    IndexType pos = invalidIndex;

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

IndexType OpenMPCSRUtils::getPosDiagonal(
    IndexType pos[],
    const IndexType numDiagonals,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const bool )
{
    IndexType numFoundDiagonals = 0;

    #pragma omp parallel
    {
        IndexType threadNumDiagonals = 0;
    
        #pragma omp for

        for ( IndexType i = 0; i < numDiagonals; ++i )
        {
            pos[i] = invalidIndex;

            for ( IndexType jj = csrIA[i]; jj < csrIA[i+1]; ++jj )
            {
                if ( csrJA[jj] == i )
                {
                    pos[i] = jj; 
                    threadNumDiagonals++;
                    break;
                }
            }
        }

        atomicAdd( numFoundDiagonals, threadNumDiagonals );
    }
 
    return numFoundDiagonals;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType OpenMPCSRUtils::shiftDiagonal(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType numDiagonals,
    const IndexType csrIA[] )
{
    IndexType numFoundDiagonals = 0;

    #pragma omp parallel
    {
        IndexType threadNumDiagonals = 0;
    
        #pragma omp for

        for ( IndexType i = 0; i < numDiagonals; ++i )
        {
            bool found = false;

            IndexType start = csrIA[i];
            IndexType end   = csrIA[i+1] - 1;

            if ( end < start )
            {
                continue;    // not found
            }

            if ( csrJA[start] == i )
            {
                threadNumDiagonals++;   // diagonal element is already first
                continue;
            }

            ValueType diagonalValue;

            // traverse reverse

            while ( end > start )
            {
                // check if it is the diagonal element, save the diagonal value

                if ( not found && csrJA[end] == i ) 
                {
                    found = true;
                    diagonalValue = csrValues[end];
                }

                // move up elements to fill the gap of diagonal element
                if ( found )
                {
                    csrJA[end] = csrJA[end - 1];
                    csrValues[end] = csrValues[end - 1];
                }

                end--;
            }

            if ( found )
            {
                // now set the first row element as the diagonal element
                csrValues[start] = diagonalValue;
                csrJA[start] = i;
                threadNumDiagonals++;
            }
        }

        atomicAdd( numFoundDiagonals, threadNumDiagonals );
    }

    return numFoundDiagonals;
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::getColumnPositions( 
    IndexType row[], 
    IndexType pos[],
    const IndexType j,
    const IndexType csrIA[], const IndexType numRows,
    const IndexType csrJA[], const IndexType )
{
    SCAI_REGION( "OpenMP.CSRUtils.getColPos" )

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
        // loop over all none zero elements of row i

        for ( IndexType jj = rIA[i]; jj < rIA[i + 1]; ++jj )
        {
            IndexType j = rJA[jj];   // inc size counter for this column
            atomicInc( cscIA[j] );   // old value not needed here
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

template<typename ValueType>
static void sumColumns(
    ValueType result[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    const IndexType numRows,
    const common::UnaryOp elemOp )
{
    #pragma omp parallel for
    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType jj = csrIA[i]; jj < csrIA[i+1]; ++jj )
        {
            IndexType j = csrJA[jj];
            ValueType v = csrValues[jj];
            v = applyUnary( elemOp, v );
            atomicAdd( result[j], v );
        }
   }
}

template<typename ValueType>
static void reduceColumns(
    ValueType result[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    const IndexType numRows,
    const common::BinaryOp reduceOp,
    const common::UnaryOp elemOp )
{
    // no parallel execution possible due to output dependencies on result

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType jj = csrIA[i]; jj < csrIA[i+1]; ++jj )
        {
            IndexType j = csrJA[jj];
            ValueType v = csrValues[jj];
            v = applyUnary( elemOp, v );
            result[j] = applyBinary( result[j], reduceOp, v );
        }
    }
}

template<typename ValueType>
void OpenMPCSRUtils::reduce(
    ValueType result[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    const IndexType numRows,
    const IndexType dim,
    const common::BinaryOp reduceOp,
    const common::UnaryOp elemOp )
{
    SCAI_REGION( "OpenMP.CSRUtils.reduce" )

    SCAI_LOG_INFO( logger, "reduce CSR[#rows = " << numRows << "], dim = 1, red = " << reduceOp << ", elem = " << elemOp )

    if ( dim == 0 )
    {
        #pragma omp parallel for
        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType jj = csrIA[i]; jj < csrIA[i+1]; ++jj )
            {
                ValueType v = csrValues[jj];
                v = applyUnary( elemOp, v );
                result[i] = applyBinary( result[i], reduceOp, v );
            }
        }
    }
    else if ( dim == 1 )
    {

        // parallel execution only if atomic updates are possible

        switch ( reduceOp )
        {
            case common::BinaryOp::ADD :

                sumColumns( result, csrIA, csrJA, csrValues, numRows, elemOp );
                break;

            default:

                reduceColumns( result, csrIA, csrJA, csrValues, numRows, reduceOp, elemOp );
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "illegal reduce dim = " << dim )
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
    const IndexType numColumns,
    const IndexType numValues,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[], 
    const common::MatrixOp op )
{
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "normalGEVM<" << TypeTraits<ValueType>::id() << ", launch it as an asynchronous task" )

        syncToken->run( std::bind( normalGEMV<ValueType>, 
                                   result, alpha, x, beta, y,
                                   numRows, numColumns, numValues,
                                   csrIA, csrJA, csrValues, op ) );
        return;
    }

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads()
                   << ">, result[" << numRows << "] = " << alpha << " * A * x + " << beta << " * y, matrix op = " << op )

    if ( op == common::MatrixOp::TRANSPOSE )
    {
        // result := alpha * x * A + beta * y -> result:= beta * y; result += alpha * x * A

        utilskernel::OpenMPUtils::binaryOpScalar( result, y, beta, numColumns, common::BinaryOp::MULT, false );

        #pragma omp parallel
        {
            // Note: region will be entered by each thread

            SCAI_REGION( "OpenMP.CSR.normalGEMV_t" )

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
    else if ( op == common::MatrixOp::NORMAL )
    {
        #pragma omp parallel
        {
            // Note: region will be entered by each thread
            SCAI_REGION( "OpenMP.CSR.normalGEMV_n" )
            #pragma omp for 

            for ( IndexType i = 0; i < numRows; ++i )
            {
                ValueType temp = 0;

                for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
                {
                    IndexType j = csrJA[jj];
                    temp += csrValues[jj] * x[j];
                }

                if ( y == NULL) 
                {
                    result[i] = alpha * temp;
                }
                else
                {
                    result[i] = alpha * temp + beta * y[i];
                }
            }
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "unsupported matrix operation: " << op );
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
    const common::MatrixOp op )
{
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        syncToken->run( std::bind( sparseGEMV<ValueType>, result, alpha,
                                   x, numNonZeroRows, rowIndexes, csrIA, csrJA, csrValues, op ) );
        return;
    }

    if ( op == common::MatrixOp::TRANSPOSE )
    {
        #pragma omp parallel
        {
            // Note: region will be entered by each thread

            SCAI_REGION( "OpenMP.CSR.sparseGEMV_t" )

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
    else if ( op == common::MatrixOp::NORMAL )
    {
        #pragma omp parallel
        {
            // Note: region will be entered by each thread
            SCAI_REGION( "OpenMP.CSR.sparseGEMV_n" )
            #pragma omp for 
    
            for ( IndexType ii = 0; ii < numNonZeroRows; ++ii )
            {
                ValueType temp = 0;
    
                IndexType i = rowIndexes[ii];
    
                for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
                {
                    IndexType j = csrJA[jj];
                    temp += csrValues[jj] * x[j];
                }

                result[i] += alpha * temp;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::gemmSD(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType nv,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    const common::MatrixOp op )
{
    SCAI_LOG_INFO( logger, "gemmSD<" << TypeTraits<ValueType>::id() << ">, " 
                           << " result " << numRows << " x " << nv << " CSR " << numRows << " x " << numColumns )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( op != common::MatrixOp::NORMAL )
    {
        COMMON_THROWEXCEPTION( "gemmSD only supported for matrix op = NORMAL" )
    }

    if ( syncToken )
    {
        SCAI_LOG_ERROR( logger, "asynchronous execution not supported here" )
    }

    #pragma omp parallel 
    {
        SCAI_REGION( "OpenMP.CSR.gemmSD" )
    
        #pragma omp for
        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType k = 0; k < nv; ++k )
            {
                ValueType temp = 0;
    
                for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
                {
                    IndexType j = csrJA[jj];
                    // SCAI_ASSERT_DEBUG( j < p , "index j = " << j << " out of range " << p )
                    // csrValues[jj] stands for CSR( i, j )
                    temp += csrValues[jj] * x[j * nv + k]; // x(j,k)
                }

                if ( y == NULL )
                {
                    result[i * nv + k] = alpha * temp;
                }
                else
                {
                    result[i * nv + k] = alpha * temp + beta * y[i * nv + k];
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::gemmDS(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType nv,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[],
    const common::MatrixOp op )
{
    SCAI_LOG_INFO( logger, "gemmDS<" << TypeTraits<ValueType>::id() << ">, " 
                           << " result " << nv << " x " << numColumns << " = " << alpha << " * "
                           << " Dense " << nv << " x " << numRows << ", CSR " << numRows << " x " << numColumns )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        syncToken->run( std::bind( gemmDS<ValueType>, result, alpha, x, beta, numRows, numColumns, nv,
                                                      csrIA, csrJA, csrValues, op ) );
        return;
    }

    if ( op != common::MatrixOp::NORMAL )
    {
        COMMON_THROWEXCEPTION( "gemmDS only supported for matrix op = NORMAL" )
    }

    if ( syncToken )
    {
        SCAI_LOG_ERROR( logger, "asynchronous execution not supported here" )
    }

    // matrix multiplicaton scatters results, so result must be set before

    utilskernel::OpenMPUtils::binaryOpScalar( result, result, beta, nv * numColumns, common::BinaryOp::MULT, false );

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.CSR.gemmDS" )

        #pragma omp for
        for ( IndexType k = 0; k < nv; ++k )
        {
            for ( IndexType i = 0; i < numRows; ++i )
            {
                for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
                {
                    IndexType j = csrJA[jj];
                    ValueType v = alpha * csrValues[jj] * x[k * numRows + i];
                    result[k * numColumns + j] += v;
                }
            }
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
        syncToken->run( std::bind( jacobi<ValueType>, solution, csrIA, csrJA, csrValues, oldSolution, rhs, omega, numRows ) );
        return;
    }

    const ValueType ZERO = 0;
    const ValueType ONE = 1;

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.CSR.jacobi" )
        #pragma omp for 

        for ( IndexType i = 0; i < numRows; i++ )
        {
            ValueType temp = rhs[i];

            ValueType diag = ZERO;

            for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; jj++ )
            {
                IndexType j = csrJA[jj];

                if ( j == i )
                {
                    diag = csrValues[jj];
                    SCAI_ASSERT_NE_DEBUG( diag, ZERO, "Diagonal element for row " << i << " is zero" )
                }
                else
                {
                    temp -= csrValues[jj] * oldSolution[j];
                }
            }

            SCAI_ASSERT_NE_DEBUG( diag, ZERO, "Diagonal element for row " << i << " not available" )

            // here we take advantange of a good branch precondiction

            if ( omega == ONE )
            {
                solution[i] = temp / diag;
            }
            else
            {
                solution[i] = omega * ( temp / diag ) + ( ONE - omega ) * oldSolution[i];
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::jacobiHalo(
    ValueType solution[],
    const ValueType localDiagonal[],
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
        #pragma omp for 

        for ( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            IndexType i = ii; // default: rowIndexes == NULL stands for identity

            if ( haloRowIndexes )
            {
                i = haloRowIndexes[ii];
            }

            ValueType temp = static_cast<ValueType>( 0.0 );

            const ValueType diag = localDiagonal[i];

            for ( IndexType j = haloIA[i]; j < haloIA[i + 1]; j++ )
            {
                temp += haloValues[j] * oldSolution[haloJA[j]];
            }

            solution[i] -= omega * ( temp / diag );
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

    std::unique_ptr<ValueType[]> denseA( new ValueType[ numRows * numRows ] );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numRows; ++j )
        {
            denseA[ i * numRows + j ] = static_cast<ValueType>( 0 );
        }
    }

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
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

    blaskernel::OpenMPBLAS2::gemv( common::MatrixOp::NORMAL,
                                   numRows, numRows, alpha, denseA.get(), numRows, rhs, inc1, beta, solution, inc1 );

}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::matrixAddSizes(
    IndexType cSizes[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType aIA[],
    const IndexType aJA[],
    const IndexType bIA[],
    const IndexType bJA[] )
{
    SCAI_LOG_INFO( logger,
                   "matrixAddSizes for " << numRows << " x " << numColumns << " matrix" )

    SCAI_REGION( "OpenMP.CSR.matrixAddSizes" )

    // determine the number of entries in output matrix

    #pragma omp parallel
    {
        BuildSparseIndexes indexList( numColumns );  // one builder for each thread

        #pragma omp for

        for ( IndexType i = 0; i < numRows; ++i )
        {
            // loop over all none zero elements of row i of input matrix a

            for ( IndexType jj = aIA[i]; jj < aIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix a
                // so we are at position a(i,j)
                IndexType j = aJA[jj];

                // element a(i,j) will generate an output element c(i,j)

                indexList.pushIndex( j );
            }

            for ( IndexType jj = bIA[i]; jj < bIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix b
                // so we are at position b(i,j)
                IndexType j = bJA[jj];

                // element a(i,j) will generate an output element c(i,j)

                indexList.pushIndex( j );
            }

            // so we have now the correct length

            cSizes[i] = indexList.getLength();

            SCAI_LOG_TRACE( logger, "row " << i << " will have " << cSizes[i] << " entries" )

            // reset indexList for next use

            indexList.reset();

        } // end loop over all rows of input matrices
    }

    IndexType numValues = sizes2offsets( cSizes, numRows );

    return numValues;
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::binaryOpSizes(
    IndexType cSizes[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType aIA[],
    const IndexType aJA[],
    const IndexType bIA[],
    const IndexType bJA[] )
{
    SCAI_LOG_INFO( logger,
                   "binaryOpSizes for " << numRows << " x " << numColumns << " storages" )

    // determine the number of entries in output matrix

    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType na = aIA[i + 1] - aIA[i];
        IndexType nb = bIA[i + 1] - bIA[i];

        // call routine for sparse vectors, column indexes are sorted for both rows

        cSizes[i] = utilskernel::OpenMPUtils::countAddSparse( &aJA[aIA[i]], na, &bJA[bIA[i]], nb );
    } 

    IndexType numValues = sizes2offsets( cSizes, numRows );

    SCAI_LOG_DEBUG( logger, "result storage: nnz = " << numValues 
                             << ", a:nnz = " << aIA[numRows] << ", b:nnz = " << bIA[numRows] )
    return numValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCSRUtils::binaryOp(
    IndexType cJA[],
    ValueType cValues[],
    const IndexType cIA[],
    const IndexType numRows,
    const IndexType,
    const IndexType aIA[],
    const IndexType aJA[],
    const ValueType aValues[],
    const IndexType bIA[],
    const IndexType bJA[],
    const ValueType bValues[],
    const common::BinaryOp op )
{
    ValueType zero = 0;

    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType na = aIA[i + 1] - aIA[i];
        IndexType nb = bIA[i + 1] - bIA[i];
        IndexType nc = cIA[i + 1] - cIA[i];

        // entries for each row are sorted, so call operation for sparse arrays

        IndexType n = utilskernel::OpenMPUtils::binopSparse( &cJA[cIA[i]], &cValues[cIA[i]],
                                                             &aJA[aIA[i]], &aValues[aIA[i]], zero, na, 
                                                             &bJA[bIA[i]], &bValues[bIA[i]], zero, nb, op );

        SCAI_ASSERT_EQ_ERROR( n, nc, "serious mismatch" )
    } 
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPCSRUtils::matrixMultiplySizes(
    IndexType cSizes[],
    const IndexType m,
    const IndexType n,
    const IndexType /* k */,
    const IndexType aIA[],
    const IndexType aJA[],
    const IndexType bIA[],
    const IndexType bJA[] )
{
    SCAI_REGION( "OpenMP.CSR.matrixMultiplySizes" )
    SCAI_LOG_INFO( logger,
                   "matrixMutliplySizes for " << m << " x " << n << " matrix" )

    // determine the number of entries in output matrix

    #pragma omp parallel
    {
        BuildSparseIndexes indexList( n );  // Each thread has its own build object, is reused

        #pragma omp for

        for ( IndexType i = 0; i < m; ++i )
        {
            indexList.reset();   // build sparse indexes for each row

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

                    indexList.pushIndex( k );
                }
            }

            // so we have now the correct length

            cSizes[i] = indexList.getLength();

            SCAI_LOG_TRACE( logger, "row " << i << " will have " << cSizes[i] << " entries" )

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
    const ValueType alpha,
    const IndexType aIA[],
    const IndexType aJA[],
    const ValueType aValues[],
    const ValueType beta,
    const IndexType bIA[],
    const IndexType bJA[],
    const ValueType bValues[] )
{
    // very close to binaryOp wito op = BinaryOp::ADD, but here also scalar mutliplication with alpha, beta

    SCAI_REGION( "OpenMP.CSR.matrixAdd" )

    SCAI_LOG_INFO( logger,
                   "matrixAdd for " << numRows << " x " << numColumns << " matrix" )

    // determine the number of entries in output matrix

    #pragma omp parallel
    {
        BuildSparseVector<ValueType> sparseRow( numColumns );

        #pragma omp for

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType jj = aIA[i]; jj < aIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix a
                // so we are at position a(i,j)

                IndexType j = aJA[jj];
                ValueType v = alpha * aValues[jj];

                sparseRow.push( j, v, common::BinaryOp::ADD );

                SCAI_LOG_TRACE( logger, "entry for [" << i << "," << j << "] by a, val = " << v )
            }

            for ( IndexType jj = bIA[i]; jj < bIA[i + 1]; ++jj )
            {
                // j is column of none zero element jj of row i of input matrix b
                // so we are at position b(i,j)
                IndexType j = bJA[jj];
                ValueType v = beta * bValues[jj];

                sparseRow.push( j, v, common::BinaryOp::ADD );

                SCAI_LOG_TRACE( logger, "entry for [" << i << "," << j << "] by b, val = " << v  )
            }

            IndexType offset = cIA[i];

            SCAI_LOG_DEBUG( logger, "fill row " << i << ", has " << sparseRow.getLength() << " entries, offset = " << offset )

            // fill in csrJA, csrValues and reset indexList, valueList for next use

            while ( !sparseRow.isEmpty() )
            {
                IndexType col;
                ValueType val;

                sparseRow.pop( col, val );

                SCAI_LOG_TRACE( logger, "row " << i << " has entry at col " << col << ", val = " << val << ", offset = " << offset )

                cJA[offset] = col;
                cValues[offset] = val;
                ++offset;
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
    const IndexType aIA[],
    const IndexType aJA[],
    const ValueType aValues[],
    const IndexType bIA[],
    const IndexType bJA[],
    const ValueType bValues[] )
{
    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.CSR.matrixMultiply" )

        BuildSparseVector<ValueType> sparseRow( n, ValueType( 0 ) ); // one for each thread

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

                    sparseRow.push( k, aValues[jj] * bValues[kk], common::BinaryOp::ADD );
                }
            }

            IndexType offset = cIA[i];

            SCAI_LOG_DEBUG( logger, "fill row " << i << ", has " << sparseRow.getLength() << " entries, offset = " << offset )

            // fill in csrJA, csrValues and reset indexList, valueList for next use

            while ( !sparseRow.isEmpty() )
            {
                IndexType col;
                ValueType val;

                sparseRow.pop( col, val );

                SCAI_LOG_TRACE( logger, "row " << i << " has entry at col " << col << ", val = " << val << ", offset = " << offset )

                cJA[offset] = col;
                cValues[offset] = alpha * val;
                ++offset;
            }

            // make sure that we have still the right offsets

            SCAI_ASSERT_EQUAL_DEBUG( offset, cIA[i + 1] )
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
RealType<ValueType> OpenMPCSRUtils::absMaxDiffRowUnsorted(
    const IndexType n1,
    const IndexType csrJA1[],
    const ValueType csrValues1[],
    const IndexType n2,
    const IndexType csrJA2[],
    const ValueType csrValues2[] )
{
    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    // No assumption about any sorting in a row

    RealType val = 0;

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

        RealType absDiff = common::Math::abs( diff );

        if ( absDiff > val )
        {
            val = absDiff;
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

        RealType absDiff = common::Math::abs( csrValues2[i2] );

        if ( absDiff > val )
        {
            val = absDiff;
        }
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> OpenMPCSRUtils::absMaxDiffRowSorted(
    const IndexType n1,
    const IndexType csrJA1[],
    const ValueType csrValues1[],
    const IndexType n2,
    const IndexType csrJA2[],
    const ValueType csrValues2[] )
{
    // Note: the implementation assumes that rows are sorted according to column indexes

    RealType<ValueType> val = 0;

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

        auto absDiff = common::Math::abs( diff );

        if ( absDiff > val )
        {
            val = absDiff;
        }
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> OpenMPCSRUtils::absMaxDiffVal(
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
    RealType<ValueType> ( *absMaxDiffRow )(
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

    RealType<ValueType> val = 0;

    #pragma omp parallel
    {
        RealType<ValueType> threadVal = 0;

        #pragma omp for 

        for ( IndexType i = 0; i < numRows; ++i )
        {
            IndexType offs1 = csrIA1[i];
            IndexType offs2 = csrIA2[i];
            IndexType n1 = csrIA1[i + 1] - offs1;
            IndexType n2 = csrIA2[i + 1] - offs2;

            RealType<ValueType> maxRow = absMaxDiffRow( n1, &csrJA1[offs1], &csrValues1[offs1], n2, &csrJA2[offs2], &csrValues2[offs2] );

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
    common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register CSRUtils OpenMP-routines for Host at kernel registry [" << flag << "]" )
    KernelRegistry::set<CSRKernelTrait::getValuePos>( getValuePos, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::getColumnPositions>( getColumnPositions, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::sizes2offsets>( sizes2offsets, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::offsets2sizes>( offsets2sizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::gatherSizes>( gatherSizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::validOffsets>( validOffsets, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::nonEmptyRows>( nonEmptyRows, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::hasDiagonalProperty>( hasDiagonalProperty, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::hasSortedRows>( hasSortedRows, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::getPosDiagonal>( getPosDiagonal, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixAddSizes>( matrixAddSizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::binaryOpSizes>( binaryOpSizes, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiplySizes>( matrixMultiplySizes, ctx, flag );
}

template<typename ValueType>
void OpenMPCSRUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register CSRUtils OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<CSRKernelTrait::convertCSR2CSC<ValueType> >( convertCSR2CSC, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::sortRows<ValueType> >( sortRows, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::shiftDiagonal<ValueType> >( shiftDiagonal, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::getDiagonal<ValueType> >( getDiagonal, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::setDiagonal<ValueType> >( setDiagonal, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::setDiagonalV<ValueType> >( setDiagonalV, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::reduce<ValueType> >( reduce, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::sparseGEMV<ValueType> >( sparseGEMV, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::gemmSD<ValueType> >( gemmSD, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::gemmDS<ValueType> >( gemmDS, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixAdd<ValueType> >( matrixAdd, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::binaryOp<ValueType> >( binaryOp, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::matrixMultiply<ValueType> >( matrixMultiply, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::absMaxDiffVal<ValueType> >( absMaxDiffVal, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::countNonZeros<ValueType> >( countNonZeros, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::compress<ValueType> >( compress, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::decomposition<ValueType> >( decomposition, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::setRows<ValueType> >( setRows, ctx, flag );
    KernelRegistry::set<CSRKernelTrait::setColumns<ValueType> >( setColumns, ctx, flag );
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
}

OpenMPCSRUtils::~OpenMPCSRUtils()
{
    SCAI_LOG_INFO( logger, "unregister CSRUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPCSRUtils OpenMPCSRUtils::guard;

} /* end namespace sparsekernel */

} /* end namespace scai */
