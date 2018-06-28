/**
 * @file OpenMPELLUtils.cpp
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
 * @brief Implementation of ELL utilities with OpenMP
 * @author Thomas Brandes
 * @date 04.07.2012
 */

// hpp
#include <scai/sparsekernel/openmp/OpenMPELLUtils.hpp>

// local library
#include <scai/sparsekernel/ELLKernelTrait.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>


#include <scai/common/macros/assert.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

// std
#include <cmath>
#include <set>
#include <map>
#include <memory>
#include <functional>

namespace scai
{

using common::TypeTraits;
using tasking::TaskSyncToken;

namespace sparsekernel
{

/* ------------------------------------------------------------------------------------------------------------------ */

SCAI_LOG_DEF_LOGGER( OpenMPELLUtils::logger, "OpenMP.ELLUtils" )

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::setRows(
    ValueType ellValues[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const ValueType values[],
    const common::BinaryOp op )
{
    SCAI_LOG_INFO( logger,
                   "setRows<" << TypeTraits<ValueType>::id() << ">" << ", #numRows = " << numRows )
    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; i++ ) //rows
    {
        for ( IndexType jj = 0; jj < ellSizes[i]; jj++ ) //elements in row
        {
            IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
            ellValues[pos] = common::applyBinary( ellValues[pos], op, values[i] );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

void OpenMPELLUtils::check(
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType numColumns,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const char* msg )
{
    SCAI_LOG_INFO( logger,
                   "check # numRows = " << numRows << ", numValuesPerRow = " << numValuesPerRow << ", numColumns = " << numColumns )

    if ( numRows > 0 )
    {
        bool integrityIA = true;
        bool integrityJA = true;
        #pragma omp parallel for reduction( && : integrityIA, integrityJA ) 

        for ( IndexType i = 0; i < numRows; i++ )
        {
            // Attention, ellSizes[i] == numValuesPerRow is also valid
            if ( common::Utils::validIndex( ellSizes[i], numValuesPerRow + 1 ) )
            {
                for ( IndexType jj = 0; jj < ellSizes[i]; jj++ )
                {
                    IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                    IndexType j = ellJA[pos];
                    integrityJA = integrityJA && common::Utils::validIndex( j, numColumns );
                }
            }
            else
            {
                integrityIA = false;
            }
        }

        SCAI_ASSERT_ERROR( integrityIA, msg << ": ellSizes: at least one value out of range" )
        SCAI_ASSERT_ERROR( integrityJA, msg << ": ellJA: at least one value out of range" )
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
ValueType OpenMPELLUtils::absMaxVal(
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const ValueType values[] )
{
    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    RealType maxValue = 0;

    #pragma omp parallel
    {
        RealType threadVal = 0;

        #pragma omp for 

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );

                RealType val = common::Math::abs( values[pos] );

                if ( val > threadVal )
                {
                    threadVal = val;
                }

                // SCAI_LOG_TRACE( logger, "absMaxVal, val[" << i << ", " << jj << "] = " << val )
            }
        }

        #pragma omp critical
        {
            SCAI_LOG_DEBUG( logger, "absMaxVal, threadVal = " << threadVal << ", maxVal = " << maxValue )

            if ( threadVal > maxValue )
            {
                maxValue = threadVal;
            }
        }
    }

    SCAI_LOG_DEBUG( logger, "absMaxVal, maxVal = " << maxValue )

    return maxValue;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::getRow(
    ValueType row[],
    const IndexType i,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType values[] )
{
    SCAI_LOG_DEBUG( logger, "get row #i = " << i )
    #pragma omp parallel for 

    for ( IndexType j = 0; j < numColumns; ++j )
    {
        row[j] = ValueType( 0 );
    }

    #pragma omp parallel for 

    for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
    {
        IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
        row[ellJA[pos]] = values[pos];
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

IndexType OpenMPELLUtils::getValuePos(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[] )
{
    IndexType vPos = invalidIndex;

    for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
    {
        IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );

        if ( ellJA[pos] == j )
        {
            vPos = pos;
            break;
        }
    }

    return vPos;
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPELLUtils::getDiagonalPositions(
    IndexType diagonalPositions[],
    const IndexType numDiagonals,
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[] )
{
    IndexType count = 0;    // counts the diagonal elements really found

    #pragma omp parallel for reduction( + : count )
    for ( IndexType i = 0; i < numDiagonals; ++i )
    {
        diagonalPositions[i] = invalidIndex;    // will be set correctly if entry found

        for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
        {
            IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );

            if ( ellJA[pos] == i )
            {
                diagonalPositions[i] = pos;
                count++;  
                break;
            }
        }
    }

    return count;
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPELLUtils::getColumnPositions( 
    IndexType row[], 
    IndexType pos[],
    const IndexType j,
    const IndexType ellIA[], 
    const IndexType numRows,
    const IndexType ellJA[], 
    const IndexType numValuesPerRow )
{
    SCAI_REGION( "OpenMP.ELLUtils.getColumnPositions" )

    IndexType cnt  = 0;   // counts number of available row entries in column j

    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType jj = 0; jj < ellIA[i]; ++jj )
        {
            IndexType p = ellindex( i, jj, numRows, numValuesPerRow );

            if ( ellJA[p] == j )
            {
                IndexType k = atomicInc( cnt );
                row[k] = i;
                pos[k] = p;
                break;
            }
        }
    }

    return cnt;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::getCSRValues(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[] )
{
    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows )
    // parallelization possible as offset array csrIA is available
    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.ELL.getCSR" )
        #pragma omp for 

        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType rowSize = ellSizes[i];
            IndexType offset = csrIA[i];
            // just make sure that csrIA and ellSizes really fit with each other
            SCAI_ASSERT_EQUAL_DEBUG( csrIA[i] + rowSize, csrIA[i + 1] )

            for ( IndexType jj = 0; jj < rowSize; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                csrJA[offset + jj] = ellJA[pos];
                csrValues[offset + jj] = ellValues[pos];
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::setCSRValues(
    IndexType ellJA[],
    ValueType ellValues[],
    const IndexType ellSizes[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    SCAI_LOG_INFO( logger,
                   "set CSRValues<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " 
                   << numRows << ", #values/row = " << numValuesPerRow )

    // parallelization possible as offset array csrIA is available

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.ELL.setCSR" )

        #pragma omp for 

        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType rowSize = ellSizes[i];
            IndexType offset = csrIA[i];
            IndexType j = 0; // will be last column index

            for ( IndexType jj = 0; jj < rowSize; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                j = csrJA[offset + jj];
                ellJA[pos] = j;
                ellValues[pos] = csrValues[offset + jj];
            }

            // fill up the remaining entries with something useful

            for ( IndexType jj = rowSize; jj < numValuesPerRow; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                ellJA[pos] = j; // last used column index
                ellValues[pos] = ValueType( 0 ); 
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::fillELLValues(
    IndexType ellJA[],
    ValueType ellValues[],
    const IndexType ellSizes[],
    const IndexType numRows,
    const IndexType numValuesPerRow )
{
    SCAI_LOG_INFO( logger, "fill ELLValues<" << TypeTraits<ValueType>::id() << ">" )

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.ELL.fillELLValues" )

        #pragma omp for 

        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType rowSize = ellSizes[i];
            IndexType j = 0; // will be last column index

            if ( rowSize > 0 && rowSize < numValuesPerRow )
            {
                IndexType pos = ellindex( i, rowSize - 1, numRows, numValuesPerRow );
                j = ellJA[pos];
            }

            // fill up the remaining entries with something useful

            for ( IndexType jj = rowSize; jj < numValuesPerRow; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                ellJA[pos] = j; // last used column index
                ellValues[pos] = static_cast<ValueType>( 0 ); // zero entry
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::compressIA(
    IndexType newIA[],
    const IndexType ellIA[],
    const IndexType[],
    const ValueType ellValues[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const RealType<ValueType> eps )
{
    SCAI_LOG_INFO( logger, "compressIA with eps = " << eps )
    #pragma omp parallel
    {
        #pragma omp for

        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType length = ellIA[i];

            for ( IndexType jj = 0; jj < ellIA[i]; jj++ )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );

                if ( common::Math::abs( ellValues[pos] ) <= eps )
                {
                    length--;    // entry is deleted
                }

            }

            newIA[i] = length;
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::compressValues(
    IndexType newJA[],
    ValueType newValues[],
    const IndexType newNumValuesPerRow,
    const IndexType IA[],
    const IndexType JA[],
    const ValueType values[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const RealType<ValueType> eps )
{
    SCAI_LOG_INFO( logger, "compressValues ( #rows = " << numRows
                   << ", values/row = " << numValuesPerRow << " / " << newNumValuesPerRow
                   << ") with eps = " << eps )

    #pragma omp parallel
    {
        #pragma omp for

        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType gap = 0;

            for ( IndexType jj = 0; jj < IA[i]; jj++ )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );

                // delete it if zero 
 
                if ( common::Math::abs( values[pos] ) > eps )
                {
                    // move entry gap positions back in this row
    
                    IndexType newpos = ellindex( i, jj - gap, numRows, newNumValuesPerRow );
                    newValues[newpos] = values[pos];
                    newJA[newpos] = JA[pos];
                }
                else
                {
                    gap++;
                }
            }

            // fill up to top

            for (  IndexType jj = IA[i] - gap; jj < newNumValuesPerRow; jj++ )
            {
                IndexType newpos = ellindex( i, jj, numRows, newNumValuesPerRow );
                newValues[newpos] = 0;
                newJA[newpos] = 0;
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

void OpenMPELLUtils::matrixMultiplySizes(
    IndexType cSizes[],
    const IndexType aNumRows,
    const IndexType SCAI_UNUSED( aNumColumns ),
    const IndexType bNumRows,
    const bool SCAI_UNUSED( diagonalProperty ),
    const IndexType aSizes[],
    const IndexType aJA[],
    const IndexType aNumValuesPerRow,
    const IndexType bSizes[],
    const IndexType bJA[],
    const IndexType bNumValuesPerRow )
{
    SCAI_LOG_INFO( logger, "matrixMultiplySizes with numRows A = " << aNumRows << " and numRows B = " << bNumRows )
    #pragma omp parallel
    {
        #pragma omp for

        for ( IndexType i = 0; i < aNumRows; i++ )
        {
            std::set<IndexType> newElements;
            std::pair<std::set<IndexType>::iterator, bool> ret;
            IndexType length = 0;

            for ( IndexType j = 0; j < aSizes[i]; j++ )
            {
                IndexType posA = ellindex( i, j, aNumRows, aNumValuesPerRow );
                IndexType jj = aJA[posA];

                for ( IndexType k = 0; k < bSizes[jj]; k++ )
                {
                    IndexType posB = ellindex( jj, k, bNumRows, bNumValuesPerRow );
                    IndexType kk = bJA[posB];
                    ret = newElements.insert( kk );

                    if ( ret.second == true )
                    {
                        length++;
                    }
                }
            }

            cSizes[i] = length;
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::matrixMultiply(
    IndexType cJA[],
    ValueType cValues[],
    const IndexType cSizes[],
    const IndexType cNumValuesPerRow,
    const IndexType aNumRows,
    const IndexType SCAI_UNUSED( aNumColumns ),
    const IndexType bNumRows,
    const bool SCAI_UNUSED( diagonalProperty ),
    const ValueType alpha,
    const IndexType aSizes[],
    const IndexType aJA[],
    const ValueType aValues[],
    const IndexType aNumValuesPerRow,
    const IndexType bSizes[],
    const IndexType bJA[],
    const ValueType bValues[],
    const IndexType bNumValuesPerRow )
{
    SCAI_LOG_INFO( logger, "matrix multiply with numRows A = " << aNumRows << " and numRows B = " << bNumRows )
    #pragma omp parallel
    {
        #pragma omp for

        for ( IndexType i = 0; i < aNumRows; i++ )
        {
            std::set<IndexType> jaRow;
            std::map<IndexType, ValueType> valuesRow;
            std::pair<std::set<IndexType>::iterator, bool> ret;

            for ( IndexType j = 0; j < aSizes[i]; j++ )
            {
                IndexType posA = ellindex( i, j, aNumRows, aNumValuesPerRow );
                IndexType jj = aJA[posA];

                for ( IndexType k = 0; k < bSizes[jj]; k++ )
                {
                    IndexType posB = ellindex( jj, k, bNumRows, bNumValuesPerRow );
                    IndexType kk = bJA[posB];
                    ret = jaRow.insert( kk );
                    ValueType mult = alpha * aValues[posA] * bValues[posB];

                    if ( ret.second == true )
                    {
                        valuesRow.insert( std::pair<IndexType, ValueType>( kk, mult ) );
                    }
                    else
                    {
                        valuesRow[kk] += mult;
                    }
                }
            }

            std::set<IndexType>::iterator jaIter;
            typename std::map<IndexType, ValueType>::iterator valuesIter;
            jaIter = jaRow.begin();
            valuesIter = valuesRow.begin();

            for ( IndexType j = 0; j < cSizes[i]; j++ )
            {
                // note: cNumRows == aNumRows
                IndexType posC = ellindex( i, j, aNumRows, cNumValuesPerRow );
                cJA[posC] = *jaIter;
                cValues[posC] = ( *valuesIter ).second;
                jaIter++;
                valuesIter++;
            }

            // fill up to top

            for ( IndexType j = cSizes[i]; j < cNumValuesPerRow; j++ )
            {
                IndexType posC = ellindex( i, j, aNumRows, cNumValuesPerRow );
                cJA[posC] = 0;
                cValues[posC] = 0;
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

void OpenMPELLUtils::matrixAddSizes(
    IndexType cSizes[],
    const IndexType m,
    const IndexType n,
    const bool SCAI_UNUSED( diagonalProperty ),
    const IndexType aSizes[],
    const IndexType aJA[],
    const IndexType aNumValuesPerRow,
    const IndexType bSizes[],
    const IndexType bJA[],
    const IndexType bNumValuesPerRow )
{
    SCAI_REGION( "OpenMP.ELL.matrixAddSizes" )

    SCAI_LOG_INFO( logger, "matrixAddSizes A + B, #rows = " << m )

    const IndexType NINIT = n + 1; // marks unused colums
    const IndexType END   = n + 2; // marks end of list

    // determine the number of entries in output matrix

    #pragma omp parallel
    {
        std::unique_ptr<IndexType[]> indexList( new IndexType[n] );

        for ( IndexType j = 0; j < n; j++ )
        {
            indexList[j] = NINIT;
        }

        #pragma omp for

        for ( IndexType i = 0; i < m; i++ )
        {
            IndexType length   = 0;
            IndexType firstCol = END;

            for ( IndexType jj = 0; jj < aSizes[i]; jj++ )
            {
                IndexType posA = ellindex( i, jj, m, aNumValuesPerRow );
                IndexType j    = aJA[posA];

                if ( indexList[j] == NINIT )
                {
                    // Add column position j to the indexList

                    indexList[j] = firstCol;
                    firstCol = j;
                    ++length;
                }
            }

            for ( IndexType jj = 0; jj < bSizes[i]; jj++ )
            {
                IndexType posB = ellindex( i, jj, m, bNumValuesPerRow );
                IndexType j    =  bJA[posB];

                if ( indexList[j] == NINIT )
                {
                    // Add column position j to the indexList

                    indexList[j] = firstCol;
                    firstCol = j;
                    ++length;
                }
            }

            cSizes[i] = length;

            // reset the indexList for next use

            while ( firstCol != END )
            {
                IndexType nextCol = indexList[firstCol];
                indexList[firstCol] = NINIT;
                firstCol = nextCol;
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::sortRowElements(
    IndexType ellJA[],
    ValueType ellValues[],
    const IndexType ellIA[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const bool diagonalFlag )
{
    SCAI_REGION( "OpenMP.ELL.sortRowElements" )

    SCAI_LOG_INFO( logger, "sort elements in each of " << numRows << " rows, diagonal flag = " << diagonalFlag )

    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        // use bubble sort as sort algorithm

        IndexType start = 0;
        IndexType end   = ellIA[i] - 1;;

        SCAI_LOG_DEBUG( logger, "row " << i << ": sort " << start << " - " << end )

        bool sorted = false;

        while ( !sorted )
        {
            sorted = true; // will be reset if any wrong order appears

            SCAI_LOG_TRACE( logger, "sort in row[" << i << "] from " << start << " - " << end )

            for ( IndexType jj = start; jj < end; ++jj )
            {
                IndexType pos  = ellindex( i, jj, numRows, numValuesPerRow );
                IndexType pos1 = ellindex( i, jj + 1, numRows, numValuesPerRow );

                bool swapIt = false;

                // if diagonalFlag is set, column i is the smallest one

                if ( diagonalFlag && ( ellJA[pos1] == i ) && ( ellJA[pos] != i ) )
                {
                    swapIt = true;
                }
                else if ( diagonalFlag && ( ellJA[pos] == i ) )
                {
                    swapIt = false;
                }
                else
                {
                    swapIt = ellJA[pos] > ellJA[pos1];
                }

                if ( swapIt )
                {
                    sorted = false;
                    std::swap( ellJA[pos], ellJA[pos1] );
                    std::swap( ellValues[pos], ellValues[pos1] );
                }
            }

            --end;
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::matrixAdd(
    IndexType cJA[],
    ValueType cValues[],
    const IndexType cSizes[],
    const IndexType cNumValuesPerRow,
    const IndexType m,
    const IndexType numColumns,
    const bool SCAI_UNUSED( diagonalProperty ),
    const ValueType alpha,
    const IndexType aSizes[],
    const IndexType aJA[],
    const ValueType aValues[],
    const IndexType aNumValuesPerRow,
    const ValueType beta,
    const IndexType bSizes[],
    const IndexType bJA[],
    const ValueType bValues[],
    const IndexType bNumValuesPerRow )
{
    SCAI_REGION( "OpenMP.ELL.matrixAdd" )

    SCAI_LOG_INFO( logger, "matrixAddELL C = " << alpha << " * A + " << beta << " * B, #rows = " << m )

    const IndexType NINIT = numColumns + 1;
    const IndexType END = numColumns + 2;

    #pragma omp parallel
    {

        std::unique_ptr<IndexType[]> indexList( new IndexType[numColumns] );
        std::unique_ptr<ValueType[]> valueList( new ValueType[numColumns] );

        for ( IndexType j = 0; j < numColumns; j++ )
        {
            indexList[j] = NINIT;
            valueList[j] = static_cast<ValueType>( 0 );
        }

        #pragma omp for

        for ( IndexType i = 0; i < m; i++ )
        {
            IndexType length = 0;
            IndexType firstCol = END;

            for ( IndexType jj = 0; jj < aSizes[i]; jj++ )
            {
                IndexType posA = ellindex( i, jj, m, aNumValuesPerRow );
                IndexType j = aJA[posA];
                valueList[j] += alpha * aValues[posA];

                // element a(i,j) will generate an output element c(i,j)

                if ( indexList[j] == NINIT )
                {
                    // Add column position j to the indexList
                    indexList[j] = firstCol;
                    firstCol = j;
                    ++length;
                }
            }

            for ( IndexType jj = 0; jj < bSizes[i]; jj++ )
            {
                IndexType posB = ellindex( i, jj, m, bNumValuesPerRow );

                IndexType j = bJA[posB];
                valueList[j] += beta * bValues[posB];

                // element b(i,j) will generate an output element c(i,j)

                if ( indexList[j] == NINIT )
                {
                    // Add column position j to the indexList
                    indexList[j] = firstCol;
                    firstCol = j;
                    ++length;
                }
            }

            // fill in cJA, cValues and reset indexList, valueList for next use

            IndexType jj = 0;

            while ( firstCol != END )
            {
                IndexType nextCol = indexList[firstCol];
                ValueType val = valueList[firstCol];
                indexList[firstCol] = NINIT;
                valueList[firstCol] = static_cast<ValueType>( 0 ); // reset for next time

                SCAI_LOG_DEBUG( logger, "entry for [" << i << "," << firstCol << "] = " << val )

                IndexType posC = ellindex( i, jj, m, cNumValuesPerRow );
                cJA[posC] = firstCol;
                cValues[posC] = val;

                firstCol = nextCol;
                ++jj;
            }

            SCAI_ASSERT_EQ_ERROR( jj, cSizes[i], "serious mismatch for cSizes[" << i << "]" )

            // fill up to top

            for (  IndexType jj = cSizes[i]; jj < cNumValuesPerRow; ++jj )
            {
                IndexType posC = ellindex( i, jj, m, cNumValuesPerRow );
                cJA[posC] = static_cast<IndexType>( 0 );
                cValues[posC] = static_cast<ValueType>( 0 );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::jacobi(
    ValueType solution[],
    const IndexType numRows,
    const IndexType ellNumValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega )
{
    SCAI_LOG_INFO( logger,
                   "jacobi<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken != NULL )
    {
        // run this method with same arguments by another thread

        syncToken->run( std::bind( jacobi<ValueType>,
                                   solution, numRows, ellNumValuesPerRow,
                                   ellSizes, ellJA, ellValues, oldSolution, rhs, omega ) );
        return;
    }

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.ELL.jacobi" )
        #pragma omp for 

        for ( IndexType i = 0; i < numRows; i++ )
        {
            ValueType temp = rhs[i];
            IndexType pos = ellindex( i, 0, numRows, ellNumValuesPerRow );
            ValueType diag = 0;

            for ( IndexType j = 0; j < ellSizes[i]; j++ )
            {
                pos = ellindex( i, j, numRows, ellNumValuesPerRow );

                if ( ellJA[pos] == i )
                {
                    diag = ellValues[pos];  // diagonal element
                }
                else
                {
                    temp -= ellValues[pos] * oldSolution[ellJA[pos]];
                }
            }

            if ( omega == scai::common::Constants::ONE )
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

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::jacobiHalo(
    ValueType solution[],
    const IndexType numRows,
    const ValueType diagonal[],
    const IndexType ellNumValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    const IndexType rowIndexes[],
    const IndexType numNonEmptyRows,
    const ValueType oldSolution[],
    const ValueType omega )
{
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken != NULL )
    {
        SCAI_LOG_WARN( logger, "jacobi called asynchronously, not supported here" )
    }

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.ELL.jacobiHalo" )
        #pragma omp for 

        for ( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            IndexType i = ii; // rowIndexes == NULL stands for all rows

            if ( rowIndexes )
            {
                i = rowIndexes[ii];
            }

            ValueType temp = static_cast<ValueType>( 0.0 );

            for ( IndexType jj = 0; jj < ellSizes[i]; jj++ )
            {
                IndexType pos = ellindex( i, jj, numRows, ellNumValuesPerRow );
                temp += ellValues[pos] * oldSolution[ellJA[pos]];
            }

            const ValueType diag = diagonal[i];

            solution[i] -= temp * ( omega / diag );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    const common::MatrixOp op )
{
    if ( numValuesPerRow == 0 )
    {
        COMMON_THROWEXCEPTION( "normalGEMV should not have been called, no entries" )
        // only compute: result = beta * y
    }

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        // combine the vectors with their scaling factors to reduze number of args
        syncToken->run( std::bind( normalGEMV<ValueType>, result,
                                   alpha, x, beta, y,
                                   numRows, numColumns, numValuesPerRow, 
                                   ellSizes, ellJA, ellValues, op ) );
        return;
    }

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads()
                   << ">, result[" << numRows << "] = " << alpha << " * A( ell, #maxNZ/row = " << numValuesPerRow 
                   << " ) * x + " << beta << " * y " )

    if ( op == common::MatrixOp::TRANSPOSE )
    {
        // result := alpha * x * A + beta * y -> result:= beta * y; result += alpha * x * A

        utilskernel::OpenMPUtils::binaryOpScalar( result, y, beta, numColumns, common::BinaryOp::MULT, false );

        #pragma  omp parallel
        {
            SCAI_REGION( "OpenMP.ELL.gemv_t" )

            #pragma omp for

            for ( IndexType i = 0; i < numRows; ++i )
            {
                for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
                {
                    IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                    IndexType j   = ellJA[pos];
                    ValueType v   = alpha * ellValues[pos] * x[i];
    
                    atomicAdd( result[j], v );
                }
            }
        }
    }
    else if ( op == common::MatrixOp::NORMAL )
    {
        #pragma omp parallel
        {
            SCAI_REGION( "OpenMP.ELL.gemv_n" )
            #pragma omp for 

            for ( IndexType i = 0; i < numRows; ++i )
            {
                ValueType temp = 0;

                for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
                {
                    IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                    IndexType j = ellJA[pos];
                    temp += ellValues[pos] * x[j];
                }
    
                if ( beta == scai::common::Constants::ZERO )
                {
                    // must be handled separately as y[i] might be uninitialized
                    result[i] = alpha * temp;
                }
                else
                {
                    result[i] = alpha * temp + beta * y[i];
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::sparseGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType numNonZeroRows,
    const IndexType rowIndexes[],
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    const common::MatrixOp op )
{
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        syncToken->run( std::bind( sparseGEMV<ValueType>, result, alpha, x,
                                   numRows, numValuesPerRow,
                                   numNonZeroRows, rowIndexes,
                                   ellSizes, ellJA, ellValues, op ) );
        return;
    }

    SCAI_LOG_INFO( logger,
                   "sparseGEMV<" << TypeTraits<ValueType>::id() << ">, n = " << numRows
                   << ", nonZeroRows = " << numNonZeroRows << ", alpha = " << alpha )

    if ( op == common::MatrixOp::TRANSPOSE )
    {
        #pragma omp parallel
        {
            SCAI_REGION( "OpenMP.ELL.spmv_t" )

            #pragma omp for

            for ( IndexType ii = 0; ii < numNonZeroRows; ++ii )
            {
                IndexType i = rowIndexes[ii];

                for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
                {
                    IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                    IndexType j   = ellJA[pos];
                    ValueType v   = alpha * ellValues[pos] * x[i];
    
                    atomicAdd( result[j], v );
                }
            }
        }
    }
    else if ( op == common::MatrixOp::NORMAL )
    {
        #pragma omp parallel
        {
            SCAI_REGION( "OpenMP.ELL.spmv_n" )

            #pragma omp for 

            for ( IndexType ii = 0; ii < numNonZeroRows; ++ii )
            {
                IndexType i = rowIndexes[ii];

                ValueType temp = 0;
    
                for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
                {
                    IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                    IndexType j = ellJA[pos];
                    temp += ellValues[pos] * x[j];
                }
    
                if ( alpha == scai::common::Constants::ONE )
                {
                    result[i] += temp;
                }
                else
                {
                    result[i] += alpha * temp;
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

void OpenMPELLUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register ELLtils OpenMP-routines for Host at kernel registry [" << flag << "]" )
    KernelRegistry::set<ELLKernelTrait::getValuePos>( getValuePos, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::getColumnPositions>( getColumnPositions, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::getDiagonalPositions>( getDiagonalPositions, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::check>( check, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::matrixMultiplySizes>( matrixMultiplySizes, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::matrixAddSizes>( matrixAddSizes, ctx, flag );
}

template<typename ValueType>
void OpenMPELLUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register ELLUtils OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<ELLKernelTrait::absMaxVal<ValueType> >( absMaxVal, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::compressIA<ValueType> >( compressIA, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::compressValues<ValueType> >( compressValues, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::sortRowElements<ValueType> >( sortRowElements, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::matrixAdd<ValueType> >( matrixAdd, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::matrixMultiply<ValueType> >( matrixMultiply, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::sparseGEMV<ValueType> >( sparseGEMV, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::getRow<ValueType> >( getRow, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::fillELLValues<ValueType> >( fillELLValues, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::setRows<ValueType> >( setRows, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::setCSRValues<ValueType> >( setCSRValues, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::getCSRValues<ValueType> >( getCSRValues, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPELLUtils::OpenMPELLUtils()
{
    SCAI_LOG_INFO( logger, "register ELLtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

OpenMPELLUtils::~OpenMPELLUtils()
{
    SCAI_LOG_INFO( logger, "unregister ELLtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPELLUtils OpenMPELLUtils::guard;

} /* end namespace sparsekernel */

} /* end namespace scai */
