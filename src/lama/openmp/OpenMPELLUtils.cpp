/**
 * @file OpenMPELLUtils.cpp
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
 * @brief Implementation of ELL utilities with OpenMP
 * @author Thomas Brandes
 * @date 04.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/openmp/OpenMPELLUtils.hpp>
#include <lama/openmp/OpenMP.hpp>

// others
#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

#include <lama/task/TaskSyncToken.hpp>

// macros
#include <lama/macros/unused.hpp>

// tracing
#include <lama/tracing.hpp>

// boost
#include <boost/bind.hpp>

// stl
#include <set>
#include <map>
#include <cmath>

namespace lama
{

/* ------------------------------------------------------------------------------------------------------------------ */

LAMA_LOG_DEF_LOGGER( OpenMPELLUtils::logger, "OpenMP.ELLUtils" )

/* ------------------------------------------------------------------------------------------------------------------ */

IndexType OpenMPELLUtils::countNonEmptyRowsBySizes( const IndexType sizes[], const IndexType numRows )
{
    IndexType counter = 0;

    #pragma omp parallel for reduction( +:counter )
    for ( IndexType i = 0; i < numRows; ++i )
    {
        if ( sizes[i] > 0 )
        {
            counter++;
        }
    }

    LAMA_LOG_INFO( logger, "#non-zero rows = " << counter << ", counted by sizes" )

    return counter;
}

/* ------------------------------------------------------------------------------------------------------------------ */

void OpenMPELLUtils::setNonEmptyRowsBySizes(
    IndexType rowIndexes[],
    const IndexType numNonEmptyRows,
    const IndexType sizes[],
    const IndexType numRows )
{
    IndexType counter = 0;

    // Note: this routine is not easy to parallelize, no offsets for rowIndexes available

    for ( IndexType i = 0; i < numRows; ++i )
    {
        if ( sizes[i] > 0 )
        {
            rowIndexes[counter] = i;
            counter++;
        }
    }

    LAMA_ASSERT_EQUAL_DEBUG( counter, numNonEmptyRows )

    LAMA_LOG_INFO( logger, "#non-zero rows = " << counter << ", set by sizes" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

bool OpenMPELLUtils::hasDiagonalProperty( const IndexType numDiagonals, const IndexType ellJA[] )
{
    LAMA_LOG_INFO( logger, "hasDiagonalProperty, #numDiagonals = " << numDiagonals )

    if ( numDiagonals == 0 )
    {
        return false;
    }

    bool diagonalProperty = true;

    #pragma omp parallel for reduction( && : diagonalProperty )
    for ( IndexType i = 0; i < numDiagonals; ++i )
    {
        if ( !diagonalProperty )
        {
            continue;
        }
        if ( ellJA[i] != i )
        {
            diagonalProperty = false;
        }
    }

    return diagonalProperty;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void OpenMPELLUtils::scaleValue(
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    ValueType ellValues[],
    const OtherValueType values[] )
{
    LAMA_LOG_INFO( logger,
                   "scaleValue, #numRows = " << numRows << ", ellSizes = " << ellSizes << ", ellValues = " << ellValues << ", values = " << values )

    #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
    for ( IndexType i = 0; i < numRows; i++ ) //rows
    {
        for ( IndexType jj = 0; jj < ellSizes[i]; jj++ ) //elements in row
        {
            IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
            ellValues[pos] *= static_cast<ValueType>( values[i] );
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
    LAMA_LOG_INFO( logger,
                   "check # numRows = " << numRows << ", numValuesPerRow = " << numValuesPerRow << ", numColumns = " << numColumns )

    if ( numRows > 0 )
    {
        bool integrityIA = true;
        bool integrityJA = true;

        #pragma omp parallel for reduction( && : integrityIA, integrityJA ) schedule( LAMA_OMP_SCHEDULE )
        for ( IndexType i = 0; i < numRows; i++ )
        {
            if ( ellSizes[i] >= 0 && ellSizes[i] <= numValuesPerRow )
            {
                for ( IndexType jj = 0; jj < ellSizes[i]; jj++ )
                {
                    IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                    IndexType j   = ellJA[pos];
                    integrityJA   = integrityJA && ( 0 <= j && j < numColumns );
                }
            }
            else
            { 
                integrityIA = false;
            }
        }

        LAMA_ASSERT_ERROR( integrityIA, msg << ": ellSizes: at least one value out of range" )
        LAMA_ASSERT_ERROR( integrityJA, msg << ": ellJA: at least one value out of range" )
    }
    else
    {
        LAMA_ASSERT_EQUAL_ERROR( 0, numValuesPerRow )
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
    ValueType maxValue = static_cast<ValueType>( 0.0 );

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>( 0.0 );

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                ValueType val = std::abs( values[pos] );

                if ( val > threadVal )
                {
                    threadVal = val;
                }

                // LAMA_LOG_TRACE( logger, "absMaxVal, val[" << i << ", " << jj << "] = " << val )
            }
        }

        #pragma omp critical
        {
            LAMA_LOG_DEBUG( logger, "absMaxVal, threadVal = " << threadVal << ", maxVal = " << maxValue )

            if ( threadVal > maxValue )
            {
                maxValue = threadVal;
            }
        }
    }

    LAMA_LOG_DEBUG( logger, "absMaxVal, maxVal = " << maxValue )

    return maxValue;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void OpenMPELLUtils::getRow(
    OtherValueType *row,
    const IndexType i,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType *values )
{
    LAMA_LOG_DEBUG( logger, "get row #i = " << i )

    #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
    for ( IndexType j = 0; j < numColumns; ++j )
    {
        row[j] = 0.0;
    }
    #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
    for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
    {
        IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
        row[ellJA[pos]] = static_cast<OtherValueType>( values[pos] );
    }
}

template<typename ValueType,typename OtherValueType>
OtherValueType OpenMPELLUtils::getValue(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[] )
{
    LAMA_LOG_TRACE( logger, "get value i = " << i << ", j = " << j )

    for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
    {
        IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );

        if ( ellJA[pos] == j )
        {
            return static_cast<OtherValueType>( ellValues[pos] );
        }
    }

    return 0.0;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ELLValueType,typename CSRValueType>
void OpenMPELLUtils::getCSRValues(
    IndexType csrJA[],
    CSRValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ELLValueType ellValues[] )
{
    LAMA_LOG_INFO( logger,
                   "get CSRValues<" << Scalar::getType<ELLValueType>()
                    << ", " << Scalar::getType<CSRValueType>() << ">"
                    << ", #rows = " << numRows )

    // parallelization possible as offset array csrIA is available

    #pragma omp parallel 
    {
        LAMA_REGION( "OpenMP.ELL->CSR_values" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType rowSize = ellSizes[i];
            IndexType offset = csrIA[i];

            // just make sure that csrIA and ellSizes really fit with each other

            LAMA_ASSERT_EQUAL_DEBUG( csrIA[i] + rowSize, csrIA[i+1] )

            for ( IndexType jj = 0; jj < rowSize; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                csrJA[offset + jj] = ellJA[pos];
                csrValues[offset + jj] = static_cast<CSRValueType>( ellValues[pos] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ELLValueType,typename CSRValueType>
void OpenMPELLUtils::setCSRValues(
    IndexType ellJA[],
    ELLValueType ellValues[],
    const IndexType ellSizes[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const CSRValueType csrValues[] )
{
    LAMA_LOG_INFO( logger,
                   "set CSRValues<" << Scalar::getType<ELLValueType>() << ", " << Scalar::getType<CSRValueType>() << ">" << ", #rows = " << numRows << ", #values/row = " << numValuesPerRow )

    // parallelization possible as offset array csrIA is available

    #pragma omp parallel 
    {
        LAMA_REGION( "OpenMP.ELL<-CSR_values" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
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
                ellValues[pos] = static_cast<ELLValueType>( csrValues[offset + jj] );
            }

            // fill up the remaining entries with something useful

            for ( IndexType jj = rowSize; jj < numValuesPerRow; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                ellJA[pos] = j; // last used column index
                ellValues[pos] = 0.0; // zero entry
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
    LAMA_LOG_INFO( logger, "fill ELLValues<" << Scalar::getType<ValueType>() )

    #pragma omp parallel 
    {
        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
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
                ellValues[pos] = 0.0; // zero entry
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::compressIA(
    const IndexType IA[],
    const IndexType JA[],
    const ValueType ellValues[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const ValueType eps,
    IndexType newIA[] )
{
    LAMA_LOG_INFO( logger, "compressIA with eps = " << eps )

    #pragma omp parallel
    {
        #pragma omp for
        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType length = IA[i];
            for ( IndexType jj = 0; jj < IA[i]; jj++ )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );

                if ( JA[pos] == i )
                {
                    continue;
                }
                if ( std::abs( ellValues[pos] ) <= eps )
                {
                    length--;
                }
            }
            newIA[i] = length;
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::compressValues(
    const IndexType IA[],
    const IndexType JA[],
    const ValueType values[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const ValueType eps,
    const IndexType newNumValuesPerRow,
    IndexType newJA[],
    ValueType newValues[] )
{
    LAMA_LOG_INFO( logger, "compressValues with eps = " << eps )

    #pragma omp parallel
    {
        #pragma omp for
        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType gap = 0;
            for ( IndexType j = 0; j < IA[i]; j++ )
            {
                IndexType pos = ellindex( i, j, numRows, numValuesPerRow );
                if ( std::abs( values[pos] ) <= eps && JA[pos] != i )
                {
                    gap++;
                    continue;
                }
                IndexType newpos = ellindex( i, j - gap, numRows, newNumValuesPerRow );
                newValues[ newpos ] = values[pos];
                newJA[ newpos ] = JA[pos];
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

void OpenMPELLUtils::matrixMultiplySizes(
    IndexType cSizes[],
    const IndexType aNumRows,
    const IndexType UNUSED( aNumColumns ),
    const IndexType bNumRows,
    const bool UNUSED( diagonalProperty ),
    const IndexType aSizes[],
    const IndexType aJA[],
    const IndexType aNumValuesPerRow,
    const IndexType bSizes[],
    const IndexType bJA[],
    const IndexType bNumValuesPerRow )
{
    LAMA_LOG_INFO( logger, "matrixMultiplySizes with numRows A = " << aNumRows << " and numRows B = " << bNumRows )

    #pragma omp parallel
    {
        #pragma omp for
        for ( IndexType i = 0; i < aNumRows; i++ )
        {
            std::set<IndexType> newElements;
            std::pair<std::set<IndexType>::iterator,bool> ret;
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
    const IndexType UNUSED( aNumColumns ),
    const IndexType bNumRows,
    const bool UNUSED( diagonalProperty ),
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
    LAMA_LOG_INFO( logger, "matrix multiply with numRows A = " << aNumRows << " and numRows B = " << bNumRows )

    #pragma omp parallel
    {
        #pragma omp for
        for ( IndexType i = 0; i < aNumRows; i++ )
        {
            std::set<IndexType> jaRow;
            std::map<IndexType,ValueType> valuesRow;
            std::pair<std::set<IndexType>::iterator,bool> ret;
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
                        valuesRow.insert( std::pair<IndexType,ValueType>( kk, mult ) );
                    }
                    else
                    {
                        valuesRow[kk] += mult;
                    }
                }
            }

            std::set<IndexType>::iterator jaIter;
            typename std::map<IndexType,ValueType>::iterator valuesIter;

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
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

void OpenMPELLUtils::matrixAddSizes(
    IndexType cSizes[],
    const IndexType m,
    const IndexType UNUSED( n ),
    const bool UNUSED( diagonalProperty ),
    const IndexType aSizes[],
    const IndexType aJA[],
    const IndexType aNumValuesPerRow,
    const IndexType bSizes[],
    const IndexType bJA[],
    const IndexType bNumValuesPerRow )
{
    LAMA_LOG_INFO( logger, "matrixAddSizes A + B, #rows = " << m )

    #pragma omp parallel
    {
        #pragma omp for
        for ( IndexType i = 0; i < m; i++ )
        {
            std::set<IndexType> iaRow;
            std::pair<std::set<IndexType>::iterator,bool> ret;
            IndexType length = 0;
            for ( IndexType j = 0; j < aSizes[i]; j++ )
            {
                IndexType posA = ellindex( i, j, m, aNumValuesPerRow );
                iaRow.insert( aJA[posA] );
                length++;
            }
            for ( IndexType j = 0; j < bSizes[i]; j++ )
            {
                IndexType posB = ellindex( i, j, m, bNumValuesPerRow );
                ret = iaRow.insert( bJA[posB] );
                if ( ret.second == true )
                {
                    length++;
                }
            }
            cSizes[i] = length;
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
    const IndexType UNUSED( n ),
    const bool UNUSED( diagonalProperty ),
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
    LAMA_LOG_INFO( logger, "matrixAdd C = " << alpha << " * A + " << beta << " * B, #rows = " << m )

    #pragma omp parallel
    {
        #pragma omp for
        for ( IndexType i = 0; i < m; i++ )
        {
            std::set<IndexType> jaRow;
            std::map<IndexType,ValueType> valuesRow;
            std::pair<std::set<IndexType>::iterator,bool> ret;

            for ( IndexType j = 0; j < aSizes[i]; j++ )
            {
                IndexType posA = ellindex( i, j, m, aNumValuesPerRow );
                jaRow.insert( aJA[posA] );
                valuesRow.insert( std::pair<IndexType,ValueType>( aJA[posA], aValues[posA] ) );
            }

            for ( IndexType j = 0; j < bSizes[i]; j++ )
            {
                IndexType posB = ellindex( i, j, m, bNumValuesPerRow );
                ret = jaRow.insert( bJA[posB] );
                if ( ret.second == true )
                {
                    valuesRow.insert( std::pair<IndexType,ValueType>( bJA[posB], beta * bValues[posB] ) );
                }
                else
                {
                    valuesRow[bJA[posB]] += beta * bValues[posB];
                }
            }

            std::set<IndexType>::iterator jaIter;
            typename std::map<IndexType,ValueType>::iterator valuesIter;

            jaIter = jaRow.begin();
            valuesIter = valuesRow.begin();

            for ( IndexType j = 0; j < cSizes[i]; j++ )
            {
                // Note: cNumRows == aNumRows
                IndexType posC = ellindex( i, j, m, cNumValuesPerRow );
                cJA[posC] = *jaIter;
                cValues[posC] = ( *valuesIter ).second;
                jaIter++;
                valuesIter++;
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
    const ValueType omega,
    class SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "jacobi<" << Scalar::getType<ValueType>() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    if ( syncToken != NULL )
    {
        LAMA_LOG_ERROR( logger, "jacobi called asynchronously, not supported here" )
    }

    const ValueType oneMinusOmega = static_cast<ValueType>( 1.0 - omega );
    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.ELL.jacobi" )
        #pragma omp for schedule(LAMA_OMP_SCHEDULE)
        for ( IndexType i = 0; i < numRows; i++ )
        {
            ValueType temp = rhs[i];
            IndexType pos  = ellindex( i, 0, numRows, ellNumValuesPerRow );
            ValueType diag = ellValues[pos]; //getDiagonal
            for ( IndexType j = 1; j < ellSizes[i]; j++ )
            {
                pos = ellindex( i, j, numRows, ellNumValuesPerRow );
                temp -= ellValues[pos] * oldSolution[ellJA[pos]];
            }

            if ( omega == 1.0 )
            {
                solution[i] = temp / diag;
            }
            else if ( omega == 0.5 )
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
    const ValueType omega,
    class SyncToken* syncToken )
{
    if ( syncToken != NULL )
    {
        LAMA_LOG_WARN( logger, "jacobi called asynchronously, not supported here" )
    }

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.ELL.jacobiHalo" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
        for ( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            IndexType i = ii; // rowIndexes == NULL stands for all rows

            if ( rowIndexes )
            {
                i = rowIndexes[ii];
            }

            ValueType temp = 0.0;

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
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "normalGEMV<" << Scalar::getType<ValueType>()
                   << ", #threads = " << omp_get_max_threads()
                   << ">, result[" << numRows << "] = " << alpha 
                   << " * A( ell, #maxNZ/row = " << numValuesPerRow << " ) * x + " << beta << " * y " )

    if ( numValuesPerRow == 0 )
    {
        LAMA_THROWEXCEPTION( "normalGEMV should not have been called, no entries" )

        // only compute: result = beta * y
    }

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "Host: asynchronous execution by task should be done at higher level" )
    }

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.ELL.normalGEMV" )

        #pragma omp for schedule(LAMA_OMP_SCHEDULE)
        for ( IndexType i = 0; i < numRows; ++i )
        {
            ValueType temp = 0.0;

            for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                IndexType j = ellJA[pos];
                LAMA_LOG_TRACE( logger,
                                "temp += dataAccess[i + jj * numRows] * xAccess[j];, jj = " << jj << ", j = " << j )
                LAMA_LOG_TRACE( logger, ", dataAccess[i + jj * numRows] = " << ellValues[ pos ] )
                LAMA_LOG_TRACE( logger, ", xAccess[j] = " << x[ j ] )
                temp += ellValues[pos] * x[j];
            }

            LAMA_LOG_TRACE( logger, "row = " << i << ", temp = " << temp )

            if ( static_cast<ValueType>( 0 ) == beta )
            {
                // must be handled separately as y[i] might be uninitialized

                result[i] = alpha * temp;
            }
            else if ( static_cast<ValueType>( 1 ) == alpha )
            {
                result[i] = temp + beta * y[i];
            }
            else
            {
                result[i] = alpha * temp + beta * y[i];
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::sparseGEMV(
    ValueType result[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const ValueType alpha,
    const ValueType x[],
    const IndexType numNonZeroRows,
    const IndexType rowIndexes[],
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    SyncToken* syncToken )
{
    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "Host: asynchronous execution by task should be done at higher level" )
    }

    LAMA_LOG_INFO( logger,
                   "sparseGEMV<" << Scalar::getType<ValueType>() << ">, n = " << numRows 
                       << ", nonZeroRows = " << numNonZeroRows << ", alpha = " << alpha )

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.ELL.sparseGEMV" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
        for ( IndexType ii = 0; ii < numNonZeroRows; ++ii )
        {
            IndexType i = rowIndexes[ii];

            //result is not initialized for performance reasons
            ValueType temp = 0.0;

            for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                IndexType j = ellJA[pos];
                temp += ellValues[pos] * x[j];
            }

            if ( 1 == alpha )
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

/* ------------------------------------------------------------------------------------------------------------------ */

void OpenMPELLUtils::setInterface( ELLUtilsInterface& ELLUtils )
{
    LAMA_INTERFACE_REGISTER( ELLUtils, countNonEmptyRowsBySizes )
    LAMA_INTERFACE_REGISTER( ELLUtils, setNonEmptyRowsBySizes )
    LAMA_INTERFACE_REGISTER( ELLUtils, hasDiagonalProperty )
    LAMA_INTERFACE_REGISTER( ELLUtils, check )

    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getRow, float, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getRow, float, double )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getRow, double, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getRow, double, double )

    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getValue, float, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getValue, float, double )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getValue, double, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getValue, double, double )

    LAMA_INTERFACE_REGISTER_TT( ELLUtils, scaleValue, float, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, scaleValue, double, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, scaleValue, float, double )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, scaleValue, double, double )

    LAMA_INTERFACE_REGISTER_TT( ELLUtils, setCSRValues, float, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, setCSRValues, float, double )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, setCSRValues, double, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, setCSRValues, double, double )

    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getCSRValues, float, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getCSRValues, float, double )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getCSRValues, double, float )
    LAMA_INTERFACE_REGISTER_TT( ELLUtils, getCSRValues, double, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, absMaxVal, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, absMaxVal, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, compressIA, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, compressIA, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, compressValues, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, compressValues, double )

    LAMA_INTERFACE_REGISTER( ELLUtils, matrixMultiplySizes )
    LAMA_INTERFACE_REGISTER( ELLUtils, matrixAddSizes )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, matrixMultiply, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, matrixMultiply, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, matrixAdd, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, matrixAdd, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, normalGEMV, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, normalGEMV, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, sparseGEMV, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, sparseGEMV, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, jacobi, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, jacobi, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, jacobiHalo, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, jacobiHalo, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, fillELLValues, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, fillELLValues, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the ELLUtils routines                             */
/* --------------------------------------------------------------------------- */

bool OpenMPELLUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::Host );
    setInterface( interface.ELLUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPELLUtils::initialized = registerInterface();

} // namespace lama
