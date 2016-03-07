/**
 * @file OpenMPELLUtils.cpp
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
 * @brief Implementation of ELL utilities with OpenMP
 * @author Thomas Brandes
 * @date 04.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/sparsekernel/openmp/OpenMPELLUtils.hpp>

// local library
#include <scai/sparsekernel/ELLKernelTrait.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>


#include <scai/common/bind.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/mepr/Container.hpp>


// std
#include <set>
#include <map>
#include <cmath>

namespace scai
{

using common::TypeTraits;
using tasking::TaskSyncToken;

namespace sparsekernel
{

/* ------------------------------------------------------------------------------------------------------------------ */

SCAI_LOG_DEF_LOGGER( OpenMPELLUtils::logger, "OpenMP.ELLUtils" )

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

    SCAI_LOG_INFO( logger, "#non-zero rows = " << counter << ", counted by sizes" )

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

    SCAI_ASSERT_EQUAL_DEBUG( counter, numNonEmptyRows )

    SCAI_LOG_INFO( logger, "#non-zero rows = " << counter << ", set by sizes" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

bool OpenMPELLUtils::hasDiagonalProperty( const IndexType numDiagonals, const IndexType ellJA[] )
{
    SCAI_LOG_INFO( logger, "hasDiagonalProperty, #numDiagonals = " << numDiagonals )

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

template<typename ValueType, typename OtherValueType>
void OpenMPELLUtils::scaleValue(
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    ValueType ellValues[],
    const OtherValueType values[] )
{
    SCAI_LOG_INFO( logger,
                   "scaleValue<" << TypeTraits<ValueType>::id() << ", " << TypeTraits<OtherValueType>::id() << ">" << ", #numRows = " << numRows )

    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

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
    SCAI_LOG_INFO( logger,
                   "check # numRows = " << numRows << ", numValuesPerRow = " << numValuesPerRow << ", numColumns = " << numColumns )

    if ( numRows > 0 )
    {
        bool integrityIA = true;
        bool integrityJA = true;

        #pragma omp parallel for reduction( && : integrityIA, integrityJA ) schedule( SCAI_OMP_SCHEDULE )

        for ( IndexType i = 0; i < numRows; i++ )
        {
            if ( ellSizes[i] >= 0 && ellSizes[i] <= numValuesPerRow )
            {
                for ( IndexType jj = 0; jj < ellSizes[i]; jj++ )
                {
                    IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                    IndexType j = ellJA[pos];
                    integrityJA = integrityJA && ( 0 <= j && j < numColumns );
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
    else
    {
        SCAI_ASSERT_EQUAL_ERROR( 0, numValuesPerRow )
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

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                ValueType val = common::Math::abs( values[pos] );

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

template<typename ValueType, typename OtherValueType>
void OpenMPELLUtils::getRow(
    OtherValueType row[],
    const IndexType i,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType values[] )
{
    SCAI_LOG_DEBUG( logger, "get row #i = " << i )

    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

    for ( IndexType j = 0; j < numColumns; ++j )
    {
        row[j] = static_cast<OtherValueType>( 0.0 );
    }

    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

    for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
    {
        IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
        row[ellJA[pos]] = static_cast<OtherValueType>( values[pos] );
    }
}

template<typename ValueType>
ValueType OpenMPELLUtils::getValue(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[] )
{
    SCAI_LOG_TRACE( logger, "get value i = " << i << ", j = " << j )

    ValueType val = 0.0;

    for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
    {
        IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );

        if ( ellJA[pos] == j )
        {
            val = ellValues[pos];
            break;
        }
    }

    return val;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ELLValueType, typename CSRValueType>
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
    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << TypeTraits<ELLValueType>::id() << ", " << TypeTraits<CSRValueType>::id() << ">" << ", #rows = " << numRows )

    // parallelization possible as offset array csrIA is available

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.ELL->CSR_values" )

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

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
                csrValues[offset + jj] = static_cast<CSRValueType>( ellValues[pos] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ELLValueType, typename CSRValueType>
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
    SCAI_LOG_INFO( logger,
                   "set CSRValues<" << TypeTraits<ELLValueType>::id() << ", " << TypeTraits<CSRValueType>::id() << ">" << ", #rows = " << numRows << ", #values/row = " << numValuesPerRow )

    // parallelization possible as offset array csrIA is available

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.ELL<-CSR_values" )

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

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
                ellValues[pos] = static_cast<ELLValueType>( 0.0 ); // zero entry
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
    SCAI_LOG_INFO( logger, "fill ELLValues<" << TypeTraits<ValueType>::id() )

    #pragma omp parallel
    {
        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

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
                ellValues[pos] = static_cast<ValueType>( 0.0 ); // zero entry
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
    SCAI_LOG_INFO( logger, "compressIA with eps = " << eps )

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

                if ( common::Math::abs( ellValues[pos] ) <= eps )
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
    SCAI_LOG_INFO( logger, "compressValues with eps = " << eps )

    #pragma omp parallel
    {
        #pragma omp for

        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType gap = 0;

            for ( IndexType j = 0; j < IA[i]; j++ )
            {
                IndexType pos = ellindex( i, j, numRows, numValuesPerRow );

                if ( common::Math::abs( values[pos] ) <= eps && JA[pos] != i )
                {
                    gap++;
                    continue;
                }

                IndexType newpos = ellindex( i, j - gap, numRows, newNumValuesPerRow );
                newValues[newpos] = values[pos];
                newJA[newpos] = JA[pos];
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
    SCAI_LOG_INFO( logger, "matrixAddSizes A + B, #rows = " << m )

    #pragma omp parallel
    {
        #pragma omp for

        for ( IndexType i = 0; i < m; i++ )
        {
            std::set<IndexType> iaRow;
            std::pair<std::set<IndexType>::iterator, bool> ret;
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
    SCAI_LOG_INFO( logger, "matrixAdd C = " << alpha << " * A + " << beta << " * B, #rows = " << m )

    #pragma omp parallel
    {
        #pragma omp for

        for ( IndexType i = 0; i < m; i++ )
        {
            std::set<IndexType> jaRow;
            std::map<IndexType, ValueType> valuesRow;
            std::pair<std::set<IndexType>::iterator, bool> ret;

            for ( IndexType j = 0; j < aSizes[i]; j++ )
            {
                IndexType posA = ellindex( i, j, m, aNumValuesPerRow );
                jaRow.insert( aJA[posA] );
                valuesRow.insert( std::pair<IndexType, ValueType>( aJA[posA], aValues[posA] ) );
            }

            for ( IndexType j = 0; j < bSizes[i]; j++ )
            {
                IndexType posB = ellindex( i, j, m, bNumValuesPerRow );
                ret = jaRow.insert( bJA[posB] );

                if ( ret.second == true )
                {
                    valuesRow.insert( std::pair<IndexType, ValueType>( bJA[posB], beta * bValues[posB] ) );
                }
                else
                {
                    valuesRow[bJA[posB]] += beta * bValues[posB];
                }
            }

            std::set<IndexType>::iterator jaIter;
            typename std::map<IndexType, ValueType>::iterator valuesIter;

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
    const ValueType omega )
{
    SCAI_LOG_INFO( logger,
                   "jacobi<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken != NULL )
    {
        SCAI_LOG_ERROR( logger, "jacobi called asynchronously, not supported here" )
    }

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.ELL.jacobi" )
        #pragma omp for schedule(SCAI_OMP_SCHEDULE)

        for ( IndexType i = 0; i < numRows; i++ )
        {
            ValueType temp = rhs[i];
            IndexType pos = ellindex( i, 0, numRows, ellNumValuesPerRow );
            ValueType diag = ellValues[pos]; //getDiagonal

            for ( IndexType j = 1; j < ellSizes[i]; j++ )
            {
                pos = ellindex( i, j, numRows, ellNumValuesPerRow );
                temp -= ellValues[pos] * oldSolution[ellJA[pos]];
            }

            if ( omega == scai::common::constants::ONE )
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

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

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
void OpenMPELLUtils::normalGEMV_a(
    ValueType result[],
    const std::pair<ValueType, const ValueType*> ax,
    const std::pair<ValueType, const ValueType*> by,
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[] )
{
    normalGEMV( result, ax.first, ax.second, by.first, by.second,
                numRows, numValuesPerRow, ellSizes, ellJA, ellValues );
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
    const ValueType ellValues[] )
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

        syncToken->run( common::bind( normalGEMV_a<ValueType>, result, 
                                      std::pair<ValueType, const ValueType*>( alpha, x ),
                                      std::pair<ValueType, const ValueType*> ( beta, y ),
                                      numRows, numValuesPerRow, ellSizes, ellJA, ellValues ) );
        return;
    }

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads() 
                   << ">, result[" << numRows << "] = " << alpha << " * A( ell, #maxNZ/row = " << numValuesPerRow << " ) * x + " << beta << " * y " )

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.ELL.normalGEMV" )

        #pragma omp for schedule(SCAI_OMP_SCHEDULE)

        for ( IndexType i = 0; i < numRows; ++i )
        {
            ValueType temp = static_cast<ValueType>( 0.0 );

            for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                IndexType j = ellJA[pos];
                SCAI_LOG_TRACE( logger,
                                "temp += dataAccess[i + jj * numRows] * xAccess[j];, jj = " << jj << ", j = " << j )
                SCAI_LOG_TRACE( logger, ", dataAccess[i + jj * numRows] = " << ellValues[ pos ] )
                SCAI_LOG_TRACE( logger, ", xAccess[j] = " << x[ j ] )
                temp += ellValues[pos] * x[j];
            }

            SCAI_LOG_TRACE( logger, "row = " << i << ", temp = " << temp )

            if ( beta == scai::common::constants::ZERO )
            {
                // must be handled separately as y[i] might be uninitialized

                result[i] = alpha * temp;
            }
            else if ( alpha == scai::common::constants::ONE )
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
void OpenMPELLUtils::sparseGEMV_a(
    ValueType result[],
    const std::pair<ValueType, const ValueType*> ax,
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const std::pair<IndexType, const IndexType*> rows,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[] )
{
    sparseGEMV( result, ax.first, ax.second, numRows, numValuesPerRow,
                rows.first, rows.second, ellSizes, ellJA, ellValues );
}

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
    const ValueType ellValues[] )
{
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        syncToken->run( common::bind( sparseGEMV_a<ValueType>, result,
                                      std::pair<ValueType, const ValueType*>( alpha, x ),
                                      numRows, numValuesPerRow,
                                      std::pair<IndexType, const IndexType*>( numNonZeroRows, rowIndexes ),
                                      ellSizes, ellJA, ellValues ) );
        return; 
    }

    SCAI_LOG_INFO( logger,
                   "sparseGEMV<" << TypeTraits<ValueType>::id() << ">, n = " << numRows 
                    << ", nonZeroRows = " << numNonZeroRows << ", alpha = " << alpha )

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.ELL.sparseGEMV" )

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for ( IndexType ii = 0; ii < numNonZeroRows; ++ii )
        {
            IndexType i = rowIndexes[ii];

            //result is not initialized for performance reasons
            ValueType temp = static_cast<ValueType>( 0.0 );

            for ( IndexType jj = 0; jj < ellSizes[i]; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                IndexType j = ellJA[pos];
                temp += ellValues[pos] * x[j];
            }

            if ( alpha == scai::common::constants::ONE )
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

template<typename ValueType>
void OpenMPELLUtils::normalGEVM(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType UNUSED( numValuesPerRow ),
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[] )
{
    SCAI_LOG_INFO( logger,
                   "normalGEVM<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads() << ">, result[" << numColumns << "] = " << alpha << " * x * A + " << beta << " * y " )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution not supported here" )
    }

    //#pragma omp parallel
    {
        SCAI_REGION( "OpenMP.ELL.normalGEVM" )

        //#pragma omp for schedule(SCAI_OMP_SCHEDULE)
        for ( IndexType i = 0; i < numColumns; ++i )
        {
            ValueType temp = static_cast<ValueType>( 0.0 );

            for ( IndexType j = 0; j < numRows; ++j )
            {
                for ( IndexType k = 0; k < ellSizes[j]; ++k )
                {
                    if ( ellJA[k * numRows + j] == i )
                    {
                        SCAI_LOG_TRACE( logger, "temp += dataAccess[k * numRows + j] * xAccess[j]; j = " << j )
                        SCAI_LOG_TRACE( logger, ", dataAccess[k * numRows + j] = " << ellValues[ k * numRows + j ] )
                        SCAI_LOG_TRACE( logger, ", xAccess[j] = " << x[ j ] )

                        temp += ellValues[k * numRows + j] * x[j];
                    }
                }
            }

            SCAI_LOG_TRACE( logger, "column = " << i << ", temp = " << temp )

            if ( beta == scai::common::constants::ZERO )
            {
                // must be handled separately as y[i] might be uninitialized

                result[i] = alpha * temp;
            }
            else if ( alpha == scai::common::constants::ONE )
            {
                result[i] = temp + beta * y[i];
            }
            else
            {
                result[i] = alpha * temp + beta * y[i];
            }
        }
    }

    if ( SCAI_LOG_TRACE_ON( logger ) )
    {
        std::cout << "NormalGEVM: result = ";

        for ( IndexType i = 0; i < numColumns; ++i )
        {
            std::cout << " " << result[i];
        }

        std::cout << std::endl;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void OpenMPELLUtils::sparseGEVM(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType UNUSED( numRows ),
    const IndexType numColumns,
    const IndexType UNUSED( numValuesPerRow ),
    const IndexType numNonZeroRows,
    const IndexType rowIndexes[],
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[] )
{
    SCAI_LOG_INFO( logger,
                   "sparseGEVM<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads() << ">, result[" << numColumns << "] = " << alpha << " * x * A " )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.ELL.sparseGEVM" )

        #pragma omp for schedule(SCAI_OMP_SCHEDULE)

        for ( IndexType i = 0; i < numColumns; ++i )
        {
            ValueType temp = static_cast<ValueType>( 0.0 );

            for ( IndexType jj = 0; jj < numNonZeroRows; ++jj )
            {
                IndexType j = rowIndexes[jj];

                for ( IndexType k = 0; k < ellSizes[j]; ++k )
                {
                    if ( ellJA[k * numNonZeroRows + j] == i )
                    {
                        SCAI_LOG_TRACE( logger, "temp += dataAccess[k * numNonZeroRows + j] * xAccess[j]; i = " << j )
                        SCAI_LOG_TRACE( logger,
                                        ", dataAccess[k * numNonZeroRows + j] = " << ellValues[ k * numNonZeroRows + j ] )
                        SCAI_LOG_TRACE( logger, ", xAccess[j] = " << x[ j ] )

                        temp += ellValues[k * numNonZeroRows + j] * x[j];
                    }
                }

            }

            if ( alpha == scai::common::constants::ONE )
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

void OpenMPELLUtils::Registrator::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    common::context::ContextType ctx = common::context::Host;

    SCAI_LOG_INFO( logger, "register ELLtils OpenMP-routines for Host at kernel registry [" << flag << "]" )

    KernelRegistry::set<ELLKernelTrait::countNonEmptyRowsBySizes>( countNonEmptyRowsBySizes, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::setNonEmptyRowsBySizes>( setNonEmptyRowsBySizes, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::hasDiagonalProperty>( hasDiagonalProperty, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::check>( check, ctx, flag );

    KernelRegistry::set<ELLKernelTrait::matrixMultiplySizes>( matrixMultiplySizes, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::matrixAddSizes>( matrixAddSizes, ctx, flag );
}

template<typename ValueType>
void OpenMPELLUtils::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    common::context::ContextType ctx = common::context::Host;

    SCAI_LOG_INFO( logger, "register ELLUtils OpenMP-routines for Host at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<ELLKernelTrait::absMaxVal<ValueType> >( absMaxVal, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::compressIA<ValueType> >( compressIA, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::compressValues<ValueType> >( compressValues, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::matrixAdd<ValueType> >( matrixAdd, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::matrixMultiply<ValueType> >( matrixMultiply, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::sparseGEMV<ValueType> >( sparseGEMV, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::normalGEVM<ValueType> >( normalGEVM, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::sparseGEVM<ValueType> >( sparseGEVM, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::getValue<ValueType> >( getValue, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::fillELLValues<ValueType> >( fillELLValues, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void OpenMPELLUtils::RegistratorVO<ValueType, OtherValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    common::context::ContextType ctx = common::context::Host;

    SCAI_LOG_INFO( logger, "register ELLUtils OpenMP-routines for Host at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )

    KernelRegistry::set<ELLKernelTrait::getRow<ValueType, OtherValueType> >( getRow, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::scaleValue<ValueType, OtherValueType> >( scaleValue, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::setCSRValues<ValueType, OtherValueType> >( setCSRValues, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::getCSRValues<ValueType, OtherValueType> >( getCSRValues, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPELLUtils::OpenMPELLUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;

    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, ARITHMETIC_HOST_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, ARITHMETIC_HOST_LIST, ARITHMETIC_HOST_LIST>::call( flag );
}

OpenMPELLUtils::~OpenMPELLUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, ARITHMETIC_HOST_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, ARITHMETIC_HOST_LIST, ARITHMETIC_HOST_LIST>::call( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPELLUtils OpenMPELLUtils::guard;

} /* end namespace sparsekernel */

} /* end namespace scai */
