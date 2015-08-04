/**
 * @file MICELLUtils.cpp
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
 * @brief Implementation of ELL utilities with MIC
 * @author Thomas Brandes
 * @date 04.07.2012
 * @since 1.1.0
 */

// hpp
#include <lama/mic/MICELLUtils.hpp>

// others
#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

#include <lama/NoSyncToken.hpp>
#include <lama/mic/MICSyncToken.hpp>
#include <lama/mic/MICContext.hpp>

// macros
#include <lama/macros/unused.hpp>

// tracing
#include <tracing.hpp>

// boost
#include <common/bind.hpp>

// stl
#include <set>
#include <map>
#include <cmath>

namespace lama
{

/* ------------------------------------------------------------------------------------------------------------------ */

LAMA_LOG_DEF_LOGGER( MICELLUtils::logger, "MIC.ELLUtils" )

/* ------------------------------------------------------------------------------------------------------------------ */

IndexType MICELLUtils::countNonEmptyRowsBySizes( const IndexType sizes[], const IndexType numRows )
{
    IndexType counter = 0;

    const void* sizesPtr = sizes;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( sizesPtr, numRows ), out( counter )
    {
        counter = 0;

        const IndexType* sizes = static_cast<const IndexType*>( sizesPtr );

        #pragma omp parallel for reduction( +:counter )

        for( IndexType i = 0; i < numRows; ++i )
        {
            if( sizes[i] > 0 )
            {
                counter++;
            }
        }
    }

    LAMA_LOG_INFO( logger, "#non-zero rows = " << counter << ", counted by sizes" )

    return counter;
}

/* ------------------------------------------------------------------------------------------------------------------ */

void MICELLUtils::setNonEmptyRowsBySizes(
    IndexType rowIndexes[],
    const IndexType numNonEmptyRows,
    const IndexType sizes[],
    const IndexType numRows )
{
    IndexType counter = 0;

    const void* sizesPtr = sizes;
    void* rowIndexesPtr = rowIndexes;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( rowIndexesPtr, sizesPtr, numRows ), out( counter )
    {
        IndexType* rowIndexes = static_cast<IndexType*>( rowIndexesPtr );
        const IndexType* sizes = static_cast<const IndexType*>( sizesPtr );

        counter = 0;

        // Note: this routine is not easy to parallelize, no offsets for rowIndexes available

        for( IndexType i = 0; i < numRows; ++i )
        {
            if( sizes[i] > 0 )
            {
                rowIndexes[counter] = i;
                counter++;
            }
        }
    }

    LAMA_ASSERT_EQUAL_DEBUG( counter, numNonEmptyRows )

    LAMA_LOG_INFO( logger, "#non-zero rows = " << counter << ", set by sizes" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

bool MICELLUtils::hasDiagonalProperty( const IndexType numDiagonals, const IndexType ellJA[] )
{
    LAMA_LOG_INFO( logger, "hasDiagonalProperty, #numDiagonals = " << numDiagonals )

    if( numDiagonals == 0 )
    {
        return false;
    }

    bool diagonalProperty = true;

    const void* ellJAPtr = ellJA;

#pragma offload target( MIC ) in( ellJAPtr, numDiagonals ), out( diagonalProperty )
    {
        const IndexType* ellJA = static_cast<const IndexType*>( ellJAPtr );

        diagonalProperty = true;

        #pragma omp parallel for reduction( && : diagonalProperty )

        for( IndexType i = 0; i < numDiagonals; ++i )
        {
            if( !diagonalProperty )
            {
                continue;
            }

            if( ellJA[i] != i )
            {
                diagonalProperty = false;
            }
        }
    }

    return diagonalProperty;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void MICELLUtils::scaleValue(
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    ValueType ellValues[],
    const OtherValueType values[] )
{
    LAMA_LOG_INFO( logger,
                   "scaleValue<" << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << ">" << ", #numRows = " << numRows )

    void* ellValuesPtr = ellValues;
    const void* ellSizesPtr = ellSizes;
    const void* valuesPtr = values;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( numRows, numValuesPerRow, ellSizesPtr, ellValuesPtr, valuesPtr )
    {
        ValueType* ellValues = static_cast<ValueType*>( ellValuesPtr );

        const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );
        const OtherValueType* values = static_cast<const OtherValueType*>( valuesPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; i++ ) //rows
        {
            for( IndexType jj = 0; jj < ellSizes[i]; jj++ ) //elements in row
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                ellValues[pos] *= static_cast<ValueType>( values[i] );
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

void MICELLUtils::check(
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType numColumns,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const char* msg )
{
    LAMA_LOG_INFO( logger,
                   "check # numRows = " << numRows << ", numValuesPerRow = " << numValuesPerRow << ", numColumns = " << numColumns )

    if( numRows > 0 )
    {
        bool integrityIA = true;
        bool integrityJA = true;

        int device = MICContext::getCurrentDevice();

        const void* ellSizesPtr = ellSizes;
        const void* ellJAPtr = ellJA;

#pragma offload target( mic : device ) in( ellSizesPtr, ellJAPtr, numRows, numValuesPerRow, numColumns ), out( integrityIA, integrityJA )
        {
            const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );
            const IndexType* ellJA = static_cast<const IndexType*>( ellJAPtr );

            integrityIA = true;
            integrityJA = true;

            #pragma omp parallel for reduction( && : integrityIA, integrityJA )

            for( IndexType i = 0; i < numRows; i++ )
            {
                if( ellSizes[i] >= 0 && ellSizes[i] <= numValuesPerRow )
                {
                    for( IndexType jj = 0; jj < ellSizes[i]; jj++ )
                    {
                        IndexType j = ellJA[jj * numRows + i];
                        integrityJA = integrityJA && ( 0 <= j && j < numColumns );
                    }
                }
                else
                {
                    integrityIA = false;
                }
            }
        }

        LAMA_ASSERT_ERROR( integrityIA, msg << ": ellSizes: at least one value out of range" );
        LAMA_ASSERT_ERROR( integrityJA, msg << ": ellJA: at least one value out of range" );
    }
    else
    {
        LAMA_ASSERT_ERROR( numValuesPerRow == 0,
                           msg << ": numRows is 0 but numValuesPerRow is: " << numValuesPerRow << ", should be 0" )
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
ValueType MICELLUtils::absMaxVal(
    const IndexType numRows,
    const IndexType /* numValuesPerRow */,
    const IndexType ia[],
    const ValueType values[] )
{
    ValueType maxValue = static_cast<ValueType>( 0.0 );

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>( 0.0 );

        #pragma omp for

        for( IndexType i = 0; i < numRows; ++i )
        {
            for( IndexType jj = 0; jj < ia[i]; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows );
                ValueType val = std::abs( values[pos] );

                if( val > threadVal )
                {
                    threadVal = val;
                }

                // LAMA_LOG_TRACE( logger, "absMaxVal, val[" << i << ", " << jj << "] = " << val )
            }
        }

        #pragma omp critical
        {
            LAMA_LOG_DEBUG( logger, "absMaxVal, threadVal = " << threadVal << ", maxVal = " << maxValue )

            if( threadVal > maxValue )
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
void MICELLUtils::getRow(
    OtherValueType row[],
    const IndexType i,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[] )
{
    LAMA_LOG_DEBUG( logger, "get row #i = " << i )

    const void* ellSizesPtr = ellSizes;
    const void* ellJAPtr = ellJA;
    const void* ellValuesPtr = ellValues;
    void* rowPtr = row;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( rowPtr, ellSizesPtr, ellJAPtr, ellValuesPtr, i, numRows, numColumns, numValuesPerRow )
    {
        const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );
        const IndexType* ellJA = static_cast<const IndexType*>( ellJAPtr );
        const ValueType* ellValues = static_cast<const ValueType*>( ellValuesPtr );

        OtherValueType* row = static_cast<OtherValueType*>( rowPtr );

        #pragma omp parallel for

        for( IndexType j = 0; j < numColumns; ++j )
        {
            row[j] = 0.0;
        }

        #pragma omp parallel for

        for( IndexType jj = 0; jj < ellSizes[i]; ++jj )
        {
            IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
            row[ellJA[pos]] = static_cast<OtherValueType>( ellValues[pos] );
        }
    }
}

template<typename ValueType,typename OtherValueType>
OtherValueType MICELLUtils::getValue(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[] )
{
    LAMA_LOG_TRACE( logger, "get value i = " << i << ", j = " << j )

    OtherValueType value = 0;

    const void* ellSizesPtr = ellSizes;
    const void* ellJAPtr = ellJA;
    const void* ellValuesPtr = ellValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), out( value ), in( ellSizesPtr, ellJAPtr, ellValuesPtr, \
                                                              i, j, numRows, numValuesPerRow )
    {
        const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );
        const IndexType* ellJA = static_cast<const IndexType*>( ellJAPtr );
        const ValueType* ellValues = static_cast<const ValueType*>( ellValuesPtr );

        value = 0; // new intialiation, has not been copied in

        for( IndexType jj = 0; jj < ellSizes[i]; ++jj )
        {
            IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );

            if( ellJA[pos] == j )
            {
                value = static_cast<OtherValueType>( ellValues[pos] );
            }
        }
    }

    return value;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ELLValueType,typename CSRValueType>
void MICELLUtils::getCSRValues(
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
                   "get CSRValues<" << common::getScalarType<ELLValueType>() << ", " << common::getScalarType<CSRValueType>() << ">" << ", #rows = " << numRows )

    LAMA_REGION( "MIC.ELL->CSR_values" )

    // parallelization possible as offset array csrIA is available

    void* csrJAPtr = csrJA;
    void* csrValuesPtr = csrValues;

    const void* ellSizesPtr = ellSizes;
    const void* ellJAPtr = ellJA;
    const void* ellValuesPtr = ellValues;
    const void* csrIAPtr = csrIA;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( ellSizesPtr, ellValuesPtr, ellJAPtr,  \
                                               csrIAPtr, csrJAPtr, csrValuesPtr, \
                                               numRows, numValuesPerRow )
    {
        const IndexType* ellJA = (IndexType*) ellJAPtr;
        const ELLValueType* ellValues = (ELLValueType*) ellValuesPtr;
        const IndexType* ellSizes = (IndexType*) ellSizesPtr;

        const IndexType* csrIA = (IndexType*) csrIAPtr;
        IndexType* csrJA = (IndexType*) csrJAPtr;
        CSRValueType* csrValues = (CSRValueType*) csrValuesPtr;

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; i++ )
        {
            IndexType rowSize = ellSizes[i];
            IndexType offset = csrIA[i];

            // just make sure that csrIA and ellSizes really fit with each other

            for( IndexType jj = 0; jj < rowSize; ++jj )
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
void MICELLUtils::setCSRValues(
    IndexType ellJA[],
    ELLValueType ellValues[],
    const IndexType ellSizes[],
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const CSRValueType csrValues[] )
{
    LAMA_REGION( "MIC.ELL<-CSR_values" )

    LAMA_LOG_INFO( logger,
                   "set CSRValues<" << common::getScalarType<ELLValueType>() << ", " << common::getScalarType<CSRValueType>() << ">" << ", #rows = " << numRows << ", #values/row = " << numValuesPerRow )

    const void* ellSizesPtr = ellSizes;
    void* ellJAPtr = ellJA;
    void* ellValuesPtr = ellValues;
    const void* csrIAPtr = csrIA;
    const void* csrJAPtr = csrJA;
    const void* csrValuesPtr = csrValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( ellSizesPtr, ellValuesPtr, ellJAPtr, \
                                               csrIAPtr, csrJAPtr, csrValuesPtr,  \
                                               numRows, numValuesPerRow )
    {
        IndexType* ellJA = static_cast<IndexType*>( ellJAPtr );
        ELLValueType* ellValues = static_cast<ELLValueType*>( ellValuesPtr );
        const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );

        const IndexType* csrIA = static_cast<const IndexType*>( csrIAPtr );
        const IndexType* csrJA = static_cast<const IndexType*>( csrJAPtr );
        const CSRValueType* csrValues = static_cast<const CSRValueType*>( csrValuesPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; i++ )
        {
            IndexType rowSize = ellSizes[i];
            IndexType offset = csrIA[i];
            IndexType j = 0; // will be last column index

            for( IndexType jj = 0; jj < rowSize; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                j = csrJA[offset + jj];
                ellJA[pos] = j;
                ellValues[pos] = static_cast<ELLValueType>( csrValues[offset + jj] );
            }

            // fill up the remaining entries with something useful

            for( IndexType jj = rowSize; jj < numValuesPerRow; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
                ellJA[pos] = j; // last used column index
                ellValues[pos] = 0; // zero entry
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void MICELLUtils::fillELLValues(
    IndexType ellJA[],
    ValueType ellValues[],
    const IndexType ellSizes[],
    const IndexType numRows,
    const IndexType numValuesPerRow )
{
    LAMA_LOG_INFO( logger, "fill ELLValues<" << typeid( ValueType ).name() )

    #pragma omp parallel
    {
        #pragma omp for

        for( IndexType i = 0; i < numRows; i++ )
        {
            IndexType rowSize = ellSizes[i];

            IndexType j = 0; // will be last column index

            if( rowSize > 0 && rowSize < numValuesPerRow )
            {
                IndexType pos = ellindex( i, rowSize - 1, numRows );
                j = ellJA[pos];
            }

            // fill up the remaining entries with something useful

            for( IndexType jj = rowSize; jj < numValuesPerRow; ++jj )
            {
                IndexType pos = ellindex( i, jj, numRows );
                ellJA[pos] = j; // last used column index
                ellValues[pos] = 0.0; // zero entry
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void MICELLUtils::compressIA(
    const IndexType IA[],
    const IndexType JA[],
    const ValueType values[],
    const IndexType numRows,
    const ValueType eps,
    IndexType newIA[] )
{
    LAMA_LOG_INFO( logger, "compressIA with eps = " << eps )

    #pragma omp parallel
    {
        #pragma omp for

        for( IndexType i = 0; i < numRows; i++ )
        {
            IndexType length = IA[i];

            for( IndexType j = 0; j < IA[i]; j++ )
            {
                if( JA[j * numRows + i] == i )
                {
                    continue;
                }

                if( std::abs( values[j * numRows + i] ) <= eps )
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
void MICELLUtils::compressValues(
    const IndexType IA[],
    const IndexType JA[],
    const ValueType values[],
    const IndexType numRows,
    const ValueType eps,
    IndexType newJA[],
    ValueType newValues[] )
{
    LAMA_LOG_INFO( logger, "compressValues with eps = " << eps )

    #pragma omp parallel
    {
        #pragma omp for

        for( IndexType i = 0; i < numRows; i++ )
        {
            IndexType gap = 0;

            for( IndexType j = 0; j < IA[i]; j++ )
            {
                if( std::abs( values[j * numRows + i] ) <= eps && JA[j * numRows + i] != i )
                {
                    gap++;
                    continue;
                }

                newValues[( j - gap ) * numRows + i] = values[j * numRows + i];
                newJA[( j - gap ) * numRows + i] = JA[j * numRows + i];
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void MICELLUtils::computeIA(
    const IndexType aIA[],
    const IndexType aJA[],
    const IndexType aNumRows,
    const IndexType bIA[],
    const IndexType bJA[],
    const IndexType bNumRows,
    IndexType cIA[] )
{
    LAMA_LOG_INFO( logger, "computeIA with numRows A = " << aNumRows << " and numRows B = " << bNumRows )

    #pragma omp parallel
    {
        #pragma omp for

        for( IndexType i = 0; i < aNumRows; i++ )
        {
            std::set<IndexType> newElements;
            std::pair<std::set<IndexType>::iterator,bool> ret;
            IndexType length = 0;

            for( IndexType j = 0; j < aIA[i]; j++ )
            {
                IndexType jj = aJA[j * aNumRows + i];

                for( IndexType k = 0; k < bIA[jj]; k++ )
                {
                    IndexType kk = bJA[k * bNumRows + jj];
                    ret = newElements.insert( kk );

                    if( ret.second == true )
                    {
                        length++;
                    }
                }
            }

            cIA[i] = length;
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void MICELLUtils::computeValues(
    const IndexType aIA[],
    const IndexType aJA[],
    const ValueType aValues[],
    const IndexType aNumRows,
    const IndexType bIA[],
    const IndexType bJA[],
    const ValueType bValues[],
    const IndexType bNumRows,
    const ValueType alpha,
    const IndexType cIA[],
    IndexType cJA[],
    ValueType cValues[] )
{
    LAMA_LOG_INFO( logger, "computeValues with numRows A = " << aNumRows << " and numRows B = " << bNumRows )

    #pragma omp parallel
    {
        #pragma omp for

        for( IndexType i = 0; i < aNumRows; i++ )
        {
            std::set<IndexType> jaRow;
            std::map<IndexType,ValueType> valuesRow;
            std::pair<std::set<IndexType>::iterator,bool> ret;

            for( IndexType j = 0; j < aIA[i]; j++ )
            {
                IndexType jj = aJA[j * aNumRows + i];

                for( IndexType k = 0; k < bIA[jj]; k++ )
                {
                    IndexType kk = bJA[k * bNumRows + jj];
                    ret = jaRow.insert( kk );
                    ValueType mult = alpha * aValues[j * aNumRows + i] * bValues[k * bNumRows + jj];

                    if( ret.second == true )
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

            for( IndexType j = 0; j < cIA[i]; j++ )
            {
                cJA[j * aNumRows + i] = *jaIter;
                cValues[j * aNumRows + i] = ( *valuesIter ).second;
                jaIter++;
                valuesIter++;
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void MICELLUtils::addComputeIA(
    const IndexType aIA[],
    const IndexType aJA[],
    const IndexType aNumRows,
    const IndexType bIA[],
    const IndexType bJA[],
    const IndexType bNumRows,
    IndexType cIA[] )
{
    LAMA_LOG_INFO( logger, "addComputeIA with numRows A = " << aNumRows << " and numRows B = " << bNumRows )

    #pragma omp parallel
    {
        #pragma omp for

        for( IndexType i = 0; i < aNumRows; i++ )
        {
            std::set<IndexType> iaRow;
            std::pair<std::set<IndexType>::iterator,bool> ret;
            IndexType length = 0;

            for( IndexType j = 0; j < aIA[i]; j++ )
            {
                iaRow.insert( aJA[j * aNumRows + i] );
                length++;
            }

            for( IndexType j = 0; j < bIA[i]; j++ )
            {
                ret = iaRow.insert( bJA[j * bNumRows + i] );

                if( ret.second == true )
                {
                    length++;
                }
            }

            cIA[i] = length;
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void MICELLUtils::addComputeValues(
    const IndexType aIA[],
    const IndexType aJA[],
    const ValueType aValues[],
    const IndexType aNumRows,
    const IndexType bIA[],
    const IndexType bJA[],
    const ValueType bValues[],
    const IndexType bNumRows,
    const ValueType beta,
    const IndexType cIA[],
    IndexType cJA[],
    ValueType cValues[] )
{
    LAMA_LOG_INFO( logger, "computeValues with numRows A = " << aNumRows << " and numRows B = " << bNumRows )

    #pragma omp parallel
    {
        #pragma omp for

        for( IndexType i = 0; i < aNumRows; i++ )
        {
            std::set<IndexType> jaRow;
            std::map<IndexType,ValueType> valuesRow;
            std::pair<std::set<IndexType>::iterator,bool> ret;

            for( IndexType j = 0; j < aIA[i]; j++ )
            {
                IndexType jj = j * aNumRows + i;
                jaRow.insert( aJA[jj] );
                valuesRow.insert( std::pair<IndexType,ValueType>( aJA[jj], aValues[jj] ) );
            }

            for( IndexType j = 0; j < bIA[i]; j++ )
            {
                IndexType jj = j * bNumRows + i;
                ret = jaRow.insert( bJA[jj] );

                if( ret.second == true )
                {
                    valuesRow.insert( std::pair<IndexType,ValueType>( bJA[jj], beta * bValues[jj] ) );
                }
                else
                {
                    valuesRow[bJA[jj]] += beta * bValues[jj];
                }
            }

            std::set<IndexType>::iterator jaIter;
            typename std::map<IndexType,ValueType>::iterator valuesIter;

            jaIter = jaRow.begin();
            valuesIter = valuesRow.begin();

            for( IndexType j = 0; j < cIA[i]; j++ )
            {
                cJA[j * aNumRows + i] = *jaIter;
                cValues[j * aNumRows + i] = ( *valuesIter ).second;
                jaIter++;
                valuesIter++;
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void MICELLUtils::jacobi(
    ValueType solution[],
    const IndexType numRows,
    const IndexType UNUSED( ellNumValuesPerRow ),
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega,
    class SyncToken* syncToken )
{
    LAMA_REGION( "MIC.ELL.jacobi" )

    LAMA_LOG_INFO( logger,
                   "jacobi<" << common::getScalarType<ValueType>() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    if( syncToken )
    {
        MICSyncToken* micSyncToken = dynamic_cast<MICSyncToken*>( syncToken );
        LAMA_ASSERT_ERROR( micSyncToken, "no MIC sync token provided" )
        LAMA_LOG_WARN( logger, "asynchronous execution on MIC not supported yet" )
    }

    void* solutionPtr = solution;
    const void* oldSolutionPtr = oldSolution;
    const void* rhsPtr = rhs;
    const void* ellSizesPtr = ellSizes;
    const void* ellJAPtr = ellJA;
    const void* ellValuesPtr = ellValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( solutionPtr, oldSolutionPtr, rhsPtr, ellSizesPtr, ellJAPtr, ellValuesPtr, omega, numRows )
    {
        ValueType* solution = static_cast<ValueType*>( solutionPtr );
        const ValueType* oldSolution = static_cast<const ValueType*>( oldSolutionPtr );
        const ValueType* rhs = static_cast<const ValueType*>( rhsPtr );
        const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );
        const IndexType* ellJA = static_cast<const IndexType*>( ellJAPtr );
        const ValueType* ellValues = static_cast<const ValueType*>( ellValuesPtr );

        const ValueType oneMinusOmega = static_cast<ValueType>( 1 ) - omega;

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; i++ )
        {
            ValueType temp = rhs[i];
            IndexType pos = i;
            ValueType diag = ellValues[pos]; //getDiagonal

            for( IndexType j = 1; j < ellSizes[i]; j++ )
            {
                pos += numRows;
                temp -= ellValues[pos] * oldSolution[ellJA[pos]];
            }

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

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void MICELLUtils::jacobiHalo(
    ValueType solution[],
    const IndexType numRows,
    const ValueType diagonal[],
    const IndexType UNUSED( ellNumValuesPerRow ),
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    const IndexType rowIndexes[],
    const IndexType numNonEmptyRows,
    const ValueType oldSolution[],
    const ValueType omega,
    class SyncToken* syncToken )
{
    if( syncToken )
    {
        MICSyncToken* micSyncToken = dynamic_cast<MICSyncToken*>( syncToken );
        LAMA_ASSERT_ERROR( micSyncToken, "no MIC sync token provided" )
        // not yet implemented: run the offload computation asynchronously
    }

    LAMA_REGION( "MIC.ELL.jacobiHalo" )

    void* solutionPtr = solution;

    const void* oldSolutionPtr = oldSolution;
    const void* diagonalPtr = diagonal;
    const void* rowIndexesPtr = rowIndexes;
    const void* ellSizesPtr = ellSizes;
    const void* ellJAPtr = ellJA;
    const void* ellValuesPtr = ellValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( solutionPtr, oldSolutionPtr, diagonalPtr, rowIndexesPtr, \
                                            ellSizesPtr, ellJAPtr, ellValuesPtr, omega, numRows, numNonEmptyRows )
    {
        ValueType* solution = static_cast<ValueType*>( solutionPtr );

        const ValueType* oldSolution = static_cast<const ValueType*>( oldSolutionPtr );
        const ValueType* diagonal = static_cast<const ValueType*>( diagonalPtr );
        const IndexType* rowIndexes = static_cast<const IndexType*>( rowIndexesPtr );
        const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );
        const IndexType* ellJA = static_cast<const IndexType*>( ellJAPtr );
        const ValueType* ellValues = static_cast<const ValueType*>( ellValuesPtr );

        #pragma omp parallel for

        for( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            IndexType i = ii; // rowIndexes == NULL stands for all rows

            if( rowIndexes )
            {
                i = rowIndexes[ii];
            }

            ValueType temp = 0.0;
            IndexType pos = i; // index for ellValues

            for( IndexType jj = 0; jj < ellSizes[i]; jj++ )
            {
                temp += ellValues[pos] * oldSolution[ellJA[pos]];
                pos += numRows; // column-major-order
            }

            const ValueType diag = diagonal[i];

            solution[i] -= temp * ( omega / diag );
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void MICELLUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numNonZerosPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "normalGEMV<" << common::getScalarType<ValueType>() << ">, result[" << numRows << "] = " << alpha << " * A( ell, #maxNZ/row = " << numNonZerosPerRow << " ) * x + " << beta << " * y " )

    if( numNonZerosPerRow == 0 )
    {
        COMMON_THROWEXCEPTION( "normalGEMV should not have been called, no entries" )

        // only compute: result = beta * y
    }

    if( syncToken )
    {
        MICSyncToken* micSyncToken = dynamic_cast<MICSyncToken*>( syncToken );
        LAMA_ASSERT_ERROR( micSyncToken, "no MIC sync token provided" )
        // not yet implemented: run the offload computation asynchronously
    }

    LAMA_REGION( "MIC.ELL.normalGEMV" )

    // conversion of pointer to size_t to cheat offload

    void* resultPtr = result;
    const void* xPtr = x;
    const void* yPtr = y;
    const void* ellSizesPtr = ellSizes;
    const void* ellJAPtr = ellJA;
    const void* ellValuesPtr = ellValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( resultPtr, xPtr, yPtr, ellSizesPtr, ellJAPtr, ellValuesPtr, alpha, numRows )
    {
        ValueType* result = static_cast<ValueType*>( resultPtr );

        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const ValueType* y = static_cast<const ValueType*>( yPtr );
        const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );
        const IndexType* ellJA = static_cast<const IndexType*>( ellJAPtr );
        const ValueType* ellValues = static_cast<const ValueType*>( ellValuesPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; ++i )
        {
            ValueType temp = 0;

            for( IndexType jj = 0; jj < ellSizes[i]; ++jj )
            {
                IndexType j = ellJA[i + jj * numRows];
                temp += ellValues[i + jj * numRows] * x[j];
            }

            if( 0 == beta )
            {
                // must be handled separately as y[i] might be uninitialized

                result[i] = alpha * temp;
            }
            else if( 1 == alpha )
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
void MICELLUtils::sparseGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const IndexType numRows,
    const IndexType UNUSED(numNonZerosPerRow),
    const IndexType numNonZeroRows,
    const IndexType rowIndexes[],
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    SyncToken* syncToken )
{
    if( syncToken )
    {
        MICSyncToken* micSyncToken = dynamic_cast<MICSyncToken*>( syncToken );
        LAMA_ASSERT_ERROR( micSyncToken, "no MIC sync token provided" )
        // not yet implemented: run the offload computation asynchronously
    }

    LAMA_REGION( "MIC.ELL.sparseGEMV" )

    LAMA_LOG_INFO( logger,
                   "sparseGEMV<" << common::getScalarType<ValueType>() << ">, n = " << numRows << ", nonZeroRows = " << numNonZeroRows << ", alpha = " << alpha )

    // conversion of pointer to void* to cheat offload

    void* resultPtr = result;

    const void* xPtr = x;
    const void* rowIndexesPtr = rowIndexes;
    const void* ellSizesPtr = ellSizes;
    const void* ellJAPtr = ellJA;
    const void* ellValuesPtr = ellValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( resultPtr, xPtr, rowIndexesPtr, ellSizesPtr, ellJAPtr, ellValuesPtr, \
                                                alpha, numRows, numNonZeroRows )
    {
        ValueType* result = static_cast<ValueType*>( resultPtr );

        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const IndexType* rowIndexes = static_cast<const IndexType*>( rowIndexesPtr );
        const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );
        const IndexType* ellJA = static_cast<const IndexType*>( ellJAPtr );
        const ValueType* ellValues = static_cast<const ValueType*>( ellValuesPtr );

        #pragma omp parallel
        {
            #pragma omp for

            for( IndexType ii = 0; ii < numNonZeroRows; ++ii )
            {
                IndexType i = rowIndexes[ii];

                //result is not initialized for performance reasons
                ValueType temp = 0.0;

                for( IndexType jj = 0; jj < ellSizes[i]; ++jj )
                {
                    IndexType pos = i + jj * numRows;
                    IndexType j = ellJA[pos];
                    temp += ellValues[pos] * x[j];
                }

                if( 1 == alpha )
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

void MICELLUtils::setInterface( ELLUtilsInterface& ELLUtils )
{
    LAMA_LOG_INFO( logger, "set ELL routines for MIC in Interface" )

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

    /*
     LAMA_INTERFACE_REGISTER_T( ELLUtils, absMaxVal, float )
     LAMA_INTERFACE_REGISTER_T( ELLUtils, absMaxVal, double )

     LAMA_INTERFACE_REGISTER_T( ELLUtils, compressIA, float )
     LAMA_INTERFACE_REGISTER_T( ELLUtils, compressIA, double )

     LAMA_INTERFACE_REGISTER_T( ELLUtils, compressValues, float )
     LAMA_INTERFACE_REGISTER_T( ELLUtils, compressValues, double )

     LAMA_INTERFACE_REGISTER_T( ELLUtils, computeIA, IndexType )
     LAMA_INTERFACE_REGISTER_T( ELLUtils, addComputeIA, IndexType )

     LAMA_INTERFACE_REGISTER_T( ELLUtils, computeValues, float )
     LAMA_INTERFACE_REGISTER_T( ELLUtils, computeValues, double )

     LAMA_INTERFACE_REGISTER_T( ELLUtils, addComputeValues, float )
     LAMA_INTERFACE_REGISTER_T( ELLUtils, addComputeValues, double )
     */

    LAMA_INTERFACE_REGISTER_T( ELLUtils, normalGEMV, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, normalGEMV, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, sparseGEMV, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, sparseGEMV, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, jacobi, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, jacobi, double )

    LAMA_INTERFACE_REGISTER_T( ELLUtils, jacobiHalo, float )
    LAMA_INTERFACE_REGISTER_T( ELLUtils, jacobiHalo, double )

    /*
     LAMA_INTERFACE_REGISTER_T( ELLUtils, fillELLValues, float )
     LAMA_INTERFACE_REGISTER_T( ELLUtils, fillELLValues, double )
     */
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the ELLUtils routines                             */
/* --------------------------------------------------------------------------- */

bool MICELLUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::MIC );
    setInterface( interface.ELLUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool MICELLUtils::initialized = registerInterface();

} // namespace lama
