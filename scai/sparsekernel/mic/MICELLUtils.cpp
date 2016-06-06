/**
 * @file MICELLUtils.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Implementation of ELL utilities with MIC
 * @author Thomas Brandes
 * @date 04.07.2012
 */

// hpp
#include <scai/sparsekernel/mic/MICELLUtils.hpp>

// local project
#include <scai/sparsekernel/ELLKernelTrait.hpp>

// internal scai projects
#include <scai/tasking/mic/MICSyncToken.hpp>
#include <scai/hmemo/mic/MICContext.hpp>
#include <scai/tasking/NoSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/OpenMP.hpp>
#include <scai/common/bind.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>

// external

// std
#include <set>
#include <map>
#include <cmath>

namespace scai
{

using tasking::SyncToken;
using tasking::MICSyncToken;

using namespace hmemo;

namespace sparsekernel
{

/* ------------------------------------------------------------------------------------------------------------------ */

SCAI_LOG_DEF_LOGGER( MICELLUtils::logger, "MIC.ELLUtils" )

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

    SCAI_LOG_INFO( logger, "#non-zero rows = " << counter << ", counted by sizes" )

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

    SCAI_ASSERT_EQUAL_DEBUG( counter, numNonEmptyRows )

    SCAI_LOG_INFO( logger, "#non-zero rows = " << counter << ", set by sizes" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

bool MICELLUtils::hasDiagonalProperty( const IndexType numDiagonals, const IndexType ellJA[] )
{
    SCAI_LOG_INFO( logger, "hasDiagonalProperty, #numDiagonals = " << numDiagonals )

    if( numDiagonals == 0 )
    {
        return false;
    }

    bool diagonalProperty = true;

    const void* ellJAPtr = ellJA;

    int device = MICContext::getCurrentDevice();

#pragma offload target( MIC : device ) in( ellJAPtr, numDiagonals ), out( diagonalProperty )
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
    SCAI_LOG_INFO( logger,
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
    SCAI_LOG_INFO( logger,
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

        SCAI_ASSERT_ERROR( integrityIA, msg << ": ellSizes: at least one value out of range" );
        SCAI_ASSERT_ERROR( integrityJA, msg << ": ellJA: at least one value out of range" );
    }
    else
    {
        SCAI_ASSERT_ERROR( numValuesPerRow == 0,
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
    ValueType maxValue = static_cast<ValueType>(0.0);

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>(0.0);

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

                // SCAI_LOG_TRACE( logger, "absMaxVal, val[" << i << ", " << jj << "] = " << val )
            }
        }

        #pragma omp critical
        {
            SCAI_LOG_DEBUG( logger, "absMaxVal, threadVal = " << threadVal << ", maxVal = " << maxValue )

            if( threadVal > maxValue )
            {
                maxValue = threadVal;
            }
        }
    }

    SCAI_LOG_DEBUG( logger, "absMaxVal, maxVal = " << maxValue )

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
    SCAI_LOG_DEBUG( logger, "get row #i = " << i )

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
            row[j] = static_cast<OtherValueType>(0.0);
        }

        #pragma omp parallel for

        for( IndexType jj = 0; jj < ellSizes[i]; ++jj )
        {
            IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );
            row[ellJA[pos]] = static_cast<OtherValueType>( ellValues[pos] );
        }
    }
}

template<typename ValueType>
ValueType MICELLUtils::getValue(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType numValuesPerRow,
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[] )
{
    SCAI_LOG_TRACE( logger, "get value i = " << i << ", j = " << j )

    ValueType value = 0;  // not really needed, just for safety

    const void* ellSizesPtr = ellSizes;
    const void* ellJAPtr = ellJA;
    const void* ellValuesPtr = ellValues;
    ValueType* valuePtr = &value;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), out( valuePtr[0:1] ), in( ellSizesPtr, ellJAPtr, ellValuesPtr, \
                                                              i, j, numRows, numValuesPerRow )
    {
        const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );
        const IndexType* ellJA = static_cast<const IndexType*>( ellJAPtr );
        const ValueType* ellValues = static_cast<const ValueType*>( ellValuesPtr );
        ValueType& valueRef = *valuePtr;

        valueRef = 0;  // new initialiation, has not been copied in

        for( IndexType jj = 0; jj < ellSizes[i]; ++jj )
        {
            IndexType pos = ellindex( i, jj, numRows, numValuesPerRow );

            if( ellJA[pos] == j )
            {
                valueRef = ellValues[pos];
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
    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << common::getScalarType<ELLValueType>() << ", " << common::getScalarType<CSRValueType>() << ">" << ", #rows = " << numRows )

    // SCAI_REGION( "MIC.ELL->CSR_values" )

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
    // SCAI_REGION( "MIC.ELL<-CSR_values" )

    SCAI_LOG_INFO( logger,
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
    SCAI_LOG_INFO( logger, "fill ELLValues<" << typeid( ValueType ).name() )

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
                ellValues[pos] = static_cast<ValueType>(0.0); // zero entry
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
    SCAI_LOG_INFO( logger, "compressIA with eps = " << eps )

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
    SCAI_LOG_INFO( logger, "compressValues with eps = " << eps )

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
    SCAI_LOG_INFO( logger, "computeIA with numRows A = " << aNumRows << " and numRows B = " << bNumRows )

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
    SCAI_LOG_INFO( logger, "computeValues with numRows A = " << aNumRows << " and numRows B = " << bNumRows )

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
    SCAI_LOG_INFO( logger, "addComputeIA with numRows A = " << aNumRows << " and numRows B = " << bNumRows )

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
    SCAI_LOG_INFO( logger, "computeValues with numRows A = " << aNumRows << " and numRows B = " << bNumRows )

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
    const IndexType SCAI_UNUSED( ellNumValuesPerRow ),
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega )
{
    // SCAI_REGION( "MIC.ELL.jacobi" )

    SCAI_LOG_INFO( logger,
                   "jacobi<" << common::getScalarType<ValueType>() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
    }

    void* solutionPtr = solution;
    const void* oldSolutionPtr = oldSolution;
    const void* rhsPtr = rhs;
    const void* ellSizesPtr = ellSizes;
    const void* ellJAPtr = ellJA;
    const void* ellValuesPtr = ellValues;

    const ValueType* omegaPtr = &omega;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( solutionPtr, oldSolutionPtr, rhsPtr, ellSizesPtr, ellJAPtr, ellValuesPtr, omegaPtr[0:1], numRows )
    {
        ValueType* solution = static_cast<ValueType*>( solutionPtr );
        const ValueType* oldSolution = static_cast<const ValueType*>( oldSolutionPtr );
        const ValueType* rhs = static_cast<const ValueType*>( rhsPtr );
        const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );
        const IndexType* ellJA = static_cast<const IndexType*>( ellJAPtr );
        const ValueType* ellValues = static_cast<const ValueType*>( ellValuesPtr );

        const ValueType& omegaRef = *omegaPtr;

        const ValueType oneMinusOmega = static_cast<ValueType>(1.0) - omegaRef;

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

            if( omegaRef == static_cast<ValueType>( 1.0 ) )
            {
                solution[i] = temp / diag;
            }
            else if( omegaRef == 0.5 )
            {
                solution[i] = omegaRef * ( temp / diag + oldSolution[i] );
            }
            else
            {
                solution[i] = omegaRef * ( temp / diag ) + oneMinusOmega * oldSolution[i];
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
    const IndexType SCAI_UNUSED( ellNumValuesPerRow ),
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[],
    const IndexType rowIndexes[],
    const IndexType numNonEmptyRows,
    const ValueType oldSolution[],
    const ValueType omega )
{
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
    }

    // SCAI_REGION( "MIC.ELL.jacobiHalo" )

    void* solutionPtr = solution;

    const void* oldSolutionPtr = oldSolution;
    const void* diagonalPtr = diagonal;
    const void* rowIndexesPtr = rowIndexes;
    const void* ellSizesPtr = ellSizes;
    const void* ellJAPtr = ellJA;
    const void* ellValuesPtr = ellValues;
    const ValueType* omegaPtr = &omega;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( solutionPtr, oldSolutionPtr, diagonalPtr, rowIndexesPtr, \
                                            ellSizesPtr, ellJAPtr, ellValuesPtr, omegaPtr[0:1], numRows, numNonEmptyRows )
    {
        ValueType* solution = static_cast<ValueType*>( solutionPtr );

        const ValueType* oldSolution = static_cast<const ValueType*>( oldSolutionPtr );
        const ValueType* diagonal = static_cast<const ValueType*>( diagonalPtr );
        const IndexType* rowIndexes = static_cast<const IndexType*>( rowIndexesPtr );
        const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );
        const IndexType* ellJA = static_cast<const IndexType*>( ellJAPtr );
        const ValueType* ellValues = static_cast<const ValueType*>( ellValuesPtr );

        const ValueType& omegaRef = *omegaPtr;

        #pragma omp parallel for

        for( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            IndexType i = ii; // rowIndexes == NULL stands for all rows

            if( rowIndexes )
            {
                i = rowIndexes[ii];
            }

            ValueType temp = static_cast<ValueType>(0.0);
            IndexType pos = i; // index for ellValues

            for( IndexType jj = 0; jj < ellSizes[i]; jj++ )
            {
                temp += ellValues[pos] * oldSolution[ellJA[pos]];
                pos += numRows; // column-major-order
            }

            const ValueType diag = diagonal[i];

            solution[i] -= temp * ( omegaRef / diag );
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
    const ValueType ellValues[] )
{
    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << common::getScalarType<ValueType>() << ">, result[" << numRows << "] = " << alpha << " * A( ell, #maxNZ/row = " << numNonZerosPerRow << " ) * x + " << beta << " * y " )

    if( numNonZerosPerRow == 0 )
    {
        COMMON_THROWEXCEPTION( "normalGEMV should not have been called, no entries" )

        // only compute: result = beta * y
    }

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution of JDS jacobi iteration for MIC not supported yet" )
    }

    // SCAI_REGION( "MIC.ELL.normalGEMV" )

    // conversion of pointer to size_t to cheat offload

    void* resultPtr = result;
    const void* xPtr = x;
    const void* yPtr = y;
    const void* ellSizesPtr = ellSizes;
    const void* ellJAPtr = ellJA;
    const void* ellValuesPtr = ellValues;
    const ValueType* alphaPtr = &alpha;
    const ValueType* betaPtr = &beta;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( resultPtr, xPtr, yPtr, ellSizesPtr, ellJAPtr, ellValuesPtr, alphaPtr[0:1], betaPtr[0:1], numRows )
    {
        ValueType* result = static_cast<ValueType*>( resultPtr );

        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const ValueType* y = static_cast<const ValueType*>( yPtr );
        const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );
        const IndexType* ellJA = static_cast<const IndexType*>( ellJAPtr );
        const ValueType* ellValues = static_cast<const ValueType*>( ellValuesPtr );

        const ValueType& alphaRef = *alphaPtr;
        const ValueType& betaRef = *betaPtr;

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; ++i )
        {
            ValueType temp = static_cast<ValueType>(0.0);

            for( IndexType jj = 0; jj < ellSizes[i]; ++jj )
            {
                IndexType j = ellJA[i + jj * numRows];
                temp += ellValues[i + jj * numRows] * x[j];
            }

            if( 0 == betaRef )
            {
                // must be handled separately as y[i] might be uninitialized

                result[i] = alphaRef * temp;
            }
            else if( alphaRef == static_cast<ValueType>( 1.0 ) )
            {
                result[i] = temp + betaRef * y[i];
            }
            else
            {
                result[i] = alphaRef * temp + betaRef * y[i];
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
    const IndexType SCAI_UNUSED(numNonZerosPerRow),
    const IndexType numNonZeroRows,
    const IndexType rowIndexes[],
    const IndexType ellSizes[],
    const IndexType ellJA[],
    const ValueType ellValues[] )
{
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution of JDS jacobi iteration for MIC not supported yet" )
    }

    // SCAI_REGION( "MIC.ELL.sparseGEMV" )

    SCAI_LOG_INFO( logger,
                   "sparseGEMV<" << common::getScalarType<ValueType>() << ">, n = " << numRows << ", nonZeroRows = " << numNonZeroRows << ", alpha = " << alpha )

    // conversion of pointer to void* to cheat offload

    void* resultPtr = result;

    const void* xPtr = x;
    const void* rowIndexesPtr = rowIndexes;
    const void* ellSizesPtr = ellSizes;
    const void* ellJAPtr = ellJA;
    const void* ellValuesPtr = ellValues;

    const ValueType* alphaPtr = &alpha;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( resultPtr, xPtr, rowIndexesPtr, ellSizesPtr, ellJAPtr, ellValuesPtr, \
                                                alphaPtr[0:1], numRows, numNonZeroRows )
    {
        ValueType* result = static_cast<ValueType*>( resultPtr );

        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const IndexType* rowIndexes = static_cast<const IndexType*>( rowIndexesPtr );
        const IndexType* ellSizes = static_cast<const IndexType*>( ellSizesPtr );
        const IndexType* ellJA = static_cast<const IndexType*>( ellJAPtr );
        const ValueType* ellValues = static_cast<const ValueType*>( ellValuesPtr );

        const ValueType& alphaRef = *alphaPtr;

        #pragma omp parallel
        {
            #pragma omp for

            for( IndexType ii = 0; ii < numNonZeroRows; ++ii )
            {
                IndexType i = rowIndexes[ii];

                //result is not initialized for performance reasons
                ValueType temp = static_cast<ValueType>(0.0);

                for( IndexType jj = 0; jj < ellSizes[i]; ++jj )
                {
                    IndexType pos = i + jj * numRows;
                    IndexType j = ellJA[pos];
                    temp += ellValues[pos] * x[j];
                }

                if( alphaRef == static_cast<ValueType>( 1.0 ) )
                {
                    result[i] += temp;
                }
                else
                {
                    result[i] += alphaRef * temp;
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

void MICELLUtils::Registrator::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_INFO( logger, "register ELLUtils OpenMP-routines for MIC at kernel registry [" << flag << "]" )

    KernelRegistry::set<ELLKernelTrait::countNonEmptyRowsBySizes>( countNonEmptyRowsBySizes, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::setNonEmptyRowsBySizes>( setNonEmptyRowsBySizes, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::hasDiagonalProperty>( hasDiagonalProperty, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::check>( check, ctx, flag );
}

template<typename ValueType>
void MICELLUtils::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_INFO( logger, "register ELLUtils OpenMP-routines for MIC at kernel registry [" << flag << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<ELLKernelTrait::getValue<ValueType> >( getValue, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::sparseGEMV<ValueType> >( sparseGEMV, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void MICELLUtils::RegistratorVO<ValueType, OtherValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_INFO( logger, "register ELLUtils OpenMP-routines for MIC at kernel registry [" << flag
                   << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )

    KernelRegistry::set<ELLKernelTrait::getRow<ValueType, OtherValueType> >( getRow, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::scaleValue<ValueType, OtherValueType> >( scaleValue, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::setCSRValues<ValueType, OtherValueType> >( setCSRValues, ctx, flag );
    KernelRegistry::set<ELLKernelTrait::getCSRValues<ValueType, OtherValueType> >( getCSRValues, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICELLUtils::RegisterGuard::RegisterGuard()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;

    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_MIC_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARITHMETIC_MIC_LIST, SCAI_ARITHMETIC_MIC_LIST>::call( flag );
}

MICELLUtils::RegisterGuard::~RegisterGuard()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_MIC_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARITHMETIC_MIC_LIST, SCAI_ARITHMETIC_MIC_LIST>::call( flag );
}

MICELLUtils::RegisterGuard MICELLUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */
