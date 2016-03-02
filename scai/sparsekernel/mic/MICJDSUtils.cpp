/**
 * @file MICJDSUtils.cpp
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
 * @brief Implementation of JDS utilities on MIC platform using offload + OpenMP
 * @author Thomas Brandes
 * @date 04.07.2012
 * @since 1.1.0
 */

// hpp
#include <scai/sparsekernel/mic/MICJDSUtils.hpp>

// local project
#include <scai/sparsekernel/JDSKernelTrait.hpp>

// other scai libraries
#include <scai/utilskernel/mic/MICUtils.hpp>
#include <scai/hmemo/mic/MICContext.hpp>
#include <scai/hmemo/mic/MICSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/mepr/Container.hpp>

#include <scai/tracing.hpp>

namespace scai
{

using utilskernel::MICUtils;

using tasking::SyncToken;
using tasking::MICSyncToken;

using common::TypeTraits;

using namespace hmemo;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( MICJDSUtils::logger, "MIC.JDSUtils" )

/* --------------------------------------------------------------------------- */
/*   Implementation of methods                                                 */
/* --------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void MICJDSUtils::getRow(
    OtherValueType row[],
    const IndexType i,
    const IndexType numColumns,
    const IndexType numRows,
    const IndexType perm[],
    const IndexType ilg[],
    const IndexType dlg[],
    const IndexType ja[],
    const ValueType values[] )
{
    SCAI_LOG_INFO( logger, "getRow with i = " << i << ", numColumns = " << numColumns << " and numRows = " << numRows )

    int device = MICContext::getCurrentDevice();

    void* rowPtr = row;
    const void* permPtr = perm;
    const void* ilgPtr = ilg;
    const void* dlgPtr = dlg;
    const void* jaPtr = ja;
    const void* valuesPtr = values;

#pragma offload target( mic : device ) in( rowPtr, permPtr, ilgPtr, dlgPtr, jaPtr, valuesPtr, \
                                               i, numColumns, numRows )
    {
        OtherValueType* row = static_cast<OtherValueType*>( rowPtr );

        const IndexType* perm = static_cast<const IndexType*>( permPtr );
        const IndexType* ilg = static_cast<const IndexType*>( ilgPtr );
        const IndexType* dlg = static_cast<const IndexType*>( dlgPtr );
        const IndexType* ja = static_cast<const IndexType*>( jaPtr );
        const ValueType* values = static_cast<const ValueType*>( valuesPtr );

        for( IndexType j = 0; j < numColumns; ++j )
        {
            row[j] = static_cast<OtherValueType>(0.0);
        }

        // one thread will set the row, but every thread searches

        #pragma omp parallel for

        for( IndexType ii = 0; ii < numRows; ii++ )
        {
            if( perm[ii] == i )
            {
                IndexType k = 0;

                for( IndexType jj = 0; jj < ilg[ii]; ++jj )
                {
                    row[ja[ii + k]] = static_cast<OtherValueType>( values[ii + k] );
                    k += dlg[jj];
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
ValueType MICJDSUtils::getValue(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType* dlg,
    const IndexType* ilg,
    const IndexType* perm,
    const IndexType* ja,
    const ValueType* values )
{
    int device = MICContext::getCurrentDevice();

    const void* permPtr = perm;
    const void* ilgPtr = ilg;
    const void* dlgPtr = dlg;
    const void* jaPtr = ja;
    const void* valuesPtr = values;

    ValueType val = 0;

#pragma offload target( mic : device ) in( permPtr, ilgPtr, dlgPtr, jaPtr, valuesPtr, \
                                               i, j, numRows ), out( val )
    {
        val = 0;

        const IndexType* perm = static_cast<const IndexType*>( permPtr );
        const IndexType* ilg = static_cast<const IndexType*>( ilgPtr );
        const IndexType* dlg = static_cast<const IndexType*>( dlgPtr );
        const IndexType* ja = static_cast<const IndexType*>( jaPtr );
        const ValueType* values = static_cast<const ValueType*>( valuesPtr );

        #pragma omp parallel for

        for( IndexType ii = 0; ii < numRows; ii++ )
        {
            if( perm[ii] == i )
            {
                // thread founds the row and now searches the column

                IndexType k = 0;

                for( IndexType jj = 0; jj < ilg[ii]; jj++ )
                {
                    if( ja[ii + k] == j )
                    {
                        val = values[ii + k];
                    }

                    k += dlg[jj];
                }
            }
        }
    }

    return val;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void MICJDSUtils::scaleValue(
    const IndexType numRows,
    const IndexType jdsPerm[],
    const IndexType ilg[],
    const IndexType dlg[],
    ValueType jdsValues[],
    const OtherValueType values[] )
{
    SCAI_LOG_INFO( logger, "scaleValue with numRows = " << numRows )

    int device = MICContext::getCurrentDevice();

    const void* jdsPermPtr = jdsPerm;
    const void* ilgPtr = ilg;
    const void* dlgPtr = dlg;
    void* jdsValuesPtr = jdsValues;
    const void* valuesPtr = values;

#pragma offload target( mic : device ) in( jdsPermPtr, ilgPtr, dlgPtr, jdsValuesPtr, \
                                               valuesPtr, numRows ),
    {
        const IndexType* jdsPerm = static_cast<const IndexType*>( jdsPermPtr );
        const IndexType* ilg = static_cast<const IndexType*>( ilgPtr );
        const IndexType* dlg = static_cast<const IndexType*>( dlgPtr );

        ValueType* jdsValues = static_cast<ValueType*>( jdsValuesPtr );

        const OtherValueType* values = static_cast<const OtherValueType*>( valuesPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; i++ )
        {
            IndexType offset = i;
            OtherValueType scalar = values[jdsPerm[i]];

            for( IndexType jj = 0; jj < ilg[i]; jj++ )
            {
                jdsValues[offset] *= static_cast<ValueType>( scalar );
                offset += dlg[jj];
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

bool MICJDSUtils::checkDiagonalProperty(
    const IndexType numDiagonals,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType jdsPerm[],
    const IndexType ja[],
    const IndexType dlg[] )
{
    SCAI_LOG_INFO( logger,
                   "checkDiagonalProperty with numDiagonals = " << numDiagonals << ", numColumns = " << numColumns << " and numRows = " << numRows )

    if( numRows > 0 )
    {
        // offload dlg[0]

        IndexType dlg0 = MICUtils::getValue( dlg, 0 );

        if( dlg0 < std::min( numDiagonals, numColumns ) )
        {
            // not even one entry for each row / column
            return false;
        }

        bool diagonalProperty = true;

        void* jdsPermPtr = (void*) jdsPerm;
        void* jaPtr = (void*) ja;

        int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( jdsPermPtr, jaPtr, numRows, numColumns, dlg0 ), out( diagonalProperty )
        {
            const IndexType* ja = static_cast<const IndexType*>( jaPtr );
            const IndexType* jdsPerm = static_cast<const IndexType*>( jdsPermPtr );

            diagonalProperty = true;

            #pragma omp parallel for

            for( IndexType ii = 0; ii < numRows; ++ii )
            {
                if( !diagonalProperty )
                {
                    continue;
                }

                const IndexType i = jdsPerm[ii];

                if( i >= numColumns )
                {
                    continue;
                }

                if( ii >= dlg0 )
                {
                    // ilg[ii] = 0, empty row

                    diagonalProperty = false;
                }
                else if( ja[ii] != i )
                {
                    diagonalProperty = false;
                }
            }
        }

        return diagonalProperty;
    }

    return false;
}

/* ------------------------------------------------------------------------------------------------------------------ */

void MICJDSUtils::setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n )
{
    SCAI_LOG_INFO( logger, "compute inverse perm, n = " << n )

    void* inversePermPtr = inversePerm;
    const void* permPtr = perm;

    // Parallel execution is safe as perm does not contain a value twice

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( permPtr, inversePermPtr, n )
    {
        IndexType* inversePerm = static_cast<IndexType*>( inversePermPtr );
        const IndexType* perm = static_cast<const IndexType*>( permPtr );

        #pragma omp parallel for

        for( IndexType ii = 0; ii < n; ii++ )
        {
            IndexType i = perm[ii];
            inversePerm[i] = ii;
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

void MICJDSUtils::sortRows( IndexType ilg[], IndexType perm[], const IndexType n )
{
    void* ilgPtr = ilg;
    void* permPtr = perm;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( ilgPtr, permPtr, n )
    {
        IndexType* ilg = static_cast<IndexType*>( ilgPtr );
        IndexType* perm = static_cast<IndexType*>( permPtr );

        // Help array needed, because bucket sort cannot be done in-place

        IndexType* input = new IndexType[n];

        // Open: can this routine be called where perm is a valid permutation as input

        #pragma omp parallel for

        for( IndexType i = 0; i < n; i++ )
        {
            input[i] = perm[i];
        }

        // The number of buckets is determined by the max value of ilg

        IndexType maxBucket = 0;

        #pragma omp parallel
        {
            IndexType threadMax = 0;
            #pragma omp for

            for( IndexType i = 0; i < n; ++i )
            {
                if( ilg[i] > threadMax )
                {
                    threadMax = ilg[i];
                }
            }

            #pragma omp critical
            {
                if( threadMax > maxBucket )
                {
                    maxBucket = threadMax;
                }
            }
        }

        // longest row = maxBucket, but rows with length 0 is possible too!

        IndexType* bucket = new IndexType[maxBucket + 1];

        for( IndexType i = 0; i <= maxBucket; i++ )
        {
            bucket[i] = 0;
        }

        // counts how many diagonals exist for each possible length
        for( IndexType i = 0; i < n; i++ )
        {
            bucket[ilg[i]]++;
        }

        // use now bucket array for finding right offsets
        // diag length:                 0   1   2   3   4   5
        // number of (now in bucket):   3   4   3   5   1   5
        // becomes (end of first for): 18  14  11   6   5   0
        // later (end of second for):  21  18  14  11   6   5

        IndexType total = 0;

        for( IndexType i = maxBucket; i >= 0; i-- )
        {
            IndexType cnt = bucket[i];
            bucket[i] = total;
            total += cnt;
        }

        // now we can build the new perm array
        // diagonals with same lengths are moved to position bucket[b] upwards
        for( IndexType i = 0; i < n; i++ )
        {
            IndexType b = ilg[i];
            perm[bucket[b]++] = input[i];
        }

        // reorganize of ilg has to wait until after filling of perm array is finished
        total = 0;

        for( IndexType i = maxBucket; i >= 0; i-- )
        {
            for( IndexType k = total; k < bucket[i]; k++ )
            {
                ilg[k] = i;
            }

            total = bucket[i];
        }

        delete[] bucket;
        delete[] input;
    }
}

/* --------------------------------------------------------------------------- */

IndexType MICJDSUtils::ilg2dlg(
    IndexType dlg[],
    const IndexType numDiagonals,
    const IndexType ilg[],
    const IndexType numRows )
{
    SCAI_LOG_INFO( logger, "ilg2dlg with numDiagonals = " << numDiagonals << ", numRows = " << numRows )

    if( numDiagonals == 0 )
    {
        return 0;
    }

    // Entries in dlg filled every time there is a change in values of consecutive elements of ilg
    //
    //   i:     0  1  2  3  4  5
    // ilg:     5  5  3  3  3  1
    // nd1:     5  5  3  3  3  1
    // nd2:     5  3  3  3  1  0
    //             x        x  x
    //             |        |  |->    6
    //             |        |---->       5  5
    //             |------------->             2   2
    // dlg:                           6  5  5  2   2

    IndexType numTotal = 0;

    void* dlgPtr = dlg;
    const void* ilgPtr = ilg;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( ilgPtr, dlgPtr, numRows ), out( numTotal )
    {
        IndexType* dlg = static_cast<IndexType*>( dlgPtr );
        const IndexType* ilg = static_cast<const IndexType*>( ilgPtr );

        numTotal = 0;

        #pragma omp parallel for reduction( +:numTotal )

        for( IndexType i = 0; i < numRows; ++i )
        {
            IndexType nd1 = ilg[i];
            IndexType nd2 = 0;

            if( i + 1 < numRows )
            {
                nd2 = ilg[i + 1];
            }

            // fill in dlg only if nd2 < nd1

            for( IndexType j = nd2; j < nd1; j++ )
            {
                dlg[j] = i + 1;
            }

            numTotal += nd1;
        }
    }

    return numTotal;
}

/* --------------------------------------------------------------------------- */

template<typename JDSValueType,typename CSRValueType>
void MICJDSUtils::getCSRValues(
    IndexType csrJA[],
    CSRValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType jdsInversePerm[],
    const IndexType jdsILG[],
    const IndexType jdsDLG[],
    const IndexType jdsJA[],
    const JDSValueType jdsValues[] )
{
    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << TypeTraits<JDSValueType>::id() << ", " << TypeTraits<CSRValueType>::id() << ">" 
                    << ", #rows = " << numRows )

    // SCAI_REGION( "MIC.JDS->CSR_values" )

    const void* jdsJAPtr = jdsJA;
    const void* jdsValuesPtr = jdsValues;
    const void* jdsInversePermPtr = jdsInversePerm;
    const void* jdsILGPtr = jdsILG;
    const void* jdsDLGPtr = jdsDLG;

    void* csrJAPtr = csrJA;
    void* csrValuesPtr = csrValues;
    const void* csrIAPtr = csrIA;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) \
                in( numRows, jdsJAPtr, jdsValuesPtr, jdsInversePermPtr, jdsILGPtr, jdsDLGPtr, \
                    csrIAPtr, csrJAPtr, csrValuesPtr )
    {
        const JDSValueType* jdsValues = static_cast<const JDSValueType*>( jdsValuesPtr );
        const IndexType* jdsJA = static_cast<const IndexType*>( jdsJAPtr );
        const IndexType* jdsInversePerm = static_cast<const IndexType*>( jdsInversePermPtr );
        const IndexType* jdsILG = static_cast<const IndexType*>( jdsILGPtr );
        const IndexType* jdsDLG = static_cast<const IndexType*>( jdsDLGPtr );

        const IndexType* csrIA = static_cast<const IndexType*>( csrIAPtr );
        IndexType* csrJA = static_cast<IndexType*>( csrJAPtr );
        CSRValueType* csrValues = static_cast<CSRValueType*>( csrValuesPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; i++ )
        {
            IndexType ii = jdsInversePerm[i]; // where to find row i in JDS storage

            const IndexType numValuesInRow = jdsILG[ii];

            IndexType jdsOffset = ii; // run through input JDS data
            IndexType offset = csrIA[i]; // run through output data

            for( IndexType jj = 0; jj < numValuesInRow; jj++ )
            {
                csrJA[offset + jj] = jdsJA[jdsOffset];
                csrValues[offset + jj] = static_cast<CSRValueType>( jdsValues[jdsOffset] );
                jdsOffset += jdsDLG[jj];
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename JDSValueType,typename CSRValueType>
void MICJDSUtils::setCSRValues(
    IndexType jdsJA[],
    JDSValueType jdsValues[],
    const IndexType numRows,
    const IndexType jdsPerm[],
    const IndexType jdsILG[],
    const IndexType /* ndlg */,
    const IndexType jdsDLG[],
    const IndexType csrIA[],
    const IndexType csrJA[],
    const CSRValueType csrValues[] )
{
    SCAI_LOG_INFO( logger,
                   "set CSRValues<" << TypeTraits<JDSValueType>::id() << ", " << TypeTraits<CSRValueType>::id() << ">"
                    << ", #rows = " << numRows )

    // SCAI_REGION( "MIC.JDS<-CSR_values" )

    void* jdsJAPtr = jdsJA;
    void* jdsValuesPtr = jdsValues;

    const void* jdsPermPtr = jdsPerm;
    const void* jdsILGPtr = jdsILG;
    const void* jdsDLGPtr = jdsDLG;
    const void* csrIAPtr = csrIA;
    const void* csrJAPtr = csrJA;
    const void* csrValuesPtr = csrValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( numRows, jdsJAPtr, jdsValuesPtr, jdsPermPtr, jdsILGPtr, jdsDLGPtr, \
                                           csrIAPtr, csrJAPtr, csrValuesPtr )
    {
        JDSValueType* jdsValues = static_cast<JDSValueType*>( jdsValuesPtr );
        IndexType* jdsJA = static_cast<IndexType*>( jdsJAPtr );

        const IndexType* jdsPerm = static_cast<const IndexType*>( jdsPermPtr );
        const IndexType* jdsILG = static_cast<const IndexType*>( jdsILGPtr );
        const IndexType* jdsDLG = static_cast<const IndexType*>( jdsDLGPtr );

        const IndexType* csrIA = static_cast<const IndexType*>( csrIAPtr );
        const IndexType* csrJA = static_cast<const IndexType*>( csrJAPtr );
        const CSRValueType* csrValues = static_cast<const CSRValueType*>( csrValuesPtr );

        #pragma omp parallel for

        for( IndexType ii = 0; ii < numRows; ii++ )
        {
            IndexType i = jdsPerm[ii];
            IndexType offset = ii;

            for( IndexType jdsJJ = 0, csrJJ = csrIA[i]; jdsJJ < jdsILG[ii]; jdsJJ++, csrJJ++ )
            {
                jdsJA[offset] = csrJA[csrJJ];
                jdsValues[offset] = static_cast<JDSValueType>( csrValues[csrJJ] );
                offset += jdsDLG[jdsJJ]; // index for next value of the row
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICJDSUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType perm[],
    const IndexType jdsILG[],
    const IndexType ndlg,
    const IndexType jdsDLG[],
    const IndexType jdsJA[],
    const ValueType jdsValues[] )
{
    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << TypeTraits<ValueType>::id() << ">, result[" << numRows << "] = " << alpha << " * A( jds, ndlg = " << ndlg << " ) * x + " << beta << " * y " )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
    }

    if( beta == common::constants::ZERO )
    {
        MICUtils::setVal( result, numRows, ValueType( 0 ), common::reduction::COPY );
    }
    else if( result == y )
    {
        // result = result * beta

        if ( beta == common::constants::ONE )
        {
            SCAI_LOG_DEBUG( logger, "result remains unchanged" )
        }
        else
        {
            MICUtils::setVal( result, numRows, beta, common::reduction::MULT );
        }
    }
    else
    {
        // result = beta * y

        MICUtils::setScale( result, beta, y, numRows );
    }

    if( ndlg == 0 )
    {
        return; // definitively empty matrix
    }

    // SCAI_REGION( "MIC.JDS.normalGEMV" )

    void* resultPtr = result;
    const void* xPtr = x;
    const void* permPtr = perm;
    const void* jdsILGPtr = jdsILG;
    const void* jdsDLGPtr = jdsDLG;
    const void* jdsJAPtr = jdsJA;
    const void* jdsValuesPtr = jdsValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( resultPtr, xPtr, jdsDLGPtr, jdsILGPtr, jdsJAPtr, jdsValuesPtr, alpha )
    {
        ValueType* result = static_cast<ValueType*>( resultPtr );
        const ValueType* x = static_cast<const ValueType*>( xPtr );
        const IndexType* perm = static_cast<const IndexType*>( permPtr );
        const IndexType* jdsILG = static_cast<const IndexType*>( jdsILGPtr );
        const IndexType* jdsDLG = static_cast<const IndexType*>( jdsDLGPtr );
        const IndexType* jdsJA = static_cast<const IndexType*>( jdsJAPtr );
        const ValueType* jdsValues = static_cast<const ValueType*>( jdsValuesPtr );

        // dlg[0] stands exactly for number of non-empty rows

        IndexType nonEmptyRows = jdsDLG[0];

        #pragma omp parallel for

        for( IndexType ii = 0; ii < nonEmptyRows; ii++ )
        {
            ValueType value = static_cast<ValueType>(0.0); // sums up final value
            IndexType offset = ii;

            for( IndexType jj = 0; jj < jdsILG[ii]; jj++ )
            {
                IndexType j = jdsJA[offset];
                value += jdsValues[offset] * x[j];
                offset += jdsDLG[jj]; // there is next value for this row
            }

            // scattering needs no synchronization as values of perm are unique

            result[perm[ii]] += alpha * value;
        }
    }
}

/* --------------------------------------------------------------------------- */
/*     Jacobi                                                                  */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICJDSUtils::jacobi(
    ValueType solution[],
    const IndexType numRows,
    const IndexType jdsPerm[],
    const IndexType jdsILG[],
    const IndexType /* jdsNumDiagonals */,
    const IndexType jdsDLG[],
    const IndexType jdsJA[],
    const ValueType jdsValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega )
{
    // SCAI_REGION( "MIC.JDS.jacobi" )

    SCAI_LOG_INFO( logger,
                   "jacobi<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution of JDS jacobi iteration for MIC not supported yet" )
    }

    void* solutionPtr = solution;
    const void* oldSolutionPtr = oldSolution;
    const void* rhsPtr = rhs;
    const void* jdsPermPtr = jdsPerm;
    const void* jdsJAPtr = jdsJA;
    const void* jdsILGPtr = jdsILG;
    const void* jdsDLGPtr = jdsDLG;
    const void* jdsValuesPtr = jdsValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( solutionPtr, oldSolutionPtr, rhsPtr, omega, numRows, \
                                       jdsPermPtr, jdsJAPtr, jdsDLGPtr, jdsValuesPtr )
    {
        ValueType* solution = static_cast<ValueType*>( solutionPtr );
        const ValueType* oldSolution = static_cast<const ValueType*>( oldSolutionPtr );
        const ValueType* rhs = static_cast<const ValueType*>( rhsPtr );
        const IndexType* jdsPerm = static_cast<const IndexType*>( jdsPermPtr );
        const IndexType* jdsJA = static_cast<const IndexType*>( jdsJAPtr );
        const IndexType* jdsDLG = static_cast<const IndexType*>( jdsDLGPtr );
        const IndexType* jdsILG = static_cast<const IndexType*>( jdsILGPtr );
        const ValueType* jdsValues = static_cast<const ValueType*>( jdsValuesPtr );

        const ValueType oneMinusOmega = static_cast<ValueType>(1.0) - omega;

        #pragma omp parallel for

        for( IndexType ii = 0; ii < numRows; ii++ )
        {
            const IndexType i = jdsPerm[ii]; // original row index

            ValueType temp = rhs[i];
            IndexType pos = jdsDLG[0] + ii; // index for jdsValues
            ValueType diag = jdsValues[ii]; // diagonal element

            for( IndexType j = 1; j < jdsILG[ii]; j++ )
            {
                temp -= jdsValues[pos] * oldSolution[jdsJA[pos]];
                pos += jdsDLG[j];
            }

            if( omega == static_cast<ValueType>( 1.0 ) )
            {
                solution[i] = temp / diag;
            }
            else if( 0.5 == omega )
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
void MICJDSUtils::jacobiHalo(
    ValueType solution[],
    const IndexType numRows,
    const ValueType localDiagonal[],
    const IndexType numDiagonals,
    const IndexType jdsHaloPerm[],
    const IndexType jdsHaloILG[],
    const IndexType jdsHaloDLG[],
    const IndexType jdsHaloJA[],
    const ValueType jdsHaloValues[],
    const ValueType oldSolution[],
    const ValueType omega )
{
    SCAI_LOG_INFO( logger,
                   "jacobiHalo<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    // SCAI_REGION( "MIC.JDS.jacobiHalo" )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
    }

    if( numRows == 0 )
    {
        return;
    }

    if( numDiagonals == 0 )
    {
        return;
    }

    void* solutionPtr = solution;

    const void* oldSolutionPtr = oldSolution;
    const void* localDiagonalPtr = localDiagonal;
    const void* jdsHaloPermPtr = jdsHaloPerm;
    const void* jdsHaloJAPtr = jdsHaloJA;
    const void* jdsHaloILGPtr = jdsHaloILG;
    const void* jdsHaloDLGPtr = jdsHaloDLG;
    const void* jdsHaloValuesPtr = jdsHaloValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in ( jdsHaloPermPtr, localDiagonalPtr, jdsHaloILGPtr, \
                                               jdsHaloDLGPtr, jdsHaloJAPtr, jdsHaloValuesPtr, \
                                               solutionPtr, oldSolutionPtr, omega )
    {
        ValueType* solution = static_cast<ValueType*>( solutionPtr );

        const ValueType* oldSolution = static_cast<const ValueType*>( oldSolutionPtr );
        const ValueType* localDiagonal = static_cast<const ValueType*>( localDiagonalPtr );
        const IndexType* jdsHaloPerm = static_cast<const IndexType*>( jdsHaloPermPtr );
        const IndexType* jdsHaloJA = static_cast<const IndexType*>( jdsHaloJAPtr );
        const IndexType* jdsHaloDLG = static_cast<const IndexType*>( jdsHaloDLGPtr );
        const IndexType* jdsHaloILG = static_cast<const IndexType*>( jdsHaloILGPtr );
        const ValueType* jdsHaloValues = static_cast<const ValueType*>( jdsHaloValuesPtr );

        // JDS has no row indexes, but number of non-zero rows is known

        const IndexType numNonEmptyRows = jdsHaloDLG[0];

        #pragma omp parallel for

        for( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            ValueType temp = static_cast<ValueType>(0.0);

            const IndexType i = jdsHaloPerm[ii];
            const ValueType diag = localDiagonal[i];

            IndexType pos = ii;

            for( IndexType j = 0; j < jdsHaloILG[ii]; j++ )
            {
                temp += jdsHaloValues[pos] * oldSolution[jdsHaloJA[pos]];
                pos += jdsHaloDLG[j];
            }

            solution[i] -= temp * omega / diag;
        }
    }
}

/* --------------------------------------------------------------------------- */

void MICJDSUtils::Registrator::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using common::context::MIC;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register JDSUtils OpenMP-routines for MIC at kernel registry [" << flag << "]" )

    KernelRegistry::set<JDSKernelTrait::sortRows>( sortRows, MIC, flag );
    KernelRegistry::set<JDSKernelTrait::setInversePerm>( setInversePerm, MIC, flag );
    KernelRegistry::set<JDSKernelTrait::ilg2dlg>( ilg2dlg, MIC, flag );
    KernelRegistry::set<JDSKernelTrait::checkDiagonalProperty>( checkDiagonalProperty, MIC, flag );
}

template<typename ValueType>
void MICJDSUtils::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using common::context::MIC;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register JDSUtils OpenMP-routines for MIC at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<JDSKernelTrait::getValue<ValueType> >( getValue, MIC, flag );
    KernelRegistry::set<JDSKernelTrait::normalGEMV<ValueType> >( normalGEMV, MIC, flag );
    KernelRegistry::set<JDSKernelTrait::jacobi<ValueType> >( jacobi, MIC, flag );
    KernelRegistry::set<JDSKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, MIC, flag );
}

template<typename ValueType, typename OtherValueType>
void MICJDSUtils::RegistratorVO<ValueType, OtherValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using common::context::MIC;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register JDSUtils OpenMP-routines for MIC at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )

    KernelRegistry::set<JDSKernelTrait::scaleValue<ValueType, OtherValueType> >( scaleValue, MIC, flag );
    KernelRegistry::set<JDSKernelTrait::getRow<ValueType, OtherValueType> >( getRow, MIC, flag );
    KernelRegistry::set<JDSKernelTrait::setCSRValues<ValueType, OtherValueType> >( setCSRValues, MIC, flag );
    KernelRegistry::set<JDSKernelTrait::getCSRValues<ValueType, OtherValueType> >( getCSRValues, MIC, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICJDSUtils::RegisterGuard::RegisterGuard()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;

    typedef common::mepr::ContainerV<RegistratorV, ARITHMETIC_MIC> ValueTypes;
    typedef common::mepr::ContainerVO<RegistratorVO, ARITHMETIC_MIC> MoreValueTypes;

    Registrator::initAndReg( flag );
    kregistry::instantiate( flag, ValueTypes() );
    kregistry::instantiate( flag, MoreValueTypes() );
}

MICJDSUtils::RegisterGuard::~RegisterGuard()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    typedef common::mepr::ContainerV<RegistratorV, ARITHMETIC_MIC> ValueTypes;
    typedef common::mepr::ContainerVO<RegistratorVO, ARITHMETIC_MIC> MoreValueTypes;

    Registrator::initAndReg( flag );
    kregistry::instantiate( flag, ValueTypes() );
    kregistry::instantiate( flag, MoreValueTypes() );
}

MICJDSUtils::RegisterGuard MICJDSUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */
