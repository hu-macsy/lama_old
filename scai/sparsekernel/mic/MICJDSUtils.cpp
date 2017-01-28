/**
 * @file MICJDSUtils.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Implementation of JDS utilities on MIC platform using offload + OpenMP
 * @author Thomas Brandes
 * @date 04.07.2012
 */

// hpp
#include <scai/sparsekernel/mic/MICJDSUtils.hpp>

// local project
#include <scai/sparsekernel/JDSKernelTrait.hpp>

// other scai libraries
#include <scai/utilskernel/mic/MICUtils.hpp>
#include <scai/utilskernel/BinaryOp.hpp>
#include <scai/hmemo/mic/MICContext.hpp>
#include <scai/tasking/mic/MICSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/OpenMP.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>

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

template<typename ValueType, typename OtherValueType>
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

#pragma offload target( mic : device ) \
    in( rowPtr, permPtr, ilgPtr, dlgPtr, jaPtr, valuesPtr, i, numColumns, numRows )
    {
        OtherValueType* row = reinterpret_cast<OtherValueType*>( rowPtr );
        const IndexType* perm = reinterpret_cast<const IndexType*>( permPtr );
        const IndexType* ilg = reinterpret_cast<const IndexType*>( ilgPtr );
        const IndexType* dlg = reinterpret_cast<const IndexType*>( dlgPtr );
        const IndexType* ja = reinterpret_cast<const IndexType*>( jaPtr );
        const ValueType* values = reinterpret_cast<const ValueType*>( valuesPtr );

        for ( IndexType j = 0; j < numColumns; ++j )
        {
            row[j] = static_cast<OtherValueType>( 0 );
        }

        // one thread will set the row, but every thread searches
        #pragma omp parallel for

        for ( IndexType ii = 0; ii < numRows; ii++ )
        {
            if ( perm[ii] == i )
            {
                IndexType k = 0;

                for ( IndexType jj = 0; jj < ilg[ii]; ++jj )
                {
                    row[ja[ii + k]] = static_cast<OtherValueType>( values[ii + k] );
                    k += dlg[jj];
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

IndexType MICJDSUtils::getValuePos(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType dlg[],
    const IndexType ilg[],
    const IndexType perm[],
    const IndexType ja[] )
{
    int device = MICContext::getCurrentDevice();

    const void* permPtr = perm;
    const void* ilgPtr = ilg;
    const void* dlgPtr = dlg;
    const void* jaPtr = ja;

    IndexType pos = nIndex;

#pragma offload target( mic : device ) in( permPtr, ilgPtr, dlgPtr, jaPtr, i, j, numRows )
    {
        const IndexType* perm = reinterpret_cast<const IndexType*>( permPtr );
        const IndexType* ilg = reinterpret_cast<const IndexType*>( ilgPtr );
        const IndexType* dlg = reinterpret_cast<const IndexType*>( dlgPtr );
        const IndexType* ja = reinterpret_cast<const IndexType*>( jaPtr );

        #pragma omp parallel for

        for ( IndexType ii = 0; ii < numRows; ii++ )
        {
            if ( perm[ii] == i )
            {
                // thread founds the row and now searches the column

                IndexType k = 0;

                for ( IndexType jj = 0; jj < ilg[ii]; jj++ )
                {
                    if ( ja[ii + k] == j )
                    {
                        pos = ii + k;
                    }

                    k += dlg[jj];
                }
            }
        }
    }

    return pos;
}

/* ------------------------------------------------------------------------------------------------------------------ */

IndexType MICJDSUtils::getValuePosCol(
    IndexType row[],
    IndexType pos[],
    const IndexType j,
    const IndexType numRows,
    const IndexType ilg[],
    const IndexType dlg[],
    const IndexType perm[],
    const IndexType ja[] )
{
    int device = MICContext::getCurrentDevice();

    const void* permPtr = perm;
    const void* ilgPtr = ilg;
    const void* dlgPtr = dlg;
    const void* jaPtr = ja;

    void* rowPtr = row;
    void* posPtr = pos;

    IndexType cnt = 0;

#pragma offload target( mic : device ) in( permPtr, ilgPtr, dlgPtr, jaPtr, j, numRows, rowPtr, posPtr )
    {
        const IndexType* perm = reinterpret_cast<const IndexType*>( permPtr );
        const IndexType* ilg = reinterpret_cast<const IndexType*>( ilgPtr );
        const IndexType* dlg = reinterpret_cast<const IndexType*>( dlgPtr );
        const IndexType* ja = reinterpret_cast<const IndexType*>( jaPtr );

        IndexType* row = reinterpret_cast<IndexType*>( rowPtr );
        IndexType* pos = reinterpret_cast<IndexType*>( posPtr );

        #pragma omp parallel for

        for ( IndexType ii = 0; ii < numRows; ii++ )
        {
            IndexType k = 0;

            for ( IndexType jj = 0; jj < ilg[ii]; jj++ )
            {
                IndexType p = ii + k;

                if ( ja[p] == j )
                {
                    IndexType n;

                    #pragma omp critical
                    {
                        n = cnt++;
                    }

                    row[n] = perm[ii];
                    pos[n] = p;

                    break;
                }

                k += dlg[jj];
            }
        }
    }

    return cnt;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType, typename OtherValueType>
void MICJDSUtils::scaleRows(
    ValueType jdsValues[],
    const IndexType numRows,
    const IndexType jdsPerm[],
    const IndexType ilg[],
    const IndexType dlg[],
    const OtherValueType values[] )
{
    SCAI_LOG_INFO( logger, "scaleRows, #rows = " << numRows )

    int device = MICContext::getCurrentDevice();
    const void* jdsPermPtr = jdsPerm;
    const void* ilgPtr = ilg;
    const void* dlgPtr = dlg;
    void* jdsValuesPtr = jdsValues;
    const void* valuesPtr = values;
#pragma offload target( mic : device ) in( jdsPermPtr, ilgPtr, dlgPtr, jdsValuesPtr, \
                                               valuesPtr, numRows ),
    {
        const IndexType* jdsPerm = reinterpret_cast<const IndexType*>( jdsPermPtr );
        const IndexType* ilg = reinterpret_cast<const IndexType*>( ilgPtr );
        const IndexType* dlg = reinterpret_cast<const IndexType*>( dlgPtr );
        ValueType* jdsValues = reinterpret_cast<ValueType*>( jdsValuesPtr );
        const OtherValueType* values = reinterpret_cast<const OtherValueType*>( valuesPtr );
        #pragma omp parallel for

        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType offset = i;
            ValueType scaleRow = static_cast<ValueType>( values[jdsPerm[i]] );

            for ( IndexType jj = 0; jj < ilg[i]; jj++ )
            {
                jdsValues[offset] *= scaleRow;
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

    if ( numRows <= 0 )
    {
        return false;
    }

    // offload dlg[0]

    IndexType dlg0 = MICUtils::getValue( dlg, 0 );

    if ( dlg0 < std::min( numDiagonals, numColumns ) )
    {
        // not even one entry for each row / column
        return false;
    }

    bool diagonalProperty = true;
    void* jdsPermPtr = ( void* ) jdsPerm;
    void* jaPtr = ( void* ) ja;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( jdsPermPtr, jaPtr, numRows, numColumns, dlg0 ), out( diagonalProperty )
    {
        const IndexType* ja = reinterpret_cast<const IndexType*>( jaPtr );
        const IndexType* jdsPerm = reinterpret_cast<const IndexType*>( jdsPermPtr );
        diagonalProperty = true;
        #pragma omp parallel for

        for ( IndexType ii = 0; ii < numRows; ++ii )
        {
            if ( !diagonalProperty )
            {
                continue;
            }

            const IndexType i = jdsPerm[ii];

            if ( i >= numColumns )
            {
                continue;
            }

            if ( ii >= dlg0 )
            {
                // ilg[ii] = 0, empty row
                diagonalProperty = false;
            }
            else if ( ja[ii] != i )
            {
                diagonalProperty = false;
            }
        }
    }

    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

IndexType MICJDSUtils::ilg2dlg(
    IndexType dlg[],
    const IndexType numDiagonals,
    const IndexType ilg[],
    const IndexType numRows )
{
    SCAI_LOG_INFO( logger, "ilg2dlg with numDiagonals = " << numDiagonals << ", numRows = " << numRows )

    if ( numDiagonals == 0 )
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
        IndexType* dlg = reinterpret_cast<IndexType*>( dlgPtr );
        const IndexType* ilg = reinterpret_cast<const IndexType*>( ilgPtr );
        numTotal = 0;
        #pragma omp parallel for reduction( +:numTotal )

        for ( IndexType i = 0; i < numRows; ++i )
        {
            IndexType nd1 = ilg[i];
            IndexType nd2 = 0;

            if ( i + 1 < numRows )
            {
                nd2 = ilg[i + 1];
            }

            // fill in dlg only if nd2 < nd1

            for ( IndexType j = nd2; j < nd1; j++ )
            {
                dlg[j] = i + 1;
            }

            numTotal += nd1;
        }
    }
    return numTotal;
}

/* --------------------------------------------------------------------------- */

template<typename JDSValueType, typename CSRValueType>
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
        const JDSValueType* jdsValues = reinterpret_cast<const JDSValueType*>( jdsValuesPtr );
        const IndexType* jdsJA = reinterpret_cast<const IndexType*>( jdsJAPtr );
        const IndexType* jdsInversePerm = reinterpret_cast<const IndexType*>( jdsInversePermPtr );
        const IndexType* jdsILG = reinterpret_cast<const IndexType*>( jdsILGPtr );
        const IndexType* jdsDLG = reinterpret_cast<const IndexType*>( jdsDLGPtr );
        const IndexType* csrIA = reinterpret_cast<const IndexType*>( csrIAPtr );
        IndexType* csrJA = reinterpret_cast<IndexType*>( csrJAPtr );
        CSRValueType* csrValues = reinterpret_cast<CSRValueType*>( csrValuesPtr );
        #pragma omp parallel for

        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType ii = jdsInversePerm[i]; // where to find row i in JDS storage
            const IndexType numValuesInRow = jdsILG[ii];
            IndexType jdsOffset = ii; // run through input JDS data
            IndexType offset = csrIA[i]; // run through output data

            for ( IndexType jj = 0; jj < numValuesInRow; jj++ )
            {
                csrJA[offset + jj] = jdsJA[jdsOffset];
                csrValues[offset + jj] = static_cast<CSRValueType>( jdsValues[jdsOffset] );
                jdsOffset += jdsDLG[jj];
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename JDSValueType, typename CSRValueType>
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
        JDSValueType* jdsValues = reinterpret_cast<JDSValueType*>( jdsValuesPtr );
        IndexType* jdsJA = reinterpret_cast<IndexType*>( jdsJAPtr );
        const IndexType* jdsPerm = reinterpret_cast<const IndexType*>( jdsPermPtr );
        const IndexType* jdsILG = reinterpret_cast<const IndexType*>( jdsILGPtr );
        const IndexType* jdsDLG = reinterpret_cast<const IndexType*>( jdsDLGPtr );
        const IndexType* csrIA = reinterpret_cast<const IndexType*>( csrIAPtr );
        const IndexType* csrJA = reinterpret_cast<const IndexType*>( csrJAPtr );
        const CSRValueType* csrValues = reinterpret_cast<const CSRValueType*>( csrValuesPtr );
        #pragma omp parallel for

        for ( IndexType ii = 0; ii < numRows; ii++ )
        {
            IndexType i = jdsPerm[ii];
            IndexType offset = ii;

            for ( IndexType jdsJJ = 0, csrJJ = csrIA[i]; jdsJJ < jdsILG[ii]; jdsJJ++, csrJJ++ )
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

    if ( beta == common::constants::ZERO )
    {
        MICUtils::setVal( result, numRows, ValueType( 0 ), utilskernel::binary::COPY );
    }
    else if ( result == y )
    {
        // result = result * beta
        if ( beta == common::constants::ONE )
        {
            SCAI_LOG_DEBUG( logger, "result remains unchanged" )
        }
        else
        {
            MICUtils::setVal( result, numRows, beta, utilskernel::binary::MULT );
        }
    }
    else
    {
        // result = beta * y
        MICUtils::binaryOpScalar1( result, beta, y, numRows, utilskernel::binary::MULT );
    }

    if ( ndlg == 0 )
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
    const ValueType* alphaPtr = &alpha;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ), in( resultPtr, xPtr, jdsDLGPtr, jdsILGPtr, jdsJAPtr, jdsValuesPtr, alphaPtr[0:1] )
    {
        ValueType* result = reinterpret_cast<ValueType*>( resultPtr );
        const ValueType* x = reinterpret_cast<const ValueType*>( xPtr );
        const IndexType* perm = reinterpret_cast<const IndexType*>( permPtr );
        const IndexType* jdsILG = reinterpret_cast<const IndexType*>( jdsILGPtr );
        const IndexType* jdsDLG = reinterpret_cast<const IndexType*>( jdsDLGPtr );
        const IndexType* jdsJA = reinterpret_cast<const IndexType*>( jdsJAPtr );
        const ValueType* jdsValues = reinterpret_cast<const ValueType*>( jdsValuesPtr );
        const ValueType& alphaRef = *alphaPtr;
        // dlg[0] stands exactly for number of non-empty rows
        IndexType nonEmptyRows = jdsDLG[0];
        #pragma omp parallel for

        for ( IndexType ii = 0; ii < nonEmptyRows; ii++ )
        {
            ValueType value = 0; // sums up final value

            IndexType offset = ii;

            for ( IndexType jj = 0; jj < jdsILG[ii]; jj++ )
            {
                IndexType j = jdsJA[offset];
                value += jdsValues[offset] * x[j];
                offset += jdsDLG[jj]; // there is next value for this row
            }

            // scattering needs no synchronization as values of perm are unique
            result[perm[ii]] += alphaRef * value;
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
    const ValueType* omegaPtr = &omega;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ), in( solutionPtr, oldSolutionPtr, rhsPtr, omegaPtr[0:1], numRows, \
                                       jdsPermPtr, jdsJAPtr, jdsDLGPtr, jdsValuesPtr )
    {
        ValueType* solution = reinterpret_cast<ValueType*>( solutionPtr );
        const ValueType* oldSolution = reinterpret_cast<const ValueType*>( oldSolutionPtr );
        const ValueType* rhs = reinterpret_cast<const ValueType*>( rhsPtr );
        const IndexType* jdsPerm = reinterpret_cast<const IndexType*>( jdsPermPtr );
        const IndexType* jdsJA = reinterpret_cast<const IndexType*>( jdsJAPtr );
        const IndexType* jdsDLG = reinterpret_cast<const IndexType*>( jdsDLGPtr );
        const IndexType* jdsILG = reinterpret_cast<const IndexType*>( jdsILGPtr );
        const ValueType* jdsValues = reinterpret_cast<const ValueType*>( jdsValuesPtr );
        const ValueType& omegaRef = *omegaPtr;
        const ValueType oneMinusOmega = static_cast<ValueType>( 1 ) - omegaRef;
        #pragma omp parallel for

        for ( IndexType ii = 0; ii < numRows; ii++ )
        {
            const IndexType i = jdsPerm[ii]; // original row index
            ValueType temp = rhs[i];
            IndexType pos = jdsDLG[0] + ii; // index for jdsValues
            ValueType diag = jdsValues[ii]; // diagonal element

            for ( IndexType j = 1; j < jdsILG[ii]; j++ )
            {
                temp -= jdsValues[pos] * oldSolution[jdsJA[pos]];
                pos += jdsDLG[j];
            }

            if ( omegaRef == static_cast<ValueType>( 1 ) )
            {
                solution[i] = temp / diag;
            }
            else if ( static_cast<ValueType>( 0.5 ) == omegaRef )
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

    if ( numRows == 0 )
    {
        return;
    }

    if ( numDiagonals == 0 )
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
    const ValueType* omegaPtr = &omega;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ) in ( jdsHaloPermPtr, localDiagonalPtr, jdsHaloILGPtr, \
                                               jdsHaloDLGPtr, jdsHaloJAPtr, jdsHaloValuesPtr, \
                                               solutionPtr, oldSolutionPtr, omegaPtr[0:1] )
    {
        ValueType* solution = reinterpret_cast<ValueType*>( solutionPtr );
        const ValueType* oldSolution = reinterpret_cast<const ValueType*>( oldSolutionPtr );
        const ValueType* localDiagonal = reinterpret_cast<const ValueType*>( localDiagonalPtr );
        const IndexType* jdsHaloPerm = reinterpret_cast<const IndexType*>( jdsHaloPermPtr );
        const IndexType* jdsHaloJA = reinterpret_cast<const IndexType*>( jdsHaloJAPtr );
        const IndexType* jdsHaloDLG = reinterpret_cast<const IndexType*>( jdsHaloDLGPtr );
        const IndexType* jdsHaloILG = reinterpret_cast<const IndexType*>( jdsHaloILGPtr );
        const ValueType* jdsHaloValues = reinterpret_cast<const ValueType*>( jdsHaloValuesPtr );
        const ValueType& omegaRef = *omegaPtr;
        // JDS has no row indexes, but number of non-zero rows is known
        const IndexType numNonEmptyRows = jdsHaloDLG[0];
        #pragma omp parallel for

        for ( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            ValueType temp = 0;
            const IndexType i = jdsHaloPerm[ii];
            const ValueType diag = localDiagonal[i];
            IndexType pos = ii;

            for ( IndexType j = 0; j < jdsHaloILG[ii]; j++ )
            {
                temp += jdsHaloValues[pos] * oldSolution[jdsHaloJA[pos]];
                pos += jdsHaloDLG[j];
            }

            solution[i] -= temp * omegaRef / diag;
        }
    }
}

/* --------------------------------------------------------------------------- */

void MICJDSUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    // context for which kernels will be added

    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_DEBUG( logger, "register[flag=" << flag << "]: untyped routines" )

    KernelRegistry::set<JDSKernelTrait::ilg2dlg>( ilg2dlg, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::checkDiagonalProperty>( checkDiagonalProperty, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::getValuePos>( getValuePos, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::getValuePosCol>( getValuePosCol, ctx, flag );
}

template<typename ValueType>
void MICJDSUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_DEBUG( logger, "register[flag=" << flag << "]: T = " << common::TypeTraits<ValueType>::id() )

    KernelRegistry::set<JDSKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void MICJDSUtils::RegistratorVO<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_DEBUG( logger, "register[flag=" << flag << "]: TT " <<
                    common::TypeTraits<ValueType>::id() << ", " << common::TypeTraits<OtherValueType>::id() )

    KernelRegistry::set<JDSKernelTrait::scaleRows<ValueType, OtherValueType> >( scaleRows, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::getRow<ValueType, OtherValueType> >( getRow, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::setCSRValues<ValueType, OtherValueType> >( setCSRValues, ctx, flag );
    KernelRegistry::set<JDSKernelTrait::getCSRValues<ValueType, OtherValueType> >( getCSRValues, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICJDSUtils::RegisterGuard::RegisterGuard()
{
    SCAI_LOG_INFO( logger, "register JDSUtils routines for MIC(OpenMP,offload) at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_MIC_LIST, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels( flag );
}

MICJDSUtils::RegisterGuard::~RegisterGuard()
{
    SCAI_LOG_INFO( logger, "unregister JDSUtils routines for MIC(OpenMP,offload) at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_MIC_LIST, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels( flag );
}

MICJDSUtils::RegisterGuard MICJDSUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */
