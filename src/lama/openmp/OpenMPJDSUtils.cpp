/**
 * @file OpenMPJDSUtils.cpp
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
 * @brief Implementation of JDS utilities with OpenMP
 * @author Thomas Brandes
 * @date 04.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/openmp/OpenMPJDSUtils.hpp>

// others
#include <lama/openmp/OpenMPUtils.hpp>
#include <lama/openmp/OpenMP.hpp>

#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

// macros
#include <lama/macros/unused.hpp>

// trace
#include <lama/tracing.hpp>

// boost
#include <boost/scoped_array.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( OpenMPJDSUtils::logger, "OpenMP.JDSUtils" )

/* --------------------------------------------------------------------------- */
/*   Implementation of methods                                                 */
/* --------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void OpenMPJDSUtils::getRow(
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
    LAMA_LOG_INFO( logger, "getRow with i = " << i << ", numColumns = " << numColumns << " and numRows = " << numRows )

    //TODO: use OpenMP
    for ( IndexType j = 0; j < numColumns; ++j )
    {
        row[j] = 0.0;
    }

    IndexType ii;

    // check the permutation of row i

    for ( ii = 0; ii < numRows; ii++ )
    {
        if ( perm[ii] == i )
        {
            break;
        }
    }

    IndexType k = 0;

    for ( IndexType jj = 0; jj < ilg[ii]; ++jj )
    {
        row[ja[ii + k]] = static_cast<OtherValueType>( values[ii + k] );
        k += dlg[jj];
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename NoType>
ValueType OpenMPJDSUtils::getValue(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType* dlg,
    const IndexType* ilg,
    const IndexType* perm,
    const IndexType* ja,
    const ValueType* values )
{
    IndexType ii;

    // check the permutation of row i

    for ( ii = 0; ii < numRows; ii++ )
    {
        if ( perm[ii] == i )
        {
            break;
        }
    }

    LAMA_LOG_TRACE( logger, "row " << i << " is now " << ii << ", has " << ilg[ii] << " elements" )
    // search in the found row
    IndexType k = 0;

    for ( IndexType jj = 0; jj < ilg[ii]; jj++ )
    {
        if ( ja[ii + k] == j )
        {
            return values[ii + k];
        }

        k += dlg[jj];
    }

    return 0.0;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType,typename OtherValueType>
void OpenMPJDSUtils::scaleValue(
    const IndexType numRows,
    const IndexType perm[],
    const IndexType ilg[],
    const IndexType dlg[],
    ValueType mValues[],
    const OtherValueType values[] )
{
    LAMA_LOG_INFO( logger, "scaleValue with numRows = " << numRows )

    //TODO: use OpenMP
    for ( IndexType i = 0; i < numRows; i++ )
    {
        IndexType offset = i;
        OtherValueType scalar = values[perm[i]];

        for ( IndexType jj = 0; jj < ilg[i]; jj++ )
        {
            mValues[offset] *= static_cast<ValueType>( scalar );
            offset += dlg[jj];
        }
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

bool OpenMPJDSUtils::checkDiagonalProperty(
    const IndexType numDiagonals,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType perm[],
    const IndexType ja[],
    const IndexType dlg[] )
{
    LAMA_LOG_INFO( logger,
                   "checkDiagonalProperty with numDiagonals = " << numDiagonals << ", numColumns = " << numColumns << " and numRows = " << numRows )
    if ( numRows > 0 )
    {
        bool diagonalProperty = true;

        if ( dlg[0] < std::min( numDiagonals, numColumns ) )
        {
            // not even one entry for each row / column
            diagonalProperty = false;
            return diagonalProperty;
        }

        #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

        for ( IndexType ii = 0; ii < numRows; ++ii )
        {
            if ( !diagonalProperty )
            {
                continue;
            }

            const IndexType i = perm[ii];

            if ( i >= numColumns )
            {
                continue;
            }

            if ( ii >= dlg[0] )
            {
                // ilg[ii] = 0, empty row
                  
                diagonalProperty = false;
            }
            else if ( ja[ii] != i )
            {
                diagonalProperty = false;
            }
        }

        return diagonalProperty;
    }

    return false;
}

/* ------------------------------------------------------------------------------------------------------------------ */

void OpenMPJDSUtils::setInversePerm( IndexType inversePerm[], const IndexType perm[], const IndexType n )
{
    LAMA_LOG_INFO( logger, "compute inverse perm, n = " << n )

    // Parallel execution is safe as perm does not contain a value twice

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
    for ( IndexType ii = 0; ii < n; ii++ )
    {
        IndexType i = perm[ii];
        LAMA_ASSERT_DEBUG( 0 <= i && i < n, "permutation value out of range, perm[" << ii << "] = " << i )
        inversePerm[i] = ii;
    }
}

void OpenMPJDSUtils::sortRows( IndexType ilg[], IndexType perm[], const IndexType n )
{
    // Help array needed, because bucket sort cannot be done in-place

    boost::scoped_array<IndexType> input( new IndexType[n] );

    // Open: can this routine be called where perm is a valid permutation as input

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
    for ( IndexType i = 0; i < n; i++ )
    {
        input[i] = perm[i];
    }

    // The number of buckets is determined by the max value of ilg

    const IndexType maxBucket = OpenMPUtils::maxval( ilg, n );

    LAMA_LOG_INFO( logger, "sort " << n << " values, number of buckets = " << maxBucket )

    // longest row = maxBucket, but rows with length 0 is possible too!

    boost::scoped_array<IndexType> bucket( new IndexType[maxBucket + 1] );

    for ( IndexType i = 0; i <= maxBucket; i++ )
    {
        bucket[i] = 0;
    }

    // counts how many diagonals exist for each possible length
    for ( IndexType i = 0; i < n; i++ )
    {
        bucket[ilg[i]]++;
    }

    for ( IndexType i = 0; i <= maxBucket; i++ )
    {
        LAMA_LOG_DEBUG( logger, "bucket " << i << " has " << bucket[i] << " entries" )
    }

    // use now bucket array for finding right offsets
    // diag length:                 0   1   2   3   4   5
    // number of (now in bucket):   3   4   3   5   1   5
    // becomes (end of first for): 18  14  11   6   5   0
    // later (end of second for):  21  18  14  11   6   5

    IndexType total = 0;
    for ( IndexType i = maxBucket; i >= 0; i-- )
    {
        IndexType cnt = bucket[i];
        bucket[i] = total;
        total += cnt;
        LAMA_LOG_TRACE( logger, "bucket " << i << " offset = " << bucket[i] << ", total = " << total )
    }

    // now we can build the new perm array
    // diagonals with same lengths are moved to position bucket[b] upwards
    for ( IndexType i = 0; i < n; i++ )
    {
        IndexType b = ilg[i];
        LAMA_LOG_TRACE( logger, "perm[" << bucket[b] << "]= " << input[i] )
        perm[bucket[b]++] = input[i];
    }

    // reorganize of ilg has to wait until after filling of perm array is finished
    total = 0;
    for ( IndexType i = maxBucket; i >= 0; i-- )
    {
        LAMA_LOG_DEBUG( logger, "set ilg[" << total << ":" << (bucket[i]-1) << "] = " << i )

        for ( IndexType k = total; k < bucket[i]; k++ )
        {
            ilg[k] = i;
        }
        total = bucket[i];
    }
}

/* --------------------------------------------------------------------------- */

IndexType OpenMPJDSUtils::ilg2dlg(
    IndexType dlg[],
    const IndexType numDiagonals,
    const IndexType ilg[],
    const IndexType numRows )
{
    LAMA_LOG_INFO( logger, "ilg2dlg with numDiagonals = " << numDiagonals << ", numRows = " << numRows )

    if ( numDiagonals == 0 )
    {
        return 0;
    }

    LAMA_ASSERT_EQUAL_DEBUG( numDiagonals, ilg[0] )

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

    #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE ) reduction( +:numTotal )
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

    return numTotal;
}

/* --------------------------------------------------------------------------- */

template<typename JDSValueType,typename CSRValueType>
void OpenMPJDSUtils::getCSRValues(
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
    LAMA_LOG_INFO( logger,
                   "get CSRValues<" << Scalar::getType<JDSValueType>() << ", " << Scalar::getType<CSRValueType>() << ">" << ", #rows = " << numRows << ", #values = " << csrIA[numRows] )

    #pragma omp parallel 
    {
        LAMA_REGION( "OpenMP.JDS->CSR_values" )

        #pragma omp for schedule(LAMA_OMP_SCHEDULE)
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

template<typename JDSValueType,typename CSRValueType>
void OpenMPJDSUtils::setCSRValues(
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
    LAMA_LOG_INFO( logger,
                   "set CSRValues<" << Scalar::getType<JDSValueType>() << ", " << Scalar::getType<CSRValueType>() << ">" << ", #rows = " << numRows << ", #values = " << csrIA[numRows] )

    // parallelization possible as offset array csrIA is available

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.JDS<-CSR_values" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
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
void OpenMPJDSUtils::normalGEMV(
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
    const ValueType jdsValues[],
    class SyncToken* /* syncToken */)
{
    LAMA_LOG_INFO( logger,
                   "normalGEMV<" << Scalar::getType<ValueType>()
                   << ", #threads = " << omp_get_max_threads()
                   << ">, result[" << numRows << "] = " << alpha 
                   << " * A( jds, ndlg = " << ndlg << " ) * x + " << beta << " * y " )

    if ( beta == 0.0 )
    {
        LAMA_LOG_DEBUG( logger, "set result = 0.0" )

        #pragma omp parallel for
        for ( IndexType i = 0; i < numRows; ++i )
        {
            result[i] = 0.0;
        }
    }
    else if ( result == y )
    {
        // result = result * beta

        if ( beta != 1.0 )
        {
            LAMA_LOG_DEBUG( logger, "set result *= beta" )

            #pragma omp parallel for
            for ( IndexType i = 0; i < numRows; ++i )
            {
                result[i] *= beta;
            }
        }
        else
        {
            LAMA_LOG_DEBUG( logger, "result remains unchanged" )
        }
    }
    else
    {
        LAMA_LOG_DEBUG( logger, "set result = beta * y" )

        #pragma omp parallel for
        for ( IndexType i = 0; i < numRows; ++i )
        {
            result[i] = beta * y[i];
        }
    }

    if ( ndlg == 0 )
    {
        return; // definitively empty matrix
    }

    // dlg[0] stands exactly for number of non-empty rows

    IndexType nonEmptyRows = jdsDLG[0];

    LAMA_LOG_DEBUG( logger, "y += alpha * A * x, #non-empty row = " << nonEmptyRows )

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.JDS.normalGEMV" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
        for ( IndexType ii = 0; ii < nonEmptyRows; ii++ )
        {
            ValueType value = 0.0; // sums up final value
            IndexType offset = ii;
            for ( IndexType jj = 0; jj < jdsILG[ii]; jj++ )
            {
                IndexType j = jdsJA[offset];
                LAMA_LOG_TRACE( logger,
                                "compute entry i = " << perm[ii] << ", j = " << j << ", val = " << jdsValues[offset] )
                value += jdsValues[offset] * x[j];
                offset += jdsDLG[jj]; // there is next value for this row
            }

            // scattering needs no synchronization as values of perm are unique

            result[perm[ii]] += alpha * value;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPJDSUtils::normalGEVM(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numColumns,
    const IndexType perm[],
    const IndexType jdsILG[],
    const IndexType ndlg,
    const IndexType jdsDLG[],
    const IndexType jdsJA[],
    const ValueType jdsValues[],
    class SyncToken* UNUSED( syncToken ) )
{
    LAMA_LOG_INFO( logger,
                   "normalGEVM<" << Scalar::getType<ValueType>()
                   << ", #threads = " << omp_get_max_threads()
                   << ">, result[" << numColumns << "] = " << alpha
                   << " * A( jds, ndlg = " << ndlg << " ) * x + " << beta << " * y " )

    if ( beta == 0.0 )
    {
        LAMA_LOG_DEBUG( logger, "set result = 0.0" )

        #pragma omp parallel for
        for ( IndexType i = 0; i < numColumns; ++i )
        {
            result[i] = 0.0;
        }
    }
    else if ( result == y )
    {
        // result = result * beta

        if ( beta != 1.0 )
        {
            LAMA_LOG_DEBUG( logger, "set result *= beta" )

            #pragma omp parallel for
            for ( IndexType i = 0; i < numColumns; ++i )
            {
                result[i] *= beta;
            }
        }
        else
        {
            LAMA_LOG_DEBUG( logger, "result remains unchanged" )
        }
    }
    else
    {
        LAMA_LOG_DEBUG( logger, "set result = beta * y" )

        #pragma omp parallel for
        for ( IndexType i = 0; i < numColumns; ++i )
        {
            result[i] = beta * y[i];
        }
    }

    if ( ndlg == 0 )
    {
        return; // definitively empty matrix
    }

    // dlg[0] stands exactly for number of non-empty rows

    IndexType nonEmptyRows = jdsDLG[0];

    LAMA_LOG_DEBUG( logger, "y += alpha * x * A, #non-empty row = " << nonEmptyRows )

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.JDS.normalGEVM" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
        for( IndexType k = 0; k < numColumns; ++k )
        {
            ValueType value = 0.0; // sums up final value

            for ( IndexType ii = 0; ii < nonEmptyRows; ii++ )
            {
                IndexType offset = ii;
                for ( IndexType jj = 0; jj < jdsILG[ii]; jj++ )
                {
                    IndexType j = jdsJA[offset];
                    if( j == k )
                    {
                        LAMA_LOG_TRACE( logger,
                                        "compute entry i = " << perm[ii] << ", j = " << j << ", matrix val = "
                                        << jdsValues[offset] << ", vector val = " << x[ perm[ii] ] )
                        value += jdsValues[offset] * x[ perm[ii] ];
                    }
                    offset += jdsDLG[jj]; // there is next value for this row
                }
            }
            result[k] += alpha * value;
        }
    }
}

/* --------------------------------------------------------------------------- */
/*     Jacobi                                                                  */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPJDSUtils::jacobi(
    ValueType solution[],
    const IndexType numRows,
    const IndexType jdsPerm[],
    const IndexType jdsILG[],
    const IndexType UNUSED( jdsNumDiagonals ),
    const IndexType jdsDLG[],
    const IndexType jdsJA[],
    const ValueType jdsValues[],
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
        LAMA_REGION( "OpenMP.JDS.jacobi" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )

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

            if ( 1.0 == omega )
            {
                solution[i] = temp / diag;
            }
            else if ( 0.5 == omega )
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
void OpenMPJDSUtils::jacobiHalo(
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
    const ValueType omega,
    class SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "jacobiHalo<" << Scalar::getType<ValueType>() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    if ( syncToken != NULL )
    {
        LAMA_LOG_ERROR( logger, "jacobi called asynchronously, not supported here" )
    }

    if ( numRows == 0 )
    {
        return;
    }

    if ( numDiagonals == 0 )
    {
        return;
    }

    // JDS has no row indexes, but number of non-zero rows is known

    const IndexType numNonEmptyRows = jdsHaloDLG[0];

    LAMA_LOG_DEBUG( logger, "#non empty rows = " << numNonEmptyRows )

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.JDS.jacobiHalo" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
        for ( IndexType ii = 0; ii < numNonEmptyRows; ++ii )
        {
            ValueType temp = 0.0;

            const IndexType i = jdsHaloPerm[ii];
            const ValueType diag = localDiagonal[i];

            IndexType pos = ii;

            for ( IndexType j = 0; j < jdsHaloILG[ii]; j++ )
            {
                temp += jdsHaloValues[pos] * oldSolution[jdsHaloJA[pos]];
                pos += jdsHaloDLG[j];
            }

            LAMA_LOG_TRACE( logger,
                            "jds row " << ii << ", is row " << i << " in halo" 
                            << ", diag = " << diag << ", temp = " << temp )

            solution[i] -= temp * omega / diag;

            LAMA_LOG_TRACE( logger, "solution[" << i << "] = " << solution[i] )
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPJDSUtils::setInterface( JDSUtilsInterface& JDSUtils )
{
    // Register all CUDA routines of this class for the LAMA interface

    LAMA_INTERFACE_REGISTER( JDSUtils, sortRows )
    LAMA_INTERFACE_REGISTER( JDSUtils, setInversePerm )
    LAMA_INTERFACE_REGISTER( JDSUtils, ilg2dlg )
    LAMA_INTERFACE_REGISTER( JDSUtils, checkDiagonalProperty )

    LAMA_INTERFACE_REGISTER_TT( JDSUtils, scaleValue, float, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, scaleValue, float, double )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, scaleValue, double, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, scaleValue, double, double )

    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getRow, float, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getRow, float, double )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getRow, double, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getRow, double, double )

    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getValue, float, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getValue, float, double )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getValue, double, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getValue, double, double )

    LAMA_INTERFACE_REGISTER_TT( JDSUtils, setCSRValues, float, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, setCSRValues, float, double )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, setCSRValues, double, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, setCSRValues, double, double )

    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getCSRValues, float, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getCSRValues, float, double )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getCSRValues, double, float )
    LAMA_INTERFACE_REGISTER_TT( JDSUtils, getCSRValues, double, double )

    LAMA_INTERFACE_REGISTER_T( JDSUtils, normalGEMV, float )
    LAMA_INTERFACE_REGISTER_T( JDSUtils, normalGEMV, double )

    LAMA_INTERFACE_REGISTER_T( JDSUtils, normalGEVM, float )
    LAMA_INTERFACE_REGISTER_T( JDSUtils, normalGEVM, double )

    LAMA_INTERFACE_REGISTER_T( JDSUtils, jacobi, float )
    LAMA_INTERFACE_REGISTER_T( JDSUtils, jacobi, double )

    LAMA_INTERFACE_REGISTER_T( JDSUtils, jacobiHalo, float )
    LAMA_INTERFACE_REGISTER_T( JDSUtils, jacobiHalo, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the JDSUtils routines                             */
/* --------------------------------------------------------------------------- */

bool OpenMPJDSUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::Host );
    setInterface( interface.JDSUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPJDSUtils::initialized = registerInterface();

} // namespace lama
