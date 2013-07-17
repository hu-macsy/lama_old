/**
 * @file OpenMPCOOUtils.cpp
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
 * @brief Implementation of COO storage utilities with OpenMP
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/openmp/OpenMPCOOUtils.hpp>
#include <lama/openmp/OpenMPUtils.hpp>

// others
#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/tracing.hpp>

#include <lama/task/TaskSyncToken.hpp>
#include <lama/openmp/OpenMP.hpp>

namespace lama
{

/* --------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( OpenMPCOOUtils::logger, "OpenMP.COOUtils" )

/* --------------------------------------------------------------------------- */
/*     Template implementations                                                */
/* --------------------------------------------------------------------------- */

void OpenMPCOOUtils::getCSRSizes(
    IndexType csrSizes[],
    const IndexType numRows,
    const IndexType numValues,
    const IndexType cooIA[] )
{
    LAMA_LOG_INFO( logger, "get CSR sizes, #rows = " << numRows << ", #values = " << numValues )

    // initialize size array for each row

    for ( IndexType i = 0; i < numRows; i++ )
    {
        csrSizes[i] = 0;
    }

    // increment size of a row for each used row value

    for ( IndexType k = 0; k < numValues; k++ )
    {
        IndexType i = cooIA[k];
        LAMA_ASSERT_DEBUG( i < numRows, "cooIA[" << k << "] = " << i << " out of range, #rows = " << numRows )
        csrSizes[i]++;
    }
}

/* --------------------------------------------------------------------------- */

template<typename COOValueType,typename CSRValueType>
void OpenMPCOOUtils::getCSRValues( IndexType csrJA[], CSRValueType csrValues[], IndexType csrIA[], // used as tmp, remains unchanged
                                   const IndexType numRows,
                                   const IndexType numValues,
                                   const IndexType cooIA[],
                                   const IndexType cooJA[],
                                   const COOValueType cooValues[] )
{
    LAMA_LOG_INFO( logger,
                   "get CSRValues<" << Scalar::getType<COOValueType>() << ", " << Scalar::getType<CSRValueType>() << ">" << ", #rows = " << numRows << ", #values = " << numValues )

    // traverse the non-zero values and put data at the right places

    for ( IndexType k = 0; k < numValues; k++ )
    {
        IndexType i = cooIA[k];

        IndexType& offset = csrIA[i];

        csrJA[offset] = cooJA[k];
        csrValues[offset] = static_cast<CSRValueType>( cooValues[k] );

        LAMA_LOG_DEBUG( logger, "row " << i << ": new offset = " << offset )

        offset++;
    }

    // set back the old offsets in csrIA

    for ( IndexType i = numRows; i > 0; --i )
    {
        csrIA[i] = csrIA[i - 1];
    }

    csrIA[0] = 0;

    LAMA_ASSERT_EQUAL_DEBUG( csrIA[numRows], numValues )
}

/* --------------------------------------------------------------------------- */

template<typename COOValueType,typename CSRValueType>
void OpenMPCOOUtils::setCSRValues(
    IndexType cooIA[],
    IndexType cooJA[],
    COOValueType cooValues[],
    const IndexType numRows,
    const IndexType numDiagonals,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const CSRValueType csrValues[],
    const bool csrDiagonalProperty )
{
    LAMA_LOG_INFO( logger,
                   "set CSRValues<" << Scalar::getType<COOValueType>() << ", " << Scalar::getType<CSRValueType>() << ">" << ", #rows = " << numRows << ", #values = " << csrIA[numRows] )

    if ( numDiagonals == 0 || csrDiagonalProperty )
    {
        LAMA_LOG_INFO( logger, "parallel fill in possible, #diagonal elements = " << numDiagonals )

        // parallel execution only possible if we have no separate diagonal elements
        // or if CSR data has diagonal property

        #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
        for ( IndexType i = 0; i < numRows; ++i )
        {
            IndexType csrOffset = csrIA[i];
            IndexType cooOffset = csrOffset;

            if ( i < numDiagonals )
            {
                // csr data must have the diagonal property, should have been checked before

                LAMA_ASSERT_DEBUG( csrOffset < csrIA[i+1],
                                   "diagonal property requires at least one entry in row " << i )
                LAMA_ASSERT_EQUAL_DEBUG( csrJA[csrOffset], i )

                // diagonal elements will be the first nrows entries

                cooIA[i] = i;
                cooJA[i] = i;
                cooValues[i] = static_cast<COOValueType>( csrValues[csrOffset] );

                csrOffset += 1; // do not fill diagonal element again
                cooOffset += numDiagonals - i; // offset in coo moves
            }

            // now fill remaining part of row i

            for ( IndexType jj = csrOffset; jj < csrIA[i + 1]; ++jj )
            {
                cooIA[cooOffset] = i;
                cooJA[cooOffset] = csrJA[jj];
                cooValues[cooOffset] = static_cast<COOValueType>( csrValues[jj] );
                cooOffset++;
            }
        }
    }
    else
    {
        // initialize diagonal elements in case of non-availablity

        #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
        for ( IndexType i = 0; i < numDiagonals; ++i )
        {
            cooIA[i] = i;
            cooJA[i] = i;
            cooValues[i] = 0.0;
        }

        LAMA_LOG_INFO( logger, "serial fill in, #diagonal elements = " << numDiagonals )

        // only serial fill-in possible as we do not now how many diagonal elements are available

        int cooOffset = numDiagonals;

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
            {
                int j = csrJA[jj];

                CSRValueType val = csrValues[jj];

                if ( i == j && i < numDiagonals )
                {
                    cooValues[i] = static_cast<COOValueType>( val );
                }
                else
                {
                    cooIA[cooOffset] = i;
                    cooJA[cooOffset] = j;
                    cooValues[cooOffset] = static_cast<COOValueType>( val );
                    cooOffset++;
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCOOUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[],
    const IndexType numValues,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "normalGEMV<" << Scalar::getType<ValueType>()
                   << ", #threads = " << omp_get_max_threads()
                   << ">, result[" << numRows << "] = " << alpha
                   << " * A( coo, #vals = " << numValues << " ) * x + " << beta << " * y " )

    if ( syncToken )
    {
        LAMA_THROWEXCEPTION( "asynchronous execution not supported here, do it by a task" )
    }

    // result := alpha * A * x + beta * y -> result:= beta * y; result += alpha * A

    OpenMPUtils::setScale( result, beta, y, numRows );

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.COO.normalGEMV" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
        for ( IndexType k = 0; k < numValues; ++k )
        {
            IndexType i = cooIA[k];
            IndexType j = cooJA[k];

            // we must use atomic updates as different threads might update same row i

            const ValueType resultUpdate = alpha * cooValues[k] * x[j];

            #pragma omp atomic
            result[i] += resultUpdate;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCOOUtils::normalGEVM(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[],
    const IndexType numValues,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "normalGEMV<" << Scalar::getType<ValueType>()
                   << ", #threads = " << omp_get_max_threads()
                   << ">, result[" << numRows << "] = " << alpha
                   << " * A( coo, #vals = " << numValues << " ) * x + " << beta << " * y " )

    if ( syncToken )
    {
        LAMA_THROWEXCEPTION( "asynchronous execution not supported here, do it by a task" )
    }

    // result := alpha * x * A + beta * y -> result:= beta * y; result += alpha * x * A

    OpenMPUtils::setScale( result, beta, y, numRows );

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.COO.normalGEMV" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
        for ( IndexType k = 0; k < numValues; ++k )
        {
            IndexType i = cooIA[k];
            IndexType j = cooJA[k];

            // we must use atomic updates as different threads might update same row i

            const ValueType resultUpdate = alpha * cooValues[k] * x[i];

            #pragma omp atomic
            result[j] += resultUpdate;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCOOUtils::jacobi(
    ValueType* solution,
    const IndexType cooNumValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega,
    const IndexType numRows,
    class SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "jacobi<" << Scalar::getType<ValueType>() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    if ( syncToken )
    {
        LAMA_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    // solution = omega * ( rhs - B * oldSolution ) * dinv + ( 1 - omega * oldSolution
    // done in two steps
    // solution = omega * rhs * dinv + ( 1 - omega * oldSolution
    // solution -= omega * B * oldSolution * dinv

    const ValueType oneMinusOmega = static_cast<ValueType>( 1.0 ) - omega;

    #pragma omp parallel for
    for ( IndexType i = 0; i < numRows; ++i )
    {
        solution[i] = omega * rhs[i] / cooValues[i] + oneMinusOmega * oldSolution[i];
    }

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.COO.jacobi" )

        #pragma omp for
        for ( IndexType k = numRows; k < cooNumValues; ++k )
        {
            IndexType i = cooIA[k];
            IndexType j = cooJA[k];

            // we must use atomic updates as different threads might update same row i

            const ValueType update = omega * cooValues[k] * oldSolution[j] / cooValues[i];

            #pragma omp atomic
            solution[i] -= update;
        }
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPCOOUtils::setInterface( COOUtilsInterface& COOUtils )
{
    LAMA_INTERFACE_REGISTER( COOUtils, getCSRSizes )

    LAMA_INTERFACE_REGISTER_TT( COOUtils, setCSRValues, float, float )
    LAMA_INTERFACE_REGISTER_TT( COOUtils, setCSRValues, float, double )
    LAMA_INTERFACE_REGISTER_TT( COOUtils, setCSRValues, double, float )
    LAMA_INTERFACE_REGISTER_TT( COOUtils, setCSRValues, double, double )

    LAMA_INTERFACE_REGISTER_TT( COOUtils, getCSRValues, float, float )
    LAMA_INTERFACE_REGISTER_TT( COOUtils, getCSRValues, float, double )
    LAMA_INTERFACE_REGISTER_TT( COOUtils, getCSRValues, double, float )
    LAMA_INTERFACE_REGISTER_TT( COOUtils, getCSRValues, double, double )

    LAMA_INTERFACE_REGISTER_T( COOUtils, normalGEMV, float )
    LAMA_INTERFACE_REGISTER_T( COOUtils, normalGEMV, double )

    LAMA_INTERFACE_REGISTER_T( COOUtils, normalGEVM, float )
    LAMA_INTERFACE_REGISTER_T( COOUtils, normalGEVM, double )

    LAMA_INTERFACE_REGISTER_T( COOUtils, jacobi, float )
    LAMA_INTERFACE_REGISTER_T( COOUtils, jacobi, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool OpenMPCOOUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::Host );
    setInterface( interface.COOUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPCOOUtils::initialized = registerInterface();

} // namespace lama
