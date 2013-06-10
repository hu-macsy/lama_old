/**
 * @file OpenMPDIAUtils.cpp
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
 * @brief Implementation of DIA utilities with OpenMP
 * @author Thomas Brandes
 * @date 04.07.2012
 * @since 1.0.0
 */

// for dll_import
#include <lama/openmp/OpenMPUtils.hpp>

// others
#include <lama/openmp/OpenMPDIAUtils.hpp>

#include <lama/LAMAInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/tracing.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

// boost
#include <boost/scoped_array.hpp>

#include <cmath>
#include <typeinfo>

namespace lama
{

LAMA_LOG_DEF_LOGGER( OpenMPDIAUtils::logger, "OpenMP.DIAUtils" )

/* --------------------------------------------------------------------------- */
/*   Implementation of methods                                                 */
/* --------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPDIAUtils::absMaxVal(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[] )
{
    ValueType maxValue = static_cast<ValueType>( 0.0 );

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>( 0.0 );

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType d = 0; d < numDiagonals; ++d )
            {
                const IndexType j = i + diaOffsets[d];

                if ( ( j < 0 ) || ( j >= numColumns ) )
                {
                    continue;
                }

                const ValueType val = std::abs( diaValues[i + d * numRows] );

                if ( val > threadVal )
                {
                    threadVal = val;
                }
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

    return maxValue;
}

/* --------------------------------------------------------------------------- */

template<typename DIAValueType,typename CSRValueType>
void OpenMPDIAUtils::getCSRValues(
    IndexType csrJA[],
    CSRValueType csrValues[],
    const IndexType csrIA[],
    const bool diagonalFlag,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const DIAValueType diaValues[],
    const DIAValueType eps )
{
    LAMA_LOG_INFO( logger,
                   "get CSRValues<" << typeid( DIAValueType ).name() << ", " << typeid( CSRValueType ).name() << ">"
                    << ", #rows = " << numRows << ", #diagonals = " << numDiagonals 
                    << ", #non-zero values = " << csrIA[numRows] << ", diagonalFlag = " << diagonalFlag )

    // we cannot check for correct sizes, but at least for valid pointers

    if ( numDiagonals == 0 )
    {
        if ( diagonalFlag )
        {
            IndexType n = std::min( numRows, numColumns );

            LAMA_ASSERT_EQUAL_DEBUG( n, csrIA[ numRows ] )

            #pragma omp parallel for schedule( LAMA_OMP_SCHEDULE )
            for ( IndexType i = 0; i < n; i++ )
            {
                csrJA[i] = i;
                csrValues[i] = 0.0;
            }
        }
        else
        {
            LAMA_ASSERT_EQUAL_DEBUG( 0, csrIA[ numRows ] )
        }
        return;
    }

    if ( numDiagonals > 0 )
    {
        LAMA_ASSERT_DEBUG( diaOffsets != NULL, "offset array of DIA data is NULL" )

        if ( numRows > 0 )
        {
            LAMA_ASSERT_DEBUG( diaValues != NULL, "value array of DIA data is NULL" )
        }
    }

    // go through the DIA the same way again and copy the non-zeros

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.DIA->CSR_values" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )

        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType offset = csrIA[i];

            IndexType ii0 = 0; // first index of diagonal

            if ( diagonalFlag && ( i < numColumns ) )
            {
                // store main diagonal at first, must be first diagonal

                LAMA_ASSERT_EQUAL_ERROR( diaOffsets[0], 0 )

                csrJA[offset] = i;
                csrValues[offset] = static_cast<CSRValueType>( diaValues[i] );

                LAMA_LOG_TRACE( logger,
                            "csrJA[" << offset << "] = " << csrJA[offset] << ", csrValues[" << offset << "] = " << csrValues[offset] )

                offset++;
                ii0 = 1;
            }

            for ( IndexType ii = ii0; ii < numDiagonals; ii++ )
            {
                IndexType j = i + diaOffsets[ii];

                if ( j < 0 )
                {
                    continue;
                }

                if ( j >= numColumns )
                {
                    break;
                }

                const DIAValueType value = diaValues[i + ii * numRows];

                bool nonZero = std::abs( value ) > eps;

                if ( nonZero )
                {
                    csrJA[offset] = j;
                    csrValues[offset] = static_cast<CSRValueType>( value );
                    LAMA_LOG_TRACE( logger,
                                    "csrJA[" << offset << "] = " << csrJA[offset] << ", csrValues[" << offset << "] = " << csrValues[offset] )

                    offset++;
                }
            }

            LAMA_ASSERT_EQUAL_DEBUG( offset, csrIA[i+1] )
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename DIAValueType>
void OpenMPDIAUtils::getCSRSizes(
    IndexType csrSizes[],
    bool diagonalFlag,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const DIAValueType diaValues[],
    const DIAValueType eps )
{
    LAMA_LOG_INFO( logger,
                   "get CSRSizes<" << typeid( DIAValueType ).name() << "> for DIA matrix " << numRows << " x " << numColumns << ", #diagonals = " << numDiagonals << ", eps = " << eps << ", diagonalFlag = " << diagonalFlag )

    #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
    for ( IndexType i = 0; i < numRows; i++ )
    {
        IndexType count = 0;

        if ( diagonalFlag && ( i < numColumns ) )
        {
            count = 1;
        }

        for ( IndexType ii = 0; ii < numDiagonals; ii++ )
        {
            IndexType j = i + diaOffsets[ii]; // column index

            if ( j < 0 )
            {
                continue;
            }

            if ( j >= numColumns )
            {
                break;
            }

            bool nonZero = std::abs( diaValues[i + ii * numRows] ) > eps;

            if ( diagonalFlag && ( i == j ) )
            {
                nonZero = false; // already counted
            }

            if ( nonZero )
            {
                count++;
            }
        }

        csrSizes[i] = count;

        LAMA_LOG_TRACE( logger, "#entries in row " << i << ": " << count )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDIAUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[] )
{
    LAMA_LOG_INFO( logger, "normalGEMV<" << typeid(ValueType).name() << ">, n = " 
                           << numRows << ", d = " << numDiagonals  )

    if ( beta == 0.0 )
    {
        #pragma omp parallel for
        for ( IndexType i = 0; i < numRows; ++i )
        {
            result[i] = 0.0;
        }
    }
    else if ( beta == 1.0 )
    {
        if ( result != y )
        {
            for ( IndexType i = 0; i < numRows; ++i )
            {
                result[i] = y[i];
            }
        }
    }
    else
    {
        #pragma omp parallel for
        for ( IndexType i = 0; i < numRows; ++i )
        {
            result[i] = beta * y[i];
        }
    }

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.DIA.normalGEMV" )

        #pragma omp for schedule ( LAMA_OMP_SCHEDULE )
        for ( IndexType i = 0; i < numRows; i++ )
        {
            ValueType accu = 0.0;

            for ( IndexType ii = 0; ii < numDiagonals; ++ii )
            {
                const IndexType j = i + diaOffsets[ii];

                if ( j >= numColumns )
                {
                    break;
                }
                if ( j >= 0 )
                {
                    accu += diaValues[ii * numRows + i] * x[j];
                }
            }
            result[i] += alpha * accu;
        }

        if ( LAMA_LOG_TRACE_ON( logger ) )
        {
            std::cout << "NormalGEMV: result = ";
            for ( IndexType i = 0; i < numRows; ++i )
            {
                std::cout << " " << result[i];
            }
            std::cout << std::endl;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDIAUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[],
    SyncToken* syncToken )
{
    if ( !syncToken )
    {
        normalGEMV( result, alpha, x, beta, y, numRows, numColumns, numDiagonals, diaOffsets, diaValues );
    }
    else
    {
        LAMA_THROWEXCEPTION( "no asynchronous support due to boost problem" )

        // asynchronous execution is done by calling an own thread at higher level
    }
}

/* --------------------------------------------------------------------------- */
/*  Jacobi                                                                     */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDIAUtils::jacobi(
    ValueType solution[],
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffset[],
    const ValueType diaValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega,
    const IndexType numRows,
    class SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "jacobi<" << typeid(ValueType).name() << ">" 
                    << ", #rows = " << numRows << ", #cols = " << numColumns 
                    << ", #diagonals = " << numDiagonals << ", omega = " << omega )

    LAMA_ASSERT_EQUAL_DEBUG( 0, diaOffset[0] )
    // main diagonal must be first

    if ( syncToken != NULL )
    {
        LAMA_LOG_ERROR( logger, "jacobi called asynchronously, not supported here" )
    }

    const ValueType oneMinusOmega = static_cast<ValueType>( 1.0 - omega );

    #pragma omp parallel
    {
        LAMA_REGION( "OpenMP.DIA.Jacobi" )

        #pragma omp for schedule( LAMA_OMP_SCHEDULE )
        for ( IndexType i = 0; i < numRows; i++ )
        {
            ValueType temp = rhs[i];
            ValueType diag = diaValues[i]; // diagonal is first
    
            for ( IndexType ii = 1; ii < numDiagonals; ++ii )
            {
                const IndexType j = i + diaOffset[ii];
    
                if ( j >= numColumns )
                {
                    break;
                }
                if ( j >= 0 )
                {
                    temp -= diaValues[ii * numRows + i] * oldSolution[j];
                }
            }

            solution[i] = omega * ( temp / diag ) + oneMinusOmega * oldSolution[i];
        }
    }
}

/* --------------------------------------------------------------------------- */

void OpenMPDIAUtils::setInterface( DIAUtilsInterface& DIAUtils )
{
    // Register all CUDA routines of this class for the LAMA interface

    LAMA_INTERFACE_REGISTER_T( DIAUtils, getCSRSizes, float )
    LAMA_INTERFACE_REGISTER_T( DIAUtils, getCSRSizes, double )

    LAMA_INTERFACE_REGISTER_TT( DIAUtils, getCSRValues, float, float )
    LAMA_INTERFACE_REGISTER_TT( DIAUtils, getCSRValues, float, double )
    LAMA_INTERFACE_REGISTER_TT( DIAUtils, getCSRValues, double, float )
    LAMA_INTERFACE_REGISTER_TT( DIAUtils, getCSRValues, double, double )

    LAMA_INTERFACE_REGISTER_T( DIAUtils, absMaxVal, float )
    LAMA_INTERFACE_REGISTER_T( DIAUtils, absMaxVal, double )

    LAMA_INTERFACE_REGISTER_T( DIAUtils, normalGEMV, float )
    LAMA_INTERFACE_REGISTER_T( DIAUtils, normalGEMV, double )

    LAMA_INTERFACE_REGISTER_T( DIAUtils, jacobi, float )
    LAMA_INTERFACE_REGISTER_T( DIAUtils, jacobi, double )
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the DIAUtils routines                             */
/* --------------------------------------------------------------------------- */

bool OpenMPDIAUtils::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::Host );
    setInterface( interface.DIAUtils );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPDIAUtils::initialized = registerInterface();

} // namespace lama
