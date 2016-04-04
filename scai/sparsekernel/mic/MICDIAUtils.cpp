/**
 * @file MICDIAUtils.cpp
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
 * @brief Implementation of DIA utilities with MIC
 * @author Thomas Brandes
 * @date 04.07.2012
 * @since 1.1.0
 */

// hpp
#include <scai/sparsekernel/mic/MICDIAUtils.hpp>

// local library
#include <scai/sparsekernel/DIAKernelTrait.hpp>

// internal scai libraries
#include <scai/utilskernel/mic/MICUtils.hpp>
#include <scai/tasking/mic/MICSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>

// std
#include <cmath>

namespace scai
{

using namespace hmemo;
using tasking::MICSyncToken;
using common::TypeTraits;
using utilskernel::MICUtils;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( MICDIAUtils::logger, "MIC.DIAUtils" )

/* --------------------------------------------------------------------------- */
/*   Implementation of methods                                                 */
/* --------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MICDIAUtils::absMaxVal(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[] )
{
    ValueType maxValue = static_cast<ValueType>(0.0);

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>(0.0);

        #pragma omp for

        for( IndexType i = 0; i < numRows; ++i )
        {
            for( IndexType d = 0; d < numDiagonals; ++d )
            {
                const IndexType j = i + diaOffsets[d];

                if( ( j < 0 ) || ( j >= numColumns ) )
                {
                    continue;
                }

                const ValueType val = std::abs( diaValues[i + d * numRows] );

                if( val > threadVal )
                {
                    threadVal = val;
                }
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

    return maxValue;
}

/* --------------------------------------------------------------------------- */

template<typename DIAValueType,typename CSRValueType>
void MICDIAUtils::getCSRValues(
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
    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << TypeTraits<DIAValueType>::id() << ", " << TypeTraits<CSRValueType>::id() 
                    << ">" << ", #rows = " << numRows << ", #diagonals = " << numDiagonals 
                    << ", #non-zero values = " << csrIA[numRows] << ", diagonalFlag = " << diagonalFlag )

    // we cannot check for correct sizes, but at least for valid pointers

    if( numDiagonals == 0 )
    {
        if( diagonalFlag )
        {
            IndexType n = std::min( numRows, numColumns );

            SCAI_ASSERT_EQUAL_DEBUG( n, csrIA[numRows] )

            #pragma omp parallel for

            for( IndexType i = 0; i < n; i++ )
            {
                csrJA[i] = i;
                csrValues[i] = static_cast<ValueType>(0.0);
            }
        }
        else
        {
            SCAI_ASSERT_EQUAL_DEBUG( 0, csrIA[numRows] )
        }

        return;
    }

    if( numDiagonals > 0 )
    {
        SCAI_ASSERT_DEBUG( diaOffsets != NULL, "offset array of DIA data is NULL" )

        if( numRows > 0 )
        {
            SCAI_ASSERT_DEBUG( diaValues != NULL, "value array of DIA data is NULL" )
        }
    }

    // go through the DIA the same way again and copy the non-zeros

    #pragma omp parallel
    {
        // SCAI_REGION( "MIC.DIA->CSR_values" )

        #pragma omp for

        for( IndexType i = 0; i < numRows; i++ )
        {
            IndexType offset = csrIA[i];

            IndexType ii0 = 0; // first index of diagonal

            if( diagonalFlag && ( i < numColumns ) )
            {
                // store main diagonal at first, must be first diagonal

                SCAI_ASSERT_EQUAL_ERROR( diaOffsets[0], 0 )

                csrJA[offset] = i;
                csrValues[offset] = static_cast<CSRValueType>( diaValues[i] );

                SCAI_LOG_TRACE( logger,
                                "csrJA[" << offset << "] = " << csrJA[offset] << ", csrValues[" << offset << "] = " << csrValues[offset] )

                offset++;
                ii0 = 1;
            }

            for( IndexType ii = ii0; ii < numDiagonals; ii++ )
            {
                IndexType j = i + diaOffsets[ii];

                if( j < 0 )
                {
                    continue;
                }

                if( j >= numColumns )
                {
                    break;
                }

                const DIAValueType value = diaValues[i + ii * numRows];

                bool nonZero = std::abs( value ) > eps;

                if( nonZero )
                {
                    csrJA[offset] = j;
                    csrValues[offset] = static_cast<CSRValueType>( value );
                    SCAI_LOG_TRACE( logger,
                                    "csrJA[" << offset << "] = " << csrJA[offset] << ", csrValues[" << offset << "] = " << csrValues[offset] )

                    offset++;
                }
            }

            SCAI_ASSERT_EQUAL_DEBUG( offset, csrIA[i + 1] )
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename DIAValueType>
void MICDIAUtils::getCSRSizes(
    IndexType csrSizes[],
    bool diagonalFlag,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const DIAValueType diaValues[],
    const DIAValueType eps )
{
    SCAI_LOG_INFO( logger,
                   "get CSRSizes<" << TypeTraits<DIAValueType>::id() << "> for DIA matrix " 
                    << numRows << " x " << numColumns << ", #diagonals = " << numDiagonals 
                    << ", eps = " << eps << ", diagonalFlag = " << diagonalFlag )

    #pragma omp parallel for

    for( IndexType i = 0; i < numRows; i++ )
    {
        IndexType count = 0;

        if( diagonalFlag && ( i < numColumns ) )
        {
            count = 1;
        }

        for( IndexType ii = 0; ii < numDiagonals; ii++ )
        {
            IndexType j = i + diaOffsets[ii]; // column index

            if( j < 0 )
            {
                continue;
            }

            if( j >= numColumns )
            {
                break;
            }

            bool nonZero = std::abs( diaValues[i + ii * numRows] ) > eps;

            if( diagonalFlag && ( i == j ) )
            {
                nonZero = false; // already counted
            }

            if( nonZero )
            {
                count++;
            }
        }

        csrSizes[i] = count;

        SCAI_LOG_TRACE( logger, "#entries in row " << i << ": " << count )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICDIAUtils::normalGEMV(
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
    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << TypeTraits<ValueType>::id() << ">, result[" << numRows << "] = " << alpha << " * A( dia, #diags = " << numDiagonals << " ) * x + " << beta << " * y " )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
    }

    // result := alpha * A * x + beta * y -> result:= beta * y; result += alpha * A

    // SCAI_REGION( "MIC.DIA.normalGEMV" )

    MICUtils::setScale( result, beta, y, numRows );

    void* resultPtr = result;
    const void* xPtr = x;
    const void* diaOffsetsPtr = diaOffsets;
    const void* diaValuesPtr = diaValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( alpha, xPtr, numRows, numColumns, numDiagonals, \
                                       diaOffsetsPtr, diaValuesPtr, resultPtr )
    {
        const IndexType* diaOffsets = static_cast<const IndexType*>( diaOffsetsPtr );
        const ValueType* diaValues = static_cast<const ValueType*>( diaValuesPtr );

        const ValueType* x = static_cast<const ValueType*>( xPtr );
        ValueType* result = static_cast<ValueType*>( resultPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; i++ )
        {
            ValueType accu = static_cast<ValueType>(0.0);

            for( IndexType ii = 0; ii < numDiagonals; ++ii )
            {
                const IndexType j = i + diaOffsets[ii];

                if( j >= numColumns )
                {
                    break;
                }

                if( j >= 0 )
                {
                    accu += diaValues[ii * numRows + i] * x[j];
                }
            }

            result[i] += alpha * accu;
        }
    }
}

/* --------------------------------------------------------------------------- */
/*  Jacobi                                                                     */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICDIAUtils::jacobi(
    ValueType solution[],
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffset[],
    const ValueType diaValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega,
    const IndexType numRows )
{
    // SCAI_REGION( "MIC.DIA.Jacobi" )

    SCAI_LOG_INFO( logger,
                   "jacobi<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", #cols = " << numColumns << ", #diagonals = " << numDiagonals << ", omega = " << omega )

    // main diagonal must be first

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
    }

    void* solutionPtr = solution;
    const void* oldSolutionPtr = oldSolution;
    const void* rhsPtr = rhs;
    const void* diaOffsetPtr = diaOffset;
    const void* diaValuesPtr = diaValues;

    bool error = false; // will be set true if diagonal is not first

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( diaOffsetPtr, diaValuesPtr, numRows, numDiagonals, numColumns, \
                                      solutionPtr, oldSolutionPtr, rhsPtr, omega ), out( error )
    {
        const IndexType* diaOffset = static_cast<const IndexType*>( diaOffsetPtr );

        if( 0 == diaOffset[0] )
        {
            error = false;

            const ValueType* diaValues = static_cast<const ValueType*>( diaValuesPtr );

            const ValueType* oldSolution = static_cast<const ValueType*>( oldSolutionPtr );
            const ValueType* rhs = static_cast<const ValueType*>( rhsPtr );
            ValueType* solution = static_cast<ValueType*>( solutionPtr );

            const ValueType oneMinusOmega = static_cast<ValueType>(1.0) - omega;

            #pragma omp parallel for

            for( IndexType i = 0; i < numRows; i++ )
            {
                ValueType temp = rhs[i];
                ValueType diag = diaValues[i]; // diagonal is first

                for( IndexType ii = 1; ii < numDiagonals; ++ii )
                {
                    const IndexType j = i + diaOffset[ii];

                    if( j >= numColumns )
                    {
                        break;
                    }

                    if( j >= 0 )
                    {
                        temp -= diaValues[ii * numRows + i] * oldSolution[j];
                    }
                }

                solution[i] = omega * ( temp / diag ) + oneMinusOmega * oldSolution[i];
            }
        }
        else
        {
            error = true;
        }
    }

    if( error )
    {
        COMMON_THROWEXCEPTION( "diagonal is not first, diaOffset[0] != 0" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICDIAUtils::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_INFO( logger, "register DIAUtils OpenMP-routines for MIC at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<DIAKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<DIAKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICDIAUtils::MICDIAUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;

    kregistry::mepr::RegistratorV<RegistratorV, ARITHMETIC_MIC_LIST>::call( flag );
}

MICDIAUtils::~MICDIAUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    kregistry::mepr::RegistratorV<RegistratorV, ARITHMETIC_MIC_LIST>::call( flag );
}

MICDIAUtils MICDIAUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */
