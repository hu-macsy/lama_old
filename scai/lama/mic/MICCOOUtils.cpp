/**
 * @file MICCOOUtils.cpp
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
 * @brief Implementation of COO storage utilities with MIC
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.1.0
 */

// hpp
#include <scai/lama/mic/MICCOOUtils.hpp>

// local library
#include <scai/lama/mic/MICUtils.hpp>

#include <scai/common/OpenMP.hpp>
#include <scai/lama/COOKernelTrait.hpp>

// internal scai libraries
#include <scai/hmemo/mic/MICSyncToken.hpp>
#include <scai/hmemo/mic/MICContext.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/ScalarType.hpp>

// std
#include <cmath>

namespace scai
{

using tasking::MICSyncToken;
using hmemo::MICContext;

namespace lama
{

__attribute__( ( target ( mic ) ) )
inline void getRange( IndexType& lb, IndexType& ub, IndexType n, PartitionId rank, PartitionId size )
{
    IndexType blockSize = ( n + size - 1 ) / size;
    lb = rank * blockSize;
    ub = ( rank + 1 ) * blockSize - 1;
    ub = std::min( ub, n - 1 );
}

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( MICCOOUtils::logger, "MIC.COOUtils" )

/* --------------------------------------------------------------------------- */
/*     Template implementations                                                */
/* --------------------------------------------------------------------------- */

void MICCOOUtils::getCSRSizes(
    IndexType csrSizes[],
    const IndexType numRows,
    const IndexType numValues,
    const IndexType cooIA[] )
{
    SCAI_LOG_INFO( logger, "get CSR sizes, #rows = " << numRows << ", #values = " << numValues )

    void* csrSizesPtr = csrSizes;
    const void* cooIAPtr = cooIA;

    // load distribution is done implicitly by block distribution of csrSizes

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( numRows, numValues, csrSizesPtr, cooIAPtr )

    {
        IndexType* csrSizes = (IndexType*) csrSizesPtr;
        const IndexType* cooIA = (const IndexType*) cooIAPtr;

        #pragma omp parallel
        {
            // initialize size array for each row

            #pragma omp for
            for( IndexType i = 0; i <= numRows; i++ )
            {
                csrSizes[i] = 0;
            }

            // increment size of a row for each used row value

            #pragma omp for

            for( IndexType k = 0; k < numValues; k++ )
            {
                IndexType i = cooIA[k];
                #pragma omp atomic
                csrSizes[i]++;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename COOValueType,typename CSRValueType>
void MICCOOUtils::getCSRValuesP(
    IndexType csrJA[],
    CSRValueType csrValues[],
    IndexType csrIA[],
    const IndexType numRows,
    const IndexType numValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const COOValueType cooValues[] )
{
    SCAI_LOG_ERROR( logger,
                    "get CSRValues<" << common::getScalarType<COOValueType>() << ", " << common::getScalarType<CSRValueType>() << ">" << ", #rows = " << numRows << ", #values = " << numValues )

    void* csrJAPtr = csrJA;
    void* csrValuesPtr = csrValues;
    void* csrIAPtr = csrIA;

    const void* cooIAPtr = cooIA;
    const void* cooJAPtr = cooJA;
    const void* cooValuesPtr = cooValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( numRows, numValues, csrJAPtr, csrValuesPtr, csrIAPtr, cooIAPtr, cooJAPtr, cooValuesPtr )
    {
        std::vector<IndexType> rowOffset( numRows ); // temp copy of csrIA

        IndexType* csrJA = (IndexType*) cooJAPtr;
        CSRValueType* csrValues = (CSRValueType*) csrValuesPtr;
        const IndexType* csrIA = (IndexType*) csrIAPtr;

        const IndexType* cooIA = (const IndexType*) cooIAPtr;
        const IndexType* cooJA = (const IndexType*) cooJAPtr;
        const COOValueType* cooValues = (const COOValueType*) cooValuesPtr;

        #pragma omp parallel
        {
            #pragma omp for

            for( IndexType i = 0; i < numRows; ++i )
            {
                rowOffset[i] = csrIA[i];
            }

            // traverse the non-zero values and put data at the right places

            #pragma omp for

            for ( IndexType k = 0; k < numValues; k++ )
            {
                IndexType i = cooIA[k];

                // better :  IndexType offset = __sync_fetch_and_add( &rowOffset[i], 1 );

                IndexType offset;

                #pragma omp critical
                {   
                    offset = rowOffset[i];
                    rowOffset[i]++;
                }

                csrJA[offset] = cooJA[k];
                csrValues[offset] = static_cast<CSRValueType>( cooValues[k] );
            }
        }
    }

    // ToDo: still some problems with diagonal property
}

/* --------------------------------------------------------------------------- */

template<typename COOValueType,typename CSRValueType>
void MICCOOUtils::getCSRValuesS(
    IndexType csrJA[],
    CSRValueType csrValues[],
    IndexType csrIA[],
    const IndexType numRows,
    const IndexType numValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const COOValueType cooValues[] )
{
    SCAI_LOG_ERROR( logger,
                    "get CSRValues<" << common::getScalarType<COOValueType>() << ", " << common::getScalarType<CSRValueType>() << ">" << ", #rows = " << numRows << ", #values = " << numValues )

    void* csrJAPtr = csrJA;
    void* csrValuesPtr = csrValues;
    void* csrIAPtr = csrIA;

    const void* cooIAPtr = cooIA;
    const void* cooJAPtr = cooJA;
    const void* cooValuesPtr = cooValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( numRows, numValues, csrJAPtr, csrValuesPtr, csrIAPtr, cooIAPtr, cooJAPtr, cooValuesPtr )
    {
        std::vector<IndexType> rowOffset( numRows ); // temp copy of csrIA

        IndexType* csrJA = (IndexType*) cooJAPtr;
        CSRValueType* csrValues = (CSRValueType*) csrValuesPtr;
        const IndexType* csrIA = (IndexType*) csrIAPtr;

        const IndexType* cooIA = (const IndexType*) cooIAPtr;
        const IndexType* cooJA = (const IndexType*) cooJAPtr;
        const COOValueType* cooValues = (const COOValueType*) cooValuesPtr;

        for( IndexType i = 0; i < numRows; ++i )
        {
            rowOffset[i] = csrIA[i];
        }

        // traverse the non-zero values and put data at the right places

        for ( IndexType k = 0; k < numValues; k++ )
        {
            IndexType i = cooIA[k];

            // better :  IndexType offset = __sync_fetch_and_add( &rowOffset[i], 1 );

            IndexType offset = rowOffset[i]++;
            csrJA[offset] = cooJA[k];
            csrValues[offset] = static_cast<CSRValueType>( cooValues[k] );
        }
    }
}

/* --------------------------------------------------------------------------- */

void MICCOOUtils::offsets2ia(
    IndexType cooIA[],
    const IndexType numValues,
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType numDiagonals )
{
    SCAI_LOG_INFO( logger,
                   "build cooIA( " << numValues << " ) from csrIA( " << ( numRows + 1 ) << " ), #diagonals = " << numDiagonals )

    const void* csrIAPtr = csrIA;
    void* cooIAPtr = cooIA;

    // parallel execution only possible if we have no separate diagonal elements
    // or if CSR data has diagonal property

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( numRows, numDiagonals, csrIAPtr, cooIAPtr )
    {
        const IndexType* csrIA = static_cast<const IndexType*>( csrIAPtr );

        IndexType* cooIA = static_cast<IndexType*>( cooIAPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; ++i )
        {
            IndexType csrOffset = csrIA[i];
            IndexType cooOffset = 0; // additional offset due to diagonals

            if( i < numDiagonals )
            {
                // diagonal elements will be the first nrows entries

                cooIA[i] = i;
                csrOffset += 1; // do not fill diagonal element again
                cooOffset = numDiagonals - i - 1; // offset in coo moves
            }

            // now fill remaining part of row i

            for( IndexType jj = csrOffset; jj < csrIA[i + 1]; ++jj )
            {
                cooIA[jj + cooOffset] = i;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename COOValueType,typename CSRValueType>
void MICCOOUtils::setCSRData(
    COOValueType cooValues[],
    const CSRValueType csrValues[],
    const IndexType numValues,
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType numDiagonals )
{
    SCAI_LOG_INFO( logger,
                   "build cooValues( << " << numValues << " from csrValues + csrIA( " << ( numRows + 1 ) << " ), #diagonals = " << numDiagonals )

    const void* csrValuesPtr = csrValues;
    const void* csrIAPtr = csrIA;

    void* cooValuesPtr = cooValues;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ) in( numRows, numDiagonals, csrIAPtr, csrValuesPtr, cooValuesPtr )
    {
        const CSRValueType* csrValues = static_cast<const CSRValueType*>( csrValuesPtr );
        const IndexType* csrIA = static_cast<const IndexType*>( csrIAPtr );

        COOValueType* cooValues = static_cast<COOValueType*>( cooValuesPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; ++i )
        {
            IndexType csrOffset = csrIA[i];
            IndexType cooOffset = 0; // additional offset due to diagonals

            if( i < numDiagonals )
            {
                // diagonal elements become the first 'numDiagonal' entries

                cooValues[i] = static_cast<COOValueType>( csrValues[csrOffset] );

                csrOffset += 1; // do not fill diagonal element again
                cooOffset = numDiagonals - i - 1; // offset in coo moves
            }

            // now fill remaining part of row i

            for( IndexType jj = csrOffset; jj < csrIA[i + 1]; ++jj )
            {
                cooValues[jj + cooOffset] = static_cast<COOValueType>( csrValues[jj] );
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICCOOUtils::normalGEMV(
    ValueType result[],
    const ValueType alpha,
    const ValueType x[],
    const ValueType beta,
    const ValueType y[],
    const IndexType numRows,
    const IndexType numValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[] )
{
    // SCAI_REGION( "MIC.COO.normalGEMV" )

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << common::getScalarType<ValueType>() << ">, result[" << numRows << "] = " 
                    << alpha << " * A( coo, #vals = " << numValues << " ) * x + " << beta << " * y " )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
    }

    // result := alpha * A * x + beta * y -> result:= beta * y; result += alpha * A

    MICUtils::setScale( result, beta, y, numRows );

    void* resultPtr = result;
    const void* cooIAPtr = cooIA;
    const void* cooJAPtr = cooJA;
    const void* cooValuesPtr = cooValues;
    const void* xPtr = x;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic : device ), in( numRows, numValues, resultPtr, alpha, cooIAPtr, cooJAPtr, cooValuesPtr, xPtr )
    {
        ValueType* result = static_cast<ValueType*>( resultPtr );
        const ValueType* x = static_cast<const ValueType*>( xPtr );

        const IndexType* cooIA = static_cast<const IndexType*>( cooIAPtr );
        const IndexType* cooJA = static_cast<const IndexType*>( cooJAPtr );
        const ValueType* cooValues = static_cast<const ValueType*>( cooValuesPtr );

        // Atomic update decreases performance by a factor of 4

        #pragma omp parallel
        {
            #pragma omp for

            for( IndexType k = 0; k < numValues; ++k )
            {
                IndexType i = cooIA[k];

                const ValueType val = alpha * cooValues[k] * x[cooJA[k]];

                #pragma omp atomic
                result[i] += val;

            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MICCOOUtils::jacobi(
    ValueType* solution,
    const IndexType cooNumValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega,
    const IndexType numRows )
{
    // SCAI_REGION( "MIC.COO.jacobi" )

    SCAI_LOG_INFO( logger,
                   "jacobi<" << common::getScalarType<ValueType>() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
    }

    // solution = omega * ( rhs - B * oldSolution ) * dinv + ( 1 - omega * oldSolution
    // done in two steps
    // solution = omega * rhs * dinv + ( 1 - omega * oldSolution
    // solution -= omega * B * oldSolution * dinv

    void* solutionPtr = solution;
    const void* cooIAPtr = cooIA;
    const void* cooJAPtr = cooJA;
    const void* cooValuesPtr = cooValues;
    const void* oldSolutionPtr = oldSolution;
    const void* rhsPtr = rhs;

    int device = MICContext::getCurrentDevice();

#pragma offload target( mic ), in( numRows, cooNumValues, solutionPtr, cooIAPtr, cooJAPtr, cooValuesPtr, \
                                       rhsPtr, oldSolutionPtr, omega )
    {
        ValueType* solution = static_cast<ValueType*>( solutionPtr );
        const ValueType* oldSolution = static_cast<const ValueType*>( oldSolutionPtr );
        const ValueType* rhs = static_cast<const ValueType*>( rhsPtr );

        const IndexType* cooIA = static_cast<const IndexType*>( cooIAPtr );
        const IndexType* cooJA = static_cast<const IndexType*>( cooJAPtr );
        const ValueType* cooValues = static_cast<const ValueType*>( cooValuesPtr );

        #pragma omp parallel for

        for( IndexType i = 0; i < numRows; ++i )
        {
            solution[i] = omega * rhs[i] / cooValues[i] + ( static_cast<ValueType>(1.0) - omega ) * oldSolution[i];
        }

        #pragma omp parallel for

        for( IndexType k = numRows; k < cooNumValues; ++k )
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

void MICCOOUtils::registerKernels( bool deleteFlag )
{
    SCAI_LOG_INFO( logger, "register COO kernels for MIC in Kernel Registry" )

    using kregistry::KernelRegistry;
    using common::context::MIC;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD ;   // add it or delete it

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

    KernelRegistry::set<COOKernelTrait::offsets2ia>( offsets2ia, MIC, flag );

    KernelRegistry::set<COOKernelTrait::setCSRData<IndexType, IndexType> >( setCSRData, MIC, flag );

    KernelRegistry::set<COOKernelTrait::setCSRData<float, float> >( setCSRData, MIC, flag );
    KernelRegistry::set<COOKernelTrait::setCSRData<float, double> >( setCSRData, MIC, flag );
    KernelRegistry::set<COOKernelTrait::setCSRData<double, float> >( setCSRData, MIC, flag );
    KernelRegistry::set<COOKernelTrait::setCSRData<double, double> >( setCSRData, MIC, flag );

    KernelRegistry::set<COOKernelTrait::getCSRSizes>( getCSRSizes, MIC, flag );

    // ToDo: routine does not work yet

    // KernelRegistry::set<COOKernelTrait::getCSRValues<float, float> >( getCSRValuesS, MIC, flag );
    // KernelRegistry::set<COOKernelTrait::getCSRValues<float, double> >( getCSRValuesS, MIC, flag );
    // KernelRegistry::set<COOKernelTrait::getCSRValues<double, float> >( getCSRValuesS, MIC, flag );
    // KernelRegistry::set<COOKernelTrait::getCSRValues<double, double> >( getCSRValuesS, MIC, flag );

    KernelRegistry::set<COOKernelTrait::normalGEMV<float> >( normalGEMV, MIC, flag );
    KernelRegistry::set<COOKernelTrait::normalGEMV<double> >( normalGEMV, MIC, flag );

    // ToDo: jacobi does not work yet

    // KernelRegistry::set<COOKernelTrait::jacobi<float> >( jacobi, MIC, flag );
    // KernelRegistry::set<COOKernelTrait::jacobi<double> >( jacobi, MIC, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICCOOUtils::RegisterGuard::RegisterGuard()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

MICCOOUtils::RegisterGuard::~RegisterGuard()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

MICCOOUtils::RegisterGuard MICCOOUtils::guard;    // guard variable for registration

} /* end namespace lama */

} /* end namespace scai */
