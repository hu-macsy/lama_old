/**
 * @file OpenMPCOOUtils.cpp
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
 * @brief Implementation of COO storage utilities with OpenMP
 * @author Thomas Brandes
 * @date 02.07.2012
 * @since 1.0.0
 */

// hpp
#include <scai/lama/openmp/OpenMPCOOUtils.hpp>

// local library
#include <scai/lama/openmp/OpenMPUtils.hpp>

#include <scai/lama/COOKernelTrait.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tracing.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/OpenMP.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

// boost
#include <boost/preprocessor.hpp>

namespace scai
{

using common::getScalarType;
using tasking::TaskSyncToken;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( OpenMPCOOUtils::logger, "OpenMP.COOUtils" )

/* --------------------------------------------------------------------------- */
/*     Template implementations                                                */
/* --------------------------------------------------------------------------- */

void OpenMPCOOUtils::getCSRSizes(
    IndexType csrSizes[],
    const IndexType numRows,
    const IndexType numValues,
    const IndexType cooIA[] )
{
    SCAI_LOG_INFO( logger, "get CSR sizes, #rows = " << numRows << ", #values = " << numValues )

    // initialize size array for each row

    for( IndexType i = 0; i < numRows; i++ )
    {
        csrSizes[i] = 0;
    }

    // increment size of a row for each used row value

    for( IndexType k = 0; k < numValues; k++ )
    {
        IndexType i = cooIA[k];
        SCAI_ASSERT_DEBUG( i < numRows, "cooIA[" << k << "] = " << i << " out of range, #rows = " << numRows )
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
    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << getScalarType<COOValueType>() << ", " 
                    << getScalarType<CSRValueType>() << ">" << ", #rows = " << numRows << ", #values = " << numValues )

    // traverse the non-zero values and put data at the right places

    for( IndexType k = 0; k < numValues; k++ )
    {
        IndexType i = cooIA[k];

        IndexType& offset = csrIA[i];

        csrJA[offset] = cooJA[k];
        csrValues[offset] = static_cast<CSRValueType>( cooValues[k] );

        SCAI_LOG_DEBUG( logger, "row " << i << ": new offset = " << offset )

        offset++;
    }

    // set back the old offsets in csrIA

    for( IndexType i = numRows; i > 0; --i )
    {
        csrIA[i] = csrIA[i - 1];
    }

    csrIA[0] = 0;

    SCAI_ASSERT_EQUAL_DEBUG( csrIA[numRows], numValues )
}

/* --------------------------------------------------------------------------- */

void OpenMPCOOUtils::offsets2ia(
    IndexType cooIA[],
    const IndexType numValues,
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType numDiagonals )
{
    SCAI_LOG_INFO( logger,
                   "build cooIA( " << numValues << " ) from csrIA( " << ( numRows + 1 ) << " ), #diagonals = " << numDiagonals )

    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

    for( IndexType i = 0; i < numRows; ++i )
    {
        IndexType csrOffset = csrIA[i];
        IndexType cooOffset = 0; // additional offset due to diagonals

        if( i < numDiagonals )
        {
            // make sure that we really have at least one element in this row

            SCAI_ASSERT_DEBUG( csrIA[i] < csrIA[i + 1], "no elem in row " << i );

            // diagonal elements will be the first nrows entries

            cooIA[i] = i;
            csrOffset += 1; // do not fill diagonal element again
            cooOffset = numDiagonals - i - 1; // offset in coo moves
        }

        // now fill remaining part of row i

        for( IndexType jj = csrOffset; jj < csrIA[i + 1]; ++jj )
        {
            SCAI_LOG_TRACE( logger,
                            "cooIA[ " << ( jj + cooOffset ) << "] = " << i << ", jj = " << jj << ", cooOffset = " << cooOffset )
            cooIA[jj + cooOffset] = i;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename COOValueType,typename CSRValueType>
void OpenMPCOOUtils::setCSRData(
    COOValueType cooValues[],
    const CSRValueType csrValues[],
    const IndexType numValues,
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType numDiagonals )
{
    SCAI_LOG_INFO( logger,
                   "build cooValues( << " << numValues << " from csrValues + csrIA( " << ( numRows + 1 ) << " ), #diagonals = " << numDiagonals )

    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

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

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPCOOUtils::normalGEMV(
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
    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << getScalarType<ValueType>() << ", #threads = " << omp_get_max_threads() << ">, result[" << numRows << "] = " << alpha << " * A( coo, #vals = " << numValues << " ) * x + " << beta << " * y " )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution not supported here, do it by a task" )
    }

    // result := alpha * A * x + beta * y -> result:= beta * y; result += alpha * A

    OpenMPUtils::setScale( result, beta, y, numRows );

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.COO.normalGEMV" )

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType k = 0; k < numValues; ++k )
        {
            IndexType i = cooIA[k];
            IndexType j = cooJA[k];

            // we must use atomic updates as different threads might update same row i

            const ValueType resultUpdate = alpha * cooValues[k] * x[j];

            atomicAdd( result[i], resultUpdate );
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
    const IndexType numColumns,
    const IndexType numValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[] )
{
    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << getScalarType<ValueType>() << ", #threads = " << omp_get_max_threads() << ">, result[" << numColumns << "] = " << alpha << " * A( coo, #vals = " << numValues << " ) * x + " << beta << " * y " )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution not supported here, do it by a task" )
    }

    // result := alpha * x * A + beta * y -> result:= beta * y; result += alpha * x * A

    OpenMPUtils::setScale( result, beta, y, numColumns );

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.COO.normalGEMV" )

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType k = 0; k < numValues; ++k )
        {
            IndexType i = cooIA[k];
            IndexType j = cooJA[k];

            // we must use atomic updates as different threads might update same row i

            const ValueType resultUpdate = alpha * cooValues[k] * x[i];

            // thread-safe atomic update

            atomicAdd( result[j], resultUpdate );
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
    const IndexType numRows )
{
    SCAI_LOG_INFO( logger,
                   "jacobi<" << getScalarType<ValueType>() << ">" << ", #rows = " << numRows << ", omega = " << omega )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    // solution = omega * ( rhs - B * oldSolution ) * dinv + ( 1 - omega * oldSolution
    // done in two steps
    // solution = omega * rhs * dinv + ( 1 - omega * oldSolution
    // solution -= omega * B * oldSolution * dinv

    #pragma omp parallel for

    for( IndexType i = 0; i < numRows; ++i )
    {
        solution[i] = omega * rhs[i] / cooValues[i] + ( static_cast<ValueType>(1.0) - omega ) * oldSolution[i];
    }

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.COO.jacobi" )

        #pragma omp for

        for( IndexType k = numRows; k < cooNumValues; ++k )
        {
            IndexType i = cooIA[k];
            IndexType j = cooJA[k];

            // we must use atomic updates as different threads might update same row i

            const ValueType update = -omega * cooValues[k] * oldSolution[j] / cooValues[i];

            atomicAdd( solution[i], update );
        }
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPCOOUtils::registerKernels()
{
    using namespace scai::kregistry;

    // ctx will contain the context for which registration is done, here Host

    common::context::ContextType ctx = common::context::Host;

    KernelRegistry::set<COOKernelTrait::offsets2ia>( offsets2ia, ctx );
    KernelRegistry::set<COOKernelTrait::getCSRSizes>( getCSRSizes, ctx );

    KernelRegistry::set<COOKernelTrait::setCSRData<IndexType, IndexType> >( setCSRData, ctx );

#define LAMA_COO_UTILS2_REGISTER(z, J, TYPE )                                                                     \
    KernelRegistry::set<COOKernelTrait::setCSRData<TYPE, ARITHMETIC_HOST_TYPE_##J> >( setCSRData, ctx );       \
    KernelRegistry::set<COOKernelTrait::getCSRValues<TYPE, ARITHMETIC_HOST_TYPE_##J> >( getCSRValues, ctx );   \

#define LAMA_COO_UTILS_REGISTER(z, I, _)                                                                \
    KernelRegistry::set<COOKernelTrait::normalGEMV<ARITHMETIC_HOST_TYPE_##I> >( normalGEMV, ctx );   \
    KernelRegistry::set<COOKernelTrait::normalGEVM<ARITHMETIC_HOST_TYPE_##I> >( normalGEVM, ctx );   \
    KernelRegistry::set<COOKernelTrait::jacobi<ARITHMETIC_HOST_TYPE_##I> >( jacobi, ctx );           \
                                                                                          \
    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT,                                            \
                     LAMA_COO_UTILS2_REGISTER,                                            \
                     ARITHMETIC_HOST_TYPE_##I )                                           \

    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_COO_UTILS_REGISTER, _ )

#undef LAMA_COO_UTILS_REGISTER
#undef LAMA_COO_UTILS2_REGISTER
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the Utils routines                                */
/* --------------------------------------------------------------------------- */

bool OpenMPCOOUtils::registerInterface()
{
    registerKernels();
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPCOOUtils::initialized = registerInterface();

} /* end namespace lama */

} /* end namespace scai */
