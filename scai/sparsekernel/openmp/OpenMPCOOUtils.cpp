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
#include <scai/sparsekernel/openmp/OpenMPCOOUtils.hpp>

// local library
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/sparsekernel/COOKernelTrait.hpp>

// internal scai libraries
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/tracing.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/OpenMP.hpp>
#include <scai/common/bind.hpp>
#include <scai/common/mepr/Container.hpp>


namespace scai
{

using common::TypeTraits;
using tasking::TaskSyncToken;

namespace sparsekernel
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

bool OpenMPCOOUtils::hasDiagonalProperty(
    const IndexType cooIA[],
    const IndexType cooJA[],
    const IndexType n )
{
    bool diagonalProperty = true;

    // The diagonal property is given if the first n entries
    // are the diagonal elements

    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

    for ( IndexType i = 0; i < n; ++i )
    {
        if ( !diagonalProperty )
        {
            continue;
        }

        if ( cooIA[i] != i || cooJA[i] != i )
        {
            diagonalProperty = false;
        }
    }

    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

template<typename COOValueType,typename CSRValueType>
void OpenMPCOOUtils::getCSRValues( IndexType csrJA[], 
                                   CSRValueType csrValues[], 
                                   IndexType csrIA[], // used as tmp, remains unchanged
                                   const IndexType numRows,
                                   const IndexType numValues,
                                   const IndexType cooIA[],
                                   const IndexType cooJA[],
                                   const COOValueType cooValues[] )
{
    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << TypeTraits<COOValueType>::id() << ", "
                    << TypeTraits<CSRValueType>::id() << ">" << ", #rows = " << numRows << ", #values = " << numValues )

    // traverse the non-zero values and put data at the right places
    // each thread reads all COO indexes but takes care for a certain range of CSR rows

    #pragma omp parallel
    {
        // get thread rank, size

        PartitionId rank = omp_get_thread_num();
        PartitionId nthreads = omp_get_num_threads();

        // compute range for which this thread is responsbile

        IndexType blockSize = ( numRows + nthreads - 1 ) / nthreads;
        IndexType lb = rank * blockSize;
        IndexType ub = ( rank + 1 ) * blockSize - 1;
        ub = std::min( ub, numRows - 1 );

        for ( IndexType k = 0; k < numValues; k++ )
        {
            IndexType i = cooIA[k];
    
            if ( i < lb || i > ub ) 
            {
                continue;
            }

            IndexType& offset = csrIA[i];
    
            csrJA[offset] = cooJA[k];
            csrValues[offset] = static_cast<CSRValueType>( cooValues[k] );
    
            SCAI_LOG_DEBUG( logger, "row " << i << ": new offset = " << offset )
    
            offset++;
        }
    }

    // set back the old offsets in csrIA

    for ( IndexType i = numRows; i > 0; --i )
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

template<typename COOValueType,typename OtherValueType>
void OpenMPCOOUtils::scaleRows( 
    COOValueType cooValues[],
    const OtherValueType rowValues[],
    const IndexType cooIA[],
    const IndexType numValues )
{
    SCAI_LOG_INFO( logger, "scaleRows in COO format" )

    #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

    for ( IndexType i = 0; i < numValues; ++i )
    {
        cooValues[i] *= static_cast<COOValueType>( rowValues[cooIA[i]] );
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
void OpenMPCOOUtils::normalGEMV_a(
    ValueType result[],
    const std::pair<ValueType, const ValueType*> ax,
    const std::pair<ValueType, const ValueType*> by,
    const IndexType numRows,
    const IndexType numValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[] )
{
    normalGEMV( result, ax.first, ax.second, by.first, by.second,
                numRows, numValues, cooIA, cooJA, cooValues );
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
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        // bind has limited number of arguments, so take help routine for call

        SCAI_LOG_INFO( logger,
                       "normalGEMV<" << TypeTraits<ValueType>::id() << "> launch it asynchronously" )

        syncToken->run( common::bind( normalGEMV_a<ValueType>,
                                      result,
                                      std::pair<ValueType, const ValueType*>( alpha, x ),
                                      std::pair<ValueType, const ValueType*>( beta, y ),
                                      numRows, numValues, cooIA, cooJA, cooValues ) );
        return;
    }

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads() << ">," 
                    << " result[" << numRows << "] = " << alpha << " * A( coo, #vals = " << numValues << " ) * x + " << beta << " * y " )

    // result := alpha * A * x + beta * y -> result:= beta * y; result += alpha * A

    utilskernel::OpenMPUtils::setScale( result, beta, y, numRows );

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
                   "normalGEMV<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads() << ">, result[" << numColumns << "] = " << alpha << " * A( coo, #vals = " << numValues << " ) * x + " << beta << " * y " )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution not supported here, do it by a task" )
    }

    // result := alpha * x * A + beta * y -> result:= beta * y; result += alpha * x * A

    utilskernel::OpenMPUtils::setScale( result, beta, y, numColumns );

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
                   "jacobi<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", omega = " << omega )

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

void OpenMPCOOUtils::Registrator::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using common::context::Host;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register COOUtils OpenMP-routines for Host at kernel registry [" << flag << "]" )

    KernelRegistry::set<COOKernelTrait::hasDiagonalProperty>( OpenMPCOOUtils::hasDiagonalProperty, Host, flag );
    KernelRegistry::set<COOKernelTrait::offsets2ia>( OpenMPCOOUtils::offsets2ia, Host, flag );
    KernelRegistry::set<COOKernelTrait::getCSRSizes>( OpenMPCOOUtils::getCSRSizes, Host, flag );
    KernelRegistry::set<COOKernelTrait::setCSRData<IndexType, IndexType> >( OpenMPCOOUtils::setCSRData, Host, flag );
}

template<typename ValueType>
void OpenMPCOOUtils::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using common::context::Host;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register COOUtils OpenMP-routines for Host at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<COOKernelTrait::normalGEMV<ValueType> >( OpenMPCOOUtils::normalGEMV, Host, flag );
    KernelRegistry::set<COOKernelTrait::normalGEVM<ValueType> >( OpenMPCOOUtils::normalGEVM, Host, flag );
    KernelRegistry::set<COOKernelTrait::jacobi<ValueType> >( OpenMPCOOUtils::jacobi, Host, flag );
}

template<typename ValueType, typename OtherValueType>
void OpenMPCOOUtils::RegistratorVO<ValueType, OtherValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using common::context::Host;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register COOUtils OpenMP-routines for Host at kernel registry [" << flag
        << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )

    KernelRegistry::set<COOKernelTrait::setCSRData<ValueType, OtherValueType> >( setCSRData, Host, flag );
    KernelRegistry::set<COOKernelTrait::getCSRValues<ValueType, OtherValueType> >( getCSRValues, Host, flag );
    KernelRegistry::set<COOKernelTrait::scaleRows<ValueType, OtherValueType> >( scaleRows, Host, flag );
}


/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPCOOUtils::OpenMPCOOUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;

    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, ARITHMETIC_HOST_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, ARITHMETIC_HOST_LIST, ARITHMETIC_HOST_LIST>::call( flag );
}

OpenMPCOOUtils::~OpenMPCOOUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    Registrator::initAndReg( flag );
    kregistry::mepr::RegistratorV<RegistratorV, ARITHMETIC_HOST_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, ARITHMETIC_HOST_LIST, ARITHMETIC_HOST_LIST>::call( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPCOOUtils OpenMPCOOUtils::guard;

} /* end namespace sparsekernel */

} /* end namespace scai */
