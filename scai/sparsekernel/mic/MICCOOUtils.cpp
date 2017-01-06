/**
 * @file MICCOOUtils.cpp
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
 * @brief Implementation of COO storage utilities with MIC
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// hpp
#include <scai/sparsekernel/mic/MICCOOUtils.hpp>

// local library
#include <scai/sparsekernel/COOKernelTrait.hpp>

// internal scai libraries
#include <scai/utilskernel/mic/MICUtils.hpp>
#include <scai/tasking/mic/MICSyncToken.hpp>
#include <scai/hmemo/mic/MICContext.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/OpenMP.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/OpenMP.hpp>

// std
#include <cmath>

namespace scai
{

using tasking::MICSyncToken;
using hmemo::MICContext;
using common::TypeTraits;
using utilskernel::MICUtils;

namespace sparsekernel
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

template<typename COOValueType, typename CSRValueType>
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
    SCAI_LOG_DEBUG( logger,
                    "get CSRValues<" << TypeTraits<COOValueType>::id() << ", " << TypeTraits<CSRValueType>::id() << ">"
                    << ", #rows = " << numRows << ", #values = " << numValues )
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
        IndexType* csrJA = ( IndexType* ) cooJAPtr;
        CSRValueType* csrValues = ( CSRValueType* ) csrValuesPtr;
        const IndexType* csrIA = ( IndexType* ) csrIAPtr;
        const IndexType* cooIA = ( const IndexType* ) cooIAPtr;
        const IndexType* cooJA = ( const IndexType* ) cooJAPtr;
        const COOValueType* cooValues = ( const COOValueType* ) cooValuesPtr;
        #pragma omp parallel
        {
            #pragma omp for

            for ( IndexType i = 0; i < numRows; ++i )
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

template<typename COOValueType, typename CSRValueType>
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
    SCAI_LOG_DEBUG( logger,
                    "get CSRValues<" << TypeTraits<COOValueType>::id() << ", " << TypeTraits<CSRValueType>::id() << ">"
                    << ", #rows = " << numRows << ", #values = " << numValues )
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

        IndexType* csrJA = reinterpret_cast<IndexType*>( cooJAPtr );
        CSRValueType* csrValues = reinterpret_cast<CSRValueType*>( csrValuesPtr );
        const IndexType* csrIA = ( IndexType* ) csrIAPtr;
        const IndexType* cooIA = ( const IndexType* ) cooIAPtr;
        const IndexType* cooJA = ( const IndexType* ) cooJAPtr;
        const COOValueType* cooValues = ( const COOValueType* ) cooValuesPtr;

        for ( IndexType i = 0; i < numRows; ++i )
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

        for ( IndexType i = 0; i < numRows; ++i )
        {
            IndexType csrOffset = csrIA[i];
            IndexType cooOffset = 0; // additional offset due to diagonals

            if ( i < numDiagonals )
            {
                // diagonal elements will be the first nrows entries
                cooIA[i] = i;
                csrOffset += 1; // do not fill diagonal element again
                cooOffset = numDiagonals - i - 1; // offset in coo moves
            }

            // now fill remaining part of row i

            for ( IndexType jj = csrOffset; jj < csrIA[i + 1]; ++jj )
            {
                cooIA[jj + cooOffset] = i;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename COOValueType, typename CSRValueType>
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

        for ( IndexType i = 0; i < numRows; ++i )
        {
            IndexType csrOffset = csrIA[i];
            IndexType cooOffset = 0; // additional offset due to diagonals

            if ( i < numDiagonals )
            {
                // diagonal elements become the first 'numDiagonal' entries
                cooValues[i] = static_cast<COOValueType>( csrValues[csrOffset] );
                csrOffset += 1; // do not fill diagonal element again
                cooOffset = numDiagonals - i - 1; // offset in coo moves
            }

            // now fill remaining part of row i

            for ( IndexType jj = csrOffset; jj < csrIA[i + 1]; ++jj )
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
                   "normalGEMV<" << TypeTraits<ValueType>::id() << ">, result[" << numRows << "] = "
                   << alpha << " * A( coo, #vals = " << numValues << " ) * x + " << beta << " * y " )
    MICSyncToken* syncToken = MICSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger, "asynchronous execution for for MIC not supported yet" )
    }

    // result := alpha * A * x + beta * y -> result:= beta * y; result += alpha * A
    MICUtils::binaryOpScalar1( result, beta, y, numRows, utilskernel::binary::MULT );
    void* resultPtr = result;
    const void* cooIAPtr = cooIA;
    const void* cooJAPtr = cooJA;
    const void* cooValuesPtr = cooValues;
    const void* xPtr = x;
    const ValueType* alphaPtr = &alpha;
    int device = MICContext::getCurrentDevice();
#pragma offload target( mic : device ), in( numRows, numValues, resultPtr, alphaPtr[0:1], cooIAPtr, cooJAPtr, cooValuesPtr, xPtr )
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

            for ( IndexType k = 0; k < numValues; ++k )
            {
                IndexType i = cooIA[k];
                const ValueType val = *alphaPtr * cooValues[k] * x[cooJA[k]];
                #pragma omp critical
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
                   "jacobi<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", omega = " << omega )
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
        for ( IndexType i = 0; i < numRows; ++i )
        {
            solution[i] = omega * rhs[i] / cooValues[i] + ( static_cast<ValueType>( 1.0 ) - omega ) * oldSolution[i];
        }

        #pragma omp parallel for
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

void MICCOOUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_DEBUG( logger, "register[flag=" << flag << "]: untyped routines" )

    KernelRegistry::set<COOKernelTrait::offsets2ia>( offsets2ia, ctx, flag );
    KernelRegistry::set<COOKernelTrait::setCSRData<IndexType, IndexType> >( setCSRData, ctx, flag );
}

template<typename ValueType>
void MICCOOUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_DEBUG( logger, "register[flag=" << flag << "]: T = " << common::TypeTraits<ValueType>::id() )

    KernelRegistry::set<COOKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void MICCOOUtils::RegistratorVO<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::MIC;

    SCAI_LOG_DEBUG( logger, "register[flag=" << flag << "]: TT " <<
                            common::TypeTraits<ValueType>::id() << ", " << common::TypeTraits<OtherValueType>::id() )

    KernelRegistry::set<COOKernelTrait::setCSRData<ValueType, OtherValueType> >( setCSRData, ctx, flag );

    // ToDo: routine does not work yet
    // KernelRegistry::set<COOKernelTrait::getCSRValues<float, float> >( getCSRValuesS, ctx, flag );
    // KernelRegistry::set<COOKernelTrait::getCSRValues<float, double> >( getCSRValuesS, ctx, flag );
    // KernelRegistry::set<COOKernelTrait::getCSRValues<double, float> >( getCSRValuesS, ctx, flag );
    // KernelRegistry::set<COOKernelTrait::getCSRValues<double, double> >( getCSRValuesS, ctx, flag );
    // ToDo: jacobi does not work yet
    // KernelRegistry::set<COOKernelTrait::jacobi<float> >( jacobi, ctx, flag );
    // KernelRegistry::set<COOKernelTrait::jacobi<double> >( jacobi, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Static initialization with registration                                  */
/* --------------------------------------------------------------------------- */

MICCOOUtils::RegisterGuard::RegisterGuard()
{
    SCAI_LOG_INFO( logger, "register COOUtils routines for MIC(OpenMP,offload) at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_MIC_LIST, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels( flag );
}

MICCOOUtils::RegisterGuard::~RegisterGuard()
{
    SCAI_LOG_INFO( logger, "unregister COOUtils routines for MIC(OpenMP,offload) at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_MIC_LIST, SCAI_NUMERIC_TYPES_MIC_LIST>::registerKernels( flag );
}

MICCOOUtils::RegisterGuard MICCOOUtils::guard;    // guard variable for registration

} /* end namespace sparsekernel */

} /* end namespace scai */
