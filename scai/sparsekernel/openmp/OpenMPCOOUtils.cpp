/**
 * @file OpenMPCOOUtils.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Implementation of COO storage utilities with OpenMP
 * @author Thomas Brandes
 * @date 02.07.2012
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

namespace scai
{

using common::TypeTraits;
using tasking::TaskSyncToken;

namespace sparsekernel
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( OpenMPCOOUtils::logger, "OpenMP.COOUtils" )

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

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType csrOffset = csrIA[i];
        IndexType cooOffset = 0; // additional offset due to diagonals

        if ( i < numDiagonals )
        {
            // make sure that we really have at least one element in this row
            SCAI_ASSERT_DEBUG( csrIA[i] < csrIA[i + 1], "no elem in row " << i );
            // diagonal elements will be the first nrows entries
            cooIA[i] = i;
            csrOffset += 1; // do not fill diagonal element again
            cooOffset = numDiagonals - i - 1; // offset in coo moves
        }

        // now fill remaining part of row i

        for ( IndexType jj = csrOffset; jj < csrIA[i + 1]; ++jj )
        {
            SCAI_LOG_TRACE( logger,
                            "cooIA[ " << ( jj + cooOffset ) << "] = " << i << ", jj = " << jj << ", cooOffset = " << cooOffset )
            cooIA[jj + cooOffset] = i;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename COOValueType, typename OtherValueType>
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

template<typename COOValueType, typename CSRValueType>
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

        for ( IndexType k = 0; k < numValues; ++k )
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
void OpenMPCOOUtils::normalGEVM_a(
    ValueType result[],
    const std::pair<ValueType, const ValueType*> ax,
    const std::pair<ValueType, const ValueType*> by,
    const IndexType numColumns,
    const IndexType numValues,
    const IndexType cooIA[],
    const IndexType cooJA[],
    const ValueType cooValues[] )
{
    normalGEVM( result, ax.first, ax.second, by.first, by.second, numColumns, numValues, cooIA, cooJA, cooValues );
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
                   "normalGEVM<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads()
                   << ">, result[" << numColumns << "] = " << alpha << " * A( coo, #vals = " << numValues << " ) * x + "
                   << beta << " * y " )
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        // bind takes maximal 9 arguments, so we put (alpha, x) and (beta, y) in pair structs
        SCAI_LOG_INFO( logger, "normalGEVM<" << TypeTraits<ValueType>::id() << ", launch it as an asynchronous task" )
        syncToken->run( common::bind( normalGEVM_a<ValueType>,
                                      result,
                                      std::pair<ValueType, const ValueType*>( alpha, x ),
                                      std::pair<ValueType, const ValueType*>( beta, y ),
                                      numColumns, numValues, cooIA, cooJA, cooValues ) );
        return;
    }

    // result := alpha * x * A + beta * y -> result:= beta * y; result += alpha * x * A
    utilskernel::OpenMPUtils::setScale( result, beta, y, numColumns );
    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.COO.normalGEMV" )
        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for ( IndexType k = 0; k < numValues; ++k )
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

    if ( syncToken )
    {
        COMMON_THROWEXCEPTION( "asynchronous execution should be done by LAMATask before" )
    }

    // solution = omega * ( rhs - B * oldSolution ) * dinv + ( 1 - omega * oldSolution
    // done in two steps
    // solution = omega * rhs * dinv + ( 1 - omega * oldSolution
    // solution -= omega * B * oldSolution * dinv
    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        solution[i] = omega * rhs[i] / cooValues[i] + ( static_cast<ValueType>( 1.0 ) - omega ) * oldSolution[i];
    }

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.COO.jacobi" )
        #pragma omp for

        for ( IndexType k = numRows; k < cooNumValues; ++k )
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

void OpenMPCOOUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::context::ContextType ctx = common::context::Host;
    SCAI_LOG_DEBUG( logger, "register COOUtils OpenMP-routines for Host at kernel registry [" << flag << "]" )
    KernelRegistry::set<COOKernelTrait::hasDiagonalProperty>( OpenMPCOOUtils::hasDiagonalProperty, ctx, flag );
    KernelRegistry::set<COOKernelTrait::offsets2ia>( OpenMPCOOUtils::offsets2ia, ctx, flag );
    KernelRegistry::set<COOKernelTrait::setCSRData<IndexType, IndexType> >( OpenMPCOOUtils::setCSRData, ctx, flag );
}

template<typename ValueType>
void OpenMPCOOUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::context::ContextType ctx = common::context::Host;
    SCAI_LOG_DEBUG( logger, "register COOUtils OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<COOKernelTrait::normalGEMV<ValueType> >( OpenMPCOOUtils::normalGEMV, ctx, flag );
    KernelRegistry::set<COOKernelTrait::normalGEVM<ValueType> >( OpenMPCOOUtils::normalGEVM, ctx, flag );
    KernelRegistry::set<COOKernelTrait::jacobi<ValueType> >( OpenMPCOOUtils::jacobi, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void OpenMPCOOUtils::RegistratorVO<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::context::ContextType ctx = common::context::Host;
    SCAI_LOG_DEBUG( logger, "register COOUtils OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )
    KernelRegistry::set<COOKernelTrait::setCSRData<ValueType, OtherValueType> >( setCSRData, ctx, flag );
    KernelRegistry::set<COOKernelTrait::scaleRows<ValueType, OtherValueType> >( scaleRows, ctx, flag );
}


/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPCOOUtils::OpenMPCOOUtils()
{
    SCAI_LOG_INFO( logger, "register COOUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_HOST_LIST, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

OpenMPCOOUtils::~OpenMPCOOUtils()
{
    SCAI_LOG_INFO( logger, "unregister COOUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_HOST_LIST, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPCOOUtils OpenMPCOOUtils::guard;

} /* end namespace sparsekernel */

} /* end namespace scai */
