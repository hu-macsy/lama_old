/**
 * @file OpenMPDIAUtils.cpp
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
 * @brief Implementation of DIA utilities with OpenMP
 * @author Thomas Brandes
 * @date 04.07.2012
 */


// local library
#include <scai/sparsekernel/openmp/OpenMPDIAUtils.hpp>
#include <scai/sparsekernel/DIAKernelTrait.hpp>

// internal scai libraries
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/tracing.hpp>
#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/OpenMP.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/Constants.hpp>

// std
#include <cmath>
#include <functional>

namespace scai
{

namespace sparsekernel
{

using common::TypeTraits;
using common::Math;
using tasking::TaskSyncToken;

SCAI_LOG_DEF_LOGGER( OpenMPDIAUtils::logger, "OpenMP.DIAUtils" )

/* --------------------------------------------------------------------------- */
/*   Implementation of methods                                                 */
/* --------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------- */

IndexType OpenMPDIAUtils::getValuePos( const IndexType i,
                                       const IndexType j,
                                       const IndexType numRows,
                                       const IndexType diaOffsets[],
                                       const IndexType numDiagonals )
{
    IndexType pos = invalidIndex;

    for ( IndexType d = 0; d < numDiagonals; ++d )
    {
        if ( i + diaOffsets[d] == j )
        {
            pos = i + d * numRows;
            break;
        }
    }

    return pos;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType OpenMPDIAUtils::absMaxVal(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[] )
{
    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    RealType maxValue = 0;

    #pragma omp parallel
    {
        RealType threadVal = 0;

        #pragma omp for 

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType d = 0; d < numDiagonals; ++d )
            {
                const IndexType j = i + diaOffsets[d];

                if ( ! common::Utils::validIndex( j, numColumns ) )
                {
                    continue;
                }

                const RealType val = common::Math::abs( diaValues[i + d * numRows] );

                if ( val > threadVal )
                {
                    threadVal = val;
                }
            }
        }

        #pragma omp critical
        {
            SCAI_LOG_DEBUG( logger, "absMaxVal, threadVal = " << threadVal << ", maxVal = " << maxValue )

            if ( threadVal > maxValue )
            {
                maxValue = threadVal;
            }
        }
    }
    return maxValue;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDIAUtils::getCSRValues(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[] )
{
    const ValueType ZERO = 0;

    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << TypeTraits<ValueType>::id() << ">" 
                   << ", #rows = " << numRows << ", #diagonals = " << numDiagonals
                   << ", #non-zero values = " << csrIA[numRows] )

    // we cannot check for correct sizes, but at least for valid pointers

    if ( numDiagonals == 0 )
    {
        SCAI_ASSERT_EQUAL_DEBUG( 0, csrIA[numRows] )
        return;
    }

    if ( numDiagonals > 0 )
    {
        SCAI_ASSERT_DEBUG( diaOffsets != NULL, "offset array of DIA data is NULL" )

        if ( numRows > 0 )
        {
            SCAI_ASSERT_DEBUG( diaValues != NULL, "value array of DIA data is NULL" )
        }
    }

    // go through the DIA the same way again and copy the non-zeros

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.DIA.getCSR" )
        #pragma omp for 
        for ( IndexType i = 0; i < numRows; i++ )
        {
            IndexType offset = csrIA[i];

            for ( IndexType ii = 0; ii < numDiagonals; ii++ )
            {
                IndexType j = i + diaOffsets[ii];

                if ( !common::Utils::validIndex( j, numColumns ) )
                {
                    continue;
                }

                const ValueType value = diaValues[i + ii * numRows];

                if ( value == ZERO )
                {
                    continue;
                }

                csrJA[offset] = j;
                csrValues[offset] = value;
                SCAI_LOG_TRACE( logger,
                                "csrJA[" << offset << "] = " << csrJA[offset] << ", csrValues[" << offset << "] = " << csrValues[offset] )
                offset++;
            }

            SCAI_ASSERT_EQUAL_DEBUG( offset, csrIA[i + 1] )
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename DIAValueType>
void OpenMPDIAUtils::getCSRSizes(
    IndexType csrSizes[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const DIAValueType diaValues[] )
{
    SCAI_LOG_INFO( logger,
                   "get CSRSizes<" << TypeTraits<DIAValueType>::id() << "> for DIA matrix " << numRows << " x " << numColumns
                   << ", #diagonals = " << numDiagonals )

    const DIAValueType ZERO = 0;

    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; i++ )
    {
        IndexType count = 0;

        for ( IndexType ii = 0; ii < numDiagonals; ii++ )
        {
            IndexType j = i + diaOffsets[ii]; // column index

            if ( !common::Utils::validIndex( j, numColumns ) )
            {
                continue;
            }

            if ( diaValues[i + ii * numRows] != ZERO )
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
    const common::MatrixOp op )
{
    IndexType nResult = common::isTranspose( op ) ? numColumns : numRows;

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads()
                   << ", op = " << op 
                   << ">, result[" << nResult << "] = " 
                   << alpha << " * A( dia, #diags = " << numDiagonals << " ) * x + " << beta << " * y ( " << y << " )" )

    // Note: launching thread gets token here, asynchronous thread will get NULL

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( beta != common::Constants::ZERO && nResult > 0 )
    {
        SCAI_ASSERT_ERROR( y != NULL, "y = " << y << " is NULL, not allowed for beta = " << beta )
    }

    if ( syncToken )
    {
        SCAI_LOG_INFO( logger,
                       "normalGEMV<" << TypeTraits<ValueType>::id() << "> launch it asynchronously" )

        syncToken->run( std::bind( normalGEMV<ValueType>, result,
                                   alpha, x, beta, y,
                                   numRows, numColumns, numDiagonals, diaOffsets, diaValues, op ) );
        return;
    }

    // result := alpha * A * x + beta * y 

    if ( op == common::MatrixOp::TRANSPOSE )
    {
        #pragma omp parallel
        {
            SCAI_REGION( "OpenMP.DIA.normalGEMV_t" )

            #pragma omp for

            for ( IndexType j = 0; j < numColumns; j++ )
            {
                ValueType accu = 0;

                for ( IndexType ii = 0; ii < numDiagonals; ++ii )
                {
                    const IndexType i = j - diaOffsets[ii];

                    if ( common::Utils::validIndex( i, numRows ) )
                    {
                        accu += diaValues[ii * numRows + i] * x[i];
                    }
                }
     
                if ( beta == 0 )
                {
                    result[j] = alpha *accu;
                }
                else
                {
                    result[j] = alpha * accu + beta * y[j];
                }
            }
        }
    }
    else if ( op == common::MatrixOp::NORMAL )
    {
        #pragma omp parallel
        {
            SCAI_REGION( "OpenMP.DIA.normalGEMV_n" )
    
            #pragma omp for 

            for ( IndexType i = 0; i < numRows; i++ )
            {
                ValueType accu = 0;
    
                for ( IndexType ii = 0; ii < numDiagonals; ++ii )
                {
                    const IndexType j = i + diaOffsets[ii];
    
                    if ( common::Utils::validIndex( j, numColumns ) )
                    {
                        accu += diaValues[ii * numRows + i] * x[j];
                    }
                }

                if ( beta == 0 )
                {
                    result[i] = alpha *accu;
                }
                else
                {
                    result[i] = alpha * accu + beta * y[i];
                }
            }
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "unsupported matrix op: " << op )
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
    const IndexType numRows )
{
    SCAI_LOG_INFO( logger,
                   "jacobi<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", #cols = " << numColumns 
                    << ", #diagonals = " << numDiagonals << ", omega = " << omega )

    // main diagonal must be first
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken != NULL )
    {
        SCAI_LOG_ERROR( logger, "jacobi called asynchronously, not supported here" )
    }

    ValueType omega1 = 1;
    omega1 -= omega;

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.DIA.Jacobi" )
        #pragma omp for 

        for ( IndexType i = 0; i < numRows; i++ )
        {
            ValueType temp = rhs[i];
            ValueType diag = 0;

            for ( IndexType ii = 0; ii < numDiagonals; ++ii )
            {
                if ( diaOffset[ii] == 0 )
                {
                    diag = diaValues[ ii * numRows + i ];
                }
                else
                {
                    const IndexType j = i + diaOffset[ii];

                    if ( common::Utils::validIndex( j, numColumns ) )
                    {
                        temp -= diaValues[ii * numRows + i] * oldSolution[j];
                    }
                }
            }

            solution[i] = omega * ( temp / diag ) + omega1 * oldSolution[i];
        }
    }
}

/* --------------------------------------------------------------------------- */
/* Registrator classes, method registerKernels                                 */
/* --------------------------------------------------------------------------- */

void OpenMPDIAUtils::Registrator::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register DIAUtils OpenMP-routines for Host at kernel registry [" << flag << "]" )
    KernelRegistry::set<DIAKernelTrait::getValuePos>( getValuePos, ctx, flag );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDIAUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register DIAUtils OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<DIAKernelTrait::getCSRSizes<ValueType> >( getCSRSizes, ctx, flag );
    KernelRegistry::set<DIAKernelTrait::absMaxVal<ValueType> >( absMaxVal, ctx, flag );
    KernelRegistry::set<DIAKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<DIAKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<DIAKernelTrait::getCSRValues<ValueType> >( getCSRValues, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPDIAUtils::OpenMPDIAUtils()
{
    SCAI_LOG_INFO( logger, "register DIAUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;

    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

OpenMPDIAUtils::~OpenMPDIAUtils()
{
    SCAI_LOG_INFO( logger, "unregister DIAUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    Registrator::registerKernels( flag );
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPDIAUtils OpenMPDIAUtils::guard;

} /* end namespace sparsekernel */

} /* end namespace scai */
