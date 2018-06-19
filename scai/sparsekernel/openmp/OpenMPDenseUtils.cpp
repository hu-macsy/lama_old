/**
 * @file OpenMPDenseUtils.cpp
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
 * @brief Implementation of Dense utilities with OpenMP
 * @author Thomas Brandes
 * @date 02.07.2012
 */

// hpp
#include <scai/sparsekernel/openmp/OpenMPDenseUtils.hpp>

// local library
#include <scai/sparsekernel/DenseKernelTrait.hpp>

#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/OpenMP.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

namespace scai
{

namespace sparsekernel
{

using common::TypeTraits;
using common::Math;
using tasking::TaskSyncToken;

SCAI_LOG_DEF_LOGGER( OpenMPDenseUtils::logger, "OpenMP.DenseUtils" )

/* --------------------------------------------------------------------------- */
/*     Template implementations                                                */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType OpenMPDenseUtils::nonZeroValues(
    const ValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const ValueType eps )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;
    RealType absEps = Math::abs( eps );

    SCAI_REGION( "OpenMP.DenseUtils.nonZeroValues" )

    IndexType count = 0;
    #pragma omp parallel for reduction( + : count )

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            if ( Math::abs( denseValues[denseindex( i, j, numRows, numColumns) ] ) > absEps )
            {
                count++;
            }
        }
    }
    return count;
}

template<typename ValueType>
void OpenMPDenseUtils::getCSRSizes(
    IndexType csrSizes[],
    const IndexType numRows,
    const IndexType numColumns,
    const ValueType denseValues[],
    const RealType<ValueType> eps )
{
    SCAI_REGION( "OpenMP.DenseUtils.getCSRSizes" )

    if ( numRows > 0 )
    {
        SCAI_ASSERT_DEBUG( csrSizes != NULL, "csrSizes is NULL" )

        if ( numColumns > 0 )
        {
            SCAI_ASSERT_DEBUG( denseValues != NULL, "denseValues is NULL" )
        }
    }

    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType nonZeros = 0; // count for each row in parallel

        for ( IndexType j = 0; j < numColumns; ++j )
        {
            const ValueType& value = denseValues[denseindex( i, j, numRows, numColumns )];

            if ( common::Math::abs( value ) > eps )
            {
                ++nonZeros;
            }
        }

        csrSizes[i] = nonZeros;
    }
}

/** Helper routine for conversion Dense to CSR */

template<typename ValueType>
void OpenMPDenseUtils::getCSRValues(
    IndexType csrJA[],
    ValueType csrValues[],
    const IndexType csrIA[],
    const IndexType numRows,
    const IndexType numColumns,
    const ValueType denseValues[],
    const RealType<ValueType> eps )
{
    typedef typename TypeTraits<ValueType>::RealType RealType;

    SCAI_REGION( "OpenMP.DenseUtils.getCSRValues" )

    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << TypeTraits<ValueType>::id() << ">" << ", size is " << numRows << " x " << numColumns )
    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType offset = csrIA[i];

        for ( IndexType j = 0; j < numColumns; ++j )
        {
            const ValueType& value = denseValues[denseindex( i, j, numRows, numColumns )];

            RealType absValue = common::Math::abs( value );

            if ( absValue > eps )
            {
                csrValues[offset] = value;
                csrJA[offset] = j;
                offset++;
            }
        }

        // verification that offset array was really a good one
        // check is not needed if non-zero values have been counted by getCSRSizes
        SCAI_ASSERT_EQUAL_DEBUG( offset, csrIA[i + 1] )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDenseUtils::setCSRValues(
    ValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const ValueType csrValues[] )
{
    SCAI_REGION( "OpenMP.DenseUtils.setCSRValues" )

    SCAI_LOG_INFO( logger,
                   "set CSRValues<" << TypeTraits<ValueType>::id() << ">" 
                   << ", size is " << numRows << " x " << numColumns )

    const ValueType ZERO = 0;

    // parallelization possible as offset array csrIA is available

    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; i++ )
    {
        // Initialize complete row with zero values
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            ValueType& elem = denseValues[denseindex( i, j, numRows, numColumns )];
            elem = ZERO;
        }

        // fill up positions for which non-zero values are given

        for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
        {
            const IndexType j = csrJA[jj];
            ValueType& elem = denseValues[denseindex( i, j, numRows, numColumns )];
            elem = csrValues[jj];
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDenseUtils::setValue(
    ValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const ValueType val,
    const common::BinaryOp op )
{
    SCAI_REGION( "OpenMP.DenseUtils.setValue" )

    using namespace common;

    switch ( op )
    {
        case BinaryOp::COPY :
        {
            // Parallel initialization very important for efficient  allocation
            #pragma omp parallel for 

            for ( IndexType i = 0; i < numRows; ++i )
            {
                for ( IndexType j = 0; j < numColumns; ++j )
                {
                    denseValues[denseindex( i, j, numRows, numColumns )] = val;
                }
            }

            break;
        }

        case BinaryOp::MULT :
        {
            // Parallel initialization very important for efficient  allocation

            #pragma omp parallel for 

            for ( IndexType i = 0; i < numRows; ++i )
            {
                for ( IndexType j = 0; j < numColumns; ++j )
                {
                    denseValues[denseindex( i, j, numRows, numColumns )] *= val;
                }
            }

            break;
        }

        default :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < numRows; ++i )
            {
                for ( IndexType j = 0; j < numColumns; ++j )
                {
                    ValueType& elem = denseValues[denseindex( i, j, numRows, numColumns )];
                    elem = applyBinary( elem, op, val );
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDenseUtils::setRows(
    ValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const ValueType rowValues[],
    common::BinaryOp op )
{
    SCAI_REGION( "OpenMP.DenseUtils.scaleRows" )

    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        const ValueType rowValue = rowValues[i];

        if ( common::BinaryOp::MULT == op )
        {
            // scale the whole row with this value

            for ( IndexType j = 0; j < numColumns; ++j )
            {
                denseValues[denseindex( i, j, numRows, numColumns )] *= rowValue;
            }
        }
        else
        {
            for ( IndexType j = 0; j < numColumns; ++j )
            {
                ValueType& updateVal = denseValues[denseindex( i, j, numRows, numColumns )];
                updateVal = applyBinary( updateVal, op, rowValue ); 
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDenseUtils::jacobi(
    ValueType solution[],
    const IndexType n,
    const ValueType denseValues[],
    const ValueType oldSolution[],
    const ValueType rhs[],
    const ValueType omega )
{
    SCAI_LOG_INFO( logger,
                   "jacobi<" << TypeTraits<ValueType>::id() << ">" << ", " << n << " x " << n
                    << ", omega = " << omega )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken != NULL )
    {
        SCAI_LOG_ERROR( logger, "jacobi called asynchronously, not supported here" )
    }

    ValueType omega1 = 1;
    omega1 -= omega;

    // computations for each element are independent

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.Dense.Jacobi" )
        #pragma omp for 

        for ( IndexType i = 0; i < n; i++ )
        {
            ValueType temp = rhs[i];
            ValueType diag = 0;

            for ( IndexType j = 0; j < n; ++j )
            {
                if ( i == j )
                {
                    // save the diagonal element for division later

                    diag = denseValues[ i * n + j ];
                }
                else
                {
                    temp -= denseValues[i * n + j] * oldSolution[j];
                }
            }

            solution[i] = omega * ( temp / diag ) + omega1 * oldSolution[i];
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDenseUtils::jacobiHalo(
    ValueType solution[],
    const ValueType diagonal[],
    const IndexType numRows,
    const IndexType numColumns,
    const ValueType denseValues[],
    const ValueType oldSolution[],
    const ValueType omega )
{
    SCAI_LOG_INFO( logger,
                   "jacobiHalo<" << TypeTraits<ValueType>::id() << ">" << ", " << numRows << " x " << numColumns
                    << ", omega = " << omega )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken != NULL )
    {
        SCAI_LOG_ERROR( logger, "jacobi called asynchronously, not supported here" )
    }

    for ( IndexType i = 0; i < numRows; i++ )
    {
        ValueType sum = 0;

        for ( IndexType j = 0; j < numColumns; ++j )
        {
            sum += denseValues[i * numColumns + j] * oldSolution[j];
        }

        solution[i] -= omega * sum / diagonal[i];
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDenseUtils::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register DenseUtils OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<DenseKernelTrait::nonZeroValues<ValueType> >( nonZeroValues, ctx, flag );
    KernelRegistry::set<DenseKernelTrait::getCSRSizes<ValueType> >( getCSRSizes, ctx, flag );
    KernelRegistry::set<DenseKernelTrait::setValue<ValueType> >( setValue, ctx, flag );
    KernelRegistry::set<DenseKernelTrait::setRows<ValueType> >( setRows, ctx, flag );
    KernelRegistry::set<DenseKernelTrait::setCSRValues<ValueType> >( setCSRValues, ctx, flag );
    KernelRegistry::set<DenseKernelTrait::getCSRValues<ValueType> >( getCSRValues, ctx, flag );
    KernelRegistry::set<DenseKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
    KernelRegistry::set<DenseKernelTrait::jacobiHalo<ValueType> >( jacobiHalo, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPDenseUtils::OpenMPDenseUtils()
{
    SCAI_LOG_INFO( logger, "register DenseUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

OpenMPDenseUtils::~OpenMPDenseUtils()
{
    SCAI_LOG_INFO( logger, "unregister DenseUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPDenseUtils OpenMPDenseUtils::guard;

/* --------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
