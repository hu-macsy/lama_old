/**
 * @file OpenMPDIAUtils.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
#include <scai/common/bind.hpp>

// std
#include <cmath>

namespace scai
{

namespace sparsekernel
{

using common::TypeTraits;
using tasking::TaskSyncToken;

SCAI_LOG_DEF_LOGGER( OpenMPDIAUtils::logger, "OpenMP.DIAUtils" )

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
    ValueType maxValue = static_cast<ValueType>(0.0);

    #pragma omp parallel
    {
        ValueType threadVal = static_cast<ValueType>(0.0);

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

        for( IndexType i = 0; i < numRows; ++i )
        {
            for( IndexType d = 0; d < numDiagonals; ++d )
            {
                const IndexType j = i + diaOffsets[d];

                if( ( j < 0 ) || ( j >= numColumns ) )
                {
                    continue;
                }

                const ValueType val = common::Math::abs( diaValues[i + d * numRows] );

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

            #pragma omp parallel for schedule( SCAI_OMP_SCHEDULE )

            for( IndexType i = 0; i < n; i++ )
            {
                csrJA[i] = i;
                csrValues[i] = static_cast<CSRValueType>(0.0);
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
        SCAI_REGION( "OpenMP.DIA->CSR_values" )

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

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

                bool nonZero = common::Math::abs( value ) > eps;

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
    SCAI_LOG_INFO( logger,
                   "get CSRSizes<" << TypeTraits<DIAValueType>::id() << "> for DIA matrix " << numRows << " x " << numColumns
                   << ", #diagonals = " << numDiagonals << ", eps = " << eps << ", diagonalFlag = " << diagonalFlag )

    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

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

            bool nonZero = common::Math::abs( diaValues[i + ii * numRows] ) > eps;

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
void OpenMPDIAUtils::normalGEMV_a(
    ValueType result[],
    const std::pair<ValueType, const ValueType*> ax,
    const std::pair<ValueType, const ValueType*> by,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[] )
{
    normalGEMV( result, ax.first, ax.second, by.first, by.second,
                numRows, numColumns, numDiagonals, diaOffsets, diaValues );
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
    // Note: launching thread gets token here, asynchronous thread will get NULL

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
                                      numRows, numColumns, numDiagonals, diaOffsets, diaValues ) );
        return;
    }

    SCAI_LOG_INFO( logger,
                   "normalGEMV<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads()
                   << ">, result[" << numRows << "] = " << alpha << " * A( dia, #diags = " << numDiagonals << " ) * x + " << beta << " * y " )

    // result := alpha * A * x + beta * y -> result:= beta * y; result += alpha * A

    utilskernel::OpenMPUtils::setScale( result, beta, y, numRows );

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.DIA.normalGEMV" )

        #pragma omp for schedule ( SCAI_OMP_SCHEDULE )

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

        if( SCAI_LOG_TRACE_ON( logger ) )
        {
            std::cout << "NormalGEMV: result = ";

            for( IndexType i = 0; i < numRows; ++i )
            {
                std::cout << " " << result[i];
            }

            std::cout << std::endl;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDIAUtils::normalGEVM_a(
    ValueType result[],
    const std::pair<ValueType, const ValueType*> ax,
    const std::pair<ValueType, const ValueType*> by,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const IndexType diaOffsets[],
    const ValueType diaValues[] )
{
    normalGEVM( result, ax.first, ax.second, by.first, by.second,
                numRows, numColumns, numDiagonals, diaOffsets, diaValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDIAUtils::normalGEVM(
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
                   "normalGEVM<" << TypeTraits<ValueType>::id() << ", #threads = " << omp_get_max_threads() << ">, result[" << numRows << "] = " << alpha << " * A( dia, #diags = " << numDiagonals << " ) * x + " << beta << " * y " )

    SCAI_LOG_INFO( logger,
                   "normalGEVM<" << TypeTraits<ValueType>::id() << ">, n = " << numRows << ", d = " << numDiagonals )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        // bind has limited number of arguments, so take help routine for call

        SCAI_LOG_INFO( logger,
                       "normalGEMV<" << TypeTraits<ValueType>::id() << "> launch it asynchronously" )

        syncToken->run( common::bind( normalGEVM_a<ValueType>,
                                      result,
                                      std::pair<ValueType, const ValueType*>( alpha, x ),
                                      std::pair<ValueType, const ValueType*>( beta, y ),
                                      numRows, numColumns, numDiagonals, diaOffsets, diaValues ) );
        return;
    }

    // result := alpha * x * A + beta * y -> result:= beta * y; result += alpha * x * A

    utilskernel::OpenMPUtils::setScale( result, beta, y, numColumns );

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.DIA.normalGEVM" )

        #pragma omp for schedule ( SCAI_OMP_SCHEDULE )

        for( IndexType k = 0; k < numColumns; ++k )
        {
            ValueType accu = static_cast<ValueType>(0.0);

            for( IndexType i = 0; i < numRows; i++ )
            {
                for( IndexType ii = 0; ii < numDiagonals; ++ii )
                {
                    const IndexType j = i + diaOffsets[ii];

                    if( j >= numColumns )
                    {
                        break;
                    }

                    if( j == k )
                    {
                        accu += diaValues[ii * numRows + i] * x[i];
                    }
                }
            }

            result[k] += alpha * accu;
        }
    }

    if( SCAI_LOG_TRACE_ON( logger ) )
    {
        std::cout << "NormalGEVM: result = ";

        for( IndexType i = 0; i < numRows; ++i )
        {
            std::cout << " " << result[i];
        }

        std::cout << std::endl;
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
                   "jacobi<" << TypeTraits<ValueType>::id() << ">" << ", #rows = " << numRows << ", #cols = " << numColumns << ", #diagonals = " << numDiagonals << ", omega = " << omega )

    SCAI_ASSERT_EQUAL_DEBUG( 0, diaOffset[0] )
    // main diagonal must be first

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if( syncToken != NULL )
    {
        SCAI_LOG_ERROR( logger, "jacobi called asynchronously, not supported here" )
    }

    #pragma omp parallel
    {
        SCAI_REGION( "OpenMP.DIA.Jacobi" )

        #pragma omp for schedule( SCAI_OMP_SCHEDULE )

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

            solution[i] = omega * ( temp / diag ) + (static_cast<ValueType>(1.0) - omega ) * oldSolution[i];
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDIAUtils::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    common::context::ContextType ctx = common::context::Host;

    SCAI_LOG_INFO( logger, "register DIAUtils OpenMP-routines for Host at kernel registry [" << flag
                   << " --> " << common::getScalarType<ValueType>() << "]" )

    KernelRegistry::set<DIAKernelTrait::getCSRSizes<ValueType> >( getCSRSizes, ctx, flag );
    KernelRegistry::set<DIAKernelTrait::absMaxVal<ValueType> >( absMaxVal, ctx, flag );
    KernelRegistry::set<DIAKernelTrait::normalGEMV<ValueType> >( normalGEMV, ctx, flag );
    KernelRegistry::set<DIAKernelTrait::normalGEVM<ValueType> >( normalGEVM, ctx, flag );
    KernelRegistry::set<DIAKernelTrait::jacobi<ValueType> >( jacobi, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void OpenMPDIAUtils::RegistratorVO<ValueType, OtherValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;

    common::context::ContextType ctx = common::context::Host;

    SCAI_LOG_INFO( logger, "register DIAUtils OpenMP-routines for Host at kernel registry [" << flag
                   << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )

    KernelRegistry::set<DIAKernelTrait::getCSRValues<ValueType, OtherValueType> >( getCSRValues, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPDIAUtils::OpenMPDIAUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_HOST_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARITHMETIC_HOST_LIST, SCAI_ARITHMETIC_HOST_LIST>::call( flag );
}

OpenMPDIAUtils::~OpenMPDIAUtils()
{
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_HOST_LIST>::call( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_ARITHMETIC_HOST_LIST, SCAI_ARITHMETIC_HOST_LIST>::call( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPDIAUtils OpenMPDIAUtils::guard;

} /* end namespace sparsekernel */

} /* end namespace scai */
