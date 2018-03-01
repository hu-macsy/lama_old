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

#include <scai/tracing.hpp>

namespace scai
{

namespace sparsekernel
{

using common::TypeTraits;
using common::Math;

SCAI_LOG_DEF_LOGGER( OpenMPDenseUtils::logger, "OpenMP.DenseUtils" )

/* --------------------------------------------------------------------------- */
/*     Template implementations                                                */
/* --------------------------------------------------------------------------- */

template<typename DenseValueType>
IndexType OpenMPDenseUtils::nonZeroValues(
    const DenseValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType eps )
{
    typedef typename TypeTraits<DenseValueType>::RealType RealType;
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

template<typename DenseValueType>
void OpenMPDenseUtils::getCSRSizes(
    IndexType csrSizes[],
    bool diagonalFlag,
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType denseValues[],
    const DenseValueType eps )
{
    typedef typename TypeTraits<DenseValueType>::RealType RealType;

    RealType absEps = Math::abs( eps );

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
            const DenseValueType& value = denseValues[denseindex( i, j, numRows, numColumns )];

            if ( common::Math::abs( value ) > absEps )
            {
                ++nonZeros;
            }
            else if ( i == j && diagonalFlag )
            {
                ++nonZeros; // count also zero elements for diagonals
            }
        }

        csrSizes[i] = nonZeros;
    }
}

/** Helper routine for conversion Dense to CSR */

template<typename CSRValueType, typename DenseValueType>
void OpenMPDenseUtils::getCSRValues(
    IndexType csrJA[],
    CSRValueType csrValues[],
    const IndexType csrIA[],
    bool diagonalFlag,
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType denseValues[],
    const DenseValueType eps )
{
    typedef typename TypeTraits<DenseValueType>::RealType RealType;

    RealType absEps = eps;

    SCAI_REGION( "OpenMP.DenseUtils.getCSRValues" )

    SCAI_LOG_INFO( logger,
                   "get CSRValues<" << TypeTraits<DenseValueType>::id() << ", " << TypeTraits<CSRValueType>::id() << ">" << ", size is " << numRows << " x " << numColumns )
    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType offset = csrIA[i];

        if ( i < numColumns && diagonalFlag )
        {
            // start with diagonal element in any case
            const DenseValueType& value = denseValues[denseindex( i, i, numRows, numColumns )];
            csrValues[offset] = static_cast<CSRValueType>( value );
            csrJA[offset] = i;
            offset++;
        }

        for ( IndexType j = 0; j < numColumns; ++j )
        {
            if ( i == j && diagonalFlag )
            {
                continue; // diagonal element is already first element in ja, values
            }

            const DenseValueType& value = denseValues[denseindex( i, j, numRows, numColumns )];

            RealType absValue = common::Math::abs( value );

            if ( absValue > absEps )
            {
                csrValues[offset] = static_cast<CSRValueType>( value );
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

template<typename DenseValueType, typename CSRValueType>
void OpenMPDenseUtils::setCSRValues(
    DenseValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType csrIA[],
    const IndexType csrJA[],
    const CSRValueType csrValues[] )
{
    SCAI_REGION( "OpenMP.DenseUtils.setCSRValues" )

    SCAI_LOG_INFO( logger,
                   "set CSRValues<" << TypeTraits<DenseValueType>::id() << ", " << TypeTraits<CSRValueType>::id() << ">" << ", size is " << numRows << " x " << numColumns )
    // parallelization possible as offset array csrIA is available
    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; i++ )
    {
        // Initialize complete row with zero values
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            DenseValueType& elem = denseValues[denseindex( i, j, numRows, numColumns )];
            elem = static_cast<DenseValueType>( 0.0 );
        }

        // fill up positions for which non-zero values are given

        for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
        {
            const IndexType j = csrJA[jj];
            DenseValueType& elem = denseValues[denseindex( i, j, numRows, numColumns )];
            elem = static_cast<DenseValueType>( csrValues[jj] );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename DenseValueType1, typename DenseValueType2>
void OpenMPDenseUtils::set(
    DenseValueType1 out[],
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType2 in[],
    const common::BinaryOp op )
{
    using namespace common;

    SCAI_REGION( "OpenMP.DenseUtils.set" )

    switch ( op )
    {
        case BinaryOp::COPY :
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < numRows; ++i )
            {
                for ( IndexType j = 0; j < numColumns; ++j )
                {
                    IndexType k = denseindex( i, j, numRows, numColumns );
                    out[k] = static_cast<DenseValueType1>( in[k] );
                }
            }

            break;
        }

        default:
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < numRows; ++i )
            {
                for ( IndexType j = 0; j < numColumns; ++j )
                {
                    IndexType k = denseindex( i, j, numRows, numColumns );
                    out[k] = applyBinary( out[k], op, static_cast<DenseValueType1>( in[k] ) );
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename DenseValueType>
void OpenMPDenseUtils::setDiagonalValue(
    DenseValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType diagonalValue )
{
    SCAI_REGION( "OpenMP.DenseUtils.setDiagonalValue" )

    IndexType numDiagonalValues = std::min( numRows, numColumns );
    #pragma omp parallel for 

    for ( IndexType i = 0; i < numDiagonalValues; ++i )
    {
        DenseValueType& elem = denseValues[denseindex( i, i, numRows, numColumns )];
        elem = diagonalValue;
    }
}

/* --------------------------------------------------------------------------- */

template<typename DenseValueType>
void OpenMPDenseUtils::setValue(
    DenseValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const DenseValueType val,
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
                    DenseValueType& elem = denseValues[denseindex( i, j, numRows, numColumns )];
                    elem = applyBinary( elem, op, val );
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPDenseUtils::scaleRows(
    ValueType denseValues[],
    const IndexType numRows,
    const IndexType numColumns,
    const ValueType rowValues[] )
{
    SCAI_REGION( "OpenMP.DenseUtils.scaleRows" )

    #pragma omp parallel for

    for ( IndexType i = 0; i < numRows; ++i )
    {
        const ValueType scaleValue = rowValues[i];

        // scale the whole row with this value

        for ( IndexType j = 0; j < numColumns; ++j )
        {
            denseValues[denseindex( i, j, numRows, numColumns )] *= scaleValue;
        }
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
    KernelRegistry::set<DenseKernelTrait::setDiagonalValue<ValueType> >( setDiagonalValue, ctx, flag );
    KernelRegistry::set<DenseKernelTrait::scaleRows<ValueType> >( scaleRows, ctx, flag );
}

template<typename ValueType, typename OtherValueType>
void OpenMPDenseUtils::RegistratorVO<ValueType, OtherValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "register DenseUtils OpenMP-routines for Host at kernel registry [" << flag
                    << " --> " << common::getScalarType<ValueType>() << ", " << common::getScalarType<OtherValueType>() << "]" )
    KernelRegistry::set<DenseKernelTrait::setCSRValues<ValueType, OtherValueType> >( setCSRValues, ctx, flag );
    KernelRegistry::set<DenseKernelTrait::getCSRValues<ValueType, OtherValueType> >( getCSRValues, ctx, flag );
    KernelRegistry::set<DenseKernelTrait::set<ValueType, OtherValueType> >( set, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPDenseUtils::OpenMPDenseUtils()
{
    SCAI_LOG_INFO( logger, "register DenseUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_HOST_LIST, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

OpenMPDenseUtils::~OpenMPDenseUtils()
{
    SCAI_LOG_INFO( logger, "unregister DenseUtils OpenMP-routines for Host at kernel registry" )

    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
    kregistry::mepr::RegistratorVO<RegistratorVO, SCAI_NUMERIC_TYPES_HOST_LIST, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPDenseUtils OpenMPDenseUtils::guard;

/* --------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
