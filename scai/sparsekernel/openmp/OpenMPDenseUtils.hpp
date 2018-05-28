/**
 * @file OpenMPDenseUtils.hpp
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
 * @brief Class defintion for OpenMP routines to be used for DenseKernelTrait.
 * @author Thomas Brandes
 * @date 03.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/kregistry/mepr/Registrator.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/BinaryOp.hpp>


// std
#include <cmath>

namespace scai
{

namespace sparsekernel
{

/** This class provides OpenMP implementations for methods in scai::lama::DenseKernelTrait */

class COMMON_DLL_IMPORTEXPORT OpenMPDenseUtils
{
public:

    /** OpenMP implementation for DenseKernelTrait::nonZeroValues */

    template<typename DenseValueType>
    static IndexType nonZeroValues(
        const DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType eps );

    /** OpenMP implementation for DenseKernelTrait::getCSRSizes */

    template<typename ValueType>
    static void getCSRSizes(
        IndexType csrSizes[],
        const IndexType numRows,
        const IndexType numColumns,
        const ValueType denseValues[],
        const RealType<ValueType> eps );

    /** OpenMP implementation for DenseKernelTrait::getCSRValues */

    template<typename CSRValueType, typename DenseValueType>
    static void getCSRValues(
        IndexType csrJA[],
        CSRValueType csrValues[],
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType denseValues[],
        const RealType<DenseValueType> eps );

    /** OpenMP implementation for DenseKernelTrait::set */

    template<typename DenseValueType1, typename DenseValueType2>
    static void set(
        DenseValueType1 newValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType2 oldValues[],
        const common::BinaryOp op );

    /** OpenMP implementation for DenseKernelTrait::setCSRValues */

    template<typename DenseValueType, typename CSRValueType>
    static void setCSRValues(
        DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const CSRValueType csrValues[] );

    /** OpenMP implementation for DenseKernelTrait::setValue */

    template<typename DenseValueType>
    static void setValue(
        DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType val,
        const common::BinaryOp op );

    /** OpenMP implementation for DenseKernelTrait::scaleRows */

    template<typename ValueType>
    static void scaleRows(
        ValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const ValueType rowValues[] );

private:

    static inline IndexType denseindex(
        const IndexType i,
        const IndexType j,
        const IndexType /* numRows */,
        const IndexType numColumns )
    {
        return i * numColumns + j;
    }

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Struct for registration of methods with one template argument.
     *
     *  Registration function is wrapped in struct/class that can be used as template
     *  argument for metaprogramming classes to expand for each supported type
     */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Struct for registration of methods with two template arguments.
     *
     *  Registration function is wrapped in struct/class that can be used as template
     *  argument for metaprogramming classes to expand for all supported types.
     */

    template<typename ValueType, typename OtherValueType>
    struct RegistratorVO
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    OpenMPDenseUtils();

    /** Destructor for unregistration. */

    ~OpenMPDenseUtils();

    /** Static variable for registration at static initialization. */

    static OpenMPDenseUtils guard;

};

/* --------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
