/**
 * @file OpenMPDenseUtils.hpp
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


// std
#include <cmath>

namespace scai
{

namespace sparsekernel
{

/** This class provides OpenMP implementations for methods in scai::lama::DenseKernelTrait
 */

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

    template<typename DenseValueType>
    static void getCSRSizes(
        IndexType csrSizes[],
        bool diagonalFlag,
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType denseValues[],
        const DenseValueType eps );

    /** OpenMP implementation for DenseKernelTrait::getCSRValues */

    template<typename DenseValueType,typename CSRValueType>
    static void getCSRValues(
        IndexType csrJA[],
        CSRValueType csrValues[],
        const IndexType csrIA[],
        const bool diagonalFlag,
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType denseValues[],
        const DenseValueType eps );

    /** OpenMP implementation for DenseKernelTrait::copyDenseValues */

    template<typename DenseValueType1,typename DenseValueType2>
    static void copyDenseValues(
        DenseValueType1 newValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType2 oldValues[] );

    /** OpenMP implementation for DenseKernelTrait::setCSRValues */

    template<typename DenseValueType,typename CSRValueType>
    static void setCSRValues(
        DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const CSRValueType csrValues[] );

    /** OpenMP implementation for DenseKernelTrait::getRow */

    template<typename RowValueType, typename DenseValueType>
    static void getRow(
        RowValueType rowValues[],
        const DenseValueType denseValues[],
        const IndexType irow, 
        const IndexType numRows,
        const IndexType numColumns );

    /** OpenMP implementation for DenseKernelTrait::getDiagonal */

    template<typename DiagonalValueType,typename DenseValueType>
    static void getDiagonal(
        DiagonalValueType diagonalValues[],
        const IndexType numDiagonalValues,
        const DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns );

    /** OpenMP implementation for DenseKernelTrait::setDiagonal */

    template<typename DenseValueType,typename DiagonalValueType>
    static void setDiagonal(
        DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const DiagonalValueType diagonalValues[],
        const IndexType numDiagonalValues );

    /** OpenMP implementation for DenseKernelTrait::setValue */

    template<typename DenseValueType>
    static void setValue(
        DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType val );

    /** OpenMP implementation for DenseKernelTrait::scaleValue */

    template<typename DenseValueType>
    static void scaleValue(
        DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType val );

    /** OpenMP implementation for DenseKernelTrait::setDiagonalValue::FuncType */

    template<typename DenseValueType>
    static void setDiagonalValue(
        DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType val );

    /** OpenMP implementation for DenseKernelTrait::scaleRows */

    template<typename DenseValueType, typename OtherType>
    static void scaleRows(
        DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const OtherType rowValues[] );

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

    /** Routine that registers all methods at the kernel registry. */

    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorV, template<typename ValueType> )
    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorVO, template<typename ValueType, typename OtherValueType> )

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
