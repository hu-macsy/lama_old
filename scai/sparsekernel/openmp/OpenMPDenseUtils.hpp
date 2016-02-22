/**
 * @file OpenMPDenseUtils.hpp
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
 * @brief Class defintion for OpenMP routines to be used for DenseKernelTrait.
 * @author Thomas Brandes
 * @date 03.07.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/assert.hpp>

#include <scai/logging.hpp>

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

    static void registerKernels( bool deleteFlag );

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
