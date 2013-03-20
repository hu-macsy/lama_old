/**
 * @file OpenMPDenseUtils.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Class defintion for OpenMP routines to be used for DenseUtilsInterface.
 * @author Thomas Brandes
 * @date 03.07.2012
 * $Id$
 */
#ifndef LAMA_OPENMP_DENSE_UTILS_HPP_
#define LAMA_OPENMP_DENSE_UTILS_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMATypes.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

// logging
#include <logging/logging.hpp>

#include <cmath>

namespace lama
{

/** This class provides OpenMP implementations for methods in lama::DenseUtilsInterface
 */

class LAMA_DLL_IMPORTEXPORT OpenMPDenseUtils
{
public:

    /** OpenMP implementation for DenseUtilsInterface::Counting::getCSRSizes */

    template<typename DenseValueType>
    static void getCSRSizes(
        IndexType csrSizes[],
        bool diagonalFlag,
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType denseValues[],
        const DenseValueType eps );

    /** OpenMP implementation for DenseUtilsInterface::Conversions::getCSRValues */

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

    /** OpenMP implementation for DenseUtilsInterface::Copy::copyDenseValues */

    template<typename DenseValueType1,typename DenseValueType2>
    static void copyDenseValues(
        DenseValueType1 newValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType2 oldValues[] );

    /** OpenMP implementation for DenseUtilsInterface::Conversions::setCSRValues */

    template<typename DenseValueType,typename CSRValueType>
    static void setCSRValues(
        DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const CSRValueType csrValues[] );

    /** OpenMP implementation for DenseUtilsInterface::Copy::getDiagonal */

    template<typename DiagonalValueType,typename DenseValueType>
    static void getDiagonal(
        DiagonalValueType diagonalValues[],
        const IndexType numDiagonalValues,
        const DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns );

    /** OpenMP implementation for DenseUtilsInterface::Copy::setDiagonal */

    template<typename DenseValueType,typename DiagonalValueType>
    static void setDiagonal(
        DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const DiagonalValueType diagonalValues[],
        const IndexType numDiagonalValues );

    /** OpenMP implementation for DenseUtilsInterface::Modify::scaleValue */

    template<typename DenseValueType>
    static void scaleValue(
        DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType val );

    /** OpenMP implementation for DenseUtilsInterface::Modify::setDiagonalValue */

    template<typename DenseValueType>
    static void setDiagonalValue(
        DenseValueType denseValues[],
        const IndexType numRows,
        const IndexType numColumns,
        const DenseValueType val );

    /** Routine that registers all routines of this class at the LAMA interface. */

    static void setInterface( struct DenseUtilsInterface& DenseUtils );

private:

    static inline IndexType denseindex(
        const IndexType i,
        const IndexType j,
        const IndexType /* numRows */,
        const IndexType numColumns )
    {
        return i * numColumns + j;
    }

    LAMA_LOG_DECL_STATIC_LOGGER( logger );

};

/* --------------------------------------------------------------------------- */

} // namespace lama

#endif  //  LAMA_DENSE_STORAGE_UTILS_HPP_
