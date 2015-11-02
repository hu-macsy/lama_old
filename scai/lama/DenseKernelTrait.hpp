/**
 * @file DenseKernelTrait.hpp
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
 * @brief Struct with traits for all LAMA utilities provided as kernels.
 * @author Thomas Brandes
 * @date 03.04.2013
 * @since 1.0.0
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

namespace scai
{

namespace lama
{

/** Traits for utility functions to be used in Dense storage.   */

struct DenseKernelTrait
{
    template <typename DenseValueType>
    struct nonZeroValues
    {
        /** Counting non-zero values in dense storage.
         *
         *  @param[in]  denseValues size is numRows x numColumns, array with all matrix elements of dense format
         *  @param[in]  numRows number of rows
         *  @param[in]  numColumns number of columns
         *  @param[in]  eps is threshold when an element is to be considered as non-zero
         *  @return     number of values whole absolute value is larger than eps
         */

        typedef IndexType ( *FuncType )(
            const DenseValueType denseValues[],
            const IndexType numRows,
            const IndexType numColumns,
            const DenseValueType eps );

        static const char* getId() { return "Dense.nonZeroValues"; }
    };

    template <typename DenseValueType>
    struct getCSRSizes
    {
        /** Counting non-zero values in dense storage for conversion to CSR
         *
         *  @param[out] csrSizes is an array that contains for each row the number of non-zero elements
         *  @param[in]  diagonalFlag if true the diagonal elements are counted in any case
         *  @param[in]  numRows number of rows
         *  @param[in]  numColumns number of columns
         *  @param[in]  denseValues size is numRows x numColumns, array with all matrix elements of dense format
         *  @param[in]  eps is threshold when an element is to be considered as non-zero
         *
         *  The matrix values are stored row-wise in denseValues.
         */

        typedef void ( *FuncType )(
            IndexType csrSizes[],
            bool diagonalFlag,
            const IndexType numRows,
            const IndexType numColumns,
            const DenseValueType denseValues[],
            const DenseValueType eps );

        static const char* getId() { return "Dense.getCSRSizes"; }
    };

    template <typename DenseValueType, typename CSRValueType>
    struct getCSRValues
    {
        /** Convesion of dense matrix to CSR storage format
         *
         *  @param[out] csrJA will contain the column indexes
         *  @param[out] csrValues will contain the matrix elements
         *  @param[in]  csrIA is the array with the offsets (must already be available before)
         *  @param[in]  diagonalFlag if true the diagonal elements are also filled
         *  @param[in]  numRows number of rows
         *  @param[in]  numColumns number of columns
         *  @param[in]  denseValues size is numRows x numColumns, array with all matrix elements of dense format
         *  @param[in]  eps is threshold when an element is to be considered as non-zero
         *
         *  Very important: the offsets in csrIA must correspond to the csrSizes computed
         *                  by getCSRSizes.
         */
        typedef void ( *FuncType ) ( IndexType csrJA[],
                        CSRValueType csrValues[],
                        const IndexType csrIA[],
                        const bool diagonalFlag,
                        const IndexType numRows,
                        const IndexType numColumns,
                        const DenseValueType denseValues[],
                        const DenseValueType eps );

        static const char* getId() { return "Dense.getCSRValues"; }
    };

    template <typename DenseValueType, typename CSRValueType>
    struct setCSRValues
    {
        /** Conversion of CSR format to dense matrix. */

        typedef void ( *FuncType ) ( DenseValueType denseValues[],
                        const IndexType numRows,
                        const IndexType numColumns,
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const CSRValueType csrValues[] );

        static const char* getId() { return "Dense.setCSRValues"; }
    };

    template<typename DenseValueType1, typename DenseValueType2>
    struct copyDenseValues
    {
        /** Copy values of dense matrix; supports also conversion. */

        typedef void ( *FuncType ) ( DenseValueType1 newValues[],
                        const IndexType numRows,
                        const IndexType numColumns,
                        const DenseValueType2 oldValues[] );

        static const char* getId() { return "Dense.copyDenseValues"; }
    };

    template<typename RowValueType, typename DenseValueType>
    struct getRow
    {
        /** Get diagonal of a dense matrix, type conversion is supported. */

        typedef void ( *FuncType ) (
            RowValueType rowValues[],
            const DenseValueType denseValues[],
            const IndexType irow,
            const IndexType numRows,
            const IndexType numColumns );

        static const char* getId() { return "Dense.getRow"; }
    };

    template<typename DenseValueType1, typename DenseValueType2>
    struct getDiagonal
    {
        /** Get diagonal of a dense matrix, type conversion is supported. */

        typedef void ( *FuncType ) ( DenseValueType1 diagonalValues[],
                        const IndexType numDiagonalValues,
                        const DenseValueType2 denseValues[],
                        const IndexType numRows,
                        const IndexType numColumns );

        static const char* getId() { return "Dense.getDiagonal"; }
    };

    template<typename DenseValueType1, typename DenseValueType2>
    struct setDiagonal
    {
        /** Set diagonal of a dense matrix, type conversion is supported. */

        typedef void ( *FuncType ) ( DenseValueType1 denseValues[],
                        const IndexType numRows,
                        const IndexType numColumns,
                        const DenseValueType2 diagonalValues[],
                        const IndexType numDiagonalValues );

        static const char* getId() { return "Dense.setDiagonal"; }
    };

    template<typename DenseValueType>
    struct setValue
    {
        /** Set all elements of the dense matrix with a value */

        typedef void ( *FuncType ) ( 
            DenseValueType denseValues[],
            const IndexType numRows,
            const IndexType numColumns,
            const DenseValueType val );

        static const char* getId() { return "Dense.setValue"; }
    };

    template<typename DenseValueType>
    struct scaleValue
    {
        /** Scale all elements of the dense matrix with a value */

        typedef void ( *FuncType ) (
            DenseValueType denseValues[],
            const IndexType numRows,
            const IndexType numColumns,
            const DenseValueType val );

        static const char* getId() { return "Dense.scaleValue"; }
    };

    template<typename DenseValueType>
    struct setDiagonalValue
    {
        /** Set diagonal elements with one and the same value. */

        typedef void ( *FuncType ) ( 
            DenseValueType denseValues[],
            const IndexType numRows,
            const IndexType numColumns,
            const DenseValueType val );

        static const char* getId() { return "Dense.setDiagonalValue"; }
    };

    template<typename DenseValueType, typename OtherValueType>
    struct scaleRows
    {
        /** Scale rows of the matrix individually.
         *
         *  @param[in,out] denseValues  data of the dense matrix, size is numRows * numColumns
         *  @param[in]     numRows      number of rows
         *  @param[in]     numColumns   number of columns
         *  @param[in]     rowValues    scale values for each row, size is numRows
         */

        typedef void ( *FuncType ) ( 
            DenseValueType denseValues[],
            const IndexType numRows,
            const IndexType numColumns,
            const OtherValueType rowValues[] );

        static const char* getId() { return "Dense.scaleRows"; }
    };
};

} /* end namespace lama */

} /* end namespace scai */
