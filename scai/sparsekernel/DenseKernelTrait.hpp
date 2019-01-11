/**
 * @file DenseKernelTrait.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Struct with traits for all LAMA utilities provided as kernels.
 * @author Thomas Brandes
 * @date 03.04.2013
 */
#pragma once

// for dll_import

#include <scai/common/config.hpp>
#include <scai/common/BinaryOp.hpp>

namespace scai
{

namespace sparsekernel
{

/** Traits for utility functions to be used in Dense storage.   */

struct DenseKernelTrait
{
    template<typename ValueType>
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
            const ValueType denseValues[],
            const IndexType numRows,
            const IndexType numColumns,
            const ValueType eps );

        static const char* getId()
        {
            return "Dense.nonZeroValues";
        }
    };

    template<typename ValueType>
    struct getCSRSizes
    {
        /** Counting non-zero values in dense storage for conversion to CSR
         *
         *  @param[out] csrSizes is an array that contains for each row the number of non-zero elements
         *  @param[in]  numRows number of rows
         *  @param[in]  numColumns number of columns
         *  @param[in]  denseValues size is numRows x numColumns, array with all matrix elements of dense format
         *  @param[in]  eps is threshold when an element is to be considered as non-zero
         *
         *  The matrix values are stored row-wise in denseValues.
         */

        typedef void ( *FuncType )(
            IndexType csrSizes[],
            const IndexType numRows,
            const IndexType numColumns,
            const ValueType denseValues[],
            const RealType<ValueType> eps );

        static const char* getId()
        {
            return "Dense.getCSRSizes";
        }
    };

    template <typename ValueType>
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
        typedef void ( *FuncType ) ( 
            IndexType csrJA[],
            ValueType csrValues[],
            const IndexType csrIA[],
            const IndexType numRows,
            const IndexType numColumns,
            const ValueType denseValues[],
            const RealType<ValueType> eps );

        static const char* getId()
        {
            return "Dense.getCSRValues";
        }
    };

    template <typename ValueType>
    struct setCSRValues
    {
        /** Conversion of CSR format to dense matrix. */

        typedef void ( *FuncType ) ( ValueType denseValues[],
                                     const IndexType numRows,
                                     const IndexType numColumns,
                                     const IndexType csrIA[],
                                     const IndexType csrJA[],
                                     const ValueType csrValues[] );

        static const char* getId()
        {
            return "Dense.setCSRValues";
        }
    };

    template<typename ValueType>
    struct setValue
    {
        /** Set all elements of the dense matrix with a value */

        typedef void ( *FuncType ) (
            ValueType denseValues[],
            const IndexType numRows,
            const IndexType numColumns,
            const ValueType val,
            const common::BinaryOp op );

        static const char* getId()
        {
            return "Dense.setValue";
        }
    };

    template<typename ValueType>
    struct setRows
    {
        /** Set/update rows of the matrix individually.
         *
         *  @param[in,out] denseValues  data of the dense matrix, size is numRows * numColumns
         *  @param[in]     numRows      number of rows
         *  @param[in]     numColumns   number of columns
         *  @param[in]     rowValues    scale values for each row, size is numRows
         *  @param[in]     op           binary operation that is applied 
         */

        typedef void ( *FuncType ) (
            ValueType denseValues[],
            const IndexType numRows,
            const IndexType numColumns,
            const ValueType rowValues[],
            const common::BinaryOp op );

        static const char* getId()
        {
            return "Dense.setRows";
        }
    };

    template<typename ValueType>
    struct setColumns
    {
        /** Set/update rows of the matrix individually.
         *
         *  @param[in,out] denseValues  data of the dense matrix, size is numRows * numColumns
         *  @param[in]     numRows      number of rows
         *  @param[in]     numColumns   number of columns
         *  @param[in]     columnValues scale values for each column, size is numColumns
         *  @param[in]     op           binary operation that is applied 
         */

        typedef void ( *FuncType ) (
            ValueType denseValues[],
            const IndexType numRows,
            const IndexType numColumns,
            const ValueType columnValues[],
            const common::BinaryOp op );

        static const char* getId()
        {
            return "Dense.setColumns";
        }
    };

    /** Structure with type definitions for solver routines */

    template<typename ValueType>
    struct jacobi
    {
        /** Method to compute one iteration step in Jacobi method
         *
         *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
         *
         */
        typedef void ( *FuncType ) (
            ValueType solution[],
            const IndexType n,
            const ValueType denseValues[],
            const ValueType oldSolution[],
            const ValueType rhs[],
            const ValueType omega );

        static const char* getId()
        {
            return "Dense.jacobi";
        }
    };

    template<typename ValueType>
    struct jacobiHalo
    {
        /** Compute one iteration step in Jacobi method for halo
         *
         *  \code
         *      solution -= omega * ( dia_halo * oldSolution ) ./ diagonal 
         *  \endcode
         *
         */
        typedef void ( *FuncType ) (
            ValueType solution[],
            const ValueType diagonal[],
            const IndexType numRows,
            const IndexType numColumns,
            const ValueType denseValues[],
            const ValueType oldSolution[],
            const ValueType omega );

        static const char* getId()
        {
            return "Dense.jacobiHalo";
        }
    };

};

} /* end namespace sparsekernel */

} /* end namespace scai */
