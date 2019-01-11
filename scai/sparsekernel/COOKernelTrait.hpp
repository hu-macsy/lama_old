/**
 * @file COOKernelTrait.hpp
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
 * @brief Struct with traits for all COO storage methods provided as kernels.
 * @author Thomas Brandes
 * @date 22.10.2015
 */
#pragma once

#include <scai/common/config.hpp>

#include <scai/common/MatrixOp.hpp>
#include <scai/common/BinaryOp.hpp>

namespace scai
{

namespace sparsekernel
{

/** Kernel traits for functions to be used in COO storage.
 *
 *  Note: routines to build CSR data from COO data are not required any more
 *        as this is done by bucket sort
 */

struct COOKernelTrait
{
    struct getValuePos
    {
        /** Returns position of element (i,j) in ja/values array
         *
         *  @param[in] i is the row of the element
         *  @param[in] j is the column of the element
         *  @param[in] cooIA is the COO ia array
         *  @param[in] cooJA is the COO ja array
         *  @returns  offset of element in values array, invalidIndex if not found
         */

        typedef IndexType ( *FuncType ) (
            const IndexType i,
            const IndexType j,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const IndexType numValues );

        static const char* getId()
        {
            return "COO.getValuePos";
        }
    };

    struct getColumn
    {
        /** Get the non-zero entries for one column
         *
         */
        typedef IndexType ( *FuncType ) (
            IndexType positions[],
            const IndexType cooJA[],
            const IndexType numValues,
            const IndexType j );

        static const char* getId()
        {
            return "COO.getColumn";
        }
    };

    struct getRow
    {
        /** Get the non-zero entries for one row
         *
         */
        typedef IndexType ( *FuncType ) (
            IndexType& offset,
            const IndexType cooIA[],
            const IndexType numValues,
            const IndexType i );

        static const char* getId()
        {
            return "COO.getRow";
        }
    };

    struct hasDiagonalProperty
    {
        /** Routine checks for diagonal property, i.e. entries for all diagonal elements are available
         *
         *  @param[in] cooIA row indexes
         *  @param[in] cooJA column indexes
         *  @param[in] numDiagonals number of diagonal elements to verify
         *  @param[in] numValues size of array cooIA and cooJA
         *  @return true if there is one entry for each diagonal element
         */

        typedef bool ( *FuncType )(
            const IndexType numDiagonals,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const IndexType numValues );

        static const char* getId()
        {
            return "COO.hasDiagonalProperty";
        }
    };

    struct offsets2ia
    {
        /** Routine for conversion of CSR offset array to COO ia array
         *
         *  @param[out] cooIA is the array with all row indexes
         *  @param[in] numValues number of non-zero values, size of cooIA
         *  @param[in] csrIA is the CSR row offset array, size is numRows+1
         *  @param[in] numRows number of rows
         */

        typedef void ( *FuncType )(
            IndexType cooIA[],
            const IndexType numValues,
            const IndexType csrIA[],
            const IndexType numRows );

        static const char* getId()
        {
            return "COO.offsets2ia";
        }
    };

    struct ia2offsets
    {
        /** Routine for conversion of COO ia array to CSR offsets array
         *
         *  @param[out] csrIA is the CSR row offset array, size is numRows+1
         *  @param[in] numRows number of rows
         *  @param[in] cooIA is the 'sorted' array with all row indexes
         *  @param[in] numValues number of non-zero values, size of cooIA
         */
        typedef void ( *FuncType )(
            IndexType csrIA[],
            const IndexType numRows,
            const IndexType cooIA[],
            const IndexType numValues );

        static const char* getId()
        {
            return "COO.ia2offsets";
        }
    };

    template<typename ValueType>
    struct setRows
    {

        /** This operation multiplies each row with an own value.
         *
         *  @param[in,out] cooValues matrix data that is scaled
         *  @param[in]     cooIA are the row indexes
         *  @param[in]     rowValues array with scale factor for each row
         *  @param[in]     numValues number of entries in cooValues and cooIA
         *
         *  This routine supports different precision for matrix values and scale values.
         */

        typedef void ( *FuncType ) (
            ValueType cooValues[],
            const IndexType cooIA[],
            const ValueType rowValues[],
            const IndexType numValues,
            common::BinaryOp op );

        static const char* getId()
        {
            return "COO.setRows";
        }
    };

    template<typename ValueType>
    struct setColumns
    {

        /** This operation updates each entry of the storage with a column-specific values
         *
         *  @param[in,out] cooValues matrix data that is scaled
         *  @param[in]     cooJA are the column indexes
         *  @param[in]     columnValues array with update value, one for each column
         *  @param[in]     numValues number of entries in cooValues and cooJA
         *
         */
        typedef void ( *FuncType ) (
            ValueType cooValues[],
            const IndexType cooJA[],
            const ValueType columnValues[],
            const IndexType numValues,
            common::BinaryOp op );

        static const char* getId()
        {
            return "COO.setColumns";
        }
    };

    template<typename ValueType>
    struct normalGEMV
    {
        /** result = alpha * CSR-Matrix * x + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for vectors result and b and rows of matrix
         *  @param numColumns is number of elements for columns of matrix, size of result, y
         *  @param cooIA, cooJA, cooValues are arrays of COO storage
         *  @param numValues is the size of the coo arrays
         */

        typedef void ( *FuncType ) (
            ValueType result[],
            const ValueType alpha,
            const ValueType x[],
            const ValueType beta,
            const ValueType y[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType nnz,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const ValueType cooValues[],
            const common::MatrixOp op );

        static const char* getId()
        {
            return "COO.normalGEMV";
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
            ValueType* const solution,
            const IndexType cooNumValues,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const ValueType cooValues[],
            const ValueType oldSolution[],
            const ValueType rhs[],
            const ValueType omega,
            const IndexType numRows );

        static const char* getId()
        {
            return "COO.jacobi";
        }
    };

    template<typename ValueType>
    struct jacobiHalo
    {
        /** Method to compute one iteration step in Jacobi method on halo storage
         */
        typedef void ( *FuncType ) (
            ValueType* const solution,
            const IndexType cooNumValues,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const ValueType cooValues[],
            const ValueType localDiagonal[],
            const ValueType oldSolution[],
            const ValueType omega,
            const IndexType numRows );

        static const char* getId()
        {
            return "COO.jacobiHalo";
        }
    };

    template<typename ValueType>
    struct getDiagonal
    {
        /** This method returns the diagonal of a coo storage.
         *
         */
        typedef void ( *FuncType ) (
            ValueType diagonal[],
            const IndexType numDiagonals,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const ValueType cooValues[],
            const IndexType numValues );

        static const char* getId()
        {
            return "COO.getDiagonal";
        }
    };

    template<typename ValueType>
    struct setDiagonalV
    {
        /** This method returns the diagonal of a coo storage.
         *
         */
        typedef void ( *FuncType ) (
            ValueType cooValues[],
            const ValueType diagonal[],
            const IndexType numDiagonals,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const IndexType numValues );

        static const char* getId()
        {
            return "COO.setDiagonalV";
        }
    };

    template<typename ValueType>
    struct setDiagonal
    {
        /** This method returns the diagonal of a coo storage.
         *
         */
        typedef void ( *FuncType ) (
            ValueType cooValues[],
            const ValueType diagonal,
            const IndexType numDiagonals,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const IndexType numValues );

        static const char* getId()
        {
            return "COO.setDiagonal";
        }
    };
};

} /* end namespace sparsekernel */

} /* end namespace scai */
