/**
 * @file ELLKernelTrait.hpp
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
 * @brief Struct with traits for all ELL storage methods provided as kernels.
 * @author Thomas Brandes
 * @date 23.10.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/MatrixOp.hpp>
#include <scai/common/TypeTraits.hpp>

namespace scai
{

namespace sparsekernel
{

/** Kernel traits for LAMA kernels used by ELL storage.   */

struct ELLKernelTrait
{
    template<typename ValueType>
    struct jacobi
    {
        /**
         *
         *  @param[out] solution is solution vector, size is numRows
         *  @param[in]  numRows is size of vectors and number of rows for matrix
         *  @param[in]  ellNumValuesPerRow is maximal number of non-zero entries
         *  @param[in]  ellSizes, ellJA, ellValues are arrays of ELL storage, numRows x numRows
         *  @param[in]  oldSolution is the old solution, size is numRows
         *  @param[in]  rhs is right hand side vector, size is numRows
         *  @param[in]  omega is scaling factor
         *
         *  The ELL storage stands for a square matrix and must have diagonal property.
         */
        typedef void ( *FuncType )(
            ValueType solution[],
            const IndexType numRows,
            const IndexType ellNumValuesPerRow,
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ValueType ellValues[],
            const ValueType oldSolution[],
            const ValueType rhs[],
            const ValueType omega );

        static const char* getId()
        {
            return "ELL.jacobi";
        }
    };

    template<typename ValueType>
    struct jacobiHalo
    {
        typedef void ( *FuncType )(
            ValueType solution[],
            const IndexType numRows,
            const ValueType diagonal[],
            const IndexType ellNumValuesPerRow,
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ValueType ellValues[],
            const IndexType rowIndexes[],
            const IndexType numNonEmptyRows,
            const ValueType oldSolution[],
            const ValueType omega );

        static const char* getId()
        {
            return "ELL.jacobiHalo";
        }
    };

    template<typename ValueType>
    struct fillELLValues
    {

        /**
         *  This method fills up the arrays ja and values of the ELL format with
         *  useful values to make matrix-vector multiplication efficient.
         *
         *  @param[in,out] ellJA is array with column indexes
         *  @param[in,out] ellValues is array with non-zero matrix values
         *  @param[in]  ellSizes array with number of entries for each row
         *  @param[in]  numRows number of rows
         *  @param[in]  numValuesPerRow number of values in each row
         */
        typedef void ( *FuncType ) (
            IndexType ellJA[],
            ValueType ellValues[],
            const IndexType ellSizes[],
            const IndexType numRows,
            const IndexType numValuesPerRow );

        static const char* getId()
        {
            return "ELL.fillELLValues";
        }
    };

    /** Conversion routines between ELL and CSR storage format. */

    template<typename ELLValueType, typename CSRValueType>
    struct getCSRValues
    {
        /** Conversion from ELL data to CSR data
         *
         *  @param[out] csrJA will contain the column indexes
         *  @param[out] csrValues will contain the matrix elements
         *  @param[in]  csrIA is the array with the offsets (must already be available before)
         *  @param[in]  numRows is number of rows
         *  @param[in]  ellSizes is the number of values in each row
         *  @param[in]  ellJA are the column indexes for ELL format
         *  @param[in]  ellValues are the stored matrix values for ELL format
         */
        typedef void ( *FuncType ) (
            IndexType csrJA[],
            CSRValueType csrValues[],
            const IndexType csrIA[],
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ELLValueType ellValues[] );

        static const char* getId()
        {
            return "ELL.getCSRValues";
        }
    };

    template<typename ELLValueType, typename CSRValueType>
    struct setCSRValues
    {
        /** Conversion from CSR data to ELL data      */

        typedef void ( *FuncType ) (
            IndexType ellJA[],
            ELLValueType ellValues[],
            const IndexType ellSizes[],
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const IndexType csrIA[],
            const IndexType csrJA[],
            const CSRValueType csrValues[] );

        static const char* getId()
        {
            return "ELL.setCSRValues";
        }
    };

    template<typename ValueType>
    struct compressIA
    {
        /** Determine the new row sizes when ELL storage data will be compressed
         *
         * @param[out] newSizes are the new sizes of the rows when storage is compressed 
         * @param[in]  ellSizes are the sizes of the rows in the original ELL data
         * @param[in]  ellJA column indexes of ELL data
         * @param[in]  ellValues are values of the ELL data
         * @param[in]  numRows number of rows
         * @param[in]  numValuesPerRow for maximal number of entries in one row
         * @param[in]  eps threshold, abs value of entry must be greater than eps to be non-zero
         * @param[in]  keepDiagonal if true do not remove diagonal elements
         */
        typedef void ( *FuncType ) (
            IndexType newSizes[],
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ValueType ellValues[],
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const RealType<ValueType> eps,
            const bool keepDiagonal );

        static const char* getId()
        {
            return "ELL.compressIA";
        }
    };

    template<typename ValueType>
    struct compressValues
    {
        /** Compute the new ja and values array for compressed ELL stroage
         *
         * @param[out] new JA array for compressed storage
         * @param[out] new values array for compressed storage
         * @param[in]  newnumValuesPerRow for maximal number of entries in one row for compressed data
         * @param[in]  ellIA is the ia array of uncompressed ELL storage
         * @param[in]  ellJA ja array of uncompressed ELL storage
         * @param[in]  ellValues is the values array of uncompressed ELL storage
         * @param[in]  numRows number of rows
         * @param[in]  numValuesPerRow for maximal number of entries in one row
         * @param[in]  eps threshold, abs value of entry must be greater than eps to be non-zero
         * @param[in]  keepDiagonal if true do not remove diagonal elements
         */
        typedef void ( *FuncType ) (
            IndexType newJA[],
            ValueType newValues[],
            const IndexType newNumValuesPerRow,
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ValueType ellValues[],
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const RealType<ValueType> eps,
            const bool keepDiagonal );

        static const char* getId()
        {
            return "ELL.compressValues";
        }
    };

    template<typename ValueType>
    struct getRow
    {
        /** Returns a row of the matrix as dense vector
         *
         *  @param[out] row as dense vector that will be returned
         *  @param[in]  i is the row that should be returned
         *  @param[in]  numRows is the number of rows of the ELL matrix
         *  @param[in]  numColums is size of ia
         *  @param[in]  ia is the ELL sizes array
         *  @param[in]  ja is the ELL ja array
         *  @param[in]  values is the ELL values array
         */

        typedef void ( *FuncType ) (
            ValueType row[],
            const IndexType i,
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numValuesPerRow,
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ValueType ellValues[] );

        static const char* getId()
        {
            return "ELL.compressValues";
        }
    };

    struct getValuePos
    {
        /** Returns one element of the matrix
         *
         *  @param[in] i is the row of the returned element
         *  @param[in] j is the column of the returned element
         *  @param[in] numRows is the number of rows of the matrix
         *  @param[in] numValuesPerRow is the maximal number of entries in one row
         *  @param[in] ellSizes is the ELL sizes array
         *  @param[in] ellJA is the ELL ja array
         */

        typedef IndexType ( *FuncType ) (
            const IndexType i,
            const IndexType j,
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const IndexType ellSizes[],
            const IndexType ellJA[] );

        static const char* getId()
        {
            return "ELL.getValuePos";
        }
    };

    struct getValuePosCol
    {
        /** This method returns for a certain column of the CSR matrix all
         *  row indexes for which elements exist and the corresponding positions
         *  in the ellJA/ellValues array
         *
         *  @param[out] row indexes of rows that have an entry for column j
         *  @param[out] pos positions of entries with col = j in csrJA,
         *  @param[in] j is the column of which positions are required
         *  @param[in] ellIA is the ELL sizes array
         *  @param[in] numRows is the number of rows
         *  @param[in] csrJA is the CSR ja array
         *  @param[in] numValuesPerRow is maximal size of one row
         *  @returns  number of entries with col index = j
         */
        typedef IndexType ( *FuncType ) (
            IndexType row[],
            IndexType pos[],
            const IndexType j,
            const IndexType ellIA[],
            const IndexType numRows,
            const IndexType ellJA[],
            const IndexType numValuesPerRow );

        static const char* getId()
        {
            return "ELL.getValuePosCol";
        }
    };

    struct hasDiagonalProperty
    {
        typedef bool ( *FuncType ) (
            const IndexType numDiagonals,
            const IndexType ellJA[] );

        static const char* getId()
        {
            return "ELL.hasDiagonalProperty";
        }
    };

    struct check
    {
        typedef void ( *FuncType ) (
            const IndexType mNumRows,
            const IndexType mNumValuesPerRow,
            const IndexType mNumColumns,
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const char* msg );

        static const char* getId()
        {
            return "ELL.check";
        }
    };

    /** Define structure for multiplication routines.  */

    template<typename ValueType>
    struct normalGEMV
    {
        /** result = alpha * ELL-Matrix * x + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param ellIA, ellJA, csrValues are arrays of ELL storage
         *  @param op specifies an operation implicitly applied to the matrix storage
         */

        typedef void ( *FuncType ) (
            ValueType result[],
            const ValueType alpha,
            const ValueType x[],
            const ValueType beta,
            const ValueType y[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numValuesPerRow,
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ValueType ellValues[],
            const common::MatrixOp op );

        static const char* getId()
        {
            return "ELL.normalGEMV";
        }
    };

    template<typename ValueType>
    struct sparseGEMV
    {
        /** result += alpha * ELL-Matrix * x, CSR matrix has only some non-zero rows
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param numNonZeroRows is size of rowIndexes
         *  @param rowIndexes are indexes of non-empty rows in matrix
         *  @param ellIA, ellJA, csrValues are arrays of ELL storage
         *
         *  Note: this routine does not provide the term 'beta * y' as it would require
         *        to run over the full result vector
         */

        typedef void ( *FuncType ) (
            ValueType result[],
            const ValueType alpha,
            const ValueType x[],
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const IndexType numNonZeroRows,
            const IndexType rowIndexes[],
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ValueType ellValues[],
            const common::MatrixOp op );

        static const char* getId()
        {
            return "ELL.sparseGEMV";
        }
    };

    /** Structure with type definitions for reduction routines */

    template<typename ValueType>
    struct absMaxVal
    {
        /** This method returns the maximal absolute value of an ELLPACK matrix. */

        typedef ValueType ( *FuncType ) (
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const IndexType ellSizes[],
            const ValueType ellValues[]
        );

        static const char* getId()
        {
            return "ELL.absMaxVal";
        }
    };

    template<typename ValueType>
    struct scaleRows
    {
        typedef void ( *FuncType ) (
            ValueType ellValues[],
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const IndexType ellSizes[],
            const ValueType values[] );

        static const char* getId()
        {
            return "ELL.scaleRows";
        }
    };

    template<typename ValueType>
    struct sortRowElements
    {
        /** This method sorts the elements of a row by increasing column indexes.
         *
         *  @param[in,out] ellJA, ellValues  the CSR matrix data and their column indexes
         *  @param[in]     ellIA             row offsets
         *  @param[in]     numRows           number of rows
         *  @param[in]     numValuesPerRow   maximal number of non-zero values per row
         *  @param[in]     diagonalFlag      if true first entry of each row will be the diagonal element if available
         *
         *  Note: This routine does not force the diagonal property, only if each diagonal element is already available
         */
        typedef void ( *FuncType )(
            IndexType ellJA[],
            ValueType ellValues[],
            const IndexType ellIA[],
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const bool diagonalFlag );

        static const char* getId()
        {
            return "ELL.sortRowElements";
        }
    };

    struct matrixMultiplySizes
    {
        /** @brief Compute the row sizes of result matrix C for matrix multiplication A x B
         *
         *  @param[out] cSizes array of length m, will contain number of entries
         *  @param[in]  m number of rows for matrix C and A
         *  @param[in]  n number of columns for matrix C and B
         *  @param[in]  k number of columns for A and number of rows for B
         *  @param[in]  diagonalProperty if true, diagonal elements will count in any case
         *  @param[in]  aSizes, aJA are the index arrays of matrix A
         *  @param[in]  aNumValuesPerRow multiplied with m is size of aJA
         *  @param[in]  bSizes, bJA are the index arrays of matrix A
         *  @param[in]  bNumValuesPerRow multiplied with k is size of bJA
         *
         *  cNumValuesPerRow can be computed afterwards as maxval( cSizes[0:m-1] )
         *
         *  Note: this routines does not need any value array as only structure is computed
         */

        typedef void ( *FuncType ) ( IndexType cSizes[],
                                     const IndexType m,
                                     const IndexType n,
                                     const IndexType k,
                                     const bool diagonalProperty,
                                     const IndexType aSizes[],
                                     const IndexType aJA[],
                                     const IndexType aNumValuesPerRow,
                                     const IndexType bSizes[],
                                     const IndexType bJA[],
                                     const IndexType bNumValuesPerRow );

        static const char* getId()
        {
            return "ELL.matrixMultiplySizes";
        }
    };

    struct matrixAddSizes
    {
        /** @brief Compute the row sizes of result matrix C for matrix addition A + B
         *
         *  @param[out] cSizes array of length m, will contain number of entries
         *  @param[in]  m number of rows for matrices A, B, and C
         *  @param[in]  n number of columns for matrices A, B, and C
         *  @param[in]  diagonalProperty if true, diagonal elements will count in any case
         *  @param[in]  aSizes, aJA are the index arrays of matrix A
         *  @param[in]  aNumValuesPerRow multiplied with m is size of aJA
         *  @param[in]  bSizes, bJA are the index arrays of matrix A
         *  @param[in]  bNumValuesPerRow multiplied with m is size of bJA
         *
         *  cNumValuesPerRow can be computed afterwards as maxval( cSizes[0:m-1] )
         *
         *  Note: this routines does not need any value array as only structure is computed
         */
        typedef void ( *FuncType ) ( IndexType cSizes[],
                                     const IndexType m,
                                     const IndexType n,
                                     const bool diagonalProperty,
                                     const IndexType aSizes[],
                                     const IndexType aJA[],
                                     const IndexType aNumValuesPerRow,
                                     const IndexType bSizes[],
                                     const IndexType bJA[],
                                     const IndexType bNumValuesPerRow );

        static const char* getId()
        {
            return "ELL.matrixAddSizes";
        }
    };

    template<typename ValueType>
    struct matrixAdd
    {
        /** @brief computes c = alpha * a + beta * b for ELL sparse matrices a, b, c
         *
         * @param[out] cJA is the column index array of c
         * @param[out] cValues is the value array of c
         * @param[in]  m number of rows for all matrices
         * @param[in]  n number of columns for all matrices
         * @param[in]  diagonalProperty if true result matrix c should have diagonal property
         * @param[in]  alpha is scaling factor of the matrix expression
         * @param[in]  aSizes, aJA, aValues, aNumValuesPerRow is data of input matrix a
         * @param[in]  beta is scaling factor of the matrix expression
         * @param[in]  bSizes, bJA, bValues, bNumValuesPerRow is data of input matrix b
         *
         * Note: the size array cValues and cNumValuePerRow must already be available.
         */

        typedef void ( *FuncType ) (
            IndexType cJA[],
            ValueType cValues[],
            const IndexType cSizes[],
            const IndexType cNumValuesPerRow,
            const IndexType m,
            const IndexType n,
            const bool diagonalProperty,
            const ValueType alpha,
            const IndexType aSizes[],
            const IndexType aJA[],
            const ValueType aValues[],
            const IndexType aNumValuesPerRow,
            const ValueType beta,
            const IndexType bSizes[],
            const IndexType bJA[],
            const ValueType bValues[],
            const IndexType bNumValuesPerRow );

        static const char* getId()
        {
            return "ELL.matrixAdd";
        }
    };

    template<typename ValueType>
    struct matrixMultiply
    {
        /** @brief computes c = alpha * a * b for ELL sparse matrices a, b, c
         *
         * @param[out] cJA is the column index array of c
         * @param[out] cValues is the value array of c
         * @param[in]  m number of rows for matrix c and a
         * @param[in]  n number of columns for matrix c and b
         * @param[in]  k number of columns for a and number of rows for b
         * @param[in]  diagonalProperty if true result matrix c should have diagonal property
         * @param[in]  alpha is scaling factor of the matrix expression
         * @param[in]  aSizes, aJA, aValues, aNumValuesPerRow is data of input matrix a
         * @param[in]  bSizes, bJA, bValues, bNumValuesPerRow is data of input matrix b
         *
         * Note: the size array cValues and cNumValuePerRow must already be available.
         */

        typedef void ( *FuncType ) (
            IndexType cJA[],
            ValueType cValues[],
            const IndexType cSizes[],
            const IndexType cNumValuesPerRow,
            const IndexType m,
            const IndexType n,
            const IndexType k,
            const bool diagonalProperty,
            const ValueType alpha,
            const IndexType aSizes[],
            const IndexType aJA[],
            const ValueType aValues[],
            const IndexType aNumValuesPerRow,
            const IndexType bSizes[],
            const IndexType bJA[],
            const ValueType bValues[],
            const IndexType bNumValuesPerRow );

        static const char* getId()
        {
            return "ELL.matrixMultiply";
        }
    };
};

} /* end namespace sparsekernel */

} /* end namespace scai */
