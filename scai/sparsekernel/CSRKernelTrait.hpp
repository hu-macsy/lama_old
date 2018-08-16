/**
 * @file CSRKernelTrait.hpp
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
 * @brief Struct with traits for all CSR utilities provided as kernels.
 * @author Thomas Brandes
 * @date 21.10.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/UnaryOp.hpp>
#include <scai/common/MatrixOp.hpp>

namespace scai
{

/** Namespace for all kernels belonging to operations for matrix storage. */

namespace sparsekernel
{

/** @brief Traits for kernel functions to be used in CSR storage.  */

struct CSRKernelTrait
{
    struct getValuePos
    {
        /** Returns position of element (i,j) in ja/values array
         *
         *  @param[in] i is the row of the element
         *  @param[in] j is the column of the element
         *  @param[in] csrIA is the CSR offset array
         *  @param[in] csrJA is the CSR ja array
         *  @returns  offset of element in values array, invalidIndex if not found
         */

        typedef IndexType ( *FuncType ) (
            const IndexType i,
            const IndexType j,
            const IndexType csrIA[],
            const IndexType csrJA[] );

        static const char* getId()
        {
            return "CSR.getValuePos";
        }
    };

    struct getColumnPositions
    {
        /** This method returns for a certain column of the CSR matrix all
         *  row indexes for which elements exist and the corresponding positions
         *  in the csrJA/csrValues array
         *
         *  @param[out] row indexes of rows that have an entry for column j
         *  @param[out] pos positions of entries in csrJA[pos[i]] == j,
         *  @param[in] j is the column of which positions are required
         *  @param[in] csrIA is the CSR offset array
         *  @param[in] numRows is the number of rows
         *  @param[in] csrJA is the CSR ja array
         *  @param[in] numValues is the number of non-zero values
         *  @returns  number of entries with col index = j
         */
        typedef IndexType ( *FuncType ) (
            IndexType row[],
            IndexType pos[],
            const IndexType j,
            const IndexType csrIA[],
            const IndexType numRows,
            const IndexType csrJA[],
            const IndexType numValues );

        static const char* getId()
        {
            return "CSR.getColumnPositions";
        }
    };

    /** Structure defining function types for operations on CSR data
     *
     *  @tparam ValueType is the value type of the matrix element, e.g. float, double
     */

    template<typename ValueType>
    struct sortRows
    {
        /** This method sorts the elements of a row by increasing column indexes.
         *
         *  @param[in,out] csrJA, csrValues  the CSR matrix data and their column indexes
         *  @param[in]     csrIA             row offsets
         *  @param[in]     numRows           number of rows
         *  @param[in]     numColumns        number of columns
         *  @param[in]     numValues         number of stored non-zero entries, size of csrJA and csrValues
         *
         *  Note: This routine does not force the diagonal property, only if each diagonal element is already available
         */
        typedef void ( *FuncType )(
            IndexType csrJA[],
            ValueType csrValues[],
            const IndexType csrIA[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numValues );

        static const char* getId()
        {
            return "CSR.sortRows";
        }
    };

    struct hasSortedRows
    {
        /** This method checks for sorted column entries:
         */
        typedef bool ( *FuncType )(
            const IndexType csrIA[],
            const IndexType csrJA[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numValues );

        static const char* getId()
        {
            return "CSR.hasSortedRows";
        }
    };

    /** Kernel definition of jacobi iteration with cSR data */

    template <typename ValueType>
    struct jacobi
    {
        /** Method to compute one iteration step in Jacobi method
         *
         *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
         *
         *  where B is the storage without diagonal and dinv the inverse of the diagonal.
         */
        typedef void ( *FuncType ) (
            ValueType solution[],
            const IndexType csrIA[],
            const IndexType csrJA[],
            const ValueType csrValues[],
            const ValueType oldSolution[],
            const ValueType rhs[],
            const ValueType omega,
            const IndexType numRows );

        static const char* getId()
        {
            return "CSR.jacobi";
        }
    };

    template <typename ValueType>
    struct jacobiHalo
    {
        /** Compute one iteration step in Jacobi method for halo
         *
         *  \code
         *      solution -= omega * ( csr_halo * oldSolution ) ./ diagonal 
         *  \endcode
         *
         */
        typedef void ( *FuncType ) (
            ValueType solution[],
            const ValueType diagonal[],
            const IndexType haloIA[],
            const IndexType haloJA[],
            const ValueType haloValues[],
            const IndexType haloRowIndexes[],
            const ValueType oldSolution[],
            const ValueType omega,
            const IndexType numNonEmptyRows );

        static const char* getId()
        {
            return "CSR.jacobiHalo";
        }
    };

    /* LU factorization with Pardiso(MKL)/cusolver */
    template <typename ValueType>
    struct decomposition
    {
        /** Direct solving of linear equations
         *
         *  @since 2.1.0
         */
        typedef void ( *FuncType ) (
            ValueType solution[],
            const IndexType csrIA[],
            const IndexType csrJA[],
            const ValueType csrValues[],
            const ValueType rhs[],
            const IndexType numRows,
            const IndexType nnz,
            const bool isSymmetic );

        static const char* getId()
        {
            return "CSR.decomposition";
        }
    };

    /** Structure with type definitions for offset routines. */

    struct sizes2offsets
    {
        /** This method makes an offset array from the sizes.
         *
         *  @param[in,out] array contains counter values and later the offsets
         *  @param[in]    n is the number of values, array must contain one additional value
         *  @returns      the total number of values
         *
         *  \code
         *    array  :    3    7   8   4   2  x
         *    array  :    0   10  15  12  16  18  -> returns 18
         *  \endcode
         *
         *  Important: sizes must have numRows + 1 allocated entries.
         *
         */

        typedef IndexType ( *FuncType ) ( IndexType array[], const IndexType n );

        static const char* getId()
        {
            return "CSR.sizes2offsets";
        }
    };

    struct offsets2sizes
    {
        /** This method computes size array from an offset array.
         *
         *  @param[out] sizes will contain the sizes (e.g. for each row ), has numRows entries
         *  @param[in] offsets contains the offsets, has numRows + 1 entries
         *  @param[in] n is size of array sizes, n + 1 is size of array offsets
         *
         *  \code
         *           offsets :   0   5  11   16  19  19  21  23
         *           sizes   :   5   6   5    3   0   2   2
         *  \endcode
         *
         */
        typedef void ( *FuncType ) ( IndexType sizes[], const IndexType offsets[], const IndexType n );

        static const char* getId()
        {
            return "CSR.offsets2sizes";
        }
    };

    struct gatherSizes
    {
        /** This method gets the sizes for a set of rows
         *
         *  @param[out] sizes will contain the sizes (e.g. for each row ), has n entries
         *  @param[in] csrIA contains the offsets, has numRows + 1 entries
         *  @param[in] numRows might be used to check for valid row indexes
         *  @param[in] rowIndexes an array of row indexes for which the size is neeted
         *  @param[in] nIndexes is size of rowIndexes and sizes
         *
         *  \code
         *           csrIA   :   0   5  11   16  19  19  21  23
         *           rows    :   0   2  4   6
         *           sizes   :   5   5  0   2
         *  \endcode
         *
         */
        typedef void ( *FuncType ) ( 
            IndexType sizes[], 
            const IndexType csrIA[], 
            const IndexType numRows,
            const IndexType rowIndexes[], 
            const IndexType nIndexes );

        static const char* getId()
        {
            return "CSR.gatherSizes";
        }
    };

    struct validOffsets
    {
        /** Check for a legal offset array.
         *
         *  @param[in] array is the array with offsets
         *  @param[in] n array has n + 1 entries
         *  @param[in] total must be last entry in array
         *  @return    true if the array is legal
         *
         *  Means: array[0] <= array[1] <= array[2] <= ... <= array[n] (= total)
         */

        typedef bool ( *FuncType ) ( const IndexType array[], const IndexType n, const IndexType total );

        static const char* getId()
        {
            return "CSR.validOffsets";
        }
    };

    struct nonEmptyRows
    {
        /** Bild a vector of indexes for non-empty rows.
         *
         *  @param[out] rowIndexes will contain all indexes of non-empty rows
         *  @param[in]  csrIA is the offset array
         *  @param[in]  n is the number of rows, offsets has size n + 1
         *  @return     the number of non-empty rowssets has size n + 1
         *  
         *  With rowIndexes == NULL the function can be used just to determine the number 
         *  of non-empty rows.
         */

        typedef IndexType ( *FuncType ) (
            IndexType rowIndexes[],
            const IndexType csrIA[],
            const IndexType numRows );

        static const char* getId()
        {
            return "CSR.nonEmptyRows";
        }
    };

    struct matrixAddSizes
    {
        /** This method computes the row sizes for result matrix C of matrix add A + B
         *
         *  @param[out] cIa array of length numRows, will contain number of entries in each row for C
         *  @param[in]  numRows number of rows for matrix A and B
         *  @param[in]  numColumns number of columns for matrix A and B
         *  @param[in]  aIA, aJA are the index arrays of matrix A
         *  @param[in]  bIA, bJA are the index arrays of matrix B
         *
         *  Note: filling the result matrix must use the same flag for diagonalProperty
         *        otherwise the row sizes/offsets will not match
         */
        typedef IndexType ( *FuncType ) (
            IndexType cIa[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType aIA[],
            const IndexType aJA[],
            const IndexType bIA[],
            const IndexType bJA[] );

        static const char* getId()
        {
            return "CSR.matrixAddSizes";
        }
    };

    struct binaryOpSizes
    {
        /** This method computes the offset array for result matrix C = matrix A + matrix B
         *
         *  @param[out] cIA is the offset array for the result matrix
         *  @param[in]  numRows number of rows for matrix A and B
         *  @param[in]  numColumns number of columns for matrix A and B
         *  @param[in]  aIA is the row offset array of matrix A
         *  @param[in]  aJA contains the column indexes of the non-zero entries
         *  @param[in]  bIA, bJA are the corresponding CSR arrays of matrix B
         *
         *  Important: this routine relies on sorted entries (by column index) for each row
         */
        typedef IndexType ( *FuncType ) (
            IndexType cIA[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType aIA[],
            const IndexType aJA[],
            const IndexType bIA[],
            const IndexType bJA[] );

        static const char* getId()
        {
            return "CSR.binaryOpSizes";
        }
    };

    struct matrixMultiplySizes
    {
        /** This method computes the row sizes for result matrix C of matrix multiplication A x B
         *
         *  @param[out] cSizes array of length numRows, will contain number of entries
         *  @param[in]  m number of rows for matrix C and A
         *  @param[in]  n number of columns for matrix C and B
         *  @param[in]  k number of columns for A and number of rows for B
         *  @param[in]  diagonalProperty if true, diagonal elements will count in any case
         *  @param[in]  aIA, aJA are the index arrays of matrix A
         *  @param[in]  bIA, bJA are the index arrays of matrix B
         */

        typedef IndexType ( *FuncType ) (
            IndexType cSizes[],
            const IndexType m,
            const IndexType n,
            const IndexType k,
            const IndexType aIA[],
            const IndexType aJA[],
            const IndexType bIA[],
            const IndexType bJA[] );

        static const char* getId()
        {
            return "CSR.matrixMultiplySizes";
        }
    };

    struct hasDiagonalProperty
    {
        /** This method checks whether the CSR structure data has the diagonal property.
         *
         *  @param[in] numDiagonals  number of first rows for which diagonal property is checked.
         *  @param[in] csrIA         offset array for the rows
         *  @param[in] csrJA         column indexes
         *  @return                  true if diagonal property is given
         *
         *  The diagonal property is given if the first column index in the row is same as the row index.
         */
        typedef bool ( *FuncType ) (
            const IndexType numDiagonals,
            const IndexType csrIA[],
            const IndexType csrJA[],
            bool isSorted );

        static const char* getId()
        {
            return "CSR.hasDiagonalProperty";
        }
    };

    template<typename ValueType>
    struct getDiagonal
    {
        /** This method returns the diagonal of a csr storage.
         *
         *  @param[in] numDiagonals number of diagonal elements, should be min( numRows, numColumns )
         *  @param[in] csrIA, csrJA, csrValues are the CSR arrays
         *  @param[in] isSorted if true the column entries for each row are sorted
         *
         *  Note: if sort flag has been set but elements were not sorted, the diagonal might have
         *        wrong zero elements.
         */
        typedef void ( *FuncType ) (
            ValueType diagonal[],
            const IndexType numDiagonals,
            const IndexType csrIA[],
            const IndexType csrJA[],
            const ValueType csrValues[],
            const bool isSorted );

        static const char* getId()
        {
            return "CSR.getDiagonal";
        }
    };

    template<typename ValueType>
    struct setDiagonalV
    {
        /** This method sets the diagonal in a CSR storage.
         *
         *  Returns true if all diagonal entries were available in the CSR storage.
         *
         *  Note: if sort flag has been set but column indexes were not sorted, there 
         *        might be diagonal elemens not set even if they were available.
         */
        typedef bool ( *FuncType ) (
            ValueType csrValues[],
            const ValueType diagonal[],
            const IndexType numDiagonals,
            const IndexType csrIA[],
            const IndexType csrJA[],
            const bool isSorted );

        static const char* getId()
        {
            return "CSR.setDiagonalV";
        }
    };

    template<typename ValueType>
    struct setDiagonal
    {
        /** This method determines sets the diagonal of a csr storage.
         *
         *  Returns true if all diagonal elements were available.
         */
        typedef bool ( *FuncType ) (
            ValueType csrValues[],
            const ValueType diagonal,
            const IndexType numDiagonals,
            const IndexType csrIA[],
            const IndexType csrJA[],
            const bool isSorted );

        static const char* getId()
        {
            return "CSR.setDiagonal";
        }
    };

    struct getPosDiagonal
    {
        /** This method determines the positions of the diagonal elements in csr arrays
         *
         *  @param[out] pos  will contain the positions in csrJA/csrValues of the i-th diagonal element
         *  @param[in] numDiagonals  number of diagonals 
         *  @param[in] csrIA         offset array for the rows
         *  @param[in] csrJA         column indexes
         *  @param[in] isSorted      if true the column indexes are sorted
         *  @return                  number of available diagonal elements
         *
         *  If a diagonal element is not available, the corresponding entry is invalidIndex.
         */
        typedef IndexType ( *FuncType ) (
            IndexType pos[],
            const IndexType numDiagonals,
            const IndexType csrIA[],
            const IndexType csrJA[],
            const bool isSorted );

        static const char* getId()
        {
            return "CSR.getPosDiagonal";
        }
    };

    template<typename ValueType>
    struct shiftDiagonal
    {
        /** This method sets the diagonal entries of a row as first entry wherever possible.
         *
         */
        typedef IndexType ( *FuncType ) (
            IndexType csrJA[], 
            ValueType csrValues[],
            const IndexType numDiagonals,
            const IndexType csrIA[] );

        static const char* getId()
        {
            return "CSR.shiftDiagonal";
        }
    };

    /** Define structure that contains type definitions for the function pointers.
     *
     *  @tparam ValueType specifies the value type used in the arrays.
     *
     *  The structure is needed as type definition templates are unsupported in C++.
     */

    template<typename ValueType>
    struct convertCSR2CSC
    {
        /** Function pointer for CSR to CSC conversion routine.
         *
         *  @param[out] cscIA, cscJA, cscValues is CSC output data
         *  @param[in]  csrIA, csrJA, csrValues is CSR input data
         *  @param numRows x numColumns is shape of input CSR matrix
         *
         *  Arrays must be big enough for the corresponing sizes.
         */

        typedef void( *FuncType ) ( IndexType cscIA[],
                                    IndexType cscJA[],
                                    ValueType cscValues[],
                                    const IndexType csrIA[],
                                    const IndexType csrJA[],
                                    const ValueType csrValues[],
                                    IndexType numRows, IndexType numColumns,
                                    IndexType numValues );

        static const char* getId()
        {
            return "CSR.convertCSR2CSC";
        }
    };

    template<typename ValueType>
    struct setRows
    {
        /** This operation multiplies each row with an own value.
         *
         *  @param[in,out] csrValues matrix data that is scaled
         *  @param[in]     csrIA offset array to identify which elements of csrValues belong to which row
         *  @param[in]     numRows number of rows in matrix and also size of values
         *  @param[in]     values array with element for each row used for scaling
         *  @param[in]     op specfies how each element is updated with row value
         *
         *  csr[i,j] *= values[i], for i = 0, ..., numRows-1
         *
         */
        typedef void ( *FuncType ) (
            ValueType csrValues[],
            const IndexType csrIA[],
            const IndexType numRows,
            const ValueType values[],
            common::BinaryOp op );

        static const char* getId()
        {
            return "CSR.setRows";
        }
    };

    template<typename ValueType>
    struct setColumns
    {
        /** This operation applies a binary operation to each element, individually for each column
         *
         *  @param[in,out] csrValues matrix data that is scaled
         *  @param[in]     csrIA offset array to identify which elements of csrValues belong to which row
         *  @param[in]     csrJA column indexes
         *  @param[in]     numRows number of rows in matrix
         *  @param[in]     values array with element for each row used for scaling
         *  @param[in]     op specfies how each element is updated with row value
         *
         *  csr[i,j] op= values[j], for i = 0, ..., numRows-1, j = 0, ..., numColumns-1
         *
         */
        typedef void ( *FuncType ) (
            ValueType csrValues[],
            const IndexType csrIA[], 
            const IndexType csrJA[], 
            const IndexType numRows,
            const ValueType values[],
            common::BinaryOp op );

        static const char* getId()
        {
            return "CSR.setColumns";
        }
    };

    /** Structure with type definitions for reduction routines. */

    template <typename ValueType>
    struct absMaxDiffVal
    {
        /** Get the maximal element-wise difference for two CSR matrices.
         *
         *  @param[in] numRows is number of rows for both matrices
         *  @param[in] sortedRows if true column indexes in JA arrays are sorted
         *  @param[in] csrIA1, csrJA1, csrValues1 is storage data of first matrix
         *  @param[in] csrIA2, csrJA2, csrValues2 is storage data of second matrix
         *  @returns maximal value of absolute difference between two matrix elements
         */

        typedef RealType<ValueType> ( *FuncType ) (
            IndexType numRows,
            bool sortedRows,
            const IndexType csrIA1[],
            const IndexType csrJA1[],
            const ValueType csrValues1[],
            const IndexType csrIA2[],
            const IndexType csrJA2[],
            const ValueType csrValues2[] );

        static const char* getId()
        {
            return "CSR.absMaxDiffVal";
        }
    };

    /** Define structure for reduce column/row distribution
     *
     *  @tparam ValueType specifies the value type used in data
     */
    template<typename ValueType>
    struct reduce
    {
        /* result = reduce( csrMatrix , dim, reduceOp, elemOp )
         *
         *  @param result is the result vector, size is numRows for dim = 0, numCols for dim = 1
         *  @param numRows is number of rows of matrix
         *  @param csrIA, csrJA, csrValues are arrays of CSR storage
         *  @param reduceOp is how elements of one row or one column are reduced
         *  @param elemOp is operator applied to the elements
         */
        typedef void ( *FuncType ) (
            ValueType result[],
            const IndexType csrIA[],
            const IndexType csrJA[],
            const ValueType csrValues[],
            const IndexType numRows,
            const IndexType dim,
            const common::BinaryOp reduceOp,
            const common::UnaryOp elemOp );

        static const char* getId()
        {
            return "CSR.reduce";
        }
    };

    /** Define structure for multiplication routines.
     *
     *  @tparam ValueType specifies the value type used in mutliplications.
     */

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
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param numColums is the number of columns of matrix
         *  @param nnz number of nonZeros, same as csrIA[ numRows ]
         *  @param csrIA, csrJA, csrValues are arrays of CSR storage
         *  @param op implicit operation applied to the matrix storage
         *
         *  Note: result, y have size numRows / numCols for normal / transpose
         *        x has size numColumns / numrows for normal / transpose
         */

        typedef void ( *FuncType ) ( ValueType result[],
                                     const ValueType alpha,
                                     const ValueType x[],
                                     const ValueType beta,
                                     const ValueType y[],
                                     const IndexType numRows,
                                     const IndexType numColumns,
                                     const IndexType nnz,
                                     const IndexType csrIA[],
                                     const IndexType csrJA[],
                                     const ValueType csrValues[],
                                     const common::MatrixOp op );

        static const char* getId()
        {
            return "CSR.normalGEMV";
        }
    };

    template<typename ValueType>
    struct sparseGEMV
    {
        /** result += alpha * CSR-Matrix * x, CSR matrix has only some non-zero rows
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param numNonZeroRows is size of rowIndexes
         *  @param rowIndexes are indexes of non-empty rows in matrix
         *  @param csrIA, csrJA, csrValues are arrays of CSR storage
         *  @param op specifies an implicit operation applied to the matrix 
         *
         *  Note: this routine does not provide the term 'beta * y' as it would require
         *        to run over the full result vector
         */

        typedef void ( *FuncType ) (
            ValueType result[],
            const ValueType alpha,
            const ValueType x[],
            const IndexType numNonZeroRows,
            const IndexType rowIndexes[],
            const IndexType csrIA[],
            const IndexType csrJA[],
            const ValueType csrValues[],
            const common::MatrixOp op );

        static const char* getId()
        {
            return "CSR.sparseGEMV";
        }
    };

    template<typename ValueType>
    struct gemm
    {
        /**  This method computes result = alpha * CSR * x + beta * y  with dense result, x, y
         *
         *   @param[out] result  has size numRows x k
         *   @param[in]  alpha scaling factor
         *   @param[in]  x has size numColumns x k 
         *   @param[in]  beta scaling factor
         *   @param[in]  y has size numRows x k
         *   @param[in]  csrIA, csrJA, csrValues are sparse matrix storage of numRows x numColumns
         */

        typedef void ( *FuncType ) (
            ValueType result[],
            const ValueType alpha,
            const ValueType x[],
            const ValueType beta,
            const ValueType y[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType k,
            const IndexType csrIA[],
            const IndexType csrJA[],
            const ValueType csrValues[],
            const common::MatrixOp op );

        static const char* getId()
        {
            return "CSR.gemm";
        }
    };

    template<typename ValueType>
    struct binaryOp
    {
        /** computes c = a 'op' b for CSR sparse storages a, b, c
         *
         *  @param[out] cJA, cValues are the matrix values of output matrix
         *  @param[in]  cIA contains already computed offsets
         *  @param[in]  numRows is number of rows for matrices a, b, c
         *  @param[in]  numColums is number of columns for matrices a, b, c
         *  @param[in]  aIA, aJA, aValues is input matrix a in CSR format
         *  @param[in]  bIA, bJA, bValues is input matrix b in CSR format
         *  @param[in]  op specifies the binary operation
         *
         *  In this routine the row offsets of C must already be determined
         *  before.
         */

        typedef void ( *FuncType ) (
            IndexType cJA[],
            ValueType cValues[],
            const IndexType cIA[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType aIA[],
            const IndexType aJA[],
            const ValueType aValues[],
            const IndexType bIA[],
            const IndexType bJA[],
            const ValueType bValues[],
            const common::BinaryOp op );

        static const char* getId()
        {
            return "CSR.binaryOp";
        }
    };

    template<typename ValueType>
    struct matrixAdd
    {
        /** computes c = alpha * a + beta * b for CSR sparse matrices a, b, c
         *
         *  @param[out] cJA, cValues are the matrix values of output matrix
         *  @param[in]  cIA contains already computed offsets
         *  @param[in]  numRows is number of rows for matrices a, b, c
         *  @param[in]  numColums is number of columns for matrices a, b, c
         *  @param[in]  alpha is a scalar scaling factor for matrix a
         *  @param[in]  aIA, aJA, aValues is input matrix a in CSR format
         *  @param[in]  beta is a scalar scaling factor for matrix b
         *  @param[in]  bIA, bJA, bValues is input matrix b in CSR format
         *
         *  In this routine the row offsets of C must already be determined
         *  before.
         */

        typedef void ( *FuncType ) (
            IndexType cJA[],
            ValueType cValues[],
            const IndexType cIA[],
            const IndexType numRows,
            const IndexType numColumns,
            const ValueType alpha,
            const IndexType aIA[],
            const IndexType aJA[],
            const ValueType aValues[],
            const ValueType beta,
            const IndexType bIA[],
            const IndexType bJA[],
            const ValueType bValues[] );

        static const char* getId()
        {
            return "CSR.matrixAdd";
        }
    };

    template<typename ValueType>
    struct matrixMultiply
    {
        /** computes c = alpha * a * b for CSR sparse matrices a, b, c
         *
         *  @param[out] cValues are the matrix values of output matrix
         *  @param[in]  cIA, cJA contain structure of output matrix
         *  @param[in]  m number of rows for matrix c and a
         *  @param[in]  n number of columns for matrix c and b
         *  @param[in]  k number of columns for a and number of rows for b
         *  @param[in]  alpha is a scalar scaling factor
         *  @param[in]  aIA, aJA, aValues is input matrix a in CSR format
         *  @param[in]  bIA, bJA, bValues is input matrix b in CSR format
         *
         *  In this routine the final structure of C must already be determined
         *  before. Only those values of C are computed for which cIA, cJA specify
         *  an available entry.
         */

        typedef void ( *FuncType ) (
            const IndexType cIa[],
            IndexType cJA[],
            ValueType cValues[],
            const IndexType m,
            const IndexType n,
            const IndexType k,
            const ValueType alpha,
            const IndexType aIA[],
            const IndexType aJA[],
            const ValueType aValues[],
            const IndexType bIA[],
            const IndexType bJA[],
            const ValueType bValues[] );

        static const char* getId()
        {
            return "CSR.matrixMultiply";
        }
    };

    template<typename ValueType>
    struct countNonZeros
    {
        /** Count non zero entries in CSR storage after matrix operation like add, mult
         *
         * @param[out] sizes are the row sizes for the compressed data
         * @param[in] ia, ja, values, numRows are the data of the current CSR storage
         * @param[in] eps     threshold value for which an element is considered to be zero
         */
        typedef void ( *FuncType )(
            IndexType sizes[],
            const IndexType ia[],
            const IndexType ja[],
            const ValueType values[],
            const IndexType numRows,
            const RealType<ValueType> eps );

        static const char* getId()
        {
            return "CSR.countNonZeros";
        }
    };

    template<typename ValueType>
    struct compress
    {
        /** Fill compressed CSR data in new data structures
         *
         * @param[out] newJA, newValues column indexes and data of the new CSR data
         * @param[in] newIA   new offsets, computed by countNonZeros + sizes2offsets
         * @param[in] ia, ja, values, numRows are the data of the current CSR storage
         * @param[in] eps     threshold value for which an element is considered to be zero
         */
        typedef void ( *FuncType )(
            IndexType newJA[],
            ValueType newValues[],
            const IndexType newIA[],
            const IndexType ia[],
            const IndexType ja[],
            const ValueType values[],
            const IndexType numRows,
            const RealType<ValueType> eps );

        static const char* getId()
        {
            return "CSR.compress";
        }
    };
};

} /* end namespace sparsekernel */

} /* end namespace scai */
