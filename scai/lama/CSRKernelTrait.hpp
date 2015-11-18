/**
 * @file CSRKernelTrait.hpp
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
 * @brief Struct with traits for all CSR utilities provided as kernels.
 * @author Thomas Brandes
 * @date 21.10.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace lama
{

/** @brief Traits for kernel functions to be used in CSR storage.  */

struct CSRKernelTrait
{
    /** Structure defining function types for operations on CSR data
     *
     *  @tparam ValueType is the value type of the matrix element, e.g. float, double
     */

    template<typename ValueType>
    struct sortRowElements
    {
        /** This metod sorts the elements of a row by increasing column indexes.
         *
         *  @param[in,out] csrJA, csrValues  the CSR matrix data and their column indexes
         *  @param[in]     csrIA             row offsets
         *  @param[in]     numRows           number of rows
         *  @param[in]     diagonalFlag      if true first entry of each row will be the diagonal element if available
         *
         *  Note: This routine does not force the diagonal property, only if each diagonal element is already available
         */
        typedef void (*FuncType)(
            IndexType csrJA[],
            ValueType csrValues[],
            const IndexType csrIA[],
            const IndexType numRows,
            const bool diagonalFlag );

        static const char* getId() { return "CSR.sortRowElements"; }
    };

    /** Structure with type definitions for solver routines */

    template <typename ValueType>
    struct jacobi
    {
        /** Method to compute one iteration step in Jacobi method
         *
         *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
         *
         */
        typedef void ( *FuncType ) ( ValueType solution[],
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const ValueType csrValues[],
                        const ValueType oldSolution[],
                        const ValueType rhs[],
                        const ValueType omega,
                        const IndexType numRows );

        static const char* getId() { return "CSR.jacobi"; }
    };

    template <typename ValueType>
    struct jacobiHalo
    {
        /** Method to compute one iteration step in Jacobi method
         *
         *  solution -= omega * ( B(halo) * oldSolution) * dinv
         *
         */
        typedef void ( *FuncType ) ( ValueType solution[],
                        const IndexType localIA[],
                        const ValueType localValues[],
                        const IndexType haloIA[],
                        const IndexType haloJA[],
                        const ValueType haloValues[],
                        const IndexType haloRowIndexes[],
                        const ValueType oldSolution[],
                        const ValueType omega,
                        const IndexType numNonEmptyRows );

        static const char* getId() { return "CSR.jacobiHalo"; }
    };

    template <typename ValueType>
    struct jacobiHaloWithDiag
    {
        /** Method to compute one iteration step in Jacobi method
         *
         *  solution -= omega * ( B(halo) * oldSolution) * dinv
         *
         *  @since 1.1.0
         */
        typedef void ( *FuncType ) ( ValueType solution[],
                        const ValueType localDiagValues[],
                        const IndexType haloIA[],
                        const IndexType haloJA[],
                        const ValueType haloValues[],
                        const IndexType haloRowIndexes[],
                        const ValueType oldSolution[],
                        const ValueType omega,
                        const IndexType numNonEmptyRows );

        static const char* getId() { return "CSR.jacobiHaloWithDiag"; }
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

        static const char* getId() { return "CSR.sizes2offsets"; }
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

        static const char* getId() { return "CSR.offsets2sizes"; }
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

        static const char* getId() { return "CSR.validOffsets"; }
    };

    struct matrixAddSizes
    {
        /** This method computes the row sizes for result matrix C of matrix add A + B
         *
         *  @param[out] cIa array of length numRows, will contain number of entries in each row for C
         *  @param[in]  numRows number of rows for matrix A and B
         *  @param[in]  numColumns number of columns for matrix A and B
         *  @param[in]  diagonalProperty if true, diagonal elements will count in any case
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
            bool diagonalProperty,
            const IndexType aIA[], 
            const IndexType aJA[],
            const IndexType bIA[], 
            const IndexType bJA[] );

        static const char* getId() { return "CSR.matrixAddSizes"; }
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

        typedef IndexType ( *FuncType ) ( IndexType cSizes[],
                        const IndexType m,
                        const IndexType n,
                        const IndexType k,
                        bool diagonalProperty,
                        const IndexType aIA[], const IndexType aJA[],
                        const IndexType bIA[], const IndexType bJA[] );

        static const char* getId() { return "CSR.matrixMultiplySizes"; }
    };

    struct matrixMultiplyJA
    {
        /** This method computes the column indexes for result matrix C of matrix multiplication A x B
         *
         *  @param[out] cJA array will contain the column indexes, size is cIA[numRows]
         *  @param[in]  CIA array with row offsets as computed by matrixMultiplySizes + sizes2offsets
         *  @param[in]  numRows number of rows for matrix C and A
         *  @param[in]  numColumns number of columns for matrix C and B
         *  @param[in]  diagonalProperty if true, diagonal elements will filled in at begin of each row
         *  @param[in]  aIA, aJA are the index arrays of matrix A
         *  @param[in]  bIA, bJA are the index arrays of matrix B
         */

        typedef void ( *FuncType ) ( IndexType cJA[], const IndexType cIA[],
                        const IndexType numRows, const IndexType numColumns,
                        bool diagonalProperty,
                        const IndexType aIA[], const IndexType aJA[],
                        const IndexType bIA[], const IndexType bJA[] );

        static const char* getId() { return "CSR.matrixMultiplyJA"; }
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
        typedef bool ( *FuncType ) ( const IndexType numDiagonals,
                        const IndexType csrIA[],
                        const IndexType csrJA[] );

        static const char* getId() { return "CSR.hasDiagonalProperty"; }
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

        static const char* getId() { return "CSR.convertCSR2CSC"; }
    };

    /** Define structure for multiplication routines.  */

    template<typename ValueType1, typename ValueType2>
    struct scaleRows
    {

        /** This operation multiplies each row with an own value.
         *
         *  @param[in,out] csrValues matrix data that is scaled
         *  @param[in]     csrIA offset array to identify which elements of csrValues belong to which row
         *  @param[in]     numRows number of rows in matrix and also size of values
         *  @param[in]     values array with element for each row used for scaling
         *
         *  csr[i,j] *= values[i], for i = 0, ..., numRows-1
         *
         *  This routine supports different precision for matrix values and scale values.
         */

        typedef void ( *FuncType ) ( ValueType1 csrValues[],
                        const IndexType csrIA[],
                        const IndexType numRows,
                        const ValueType2 values[] );

        static const char* getId() { return "CSR.scaleRows"; }
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

        typedef ValueType ( *FuncType ) ( IndexType numRows, bool sortedRows,
                        const IndexType csrIA1[], const IndexType csrJA1[], const ValueType csrValues1[],
                        const IndexType csrIA2[], const IndexType csrJA2[], const ValueType csrValues2[] );

        static const char* getId() { return "CSR.absMaxDiffVal"; }
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
         *
         *  The number of columns is not really needed to implement the operation but
         *  might be helpful to choose between different implementations.
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
                        const ValueType csrValues[] );

        static const char* getId() { return "CSR.normalGEMV"; }
    };

    template<typename ValueType>
    struct normalGEVM
    {
        typedef void ( *FuncType ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const ValueType beta,
                        const ValueType y[],
                        const IndexType numRows,
                        const IndexType numColumns,
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const ValueType csrValues[] );

        static const char* getId() { return "CSR.normalGEVM"; }
    };

    template<typename ValueType>
    struct sparseGEMV
    {
        /** result = alpha * CSR-Matrix * x, CSR matrix has only some non-zero rows
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param numNonZeroRows is size of rowIndexes
         *  @param rowIndexes are indexes of non-empty rows in matrix
         *  @param csrIA, csrJA, csrValues are arrays of CSR storage
         *
         *  Note: this routine does not provide the term 'beta * y' as it would require
         *        to run over the full result vector
         */

        typedef void ( *FuncType ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const IndexType numNonZeroRows,
                        const IndexType rowIndexes[],
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const ValueType csrValues[] );

        static const char* getId() { return "CSR.sparseGEMV"; }
    };

    template<typename ValueType>
    struct sparseGEVM
    {
        typedef void ( *FuncType ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const IndexType numColumns,
                        const IndexType numNonZeroRows,
                        const IndexType rowIndexes[],
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const ValueType csrValues[] );

        static const char* getId() { return "CSR.sparseGEMV"; }
    };

    template<typename ValueType>
    struct gemm
    {
        /**  This method computes result = alpha * CSR * x + beta * y  with dense result, x, y
         *
         *   @param[out] result  has size m x n
         *   @param[in]  alpha scaling factor
         *   @param[in]  x has size p x n
         *   @param[in]  beta scaling factor
         *   @param[in]  y has size m x n
         *   @param[in]  csrIA offset array of CSR matrix, has size m + 1
         *   @param[in]  csrJA has size csrIA[m], values between 0 and p-1
         *   @param[in]  csrVaues is value array of CSR matrix
         */

        typedef void ( *FuncType ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const ValueType beta,
                        const ValueType y[],
                        const IndexType m,
                        const IndexType n,
                        const IndexType p,
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const ValueType csrValues[] );

        static const char* getId() { return "CSR.gemm"; }
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
         *  @param[in]  diagonalProperty if true result matrix gets diagonal property
         *  @param[in]  alpha is a scalar scaling factor for matrix a
         *  @param[in]  aIA, aJA, aValues is input matrix a in CSR format
         *  @param[in]  beta is a scalar scaling factor for matrix b
         *  @param[in]  bIA, bJA, bValues is input matrix b in CSR format
         *
         *  In this routine the row offsets of C must already be determined
         *  before.
         */

        typedef void ( *FuncType ) ( IndexType cJA[],
                        ValueType cValues[],
                        const IndexType cIA[],
                        const IndexType numRows,
                        const IndexType numColumns,
                        const bool diagonalProperty,
                        const ValueType alpha,
                        const IndexType aIA[],
                        const IndexType aJA[],
                        const ValueType aValues[],
                        const ValueType beta,
                        const IndexType bIA[],
                        const IndexType bJA[],
                        const ValueType bValues[] );

        static const char* getId() { return "CSR.matrixAdd"; }
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

        typedef void ( *FuncType ) ( const IndexType cIa[],
                        IndexType cJA[],
                        ValueType cValues[],
                        const IndexType m,
                        const IndexType n,
                        const IndexType k,
                        const ValueType alpha,
                        bool diagonalProperty,
                        const IndexType aIA[],
                        const IndexType aJA[],
                        const ValueType aValues[],
                        const IndexType bIA[],
                        const IndexType bJA[],
                        const ValueType bValues[] );

        static const char* getId() { return "CSR.matrixMultiply"; }
    };
};

    template<typename ValueType>
    struct scale
    {
        /** @brief scale array of values with a value in place
         *
         *  @param[in,out]  values is the array with entries to scale
         *  @param[in]      value  is the scaling factor
         *  @param[in]      n      is the number of entries in values
         */
        typedef void ( *FuncType ) ( ValueType values[],
                        const ValueType value,
                        const IndexType n );

        static const char* getId() { return "scale"; }
    };

} /* end namespace lama */

} /* end namespace scai */
