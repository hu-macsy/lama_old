/**
 * @file UtilsInterface.hpp
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
 * @brief Interface classes for utilities and matrix storage routines
 * @author Thomas Brandes
 * @date 03.04.2013
 * @since 1.0.0
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/lama/BLASInterface.hpp>

// internal scai libraries
#include <scai/common/macros/interface.hpp>

namespace scai
{

namespace tasking
{
    class SyncToken;
}

// forward declaration
namespace lama
{


/** Structure with pointers for all Utils methods. */

struct UtilsInterface
{
    /** @brief Trait for register kernel function validIndexes */

    struct validIndexes
    {
        /** Check that all values in array are in a certain range.
         *
         *  @param array is an array of index values
         *  @param n is the size of array
         *  @param size specifies the range in which array values must fit
         *  @return true if \f$ 0 \le array[i] < size \forall i = 0, ..., n-1\f$
         */

        typedef bool ( *FuncType )( const IndexType array[], const IndexType n, const IndexType size );
        static const char* getId() { return "validIndexes"; }
    };

    /** @brief Trait for register kernel function sum that sums elements of an array
     *
     *  @tparam ValueType specifies the value type used in the sum reduction.
     */
    template <typename ValueType>
    struct sum
    {
        /** @brief Sum n contiguously stored values.
         *
         *  @param[in] array is an array of values
         *  @param[in] n is the size of array
         *  @return sum of all values in array
         */
        typedef ValueType ( *FuncType ) ( const ValueType array[], const IndexType n );
        static const char* getId () { return "sum"; }
    };

    template <typename ValueType>
    struct maxval
    {
        /** @brief Find maximal value of n contiguously stored values.
         *
         *  @param[in] array is an array of values
         *  @param[in] n is the size of array
         *  @return maximum of all values in array
         */

        typedef ValueType ( *FuncType ) ( const ValueType array[], const IndexType n );
        static const char* getId() { return "maxval"; }
    };

    template <typename ValueType>
    struct absMaxVal
    {
        /** @brief Find absolute maximal value of n contiguously stored values. */

        typedef ValueType ( *FuncType ) ( const ValueType array[], const IndexType n );
        static const char* getId() { return "absMaxVal"; }
    };

    template <typename ValueType>
    struct absMaxDiffVal
    {
        /** @brief Building absolute maximum of element-wise difference of vector elements.
         *
         *  @param array1i[in] first array
         *  @param array2i[in] second array
         *  @param n           size of array1 and array2
         *  @returns           max( abs( array1[i] - array2[i] ) ), \f$ 0 \le i < n \f$
         *
         *  Function is helpful to compute maximum norm for vectors and matrices
         */

        typedef ValueType ( *FuncType ) ( const ValueType array1[], const ValueType array2[], const IndexType n );
        static const char* getId() { return "absMaxDiffVal"; }
    };

    template <typename ValueType>
    struct isSorted
    {
        /** @brief Predicate that tests whether a sequene is sorted.
         *
         *  @param[in] array values to be checked
         *  @param[in] n number of values to check
         *  @param[in] ascending if true check for ascending order, otherwise for descending
         */

        typedef bool ( *FuncType ) ( const ValueType array[], const IndexType n, bool ascending );
        static const char* getId() { return "isSorted"; }
    };

    /** @brief Structure with functiooń pointer type defintions for setter methods.
     *
     *  @tparam ValueType specifies the value type used in the set operations.
     */

    template<typename ValueType>
    struct setVal
    {
        /** Set all elements of a contiguous array with a value. */

        typedef void ( *FuncType ) ( ValueType array[], const IndexType n, const ValueType val );
        static const char* getId() { return "setVal"; }
    };

    template<typename ValueType>
    struct setOrder
    {
        /** Set all elements of a contiguous array with its order number 0, 1, 2, ... */

        typedef void ( *FuncType ) ( ValueType array[], const IndexType n );
        static const char* getId() { return "setOrder"; }
    };

    template<typename ValueType>
    struct getValue
    {
        typedef ValueType ( *FuncType ) ( const ValueType* array, const IndexType i );
        static const char* getId() { return "getValue"; }
    };

    template<typename ValueType1, typename ValueType2>
    struct set
    {
        /** Set out[i] = in[i],  0 <= i < n */

        typedef void ( *FuncType ) ( ValueType1 out[], const ValueType2 in[], const IndexType n );
        static const char* getId() { return "set"; }
    };

    template<typename ValueType1, typename ValueType2>
    struct setScale
    {
        /** @brief scaled array assignment, out = in * value
         *
         *  Set out[i] = scale * in[i],  0 <= i < n
         *
         *  @param[in,out]  outValues  is the output array
         *  @param[in]      scaleValue scaling factor
         *  @param[in,out]  inValues   is the array with entries to scale
         *  @param[in]      n          is the number of entries
         */
        typedef void ( *FuncType ) ( ValueType1 outValues[],
                        const ValueType1 scaleValue,
                        const ValueType2 inValues[],
                        const IndexType n );
 
         static const char* getId() { return "setScale"; } 
    };

    template<typename ValueType1, typename ValueType2>
    struct setGather
    {
        /** Set out[i] = in[ indexes[i] ],  \f$0 \le i < n\f$ */

        typedef void ( *FuncType ) ( ValueType1 out[],
                        const ValueType2 in[],
                        const IndexType indexes[],
                        const IndexType n );

        static const char* getId() { return "setGather"; } 
    };

    template<typename ValueType1, typename ValueType2>
    struct setScatter
    {
        /** Set out[ indexes[i] ] = in [i] */

        typedef void ( *FuncType ) ( ValueType1 out[],
                        const IndexType indexes[],
                        const ValueType2 in[],
                        const IndexType n );

        static const char* getId() { return "setScatter"; }
    };

    template<typename ValueType>
    struct invert
    {
        /** @brief Set array[i] = 1.0 / array[i],  0 <= i < n
         *
         *  @param[in,out] array is the array to invert
         *  @param         n     is the number of entries to invert
         */

        typedef void ( *FuncType ) ( ValueType array[], const IndexType n );

        static const char* getId() { return "invert"; }
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

    /** Constructor initializes all function pointers with nullPtr */

    UtilsInterface ();
};

/** @brief Interface for utility functions to be used in CSR storage.
 *
 *  This interface contains function pointer type definitions for all used routines
 *  and tables with actual values for the functions.
 */

struct CSRUtilsInterface
{
    /** Structure defining function types for operations on CSR data
     *
     *  @tparam ValueType is the value type of the matrix element, e.g. float, double
     */

    template<typename ValueType>
    struct Operations
    {
        /** This metod sorts the elemnts of a row by increasing column indexes.
         *
         *  @param[in,out] csrJA, csrValues  the CSR matrix data and their column indexes
         *  @param[in]     csrIA             row offsets
         *  @param[in]     numRows           number of rows
         *  @param[in]     diagonalFlag      if true first entry of each row will be the diagonal element if available
         *
         *  Note: This routine does not force the diagonal property, only if each diagonal element is already available
         */
        typedef void (*sortRowElements)(
            IndexType csrJA[],
            ValueType csrValues[],
            const IndexType csrIA[],
            const IndexType numRows,
            const bool diagonalFlag );
    };

    LAMA_INTERFACE_DEFINE_T( Operations, sortRowElements )

    /** Structure with type definitions for solver routines */

template    <typename ValueType>
    struct Solver
    {
        /** Method to compute one iteration step in Jacobi method
         *
         *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
         *
         */
        typedef void ( *jacobi ) ( ValueType solution[],
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const ValueType csrValues[],
                        const ValueType oldSolution[],
                        const ValueType rhs[],
                        const ValueType omega,
                        const IndexType numRows,
                        SyncToken* syncToken );

        /** Method to compute one iteration step in Jacobi method
         *
         *  solution -= omega * ( B(halo) * oldSolution) * dinv
         *
         */
        typedef void ( *jacobiHalo ) ( ValueType solution[],
                        const IndexType localIA[],
                        const ValueType localValues[],
                        const IndexType haloIA[],
                        const IndexType haloJA[],
                        const ValueType haloValues[],
                        const IndexType haloRowIndexes[],
                        const ValueType oldSolution[],
                        const ValueType omega,
                        const IndexType numNonEmptyRows );

        /** Method to compute one iteration step in Jacobi method
         *
         *  solution -= omega * ( B(halo) * oldSolution) * dinv
         *
         *  @since 1.1.0
         */
        typedef void ( *jacobiHaloWithDiag ) ( ValueType solution[],
                        const ValueType localDiagValues[],
                        const IndexType haloIA[],
                        const IndexType haloJA[],
                        const ValueType haloValues[],
                        const IndexType haloRowIndexes[],
                        const ValueType oldSolution[],
                        const ValueType omega,
                        const IndexType numNonEmptyRows );
    };

    LAMA_INTERFACE_DEFINE_T( Solver, jacobi )
    LAMA_INTERFACE_DEFINE_T( Solver, jacobiHalo )
    LAMA_INTERFACE_DEFINE_T( Solver, jacobiHaloWithDiag )

    /** Structure with type definitions for offset routines. */

    struct Offsets
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

        typedef IndexType ( *sizes2offsets ) ( IndexType array[], const IndexType n );

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

        typedef void ( *offsets2sizes ) ( IndexType sizes[], const IndexType offsets[], const IndexType n );

        /** Check for a legal offset array.
         *
         *  @param[in] array is the array with offsets
         *  @param[in] n array has n + 1 entries
         *  @param[in] total must be last entry in array
         *  @return    true if the array is legal
         *
         *  Means: array[0] <= array[1] <= array[2] <= ... <= array[n] (= total)
         */

        typedef bool ( *validOffsets ) ( const IndexType array[], const IndexType n, const IndexType total );

        /** This method computes the row sizes for result matrix C of matrix add A + B
         *
         *  @param[out] cSizes array of length numRows, will contain number of entries
         *  @param[in]  numRows number of rows for matrix A and B
         *  @param[in]  numColumns number of columns for matrix A and B
         *  @param[in]  diagonalProperty if true, diagonal elements will count in any case
         *  @param[in]  aIA, aJA are the index arrays of matrix A
         *  @param[in]  bIA, bJA are the index arrays of matrix B
         *
         *  Note: filling the result matrix must use the same flag for diagonalProperty
         *        otherwise the row sizes/offsets will not match
         */

        typedef IndexType ( *matrixAddSizes ) ( IndexType cIa[], const IndexType numRows,
                        const IndexType numColumns, bool diagonalProperty,
                        const IndexType aIA[], const IndexType aJA[],
                        const IndexType bIA[], const IndexType bJA[] );

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

        typedef IndexType ( *matrixMultiplySizes ) ( IndexType cSizes[],
                        const IndexType m,
                        const IndexType n,
                        const IndexType k,
                        bool diagonalProperty,
                        const IndexType aIA[], const IndexType aJA[],
                        const IndexType bIA[], const IndexType bJA[] );

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

        typedef void ( *matrixMultiplyJA ) ( IndexType cJA[], const IndexType cIA[],
                        const IndexType numRows, const IndexType numColumns,
                        bool diagonalProperty,
                        const IndexType aIA[], const IndexType aJA[],
                        const IndexType bIA[], const IndexType bJA[] );

        /** This method checks whether the CSR structure data has the diagonal property.
         *
         *  @param[in] numDiagonals  number of first rows for which diagonal property is checked.
         *  @param[in] csrIA         offset array for the rows
         *  @param[in] csrJA         column indexes
         *  @return                  true if diagonal property is given
         *
         *  The diagonal property is given if the first column index in the row is same as the row index.
         */
        typedef bool ( *hasDiagonalProperty ) ( const IndexType numDiagonals,
                        const IndexType csrIA[],
                        const IndexType csrJA[] );

    };

    LAMA_INTERFACE_DEFINE( Offsets, sizes2offsets )
    LAMA_INTERFACE_DEFINE( Offsets, offsets2sizes )
    LAMA_INTERFACE_DEFINE( Offsets, validOffsets )
    LAMA_INTERFACE_DEFINE( Offsets, matrixAddSizes )
    LAMA_INTERFACE_DEFINE( Offsets, matrixMultiplySizes )
    LAMA_INTERFACE_DEFINE( Offsets, matrixMultiplyJA )
    LAMA_INTERFACE_DEFINE( Offsets, hasDiagonalProperty )

    /** Define structure that contains type definitions for the function pointers.
     *
     *  @tparam ValueType specifies the value type used in the arrays.
     *
     *  The structure is needed as type definition templates are unsupported in C++.
     */

    template<typename ValueType>
    struct Transpose
    {
        /** Function pointer for CSR to CSC conversion routine.
         *
         *  @param[out] cscIA, cscJA, cscValues is CSC output data
         *  @param[in]  csrIA, csrJA, csrValues is CSR input data
         *  @param numRows x numColumns is shape of input CSR matrix
         *
         *  Arrays must be big enough for the corresponing sizes.
         */

        typedef void( *convertCSR2CSC ) ( IndexType cscIA[],
                        IndexType cscJA[],
                        ValueType cscValues[],
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const ValueType csrValues[],
                        IndexType numRows, IndexType numColumns,
                        IndexType numValues );
    };

    LAMA_INTERFACE_DEFINE_T( Transpose, convertCSR2CSC )

    /** Define structure for multiplication routines.  */

    template<typename ValueType1, typename ValueType2>
    struct Scale
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

        typedef void ( *scaleRows ) ( ValueType1 csrValues[],
                        const IndexType csrIA[],
                        const IndexType numRows,
                        const ValueType2 values[] );
    };

    LAMA_INTERFACE_DEFINE_TT( Scale, scaleRows )

    /** Structure with type definitions for reduction routines. */

    template <typename ValueType>
    struct Reductions
    {
        /** Get the maximal element-wise difference for two CSR matrices.
         *
         *  @param[in] numRows is number of rows for both matrices
         *  @param[in] sortedRows if true column indexes in JA arrays are sorted
         *  @param[in] csrIA1, csrJA1, csrValues1 is storage data of first matrix
         *  @param[in] csrIA2, csrJA2, csrValues2 is storage data of second matrix
         *  @returns maximal value of absolute difference between two matrix elements
         */

        typedef ValueType ( *absMaxDiffVal ) ( IndexType numRows, bool sortedRows,
                        const IndexType csrIA1[], const IndexType csrJA1[], const ValueType csrValues1[],
                        const IndexType csrIA2[], const IndexType csrJA2[], const ValueType csrValues2[] );
    };

    LAMA_INTERFACE_DEFINE_T( Reductions, absMaxDiffVal )

    /** Define structure for multiplication routines.
     *
     *  @tparam ValueType specifies the value type used in mutliplications.
     */

    template<typename ValueType>
    struct Mult
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
         *  @param syncToken optional, if available starts asynchronous computation
         *
         *  The number of columns is not really needed to implement the operation but
         *  might be helpful to choose between different implementations.
         */

        typedef void ( *normalGEMV ) ( ValueType result[],
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
                        SyncToken* syncToken );

        typedef void ( *normalGEVM ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const ValueType beta,
                        const ValueType y[],
                        const IndexType numRows,
                        const IndexType numColumns,
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const ValueType csrValues[],
                        SyncToken* syncToken );

        /** result = alpha * CSR-Matrix * x, CSR matrix has only some non-zero rows
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param numNonZeroRows is size of rowIndexes
         *  @param rowIndexes are indexes of non-empty rows in matrix
         *  @param csrIA, csrJA, csrValues are arrays of CSR storage
         *  @param syncToken optional, if available starts asynchronous computation
         *
         *  Note: this routine does not provide the term 'beta * y' as it would require
         *        to run over the full result vector
         */

        typedef void ( *sparseGEMV ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const IndexType numNonZeroRows,
                        const IndexType rowIndexes[],
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const ValueType csrValues[],
                        SyncToken* syncToken );

        typedef void ( *sparseGEVM ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const IndexType numColumns,
                        const IndexType numNonZeroRows,
                        const IndexType rowIndexes[],
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const ValueType csrValues[],
                        SyncToken* syncToken );

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
         *   @param[in,out] syncToken might be used for asynchronous execution
         */

        typedef void ( *gemm ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const ValueType beta,
                        const ValueType y[],
                        const IndexType m,
                        const IndexType n,
                        const IndexType p,
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const ValueType csrValues[],
                        SyncToken* syncToken );

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

        typedef void ( *matrixAdd ) ( IndexType cJA[],
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

        typedef void ( *matrixMultiply ) ( const IndexType cIa[],
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
    };

    LAMA_INTERFACE_DEFINE_T( Mult, normalGEMV )
    LAMA_INTERFACE_DEFINE_T( Mult, normalGEVM )
    LAMA_INTERFACE_DEFINE_T( Mult, gemm )
    LAMA_INTERFACE_DEFINE_T( Mult, sparseGEMV )
    LAMA_INTERFACE_DEFINE_T( Mult, sparseGEVM )
    LAMA_INTERFACE_DEFINE_T( Mult, matrixAdd )
    LAMA_INTERFACE_DEFINE_T( Mult, matrixMultiply )

    CSRUtilsInterface ();
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

/** Traits for utility functions to be used in Dense storage.   */

struct DenseUtilsInterface
{
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

    /** Function pointer type definitions for modification of dense storage. */

    template<typename DenseValueType>
    struct scaleValue
    {
        /** Scale all elements of the dense matrix with a value */

        typedef void ( *FuncType ) ( DenseValueType denseValues[],
                        const IndexType numRows,
                        const IndexType numColumns,
                        const DenseValueType val );

        static const char* getId() { return "Dense.scaleValue"; }
    };

    template<typename DenseValueType>
    struct setDiagonalValue
    {
        /** Set diagonal elements with one and the same value. */

        typedef void ( *FuncType ) ( DenseValueType denseValues[],
                        const IndexType numRows,
                        const IndexType numColumns,
                        const DenseValueType val );

        static const char* getId() { return "Dense.setDiagonalValue"; }
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
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param numColumns is the number of columns, size of x
         *  @param numValues, array with the dense matrix values
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *FuncType ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const ValueType beta,
                        const ValueType y[],
                        const IndexType numRows,
                        const IndexType numColumns,
                        const ValueType denseValues[],
                        SyncToken* syncToken );

        static const char* getId() { return "Dense.normalGEMV"; }
    };
};

/** Interface for utility functions to be used in ELL storage.
 *
 *  This interface contains function pointer type definitions for all used routines
 *  and tables with actual values for the functions.
 */

struct ELLUtilsInterface
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
         *  @param[in]  syncToken optional, NULL synchronous execution, otherwise asynchronous
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
            const ValueType omega,
            SyncToken* syncToken );

        static const char* getId() { return "ELL.jacobi"; }
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
            const ValueType omega,
            SyncToken* syncToken );

        static const char* getId() { return "ELL.jacobiHalo"; }
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

        static const char* getId() { return "ELL.fillELLValues"; }
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

        static const char* getId() { return "ELL.getCSRValues"; }
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

        static const char* getId() { return "ELL.setCSRValues"; }
    };

    template<typename ValueType>
    struct compressIA
    {
        /** Compresses the given IA array using the values array and epsilon
         *
         * @param[in]  IA that should be compressed
         * @param[in]  related JA array
         * @param[in]  related values array
         * @param[in]  number of rows
         * @param[in]  epsilon
         * @param[out] new created IA
         */
        typedef void ( *FuncType ) ( 
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ValueType ellValues[],
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const ValueType eps,
            IndexType newIA[] );

        static const char* getId() { return "ELL.compressIA"; }
    };

    template<typename ValueType>
    struct compressValues
    {
        /** Compresses the given JA and values array using epsilon
         *
         * @param[in]  IA that should be compressed
         * @param[in]  related JA array
         * @param[in]  related values array
         * @param[in]  number of rows
         * @param[in]  epsilon
         * @param[out] new created JA
         * @param[out] new created values
         */
        typedef void ( *FuncType ) ( 
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ValueType ellValues[],
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const ValueType eps,
            const IndexType newNumValuesPerRow,
            IndexType newJA[],
            ValueType newValues[] );

        static const char* getId() { return "ELL.compressValues"; }
    };

    template<typename ValueType, typename OtherValueType>
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
            OtherValueType row[],
            const IndexType i,
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numValuesPerRow,
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ValueType ellValues[] );

        static const char* getId() { return "ELL.compressValues"; }
    };

    template<typename ValueType, typename OtherValueType>
    struct getValue
    {
        /** Returns one element of the matrix
         *
         *  @param[in] i is the row of the returned element
         *  @param[in] j is the column of the returned element
         *  @param[in] numRows is the number of rows of the matrix
         *  @param[in] ellSizes is the ELL sizes array
         *  @param[in] ellJA is the ELL ja array
         *  @param[in] ellValues is the ELL values array
         */

        typedef OtherValueType ( *FuncType ) ( 
            const IndexType i,
            const IndexType j,
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ValueType ellValues[] );
        
        static const char* getId() { return "ELL.getValue"; }
    };

    struct countNonEmptyRowsBySizes
    {
        typedef IndexType ( *FuncType ) ( 
            const IndexType ellSizes[],
            const IndexType numRows );

        static const char* getId() { return "ELL.countNonEmptyRowsBySizes"; }
    };

    struct setNonEmptyRowsBySizes
    {
        typedef void ( *FuncType ) ( 
            IndexType rowIndexes[],
            const IndexType numNonEmptyRows,
            const IndexType ellSizes[],
            const IndexType numRows );

        static const char* getId() { return "ELL.setNonEmptyRowsBySizes"; }
    };

    struct hasDiagonalProperty
    {
        typedef bool ( *FuncType ) ( 
            const IndexType numDiagonals,
            const IndexType ellJA[] );

        static const char* getId() { return "ELL.hasDiagonalProperty"; }
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

        static const char* getId() { return "ELL.check"; }
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
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *FuncType ) ( 
            ValueType result[],
            const ValueType alpha,
            const ValueType x[],
            const ValueType beta,
            const ValueType y[],
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ValueType ellValues[],
            SyncToken* syncToken );

        static const char* getId() { return "ELL.normalGEMV"; }
    };

    template<typename ValueType>
    struct sparseGEMV
    {
        /** result = alpha * ELL-Matrix * x, CSR matrix has only some non-zero rows
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param numNonZeroRows is size of rowIndexes
         *  @param rowIndexes are indexes of non-empty rows in matrix
         *  @param ellIA, ellJA, csrValues are arrays of ELL storage
         *  @param syncToken optional, if available starts asynchronous computation
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
            SyncToken* syncToken );

        static const char* getId() { return "ELL.sparseGEMV"; }
    };

    template<typename ValueType>
    struct normalGEVM
    {
        /** Implementation for ELLUtilsInterface::Mult::normalGEVM  */

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
            SyncToken* syncToken );

        static const char* getId() { return "ELL.normalGEVM"; }
    };

    template<typename ValueType>
    struct sparseGEVM
    {
        /** Implementation for ELLUtilsInterface::Mult::sparseGEVM  */

        typedef void ( *FuncType ) ( 
            ValueType result[],
            const ValueType alpha,
            const ValueType x[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numValuesPerRow,
            const IndexType numNonZeroRows,
            const IndexType rowIndexes[],
            const IndexType ellSizes[],
            const IndexType ellJA[],
            const ValueType ellValues[],
            SyncToken* syncToken );

        static const char* getId() { return "ELL.sparseGEVM"; }
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

        static const char* getId() { return "ELL.absMaxVal"; }
    };

    template<typename ValueType, typename OtherValueType>
    struct scaleValue
    {
        typedef void ( *FuncType ) ( 
            const IndexType numRows,
            const IndexType numValuesPerRow,
            const IndexType ellSizes[],
            ValueType ellValues[],
            const OtherValueType values[] );

        static const char* getId() { return "ELL.scaleValue"; }
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

        static const char* getId() { return "ELL.matrixMultiplySizes"; }
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

        static const char* getId() { return "ELL.matrixAddSizes"; }
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

        static const char* getId() { return "ELL.matrixAdd"; }
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

        static const char* getId() { return "ELL.matrixMultiply"; }
    };
};

/** Interface for utility functions to be used in JDS storage.
 *
 *  This interface contains function pointer type definitions for all used routines
 *  and tables with actual values for the functions.
 */

struct JDSUtilsInterface
{
    /** Structure with type definitions for solver routines */
    template<typename ValueType>
    struct Solver
    {
        /** Method to compute one iteration step in Jacobi method
         *
         *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
         *
         */
        typedef void (*jacobi)(
            ValueType* const solution,
            const IndexType numRows,
            const IndexType jdsPerm[],
            const IndexType jdsIlg[],
            const IndexType jdsNumDiagonals,
            const IndexType jdsDlg[],
            const IndexType jdsJA[],
            const ValueType jdsValues[],
            const ValueType oldSolution[],
            const ValueType rhs[],
            const ValueType omega,
            SyncToken* syncToken );

        /** Method to compute one iteration step in Jacobi method with halo.  */

        typedef void (*jacobiHalo)(
            ValueType solution[],
            const IndexType numRows,
            const ValueType invDiagonal[],
            const IndexType numDiagonals,
            const IndexType jdsHaloPerm[],
            const IndexType jdsHaloIlg[],
            const IndexType jdsHaloDlg[],
            const IndexType jdsHaloJA[],
            const ValueType jdsHaloValues[],
            const ValueType oldSolution[],
            const ValueType omega,
            SyncToken* syncToken );
    };

    LAMA_INTERFACE_DEFINE_T( Solver, jacobi )LAMA_INTERFACE_DEFINE_T( Solver, jacobiHalo )

    struct Sort
    {
        /** Stable sorting of values in array in descending order.
         *
         *  @param[in,out] array are the values to be sorted
         *  @param[in,out] perm, where perm[i] has the value of the original position
         *  @param[in]    n is the number of values to be sorted
         *
         *  \code
         *           array =   1  4   1  8  5  7
         *           perm  =   0  1   2  3  4  5
         +
         *           array =   8  7   5  4  1  1
         *           perm  =   3  5   4  1  0  2
         *  \endcode
         */

        typedef void ( *sortRows ) ( IndexType array[],
                        IndexType perm[],
                        const IndexType n );

        /** Compute the inverse permutation for a given permutation.
         *
         *  inversePerm [ perm [i] ] == i , 0 <= i < n
         *
         *  @param[out] inversePerm, size = n, will contain the inverse permutation
         *  @param[in] perm, size = n, is input permuation of 0, ..., n-1
         *  @param[in] n specifies the size of perm and inversePerm
         *
         *  /code
         *       perm      2  5  1  4  6  3  0
         *     inperm      6  2  0  5  3  1  4
         *  /endcode
         */

        typedef void ( *setInversePerm ) ( IndexType inversePerm[],
                        const IndexType perm[],
                        const IndexType n );

        /** Compute dlg array from ilg array.
         *
         *  @param[out] dlg is the array with sizes of the columns
         *  @param[in]  numDiagonals is size of dlg
         *  @param[in]  ilg is the array with sizes of the rows
         *  @param[in]  numRows is the number of rows, size of ilg
         *
         *  The values in ilg must be descending. The same will be true
         *  for the output array dlg.
         *
         *  /code
         *       ilg       4  3  2  2  1   dlg
         *       5         x  x  x  x  x
         *       4         x  x  x  x
         *       2         x  x
         *       1         x
         *  /endcode
         */

        typedef IndexType ( *ilg2dlg ) ( IndexType dlg[], const IndexType numDiagonals,
                        const IndexType ilg[], const IndexType numRows );

    };

    LAMA_INTERFACE_DEFINE( Sort, sortRows )
    LAMA_INTERFACE_DEFINE( Sort, setInversePerm )
    LAMA_INTERFACE_DEFINE( Sort, ilg2dlg )

    template<typename JDSValueType, typename CSRValueType>
    struct Conversions
    {
        /** Conversion of JDS storage data to CSR data
         *
         *  @param[out]  csrJA is array with column indexes
         *  @param[out]  csrValues is array with non-zero values
         *  @param[in]   csrIA is offset array (must be computed before)
         *  @param[in]   numRows number of rows in matrix
         *  @param[in]   jdsPerm with jdsPerm[ii] is original index of row i
         *  @param[in]   jdsILG with size of entries in row i
         *  @param[in]   jdsDLG distances of columns
         *  @param[in]   jdsJA column indexes
         *  @param[in]   jdsValues matrix values
         */

        typedef void ( *getCSRValues ) ( IndexType csrJA[],
                        CSRValueType csrValues[],
                        const IndexType csrIA[],
                        const IndexType numRows,
                        const IndexType jdsPerm[],
                        const IndexType jdsILG[],
                        const IndexType jdsDLG[],
                        const IndexType jdsJA[],
                        const JDSValueType jdsValues[] );

        /** Conversion of CSR storage data to JDS data
         *
         *  @param[out]  jdsJA column indexes
         *  @param[out]  jdsValues matrix values
         *  @param[in]   numRows number of rows in matrix
         *  @param[in]   jdsPerm with jdsPerm[ii] is original index of row i
         *  @param[in]   jdsILG with size of entries in row i
         *  @param[in]   numDiagonals size of array jdsDLG
         *  @param[in]   jdsDLG distances of columns
         *  @param[in]   csrIA is offset array (must be computed before)
         *  @param[in]   csrJA is array with column indexes
         *  @param[in]   csrValues is array with non-zero values
         */

        typedef void( *setCSRValues ) ( IndexType jdsJA[],
                        JDSValueType jdsValues[],
                        const IndexType numRows,
                        const IndexType jdsPerm[],
                        const IndexType jdsILG[],
                        const IndexType numDiagonals,
                        const IndexType jdsDLG[],
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const CSRValueType csrValues[] );
    };

    // Define tables ( indexed by template value types) for all methods

    LAMA_INTERFACE_DEFINE_TT( Conversions, setCSRValues )
    LAMA_INTERFACE_DEFINE_TT( Conversions, getCSRValues )

    template<typename ValueType>
    struct Mult
    {
        /** result = alpha * CSR-Matrix * x + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param numDiagonals are number of diagonals, is size of jdsDLG
         *  @param jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues are arrays of JDS storage
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *normalGEMV ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const ValueType beta,
                        const ValueType y[],
                        const IndexType numRows,
                        const IndexType jdsPerm[],
                        const IndexType jdsILG[],
                        const IndexType ndlg,
                        const IndexType jdsDLG[],
                        const IndexType jdsJA[],
                        const ValueType jdsValues[],
                        SyncToken* syncToken );

        typedef void ( *normalGEVM ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const ValueType beta,
                        const ValueType y[],
                        const IndexType numColumns,
                        const IndexType jdsPerm[],
                        const IndexType jdsILG[],
                        const IndexType ndlg,
                        const IndexType jdsDLG[],
                        const IndexType jdsJA[],
                        const ValueType jdsValues[],
                        SyncToken* syncToken );
    };

    LAMA_INTERFACE_DEFINE_T( Mult, normalGEMV )
    LAMA_INTERFACE_DEFINE_T( Mult, normalGEVM )

    template<typename ValueType, typename OtherValueType>
    struct Getter
    {
        typedef void ( *getRow ) ( OtherValueType row[],
                        const IndexType i,
                        const IndexType numColumns,
                        const IndexType numRows,
                        const IndexType perm[],
                        const IndexType ilg[],
                        const IndexType dlg[],
                        const IndexType ja[],
                        const ValueType values[] );

        typedef ValueType (*getValue ) ( const IndexType i,
                        const IndexType j,
                        const IndexType numRows,
                        const IndexType* dlg,
                        const IndexType* ilg,
                        const IndexType* perm,
                        const IndexType* ja,
                        const ValueType* values );

    };

    LAMA_INTERFACE_DEFINE_TT( Getter, getRow )
    LAMA_INTERFACE_DEFINE_TT( Getter, getValue )

    template<typename ValueType, typename OtherValueType>
    struct Scale
    {
        typedef void ( *scaleValue ) ( const IndexType numRows,
                        const IndexType perm[],
                        const IndexType ilg[],
                        const IndexType dlg[],
                        ValueType mValues[],
                        const OtherValueType values[] );
    };

    LAMA_INTERFACE_DEFINE_TT( Scale, scaleValue )

    struct Helper
    {
        typedef bool ( *checkDiagonalProperty ) ( const IndexType numDiagonals,
                        const IndexType numRows,
                        const IndexType numColumns,
                        const IndexType perm[],
                        const IndexType ja[],
                        const IndexType dlg[] );
    };

    LAMA_INTERFACE_DEFINE( Helper, checkDiagonalProperty )

    JDSUtilsInterface ();
};

/** Interface for utility functions to be used in DIA storage.
 *
 *  This interface contains function pointer type definitions for all used routines
 *  and tables with actual values for the functions.
 */

struct DIAUtilsInterface
{
    template<typename ValueType>
    struct Counting
    {
        /** Type definition of function pointer for counting sparse values in DIA storage
         *
         *  @param[out] csrSizes array with number of non-zero entries in each row
         *  @param[in]  diagonalFlag if true count also zero diagonal elements
         *  @param[in] numRows is the number of rows
         *  @param[in] numColumns is the number of columns
         *  @param[in] numDiagonals number of diagonals used in the DIA format
         *  @param[in] diaOffsets diagonal offsets, size is numDiagonals
         *  @param[in] diaValues are stored values of the diagonals
         *  @param[in] eps threshold value when an element should be considered as zero
         *
         *  Note: the diagonals might contain zero entries so the number of non-zero
         *        elements might be less than number of stored elements
         *
         *  - csrSizes must have been allocated with at least numRows entries
         *  - diaOffsets has at least numDiagonals entries
         *  - diaValues has numDiagonals x max(numRows, numColumns) entries
         */

        typedef void (*getCSRSizes)(
            IndexType csrSizes[],
            bool diagonalFlag,
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numDiagonals,
            const IndexType diaOffsets[],
            const ValueType diaValues[],
            const ValueType eps );

    };

    LAMA_INTERFACE_DEFINE_T( Counting, getCSRSizes )

template    <typename DIAValueType, typename CSRValueType>
    struct Conversions
    {
        /** Type definition of function pointer for conversion of DIA storage data to CSR data.
         *
         *  @param[out] csrJA will contain the column indexes
         *  @param[out] csrValues will contain the matrix elements
         *  @param[in] csrIA is the array with the offsets (must already be available before)
         *  @param[in] numRows is the number of rows
         *  @param[in] numColumns is the number of columns
         *  @param[in] numDiagonals number of diagonals used in the DIA format
         *  @param[in] diaOffsets diagonal offsets, size is numDiagonals
         *  @param[in] diaValues are stored values of the diagonals
         *  @param[in] eps threshold value when an element should be considered as zero
         *
         *   - csrIA has numRows + 1 entries
         *   - csrJA and csrValues must have at least numValues entries, numValues = csrIA[numRows]
         */

        typedef void ( *getCSRValues ) ( IndexType csrJA[],
                        CSRValueType csrValues[],
                        const IndexType csrIA[],
                        const bool diagonalFlag,
                        const IndexType numRows,
                        const IndexType numColumns,
                        const IndexType numDiagonals,
                        const IndexType diaOffsets[],
                        const DIAValueType diaValues[],
                        const DIAValueType eps );
    };

    LAMA_INTERFACE_DEFINE_TT( Conversions, getCSRValues )

    template<typename ValueType>
    struct Mult
    {
        /** result = alpha * CSR-Matrix * x + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param numValues is the number of diagonals in DIA storage
         *  @param diaOffsets, diaValues are arrays of DIA storage
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *normalGEMV ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const ValueType beta,
                        const ValueType y[],
                        const IndexType numRows,
                        const IndexType numColumns,
                        const IndexType numDiagonals,
                        const IndexType diaOffsets[],
                        const ValueType diaValues[],
                        SyncToken* syncToken );

        /** result = alpha * x * CSR-Matrix + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param numValues is the number of diagonals in DIA storage
         *  @param diaOffsets, diaValues are arrays of DIA storage
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *normalGEVM ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const ValueType beta,
                        const ValueType y[],
                        const IndexType numRows,
                        const IndexType numColumns,
                        const IndexType numDiagonals,
                        const IndexType diaOffsets[],
                        const ValueType diaValues[],
                        SyncToken* syncToken );
    };

    LAMA_INTERFACE_DEFINE_T( Mult, normalGEMV )
    LAMA_INTERFACE_DEFINE_T( Mult, normalGEVM )

    /** Structure with type definitions for solver routines */

    template<typename ValueType>
    struct Solver
    {
        /** Method to compute one iteration step in Jacobi method
         *
         *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
         *
         */
        typedef void ( *jacobi ) ( ValueType* const solution,
                        const IndexType numColumns,
                        const IndexType numDiagonals,
                        const IndexType diaOffset[],
                        const ValueType diaValues[],
                        const ValueType oldSolution[],
                        const ValueType rhs[],
                        const ValueType omega,
                        const IndexType numRows,
                        SyncToken* syncToken );
    };

    LAMA_INTERFACE_DEFINE_T( Solver, jacobi )

    /** Structure with type definitions for reduction routines */

    template<typename ValueType>
    struct Reductions
    {
        /** This method returns the maximal absolute value of an ELLPACK matrix. */

        typedef ValueType ( *absMaxVal ) ( const IndexType numRows,
                        const IndexType numColumns,
                        const IndexType numDiagonals,
                        const IndexType diaOffsets[],
                        const ValueType diaValues[]
        );
    };

    LAMA_INTERFACE_DEFINE_T( Reductions, absMaxVal )

    /** Constructur of interface sets all function pointers to nullPtr. */

    DIAUtilsInterface ();
};

/** Interface for utility functions to be used in COO storage.
 *
 *  This interface contains function pointer type definitions for all used routines
 *  and tables with actual values for the functions.
 */

struct COOUtilsInterface
{
    struct Counting
    {
        /** Helper routine for conversion of COO format to CSR format to get sparse row sizes.
         *
         *  @param[out] csrSizes array with number of non-zero entries in each row
         *  @param[in] numRows number of rows
         *  @param[in] numValues number of non-zero values
         *  @param[in] array with row indexes of COO storage (size is numValues)
         */

        typedef void (*getCSRSizes)(
            IndexType csrSizes[],
            const IndexType numRows,
            const IndexType numValues,
            const IndexType cooIA[] );

        /** Routine for conversion of CSR offset array to COO ia array
         *
         *  @param[out] cooIA is the array with all row indexes
         *  @param[in] numValues number of non-zero values, size of cooIA
         *  @param[in] csrIA is the CSR row offset array, size is numRows+1
         *  @param[in] numRows number of rows
         *  @param[in] numDiagonals is number of diagonals
         *
         *  The diagonal values will be stored at the beginning of the array cooIA.
         */

        typedef void (*offsets2ia)(
            IndexType cooIA[],
            const IndexType numValues,
            const IndexType csrIA[],
            const IndexType numRows,
            const IndexType numDiagonals );

    };

    LAMA_INTERFACE_DEFINE( Counting, getCSRSizes )LAMA_INTERFACE_DEFINE( Counting, offsets2ia )

    template<typename COOValueType, typename CSRValueType>
    struct Conversions
    {
        /** Helper routine for conversion COO to CSR
         *
         *  @param[out] csrJA will contain the column indexes
         *  @param[out] csrValues will contain the matrix elements
         *  @param[in] csrIA is the array with the offsets (must already be available before)
         *
         *   - csrIA has numRows + 1 entries
         *   - csrJA and csrValues must have at least numValues entries, numValues = csrIA[numRows]
         *
         *  Note: this routine preserves the diagonal property of the COO format
         */

        typedef void ( *getCSRValues )( IndexType csrJA[],
                        CSRValueType csrValues[],
                        IndexType csrIA[],
                        const IndexType numRow,
                        const IndexType numValues,
                        const IndexType cooIA[],
                        const IndexType cooJA[],
                        const COOValueType cooValues[] );

        /** Conversion of CSR data (ja, values) to COO data.
         *
         *  @param[out] cooValues is new COO data with diagonal elements at first
         *  @param[in] csrValues is given CSR data with diagonal elements first in each row
         *  @param[in] numValues is size of arrays cooValues and csrValues
         *  @param[in] csrIA is CSR offset array
         *  @param[in] numRows is number of rows
         *  @param[in] numDiagonals is number of diagonal elements to be stored at the beginning
         *  @param[in] csrDiagonalProperty is true if CSR data has diagonal property
         *
         *  Note: Diagonal elements must be first in each row for CSR data, no resort done for this
         *  Note: For numDiagonals == 0, this routine can be replaced with Utils::set.
         */

        typedef void ( *setCSRData ) ( COOValueType cooValues[],
                        const CSRValueType csrValues[],
                        const IndexType numValues,
                        const IndexType csrIA[],
                        const IndexType numRows,
                        const IndexType numDiagonals );
    };

    LAMA_INTERFACE_DEFINE_TT( Conversions, getCSRValues )
    LAMA_INTERFACE_DEFINE_TT( Conversions, setCSRData )

    template<typename ValueType>
    struct Mult
    {
        /** result = alpha * CSR-Matrix * x + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param cooIA, cooJA, cooValues are arrays of COO storage
         *  @param numValues is the size of the coo arrays
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *normalGEMV ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const ValueType beta,
                        const ValueType y[],
                        const IndexType numRows,
                        const IndexType nnz,
                        const IndexType cooIA[],
                        const IndexType cooJA[],
                        const ValueType cooValues[],
                        SyncToken* syncToken );

        typedef void ( *normalGEVM ) ( ValueType result[],
                        const ValueType alpha,
                        const ValueType x[],
                        const ValueType beta,
                        const ValueType y[],
                        const IndexType numRows,
                        const IndexType nnz,
                        const IndexType cooIA[],
                        const IndexType cooJA[],
                        const ValueType cooValues[],
                        SyncToken* syncToken );
    };

    LAMA_INTERFACE_DEFINE_T( Mult, normalGEMV )
    LAMA_INTERFACE_DEFINE_T( Mult, normalGEVM )

    /** Structure with type definitions for solver routines */
    template<typename ValueType>
    struct Solver
    {
        /** Method to compute one iteration step in Jacobi method
         *
         *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
         *
         */
        typedef void ( *jacobi ) ( ValueType* const solution,
                        const IndexType cooNumValues,
                        const IndexType cooIA[],
                        const IndexType cooJA[],
                        const ValueType cooValues[],
                        const ValueType oldSolution[],
                        const ValueType rhs[],
                        const ValueType omega,
                        const IndexType numRows,
                        SyncToken* syncToken );

        /** Method to compute one iteration step in Jacobi method
         *
         *  solution -= omega * ( B(halo) * oldSolution) * dinv
         *
         */
        typedef void ( *jacobiHalo ) ( ValueType solution[],
                        const ValueType diaValues[],
                        const IndexType haloIA[],
                        const IndexType haloJA[],
                        const ValueType haloValues[],
                        const IndexType haloRowIndexes[],
                        const ValueType oldSolution[],
                        const ValueType omega,
                        const IndexType numNonEmptyRows );
    };

    LAMA_INTERFACE_DEFINE_T( Solver, jacobi )
    LAMA_INTERFACE_DEFINE_T( Solver, jacobiHalo )

    /** Default constructor initializes all function pointers with NULL. */

    COOUtilsInterface ();
};

} /* end namespace lama */

} /* end namespace scai */
