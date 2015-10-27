/**
 * @file JDSKernelTrait.hpp
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
 * @brief Struct with traits for all JDS storage methods provided as kernels.
 * @author Thomas Brandes
 * @date 23.10.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

namespace scai
{

// forward declaration

namespace tasking
{
    class SyncToken;
}

namespace lama
{

/** Kernel traits for functions to be used in JDS storage.  */

struct JDSKernelTrait
{
    template<typename ValueType>
    struct jacobi
    {
        /** Method to compute one iteration step in Jacobi method
         *
         *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
         *
         */
        typedef void ( *FuncType )(
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
            tasking::SyncToken* syncToken );

        static const char* getId() { return "JDS.jacobi"; }
    };

    /** Structure with type definitions for solver routines */
    template<typename ValueType>
    struct jacobiHalo
    {
        /** Method to compute one iteration step in Jacobi method with halo.  */

        typedef void ( *FuncType )(
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
            tasking::SyncToken* syncToken );

        static const char* getId() { return "JDS.jacobiHalo"; }
    };

    struct sortRows
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

        typedef void ( *FuncType ) ( IndexType array[],
                        IndexType perm[],
                        const IndexType n );

        static const char* getId() { return "JDS.sortRows"; }
    };

    struct setInversePerm
    {
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

        typedef void ( *FuncType ) ( IndexType inversePerm[],
                        const IndexType perm[],
                        const IndexType n );

        static const char* getId() { return "JDS.setInversePerm"; }
    };

    struct ilg2dlg
    {

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

        typedef IndexType ( *FuncType ) ( IndexType dlg[], const IndexType numDiagonals,
                        const IndexType ilg[], const IndexType numRows );

        static const char* getId() { return "JDS.ilg2dlg"; }
    };

    template<typename JDSValueType, typename CSRValueType>
    struct getCSRValues
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

        typedef void ( *FuncType ) ( IndexType csrJA[],
                        CSRValueType csrValues[],
                        const IndexType csrIA[],
                        const IndexType numRows,
                        const IndexType jdsPerm[],
                        const IndexType jdsILG[],
                        const IndexType jdsDLG[],
                        const IndexType jdsJA[],
                        const JDSValueType jdsValues[] );

        static const char* getId() { return "JDS.getCSRValues"; }
    };

    template<typename JDSValueType, typename CSRValueType>
    struct setCSRValues
    {

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

        typedef void( *FuncType ) ( IndexType jdsJA[],
                        JDSValueType jdsValues[],
                        const IndexType numRows,
                        const IndexType jdsPerm[],
                        const IndexType jdsILG[],
                        const IndexType numDiagonals,
                        const IndexType jdsDLG[],
                        const IndexType csrIA[],
                        const IndexType csrJA[],
                        const CSRValueType csrValues[] );

        static const char* getId() { return "JDS.setCSRValues"; }
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
         *  @param numDiagonals are number of diagonals, is size of jdsDLG
         *  @param jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues are arrays of JDS storage
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *FuncType ) ( ValueType result[],
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
                        tasking::SyncToken* syncToken );

        static const char* getId() { return "JDS.normalGEMV"; }
    };

    template<typename ValueType>
    struct normalGEVM
    {
        typedef void ( *FuncType ) ( ValueType result[],
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
                        tasking::SyncToken* syncToken );

        static const char* getId() { return "JDS.normalGEVM"; }
    };

    template<typename ValueType, typename OtherValueType>
    struct getRow
    {
        typedef void ( *FuncType ) ( OtherValueType row[],
                        const IndexType i,
                        const IndexType numColumns,
                        const IndexType numRows,
                        const IndexType perm[],
                        const IndexType ilg[],
                        const IndexType dlg[],
                        const IndexType ja[],
                        const ValueType values[] );

        static const char* getId() { return "JDS.getRow"; }
    };

    template<typename ValueType>
    struct getValue
    {
        typedef ValueType ( *FuncType ) ( const IndexType i,
                        const IndexType j,
                        const IndexType numRows,
                        const IndexType* dlg,
                        const IndexType* ilg,
                        const IndexType* perm,
                        const IndexType* ja,
                        const ValueType* values );

        static const char* getId() { return "JDS.getValue"; }
    };

    template<typename ValueType, typename OtherValueType>
    struct scaleValue
    {
        typedef void ( *FuncType ) ( const IndexType numRows,
                        const IndexType perm[],
                        const IndexType ilg[],
                        const IndexType dlg[],
                        ValueType mValues[],
                        const OtherValueType values[] );

        static const char* getId() { return "JDS.scaleValue"; }
    };

    struct checkDiagonalProperty
    {
        typedef bool ( *FuncType ) ( const IndexType numDiagonals,
                        const IndexType numRows,
                        const IndexType numColumns,
                        const IndexType perm[],
                        const IndexType ja[],
                        const IndexType dlg[] );

        static const char* getId() { return "JDS.checkDiagonalProperty"; }
    };
};

} /* end namespace lama */

} /* end namespace scai */
