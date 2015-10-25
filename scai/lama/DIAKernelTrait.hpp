/**
 * @file DIAKernelTrait.hpp
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
 * @brief Struct with traits for all LAMA kernels needed by DIA storage 
 * @author Thomas Brandes
 * @date 03.04.2013
 * @since 1.0.0
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

/** Traits for LAMA kernels used in DIA storage.  */

struct DIAKernelTrait
{
    template<typename ValueType>
    struct getCSRSizes
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

        typedef void ( *FuncType )(
            IndexType csrSizes[],
            bool diagonalFlag,
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numDiagonals,
            const IndexType diaOffsets[],
            const ValueType diaValues[],
            const ValueType eps );

        static const char* getId() { return "DIA.getCSRSizes"; }
    };

    template<typename DIAValueType, typename CSRValueType>
    struct getCSRValues
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

        typedef void ( *FuncType ) ( 
            IndexType csrJA[],
            CSRValueType csrValues[],
            const IndexType csrIA[],
            const bool diagonalFlag,
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numDiagonals,
            const IndexType diaOffsets[],
            const DIAValueType diaValues[],
            const DIAValueType eps );

        static const char* getId() { return "DIA.getCSRValues"; }
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
         *  @param numValues is the number of diagonals in DIA storage
         *  @param diaOffsets, diaValues are arrays of DIA storage
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *FuncType ) ( 
            ValueType result[],
            const ValueType alpha,
            const ValueType x[],
            const ValueType beta,
            const ValueType y[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numDiagonals,
            const IndexType diaOffsets[],
            const ValueType diaValues[],
            tasking::SyncToken* syncToken );

        static const char* getId() { return "DIA.normalGEMV"; }
    };

    template<typename ValueType>
    struct normalGEVM
    {
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

        typedef void ( *FuncType ) ( 
            ValueType result[],
            const ValueType alpha,
            const ValueType x[],
            const ValueType beta,
            const ValueType y[],
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numDiagonals,
            const IndexType diaOffsets[],
            const ValueType diaValues[],
            tasking::SyncToken* syncToken );

        static const char* getId() { return "DIA.normalGEVM"; }
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
            ValueType* solution,
            const IndexType numColumns,
            const IndexType numDiagonals,
            const IndexType diaOffset[],
            const ValueType diaValues[],
            const ValueType oldSolution[],
            const ValueType rhs[],
            const ValueType omega,
            const IndexType numRows,
            tasking::SyncToken* syncToken );

        static const char* getId() { return "DIA.jacobi"; }
    };

    /** Structure with type definitions for reduction routines */

    template<typename ValueType>
    struct absMaxVal
    {
        /** This method returns the maximal absolute value of an ELLPACK matrix. */

        typedef ValueType (  *FuncType ) ( 
            const IndexType numRows,
            const IndexType numColumns,
            const IndexType numDiagonals,
            const IndexType diaOffsets[],
            const ValueType diaValues[]
        );

        static const char* getId() { return "DIA.absMaxVal"; }
    };
};

} /* end namespace lama */

} /* end namespace scai */
