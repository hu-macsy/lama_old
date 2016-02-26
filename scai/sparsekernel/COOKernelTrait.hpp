/**
 * @file COOKernelTrait.hpp
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
 * @brief Struct with traits for all COO storage methods provided as kernels.
 * @author Thomas Brandes
 * @date 22.10.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

namespace scai
{

namespace sparsekernel
{

/** Kernel traits for functions to be used in COO storage. */

struct COOKernelTrait
{
    struct hasDiagonalProperty
    {
        /** Routine checks for diagonal property, first n entries are the diagonal elements.
         *
         *  @param[in] cooIA row indexes 
         *  @param[in] cooJA column indexes
         *  @param[in] n number of diagonal elements
         *  @return true if first n entries stand for the diagonal elements
         */

        typedef bool ( *FuncType )(
            const IndexType cooIA[],
            const IndexType cooJA[],
            const IndexType n ); 

        static const char* getId() { return "COO.hasDiagonalProperty"; }
    };

    struct getCSRSizes
    {
        /** Helper routine for conversion of COO format to CSR format to get sparse row sizes.
         *
         *  @param[out] csrSizes array with number of non-zero entries in each row
         *  @param[in] numRows number of rows
         *  @param[in] numValues number of non-zero values
         *  @param[in] array with row indexes of COO storage (size is numValues)
         */

        typedef void ( *FuncType )(
            IndexType csrSizes[],
            const IndexType numRows,
            const IndexType numValues,
            const IndexType cooIA[] );

        static const char* getId() { return "COO.getCSRSizes"; }
    };

    struct offsets2ia
    {
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

        typedef void ( *FuncType )(
            IndexType cooIA[],
            const IndexType numValues,
            const IndexType csrIA[],
            const IndexType numRows,
            const IndexType numDiagonals );

        static const char* getId() { return "COO.offsets2ia"; }
    };

    template<typename COOValueType, typename OtherValueType>
    struct scaleRows
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
            COOValueType cooValues[],
            const OtherValueType rowValues[],
            const IndexType cooIA[],
            const IndexType numValues );

        static const char* getId() { return "COO.scaleRows"; }
    };

    template<typename COOValueType, typename CSRValueType>
    struct getCSRValues
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

        typedef void ( *FuncType )( 
            IndexType csrJA[],
            CSRValueType csrValues[],
            IndexType csrIA[],
            const IndexType numRow,
            const IndexType numValues,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const COOValueType cooValues[] );

        static const char* getId() { return "COO.getCSRValues"; }
    };

    template<typename COOValueType, typename CSRValueType>
    struct setCSRData
    {

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

        typedef void ( *FuncType ) ( 
            COOValueType cooValues[],
            const CSRValueType csrValues[],
            const IndexType numValues,
            const IndexType csrIA[],
            const IndexType numRows,
            const IndexType numDiagonals );

        static const char* getId() { return "COO.setCSRData"; }
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
            const IndexType nnz,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const ValueType cooValues[] );

        static const char* getId() { return "COO.normalGEMV"; }
    };

    template<typename ValueType>
    struct normalGEVM
    {
        typedef void ( *FuncType ) ( 
            ValueType result[],
            const ValueType alpha,
            const ValueType x[],
            const ValueType beta,
            const ValueType y[],
            const IndexType numRows,
            const IndexType nnz,
            const IndexType cooIA[],
            const IndexType cooJA[],
            const ValueType cooValues[] );

        static const char* getId() { return "COO.normalGEVM"; }
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

        static const char* getId() { return "COO.jacobi"; }
    };

    template<typename ValueType>
    struct jacobiHalo
    {
        /** Method to compute one iteration step in Jacobi method
         *
         *  solution -= omega * ( B(halo) * oldSolution) * dinv
         *
         */
        typedef void ( *FuncType ) ( 
            ValueType solution[],
            const ValueType diaValues[],
            const IndexType haloIA[],
            const IndexType haloJA[],
            const ValueType haloValues[],
            const IndexType haloRowIndexes[],
            const ValueType oldSolution[],
            const ValueType omega,
            const IndexType numNonEmptyRows );

        static const char* getId() { return "COO.jacobiHalo"; }
    };
};

} /* end namespace sparsekernel */

} /* end namespace scai */
