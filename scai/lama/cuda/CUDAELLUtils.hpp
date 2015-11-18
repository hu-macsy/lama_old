/**
 * @file CUDAELLUtils.hpp
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
 * @brief General conversion routines for ELL sparse matrices.
 * @author Thomas Brandes
 * @date 03.07.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai library
#include <scai/common/SCAITypes.hpp>

#include <scai/logging.hpp>

namespace scai
{

namespace lama
{

/** This class provides routines to converse ELL storage data to CSR storage data and vice versa.
 *
 *  All routines work on already allocated data and utilize CUDA for their parallelization.
 */

class COMMON_DLL_IMPORTEXPORT CUDAELLUtils
{
public:

    /** Addressing function for the arrays ia and ja: column-wise */

    static inline IndexType ellindex( const IndexType i, const IndexType jj, const IndexType numRows )
    {
        return jj * numRows + i;
    }

    /** This method computes the total number of non-zero rows by the size array  */

    static IndexType countNonEmptyRowsBySizes( const IndexType sizes[], const IndexType numRows );

    /** check diagonal property. ELL format with diagonal property: diagonal is just the first column in mValues */

    static bool hasDiagonalProperty( const IndexType numDiagonals, const IndexType ellJA[] );

    /** Build a vector of indexes for non-empty rows. */

    static void setNonEmptyRowsBySizes(
        IndexType rowIndexes[],
        const IndexType numNonEmptyRows,
        const IndexType sizes[],
        const IndexType numRows );

    static void check(
        const IndexType mNumRows,
        const IndexType mNumValuesPerRow,
        const IndexType mNumColumns,
        const IndexType *ia,
        const IndexType *ja,
        const char* msg );

    /** Returns one row of the matrix */

    template<typename ValueType,typename OtherValueType>
    static void getRow(
        OtherValueType row[],
        const IndexType i,
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    /** Returns one value of the matrix */

    template<typename ValueType>
    static ValueType getValue(
        const IndexType i,
        const IndexType j,
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    /** Scales matrix using an vector */

    template<typename ValueType,typename OtherValueType>
    static void scaleValue(
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        ValueType ellValues[],
        const OtherValueType values[] );

    /** Implementation for ELLKernelTrait::Conversions::getCSRValues */

    template<typename ELLValueType,typename CSRValueType>
    static void getCSRValues(
        IndexType csrJA[],
        CSRValueType csrValues[],
        const IndexType csrIA[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ELLValueType ellValues[] );

    template<typename ValueType>
    static void fillELLValues(
        IndexType ellJA[],
        ValueType ellValues[],
        const IndexType ellSizes[],
        const IndexType numRows,
        const IndexType numValuesPerRow );

    /** Helper routine for conversion CSR to ELL format.  */

    template<typename ELLValueType,typename CSRValueType>
    static void setCSRValues(
        IndexType ellJA[],
        ELLValueType ellValues[],
        const IndexType ellSizes[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const CSRValueType csrValues[] );

    /** Implementation for ELLKernelTrait::normalGEMV on CUDA devices. */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellIA[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    /** Implementation for ELLKernelTrait::normalGEVM  */

    template<typename ValueType>
    static void normalGEVM(
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
        const ValueType ellValues[] );

    /** Implementation for ELLKernelTrait::sparseGEMV  */

    template<typename ValueType>
    static void sparseGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType numNonZeroRows,
        const IndexType rowIndexes[],
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    /** Implementation for ELLKernelTrait::sparseGEVM  */

    template<typename ValueType>
    static void sparseGEVM(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numNonZerosPerRow,
        const IndexType numNonZeroRows,
        const IndexType rowIndexes[],
        const IndexType ellIA[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    /** Implementation for ELLKernelTrait::jacobi  */

    template<typename ValueType>
    static void jacobi(
        ValueType solution[],
        const IndexType numRows,
        const IndexType ellNumValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType ellValues[],
        const ValueType oldSolution[],
        const ValueType rhs[],
        const ValueType omega );

    /** Implementation for ELLKernelTrait::Solver::jacobiHalo  */

    template<typename ValueType>
    static void jacobiHalo(
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

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Routine that registers all methods at the kernel registry. */

    static void registerKernels( bool deleteFlag );

    /** Constructor for registration. */

    CUDAELLUtils();

    /** Destructor for unregistration. */

    ~CUDAELLUtils();

    /** Static variable for registration at static initialization. */

    static CUDAELLUtils guard;
};

/* --------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
