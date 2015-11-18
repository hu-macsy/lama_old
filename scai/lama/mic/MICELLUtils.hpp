/**
 * @file MICELLUtils.hpp
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
 * @since 1.1.0
 */

#pragma once

// for dll import
#include <scai/common/config.hpp>

// others
#include <scai/common/SCAITypes.hpp>

// other scai projects
#include <scai/logging.hpp>

namespace scai
{

namespace lama
{

/** This class provides routines to converse ELL storage data to CSR storage data and vice versa.
 *
 *  All routines work on already allocated data and utilize MIC for their parallelization.
 */

class COMMON_DLL_IMPORTEXPORT MICELLUtils
{

public:

    /** This method computes the total number of non-zero rows by the size array  */

    static IndexType countNonEmptyRowsBySizes( const IndexType sizes[], const IndexType numRows );

    /** Build a vector of indexes for non-empty rows. */

    static void setNonEmptyRowsBySizes(
        IndexType rowIndexes[],
        const IndexType numNonEmptyRows,
        const IndexType sizes[],
        const IndexType numRows );

    /** Addressing function for the arrays ellJA[numRows*numValuesPerRow] and ellValues: column-major order */

    __attribute__( ( target( mic ) ))
    static inline IndexType ellindex(
        const IndexType i,
        const IndexType jj,
        const IndexType numRows,
        const IndexType /* numValuesPerRow */)
    {
        return jj * numRows + i; // column major-order
        // return i * numValuesPerRow + jj;    // row major-order
    }

    /** Returns the maximal absolute value of the ELL storage. */

    template<typename ValueType>
    static ValueType absMaxVal(
        const IndexType mNumRows,
        const IndexType mNumValuesPerRow,
        const IndexType ia[],
        const ValueType values[] );

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

    /** check diagonal property */

    static bool hasDiagonalProperty( const IndexType numDiagonals, const IndexType csrJA[] );

    /** Implementation for ELLKernelTrait::scaleValue */

    template<typename ValueType,typename OtherValueType>
    static void scaleValue(
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        ValueType ellValues[],
        const OtherValueType values[] );

    /** Implementation for ELLKernelTrait::compressIA */

    template<typename ValueType>
    static void compressIA(
        const IndexType IA[],
        const IndexType JA[],
        const ValueType values[],
        const IndexType numRows,
        const ValueType eps,
        IndexType newIA[] );

    /** Implementation for ELLKernelTrait::compressValues */

    template<typename ValueType>
    static void compressValues(
        const IndexType IA[],
        const IndexType JA[],
        const ValueType values[],
        const IndexType numRows,
        const ValueType eps,
        IndexType newJA[],
        ValueType newValues[] );

    /** Implementation for ELLKernelTrait::getCSRValues */

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

    /** Implementation for ELLKernelTrait::setCSRValues */

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

    /** Implementation for ELLKernelTrait::computeIA */

    template<typename ValueType>
    static void computeIA(
        const IndexType aIA[],
        const IndexType aJA[],
        const IndexType aNumRows,
        const IndexType bIA[],
        const IndexType bJA[],
        const IndexType bNumRows,
        IndexType cIA[] );

    /** Implementation for ELLKernelTrait::computeValues */

    template<typename ValueType>
    static void computeValues(
        const IndexType aIA[],
        const IndexType aJA[],
        const ValueType aValues[],
        const IndexType aNumRows,
        const IndexType bIA[],
        const IndexType bJA[],
        const ValueType bValues[],
        const IndexType bNumRows,
        const ValueType alpha,
        const IndexType cIA[],
        IndexType cJA[],
        ValueType cValues[] );

    /** Implementation for ELLKernelTrait::addComputeIA */

    template<typename ValueType>
    static void addComputeIA(
        const IndexType aIA[],
        const IndexType aJA[],
        const IndexType aNumRows,
        const IndexType bIA[],
        const IndexType bJA[],
        const IndexType bNumRows,
        IndexType cIA[] );

    /** Implementation for ELLKernelTrait::addComputeValues */

    template<typename ValueType>
    static void addComputeValues(
        const IndexType aIA[],
        const IndexType aJA[],
        const ValueType aValues[],
        const IndexType aNumRows,
        const IndexType bIA[],
        const IndexType bJA[],
        const ValueType bValues[],
        const IndexType bNumRows,
        const ValueType beta,
        const IndexType cIA[],
        IndexType cJA[],
        ValueType cValues[] );

    /** Implementation for ELLKernelTrait::jacobi */

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

    /** Implementation for ELLKernelTrait::jacobiHalo */

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

    /** Implementation for ELLKernelTrait::normalGEMV  */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const ValueType beta,
        const ValueType y[],
        const IndexType numRows,
        const IndexType numNonZerosPerRow,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[] );

    /** Implementation for ELLKernelTrait::sparseGEMV  */

    template<typename ValueType>
    static void sparseGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType numRows,
        const IndexType numNonZerosPerRow,
        const IndexType numNonZeroRows,
        const IndexType rowIndexes[],
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[] );

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    /** Routine that registers all methods at the kernel registry. */

    static void registerKernels( bool deleteFlag );

    /** Helper class for (un) registration of kernel routines at static initialization. */

    class RegisterGuard
    {
    public:
        RegisterGuard();
        ~RegisterGuard();
    };

    static RegisterGuard guard;  // registration of kernels @ static initialization
};

/* --------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
