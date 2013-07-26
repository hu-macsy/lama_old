/**
 * @file OpenMPELLUtils.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
#ifndef LAMA_OPENMP_ELL_UTILS_HPP_
#define LAMA_OPENMP_ELL_UTILS_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMATypes.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

class SyncToken;  // forward declaration 

/** This class provides routines to converse ELL storage data to CSR storage data and vice versa.
 *
 *  All routines work on already allocated data and utilize OpenMP for their parallelization.
 */

class LAMA_DLL_IMPORTEXPORT OpenMPELLUtils
{

public:

    /** Routine that registers all routines of this class at the LAMA interface. */

    static void setInterface( struct ELLUtilsInterface& ELLUtils );

private:

    /** This method computes the total number of non-zero rows by the size array  */

    static IndexType countNonEmptyRowsBySizes( const IndexType sizes[], const IndexType numRows );

    /** Build a vector of indexes for non-empty rows. */

    static void setNonEmptyRowsBySizes(
        IndexType rowIndexes[],
        const IndexType numNonEmptyRows,
        const IndexType sizes[],
        const IndexType numRows );

    /** Addressing function for the arrays ellJA[numRows*numValuesPerRow] and ellValues: column-major order */

    static inline IndexType ellindex( const IndexType i, const IndexType jj, const IndexType numRows, const IndexType /* numValuesPerRow */ )
    {
        return jj * numRows + i;            // column major-order
        // return i * numValuesPerRow + jj;    // row major-order
    }

    /** Returns the maximal absolute value of the ELL storage. */

    template<typename ValueType>
    static ValueType absMaxVal(
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const ValueType ellValues[] );

    static void check(
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType numColumns,
        const IndexType ellSizes[],
        const IndexType ellJA[],
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

    template<typename ValueType,typename OtherValueType>
    static OtherValueType getValue(
        const IndexType i,
        const IndexType j,
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        const IndexType ellJA[],
        const ValueType ellValues[] );

    /** check diagonal property */

    static bool hasDiagonalProperty( const IndexType numDiagonals, const IndexType csrJA[] );

    template<typename ValueType,typename OtherValueType>
    static void scaleValue(
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const IndexType ellSizes[],
        ValueType ellValues[],
        const OtherValueType values[] );

    /** Implementation for ELLUtilsInterface::Conversions::compressIA */

    template<typename ValueType>
    static void compressIA(
        const IndexType ellIA[],
        const IndexType ellJA[],
        const ValueType ellValues[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const ValueType eps,
        IndexType newIA[] );

    /** Implementation for ELLUtilsInterface::Conversions::compressValues */

    template<typename ValueType>
    static void compressValues(
        const IndexType ellIA[],
        const IndexType ellJA[],
        const ValueType ellValues[],
        const IndexType numRows,
        const IndexType numValuesPerRow,
        const ValueType eps,
        const IndexType newNumValuesPerRow,
        IndexType newJA[],
        ValueType newValues[] );

    /** Implementation for ELLUtilsInterface::Conversions::getCSRValues */

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

    /** Implementation for ELLUtilsInterface::Conversions::setCSRValues */

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

    /** Implementation for ELLUtilsInterface::MatrixExpBuild::matrixMultiplySizes */

    static void matrixMultiplySizes(
        IndexType cSizes[],
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

    /** Implementation for ELLUtilsInterface::MatrixExp::matrixMultiply */

    template<typename ValueType>
    static void matrixMultiply(
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

    /** Implementation for ELLUtilsInterface::MatrixExpBuild::matrixAddSizes */

    static void matrixAddSizes(
        IndexType csizes[],
        const IndexType m,
        const IndexType n,
        const bool diagonalProperty,
        const IndexType aSizes[],
        const IndexType aJA[],
        const IndexType aNumValuesPerRow,
        const IndexType bSizes[],
        const IndexType bJA[],
        const IndexType bNumValuesPerRow );

    /** Implementation for ELLUtilsInterface::MatrixExp::matrixAdd */

    template<typename ValueType>
    static void matrixAdd(
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

    /** Implementation for ELLUtilsInterface::Solver::jacobi */

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
        const ValueType omega,
        SyncToken* syncToken );

    /** Implementation for ELLUtilsInterface::Solver::jacobiHalo */

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
        const ValueType omega,
        SyncToken* syncToken );

    /** Implementation for ELLUtilsInterface::Mult::normalGEMV  */

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
        const ValueType csrValues[],
        SyncToken* syncToken );

    /** Implementation for ELLUtilsInterface::Mult::sparseGEMV  */

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
        const ValueType csrValues[],
        SyncToken* syncToken );

    /** Implementation for CSRUtilsInterface::Mult:normalGEVM */

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
        const ValueType ellValues[],
        SyncToken* syncToken );

    /** Implementation for CSRUtilsInterface::Mult::sparseGEVM  */

    template<typename ValueType>
    static void sparseGEVM(
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

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    static bool initialized;

    static bool registerInterface();
};

/* --------------------------------------------------------------------------- */

} // namespace lama

#endif  //  LAMA_ELL_STORAGE_UTILS_HPP_
