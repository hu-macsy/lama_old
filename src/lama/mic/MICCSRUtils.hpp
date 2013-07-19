/**
 * @file MICCSRUtils.hpp
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
 * @brief Implementation of CSR utilities with MIC
 * @author Thomas Brandes
 * @date 02.07.2013
 * @since 1.1.0
 */
#ifndef LAMA_MIC_CSR_UTILS_HPP_
#define LAMA_MIC_CSR_UTILS_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMATypes.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

class SyncToken;   // forward declaration 

/** This class provides routines on compressed sparse row data
 */

class LAMA_DLL_IMPORTEXPORT MICCSRUtils
{
public:

    /** This method computes the total number of non-zero rows by the offset array
     *
     */

    static IndexType countNonEmptyRowsByOffsets( const IndexType offsets[], const IndexType numRows );

    /** Build a vector of indexes for non-empty rows. */

    static void setNonEmptyRowsByOffsets(
        IndexType rowIndexes[],
        const IndexType numNonEmptyRows,
        const IndexType offsets[],
        const IndexType numRows );

    /** This function build an offset array from a counter array.
     *
     *  @param[in,out] array contains counter values and later the offsets
     *  @param[in]     numValues is the number of values, array must contain one additional value
     *  @returns       the total number of values
     *
     *  \code
     *    array  :    3    7   8   4   2
     *    array  :    0   10  15  12  16    -> returns 18
     *  \endcode
     *
     *  CSR sparse representation of matrices store the sum of all values at an additional
     *  position at the end of the array.
     */

    static IndexType scan( IndexType array[], const IndexType numValues );

    /** Implementation for CSRUtilsInterface::Offsets::sizes2offsets */

    static IndexType sizes2offsets( IndexType sizes[], const IndexType numRows );

    /** Implementation for CSRUtilsInterface::Offsets::offsets2sizes */

    static void offsets2sizes( IndexType sizes[], const IndexType offsets[], const IndexType n );

    /** offset2sizes for indexed rows */

    static void offsets2sizesGather(
        IndexType sizes[],
        const IndexType offsets[],
        const IndexType rowIndexes[],
        const IndexType numRows );

    /** Implementation for CSRUtilsInterface::Offsets::validOffsets  */

    static bool validOffsets( const IndexType array[], const IndexType n, const IndexType total );

    /** This function checks if csr data has diagonal property.
     *
     *  @param[in] numDiagonals is min( numRows, numColumns )
     *  @param[in] csrIA is the csr offset array
     *  @param[in] csrJA is array with column indexes
     */

    static bool hasDiagonalProperty( const IndexType numDiagonals, const IndexType csrIA[], const IndexType csrJA[] );

    /** This function sorts the column indexes of each row in ascending order.
     *
     *  If the diagonal flag is set, the first entry in the row will be the
     *  diagonal element ( if available ).
     */

    template<typename ValueType>
    static void sortRowElements(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType numRows,
        const bool diagonalFlag );

    /** Implementation for CSRUtilsInterface::Transpose::convertCSR2CSC  */

    template<typename ValueType>
    static void convertCSR2CSC(
        IndexType cscIA[],
        IndexType cscJA[],
        ValueType cscValues[],
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        IndexType numRows,
        IndexType numColumns,
        IndexType numValues );

    /** Implementation for CSRUtilsInterface::Mult::scaleRows  */

    template<typename ValueType1,typename ValueType2>
    static void scaleRows(
        ValueType1 csrValues[],
        const IndexType csrIA[],
        const IndexType numRows,
        const ValueType2 values[] );

    /** Implementation for CSRUtilsInterface::Mult::normalGEMV  */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[],
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

    /** Implementation for CSRUtilsInterface::Mult::sparseGEMV  */

    template<typename ValueType>
    static void sparseGEMV(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType numNonZeroRows,
        const IndexType rowIndexes[],
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        SyncToken* syncToken );

    /** Implementation for CSRUtilsInterface::Mult::gemm  */

    template<typename ValueType>
    static void gemm(
        ValueType result[],
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

    /** Implementation for CSRUtilsInterface::Jacobi::jacobi(Async/Halo) */

    template<typename ValueType>
    static void jacobi(
        ValueType* const solution,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        const ValueType rhs[],
        const ValueType oldSolution[],
        const ValueType omega,
        const IndexType numRows,
        SyncToken* syncToken );

    template<typename ValueType>
    static void jacobiHalo(
        ValueType solution[],
        const IndexType localIA[],
        const ValueType localValues[],
        const IndexType haloIA[],
        const IndexType haloJA[],
        const ValueType haloValues[],
        const IndexType haloRowIndexes[],
        const ValueType oldSolution[],
        const ValueType omega,
        const IndexType numNonEmptyRows );

    /** Implementation for CSRUtilsInterface::Jacobi::jacobiHaloWithDiag
     *  @since 1.1.0
     */

    template<typename ValueType>
    static void jacobiHaloWithDiag(
        ValueType solution[],
        const ValueType localDiagValues[],
        const IndexType haloIA[],
        const IndexType haloJA[],
        const ValueType haloValues[],
        const IndexType haloRowIndexes[],
        const ValueType oldSolution[],
        const ValueType omega,
        const IndexType numNonEmptyRows );

    /** Implementation for CSRUtilsInterface::Offsets::matrixAddSizes  */

    static IndexType matrixAddSizes(
        IndexType cSizes[],
        const IndexType numRows,
        const IndexType numColumns,
        bool diagonalProperty,
        const IndexType aIA[],
        const IndexType aJA[],
        const IndexType bIA[],
        const IndexType bJA[] );

    /** Implementation for CSRUtilsInterface::Offsets::matrixMultiplySizes  */

    static IndexType matrixMultiplySizes(
        IndexType cSizes[],
        const IndexType m,
        const IndexType n,
        const IndexType k,
        bool diagonalProperty,
        const IndexType aIA[],
        const IndexType aJA[],
        const IndexType bIA[],
        const IndexType bJA[] );

    /** Implementation for CSRUtilsInterface::Mult::matrixAdd */

    template<typename ValueType>
    static void matrixAdd(
        IndexType cJA[],
        ValueType cValues[],
        const IndexType cIA[],
        const IndexType numRows,
        const IndexType numColumns,
        bool diagonalProperty,
        const ValueType alpha,
        const IndexType aIA[],
        const IndexType aJA[],
        const ValueType aValues[],
        const ValueType beta,
        const IndexType bIA[],
        const IndexType bJA[],
        const ValueType bValues[] );

    /** Implementation for CSRUtilsInterface::Mult::matrixMultiply */

    template<typename ValueType>
    static void matrixMultiply(
        const IndexType cIa[],
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

    /** Implementation for CSRUtilsInterface::Reductions::absMaxDiffVal */

    template<typename ValueType>
    static ValueType absMaxDiffVal(
        IndexType numRows,
        bool sortedRows,
        const IndexType csrIA1[],
        const IndexType csrJA1[],
        const ValueType csrValues1[],
        const IndexType csrIA2[],
        const IndexType csrJA2[],
        const ValueType csrValues2[] );

    /** Routine that registers all routines of this class at the LAMA interface. */

    static void setInterface( struct CSRUtilsInterface& CSRUtils );

protected:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

private:

    static bool initialized;

    static bool registerInterface();

    static IndexType scanSerial( IndexType array[], const IndexType numValues );

    static IndexType scanParallel( PartitionId numThreads, IndexType array[], const IndexType numValues );

    template<typename ValueType>
    __attribute__( ( target( mic ) ) )
    static ValueType absMaxDiffRowSorted(
        const IndexType n1, const IndexType csrJA1[], const ValueType csrValues1[],
        const IndexType n2, const IndexType csrJA2[], const ValueType csrValues2[] );

    template<typename ValueType>
    __attribute__( ( target( mic ) ) )
    static ValueType absMaxDiffRowUnsorted(
        const IndexType n1, const IndexType csrJA1[], const ValueType csrValues1[],
        const IndexType n2, const IndexType csrJA2[], const ValueType csrValues2[] );
};

} // namespace lama

#endif //  LAMA_MIC_CSR_UTILS_HPP_
