/**
 * @file CUDACSRUtils.hpp
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
 * @brief General conversion routines for CSR sparse matrices.
 * @author Thomas Brandes
 * @date 03.07.2012
 * @since 1.0.0
 */
#ifndef LAMA_CUDA_CSR_UTILS_HPP_
#define LAMA_CUDA_CSR_UTILS_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/SyncToken.hpp>
#include <lama/LAMATypes.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

class LAMA_DLL_IMPORTEXPORT CUDACSRUtils
{
public:
    /** Helper routine for conversion CSR to CSR format.  */

    static IndexType sizes2offsets( IndexType sizes[], const IndexType numRows );

    static void offsets2sizes( IndexType sizes[], const IndexType offsets[], const IndexType n );

    static bool hasDiagonalProperty( const IndexType numDiagonals, const IndexType csrIA[], const IndexType csrJA[] );

    /** Matrix transpose for CSR matrices on CUDA device. */

    template<typename ValueType>
    static void convertCSR2CSC(
        IndexType cscIA[],
        IndexType cscJA[],
        ValueType cscValues[],
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        int numRows,
        int numColumns,
        int numValues );

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

    /** Implementation for CSRUtilsInterface::Mult::normalGEVM  */

    template<typename ValueType>
    static void normalGEVM(
        ValueType result[],
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

    /** Implementation for CSRUtilsInterface::Mult::sparseGEVM  */

    template<typename ValueType>
    static void sparseGEVM(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType numColumns,
        const IndexType numNonZeroRows,
        const IndexType rowIndexes[],
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        SyncToken* syncToken );

    /** Implementation for CSRUtilsInterface::Solver::jacobi  */

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

    /** Implementation for CSRUtilsInterface::Solver::jacobiHalo  */

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

    /** Implementation for CSRUtilsInterface::Solver::jacobiHaloWithDiag
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

    /** Routine that registers all routines of this class at the LAMA interface. */

    static void setInterface( struct CSRUtilsInterface& CSRUtils );

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    static    bool initialized; //!< static initialization used for registration

    static unsigned int lastHashTableSize;// local variable to handhover hash table size for multiply

    static bool registerInterface();//!< registration
};

/* --------------------------------------------------------------------------- */

}
// namespace lama

#endif  //  LAMA_CUDA_CSR_UTILS_HPP_
