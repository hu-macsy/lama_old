/**
 * @file CSRUtils.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Utility functions for CSR data
 * @author Thomas Brandes
 * @date 14.02.2018
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/BinaryOp.hpp>

namespace scai
{

namespace sparsekernel
{

class COMMON_DLL_IMPORTEXPORT CSRUtils
{
public:

    /** 
     *  @brief check that csrIA is a valid row offset array for the CSR format.
     *
     *  @param[in] csrIA array to check
     *  @param[in] numValues is the total number of non-zeros 
     *  @param[in] prefLoc specfies the context where to execute it
     */
    static bool validOffsets(
        const hmemo::HArray<IndexType>& csrIA,
        const IndexType numValues,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief compute the offset array from size array
     *
     *  @param[in]  sizes contains number of entries for each row
     *  @param[out] offsets is the offset array, size will be sizes.size()
     *  @param[in] prefLoc specfies the context where to execute it
     * 
     *  \code
     *     HArray<IndexType> sizes( { 2, 5, 7 } );
     *     CSRUtils::sizes2offsets( sizes, sizes, prefLoc );
     *     -> sizes == { 0, 2, 7, 14 }, returns 14
     *  \endcode
     */
    static IndexType sizes2offsets( 
        hmemo::HArray<IndexType>& offsets, 
        const hmemo::HArray<IndexType>& sizes, 
        hmemo::ContextPtr prefLoc );

    static void offsets2sizes( 
        hmemo::HArray<IndexType>& sizes, 
        const hmemo::HArray<IndexType>& offsets, 
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Get the sizes for a set of rows 
     *
     *  @param[in] csrIA is the offset array of the CSR storage
     *  @param[in] rowIndexes is an array with the row indexes for which sizes are needed
     *  @param[out] sizes contains the sizes of the rows as specified in rowIndexes
     *  @param[in] prefLoc is the context where operation should be executed.
     */
    static void gatherSizes( 
        hmemo::HArray<IndexType>& sizes, 
        const hmemo::HArray<IndexType>& csrIA, 
        const hmemo::HArray<IndexType>& rowIndexes, 
        hmemo::ContextPtr prefLoc );

    /** @brief get the indexes of non-empty rows 
     *
     *  @param[out] rowIndexes contains the indexes of non-zero rows
     *  @param[in]  csrIA the CSR row offset array
     *  @param[in]  threshold builds rowIndexes only if #nonZeroRows/#numRows < threshhold
     *  @param[in]  prefLoc is the context where operation is executed
     *  @returns the number of non-zero rows, will also be the size of rowIndexes if built
     *
     *  If threshhold is 0, the row indexes are never built.
     */
    static IndexType nonEmptyRows( 
        hmemo::HArray<IndexType>& rowIndexes, 
        const hmemo::HArray<IndexType>& csrIA, 
        float threshold,
        hmemo::ContextPtr prefLoc );

    /** @brief This method generates new CSR data where all zero elemens are removed.
     *
     *  @param[in,out] csrIA, csrJA, csrValues is the CSR data that is compressed
     *  @param[in] eps a value is considered to be zero if abs( value ) <= eps
     *  @param[in] prefLoc specficies the context where compression should be done
     */  
    template<typename ValueType>
    static void compress( 
        hmemo::HArray<IndexType>& csrIA, 
        hmemo::HArray<IndexType>& csrJA, 
        hmemo::HArray<ValueType>& csrValues,
        const RealType<ValueType> eps, 
        hmemo::ContextPtr prefLoc );

    /** 
     *  @brief sort column entries in each row of CSR data
     *
     *  This routine sorts the column indexes for each row in ascending order.
     *
     *  @param[in,out] ja array with the column indexes
     *  @param[in,out] values array with the matrix values
     *  @param[in] ia offset array, size is numRows + 1
     *  @param[in] numRows the number of rows, only needed for convenience
     *  @param[in] numColumns the number of columns, only needed for convenience
     *  @param[in] prefLoc specifies the context where the operation should be executed
     */
    template<typename ValueType>
    static void sortRows(
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>& values,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Check if the column entries of CSR data are sorted 
     *
     *  @param[in] ia is the CSR offset array
     *  @param[in] ja are the column indexes
     *  @param[in] numRows needed for convenience, same as ia.size() - 1
     *  @param[in] numColumns needed for convenience
     *  @param[prefLoc] specifies the context where the operation should be executed
     */
    static bool hasSortedRows(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        hmemo::ContextPtr prefLoc );
 
    /**
     *  @brief This routine moves the diagonal entries to the first entry entry for each row
     * 
     *  @param[in,out] ja array with the column indexes
     *  @param[in,out] values array with the matrix values
     *  @param[in] ia offset array, size is numRows + 1
     *  @param[in] numRows the number of rows, only needed for convenience
     *  @param[in] numColumns the number of columns, only needed for convenience
     *  @param[in] prefLoc specifies the context where the operation should be executed
     *  @return the number of rows where diagonal element is first element
     * 
     *  Only the diagonal element is moved to the first entry, all other elements
     *  remain in the original order. This operation might be used for optimization
     *  of CSR storage to access the diagonal (set/get diagonal).
     */
    template<typename ValueType>
    static IndexType shiftDiagonalFirst(
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>& values,
        const IndexType numRows,   
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        hmemo::ContextPtr prefLoc );

    /** Conversion routine of compressed sparse row data to compressed sparse column.
     *
     *  This method does not keep diagonalProperty.
     *  But the data might be sorted in the columns.
     */
    template<typename ValueType>
    static void convertCSR2CSC(
        hmemo::HArray<IndexType>& colIA,
        hmemo::HArray<IndexType>& colJA,
        hmemo::HArray<ValueType>& colValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& rowIA,
        const hmemo::HArray<IndexType>& rowJA,
        const hmemo::HArray<ValueType>& rowValues,
        hmemo::ContextPtr prefLoc );

    /** Matrix multiplication 
     *
     *  This method multiplies two matrices A and B (both CSR format) and 
     *  delivers the result also in CSR format.
     *
     *  @param[out] cIA, cJA, cValues for output storage in CSR format
     *  @param[in] aIA, aJA, aValues for first input storage
     *  @param[in] bIA, bJA, bValues for second input storage
     *  @param[in] m, n  are the sizes of the output matrix
     *  @param[in] k    are the number of columns of a and number of rows of b
     *  @param[in] alpha scaling factor 
     *  @param[in] prefLoc specifies the preferred context for execution.
     */
    template<typename ValueType>
    static void matrixMultiply(
        hmemo::HArray<IndexType>& cIA,
        hmemo::HArray<IndexType>& cJA,
        hmemo::HArray<ValueType>& cValues,
        const ValueType alpha,
        const hmemo::HArray<IndexType>& aIA,
        const hmemo::HArray<IndexType>& aJA,
        const hmemo::HArray<ValueType>& aValues,
        const hmemo::HArray<IndexType>& bIA,
        const hmemo::HArray<IndexType>& bJA,
        const hmemo::HArray<ValueType>& bValues,
        const IndexType m,
        const IndexType n,
        const IndexType k ,
        hmemo::ContextPtr prefLoc );

    template<typename ValueType>
    static void matrixAdd(
        hmemo::HArray<IndexType>& cIA,
        hmemo::HArray<IndexType>& cJA,
        hmemo::HArray<ValueType>& cValues,
        const ValueType alpha,
        const hmemo::HArray<IndexType>& aIA,
        const hmemo::HArray<IndexType>& aJA,
        const hmemo::HArray<ValueType>& aValues,
        const ValueType beta,
        const hmemo::HArray<IndexType>& bIA,
        const hmemo::HArray<IndexType>& bJA,
        const hmemo::HArray<ValueType>& bValues,
        const IndexType m,
        const IndexType n,
        hmemo::ContextPtr prefLoc );

    template<typename ValueType>
    static void binaryOp(
        hmemo::HArray<IndexType>& cIA,
        hmemo::HArray<IndexType>& cJA,
        hmemo::HArray<ValueType>& cValues,
        const hmemo::HArray<IndexType>& aIA,
        const hmemo::HArray<IndexType>& aJA,
        const hmemo::HArray<ValueType>& aValues,
        const hmemo::HArray<IndexType>& bIA,
        const hmemo::HArray<IndexType>& bJA,
        const hmemo::HArray<ValueType>& bValues,
        const IndexType m,
        const IndexType n,
        const common::BinaryOp op,
        hmemo::ContextPtr prefLoc );

    /** 
     *  @brief Check if the CSR data has for each diagonal element an entry.
     *
     *  Some operations on matrices require that a full diagonal is available
     *  (setDiagonal, jacobi).
     */
    static bool hasDiagonalProperty(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const bool isSorted,
        hmemo::ContextPtr prefLoc );

    /** @brief Get the diagonal of CSR storage
     *
     *  This routine is very efficient if diagonal elements are stored first.
     *  It is also significantly faster if the column indexes are sorted for each row.
     */
    template<typename ValueType>
    static void getDiagonal(
        hmemo::HArray<ValueType>& diagonal,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::HArray<ValueType>& values,
        const bool isSorted,
        hmemo::ContextPtr prefLoc );

    /** @brief set the diagonal */

    template<typename ValueType>
    static bool setDiagonalV(
        hmemo::HArray<ValueType>& values,
        const hmemo::HArray<ValueType>& diagonal,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const bool isSorted,
        hmemo::ContextPtr prefLoc );

    /** @brief set the diagonal */

    template<typename ValueType>
    static bool setDiagonal(
        hmemo::HArray<ValueType>& values,
        const ValueType diagonalValue,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const bool isSorted,
        hmemo::ContextPtr prefLoc );

    /** 
     *  @brief return the position for an entry (i,j) in the CSR data
     *
     *  @param[in] i, j are the row and column index for the searched entry
     *  @param[in] csrIA, csrJA are the row offset array and the column indexes
     *  @param[in] prefLoc specifies the location where operation is done
     *  @return invalidIndex if not found, otherwise k with ja[k] == j, ia[i] <= k < ia[i+1]
     *
     *  The corresponding matrix value can be found via csrValues[k] if k is the not invalid.
     */
    static IndexType getValuePos(
        const IndexType i,
        const IndexType j,
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Identify all entries of a given column in the CSR data
     *
     *  @param[out] ia contains the row indexes of the entries 
     *  @param[out] positions contains the positions of the entries in the arrays csrJA
     *  @param[in] csrIA the row offset array
     *  @param[in] csrJA is the array of column indexes
     *  @param[in] j the column for which entries are needed
     *  @param[in] prefLoc specifies the context where operation should be executed.
     *
     *  Be careful, the entries might be in any order.
     */
    static void getColumnPositions(
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>& positions,
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        const IndexType j,
        const hmemo::ContextPtr prefLoc );

    /**
     *  @brief Jacobi iteration step with CSR storage
     *
     *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
     */
    template<typename ValueType>
    static void jacobi(
        hmemo::HArray<ValueType>& solution,
        const ValueType omega,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::HArray<ValueType>& csrValues,
        hmemo::ContextPtr prefLoc );

    /** Jacobi iteration step using a halo storage.
     *
     *  solution -= omega * ( B(halo) * oldSolution) ./ localDiagonal
     *
     *  @param[in,out] localSolution is the solution vector that is updated
     *  @param[in]     localDiagonal pointer to the diagonal of local storage
     *  @param[in]     oldSolution is the old solution vector of halo part
     *  @param[in]     omega is the scaling factor.
     *  @param[in]     csrIA, csrJA, csrValues are the CSR containers
     *  @param[in]     rowIndexes if not empty it contains row indexes of non-empty rows
     */
    template<typename ValueType>
    static void jacobiHalo(
        hmemo::HArray<ValueType>& localSolution,
        const ValueType omega,
        const hmemo::HArray<ValueType>& localDiagonal,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::HArray<ValueType>& csrValues,
        const hmemo::HArray<IndexType>& rowIndexes,
        hmemo::ContextPtr prefLoc );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
