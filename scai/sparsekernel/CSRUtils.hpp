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

    /** @brief get the indexes of non-empty rows 
     *
     *  @param[out] rowIndexes contains the indexes of non-zero rows
     *  @param[in]  the CSR row offset array
     *  @param[in]  threshold builds rowIndexes only if #nonZeroRows/#numRows < threshhold
     *  @param[in]  loc is the context where operation is executed
     *  @returns the number of non-zero rows, will also be the size of rowIndexes if built
     *
     *  If threshhold is 0, the row indexes are never built.
     */
    static IndexType nonEmptyRows( 
        hmemo::HArray<IndexType>& rowIndexes, 
        const hmemo::HArray<IndexType>& csrIA, 
        float threshhold,
        hmemo::ContextPtr loc );

    /** @brief This method generates new CSR data where all zero elemens are removed.
     *
     *  @param[in,out] csrIA, csrJA, csrValues is the CSR data that is compressed
     *  @param[in] eps a value is considered to be zero if abs( value ) <= eps
     *  @param[in] loc specficies the context where compression should be done
     */  
    template<typename ValueType>
    static void compress( 
        hmemo::HArray<IndexType>& csrIA, 
        hmemo::HArray<IndexType>& csrJA, 
        hmemo::HArray<ValueType>& csrValues,
        const RealType<ValueType> eps, 
        hmemo::ContextPtr loc );

    /** 
     *  @brief sort column entries in each row of CSR data
     *
     *  This routine sorts the column indexes for each row in ascending order.
     *
     *  @param[in,out] ja array with the column indexes
     *  @param[in,out] values array with the matrix values
     *  @param[in] ia offset array, size is numRows + 1
     *  @param[in] numRows the number of rows, only needed for convenience
     *  @param[in] numColums the number of columns, only needed for convenience
     *  @param[in] loc specifies the context where the operation should be executed
     */
    template<typename ValueType>
    static void sortRows(
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>& values,
        const hmemo::HArray<IndexType>& ia,
        const IndexType numRows,
        const IndexType numColumns,
        hmemo::ContextPtr loc );

    /**
     *  @brief Check if the column entries of CSR data are sorted 
     *
     *  @param[in] ia is the CSR offset array
     *  @param[in] ja are the column indexes
     *  @param[in] numRows needed for convenience, same as ia.size() - 1
     *  @param[in] numColumns needed for convenience
     *  @param[loc] specifies the context where the operation should be executed
     */
    static bool hasSortedRows(
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const IndexType numRows,
        const IndexType numColumns,
        hmemo::ContextPtr loc );
 
    /**
     *  @brief This routine moves the diagonal entries to the first entry entry for each row
     * 
     *  @param[in,out] ja array with the column indexes
     *  @param[in,out] values array with the matrix values
     *  @param[in] ia offset array, size is numRows + 1
     *  @param[in] numColums the number of columns, only needed for convenience
     *  @param[in] loc specifies the context where the operation should be executed
     *  @return the number of rows where diagonal element is first element
     * 
     *  Only the diagonal element is moved to the first entry, all other elements
     *  remain in the original order. So this routine will not destroy the sorting of
     *  the non-diagonal entries.
     */
    template<typename ValueType>
    static IndexType setDiagonalFirst(
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>& values,
        const hmemo::HArray<IndexType>& ia,
        const IndexType numColumns,
        hmemo::ContextPtr loc );

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
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& rowIA,
        const hmemo::HArray<IndexType>& rowJA,
        const hmemo::HArray<ValueType>& rowValues,
        hmemo::ContextPtr loc );

    /** Matrix multiplication 
     *
     *  @param[out] cIA, cJA, cValues for output storage
     *  @param[in] aIA, aJA, aValues for first input storage
     *  @param[in] bIA, bJA, bValues for second input storage
     *  @param[in] m, n  are the sizes of the output matrix
     *  @param[in] k    are the number of columns of a and number of rows of b
     *  @param[in] diagonalProperty  if true diagonal elements in output storage will be available
     *  @param[in] loc specifies the preferred context for execution.
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
        hmemo::ContextPtr loc );

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
        hmemo::ContextPtr loc );

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
        hmemo::ContextPtr loc );

    /** 
     *  @brief Check if all diagonal elements are available
     */
    static bool hasDiagonalProperty(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const bool isSorted,
        hmemo::ContextPtr loc );

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
        hmemo::ContextPtr loc );

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
        hmemo::ContextPtr loc );

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
        hmemo::ContextPtr loc );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
