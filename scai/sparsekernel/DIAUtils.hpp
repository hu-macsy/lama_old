/**
 * @file DIAUtils.hpp
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
 * @brief Utility functions for the DIAGONAL storage format.
 * @author Thomas Brandes
 * @date 28.05.2018
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/tasking/SyncToken.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/MatrixOp.hpp>

namespace scai
{

namespace sparsekernel
{

// type for diagonal offsets should be a signed type.

typedef IndexType OffsetType;

/**
 *  This class provides a lot of utility functions regarding the DIA storage format
 *  that contains the following arrays:
 *
 *   - array offset contains the offsets for the stored diagonals, its size stands for 
 *     the number of diagonals.
 *   - array values contains the entries for the stored diagonals. Therefore it is a
 *     two-dimensional array, the diagonals are stored column-wise.
 *
 *   An entry values( i, d )  stands for entry ( i, i + offsets[d] ) in the matrix.
 */
class COMMON_DLL_IMPORTEXPORT DIAUtils
{
public:

    /** 
     *  @brief Conversion of DIA storage data to the CSR storage format.
     */
    template<typename ValueType>
    static void convertDIA2CSR( 
        hmemo::HArray<IndexType>& csrIA, 
        hmemo::HArray<IndexType>& csrJA, 
        hmemo::HArray<ValueType>& csrValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<OffsetType>& diaOffset,
        const hmemo::HArray<ValueType>& diaValues,
        hmemo::ContextPtr loc );

    /** 
     *  @brief Get the number of non-zero entries for each row
     */
    template<typename ValueType>
    static void getCSRSizes( 
        hmemo::HArray<IndexType>& csrIA, 
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<OffsetType>& diaOffset,
        const hmemo::HArray<ValueType>& diaValues,
        hmemo::ContextPtr loc );

    /** 
     *  @brief Convert CSR storage data to DIA storage data
     */
    template<typename ValueType>
    static void convertCSR2DIA(
        hmemo::HArray<OffsetType>& diaOffset,
        hmemo::HArray<ValueType>& diaValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::HArray<ValueType>& csrValues,
        hmemo::ContextPtr loc );

    /**
     *  @brief Determine the offsets of all used diagonals for CSR data
     */
    template<typename ValueType>
    static void getDIAOffset(
        hmemo::HArray<IndexType>& diaOffset,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::HArray<ValueType>& csrValues,
        hmemo::ContextPtr );

    /**
     *  @brief Get the position for the matrix entry (i, j) in the values array.
     */
    static IndexType getValuePos(
        const IndexType i,
        const IndexType j,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<OffsetType>& diaOffset,
        hmemo::ContextPtr );

    /**
     *  @brief Get the positions for the matrix row (i, : ) 
     * 
     *  @param[out] indexes contains all column indexes j where entry (i, j ) is available
     *  @param[out] positions contains the positions in values array where to find the matrix value
     *  @param[in]  i is the index of the row for which positions are required
     *  @param[in]  numRows is the number of row for the storage
     *  @param[in]  numColumns is the number of columns for the storage
     *  @param[in]  diaOffset are the offsets of the available diagonals
     *  @param[in]  prefLoc specifies context where operation should be done
     */
    static void getRowPositions(
        hmemo::HArray<IndexType>& indexes,
        hmemo::HArray<IndexType>& positions,
        const IndexType i,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<OffsetType>& diaOffset,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Get the positions for the matrix columns (:, j ) 
     * 
     *  @param[out] indexes contains all row indexes i where entry (i, j ) is available
     *  @param[out] positions contains the positions in values array where to find the matrix value
     *  @param[in]  j is the index of the columns for which positions are required
     *  @param[in]  numRows is the number of row for the storage
     *  @param[in]  numColumns is the number of columns for the storage
     *  @param[in]  diaOffset are the offsets of the available diagonals
     *  @param[in]  prefLoc specifies context where operation should be done
     */
    static void getColPositions(
        hmemo::HArray<IndexType>& indexes,
        hmemo::HArray<IndexType>& positions,
        const IndexType j,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<OffsetType>& diaOffset,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Jacobi iteration step with DIA storage
     *
     *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
     */
    template<typename ValueType>
    static void jacobi(
        hmemo::HArray<ValueType>& solution,
        const ValueType omega,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const IndexType n,
        const hmemo::HArray<OffsetType>& diaOffset,
        const hmemo::HArray<ValueType>& diaValues,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Jacobi halo iteration step with DIA storage
     *
     *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
     */
    template<typename ValueType>
    static void jacobiHalo(
        hmemo::HArray<ValueType>& solution,
        const ValueType omega,
        const hmemo::HArray<ValueType>& diagonal,
        const hmemo::HArray<ValueType>& oldSolution,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<OffsetType>& diaOffset,
        const hmemo::HArray<ValueType>& diaValues,
        hmemo::ContextPtr prefLoc );

    template<typename ValueType>
    static RealType<ValueType> maxNorm(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<OffsetType>& diaOffset,
        const hmemo::HArray<ValueType>& diaValues,
        hmemo::ContextPtr prefLoc );

    template<typename ValueType>
    static tasking::SyncToken* gemv(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<OffsetType>& diaOffset,
        const hmemo::HArray<ValueType>& diaValues,
        const common::MatrixOp op,
        bool async,
        hmemo::ContextPtr prefLoc );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
