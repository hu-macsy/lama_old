/**
 * @file JDSUtils.hpp
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
 * @brief Utility functions for JDS data
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

/** 
 *  @brief Class that provides typical operations for JDS storage arrays.
 *
 */
class COMMON_DLL_IMPORTEXPORT JDSUtils
{
public:

    /** @brief Get the diagonal entries for an JDS storage
     *
     *  This routine is very efficient if diagonal elements are stored first.
     *
     *  Note: if one entry in diagonalPositions is invalidIndex, the diagonal property is not given
     */
    static IndexType getDiagonalPositions(
        hmemo::HArray<IndexType>& diagonalPositions,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& jdsILG,
        const hmemo::HArray<IndexType>& jdsDLG,
        const hmemo::HArray<IndexType>& jdsPerm,
        const hmemo::HArray<IndexType>& jdsJA,
        hmemo::ContextPtr loc );

    /** @brief Get the diagonal of JDS storage
     *
     *  This routine is very efficient if diagonal elements are stored first.
     *  It is also significantly faster if the column indexes are sorted for each row.
     */
    template<typename ValueType>
    static void getDiagonal(
        hmemo::HArray<ValueType>& diagonal,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& jdsILG,
        const hmemo::HArray<IndexType>& jdsDLG,
        const hmemo::HArray<IndexType>& jdsPerm,
        const hmemo::HArray<IndexType>& jdsJA,
        const hmemo::HArray<ValueType>& jdsValues,
        hmemo::ContextPtr loc );

    /** 
     *  @brief return the position for an entry (i,j) in the JDS data
     *
     *  @param[in] i, j are the row and column index for the searched entry
     *  @param[in] jdsILG, jdsDLG, jdsPerm, jdsJA are the corresponding arrays of JDS format
     *  @return invalidIndex if not found, otherwise k with ja[k] == j, k % numRows = i
     *
     *  The corresponding matrix value can be found via csrValues[k] if k is the not invalid.
     */
    static IndexType getValuePos(
        const IndexType i,
        const IndexType j,
        const hmemo::HArray<IndexType>& jdsILG,
        const hmemo::HArray<IndexType>& jdsDLG,
        const hmemo::HArray<IndexType>& jdsPerm,
        const hmemo::HArray<IndexType>& jdsJA,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Identify all entries of a given column in the CSR data
     *
     *  @param[out] ia contains the row indexes of the entries 
     *  @param[out] positions contains the positions of the entries in the arrays jdsJA, jdsValues
     *  @param[in] jdsILG, jdsDLG, jdsPerm, jdsJA are the corresponding arrays of the JDS format
     *  @param[in] j the column for which entries are needed
     *  @param[in] prefLoc specifies the context where operation should be executed.
     *
     *  Be careful, the entries might be in any order.
     */
    static void getColumnPositions(
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>& positions,
        const hmemo::HArray<IndexType>& jdsILG,
        const hmemo::HArray<IndexType>& jdsDLG,
        const hmemo::HArray<IndexType>& jdsPerm,
        const hmemo::HArray<IndexType>& jdsJA,
        const IndexType j,
        hmemo::ContextPtr prefLoc );

    /** 
     *  @brief Compute the dlg array by the ilg array for the JDS format.
     */
    static void ilg2dlg(
        hmemo::HArray<IndexType>& jdsDLG,
        const IndexType numDiagonals,
        const hmemo::HArray<IndexType>& jdsILG,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Jacobi iteration step with JDS storage
     *
     *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
     */
    template<typename ValueType>
    static void jacobi(
        hmemo::HArray<ValueType>& solution,
        const ValueType omega,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const hmemo::HArray<IndexType>& jdsILG,
        const hmemo::HArray<IndexType>& jdsDLG,
        const hmemo::HArray<IndexType>& jdsPerm,
        const hmemo::HArray<IndexType>& jdsJA,
        const hmemo::HArray<ValueType>& jdsValues,
        hmemo::ContextPtr prefLoc );
    /**
     *  @brief Jacobi halo iteration step with JDS storage
     *
     *  solution -= omega * ( B(halo) * oldSolution) ./ diagonal
     */
    template<typename ValueType>
    static void jacobiHalo(
        hmemo::HArray<ValueType>& solution,
        const ValueType omega,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& diagonal,
        const hmemo::HArray<IndexType>& jdsILG,
        const hmemo::HArray<IndexType>& jdsDLG,
        const hmemo::HArray<IndexType>& jdsPerm,
        const hmemo::HArray<IndexType>& jdsJA,
        const hmemo::HArray<ValueType>& jdsValues,
        hmemo::ContextPtr prefLoc );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
