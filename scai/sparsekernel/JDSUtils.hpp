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

#include <scai/tasking/SyncToken.hpp>

#include <scai/tasking/SyncToken.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/MatrixOp.hpp>

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

    /** Build the row sizes by using JDS data */

    static void buildRowSizes(
        hmemo::HArray<IndexType>& rowSizes,
        const hmemo::HArray<IndexType>& jdsIlg,
        const hmemo::HArray<IndexType>& jdsPerm,
        hmemo::ContextPtr loc );

    /** 
     *  @brief Conversion of JDS storage data to the CSR storage format.
     *
     *  Note: the order of the entries in each row remains unchanged.
     */
    template<typename ValueType>
    static void convertJDS2CSR(
        hmemo::HArray<IndexType>& csrIA,
        hmemo::HArray<IndexType>& csrJA,
        hmemo::HArray<ValueType>& csrValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& jdsIlg,
        const hmemo::HArray<IndexType>& jdsDlg,
        const hmemo::HArray<IndexType>& jdsPerm,
        const hmemo::HArray<IndexType>& jdsJA,
        const hmemo::HArray<ValueType>& jdsValues,
        hmemo::ContextPtr loc );

    /** 
     *  @brief Conversion of CSR storage data to the JDS storage format.
     *
     *  Note: the order of the entries in each row remains unchanged.
     */
    template<typename ValueType>
    static void convertCSR2JDS(
        hmemo::HArray<IndexType>& jdsIlg,
        hmemo::HArray<IndexType>& jdsDlg,
        hmemo::HArray<IndexType>& jdsPerm,
        hmemo::HArray<IndexType>& jdsJA,
        hmemo::HArray<ValueType>& jdsValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::HArray<ValueType>& csrValues,
        hmemo::ContextPtr loc );

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
        hmemo::ContextPtr prefLoc );

    /** @brief Get the diagonal of JDS storage
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
        hmemo::ContextPtr prefLoc );

    /** 
     *  @brief return the position for an entry (i,j) in the JDS data
     *
     *  @param[in] i, j are the row and column index for the searched entry
     *  @param[in] jdsILG, jdsDLG, jdsPerm, jdsJA are the corresponding arrays of JDS format
     *  @param[in] prefLoc specifies the context where the operation should be executed
     *  @return invalidIndex if not found, otherwise k with ja[k] == j, k % numRows = i
     *
     *  The corresponding matrix value can be found via csrValues[k] if k is not invalid.
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
     *  @brief Identify all entries of a given row in the JDS data
     *
     *  @param[out] positions contains the positions of the entries in the arrays jdsJA, jdsValues
     *  @param[in] jdsILG, jdsDLG, jdsPerm are the corresponding arrays of the JDS format
     *  @param[in] i the row for which entries are needed
     *  @param[in] prefLoc specifies the context where operation should be executed.
     *
     *  The following code shows how to get sparse arrays for a row i of the JDS storage.
     *
     *  \code
     *      HArray<IndexType> pos;   // positions of row entries
     *
     *      JDSUtils::getRowPositions( pos, jdsIlg, jdsDlg, jdsPerm, i, getContextPtr() );
     *
     *      HArray<IndexType> jA;      // available column indexes for row i
     *      HArray<IndexType> values;  // corresponding values of row i
     *
     *      HArrayUtils::gather( jA, jdsJA, pos, BinaryOp::COPY, loc );
     *      HArrayUtils::gather( values, jdsValues, pos, BinaryOp::COPY, loc );
     *  \endcode
     */
    static void getRowPositions(
        hmemo::HArray<IndexType>& positions,
        const hmemo::HArray<IndexType>& jdsILG,
        const hmemo::HArray<IndexType>& jdsDLG,
        const hmemo::HArray<IndexType>& jdsPerm,
        const IndexType i,
        hmemo::ContextPtr prefLoc );

    /** 
     *  @brief Compute the dlg array by the ilg array for the JDS format.
     *  
     *  @param[in] jdsILG contains sizes of all rows (weak descending, as rows are sorted by size)
     *  @param[out] jdsDLG where jdsDLG[i] contains number of rows that have more than i entries
     *  @param[in] prefLoc is the context where operation should be executed
     *  @return total number of entries, is sum(jdsDLG), same as sum(jdsILG)
     */
    static IndexType ilg2dlg(
        hmemo::HArray<IndexType>& jdsDLG,
        const hmemo::HArray<IndexType>& jdsILG,
        hmemo::ContextPtr prefLoc );

    template<typename ValueType>
    static void getRow(
        hmemo::HArray<ValueType>& row,
        const IndexType numColumns,
        const IndexType i,
        const hmemo::HArray<IndexType>& jdsIlg,
        const hmemo::HArray<IndexType>& jdsDlg,
        const hmemo::HArray<IndexType>& jdsPerm,
        const hmemo::HArray<IndexType>& jdsJA,
        const hmemo::HArray<ValueType>& jdsValues,
        hmemo::ContextPtr prefLoc );

    template<typename ValueType>
    static void setRow(
        hmemo::HArray<ValueType>& jdsValues,
        const IndexType i,
        const hmemo::HArray<ValueType>& row,
        const hmemo::HArray<IndexType>& jdsIlg,
        const hmemo::HArray<IndexType>& jdsDlg,
        const hmemo::HArray<IndexType>& jdsPerm,
        const hmemo::HArray<IndexType>& jdsJA,
        const common::BinaryOp op,
        hmemo::ContextPtr prefLoc );

    template<typename ValueType>
    static tasking::SyncToken* gemv0(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& jdsIlg,
        const hmemo::HArray<IndexType>& jdsDlg,
        const hmemo::HArray<IndexType>& jdsPerm,
        const hmemo::HArray<IndexType>& jdsJA,
        const hmemo::HArray<ValueType>& jdsValues,
        const common::MatrixOp op,
        const bool async,
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
        const hmemo::HArray<IndexType>& jdsIlg,
        const hmemo::HArray<IndexType>& jdsDlg,
        const hmemo::HArray<IndexType>& jdsPerm,
        const hmemo::HArray<IndexType>& jdsJA,
        const hmemo::HArray<ValueType>& jdsValues,
        const common::MatrixOp op,
        const bool async,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Jacobi iteration step with JDS storage
     *
     *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
     */
    template<typename ValueType>
    static tasking::SyncToken* jacobi(
        hmemo::HArray<ValueType>& solution,
        const ValueType omega,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const hmemo::HArray<IndexType>& jdsILG,
        const hmemo::HArray<IndexType>& jdsDLG,
        const hmemo::HArray<IndexType>& jdsPerm,
        const hmemo::HArray<IndexType>& jdsJA,
        const hmemo::HArray<ValueType>& jdsValues,
        bool async,
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

    /**
     * @brief Set/update each row in a JDS storage with an individual value
     *
     *  \code
     *      for all i = 0, ..., n-1; j = 0, ..., m-1
     *      jdsValues( i, j ) = jdsValues( i, j ) <op> rowValues( i )
     *  \endcode
     */
    template<typename ValueType>
    static void setRows(
        hmemo::HArray<ValueType>& jdsValues,
        const hmemo::HArray<IndexType>& jdsILG,
        const hmemo::HArray<IndexType>& jdsDLG,
        const hmemo::HArray<IndexType>& jdsPerm,
        const hmemo::HArray<ValueType>& rowValues,
        const common::BinaryOp op,
        hmemo::ContextPtr prefLoc );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
