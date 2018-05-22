/**
 * @file COOUtils.hpp
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
 * @brief Utility functions for COO data
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

class COMMON_DLL_IMPORTEXPORT COOUtils
{
public:

    /**
     *  @brief Convert the COO ia array to CSR ia offset array
     *
     *  \code
     *    cooIA = { 0, 0, 1, 1, 1, 2, 4, 4, 5, 6, 6 } -> csrIA = { 0, 2, 5, 6, 6, 8, 9, 11 } 
     *  \endcode
     */
    static void convertCOO2CSR( 
        hmemo::HArray<IndexType>& csrIA, 
	const hmemo::HArray<IndexType>& cooIA, 
	const IndexType numRows,
	hmemo::ContextPtr loc );

    /**
     *  @brief convert the CSR offset array to an COO ia array
     *
     *  @param[in]  csrIA contains the offset array, number of rows is csrIA.size() - 1
     *  @param[in]  nnz are the number of non-zero elements, must be same as csrIA[ numRows ]
     *  @param[out] cooIA contains the row positions, size will be nnz 
     *
     *  \code
     *    csrIA = { 0, 2, 5, 6, 6, 8, 9, 11 } -> cooIA = { 0, 0, 1, 1, 1, 2, 4, 4, 5, 6, 6 } 
     *  \endcode
     */
    static void convertCSR2COO( 
        hmemo::HArray<IndexType>& cooIA, 
	const hmemo::HArray<IndexType>& csrIA, 
        const IndexType nnz,
	hmemo::ContextPtr loc );

    /**
     *  @brief check if COO arrays ia and ja are sorted.
     */
    static bool isSorted(
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        hmemo::ContextPtr prefContext );

    /** 
     *  @brief sort COO data
     *
     *  The sorting is stable, i.e. entries with same coordinates keep their order.
     */
    template<typename ValueType>
    static void sort(
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>& values );

    /**
     *  @brief This method eliminates consecutive entries with same coordinates.
     *
     *  If op is COPY, only the latest element is taken. Otherwise the elements are 
     *  combined corresponding to the specified operation.
     *
     *  Only if the COO data is sorted, it can be
     *  guaranteed that COO data does not contain any two elements with same coordinates.
     */
    template<typename ValueType>
    static void unique(
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>& values,
        common::BinaryOp op = common::BinaryOp::COPY );

    static IndexType getValuePos( 
        const IndexType i, 
        const IndexType j,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief This method checks if COO indexes have an entry for each diagonal element
     *
     *  If diagonal property is not given, some operations on COO storage might fail, e.g.
     *  jacobi, setDiagonal.
     */
    static bool hasDiagonalProperty(
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        hmemo::ContextPtr prefLoc );

    /** @brief Get the diagonal of COO storage
     *
     *  This routine might not work correctly if the arrays (ia, ja) are not sorted
     */
    template<typename ValueType>
    static void getDiagonal(
        hmemo::HArray<ValueType>& diagonal,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        const hmemo::HArray<ValueType>& values,
        hmemo::ContextPtr prefLoc );

    /** @brief Set the diagonal of COO storage
     *
     */
    template<typename ValueType>
    static void setDiagonalV(
        hmemo::HArray<ValueType>& values,
        const hmemo::HArray<ValueType>& diagonal,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        hmemo::ContextPtr prefLoc );

    template<typename ValueType>
    static void setDiagonal(
        hmemo::HArray<ValueType>& values,
        const ValueType diagonal,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        hmemo::ContextPtr prefLoc );
};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
