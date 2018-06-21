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

#include <scai/tasking/SyncToken.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/BinaryOp.hpp>
#include <scai/common/MatrixOp.hpp>

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
     *  @param[out] csrIA is the offset array
     *  @param[in]  cooIA is the 'sorted' array with the row indexes
     *  @param[in]  numRows number of rows
     *  @param[in]  prefLoc specifies the context where operation should be executed
     *
     *  \code
     *    cooIA = { 0, 0, 1, 1, 1, 2, 4, 4, 5, 6, 6 } -> csrIA = { 0, 2, 5, 6, 6, 8, 9, 11 } 
     *  \endcode
     */
    static void convertCOO2CSR( 
        hmemo::HArray<IndexType>& csrIA, 
        const hmemo::HArray<IndexType>& cooIA, 
        const IndexType numRows,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief convert the CSR offset array to an COO ia array
     *
     *  @param[in]  csrIA contains the offset array, number of rows is csrIA.size() - 1
     *  @param[in]  nnz are the number of non-zero elements, must be same as csrIA[ numRows ]
     *  @param[out] cooIA contains the row positions, size will be nnz 
     *  @param[in]  prefLoc specifies the context where operation should be executed
     *
     *  \code
     *    csrIA = { 0, 2, 5, 6, 6, 8, 9, 11 } -> cooIA = { 0, 0, 1, 1, 1, 2, 4, 4, 5, 6, 6 } 
     *  \endcode
     */
    static void convertCSR2COO( 
        hmemo::HArray<IndexType>& cooIA, 
        const hmemo::HArray<IndexType>& csrIA, 
        const IndexType nnz,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief check if COO arrays ia and ja are sorted.
     */
    static bool isSorted(
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<IndexType>& ja,
        hmemo::ContextPtr prefContext );

    /** 
     *  @brief normalize COO data
     *
     *  Normalizes COO data, sorts it and eliminates double entries.
     */
    template<typename ValueType>
    static void normalize(
        hmemo::HArray<IndexType>& cooIA,
        hmemo::HArray<IndexType>& cooJA,
        hmemo::HArray<ValueType>& cooValues,
        common::BinaryOp,
        hmemo::ContextPtr prefLoc );

    /** 
     *  @brief sort COO data
     *
     *  The sorting is stable, i.e. entries with same coordinates keep their order.
     */
    template<typename ValueType>
    static void sort(
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>& values,
        hmemo::ContextPtr prefLoc );

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
        common::BinaryOp op,
        hmemo::ContextPtr prefLoc );

    /** 
     *  @brief return the position for an entry (i,j) in the COO data
     *
     *  @param[in] i, j are the row and column index for the searched entry
     *  @param[in] cooIA, cooJA are the sorted row/column indexes of the non-zero entries
     *  @param[in]  prefLoc specifies the context where operation should be executed
     *  @return invalidIndex if not found, otherwise k with ia[k] == i & ja[k] == j
     *
     *  When the position is known, the corresponding matrix value entry can be read or updated.
     *
     *  \code
     *    IndexType pos = COOUtils::getValuePos( i, j, cooIA, cooJA, loc );
     *    if ( pos != invalidIndex )
     *    {
     *        ValueType v = cooValues[pos];
     *        std::cout << "Value at (" << i << ", " << j << ") is " << v << std::endl;
     *    }
     *  \endcode
     */
    static IndexType getValuePos( 
        const IndexType i, 
        const IndexType j,
        const hmemo::HArray<IndexType>& cooIA,
        const hmemo::HArray<IndexType>& cooJA,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Get all positions with entries for a given column j
     *
     *  @param[out] positions will contain all positions for entries with colJA[pos[k]] == j
     *  @param[in]  cooJA are the column indexes of the non-zero entries
     *  @param[in]  j is the column for which non-zero entries are selected
     *  @param[in]  prefLoc specifies the context where operation should be executed
     *
     *  The row indexes and the values might be gathered from cooIA and cooValues via the positions.
     */
    static void getColumnPositions(
        hmemo::HArray<IndexType>& positions,
        const hmemo::HArray<IndexType>& cooJA,
        const IndexType j,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Get the positions of all non-zero entries for given row i
     *
     *  @param[out] offset is the first position, invalidIndex if not found
     *  @param[out] n      is the number of available entries for row i
     *  @param[in]  cooIA   the (sorted) array with the row indexes
     *  @param[in]  i      is the index of row for which entries are queried
     *  @param[in]  prefLoc specifies the context where operation should be executed
     *
     *  As the entries of one row are stored contiguously, here is no need for a position array
     */
    static void getRowPositions(
        IndexType& offset,
        IndexType& n,
        const hmemo::HArray<IndexType>& cooIA,
        const IndexType i,
        hmemo::ContextPtr prefLoc );
    /**
     *  @brief This method checks if COO arrays have an entry for each diagonal element
     *
     *  If not all diagonal elements are available, some operations on COO storage might fail, e.g.
     *  jacobi, setDiagonal.
     */
    static bool hasDiagonal(
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

    /**
     * @brief Set the diagonal of COO storage
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

    /** 
     *  @brief Scale each row of a COO storage with an individual values
     * 
     *  @param[in,out] values are the non-zero entries of the storage
     *  @param[in]     ia are the row indexes of the non-zero entries
     *  @param[in]     scale are the scaling factors, size is number of rows
     *  @param[in]     prefLoc is context where the operation should be done
     */
    template<typename ValueType>
    static void scaleRows(
        hmemo::HArray<ValueType>& values,
        const hmemo::HArray<IndexType>& ia,
        const hmemo::HArray<ValueType>& scale,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Apply a binary operation element-wise for sparse storage
     */
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
     *  @brief Jacobi iteration step with COO storage
     *
     *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
     */
    template<typename ValueType>
    static void jacobi(
        hmemo::HArray<ValueType>& solution,
        const ValueType omega,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const hmemo::HArray<IndexType>& cooIA,
        const hmemo::HArray<IndexType>& cooJA,
        const hmemo::HArray<ValueType>& cooValues,
        hmemo::ContextPtr prefLoc );

    /** Jacobi iteration step using a COO halo storage.
     *
     *  \f[ solution -= omega * ( B(halo) * oldSolution) ./ localDiagonal \f]
     *
     *  @param[in,out] localSolution is the solution vector that is updated
     *  @param[in]     localDiagonal pointer to the diagonal of local storage
     *  @param[in]     oldSolution is the old solution vector of halo part
     *  @param[in]     omega is the scaling factor.
     *  @param[in]     cooIA, cooJA, cooValues are the COO containers
     *  @param[in]     prefLoc is context where the operation should be done
     */
    template<typename ValueType>
    static void jacobiHalo(
        hmemo::HArray<ValueType>& localSolution,
        const ValueType omega,
        const hmemo::HArray<ValueType>& localDiagonal,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<IndexType>& cooIA,
        const hmemo::HArray<IndexType>& cooJA,
        const hmemo::HArray<ValueType>& cooValues,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief matrix-vector multiplication 
     */
    template<typename ValueType>
    static tasking::SyncToken* gemv(
        hmemo::HArray<ValueType>& result,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& cooIA,
        const hmemo::HArray<IndexType>& cooJA,
        const hmemo::HArray<ValueType>& cooValues,
        const common::MatrixOp op,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        bool async,
        hmemo::ContextPtr prefLoc );

};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
