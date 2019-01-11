/**
 * @file ELLUtils.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Utility functions for ELL data
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

class COMMON_DLL_IMPORTEXPORT ELLUtils
{
public:

    /** 
     *  @brief Conversion of ELL storage data to the CSR storage format.
     *
     *  Note: the order of the entries in each row remains unchanged.
     */
    template<typename ValueType>
    static void convertELL2CSR(
        hmemo::HArray<IndexType>& csrIA,
        hmemo::HArray<IndexType>& csrJA,
        hmemo::HArray<ValueType>& csrValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ellIA,
        const hmemo::HArray<IndexType>& ellJA,
        const hmemo::HArray<ValueType>& ellValues,
        hmemo::ContextPtr loc );

    /** 
     *  @brief Conversion of ELL storage data to the CSR storage format.
     *
     *  Note: the order of the entries in each row remains unchanged.
     */
    template<typename ValueType>
    static void convertCSR2ELL(
        hmemo::HArray<IndexType>& ellIA,
        hmemo::HArray<IndexType>& ellJA,
        hmemo::HArray<ValueType>& ellValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::HArray<ValueType>& csrValues,
        hmemo::ContextPtr loc );

    /** @brief get the indexes of non-empty rows 
     *
     *  @param[out] rowIndexes contains the indexes of non-zero rows
     *  @param[in]  ellSizes continas the number of non-zero entries of each row
     *  @param[in]  threshold builds rowIndexes only if nonZeroRows / numRows < threshhold
     *  @param[in]  prefLoc is the context where operation is executed
     *  @returns the number of non-zero rows, will also be the size of rowIndexes if built
     *
     *  If threshhold is 0, the row indexes are never built.
     */
    static IndexType nonEmptyRows( 
        hmemo::HArray<IndexType>& rowIndexes, 
        const hmemo::HArray<IndexType>& ellSizes,
        float threshold,
        hmemo::ContextPtr prefLoc );

    /** @brief Get the diagonal entries for an ELL storage
     *
     *  This routine is very efficient if diagonal elements are stored first.
     *
     *  Note: if one entry in diagonalPositions is invalidIndex, the diagonal property is not given
     */
    static IndexType getDiagonalPositions(
        hmemo::HArray<IndexType>& diagonalPositions,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ellIA,
        const hmemo::HArray<IndexType>& ellJA,
        hmemo::ContextPtr prefLoc );

    /** @brief Get the diagonal of ELL storage
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
        hmemo::ContextPtr prefLoc );

    /** @brief set the diagonal */

    template<typename ValueType>
    static void setDiagonalV(
        hmemo::HArray<ValueType>& ellValues,
        const hmemo::HArray<ValueType>& diagonal,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ellIA,
        const hmemo::HArray<IndexType>& ellJA,
        hmemo::ContextPtr prefLoc );

    /** @brief set the diagonal */

    template<typename ValueType>
    static void setDiagonal(
        hmemo::HArray<ValueType>& ellValues,
        const ValueType diagonalValue,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ellIA,
        const hmemo::HArray<IndexType>& ellJA,
        hmemo::ContextPtr prefLoc );

    /** 
     *  @brief return the position for an entry (i,j) in the ELL data
     *
     *  @param[in] i, j are the row and column index for the searched entry
     *  @param[in] ellIA, ellJA are the sizes array and the column indexes
     *  @param[in] prefLoc specifies the context where the operation should be executed
     *  @return invalidIndex if not found, otherwise k with ja[k] == j, k % numRows = i
     *
     *  The corresponding matrix value can be found via ellValues[k] if k is not invalid.
     */
    static IndexType getValuePos(
        const IndexType i,
        const IndexType j,
        const hmemo::HArray<IndexType>& ellIA,
        const hmemo::HArray<IndexType>& ellJA,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Identify all entries of a given column in the ELL data
     *
     *  @param[out] ia contains the row indexes of the entries 
     *  @param[out] positions contains the positions of the entries in the arrays ellJA
     *  @param[in] ellIA the sizes of each row
     *  @param[in] ellJA is the array of column indexes
     *  @param[in] j the column for which entries are needed
     *  @param[in] prefLoc specifies the context where operation should be executed.
     *
     *  The array ia is not really needed as the rows can be identified by the positions.
     *  ia[i] is same as positions[i] % numRows.
     *  
     *  Be careful, the entries might be in any order.
     */
    static void getColumnPositions(
        hmemo::HArray<IndexType>& ia,
        hmemo::HArray<IndexType>& positions,
        const hmemo::HArray<IndexType>& ellIA, 
        const hmemo::HArray<IndexType>& ellJA,
        const IndexType j,
        const hmemo::ContextPtr prefLoc );

    /** @brief This method generates new ELL data where all zero elemens are removed.
     *
     *  @param[in,out] ellIA, ellJA, ellValues is the ELL data that is compressed
     *  @param[in,out] numValuesPerRow size for each row in the arrays ellJA and ellValues
     *  @param[in] eps a value is considered to be zero if abs( value ) <= eps
     *  @param[in] prefLoc specficies the context where compression should be done
     */
    template<typename ValueType>
    static void compress(
        hmemo::HArray<IndexType>& ellIA,
        hmemo::HArray<IndexType>& ellJA,
        hmemo::HArray<ValueType>& ellValues,
        IndexType& numValuesPerRow,
        const RealType<ValueType> eps,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief matrix-vector multiplication, result = alpha * ELLstorage * x + beta * y
     */
    template<typename ValueType>
    static tasking::SyncToken* gemv(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValuesPerRow,
        const hmemo::HArray<IndexType>& ellIA,
        const hmemo::HArray<IndexType>& ellJA,
        const hmemo::HArray<ValueType>& ellValues,
        const common::MatrixOp op,
        bool async,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief matrix-vector multiplication, result = alpha * ELLstorage * x
     */
    template<typename ValueType>
    static tasking::SyncToken* gemv0(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValuesPerRow,
        const hmemo::HArray<IndexType>& ellIA,
        const hmemo::HArray<IndexType>& ellJA,
        const hmemo::HArray<ValueType>& ellValues,
        const common::MatrixOp op,
        bool async,
        hmemo::ContextPtr prefLoc );
    /**
     *  @brief Computes result += alpha * ELLStorage * x
     *
     *  This operation touches only elements of result, for which entries are available.
     *  This routine can be optimized if the array with number of nonZeroRowIndexes is
     *  also availabe.
     */
    template<typename ValueType>
    static tasking::SyncToken* gemvSp(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const IndexType numRows,
        const IndexType numColumns,
        const IndexType numValuesPerRow,
        const hmemo::HArray<IndexType>& ellIA,
        const hmemo::HArray<IndexType>& ellJA,
        const hmemo::HArray<ValueType>& ellValues,
        const common::MatrixOp op,
        const hmemo::HArray<IndexType>& nonZeroRowIndexes,
        bool async,
        hmemo::ContextPtr prefLoc );

    /**
     *  @brief Jacobi iteration step with ELL storage
     *
     *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
     */
    template<typename ValueType>
    static tasking::SyncToken* jacobi(
        hmemo::HArray<ValueType>& solution,
        const ValueType omega,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<ValueType>& rhs,
        const hmemo::HArray<IndexType>& ellIA,
        const hmemo::HArray<IndexType>& ellJA,
        const hmemo::HArray<ValueType>& ellValues,
        bool async,
        hmemo::ContextPtr prefLoc );

    /** Jacobi iteration step using a halo storage.
     *
     *  solution -= omega * ( B(halo) * oldSolution) ./ localDiagonal
     *
     *  @param[in,out] localSolution is the solution vector that is updated
     *  @param[in]     localDiagonal pointer to the diagonal of local storage
     *  @param[in]     oldSolution is the old solution vector of halo part
     *  @param[in]     omega is the scaling factor.
     *  @param[in]     ellIA, ellJA, ellValues are the ELL containers
     *  @param[in]     rowIndexes if not empty it contains row indexes of non-empty rows
     *  @param[in]     prefLoc specifies the context where iteration should be done 
     */
    template<typename ValueType>
    static void jacobiHalo(
        hmemo::HArray<ValueType>& localSolution,
        const ValueType omega,
        const hmemo::HArray<ValueType>& localDiagonal,
        const hmemo::HArray<ValueType>& oldSolution,
        const hmemo::HArray<IndexType>& ellIA,
        const hmemo::HArray<IndexType>& ellJA,
        const hmemo::HArray<ValueType>& ellValues,
        const hmemo::HArray<IndexType>& rowIndexes,
        hmemo::ContextPtr prefLoc );

    /**
     * @brief Set/update each row in an ELL storage with an individual value
     *
     *  \code
     *      for all i = 0, ..., n-1; j = 0, ..., m-1
     *      ellValues( i, j ) = ellValues( i, j ) <op> rowValues( i )
     *  \endcode
     */
    template<typename ValueType>
    static void setRows(
        hmemo::HArray<ValueType>& ellValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ellIA,
        const hmemo::HArray<IndexType>& ellJA,
        const hmemo::HArray<ValueType>& rowValues,
        const common::BinaryOp op,
        hmemo::ContextPtr prefLoc );

    /**
     * @brief Set/update each column in an ELL storage with an individual value
     *
     *  \code
     *      for all i = 0, ..., n-1; j = 0, ..., m-1
     *      ellValues( i, j ) = ellValues( i, j ) <op> columnValues( j )
     *  \endcode
     */
    template<typename ValueType>
    static void setColumns(
        hmemo::HArray<ValueType>& ellValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& ellIA,
        const hmemo::HArray<IndexType>& ellJA,
        const hmemo::HArray<ValueType>& columnValues,
        const common::BinaryOp op,
        hmemo::ContextPtr prefLoc );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
