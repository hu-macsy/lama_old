/**
 * @file ELLUtils.hpp
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
 * @brief Utility functions for ELL data
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

class COMMON_DLL_IMPORTEXPORT ELLUtils
{
public:

    /** @brief get the indexes of non-empty rows 
     *
     *  @param[out] rowIndexes contains the indexes of non-zero rows
     *  @param[in]  ellSizes continas the number of non-zero entries of each row
       @param[in]  threshold builds rowIndexes only if #nonZeroRows/#numRows < threshhold
     *  @param[in]  loc is the context where operation is executed
     *  @returns the number of non-zero rows, will also be the size of rowIndexes if built
     *
     *  If threshhold is 0, the row indexes are never built.
     */
    static IndexType nonEmptyRows( 
        hmemo::HArray<IndexType>& rowIndexes, 
        const hmemo::HArray<IndexType>& ellSizes,
        float threshhold,
        hmemo::ContextPtr loc );

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
        hmemo::ContextPtr loc );

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
        hmemo::ContextPtr loc );

    /** @brief This method generates new ELL data where all zero elemens are removed.
     *
     *  @param[in,out] ellIA, ellJA, ellValues is the ELL data that is compressed
     *  @param[in,out] numValuesPerRow size for each row in the arrays ellJA and ellValues
     *  @param[in] eps a value is considered to be zero if abs( value ) <= eps
     *  @param[in] loc specficies the context where compression should be done
     */
    template<typename ValueType>
    static void compress(
        hmemo::HArray<IndexType>& ellIA,
        hmemo::HArray<IndexType>& ellJA,
        hmemo::HArray<ValueType>& ellValues,
        IndexType& numValuesPerRow,
        const RealType<ValueType> eps,
        hmemo::ContextPtr loc );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
