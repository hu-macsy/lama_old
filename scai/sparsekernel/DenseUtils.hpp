/**
 * @file DenseUtils.hpp
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
 * @brief Utility functions for dense storage data
 * @author Thomas Brandes
 * @date 28.05.2018
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

class COMMON_DLL_IMPORTEXPORT DenseUtils
{
public:

    /** 
     *  @brief Conversion of dense matrix storage data to the CSR storage format.
     *
     *  \param{in] denseValues is an array of size numRows * numColumns, values are stored row-wise
     *  \param[in] numRows is the number of rows
     *  \param[out] csrIA will be the offset array for the row sizes
     *  \param[out] csrJA contains the columns of the nonz-zero entries
     *  \param[out] csrValues contains the non-zero values
     *  \param[in] loc is the context where conversion should be executed.
     */
    template<typename ValueType>
    static void convertDense2CSR( 
        hmemo::HArray<IndexType>& csrIA, 
        hmemo::HArray<IndexType>& csrJA, 
        hmemo::HArray<ValueType>& csrValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<ValueType>& denseValues,
        hmemo::ContextPtr loc );

    /** 
     *  @brief Determine for each row the number of non-zero entries
     *
     *  \param{in] denseValues is an array of size numRows * numColumns, values are stored row-wise
     *  \param[in] numRows is the number of rows
     *  \param[out] rowSizes will contain for each row the number of non-zero elements
     *  \param[in] loc is the context where conversion should be executed.
     */
    template<typename ValueType>
    static void getSparseRowSizes( 
        hmemo::HArray<IndexType>& rowSizes,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<ValueType>& denseValues,
        hmemo::ContextPtr loc );

    /** 
     *  @brief Determine for each row the number of non-zero entries
     *
     *  \param{out] denseValues is an array of size numRows * numColumns, values are stored row-wise
     *  \param[in] numRows is the number of rows
     */
    template<typename ValueType>
    static void convertCSR2Dense(
        hmemo::HArray<ValueType>& denseValues,
        const IndexType numRows,
        const IndexType numColumns,
        const hmemo::HArray<IndexType>& csrIA,
        const hmemo::HArray<IndexType>& csrJA,
        const hmemo::HArray<ValueType>& csrValues,
        hmemo::ContextPtr loc );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
