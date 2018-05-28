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

// type for diagonal offsets should be a signed type.

typedef IndexType OffsetType;

class COMMON_DLL_IMPORTEXPORT DIAUtils
{
public:

    /** 
     *  @brief Conversion of DIA storage data to the CSR storage format.
     *
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

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
