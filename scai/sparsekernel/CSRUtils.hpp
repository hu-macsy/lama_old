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

    /** 
     *  @brief sort column entries in each row of CSR data
     *
     *  This routine sorts the column indexes for each row in ascending order.
     *
     *  @param[in,out] ja array with the column indexes
     *  @param[in,out] values array with the matrix values
     *  @param[in] ia offset array, size is numRows + 1
     *  @param[in] numColums the number of columns, only needed for convenience
     *  @param[in] keepDiagonalFirst if true the first element remains diagonal element
     *  @param[in] loc specifies the context where the operation should be executed
     */
    template<typename ValueType>
    static void sort(
        hmemo::HArray<IndexType>& ja,
        hmemo::HArray<ValueType>& values,
        const hmemo::HArray<IndexType>& ia,
        const IndexType numColumns,
        bool keepDiagonalFirst,
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

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
