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

namespace lama
{

class COMMON_DLL_IMPORTEXPORT COOUtils
{
public:

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
};

/* -------------------------------------------------------------------------- */

} /* end namespace lama */

} /* end namespace scai */
