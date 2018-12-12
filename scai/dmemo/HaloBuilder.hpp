/**
 * @file HaloBuilder.hpp
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
 * @brief HaloBuilder.hpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 24.03.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/dmemo/Halo.hpp>
#include <scai/dmemo/Distribution.hpp>

namespace scai
{

namespace dmemo
{

class COMMON_DLL_IMPORTEXPORT HaloBuilder
{
public:

    /** Build a halo by an array of required indexes 
     *
     *  @param[out] halo is the Halo object that contains (sorted) required and provides indexes and exchange plans
     *  @param[in]  distribution stands for the mapping of the global indexes
     *  @param[in]  requiredIndexes are global indexes for required values from other processors
     *
     *  Note: requiredIndexes should not contain global indexes that are owned by this processor (isLocal)
     *        and there should be no double values in it to reduce communication volume.
     */
    static void buildFromRequired( Halo& halo, const Distribution& distribution, const hmemo::HArray<IndexType>& requiredIndexes );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace dmemo */

} /* end namespace scai */
