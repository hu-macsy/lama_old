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

    // TODO: Rename to buildFromRequired
    static void build( const Distribution& distribution, const hmemo::HArray<IndexType>& requiredIndexes, Halo& halo );

    /**
     * Build a Halo data structure from prior knowledge about the owners of provided indexes.
     *
     * Because you might want to build a Halo data structure for only a subset of the local indices
     * of a distribution, this routine works with an arbitrary "halo set". That is, the provided indexes
     * of the returned Halo will be local indexes in this halo set. The halo set is related to the set
     * of global indexes by the array halo2global, which maps local halo indexes to global indexes.
     * The (global) required indexes of the Halo are thus defined by this halo2global mapping.
     *
     * As an example, if halo2global[i] == distribution.local2global(i), then the halo set would correspond
     * to the set of local indexes in the distribution.
     *
     * @param[in]  comm             The communicator through which communication with other processors will take place.
     * @param[in]  halo2global      An array mapping halo indexes to global indexes, i.e. halo2global[i] is the global index
     *                              corresponding to halo index i.
     * @param[in]  ownersOfProvided An array mapping local (halo), provided indices to owners. More precisely, ownersOfProvided[i] must contain
     *                              the owner of halo index i. Note that ownersOfProvided must have the same length as halo2global.
     * @param[out] halo             The Halo data structure to store the result in.
     */
    static void buildFromProvidedOwners( const Communicator& comm,
                                         const hmemo::HArray<IndexType>& halo2global,
                                         const hmemo::HArray<PartitionId>& ownersOfProvided,
                                         Halo& halo );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace dmemo */

} /* end namespace scai */
