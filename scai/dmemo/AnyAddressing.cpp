/**
 * @file AnyAddressing.cpp
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
 * @brief Implementations of methods for any addressing.
 * @author Thomas Brandes
 * @date 17.12.2018
 */

// hpp
#include <scai/dmemo/AnyAddressing.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

namespace scai
{

using utilskernel::HArrayUtils;

namespace dmemo
{

AnyAddressing::AnyAddressing( const Distribution& dist )
{
    // compute mAllOwners

    hmemo::HArray<IndexType> indexes;   // will contain all column indexes to get all owners

    HArrayUtils::setOrder( indexes, dist.getGlobalSize() );

    dist.Distribution::computeOwners( allOwners, indexes );

    // bucket sort the owners, gives offsets and permutation to block values according to owners

    HArrayUtils::bucketSortOffsets( allLocalOffsets, allLocal2Global, allOwners, dist.getCommunicator().getSize() );

    HArrayUtils::inversePerm( allGlobal2Local, allLocal2Global ); // global2Local
}

} /* end namespace dmemo */

} /* end namespace scai */
