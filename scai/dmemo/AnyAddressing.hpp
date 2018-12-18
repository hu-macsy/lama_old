/**
 * @file AnyAddressing.hpp
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
 * @brief GeneralDistribution.hpp
 * @author brandes
 * @date 25.02.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/common/SCAITypes.hpp>

#include <scai/hmemo/HArray.hpp>

namespace scai
{

namespace dmemo
{

class Distribution;

/** 
 *   Data structure that contains owner and local offsets for all 'global' indexes

    // the following arrays will only be available if enableAnyAddressing has been called
    // Note: if set the array mGlobal2Local is no more needed

    // mutable hmemo::HArray<PartitionId> mAllOwners;         // will have globalSize entries on each processor
    // mutable hmemo::HArray<IndexType> mAllLocalOffsets;     // local size on each partition
    // mutable hmemo::HArray<IndexType> mAllLocal2Global;     // sorts elements into buckets
    // mutable hmemo::HArray<IndexType> mAllGlobal2Local;     // sorts elements into buckets

    // Example
    // index       0    1    2    3   4    5    6    7   8   9   10   11   12
    // mOwners:    0    1    2    0   2    0    1    0   0   1    1    2    2
    // Offsets:    0                       5                 9                    13
    // perm   :    0    3    5    7   8    1    6    9  10   2    4   11   12     local2Global
    // perm'  :    0    5    9    1  10    2    6    3   4   7    8   11   12     global2Local
    //
    // Note: perm is identity iff we have a block distribution
 */

struct COMMON_DLL_IMPORTEXPORT AnyAddressing
{
    hmemo::HArray<PartitionId> allOwners;         // will have globalSize entries on each processor
    hmemo::HArray<IndexType> allLocalOffsets;     // local size on each partition
    hmemo::HArray<IndexType> allLocal2Global;     // sorts elements into buckets
    hmemo::HArray<IndexType> allGlobal2Local;     // inverse to allLocal2Global

    AnyAddressing( const Distribution& dist );

    IndexType localSize( const PartitionId rank ) const
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( rank, allLocalOffsets.size() - 1, "illegal" )
        return allLocalOffsets[ rank + 1] - allLocalOffsets[rank];
    }

    PartitionId owner( const IndexType globalIndex ) const
    {
        return allOwners[ globalIndex ];
    }

    IndexType localIndex( const IndexType globalIndex, const PartitionId owner ) const
    {
        // here the owner is important as local index  requires size offsets
        return allGlobal2Local[ globalIndex ] - allLocalOffsets[ owner ];
    }

    IndexType globalIndex( const IndexType localIndex, const PartitionId owner ) const
    {
        return allLocal2Global[ localIndex + allLocalOffsets[ owner ] ];
    }
};

} /* end namespace dmemo */

} /* end namespace scai */
