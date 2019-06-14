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
 * @brief Data structure that can be used for global addressing without any communication
 * @author brandes
 * @date 14.11.2018
 */

#pragma once

#include <scai/common/config.hpp>


#include <scai/hmemo/HArray.hpp>
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace dmemo
{

class Distribution;

/** 
 *   Data structure that provides owner and local offsets for all 'global' indexes
 */
struct COMMON_DLL_IMPORTEXPORT AnyAddressing
{
    /** Construct arrays for any addressing from an aribtrary distribution. */

    AnyAddressing( const Distribution& dist );

    /** 
     * @brief Get the local size for a given processor.
     */
    inline IndexType localSize( const PartitionId rank ) const;

    /** 
     * @brief Get the owner for a global index.
     */
    inline PartitionId owner( const IndexType globalIndex ) const;

    /** 
     * @brief Get the local index for a global index.
     *
     * @param[in] globalIndex is the queried global index
     * @param[in] owner must be owning processor, i.e. owner( globalIndex )
     * @returns   the local index on the owner processor
     */
    inline IndexType localIndex( const IndexType globalIndex, const PartitionId owner ) const;

    /** 
     * @brief Get the global index for a local index on any processor
     *
     * @param[in] localIndex is the local index 
     * @param[in] owner is the processor for which the query is done
     * @returns   the global index 
     */
    inline IndexType globalIndex( const IndexType localIndex, const PartitionId owner ) const;

private:

    /** Array that contains for each global index the owner 
     *
     *  \code
     *     allOwners = {  0,   1,   2,   0,  2,   0,   1,   0,  0,  1,   1,   2,   2 }
     *  \endcode
     */
    hmemo::HArray<PartitionId> allOwners;         // will have globalSize entries on each processor

    /** Running offsets for all local sizes on each processor */

    hmemo::HArray<IndexType> allLocalOffsets;     // local size on each partition

    /** Array that contains for each processor contiguously the owned indexes
     *
     *  \code
     *     allOwners       = {  0,   1,   2,   0,  2,   0,   1,   0,  0,   1,   1,   2,   2  }
     *     allLocal2Global = {  0,   3,   5,   7,  8 |  1,   6,   9,  10 | 2,   4,   11,  12 }
     *     allLocalOffsets = {  0,                   ,  5,               , 9,               , 13 }
     *  \endcode
     */
    hmemo::HArray<IndexType> allLocal2Global; 

    /**
     *   The inverse permutation to allLocal2Global.
     *
     *     allOwners       = {  0,   1,   2,   0,  2,   0,   1,   0,  0,   1,   1,   2,   2  }
     *     allGlobal2Local = {  0,   5,   9,   1,  10,  2,   6,   3,  4,   7,   8,  11,   12 }
     *
     *     So for a global index k we get the owner by allOwners[k] and the local offset there
     *     by allGlobal2Local[k] - allLocalOffsets[owner].
     */
    hmemo::HArray<IndexType> allGlobal2Local;
};

/* ------------------------------------------------------------------------ */
/* Inline methods                                                           */
/* ------------------------------------------------------------------------ */

IndexType AnyAddressing::localSize( const PartitionId rank ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( rank, allLocalOffsets.size() - 1, "illegal" )
    return allLocalOffsets[ rank + 1] - allLocalOffsets[rank];
}

PartitionId AnyAddressing::owner( const IndexType globalIndex ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( globalIndex, allOwners.size(), "illegal" )
    return allOwners[ globalIndex ];
}

IndexType AnyAddressing::localIndex( const IndexType globalIndex, const PartitionId owner ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( globalIndex, allOwners.size(), "illegal" )
    // here the owner is important as local index  requires size offsets
    return allGlobal2Local[ globalIndex ] - allLocalOffsets[ owner ];
}

IndexType AnyAddressing::globalIndex( const IndexType localIndex, const PartitionId owner ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( localIndex, localSize( owner ), "illegal" )
    return allLocal2Global[ localIndex + allLocalOffsets[ owner ] ];
}

} /* end namespace dmemo */

} /* end namespace scai */
