/**
 * @file Halo.cpp
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
 * @brief Halo.cpp
 * @author Thomas Brandes
 * @date 23.02.2011
 */

// hpp
#include <scai/dmemo/Halo.hpp>

#include <scai/hmemo/HostReadAccess.hpp>

namespace scai
{

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( Halo::logger, "Halo" )

Halo::Halo()
{
}

Halo::~Halo()
{
}

void Halo::clear()
{
    mRequiredPlan.clear();
    mProvidesPlan.clear();
    mRequiredIndexes.clear();
    mProvidesIndexes.clear();
    mGlobal2Halo.clear();
}

void Halo::purge()
{
    mRequiredPlan.purge();
    mProvidesPlan.purge();
    mRequiredIndexes.purge();
    mProvidesIndexes.purge();
    // free memory of map by reallocation
    std::map<IndexType, IndexType>().swap( mGlobal2Halo );
}

/* ---------------------------------------------------------------------- */

void Halo::halo2Global( hmemo::HArray<IndexType>& indexes ) const
{
    IndexType numValues = indexes.size();
    IndexType haloSize  = mRequiredIndexes.size();

    hmemo::WriteAccess<IndexType> wIndexes( indexes );
    hmemo::ReadAccess<IndexType> halo2global( mRequiredIndexes );

    for ( IndexType i = 0; i < numValues; i++ )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( wIndexes[i], haloSize, "illegal halo index" )
        wIndexes[i] = halo2global[wIndexes[i]];
    }
}

/* ---------------------------------------------------------------------- */

void Halo::global2Halo( hmemo::HArray<IndexType>& indexes ) const
{
    IndexType numValues = indexes.size();

    hmemo::WriteAccess<IndexType> wIndexes( indexes );

    for ( IndexType i = 0; i < numValues; i++ )
    {
        const std::map<IndexType, IndexType>::const_iterator elem = mGlobal2Halo.find( wIndexes[i] );
        SCAI_ASSERT_ERROR( elem != mGlobal2Halo.end(), "global Index " << wIndexes[i] << " no halo index, never required" )
        wIndexes[i] = elem->second;
    }
}

/* ---------------------------------------------------------------------- */

void Halo::writeAt( std::ostream& stream ) const
{
    // write info this object
    stream << "Halo( size = " << getHaloSize() << ", required plan = " << mRequiredPlan << ", provides plan = "
           << mProvidesPlan << ")";
}

} /* end namespace dmemo */

} /* end namespace scai */
