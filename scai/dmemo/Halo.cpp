/**
 * @file Halo.cpp
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
 * @brief Halo.cpp
 * @author Thomas Brandes
 * @date 23.02.2011
 */

// hpp
#include <scai/dmemo/Halo.hpp>

namespace scai
{

namespace dmemo
{

SCAI_LOG_DEF_LOGGER( Halo::logger, "Halo" )

Halo::Halo()
{
}

Halo::Halo( const Halo& other )
{
    operator=( other );
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

Halo& Halo::operator=( const Halo& other )
{
    if ( this != &other )
    {
        SCAI_LOG_DEBUG( logger, "make deep copy of Halo" )
        mRequiredPlan = other.mRequiredPlan;
        mProvidesPlan = other.mProvidesPlan;
        mRequiredIndexes = other.mRequiredIndexes;
        mProvidesIndexes = other.mProvidesIndexes;
        mGlobal2Halo = other.mGlobal2Halo;
    }

    return *this;
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
