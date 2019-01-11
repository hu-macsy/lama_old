/**
 * @file VTInterface.cpp
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
 * @brief Implementation of interface between LAMA and VampirTrace
 * @author Thomas Brandes
 * @date 11.06.2015
 */

// hpp
#include <scai/tracing/VTInterface.hpp>

// local library
#include <scai/tracing/RegionEntry.hpp>

// std
#include <string>

namespace scai
{

namespace tracing
{

/** Macro that should be defined if a version of VampirTrace is available that also
 *  supports the definition of a group for a region.
 */

#undef NEW_VT

#ifdef USE_VAMPIRTRACE

#ifdef NEW_VT
extern "C" unsigned int VT_User_def__( const char* name, const char* group, const char* file, int lno );
extern "C" void VT_User_start2__( unsigned int rid );
extern "C" void VT_User_end2__( unsigned int rid );
#else
extern "C" void VT_User_start__( const char* name, const char* file, int lno );
extern "C" void VT_User_end__( const char* name );
#endif

extern "C" void VT_User_trace_on__();
extern "C" void VT_User_trace_off__();

#endif

void VTInterface::define( RegionEntry& region )
{
    region.mVTId = 0;
#if defined( USE_VAMPIRTRACE ) && defined( NEW_VT )
    const std::string& fullName = region.mName;
    //  check region name for <group_name>.<method_name>
    size_t pindex = fullName.find_first_of( "." );

    if ( pindex == std::string::npos )
    {
        // region id without "." belong to group LAMA
        region.mVTId = VT_User_def__( region.getRegionName(), "LAMA", region.getFileName(), region.mLine );
    }
    else
    {
        // split the region name and define it
        std::string groupName = fullName.substr( 0, pindex );
        std::string regionName = fullName.substr( pindex + 1 );
        region.mVTId = VT_User_def__( regionName.c_str(), groupName.c_str(), region.getFileName(), region.mLine );
    }

#endif
}

#if defined( USE_VAMPIRTRACE )

void VTInterface::enter( const RegionEntry& region )
{
#ifdef NEW_VT
    VT_User_start2__( region.mVTId );
#else
    VT_User_start__( region.getRegionName(),
                     region.getFileName(),
                     region.getLine() );
#endif
}

#else

// empty routine if VT tracing is disabled at compile time

void VTInterface::enter( const RegionEntry& )
{
}

#endif

#if defined( USE_VAMPIRTRACE )

void VTInterface::leave( const RegionEntry& region )
{
#ifdef NEW_VT
    VT_User_end2__( region.mVTId );
#else
    VT_User_end__( region.getRegionName() );
#endif
}

#else

// empty routine if VT tracing is disabled at compile time

void VTInterface::leave( const RegionEntry& )
{
}

#endif

#if defined( USE_VAMPIRTRACE )

void VTInterface::enable( const bool flag )
{
    if ( flag )
    {
        VT_User_trace_on__();
    }
    else
    {
        VT_User_trace_off__();
    }
}

#else

void VTInterface::enable( const bool )
{
}

#endif

} /* end namespace tracing */

} /* end namespace scai */
