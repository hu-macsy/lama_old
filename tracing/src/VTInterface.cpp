/**
 * @file VTInterface.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Implementation of interface between LAMA and VampirTrace
 * @author Thomas Brandes
 * @date 11.06.2015
 */

// hpp
#include "tracing/VTInterface.hpp"

// others
#include "tracing/RegionEntry.hpp"

#include <string>

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

} // namespace
