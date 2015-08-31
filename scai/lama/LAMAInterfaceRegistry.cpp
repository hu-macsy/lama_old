/**
 * @file LAMAInterfaceRegistry.cpp
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
 * @brief LAMAInterfaceRegistry.cpp
 * @author Thomas Brandes
 * @date 12.04.2013
 * @since 1.0.0
 */

// hpp
#include <scai/lama/LAMAInterfaceRegistry.hpp>

// others
#include <scai/common/exception/Exception.hpp>

using namespace std;

namespace scai
{

namespace lama
{

LAMAInterfaceRegistry* LAMAInterfaceRegistry::instance = 0;

LAMAInterfaceRegistry::LAMAInterfaceRegistry()
{
}

LAMAInterfaceRegistry::~LAMAInterfaceRegistry()
{
    // free all allocated interfaces

    while( !mInterfaceMap.empty() )
    {
        InterfaceMapType::iterator begin = mInterfaceMap.begin();
        LAMAInterface* ptr = begin->second;
        mInterfaceMap.erase( begin );
        delete ptr;
    }
}

const LAMAInterface* LAMAInterfaceRegistry::getInterface( const hmemo::ContextType location ) const
{
    InterfaceMapType::const_iterator loc = mInterfaceMap.find( location );

    if( loc == mInterfaceMap.end() )
    {
        COMMON_THROWEXCEPTION( "No interface on location " << location << " available." )
    }

    return loc->second;
}

LAMAInterface& LAMAInterfaceRegistry::modifyInterface( const hmemo::ContextType location )
{
    InterfaceMapType::const_iterator loc = mInterfaceMap.find( location );

    if( loc == mInterfaceMap.end() )
    {
        // create a new default interface that can be filled for this context type

        LAMAInterface* lamaInterface = new LAMAInterface();
        mInterfaceMap[location] = lamaInterface;
        return *lamaInterface;
    }

    return *loc->second;
}

bool LAMAInterfaceRegistry::hasInterface( const hmemo::ContextType location ) const
{
    bool hasInterface = false;

    InterfaceMapType::const_iterator loc = mInterfaceMap.find( location );

    if( loc != mInterfaceMap.end() )
    {
        hasInterface = true;
    }

    return hasInterface;
}

LAMAInterfaceRegistry& LAMAInterfaceRegistry::getRegistry()
{
    static CGuard g;

    if ( instance == 0 )
    {
        instance = new LAMAInterfaceRegistry();
    }

    return *instance;
}

LAMAInterfaceRegistry::CGuard::CGuard()
{
}

LAMAInterfaceRegistry::CGuard::~CGuard()
{
    if( LAMAInterfaceRegistry::instance != 0 )
    {
        delete LAMAInterfaceRegistry::instance;
        LAMAInterfaceRegistry::instance = 0;
    }
}

}

/* -----------------------------------------------------------------------------*/

namespace hmemo
{

using namespace lama;

const LAMAInterface& Context::getInterface() const
{
    const LAMAInterface* lamaInterface = LAMAInterfaceRegistry::getRegistry().getInterface( getType() );

    // Registry throws an exception if no interface is available

    SCAI_ASSERT( lamaInterface, "No lama interface available on " << *this )

    return *lamaInterface;
}

} /* end namespace lama */

} /* end namespace scai */
