/**
 * @file KernelRegistry.cpp
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
 * @brief Implementation of routines for KernelRegistry
 * @author Thomas Brandes
 * @date 13.10.2015
 */

#include "KernelRegistry.hpp"

namespace scai
{

using common::ContextType;

namespace kregistry
{

/* -----------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( KernelRegistry::logger, "KernelRegistry" )

/* -----------------------------------------------------------------------------*/

// define static variable for KernelMap here

KernelRegistry::KernelMap KernelRegistry::theKernelMap;

/* -----------------------------------------------------------------------------*/

void KernelRegistry::registerContextFunction( const KernelRegistryKey& key, ContextType ctx, VoidFunction fn, bool replace )
{
    SCAI_LOG_INFO( logger, "register ctx = " << ctx << " with " << key )

    KernelMap::iterator it = theKernelMap.find( key );

    if ( it == theKernelMap.end() )
    {
        SCAI_LOG_DEBUG( logger, "register: no entry yet, will add it" )

        _ContextFunction routine;

        routine.set( ctx, fn );

        theKernelMap.insert( std::pair<KernelRegistryKey, _ContextFunction>( key, routine ) );

        SCAI_LOG_DEBUG( logger, "added" )
    }
    else
    {
        VoidFunction old_fn = it->second.get( ctx );

        if ( old_fn != NULL )
        {
            if ( fn == old_fn )
            {
                SCAI_LOG_WARN( logger, "same function registered again, " << key )
            }
            else if ( replace )
            {
                SCAI_LOG_INFO( logger, "kernel function replaced, " << key )
            }
        }

        SCAI_LOG_DEBUG( logger, "register: entry available, set it for ctx = " << ctx )

        if ( old_fn == NULL || replace )
        {
            it->second.set( ctx, fn );
        }
    }
}

/* -----------------------------------------------------------------------------*/

bool KernelRegistry::Compare::operator()( const KernelRegistryKey& x, const KernelRegistryKey& y )
{
    // first compare the id of the routine (is second key argument)

    int compareName = std::strcmp( x.second, y.second );

    if ( compareName < 0 )
    {
         return true;
    }

    if ( compareName > 0 )
    {
         return false;
    }

    // both have same id, so take typename to distinguish

    return x.first.name() > y.first.name();
}

} /* end namespace kregistry */

} /* end namespace scai */

