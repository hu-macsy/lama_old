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

#include <scai/kregistry/KernelRegistry.hpp>

namespace scai
{

using common::context;

namespace kregistry
{

/* -----------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( KernelRegistry::logger, "KernelRegistry" )

/* -----------------------------------------------------------------------------*/

// define static variable for KernelMap here

//KernelRegistry::KernelMap KernelRegistry::theKernelMap;
KernelRegistry* KernelRegistry::mInstance;

/* -----------------------------------------------------------------------------*/

void KernelRegistry::registerContextFunction( const KernelRegistryKey& key, context::ContextType ctx, VoidFunction fn, bool replace )
{
    SCAI_LOG_INFO( logger, "register ctx = " << ctx << " with " << key )

	KernelRegistry& kreg = getInstance();

    KernelMap::iterator it = kreg.theKernelMap.find( key );

    if ( it == kreg.theKernelMap.end() )
    {
        SCAI_LOG_DEBUG( logger, "register: no entry yet, will add it" )

        _ContextFunction routine;

        routine.set( ctx, fn );

        kreg.theKernelMap.insert( std::pair<KernelRegistryKey, _ContextFunction>( key, routine ) );

        SCAI_LOG_DEBUG( logger, "added" )
    }
    else
    {
        VoidFunction old_fn = it->second.get( ctx );

        if ( old_fn != NULL )
        {
            if ( fn == old_fn )
            {
                SCAI_LOG_WARN( logger, "same function registered again, " << key << ", context = " << ctx )
            }
            else if ( replace )
            {
                SCAI_LOG_INFO( logger, "kernel function replaced, " << key << ", context = " << ctx )
            }
            else 
            {
                SCAI_LOG_INFO( logger, "kernel function not replaced, " << key << ", context = " << ctx )
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

void KernelRegistry::unregisterContextFunction( const KernelRegistryKey& key, context::ContextType ctx, VoidFunction fn )
{
    // can cause serious segmentation violoation in case of unregister after destructo of KernelMap

    return;

    SCAI_LOG_INFO( logger, "unregister ctx = " << ctx << " with " << key )

	KernelRegistry& kreg = getInstance();

    KernelMap::iterator it = kreg.theKernelMap.find( key );

    if ( it == kreg.theKernelMap.end() )
    {
//        SCAI_LOG_ERROR( logger, "unregister: no entry for key = " << key )
    	SCAI_LOG_INFO( logger, "unregister: entry for key = " << key << " not found")
    }
    else
    {
        VoidFunction old_fn = it->second.get( ctx );

        if ( old_fn == NULL )
        {
            // rather likely to happen for multiple versions of kernel routines

            SCAI_LOG_INFO( logger, "unregister: entry for key = " << key << " available, but not for " << ctx );
        }
        else if ( old_fn != fn )
        {
            // rather likely to happen if kernel routines have been replaced

            SCAI_LOG_INFO( logger, "unregister: entry for key = " << key << " and ctx = " << ctx << " available, but does not match" )
        }
        else
        {
        	SCAI_LOG_INFO( logger, "unregister: setting ctx to NULL" )
            it->second.set( ctx, NULL );
        }

        if( it->second.isEmpty() )
        {
        	SCAI_LOG_INFO( logger, "erasing complete entry" )
        	kreg.theKernelMap.erase( it );
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

/* -----------------------------------------------------------------------------*/

void KernelRegistry::printAll()
{
	KernelRegistry& kreg = getInstance();

    KernelMap::const_iterator it;

    std::cout << "KernelRegistry:" << std::endl;
    std::cout << "================" << std::endl;

    for ( it = kreg.theKernelMap.begin(); it != kreg.theKernelMap.end(); ++it )
    {
        std::cout << "Entry: key = " << it->first;
        std::cout << ", ctx = " << it->second.printIt();
        std::cout << std::endl;
    }

    std::cout << "================" << std::endl;
}

} /* end namespace kregistry */

} /* end namespace scai */

