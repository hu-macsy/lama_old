/**
 * @file KernelRegistry.cpp
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
 * @brief Implementation of routines for KernelRegistry
 * @author Thomas Brandes
 * @date 13.10.2015
 */

#include <scai/kregistry/KernelRegistry.hpp>

namespace scai
{

using common::ContextType;

namespace kregistry
{

/* -----------------------------------------------------------------------------*/

SCAI_LOG_DEF_LOGGER( KernelRegistry::logger, "KernelRegistry" )

/* -----------------------------------------------------------------------------*/

// define static variable for KernelMap here

KernelRegistry* KernelRegistry::mInstance = NULL;

/* -----------------------------------------------------------------------------*/

void KernelRegistry::registerContextFunction( const KernelRegistryKey& key, ContextType ctx, VoidFunction fn, bool replace )
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

void KernelRegistry::unregisterContextFunction( const KernelRegistryKey& key, ContextType ctx, VoidFunction fn )
{
    // this is safe at program exit as registry is still alive
    SCAI_LOG_INFO( logger, "unregister ctx = " << ctx << " with " << key )
    KernelRegistry& kreg = getInstance();
    KernelMap::iterator it = kreg.theKernelMap.find( key );

    if ( it == kreg.theKernelMap.end() )
    {
        SCAI_LOG_INFO( logger, "unregister: entry for key = " << key << " not found" )
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

        if ( it->second.isEmpty() )
        {
            kreg.theKernelMap.erase( it );
            SCAI_LOG_INFO( logger, "erased complete entry, #entries = " << kreg.theKernelMap.size() )
        }
    }
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

