/**
 * @file thread.cpp
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
 * @brief Implementation of methods for naming threads
 * @author Thomas Brandes
 * @date 10.06.2015
 */

// hpp
#include <scai/common/thread.hpp>

// local library
#include <scai/common/macros/system_call.hpp>

// std
#include <map>
#include <string>
#include <iostream>
#include <string.h>

// define LOCAL_DEBUG for debugging this source code

#undef LOCAL_DEBUG

using namespace std;

namespace scai
{

namespace common
{

namespace thread
{

// Map that defines mapping thread ids -> thread names (as strings)

typedef map<std::thread::id, string> MapThreads;

static MapThreads& getMapThreads()
{
    static MapThreads* mapThreads;

    if ( mapThreads == NULL )
    {
        // allocate it for when its used the first time
        mapThreads = new MapThreads;
    }

    return *mapThreads;
}

std::mutex map_mutex; // Make access to map thread safe

void defineCurrentThreadName( const char* name )
{
    std::unique_lock<std::mutex> lock( map_mutex );
    std::thread::id id = std::this_thread::get_id();

#ifdef LOCAL_DEBUG
    cout << "defineCurrentThreadName, id = " << id << ", name = " << name << endl;
#endif
    MapThreads& mapThreads = getMapThreads();
    map<Id, string>::iterator it = mapThreads.find( id );

    if ( it == mapThreads.end() )
    {
        // name not defined yet
        mapThreads.insert( std::pair<Id, string>( id, name ) );
#ifdef LOCAL_DEBUG
        cout << "Thread " << id << " defines name " << name << endl;
#endif
    }
    else
    {
        // already defined, but probably on purporse
        // cout << "Redefine Thread " << id << " = " << it->second << " as " << name << endl;
#ifdef LOCAL_DEBUG
        cout << "Thread " << id << " named " << it->second << ", renamed to " << name << endl;
#endif
        it->second = name;
    }
}

const char* getThreadName( Id id )
{
    std::unique_lock<std::mutex> lock( map_mutex );
    MapThreads& mapThreads = getMapThreads();
    MapThreads::iterator it = mapThreads.find( id );

    if ( it == mapThreads.end() )
    {
        // No name defined yet, give it one, use internal numbering
        // Tracing requires unique name
        ostringstream thread_name;
        thread_name << "thread_" << mapThreads.size();
        // Attention: This would not possible if mapThreads is not statically initialized
        mapThreads.insert( std::pair<Id, string>( id, thread_name.str() ) );
        it = mapThreads.find( id );
    }

    // return the defined name
    return it->second.c_str();
}

const char* getCurrentThreadName()
{
    return getThreadName( std::this_thread::get_id() );
}

} /* end namespace thread */

} /* end namespace common */

} /* end namespace scai */
