/**
 * @file Thread.cpp
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
 * @brief Implementation of methods for threads using C++11 standard
 * @author Thomas Brandes
 * @date 10.06.2015
 */

// hpp
#include <scai/common/Thread.hpp>

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

Thread::Id Thread::getSelf()
{
    return std::this_thread::get_id();
}

// Map that defines mapping thread ids -> thread names (as strings)

typedef map<Thread::Id, string> MapThreads;

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

Thread::Mutex::Mutex( bool isRecursive ) : mIsRecursive( isRecursive )
{
    // Only one of the member objects is allocated
    if ( mIsRecursive )
    {
        mRecursiveMutex.reset( new std::recursive_mutex() );
    }
    else
    {
        mMutex.reset( new std::mutex() );
    }
}

Thread::Condition::Condition() : std::condition_variable_any()
{
}

Thread::Mutex::~Mutex()
{
}

Thread::Condition::~Condition()
{
}

void Thread::Mutex::lock()
{
    if ( mIsRecursive )
    {
        mRecursiveMutex->lock();
    }
    else
    {
        mMutex->lock();
    }
}

void Thread::Mutex::unlock()
{
    if ( mIsRecursive )
    {
        mRecursiveMutex->unlock();
    }
    else
    {
        mMutex->unlock();
    }
}

void Thread::Condition::notifyOne()
{
    std::condition_variable_any::notify_one();
}

void Thread::Condition::notifyAll()
{
    std::condition_variable_any::notify_all();
}

void Thread::Condition::wait( ScopedLock& lock )
{
    if ( lock.mMutex.mIsRecursive )
    {
        std::condition_variable_any::wait( *lock.mMutex.mRecursiveMutex );
    }
    else
    {
        std::condition_variable_any::wait( *lock.mMutex.mMutex );
    }
}

Thread::ScopedLock::ScopedLock( Mutex& mutex ) :

    mMutex( mutex )
{
    mMutex.lock();
}

Thread::ScopedLock::~ScopedLock( )
{
    mMutex.unlock();
}

Thread::Mutex map_mutex; // Make access to map thread safe

void Thread::defineCurrentThreadName( const char* name )
{
    ScopedLock lock( map_mutex );
    Thread::Id id = getSelf();
#ifdef LOCAL_DEBUG
    cout << "defineCurrentThreadName, id = " << id << ", name = " << name << endl;
#endif
    MapThreads& mapThreads = getMapThreads();
    map<Thread::Id, string>::iterator it = mapThreads.find( id );

    if ( it == mapThreads.end() )
    {
        // name not defined yet
        mapThreads.insert( std::pair<Thread::Id, string>( id, name ) );
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

const char* Thread::getThreadName( Thread::Id id )
{
    Thread::ScopedLock lock( map_mutex );
    MapThreads& mapThreads = getMapThreads();
    MapThreads::iterator it = mapThreads.find( id );

    if ( it == mapThreads.end() )
    {
        // No name defined yet, give it one, use internal numbering
        // Tracing requires unique name
        ostringstream thread_name;
        thread_name << "thread_" << mapThreads.size();
        // Attention: This would not possible if mapThreads is not statically initialized
        mapThreads.insert( std::pair<Thread::Id, string>( id, thread_name.str() ) );
        it = mapThreads.find( id );
    }

    // return the defined name
    return it->second.c_str();
}

const char* Thread::getCurrentThreadName()
{
    return getThreadName( getSelf() );
}

void Thread::start( ThreadFunction start_routine, void* arg )
{
    mThread = new std::thread( start_routine, arg );
}

void Thread::join()
{
    if ( mThread != NULL )
    {
        mThread->join();
        delete mThread;
        mThread = NULL;
    }
}

Thread::~Thread()
{
    join();
}

} /* end namespace common */

} /* end namespace scai */
