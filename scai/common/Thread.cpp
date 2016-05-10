/**
 * @file Thread.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Implementation of methods for handling of thread names
 * @author Thomas Brandes
 * @date 10.06.2015
 */

// hpp
#include <scai/common/Thread.hpp>

// local library
#include <scai/common/macros/assert.hpp>
#include <scai/common/macros/throw.hpp>

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
    return pthread_self();
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

Thread::Mutex::Mutex( bool isRecursive )
{
    pthread_mutexattr_init( &p_mutexattr );

    if ( isRecursive )
    {
        pthread_mutexattr_settype( &p_mutexattr, PTHREAD_MUTEX_RECURSIVE );
    }

    int rc = pthread_mutex_init( &p_mutex, &p_mutexattr );

    if ( rc != 0 )
    {
        COMMON_THROWEXCEPTION( "mutex init failed, rc = " << rc )
    }
}

Thread::Condition::Condition()
{
    int rc = pthread_cond_init( &p_condition, NULL );

    if ( rc != 0 )
    {
        COMMON_THROWEXCEPTION( "condition init failed, rc = " << rc )
    }
}

Thread::Mutex::~Mutex()
{
    int rc = pthread_mutex_destroy( &p_mutex );

    if ( rc != 0 )
    {
        COMMON_THROWEXCEPTION( "mutex destroy failed, rc = " << rc )
    }
}

Thread::Condition::~Condition()
{
    pthread_cond_destroy( &p_condition );
}

void Thread::Mutex::lock()
{
    pthread_mutex_lock( &p_mutex );
}

void Thread::Mutex::unlock()
{
    pthread_mutex_unlock( &p_mutex );
}

void Thread::Condition::notifyOne()
{
    pthread_cond_signal ( &p_condition );
}

void Thread::Condition::notifyAll()
{
    pthread_cond_broadcast ( &p_condition );
}

void Thread::Condition::wait( ScopedLock& lock )
{
    pthread_cond_wait( &p_condition, &lock.mMutex.p_mutex );
}

Thread::ScopedLock::ScopedLock( Mutex& mutex ) : mMutex( mutex )
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

void Thread::start( pthread_routine start_routine, void* arg )
{
    int rc = pthread_create( &mTId, NULL, start_routine, arg );

    if ( rc != 0 )
    {
        COMMON_THROWEXCEPTION( "pthread_create failed, err = " << rc << ", " << strerror( rc ) )
    }

    running = true;
}

void Thread::join()
{
    if ( running )
    {
        int rc = pthread_join( mTId, NULL );
        SCAI_ASSERT_EQUAL( 0, rc, "pthread_join failed" )

        // std::cout << "PThread " << std::hex << mTId << " terminated." << std::endl;
    }

    running = false;
}

Thread::~Thread()
{
    join();
}

} /* end namespace common */

} /* end namespace scai */
