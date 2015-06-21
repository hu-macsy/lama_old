/**
 * @file Thread.cpp
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
 * @brief Implementation of methods for handling of thread names
 * @author Thomas Brandes
 * @date   10.06.2015
 */

#include <common/Thread.hpp>
#include <common/Exception.hpp>

#include <map>
#include <string>
#include <iostream>

using namespace std;

namespace common
{

Thread::Id Thread::getSelf()
{
    return pthread_self();
}

// Map that defines mapping thread ids -> thread names (as strings)

static map<Thread::Id, string> mapThreads;

Thread::Mutex::Mutex( bool isRecursive )
{
    pthread_mutexattr_init( &p_mutexattr);

    if ( isRecursive )
    {
        pthread_mutexattr_settype( &p_mutexattr, PTHREAD_MUTEX_RECURSIVE);
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

void Thread::Condition::notify_one()
{
    pthread_cond_signal ( &p_condition );
}

void Thread::Condition::notify_all()
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

void Thread::defineCurrentThreadId( const char* name )
{
    ScopedLock lock( map_mutex );

    Thread::Id id = getSelf();

    cout << "defineCurrentThreadId, id = " << id << ", name = " << name << endl;

    map<Thread::Id, string>::iterator it = mapThreads.find( id );

    if ( it == mapThreads.end() )
    {
        // name not defined yet

        mapThreads[ id ] = string( name );

        cout << "Thread " << id << " defines name " << name << endl;
    }
    else
    {
        // already defined, but probably on purporse

        // cout << "Redefine Thread " << id << " = " << it->second << " as " << name << endl;

        cout << "Thread " << id << " named " << it->second << ", renamed to " << name << endl;

        it->second = name;
    }
}

const char* Thread::getCurrentThreadId()
{
    Thread::ScopedLock lock( map_mutex );

    Thread::Id id = getSelf();

    map<Thread::Id, string>::iterator it = mapThreads.find( id );

    if ( it == mapThreads.end() )
    {
        // No name defined yet

        ostringstream thread_name;

        thread_name << "t_" << id;

        /* Attention: This fails when called before program start:

           mapThreads[ id ] = thread_name.str();

           it = mapThreads.find( id );

           return it->second.c_str();
        */

        return thread_name.str().c_str();
    }
    else
    {
        // return the defined name

        return it->second.c_str();
    }
}

} // namespace
