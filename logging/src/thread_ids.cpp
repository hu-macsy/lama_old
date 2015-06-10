/**
 * @file thread_ids.cpp
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

#include <logging/thread_ids.hpp>

#include <boost/thread.hpp>
#include <map>

using namespace boost;
using namespace std;

namespace log4lama
{

typedef thread::id thread_id;

// Map that defines mapping thread ids -> thread names (as strings) 

static map<thread_id, string> mapThreads;

mutex map_mutex; // Make access to map thread safe

void defineCurrentThreadId( const char* name )
{
    mutex::scoped_lock scoped_lock( map_mutex ); 

    thread_id id = this_thread::get_id();
    
    map<thread_id, string>::iterator it = mapThreads.find( id );

    if ( it == mapThreads.end() )
    {
        // name not defined yet

        mapThreads[ id ] = string( name );
    }
    else
    {
        // already defined, but probably on purporse

        // cout << "Redefine Thread " << id << " = " << it->second << " as " << name << endl;

        it->second = name;
    }
}

const char* getCurrentThreadId()
{
    mutex::scoped_lock scoped_lock( map_mutex ); 

    thread_id id = this_thread::get_id();

    map<thread_id, string>::const_iterator it = mapThreads.find( id );

    if ( it == mapThreads.end() )
    {
        // No name define yet

        return "<unk_thread>";

    }
    else
    {
        // return the defined name

        return it->second.c_str();
    }
}

} // namespace
