/**
 * @file LogThreads.hpp
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
 * @brief Logging with multiple threads.
 *
 * @author Thomas Brandes
 * @date 08.07.2015
 */

#include "logging.hpp"

#include <cstdlib>

#include <boost/thread.hpp>

LAMA_LOG_DEF_LOGGER( myLogger, "LogTest" )

int threadRoutine( int id, int param )
{
    LAMA_LOG_THREAD( "thread_" << id )
 
    LAMA_LOG_INFO( myLogger, "starts, param = " << param )

    sleep( param );

    LAMA_LOG_INFO( myLogger, "stops, param = " << param )

    return 0;
}

int main( int, char** )
{
    // macro to give the current thread a name that appears in further logs

    LAMA_LOG_THREAD( "main" )
    
    int params[] = { 1, 2, 3, 5 };

    const int N = sizeof( params ) / sizeof( int );

    // log macros handle arguments like streams do

    LAMA_LOG_INFO( myLogger, "start " << N << " threads" )

    std::vector<boost::thread*> threads;

    threads.resize( N );

    for ( int i = 0; i < N; ++i )
    {
        LAMA_LOG_INFO( myLogger, "create thread " << i )
        threads[i] = new boost::thread( threadRoutine, i, params[i] );
    }

    LAMA_LOG_INFO( myLogger, "go sleep for 5 seconds" )

    sleep( 4 );

    LAMA_LOG_INFO( myLogger, "wait for threads" )

    for ( int i = 0; i < N; ++i )
    {
        threads[i]->join();
        delete threads[i];
    }

    LAMA_LOG_INFO( myLogger, "Threads finished, terminate" )
}
