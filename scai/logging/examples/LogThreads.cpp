/**
 * @file LogThreads.cpp
 * @brief Logging with multiple threads.
 */

#define SCAI_LOG_LEVEL_WARN

#include <scai/logging.hpp>

#include <cstdlib>

#include <boost/thread.hpp>

SCAI_LOG_DEF_LOGGER( myLogger, "LogTest" )

int threadRoutine( int id, int param )
{
    SCAI_LOG_THREAD( "thread_" << id )
 
    SCAI_LOG_INFO( myLogger, "starts, param = " << param )

    sleep( param );

    SCAI_LOG_INFO( myLogger, "stops, param = " << param )

    return 0;
}

int main( int, char** )
{
    // macro to give the current thread a name that appears in further logs

    SCAI_LOG_THREAD( "main" )
    
    int params[] = { 1, 2, 3, 5 };

    const int N = sizeof( params ) / sizeof( int );

    // log macros handle arguments like streams do

    SCAI_LOG_INFO( myLogger, "start " << N << " threads" )

    std::vector<boost::thread*> threads;

    threads.resize( N );

    for ( int i = 0; i < N; ++i )
    {
        SCAI_LOG_INFO( myLogger, "create thread " << i )
        threads[i] = new boost::thread( threadRoutine, i, params[i] );
    }

    SCAI_LOG_INFO( myLogger, "go sleep for 5 seconds" )

    sleep( 5 );

    SCAI_LOG_INFO( myLogger, "wait for threads" )

    for ( int i = 0; i < N; ++i )
    {
        threads[i]->join();
        delete threads[i];
    }

    SCAI_LOG_INFO( myLogger, "Threads finished, terminate" )
}
