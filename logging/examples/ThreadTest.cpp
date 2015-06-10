
#include "logging/logging.hpp"

#include <cstdlib>

#include <boost/thread.hpp>

LAMA_LOG_DEF_LOGGER( myLogger, "LogTest" )

int threadRoutine( int param )
{
    LAMA_LOG_THREAD( "thread_" << param )
 
    LAMA_LOG_INFO( myLogger, "starts" )

    sleep( param );

    LAMA_LOG_INFO( myLogger, "stops" )
}

int main( int argc, char** argv )
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
        threads[i] = new boost::thread( threadRoutine, params[i] );
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
