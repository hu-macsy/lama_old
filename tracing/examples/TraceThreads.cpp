
#include <common/Thread.hpp>
#include <tracing/tracing.hpp>
#include <unistd.h>

using common::Thread;

void subA( int& )
{
    LAMA_LOG_THREAD( "Thread1" )
    LAMA_REGION( "A" )
    sleep( 2 );
}

void subB( int& )
{
    LAMA_LOG_THREAD( "Thread2" )
    LAMA_REGION( "B" )
    sleep( 2 );
}

int main()
{
    LAMA_LOG_THREAD( "master" )
    LAMA_REGION( "main" )
    int dummy = 0;
    Thread t1( subA, dummy );
    Thread t2( subB, dummy );
}

