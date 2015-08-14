
#include <scai/common/Thread.hpp>
#include <scai/tracing.hpp>
#include <unistd.h>

using common::Thread;

void subA( int& )
{
    SCAI_LOG_THREAD( "Thread1" )
    LAMA_REGION( "A" )
    sleep( 2 );
}

void subB( int& )
{
    SCAI_LOG_THREAD( "Thread2" )
    LAMA_REGION( "B" )
    sleep( 2 );
}

int main()
{
    SCAI_LOG_THREAD( "master" )
    LAMA_REGION( "main" )
    int dummy = 0;
    Thread t1( subA, dummy );
    Thread t2( subB, dummy );
}

