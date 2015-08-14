
#include <scai/common/Thread.hpp>
#include <scai/tracing.hpp>
#include <unistd.h>

using scai::common::Thread;

void subA( int& )
{
    SCAI_LOG_THREAD( "Thread1" )
    SCAI_REGION( "A" )
    sleep( 2 );
}

void subB( int& )
{
    SCAI_LOG_THREAD( "Thread2" )
    SCAI_REGION( "B" )
    sleep( 2 );
}

int main()
{
    SCAI_LOG_THREAD( "master" )
    SCAI_REGION( "main" )
    int dummy = 0;
    Thread t1( subA, dummy );
    Thread t2( subB, dummy );
}

