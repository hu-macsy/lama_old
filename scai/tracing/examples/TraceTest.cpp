
#include <scai/tracing.hpp>
#include <unistd.h>

void subA()
{
    SCAI_REGION( "A" )
    sleep( 1 );
}

void subB()
{
    SCAI_REGION( "B" )
    subA();
    sleep( 2 );
}

int main()
{
    SCAI_LOG_THREAD( "master" )
    SCAI_REGION( "main" )
    subA();
    subB();
    sleep( 3 );
}

