
#include <tracing.hpp>
#include <unistd.h>

void subA()
{
    LAMA_REGION( "A" )
    sleep( 1 );
}

void subB()
{
    LAMA_REGION( "B" )
    subA();
    sleep( 2 );
}

int main()
{
    LAMA_LOG_THREAD( "master" )
    LAMA_REGION( "main" )
    subA();
    subB();
    sleep( 3 );
}

