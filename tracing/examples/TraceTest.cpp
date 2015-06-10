
#include <tracing/tracing.hpp>

void subA()
{
    LAMA_REGION( "A" )
    sleep(1);
}

void subB()
{
    LAMA_REGION( "B" )
    sleep(1);
    subA();
}

int main()
{
    LAMA_LOG_THREAD( "master" )

    LAMA_REGION( "main" )

    sleep(1);

    subA();
    subB();
}

