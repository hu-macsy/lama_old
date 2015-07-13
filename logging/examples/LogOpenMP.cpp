#include "logging/logging.hpp"
#include "omp.h"

LAMA_LOG_DEF_LOGGER( myLogger, "LogOpenMP" )

int main( int argc, char** argv )
{
    LAMA_LOG_THREAD( "main" )

    // Here we name the OpenMP Threads

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();

        if ( thread_id > 0 ) 
        {
            LAMA_LOG_THREAD( "OMP_Thread_" << thread_id )
        }
    }

    const int N = 20;

    #pragma omp parallel for schedule( static, 2 )
    for ( int i = 0; i < N; ++i)
    {
        // logging per thread shows exactly which thread executes which iteration

        LAMA_LOG_INFO( myLogger, "executes iteration " << i << " of " << N )
    }
}
