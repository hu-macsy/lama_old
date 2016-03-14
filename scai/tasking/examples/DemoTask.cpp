/**
 * @file DemoTask.cpp
 *
 * @brief Use of class Task
 * @author Thomas Brandes
 * @date 13.07.2015
 */

#include <scai/tasking/Task.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

#include <unistd.h>

using namespace scai::common;
using namespace scai::tasking;

/** For demo purpose we take a work routine that sleeps for a certain time. */

void work( const int in )
{
    sleep( in );
}

int main( int argc, const char** argv )
{
    Settings::parseArgs( argc, argv );

    // if no value specified take 2 as default

    bool replace = false;  // do not override if already set

    Settings::putEnvironment( "SCAI_THREADPOOL_SIZE", 2, replace );

    double t0 = Walltime::get();

    {
        // Note: Task use own thread pool, size specified by SCAI_THREADPOOL_SIZE

        Task t1( bind( &work, 3 ) );
        Task t2( bind( &work, 2 ) );
        Task t3( bind( &work, 1 ) );

        std::cout << "Bundle: all tasks scheduled" << std::endl;
    }

    double t1 = Walltime::get() - t0;

    std::cout << "All tasks finished, time = " << t1 << ", expected ~3 seconds" << std::endl;
}
