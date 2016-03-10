/**
 * @file DemoThreadPool.cpp
 *
 * @brief Use of class ThreadPool
 * @author Thomas Brandes
 * @date 13.07.2015
 */

#include <scai/tasking/ThreadPool.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

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

    int n = 2;

    Settings::getEnvironment( n, "SCAI_THREADPOOL_SIZE" );

    ThreadPool pool( n );   // creates a pool with n threads

    std::cout << "ThreadPool, size = " << pool.size() << " created" << std::endl;

    double t0 = Walltime::get();

    // schedule 3 tasks, task3 must wait for completion of task1 or task2 

    shared_ptr<ThreadPoolTask> task1 = pool.schedule( bind( &work, 1 ) );
    shared_ptr<ThreadPoolTask> task2 = pool.schedule( bind( &work, 2 ) );
    shared_ptr<ThreadPoolTask> task3 = pool.schedule( bind( &work, 3 ) );

    std::cout << "Bundle1: all tasks scheduled" << std::endl;

    // wait for completion

    pool.wait( task1 );
    pool.wait( task2 );
    pool.wait( task3 );

    double t1 = Walltime::get() - t0;

    std::cout << "All tasks finished, time = " << t1 << ", expected ~4 seconds" << std::endl;
}
