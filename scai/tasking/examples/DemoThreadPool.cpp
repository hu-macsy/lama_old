/**
 * @file DemoThreadPool.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Use of class ThreadPool
 * @author Thomas Brandes
 * @date 13.07.2015
 */

#include <scai/tasking/ThreadPool.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

#include <memory>
#include <functional>

using std::shared_ptr;
using std::bind;

using namespace scai::common;
using namespace scai::tasking;

/** For demo purpose we take a work routine that sleeps for a certain time. */

void work( const int in )
{
    Walltime::sleep( in * 1000 );
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
    shared_ptr<ThreadPoolTask> task1 = pool.schedule( std::bind( &work, 1 ) );
    shared_ptr<ThreadPoolTask> task2 = pool.schedule( std::bind( &work, 2 ) );
    shared_ptr<ThreadPoolTask> task3 = pool.schedule( std::bind( &work, 3 ) );
    std::cout << "Bundle1: all tasks scheduled" << std::endl;
    // wait for completion
    pool.wait( task1 );
    pool.wait( task2 );
    pool.wait( task3 );
    double t1 = Walltime::get() - t0;
    std::cout << "All tasks finished, time = " << t1 << ", expected ~4 seconds" << std::endl;
}
