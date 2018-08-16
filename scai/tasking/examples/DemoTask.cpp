/**
 * @file DemoTask.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Use of class Task
 * @author Thomas Brandes
 * @date 13.07.2015
 */

#include <scai/tasking/Task.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/Settings.hpp>

#include <functional>

using namespace scai::tasking;
using scai::common::Walltime;
using scai::common::Settings;
using std::bind;

/** For demo purpose we take a work routine that sleeps for a certain time. */

void work( const int in )
{
    Walltime::sleep( in * 1000 );
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
