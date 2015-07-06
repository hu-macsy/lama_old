/**
 * @file Aliasing.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Demo on aliasing problem with LAMA arrays
 * @author: Thomas Brandes
 * @date 03.07.2015
 **/

#include <memory/Context.hpp>
#include <memory/HostReadAccess.hpp>
#include <memory/HostWriteAccess.hpp>
#include <common/Exception.hpp>
#include <tasking/ThreadPool.hpp>
#include <tasking/Task.hpp>
#include <boost/bind.hpp>

using namespace memory;
using namespace tasking;

/* --------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( logger, "Threading" )

using namespace memory;


void readJob( LAMAArray<double>& X )
{
    HostReadAccess<double> read( X );

    double s = read[0];

    LAMA_LOG_INFO( logger, "Do Read job, size = " << read.size() << ", val = " << s )

    int error = 0;

    for ( int k = 0; k < 100; ++k )
    {
        for ( int i = 1; i < read.size(); ++i )
        {
            double v = read[i];

            if ( v !=  s )
            {
                if ( error < 5 )
                {
                    LAMA_LOG_ERROR( logger, "Read error at i = " << i << ", is " << v << ", expected " << s )
                }
                ++error;
            }
        }
    }
}


void writeJob( LAMAArray<double>& X )
{
    HostWriteAccess<double> write( X );

    LAMA_LOG_INFO( logger, "Do Write job, size = " << write.size() << ", val = " << write[0] )

    for ( int i = 0; i < write.size(); ++i )
    {
        write[i] = write[i] + 1.0;
    }
}


void job( LAMAArray<double>* X )
{
    int r = rand();

    int kind = r % 2;

    LAMA_LOG_INFO( logger, "job, r = " << r << ", kind = " << kind << ", X = " << X )

    if ( kind > 0 )
    {
        readJob( *X );
    }
    else
    {
        writeJob( *X );
    }

    LAMA_LOG_INFO( logger, "job, r = " << r << ", kind = " << kind << ", finished" )
}

using namespace tasking;

int main()
{
    std::vector<boost::shared_ptr<tasking::Task> > tasks;

    LAMAArray<double> X( 100000, 10 );

    LAMA_LOG_INFO( logger, "X = " << X << " at " << ( &X ) )

    tasking::ThreadPool pool( 4 );

    for ( int k = 0; k < 16; ++k )
    {
        boost::shared_ptr<tasking::Task> task( new Task( boost::bind( &job, &X )  ) );
        tasks.push_back( task );
    }

    LAMA_LOG_INFO( logger, "synchronize" )

    for ( int k = 0; k < 16; ++k )
    {
        tasks[k]->synchronize();
    }
    // pool.wait()  // otherwise
}

