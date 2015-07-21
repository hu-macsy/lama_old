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
#include <memory/ReadAccess.hpp>
#include <memory/WriteAccess.hpp>
#include <common/Exception.hpp>
#include <tasking/ThreadPool.hpp>
#include <tasking/Task.hpp>
#include <common/bind.hpp>

using namespace memory;
using namespace tasking;

/* --------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( logger, "Threading" )

using namespace memory;

void readJob( LAMAArray<double>& X )
{
    ReadAccess<double> read( X );
    const double* data = read.get();
    double s = data[0];
    LAMA_LOG_INFO( logger, "Do Read job, size = " << X.size() << ", val = " << s )
    int error = 0;

    for ( int k = 0; k < 100; ++k )
    {
        for ( int i = 1; i < X.size(); ++i )
        {
            double v = data[i];

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
    // Note: different thread on same context will wait until other access is released

    WriteAccess<double> write( X );
 
    double* data = write.get();

    LAMA_LOG_INFO( logger, "Do Write job, size = " << write.size() << ", val = " << data[0] )

    for ( int i = 0; i < write.size(); ++i )
    {
        data[i] = data[i] + 1.0;
    }
}

void job( LAMAArray<double>* X )
{
    int r = rand();
    int kind = r % 5;

    try
    {
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
    catch ( common::Exception ex )
    {
        LAMA_LOG_ERROR( logger, "job, r = " << r << ", kind = " << kind << ", caught exception: " << ex.what() )
    }
}
using namespace tasking;
int main()
{
    LAMA_LOG_THREAD( "main" )
    LAMAArray<double> X( 100000, 10 );
    LAMA_LOG_INFO( logger, "X = " << X << " at " << ( &X ) )
    tasking::ThreadPool pool( 10 );

    for ( int k = 0; k < 100; ++k )
    {
        pool.schedule( common::bind( &job, &X )  );
    }

    LAMA_LOG_INFO( logger, "synchronize" )

    // Wait on pool until all tasks are finished

    pool.shutdown();  // otherwise
}
