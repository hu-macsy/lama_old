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

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/common/exception/Exception.hpp>
#include <scai/tasking/ThreadPool.hpp>
#include <scai/tasking/Task.hpp>
#include <scai/common/bind.hpp>

using namespace scai::hmemo;
using namespace scai::tasking;

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Threading" )

using namespace scai::hmemo;

void readJob( LAMAArray<double>& X )
{
    ContextPtr contextPtr = Context::getHostPtr();

    ReadAccess<double> read( X, contextPtr );
    const double* data = read.get();
    double s = data[0];
    SCAI_LOG_INFO( logger, "Do Read job, size = " << X.size() << ", val = " << s )
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
                    SCAI_LOG_ERROR( logger, "Read error at i = " << i << ", is " << v << ", expected " << s )
                }

                ++error;
            }
        }
    }
}

void writeJob( LAMAArray<double>& X )
{
    // Note: different thread on same context will wait until other access is released

    ContextPtr contextPtr = Context::getHostPtr();

    WriteAccess<double> write( X, contextPtr );
 
    double* data = write.get();

    SCAI_LOG_INFO( logger, "Do Write job, size = " << write.size() << ", val = " << data[0] )

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
        SCAI_LOG_INFO( logger, "job, r = " << r << ", kind = " << kind << ", X = " << X )

        if ( kind > 0 )
        {
            readJob( *X );
        }
        else
        {
            writeJob( *X );
        }

        SCAI_LOG_INFO( logger, "job, r = " << r << ", kind = " << kind << ", finished" )
    }
    catch ( scai::common::Exception ex )
    {
        SCAI_LOG_ERROR( logger, "job, r = " << r << ", kind = " << kind << ", caught exception: " << ex.what() )
    }
}
using namespace scai::tasking;
int main()
{
    SCAI_LOG_THREAD( "main" )
    LAMAArray<double> X( 100000, 10 );
    SCAI_LOG_INFO( logger, "X = " << X << " at " << ( &X ) )
    scai::tasking::ThreadPool pool( 10 );

    for ( int k = 0; k < 100; ++k )
    {
        pool.schedule( scai::common::bind( &job, &X )  );
    }

    SCAI_LOG_INFO( logger, "synchronize" )

    // Wait on pool until all tasks are finished

    pool.shutdown();  // otherwise
}
