/**
 * @file hmemo/examples/Threading.cpp
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
 * @brief Demo on aliasing problem with LAMA arrays
 * @author Thomas Brandes
 * @date 03.07.2015
 */

#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/ReadAccess.hpp>
#include <scai/hmemo/WriteAccess.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/tasking/ThreadPool.hpp>
#include <scai/tasking/Task.hpp>

#include <functional>

using namespace scai::hmemo;
using namespace scai::tasking;

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Threading" )

using namespace scai;
using namespace hmemo;

void readJob( HArray<double>& X )
{
    ContextPtr contextPtr = Context::getHostPtr();
    ReadAccess<double> read( X, contextPtr );
    const double* data = read.get();
    double s = data[0];
    SCAI_LOG_INFO( logger, "Do Read job, size = " << X.size() << ", val = " << s )
    int error = 0;

    for ( int k = 0; k < 100; ++k )
    {
        for ( IndexType i = 1; i < X.size(); ++i )
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

void writeJob( HArray<double>& X )
{
    // Note: different thread on same context will wait until other access is released
    ContextPtr contextPtr = Context::getHostPtr();
    WriteAccess<double> write( X, contextPtr );
    double* data = write.get();
    SCAI_LOG_INFO( logger, "Do Write job, size = " << write.size() << ", val = " << data[0] )

    for ( IndexType i = 0; i < write.size(); ++i )
    {
        data[i] = data[i] + 1.0;
    }
}

void job( HArray<double>* X )
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
    catch ( scai::common::Exception& ex )
    {
        SCAI_LOG_ERROR( logger, "job, r = " << r << ", kind = " << kind << ", caught exception: " << ex.what() )
    }
}
using namespace scai::tasking;
int main()
{
    SCAI_LOG_THREAD( "main" )
    HArray<double> X( 100000, 10 );
    SCAI_LOG_INFO( logger, "X = " << X << " at " << ( &X ) )
    scai::tasking::ThreadPool pool( 10 );

    for ( int k = 0; k < 100; ++k )
    {
        pool.schedule( std::bind( &job, &X )  );
    }

    SCAI_LOG_INFO( logger, "synchronize" )
    // Wait on pool until all tasks are finished
    pool.shutdown();  // otherwise
}
