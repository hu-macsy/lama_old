/**
 * @file MPISyncToken.cpp
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
 * @brief Implementation of methods for class MPISyncToken.
 * @author Thomas Brandes, Jiri Kraus
 * @date 23.03.2011
 * @since 1.0.0
 */

// hpp
#include <scai/dmemo/mpi/MPISyncToken.hpp>

// local library
#include <scai/dmemo/mpi/MPIUtils.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>

using namespace scai::tasking;

namespace scai
{

namespace dmemo
{

MPISyncToken::MPISyncToken( PartitionId nRequests )
    : mNRequests( nRequests ), mUsedRequests( 0 ), mRequests( new MPI_Request[mNRequests] ), mStatuses(
          new MPI_Status[mNRequests] )
{
    SCAI_LOG_INFO( logger, "MPISyncToken for " << mNRequests << " requests constructed" )
}

MPISyncToken::~MPISyncToken()
{
    if( !isSynchronized() )
    {
        SCAI_LOG_DEBUG( logger, *this << ": synchnronized at destructor" )
        wait();
    }
}

void MPISyncToken::writeAt( std::ostream& stream ) const
{
    stream << "MPISyncToken( requests = " << mUsedRequests << ", synchronized = " << isSynchronized() << " )";
}

void MPISyncToken::wait()
{
    SCAI_REGION( "SyncToken.MPI.wait" )

    if( isSynchronized() )
    {
        SCAI_LOG_WARN( logger, *this << ": waiting twice" )

        return; // do not wait twice, especially do not clean-up twice
    }

    SCAI_LOG_INFO( logger, *this << ": wait" )

    LAMA_MPICALL( logger, MPI_Waitall( mUsedRequests, mRequests.get(), mStatuses.get() ), "MPI_Waitall" );

    SCAI_LOG_INFO( logger, *this << ": synchronized, clean up and free accesses" )

    setSynchronized();
}

bool MPISyncToken::probe() const
{
    // not well defined yet

    int flag = 0;
    LAMA_MPICALL( logger, MPI_Testall( mUsedRequests, mRequests.get(), &flag, mStatuses.get() ), "MPI_Testall" );
    return flag != 0;
}

void MPISyncToken::pushRequest( const MPI_Request& request )
{
    SCAI_ASSERT_ERROR( mUsedRequests < mNRequests,
                       "too many requests" << ", MPISyncToken allocated for maximal " << mNRequests << " requests" )

    mRequests[mUsedRequests++] = request;

    SCAI_LOG_INFO( logger, "# used requests = " << mUsedRequests )
}

} /* end namespace dmemo */

} /* end namespace scai */
