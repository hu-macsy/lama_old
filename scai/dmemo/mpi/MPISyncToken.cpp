/**
 * @file MPISyncToken.cpp
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
 * @brief Implementation of methods for class MPISyncToken.
 * @author Thomas Brandes, Jiri Kraus
 * @date 23.03.2011
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
    if ( !isSynchronized() )
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

    if ( isSynchronized() )
    {
        SCAI_LOG_WARN( logger, *this << ": waiting twice" )
        return; // do not wait twice, especially do not clean-up twice
    }

    SCAI_LOG_INFO( logger, *this << ": wait" )
    SCAI_MPICALL( logger, MPI_Waitall( mUsedRequests, mRequests.get(), mStatuses.get() ), "MPI_Waitall" );
    SCAI_LOG_INFO( logger, *this << ": synchronized, clean up and free accesses" )
    setSynchronized();
}

bool MPISyncToken::probe() const
{
    // not well defined yet
    int flag = 0;
    SCAI_MPICALL( logger, MPI_Testall( mUsedRequests, mRequests.get(), &flag, mStatuses.get() ), "MPI_Testall" );
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
