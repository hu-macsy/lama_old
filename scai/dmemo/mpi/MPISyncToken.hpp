/**
 * @file MPISyncToken.hpp
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
 * @brief MPISyncToken.hpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 23.03.2011
 */
#pragma once

#include <mpi.h> //Intel MPI need mpi.h to be included before stdio.h so this header comes first

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/tasking/SyncToken.hpp>

//Intel MPI need mpi.h to be included before stdio.h so this header comes first
#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace dmemo
{

/** Class for MPI synchronization that waits on pending messages. */

class COMMON_DLL_IMPORTEXPORT MPISyncToken: public tasking::SyncToken
{
public:

    /** Constructor for an MPI synchronization token.
     *
     *  \param nRequests is the number of requests to be allocated
     *
     *  Each pending send and receive will have its own request.
     */

    MPISyncToken( PartitionId nRequests );

    /** Destructor, will also wait for synchronization and cleanup. */

    virtual ~MPISyncToken();

    void pushRequest( const MPI_Request& request );

    /** This method waits for all requests. */

    virtual void wait();

    /** Method to check wheter communications have already been finished. */

    virtual bool probe() const;

    virtual void writeAt( std::ostream& stream ) const;

private:

    MPISyncToken();

    PartitionId mNRequests; // allocated size
    PartitionId mUsedRequests; // used size

    std::unique_ptr<MPI_Request[]> mRequests;
    std::unique_ptr<MPI_Status[]> mStatuses;
};

} /* end namespace dmemo */

} /* end namespace scai */
