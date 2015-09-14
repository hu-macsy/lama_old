/**
 * @file MPISyncToken.hpp
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
 * @brief MPISyncToken.hpp
 * @author Thomas Brandes, Jiri Kraus
 * @date 23.03.2011
 */
#pragma once

#include <mpi.h>

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/tasking/SyncToken.hpp>

// others
#include <scai/common/SCAITypes.hpp>

#include <scai/common/unique_ptr.hpp>

namespace scai
{

namespace lama
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

    common::scoped_array<MPI_Request> mRequests;
    common::scoped_array<MPI_Status> mStatuses;
};

} /* end namespace lama */

} /* end namespace scai */
