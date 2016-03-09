/**
 * @file GPISyncToken.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Definition of derived class GPISyncToken that defines a SyncToken for GPI
 * @author Thomas Brandes    
 * @date 25.02.2014
 * @since 1.1.0
 */

#pragma once

#include <GASPI.h>

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/tasking/SyncToken.hpp>

// others
#include <scai/common/SCAITypes.hpp>
#include <scai/dmemo/gpi/GPICommunicator.hpp>

// boost
#include <vector>

namespace scai
{

namespace dmemo
{

/** Class for GPI synchronization that waits for a certain number of
 *  notifications at a given GASPI segment.
 */

class COMMON_DLL_IMPORTEXPORT GPISyncToken: public tasking::SyncToken
{
public:

    // TODO: doxy docu
    /** Constructor for an GPI synchronization token.
     *
     *  @param segId 
     *  @param numNotifications 
     *
     *  Each pending send and receive will have its own request.
     */

    GPISyncToken( const gaspi_segment_id_t segId, PartitionId numNotifications );

    /** Destructor, will also wait for synchronization and cleanup. */

    virtual ~GPISyncToken();

    /** This method adds a notification id that token will wait for. */

    void pushNotification( const PartitionId notification );

    /** This method waits for all requests. */

    virtual void wait();

    /** Method to check wheter communications have already been finished. */

    virtual bool probe() const;

    virtual void writeAt( std::ostream& stream ) const;

private:

    gaspi_segment_id_t mSegId;

    GPISyncToken();

    std::vector<PartitionId> mNotifications;
};

}  // namespace dmemo

}  // namespace scai
