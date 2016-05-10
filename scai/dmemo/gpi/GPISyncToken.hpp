/**
 * @file GPISyncToken.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Definition of derived class GPISyncToken that defines a SyncToken for GPI
 * @author Thomas Brandes    
 * @date 25.02.2014
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
