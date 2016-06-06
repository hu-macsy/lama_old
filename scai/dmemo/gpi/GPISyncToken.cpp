/**
 * @file GPISyncToken.cpp
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
 * @brief GPISyncToken.cpp
 * @author Lauretta Schubert
 * @date 25.02.2014
 */

// hpp
#include <scai/dmemo/gpi/GPISyncToken.hpp>
#include <scai/dmemo/gpi/GPIUtils.hpp>

// others
#include <scai/tracing.hpp>
#include <GASPI.h>

#include <scai/common/macros/assert.hpp>

namespace scai
{

namespace dmemo
{

GPISyncToken::GPISyncToken( const gaspi_segment_id_t segId,
                            const PartitionId numNotifications )
    : mSegId( segId )

{
    SCAI_LOG_INFO( logger, "GPISyncToken for constructed" )
    mNotifications.reserve( numNotifications );
}

GPISyncToken::~GPISyncToken()
{
    if ( !isSynchronized() )
    {
        SCAI_LOG_DEBUG( logger, *this << ": synchnronized at destructor" )
        wait();
    }
}

void GPISyncToken::pushNotification( const PartitionId notification )
{
    mNotifications.push_back( notification );
}

void GPISyncToken::writeAt( std::ostream& stream ) const
{
    if ( isSynchronized() )
    {
        stream << "GPISyncToken( synchronized )";
    }
    else
    {
        stream << "GPISyncToken( wait for " << mNotifications.size() << " notifications )";
    }
}

void GPISyncToken::wait()
{
    SCAI_REGION( "SyncToken.GPI.wait" )

    if ( isSynchronized() )
    {
        SCAI_LOG_WARN( logger, *this << ": waiting twice" )

        return; // do not wait twice, especially do not clean-up twice
    }

    for ( size_t i = 0; i < mNotifications.size(); ++i )
    {
        gaspi_notification_id_t waitNotId = mNotifications[i];
        gaspi_notification_id_t gotNotId;
        SCAI_GASPI_CALL ( gaspi_notify_waitsome( mSegId, waitNotId, 1 , &gotNotId, GASPI_BLOCK) );
        gaspi_notification_t notifyVal;
        SCAI_GASPI_CALL( gaspi_notify_reset( mSegId, gotNotId, &notifyVal) );
    }

    mNotifications.clear();

    SCAI_LOG_INFO( logger, *this << ": synchronized, clean up and free accesses" )

    setSynchronized();
}

bool GPISyncToken::probe() const
{
    // no accurate functionality in GPI

    return false;
}

}

}
