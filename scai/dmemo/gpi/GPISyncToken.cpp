/**
 * @file GPISyncToken.cpp
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
 * @brief GPISyncToken.cpp
 * @author Lauretta Schubert
 * @date 25.02.2014
 * @since 1.1.0
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
