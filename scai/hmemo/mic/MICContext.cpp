/**
 * @file MICContext.cpp
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
 * @brief Contains the implementation of the class MICContext.
 * @author Thomas Brandes
 * @date 15.07.2011
 * @since 1.1.0
 */

// hpp
#include <scai/hmemo/mic/MICContext.hpp>

// local library
#include <scai/hmemo/mic/MICMemory.hpp>
#include <scai/hmemo/mic/MICSyncToken.hpp>

// internal scai libraries
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/common/macros/assert.hpp>

// external
#include <omp.h>
#include <mkl.h>

#define SCAI_MAX_MIC_DEVICES 4
#define SCAI_DEFAULT_DEVICE_NUMBER -1

namespace scai
{

namespace hmemo
{

/**  static variables *****************************************************/

SCAI_LOG_DEF_LOGGER( MICContext::logger, "MIC.Context" )

int MICContext::currentDeviceNr = -1;

int MICContext::numUsedDevices = 0;

/**  constructor  *********************************************************/


MICContext::MICContext( int deviceNr )
    : Context( common::context::MIC ), mDeviceNr( deviceNr )
{
    SCAI_LOG_INFO( logger, "construct MICContext, device nr = = " << deviceNr )

    int numDevices = 0;

    SCAI_LOG_INFO( logger, "Checking for Intel(R) MIC Architecture (Target CPU) devices" )

#ifdef __INTEL_OFFLOAD
    numDevices = _Offload_number_of_devices();
#endif

    if( numDevices < 1 )
    {
        COMMON_THROWEXCEPTION( "No mic devices available" )
    }

    // ToDo: allow for any device

    mDeviceNr = deviceNr;

    SCAI_ASSERT_LT( deviceNr, numDevices, "Illegal deviceNr" )

    SCAI_LOG_INFO( logger, "Using device " << mDeviceNr << " of " << numDevices << " installed devices" )

    bool targetOK = false;

    int numCores;

    #pragma offload target( mic: mDeviceNr ) out( numCores )
    {
        #pragma omp parallel
        {
            #pragma omp single
            {
                numCores = omp_get_num_threads();
            }
        }
    }

    SCAI_LOG_INFO( logger, "Uses " << numCores << " threads for parallel execution" );

    mNumThreads = numCores;
}

/**  destructor   *********************************************************/

MICContext::~MICContext()
{
    SCAI_LOG_INFO( logger, "~MICContext: " << *this )
}

/* ----------------------------------------------------------------------------- */

void MICContext::writeAt( std::ostream& stream ) const
{
    stream << "MICContext(" << mDeviceNr << ": " << mNumThreads << " Threads)";
}

/* ----------------------------------------------------------------------------- */

void MICContext::disable( const char* file, int line ) const
{
    currentDeviceNr = -1;
}

/* ----------------------------------------------------------------------------- */

void MICContext::enable( const char* file, int line ) const
{
    SCAI_LOG_DEBUG( logger, *this << ": enable" )
    currentDeviceNr = mDeviceNr;
}

/* ----------------------------------------------------------------------------- */

bool MICContext::canUseMemory( const Memory& other ) const
{
    bool canUse = false;

    // CUDA device can use only data on same CUDA device

    if ( other.getType() == memtype::MICMemory )
    {
        const MICMemory* otherMICMem = dynamic_cast<const MICMemory*>( &other );

        SCAI_ASSERT( otherMICMem, "serious type mismatch" )

        canUse = otherMICMem->getDeviceNr() == mDeviceNr;
    }

    SCAI_LOG_DEBUG( logger, *this << ": " << ( canUse ? "can use " : "can't use " )
                            << other )

    return canUse;
}

/* ----------------------------------------------------------------------------- */

MemoryPtr MICContext::getLocalMemoryPtr() const
{
    MemoryPtr memory;

    if ( mMemory.expired() )
    {
        memory.reset( new MICMemory( shared_from_this() ) );
        mMemory = memory; // save it here as a weak pointer to avoid cycles
    }
    else
    {
        // the last memory instance is still valid, so we return just shared pointer to it
        memory = mMemory.lock();
    }

    return memory;
}

/* ----------------------------------------------------------------------------- */

tasking::SyncToken* MICContext::getSyncToken() const
{
    return new tasking::MICSyncToken( shared_from_this() );
}

/* ----------------------------------------------------------------------------- */

int MICContext::getCurrentDevice()
{
    // ToDo: get current device, make sure that access has been enabled

    SCAI_ASSERT_LE( 0, currentDeviceNr, "MICContext not enabled" )

    return currentDeviceNr;
}

/* ----------------------------------------------------------------------------- */

static int getDefaultDeviceNr()
{
    return 0;
}

/* ----------------------------------------------------------------------------- */
/*      Factory::Register - create( int )                                        */
/* ----------------------------------------------------------------------------- */

static common::weak_ptr<MICContext> mMICContext[SCAI_MAX_MIC_DEVICES];

ContextPtr MICContext::create( int deviceNr )
{
    int micDeviceNr = deviceNr;

    if( micDeviceNr == SCAI_DEFAULT_DEVICE_NUMBER )
    {
        micDeviceNr = getDefaultDeviceNr();

        // no need here to check for a good value
    }
    else
    {
        SCAI_ASSERT(
            0 <= micDeviceNr && micDeviceNr < SCAI_MAX_MIC_DEVICES,
            "device = " << micDeviceNr << " out of range" << ", max supported device = " << SCAI_MAX_MIC_DEVICES )
    }

    common::shared_ptr<MICContext> context = common::shared_ptr<MICContext>();

    if( mMICContext[micDeviceNr].expired() )
    {
        // create a new context for the device and return the shared pointer

        context = common::shared_ptr<MICContext>( new MICContext( micDeviceNr ) );

        // we keep a weak pointer so that we can return

        mMICContext[micDeviceNr] = context;
    }
    else
    {
        // the weak pointer to the device is still okay, so return a shared pointer for it

        context = mMICContext[micDeviceNr].lock();
    }

    return context;
}

} /* end namespace lama */

} /* end namespace scai */
