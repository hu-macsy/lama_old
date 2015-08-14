/**
 * @file MICContextManager.cpp
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
 * @brief Contains the implementation of the singleton class MICContextManager.
 * @author Thomas Brandes
 * @date 15.07.2011
 * @since 1.1.0
 */

// hpp
#include <scai/lama/mic/MICContextManager.hpp>

// others
#include <scai/lama/mic/MICContext.hpp>

#include <scai/lama/ContextFactory.hpp>

// assert
#include <scai/lama/exception/LAMAAssert.hpp>

namespace scai
{

namespace lama
{

/* ----------------------------------------------------------------------------- */

// definition of array with weak pointers so that we can return shared pointers without allocating again
boost::weak_ptr<MICContext> MICContextManager::mMICContext[LAMA_MAX_MIC_DEVICES];

MICContextManager MICContextManager::theInstance;

int MICContextManager::defaultDeviceNr = LAMA_DEFAULT_DEVICE_NUMBER;

/* ----------------------------------------------------------------------------- */

int MICContextManager::getDefaultDeviceNr()
{
    if( defaultDeviceNr == LAMA_DEFAULT_DEVICE_NUMBER )
    {

        // not yet set, so do it now exactly once

        if( getenv( LAMA_MIC_ENV_FOR_DEVICE ) )
        {
            std::string devNumber( getenv( LAMA_MIC_ENV_FOR_DEVICE ) );
            std::istringstream devNumberReader( devNumber );
            devNumberReader >> defaultDeviceNr;

            SCAI_LOG_INFO( logger, LAMA_MIC_ENV_FOR_DEVICE << " = " << defaultDeviceNr << " set, take it" )
        }
        else
        {
            SCAI_LOG_WARN( logger, LAMA_MIC_ENV_FOR_DEVICE << " not set, take device 0" )
            defaultDeviceNr = 0;
        }
    }

    return defaultDeviceNr;
}

/* ----------------------------------------------------------------------------- */

MICContextManager::MICContextManager()
    : ContextManager( Context::MIC )
{
    // Note: do not any logging here as the only one MICContextManager is created
    //       during static initialization, logger might not be available

    registerFactory();

    // initialize the weak pointers for different devices ( probably not necessary )

    for( int i = 0; i < LAMA_MAX_MIC_DEVICES; i++ )
    {
        mMICContext[i] = boost::weak_ptr<MICContext>();
    }
}

/* ----------------------------------------------------------------------------- */

MICContextManager::~MICContextManager()
{
    SCAI_LOG_DEBUG( logger, "~MICContextManager" )

    for( int i = 0; i < LAMA_MAX_MIC_DEVICES; i++ )
    {
        if( mMICContext[i].expired() )
        {
            SCAI_LOG_DEBUG( logger, "expired MICContext for device " << i )
        }
        else
        {
            SCAI_LOG_DEBUG( logger, "available MICContext for device " << i )
        }
    }
}

/* ----------------------------------------------------------------------------- */

ContextPtr MICContextManager::getInstance( int deviceNr )
{
    int micDeviceNr = deviceNr;

    if( micDeviceNr == LAMA_DEFAULT_DEVICE_NUMBER )
    {
        micDeviceNr = getDefaultDeviceNr();

        // no need here to check for a good value
    }
    else
    {
        SCAI_ASSERT_ERROR(
            0 <= micDeviceNr && micDeviceNr < LAMA_MAX_MIC_DEVICES,
            "device = " << micDeviceNr << " out of range" << ", max supported device = " << LAMA_MAX_MIC_DEVICES )
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
