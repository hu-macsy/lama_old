/**
 * @file CUDAContextManager.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Contains the implementation of the singleton class CUDAContextManager.
 * @author Thomas Brandes
 * Created on: 15.07.2011
 */

// hpp
#include <lama/cuda/CUDAContextManager.hpp>

// others
#include <lama/cuda/CUDAContext.hpp>

#include <lama/ContextFactory.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

namespace lama

{

/* ----------------------------------------------------------------------------- */

// definition of array with weak pointers so that we can return shared pointers without allocating again
boost::weak_ptr<CUDAContext> CUDAContextManager::mCUDAContext[LAMA_MAX_CUDA_DEVICES];

CUDAContextManager CUDAContextManager::theInstance;

CUDAContextManager::CUDAContextManager()
    : ContextManager( Context::CUDA )
{
    // Note: do not any logging here as the only one CUDAContextManager is created
    //       during static initialization, logger might not be available

    registerFactory();

    // initialize the weak pointers for different devices ( probably not necessary )

    for ( int i = 0; i < LAMA_MAX_CUDA_DEVICES; i++ )
    {
        mCUDAContext[i] = boost::weak_ptr<CUDAContext>();
    }
}

/* ----------------------------------------------------------------------------- */

CUDAContextManager::~CUDAContextManager()
{
    LAMA_LOG_DEBUG( logger, "~CUDAContextManager" )

    for ( int i = 0; i < LAMA_MAX_CUDA_DEVICES; i++ )
    {
        if ( mCUDAContext[i].expired() )
        {
            LAMA_LOG_DEBUG( logger, "expired CUDAContext for device " << i )
        }
        else
        {
            LAMA_LOG_DEBUG( logger, "available CUDAContext for device " << i )
        }
    }
}

/* ----------------------------------------------------------------------------- */

ContextPtr CUDAContextManager::getInstance( int deviceNr )
{
    int cudaDeviceNr = deviceNr;

    if ( cudaDeviceNr == LAMA_DEFAULT_DEVICE_NUMBER )
    {
        // if no device has been specified we take the value of the environment variable

        if ( getenv( LAMA_CUDA_ENV_FOR_DEVICE ) )
        {
            std::string devNumber( getenv( LAMA_CUDA_ENV_FOR_DEVICE ) );
            std::istringstream devNumberReader( devNumber );
            devNumberReader >> cudaDeviceNr;

            LAMA_LOG_INFO( logger, LAMA_CUDA_ENV_FOR_DEVICE << " = " << cudaDeviceNr << " set, take it" )
        }
        else
        {
            LAMA_LOG_WARN( logger, LAMA_CUDA_ENV_FOR_DEVICE << " not set, take device 0" )
            cudaDeviceNr = 0;
        }
    }

    LAMA_ASSERT_ERROR(
        0 <= cudaDeviceNr && cudaDeviceNr < LAMA_MAX_CUDA_DEVICES,
        "device = " << cudaDeviceNr << " out of range, max supported device = " << LAMA_MAX_CUDA_DEVICES )

    boost::shared_ptr<CUDAContext> context = boost::shared_ptr<CUDAContext>();

    if ( mCUDAContext[cudaDeviceNr].expired() )
    {
        // create a new context for the device and return the shared pointer

        context = boost::shared_ptr<CUDAContext>( new CUDAContext( cudaDeviceNr ) );

        // we keep a weak pointer so that we can return

        mCUDAContext[cudaDeviceNr] = context;
    }
    else
    {
        // the weak pointer to the device is still okay, so return a shared pointer for it

        context = mCUDAContext[cudaDeviceNr].lock();
    }

    return context;
}

} //namespace lama

