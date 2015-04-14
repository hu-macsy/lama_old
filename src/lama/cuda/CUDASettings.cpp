/**
 * @file CUDASettings.cpp
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
 * @brief Implementation of utilties for CUDA Settings
 * @author Thomas Brandes
 * @date 24.05.2013
 */

#include <lama/Settings.hpp>
#include <lama/cuda/CUDASettings.hpp>

#include <lama/cuda/CUDAError.hpp>

// Some calls are needed of the CUDA driver API

#include <cuda.h>

#include <cstdlib>

namespace lama
{

LAMA_LOG_DEF_LOGGER( CUDASettings::logger, "CUDA.Settings" )

/* ----------------------------------------------------------------------------- */

bool CUDASettings::initialized = false;

/* ----------------------------------------------------------------------------- */

bool CUDASettings::theUseTextureFlag = false; // good value available after initialized

/* ----------------------------------------------------------------------------- */

// flag for using SharedMem will be set at initialization
bool CUDASettings::theUseSharedMemFlag = false; // good value available after initialized

/* ----------------------------------------------------------------------------- */

int CUDASettings::theBlockSize = 128; // good value available after initialized

/* ----------------------------------------------------------------------------- */

int CUDASettings::getComputeCapability()
{
    CUdevice dev; // curent device

    LAMA_CUDA_DRV_CALL( cuCtxGetDevice( &dev ), "get current device" )

    int major = 0;
    int minor = 0;

    // newer CUDA versions ( not tested yet )
    // CUdevice_attribute attr = CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR;
    // cuDeviceGetAttribute ( &major, attr, dev ),

    LAMA_CUDA_DRV_CALL( cuDeviceComputeCapability( &major, &minor, dev ), "getComputeCapability" )

    LAMA_LOG_INFO( logger, "compute capabilty = " << major << "." << minor );

    return major;
}

/* ----------------------------------------------------------------------------- */

void CUDASettings::initialize()
{
    // check environment variables for settings

    bool setTexture = Settings::getEnvironment( theUseTextureFlag, "LAMA_CUDA_USE_TEXTURE" );
    bool setSharedMem = Settings::getEnvironment( theUseSharedMemFlag, "LAMA_CUDA_USE_SHARED_MEM" );
    bool setBlockSize = Settings::getEnvironment( theBlockSize, "LAMA_CUDA_BLOCK_SIZE" );

    if( !setTexture || !setSharedMem || !setBlockSize )
    {
        // at least one environment variable not set, so define by compute capability

        int computeCapability = getComputeCapability();

        if( !setTexture )
        {
            theUseTextureFlag = computeCapability == 1;

            LAMA_LOG_INFO( logger,
                           "useTexture = " << theUseTextureFlag << ", due to compute capability = " << computeCapability );
        }

        if( !setSharedMem )
        {
            theUseSharedMemFlag = computeCapability == 1;

            LAMA_LOG_INFO( logger,
                           "useSharedMem = " << theUseSharedMemFlag << ", due to compute capability = " << computeCapability );
        }

        if( !setBlockSize )
        {
            if( computeCapability == 1 )
            {
                theBlockSize = 128;
            }
            else
            {
                theBlockSize = 256;
            }

            LAMA_LOG_INFO( logger,
                           "blockSize = " << theBlockSize << ", due to compute capability = " << computeCapability );
        }
    }

    initialized = true;
}

/* ----------------------------------------------------------------------------- */

bool CUDASettings::useTexture()
{
    if( !initialized )
    {
        initialize();
    }

    return theUseTextureFlag;
}

/* ----------------------------------------------------------------------------- */

bool CUDASettings::useSharedMem()
{
    if( !initialized )
    {
        initialize();
    }

    return theUseSharedMemFlag;
}

/* ----------------------------------------------------------------------------- */

int CUDASettings::getBlockSize()
{
    if( !initialized )
    {
        initialize();
    }

    return theBlockSize;
}

/* ----------------------------------------------------------------------------- */

int CUDASettings::getBlockSize( int n )
{
    if( !initialized )
    {
        initialize();
    }

    static int mp = 8; // multiprocessors available

    if( n >= mp * theBlockSize )
    {
        return theBlockSize;
    }
    else
    {
        // take smaller block size to get all multiprocessors occupied

        return theBlockSize > 1; // mod 2
    }
}

/* ----------------------------------------------------------------------------- */

void CUDASettings::set( bool useSharedMemFlag, bool useTextureFlag )
{
    theUseTextureFlag = useTextureFlag;
    theUseSharedMemFlag = useSharedMemFlag;

    initialized = true;

    LAMA_LOG_INFO( logger, "set: useSharedMem = " << theUseSharedMemFlag << ", useTexture = " << theUseTextureFlag )
}

/* ----------------------------------------------------------------------------- */

} // namespace

