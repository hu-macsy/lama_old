/**
 * @file CUDASettings.cpp
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
 * @brief Implementation of utilties for CUDA Settings
 * @author Thomas Brandes
 * @date 24.05.2013
 */

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

bool CUDASettings::theUseTextureFlag = false;  // good value available after initialized

/* ----------------------------------------------------------------------------- */

// flag for using SharedMem will be set at initialization

bool CUDASettings::theUseSharedMemFlag = false;  // good value available after initialized

/* ----------------------------------------------------------------------------- */

bool CUDASettings::convertYesNoString( bool& flag, const char* stringVal )
{
    char key = toupper( stringVal[0] );

    bool done = true;  // becomes false if no legal value has been found

    // to upper

    if ( key == '0' ) 
    {
        flag = false;
    }
    else if ( key == '1' ) 
    {
        flag = true;
    }
    else if ( key == 'J' ) 
    {
        flag = true;
    }
    else if ( key == 'Y' ) 
    {
        flag = true;
    }
    else if ( key == 'T' ) 
    {
        flag = true;
    }
    else if ( key == 'N' ) 
    {
        flag = false;
    }
    else if ( key == 'F' ) 
    {
        flag = false;
    }
    else
    {
        // could not identify meaning
        done = false;
    }
 
    return done;
}

/* ----------------------------------------------------------------------------- */

int CUDASettings::getComputeCapability()
{
    CUdevice dev;   // curent device

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

bool CUDASettings::getEnvironmentSetting( bool& flag, const char* envVarName )
{
    const char* env = getenv ( envVarName );

    if ( !env ) 
    {
        LAMA_LOG_INFO( logger, envVarName << " not set, will select by compute capability" )

        return false;   // no initialization by environment
    }

    bool done = convertYesNoString( flag, env );

    if ( !done )
    {
        LAMA_LOG_ERROR( logger, "Environment variable " << envVarName << "=" << env 
                                << ", is illegal setting, assume FALSE" )

        flag = false;
    }

    return true;   // environment variable was available
}

/* ----------------------------------------------------------------------------- */

void CUDASettings::initialize()
{
    // check environment variables for settings

    bool setTexture   = getEnvironmentSetting ( theUseTextureFlag, "LAMA_CUDA_USE_TEXTURE" );
    bool setSharedMem = getEnvironmentSetting ( theUseSharedMemFlag, "LAMA_CUDA_USE_SHARED_MEM" );

    if ( !setTexture || !setSharedMem )
    {
        // at least one environment variable not set, so define by compute capability

        int computeCapability = getComputeCapability();

        if ( !setTexture )
        {
            theUseTextureFlag = computeCapability == 1;

            LAMA_LOG_INFO( logger, "useTexture = " << theUseTextureFlag 
                                   << ", due to compute capability = " << computeCapability );
        }

        if ( !setSharedMem )
        {
            theUseSharedMemFlag = computeCapability == 1;

            LAMA_LOG_INFO( logger, "useSharedMem = " << theUseSharedMemFlag 
                                   << ", due to compute capability = " << computeCapability );
        }
    }

    initialized = true;
}

/* ----------------------------------------------------------------------------- */

bool CUDASettings::useTexture()
{
    if ( !initialized )
    {
        initialize();
    }

    return theUseTextureFlag;
}

/* ----------------------------------------------------------------------------- */

bool CUDASettings::useSharedMem()
{
    if ( !initialized )
    {
        initialize();
    }

    return theUseSharedMemFlag;
}

/* ----------------------------------------------------------------------------- */

void CUDASettings::set( bool useSharedMemFlag, bool useTextureFlag )
{
    theUseTextureFlag   = useTextureFlag;
    theUseSharedMemFlag = useSharedMemFlag;

    initialized = true;

    LAMA_LOG_INFO( logger, "set: useSharedMem = " << theUseSharedMemFlag 
                           << ", useTexture = " << theUseTextureFlag )
}

/* ----------------------------------------------------------------------------- */

}  // namespace

