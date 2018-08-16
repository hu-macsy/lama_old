/**
 * @file CUDASettings.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Implementation of utilties for CUDA Settings
 * @author Thomas Brandes
 * @date 24.05.2013
 */

// hpp
#include <scai/common/cuda/CUDASettings.hpp>

// internal scai library
#include <scai/common/Settings.hpp>
#include <scai/common/cuda/CUDAError.hpp>

// CUDA
#include <cuda.h> // Some calls are needed of the CUDA driver API

// std
#include <cstdlib>

namespace scai
{

namespace common
{

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
    SCAI_CUDA_DRV_CALL( cuCtxGetDevice( &dev ), "get current device" )
    int major = 0;
    int minor = 0;
    // newer CUDA versions ( not tested yet )
    // CUdevice_attribute attr = CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR;
    // cuDeviceGetAttribute ( &major, attr, dev ),
    SCAI_CUDA_DRV_CALL( cuDeviceComputeCapability( &major, &minor, dev ), "getComputeCapability" )
    return major;
}

/* ----------------------------------------------------------------------------- */

void CUDASettings::initialize()
{
    // check environment variables for settings
    bool setTexture = common::Settings::getEnvironment( theUseTextureFlag, "SCAI_CUDA_USE_TEXTURE" );
    bool setSharedMem = common::Settings::getEnvironment( theUseSharedMemFlag, "SCAI_CUDA_USE_SHARED_MEM" );
    bool setBlockSize = common::Settings::getEnvironment( theBlockSize, "SCAI_CUDA_BLOCK_SIZE" );

    if ( !setTexture || !setSharedMem || !setBlockSize )
    {
        // at least one environment variable not set, so define by compute capability
        int computeCapability = getComputeCapability();

        if ( !setTexture )
        {
            theUseTextureFlag = computeCapability == 1;
        }

        if ( !setSharedMem )
        {
            theUseSharedMemFlag = computeCapability == 1;
        }

        if ( !setBlockSize )
        {
            if ( computeCapability == 1 )
            {
                theBlockSize = 128;
            }
            else
            {
                theBlockSize = 256;
            }
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

int CUDASettings::getBlockSize()
{
    if ( !initialized )
    {
        initialize();
    }

    return theBlockSize;
}

/* ----------------------------------------------------------------------------- */

int CUDASettings::getBlockSize( int n )
{
    if ( !initialized )
    {
        initialize();
    }

    static int mp = 8; // multiprocessors available

    if ( n >= mp * theBlockSize )
    {
        return theBlockSize;
    }
    else
    {
        // take smaller block size to get all multiprocessors occupied
        return theBlockSize >> 1; // mod 2
    }
}

/* ----------------------------------------------------------------------------- */

void CUDASettings::set( bool useSharedMemFlag, bool useTextureFlag )
{
    theUseTextureFlag = useTextureFlag;
    theUseSharedMemFlag = useSharedMemFlag;
    initialized = true;
}

/* ----------------------------------------------------------------------------- */

} /* end namespace common */

} /* end namespace scai */
