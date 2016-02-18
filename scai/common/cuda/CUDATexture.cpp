/**
 * @file CUDATexture.cpp
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
 * @brief Implementation of utilties for CUDA Texture
 * @author Thomas Brandes
 * @date 04.05.2013
 * $Id$
 */

// hpp
#include <scai/common/cuda/CUDATexture.hpp>

// internal scai libraries
#include <scai/common/cuda/CUDAError.hpp>

// CUDA
#include <cuda.h> // Some calls are needed of the CUDA driver API

// std
#include <cstdlib>

namespace scai
{

namespace common
{

bool CUDATexture::initialized = false;

// flag for using Texture will be set at initialization

bool CUDATexture::theUseTextureFlag = getUseTextureByEnv();

bool CUDATexture::getUseTextureByEnv()
{
    const char* env = getenv( "SCAI_CUDA_USE_TEXTURE" );

    if( !env )
        return false; // no initialization by environment

    char key = static_cast<char>( toupper( env[0] ) );

    // to upper

    if( key == '0' )
    {
        theUseTextureFlag = false;
    }
    else if( key == '1' )
    {
        theUseTextureFlag = true;
    }
    else if( key == 'J' )
    {
        theUseTextureFlag = true;
    }
    else if( key == 'Y' )
    {
        theUseTextureFlag = true;
    }
    else if( key == 'N' )
    {
        theUseTextureFlag = false;
    }
    else
    {
        return false;
    }

    return true; // intialization done
}

void CUDATexture::setUseTextureByDevice()
{
    CUdevice dev; // curent device

    SCAI_CUDA_DRV_CALL( cuCtxGetDevice( &dev ), "get current device" )

    int major = 0;
    int minor = 0;

    // newer CUDA versions ( not tested yet )
    // CUdevice_attribute attr = CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR;
    // cuDeviceGetAttribute ( &major, attr, dev ),

    SCAI_CUDA_DRV_CALL( cuDeviceComputeCapability( &major, &minor, dev ), "getComputeCapability" )

    theUseTextureFlag = major == 1; // use texture on devices with compute capability 1.x

    initialized = true;
}

bool CUDATexture::useTexture()
{
    if( !initialized )
    {
        // initialization not done by enviroment variable, so do it by device

        setUseTextureByDevice();
    }

    return theUseTextureFlag;
}

} /* end namespace common */

} /* end namespace scai */

