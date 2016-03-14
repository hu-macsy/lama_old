/**
 * @file CUDADevice.hpp
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
 * @brief Class for a CUDA device
 * @author: Thomas Brandes
 * @date 08.03.2016
 **/

#pragma once

// for dll_import
#include <scai/common/config.hpp>


#include <scai/common/NonCopyable.hpp>
#include <scai/common/function.hpp>

#include <cuda.h>
#include <cublas_v2.h>

#include <vector>

namespace scai
{

namespace common
{

/* --------------------------------------------------------------------- */

class COMMON_DLL_IMPORTEXPORT CUDADevice : private NonCopyable
{   

public:

    /** Constructor will check for a CUDA device and create a context for it 
     *
     *  @param[in] deviceNr is the number of the device to access
     *  @throws Exception if the device cannot be accessed
     */

    CUDADevice( int deviceNr );

    /** Destructor destroys the CUDA context used for the device.
     *
     *  Be careful: all allocated resources on the device should have been freed before.
     */

    ~CUDADevice();

    /** Getter for the CUcontext to be used for CUDA driver API operations. */

    CUcontext getCUcontext() const
    {
        return mCUcontext;
    }

    /** Getter for the CUdevice to be used for CUDA driver API operations. */

    CUdevice getCUdevice() const
    {
        return mCUdevice;
    }

    /** Getter for the cuBLAS handle, will be created with first use */

    cublasHandle_t getcuBLASHandle() const;

    /** Getter for the number of the device. */

    int getDeviceNr() const
    {
        return mDeviceNr;
    }
    
    /** Add a routine to be called with destructor */

    void addShutdown( common::function<void()> routine );

private:

    CUcontext mCUcontext;

    CUdevice mCUdevice;

    int mDeviceNr;

    cublasHandle_t mcuBLASHandle;

    std::vector< common::function<void()> > mShutdownFunctions;
};

}  // namespace common

}  // namespace scai
