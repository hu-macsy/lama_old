/**
 * @file CUDAContext.hpp
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
 * @brief Definition of context class for CUDA devices.
 * @author Thomas Brandes, Jiri Kraus
 * @date 08.07.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/Context.hpp>

// others
#include <scai/common/Thread.hpp>
#include <scai/common/weak_ptr.hpp>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>

#include <string>

namespace scai
{

namespace tasking
{
    class CUDAStreamSyncToken;  // forward 
}

namespace hmemo
{

/**
 *  Define here the name of the environment variable that can be used
 *  to specify the default CUDA device.
 */

#define SCAI_CUDA_ENV_FOR_DEVICE "SCAI_DEVICE"

/**
 * @brief CUDAContext initializes the CUDA device with the given number.
 *
 * CUDAContext initializes the CUDA device with the given number according to the
 * RAII Programming Idiom (http://www.hackcraft.net/raii/). The device is
 * initialized in the constructor and cleaned up in the destructor. CUDAContext
 * uses a static counter to avoid multiple device initializations.
 */
class COMMON_DLL_IMPORTEXPORT CUDAContext: 

    public Context, 
    public Context::Register<CUDAContext>,
    public common::enable_shared_from_this<CUDAContext>
{

public:

    /**
     * @brief The destructor destroys this CUDA device, and frees the initialized
     *        CUDA device if needed.
     */
    virtual ~CUDAContext();

    /**
     *  Override Context::getHostMemoryPtr
     *
     *  CUDAHostMemory is used as host memory to support faster and asynchronous memory transfer.
     */
    virtual MemoryPtr getHostMemoryPtr() const;

    /**
     *  @brief Implementation of Context::getMemory for this class.
     */
    virtual MemoryPtr getMemoryPtr() const;

    int getDeviceNr() const
    {
        return mDeviceNr;
    }

    /** 
     *  @brief Implementation of Context::canUseMemory for this class. 
     */

    virtual bool canUseMemory( const Memory& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    /** The CUDA interface used for the implementation requires that the device
     *  must be set via a setDevice routine. This method takes care of it if this
     *  context device is not the current one. So this method must be called before
     *  any CUDA code is executed (includes also memory transfer routines).
     */

    virtual void enable( const char* filename, int line ) const;

    virtual void disable( const char* filename, int line ) const;

    tasking::CUDAStreamSyncToken* getComputeSyncToken() const;

    tasking::CUDAStreamSyncToken* getTransferSyncToken() const;

    /** 
     *  @brief Implementation for Context::getSyncToken 
     */

    virtual tasking::SyncToken* getSyncToken() const;

    // Using this context might require accessing the context

    CUcontext getCUcontext() const
    {
        return mCUcontext;
    }

    /** This routine is required for Register in Context Factory. */

    static ContextType createValue() 
    {
        return context::CUDA;
    }

    /** This routine is required for Register in Context Factory. */

    static ContextPtr create( int deviceNr );

protected:

    /**
     * @brief Constructor for a context on a certain device.
     *
     * @param device    the number of the CUDA device to initialize. [IN]
     * @throw           Exception if the device initialization fails.
     *
     * If device is DEFAULT_DEVICE_NUMBER, the device number is taken by
     * the environment variable SCAI_DEVICE. If this variabe is not set,
     * it will take device 0 as default.
     */
    CUDAContext( int device );

private:

    CUcontext mCUcontext; //!< data structure for context

    mutable common::weak_ptr<class Memory> mMemory;     //!< memory management for this devie
    mutable common::weak_ptr<class Memory> mHostMemory; //!< preferred host memory

    int mDeviceNr; //!< number of device for this context

    CUdevice mCUdevice; //!< data structure for device

    CUstream mTransferStream; //!< stream for memory transfers
    CUstream mComputeStream; //!< stream for asynchronous computations

//    cublasHandle_t   mCublasHandle;   //!< handle to cublas library
//    cusparseHandle_t mCusparseHandle; //!< handle to cusparse library

    std::string mDeviceName; //!< name set during initialization

    common::Thread::Id mOwnerThread;

    static int currentDeviceNr; //!< number of device currently set for CUDA

    static int numUsedDevices; //!< total number of used devices

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace hmemo */

} /* end namespace scai */
