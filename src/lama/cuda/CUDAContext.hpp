/**
 * @file CUDAContext.hpp
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
 * @brief Definition of context class for CUDA devices and a context manager class.
 * @author Thomas Brandes, Jiri Kraus
 * @date 15.07.2011
 * @since 1.0.0
 */
#ifndef LAMA_CUDA_CONTEXT_HPP_
#define LAMA_CUDA_CONTEXT_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/Context.hpp>
#include <boost/enable_shared_from_this.hpp>

// others
#include <lama/task/Thread.hpp>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>

// boost
#include <boost/weak_ptr.hpp>
#include <boost/version.hpp>

#include <string>

namespace lama
{

class LAMA_DLL_IMPORTEXPORT CUDAStreamSyncToken;

/**
 *  Define here the name of the environment variable that can be used
 *  to specify the default CUDA device.
 */

#define LAMA_CUDA_ENV_FOR_DEVICE "LAMA_DEVICE"

/**
 * @brief CUDAContext initializes the CUDA device with the given number.
 *
 * CUDAContext initializes the CUDA device with the given number according to the
 * RAII Programming Idiom (http://www.hackcraft.net/raii/). The device is
 * initialized in the constructor and cleaned up in the destructor. CUDAContext
 * uses a static counter to avoid multiple device initializations.
 */
class LAMA_DLL_IMPORTEXPORT CUDAContext: public Context, public boost::enable_shared_from_this<CUDAContext>
{

    friend class CUDAContextManager;

public:

    /**
     * @brief The destructor destroys this CUDA device, and frees the initialized
     *        CUDA device if needed.
     */
    virtual ~CUDAContext();

    /**
     *  Override Context::getHostContext
     *
     *  CUDAHostContext is used as host context to support faster and asynchronous memory transfer.
     */
    virtual ContextPtr getHostContext() const;

    int getDeviceNr() const
    {
        return mDeviceNr;
    }

    /** Each host context can use data of each other host context.  */

    virtual bool canUseData( const Context& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    virtual void* allocate( const size_t size ) const;

    virtual void allocate( ContextData& contextData, const size_t size ) const;

    virtual void free( void* pointer, const size_t size ) const;

    virtual void free( ContextData& contextData ) const;

    virtual void memcpy( void* dst, const void* src, const size_t size ) const;

    virtual SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const;

    virtual bool cancpy( const ContextData& dst, const ContextData& src ) const;

    virtual void memcpy( ContextData& dst, const ContextData& src, const size_t size ) const;

    /** This method implements Context::memcpyAsync */

    virtual SyncToken* memcpyAsync( ContextData& dst, const ContextData& src, const size_t size ) const;

    /** The CUDA interface used for the implementation requires that the device
     *  must be set via a setDevice routine. This method takes care of it if this
     *  context device is not the current one. So this method must be called before
     *  any CUDA code is executed (includes also memory transfer routines).
     */

    virtual void enable( const char* filename, int line ) const;

    virtual void disable( const char* filename, int line ) const;

    CUDAStreamSyncToken* getComputeSyncToken() const;

    CUDAStreamSyncToken* getTransferSyncToken() const;

    /** Implementation for Context::getSyncToken */

    virtual SyncToken* getSyncToken() const;

protected:

    /**
     * @brief Constructor for a context on a certain device.
     *
     * @param device    the number of the CUDA device to initialize. [IN]
     * @throw           Exception if the device initialization fails.
     *
     * If device is DEFAULT_DEVICE_NUMBER, the device number is taken by
     * the environment variable LAMA_DEVICE. If this variabe is not set,
     * it will take device 0 as default.
     */
    CUDAContext( int device );

private:

    void memcpyFromHost( void* dst, const void* src, const size_t size ) const;
    void memcpyToHost( void* dst, const void* src, const size_t size ) const;
    void memcpyFromCUDAHost( void* dst, const void* src, const size_t size ) const;
    void memcpyToCUDAHost( void* dst, const void* src, const size_t size ) const;

    SyncToken* memcpyAsyncFromHost( void* dst, const void* src, const size_t size ) const;
    SyncToken* memcpyAsyncToHost( void* dst, const void* src, const size_t size ) const;
    SyncToken* memcpyAsyncFromCUDAHost( void* dst, const void* src, const size_t size ) const;
    SyncToken* memcpyAsyncToCUDAHost( void* dst, const void* src, const size_t size ) const;

    mutable boost::weak_ptr<const class Context> mHostContext;  //!< preferred host context

    int mDeviceNr; //!< number of device for this context

    CUdevice mCUdevice; //!< data structure for device

    CUcontext mCUcontext; //!< data structure for context

    CUstream mTransferStream; //!< stream for memory transfers
    CUstream mComputeStream; //!< stream for asynchronous computations

    cublasHandle_t mCublasHandle; //!< handle to cublas library
    cusparseHandle_t mCusparseHandle; //!< handle to cusparse library

    std::string mDeviceName; //!< name set during initialization

    Thread::Id mOwnerThread;

    mutable int mNumberOfAllocates; //!< variable counts allocates
    mutable long long mNumberOfAllocatedBytes; //!< variable counts allocated bytes on device
    mutable long long mMaxNumberOfAllocatedBytes; //!< variable counts the maximum allocated bytes

    static int currentDeviceNr; //!< number of device currently set for CUDA

    static int numUsedDevices; //!< total number of used devices

    static size_t minPinnedSize;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} //namespace lama

#endif // LAMA_CUDA_CONTEXT_HPP_
