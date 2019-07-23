/**
 * @file CUDAContext.hpp
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
 * @brief Definition of context class for CUDA devices.
 * @author Thomas Brandes, Jiri Kraus
 * @date 08.07.2015
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/Context.hpp>
#include <scai/common/cuda/CUDACtx.hpp>

// CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>

// std
#include <string>
#include <stack>
#include <memory>

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
    public std::enable_shared_from_this<CUDAContext>,
    public common::CUDACtx
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
     *  @brief Implementation of Context::getLocalMemory for this class.
     */
    virtual MemoryPtr getLocalMemoryPtr() const;

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

    /** This routine is required for Register in Context Factory. */

    static common::ContextType createValue()
    {
        return common::ContextType::CUDA;
    }

    /** Implementation of pure method Context::isEqual */

    bool isEqual( const Context& other ) const;

    /** This routine is required for Register in Context Factory. */

    static ContextPtr create( int deviceNr );

    using common::CUDACtx::getcuBLASHandle;

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

    mutable std::weak_ptr<class Memory> mMemory;     //!< memory management for this devie
    mutable std::weak_ptr<class Memory> mHostMemory; //!< preferred host memory

    //    cublasHandle_t   mCublasHandle;   //!< handle to cublas library
    //    cusparseHandle_t mCusparseHandle; //!< handle to cusparse library

    std::string mDeviceName; //!< name set during initialization

    // static int numUsedDevices; //!< total number of used devices

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace hmemo */

} /* end namespace scai */
