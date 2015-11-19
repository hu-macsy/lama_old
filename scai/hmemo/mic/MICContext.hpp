/**
 * @file MICContext.hpp
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
 * @brief Definition of context class for MIC devices and a context manager class.
 * @author Thomas Brandes
 * @date 01.07.2013
 * @since 1.1.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/Context.hpp>

// internal scai libraries
#include <scai/tasking/SyncToken.hpp>

#include <scai/logging.hpp>

#include <scai/common/shared_ptr.hpp>
#include <scai/common/weak_ptr.hpp>

// std
#include <string>

namespace scai
{

namespace hmemo
{

/**
 * @brief MICContext initializes the MIC device with the given number.
 *
 */
class COMMON_DLL_IMPORTEXPORT MICContext: 

    public Context, 
    public Context::Register<MICContext>,
    public common::enable_shared_from_this<MICContext>
{

public:

    /**
     * @brief The destructor destroys this MIC device, and frees the initialized
     *        MIC device if needed.
     */
    virtual ~MICContext();

    int getDeviceNr() const
    {
        return mDeviceNr;
    }

    virtual void writeAt( std::ostream& stream ) const;

    virtual void enable( const char* filename, int line ) const;

    virtual void disable( const char* filename, int line ) const;

    /** Getter routine for a new sync token that allows to asynchronous computations on the context. */

    virtual tasking::SyncToken* getSyncToken() const;

    static int getCurrentDevice();

    /** This routine is required for Register in Context Factory. */

    static common::context::ContextType createValue()
    {
        return common::context::MIC;
    }

    /** This routine is required for Register in Context Factory. */

    static ContextPtr create( int deviceNr );

    /**
     *  No Override for Context::getHostMemoryPtr
     */

    /**
     *  @brief Implementation of Context::getLocalMemoryPtr for this class.
     */
    virtual MemoryPtr getLocalMemoryPtr() const;

    /** 
     *  @brief Implementation of Context::canUseMemory for this class. 
     */

    virtual bool canUseMemory( const Memory& other ) const;

protected:

    /**
     * @brief Constructor for a context on a certain device.
     *
     * @param device    the number of the MIC device to initialize. [IN]
     * @throw           Exception if the device initialization fails.
     *
     * If device is DEFAULT_DEVICE_NUMBER, the device number is taken by
     * the environment variable SCAI_DEVICE. If this variabe is not set,
     * it will take device 0 as default.
     */
    MICContext( int device );

private:

    mutable common::weak_ptr<class Memory> mMemory;     //!< memory management for this devie

    static int currentDeviceNr; //!< number of device currently set for MIC

    static int numUsedDevices; //!< total number of used MIC devices

    int mDeviceNr; //!< number of device for this context

    int mNumThreads; //!< number of threads used on this device

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace hmemo */

} /* end namespace scai */
