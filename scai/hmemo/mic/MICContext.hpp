/**
 * @file MICContext.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Definition of context class for MIC devices and a context manager class.
 * @author Thomas Brandes
 * @date 01.07.2013
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
