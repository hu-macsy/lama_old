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
#ifndef LAMA_MIC_CONTEXT_HPP_
#define LAMA_MIC_CONTEXT_HPP_

// for dll_import
#include <common/config.hpp>

// base classes
#include <lama/Context.hpp>
#include <boost/enable_shared_from_this.hpp>

#include <logging/logging.hpp>

// others
#include <string>

namespace lama
{

/**
 * @brief MICContext initializes the MIC device with the given number.
 *
 */
class COMMON_DLL_IMPORTEXPORT MICContext: public Context, public boost::enable_shared_from_this<MICContext>
{

    friend class MICContextManager;

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

    /** Each host contxt can data of each other host contxt.  */

    virtual bool canUseData( const Context& other ) const;

    virtual void writeAt( std::ostream& stream ) const;

    virtual void* allocate( const size_t size ) const;

    virtual void allocate( ContextData& contextData, const size_t size ) const;

    virtual void free( void* pointer, const size_t size ) const;

    virtual void free( ContextData& contextData ) const;

    virtual void memcpy( void* dst, const void* src, const size_t size ) const;

    virtual void memcpyToHost( void* dst, const void* src, const size_t size ) const;

    virtual void memcpyFromHost( void* dst, const void* src, const size_t size ) const;

    virtual SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const;

    virtual bool cancpy( const ContextData& dst, const ContextData& src ) const;

    virtual void memcpy( ContextData& dst, const ContextData& src, const size_t size ) const;

    /** This method implements Context::memcpyAsync */

    virtual SyncToken* memcpyAsync( ContextData& dst, const ContextData& src, const size_t size ) const;

    /** The MIC interface used for the implementation requires that the device
     *  must be set via a setDevice routine. This method takes care of it if this
     *  context device is not the current one. So this method must be called before
     *  any CUDA code is executed (includes also memory transfer routines).
     */

    virtual void enable( const char* filename, int line ) const;

    virtual void disable( const char* filename, int line ) const;

    /** Getter routine for a new sync token that allows to asynchronous computations on the context. */

    virtual SyncToken* getSyncToken() const;

    static int getCurrentDevice();

protected:

    /**
     * @brief Constructor for a context on a certain device.
     *
     * @param device    the number of the MIC device to initialize. [IN]
     * @throw           Exception if the device initialization fails.
     *
     * If device is DEFAULT_DEVICE_NUMBER, the device number is taken by
     * the environment variable LAMA_DEVICE. If this variabe is not set,
     * it will take device 0 as default.
     */
    MICContext( int device );

private:

    static int currentDeviceNr; //!< number of device currently set for MIC

    static int numUsedDevices; //!< total number of used MIC devices

    int mDeviceNr; //!< number of device for this context

    int mNumThreads; //!< number of threads used on this device

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

}
//namespace lama

#endif // LAMA_MIC_CONTEXT_HPP_
