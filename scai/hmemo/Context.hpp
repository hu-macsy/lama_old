/**

 * @file Context.hpp
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
 * @brief Definition of base class for a context that specifies where data is
 *        allocated and where computation takes place.
 *
 * @author Thomas Brandes
 * @date 14.07.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Factory1.hpp>
#include <scai/common/Printable.hpp>
#include <scai/common/NonCopyable.hpp>
#include <scai/common/ContextType.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/shared_ptr.hpp>
#include <scai/common/function.hpp>

namespace scai
{

namespace tasking
{
    class SyncToken;    // forward declaration
}

/** Namespace for all data structures of the context memory management. */

namespace hmemo
{

class Memory;  // forward declaration

typedef common::shared_ptr<Memory> MemoryPtr;

class Context;   // forward declaration

/** Context pointers will be always const, so context can never be modified. */

typedef common::shared_ptr<const Context> ContextPtr;

/** @brief This class is a common base class for all possible contexts.
 *
 *  A context stands for a compute and data environment to run certain code.
 *
 *  Each context must provide routines to allocate and free data for the context as
 *  well as copy data in suitable memory locations accessible for the context.
 *
 *  Furthermore, each context should provide routines to transfer memory from the HOST to
 *  to the context and vice versa.
 *
 *  A copy constructor for a context is not provided.
 */
class COMMON_DLL_IMPORTEXPORT Context: 
  
    public  common::Factory1<common::context::ContextType, int, ContextPtr>,
    public  common::Printable,
    private common::NonCopyable
{
public:

    virtual ~Context();

    /** Method to get the type of the context. */

    common::context::ContextType getType() const;

    /** @brief  Predicate to check in a context whether a certain memory class can be used.
     *
     *  If an incarnation of a LAMAArray has a valid copy that can be used, no additional memory
     *  or memcopy is required. But be careful: it might be faster to use other memory at this
     *  context.
     *
     *  This pure routine must be implemented by derived classes.
     *
     *  @param[in] memory is the memory against which the check is done
     */
    virtual bool canUseMemory( const Memory& memory ) const = 0;

    virtual void writeAt( std::ostream& stream ) const;

    /** Getter routine for a new sync token that allows to asynchronous computations on the context.
     *
     *  @returns new SyncToken object
     */

    virtual tasking::SyncToken* getSyncToken() const = 0;

    /**
     * @brief Enable computations in the context.
     *
     * Operations on the context like allocation, deallocation of data, computations can
     * only be done if the context is enabled.
     *
     * Different threads can enable a context, but only one at a time.
     */
    virtual void enable( const char* file, int line ) const;

    /**
     * @brief Disable computations in the context.
     */
    virtual void disable( const char* file, int line ) const;

    /** This method returns the memory that can be used at this context. 
     *
     *  Note: canUseMemory( *getMemory() ) must be true.
     *
     *  It might be possible that other memories can also be used. This
     *  method should return the memory that gives the best performance.
     */

    virtual MemoryPtr getMemoryPtr() const = 0;

    /** Get the preferred host memory for a context.
     *
     *  Depending on the device, the host data might be allocated
     *  in such a way that faster and/or asynchronous memory transfer is supported.
     *
     *  As default, the standard Host memory is returned. Derived classes of
     *  Context should override this method to provide more efficient solutions.
     */

    virtual MemoryPtr getHostMemoryPtr() const;
 
    /** @brief Get a context of a certain type from the Context factory.
     *
     *  Note: This is the same as Factory::create but with default values.
     *
     *  The shared pointer guarantees that a context will live at least
     *  as long as arrays are allocated on it.
     *
     *  @param[in] type     is the type of context that is wanted
     *  @param[in] deviceNr is used for multiple devices of the same type (default is -1 for not set)
     *
     *  @return             a context of the requested type, if available.
     *
     *  @throws Exception if the context of the requested type is not available
     */
    static ContextPtr getContextPtr( const common::context::ContextType type = common::context::Host, int deviceNr = -1 );

    /** @brief getHostPtr() as abbreviation of getContextPtr( context::Host ) */

    static ContextPtr getHostPtr()
    {
        return getContextPtr( common::context::Host );
    }

    /** Checks if a context of the passed type is available.
     *
     * @param[in] type  is the type of context that is wanted
     * @return          if a context of the passed type is available
     */
    static bool hasContext( const common::context::ContextType type );

protected:

    /** Default constructor, can only be called by base classes. */

    Context( common::context::ContextType type );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    common::context::ContextType mContextType;

    mutable bool mEnabled; //!<  if true the context is currently accessed

    mutable const char* mFile;//!< File name where context has been enabled

    mutable int mLine;//!< Line number where context has been enabled
};

inline common::context::ContextType Context::getType() const
{
    return mContextType;
}

inline bool Context::hasContext( const common::context::ContextType type )
{
    return canCreate( type );
}

/** Output of context type in stream. */


/** Output of AccessKind in stream is supported and very useful.  */

} /* end namespace hmemo */

} /* end namespace scai */
