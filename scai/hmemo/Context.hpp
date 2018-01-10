/**
 * @file Context.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Definition of base class for a context that specifies where data is
 *        allocated and where computation takes place.
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

#include <stack>
#include <memory>

namespace scai
{

namespace tasking
{
class SyncToken;    // forward declaration
}

/** Namespace for all data structures of the the heterogeneous memory management. */

namespace hmemo
{

class Memory;  // forward declaration

typedef std::shared_ptr<Memory> MemoryPtr;

class Context;   // forward declaration

/** Context pointers will be always const, so context can never be modified. */

typedef std::shared_ptr<const Context> ContextPtr;

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

    public  common::Factory1<common::ContextType, int, ContextPtr>,
    public  common::Printable,
    private common::NonCopyable
{
public:

    virtual ~Context();

    /** Method to get the type of the context. */

    inline common::ContextType getType() const;

    /** @brief  Predicate to check in a context whether a certain memory class can be used.
     *
     *  If an incarnation of a HArray has a valid copy that can be used, no additional memory
     *  or memcopy is required. But be careful: it might be faster to use other memory at this
     *  context.
     *
     *  This pure routine must be implemented by derived classes.
     *
     *  @param[in] memory is the memory against which the check is done
     */
    virtual bool canUseMemory( const Memory& memory ) const = 0;

    virtual void writeAt( std::ostream& stream ) const;

    /** Getter routine for a new sync token that allows to asynchronous computations on the context.  //!< if true getMemoryPtr() returns HostMemory
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

    virtual MemoryPtr getLocalMemoryPtr() const = 0;

    /** Get the preferred host memory for a context.
     *
     *  Depending on the device, the host data might be allocated
     *  in such a way that faster and/or asynchronous memory transfer is supported.
     *
     *  As default, the standard Host memory is returned. Derived classes of
     *  Context should override this method to provide more efficient solutions.
     */

    virtual MemoryPtr getHostMemoryPtr() const;

    /** Get the memory on which the context will work. */

    MemoryPtr getMemoryPtr() const;

    /** If zero copy is enabled, the default memory is not the local memory
     *  but the host memory so Host and this context can work on the same memory.
     */

    virtual void enableZeroCopy( bool flag ) const;

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
    static ContextPtr getContextPtr( const common::ContextType type, int deviceNr = -1 );

    /** @brief get context as set by SCAI_CONTEXT and SCAI_DEVICE
     *
     *  @return             a context of the type as set by environment variable SCAI_CONTEXT (Host if not set)
     *  @throws Exception if the context set by SCAI_CONTEXT is not available or unknown
     */

    static ContextPtr getContextPtr();

    /** @brief getHostPtr() as abbreviation of getContextPtr( common::ContextType::Host ) */

    static inline ContextPtr getHostPtr();

    /** Checks if a context of the passed type is available.
     *
     * @param[in] type  is the type of context that is wanted
     * @return          if a context of the passed type is available
     */
    static inline bool hasContext( const common::ContextType type );

    /** Get the currently accessed context of this thread */

    static const Context* getCurrentContext();

protected:

    /** Default constructor, can only be called by base classes. */

    Context( common::ContextType type );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    common::ContextType mContextType;

    mutable bool mUseZeroCopy;   //!< if true getMemoryPtr() returns HostMemory

    mutable bool mEnabled; //!<  if true the context is currently accessed

    mutable const char* mFile;//!< File name where context has been enabled

    mutable int mLine;//!< Line number where context has been enabled

private:

    /**
     *  Set this Context as the current one of this thread.
     *
     *  Only one Context can be the current one.
     *
     *  Global access to the current context makes design easier
     *  as it can be decided locally. 
     */

    void setCurrent() const;

    /**
     *  Current Context will be no more current one.
     */
    void unsetCurrent() const;

    /** thread-private variable where the current context of a thread can be asked for */

    typedef std::stack<const Context*> ContextStack;

    static thread_local ContextStack contextStack;
};

/* ======================================================================== */
/*             Inline methods                                               */
/* ======================================================================== */

common::ContextType Context::getType() const
{
    return mContextType;
}

bool Context::hasContext( const common::ContextType type )
{
    return canCreate( type );
}

ContextPtr Context::getHostPtr()
{
    return getContextPtr( common::ContextType::Host );
}

} /* end namespace hmemo */

} /* end namespace scai */
