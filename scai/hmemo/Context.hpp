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

namespace lama
{
    class LAMAInterface;
}

/** Namespace for all data structures of the context memory management. */

namespace hmemo
{

class Memory;  // forward declaration

typedef scai::common::shared_ptr<Memory> MemoryPtr;

class Context;   // forward declaration

/** Context pointers will be always const, so context can never be modified. */

typedef scai::common::shared_ptr<const Context> ContextPtr;

/** Namespace for enumeration of context types and access kinds. */

namespace context
{
    /** Enumeration type for the supported contexts. The type is used to select
     *  the appropriate code that will be used for the computations in the context.
     *
     *  The same context type does not imply that two different contexts can use
     *  the same data. Two CUDA contexts might allocate their own data where data
     *  must be transfered explicitly.
     */
    enum ContextType
    {
        Host,          //!< context for cpu + main memory
        CUDA,          //!< CUDA GPU device
        OpenCL,        //!< OpenCL GPU device, currently not supported
        MIC,           //!< Intel MIC
        UserContext,   //!< can be used for a new derived Context class
        MaxContext     //!< used for dimension of ContextType arrays
    };

    /** Enumeration type for access kind, may be read or write */

    enum AccessKind
    {
        Read, //!<  read access to the array, can be multiple
        Write, //!<  write access to the array, only one at a time
        MaxAccessKind //!<  internal use for dimension of arrays
    };

    COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const ContextType& type );

    COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const AccessKind& kind );
}

// Make ContexType and AccessKind visible, but not enum values. 

using context::ContextType;
using context::AccessKind;

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
  
    public  scai::common::Factory1<ContextType, int, ContextPtr>,
    public  scai::common::Printable,
    private scai::common::NonCopyable
{
public:

    virtual ~Context();

    /** Method to get the type of the context. */

    ContextType getType() const;

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

    /** This method returns the LAMA interface for a given context. */

    const lama::LAMAInterface& getInterface() const;

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
    static ContextPtr getContextPtr( const ContextType type = context::Host, int deviceNr = -1 );

    /** @brief getHostPtr() as abbreviation of getContextPtr( context::Host ) */

    static ContextPtr getHostPtr()
    {
        return getContextPtr( context::Host );
    }

    /** Checks if a context of the passed type is available.
     *
     * @param[in] type  is the type of context that is wanted
     * @return          if a context of the passed type is available
     */
    static bool hasContext( const ContextType type );

protected:

    /** Default constructor, can only be called by base classes. */

    Context( ContextType type );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    ContextType mContextType;

    mutable bool mEnabled; //!<  if true the context is currently accessed

    mutable const char* mFile;//!< File name where context has been enabled

    mutable int mLine;//!< Line number where context has been enabled
};

inline ContextType Context::getType() const
{
    return mContextType;
}

inline bool Context::hasContext( const ContextType type )
{
    return canCreate( type );
}

/** Output of context type in stream. */


/** Output of AccessKind in stream is supported and very useful.  */

} /* end namespace hmemo */

} /* end namespace scai */