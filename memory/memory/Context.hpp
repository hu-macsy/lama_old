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

#include <memory/ContextData.hpp>

// for dll_import
#include <common/config.hpp>

// base classes
#include <common/Printable.hpp>
#include <common/NonCopyable.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

namespace memory
{

class SyncToken;
// forward declaration

class Context;
// forward declaration

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
class COMMON_DLL_IMPORTEXPORT Context: public Printable, private common::NonCopyable
{
public:

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
        MIC,           //!< Intel Many-Integrated-Core Architecture
        OpenCL,        //!< OpenCL GPU device, currently not supported
        UserContext,   //!< can be used for a new derived Context class
        UserContext1,  //!< ucan be used for a new derived Context class
        MaxContext     //!< used for dimension of ContextType arrays
    };

    virtual ~Context();

    /** Method to get the type of the context. */

    ContextType getType() const
    {
        return mContextType;
    }

    /** This predicate returns true if two contexts can use same data and
     *  have the same type.
     */
    bool operator==( const Context& other ) const;

    /** The inequality operator is just the inverse to operator==. */

    bool operator!=( const Context& other ) const;

    /** This predicate returns true if one context can use data of another context.
     *  This pure routine must be implemented by base classes.
     *
     *  @param[in] other is the context against which the check is done
     */
    virtual bool canUseData( const Context& other ) const = 0;

    virtual bool canCopyFrom( const Context& other ) const { return false; }

    virtual bool canCopyTo( const Context& other ) const { return false; }

    virtual void memcpyFrom( void* dst, const Context& srcContext, const void* src, size_t size ) const {}

    virtual void memcpyTo( const Context& dstContext, void* dst, const void* src, size_t size ) const {}

    virtual SyncToken* memcpyFromAsync( void* dst, const Context& srcContext, const void* src, size_t size ) const
    {
        memcpyFrom( dst, srcContext, src, size );
        return NULL;
    }

    virtual SyncToken* memcpyToAsync( const Context& dstContext, void* dst, const void* src, size_t size ) const 
    {
        memcpyTo( dstContext, dst, src, size );
        return NULL;
    }

    virtual void writeAt( std::ostream& stream ) const;

    /** This method allocates memory and must be implemented by
     *  each Context.
     *
     *  @param[in] size is the number of bytes needed
     *  @return pointer to the allocated data, NULL if not enough data is available
     */

    virtual void* allocate( const size_t size ) const = 0;

    /** This method allocates memory and must be implemented by
     *  each Context. The pointer to the allocated memory is stored in contextData.pointer.
     *
     *  @param[out] contextData   the ContextData the memory should be allocated for
     *  @param[in] size           is the number of bytes needed
     *  @return                   pointer to the allocated data,
     *                            NULL if not enough data is available
     */
    virtual void allocate( ContextData& contextData, const size_t size ) const = 0;

    /**
     * This method free's allocated data allocated by this allocator.
     *
     * @param[in] pointer is the pointer to the allocated data.
     * @param[in] size is the number of bytes that have been allocated with pointer.
     *
     * The pointer must have been allocated by the same allocator with the given size.
     */
    virtual void free( void* pointer, const size_t size ) const = 0;

    /**
     * This method free's allocated data allocated by this allocator.
     *
     * @param[in] contextData which holds the pointer to the allocated data.
     *
     * The pointer must have been allocated by the same allocator.
     */
    virtual void free( ContextData& contextData ) const = 0;

    /** Memory copy within the same context.
     *
     * param[in] dst pointer to the destination
     * param[in] src pointer to the source
     * param[in] size is the number of bytes to be copied.
     *
     * This memory copies size values.
     *
     */
    virtual void memcpy( void* dst, const void* src, const size_t size ) const = 0;

    /** Asynchronous memory copy within the same context.
     *
     * param[in] dst pointer to the destination
     * param[in] src pointer to the source
     * param[in] size is the number of bytes to be copied.
     * return new SyncToken object for the asynchronous operation.
     *
     * This memory copies size values.
     */
    virtual SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const = 0;

    /** Checks if this Context can copy from src to dst
     *
     * @param[in] dst   the dst ContextData
     * @param[in] src   the src ContextData
     * @return          if this Context can copy from src to dst.
     */
    virtual bool cancpy( const ContextData& dst, const ContextData& src ) const =0;

    /** Memory copy.
     *
     * Copies the first size bytes from src->pointer to dst->pointer. the memory pointed to by src
     * might be registered in a certain way to allow faster memory transfers.
     * However the memory pointed to by src->pointer will not be altered by this function.
     *
     * param[in] dst pointer to the destination
     * param[in] src pointer to the source
     * param[in] size is the number of bytes to be copied.
     *
     * This memory copies size values
     *
     */
    virtual void memcpy( ContextData& dst, const ContextData& src, const size_t size ) const = 0;

    virtual SyncToken* memcpyAsync( ContextData& dst, const ContextData& src, const size_t size ) const = 0;

    /** Getter routine for a new sync token that allows to asynchronous computations on the context.
     *
     *  @returns new SyncToken object
     */

    virtual SyncToken* getSyncToken() const = 0;

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

    const class LAMAInterface& getInterface() const;

    /** Get the preferred host context for a context.
     *
     *  Depending on the device, the host data might be allocated
     *  in such a way that faster and/or asynchronous memory transfer is supported.
     *
     *  As default, the default Host context is returned. Derived classes of
     *  Context should override this method to provide more efficient solutions.
     */

    virtual ContextPtr getHostContext() const;
 
    /** This method can be used to get a context of a certain type.
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
    static ContextPtr getContext( const Context::ContextType type, int deviceNr = -1 );

    /** Checks if a context of the passed type is available.
     *
     * @param[in] type  is the type of context that is wanted
     * @return          if a context of the passed type is available
     */
    static bool hasContext( const Context::ContextType type );

    /** Type definition of a function that delivers a context by a device id */

    typedef ContextPtr (*CreateFn)( int );

protected:

    /** Default constructor, can only be called by base classes. */

    Context( ContextType type );

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    ContextType mContextType;

    mutable bool mEnabled; //!<  if true the context is currently accessed

    mutable const char* mFile;//!< File name where context has been enabled

    mutable int mLine;//!< Line number where context has been enabled

    /** This method adds a Context creator routine for a certain type.
     *  By this way it is possible to create an object of a derived class 
     *  that might be added later dynamically.
     */

    static void addCreator( const Context::ContextType type,  CreateFn create );
};

/** Make ContextType visible in namespace, but not the different enumeration values. */

typedef Context::ContextType ContextType;

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const ContextType& type );

}
