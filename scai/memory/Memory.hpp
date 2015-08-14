/**
 * @file Memory.hpp
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
 * @brief Definition of base class for memory management.
 *
 * @author Thomas Brandes
 * @date 14.07.2015
 */
#pragma once

#include <scai/memory/Context.hpp>

#include <scai/common/config.hpp>
#include <scai/common/Printable.hpp>
#include <scai/common/NonCopyable.hpp>
#include <scai/common/shared_ptr.hpp>

// logging
#include <scai/logging.hpp>

namespace tasking
{
    class SyncToken;    // forward declaration
}

/** Namespace for all data structures of the context memory management. */

namespace memory
{

/** Namespace for enumeration of memory types */

namespace memtype
{

/** Enumeration type for the supported memory locations. 
 *
 *  The same memory type does not imply that two Memory objects are the same.
 *  E.g. the memory of two different CUDA devices are not the same.
 */
enum MemoryType
{
    HostMemory,       //!< memory for CPU as host, is main memory
    CUDAMemory,       //!< CUDA GPU memory on a device
    CUDAHostMemory,   //!< pinned memory that allows faster transfer to a certain CUDA Device
    UserMemory        //!< can be used for a new derived Context class
};
 
}

// make enumeration type visible but not the values

using memtype::MemoryType;

/** @brief This class is a common base class for all memory classes.
 *
 *  A memorys stands for memory management at a certain place.
 *
 *  Each memory must provide routines to allocate and free data as
 *  well as copy data within one memory class but also between different memory classes.
 *
 *  Especially each memory should provide routines to transfer data into the HostMem.
 *
 *  A copy constructor for a memory is not provided (singleton class).
 */
class COMMON_DLL_IMPORTEXPORT Memory: 
  
    public Printable, 
    private common::NonCopyable
{
public:

    virtual ~Memory();

    /** Method to get the type of the memory. */

    MemoryType getType() const;

    /** Predicate to check whether copy from other memory to this memory is supported. 
     *  If the method returns true, a call of memcpyFrom with srcMemory is safe. 
     *
     *  Note:  dstMemory.canCopyFrom( srcMemory ) and srcMemory.canCopyTo( dstMemory )
     *         can have different values, i.e. the corresponding memory transfer is only
     *         implemented by one memory.
     */

    virtual bool canCopyFrom( const Memory& srcMemory ) const;

    /** Predicate to check whether copy to other memory from this memory is supported. 
     *  If the method returns true, a call of memcpyTo with dstMemory is safe. 
     *
     *  Note:  dstMemory.canCopyFrom( srcMemory ) and srcMemory.canCopyTo( dstMemory )
     *         can have different values, i.e. the corresponding memory transfer is only
     *         implemented by one memory.
     */

    virtual bool canCopyTo( const Memory& dstMemory ) const;

    /** Copy from other memory to this memory. 
     *
     *  If canCopyFrom( srcMemory ) is false, this method throws an exception.
     */
    virtual void memcpyFrom( void* dst, const Memory& srcMemory, const void* src, size_t size ) const;

    /** Copy to other memory from this memory. 
     *
     *  If canCopyTo( dstMemory ) is false, this method throws an exception.
     */
    virtual void memcpyTo( const Memory& dstMemory, void* dst, const void* src, size_t size ) const;

    /** @brief Asynchronous version of memcpyFrom.
     *
     *  In this base class a default implementation is provided that does the memory transfer synchronously.
     */
    virtual tasking::SyncToken* memcpyFromAsync( void* dst, const Memory& srcMemory, const void* src, size_t size ) const;

    /** @brief Asynchronous version of memcpyTo.
     *
     *  In this base class a default implementation is provided that does the memory transfer synchronously.
     */
    virtual tasking::SyncToken* memcpyToAsync( const Memory& dstMemory, void* dst, const void* src, size_t size ) const;

    virtual void writeAt( std::ostream& stream ) const;

    /** This method allocates memory and must be implemented by
     *  each Memory.
     *
     *  @param[in] size is the number of bytes needed
     *  @return pointer to the allocated data, NULL if not enough data is available
     */

    virtual void* allocate( const size_t size ) const = 0;

    /**
     * This method free's allocated data allocated by this allocator.
     *
     * @param[in] pointer is the pointer to the allocated data.
     * @param[in] size is the number of bytes that have been allocated with pointer.
     *
     * The pointer must have been allocated by the same allocator with the given size.
     */
    virtual void free( void* pointer, const size_t size ) const = 0;

    /** Copy within the same memory.
     *
     * param[in] dst pointer to the destination
     * param[in] src pointer to the source
     * param[in] size is the number of bytes to be copied.
     *
     * This memory copies size values.
     *
     */
    virtual void memcpy( void* dst, const void* src, const size_t size ) const = 0;

    /** Asynchronous copy within the same memory.
     *
     * param[in] dst pointer to the destination
     * param[in] src pointer to the source
     * param[in] size is the number of bytes to be copied.
     * return new SyncToken object for the asynchronous operation.
     *
     * Base class provides default implementation where asynchronous copy is done synchronously.
     */
    virtual tasking::SyncToken* memcpyAsync( void* dst, const void* src, const size_t size ) const;

    /** Return a context at which memory can be used, e.g. to be initialized. 
     *
     *  This method must be implemented by all derived classes.
     *
     *  Hint: Often, Context and Memory come along together. Cyclic references with shared pointers
     *        should be avoided. 
     */

    virtual ContextPtr getContextPtr() const = 0;

    const Context& getContext() const 
    {
        return *getContextPtr();
    }

protected:

    /** Constructor, can only be called by derived classes. */

    Memory( MemoryType type );

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    MemoryType mMemoryType;

private:

    Memory();  // disable default constructor
};

inline MemoryType Memory::getType() const
{
    return mMemoryType;
}

}  // namespace

/** This method make is possible to use enum values of MemoryType in output streams. 
 *
 *  Note: It should not be defined in a namespace.
 */

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const memory::MemoryType& type );

