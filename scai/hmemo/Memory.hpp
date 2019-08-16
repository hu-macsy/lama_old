/**
 * @file Memory.hpp
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
 * @brief Definition of base class for memory management.
 * @author Thomas Brandes
 * @date 14.07.2015
 */
#pragma once

// local library
#include <scai/hmemo/Context.hpp>


// internal scai libraries

#include <scai/logging.hpp>

#include <scai/common/config.hpp>
#include <scai/common/Printable.hpp>
#include <scai/common/NonCopyable.hpp>

namespace scai
{

namespace tasking
{
class SyncToken;    // forward declaration
}

namespace hmemo
{

/** Enumeration type for the supported memory locations.
 *
 *  The same memory type does not imply that two Memory objects are the same.
 *  E.g. the memory of two different CUDA devices are not the same.
 */
enum class MemoryType
{
    HostMemory,       //!< memory for CPU as host, is main memory
    CUDAMemory,       //!< CUDA GPU memory on a device
    CUDAHostMemory,   //!< pinned memory that allows faster transfer to a certain CUDA Device
    UserMemory        //!< can be used for a new derived memory class
};

/**
 * This method make is possible to use enum values of MemoryType in output streams.
 */
COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const MemoryType& type );

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

    public  common::Printable,
    private common::NonCopyable
{
public:

    virtual ~Memory();

    /** Method to get the type of the memory. */

    inline MemoryType getType() const;

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

    /** This method allocates memory and must be implemented by each derived class.
     *
     *  @param[in] size is the number of bytes needed
     *  @return pointer to the allocated data, NULL if not enough data is available
     */

    virtual void* allocate( const size_t size ) = 0;

    /** This method returns the maximal number of allocated bytes during the lifetime of the memory. */

    inline size_t maxAllocatedBytes() const;

    /** This method returns the number of allocated bytes currently allocated by this memory. */

    inline size_t allocatedBytes() const;

    /** This method returns the number of allocates currently done on this memory */

    inline size_t allocates() const;

    /**
     * This method free's allocated data allocated by this allocator.
     *
     * @param[in] pointer is the pointer to the allocated data.
     * @param[in] size is the number of bytes that have been allocated with pointer.
     *
     * The pointer must have been allocated by the same allocator with the given size.
     */
    virtual void free( void* pointer, const size_t size ) = 0;

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

    /** Fill the memory with some byte value
     *
     * param[in] dst pointer to the destination
     * param[in] val will be interpreted as unsigned_char
     * param[in] size is the number of bytes to be filled
     *
     */
    virtual void memset( void* dst, const int val, const size_t size ) const = 0;

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

    /** @brief return the context to use for this memory (as reference) */

    inline const Context& getContext() const;

protected:

    /** Constructor, can only be called by derived classes. */

    Memory( MemoryType type );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Derived classes should call this method when they allocate data for correct statistics */

    void setAllocated( size_t nBytes );

    /** Derived classes should call this method when they allocate data for correct statistics */

    void setFreed( size_t nBytes );

    void checkAllFreed();

    Memory() = delete;  // disable default constructor

    Memory( const Memory& ) = delete;   // disable default copy constructor

    Memory& operator= ( const Memory& ) = delete;   // disable default assignment operator

    MemoryType mMemoryType;

    size_t mAllocates; //!< variable counts allocates

    size_t mAllocatedBytes;//!< variable counts allocated bytes

    size_t mMaxAllocatedBytes;//!< variable counts max allocated bytes
};

/* ------------------------------------------------------------------------- */
/*  Implementation of inline methods                                         */
/* ------------------------------------------------------------------------- */

MemoryType Memory::getType() const
{
    return mMemoryType;
}

size_t Memory::maxAllocatedBytes() const
{
    return mMaxAllocatedBytes;
}

size_t Memory::allocatedBytes() const
{
    return mAllocatedBytes;
}

size_t Memory::allocates() const
{
    return mAllocates;
}

const Context& Memory::getContext() const
{
    return *getContextPtr();
}

} /* end namespace hmemo */

} /* end namespace scai */
