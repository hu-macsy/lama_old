/**
 * @file _HArray.hpp
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
 * @brief Common base class for the heterogeneous array where the data can be
 *        used in different contexts and where the data is moved implicitly
 * @author Thomas Brandes
 * @date 03.07.2015
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>
#include <scai/tasking/SyncToken.hpp>
#include <scai/common/Factory.hpp>

// local library
#include <scai/hmemo/Context.hpp>
#include <scai/hmemo/ContextDataManager.hpp>

// internal scai libraries
#include <scai/common/ScalarType.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/assert.hpp>

// std
#include <vector>
#include <map>

/** Number of contexts that might be used in maximum. This number
 *  is used for reservation of entries but does not imply any restrictions.
 */

#define SCAI_MAX_CONTEXTS 4

namespace scai
{

namespace hmemo
{

// Forward declaration of friend classes.

template<typename ValueType>
class ReadAccess;

template<typename ValueType>
class WriteAccess;

/** Common base class for typed HArray class.
 *
 *  This base class provides all methods for heterogeneous array in an untyped fashion,
 *  i.e. it knowns the size for one value entry but allocate, free, memory transfers are
 *  just in terms of bytes.
 *
 *  Base class provides also a factory for creating arrays.
 *  Base class derives also from SyncTokenMember, i.e. these arrays might be freed
 *  only after an asynchronous operation on it has finished.
 */
class COMMON_DLL_IMPORTEXPORT _HArray:

    public common::Printable,
    public tasking::SyncTokenMember,
    public common::Factory<common::ScalarType, _HArray*>
{

private:

    /* ==================================================================== */
    /*    member variables                                                  */
    /* ==================================================================== */

    IndexType mSize;              //!< number of entries for the context array, common for all contexts

    const IndexType mValueSize;   //!< number of bytes needed for one data element, can never be changed

    bool constFlag;               //!< if true the array cannot be written

    mutable ContextDataManager mContextDataManager;  //!< takes control of accesses and allocations

public:

    /** Virtual destructor required. */

    virtual ~_HArray()
    {
    }

    /** Just create an entry for managed locations. */

    void touch( ContextPtr ctx )
    {
        mContextDataManager.getContextData( ctx );
    }

    /** Just create an entry for certain memory. */

    void touch( MemoryPtr memory )
    {
        mContextDataManager.getMemoryData( memory );
    }

    /**
     * @brief reserve a certain amount of data at a specific context
     *
     * @param[in] context where a certain amount of data should be reserved
     * @param[in] capacity amount of data to be allocated
     *
     */
    inline void reserve( ContextPtr context, const IndexType capacity );

    void _setRawData( const IndexType size, const void* src )
    {
        mSize = size;

        // context manager copies the data to the first touch location

        mContextDataManager.init( src, mValueSize * size );
    }

    /** 
     *  @brief Get one value from the array at a certain pos
     *
     *  This operation uses memcpy to a valid context.
     */
    inline void _getValue( void* data, const IndexType pos ) const;

    /** 
     *  @brief Set one value from the array at a certain pos
     *
     *  This operation uses memcpy to a valid context.
     */
    inline void _setValue( const void* data, const IndexType pos );

    /**
     * @brief sets the size of this array to 0 an frees all memory
     */
    inline void purge();


    bool isConst() const { return constFlag; } 

    /**
     * @brief Query the value type of the array elements, e.g. DOUBLE or FLOAT.
     */
    virtual common::ScalarType getValueType() const = 0;

    using common::Factory<common::ScalarType, _HArray*>::create;

    /**
     *  @brief Create a new empty array with same value type.
     *
     *  Each derived class must provide a constructor function to create
     *  a new array of the same type. This allows writing general routines that 
     *  require temporary arrays of the same type.
     *
     *  Note: derived class might implement this routine by using covariant return types.
     *  Note: will be the same as _HArray::create( this->getValueType() )
     */

    virtual _HArray* newArray() const = 0;

    /**
     *  @brief Create a copy of this array.
     *
     *  Similiar to newArray, but here the copy constructor is called, i.e. the
     *  new array will have the same size and same values as the calling object.
     */

    virtual _HArray* copy() const = 0;

    /**
     * @brief Query the current size of the LAMA array, i.e. number of entries.
     *
     * @return the number of entries of the array.
     */
    inline IndexType size() const;

    /**
     * @brief Gets the first context where the data of this HArray is available.
     *
     * An argument can be passed to give a preferred context if the data is available
     * at multiple locations.
     *
     * @param[in] prefContext the preferred context to look for valid data
     * @return                a context there the data of this HArray is valid.
     *
     * Note: if the array has never been written to, no valid context is available.
     *       In this case this method returns getFirstTouchContextPtr()
     */
    inline ContextPtr getValidContext( const ContextPtr prefContext = ContextPtr() ) const;

    /**
     * @brief Get the context where the array has been touched the first time.
     *
     * Returns the context of the memory where the array has been used the first time.
     * For new arrays without entry, the Host context will be returned.
     *
     * This method can be used for Write only accesses where valid data is not required.
     */
    inline ContextPtr getFirstTouchContextPtr() const;

    /**
     * @brief Prefetches the content of the container to a certain location.
     *
     * @param[in] context  the location to prefetch to
     *
     * This method prefetches the content of the container to context.
     * If it is already valid at location nothing happens, otherwise a transfer from a valid location
     * to the passed location is started. Because the transfer is handled by HArray and to
     * maintain the consistency of the container only one running transfer can exist at any
     * point in time. There for if two prefetches to two different invalid locations are
     * started one after the other the second transfer does not start before the first one is finished.
     */
    void prefetch( ContextPtr context ) const;

    /** Waits for completion of prefetches. */

    void wait() const
    {
        mContextDataManager.wait();
    }

    /**
     * @brief Query the capacity ( in number of elements ) at a certain context.
     */
    inline IndexType capacity( ContextPtr context ) const;

    /**
     * @brief Query if data is valid in a certain context
     */
    bool isValid( ContextPtr context ) const;

    /**
     * @brief resize the array
     *
     * Changes the size of the array. No memory is freed if the size becomes smaller.
     * Reserves additional memory on all valid locations.
     */
    void resize( IndexType size );

    /**
     * @brief clear of an array is the same as resize 0
     *
     * \code
     *   _HArray& array = ...
     *   // this sequence of statements just invalidates all data
     *   IndexType size = array.size();
     *   array.clear();
     *   array.resize( 0 );
     * \endcode
     */

    void clear();

    /** Swap data with other array to avoid additional memory allocation.
     *
     *  @param[in,out] other array for swapping, must have same value type.
     *
     *  This method allows swapping for heterogeneous arrays where the value
     *  type is not known at compile time.
     *
     *  \code
     *  std::unique_ptr<_Harray> arr1( _HArray::create( type ) );
     *  std::unique_ptr<_Harray> arr2( _HArray::create( type ) );
     *  ...
     *  arr1.swap( arr2 );   // is okay as they have same 'unknown' type
     *  \endcode
     */
    void swap( _HArray& other );

    bool isInitialized() 
    {
        return mContextDataManager.isInitialized();
    }

    /** Method to override Printable::writeAt */

    void writeAt( std::ostream& stream ) const;

    /** Help method for writeAt of derived classes */

    void writeAtTyped( std::ostream& stream, const char* typeName ) const;

    /** This method initializes for an HArray reference.
     *
     *  @param[in] size is the number of entries 
     *  @param[in] pointer is the reference to the data
     */

    void setHostRef( const IndexType size, void* pointer );

    /** This method initializes for an HArray const reference.
     *
     *  @param[in] size is the number of entries 
     *  @param[in] pointer is the reference to the data
     */
    void setHostRef( const IndexType size, const void* pointer );

protected:

    /** copy constructor, only visible for derived classes. */

    _HArray( const _HArray& other );

    /** move constructor, only visible for derived classes. */

    _HArray( _HArray&& other ) noexcept;

    /** Assignment operator, only visible for derived classes. */

    _HArray& operator=( const _HArray& other );

    /** move assignment operator, only visible for derived classes
     *
     *  Very important: this move assignment must only be called between arrays
     *                  of the same value type, no type conversion supported here.
     */
    _HArray& operator=( _HArray&& other ) noexcept;


    /** Assignment, but here the target array will have only valid data at the 
     *  specified context. 
     */
    void assign( const _HArray& other, ContextPtr context );

    explicit _HArray( const IndexType n, const IndexType size ) :

        mSize( n ),
        mValueSize( size ),
        constFlag( false )
    {
    }

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** reserve of memory at a certain context/memory, i.e. pos is known 
     *
     *  @param[in] index must be legal index in the array on context data entries
     *  @param[in] capacity is number of data entries for which memory has to be reserved
     *
     *  Note: position in array is used instead of reference as array might be resized
     */
    void reserveWithIndex( ContextDataIndex index, const IndexType capacity ) const;

    /** resize of memory at a certain context/memory, i.e. pos is known */

    void resizeWithIndex( ContextDataIndex index, const IndexType capacity );

    /** Clear during a write access where index specifies current location. */

    void clearWithIndex( const ContextDataIndex index );

    /** Query the capacity for a certain access.
     *
     *  @param[in] index is the reference to the context data as the result of an acquired access.
     */

    inline IndexType capacityWithIndex( ContextDataIndex index ) const;

    /** Get pointer to the data at a certain location */

    void* get( ContextDataIndex index )
    {
        return mContextDataManager[index].get();
    }

    /** Get (read-only) pointer to the data at a certain location */

    const void* get( ContextDataIndex index ) const
    {
        return mContextDataManager[index].get();
    }

public:

    /** Complete handling to get read access for a certain context.
     *
     *  @param[in] context is the context where the read access is needed
     *
     *  For a read access it will be always guaranteed that a valid copy of the
     *  data is available at the context.
     *
     * \returns index of context data array that contains the valid entry.
     */
    ContextDataIndex acquireReadAccess( ContextPtr context ) const;

    /** Release an acquired read access. */

    void releaseReadAccess( ContextDataIndex ref ) const;

    /** Complete handling to get write access for a certain context.
     *
     *  @param[in] context is the context where the write access is needed
     *  @param[in] keepFlag if true it will be guaranteed that a valid copy is at the context.
     *
     *  A write only access will set keepFlag = false.
     *
     * \returns index of context data array that contains the valid entry.
     */
    ContextDataIndex acquireWriteAccess( ContextPtr context, bool keepFlag );

    /** Release an acquired write access. */

    inline void releaseWriteAccess( ContextDataIndex );

    /** Query the memory for a certain access. */

    const Memory& getMemory( ContextDataIndex index ) const;
};

/* ---------------------------------------------------------------------------------*/

typedef _HArray _HArray;

/* ---------------------------------------------------------------------------------*/

inline IndexType _HArray::size() const
{
    return mSize;
}

/* ---------------------------------------------------------------------------------*/

inline void _HArray::prefetch( ContextPtr context ) const
{
    size_t memSize = static_cast<size_t>( mSize ) * static_cast<size_t>( mValueSize );

    mContextDataManager.prefetch( context, memSize );
}

/* ---------------------------------------------------------------------------------*/

inline bool _HArray::isValid( ContextPtr context ) const
{
    return mContextDataManager.isValid( context );
}

/* ---------------------------------------------------------------------------------*/

inline void _HArray::resize( IndexType size )
{
    IndexType newSize = size;

    // resize on all valid locations

    if ( newSize  > mSize )
    {
        size_t oldMemSize = static_cast<size_t>( mSize ) * static_cast<size_t>( mValueSize );
        size_t newMemSize = static_cast<size_t>( newSize ) * static_cast<size_t>( mValueSize );

        mContextDataManager.resize( newMemSize, oldMemSize );
    }

    mSize = newSize;
}

/* ---------------------------------------------------------------------------------*/

void _HArray::reserve( ContextPtr context, const IndexType capacity )
{
    mContextDataManager.reserve( context, capacity * mValueSize, mSize * mValueSize );
}

/* ---------------------------------------------------------------------------------*/

void _HArray::purge()
{
    mContextDataManager.purge();
    mSize = 0;
}

/* ---------------------------------------------------------------------------------*/

inline void _HArray::clear()
{
    SCAI_ASSERT( !mContextDataManager.locked(), "Tried to clear a locked HArray " << *this )
    mSize = 0;
}

/* ---------------------------------------------------------------------------------*/

inline IndexType _HArray::capacity( ContextPtr context ) const
{
    // will return 0 if no data is available at the specified context
    return mContextDataManager.capacity( context ) / mValueSize;
}

/* ---------------------------------------------------------------------------------*/

IndexType _HArray::capacityWithIndex( ContextDataIndex index ) const
{
    const ContextData& entry = mContextDataManager[index];
    return static_cast<IndexType>( entry.capacity() / mValueSize );
}

/* ---------------------------------------------------------------------------------*/

inline const Memory& _HArray::getMemory( ContextDataIndex index ) const
{
    return mContextDataManager[ index ].getMemory();
}

/* ---------------------------------------------------------------------------------*/

inline ContextPtr _HArray::getValidContext( const ContextPtr prefContext ) const
{
    return mContextDataManager.getValidContext( prefContext );
}

/* ---------------------------------------------------------------------------------*/

inline ContextPtr _HArray::getFirstTouchContextPtr() const
{
    return mContextDataManager.getFirstTouchContextPtr();
}

/* ---------------------------------------------------------------------------------*/

inline ContextDataIndex _HArray::acquireReadAccess( ContextPtr context ) const
{
    // use of size_t for bytes allows to allocate larger memory sizes

    size_t allocSize = static_cast<size_t>( mSize ) * static_cast<size_t>( mValueSize );
    size_t validSize = allocSize;                   // read access needs valid data in any case

    return mContextDataManager.acquireAccess( context, common::AccessKind::Read, allocSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

inline void _HArray::releaseReadAccess( ContextDataIndex index ) const
{
    mContextDataManager.releaseAccess( index, common::AccessKind::Read );
}

/* ---------------------------------------------------------------------------------*/

inline ContextDataIndex _HArray::acquireWriteAccess( ContextPtr context, bool keepFlag )
{
    // use of size_t for bytes allows to allocate larger memory sizes

    size_t allocSize = static_cast<size_t>( mSize ) * static_cast<size_t>( mValueSize );

    size_t validSize = keepFlag ? allocSize : 0 ;    // valid data only if keepFlag is set
    return mContextDataManager.acquireAccess( context, common::AccessKind::Write, allocSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

void _HArray::releaseWriteAccess( ContextDataIndex index )
{
    mContextDataManager.releaseAccess( index, common::AccessKind::Write );
}

/* ---------------------------------------------------------------------------------*/

void _HArray::_getValue( void* data, const IndexType pos ) const
{
    mContextDataManager.getData( data, pos * mValueSize, mValueSize );
}

/* ---------------------------------------------------------------------------------*/

void _HArray::_setValue( const void* data, const IndexType pos )
{
    size_t allocSize = static_cast<size_t>( mSize ) * static_cast<size_t>( mValueSize );
    size_t offset    = static_cast<size_t>( pos ) * static_cast<size_t>( mValueSize );

    mContextDataManager.setData( data, offset, mValueSize, allocSize );
}

} /* end namespace hmemo */

} /* end namespace scai */
