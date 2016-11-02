/**
 * @file _HArray.hpp
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

/** Common base class for typed HArray.
 *
 *  Base class provides also a factory for creating arrays.
 */

class COMMON_DLL_IMPORTEXPORT _HArray:

    public common::Printable,
    public tasking::SyncTokenMember,
    public common::Factory<common::scalar::ScalarType, _HArray*>
{
    // Member variables of this class

protected:

    IndexType mSize;        //!< number of entries for the context array, common for all contexts
    IndexType mValueSize;   //!< number of bytes needed for one data element

    bool constFlag;         //!< if true the array cannot be written

    mutable ContextDataManager mContextDataManager;  //!< takes control of accesses and allocations

public:

    /** Virtual destructor required. */

    virtual ~_HArray()
    {
    }

    /**
     * @brief Query the value type of the array elements, e.g. DOUBLE or FLOAT.
     */
    virtual common::scalar::ScalarType getValueType() const = 0;

    using common::Factory<common::scalar::ScalarType, _HArray*>::create;

    /**
     *  Each derived class must provide a copy function. This will
     *  allow writing general routines that require temporary data.
     *
     *  Note: derived class might implement this routine by using covariant return types.
     *  Note: will be the same as _HArray::create( this->getValueType() )
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
    IndexType capacity( ContextPtr context ) const;

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
     *  common::unique_ptr<_Harray> arr1( _HArray::create( type ) );
     *  common::unique_ptr<_Harray> arr2( _HArray::create( type ) );
     *  ...
     *  arr1.swap( arr2 );   // is okay as they have same 'unknown' type
     *  \endcode
     */
    virtual void swap( _HArray& other ) = 0;

protected:

    explicit _HArray( const IndexType n, const IndexType size ) :

        mSize( n ),
        mValueSize( size ),
        constFlag( false )
    {
    }

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

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

    void releaseWriteAccess( ContextDataIndex );

    /** Query the capacity for a certain access.
     *
     *  @param[in] index is the reference to the context data as the result of an acquired access.
     */

    IndexType capacity( ContextDataIndex index ) const;

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

inline IndexType _HArray::capacity( ContextDataIndex index ) const
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

    return mContextDataManager.acquireAccess( context, common::context::Read, allocSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

inline void _HArray::releaseReadAccess( ContextDataIndex index ) const
{
    mContextDataManager.releaseAccess( index, common::context::Read );
}

/* ---------------------------------------------------------------------------------*/

inline ContextDataIndex _HArray::acquireWriteAccess( ContextPtr context, bool keepFlag )
{
    // use of size_t for bytes allows to allocate larger memory sizes

    size_t allocSize = static_cast<size_t>( mSize ) * static_cast<size_t>( mValueSize );

    size_t validSize = keepFlag ? allocSize : 0 ;    // valid data only if keepFlag is set
    return mContextDataManager.acquireAccess( context, common::context::Write, allocSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

inline void _HArray::releaseWriteAccess( ContextDataIndex index )
{
    mContextDataManager.releaseAccess( index, common::context::Write );
}

} /* end namespace hmemo */

} /* end namespace scai */
