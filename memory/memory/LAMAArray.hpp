/**
 * @file LAMAArray.hpp
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
 * @brief Definition of a new dynamic array class where the array data can be
 *        used in different contexts and where the data is moved implicitly
 *        when corresponding read/write accesses are required
 *
 * @author Thomas Brandes, Jiri Krause
 * @date 14.03.2011
 * @revised 03.07.2015
 */

#pragma once

// for dll_import
#include <common/config.hpp>

// base classes

// others
#include <memory/Context.hpp>
#include <memory/ContextManager.hpp>
#include <memory/SyncToken.hpp>
#include <memory/Scalar.hpp>

// common
#include <common/Printable.hpp>

#include <vector>

/** Number of contexts that might be used in maximum. This number
 *  is used for reservation of entries but does not imply any restrictions.
 */

#define LAMA_MAX_CONTEXTS 4

namespace memory
{

// Forward declaration of friend classes.

template<typename ValueType>
class ReadAccess;

template<typename ValueType>
class WriteAccess;

/** Common base class for typed LAMAArray. */

class COMMON_DLL_IMPORTEXPORT _LAMAArray: public Printable
{
    // Member variables of this class

protected:

    IndexType mSize;        //!< number of entries for the context array, common for all contexts
    IndexType mValueSize;   //!< number of bytes needed for one data element

    bool constFlag;         //!< if true the array cannot be written

    mutable ContextManager mContextManager;  //!< takes control of accesses and allocations

public:

    virtual ~_LAMAArray()
    {
    }

    /**
     * @brief Query the value type of the array elements, e.g. DOUBLE or FLOAT.
     */
    virtual Scalar::ScalarType getValueType() const = 0;

    /**
     * @brief Create an empty array of a certain type.
     *
     * @param type specifies the type of the array to be created.
     *
     */

    static _LAMAArray* create( const Scalar::ScalarType type );

    /**
     * @brief Query the current size of the LAMA array, i.e. number of entries.
     *
     * @return the number of entries of the array.
     */
    inline IndexType size() const;

    /**
     * @brief Gets the first context where the data of this LAMAArray is available.
     *
     * If possible a context of the passed preferred type is returned.
     *
     * @param[in] preferredType the preferred type for the valid context.
     * @return                  a context there the data of this LAMAArray is available.
     * 
     * Note: NULL pointer is returned if no valid data is available
     */
    ContextPtr getValidContext( const Context::ContextType preferredType = Context::Host ) const;

    /**
     * @brief Prefetches the contents of the container to the passed context.
     *
     * @param[in] context  the context to prefetch to
     *
     * This method prefetches the contents of the container to the context.
     * If this valid at location nothing happens,if not a transfer from a valid location
     * to the passed location is started. Because the transfer is handled by LAMAArray and to
     * maintain the consistency of the container only one running transfer can exist at any
     * point in time. There for if two prefetches to two different invalid locations are
     * started one after the other the second transfer does not start before the first one is finished.
     */
    void prefetch( ContextPtr context ) const;

    /**
     * @brief Query the capacity ( in number of elements ) at a certain context.
     */
    IndexType capacity( ContextPtr context ) const;
    
    /**
     * @brief Query if data is valid in a certain context
     */
    bool isValid( ContextPtr context ) const;
    
protected:

    explicit _LAMAArray( const IndexType n, const IndexType size )
                    : mSize( n ), mValueSize( size ), constFlag( false )
    {
        LAMA_LOG_DEBUG( logger, "construct LAMAArray, mSize = " << mSize << ", mValueSize " << mValueSize )
    }

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

/**
 * @brief LAMAArray is the base array container for all compute relevant data within LAMA.
 *
 * @tparam ValueType is the type stored in this container.
 *
 * LAMAArray its contents on all supported locations, e.g. Host, CUDA and OpenCL. It transparently handles
 * synchronization between the locations. To enforce the consistency of the data a LAMAArray can be only
 * indirectly accessed via a ReadAccess or a WriteAccess.
 *
 * Compared to a C++ container like std::vector some differences must be taken into account:
 * 
 *  - There is never a call of the default constructor, destructor or copy constructor of ValueType
 *    (so data is always handled bytewise in case of context transfers or reallocation)
 *  - Iterators are not provided
 *
 *  Even if ValueType is usally float or double, other data types might also be used
 *  (especially structures of such data). Do not use any ValueType that contains pointers
 *  or references; these might be invalid when data is moved to another context.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT LAMAArray: public _LAMAArray
{
public:

    // WriteAccess and ReadAccess should be allowed to call some methods that
    // use ContextDataIndex for more efficient usage

    friend class ReadAccess<ValueType> ;
    friend class WriteAccess<ValueType> ;

    /**
     * @brief LAMAArray() creates an empty LAMAArray with size 0
     */
    LAMAArray();

    /**
     * @brief LAMAArray( const IndexType n ) creates a LAMAArray of size n
     *
     * @param[in] n the size of the LAMAArray to create
     *
     * LAMAArray( const IndexType n ) creates a LAMAArray of size n and allocates uninitialized Host memory.
     */
    explicit LAMAArray( const IndexType n );

    /**
     * @brief Creates a LAMAArray of size n.
     *
     * @param[in] n     the size of the LAMAArray to create
     * @param[in] value the value to initialize the container contens with
     *
     * LAMAArray( const IndexType n ) creates a LAMAArray of size n, allocates Host memory and fills the Host memory with
     * the passed value.
     */
    LAMAArray( const IndexType n, const ValueType& value );

    /**
     * @brief Creates a LAMAArray of size n.
     *
     * @param[in] n         the size of the LAMAArray to create
     * @param[in] values    the values to initialize the container contens with
     *
     * LAMAArray( const IndexType n ) creates a LAMAArray of size n, allocates Host memory and fills the Host memory with
     * the passed values.
     */
    template<typename OtherValueType>
    LAMAArray( const IndexType n, const OtherValueType* const values );

    /**
     * @brief Creates a copy of the passed LAMAArray.
     *
     * @param[in] other the LAMAArray to copy
     *
     * LAMAArray(const LAMAArray<ValueType>& other) copies the passed LAMAArray. The container contens is copied for all currently valid
     * Locations.
     */
    LAMAArray( const LAMAArray<ValueType>& other );

    /**
     * @brief Destructor, releases all used resources.
     */
    virtual ~LAMAArray();

    /**
     * @brief Copies the passed LAMAArray into this.
     *
     * @param[in] other the LAMAArray to copy
     *
     * The assignment operator copies the passed LAMAArray.
     * The container content is copied for all contexts where a
     * valid copy is available (at least one).
     */
    LAMAArray<ValueType>& operator=( const LAMAArray<ValueType>& other );

    /**
     * @brief Assignment of array values with valid values at a given other.
     *
     * @param[in] other     the LAMAArray whose values are copied
     * @param[in] context   the context where the assignment should be carried out
     *
     * The assign method copies the passed LAMAArray.
     * The container content is copied for the passed contexts. If necessary other is
     * copied to context to carry this out.
     */
    void assign( const LAMAArray<ValueType>& other, ContextPtr context );

    /**
     * @brief Swaps the contens of this with other.
     *
     * @param[in] other the LAMAArray to swap the contens with.
     */
    void swap( LAMAArray<ValueType>& other );

    /**
     * @brief sets the size of this to 0 but does not free any memory
     *
     * Clear also invalidates all data. So this operation
     * should be used before a write-only access that is
     * followed by a resize of the array. It avoids data transfer between
     * two contextes.
     */
    void clear();

    /**
     * @brief sets the size of this to 0 an frees all memory
     */
    void purge();

    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Implementation of pure method.
     */
    virtual Scalar::ScalarType getValueType() const;

    /**
     * @brief reserve a certain amount of data at a specific context
     *
     * @param[in] context where a certain amount of data should be reserved
     * @param[in] capacity amount of data to be allocated
     *
     */
    void reserve( ContextPtr context, const IndexType capacity );

    using _LAMAArray::capacity;

protected:

    using _LAMAArray::mSize;
    using _LAMAArray::mValueSize;
    using _LAMAArray::constFlag;

    ValueType* get( ContextDataIndex index );

    const ValueType* get( ContextDataIndex index ) const;

    void clear( ContextDataIndex index );

    void resize( ContextDataIndex index, const IndexType newSize );

    void reserve( ContextDataIndex index, const IndexType capacity ) const;

    IndexType capacity( ContextDataIndex index ) const;

    /** Complete handling to get read access for a certain context.
     *
     * \returns index of context data array that contains the valid entry.
     */
    ContextDataIndex acquireReadAccess( ContextPtr context ) const;

    void releaseReadAccess( ContextDataIndex ref ) const;

    /** Complete handling to get write access for a certain context.
     *
     * \returns index of context data array that contains the valid entry.
     */
    ContextDataIndex acquireWriteAccess( ContextPtr context, bool keepFlag );

    void releaseWriteAccess( ContextDataIndex );

    void fetch( ContextData& target, const ContextData& source ) const;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

/**
 * @brief LAMAArrayRef is a container that uses external data.
 *
 * @tparam ValueType is the type stored in this container.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT LAMAArrayRef: public LAMAArray<ValueType>
{
public:

    /** Contruct a container for a host array. */

    LAMAArrayRef( ValueType* pointer, IndexType size );

    /** Contruct a container for a const host array.
     *  Due to the const pointer it is guaranteed that the array cannot be modified
     */

    LAMAArrayRef( const ValueType* pointer, IndexType size );

protected:

    using LAMAArray<ValueType>::mSize;

    using LAMAArray<ValueType>::mContextManager;
    using LAMAArray<ValueType>::constFlag;
};

/* ---------------------------------------------------------------------------------*/

inline IndexType _LAMAArray::size() const
{
    return mSize;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
template<typename OtherValueType>
LAMAArray<ValueType>::LAMAArray( const IndexType n, const OtherValueType* const values )
                : _LAMAArray( n, sizeof( ValueType ) )
{
    ContextPtr hostContext = Context::getContext( Context::Host );

    ContextData& host = mContextManager[ hostContext ];

    if( n <= 0 )
    {
        LAMA_LOG_DEBUG( logger, "Zero-sized array with value constructed: " << *this )
        return;
    }

    host.allocate( mSize * sizeof(ValueType) );

    LAMA_LOG_DEBUG( logger, "constructed: " << *this )

    ValueType* hostData = static_cast<ValueType*>( host.get() );

#pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

    for( int i = 0; i < mSize; ++i )
    {
        hostData[i] = static_cast<ValueType>( values[i] );
    }

    host.setValid( true );

    LAMA_LOG_DEBUG( logger, "constructed: " << *this )
}

/* ---------------------------------------------------------------------------------*/

void _LAMAArray::prefetch( ContextPtr context ) const
{
    mContextManager.prefetch( context, mSize * mValueSize );
}

/* ---------------------------------------------------------------------------------*/

bool _LAMAArray::isValid( ContextPtr context ) const
{
    return mContextManager.isValid( context );
}

/* ---------------------------------------------------------------------------------*/

IndexType _LAMAArray::capacity( ContextPtr context ) const
{
    return mContextManager.capacity( context ) / mValueSize;
}

/* ---------------------------------------------------------------------------------*/

ContextPtr _LAMAArray::getValidContext( const Context::ContextType preferredType ) const
{
    return mContextManager.getValidContext( preferredType );
}

typedef ContextData::AccessKind AccessKind;

/* ---------------------------------------------------------------------------------*/

LAMA_LOG_DEF_LOGGER( _LAMAArray::logger, "LAMAArray" )

/* -------------------------------------------------------------------------- */

_LAMAArray* _LAMAArray::create( const Scalar::ScalarType valueType )
{
    switch( valueType )
    {
        case Scalar::FLOAT:
        {
            return new LAMAArray<float>();
        }

        case Scalar::DOUBLE:
        {
            return new LAMAArray<double>();
        }

        default:
        {
            COMMON_THROWEXCEPTION( "Unsupported ValueType " << valueType )
        }
    }
}

/* ---------------------------------------------------------------------------------*/

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, LAMAArray<ValueType>::logger, "LAMAArray" )

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray()
                : _LAMAArray( 0, sizeof( ValueType ) )
{
    // just make an entry for host context

    ContextPtr host = Context::getContext( Context::Host );

    ContextDataIndex data = mContextManager.getContextData( host );

    LAMA_LOG_DEBUG( logger, "created new LAMA array, mSize = " << mSize )
    LAMA_LOG_DEBUG( logger, "created new LAMA array, mValueSize = " << mValueSize )
    
    // LAMA_LOG_DEBUG( logger, "created new context array: " << *this )
}

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const IndexType n ) : _LAMAArray( n, sizeof( ValueType) )
{
    ContextPtr host = Context::getContext( Context::Host );

    // just allocate the data at host context with a corresponding write access

    size_t validSize = 0;   // no valid data availalbe, so even don't search for it

    // Use of acquireAccess guarantees allocation of data

    ContextDataIndex index = mContextManager.acquireAccess( host, ContextData::Write, mSize * mValueSize, validSize );

    releaseWriteAccess( index );
}

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const IndexType n, const ValueType& value ) : _LAMAArray( n, sizeof( ValueType ) )

{
    // In constructor of the LAMA array lock of accesses is not required 

    ContextPtr host = Context::getContext( Context::Host );

    size_t validSize = 0;   // no valid data availalbe, so even don't search for it

    // Use of acquireAccess guarantees allocation of data

    ContextDataIndex index = mContextManager.acquireAccess( host, ContextData::Write, mSize * mValueSize, validSize );

    ContextData& data = mContextManager[index];

    if ( n > 0 )
    {
        ValueType* hostData = static_cast<ValueType*>( data.get() );

#pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)

        for ( int i = 0; i < mSize; ++i )
        {
            hostData[i] = value;
        }
    }

    releaseWriteAccess( index );

    LAMA_LOG_DEBUG( logger, "constructed: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const LAMAArray<ValueType>& other )

    : _LAMAArray( other.mSize, sizeof( ValueType ) )
{
    mContextManager.copyAllValidEntries( other.mContextManager, mSize * mValueSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::~LAMAArray()
{
    // destructor of ContextManager does all the release/check stuff

    LAMA_LOG_DEBUG( logger, "~LAMAArray = " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
Scalar::ScalarType LAMAArray<ValueType>::getValueType() const
{
    // Note: this is implementation of the pure method of base class _LAMAArray.

    return Scalar::getType<ValueType>();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>& LAMAArray<ValueType>::operator=( const LAMAArray<ValueType>& other )
{
    LAMA_LOG_DEBUG( logger, other << " will be assigned to " << *this )

    if ( other == *this )
    {
        return *this;
    }

    mSize      = other.mSize;
    mValueSize = other.mValueSize;

    COMMON_ASSERT( !mContextManager.locked(), "assign to a locked array (read/write access)" )

    // ToDo: we might add an exception on same thread: only valid write location is copied

    COMMON_ASSERT( !other.mContextManager.writeLocked(), "assign of a write locked array" )

    mContextManager.invalidateAll();

    // Now the same stuff as in copy constructor

    mContextManager.copyAllValidEntries( other.mContextManager, mSize * mValueSize );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::assign( const LAMAArray<ValueType>& other, ContextPtr context )
{
    LAMA_LOG_DEBUG( logger, other << " will be assigned to " << *this )

    if ( other == *this )
    {
         mContextManager.setValidData( context, mContextManager, mSize * mValueSize );
         return;
    }

    mSize      = other.mSize;
    mValueSize = other.mValueSize;

    COMMON_ASSERT( !mContextManager.locked(), "assign to a locked array (read/write access)" )

    // ToDo: we might add an exception on same thread: only valid write location is copied

    COMMON_ASSERT( !other.mContextManager.writeLocked(), "assign of a write locked array" )

    mContextManager.invalidateAll();

    mContextManager.setValidData( context, other.mContextManager, mSize * mValueSize );
 
    LAMA_LOG_DEBUG( logger, *this << " has now been assigned at " << *context )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::swap( LAMAArray<ValueType>& other )
{
    LAMA_LOG_DEBUG( logger, *this << ": swap with other = " << other )

    // we cannot swap if there is any access for any array

    COMMON_ASSERT( !other.mContextManager.locked(), "swap: other array locked" )
    COMMON_ASSERT( mContextManager.locked(), "this array locked" )

    COMMON_ASSERT_EQUAL( mValueSize, other.mValueSize, "serious size mismatch" )

    std::swap( mSize, other.mSize );

    mContextManager.swap( other.mContextManager );

    LAMA_LOG_DEBUG( logger, *this << ": has been swapped with other = " << other )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::reserve( ContextPtr context, const IndexType capacity )
{
    mContextManager.reserve( context, capacity * mValueSize, mSize * mValueSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::clear()
{
    COMMON_ASSERT( !mContextManager.locked(), "Tried to clear a locked LAMAArray " << *this )

    mSize = 0;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::purge()
{
    mContextManager.purge();

    mSize = 0;
 
    LAMA_LOG_DEBUG( logger, *this << " purged" )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType* LAMAArray<ValueType>::get( ContextDataIndex index )
{
    return static_cast<ValueType*>( mContextManager[index].get() );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
const ValueType* LAMAArray<ValueType>::get( ContextDataIndex index ) const
{
    return static_cast<const ValueType*>( mContextManager[index].get() );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::clear( const ContextDataIndex index )
{
    ContextData& data = mContextManager[index];

    // ToDo: COMMON_ASSERT( data.locked( ContextData::Write ), "clear illegal here " << data )

    mSize = 0;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::resize( ContextDataIndex index, const IndexType size )
{
    ContextData& entry = mContextManager[index];

    // COMMON_ASSERT( entry.locked( ContextData::Write ), "resize illegal here " << entry )

    size_t allocSize = size * mValueSize;

    size_t validSize = mSize * mValueSize;

    if ( validSize > allocSize )
    {
        validSize = allocSize;   // some entries are no more needed
    }

    LAMA_LOG_INFO( logger, *this << ": resize, needed = " << allocSize << " bytes, used = " 
                         << validSize << " bytes, capacity = " << entry.capacity() << " bytes" )

    entry.reserve( allocSize, validSize );

    // capacity is now sufficient for size elements

    mSize = size;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::reserve( ContextDataIndex index, const IndexType size ) const
{
    COMMON_THROWEXCEPTION( "not available yet" )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
IndexType LAMAArray<ValueType>::capacity( ContextDataIndex index ) const
{
    ContextData& entry = mContextManager[index];

    return static_cast<IndexType>( entry.capacity() / sizeof( ValueType ) );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ContextDataIndex LAMAArray<ValueType>::acquireReadAccess( ContextPtr context ) const
{
    LAMA_LOG_DEBUG( logger, "acquireReadAccess for " << *this );

    size_t allocSize = mSize * mValueSize;
    size_t validSize = allocSize;                   // read access needs valid data in any case

    return mContextManager.acquireAccess( context, ContextData::Read, allocSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::releaseReadAccess( ContextDataIndex index ) const
{
    mContextManager.releaseAccess( index, ContextData::Read );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ContextDataIndex LAMAArray<ValueType>::acquireWriteAccess( ContextPtr context, bool keepFlag )
{
    size_t allocSize = mSize * sizeof( ValueType );
    size_t validSize = keepFlag ? allocSize : 0 ;

    return mContextManager.acquireAccess( context, ContextData::Write, allocSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::releaseWriteAccess( ContextDataIndex index )
{
    ContextData& data = mContextManager[index];

    mContextManager.releaseAccess( index, ContextData::Write );

    LAMA_LOG_INFO( logger, "releaseWriteAccess: " << data );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "LAMAArray<";
    stream << Scalar::getType<ValueType>();
    stream << ">(" << mSize; stream << ")";
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArrayRef<ValueType>::LAMAArrayRef( ValueType* pointer, IndexType size )
                : LAMAArray<ValueType>()
{
    // Important: context must be set to the DefaultHostContext

    if( size != 0 && pointer == NULL )
    {
        COMMON_THROWEXCEPTION( "LAMAArryRef with NULL pointer" )
    }

/*
    ContextData& host = *mContextData[0];
    host.setRef( pointer, size * sizeof(ValueType) );
*/

    mSize = size;
}

template<typename ValueType>
LAMAArrayRef<ValueType>::LAMAArrayRef( const ValueType* pointer, IndexType size )
                : LAMAArray<ValueType>()
{
    // Important: context must be set to the DefaultHostContext

    if( size != 0 && pointer == NULL )
    {
        COMMON_THROWEXCEPTION( "LAMAArryRef with NULL pointer" )
    }

/*
    ContextData& host = *mContextData[0];
    host.setRef( const_cast<ValueType*>( pointer ), size * sizeof(ValueType) );

    constFlag = true; // makes sure that we cannot have a WriteAccess

    mSize = size;
*/

}

}  // namespace 
