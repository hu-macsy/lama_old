/**
 * @file HArray.hpp
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
 * @brief Definition of a dynamic array class where the array data can be
 *        used at different locations and where the data is moved implicitly
 * @author Thomas Brandes, Jiri Krause
 * @date 03.07.2015
 */

#pragma once

// base classes
#include <scai/hmemo/_HArray.hpp>

// common library
#include <scai/common/TypeTraits.hpp>
#include <scai/common/unique_ptr.hpp>

namespace scai
{

namespace hmemo
{

/**
 * @brief HArray is the typed version of _HArray.
 *
 * @tparam ValueType is the type stored in this container.
 *
 * HArray its contents on all supported locations, e.g. Host, CUDA and OpenCL. It transparently handles
 * synchronization between the locations. To enforce the consistency of the data a HArray can be only
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
class COMMON_DLL_IMPORTEXPORT HArray:

    public _HArray,
    public _HArray::Register<HArray<ValueType> >

{
public:

    typedef ValueType type;

    // WriteAccess and ReadAccess should be allowed to call some methods that
    // use ContextDataIndex for more efficient usage

    friend class ReadAccess<ValueType> ;
    friend class WriteAccess<ValueType> ;

    /**
     * @brief HArray() creates an empty HArray with size 0
     */
    HArray();

    /**
     * @brief Initialize an array with values from host
     */
    void init( const ValueType src[], const IndexType size );

    /**
     * @brief Initialize an array with same value for each entry
     */
    void init( const ValueType src, const IndexType size );

    /**
     * @brief Create a Heterogeneous array and give it a first touch on a context
     *
     * The first context decides about the default context and the default memory
     * used for the array. This is only a performance issue in some situations.
     */
    explicit HArray( ContextPtr context );

    /**
     * @brief Create a Heterogeneous array and give it a first touch on a memory.
     *
     * The first memory decides which data should be used whenever possible.
     */
    explicit HArray( MemoryPtr context );

    /**
     * @brief HArray( const IndexType n ) creates a HArray of size n
     *
     * @param[in] n the size of the HArray to create
     *
     * HArray( const IndexType n ) creates a HArray of size n and allocates uninitialized Host memory.
     */
    explicit HArray( const IndexType n );

    /**
     * @brief Create a HArray of size n and initialize it with one value.
     *
     * @param[in] n       the size of the HArray to create
     * @param[in] value   the value to initialize the container contens with
     * @param[in] context the location the HArray should be located first
     *
     * HArray( const IndexType n ) creates a HArray of size n, allocates Host memory and fills the Host memory with
     * the passed value.
     */
    HArray( const IndexType n, const ValueType& value, ContextPtr context = Context::getHostPtr() );

    /**
     * @brief Create a HArray of size n and initialize it with values from Host
     *
     * @param[in] n      the size of the HArray to create
     * @param[in] values host values
     * @param[in] context the location the HArray should be located first
     *
     * Note: the data is directly copied and the array can have no incarnation on Host
     */
    HArray( const IndexType n, const ValueType[], ContextPtr context = Context::getHostPtr() );

    /**
     * @brief Override the default copy constructor with appropriate version.
     *
     * @param[in] other the HArray to copy
     *
     * HArray(const HArray<ValueType>& other) copies the passed HArray. The container contens is copied for all currently valid
     * Locations.
     */
    HArray( const HArray<ValueType>& other );

    /**
     * @brief Destructor, releases all used resources.
     */
    virtual ~HArray();

    /**
     *  The method copy is a function that returns a new object of the
     *  same class as the object for which it is called. The copy
     *  constructor is called.
     */

    HArray<ValueType>* copy() const;

    /**
     *  Static create routine that is used for the _HArray factory.
     */

    static _HArray* create();

    /**
     * @brief Assignment operator for Heterogeneous arrays.
     *
     * @param[in] other the HArray to assign
     * @return this array as a copy of the other array
     *
     * The assignment operator copies the passed HArray.
     * The container content is copied for all contexts where a
     * valid copy is available (at least one).
     */
    HArray<ValueType>& operator=( const HArray<ValueType>& other );

    /**
     * @brief Assignment of array values with valid values at a given other.
     *
     * @param[in] other     the HArray whose values are copied
     * @param[in] context   the context where the assignment should be carried out
     *
     * The assign method copies the passed HArray.
     * The container content is copied for the passed contexts. If necessary other is
     * copied to context to carry this out.
     */
    void assign( const HArray<ValueType>& other, ContextPtr context );

    /**
     * @brief Swaps the contens of this with other.
     *
     * @param[in] other the HArray to swap the contens with.
     */
    void swapImpl( HArray<ValueType>& other );

    void swap( _HArray& other );

    /**
     * @brief sets the size of this to 0 an frees all memory
     */
    void purge();

    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Implementation of pure method.
     */
    virtual common::scalar::ScalarType getValueType() const;

    /**
     * @brief reserve a certain amount of data at a specific context
     *
     * @param[in] context where a certain amount of data should be reserved
     * @param[in] capacity amount of data to be allocated
     *
     */
    void reserve( ContextPtr context, const IndexType capacity );

    using _HArray::capacity;
    using _HArray::clear;

    /** Method that must be provided to guarantee a correct registration in
     *  the _HArray factory.
     */

    static common::scalar::ScalarType createValue()
    {
        return common::TypeTraits<ValueType>::stype;
    }

    static HArray<ValueType>* create( common::scalar::ScalarType key );

    using _HArray::resize;

protected:

    using _HArray::mSize;
    using _HArray::mValueSize;
    using _HArray::constFlag;

    ValueType* get( ContextDataIndex index );

    const ValueType* get( ContextDataIndex index ) const;

    void clear( ContextDataIndex index );

    void resize( ContextDataIndex index, const IndexType newSize );

    void reserve( ContextDataIndex index, const IndexType capacity ) const;

    IndexType capacity( ContextDataIndex index ) const;

};

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::HArray() :

    _HArray( 0, sizeof( ValueType ) )

{
    SCAI_LOG_DEBUG( logger, "created new HArray: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::HArray( ContextPtr context ) :

    _HArray( 0, sizeof( ValueType ) )

{
    // just make the first entry for the context
    /* ContextDataIndex data = */  mContextDataManager.getContextData( context );
    SCAI_LOG_DEBUG( logger, "created new HArray: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::HArray( MemoryPtr memory ) :
    _HArray( 0, sizeof( ValueType ) )
{
    // just make the first entry for the memory
    /* ContextDataIndex data = */  mContextDataManager.getMemoryData( memory );
    SCAI_LOG_DEBUG( logger, "created new HArray: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::HArray( const IndexType n ) :
    _HArray( n, sizeof( ValueType ) )
{
    // reserves already memory on the host, but this data is not valid
    ContextPtr hostPtr = Context::getHostPtr();
    mContextDataManager.reserve( hostPtr, n * mValueSize, 0 );
    SCAI_LOG_DEBUG( logger, "created new HArray: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::HArray( const IndexType n, const ValueType& value, ContextPtr context ) :

    _HArray( n, sizeof( ValueType ) )

{
    mContextDataManager.getContextData( context );  // first touch here
    // In constructor of the HArray lock of accesses is not required
    ContextPtr host = Context::getHostPtr();
    size_t validSize = 0;   // no valid data availalbe, so even don't search for it
    // Use of acquireAccess guarantees allocation of data
    ContextDataIndex index = mContextDataManager.acquireAccess( host, common::context::Write, mSize * mValueSize, validSize );
    ContextData& data = mContextDataManager[index];

    if ( n > 0 )
    {
        ValueType* hostData = static_cast<ValueType*>( data.get() );
        // ToDo: iterator changed from size_t to IndexType due to errors using compilers just supporting OpenMP 2.5
        #pragma omp parallel for
        for ( IndexType i = 0; i < mSize; ++i )
        {
            hostData[i] = value;
        }
    }

    releaseWriteAccess( index );

    SCAI_LOG_DEBUG( logger, "constructed: " << *this )

    // if array is constructed for other context than host, we prefetch the data to it

    if ( context->getType() != host->getType() )
    {
        prefetch( context );
    }
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::HArray( const IndexType n, const ValueType values[], ContextPtr context ) :

    _HArray( 0, sizeof( ValueType ) )

{
    /* ContextDataIndex data = */  mContextDataManager.getContextData( context );
    init( values, n );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::HArray( const HArray<ValueType>& other ):

    _HArray( other.mSize, sizeof( ValueType ) )

{
    mContextDataManager.copyAllValidEntries( other.mContextDataManager, mSize * mValueSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::~HArray()
{
    // destructor of ContextDataManager does all the release/check stuff
    SCAI_LOG_DEBUG( logger, "~HArray = " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
_HArray* HArray<ValueType>::create()
{
    return new HArray<ValueType>();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>* HArray<ValueType>::copy() const
{
    return new HArray<ValueType>( *this );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
common::scalar::ScalarType HArray<ValueType>::getValueType() const
{
    // Note: this is implementation of the pure method of base class _HArray.
    return common::TypeTraits<ValueType>::stype;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::init( const ValueType src[], const IndexType size )
{
    // context manager copies the data to the first touch location
    mContextDataManager.init( src, sizeof( ValueType ) * size );
    mSize = size;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::init( const ValueType value, const IndexType size )
{
    common::scoped_array<ValueType> data( new ValueType[size] );

    for ( IndexType i = 0; i < size; ++i )
    {
        data[i] = value;
    }

    init( data.get(), size );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>& HArray<ValueType>::operator=( const HArray<ValueType>& other )
{
    SCAI_LOG_DEBUG( logger, other << " will be assigned to " << *this )

    if ( &other == this )
    {
        return *this;
    }

    mSize      = other.mSize;
    mValueSize = other.mValueSize;
    SCAI_ASSERT( !mContextDataManager.locked(), "assign to a locked array (read/write access)" )
    // ToDo: we might add an exception on same thread: only valid write location is copied
    SCAI_ASSERT( !other.mContextDataManager.locked( common::context::Write ), "assign of a write locked array" )
    mContextDataManager.invalidateAll();
    // Now the same stuff as in copy constructor
    mContextDataManager.copyAllValidEntries( other.mContextDataManager, mSize * mValueSize );
    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::assign( const HArray<ValueType>& other, ContextPtr context )
{
    SCAI_LOG_DEBUG( logger, other << " will be assigned to " << *this )

    if ( &other == this )
    {
        mContextDataManager.setValidData( context, mContextDataManager, mSize * mValueSize );
        return;
    }

    mSize      = other.mSize;
    mValueSize = other.mValueSize;
    SCAI_ASSERT( !mContextDataManager.locked(), "assign to a locked array (read/write access)" )
    // ToDo: we might add an exception on same thread: only valid write location is copied
    SCAI_ASSERT( !other.mContextDataManager.locked( common::context::Write ), "assign of a write locked array" )
    mContextDataManager.invalidateAll();
    mContextDataManager.setValidData( context, other.mContextDataManager, mSize * mValueSize );
    SCAI_LOG_DEBUG( logger, *this << " has now been assigned at " << *context )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::swap( _HArray& other )
{
    SCAI_ASSERT_EQUAL( getValueType(), other.getValueType(), "swap only for same value type" )

    swapImpl( reinterpret_cast<HArray<ValueType>& >( other ) );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::swapImpl( HArray<ValueType>& other )
{
    SCAI_LOG_DEBUG( logger, *this << ": swap with other = " << other )
    // we cannot swap if there is any access for any array
    SCAI_ASSERT_EQUAL( 0, other.mContextDataManager.locked(), "swap: other array locked: " << other )
    SCAI_ASSERT_EQUAL( 0, mContextDataManager.locked(), "this array locked: " << *this )
    SCAI_ASSERT_EQUAL( mValueSize, other.mValueSize, "serious size mismatch" )
    mContextDataManager.swap( other.mContextDataManager );
    std::swap( mSize, other.mSize );
    SCAI_LOG_DEBUG( logger, *this << ": has been swapped with other = " << other )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::reserve( ContextPtr context, const IndexType capacity )
{
    mContextDataManager.reserve( context, capacity * mValueSize, mSize * mValueSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::purge()
{
    mContextDataManager.purge();
    mSize = 0;
    SCAI_LOG_DEBUG( logger, *this << " purged" )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType* HArray<ValueType>::get( ContextDataIndex index )
{
    return static_cast<ValueType*>( mContextDataManager[index].get() );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
const ValueType* HArray<ValueType>::get( ContextDataIndex index ) const
{
    return static_cast<const ValueType*>( mContextDataManager[index].get() );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::clear( const ContextDataIndex index )
{
    // make sure that we have exactly one write access at this context
    ContextData& data = mContextDataManager[index];
    SCAI_ASSERT_EQUAL( 1, mContextDataManager.locked( common::context::Write ), "multiple write access for clear" << data )
    SCAI_ASSERT_EQUAL( 0, mContextDataManager.locked( common::context::Read ), "further read access, cannot clear " << data )
    mSize = 0;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::resize( ContextDataIndex index, const IndexType size )
{
    ContextData& entry = mContextDataManager[index];
    bool inUse =  mContextDataManager.locked() > 1;   // further accesses on this array

    // SCAI_ASSERT( entry.locked( common::context::Write ), "resize illegal here " << entry )
    // static cast to have multiplication with 64 bit values

    size_t allocSize = static_cast<size_t>( size )  * static_cast<size_t>( mValueSize );
    size_t validSize = static_cast<size_t>( mSize ) * static_cast<size_t>( mValueSize );

    if ( validSize > allocSize )
    {
        validSize = allocSize;   // some entries are no more needed
    }

    SCAI_LOG_INFO( logger, *this << ": resize, needed = " << allocSize << " bytes, used = "
                   << validSize << " bytes, capacity = " << entry.capacity() << " bytes" )
    entry.reserve( allocSize, validSize, inUse );
    // capacity is now sufficient for size elements
    mSize = size;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::reserve( ContextDataIndex index, const IndexType size ) const
{
    if ( size <= mSize )
    {
        return;   // nothing to do
    }

    bool inUse =  mContextDataManager.locked() > 1;   // further accesses on this array
    ContextData& entry = mContextDataManager[index];

    size_t allocSize = static_cast<size_t>( size ) * static_cast<size_t>( mValueSize );
    size_t validSize = static_cast<size_t>( mSize ) * static_cast<size_t>( mValueSize );

    entry.reserve( allocSize, validSize, inUse );

    // Note: mSize does not change by the reserve
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
IndexType HArray<ValueType>::capacity( ContextDataIndex index ) const
{
    const ContextData& entry = mContextDataManager[index];
    return entry.capacity() / mValueSize;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "HArray<";
    stream << common::TypeTraits<ValueType>::id();
    stream << ">(" << mSize;
    stream << ") ";
    stream << mContextDataManager;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>* HArray<ValueType>::create( common::scalar::ScalarType key )
{
    if ( key == createValue() )
    {
        return reinterpret_cast<HArray<ValueType>* >( _HArray::create( key ) );
    }
    else
    {
        COMMON_THROWEXCEPTION( "creation not possible here" )
    }
}

/* ---------------------------------------------------------------------------------*/

} /* end namespace hmemo */

} /* end namespace scai */
