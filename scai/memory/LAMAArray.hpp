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
 * @brief Definition of a dynamic array class where the array data can be
 *        used at different locations and where the data is moved implicitly
 *        when corresponding read/write accesses are required
 *
 * @author Thomas Brandes, Jiri Krause
 * @date 03.07.2015
 */

#pragma once

#include <scai/memory/ContextArray.hpp>

namespace scai
{

namespace memory
{

using common::IndexType;   // for convenience

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
class COMMON_DLL_IMPORTEXPORT LAMAArray: 

    public ContextArray,
    public ContextArray::Register<LAMAArray<ValueType> >

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
     * @brief Create a LAMA array and give it a first touch on a context
     *
     * The first context decides about the default context and the default memory          
     * used for the array. This is only a performance issue in some situations.
     */
    explicit LAMAArray( ContextPtr context );

    /**
     * @brief Create a LAMA array and give it a first touch on a memory.
     *
     * The first memory decides which data should be used whenever possible.
     */
    explicit LAMAArray( MemoryPtr context );

    /**
     * @brief LAMAArray( const IndexType n ) creates a LAMAArray of size n
     *
     * @param[in] n the size of the LAMAArray to create
     *
     * LAMAArray( const IndexType n ) creates a LAMAArray of size n and allocates uninitialized Host memory.
     */
    explicit LAMAArray( const common::IndexType n );

    /**
     * @brief Creates a LAMAArray of size n.
     *
     * @param[in] n     the size of the LAMAArray to create
     * @param[in] value the value to initialize the container contens with
     *
     * LAMAArray( const IndexType n ) creates a LAMAArray of size n, allocates Host memory and fills the Host memory with
     * the passed value.
     */
    LAMAArray( const common::IndexType n, const ValueType& value );

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
    LAMAArray( const common::IndexType n, const OtherValueType* const values );

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
     *  The method clone is a function that returns a new object of the
     *  same class as the object for which it is called. The default 
     *  constructor is always called.
     */

    LAMAArray<ValueType>* clone();

    /**
     *  Similiar to clone but here the copy constructor is called.
     */

    LAMAArray<ValueType>* copy();

    /** 
     *  Static create routine that is used for the ContextArray factory.
     */

    static ContextArray* create();

    /**
     * @brief Assignment operator for LAMA arrays.  
     *
     * @param[in] other the LAMAArray to assign
     * @return this array as a copy of the other array
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
     * @brief sets the size of this to 0 an frees all memory
     */
    void purge();

    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Implementation of pure method.
     */
    virtual common::ScalarType getValueType() const;

    /**
     * @brief reserve a certain amount of data at a specific context
     *
     * @param[in] context where a certain amount of data should be reserved
     * @param[in] capacity amount of data to be allocated
     *
     */
    void reserve( ContextPtr context, const common::IndexType capacity );

    using ContextArray::capacity;
    using ContextArray::clear;

    /** Method that must be provided to guarantee a correct registration in
     *  the ContextArray factory. 
     */

    static common::ScalarType createValue() 
    {
        return common::getScalarType<ValueType>();
    }

protected:

    using ContextArray::mSize;
    using ContextArray::mValueSize;
    using ContextArray::constFlag;

    ValueType* get( ContextDataIndex index );

    const ValueType* get( ContextDataIndex index ) const;

    void clear( ContextDataIndex index );

    void resize( ContextDataIndex index, const common::IndexType newSize );

    void reserve( ContextDataIndex index, const common::IndexType capacity ) const;

    common::IndexType capacity( ContextDataIndex index ) const;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

};

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
template<typename OtherValueType>
LAMAArray<ValueType>::LAMAArray( const common::IndexType n, const OtherValueType* const values )
                : ContextArray( n, sizeof( ValueType ) )
{
    ContextPtr hostContextPtr = Context::getContextPtr( context::Host );

    ContextData& host = mContextDataManager[ hostContextPtr ];

    if ( n <= 0 )
    {
        SCAI_LOG_DEBUG( logger, "Zero-sized array with value constructed: " << *this )
        return;
    }

    host.allocate( mSize * sizeof(ValueType) );

    SCAI_LOG_DEBUG( logger, "constructed: " << *this )

    ValueType* hostData = static_cast<ValueType*>( host.get() );

#pragma omp parallel for

    for( int i = 0; i < mSize; ++i )
    {
        hostData[i] = static_cast<ValueType>( values[i] );
    }

    host.setValid( true );

    SCAI_LOG_DEBUG( logger, "constructed: " << *this )
}

/* ---------------------------------------------------------------------------------*/

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, LAMAArray<ValueType>::logger, "LAMAArray" )

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray() : 

    ContextArray( 0, sizeof( ValueType ) )

{
    SCAI_LOG_DEBUG( logger, "created new LAMA array: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( ContextPtr context ) :

    ContextArray( 0, sizeof( ValueType ) )

{
    // just make the first entry for the context

    /* ContextDataIndex data = */  mContextDataManager.getContextData( context );

    SCAI_LOG_DEBUG( logger, "created new LAMA array: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( MemoryPtr memory ) :

    ContextArray( 0, sizeof( ValueType ) )

{
    // just make the first entry for the memory

    /* ContextDataIndex data = */  mContextDataManager.getMemoryData( memory );

    SCAI_LOG_DEBUG( logger, "created new LAMA array: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const common::IndexType n ) : 

    ContextArray( n, sizeof( ValueType) )

{
    // reserves already memory on the host, but this data is not valid

    ContextPtr hostPtr = Context::getContextPtr( context::Host );
    mContextDataManager.reserve( hostPtr, n * mValueSize, 0 );

    SCAI_LOG_DEBUG( logger, "created new LAMA array: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const common::IndexType n, const ValueType& value ) : ContextArray( n, sizeof( ValueType ) )

{
    // In constructor of the LAMA array lock of accesses is not required 

    ContextPtr host = Context::getContextPtr( context::Host );

    size_t validSize = 0;   // no valid data availalbe, so even don't search for it

    // Use of acquireAccess guarantees allocation of data

    ContextDataIndex index = mContextDataManager.acquireAccess( host, context::Write, mSize * mValueSize, validSize );

    ContextData& data = mContextDataManager[index];

    if ( n > 0 )
    {
        ValueType* hostData = static_cast<ValueType*>( data.get() );

#pragma omp parallel for 

        for ( int i = 0; i < mSize; ++i )
        {
            hostData[i] = value;
        }
    }

    releaseWriteAccess( index );

    SCAI_LOG_DEBUG( logger, "constructed: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::LAMAArray( const LAMAArray<ValueType>& other ): 

    ContextArray( other.mSize, sizeof( ValueType ) )

{
    mContextDataManager.copyAllValidEntries( other.mContextDataManager, mSize * mValueSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>::~LAMAArray()
{
    // destructor of ContextDataManager does all the release/check stuff

    SCAI_LOG_DEBUG( logger, "~LAMAArray = " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ContextArray* LAMAArray<ValueType>::create()
{
    return new LAMAArray<ValueType>();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>* LAMAArray<ValueType>::clone()
{
    return new LAMAArray<ValueType>();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>* LAMAArray<ValueType>::copy()
{
    return new LAMAArray<ValueType>( *this );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
common::ScalarType LAMAArray<ValueType>::getValueType() const
{
    // Note: this is implementation of the pure method of base class ContextArray.

    return common::getScalarType<ValueType>();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArray<ValueType>& LAMAArray<ValueType>::operator=( const LAMAArray<ValueType>& other )
{
    SCAI_LOG_DEBUG( logger, other << " will be assigned to " << *this )

    if ( &other == this )
    {
        return *this;
    }

    mSize      = other.mSize;
    mValueSize = other.mValueSize;

    COMMON_ASSERT( !mContextDataManager.locked(), "assign to a locked array (read/write access)" )

    // ToDo: we might add an exception on same thread: only valid write location is copied

    COMMON_ASSERT( !other.mContextDataManager.locked( context::Write ), "assign of a write locked array" )

    mContextDataManager.invalidateAll();

    // Now the same stuff as in copy constructor

    mContextDataManager.copyAllValidEntries( other.mContextDataManager, mSize * mValueSize );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::assign( const LAMAArray<ValueType>& other, ContextPtr context )
{
    SCAI_LOG_DEBUG( logger, other << " will be assigned to " << *this )

    if ( &other == this )
    {
         mContextDataManager.setValidData( context, mContextDataManager, mSize * mValueSize );
         return;
    }

    mSize      = other.mSize;
    mValueSize = other.mValueSize;

    COMMON_ASSERT( !mContextDataManager.locked(), "assign to a locked array (read/write access)" )

    // ToDo: we might add an exception on same thread: only valid write location is copied

    COMMON_ASSERT( !other.mContextDataManager.locked( context::Write ), "assign of a write locked array" )

    mContextDataManager.invalidateAll();

    mContextDataManager.setValidData( context, other.mContextDataManager, mSize * mValueSize );
 
    SCAI_LOG_DEBUG( logger, *this << " has now been assigned at " << *context )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::swap( LAMAArray<ValueType>& other )
{
    SCAI_LOG_DEBUG( logger, *this << ": swap with other = " << other )

    // we cannot swap if there is any access for any array

    COMMON_ASSERT_EQUAL( 0, other.mContextDataManager.locked(), "swap: other array locked: " << other )
    COMMON_ASSERT_EQUAL( 0, mContextDataManager.locked(), "this array locked: " << *this )

    COMMON_ASSERT_EQUAL( mValueSize, other.mValueSize, "serious size mismatch" )

    mContextDataManager.swap( other.mContextDataManager );

    std::swap( mSize, other.mSize );

    SCAI_LOG_DEBUG( logger, *this << ": has been swapped with other = " << other )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::reserve( ContextPtr context, const common::IndexType capacity )
{
    mContextDataManager.reserve( context, capacity * mValueSize, mSize * mValueSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::purge()
{
    mContextDataManager.purge();

    mSize = 0;
 
    SCAI_LOG_DEBUG( logger, *this << " purged" )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType* LAMAArray<ValueType>::get( ContextDataIndex index )
{
    return static_cast<ValueType*>( mContextDataManager[index].get() );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
const ValueType* LAMAArray<ValueType>::get( ContextDataIndex index ) const
{
    return static_cast<const ValueType*>( mContextDataManager[index].get() );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::clear( const ContextDataIndex index )
{
    // make sure that we have exactly one write access at this context

    ContextData& data = mContextDataManager[index];

    COMMON_ASSERT_EQUAL( 1, mContextDataManager.locked( context::Write ), "multiple write access for clear" << data )
    COMMON_ASSERT_EQUAL( 0, mContextDataManager.locked( context::Read ), "further read access, cannot clear " << data )

    mSize = 0;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::resize( ContextDataIndex index, const common::IndexType size )
{
    ContextData& data = mContextDataManager[index];

    COMMON_ASSERT_EQUAL( 1, mContextDataManager.locked( context::Write ), "multiple write access for resize" << data )
    COMMON_ASSERT_EQUAL( 0, mContextDataManager.locked( context::Read ), "further read access, cannot resize" << data )

    ContextData& entry = mContextDataManager[index];

    // COMMON_ASSERT( entry.locked( context::Write ), "resize illegal here " << entry )

    size_t allocSize = size * mValueSize;

    size_t validSize = mSize * mValueSize;

    if ( validSize > allocSize )
    {
        validSize = allocSize;   // some entries are no more needed
    }

    SCAI_LOG_INFO( logger, *this << ": resize, needed = " << allocSize << " bytes, used = " 
                         << validSize << " bytes, capacity = " << entry.capacity() << " bytes" )

    entry.reserve( allocSize, validSize );

    // capacity is now sufficient for size elements

    mSize = size;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::reserve( ContextDataIndex index, const common::IndexType size ) const
{
    if ( size <= mSize )
    {
        return;   // nothing to do
    }

    ContextData& entry = mContextDataManager[index];
   
    size_t allocSize = size * mValueSize;
    size_t validSize = mSize * mValueSize;

    entry.reserve( allocSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
common::IndexType LAMAArray<ValueType>::capacity( ContextDataIndex index ) const
{
    const ContextData& entry = mContextDataManager[index];
    return entry.capacity();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "LAMAArray<";
    stream << common::getScalarType<ValueType>();
    stream << ">(" << mSize; stream << ") ";
    stream << mContextDataManager;
}

/* ---------------------------------------------------------------------------------*/

} /* end namespace memory */

} /* end namespace scai */
