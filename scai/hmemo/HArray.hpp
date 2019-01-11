/**
 * @file HArray.hpp
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

#include <memory>
#include <initializer_list>

namespace scai
{

namespace hmemo
{

/**
 * @brief HArray is the typed version of the container class _HArray (heterogeneous array).
 *
 * @tparam ValueType is the type stored in this container.
 *
 * HArray can have incarnations on multiple locations, e.g. Host, CUDA and OpenCL. It transparently handles
 * synchronization between the locations. To enforce the consistency of the data a HArray can be only
 * indirectly accessed via a ReadAccess or a WriteAccess.
 *
 * Compared to a C++ container like std::vector some differences must be taken into account:
 *
 *  - There is never a call of the default constructor, destructor or copy constructor of the 
 *    array elements (so data is always handled bytewise in case of context transfers or reallocation)
 *  - Iterators are not provided
 *
 *  Even if ValueType is usally float or double, other data types might also be used
 *  (especially structures of such data). Do not use any ValueType that contains pointers
 *  or references; these might be invalid when data is moved to another context.
 *
 *  Move constructor, move assignment are also provided for this (template) class. This allows
 *  using the arrays in the new C++ way to avoid copies. Especially, heterogeneous 
 *  arrays might be used in C++ container classes without allocation/reallocation and/or copying of data.
 *
 *  In contrary to other LAMA classes, this template class is a header-only class, so it
 *  can be used in applications for any value type.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT HArray:

    public _HArray,
    public _HArray::Register<HArray<ValueType> >

{

public:

    /** Help class to observe the further use of operator[] in HArray */

    class IndexProxy
    {
    public:

        /** Proxy constructed by ref to the array and the index value. */

        IndexProxy( HArray<ValueType>& array, const IndexType i ) :

            mArray( array ),
            mIndex( i )
        {
        }

        /** indexed value proxy can be used to get its value */

        operator ValueType() const
        {
            ValueType val;
            mArray._getValue( &val, mIndex );
            return val;
        }

        /** indexed value proxy can be assigned a value */

        IndexProxy& operator= ( ValueType val )
        {
            mArray._setValue( &val, mIndex );
            return *this;
        }

        /** Override the default assignment operator to avoid ambiguous interpretation of a[i] = b[i] */

        IndexProxy& operator= ( const IndexProxy& other )
        {
            ValueType tmp;
            other.mArray._getValue( &tmp, other.mIndex );
            mArray._setValue( &tmp, mIndex );
            return *this;
        }

    private:

        HArray<ValueType>& mArray;
        IndexType mIndex;

    };

    // WriteAccess and ReadAccess should be allowed to call some methods that
    // use ContextDataIndex for more efficient usage

    friend class ReadAccess<ValueType> ;
    friend class WriteAccess<ValueType> ;

    /* ======================================================== */
    /*   Constructor / Destructor                               */
    /* ======================================================== */

    /**
     * @brief HArray() creates an empty HArray with size 0
     */
    HArray();

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
     * @brief Create an (unitialized) HArray of size n 
     *
     * @param[in] n       the size of the HArray to create
     * @param[in] context the location the HArray should be located first
     *
     */
    HArray( const IndexType n, ContextPtr context );

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
     * @brief Override the default move constructor with appropriate version.
     *
     * @param[in] other the HArray from which data is moved
     *
     * The new created object gets its resources from the passed array.
     */
    HArray( HArray<ValueType>&& other ) noexcept;

    /**
     * @brief Construct an instance of HArray from an initializer list.
     *
     * @param[in] init an initializer list of data
     * @param[in] context the context for which to allocate memory to hold the data (optional)
     *
     */
    HArray( std::initializer_list<ValueType> init, ContextPtr context = Context::getHostPtr() );

    /**
     * @brief Construct an instance of HArray form a vector 
     */
    HArray( const std::vector<ValueType>& vector, ContextPtr context = Context::getHostPtr() );

    /**
     * @brief Destructor, releases all used resources.
     */
    virtual ~HArray();

    /* ======================================================== */
    /*   operator []                                            */
    /* ======================================================== */

    IndexProxy operator[] ( const IndexType i )
    {
        return IndexProxy( *this, i );
    }

    ValueType operator[] ( const IndexType i ) const
    {
        ValueType val;
        _HArray::_getValue( &val, i );
        return val;
    }

    /* ======================================================== */
    /*   Setter methods                                         */
    /* ======================================================== */

    /**
     * @brief Initialize an array with raw data values from host 
     *
     * @param[in] size is the number of elements to set
     * @param[in] data is the value array with at least size values
     */
    void setRawData( const IndexType size, const ValueType data[] );

    /**
     * @brief Initialize an array with same value for each entry
     */
    void setSameValue( const IndexType size, const ValueType value );

    /* ======================================================== */
    /*   Dynamic creators                                       */
    /* ======================================================== */

    /**
     *  The method copy is a function that returns a new object of the
     *  same class as the object for which it is called. The copy
     *  constructor is called.
     */
    HArray<ValueType>* newArray() const;

    /**
     *  Implementation of pure methods _HArray::copy(), returns covariant type.
     */
    HArray<ValueType>* copy() const;

    /**
     *  Static create routine that is used for the _HArray factory.
     */
    static _HArray* create();

    /** Get the key to create an heterogeneous array of this type via the 
     *  factory. 
     *
     *  Method is mandatory to guarantee a correct registration in
     *  the _HArray factory.
     */
    static common::ScalarType createValue()
    {
        return common::TypeTraits<ValueType>::stype;
    }

    /* ======================================================== */
    /*   Assignment operators                                   */
    /* ======================================================== */

    /**
     * @brief Assignment operator for Heterogeneous arrays.
     *
     * @param[in] other the HArray to assign
     * @return this array as a copy of the other array
     *
     * The assignment operator copies the passed array, i.e. as a deep copy.
     * The container content is copied for all contexts where a
     * valid copy is available (at least one).
     */
    HArray<ValueType>& operator=( const HArray<ValueType>& other );

    /**
     * @brief Move assignment operator for heterogeneous arrays.
     *
     * @param[in] other the heterogeneous array to move
     * @return reference to this array for further assignments
     *
     * The other array will not contain any valid data afterwards.
     */
    HArray<ValueType>& operator=( HArray<ValueType>&& other ) noexcept;

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
    void swap( HArray<ValueType>& other );

    /** Override the method Printable::writeAt */

    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Implementation of pure method.
     */
    virtual common::ScalarType getValueType() const;

    // Make methods visible of _HArray (independent of the value type)

    using _HArray::size;
    using _HArray::reserve;
    using _HArray::purge;
    using _HArray::capacity;
    using _HArray::clear;
    using _HArray::resize;


protected:

    using _HArray::touch;

private:

    // Provide methods for Write and Read access where context is already found at a certain index

    using _HArray::reserveWithIndex;
    using _HArray::resizeWithIndex;
    using _HArray::capacityWithIndex;
    using _HArray::clearWithIndex;

    ValueType* get( ContextDataIndex index );

    const ValueType* get( ContextDataIndex index ) const;
};

/* ==================================================================== */
/*    Implementation of methods                                         */
/* ==================================================================== */

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
    touch( context );   // first context
    SCAI_LOG_DEBUG( logger, "created new HArray: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::HArray( MemoryPtr memory ) :
    _HArray( 0, sizeof( ValueType ) )
{
    touch( memory );   // first context
    SCAI_LOG_DEBUG( logger, "created new HArray: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::HArray( const IndexType n ) :

    _HArray( n, sizeof( ValueType ) )

{
    // reserve memory on the host

    ContextPtr hostPtr = Context::getHostPtr();

    _HArray::reserve( hostPtr, n );

    SCAI_LOG_DEBUG( logger, "created new HArray: " << *this )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::HArray( const IndexType n, ContextPtr context ) :

    _HArray( n, sizeof( ValueType ) )

{
    _HArray::reserve( context, n );

    SCAI_LOG_DEBUG( logger, "created new HArray: " << *this)
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::HArray( const IndexType n, const ValueType& value, ContextPtr context ) :

    _HArray( n, sizeof( ValueType ) )

{
    touch( context );  // first touch here

    // Note: in constructor of the HArray lock of accesses is not required but done for consistency

    ContextPtr host = Context::getHostPtr();

    // Use of acquireAccess guarantees allocation of data

    ContextDataIndex index = acquireWriteAccess( host, false );

    if ( n > 0 )
    {
        ValueType* hostData = get( index );

        #pragma omp parallel for

        for ( IndexType i = 0; i < n; ++i )
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
    touch( context );         // make sure that raw data is copied to this context
    _HArray::_setRawData( n, values );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::HArray( const HArray<ValueType>& other ):

    _HArray( other )

{
    SCAI_LOG_DEBUG( logger, "copy constructor, other = " << other )
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>::HArray( HArray<ValueType>&& other ) noexcept:

    _HArray( std::move( other ) )

{
    SCAI_LOG_DEBUG( logger, "move constructor of HArray, this = " << *this << ", other = " << other )
}

template <typename ValueType>
HArray<ValueType>::HArray( std::initializer_list<ValueType> init, ContextPtr context )
    : _HArray ( 0, sizeof( ValueType ) )
{
    touch( context );
    // std::initializer_list actually guarantees that begin() returns a pointer (to const)
    // to the data, so this should be safe!
    const ValueType* data = init.begin();
    _HArray::_setRawData( init.size(), data );
}

/* ---------------------------------------------------------------------------------*/

template <typename ValueType>
HArray<ValueType>::HArray( const std::vector<ValueType>& values, ContextPtr context ) :

    _HArray( 0, sizeof( ValueType ) )
{
    touch( context );
    _HArray::_setRawData( values.size(), values.data() );
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
HArray<ValueType>* HArray<ValueType>::newArray() const
{
    return new HArray<ValueType>();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
common::ScalarType HArray<ValueType>::getValueType() const
{
    // Note: this is implementation of the pure method of base class _HArray.
    return common::TypeTraits<ValueType>::stype;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::setRawData( const IndexType size, const ValueType src[] )
{
    _HArray::_setRawData( size, src );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::setSameValue( const IndexType size, const ValueType value )
{
    std::unique_ptr<ValueType[]> data( new ValueType[size] );

    for ( IndexType i = 0; i < size; ++i )
    {
        data[i] = value;
    }

    _HArray::_setRawData( size, data.get() );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>& HArray<ValueType>::operator=( const HArray<ValueType>& other )
{
    SCAI_LOG_DEBUG( logger, "copy assignment, other = " << other )

    _HArray::operator=( other );

    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArray<ValueType>& HArray<ValueType>::operator=( HArray<ValueType>&& other ) noexcept
{
    // Note: HArray has no own member variables

    SCAI_LOG_DEBUG( logger, "move assignment: other = " << other << " will be moved to this = " << *this )
    _HArray::operator=( std::move( other ) );
    SCAI_LOG_DEBUG( logger, "move assignment: other = " << other << " has been moved to this = " << *this )
    return *this;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::assign( const HArray<ValueType>& other, ContextPtr context )
{
    SCAI_LOG_DEBUG( logger, other << " will be assigned to " << *this )

    _HArray::assign( other, context );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::swap( HArray<ValueType>& other )
{
    _HArray::swap( other );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ValueType* HArray<ValueType>::get( ContextDataIndex index )
{
    return static_cast<ValueType*>( _HArray::get( index ) );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
const ValueType* HArray<ValueType>::get( ContextDataIndex index ) const
{
    return static_cast<const ValueType*>( _HArray::get( index ) );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void HArray<ValueType>::writeAt( std::ostream& stream ) const
{
    _HArray::writeAtTyped( stream, common::TypeTraits<ValueType>::id() );
}

/* ---------------------------------------------------------------------------------*/

} /* end namespace hmemo */

} /* end namespace scai */
