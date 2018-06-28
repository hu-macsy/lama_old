/**
 * @file WriteAccess.hpp
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
 * @brief Definition of a template class for writing to a LAMA array.
 * @author Thomas Brandes
 * @date 29.04.2011
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/hmemo/Access.hpp>

// local library
#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/Context.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/macros/throw.hpp>

#include <functional>

namespace scai
{

namespace hmemo
{

/**
 * @brief The Template WriteAccess is used to enforce the consistency of the template HArray.
 *
 * WriteAccess enforces the consistency of the template HArray by following the RAII Idiom. This is
 * done by acquiring a write lock on a HArray in the constructor and releasing this write lock in
 * the destructor. There for a WriteAccess should be only used as a stack object.
 *
 * @tparam ValueType is the value type stored in the wrapped container.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT WriteAccess: public Access
{
    // member variables

protected:

    HArray <ValueType>* mArray;

    ValueType* mData;   // current pointer to the data, cache value for mArray->get( mContextDataIndex )

    ContextDataIndex mContextDataIndex; // index for this context access

public:

    /**
     * @brief Acquires a WriteAccess to the passed HArray for the passed location.
     *
     * @param[in] array       the HArray to acquire a WriteAccess for
     * @param[in] contextPtr  the context to acquire a WriteAccess for
     * @param[in] keep        if false, implicit clear, old values of the array are no more needed
     * @throws Exception      if the WriteAccess can not be acquired, e.g. because another WriteAccess exists.
     *
     * This constructor acquires a WriteAccess to the passed HArray for the passed location. Depending on the flag
     * keep the HArray is updated on location or not. This flag can be used to acquire write only access, e.g. if
     * one want to write only to the passed wrapped HArray. If the passed HArray is not valid at location and keep
     * is true the HArray will be made valid at location before this constructor returns.
     */
    WriteAccess( HArray<ValueType>& array, ContextPtr contextPtr, const bool keep = true );

    /**
     * @brief Acquires a WriteAccess to the passed HArray for the passed location.
     *
     * @param[in] array     the HArray to acquire a WriteAccess for
     * @param[in] keep      if false, implicit clear, old values of the array are no more needed
     * @throws Exception    if the WriteAccess can not be acquired, e.g. because another WriteAccess exists.
     *
     * This constructor acquires a WriteAccess to the passed HArray for the passed location. Depending on the flag
     * keep the HArray is updated on location or not. This flag can be used to acquire write only access, e.g. if
     * one want to write only to the passed wrapped HArray. If the passed HArray is not valid at location and keep
     * is true the HArray will be made valid at location before this constructor returns.
     */
    WriteAccess( HArray<ValueType>& array, const bool keep = true );

    /**
     * @brief Move constructor for WriteAccess.
     */
    WriteAccess( WriteAccess<ValueType>&& other ) noexcept;

    /**
     * @brief Releases the WriteAccess on the associated HArray.
     */
    virtual ~WriteAccess();

    /**
     * @brief Returns a pointer to the data of the wrapped HArray.
     *
     * @return a pointer to the data of the wrapped HArray.
     */
    ValueType* get();

    /**
     * @brief Returns a pointer to the data of the wrapped HArray.
     *
     * @return a pointer to the data of the wrapped HArray.
     */
    const ValueType* get() const;

    /**
     * @brief Support implicit type conversion to pointer of the data.
     *
     * @return a pointer to the wrapped HArray.
     */
    operator ValueType* ();

    /**
     * @brief Clear of the LAMA array.
     *
     * This operation is the same as calling it directly for the LAMA array.
     * It does not free any data but all data is no more valid.
     */
    void clear();

    /**
     * @brief Resizes the wrapped HArray.
     *
     * @param[in] newSize   the new size of the wrapped HArray
     *
     * If a reallocation is necessary it is only done at the associated location
     * for all other locations this is done lazy.
     *
     */
    void resize( const IndexType newSize );

    /**
     * @brief Reserves storage for the wrapped HArray at the associated location.
     *
     * @param[in] capacity  the number of elements that should fit into the new storage
     */
    void reserve( const IndexType capacity );

    /**
     *  @brief Queries the capacity of the array on the reserved context.
     */
    IndexType capacity() const;

    /**
     * @brief return the memory where data has been allocated
     */

    const Memory& getMemory() const;

    /**
     * @brief This method writes just one value at a certain position
     *
     * @param[in] val is the value to write
     * @param[in] pos is the position where to write 0 <= pos < capacity()
     */

    void setValue( const ValueType val, const IndexType pos );

    /**
     * @brief Return a single value to host memory.
     *
     * @param[out] val will contain value from array[pos]
     * @param[in] pos is the position of array to read from, 0 <= pos < size()
     */
    void getValue( ValueType& val, const IndexType pos ) const;

    /**
     * @brief Releases the WriteAccess on the associated HArray.
     */
    virtual void release();

    /**
     * @brief Return function that does the release later.
     *
     * This method can be used for asynchronous operations so the array is kept
     * locked even if the ~WriteAccess has been called.
     */

    std::function<void()> releaseDelayed();

    /**
     * @brief Override method of base class Printable
     */
    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Returns the size of the wrapped HArray.
     *
     * @return  the size of the wrapped HArray
     */
    inline IndexType size() const;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/**
 * @brief Obtain a WriteAccess to the given array.
 *
 * Analogous to readAccess(const HArray & array). See its documentation for
 * motivation and intended usage of this function.
 */
template <typename ValueType>
WriteAccess<ValueType> writeAccess( HArray<ValueType>& array )
{
    return WriteAccess<ValueType>( array );
}

/**
 * @brief Obtain a WriteAccess for the supplied context to the given array.
 *
 * Analogous to readAccess(const HArray & array, ContextPtr). See its documentation for
 * motivation and intended usage of this function.
 */
template <typename ValueType>
WriteAccess<ValueType> writeAccess( HArray<ValueType>& array, ContextPtr context )
{
    return WriteAccess<ValueType>( array, context );
}

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, WriteAccess<ValueType>::logger, "WriteAccess" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
WriteAccess<ValueType>::WriteAccess( HArray<ValueType>& array, ContextPtr contextPtr, const bool keep )
    : mArray( &array )
{
    SCAI_ASSERT( !array.isConst(), "WriteAccess on const array not allowed: " << array )
    SCAI_LOG_DEBUG( logger, "acquire write access for " << *mArray << " at " << *contextPtr << ", keep = " << keep )
    mContextDataIndex = mArray->acquireWriteAccess( contextPtr, keep );
    mData = mArray->get( mContextDataIndex );     // cache the data pointer
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
WriteAccess<ValueType>::WriteAccess( HArray<ValueType>& array, const bool keep /* = true*/ )
    : mArray( &array )
{
    SCAI_ASSERT( !array.isConst(), "WriteAccess on const array not allowed: " << array )
    ContextPtr contextPtr = Context::getContextPtr( common::ContextType::Host );
    SCAI_LOG_DEBUG( logger, "acquire write access for " << *mArray << " at " << *contextPtr << ", keep = " << keep )
    mContextDataIndex = mArray->acquireWriteAccess( contextPtr, keep );
    mData = mArray->get( mContextDataIndex );     // cache the data pointer
}

template <typename ValueType>
WriteAccess<ValueType>::WriteAccess( WriteAccess<ValueType>&& other ) noexcept
    :   mArray( other.mArray ),
        mData( other.mData ),
        mContextDataIndex( other.mContextDataIndex )
{
    other.mArray = nullptr;
    other.mData = nullptr;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
WriteAccess<ValueType>::~WriteAccess()
{
    SCAI_LOG_DEBUG( logger, "~WriteAccess: release" )
    release();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
const Memory& WriteAccess<ValueType>::getMemory() const
{
    SCAI_ASSERT( mArray, "WriteAccess has already been released." )
    return mArray->getMemory( mContextDataIndex );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void WriteAccess<ValueType>::setValue( const ValueType val, const IndexType pos )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mArray->size(), "Index out of range" )

    const Memory& mem = mArray->getMemory( mContextDataIndex );
    const Memory& hostMem = *Context::getHostPtr()->getLocalMemoryPtr();
    mem.memcpyFrom( mData + pos, hostMem, &val, sizeof( ValueType ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void WriteAccess<ValueType>::getValue( ValueType& val, const IndexType pos ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mArray->size(), "Index out of range" )

    const Memory& mem = mArray->getMemory( mContextDataIndex );
    const Memory& hostMem = *Context::getHostPtr()->getLocalMemoryPtr();
    mem.memcpyTo( hostMem, &val, mData + pos, sizeof( ValueType ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType* WriteAccess<ValueType>::get()
{
    SCAI_ASSERT( mArray, "illegal get(): access has already been released." )
    return mData;    // mData might be NULL if size of array is 0
}

template<typename ValueType>
const ValueType* WriteAccess<ValueType>::get() const
{
    SCAI_ASSERT( mArray, "illegal get(): access has already been released." )
    return mData;    // mData might be NULL if size of array is 0
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
WriteAccess<ValueType>::operator ValueType* ()
{
    SCAI_ASSERT( mArray, "illegal get(): access has already been released." )
    return mData;    // mData might be NULL if size of array is 0
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
inline IndexType WriteAccess<ValueType>::size() const
{
    if ( mArray )
    {
        return mArray->size();
    }
    else
    {
        COMMON_THROWEXCEPTION( "cannot call size on released array" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void WriteAccess<ValueType>::clear()
{
    SCAI_ASSERT( mArray, "WriteAccess has already been released." )
    mArray->clearWithIndex( mContextDataIndex );
    SCAI_LOG_DEBUG( logger, "cleared " << *mArray )
    mData = mArray->get( mContextDataIndex );     // not really needed
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void WriteAccess<ValueType>::resize( const IndexType newSize )
{
    SCAI_ASSERT( mArray, "WriteAccess has already been released." )
    // do not log before check of mArray
    SCAI_LOG_DEBUG( logger, "resize " << *mArray << " to new size " << newSize )
    mArray->resizeWithIndex( mContextDataIndex, newSize );
    mData = mArray->get( mContextDataIndex );     // data might be reallocated
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void WriteAccess<ValueType>::reserve( const IndexType capacity )
{
    SCAI_ASSERT( mArray, "WriteAccess has already been released." )
    SCAI_LOG_DEBUG( logger, "reserve " << *mArray << " to new capacity " << capacity )
    mArray->reserveWithIndex( mContextDataIndex, capacity ); // copy = true for old data
    mData = mArray->get( mContextDataIndex );     // data might be reallocated
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType WriteAccess<ValueType>::capacity() const
{
    SCAI_ASSERT( mArray, "WriteAccess has already been released." )
    return mArray->capacityWithIndex( mContextDataIndex );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void WriteAccess<ValueType>::release()
{
    if ( mArray )
    {
        SCAI_LOG_DEBUG( logger, "release write access for " << *mArray )
        mArray->releaseWriteAccess( mContextDataIndex );
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "release write access for an already released HArray" )
    }

    mArray = 0;
    mData = 0;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
std::function<void()> WriteAccess<ValueType>::releaseDelayed()
{
    SCAI_ASSERT( mArray, "releaseDelay not possible on released access" )
    void ( _HArray::*releaseAccess ) ( ContextDataIndex ) = &_HArray::releaseWriteAccess;
    _HArray* ctxArray = mArray;
    // This access itself is treated as released
    mArray = 0;
    mData  = 0;
    return std::bind( releaseAccess, ctxArray, mContextDataIndex );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void WriteAccess<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "WriteAccess to ";

    if ( mArray )
    {
        stream << *mArray;
    }
    else
    {
        stream << "already releases array view.";
    }
}

} /* end namespace hmemo */

} /* end namespace scai */
