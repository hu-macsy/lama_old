/**
 * @file WriteAccess.hpp
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
#include <scai/common/function.hpp>
#include <scai/common/bind.hpp>

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
     * @brief Releases the WriteAccess on the associated HArray.
     */
    virtual ~WriteAccess();

    /**
     * @brief Returns a pointer to the data of the wrapped HArray.
     *
     * @return a pointer to the wrapped HArray.
     */
    ValueType* get();

    /**
     * @brief Support implicit type conversion to pointer of the data.
     *
     * @return a pointer to the wrapped HArray.
     */
    operator ValueType*();

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
     * @brief Releases the WriteAccess on the associated HArray.
     */
    virtual void release();

    /**
     * @brief Return function that does the release later.
     *
     * This method can be used for asynchronous operations so the array is kept 
     * locked even if the ~WriteAccess has been called.
     */

    common::function<void()> releaseDelayed();

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

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, WriteAccess<ValueType>::logger, "WriteAccess" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
WriteAccess<ValueType>::WriteAccess( HArray<ValueType>& array, ContextPtr contextPtr, const bool keep )
    : mArray( &array )
{
    SCAI_ASSERT( !array.constFlag, "WriteAccess on const array not allowed: " << array )

    SCAI_LOG_DEBUG( logger, "acquire write access for " << *mArray << " at " << *contextPtr << ", keep = " << keep )
    mContextDataIndex = mArray->acquireWriteAccess( contextPtr, keep );
    mData = mArray->get( mContextDataIndex );     // cache the data pointer
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
WriteAccess<ValueType>::WriteAccess( HArray<ValueType>& array, const bool keep /* = true*/ )
    : mArray( &array )
{
    SCAI_ASSERT( !array.constFlag, "WriteAccess on const array not allowed: " << array )

    ContextPtr contextPtr = Context::getContextPtr( common::context::Host );

    SCAI_LOG_DEBUG( logger, "acquire write access for " << *mArray << " at " << *contextPtr << ", keep = " << keep )
    mContextDataIndex = mArray->acquireWriteAccess( contextPtr, keep );
    mData = mArray->get( mContextDataIndex );     // cache the data pointer
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
WriteAccess<ValueType>::~WriteAccess()
{
    SCAI_LOG_DEBUG( logger, "~WriteAccess: release" )
    release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType* WriteAccess<ValueType>::get()
{
    SCAI_ASSERT( mArray, "illegal get(): access has already been released." )

    return mData;    // mData might be NULL if size of array is 0
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
WriteAccess<ValueType>::operator ValueType*()
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
    mArray->clear( mContextDataIndex );
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
    mArray->resize( mContextDataIndex, newSize );
    mData = mArray->get( mContextDataIndex );     // data might be reallocated
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void WriteAccess<ValueType>::reserve( const IndexType capacity )
{
    SCAI_ASSERT( mArray, "WriteAccess has already been released." )
    SCAI_LOG_DEBUG( logger, "reserve " << *mArray << " to new capacity " << capacity )
    mArray->reserve( mContextDataIndex, capacity ); // copy = true for old data
    mData = mArray->get( mContextDataIndex );     // data might be reallocated
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType WriteAccess<ValueType>::capacity() const
{
    SCAI_ASSERT( mArray, "WriteAccess has already been released." )
    return mArray->capacity( mContextDataIndex );
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
common::function<void()> WriteAccess<ValueType>::releaseDelayed()
{
    SCAI_ASSERT( mArray, "releaseDelay not possible on released access" )

    void ( _HArray::*releaseAccess ) ( ContextDataIndex ) = &_HArray::releaseWriteAccess;

    _HArray* ctxArray = mArray; 

    // This access itself is treated as released

    mArray = 0;
    mData  = 0;

    return common::bind( releaseAccess, ctxArray, mContextDataIndex ); 
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
