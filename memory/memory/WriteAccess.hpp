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
 * @brief WriteAccess.hpp
 * @author Thomas Brandes
 * @date 29.04.2011
 */

#pragma once

// for dll_import
#include <common/config.hpp>

// base classes
#include <memory/Access.hpp>

// others
#include <memory/LAMAArray.hpp>
#include <memory/Context.hpp>

#include <common/Exception.hpp>

// logging
#include <logging/logging.hpp>

namespace memory
{

/**
 * @brief The Template WriteAccess is used to enforce the consistency of the template LAMAArray.
 *
 * WriteAccess enforces the consistency of the template LAMAArray by following the RAII Idiom. This is
 * done by acquiring a write lock on a LAMAArray in the constructor and releasing this write lock in
 * the destructor. There for a WriteAccess should be only used as a stack object.
 *
 * @tparam ValueType is the value type stored in the wrapped container.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT WriteAccess: public Access
{
    // member variables

protected:

    LAMAArray <ValueType>* mArray;

    ValueType* mData;   // current pointer to the data, cache value for mArray->get( mContextDataIndex )

    ContextDataIndex mContextDataIndex; // index for this context access

public:

    /**
     * @brief Acquires a WriteAccess to the passed LAMAArray for the passed location.
     *
     * @param[in] array     the LAMAArray to acquire a WriteAccess for
     * @param[in] context   the context to acquire a WriteAccess for
     * @param[in] keep      if false, implicit clear, old values of the array are no more needed
     * @throws Exception    if the WriteAccess can not be acquired, e.g. because another WriteAccess exists.
     *
     * This constructor acquires a WriteAccess to the passed LAMAArray for the passed location. Depending on the flag
     * keep the LAMAArray is updated on location or not. This flag can be used to acquire write only access, e.g. if
     * one want to write only to the passed wrapped LAMAArray. If the passed LAMAArray is not valid at location and keep
     * is true the LAMAArray will be made valid at location before this constructor returns.
     */
    WriteAccess( LAMAArray<ValueType>& array, ContextPtr context, const bool keep = true );

    /**
     * @brief Acquires a WriteAccess to the passed LAMAArray for a valid context.
     *
     * @param[in] array     the LAMAArray to acquire a WriteAccess for
     * @throws Exception    if the WriteAccess can not be acquired, e.g. because another WriteAccess exists.
     *
     * This constructor acquires a WriteAccess to the passed LAMAArray for a valid location. The chosen Location
     * will be the first valid location. This Access can be used to call, e.g. resize. If it does not matter where
     * the data currently resides.
     */
    explicit WriteAccess( LAMAArray<ValueType>& array );

    /**
     * @brief Releases the WriteAccess on the associated LAMAArray.
     */
    virtual ~WriteAccess();

    /**
     * @brief Returns a pointer to the data of the wrapped LAMAArray.
     *
     * @return a pointer to the wrapped LAMAArray.
     */
    ValueType* get();

    /**
     * @brief Clear of the LAMA array.
     *
     * This operation is the same as calling it directly for the LAMA array.
     * It does not free any data but all data is no more valid.
     */
    void clear();

    /**
     * @brief Resizes the wrapped LAMAArray.
     *
     * @param[in] newSize   the new size of the wrapped LAMAArray
     *
     * If a reallocation is necessary it is only done at the associated location
     * for all other locations this is done lazy.
     *
     */
    void resize( const IndexType newSize );

    /**
     * @brief Reserves storage for the wrapped LAMAArray at the associated location.
     *
     * @param[in] capacity  the number of elements that should fit into the new storage
     */
    void reserve( const IndexType capacity );

    /**
     *  @brief Queries the capacity of the array on the reserved context.
     */
    IndexType capacity() const;

    /**
     * @brief Releases the WriteAccess on the associated LAMAArray.
     */
    virtual void release();

    /**
     * @brief Override method of base class Printable
     */
    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Returns the size of the wrapped LAMAArray.
     *
     * @return  the size of the wrapped LAMAArray
     */
    inline IndexType size() const;

protected:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

/**
 * @brief WriteOnlyAccess is a write access where no existing values of the array are needed (keepFlag).
 *
 * This derived class has been added for more convenience as it avoids the use of the keepFlag param.
 *
 * A WriteOnlyAccess should be used whenever possible. It avoids any memory transfer of no more
 * needed values between devices and in case of a reallocation it avoids copying of old values.
 *
 * @tparam ValueType is the value type stored in the wrapped container.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT WriteOnlyAccess: public WriteAccess<ValueType>
{
public:

    /** Create a write access with keep flag = false.
     *
     *  Attention: this kind of write access assumes that the array is completely new written.
     */

    WriteOnlyAccess( LAMAArray<ValueType>& array, ContextPtr context )
        : WriteAccess<ValueType>( array, context, false )
    {
        LAMA_LOG_DEBUG( logger, "WriteOnlyAccess<" << Scalar::getType<ValueType>() << ">" )
    }

    /**
     *  @brief Acquire a write only access with a resize.
     *
     * @param[in] array     the LAMAArray for which access is required
     * @param[in] context   the context where data will be written later
     * @param[in] size      the new size of the LAMA array.
     *
     * Attention: this kind of write access assumes that the array is completely new written.
     */
    WriteOnlyAccess( LAMAArray<ValueType>& array, ContextPtr context, const IndexType size )
        : WriteAccess<ValueType>( array, context, size, false )
    {
        this->clear( mContextDataIndex );           // invalides all data before resize
        this->resize( mContextDataIndex, size );
    }

    ~WriteOnlyAccess()
    {
        LAMA_LOG_DEBUG( WriteAccess<ValueType>::logger, "~WriteOnlyAccess<" << Scalar::getType<ValueType>() << ">" )
    }

protected:
    using WriteAccess<ValueType>::logger;   // no additinal logger for this derived class

private:

    using WriteAccess<ValueType>::mArray;
    using WriteAccess<ValueType>::mContextDataIndex;
};

/* --------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, WriteAccess<ValueType>::logger, "WriteAccess" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
WriteAccess<ValueType>::WriteAccess( LAMAArray<ValueType>& view, ContextPtr context, const bool keep /* = true*/ )
    : mArray( &view )
{
    if ( view.constFlag )
    {
        COMMON_THROWEXCEPTION( "write on const array not allowed" )
    }

    LAMA_LOG_DEBUG( logger, "acquire write access for " << *mArray << " at " << *context << ", keep = " << keep )
    mContextDataIndex = mArray->acquireWriteAccess( context, keep );
    mData = mArray->get( mContextDataIndex );     // cache the data pointer
}

template<typename ValueType>
WriteAccess<ValueType>::WriteAccess( LAMAArray<ValueType>& array )
    : mArray( &array )
{
    if ( array.constFlag )
    {
        COMMON_THROWEXCEPTION( "write on const array not allowed" )
    }

    LAMA_LOG_DEBUG( logger, "acquire write access for " << *mArray << " at first valid context " )
    const bool keepFlag = true;
    mContextDataIndex = mArray->acquireWriteAccess( Context::getContext( context::Host ), keepFlag );
    mData = mArray->get( mContextDataIndex );     // cache the data pointer
}

template<typename ValueType>
WriteAccess<ValueType>::~WriteAccess()
{
    LAMA_LOG_DEBUG( logger, "~WriteAccess: release" )
    release();
}

template<typename ValueType>
ValueType* WriteAccess<ValueType>::get()
{
    COMMON_ASSERT( mArray, "illegal get(): access has already been released." )

    return mData;    // mData might be NULL if size of array is 0
}

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

template<typename ValueType>
void WriteAccess<ValueType>::clear()
{
    COMMON_ASSERT( mArray, "WriteAccess has already been released." )
    mArray->clear( mContextDataIndex );
    LAMA_LOG_DEBUG( logger, "cleared " << *mArray )
    mData = mArray->get( mContextDataIndex );     // not really needed
}

template<typename ValueType>
void WriteAccess<ValueType>::resize( const IndexType newSize )
{
    COMMON_ASSERT( mArray, "WriteAccess has already been released." )
    // do not log before check of mArray
    LAMA_LOG_DEBUG( logger, "resize " << *mArray << " to new size " << newSize )
    mArray->resize( mContextDataIndex, newSize );
    mData = mArray->get( mContextDataIndex );     // data might be reallocated
}

template<typename ValueType>
void WriteAccess<ValueType>::reserve( const IndexType capacity )
{
    COMMON_ASSERT( mArray, "WriteAccess has already been released." )
    LAMA_LOG_DEBUG( logger, "reserve " << *mArray << " to new capacity " << capacity )
    mArray->reserve( mContextDataIndex, capacity ); // copy = true for old data
    mData = mArray->get( mContextDataIndex );     // data might be reallocated
}

template<typename ValueType>
IndexType WriteAccess<ValueType>::capacity() const
{
    COMMON_ASSERT( mArray, "WriteAccess has already been released." )
    return mArray->capacity( mContextDataIndex );
}

template<typename ValueType>
void WriteAccess<ValueType>::release()
{
    if ( mArray )
    {
        LAMA_LOG_DEBUG( logger, "release write access for " << *mArray )
        mArray->releaseWriteAccess( mContextDataIndex );
    }
    else
    {
        LAMA_LOG_DEBUG( logger, "release write access for an already released LAMAArray" )
    }

    mArray = 0;
    mData = 0;
}

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

} // namespace

