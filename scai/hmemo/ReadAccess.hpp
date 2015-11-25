/**
 * @file ReadAccess.hpp
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
 * @brief Definition of a template class ReadAccess for reading a HArray.
 * @author Thomas Brandes, Jiri Kraus
 * @date 29.04.2011
 */

#pragma once

// local library
#include <scai/hmemo/Access.hpp>
#include <scai/hmemo/HArray.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/config.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/function.hpp>
#include <scai/common/bind.hpp>

namespace scai
{

namespace hmemo
{

/**
 * @brief The template ReadAccess is used to enforce the consistency of the template HArray.
 *
 * ReadAccess enforces the consistency of the template HArray by following the RAII Idiom. This is
 * done by acquiring a read lock on a HArray in the constructor and releasing this read lock in
 * the destructor. Therefore a ReadAccess should be only used as a stack object.
 *
 * @tparam ValueType is the value type stored in the wrapped container.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT ReadAccess: public Access
{
    // Member variables

private:

    const HArray<ValueType>* mArray;   // read access to this associated LAMA array

    const ValueType* mData;               // pointer to the data used by the access

    ContextDataIndex mContextDataIndex;   // reference to the incarnation of the array

public:

    /**
     * @brief Acquires a ReadAccess to the passed HArray for a given context.
     *
     * @param[in] array      the HArray to acquire a ReadAccess for
     * @param[in] contextPtr the context that needs a read acess
     * @throws Exception     if the ReadAccess can not be acquired, e.g. because a WriteAccess exists.
     */

    ReadAccess( const HArray<ValueType>& array, ContextPtr contextPtr );

    /**
     * @brief Acquires a ReadAccess to the passed HArray for the host context.
     *
     * @param[in] array     the HArray to acquire a ReadAccess for
     * @throws Exception    if the ReadAccess can not be acquired, e.g. because a WriteAccess exists.
     */
    ReadAccess( const HArray<ValueType>& array );

    /**
     * @brief Releases the ReadAccess on the associated HArray.
     */
    virtual ~ReadAccess();

    /**
     * @brief Returns a valid pointer to the data usable for the context.
     *
     * @return a pointer to the wrapped HArray.
     */
    const ValueType* get() const;

    /**
     * @brief Allow type conversion.
     */
    operator const ValueType*() const;

    /**
     * @brief Releases the acquired ReadAccess.
     *
     * Release is mandatory to unlock the array so that it might be
     * used for further write accesses.
     */
    virtual void release();

    /**
     *  @brief Delay the release of the access in an own function
     *
     *  The access can no more be used afterwards but there is no
     *  release done with the destructor.
     */

    common::function<void()> releaseDelayed();

    /** 
     * @brief Output of this object in a stream. 
     */
    virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Returns the size of the array
     *
     * @return  the size of the wrapped HArray
     */
    IndexType size() const;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* ---------------------------------------------------------------------------------*/

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, ReadAccess<ValueType>::logger, "ReadAccess" )

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ReadAccess<ValueType>::ReadAccess( const HArray<ValueType>& array, ContextPtr contextPtr ) : mArray( &array )
{
    SCAI_ASSERT( contextPtr.get(), "NULL context for read access not allowed" )

    SCAI_LOG_DEBUG( logger, "ReadAccess<" << common::TypeTraits<ValueType>::id()
                    << "> : create for " << array << " @ " << *contextPtr )

    mContextDataIndex = mArray->acquireReadAccess( contextPtr );

    mData = mArray->get( mContextDataIndex );
}

template<typename ValueType>
ReadAccess<ValueType>::ReadAccess( const HArray<ValueType>& array ) : mArray( &array )
{
    ContextPtr contextPtr = Context::getContextPtr( common::context::Host );

    SCAI_LOG_DEBUG( logger, "ReadAccess<" << common::TypeTraits<ValueType>::id()
                    << "> : create for " << array << " @ " << *contextPtr )

    mContextDataIndex = mArray->acquireReadAccess( contextPtr );

    mData = mArray->get( mContextDataIndex );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ReadAccess<ValueType>::~ReadAccess()
{
    SCAI_LOG_DEBUG( logger, "~ReadAccess<" << common::TypeTraits<ValueType>::id() << ">" )
    release();
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void ReadAccess<ValueType>::release()
{
    if ( mArray )
    {
        SCAI_LOG_DEBUG( logger, "ReadAccess<" << common::TypeTraits<ValueType>::id() << ">: realase for " << *mArray )
        mArray->releaseReadAccess( mContextDataIndex );
    }

    mArray = NULL;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
common::function<void()> ReadAccess<ValueType>::releaseDelayed()
{
    SCAI_ASSERT( mArray, "releaseDelay not possible on released access" )

    void ( ContextArray::*releaseAccess ) ( ContextDataIndex ) const = &ContextArray::releaseReadAccess;

    const ContextArray* ctxArray = mArray;

    // This access itself is treated as released

    mArray = NULL;

    return common::bind( releaseAccess, ctxArray, mContextDataIndex );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void ReadAccess<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "ReadAccess<" << common::TypeTraits<ValueType>::id() << "> ";

    if ( mArray )
    {
        stream << "to " << *mArray;
    }
    else
    {
        stream << "released.";
    }
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
IndexType ReadAccess<ValueType>::size() const
{
    if ( mArray )
    {
        return mArray->size();
    }
    else
    {
        return 0;
    }
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
const ValueType* ReadAccess<ValueType>::get() const
{
    SCAI_ASSERT( mArray, "ReadAccess::get fails, has already been released." )

    return mData;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
ReadAccess<ValueType>::operator const ValueType*() const
{
    SCAI_ASSERT( mArray, "ReadAccess has already been released." )

    return mData;
}

} /* end namespace hmemo */

} /* end namespace scai */
