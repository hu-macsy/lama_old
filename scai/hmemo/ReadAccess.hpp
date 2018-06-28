/**
 * @file ReadAccess.hpp
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

#include <functional>

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

    const HArray<ValueType>* mArray;      // read access to this associated LAMA array

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
     * @brief Move constructor for ReadAccess.
     */
    ReadAccess( ReadAccess<ValueType>&& other ) noexcept;

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
    operator const ValueType* () const;

    /**
     * @brief return the memory where data has been allocated
     */

    const Memory& getMemory() const;

    /**
     * @brief Return a single value to host memory.
     *
     * @param[out] val will contain value from array[pos]
     * @param[in] pos is the position of array to read from, 0 <= pos < size()
     */
    void getValue( ValueType& val, const IndexType pos ) const;

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

    std::function<void()> releaseDelayed();

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

/**
 * @brief Return a ReadAccess for the provided array.
 *
 * This method is just a convenience function which allows acquiring a ReadAccess without
 * redundantly specifying the valuetype. See the example below for an idea of intended usage.
 *
 * @code{.cpp}
 * const auto array = HArray<double>({ 2.0, 3.0 });
 *
 * // Call ReadAccess constructor, must redundantly specify type
 * const ReadAccess<double> read1(array);
 *
 * // Call readAccess, type is deduced by compiler.
 * const auto read2 = readAccess(array);
 * @endcode
 */
template <typename ValueType>
inline ReadAccess<ValueType> readAccess( const HArray<ValueType>& array )
{
    return ReadAccess<ValueType>( array );
}

/**
 * @brief Return a ReadAccess on the given context for the provided array.
 *
 * Same as readAccess(const HArray<ValueType> &), but data is copied directly
 * into memory provided by the given context.
 */
template <typename ValueType>
inline ReadAccess<ValueType> readAccess( const HArray<ValueType>& array, ContextPtr context )
{
    return ReadAccess<ValueType>( array, context );
}

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
    ContextPtr contextPtr = Context::getContextPtr( common::ContextType::Host );
    SCAI_LOG_DEBUG( logger, "ReadAccess<" << common::TypeTraits<ValueType>::id()
                    << "> : create for " << array << " @ " << *contextPtr )
    mContextDataIndex = mArray->acquireReadAccess( contextPtr );
    mData = mArray->get( mContextDataIndex );
}

template <typename ValueType>
ReadAccess<ValueType>::ReadAccess( ReadAccess&& other ) noexcept
    : mArray( other.mArray ),
      mData( other.mData ),
      mContextDataIndex( other.mContextDataIndex )
{
    other.mArray = nullptr;
    other.mData = nullptr;
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
const Memory& ReadAccess<ValueType>::getMemory() const
{
    SCAI_ASSERT_ERROR( mArray, "ReadAccess has already been released." )
    return mArray->getMemory( mContextDataIndex );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ReadAccess<ValueType>::getValue( ValueType& val, const IndexType pos ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mArray->size(), "Index out of range" )

    const Memory& mem = mArray->getMemory( mContextDataIndex );
    const Memory& hostMem = *Context::getHostPtr()->getLocalMemoryPtr();
    mem.memcpyTo( hostMem, &val, mData + pos, sizeof( ValueType ) );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
std::function<void()> ReadAccess<ValueType>::releaseDelayed()
{
    SCAI_ASSERT( mArray, "releaseDelay not possible on released access" )
    void ( _HArray::*releaseAccess ) ( ContextDataIndex ) const = &_HArray::releaseReadAccess;
    const _HArray* ctxArray = mArray;
    // This access itself is treated as released
    mArray = NULL;
    return std::bind( releaseAccess, ctxArray, mContextDataIndex );
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
ReadAccess<ValueType>::operator const ValueType* () const
{
    SCAI_ASSERT( mArray, "ReadAccess has already been released." )
    return mData;
}

} /* end namespace hmemo */

} /* end namespace scai */
