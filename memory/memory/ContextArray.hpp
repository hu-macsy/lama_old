/**
 * @file ContextArray.hpp
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
 * @brief Common base class for dynamic array classes where the array data can be
 *        used in different contexts and where the data is moved implicitly
 *        when corresponding read/write accesses are required.
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

class COMMON_DLL_IMPORTEXPORT ContextArray: public Printable
{
    // Member variables of this class

protected:

    IndexType mSize;        //!< number of entries for the context array, common for all contexts
    IndexType mValueSize;   //!< number of bytes needed for one data element

    bool constFlag;         //!< if true the array cannot be written

    mutable ContextManager mContextManager;  //!< takes control of accesses and allocations

public:

    virtual ~ContextArray()
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

    static ContextArray* create( const Scalar::ScalarType type );

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

    explicit ContextArray( const IndexType n, const IndexType size )
                    : mSize( n ), mValueSize( size ), constFlag( false )
    {
        LAMA_LOG_DEBUG( logger, "construct LAMAArray, mSize = " << mSize << ", mValueSize " << mValueSize )
    }

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

    void releaseWriteAccess( ContextDataIndex );

    /** Static class variable for logger. */

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

/* ---------------------------------------------------------------------------------*/

inline IndexType ContextArray::size() const
{
    return mSize;
}

/* ---------------------------------------------------------------------------------*/

inline void ContextArray::prefetch( ContextPtr context ) const
{
    mContextManager.prefetch( context, mSize * mValueSize );
}

/* ---------------------------------------------------------------------------------*/

inline bool ContextArray::isValid( ContextPtr context ) const
{
    return mContextManager.isValid( context );
}

/* ---------------------------------------------------------------------------------*/

inline IndexType ContextArray::capacity( ContextPtr context ) const
{
    return mContextManager.capacity( context ) / mValueSize;
}

/* ---------------------------------------------------------------------------------*/

IndexType ContextArray::capacity( ContextDataIndex index ) const
{
    ContextData& entry = mContextManager[index];

    return static_cast<IndexType>( entry.capacity() / sizeof( ValueType ) );
}

/* ---------------------------------------------------------------------------------*/

inline ContextPtr ContextArray::getValidContext( const Context::ContextType preferredType ) const
{
    return mContextManager.getValidContext( preferredType );
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ConextArray::acquireReadAccess( ContextPtr context ) const
{
    LAMA_LOG_DEBUG( logger, "acquireReadAccess for " << *this );

    size_t allocSize = mSize * mValueSize;
    size_t validSize = allocSize;                   // read access needs valid data in any case

    return mContextManager.acquireAccess( context, ContextData::Read, allocSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

void ContextArray::releaseReadAccess( ContextDataIndex index ) const
{
    mContextManager.releaseAccess( index, ContextData::Read );
}

/* ---------------------------------------------------------------------------------*/

ContextDataIndex ContextArray::acquireWriteAccess( ContextPtr context, bool keepFlag )
{
    size_t allocSize = mSize * sizeof( ValueType );
    size_t validSize = keepFlag ? allocSize : 0 ;    // valid data only if keepFlag is set

    return mContextManager.acquireAccess( context, ContextData::Write, allocSize, validSize );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
void LAMAArray<ValueType>::releaseWriteAccess( ContextDataIndex index )
{
    mContextManager.releaseAccess( index, ContextData::Write );
}

/* ---------------------------------------------------------------------------------*/

LAMA_LOG_DEF_LOGGER( ContextArray::logger, "ContextArray" )

}  // namespace 
