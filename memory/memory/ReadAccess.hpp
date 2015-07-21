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
 * @brief Definition of a template class ReadAccess for reading a LAMAArray.
 * @author Thomas Brandes, Jiri Kraus
 * @date 29.04.2011
 */

#pragma once

// for dll_import
#include <common/config.hpp>

#include <memory/Access.hpp>
#include <memory/LAMAArray.hpp>

#include <logging/logging.hpp>

namespace memory
{

/**
 * @brief The template ReadAccess is used to enforce the consistency of the template LAMAArray.
 *
 * ReadAccess enforces the consistency of the template LAMAArray by following the RAII Idiom. This is
 * done by acquiring a read lock on a LAMAArray in the constructor and releasing this read lock in
 * the destructor. Therefore a ReadAccess should be only used as a stack object.
 *
 * @tparam ValueType is the value type stored in the wrapped container.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT ReadAccess: public Access
{
    // Member variables

private:

    const LAMAArray<ValueType>* mArray;   // read access to this associated LAMA array

    ContextDataIndex mContextDataIndex;   // reference to the context data of the read

public:

    /**
     * @brief Acquires a ReadAccess to the passed LAMAArray for the passed context.
     *
     * @param[in] array     the LAMAArray to acquire a ReadAccess for
     * @param[in] context   the context that needs a read acess
     * @throws Exception    if the ReadAccess can not be acquired, e.g. because a WriteAccess exists.
     */
    ReadAccess( const LAMAArray<ValueType>& array, ContextPtr context ) : mArray( &array )
    {
        COMMON_ASSERT( context, "NULL context for read access not allowed" );

        LAMA_LOG_DEBUG( logger, "ReadAccess<" << getScalarType<ValueType>()
                        << "> : create for " << array << " @ " << *context )

        mContextDataIndex = mArray->acquireReadAccess( context );
    }

    /**
     * @brief Acquires a ReadAccess to the passed LAMAArray for the host context.
     *
     * @param[in] context   the context that needs a read acess
     * @throws Exception    if the ReadAccess can not be acquired, e.g. because a WriteAccess exists.
     */
    ReadAccess( const LAMAArray<ValueType>& array ) : mArray( &array )
    {
        ContextPtr contextPtr = Context::getContextPtr( context::Host );

        LAMA_LOG_DEBUG( logger, "ReadAccess<" << getScalarType<ValueType>()
                        << "> : create for " << array << " @ " << *contextPtr )

        mContextDataIndex = mArray->acquireReadAccess( contextPtr );
    }

    /**
     * @brief Releases the ReadAccess on the associated LAMAArray.
     */
    virtual ~ReadAccess()
    {
        LAMA_LOG_DEBUG( logger, "~ReadAccess<" << getScalarType<ValueType>() << ">" )
        release();
    }

    /**
     * @brief Returns a valid pointer to the data usable for the context.
     *
     * @return a pointer to the wrapped LAMAArray.
     */
    const ValueType* get() const
    {
        COMMON_ASSERT( mArray, "ReadAccess::get fails, has already been released." )

        return mArray->get( mContextDataIndex );
    }

    /**
     * @brief Releases the acquired ReadAccess.
     *
     * Release is mandatory to unlock the array so that it might be
     * used for further write accesses.
     */
    virtual void release()
    {
        if ( mArray )
        {
            LAMA_LOG_DEBUG( logger, "ReadAccess<" << getScalarType<ValueType>() << ">: realase for " << *mArray )
            mArray->releaseReadAccess( mContextDataIndex );
        }

        mArray = NULL;
    }

    virtual void writeAt( std::ostream& stream ) const
    {
        stream << "ReadAccess<" << getScalarType<ValueType>() << "> ";

        if ( mArray )
        {
            stream << "to " << *mArray;
        }
        else
        {
            stream << "released.";
        }
    }

    /**
     * @brief Returns the size of the array
     *
     * @return  the size of the wrapped LAMAArray
     */
    inline IndexType size() const
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

protected:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, ReadAccess<ValueType>::logger, "ReadAccess" )

} // namespace
