/**
 * @file HostReadAccess.hpp
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
 * @brief HostReadAccess.hpp
 * @author Thomas Brandes
 * @date 02.05.2011
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <common/config.hpp>

// base classes
#include <memory/ReadAccess.hpp>

namespace memory
{

/**
 * @brief HostReadAccess is a specialization of ReadAccess for a host context with an extended interface.
 *
 * @tparam ValueType is the type stored in the wrapped container.
 */
template<typename ValueType>
class HostReadAccess: public ReadAccess<ValueType>
{
public:
    /**
     * @brief Acquire a ReadAccess to the passed LAMAArray for the host location.
     *
     * @param[in] array     the LAMAArray to acquire a ReadAccess for
     * @throws Exception    if the ReadAccess can not be acquired, e.g. because a WriteContext exists.
     */
    HostReadAccess( const LAMAArray<ValueType>& array );

    /**
     * @brief Releases the ReadAccess on the associated LAMAArray.
     */
    virtual ~HostReadAccess();

    /**
     * @brief Constant access to the element i of the wrapped LAMAArray.
     *
     * @param[in] i   the index of the element of the wrapped LAMAArray to access
     * @return        a constant reference to the element i of the wrapped LAMAArray
     */
    inline const ValueType& operator[]( const IndexType i ) const;

    /**
     * @brief Conversion to a constant ValueType pointer (shortcut for get() ).
     *
     * @return  a constant pointer to the data of the wrapped LAMAArray
     */
    inline operator const ValueType* () const;

    using ReadAccess<ValueType>::get;

    using ReadAccess<ValueType>::size;

private:
    const ValueType* const mData;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, HostReadAccess<ValueType>::logger,
                              "ReadAccess.HostReadAccess" )

template<typename ValueType>
HostReadAccess<ValueType>::HostReadAccess( const LAMAArray<ValueType>& array ) :

    ReadAccess<ValueType>( array, Context::getContext( context::Host ) ), mData( get() )

{
    LAMA_LOG_DEBUG( logger, "read access on host, mData = " << mData );
}

template<typename ValueType>
HostReadAccess<ValueType>::~HostReadAccess()
{
    LAMA_LOG_DEBUG( logger, "~HostReadAccess" )
}

template<typename ValueType>
inline const ValueType& HostReadAccess<ValueType>::operator[]( const IndexType i ) const
{
    return mData[i];
}

template<typename ValueType>
inline HostReadAccess<ValueType>::operator const ValueType* () const
{
    LAMA_LOG_TRACE( logger, "mData = " << mData )
    return mData;
}

}
