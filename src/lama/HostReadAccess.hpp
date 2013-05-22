/**
 * @file HostReadAccess.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * $Id$
 */
#ifndef LAMA_HOSTREADACCESS_HPP_
#define LAMA_HOSTREADACCESS_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/ReadAccess.hpp>

// others
#include <lama/ContextFactory.hpp>

namespace lama
{

/**
 * @brief HostReadAccess is a specialization of ReadAccess for a host context with an extended interface.
 *
 * @tparam T is the type stored in the wrapped container.
 */
template<typename T>
class HostReadAccess: public ReadAccess<T>
{
public:
    typedef T ValueType; //!< This is the type stored in the wrapped container.

    /**
     * @brief Acquire a ReadAccess to the passed LAMAArray for the host location.
     *
     * @param[in] array     the LAMAArray to acquire a ReadAccess for
     * @throws Exception    if the ReadAccess can not be acquired, e.g. because a WriteContext exists.
     */
    HostReadAccess( const LAMAArray<ValueType>& array );

    /**
     * @brief Acquires a ReadAccess to the passed LAMAArrayConstView for the host location.
     *
     * @param[in] view      the LAMAArrayConstView to acquire a ReadAccess for
     * @throws Exception    if the ReadAccess can not be acquired, e.g. because a WriteContext exists.
     */
    HostReadAccess( const LAMAArrayConstView<ValueType>& view );

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
    inline operator const T*() const;

    using ReadAccess<ValueType>::get;

    using ReadAccess<ValueType>::size;

private:
    const ValueType* const mData;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename T>, HostReadAccess<T>::logger, "ReadAccess.HostReadAccess" )

template<typename T>
HostReadAccess<T>::HostReadAccess( const LAMAArray<ValueType>& array )
    : ReadAccess<T>( array, ContextFactory::getContext( Context::Host ) ),
      mData( get() )
{
    LAMA_LOG_DEBUG(logger, "read access on host, mData = " << mData);
}

template<typename T>
HostReadAccess<T>::HostReadAccess( const LAMAArrayConstView<ValueType>& view )
    : ReadAccess<T>( view, ContextFactory::getContext( Context::Host ) ), mData( get() )
{
    LAMA_LOG_DEBUG( logger, "read access on host, mData = " << mData )
}

template<typename T>
HostReadAccess<T>::~HostReadAccess()
{
    LAMA_LOG_DEBUG( logger, "~HostReadAccess" )
}

template<typename T>
inline const T& HostReadAccess<T>::operator[]( const IndexType i ) const
{
    return mData[i];
}

template<typename T>
inline HostReadAccess<T>::operator const T*() const
{
    LAMA_LOG_TRACE( logger, "mData = " << mData )
    return mData;
}

}

#endif // LAMA_HOSTREADACCESS_HPP_
