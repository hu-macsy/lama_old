/**
 * @file HostWriteAccess.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief HostWriteAccess.hpp
 * @author brandes
 * @date 02.05.2011
 * $Id$
 */
#ifndef LAMA_HOSTWRITEACCESS_HPP_
#define LAMA_HOSTWRITEACCESS_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/WriteAccess.hpp>

// others
#include <lama/ContextFactory.hpp>

namespace lama
{

/**
 * @brief HostWriteAccess is a specialization of WriteAccess for the Host Location with an extended interface.
 */
template<typename T>
class HostWriteAccess: public WriteAccess<T>
{
public:

    /**
     * @brief ValueType is the type stored in the wrapped container.
     */
    typedef T ValueType;

    /**
     * @brief acquire a WriteAccess to the passed LAMAArray for the Host location
     *
     * @param[in] array     the LAMAArray to acquire a WriteAccess for
     * @throws Exception    if the WriteAccess can not be acquired, e.g. because another WriteAccess exists.
     *
     * @see WriteAccess for more details
     */
    HostWriteAccess( LAMAArray<ValueType>& array );

    /**
     * @brief acquire a WriteAccess to the passed LAMAArray for the Host location
     *
     * @param[in] array     the LAMAArray to acquire a WriteAccess for
     * @param[in] size      TODO[doxy] Complete Description.
     * @param[in] keep      if the contents of the LAMAArray should be kept or not (default: true)
     * @throws Exception    if the WriteAccess can not be acquired, e.g. because another WriteAccess exists.
     *
     * @see WriteAccess for more details
     */
    HostWriteAccess( LAMAArray<ValueType>& array, const IndexType size, const bool keep );

    /**
     * @brief acquire a WriteAccess to the passed LAMAArray for the Host location
     *
     * @param[in] view      TODO[doxy] Complete Description.
     * @param[in] keep      if the contents of the LAMAArray should be kept or not (default: true)
     * @throws Exception    if the WriteAccess can not be acquired, e.g. because another WriteAccess exists.
     *
     * @see WriteAccess for more details
     */
    HostWriteAccess( LAMAArrayView<ValueType>& view, const bool keep = true );

    /**
     * @brief Releases the WriteAccess on the associated LAMAArray.
     */
    virtual ~HostWriteAccess();

    /**
     * @brief access to the element i of the wrapped LAMAArray
     *
     * @param[in] i   the index of the element of the wrapped LAMAArray to access
     * @return        a reference to the element i of the wrapped LAMAArray
     */
    inline ValueType& operator[]( const IndexType i );

    /**
     * @brief conversion to a ValueType pointer (shortcut for get() )
     *
     * @return  a pointer to the data of the wrapped LAMAArray
     */
    inline operator ValueType*();

    /**
     * @brief appends val to the end of the wrapped LAMAArray.
     *
     * @param[in] val   the value to append
     */
    void push_back( const ValueType val );

    using WriteAccess<ValueType>::size;

    using WriteAccess<ValueType>::resize;

    using WriteAccess<ValueType>::reserve;

protected:

    using WriteAccess<ValueType>::mData;
};

/**
 * @brief HostWriteOnlyAccess is a write access where no existing values of the array are needed (keepFlag = false).
 *
 * This derived class has been added for more convenience as it avoids the use of the keepFlag param.
 *
 * A HostWriteOnlyAccess should be used whenever possible. It avoids any memory transfer of no more
 * needed values between devices and in case of a reallocation it avoids copying of old values.
 */
template<typename T>
class HostWriteOnlyAccess: public HostWriteAccess<T>
{
public:

    /**
     * @brief ValueType is the type stored in the wrapped container.
     */
    typedef T ValueType;

    /** Create a write access with keep flag = false. */

    HostWriteOnlyAccess( LAMAArray<ValueType>& array )

        : HostWriteAccess<ValueType>( array, false )

    {
    }

    /** Create a write access with keep flag = false and do also a resize. */

    HostWriteOnlyAccess( LAMAArray<ValueType>& array, const IndexType size )

        : HostWriteAccess<ValueType>( array, size, false )

    {
    }

    ~HostWriteOnlyAccess()
    {
    }
};

template<typename T>
HostWriteAccess<T>::HostWriteAccess( LAMAArray<ValueType>& array )
    : WriteAccess<T>( array, ContextFactory::getContext( Context::Host ) )
{
}

template<typename T>
HostWriteAccess<T>::HostWriteAccess( LAMAArray<ValueType>& array, const IndexType size, const bool keep )
    : WriteAccess<T>( array, ContextFactory::getContext( Context::Host ), size, keep )
{
}

template<typename T>
HostWriteAccess<T>::HostWriteAccess( LAMAArrayView<ValueType>& view, const bool keep /* = true */)
    : WriteAccess<T>( view, ContextFactory::getContext( Context::Host ), keep )
{
}

template<typename T>
HostWriteAccess<T>::~HostWriteAccess()
{
}

template<typename T>
T& HostWriteAccess<T>::operator[]( const IndexType i )
{
    LAMA_ASSERT_ERROR( mData,
                       "[" << i << "]: HostWriteAccess has already" << " been released or has not been allocated." );

    return mData[i];
}

template<typename T>
inline HostWriteAccess<T>::operator ValueType*()
{
    return mData;
}

template<typename T>
void HostWriteAccess<T>::push_back( const ValueType val )
{
    const IndexType currentSize = size();
    resize( currentSize + 1 );
    mData[currentSize] = val;
}

}

#endif // LAMA_HOSTWRITEACCESS_HPP_
