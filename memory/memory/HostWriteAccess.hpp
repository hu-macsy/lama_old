/**
 * @file HostWriteAccess.hpp
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
 * @brief Definition of a write access that uses directly the host context.
 * @author Thomas Brandes
 * @date 02.05.2011
 */

#pragma once

// for dll_import
#include <common/config.hpp>

// base classes
#include <memory/WriteAccess.hpp>

namespace memory
{

/**
 * @brief HostWriteAccess is a specialization of WriteAccess for the Host Location with an extended interface.
 *
 * @tparam ValueType is the type stored in the wrapped container.
 */
template<typename ValueType>
class HostWriteAccess: public WriteAccess<ValueType>
{
public:

    /**
     * @brief acquire a WriteAccess to the passed LAMAArray for the host location
     *
     * @param[in] array     the LAMAArray to acquire a WriteAccess for
     * @param[in] keep      if the contents of the LAMAArray should be kept or not (default: true)
     * @throws Exception    if the WriteAccess can not be acquired, e.g. because another WriteAccess exists.
     *
     * @see WriteAccess for more details
     */
    HostWriteAccess( LAMAArray<ValueType>& array, const bool keep = true );

    /**
     * @brief acquire a WriteAccess to the passed LAMAArray for the host location
     *
     * @param[in] array     the LAMAArray to acquire a WriteAccess for
     * @param[in] size      the new size of the LAMA array
     * @param[in] keep      if the contents of the LAMAArray should be kept or not (default: true)
     * @throws Exception    if the WriteAccess can not be acquired, e.g. because another WriteAccess exists.
     *
     * @see WriteAccess for more details
     */
    HostWriteAccess( LAMAArray<ValueType>& array, const IndexType size, const bool keep );

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
    inline operator ValueType* ();

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
 *
 * @tparam ValueType is the value type for an element of this.
 */
template<typename ValueType>
class HostWriteOnlyAccess: public HostWriteAccess<ValueType>
{
public:

    /** Creates a write access with keep flag = false. */

    explicit HostWriteOnlyAccess( LAMAArray<ValueType>& array );

    /** Creates a write access with keep flag = false and do also a resize. */

    HostWriteOnlyAccess( LAMAArray<ValueType>& array, const IndexType size );

    /** Destructor. */

    ~HostWriteOnlyAccess();
};

template<typename ValueType>
HostWriteAccess<ValueType>::HostWriteAccess( LAMAArray<ValueType>& array, const IndexType size, const bool keep ) :

    WriteAccess<ValueType>( array, Context::getContext( context::Host ), size, keep )

{
}

template<typename ValueType>
HostWriteAccess<ValueType>::HostWriteAccess( LAMAArray<ValueType>& array, const bool keep ) :

    WriteAccess<ValueType>( array, Context::getContext( context::Host ), keep )

{
}

template<typename ValueType>
HostWriteAccess<ValueType>::~HostWriteAccess()
{
}

template<typename ValueType>
ValueType& HostWriteAccess<ValueType>::operator[]( const IndexType i )
{
    COMMON_ASSERT( mData, "[" << i << "]: HostWriteAccess has already been released or has not been allocated." )
    return mData[i];
}

template<typename ValueType>
inline HostWriteAccess<ValueType>::operator ValueType* ()
{
    return mData;
}

template<typename ValueType>
inline HostWriteOnlyAccess<ValueType>::HostWriteOnlyAccess( LAMAArray<ValueType>& array ) :

    HostWriteAccess<ValueType>( array, false )

{
}

template<typename ValueType>
inline HostWriteOnlyAccess<ValueType>::HostWriteOnlyAccess( LAMAArray<ValueType>& array, const IndexType size ) :

    HostWriteAccess<ValueType>( array, false )

{
    this->resize( 0 );
    this->resize( size );
}

template<typename ValueType>
inline HostWriteOnlyAccess<ValueType>::~HostWriteOnlyAccess()
{
}

}  // namespace

