/**
 * @file HostWriteAccess.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 *
 * @endlicense
 *
 * @brief Definition of a template class HostWriteAccess for writing to an HArray.
 * @author Andreas Longva
 * @date 15.01.2018
 */

#pragma once

#include <scai/hmemo/WriteAccess.hpp>

#include <iterator>

namespace scai
{

namespace hmemo
{

/**
 * @brief  Acquire a WriteAccess on the host context.
 *
 * HostWriteAccess is exactly equivalent to WriteAccess, except that it can only
 * acquire a WriteAccess to the Host context. This guarantees that the data pointer
 * resides in host memory, which enables us to provide some functionality which we
 * could not safely do if the memory were allowed to reside anywhere else.
 * In particular, `begin`/`cbegin` and `end`/`cend` methods are provided
 * for interoperability with the C++ standard library.
 *
 * This enables users to more easily use HArray with algorithms from the STL,
 * as well as enabling features such as range-based for loops when iterating over
 * arrays. See the documentation for HostReadAccess for more information.
 */

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT HostWriteAccess: public WriteAccess<ValueType>
{
public:

    // Typedefs for the Container STL concept. We can only partially
    // satisfy the concept as we offer only a read-only view.
    typedef ValueType                               value_type;
    typedef const ValueType&                        const_reference;
    typedef ValueType*                              iterator;
    typedef const ValueType*                        const_iterator;
    typedef std::iterator_traits<const_iterator>    difference_type;

    /**
     * @brief Obtain write access to the given array on the host context.
     *
     * This is exactly equivalent to WriteAccess(array, context::getHostPtr());
     */
    explicit HostWriteAccess( HArray<ValueType>& array )
        :   WriteAccess<ValueType>( array, Context::getHostPtr(), true )
    { }

    /**
     * @brief Move constructor for HostWriteAccess.
     */
    HostWriteAccess( HostWriteAccess<ValueType>&& other )
        : WriteAccess<ValueType>( std::move( other ) )
    { }

    iterator begin()
    {
        return WriteAccess<ValueType>::get();
    }

    iterator end()
    {
        return WriteAccess<ValueType>::get() + WriteAccess<ValueType>::size();
    }

    const_iterator begin() const
    {
        return cbegin();
    }

    const_iterator end() const
    {
        return cend();
    }

    const_iterator cbegin() const
    {
        return WriteAccess<ValueType>::get();
    }
    const_iterator cend() const
    {
        return WriteAccess<ValueType>::get() + WriteAccess<ValueType>::size();
    }
};

/**
 * @return Return a HostWriteAccess for the provided array.
 *
 * This is exactly equivalent to writeAccess(const HArray<ValueType> & array), except
 * that a HostWriteAccess is returned instead of WriteAccess.
 */
template <typename ValueType>
HostWriteAccess<ValueType> hostWriteAccess( HArray<ValueType>& array )
{
    return HostWriteAccess<ValueType>( array );
}

}

}
