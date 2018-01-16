/**
 * @file HostReadAccess.hpp
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
 * @brief Definition of a template class HostReadAccess for reading a HArray.
 * @author Andreas Longva
 * @date 15.01.2018
 */

#include <scai/hmemo/ReadAccess.hpp>

#include <iterator>

namespace scai
{

namespace hmemo
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT HostReadAccess: public ReadAccess<ValueType>
{
public:

    // Typedefs for the Container STL concept. We can only partially
    // satisfy the concept as we offer only a read-only view.
    typedef ValueType                               value_type;
    typedef const ValueType&                        const_reference;
    typedef const ValueType *                       const_iterator;
    typedef std::iterator_traits<const_iterator>    difference_type;

    explicit HostReadAccess( const HArray<ValueType> & array)
        :   ReadAccess<ValueType>(array, Context::getHostPtr())
    { }

    HostReadAccess( HostReadAccess<ValueType> && other )
        : ReadAccess<ValueType>( std::move(other) )
    { }

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
        return ReadAccess<ValueType>::get();
    }
    const_iterator cend() const
    {
        return ReadAccess<ValueType>::get() + ReadAccess<ValueType>::size();
    }
};

template <typename ValueType>
HostReadAccess<ValueType> hostReadAccess(const HArray<ValueType> & array)
{
    return HostReadAccess<ValueType>(array);
}

}

}

#pragma once
