/**
 * @file HArrayRef.hpp
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
 * @brief Definition of a new dynamic array class where the array data can be
 *        used in different contexts and where the data is moved implicitly
 * @author Thomas Brandes, Jiri Krause
 * @date 03.07.2015
 */

#pragma once

// base classes
#include <scai/hmemo/HArray.hpp>

// local library
#include <scai/hmemo/HostMemory.hpp>

namespace scai
{

namespace hmemo
{

/**
 * @brief HArrayRef is a container that uses already allocated Host memory
 *
 * @tparam ValueType is the type stored in this container.
 *
 * In some situations data is already available at the host and copying the
 * data would cause some loss of performance. This class allows construction
 * of a LAMA array that uses this allocated data.
 *
 * An object of this class is restricted in its use. Any resize operation that
 * would cause reallocation of data at the host throws an exception.
 *
 * If a const pointer/reference is passed, only read accesses will be possible.
 *
 * Instead of using an additional class like ConstHArrayRef this approach has been
 * chosen to avoid additional class definitions. The penalty is that illegal access
 * can only be identified at runtime and not at compile time.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT HArrayRef: public HArray<ValueType>
{
public:

    /** Contruct a container for a host array. */

    HArrayRef( IndexType size, ValueType* pointer );

    /** Contruct a container for a const host array.
     *  Due to the const pointer it is guaranteed that the array cannot be modified
     */

    HArrayRef( IndexType size, const ValueType* pointer );

    HArrayRef( const std::vector<ValueType>& data );

    HArrayRef( std::vector<ValueType>& data );

};

/* ---------------------------------------------------------------------------------*/
/*   Implementation of template methods                                             */
/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArrayRef<ValueType>::HArrayRef( IndexType size, ValueType* pointer )
    : HArray<ValueType>()
{
    // Important: context must be set to the DefaultHostContext
    if ( size != 0 && pointer == NULL )
    {
        COMMON_THROWEXCEPTION( "LAMAArryRef with NULL pointer" )
    }

    _HArray::setHostRef( size, pointer );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArrayRef<ValueType>::HArrayRef( std::vector<ValueType>& data ) : HArray<ValueType>()
{
    IndexType size = data.size();

    if ( size == 0 )
    {
        // in this case the HArray behaves like a usual array
        return;
    }

    _HArray::setHostRef( size, &data[0] );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArrayRef<ValueType>::HArrayRef( IndexType size, const ValueType* pointer ) : 

    HArray<ValueType>()

{
    if ( size != 0 && pointer == NULL )
    {
        COMMON_THROWEXCEPTION( "LAMAArryRef with NULL pointer" )
    }

    // const cast required but const flag is set true

    _HArray::setHostRef( size, pointer );
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArrayRef<ValueType>::HArrayRef( const std::vector<ValueType>& data ) : HArray<ValueType>()
{
    IndexType size = data.size();

    if ( size == 0 )
    {
        // in this case the HArray behaves like a usual array
        return;
    }

    _HArray::setHostRef( size, &data[0] );
}

} /* end namespace hmemo */

} /* end namespace scai */
