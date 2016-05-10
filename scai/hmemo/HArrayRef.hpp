/**
 * @file HArrayRef.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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

protected:

    using HArray<ValueType>::mSize;
    using HArray<ValueType>::mValueSize;

    using HArray<ValueType>::mContextDataManager;
    using HArray<ValueType>::constFlag;
};

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArrayRef<ValueType>::HArrayRef( IndexType size, ValueType* pointer )
                : HArray<ValueType>()
{
    // Important: context must be set to the DefaultHostContext

    if( size != 0 && pointer == NULL )
    {
        COMMON_THROWEXCEPTION( "LAMAArryRef with NULL pointer" )
    }

    ContextData& host = mContextDataManager[ HostMemory::getIt() ];
    host.setRef( pointer, size * mValueSize );

    mSize = size;
}

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
HArrayRef<ValueType>::HArrayRef( IndexType size, const ValueType* pointer )
                : HArray<ValueType>()
{
    // Important: context must be set to the DefaultHostContext

    if( size != 0 && pointer == NULL )
    {
        COMMON_THROWEXCEPTION( "LAMAArryRef with NULL pointer" )
    }

    ContextData& host = mContextDataManager[ HostMemory::getIt() ];

    // dealing with const references in ContextData is not supported

    host.setRef( const_cast<ValueType*>( pointer ), size * sizeof(ValueType) );

    // Take care of const awareness by setting a flag

    constFlag = true; 

    mSize = size;
}

} /* end namespace hmemo */

} /* end namespace scai */
