/**
 * @file LAMAArrayRef.hpp
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
 * @brief Definition of a new dynamic array class where the array data can be
 *        used in different contexts and where the data is moved implicitly
 *        when corresponding read/write accesses are required
 *
 * @author Thomas Brandes, Jiri Krause
 * @date 03.07.2015
 */

#pragma once

#include <scai/memory/LAMAArray.hpp>
#include <scai/memory/HostMemory.hpp>

namespace memory
{

/**
 * @brief LAMAArrayRef is a container that uses already allocated Host memory
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
class COMMON_DLL_IMPORTEXPORT LAMAArrayRef: public LAMAArray<ValueType>
{
public:

    /** Contruct a container for a host array. */

    LAMAArrayRef( IndexType size, ValueType* pointer );

    /** Contruct a container for a const host array.
     *  Due to the const pointer it is guaranteed that the array cannot be modified
     */

    LAMAArrayRef( IndexType size, const ValueType* pointer );

protected:

    using LAMAArray<ValueType>::mSize;
    using LAMAArray<ValueType>::mValueSize;

    using LAMAArray<ValueType>::mContextDataManager;
    using LAMAArray<ValueType>::constFlag;
};

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArrayRef<ValueType>::LAMAArrayRef( IndexType size, ValueType* pointer )
                : LAMAArray<ValueType>()
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
LAMAArrayRef<ValueType>::LAMAArrayRef( IndexType size, const ValueType* pointer )
                : LAMAArray<ValueType>()
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

}  // namespace 
