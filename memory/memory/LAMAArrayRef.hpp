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
 * @date 14.03.2011
 * @revised 03.07.2015
 */

#pragma once

#include <memory/LAMAArray.hpp>

/** Number of contexts that might be used in maximum. This number
 *  is used for reservation of entries but does not imply any restrictions.
 */

#define LAMA_MAX_CONTEXTS 4

namespace memory
{

/**
 * @brief LAMAArrayRef is a container that uses external data.
 *
 * @tparam ValueType is the type stored in this container.
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT LAMAArrayRef: public LAMAArray<ValueType>
{
public:

    /** Contruct a container for a host array. */

    LAMAArrayRef( ValueType* pointer, IndexType size );

    /** Contruct a container for a const host array.
     *  Due to the const pointer it is guaranteed that the array cannot be modified
     */

    LAMAArrayRef( const ValueType* pointer, IndexType size );

protected:

    using LAMAArray<ValueType>::mSize;

    using LAMAArray<ValueType>::mContextManager;
    using LAMAArray<ValueType>::constFlag;
};

/* ---------------------------------------------------------------------------------*/

template<typename ValueType>
LAMAArrayRef<ValueType>::LAMAArrayRef( ValueType* pointer, IndexType size )
                : LAMAArray<ValueType>()
{
    // Important: context must be set to the DefaultHostContext

    if( size != 0 && pointer == NULL )
    {
        COMMON_THROWEXCEPTION( "LAMAArryRef with NULL pointer" )
    }

/*
    ContextData& host = *mContextData[0];
    host.setRef( pointer, size * sizeof(ValueType) );
*/

    mSize = size;
}

template<typename ValueType>
LAMAArrayRef<ValueType>::LAMAArrayRef( const ValueType* pointer, IndexType size )
                : LAMAArray<ValueType>()
{
    // Important: context must be set to the DefaultHostContext

    if( size != 0 && pointer == NULL )
    {
        COMMON_THROWEXCEPTION( "LAMAArryRef with NULL pointer" )
    }

/*
    ContextData& host = *mContextData[0];
    host.setRef( const_cast<ValueType*>( pointer ), size * sizeof(ValueType) );

    constFlag = true; // makes sure that we cannot have a WriteAccess

    mSize = size;
*/

}

}  // namespace 
