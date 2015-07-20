/**
 * @file HostWriteOnlyAccess.hpp
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
 * @brief Definition of a write only access that uses directly the host context.
 * @author Thomas Brandes
 * @date 20.07.2015
 */

#pragma once

// for dll_import
#include <common/config.hpp>

// base classes
#include <memory/HostWriteAccess.hpp>

namespace memory
{

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
inline HostWriteOnlyAccess<ValueType>::HostWriteOnlyAccess( LAMAArray<ValueType>& array ) :

    HostWriteAccess<ValueType>( array, false )

{
}

template<typename ValueType>
inline HostWriteOnlyAccess<ValueType>::HostWriteOnlyAccess( LAMAArray<ValueType>& array, const IndexType size ) :

    HostWriteAccess<ValueType>( array, false )

{
    // keep flag avoids transfer of invalid data, but also it should not be copied for resize

    this->resize( 0 );
    this->resize( size );
}

template<typename ValueType>
inline HostWriteOnlyAccess<ValueType>::~HostWriteOnlyAccess()
{
}

}  // namespace

