/**
 * @file Constants.hpp
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
 * @brief Definition of the Constants One, Zero, Eps for comparisons
 *
 * @author Eric Schricker
 * @date 07.10.2015
 */

#pragma once

#include <scai/common/ScalarType.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

#include <cmath>
#include <limits>

namespace scai
{

namespace common
{

/** Own namespace for the enum type ConstantType and its values. */

namespace constants
{

/** Enumeration type for constants for which type-specific values are provided */

enum ConstantType
{
    ONE,   //!< stands for value 1
    ZERO   //!< stands for value 0
};

/** This method returns the type specific value for each constant */

template<typename ValueType>
inline ValueType getConstant( const enum ConstantType& c )
{
    ValueType val( 0 );

    switch( c )
    {
        case ONE:
            val = ValueType( 1 );
            break;
        case ZERO:
            val = ValueType( 0 );
            break;
    }
    return val;
}

/** Comparison against constant ZERO or ONE uses machine-specific EPS */

template<typename ValueType>
bool operator==( const ValueType& x, const enum ConstantType& c )
{
    return Math::abs( x - getConstant<ValueType>( c ) ) < TypeTraits<ValueType>::getEps();
}

/** Operator not equal also provided for convenience */

template<typename ValueType>
bool operator!=( const ValueType& x, const enum ConstantType& c )
{
    return ! ( x == c );
}

} /* end namespace constants */

} /* end namespace common */

} /* end namespace scai */
