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

#include <cmath>
#include <limits>

namespace scai
{

namespace common
{

namespace constants
{

enum ConstantType
{
    ONE,
    ZERO,
    EPS
};

template<typename ValueType>
inline ValueType getEps()
{
	return std::numeric_limits<ValueType>::epsilon();
}

template<>
inline ComplexFloat getEps()
{
	return std::numeric_limits<float>::epsilon();
}

template<>
inline ComplexDouble getEps()
{
	return std::numeric_limits<double>::epsilon();
}

template<>
inline ComplexLongDouble getEps()
{
	return std::numeric_limits<long double>::epsilon();
}

template<typename ValueType>
inline ValueType getMax()
{
	return std::numeric_limits<ValueType>::max();
}

template<>
inline ComplexFloat getMax()
{
	return std::numeric_limits<float>::max();
}

template<>
inline ComplexDouble getMax()
{
	return std::numeric_limits<double>::max();
}

template<>
inline ComplexLongDouble getMax()
{
	return std::numeric_limits<long double>::max();
}

template<typename ValueType>
inline ValueType getConstant( const enum ConstantType& c )
{
    ValueType val = getMax<ValueType>();
    switch( c )
    {
        case ONE:
            val = ValueType( 1.0 );
            break;
        case ZERO:
            val = ValueType( 0.0 );
            break;
        case EPS:
        	val = getEps<ValueType>();
            break;
    }
    return val;
}

template<typename ValueType>
bool operator==( const ValueType& x, const enum ConstantType& c )
{
    using std::fabs;
    return fabs( x - getConstant<ValueType>( c ) ) < getConstant<ValueType>( EPS );
}

template<typename ValueType>
bool operator!=( const ValueType& x, const enum ConstantType& c )
{
    using std::fabs;
    return !(fabs( x - getConstant<ValueType>( c ) ) < getConstant<ValueType>( EPS ));
}

} /* end namespace constants */

} /* end namespace common */

} /* end namespace scai */
