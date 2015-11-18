/**
 * @file TypeTraits.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Type specific implementations for arithmetic value types.
 * @author Thomas Brandes
 * @date 17.11.2015
 */

#pragma once

#include <scai/common/Complex.hpp>
#include <scai/common/ScalarType.hpp>

#include <cmath>
#include <Complex.hpp>

namespace scai
{

namespace common
{

/**
 * @brief The template class TypeTraits determines type-specific implementations.
 *
 * @tparam T The type of any arithmetic value type used
 */
template<typename T>
class TypeTraits
{
public:

    static inline T sqrt( T val );
    static inline T abs( T val );

    static inline T getEps()
    {
        return std::numeric_limits<T>::epsilon();
    }

    static inline scalar::ScalarType getValueType()
    {
        return scalar::UNKNOWN;
    }

    static inline const char* id()
    {
        return "UNKNOWN";
    }
};

template<>
class TypeTraits<IndexType>
{
public:

    static inline IndexType abs( IndexType x )
    {
        return ::abs( x );
    }
    static inline IndexType getEps()
    {
        return 0;
    }

    static inline scalar::ScalarType getValueType()
    {
        return scalar::INDEX_TYPE;
    }

    static inline const char* id()
    {
        return "IndexType";
    }
};

template<>
class TypeTraits<long double>
{
public:

    static inline long double sqrt( long double x )
    {
        return ::sqrtl( x );
    }
    static inline long double abs( long double x )
    {
        return ::fabsl( x );
    }
    static inline long double getEps()
    {
        return std::numeric_limits<long double>::epsilon();
    }
    static inline scalar::ScalarType getValueType()
    {
        return scalar::LONG_DOUBLE;
    }

    static inline const char* id()
    {
        return "LongDouble";
    }
};

template<>
class TypeTraits<double>
{
public:

    static inline double sqrt( double x )
    {
        return ::sqrt( x );
    }
    static inline double abs( double x )
    {
        return ::fabs( x );
    }
    static inline double getEps()
    {
        return std::numeric_limits<double>::epsilon();
    }
    static inline scalar::ScalarType getValueType()
    {
        return scalar::DOUBLE;
    }
    static inline const char* id()
    {
        return "double";
    }
};

template<>
class TypeTraits<float>
{
public:

    static inline float sqrt( float x )
    {
        return ::sqrtf( x );
    }
    static inline float abs( float x )
    {
        return ::fabsf( x );
    }
    static inline float getEps()
    {
        return std::numeric_limits<float>::epsilon();
    }
    static inline scalar::ScalarType getValueType()
    {
        return scalar::FLOAT;
    }
    static inline const char* id()
    {
        return "float";
    }
};

template<>
class TypeTraits<ComplexFloat>
{
public:

    static inline ComplexFloat sqrt( ComplexFloat x )
    {
        return scai::common::sqrt( x );
    }
    static inline ComplexFloat abs( ComplexFloat x )
    {
        return scai::common::abs( x );
    }

    static inline ComplexFloat getEps()
    {
        return std::numeric_limits<float>::epsilon();
    }
    static inline scalar::ScalarType getValueType()
    {
        return scalar::COMPLEX;
    }
    static inline const char* id()
    {
        return "ComplexFloat";
    }
};

template<>
class TypeTraits<ComplexDouble>
{
public:

    static inline ComplexDouble sqrt( ComplexDouble x )
    {
        return scai::common::sqrt( x );
    }
    static inline ComplexDouble abs( ComplexDouble x )
    {
        return scai::common::abs( x );
    }
    static inline ComplexDouble getEps()
    {
        return std::numeric_limits<double>::epsilon();
    }
    static inline scalar::ScalarType getValueType()
    {
        return scalar::DOUBLE_COMPLEX;
    }
    static inline const char* id()
    {
        return "ComplexDouble";
    }
};

template<>
class TypeTraits<ComplexLongDouble>
{
public:

    static inline ComplexLongDouble sqrt( ComplexLongDouble x )
    {
        return scai::common::sqrt( x );
    }
    static inline ComplexLongDouble abs( ComplexLongDouble x )
    {
        return scai::common::abs( x );
    }
    static inline ComplexLongDouble getEps()
    {
        return std::numeric_limits<long double>::epsilon();
    }
    static inline scalar::ScalarType getValueType()
    {
        return scalar::LONG_DOUBLE_COMPLEX;
    }
    static inline const char* id()
    {
        return "ComplexLong";
    }
};

}  // namespace common

}  // namespace scai
