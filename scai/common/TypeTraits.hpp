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

#include <scai/common/ScalarType.hpp>
#include <scai/common/Math.hpp>

#include <cmath>
#include <cstdlib>

namespace scai
{

namespace common
{

/**
 * @brief The template class TypeTraits determines type-specific implementations.
 *
 * @tparam T The type of any arithmetic value type used
 */
template<typename ValueType>
class TypeTraits
{
public:

	typedef ValueType PartsType;

    /** Get value-specific epsilon for comparison. */

    static inline ValueType getEps()
    {
        return std::numeric_limits<ValueType>::epsilon();
    }

    /** Get maximal value of a ValueType, used for min, max reductions on arrays. */

    static inline ValueType getMax()
    {
        return std::numeric_limits<ValueType>::max();
    }

    /** Get maximal value of a ValueType, used for min, max reductions on arrays. */

    static inline ValueType getMin()
    {
        return - std::numeric_limits<ValueType>::max();
    }

    /**
     * @brief Corresponding type value of enum ScalarType.
    */
    static const scalar::ScalarType stype = scalar::UNKNOWN;

    /**
     * @brief Return name of the type for identification.
     *
     * Note: this routine returns more meaningful names than typeinfo( T ).name()
     */
    static inline const char* id()
    {
        return scalar2str( scalar::UNKNOWN );
    }

    static inline bool isComplex()
    {
    	return false;
    }
};

/** Type specific traits for IndexType */

template<>
class TypeTraits<IndexType>
{
public:
	typedef IndexType PartsType;

    static inline IndexType getEps()
    {
        return 0;
    }

    static inline IndexType getMax()
    {
        return std::numeric_limits<IndexType>::max();
    }

    static inline IndexType getMin()
    {
        return - std::numeric_limits<IndexType>::max();
    }

    static const scalar::ScalarType stype = scalar::INDEX_TYPE;

    static inline const char* id()
    {
        return scalar2str( scalar::INDEX_TYPE );
    }

    static inline bool isComplex()
	{
		return false;
	}
};

/** Type specific traits for long double */

template<>
class TypeTraits<long double>
{
public:
	typedef long double PartsType;

    static inline long double getEps()
    {
        return std::numeric_limits<long double>::epsilon();
    }

    static inline long double getMax()
    {
        return std::numeric_limits<long double>::max();
    }

    static inline long double getMin()
    {
        return - std::numeric_limits<long double>::max();
    }

    static const scalar::ScalarType stype = scalar::LONG_DOUBLE;

    static inline const char* id()
    {
        return scalar2str( scalar::LONG_DOUBLE );
    }

    static inline bool isComplex()
	{
		return false;
	}
};

/** Type specific traits for double */

template<>
class TypeTraits<double>
{
public:
	typedef double PartsType;

    static inline double getEps()
    {
        return std::numeric_limits<double>::epsilon();
    }
    static inline double getMax()
    {
        return std::numeric_limits<double>::max();
    }
    static inline double getMin()
    {
        return - std::numeric_limits<double>::max();
    }

    static const scalar::ScalarType stype = scalar::DOUBLE;

    static inline const char* id()
    {
        return scalar2str( scalar::DOUBLE );
    }

    static inline bool isComplex()
	{
		return false;
	}
};

/** Type specific traits for float */

template<>
class TypeTraits<float>
{
public:
	typedef float PartsType;

    static inline float getEps()
    {
        return std::numeric_limits<float>::epsilon();
    }
    static inline float getMax()
    {
        return std::numeric_limits<float>::max();
    }
    static inline float getMin()
    {
        return - std::numeric_limits<float>::max();
    }

    static const scalar::ScalarType stype = scalar::FLOAT;

    static inline const char* id()
    {
        return scalar2str( stype );
    }

    static inline bool isComplex()
	{
		return false;
	}
};

#ifdef SCAI_COMPLEX_SUPPORTED

/** Type specific traits for complex(float) */

template<>
class TypeTraits<ComplexFloat>
{
public:
	typedef float PartsType;

    static inline ComplexFloat getEps()
    {
        return std::numeric_limits<float>::epsilon();
    }
    static inline ComplexFloat getMax()
    {
        return std::numeric_limits<float>::max();
    }
    static inline ComplexFloat getMin()
    {
        return 0;
    }

    static const scalar::ScalarType stype = scalar::COMPLEX;

    static inline const char* id()
    {
        return "ComplexFloat";
    }

    static inline bool isComplex()
	{
		return true;
	}
};

/** Type specific traits for complex(double) */

template<>
class TypeTraits<ComplexDouble>
{
public:
	typedef double PartsType;

    static inline ComplexDouble getEps()
    {
        return std::numeric_limits<double>::epsilon();
    }
    static inline ComplexDouble getMax()
    {
        return std::numeric_limits<double>::epsilon();
    }
    static inline ComplexDouble getMin()
    {
        return 0;
    }

    static const scalar::ScalarType stype = scalar::DOUBLE_COMPLEX;

    static inline const char* id()
    {
        return scalar2str( scalar::DOUBLE_COMPLEX );
    }

    static inline bool isComplex()
	{
		return true;
	}
};

/** Type specific traits for complex(long double) */

template<>
class TypeTraits<ComplexLongDouble>
{
public:
	typedef long double PartsType;

    static inline ComplexLongDouble getEps()
    {
        return std::numeric_limits<long double>::epsilon();
    }

    static inline ComplexLongDouble getMax()
    {
        return std::numeric_limits<long double>::max();
    }

    static inline ComplexLongDouble getMin()
    {
        return 0;
    }

    static const scalar::ScalarType stype = scalar::LONG_DOUBLE_COMPLEX;

    static inline const char* id()
    {
        return scalar2str( stype );
    }

    static inline bool isComplex()
	{
		return true;
	}
};

#endif

/** For convenience and for compatibility make own routine of getScalarType */

template<typename ValueType> inline scalar::ScalarType getScalarType()
{
    return TypeTraits<ValueType>::stype;
}

}  // namespace common

}  // namespace scai

