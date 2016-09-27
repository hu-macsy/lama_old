/**
 * @file TypeTraits.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @endlicense
 *
 * @brief Type specific implementations for arithmetic value types.
 * @author Thomas Brandes
 * @date 17.11.2015
 */

#pragma once

#include <scai/common/ScalarType.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/Math.hpp>

#include <cmath>
#include <cstdlib>
#include <limits>
#include <typeinfo>

namespace scai
{

namespace common
{

/**
 * @brief The template class TypeTraits determines type-specific implementations.
 *
 * @tparam T The type of any arithmetic value type used
 *
 * /code
 *  typedef ResultType AbsType;   // return type of function Math::abs( ValueType )
 *  AbsType eps0();               // closest value to 0
 *  AbsType eps1();               // closest value to 1
 *  AbsType small();              // accepted tolerance for comparisons of results
 *  ValueType min();              // used as neutral element for max reduction
 *  ValueType max();              // used as neutral element for min reduction
 *  int precision();              // number of relevant digits if printed 
 * /endcode
 */
template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT TypeTraits
{
public:
    /**
     * @brief Return name of the type for identification.
     *
     * Note: this routine returns more meaningful names than typeinfo( T ).name()
     */
    static inline const char* id()
    {
        return typeid( ValueType ).name();
    }

    /**
     * @brief Corresponding type value of enum ScalarType.
    */
    static const scalar::ScalarType stype = scalar::UNKNOWN;
};

/** Type specific traits for int */

template<>
class TypeTraits<int>
{
public:

    typedef int AbsType;

    static inline int eps1()
    {
        return 0;
    }

    static inline int small()
    {
        return 0;
    }

    static inline int eps0()
    {
        return std::numeric_limits<int>::min();
    }

    static inline int getMax()
    {
        return std::numeric_limits<int>::max();
    }

    static inline int precision()
    {
        return 0;
    }

    static inline int getMin() // with C++11: use lowest instead
    {
        return - std::numeric_limits<int>::max();
    }

    static const scalar::ScalarType stype = scalar::INT;

    static inline const char* id()
    {
        return scalarType2str( scalar::INT );
    }
};

/** Type specific traits for char */

template<>
class TypeTraits<char>
{
public:

    typedef char AbsType;

    static inline char eps1()
    {
        return 0;
    }

    static inline char small()
    {
        return 0;
    }

    static inline char eps0()
    {
        return std::numeric_limits<char>::min();
    }

    static inline char getMax()
    {
        return std::numeric_limits<char>::max();
    }

    static inline int precision()
    {
        return 0;
    }

    static inline char getMin() // with C++11: use lowest instead
    {
        return static_cast<char> (- std::numeric_limits<char>::max() );
    }

    static const scalar::ScalarType stype = scalar::CHAR;

    static inline const char* id()
    {
        return scalarType2str( scalar::CHAR );
    }
};

/** Type specific traits for long */

template<>
class TypeTraits<long>
{
public:
    typedef long AbsType;

    static inline long eps1()
    {
        return 0;
    }

    static inline long small()
    {
        return 0;
    }

    static inline long eps0()
    {
        return std::numeric_limits<long>::min();
    }

    static inline long getMax()
    {
        return std::numeric_limits<long>::max();
    }

    static inline int precision()
    {
        return 0;
    }

    static inline long getMin() // with C++11: use lowest instead
    {
        return - std::numeric_limits<long>::max();
    }

    static const scalar::ScalarType stype = scalar::LONG;

    static inline const char* id()
    {
        return scalarType2str( scalar::LONG );
    }
};

/** Type specific traits for unsigned int */

template<>
class TypeTraits<unsigned int>
{
public:
    typedef unsigned int AbsType;

    static inline unsigned int eps1()
    {
        return 0;
    }

    static inline unsigned int small()
    {
        return 0;
    }

    static inline unsigned int eps0()
    {
        return std::numeric_limits<unsigned int>::min();
    }

    static inline unsigned int getMax()
    {
        return std::numeric_limits<unsigned int>::max();
    }

    static inline int precision()
    {
        return 0;
    }

    static inline unsigned int getMin() // with C++11: use lowest instead
    {
        return 0;
    }

    static const scalar::ScalarType stype = scalar::UNSIGNED_INT;

    static inline const char* id()
    {
        return scalarType2str( scalar::UNSIGNED_INT );
    }
};

/** Type specific traits for unsigned long */

template<>
class TypeTraits<unsigned long>
{
public:
    typedef unsigned long AbsType;

    static inline unsigned long eps1()
    {
        return 0;
    }

    static inline unsigned long small()
    {
        return 0;
    }

    static inline unsigned long eps0()
    {
        return std::numeric_limits<unsigned long>::min();
    }

    static inline unsigned long getMax()
    {
        return std::numeric_limits<unsigned long>::max();
    }

    static inline int precision()
    {
        return 0;
    }

    static inline unsigned long getMin() // with C++11: use lowest instead
    {
        return 0;
    }

    static const scalar::ScalarType stype = scalar::UNSIGNED_LONG;

    static inline const char* id()
    {
        return scalarType2str( scalar::UNSIGNED_LONG );
    }
};

/** Type specific traits for long double */

template<>
class TypeTraits<long double>
{
public:
    typedef long double AbsType;

    static inline long double eps1()
    {
        return std::numeric_limits<long double>::epsilon();
    }

    static inline long double small()
    {
        return 1e-8L;
    }

    static inline long double eps0()
    {
        return std::numeric_limits<long double>::min();
    }

    static inline long double getMax()
    {
        return std::numeric_limits<long double>::max();
    }

    static inline int precision()
    {
        return std::numeric_limits<long double>::digits10;
    }

    static inline long double getMin() // with C++11: use lowest instead
    {
        return - std::numeric_limits<long double>::max();
    }

    static const scalar::ScalarType stype = scalar::LONG_DOUBLE;

    static inline const char* id()
    {
        return scalarType2str( scalar::LONG_DOUBLE );
    }
};

/** Type specific traits for double */

template<>
class TypeTraits<double>
{
public:
    typedef double AbsType;

    static inline double eps1()
    {
        return std::numeric_limits<double>::epsilon();
    }

    static inline double small()
    {
        return 1e-5;
    }

    static inline double eps0()
    {
        return std::numeric_limits<double>::min();
    }

    static inline double getMax()
    {
        return std::numeric_limits<double>::max();
    }

    static inline int precision()
    {
        return std::numeric_limits<double>::digits10;
    }

    static inline double getMin() // with C++11: use lowest instead
    {
        return - std::numeric_limits<double>::max();
    }

    static const scalar::ScalarType stype = scalar::DOUBLE;

    static inline const char* id()
    {
        return scalarType2str( scalar::DOUBLE );
    }
};

/** Type specific traits for float */

template<>
class TypeTraits<float>
{
public:
    typedef float AbsType;

    static inline float eps1()
    {
        return std::numeric_limits<float>::epsilon();
    }

    static inline float small()
    {
        return 1e-3f;
    }

    static inline float eps0()
    {
        return std::numeric_limits<float>::min();
    }

    static inline int precision()
    {
        return std::numeric_limits<float>::digits10;
    }

    static inline float getMax()
    {
        return std::numeric_limits<float>::max();
    }
    static inline float getMin() // with C++11: use lowest instead
    {
        return - std::numeric_limits<float>::max();
    }

    static const scalar::ScalarType stype = scalar::FLOAT;

    static inline const char* id()
    {
        return scalarType2str( stype );
    }
};

#ifdef SCAI_COMPLEX_SUPPORTED

/** Type specific traits for complex(float) */

template<>
class TypeTraits<ComplexFloat>
{
public:
    typedef float AbsType;

    static inline float eps1()
    {
        return std::numeric_limits<float>::epsilon();
    }

    static inline float small()
    {
        return ComplexFloat( 1e-3f );
    }

    static inline float eps0()
    {
        return std::numeric_limits<float>::min();
    }

    static inline int precision()
    {
        return std::numeric_limits<float>::digits10;
    }

    static inline ComplexFloat getMax()
    {
        return std::numeric_limits<float>::max();
    }

    static inline ComplexFloat getMin() // with C++11: use lowest instead
    {
        return 0;
    }

    static const scalar::ScalarType stype = scalar::COMPLEX;

    static inline const char* id()
    {
        return "ComplexFloat";
    }
};

/** Type specific traits for complex(double) */

template<>
class TypeTraits<ComplexDouble>
{
public:
    typedef double AbsType;

    static inline double eps1()
    {
        return std::numeric_limits<double>::epsilon();
    }

    static inline double small()
    {
        return double( 1e-5 );
    }

    static inline double eps0()
    {
        return std::numeric_limits<double>::min();
    }

    static inline int precision()
    {
        return std::numeric_limits<double>::digits10;
    }

    static inline ComplexDouble getMax()
    {
        return std::numeric_limits<double>::max();
    }

    static inline ComplexDouble getMin() // with C++11: use lowest instead
    {
        return 0;
    }

    static const scalar::ScalarType stype = scalar::DOUBLE_COMPLEX;

    static inline const char* id()
    {
        return scalarType2str( scalar::DOUBLE_COMPLEX );
    }
};

/** Type specific traits for complex(long double) */

template<>
class TypeTraits<ComplexLongDouble>
{
public:
    typedef long double AbsType;

    static inline LongDouble eps1()
    {
        return std::numeric_limits<long double>::epsilon();
    }

    static inline LongDouble small()
    {
        return LongDouble( 1e-8L );
    }

    static inline LongDouble eps0()
    {
        return std::numeric_limits<long double>::min();
    }

    static inline int precision()
    {
        return std::numeric_limits<long double>::digits10;
    }

    static inline ComplexLongDouble getMax()
    {
        return std::numeric_limits<long double>::max();
    }

    static inline ComplexLongDouble getMin() // with C++11: use lowest instead
    {
        return 0;
    }

    static const scalar::ScalarType stype = scalar::LONG_DOUBLE_COMPLEX;

    static inline const char* id()
    {
        return scalarType2str( stype );
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

