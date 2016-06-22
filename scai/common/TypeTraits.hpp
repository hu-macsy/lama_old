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

    /** Return type of function abs( ValueType )
     *
     *  Note: ValueType != TypeTraits<ValueType>::AbsType <=> ValueType is complex
     */

    typedef ValueType AbsType;

    /** Get value-specific epsilon for comparison. */

    static inline ValueType getEps()
    {
        return std::numeric_limits<ValueType>::epsilon();
    }

    /** Get type-specific precision to be used for comparison in matrix-vector operations */

    static inline ValueType small()
    {
        return getEps();
    }

    /** Get maximal value of a ValueType, used for min, max reductions on arrays. */

    static inline ValueType getMax()
    {
        return std::numeric_limits<ValueType>::max();
    }

    /** Get number of digits to be used as default precision for formatted ouptut */

    static inline int precision()
    {
        return std::numeric_limits<ValueType>::digits * std::log10( 2 ) + 2;
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
        return scalarType2str( scalar::UNKNOWN );
    }
};

/** Type specific traits for IndexType */

template<>
class TypeTraits<IndexType>
{
public:
    typedef IndexType AbsType;

    static inline IndexType getEps()
    {
        return 0;
    }

    static inline IndexType small()
    {
        return 0;
    }

    static inline IndexType getMax()
    {
        return std::numeric_limits<IndexType>::max();
    }

    static inline int precision()
    {
        return 0;
    }

    static inline IndexType getMin()
    {
        return - std::numeric_limits<IndexType>::max();
    }

    static const scalar::ScalarType stype = scalar::INDEX_TYPE;

    static inline const char* id()
    {
        return scalarType2str( scalar::INDEX_TYPE );
    }
};

/** Type specific traits for long double */

template<>
class TypeTraits<long double>
{
public:
    typedef long double AbsType;

    static inline long double getEps()
    {
        return std::numeric_limits<long double>::epsilon();
    }

    static inline long double small()
    {
        return 1e-8L;
    }

    static inline long double getMax()
    {
        return std::numeric_limits<long double>::max();
    }

    static inline int precision()
    {
        return std::numeric_limits<long double>::digits * std::log10( 2 ) + 2;
    }

    static inline long double getMin()
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

    static inline double getEps()
    {
        return std::numeric_limits<double>::epsilon();
    }

    static inline double small()
    {
        return 1e-5;
    }

    static inline double getMax()
    {
        return std::numeric_limits<double>::max();
    }

    static inline int precision()
    {
        return std::numeric_limits<double>::digits * std::log10( 2 ) + 2;
    }

    static inline double getMin()
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

    static inline float getEps()
    {
        return std::numeric_limits<float>::epsilon();
    }

    static inline float small()
    {
        return 1e-3f;
    }

    static inline int precision()
    {
        return std::numeric_limits<float>::digits * std::log10( 2 ) + 2;
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

    static inline ComplexFloat getEps()
    {
        return std::numeric_limits<float>::epsilon();
    }

    static inline ComplexFloat small()
    {
        return ComplexFloat( 1e-3f );
    }

    static inline int precision()
    {
        return std::numeric_limits<float>::digits * std::log10( 2 ) + 2;
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
};

/** Type specific traits for complex(double) */

template<>
class TypeTraits<ComplexDouble>
{
public:
    typedef double AbsType;

    static inline ComplexDouble getEps()
    {
        return std::numeric_limits<double>::epsilon();
    }

    static inline ComplexDouble small()
    {
        return ComplexDouble( 1e-5 );
    }

    static inline int precision()
    {
        return std::numeric_limits<double>::digits * std::log10( 2 ) + 2;
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
        return scalarType2str( scalar::DOUBLE_COMPLEX );
    }
};

/** Type specific traits for complex(long double) */

template<>
class TypeTraits<ComplexLongDouble>
{
public:
    typedef long double AbsType;

    static inline ComplexLongDouble getEps()
    {
        return std::numeric_limits<long double>::epsilon();
    }

    static inline ComplexLongDouble small()
    {
        return ComplexLongDouble( 1e-8L );
    }

    static inline int precision()
    {
        return std::numeric_limits<long double>::digits * std::log10( 2 ) + 2;
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

