/**
 * @file Scalar.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Scalar.hpp
 * @author Jiri Kraus
 * @date 22.02.2011
 * $Id$
 */
#ifndef LAMA_SCALAR_HPP_
#define LAMA_SCALAR_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/Printable.hpp>

// others
#include <lama/LAMATypes.hpp>

#include <lama/exception/LAMAAssert.hpp>

// logging
#include <logging/logging.hpp>

#include <complex>
#include <cmath>
#include <cstdio>

namespace lama
{

/**
 * @brief The class Scalar represents a multi precision scalar.
 */
class LAMA_DLL_IMPORTEXPORT Scalar: public Printable
{
public:

    enum ScalarType
    {
        INDEX_TYPE, FLOAT, DOUBLE, LONG_DOUBLE, COMPLEX, DOUBLE_COMPLEX, LONG_DOUBLE_COMPLEX, UNKNOWN
    };

    /**
     * @brief ExpressionMemberType is the type that is used the template Expression to store a Scalar.
     */
    typedef const Scalar ExpressionMemberType;

    /**
     * @brief Constructs a scalar representing 0.
     */
    inline Scalar();

    /**
     * @brief Constructs a scalar representing the passed real value.
     *
     * The templated converstion constructor needs to be explicit, because the operator==(Scalar,Scalar) can lead to ambiguities.
     *
     * @tparam T          TODO[doxy] Complete Description.
     * @param[in] value   the value this scalar should represent
     */
    template<typename T>
    explicit inline Scalar( const T value );

    /**
     * @brief Constructs a scalar representing the passed real value.
     *
     * @param[in] value the value this scalar should represent
     */
    inline Scalar( const float value );

    /**
     * @brief Constructs a scalar representing the passed real value.
     *
     * @param[in] value the value this scalar should represent
     */
    inline Scalar( const double value );

    /**
     * @brief Constructs a scalar representing the passed real value.
     *
     * @param[in] value the value this scalar should represent
     */
    inline Scalar( const long double value );

    /**
     * @brief Constructs a scalar representing the passed complex value.
     *
     * @tparam T    TODO[doxy] Complete Description.
     * @param[in]   value the value this scalar should represent
     */
    template<typename T>
    inline Scalar( const std::complex<T> value );

    /**
     * @brief Releases all allocated resources.
     */
    inline virtual ~Scalar();

    /**
     * @brief Returns the value this Scalar represents as type T.
     *
     * @tparam T    TODO[doxy] Complete Description.
     * @return      the value this Scalar represents.
     */
    template<typename T>
    inline T getValue() const;

    // Removed: template<typename T> operator T () const  for type conversions
    // Might cause problems due to implicit conversions, should only be used explicitly
    // Now should be done as: cast<T>( scalar )

    /**
     * @brief Unary minus operator for Scalar.
     */
    inline Scalar operator-() const;

    /**
     * @brief Binary operator
     */
    Scalar& operator+=( Scalar& other );
    Scalar& operator-=( Scalar& other );
    Scalar& operator*=( Scalar& other );
    Scalar& operator/=( Scalar& other );

    /**
     *  @brief Query that scalar values has no imaginary part.
     */
    inline bool isReal() const;

    inline virtual void writeAt( std::ostream& stream ) const;

    /**
     * @brief Returns the type this Scalar represents.
     *
     * @tparam T    TODO[doxy] Complete Description.
     * @return      the ScalarType
     */
    template<typename T>
    inline static ScalarType getType();

    /**
     * @brief Returns the size of the given ScalarType.
     *
     * @param[in] type    the given ScalarType.
     * @return            the size of the given ScalarType.
     */
    inline static size_t getTypeSize( const ScalarType type );

protected:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

private:

    std::complex<long double> mValue;
};

LAMA_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const Scalar::ScalarType& object );

/** @brief Cast operator to convert a Scalar into corresponding basic type.
 *
 *  This solutions avoids an implicit conversion of a Scalar to a basic type.
 */

template<typename T>
T cast( const Scalar& scalar )
{
    return scalar.getValue<T>();
}

const Scalar zero;

inline Scalar::Scalar()
    : mValue( 0.0, 0.0 )
{
}

template<typename T>
inline Scalar::Scalar( const T value )
    : mValue( value, 0.0 )
{
}

inline Scalar::Scalar( const float value )
    : mValue( value, 0.0 )
{
}

inline Scalar::Scalar( const double value )
    : mValue( value, 0.0 )
{
}

inline Scalar::Scalar( const long double value )
    : mValue( value, 0.0 )
{
}

template<typename T>
inline Scalar::Scalar( const std::complex<T> value )
    : mValue( value.real(), value.imag() )
{
}

inline Scalar::~Scalar()
{
}

template<typename T>
inline T Scalar::getValue() const
{
    return static_cast<T>( mValue.real() );
}

template<>
inline std::complex<float> Scalar::getValue() const
{
    return std::complex<float>( static_cast<float>( mValue.real() ), static_cast<float>( mValue.imag() ) );
}

template<>
inline std::complex<double> Scalar::getValue() const
{
    return std::complex<double>( static_cast<double>( mValue.real() ), static_cast<double>( mValue.imag() ) );
}

template<>
inline std::complex<long double> Scalar::getValue() const
{
    return std::complex<long double>( mValue.real(), mValue.imag() );
}

inline Scalar Scalar::operator-() const
{
    return Scalar( -mValue );
}

inline bool Scalar::isReal() const
{
    return mValue.imag() == 0.0;
}

inline void Scalar::writeAt( std::ostream& stream ) const
{
    stream << "Scalar(" << mValue.real() << ")";
}

template<typename T>
inline Scalar::ScalarType Scalar::getType()
{
    return UNKNOWN;
}

template<>
inline Scalar::ScalarType Scalar::getType<IndexType>()
{
    return INDEX_TYPE;
}

template<>
inline Scalar::ScalarType Scalar::getType<float>()
{
    return FLOAT;
}

template<>
inline Scalar::ScalarType Scalar::getType<double>()
{
    return DOUBLE;
}

template<>
inline Scalar::ScalarType Scalar::getType<long double>()
{
    return LONG_DOUBLE;
}

template<>
inline Scalar::ScalarType Scalar::getType<std::complex<float> >()
{
    return COMPLEX;
}

template<>
inline Scalar::ScalarType Scalar::getType<std::complex<double> >()
{
    return DOUBLE_COMPLEX;
}

template<>
inline Scalar::ScalarType Scalar::getType<std::complex<long double> >()
{
    return LONG_DOUBLE_COMPLEX;
}

inline size_t Scalar::getTypeSize( const ScalarType type )
{
    size_t typeSize = 0;

    switch ( type )
    {
    case FLOAT:
        typeSize = 4;
        break;

    case DOUBLE:
        typeSize = 8;
        break;

    case LONG_DOUBLE:
        typeSize = 16;
        break;

    case COMPLEX:
        typeSize = 8;
        break;

    case DOUBLE_COMPLEX:
        typeSize = 16;
        break;

    case LONG_DOUBLE_COMPLEX:
        typeSize = 32;
        break;

    default:
        typeSize = 0;
    }

    return typeSize;
}

/**
 * @brief Add Scalar a with Scalar b
 *
 * @param[in] a     1st Scalar.
 * @param[in] b     2nd Scalar.
 * @return          sum
 */
inline Scalar operator+( const Scalar& a, const Scalar& b )
{
    return Scalar( a.getValue<long double>() + b.getValue<long double>() );
}

/**
 * @brief Subtract Scalar a with Scalar b
 *
 * @param[in] a     1st Scalar.
 * @param[in] b     2nd Scalar.
 * @return          difference
 */
inline Scalar operator-( const Scalar& a, const Scalar& b )
{
    return Scalar( a.getValue<long double>() - b.getValue<long double>() );
}

/**
 * @brief Multiply Scalar a with Scalar b
 *
 * @param[in] a     1st Scalar.
 * @param[in] b     2nd Scalar.
 * @return          product
 */
inline Scalar operator*( const Scalar& a, const Scalar& b )
{
    return Scalar( a.getValue<long double>() * b.getValue<long double>() );
}

/**
 * @brief Divide Scalar a with Scalar b
 *
 * @param[in] a     1st Scalar.
 * @param[in] b     2nd Scalar.
 * @return          quotient
 */
inline Scalar operator/( const Scalar& a, const Scalar& b )
{
    return Scalar( a.getValue<long double>() / b.getValue<long double>() );
}

/**
 * @brief Check equality of a and b.
 *
 * @param[in] a     the 1st Scalar to compare this to.
 * @param[in] b     the 2nd Scalar to compare this to.
 * @return          if a is equal to b
 */
inline bool operator==( const Scalar& a, const Scalar& b )
{
    return a.getValue<long double>() == b.getValue<long double>();
}

/**
 * @brief Check inequality of a and b.
 *
 * @param[in] a     the 1st Scalar to compare this to.
 * @param[in] b     the 2nd Scalar to compare this to.
 * @return          if a is unequal to b
 */
inline bool operator!=( const Scalar& a, const Scalar& b )
{
    return !( a == b );
}

inline bool operator<( const Scalar& a, const Scalar& b )
{
    if ( !a.isReal() || !b.isReal() )
    {
        LAMA_THROWEXCEPTION(
            "Could not call operator< for Scalar a = " << a << ", b = " << b << ", because one of them is not real." );
    }

    return a.getValue<long double>() < b.getValue<long double>();
}

inline bool operator>( const Scalar& a, const Scalar& b )
{
    if ( !a.isReal() || !b.isReal() )
    {
        LAMA_THROWEXCEPTION(
            "Could not call operator> for Scalar a = " << a << ", b = " << b << ", because one of them is not real." );
    }

    return a.getValue<long double>() > b.getValue<long double>();
}

inline bool operator<=( const Scalar& a, const Scalar& b )
{
    return !( a > b );
}

inline bool operator>=( const Scalar& a, const Scalar& b )
{
    return !( a < b );
}

inline Scalar sqrt( const Scalar scalar )
{
    return std::sqrt( scalar.getValue<std::complex<long double> >() );
}

inline Scalar abs( const Scalar scalar )
{
    return std::abs( scalar.getValue<std::complex<long double> >() );
}

inline Scalar max( const Scalar a, const Scalar b )
{
    LAMA_ASSERT_DEBUG( a.isReal(), "Non-real value in max : " << a )
    LAMA_ASSERT_DEBUG( b.isReal(), "Non-real value in max : " << b )

    return std::max( a.getValue<long double>(), b.getValue<long double>() );
}

inline Scalar min( const Scalar a, const Scalar b )
{
    LAMA_ASSERT_DEBUG( a.isReal(), "Non-real value in max : " << a )
    LAMA_ASSERT_DEBUG( b.isReal(), "Non-real value in max : " << b )

    return std::min( a.getValue<long double>(), b.getValue<long double>() );
}

} //namespace lama

#endif // LAMA_SCALAR_HPP_
