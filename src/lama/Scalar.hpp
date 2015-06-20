/**
 * @file Scalar.hpp
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
 * @brief Scalar.hpp
 * @author Jiri Kraus
 * @date 22.02.2011
 * @since 1.0.0
 */
#ifndef LAMA_SCALAR_HPP_
#define LAMA_SCALAR_HPP_

// for dll_import
#include <common/config.hpp>

// base classes
#include <common/Printable.hpp>

// others
#include <lama/LAMATypes.hpp>

#include <lama/exception/LAMAAssert.hpp>

// logging
#include <logging/logging.hpp>

#include <cstdio>
#include <lama/Complex.hpp>

namespace lama
{

/**
 * @brief The class Scalar represents a multi precision scalar.
 *
 * An object of the class Scalar is used in LAMA in all code parts
 * that are universal for all arithmetic types, especially for code
 * parts that use book syntax.
 *
 * For a Scalar the arithmetic operations +, -, *, / etc. are
 * also supported to allow a high flexibility. But for efficiency
 * these operations should be avoided in all critical code parts.
 */
class COMMON_DLL_IMPORTEXPORT Scalar: public Printable
{
public:

    /** Enumeration type for supported value types in LAMA.
     *
     *  This enumeration type is result of many query operations for LAMA classes
     *  and avoids expensive calls of the typeid operator.
     *
     *  \code
     *    CSRSparseMatrix<double> a;
     *    a.getValueType()  // returns ScalarType::DOUBLE
     *    Scalar::getType<double>()  // return ScalarType DOUBLE
     *  \endcode
     *
     *  It is especially useful when casting variables of base classes to derived classes.
     */
    enum ScalarType
    {
        INDEX_TYPE, //!<  synonymous for IndexType
        FLOAT, //!<  synonymous for float
        DOUBLE, //!<  synonymous for double
        LONG_DOUBLE, //!<  synonymous for long double
        COMPLEX, //!<  synonymous for complex
        DOUBLE_COMPLEX, //!<  synonymous for double complex
        LONG_DOUBLE_COMPLEX, //!<  synonymous for long double complex
        PATTERN, //!<  dummy type of size 0
        INTERNAL, //!<  take the type currently in use, getType<ValueType>()
        UNKNOWN
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
     * The templated conversion constructor needs to be explicit,
     * because the operator==(Scalar,Scalar) can lead to ambiguities.
     *
     * @tparam ValueType          type of the input argument value for constructor of Scalar
     * @param[in] value   the value this scalar should represent
     */
    template<typename ValueType>
    explicit inline Scalar( const ValueType value );

    /**
     * @brief Constructs a scalar representing the passed real value.
     *
     * @param[in] value the value this scalar should represent
     */
    inline Scalar( const float value );

    /**
     * @brief Constructs a scalar representing the passed real and imaginary value.
     *
     * @param[in] real the real part this scalar should represent
     * @param[in] imag the imaginary part this scalar should represent
     */
    inline Scalar( const float real, const float imag );

    /**
     * @brief Constructs a scalar representing the passed real value.
     *
     * @param[in] value the value this scalar should represent
     */
    inline Scalar( const double value );

    /**
     * @brief Constructs a scalar representing the passed real and imaginary value.
     *
     * @param[in] real the real part this scalar should represent
     * @param[in] imag the imaginary part this scalar should represent
     */
    inline Scalar( const double real, const double imag );

    /**
     * @brief Constructs a scalar representing the passed real value.
     *
     * @param[in] value the value this scalar should represent
     */
    inline Scalar( const LongDouble value );

    /**
     * @brief Constructs a scalar representing the passed real and imaginary value.
     *
     * @param[in] real the real part this scalar should represent
     * @param[in] imag the imaginary part this scalar should represent
     */
    inline Scalar( const LongDouble real, const LongDouble imag );

    /**
     * @brief Constructs a scalar representing the passed complex value.
     *
     * @tparam ValueType    base type of the complex type, e.g. float, double or long double
     * @param[in]   value the value this scalar should represent
     */
    template<typename ValueType>
    inline Scalar( const Complex<ValueType> value );

    /**
     * @brief Releases all allocated resources.
     */
    inline virtual ~Scalar();

    /**
     * @brief Returns the value this Scalar represents as type ValueType.
     *
     * @tparam ValueType    arithmetic type of the return argument
     * @return      the value this Scalar represents as type ValueType
     */
    template<typename ValueType>
    inline ValueType getValue() const;

    // Removed: template<typename ValueType> operator ValueType () const  for type conversions
    // Might cause problems due to implicit conversions, should only be used explicitly
    // Now should be done as: cast<ValueType>( scalar )

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
     * @brief Conversion of a C type into value of enum ScalarType.
     *
     * @tparam ValueType    C++ type that should be converted
     * @return      value of enum type ScalarType that represents the C++ type.
     */
    template<typename ValueType>
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

private    :

    ComplexLongDouble mValue; //!< use highest precision for presentation
};

/** Output of ScalarType in stream is supported and very useful.
 */

COMMON_DLL_IMPORTEXPORT std::ostream& operator<<( std::ostream& stream, const Scalar::ScalarType& object );

/** @brief Cast operator to convert a Scalar into corresponding basic type.
 *
 *  This solutions avoids an implicit conversion of a Scalar to a basic type.
 */

template<typename ValueType>
ValueType cast( const Scalar& scalar )
{
    return scalar.getValue<ValueType>();
}

const Scalar zero;

inline Scalar::Scalar()
                : mValue( 0.0, 0.0 )
{
}

template<typename ValueType>
inline Scalar::Scalar( const ValueType value )
                : mValue( value, 0.0 )
{
}

inline Scalar::Scalar( const float value )
                : mValue( value, 0.0 )
{
}

inline Scalar::Scalar( const float real, const float imag )
                : mValue( real, imag )
{
}

inline Scalar::Scalar( const double value )
                : mValue( value, 0.0 )
{
}

inline Scalar::Scalar( const double real, const double imag )
                : mValue( real, imag )
{
}

inline Scalar::Scalar( const LongDouble value )
                : mValue( value, 0.0 )
{
}

inline Scalar::Scalar( const LongDouble real, const LongDouble imag )
                : mValue( real, imag )
{
}

template<typename ValueType>
inline Scalar::Scalar( const Complex<ValueType> value )
                : mValue( value.real(), value.imag() )
{
}

inline Scalar::~Scalar()
{
}

template<typename ValueType>
inline ValueType Scalar::getValue() const
{
    return static_cast<ValueType>( mValue.real() );
}

template<>
inline ComplexFloat Scalar::getValue() const
{
    return ComplexFloat( static_cast<float>( mValue.real() ), static_cast<float>( mValue.imag() ) );
}

template<>
inline ComplexDouble Scalar::getValue() const
{
    return ComplexDouble( static_cast<double>( mValue.real() ), static_cast<double>( mValue.imag() ) );
}

template<>
inline ComplexLongDouble Scalar::getValue() const
{
    return ComplexLongDouble( mValue.real(), mValue.imag() );
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
    if( isReal() )
    {
        stream << "Scalar(" << mValue.real() << ")";
    }
    else
    {
        stream << "Scalar(" << mValue.real() << "," << mValue.imag() << ")";
    }
}

template<typename ValueType>
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
inline Scalar::ScalarType Scalar::getType<LongDouble>()
{
    return LONG_DOUBLE;
}

template<>
inline Scalar::ScalarType Scalar::getType<ComplexFloat>()
{
    return COMPLEX;
}

template<>
inline Scalar::ScalarType Scalar::getType<ComplexDouble>()
{
    return DOUBLE_COMPLEX;
}

template<>
inline Scalar::ScalarType Scalar::getType<ComplexLongDouble>()
{
    return LONG_DOUBLE_COMPLEX;
}

inline size_t Scalar::getTypeSize( const ScalarType type )
{
    size_t typeSize = 0;

    switch( type )
    {
        case FLOAT:
            typeSize = sizeof(float);
            break;

        case DOUBLE:
            typeSize = sizeof(double);
            break;

        case LONG_DOUBLE:
            typeSize = sizeof(LongDouble);
            break;

        case COMPLEX:
            typeSize = sizeof(ComplexFloat);
            break;

        case DOUBLE_COMPLEX:
            typeSize = sizeof(ComplexDouble);
            break;

        case LONG_DOUBLE_COMPLEX:
            typeSize = sizeof(ComplexLongDouble);
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
    return Scalar( a.getValue<ComplexLongDouble>() + b.getValue<ComplexLongDouble>() );
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
    return Scalar( a.getValue<ComplexLongDouble>() - b.getValue<ComplexLongDouble>() );
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
    return Scalar( a.getValue<ComplexLongDouble>() * b.getValue<ComplexLongDouble>() );
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
    return Scalar( a.getValue<ComplexLongDouble>() / b.getValue<ComplexLongDouble>() );
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
    return a.getValue<ComplexLongDouble>() == b.getValue<ComplexLongDouble>();
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
    return a.getValue<ComplexLongDouble>() != b.getValue<ComplexLongDouble>();
}

inline bool operator<( const Scalar& a, const Scalar& b )
{
    return a.getValue<ComplexLongDouble>() < b.getValue<ComplexLongDouble>();
}

inline bool operator>( const Scalar& a, const Scalar& b )
{
    return a.getValue<ComplexLongDouble>() > b.getValue<ComplexLongDouble>();
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
    // note: uses sqrt for Complex<ValueType> with ValueType == long double
    return Scalar( sqrt( scalar.getValue<ComplexLongDouble>() ) );
}

inline Scalar abs( const Scalar scalar )
{
    // note: uses abs for Complex<ValueType> with ValueType == long double
    return Scalar( abs( scalar.getValue<ComplexLongDouble>() ) );
}

inline Scalar max( const Scalar a, const Scalar b )
{
    // note: must use std::max, otherwise infinite recursion
    return Scalar( std::max( a.getValue<ComplexLongDouble>(), b.getValue<ComplexLongDouble>() ) );
}

inline Scalar min( const Scalar a, const Scalar b )
{
    // note: must use std::min, otherwise infinite recursion
    return Scalar( std::min( a.getValue<ComplexLongDouble>(), b.getValue<ComplexLongDouble>() ) );
}

} //namespace lama

#endif // LAMA_SCALAR_HPP_
