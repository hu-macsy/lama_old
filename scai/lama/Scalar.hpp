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
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>
#include <scai/common/TypeTraits.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/preprocessor.hpp>

// std
#include <cstdio>

namespace scai
{

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
 *
 * ScalarRepType is used internally for the representation of
 * the value. For each supported arithmetic type ARITHMETIC_TYPE the following
 * conversions must be supported:
 *
 *    - ScalarRepType( ARITHMETIC_TYPE v )
 *    - ARITHMETIC_TYPE( ScalarRepType v )
 *
 * Conversion into the representation type and back should be lossless, i. e. the
 * following relation must / should  hold:
 *
 * \code
 *    ARITHEMTIC_TYPE( ScalarRepType( x ) ) == x  
 * \endcode
 */
class COMMON_DLL_IMPORTEXPORT Scalar: public common::Printable
{
public:

    /**
     * @brief ExpressionMemberType is the type that is used the template Expression to store a Scalar.
     */
    typedef const Scalar ExpressionMemberType;

    /**
     * @brief Constructs a scalar representing 0.
     */
    inline Scalar() : mValue( 0 )
    {
    }

    /**
     * @brief Constructs a scalar representing the passed real value.
     *
     * The templated conversion constructor needs to be explicit,
     * because the operator==(Scalar,Scalar) can lead to ambiguities.
     *
     * @tparam ValueType  type of the input argument value for constructor of Scalar
     * @param[in] value   the value this scalar should represent
     */
    template<typename ValueType>
    explicit inline Scalar( const ValueType value ) : mValue( value )
    {
    }

    /**
     * @brief Constructor of scalar for each supported arithmetic type.
     */

#define LAMA_SCALAR_METHODS(z, I, _ )                                        \
    inline Scalar( const ARITHMETIC_HOST_TYPE_##I value ) : mValue( value )  \
    {                                                                        \
    }                                                                        \

    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_SCALAR_METHODS, _ )

#undef LAMA_SCALAR_METHODS

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
     * @brief Overload assignment operator +=
     */
    Scalar& operator+=( Scalar& other )
    {
        mValue += other.mValue;
        return *this;
    }

    /**
     * @brief Overload assignment operator -=
     */
    Scalar& operator-=( Scalar& other )
    {
        mValue -= other.mValue;
        return *this;
    }

    /**
     * @brief Overload assignment operator *=
     */
    Scalar& operator*=( Scalar& other )
    {
        mValue *= other.mValue;
        return *this;
    }

    /**
     * @brief Overload assignment operator /=
     */
    Scalar& operator/=( Scalar& other )
    {
        mValue /= other.mValue;
        return *this;
    }

    /**
     *  @brief Query that scalar values has no imaginary part.
     */

    inline virtual void writeAt( std::ostream& stream ) const;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    ScalarRepType mValue;  //!< use highest precision for representation
};

/* --------------------------------------------------------------------------- *
 *  Implementation of methods for Scalar                                       *
 * --------------------------------------------------------------------------- */

inline Scalar::~Scalar()
{
}

template<typename ValueType>
inline ValueType Scalar::getValue() const
{
    return static_cast<ValueType>( mValue );
}

inline Scalar Scalar::operator-() const
{
    return Scalar( -mValue );
}

inline void Scalar::writeAt( std::ostream& stream ) const
{
    stream << "Scalar(" << mValue << ")";
}

/* --------------------------------------------------------------------------- *
 *  Binaray operators for Scalar: +, -, +, /                                   *
 * --------------------------------------------------------------------------- */

/* Instead of definining the binary operators as class methods we use this
 * global syntax to support also:
 *
 *  Scalar a = 1.0;
 *  Scalar b = 2 * a;   // becomes Scalar( 2) * Scalar( a )
 */

/**
 * @brief Add Scalar a with Scalar b
 *
 * @param[in] a     1st Scalar.
 * @param[in] b     2nd Scalar.
 * @return          sum
 */
inline Scalar operator+( const Scalar& a, const Scalar& b )
{
    return Scalar( a.getValue<ScalarRepType>() + b.getValue<ScalarRepType>() );
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
    return Scalar( a.getValue<ScalarRepType>() - b.getValue<ScalarRepType>() );
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
    return Scalar( a.getValue<ScalarRepType>() * b.getValue<ScalarRepType>() );
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
    return Scalar( a.getValue<ScalarRepType>() / b.getValue<ScalarRepType>() );
}

/* --------------------------------------------------------------------------- *
 *  Compare operators for Scalar: ==, !=, <, >, <=, >=                         *
 * --------------------------------------------------------------------------- */

/**
 * @brief Check equality of a and b.
 *
 * @param[in] a     the 1st Scalar to compare this to.
 * @param[in] b     the 2nd Scalar to compare this to.
 * @return          if a is equal to b
 */

inline bool operator==( const Scalar& a, const Scalar& b )
{
    return a.getValue<ScalarRepType>() == b.getValue<ScalarRepType>();
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
    return a.getValue<ScalarRepType>() != b.getValue<ScalarRepType>();
}

inline bool operator<( const Scalar& a, const Scalar& b )
{
    return a.getValue<ScalarRepType>() < b.getValue<ScalarRepType>();
}

inline bool operator>( const Scalar& a, const Scalar& b )
{
    return a.getValue<ScalarRepType>() > b.getValue<ScalarRepType>();
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
    // call sqrt for ScalarRepType

    return Scalar( common::TypeTraits<ScalarRepType>::sqrt( scalar.getValue<ScalarRepType>() ) );
}

/* --------------------------------------------------------------------------- *
 *  Unary functions for Scalar: abs, conj                                      *
 * --------------------------------------------------------------------------- */

inline Scalar abs( const Scalar scalar )
{
    // call abs for ScalarRepType

    return Scalar( common::TypeTraits<ScalarRepType>::abs( scalar.getValue<ScalarRepType>() ) );
}

inline Scalar conj( const Scalar scalar )
{
    return Scalar( common::TypeTraits<ScalarRepType>::conj( scalar.getValue<ScalarRepType>() ) );
}

/* --------------------------------------------------------------------------- *
 *  Binary functions for Scalar: min, max                                      *
 * --------------------------------------------------------------------------- */

inline Scalar max( const Scalar a, const Scalar b )
{
    return Scalar(   a.getValue<ScalarRepType>() >= b.getValue<ScalarRepType>() 
                   ? a.getValue<ScalarRepType>() 
                   : b.getValue<ScalarRepType>() );
}

inline Scalar min( const Scalar a, const Scalar b )
{
    return Scalar(   a.getValue<ScalarRepType>() <= b.getValue<ScalarRepType>() 
                   ? a.getValue<ScalarRepType>() 
                   : b.getValue<ScalarRepType>() );
}

} /* end namespace lama */

} /* end namespace scai */
