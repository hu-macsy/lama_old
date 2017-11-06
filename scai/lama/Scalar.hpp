/**
 * @file Scalar.hpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief Definition of dummy class whose objects stand for arbitrary arithmetic types.
 * @author Jiri Kraus
 * @date 22.02.2011
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

// base classes
#include <scai/common/Printable.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/macros/loop.hpp>

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
 * the value. For each supported arithmetic type SCAI_NUMERIC_TYPES_TYPE the following
 * conversions must be supported:
 *
 *    - ScalarRepType( SCAI_NUMERIC_TYPES_TYPE v )
 *    - SCAI_NUMERIC_TYPES_TYPE( ScalarRepType v )
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

    inline Scalar( const Scalar& x ) : mValue( x.mValue )
    {
    }

    /**
     * @brief Constructor of scalar for each supported arithmetic type.
     */
#define SCAI_LAMA_SCALAR_CONSTRUCTORS( type )                \
    inline Scalar( const type value ) : mValue( value ) \
    { }

    SCAI_COMMON_LOOP( SCAI_LAMA_SCALAR_CONSTRUCTORS, SCAI_ALL_TYPES )

#undef SCAI_LAMA_SCALAR_CONSTRUCTORS

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
    Scalar& operator+=( const Scalar& other )
    {
        mValue += other.mValue;
        return *this;
    }

    /**
     * @brief Overload assignment operator -=
     */
    Scalar& operator-=( const Scalar& other )
    {
        mValue -= other.mValue;
        return *this;
    }

    /**
     * @brief Overload assignment operator *=
     */
    Scalar& operator*=( const Scalar& other )
    {
        mValue *= other.mValue;
        return *this;
    }

    /**
     * @brief Overload assignment operator /=
     */
    Scalar& operator/=( const Scalar& other )
    {
        mValue /= other.mValue;
        return *this;
    }

    /**
     *  @brief Query that scalar values has no imaginary part.
     */

    inline virtual void writeAt( std::ostream& stream ) const;

    inline bool hasComplexValue() const
    {
        return common::Math::imag( mValue ) != common::Constants::ZERO;
    }

    /** Return a Scalar with the corresponding eps0 value of a type.
     *
     *  @param[in] type is the enum value of the required type
     *
     *  @returns TypeTraits<ValueType>::eps0() for ValueType with TypeTraits<ValueType>::sid == type
     */

    static inline Scalar eps0( const common::ScalarType type );

    /** Return a Scalar with the corresponding eps1 value of a type.
     *
     *  @param[in] type is the enum value of the required type
     *
     *  @returns TypeTraits<ValueType>::eps1() for ValueType with TypeTraits<ValueType>::sid == type
     */

    static inline Scalar eps1( const common::ScalarType type );

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
    // comparison only defined for real types, not for complex types

    if ( common::isComplex( common::TypeTraits<ScalarRepType>::stype ) )
    {
        SCAI_ASSERT_EQ_ERROR( common::Math::imag( a.getValue<ScalarRepType>() ), 0, "complex value in a < b, a = " << a )
        SCAI_ASSERT_EQ_ERROR( common::Math::imag( b.getValue<ScalarRepType>() ), 0, "complex value in a < b, b = " << a )
    }

    typedef common::TypeTraits<ScalarRepType>::AbsType RealType;

    return a.getValue<RealType>() < b.getValue<RealType>();
}

inline bool operator<=( const Scalar& a, const Scalar& b )
{
    // comparison only defined for real types, not for complex types

    if ( common::isComplex( common::TypeTraits<ScalarRepType>::stype ) )
    {
        SCAI_ASSERT_EQ_ERROR( common::Math::imag( a.getValue<ScalarRepType>() ), 0, "complex value in a < b, a = " << a )
        SCAI_ASSERT_EQ_ERROR( common::Math::imag( b.getValue<ScalarRepType>() ), 0, "complex value in a < b, b = " << a )
    }

    typedef common::TypeTraits<ScalarRepType>::AbsType RealType;

    return a.getValue<RealType>() <= b.getValue<RealType>();
}

inline bool operator>( const Scalar& a, const Scalar& b )
{
    // comparison only defined for real types, not for complex types

    if ( common::isComplex( common::TypeTraits<ScalarRepType>::stype ) )
    {
        SCAI_ASSERT_EQ_ERROR( common::Math::imag( a.getValue<ScalarRepType>() ), 0, "complex value in a < b, a = " << a )
        SCAI_ASSERT_EQ_ERROR( common::Math::imag( b.getValue<ScalarRepType>() ), 0, "complex value in a < b, b = " << a )
    }

    typedef common::TypeTraits<ScalarRepType>::AbsType RealType;

    return a.getValue<RealType>() > b.getValue<RealType>();
}

inline bool operator>=( const Scalar& a, const Scalar& b )
{
    // comparison only defined for real types, not for complex types

    if ( common::isComplex( common::TypeTraits<ScalarRepType>::stype ) )
    {
        SCAI_ASSERT_EQ_ERROR( common::Math::imag( a.getValue<ScalarRepType>() ), 0, "complex value in a < b, a = " << a )
        SCAI_ASSERT_EQ_ERROR( common::Math::imag( b.getValue<ScalarRepType>() ), 0, "complex value in a < b, b = " << a )
    }

    typedef common::TypeTraits<ScalarRepType>::AbsType RealType;

    return a.getValue<RealType>() >= b.getValue<RealType>();
}

inline Scalar sqrt( const Scalar scalar )
{
    // call sqrt for ScalarRepType
    return Scalar( common::Math::sqrt( scalar.getValue<ScalarRepType>() ) );
}

/* --------------------------------------------------------------------------- *
 *  Unary functions for Scalar: abs, conj                                      *
 * --------------------------------------------------------------------------- */

inline Scalar abs( const Scalar scalar )
{
    // call abs for ScalarRepType
    return Scalar( common::Math::abs( scalar.getValue<ScalarRepType>() ) );
}

inline Scalar conj( const Scalar scalar )
{
    return Scalar( common::Math::conj( scalar.getValue<ScalarRepType>() ) );
}

/* --------------------------------------------------------------------------- *
 *  Binary functions for Scalar: min, max                                      *
 * --------------------------------------------------------------------------- */

inline Scalar max( const Scalar a, const Scalar b )
{
    SCAI_ASSERT_ERROR( !a.hasComplexValue(), "complex value in max: " << a )
    SCAI_ASSERT_ERROR( !b.hasComplexValue(), "complex value in max: " << b )

    typedef common::TypeTraits<ScalarRepType>::AbsType RealType;

    return Scalar( common::Math::max(
                       common::Math::real( a.getValue<RealType>() ),
                       common::Math::real( b.getValue<RealType>() ) ) );
}

inline Scalar min( const Scalar a, const Scalar b )
{
    SCAI_ASSERT_ERROR( !a.hasComplexValue(), "complex value in min: " << a )
    SCAI_ASSERT_ERROR( !b.hasComplexValue(), "complex value in min: " << b )

    typedef common::TypeTraits<ScalarRepType>::AbsType RealType;

    return Scalar( common::Math::min(
                       common::Math::real( a.getValue<RealType>() ),
                       common::Math::real( b.getValue<RealType>() ) ) );
}

template<typename TList>
struct TypeTraitAccess;

template<>
struct TypeTraitAccess<common::mepr::NullType>
{
    static Scalar eps1( const common::ScalarType& )
    {
        return Scalar( 0 );
    }

    static Scalar eps0( const common::ScalarType& )
    {
        return Scalar( 0 );
    }
};

template<typename H, typename T>
struct TypeTraitAccess<common::mepr::TypeList<H, T> >
{
    static Scalar eps1( const common::ScalarType& type )
    {
        if ( common::TypeTraits<H>::stype == type )
        {
            return Scalar( common::TypeTraits<H>::eps1() );
        }
        else
        {
            return TypeTraitAccess<T>::eps1( type );
        }
    }

    static Scalar eps0( const common::ScalarType& type )
    {
        if ( common::TypeTraits<H>::stype == type )
        {
            return Scalar( common::TypeTraits<H>::eps0() );
        }
        else
        {
            return TypeTraitAccess<T>::eps0( type );
        }
    }
};

Scalar Scalar::eps0( const common::ScalarType type )
{
    return TypeTraitAccess<SCAI_NUMERIC_TYPES_HOST_LIST>::eps0( type );
}

Scalar Scalar::eps1( const common::ScalarType type )
{
    return TypeTraitAccess<SCAI_NUMERIC_TYPES_HOST_LIST>::eps1( type );
}

} /* end namespace lama */

} /* end namespace scai */
