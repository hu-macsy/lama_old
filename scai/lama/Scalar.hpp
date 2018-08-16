/**
 * @file Scalar.hpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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

namespace intern
{

/**
 * @brief The class Scalar represents a container for all kind of ValueType used for Matrix and/or Vector.
 *
 * An object of the class Scalar is needed in template expressions.
 *
 * For a Scalar the arithmetic operations +, -, *, / etc. are
 * also supported to allow a high flexibility. 
 *
 * ScalarRepType is used internally for the representation of
 * the value. For each supported arithmetic type SCAI_NUMERIC_TYPES the following
 * conversions must be supported:
 *
 *    - ScalarRepType( SCAI_NUMERIC_TYPES v )
 *    - SCAI_NUMERIC_TYPE( ScalarRepType v )
 *
 * Conversion into the representation type and back should be lossless, i. e. the
 * following relation must / should  hold:
 *
 * \code
 *    SCAI_NUMERIC_TYPE( ScalarRepType( x ) ) == x
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
    inline Scalar();

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
    inline Scalar( const ValueType value );

    /**
     *  Copy constructor
     */
    inline Scalar( const Scalar& x );

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

    inline virtual void writeAt( std::ostream& stream ) const;

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    ScalarRepType mValue;  //!< use highest precision for representation
};

/* --------------------------------------------------------------------------- *
 *  Implementation of methods for Scalar                                       *
 * --------------------------------------------------------------------------- */

inline Scalar::Scalar() : mValue( 0 )
{
}

template<typename ValueType>
inline Scalar::Scalar( const ValueType value ) : mValue( value )
{
}

inline Scalar::Scalar( const Scalar& x ) : mValue( x.mValue )
{
}

inline Scalar::~Scalar()
{
}

template<typename ValueType>
inline ValueType Scalar::getValue() const
{
    return static_cast<ValueType>( mValue );
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

inline Scalar operator-( const Scalar& a )
{
    return Scalar( -a.getValue<ScalarRepType>() );
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

} /* end namespace intern */

} /* end namespace lama */

} /* end namespace scai */
