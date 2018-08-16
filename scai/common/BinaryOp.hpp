/**
 * @file BinaryOp.hpp
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
 * @brief Enumeration type for the different binary operators
 * @author Thomas Brandes
 * @date 15.01.2016
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/macros/throw.hpp>

// std
#include <iostream>

namespace scai
{

namespace common
{

/** Enumeration type for binary operators used in set/scatter ops
 *
 *  The binary operator specifies in different kernel routines what kind
 *  of operation is applied to combine two elements, especially old and new values.
 *
 *  The associative binary operations (MIN, MAX, ADD, MULT) can also be used in
 *  redcution operations.
 */
enum class BinaryOp
{
    COPY,         //!< for assign   just take y
    ADD,          //!< for operator x + y
    SUB,          //!< for operator x - y
    MULT,         //!< for operator x * y
    DIVIDE,       //!< for operator x / y
    MODULO,       //!< for operator x % y
    MIN,          //!< for operator min( x, y )
    MAX,          //!< for operator max( x, y )
    ABS_MAX,      //!< for operator max( abs(x), abs(y) )
    ABS_DIFF,     //!< for operator abs( x - y )
    POW,          //!< for operator pow( x, y )
    COPY_SIGN,    //!< for operator magnitude(x) * sign(y)
    MAX_BINARY_OP //!< for internal use only
};

/** Method that applies a binary operation for two operands.
 *
 *  @param[in] x1, x2 the operands
 *  @param[in] op     specifies the binary operator
 *  @return    x1 op x2
 *
 *  Note: keep in mind that there is no exception handling in device code
 */
template <typename ValueType>
CUDA_CALLABLE_MEMBER
inline ValueType applyBinary( const ValueType& x1, const BinaryOp op, const ValueType& x2 )
{
    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    switch ( op )
    {
        case BinaryOp::COPY:
            return x2;
        case BinaryOp::ADD:
            return x1 + x2;
        case BinaryOp::SUB:
            return x1 - x2;
        case BinaryOp::MULT:
            return x1 * x2;
        case BinaryOp::DIVIDE:
            return x1 / x2;
        case BinaryOp::MODULO:
            return Math::mod( x1, x2 );
        case BinaryOp::COPY_SIGN:
            return Math::copysign( x1, x2 );
        case BinaryOp::POW:
            return Math::pow( x1, x2 );
        case BinaryOp::MIN:
            return Math::min( RealType( x1 ), RealType( x2 ) );
        case BinaryOp::MAX:
            return Math::max( RealType( x1 ), RealType( x2 ) );
        case BinaryOp::ABS_MAX:
            return Math::max( Math::abs( x1 ), Math::abs( x2 ) );
        case BinaryOp::ABS_DIFF:
            return Math::abs( x1 - x2 );
        default:
            return ValueType( 0 );
    }
}

// Template specialization for IndexType required, as some math operations, e.g pow, copysign are not defined

template <>
inline IndexType applyBinary( const IndexType& x1, const BinaryOp op, const IndexType& x2 )
{
    switch ( op )
    {
        case BinaryOp::COPY:
            return x2;
        case BinaryOp::ADD:
            return x1 + x2;
        case BinaryOp::SUB:
            return x1 - x2;
        case BinaryOp::MULT:
            return x1 * x2;
        case BinaryOp::DIVIDE:
            return x1 / x2;
        case BinaryOp::MODULO:
            return x1 % x2;
        case BinaryOp::MIN:
            return Math::min( x1, x2 );
        case BinaryOp::MAX:
            return Math::max( x1, x2 );
        case BinaryOp::ABS_MAX:
            return Math::max( Math::abs( x1 ), Math::abs( x2 ) );
        case BinaryOp::ABS_DIFF:
            // abs( x1 - x2 ) does not work for unsigned int
            if ( x1 < x2 ) 
            {
                return x2 - x1;
            }
            else
            {
                return x1 - x2;
            }
        default:
            return IndexType( 0 );
    }
}

/**
 * Predicate to check whether a binary op is supported for a certain value type.
 *
 * @tparam   ValueType specifies the type for which check is done
 * @param    op        is the binary operation for which the query is made
 * @return   true      if binary op is supported for the ValueType.
 *
 * Note: This predicate is helpful as applyBinary has no error handling.
 */
template <typename ValueType>
inline bool isBinarySupported( const BinaryOp op )
{
    return op < BinaryOp::MAX_BINARY_OP;
}

template <>
inline bool isBinarySupported<IndexType>( const BinaryOp op )
{
    return op <= BinaryOp::ABS_MAX;
}

/*
 * Output of BinaryOp in stream by writing strings instead of numbers
 */

inline std::ostream& operator<<( std::ostream& stream, const BinaryOp& op )
{
    switch ( op )
    {
        case BinaryOp::COPY:
            stream << "COPY";
            break;

        case BinaryOp::ADD:
            stream << "ADD";
            break;

        case BinaryOp::SUB:
            stream << "SUB";
            break;

        case BinaryOp::MULT:
            stream << "MULT";
            break;

        case BinaryOp::DIVIDE:
            stream << "DIVIDE";
            break;

        case BinaryOp::MODULO:
            stream << "MODULO";
            break;

        case BinaryOp::MIN:
            stream << "MIN";
            break;

        case BinaryOp::MAX:
            stream << "MAX";
            break;

        case BinaryOp::ABS_MAX:
            stream << "ABS_MAX";
            break;

        case BinaryOp::ABS_DIFF:
            stream << "ABS_DIFF";
            break;

        case BinaryOp::POW:
            stream << "POW";
            break;

        case BinaryOp::COPY_SIGN:
            stream << "COPY_SIGN";
            break;

        default:
            stream << "<unknown_binary_op>";
            break;
    }

    return stream;
}

/** This method returns for a binary operation reduced in a reduction the corresponding zero element.
 *
 *  @tparam ValueType specifies the type for which the operation is used
 *  @param[in] op must be a binary reduction operation
 *  @return the zero elemen for the reduction
 */

template <typename ValueType>
inline ValueType zeroBinary( const BinaryOp op )
{
    switch ( op )
    {
        case BinaryOp::COPY:
        case BinaryOp::ADD:
        case BinaryOp::SUB:
        case BinaryOp::ABS_MAX:
            return ValueType( 0 );
        case BinaryOp::MULT:
        case BinaryOp::DIVIDE:
            return ValueType( 1 );
        case BinaryOp::MIN:
            return TypeTraits<ValueType>::getMax();
        case BinaryOp::MAX:
            return TypeTraits<ValueType>::getMin();
        default:
        {
            COMMON_THROWEXCEPTION( "Illegal reduction operator: " << op )
        }
    }
}


} /* end namespace common */

} /* end namespace scai */
