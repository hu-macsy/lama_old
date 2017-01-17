/**
 * @file BinaryOp.hpp
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

namespace utilskernel
{

/** Own struct for enum type of binary operators */

struct binary
{
    /** Enumeration type for binary operators used in set/scatter ops
     *
     *  The binary operator specifies in different kernel routines what kind
     *  of operation is applied to combine two elements, especially old and new values.
     *
     *  The associative binary operations (MIN, MAX, ADD, MULT) can also be used in
     *  redcution operations.
     */
    typedef enum
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
        POW,          //!< for operator pow( x, y )
        COPY_SIGN,    //!< for operator magnitude(x) * sign(y)
        MAX_BINARY_OP //!< for internal use only
    } BinaryOp;
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
MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER
inline ValueType applyBinary( const ValueType& x1, const binary::BinaryOp op, const ValueType& x2 )
{
    switch ( op )
    {
        case binary::COPY:
            return x2;
        case binary::ADD:
            return x1 + x2;
        case binary::SUB:
            return x1 - x2;
        case binary::MULT:
            return x1 * x2;
        case binary::DIVIDE:
            return x1 / x2;
        case binary::MODULO:
            return common::Math::mod( x1, x2 );
        case binary::COPY_SIGN:
            return common::Math::copysign( x1, x2 );
        case binary::POW:
            return common::Math::pow( x1, x2 );
        case binary::MIN:
            return common::Math::min( x1, x2 );
        case binary::MAX:
            return common::Math::max( x1, x2 );
        case binary::ABS_MAX:
            return common::Math::max( common::Math::abs( x1 ), common::Math::abs( x2 ) );
        default:
            return ValueType( 0 );
    }
}

// Template specialization for IndexType required, as some math operations, e.g pow, copysign are not defined

template <>
inline IndexType applyBinary( const IndexType& x1, const binary::BinaryOp op, const IndexType& x2 )
{
    switch ( op )
    {
        case binary::COPY:
            return x2;
        case binary::ADD:
            return x1 + x2;
        case binary::SUB:
            return x1 - x2;
        case binary::MULT:
            return x1 * x2;
        case binary::DIVIDE:
            return x1 / x2;
        case binary::MODULO:
            return x1 % x2;
        case binary::MIN:
            return common::Math::min( x1, x2 );
        case binary::MAX:
            return common::Math::max( x1, x2 );
        case binary::ABS_MAX:
            return common::Math::max( common::Math::abs( x1 ), common::Math::abs( x2 ) );
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
inline bool isBinarySupported( const binary::BinaryOp op )
{
    return op < binary::MAX_BINARY_OP;
}

template <>
inline bool isBinarySupported<IndexType>( const binary::BinaryOp op )
{
    return op <= binary::ABS_MAX;
}

/** This method returns for a binary operation reduced in a reduction the corresponding zero element.
 *
 *  @tparam ValueType specifies the type for which the operation is used
 *  @param[in] op must be a binary reduction operation
 *  @return the zero elemen for the reduction
 */

template <typename ValueType>
inline ValueType zeroBinary( const binary::BinaryOp op )
{
    switch ( op )
    {
        case binary::ADD:
        case binary::ABS_MAX:
            return ValueType( 0 );
        case binary::MIN:
            return common::TypeTraits<ValueType>::getMax();
        case binary::MAX:
            return common::TypeTraits<ValueType>::getMin();
        default:
        {
            COMMON_THROWEXCEPTION( "Illegal reduction operator: " << op )
        }
    }
}

/*
 * Output of BinaryOp in stream by writing strings instead of numbers
 */

inline std::ostream& operator<<( std::ostream& stream, const binary::BinaryOp& op )
{
    switch ( op )
    {
        case binary::COPY:
            stream << "COPY";
            break;

        case binary::ADD:
            stream << "ADD";
            break;

        case binary::SUB:
            stream << "SUB";
            break;

        case binary::MULT:
            stream << "MULT";
            break;

        case binary::DIVIDE:
            stream << "DIVIDE";
            break;

        case binary::MODULO:
            stream << "MODULO";
            break;

        case binary::MIN:
            stream << "MIN";
            break;

        case binary::MAX:
            stream << "MAX";
            break;

        case binary::ABS_MAX:
            stream << "ABS_MAX";
            break;

        case binary::POW:
            stream << "POW";
            break;

        case binary::COPY_SIGN:
            stream << "COPY_SIGN";
            break;

        default:
            stream << "<unknown_binary_op>";
            break;
    }

    return stream;
}

} /* end namespace utilskernel */

} /* end namespace scai */
