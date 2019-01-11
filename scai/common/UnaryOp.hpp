/**
 * @file UnaryOp.hpp
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
 * @brief Enum typ for the different elementwise functions.
 * @author Lauretta Schubert
 * @date 05.10.2016
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/Math.hpp>

// std
#include <iostream>

namespace scai
{

namespace common
{

/** Enumeration type for unary operators used in elemental array operation
 *
 *  The unary operator specifies the function to be applied for each (array) element
 *
 *  \code
 *  A[i] = unary( A[i] )
 *  \endcode
 *
 */
enum class UnaryOp
{
    COPY,       //!< just identity
    CONJ,       //!< for conjugate of a vector
    ABS,        //!< for absolute value
    ASUM,       //!< sum of magnitudes of real and imaginary part, same as abs for real
    SQR,        //!< square the value
    MINUS,      //!< for negative value
    EXP,        //!< call exp on each vector element
    SQRT,       //!< call sqare root on each vector element
    SIN,        //!< call sin on each vector element
    COS,        //!< trigonometric function cos for each vector element
    TAN,        //!< trigonometric function tan on each vector element
    ATAN,       //!< call atan on each vector element
    LOG,        //!< call log on each vector element
    FLOOR,      //!< rounds downward
    CEIL,       //!< rounds upward
    SIGN,       //!< signum function
    RECIPROCAL, //!< builds multiplicate inverse, x^-1


    MAX_UNARY_OP //!< internal use only
};

/** This method provides a general routine for applying an unary operator.
 *
 *  @param[in] op specifies the operator
 *  @param[in] x  is the input value
 *  @returns   op( x )
 *
 *  Due to inlining this method can also be used within loops where op
 *  is loop invariant without losing performance. But therefore good
 *  compiler optimization must be switched on.
 */
template <typename ValueType>
CUDA_CALLABLE_MEMBER
inline ValueType applyUnary( const UnaryOp op, const ValueType& x )
{
    switch ( op )
    {
        case UnaryOp::COPY:
            return x;
        case UnaryOp::CONJ:
            return common::Math::conj( x );
        case UnaryOp::ABS:
            return common::Math::abs( x );
        case UnaryOp::ASUM:
            return common::Math::abs( common::Math::real( x ) ) + common::Math::abs( common::Math::imag( x ));
        case UnaryOp::SQR:
            return x * x;
        case UnaryOp::MINUS:
            return -x;
        case UnaryOp::EXP:
            return common::Math::exp( x );
        case UnaryOp::SQRT:
            return common::Math::sqrt( x );
        case UnaryOp::SIN:
            return common::Math::sin( x );
        case UnaryOp::COS:
            return common::Math::cos( x );
        case UnaryOp::TAN:
            return common::Math::tan( x );
        case UnaryOp::ATAN:
            return common::Math::atan( x );
        case UnaryOp::LOG:
            return common::Math::log( x );
        case UnaryOp::FLOOR:
            return common::Math::floor( x );
        case UnaryOp::CEIL:
            return common::Math::ceil( x );
        case UnaryOp::SIGN:
            return common::Math::sign( x );
        case UnaryOp::RECIPROCAL:
            return ValueType( 1 ) / x;
        default:
            return ValueType( 0 );
    }
}

template <>
inline IndexType applyUnary( const UnaryOp op, const IndexType& x )
{
    switch ( op )
    {
        case UnaryOp::COPY:
            return x;
        case UnaryOp::CONJ:
            return x;
        case UnaryOp::ABS:
            return common::Math::abs( x );
        case UnaryOp::ASUM:
            return common::Math::abs( x );
        case UnaryOp::MINUS:
            return -x;
        case UnaryOp::SQR:
            return x * x;
        default:
            return IndexType( 0 );
    }
}

template <typename ValueType>
inline bool isUnarySupported( const UnaryOp op )
{
    return op < UnaryOp::MAX_UNARY_OP;
}

template <>
inline bool isUnarySupported<IndexType>( const UnaryOp op )
{
    return op <= UnaryOp::MINUS;
}

/*
 * Output of UnaryOp in stream by writing strings instead of numbers
 */

inline std::ostream& operator<<( std::ostream& stream, const UnaryOp& op )
{
    switch ( op )
    {
        case UnaryOp::COPY:
            stream << "COPY";
            break;

        case UnaryOp::CONJ:
            stream << "CONJ";
            break;

        case UnaryOp::ABS:
            stream << "ABS";
            break;

        case UnaryOp::ASUM:
            stream << "ASUM";
            break;

        case UnaryOp::SQR:
            stream << "SQR";
            break;

        case UnaryOp::MINUS:
            stream << "MINUS";
            break;

        case UnaryOp::EXP:
            stream << "EXP";
            break;

        case UnaryOp::SQRT:
            stream << "SQRT";
            break;

        case UnaryOp::SIN:
            stream << "SIN";
            break;

        case UnaryOp::COS:
            stream << "COS";
            break;

        case UnaryOp::TAN:
            stream << "TAN";
            break;

        case UnaryOp::ATAN:
            stream << "ATAN";
            break;

        case UnaryOp::LOG:
            stream << "LOG";
            break;

        case UnaryOp::FLOOR:
            stream << "FLOOR";
            break;

        case UnaryOp::CEIL:
            stream << "CEIL";
            break;

        case UnaryOp::SIGN:
            stream << "SIGN";
            break;

        case UnaryOp::RECIPROCAL:
            stream << "RECIPROCAL";
            break;

        case UnaryOp::MAX_UNARY_OP:
            stream << "MAX_UNARY_OP for tests";
            break;

        default:
            stream << "<unknown_unary_op>";
            break;
    }

    return stream;
}

} /* end namespace common */

} /* end namespace scai */
