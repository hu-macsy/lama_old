/**
 * @file UnaryOp.hpp
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

namespace utilskernel
{

/** Own struct for enum type of elementwise functions */

struct unary
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

    typedef enum
    {
        CONJ,    //!< for conjugate of a vector
        ABS,     //!< for absolute value
        MINUS,   //!< for negative value
        EXP,     //!< call exp on each vector element
        SQRT,    //!< call sqrt on each vector element
        SIN,     //!< call sin on each vector element
        COS,     //!< trigonometric function cos for each vector element
        TAN,     //!< trigonometric function tan on each vector element
        ATAN,    //!< call atan on each vector element
        LOG,     //!< call log on each vector element
        FLOOR,   //!< rounds downward
        CEIL,    //!< rounds upward


        MAX_UNARY_OP //!< internal use only

    } UnaryOp;
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
MIC_CALLABLE_MEMBER CUDA_CALLABLE_MEMBER
inline ValueType applyUnary( const unary::UnaryOp op, const ValueType& x )
{
    switch ( op )
    {
        case unary::CONJ:   return common::Math::conj( x );
        case unary::ABS:    return common::Math::abs( x );
        case unary::MINUS:  return -x;
        case unary::EXP:    return common::Math::exp( x );
        case unary::SQRT:   return common::Math::sqrt( x );
        case unary::SIN:    return common::Math::sin( x );
        case unary::COS:    return common::Math::cos( x );
        case unary::TAN:    return common::Math::tan( x );
        case unary::ATAN:   return common::Math::atan( x );
        case unary::LOG:    return common::Math::log( x );
        case unary::FLOOR:  return common::Math::floor( x );
        case unary::CEIL:   return common::Math::ceil( x );
        default:            return ValueType( 0 );
    }
}

template <>
inline IndexType applyUnary( const unary::UnaryOp op, const IndexType& x )
{
    switch ( op )
    {
        case unary::CONJ:   return x;
        case unary::ABS:    return common::Math::abs( x );
        case unary::MINUS:  return -x;
        default:            return IndexType( 0 );
    }
}

template <typename ValueType>
inline bool isUnarySupported( const unary::UnaryOp op )
{
    return op < unary::MAX_UNARY_OP;
}

template <>
inline bool isUnarySupported<IndexType>( const unary::UnaryOp op )
{
    return op <= unary::MINUS;
}

/*
 * Output of UnaryOp in stream by writing strings instead of numbers
 */

inline std::ostream& operator<<( std::ostream& stream, const unary::UnaryOp& op )
{
    switch ( op )
    {
        case unary::unary::CONJ:
            stream << "CONJ";
            break;

        case unary::ABS:
            stream << "ABS";
            break;

        case unary::MINUS:
            stream << "MINUS";
            break;

        case unary::EXP:
            stream << "EXP";
            break;

        case unary::SQRT:
            stream << "SQRT";
            break;

        case unary::SIN:
            stream << "SIN";
            break;

        case unary::COS:
            stream << "COS";
            break;

        case unary::TAN:
            stream << "TAN";
            break;

        case unary::ATAN:
            stream << "ATAN";
            break;

        case unary::LOG:
            stream << "LOG";
            break;

        case unary::FLOOR:
            stream << "FLOOR";
            break;

        case unary::CEIL:
            stream << "CEIL";
            break;

        case unary::MAX_UNARY_OP:
            stream << "MAX_UNARY_OP for tests";
            break;

        default:
            stream << "<unknown_unary_op>";
            break;
    }

    return stream;
}

} /* end namespace utilskernel */

} /* end namespace scai */
