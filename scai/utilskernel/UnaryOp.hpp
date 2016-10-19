/**
 * @file UnaryOp.hpp
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
 * @brief Enum typ for the different elementwise functions.
 * @author Lauretta Schubert
 * @date 05.10.2016
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

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

/*
 * Output of UnaryOp in stream by writing strings instead of numbers
 */

inline std::ostream& operator<<( std::ostream& stream, const unary::UnaryOp& op )
{
    switch ( op )
    {
        case unary::CONJ:
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
