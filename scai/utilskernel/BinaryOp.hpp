/**
 * @file BinaryOp.hpp
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
 * @brief Enumeration type for the different binary operators 
 * @author Thomas Brandes
 * @date 15.01.2016
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

/** Own struct for enum type of binary operators */

struct binary
{
    /** Enumeration type for binary operators used in set/scatter ops
     *
     *  The binary operator specifies different kernel routines what kind
     *  of operator is applied to combine two elements.
     *
     *  \code
     *  A[i] = A[i] binary_op B[i]
     *  A[i] = A[i] binary_op val
     *  \endcode
     *
     *  The associative binary operations (MIN, MAX, ADD, MULT) can also be used in
     *  redcution operations.
     */

    typedef enum
    {
        COPY,         //!< for assign   x = y
        ADD,          //!< for operator x += y
        SUB,          //!< for operator x -= y
        MULT,         //!< for operator x *= y
        DIVIDE,       //!< for operator x /= y
        MIN,          //!< for operator x = min( x, y )
        MAX,          //!< for operator x = max( x, y )
        ABS_MAX,      //!< for operator x = max( x, abs(y) )
        POW,          //!< for operator x = pow( x, y )
        COPY_SIGN,    //!< for operator x = magnitude(x), sign(y)
        MAX_BINARY_OP //!< for internal use only
    } BinaryOp;

};

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
