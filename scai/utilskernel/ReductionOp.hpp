/**
 * @file ReductionOp.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Enum typ for the different reduction operator.s
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

/** Own struct for enum type of reduction operators */

struct reduction
{
    /** Enumeration type for reduction operators used in set/scatter ops  
     *
     *  The reduction operator specifies for typical reductions what kind
     *  of operator is applied to combine two elements.
     *
     *  \code
     *  A[i] = A[i] reduction_op B[i]
     *  \endcode
     *  
     */

    typedef enum
    {
        COPY,     // for assign   x = y
        ADD,      // for operator x += y
        SUB,      // for operator x -= y
        MULT,     // for operator x *= y
        DIVIDE,   // for operator x /= y
        MIN,      // for operator x = min( x, y )
        MAX,      // for operator x = max( x, y )
        ABS_MAX   // for operator x = max( x, abs(y) )
    } ReductionOp;

};

/*
 * Output of ReductionOp in stream by writing strings instead of numbers
 */

inline std::ostream& operator<<( std::ostream& stream, const reduction::ReductionOp& op )
{
    switch ( op )
    {
        case reduction::COPY:
            stream << "COPY";
            break;
        case reduction::ADD:
            stream << "ADD";
            break;
        case reduction::SUB:
            stream << "SUB";
            break;
        case reduction::MULT:
            stream << "MULT";
            break;
        case reduction::DIVIDE:
            stream << "DIVIDE";
            break;
        case reduction::MIN:
            stream << "MIN";
            break;
        case reduction::MAX:
            stream << "MAX";
            break;
        case reduction::ABS_MAX:
            stream << "ABS_MAX";
            break;
        default:
            stream << "<unknown_reduction_op>";
            break;
    }
    return stream;
}

} /* end namespace utilskernel */

} /* end namespace scai */
