/**
 * @file ElementwiseOp.hpp
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

struct elementwise
{
    /** Enumeration type for elementwise operators used in set/scatter ops
     *
     *  The elementwise operator specifies for typical elementwises what kind
     *  of operator is applied to combine two elements.
     *
     *  \code
     *  A[i] = A[i] elementwise_op B[i]
     *  \endcode
     *
     */

    typedef enum
    {
        INVERT,  // for inverse/reciprocal of a vector
        CONJ,    // for conjugate of a vector
        EXP,
        SQRT,    // call sqrt on each vector element
        SIN,     // call sin on each vector element
        COS,     // call cos on each vector element
        TAN,     // call tan on each vector element
        ATAN,    // call atan on each vector element
        LOG,     // call log on each vector element

        MAX_ELEMENTWISE_OP // only for tests, leave this at the end
    } ElementwiseOp;

};

/*
 * Output of ElementwiseOp in stream by writing strings instead of numbers
 */

inline std::ostream& operator<<( std::ostream& stream, const elementwise::ElementwiseOp& op )
{
    switch ( op )
    {
        case elementwise::INVERT:
            stream << "INVERT";
            break;

        case elementwise::CONJ:
            stream << "CONJ";
            break;

        case elementwise::EXP:
            stream << "EXP";
            break;

        case elementwise::SQRT:
            stream << "SQRT";
            break;

        case elementwise::SIN:
            stream << "SIN";
            break;

        case elementwise::COS:
            stream << "COS";
            break;

        case elementwise::TAN:
            stream << "TAN";
            break;

        case elementwise::ATAN:
            stream << "ATAN";
            break;

        case elementwise::LOG:
            stream << "LOG";
            break;

        case elementwise::MAX_ELEMENTWISE_OP:
            stream << "MAX_ELEMENTWISE_OP-only for tests";
            break;

        default:
            stream << "<unknown_elementwise_op>";
            break;
    }

    return stream;
}

} /* end namespace utilskernel */

} /* end namespace scai */
