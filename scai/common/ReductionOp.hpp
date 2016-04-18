/**
 * @file ReductionOp.hpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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

namespace common
{

/** Own namespace for enum type of reduction operators */

namespace reduction
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

    /*
     * Output of ReductionOp in stream by writing strings instead of numbers
     */

    inline std::ostream& operator<<( std::ostream& stream, const ReductionOp& op )
    {
        switch ( op )
        {
            case COPY:
                stream << "COPY";
                break;
            case ADD:
                stream << "ADD";
                break;
            case SUB:
                stream << "SUB";
                break;
            case MULT:
                stream << "MULT";
                break;
            case DIVIDE:
                stream << "DIVIDE";
                break;
            case MIN:
                stream << "MIN";
                break;
            case MAX:
                stream << "MAX";
                break;
            case ABS_MAX:
                stream << "ABS_MAX";
                break;
            default:
                stream << "<unknown_reduction_op>";
                break;
        }
        return stream;
    }

} /* end namespace reduction */

} /* end namespace common */

} /* end namespace scai */
