/**
 * @file ScalarType.cpp
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
 * @brief Implementation of operations on ScalarType.
 * @author Thomas Brandes
 * @date 07.11.2011
 */

// hpp
#include <scai/common/ScalarType.hpp>

namespace scai
{

namespace common
{

MIC_CALLABLE_MEMBER const char* scalar2str( const scalar::ScalarType stype )
{
    switch ( stype )
    {
        case scalar::FLOAT:
            return "float";

        case scalar::DOUBLE:
            return "double";

        case scalar::INDEX_TYPE:
            return "IndexType";

        case scalar::LONG_DOUBLE:
            return "LongDouble";

        case scalar::COMPLEX:
            return "ComplexFloat";

        case scalar::DOUBLE_COMPLEX:
            return "ComplexDouble";

        case scalar::LONG_DOUBLE_COMPLEX:
            return "ComplexLongDouble";

        case scalar::INTERNAL:
            return "_Internal";

        case scalar::PATTERN:
            return "_Pattern";

        default:
            return "Unknown";
    }
}

namespace scalar
{

bool isComplex( const ScalarType t )
{
    bool is = false;

    switch ( t )
    {
        case DOUBLE_COMPLEX:
        case COMPLEX:
        case LONG_DOUBLE_COMPLEX :
            is = true;
            break;

        default:
            is = false;
    }

    return is;
}

std::ostream& operator<<( std::ostream& stream, const ScalarType& object )
{
    stream << scalar2str( object );
    return stream;
}

} /* end namespace scalar */

} /* end namespace common */

} /* end namespace scai */
