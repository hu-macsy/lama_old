/**
 * @file ContextType.cpp
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
 * @author Jiri Kraus
 * @date 07.11.2011
 */

// hpp
#include <scai/common/ContextType.hpp>

namespace scai
{

namespace common
{

namespace context
{

std::ostream& operator<<( std::ostream& stream, const ContextType& type )
{
    switch ( type )
    {
        case Host :
            stream << "Host";
            break;

        case CUDA :
            stream << "CUDA";
            break;

        case MIC :
            stream << "MIC";
            break;

        case OpenCL :
            stream << "OpenCL";
            break;

        case UserContext :
            stream << "UserContext";
            break;

        default:
            stream << "ContextType_" << (int) type;
    }

    return stream;
}

/* -----------------------------------------------------------------------------*/

std::ostream& operator<<( std::ostream& stream, const AccessKind& kind )
{
    switch ( kind )
    {
        case Write :
            stream << "Write";
            break;

        case Read :
            stream << "Read";
            break;

        default:
            stream << "AccessKind_" << (int) kind;
    }

    return stream;
}

} /* end namespace context */

} /* end namespace common */

} /* end namespace scai */
