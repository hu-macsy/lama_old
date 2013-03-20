/**
 * @file L2Norm.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief L2Norm.cpp
 * @author Jiri Kraus
 * @date 01.06.2011
 * $Id$
 */

// hpp
#include <lama/norm/L2Norm.hpp>

namespace lama
{

L2Norm::L2Norm()
{
}

L2Norm::~L2Norm()
{
}

Scalar L2Norm::apply( const Scalar& scalar ) const
{
    return l2Norm( scalar );
}

Scalar L2Norm::apply( const Vector& vector ) const
{
    return l2Norm( vector );
}

Scalar l2Norm( const Scalar& scalar )
{
    return abs( scalar );
}

Scalar l2Norm( const Vector& vector )
{
    return vector.l2Norm();
}

}
