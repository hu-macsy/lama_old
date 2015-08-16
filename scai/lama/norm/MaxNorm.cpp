/**
 * @file MaxNorm.cpp
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
 * @brief MaxNorm.cpp
 * @author Jiri Kraus
 * @date 14.06.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/norm/MaxNorm.hpp>

namespace scai
{

namespace lama
{

MaxNorm::MaxNorm()
{
}

MaxNorm::~MaxNorm()
{
}

Scalar MaxNorm::apply( const Scalar& scalar ) const
{
    return maxNorm( scalar );
}

Scalar MaxNorm::apply( const Vector& vector ) const
{
    return maxNorm( vector );
}

Scalar MaxNorm::apply( const Matrix& matrix ) const
{
    return maxNorm( matrix );
}

Scalar maxNorm( const Scalar& scalar )
{
    return abs( scalar );
}

Scalar maxNorm( const Vector& vector )
{
    return vector.maxNorm();
}

Scalar maxNorm( const Matrix& matrix )
{
    return matrix.maxNorm();
}

} /* end namespace lama */

} /* end namespace scai */
