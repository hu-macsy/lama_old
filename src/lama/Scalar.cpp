/**
 * @file Scalar.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Scalar.cpp
 * @author Jiri Kraus
 * @date 07.11.2011
 * @since 1.0.0
 */

// hpp
#include <lama/Scalar.hpp>

namespace lama
{

Scalar& Scalar::operator+=( Scalar& other )
{
    mValue += other.mValue;
    return *this;
}

Scalar& Scalar::operator-=( Scalar& other )
{
    mValue -= other.mValue;
    return *this;
}

Scalar& Scalar::operator*=( Scalar& other )
{
    mValue *= other.mValue;
    return *this;
}

Scalar& Scalar::operator/=( Scalar& other )
{
    mValue /= other.mValue;
    return *this;
}

std::ostream& operator<<( std::ostream& stream, const Scalar::ScalarType& object )
{
    switch ( object )
    {
    case Scalar::FLOAT:
        stream << "float";
        break;
    case Scalar::DOUBLE:
        stream << "double";
        break;
    case Scalar::LONG_DOUBLE:
        stream << "long double";
        break;
    case Scalar::COMPLEX:
        stream << "complex<float>";
        break;
    case Scalar::DOUBLE_COMPLEX:
        stream << "complex<double>";
        break;
    case Scalar::LONG_DOUBLE_COMPLEX:
        stream << "complex<long double>";
        break;
    default:
        stream << "Unknown";
    }
    return stream;
}

} //namespace lama
