/**
 * @file SolutionProxy.cpp
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
 * @brief SolutionProxy.cpp
 * @author Jiri Kraus
 * @date 07.06.2011
 * @since 1.0.0
 */

// hpp
#include <lama/solver/SolutionProxy.hpp>

namespace lama
{

SolutionProxy::SolutionProxy()
    : mIsDirty( true )
{
}

SolutionProxy::SolutionProxy( Vector* const solution )
    : mSolution( solution ), mIsDirty( true )
{
}

SolutionProxy::~SolutionProxy()
{
}

const Vector& SolutionProxy::getConstReference() const
{
    return ( *mSolution );
}

Vector& SolutionProxy::operator*()
{
    return getReference();
}

void SolutionProxy::operator=( Vector* const newVector )
{
    setDirty( true );
    mSolution = newVector;
}

bool SolutionProxy::isDirty() const
{
    return mIsDirty;
}

void SolutionProxy::setDirty( bool isDirty )
{
    mIsDirty = isDirty;
}

Vector& SolutionProxy::getReference()
{
    setDirty( true );
    return ( *mSolution );
}

std::auto_ptr<Vector> SolutionProxy::create()
{
    return std::auto_ptr<Vector>( mSolution->create() );
}

void SolutionProxy::swap( Vector*& other )
{
    setDirty( true );
    std::swap( other, mSolution );
}

} // namespace lama
