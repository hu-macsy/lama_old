/**
 * @file SolutionProxy.cpp
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
 * @brief SolutionProxy.cpp
 * @author Jiri Kraus
 * @date 07.06.2011
 * @since 1.0.0
 */

// hpp
#include <scai/solver/SolutionProxy.hpp>

namespace scai
{

namespace solver
{

SolutionProxy::SolutionProxy()
    : mIsDirty( true )
{
}

SolutionProxy::SolutionProxy( lama::Vector* const solution )
    : mSolution( solution ), mIsDirty( true )
{
}

SolutionProxy::~SolutionProxy()
{
}

const lama::Vector& SolutionProxy::getConstReference() const
{
    return ( *mSolution );
}

lama::Vector& SolutionProxy::operator*()
{
    return getReference();
}

void SolutionProxy::operator=( lama::Vector* const newVector )
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

lama::Vector& SolutionProxy::getReference()
{
    setDirty( true );
    return *mSolution;
}

lama::Vector* SolutionProxy::create()
{
    return lama::Vector::create( mSolution->getCreateValue() );
}

void SolutionProxy::swap( lama::Vector*& other )
{
    setDirty( true );
    std::swap( other, mSolution );
}

} /* end namespace solver */

} /* end namespace scai */
