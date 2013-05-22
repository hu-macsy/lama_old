/**
 * @file TrivialPreconditioner.cpp
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
 * @brief TrivialPreconditioner.cpp
 * @author Matthias Makulla
 * @date 06.04.2011
 * $Id$
 */

// hpp
#include <lama/solver/TrivialPreconditioner.hpp>

namespace lama
{

TrivialPreconditioner::TrivialPreconditioner( const std::string& id )
    : Solver( id )
{
}

TrivialPreconditioner::TrivialPreconditioner( const std::string& id, LoggerPtr logger )
    : Solver( id, logger )
{
}

TrivialPreconditioner::TrivialPreconditioner( const TrivialPreconditioner& other )
    : Solver( other )
{
}

TrivialPreconditioner::TrivialPreconditionerRuntime::TrivialPreconditionerRuntime()
    : SolverRuntime()
{
}

TrivialPreconditioner::~TrivialPreconditioner()
{
}

TrivialPreconditioner::TrivialPreconditionerRuntime::~TrivialPreconditionerRuntime()
{
}

void TrivialPreconditioner::initialize( const Matrix& coefficients )
{
    Solver::initialize( coefficients );
}

void TrivialPreconditioner::solveImpl()
{
    *( getRuntime().mSolution ) = *( getRuntime().mRhs );
}

SolverPtr TrivialPreconditioner::copy()
{
    return SolverPtr( new TrivialPreconditioner( *this ) );
}

TrivialPreconditioner::TrivialPreconditionerRuntime& TrivialPreconditioner::getRuntime()
{
    return mTrivialPreconditionerRuntime;
}

const TrivialPreconditioner::TrivialPreconditionerRuntime& TrivialPreconditioner::getConstRuntime() const
{
    return mTrivialPreconditionerRuntime;
}

} // namespace lama
