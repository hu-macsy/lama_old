/**
 * @file OmegaSolver.hpp
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
 * @brief OmegaSolver.hpp
 * @author Kai Buschulte
 * @date 10.08.2011
 * @since 1.0.0
 */

// hpp
#include <scai/solver/OmegaSolver.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( OmegaSolver::logger, "Solver.IterativeSolver.OmegaSolver" )

using lama::Matrix;
using lama::Vector;

OmegaSolver::OmegaSolver( const std::string& id )
    : IterativeSolver( id ), mOmega( 0.5 )
{
}

OmegaSolver::OmegaSolver( const std::string& id, const lama::Scalar omega )
    : IterativeSolver( id ), mOmega( omega )
{
}

OmegaSolver::OmegaSolver( const std::string& id, LoggerPtr logger )
    : IterativeSolver( id, logger ), mOmega( 0.5 )
{
}

OmegaSolver::OmegaSolver( const std::string& id, const lama::Scalar omega, LoggerPtr logger )
    : IterativeSolver( id, logger ), mOmega( omega )
{
}

OmegaSolver::OmegaSolver( const OmegaSolver& other )
    : IterativeSolver( other ), mOmega( other.mOmega )
{
}

OmegaSolver::OmegaSolverRuntime::OmegaSolverRuntime()
    : IterativeSolverRuntime()
{
}

OmegaSolver::~OmegaSolver()
{
}

OmegaSolver::OmegaSolverRuntime::~OmegaSolverRuntime()
{
}

void OmegaSolver::initialize( const Matrix& coefficients )
{
    IterativeSolver::initialize( coefficients );
}

void OmegaSolver::setOmega( const lama::Scalar omega )
{
    mOmega = omega;
}

lama::Scalar OmegaSolver::getOmega() const
{
    return mOmega;
}

} /* end namespace solver */

} /* end namespace scai */
