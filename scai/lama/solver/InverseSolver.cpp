/**
 * @file InverseSolver.cpp
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
 * @brief InverseSolver.cpp
 * @author Jiri Kraus
 * @date 08.06.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/solver/InverseSolver.hpp>

// local library
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/LAMAInterface.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

// std
#include <sstream>

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( InverseSolver::logger, "Solver.InverseSolver" )

InverseSolver::InverseSolver( const std::string & id )
    : Solver( id )
{
    SCAI_LOG_INFO( InverseSolver::logger, "InverseSolver, id = " << id )
}

InverseSolver::InverseSolver( const std::string & id, LoggerPtr logger )
    : Solver( id, logger )
{
    SCAI_LOG_INFO( InverseSolver::logger, "InverseSolver, id = " << id )
}

InverseSolver::InverseSolver( const InverseSolver& other )
    : Solver( other )
{
    SCAI_LOG_INFO( InverseSolver::logger, "InverseSolver, id = " << other.mId )
}

InverseSolver::InverseSolverRuntime::InverseSolverRuntime()
    : SolverRuntime()
{
}

InverseSolver::~InverseSolver()
{
    SCAI_LOG_INFO( logger, "~InverseSolver" )
}

InverseSolver::InverseSolverRuntime::~InverseSolverRuntime()
{
}

/* --------------------------------------------------------------------------- */

void InverseSolver::initialize( const Matrix& coefficients )
{
    SCAI_REGION( "Solver.Inverse.intialize" )

    SCAI_LOG_INFO( logger, "Initializing with " << coefficients )

    getRuntime().mInverse = MatrixPtr( coefficients.clone() );

    getRuntime().mInverse->invert( coefficients );

    getRuntime().mInverse->setContext( coefficients.getContextPtr() );

    getRuntime().mInverse->prefetch();

    Solver::initialize( coefficients );
}

/* --------------------------------------------------------------------------- */

const Matrix& InverseSolver::getInverse() const
{
    SCAI_ASSERT_ERROR( getConstRuntime().mInverse, "inverse not available (no call of initialize before)" );

    return *getConstRuntime().mInverse;
}

/* --------------------------------------------------------------------------- */

void InverseSolver::solveImpl()
{
    SCAI_REGION( "Solver.Inverse.solve" )

    InverseSolverRuntime& runtime = getRuntime();

    SCAI_ASSERT_ERROR( runtime.mInverse.get(), "solve, but mInverse is NULL" )

    logStartSolve();
    *runtime.mSolution = ( *runtime.mInverse ) * ( *runtime.mRhs );
    logEndSolve();
}

/* --------------------------------------------------------------------------- */

void InverseSolver::setContext( ContextPtr context )
{
    Solver::setContext( context );

    if( getRuntime().mInverse )
    {
        getRuntime().mInverse->setContext( mContext );
    }
    else
    {
        SCAI_LOG_WARN( logger, "setContext on uninitialized solver" )
    }
}

/* --------------------------------------------------------------------------- */

void InverseSolver::logStartSolve()
{
    mLogger->startTimer( "SolutionTimer" );
}

/* --------------------------------------------------------------------------- */

void InverseSolver::logEndSolve()
{
    L2Norm l2Norm;
    mLogger->logResidual( LogLevel::convergenceHistory, *this, l2Norm, "Final " );
    mLogger->logTime( "SolutionTimer", LogLevel::solverInformation, "Total Runtime [s]: " );
    mLogger->stopAndResetTimer( "SolutionTimer" );
    mLogger->logNewLine( LogLevel::solverInformation );

}

/* --------------------------------------------------------------------------- */

InverseSolver::InverseSolverRuntime& InverseSolver::getRuntime()
{
    return mInverseSolverRuntime;
}

/* --------------------------------------------------------------------------- */

const InverseSolver::InverseSolverRuntime& InverseSolver::getConstRuntime() const
{
    return mInverseSolverRuntime;
}

/* --------------------------------------------------------------------------- */

SolverPtr InverseSolver::copy()
{
    return SolverPtr( new InverseSolver( *this ) );
}

} /* end namespace lama */

} /* end namespace scai */
