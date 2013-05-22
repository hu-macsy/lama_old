/**
 * @file IterativeSolver.cpp
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
 * @brief IterativeSolver.cpp
 * @author Kai Buschulte
 * @date 19.07.2011
 * @since 1.0.0
 */

// hpp
#include <lama/solver/IterativeSolver.hpp>

// others
#include <lama/norm/L2Norm.hpp>

#include <lama/solver/criteria/IterationCount.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( IterativeSolver::logger, "Solver.IterativeSolver" )

IterativeSolver::IterativeSolver( const std::string& id )
    : lama::IterativeSolver::Solver( id ), mCriterionRootComponent( new IterationCount( 1 ) )
{
}

IterativeSolver::IterativeSolver( const std::string& id, LoggerPtr logger )
    : lama::IterativeSolver::Solver( id, logger ), mCriterionRootComponent( new IterationCount( 1 ) )
{
}

IterativeSolver::IterativeSolver( const IterativeSolver& other )
    : Solver( other )
{
    if ( mCriterionRootComponent )
    {
        mCriterionRootComponent.reset( new Criterion( *other.mCriterionRootComponent ) );
    }
    else
    {
        LAMA_LOG_INFO( logger, other <<" has no conditions, that can be copied." )
        mCriterionRootComponent.reset( new IterationCount( 1 ) );
    }

    if ( other.getPreconditioner() )
    {
        mPreconditioner = other.getPreconditioner()->copy();
    }
    else
    {
        LAMA_LOG_INFO( logger, other << " has no preconditioner." )
    }
}

IterativeSolver::IterativeSolverRuntime::IterativeSolverRuntime()
    : SolverRuntime(), mIterations( 0 )
{
}

IterativeSolver::~IterativeSolver()
{
}

IterativeSolver::IterativeSolverRuntime::~IterativeSolverRuntime()
{
}

void IterativeSolver::initialize( const Matrix& coefficients )
{
    Solver::initialize( coefficients );

    if ( mPreconditioner )
    {
        if ( mPreconditioner->getConstRuntime().mInitialized )
        {
            LAMA_LOG_INFO( logger, "Preconditioner already initialized, skipping recursive init." )
            mLogger->logMessage( LogLevel::solverInformation,
                                 "Preconditioner already initialized, skipping recursive init\n" );
        }
        else
        {
            mPreconditioner->initialize( coefficients );
            mLogger->logMessage( LogLevel::solverInformation, "Preconditioner initialized\n" );
        }
    }
}

void IterativeSolver::solveImpl()
{
    getRuntime().mIterations = 0;

    if ( !getConstRuntime().mSolveInit )
    {
        LAMA_THROWEXCEPTION(
            "Solver " + this->getId() + " has not been initialized. Call solveInit( Vector& solution, const Vector& rhs ) before solving " + this->getId() )
    }

    logStartSolve();

    while ( !criteriaAreSatisfied() )
    {
        logIterationStart();

        iterate();
        getRuntime().mIterations++;

        logIterationEndAndResidual();
    }

    logEndSolve();
}

void IterativeSolver::setStoppingCriterion( const CriterionPtr criterion )
{
    LAMA_ASSERT_ERROR( criterion, "Criterion defined is NULL." )

    LAMA_LOG_INFO( logger, "Criteria " << *criterion << " defined." )

    mCriterionRootComponent = criterion;
}

bool IterativeSolver::criteriaAreSatisfied() const
{
    if ( mCriterionRootComponent.get() == 0 )
    {
        LAMA_THROWEXCEPTION( this->getId() + ": No stopping criterion set." )
    }

    return mCriterionRootComponent->isSatisfied( *this );
}

void IterativeSolver::setPreconditioner( SolverPtr const conditioner )
{
    LAMA_LOG_INFO( logger, "Preconditioner " << conditioner->getId() << " defined." )
    mPreconditioner = conditioner;
}

const SolverPtr IterativeSolver::getPreconditioner() const
{
    return mPreconditioner;
}

int IterativeSolver::getIterationCount() const
{
    return getConstRuntime().mIterations;
}

void IterativeSolver::logStartSolve()
{
    mLogger->logNewLine( LogLevel::solverInformation );
    mLogger->logNewLine( LogLevel::solverInformation );
    mLogger->logMessage( LogLevel::solverInformation, "Beginning solve.\n" );
    L2Norm l2Norm;
    mLogger->logResidual( LogLevel::convergenceHistory, *this, l2Norm, "Start " );
    mLogger->startTimer( "SolutionTimer" );
}

void IterativeSolver::logIterationEndAndResidual()
{
    IterativeSolverRuntime& runtime = getRuntime();
    mLogger->stopTimer( "IterationTimer" );
    std::stringstream iterationPrefix;
    iterationPrefix << "Iteration #" << runtime.mIterations << " Duration [s]: ";
    mLogger->logTime( "IterationTimer", LogLevel::advancedInformation, iterationPrefix.str() );
    mLogger->stopAndResetTimer( "IterationTimer" );
    L2Norm l2Norm;
    iterationPrefix.str( "" );
    iterationPrefix << "Iteration #" << runtime.mIterations << " ";
    mLogger->logResidual( LogLevel::convergenceHistory, *this, l2Norm, iterationPrefix.str() );
}

void IterativeSolver::logEndSolve()
{
    mLogger->logNewLine( LogLevel::solverInformation );
    mLogger->logMessage( LogLevel::solverInformation, "-------------------------\n" );
    mLogger->logType( LogLevel::solverInformation, "Solve completed. Iterations: ", getRuntime().mIterations );
    mLogger->logTime( "SolutionTimer", LogLevel::solverInformation, "Total Runtime [s]: " );
    mLogger->stopAndResetTimer( "SolutionTimer" );
    mLogger->logNewLine( LogLevel::solverInformation );

}

void IterativeSolver::logIterationStart()
{
    mLogger->logNewLine( LogLevel::completeInformation );
    mLogger->logType( LogLevel::completeInformation, "Beginning iteration #", getRuntime().mIterations + 1 );
    mLogger->startTimer( "IterationTimer" );
}

} // namespace lama
