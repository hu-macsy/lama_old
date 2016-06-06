/**
 * @file IterativeSolver.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief IterativeSolver.cpp
 * @author Kai Buschulte
 * @date 19.07.2011
 */

// hpp
#include <scai/solver/IterativeSolver.hpp>

// local library
#include <scai/solver/criteria/IterationCount.hpp>

// internal scai libraries
#include <scai/lama/norm/L2Norm.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( IterativeSolver::logger, "Solver.IterativeSolver" )

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

IterativeSolver::IterativeSolver( const std::string& id )
    : Solver( id ), mCriterionRootComponent( new IterationCount( 1 ) )
{
}

IterativeSolver::IterativeSolver( const std::string& id, LoggerPtr logger )
    : Solver( id, logger ), mCriterionRootComponent( new IterationCount( 1 ) )
{
}

IterativeSolver::IterativeSolver( const IterativeSolver& other )
    : Solver( other )
{
    if( mCriterionRootComponent )
    {
        mCriterionRootComponent.reset( new Criterion( *other.mCriterionRootComponent ) );
    }
    else
    {
        SCAI_LOG_INFO( logger, other <<" has no conditions, that can be copied." )
        mCriterionRootComponent.reset( new IterationCount( 1 ) );
    }

    if( other.getPreconditioner() )
    {
        mPreconditioner = other.getPreconditioner()->copy();
    }
    else
    {
        SCAI_LOG_INFO( logger, other << " has no preconditioner." )
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

    if( mPreconditioner )
    {
        if( mPreconditioner->getConstRuntime().mInitialized )
        {
            SCAI_LOG_INFO( logger, "Preconditioner already initialized, skipping recursive init." )
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

    if( !getConstRuntime().mSolveInit )
    {
        COMMON_THROWEXCEPTION(
            "Solver " + this->getId()
            + " has not been initialized. Call solveInit( Vector& solution, const Vector& rhs ) before solving "
            + this->getId() )
    }

    logStartSolve();

    while( !criteriaAreSatisfied() )
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
    SCAI_ASSERT_ERROR( criterion, "Criterion defined is NULL." )

    SCAI_LOG_INFO( logger, "Criteria " << *criterion << " defined." )

    mCriterionRootComponent = criterion;
}

bool IterativeSolver::criteriaAreSatisfied() const
{
    if( mCriterionRootComponent.get() == 0 )
    {
        COMMON_THROWEXCEPTION( this->getId() + ": No stopping criterion set." )
    }

    return mCriterionRootComponent->isSatisfied( *this );
}

void IterativeSolver::setPreconditioner( SolverPtr const conditioner )
{
    SCAI_LOG_INFO( logger, "Preconditioner " << conditioner->getId() << " defined." )
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
    lama::L2Norm l2Norm;
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
    lama::L2Norm l2Norm;
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

IterativeSolver* IterativeSolver::create( const std::string type, const std::string name )
{
    IterativeSolver* sov = dynamic_cast<IterativeSolver*>( Solver::create( type, name ) );

    if( !sov )
    {
        COMMON_THROWEXCEPTION( "requested Solver is not inherited from IterativeSolver" )
    }

    return sov;
}

void IterativeSolver::writeAt( std::ostream& stream ) const
{
    stream << "IterativeSolver ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

} /* end namespace solver */

} /* end namespace scai */
