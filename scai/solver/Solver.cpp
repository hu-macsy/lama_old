/**
 * @file Solver.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Implementation of methods for the base class Solver.
 * @author Jiri Kraus
 * @date 08.06.2011
 */

// hpp
#include <scai/solver/Solver.hpp>

// local library
#include <scai/solver/logger/CommonLogger.hpp>
#include <scai/solver/logger/Timer.hpp>

#include <scai/lama/expression/MatrixVectorExpressions.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( Solver::logger, "Solver" )

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

Solver::Solver( const std::string& id )
    : mId( id ), mLogger(
        new CommonLogger( "dummyLog", LogLevel::noLogging,
                          LoggerWriteBehaviour::toConsoleOnly,
                          common::shared_ptr<Timer>( new Timer() ) ) )
{
    SCAI_LOG_INFO( Solver::logger, "Solver id = " << mId << " created, dummy log" )
}

Solver::Solver( const std::string& id, LoggerPtr logger )
    : mId( id ), mLogger( logger )
{
    SCAI_LOG_INFO( Solver::logger, "Solver id = " << mId << " created, with logger" )
}

Solver::Solver( const Solver& other )
    : mId( other.mId ), mLogger( other.mLogger )
{
}

Solver::SolverRuntime::SolverRuntime()
    : mCoefficients( 0 ), mRhs( 0 ), mResidual(), mInitialized( false ), mSolveInit( false )
{
}

Solver::~Solver()
{
    // mRhs, mCoefficents are used as 'dynamic' references, no free
    SCAI_LOG_INFO( Solver::logger, "~Solver " << mId )
}

Solver::SolverRuntime::~SolverRuntime()
{
    SCAI_LOG_INFO( logger, "~SolverRuntime" )
}

void Solver::initialize( const Matrix& coefficients )
{
    if ( getConstRuntime().mInitialized )
    {
        SCAI_LOG_DEBUG( logger, "Previous initialization of solver found! Will be overridden!" )
        mLogger->logMessage( LogLevel::solverInformation, "Solver already initialized, will be overridden\n" );
    }

    getRuntime().mCoefficients = &coefficients;
    getRuntime().mInitialized = true;
    mLogger->logMessage( LogLevel::solverInformation, "Solver initialized\n" );
}

void Solver::solve( Vector& solution, const Vector& rhs )
{
    SCAI_REGION( "Solver.solve" )
    SCAI_ASSERT( getConstRuntime().mInitialized, "Solver not initialized, solve cannot be called" )

    if ( getConstRuntime().mSolveInit )
    {
        SCAI_LOG_DEBUG( logger, "Previous initialization of 'solve'-process found! Will be overridden!" )
    }

    solveInit( solution, rhs );
    solveImpl();
    solveFinalize();
}

void Solver::solveInit( Vector& solution, const Vector& rhs )
{
    SolverRuntime& runtime = getRuntime();
    runtime.mRhs = &rhs;
    runtime.mSolution = &solution;
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumRows(), rhs.size(), "mismatch: #rows of matrix, rhs" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumColumns(), solution.size(), "mismatch: #cols of matrix, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getColDistribution(), solution.getDistribution(), "mismatch: matrix col dist, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getRowDistribution(), rhs.getDistribution(), "mismatch: matrix row dist, rhs dist" )
    runtime.mSolveInit = true;
}

void Solver::solveFinalize()
{
}

const std::string& Solver::getId() const
{
    SCAI_LOG_TRACE( logger, "Returning Solver Id " << mId )
    return mId;
}

const Vector& Solver::getResidual() const
{
    const SolverRuntime& runtime = getConstRuntime();

    if ( runtime.mResidual.get() )
    {
        SCAI_LOG_DEBUG( logger, "getResidual of solver " << mId << ", is dirty = " << runtime.mSolution.isDirty()
                        << ", runtime.mResidual = " << *runtime.mResidual )
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "getResidual of solver " << mId << ", residual not available yet" )
    }

    // initialize and solveInit must have been called before

    SCAI_ASSERT_DEBUG( runtime.mCoefficients, "mCoefficients == NULL" )
    SCAI_ASSERT_DEBUG( runtime.mRhs, "mRhs == NULL" )

    //mLogger->logMessage(LogLevel::completeInformation,"Request for residual received.\n");

    if ( runtime.mSolution.isDirty() || !runtime.mResidual.get() )
    {
        SCAI_REGION( "Solver.computeResidual" )
        SCAI_LOG_DEBUG( logger, "calculating residual of = " << runtime.mSolution.getConstReference() )

        if ( !runtime.mResidual.get() )
        {
            // VERY IMPORTANT: newVector makes sure that residual has same context
            //                 otherwise: many unnecessary data movements !!!

            runtime.mResidual.reset( runtime.mRhs->newVector() );
        }

        //mLogger->logMessage(LogLevel::completeInformation,"Residual needs revaluation.\n");
        mLogger->startTimer( "ResidualTimer" );
        *runtime.mResidual = ( *runtime.mRhs ) - ( *runtime.mCoefficients ) * ( runtime.mSolution.getConstReference() );
        mLogger->stopTimer( "ResidualTimer" );
        mLogger->logTime( "ResidualTimer", LogLevel::completeInformation, "Revaluation of residual took [s]: " );
        mLogger->stopAndResetTimer( "ResidualTimer" );
        runtime.mSolution.setDirty( false );
    }

    return *runtime.mResidual;
}

const Matrix& Solver::getCoefficients() const
{
    SCAI_ASSERT_DEBUG( getConstRuntime().mCoefficients, "mCoefficents == NULL" )
    return *getConstRuntime().mCoefficients;
}

void Solver::setLogger( LoggerPtr logger )
{
    mLogger = logger;
}

void Solver::setLogLevel( LogLevel::LogLevel level )
{
    mLogger->setLogLevel( level );
}

void Solver::writeAt( std::ostream& stream ) const
{
    stream << "Solver ( id = " << mId << " )";
}

} /* end namespace solver */

} /* end namespace scai */
