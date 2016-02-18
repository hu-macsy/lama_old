/**
 * @file Solver.cpp
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
 * @brief Solver.cpp
 * @author Jiri Kraus
 * @date 08.06.2011
 * @since 1.0.0
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
// TODO mContext
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
    if( getConstRuntime().mInitialized )
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

    if( getConstRuntime().mSolveInit )
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

    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumRows(), runtime.mRhs->size(), "mismatch: #rows of matrix, rhs" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumColumns(), solution.size(), "mismatch: #cols of matrix, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getColDistribution(), solution.getDistribution(), "mismatch: matrix col dist, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getDistribution(), runtime.mRhs->getDistribution(), "mismatch: matrix row dist, rhs dist" )

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

    SCAI_LOG_DEBUG( logger, "getResidual of solver " << mId << ", is dirty = " << runtime.mSolution.isDirty() 
                            << ", runtime.mResidual.get() = " << runtime.mResidual.get() )

    SCAI_ASSERT_DEBUG( runtime.mCoefficients, "mCoefficients == NULL" )
    SCAI_ASSERT_DEBUG( runtime.mRhs, "mRhs == NULL" )

    //mLogger->logMessage(LogLevel::completeInformation,"Request for residual received.\n");

    if( runtime.mSolution.isDirty() || !runtime.mResidual.get() )
    {
        SCAI_REGION( "Solver.computeResidual" )

        SCAI_LOG_DEBUG( logger, "calculating residual of = " << runtime.mSolution.getConstReference() )

        if( !runtime.mResidual.get() )
        {
            runtime.mResidual.reset( Vector::create( runtime.mRhs->getCreateValue() ) );
        }

        //mLogger->logMessage(LogLevel::completeInformation,"Residual needs revaluation.\n");

        mLogger->startTimer( "ResidualTimer" );

        *runtime.mResidual = ( *runtime.mRhs ) - ( *runtime.mCoefficients ) * ( runtime.mSolution.getConstReference() );

        mLogger->stopTimer( "ResidualTimer" );
        mLogger->logTime( "ResidualTimer", LogLevel::completeInformation, "Revaluation of residual took [s]: " );
        mLogger->stopAndResetTimer( "ResidualTimer" );

        runtime.mSolution.setDirty( false );
    }

    return ( *runtime.mResidual );
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

void Solver::setContextPtr( hmemo::ContextPtr context )
{
    SCAI_LOG_DEBUG( logger, "Set context to " << *context )
    mContext = context;
}

void Solver::writeAt( std::ostream& stream ) const
{
    stream << "Solver ( id = " << mId << " )";
}

} /* end namespace solver */

} /* end namespace scai */
