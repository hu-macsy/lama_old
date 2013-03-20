/**
 * @file Solver.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * $Id$
 */

// hpp
#include <lama/solver/Solver.hpp>

// others
#include <lama/solver/logger/CommonLogger.hpp>
#include <lama/solver/logger/OpenMPTimer.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>

// assert
#include <lama/tracing.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( Solver::logger, "Solver" );

Solver::Solver( const std::string& id )
    : mId( id ), mLogger(
        new CommonLogger( "dummyLog", lama::LogLevel::noLogging,
                          lama::LoggerWriteBehaviour::toConsoleOnly,
                          std::auto_ptr<Timer>( new OpenMPTimer() ) ) )
{
    LAMA_LOG_INFO( Solver::logger, "Solver id = " << mId << " created, dummy log" );
}

Solver::Solver( const std::string& id, LoggerPtr logger )
    : mId( id ), mLogger( logger )
{
    LAMA_LOG_INFO( Solver::logger, "Solver id = " << mId << " created, with logger" );
}

Solver::Solver( const Solver& other )
    : mId( other.mId ), mLogger( other.mLogger )
// TODO mContext
{
}

Solver::SolverRuntime::SolverRuntime()
    : mCoefficients( 0 ), mRhs( 0 ), mResidual( 0 ), mInitialized( false ), mSolveInit( false )
{
}

Solver::~Solver()
{
    // mRhs, mCoefficents are used as 'dynamic' references, no free

    LAMA_LOG_INFO( Solver::logger, "~Solver " << mId );
}

Solver::SolverRuntime::~SolverRuntime()
{
}

void Solver::initialize( const Matrix& coefficients )
{
    if ( getConstRuntime().mInitialized )
    {
        LAMA_LOG_DEBUG( logger, "Previous initialization of solver found! Will be overridden!" );
        mLogger->logMessage( LogLevel::solverInformation, "Solver already initialized, will be overridden\n" );
    }
    getRuntime().mCoefficients = &coefficients;
    getRuntime().mInitialized = true;
    mLogger->logMessage( LogLevel::solverInformation, "Solver initialized\n" );
}

void Solver::solve( Vector& solution, const Vector& rhs )
{
    LAMA_REGION( "Solver.solve" );
    if ( getConstRuntime().mSolveInit )
    {
        LAMA_LOG_DEBUG( logger, "Previous initialization of 'solve'-process found! Will be overridden!" );
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

    if ( runtime.mCoefficients->getNumRows() != runtime.mRhs->size() )
    {
        LAMA_THROWEXCEPTION(
            "Size of rhs vector " << *runtime.mRhs << " does not match column size of matrix " << *runtime.mCoefficients );
    }

    if ( runtime.mCoefficients->getNumColumns() != solution.size() )
    {
        LAMA_THROWEXCEPTION(
            "Size of solution vector " << solution << " does not match row size of matrix " << *runtime.mCoefficients );
    }

    if ( runtime.mCoefficients->getColDistribution() != solution.getDistribution() )
    {
        LAMA_THROWEXCEPTION(
            "Distribution of lhs " << solution << " = " << solution.getDistribution() << " does not match (row) distribution of " << *runtime.mCoefficients << " = " << runtime.mCoefficients->getColDistribution() );
    }

    if ( runtime.mCoefficients->getDistribution() != runtime.mRhs->getDistribution() )
    {
        LAMA_THROWEXCEPTION(
            "Distribution of old Solution " << *runtime.mRhs << " = " << runtime.mRhs->getDistribution() << " does not match (row) distribution of " << *runtime.mCoefficients << " = " << runtime.mCoefficients->getDistribution() );
    }

    runtime.mSolveInit = true;
}

void Solver::solveFinalize()
{
}

const std::string& Solver::getId() const
{
    LAMA_LOG_TRACE( logger, "Returning Solver Id " << mId );
    return mId;
}

const Vector& Solver::getResidual() const
{
    const SolverRuntime& runtime = getConstRuntime();
    LAMA_ASSERT_DEBUG( runtime.mCoefficients, "mCoefficients == NULL" );
    LAMA_ASSERT_DEBUG( runtime.mRhs, "mRhs == NULL" );

    //mLogger->logMessage(LogLevel::completeInformation,"Request for residual received.\n");

    if ( runtime.mSolution.isDirty() || !runtime.mResidual.get() )
    {
        LAMA_LOG_DEBUG( logger, "calculating residual of = " << &(runtime.mSolution.getConstReference()) );
        if ( !runtime.mResidual.get() )
        {
            runtime.mResidual = runtime.mRhs->create();
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
    LAMA_ASSERT_DEBUG( getConstRuntime().mCoefficients, "mCoefficents == NULL" );

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

void Solver::setContext( ContextPtr context )
{
    LAMA_LOG_DEBUG( logger, "Set context to " << *context );
    mContext = context;
}

void Solver::writeAt( std::ostream& stream ) const
{
    stream << "Solver ( id = " << mId << " )";
}

} //namespace lama
