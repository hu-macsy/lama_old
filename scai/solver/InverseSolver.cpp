/**
 * @file InverseSolver.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @endlicense
 *
 * @brief InverseSolver.cpp
 * @author Jiri Kraus
 * @date 08.06.2011
 */

// hpp
#include <scai/solver/InverseSolver.hpp>

// local library
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/expression/MatrixVectorExpressions.hpp>

// internal scai libraries
#include <scai/tracing.hpp>

// std
#include <sstream>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( InverseSolver::logger, "Solver.InverseSolver" )

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

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

    getRuntime().mInverse = lama::MatrixPtr( coefficients.newMatrix() );

    getRuntime().mInverse->invert( coefficients );

    getRuntime().mInverse->setContextPtr( coefficients.getContextPtr() );

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

void InverseSolver::setContextPtr( hmemo::ContextPtr context )
{
    Solver::setContextPtr( context );

    if( getRuntime().mInverse )
    {
        getRuntime().mInverse->setContextPtr( mContext );
    }
    else
    {
        SCAI_LOG_WARN( logger, "setContextPtr on uninitialized solver" )
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
    lama::L2Norm l2Norm;
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

/* --------------------------------------------------------------------------- */

void InverseSolver::writeAt( std::ostream& stream ) const
{
    stream << "InverseSolver ( id = " << mId << " )";
}

/* --------------------------------------------------------------------------- */

std::string InverseSolver::createValue()
{
    return "InverseSolver";
}

Solver* InverseSolver::create( const std::string name )
{
    return new InverseSolver( name );
}

} /* end namespace solver */

} /* end namespace scai */
