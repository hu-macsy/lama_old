/**
 * @file IterativeSolver.cpp
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

#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, IterativeSolver<ValueType>::logger, "Solver.IterativeSolver" )

using lama::Matrix;

template<typename ValueType>
IterativeSolver<ValueType>::IterativeSolver( const std::string& id ) : 

    Solver<ValueType>( id ), 
    mCriterionRootComponent( new IterationCount<ValueType>( 1 ) )
{
}

template<typename ValueType>
IterativeSolver<ValueType>::IterativeSolver( const std::string& id, LoggerPtr logger ) : 

    Solver<ValueType>( id, logger ), 
    mCriterionRootComponent( new IterationCount<ValueType>( 1 ) )
{
}

template<typename ValueType>
IterativeSolver<ValueType>::IterativeSolver( const IterativeSolver& other ) : 

    Solver<ValueType>( other )

{
    if ( mCriterionRootComponent )
    {
        mCriterionRootComponent.reset( new Criterion<ValueType>( *other.mCriterionRootComponent ) );
    }
    else
    {
        SCAI_LOG_INFO( logger, other << " has no conditions, that can be copied." )
        mCriterionRootComponent.reset( new IterationCount<ValueType>( 1 ) );
    }

    if ( other.getPreconditioner() )
    {
        mPreconditioner.reset( other.getPreconditioner()->copy() );
    }
    else
    {
        SCAI_LOG_INFO( logger, other << " has no preconditioner." )
    }
}

template<typename ValueType>
IterativeSolver<ValueType>::IterativeSolverRuntime::IterativeSolverRuntime() : 
 
    mIterations( 0 )
{
}

template<typename ValueType>
IterativeSolver<ValueType>::~IterativeSolver()
{
}

template<typename ValueType>
IterativeSolver<ValueType>::IterativeSolverRuntime::~IterativeSolverRuntime()
{
}

template<typename ValueType>
void IterativeSolver<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    Solver<ValueType>::initialize( coefficients );

    if ( mPreconditioner )
    {
        if ( mPreconditioner->getRuntime().mInitialized )
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

template<typename ValueType>
void IterativeSolver<ValueType>::solveImpl()
{
    getRuntime().mIterations = 0;

    if ( !getRuntime().mSolveInit )
    {
        COMMON_THROWEXCEPTION(
            "Solver " + this->getId()
            + " has not been initialized. Call solveInit( Vector<T>& solution, const Vector<T>& rhs ) before solving "
            + this->getId() )
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

template<typename ValueType>
void IterativeSolver<ValueType>::setStoppingCriterion( const CriterionPtr<ValueType> criterion )
{
    SCAI_ASSERT_ERROR( criterion, "Criterion defined is NULL." )
    SCAI_LOG_INFO( logger, "Criteria " << *criterion << " defined." )
    mCriterionRootComponent = criterion;
}

template<typename ValueType>
bool IterativeSolver<ValueType>::criteriaAreSatisfied() const
{
    if ( mCriterionRootComponent.get() == 0 )
    {
        COMMON_THROWEXCEPTION( this->getId() + ": No stopping criterion set." )
    }

    return mCriterionRootComponent->isSatisfied( *this );
}

template<typename ValueType>
void IterativeSolver<ValueType>::setPreconditioner( SolverPtr<ValueType> const conditioner )
{
    SCAI_LOG_INFO( logger, "Preconditioner " << conditioner->getId() << " defined." )
    mPreconditioner = conditioner;
}

template<typename ValueType>
const SolverPtr<ValueType> IterativeSolver<ValueType>::getPreconditioner() const
{
    return mPreconditioner;
}

template<typename ValueType>
IndexType IterativeSolver<ValueType>::getIterationCount() const
{
    return this->getRuntime().mIterations;
}

template<typename ValueType>
void IterativeSolver<ValueType>::logStartSolve()
{
    mLogger->logNewLine( LogLevel::solverInformation );
    mLogger->logNewLine( LogLevel::solverInformation );
    mLogger->logMessage( LogLevel::solverInformation, "Beginning solve.\n" );
    lama::L2Norm<ValueType> l2Norm;
    mLogger->logResidual( LogLevel::convergenceHistory, *this, l2Norm, "Start " );
    mLogger->startTimer( "SolutionTimer" );
}

template<typename ValueType>
void IterativeSolver<ValueType>::logIterationEndAndResidual()
{
    IterativeSolverRuntime& runtime = getRuntime();
    mLogger->stopTimer( "IterationTimer" );
    std::stringstream iterationPrefix;
    iterationPrefix << "Iteration #" << runtime.mIterations << " Duration [s]: ";
    mLogger->logTime( "IterationTimer", LogLevel::advancedInformation, iterationPrefix.str() );
    mLogger->stopAndResetTimer( "IterationTimer" );
    lama::L2Norm<ValueType> l2Norm;
    iterationPrefix.str( "" );
    iterationPrefix << "Iteration #" << runtime.mIterations << " ";
    mLogger->logResidual( LogLevel::convergenceHistory, *this, l2Norm, iterationPrefix.str() );
}

template<typename ValueType>
void IterativeSolver<ValueType>::logEndSolve()
{
    mLogger->logNewLine( LogLevel::solverInformation );
    mLogger->logMessage( LogLevel::solverInformation, "-------------------------\n" );
    mLogger->logType( LogLevel::solverInformation, "Solve completed. Iterations: ", getRuntime().mIterations );
    mLogger->logTime( "SolutionTimer", LogLevel::solverInformation, "Total Runtime [s]: " );
    mLogger->stopAndResetTimer( "SolutionTimer" );
    mLogger->logNewLine( LogLevel::solverInformation );
}

template<typename ValueType>
void IterativeSolver<ValueType>::logIterationStart()
{
    mLogger->logNewLine( LogLevel::completeInformation );
    mLogger->logType( LogLevel::completeInformation, "Beginning iteration #", getRuntime().mIterations + 1 );
    mLogger->startTimer( "IterationTimer" );
}

template<typename ValueType>
IterativeSolver<ValueType>* IterativeSolver<ValueType>::getSolver( const std::string& solverType )
{
    // use of unique_ptr makes sure that the generated solver will be deleted if an exception is thrown

    std::unique_ptr<Solver<ValueType> > newSolver( Solver<ValueType>::getSolver( solverType ) );
    
    IterativeSolver<ValueType>* itSolver = dynamic_cast<IterativeSolver<ValueType>*>( newSolver.get() );

    if ( !itSolver )
    {
        COMMON_THROWEXCEPTION( "requested Solver is not inherited from IterativeSolver" )
    }

    return reinterpret_cast<IterativeSolver<ValueType>*>( newSolver.release() );
}

template<typename ValueType>
void IterativeSolver<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "IterativeSolver ( id = " << this->getId() << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( IterativeSolver, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
