/**
 * @file Solver.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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

// internal scai libraries

#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/tracing.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace solver
{

using lama::Matrix;
using lama::Vector;

/* ========================================================================= */
/*    Solver<ValueType>:  static methods for factory                         */
/* ========================================================================= */

template<typename ValueType>
void Solver<ValueType>::getCreateValues( std::vector<std::string>& values )
{
    _Solver::getTypedCreateValues( values, common::TypeTraits<ValueType>::stype );
}

template<typename ValueType>
Solver<ValueType>* Solver<ValueType>::getSolver( const std::string& solverType )
{
    _Solver* solver = _Solver::getSolver( common::TypeTraits<ValueType>::stype, solverType );

    SCAI_ASSERT_DEBUG( dynamic_cast<Solver<ValueType>*>( solver ), "Illegal solver" )

    return reinterpret_cast<Solver<ValueType>*>( solver );
}

/* ========================================================================= */
/*    Constructor / Destructor of Solver<ValueType>                          */
/* ========================================================================= */

template<typename ValueType>
Solver<ValueType>::Solver( const std::string& id ) : _Solver( id )
{
    // Note: solver runtime calls his own constructor
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Solver<ValueType>::Solver( const std::string& id, LoggerPtr logger ) : _Solver( id, logger )
{
    // Note: solver runtime calls his own constructor
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Solver<ValueType>::Solver( const Solver& other ) : _Solver( other )
{
    // does not copy any runtime stuff here
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
common::ScalarType Solver<ValueType>::getValueType() const
{
    return common::TypeTraits<ValueType>::stype;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
bool Solver<ValueType>::canCreate( const std::string& solverType )
{
    common::ScalarType valType = common::TypeTraits<ValueType>::stype;

    return _Solver::canCreate( SolverCreateKeyType( valType, solverType ) );
}

template<typename ValueType>
Solver<ValueType>::SolverRuntime::SolverRuntime() : 

    mCoefficients( NULL ), 
    mRhs(), 
    mResidual(), 
    mInitialized( false ), 
    mSolveInit( false )
{
}

template<typename ValueType>
Solver<ValueType>::~Solver()
{
    // mRhs, mCoefficents are used as 'dynamic' references, no ownership here
    SCAI_LOG_INFO( _Solver::logger, "~Solver " << this->getId() )
}

template<typename ValueType>
Solver<ValueType>::SolverRuntime::~SolverRuntime()
{
    SCAI_LOG_INFO( logger, "~SolverRuntime" )
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void Solver<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    if ( getRuntime().mInitialized )
    {
        SCAI_LOG_DEBUG( logger, "Previous initialization of solver found! Will be overridden!" )
        mLogger->logMessage( LogLevel::solverInformation, "Solver already initialized, will be overridden\n" );
    }

    getRuntime().mCoefficients = &coefficients;
    getRuntime().mInitialized = true;
    getRuntime().mSolveInit = false;
    getRuntime().mResidual.reset( coefficients.newTargetVector() );
    mLogger->logMessage( LogLevel::solverInformation, "Solver initialized\n" );
}

/* ========================================================================= */
/*    solve( solution, rhs )                                                 */
/* ========================================================================= */

template<typename ValueType>
void Solver<ValueType>::solve( Vector<ValueType>& solution, const Vector<ValueType>& rhs )
{
    SCAI_REGION( "Solver.solve" )

    const SolverRuntime& runtime = getRuntime();

    SCAI_ASSERT_ERROR( runtime.mInitialized, "Solver not initialized, solve cannot be called" )

    solveInit( solution, rhs );

    // check to guarantee that derived classes have called also solveInit of this base class

    SCAI_ASSERT_ERROR( runtime.mSolveInit, "Solver::solveInit not called by derived class." )

    solveImpl();

    solveFinalize();
}

/* ========================================================================= */
/*    solve: init( solution, rhs )                                           */
/* ========================================================================= */

template<typename ValueType>
void Solver<ValueType>::solveInit( Vector<ValueType>& solution, const Vector<ValueType>& rhs )
{
    SolverRuntime& runtime = getRuntime();

    const Matrix<ValueType>& m = *runtime.mCoefficients;

    SCAI_ASSERT_EQ_ERROR( m.getColDistribution(), solution.getDistribution(), "mismatch source space" )
    SCAI_ASSERT_EQ_ERROR( m.getRowDistribution(), rhs.getDistribution(), "mismatch target space" )

    runtime.mRhs = &rhs;
    runtime.mSolution = &solution;

    runtime.mSolveInit = true;
}

template<typename ValueType>
void Solver<ValueType>::solveFinalize()
{
}

/* ========================================================================= */
/*    getResidual() : only recompute if solution has been modified           */
/* ========================================================================= */

template<typename ValueType>
const Vector<ValueType>& Solver<ValueType>::getResidual() const
{
    const SolverRuntime& runtime = getRuntime();

    // initialize and solveInit must have been called before

    SCAI_ASSERT_DEBUG( runtime.mCoefficients, "mCoefficients == NULL" )
    SCAI_ASSERT_DEBUG( runtime.mRhs, "mRhs == NULL" )

    // residual only computed if solution has changed in the meantime

    if ( runtime.mSolution.isDirty() )
    {
        SCAI_REGION( "Solver.computeResidual" )

        const Vector<ValueType>& x = runtime.mSolution.getConstReference();  // solution not modified here

        SCAI_LOG_DEBUG( logger, "calculating residual current solution = " << x )

        const Matrix<ValueType>& A = *runtime.mCoefficients;
        const Vector<ValueType>& y = *runtime.mRhs;

        mLogger->startTimer( "ResidualTimer" );
        *runtime.mResidual = y - A  * x;
        mLogger->stopTimer( "ResidualTimer" );
        mLogger->logTime( "ResidualTimer", LogLevel::completeInformation, "Revaluation of residual took [s]: " );
        mLogger->stopAndResetTimer( "ResidualTimer" );
        runtime.mSolution.setDirty( false );
    }

    return *runtime.mResidual;
}

template<typename ValueType>
const Matrix<ValueType>& Solver<ValueType>::getCoefficients() const
{
    SCAI_ASSERT_DEBUG( getRuntime().mCoefficients, "mCoefficents == NULL" )
    return *getRuntime().mCoefficients;
}

template<typename ValueType>
void Solver<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "Solver ( id = " << getId() << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( Solver, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
