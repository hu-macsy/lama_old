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

// internal scai libraries

#include <scai/lama/DenseVector.hpp>
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
/*    Implementation of methods for _Solver                                  */
/* ========================================================================= */

SCAI_LOG_DEF_LOGGER( _Solver::logger, "Solver" )

_Solver* _Solver::getSolver( const common::ScalarType scalarType, const std::string& solverType )
{
    return create( SolverCreateKeyType( scalarType, solverType ) );
}

template<typename ValueType>
void Solver<ValueType>::getCreateValues( std::vector<std::string>& values )
{
    std::vector<SolverCreateKeyType> createValues;

    _Solver::getCreateValues( createValues );  // all solvers ( valueType, solvertype )

    values.clear();
 
    for ( size_t i = 0; i < createValues.size(); ++i )
    {
        if ( createValues[i].first == common::TypeTraits<ValueType>::stype )
        {
            // Solver for this value type
            values.push_back( createValues[i].second );
        }
    }
}

/* ========================================================================= */
/*    Constructor / Destructor of Solver<ValueType>                          */
/* ========================================================================= */

template<typename ValueType>
Solver<ValueType>::Solver( const std::string& id ) : 

    mLogger( new CommonLogger( "dummyLog", LogLevel::noLogging,
                               LoggerWriteBehaviour::toConsoleOnly,
                               std::shared_ptr<Timer>( new Timer() ) ) ),
    mId( id )

{
    SCAI_LOG_INFO( _Solver::logger, "Solver id = " << mId << " created, dummy log" )
}

template<typename ValueType>
Solver<ValueType>::Solver( const std::string& id, LoggerPtr logger ) : 

    mLogger( logger ),
    mId( id )

{
    SCAI_LOG_INFO( _Solver::logger, "Solver id = " << mId << " created, with logger" )
}

template<typename ValueType>
Solver<ValueType>::Solver( const Solver& other ) : 

    mLogger( other.mLogger ),
    mId( other.mId )
{
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
common::ScalarType Solver<ValueType>::getValueType() const
{
    return common::TypeTraits<ValueType>::stype;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Solver<ValueType>* Solver<ValueType>::getSolver( const std::string& solverType )
{
    _Solver* solver = _Solver::getSolver( common::TypeTraits<ValueType>::stype, solverType );

    SCAI_ASSERT_DEBUG( dynamic_cast<Solver<ValueType>*>( solver ), "Illegal solver" )

    return reinterpret_cast<Solver<ValueType>*>( solver );
}


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
    // mRhs, mCoefficents are used as 'dynamic' references, no free
    SCAI_LOG_INFO( _Solver::logger, "~Solver " << mId )
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
    getRuntime().mResidual.setContextPtr( coefficients.getContextPtr() );           // IMPORTANT
    getRuntime().mResidual.allocate( coefficients.getRowDistributionPtr() );
    mLogger->logMessage( LogLevel::solverInformation, "Solver initialized\n" );
}

template<typename ValueType>
void Solver<ValueType>::solve( Vector<ValueType>& solution, const Vector<ValueType>& rhs )
{
    SCAI_REGION( "Solver.solve" )

    const SolverRuntime& runtime = getRuntime();

    SCAI_ASSERT( runtime.mInitialized, "Solver not initialized, solve cannot be called" )

    solveInit( solution, rhs );

    // check to guarantee that derived classes have called also solveInit of this base class

    SCAI_ASSERT_ERROR( runtime.mSolveInit, "Solver::solveInit not called by derived class." )

    solveImpl();

    solveFinalize();
}

template<typename ValueType>
void Solver<ValueType>::solveInit( Vector<ValueType>& solution, const Vector<ValueType>& rhs )
{
    SolverRuntime& runtime = getRuntime();

    const Matrix<ValueType>& m = *runtime.mCoefficients;

    SCAI_ASSERT_EQ_ERROR( m.getColDistribution(), solution.getDistribution(), "mismatch: matrix col dist, solution" )
    SCAI_ASSERT_EQ_ERROR( m.getRowDistribution(), rhs.getDistribution(), "mismatch: matrix row dist, rhs dist" )

    runtime.mRhs = &rhs;
    runtime.mSolution = &solution;

    runtime.mSolveInit = true;
}

template<typename ValueType>
void Solver<ValueType>::solveFinalize()
{
}

template<typename ValueType>
const std::string& Solver<ValueType>::getId() const
{
    SCAI_LOG_TRACE( logger, "Returning Solver Id " << mId )
    return mId;
}

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
        runtime.mResidual = y - A  * x;
        mLogger->stopTimer( "ResidualTimer" );
        mLogger->logTime( "ResidualTimer", LogLevel::completeInformation, "Revaluation of residual took [s]: " );
        mLogger->stopAndResetTimer( "ResidualTimer" );
        runtime.mSolution.setDirty( false );
    }

    return runtime.mResidual;
}

template<typename ValueType>
const Matrix<ValueType>& Solver<ValueType>::getCoefficients() const
{
    SCAI_ASSERT_DEBUG( getRuntime().mCoefficients, "mCoefficents == NULL" )
    return *getRuntime().mCoefficients;
}

template<typename ValueType>
void Solver<ValueType>::setLogger( LoggerPtr logger )
{
    mLogger = logger;
}

template<typename ValueType>
void Solver<ValueType>::setLogLevel( LogLevel::LogLevel level )
{
    mLogger->setLogLevel( level );
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
