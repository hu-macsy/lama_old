/**
 * @file InverseSolver.cpp
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
#include <scai/common/macros/instantiate.hpp>

// std
#include <sstream>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, InverseSolver<ValueType>::logger, "Solver.InverseSolver" )

using lama::Matrix;
using lama::Vector;

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_Solver* InverseSolver<ValueType>::create()
{
    return new InverseSolver<ValueType>( "_genByFactory" );
}

template<typename ValueType>
SolverCreateKeyType InverseSolver<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "InverseSolver" );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
InverseSolver<ValueType>::InverseSolver( const std::string& id ) : 

    Solver<ValueType>( id )

{
    SCAI_LOG_INFO( InverseSolver<ValueType>::logger, "InverseSolver, id = " << id )
}

template<typename ValueType>
InverseSolver<ValueType>::InverseSolver( const std::string& id, LoggerPtr logger ) : 

    Solver<ValueType>( id, logger )

{
    SCAI_LOG_INFO( InverseSolver<ValueType>::logger, "InverseSolver, id = " << id )
}

template<typename ValueType>
InverseSolver<ValueType>::InverseSolver( const InverseSolver& other ) : 

    Solver<ValueType>( other )

{
    SCAI_LOG_INFO( InverseSolver<ValueType>::logger, "InverseSolver, id = " << other.getId() )
}

template<typename ValueType>
InverseSolver<ValueType>::~InverseSolver()
{
    SCAI_LOG_INFO( logger, "~InverseSolver" )
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void InverseSolver<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    SCAI_REGION( "Solver.Inverse.intialize" )

    Solver<ValueType>::initialize( coefficients );

    // allocate and build the inverse matrix, might be very time consuming

    SCAI_LOG_INFO( logger, "Initializing with " << coefficients )

    getRuntime().mInverse.reset( coefficients.newMatrix() );

    Matrix<ValueType>& inverse = *getRuntime().mInverse;

    inverse.invert( coefficients );
    inverse.setContextPtr( coefficients.getContextPtr() );
    inverse.prefetch();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const Matrix<ValueType>& InverseSolver<ValueType>::getInverse() const
{
    SCAI_ASSERT_ERROR( getRuntime().mInverse, "inverse not available (no call of initialize before)" );
    return *getRuntime().mInverse;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void InverseSolver<ValueType>::solveImpl()
{
    SCAI_REGION( "Solver.Inverse.solve" )
    InverseSolverRuntime& runtime = getRuntime();
    SCAI_ASSERT_ERROR( runtime.mInverse.get(), "solve, but mInverse is NULL" )
    logStartSolve();

    const Matrix<ValueType>& inverse = *runtime.mInverse;
    const Vector<ValueType>& rhs      = *runtime.mRhs;

    Vector<ValueType>& solution = runtime.mSolution.getReference(); // dirty

    solution = inverse * rhs;

    logEndSolve();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void InverseSolver<ValueType>::logStartSolve()
{
    mLogger->startTimer( "SolutionTimer" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void InverseSolver<ValueType>::logEndSolve()
{
    lama::L2Norm<ValueType> l2Norm;
    mLogger->logResidual( LogLevel::convergenceHistory, *this, l2Norm, "Final " );
    mLogger->logTime( "SolutionTimer", LogLevel::solverInformation, "Total Runtime [s]: " );
    mLogger->stopAndResetTimer( "SolutionTimer" );
    mLogger->logNewLine( LogLevel::solverInformation );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
typename InverseSolver<ValueType>::InverseSolverRuntime& InverseSolver<ValueType>::getRuntime()
{
    return mInverseSolverRuntime;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const typename InverseSolver<ValueType>::InverseSolverRuntime& InverseSolver<ValueType>::getRuntime() const
{
    return mInverseSolverRuntime;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
InverseSolver<ValueType>* InverseSolver<ValueType>::copy()
{
    return new InverseSolver( *this );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void InverseSolver<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "InverseSolver<" << common::TypeTraits<ValueType>::id()
           << "> ( id = " << Solver<ValueType>::getId() << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( InverseSolver, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
