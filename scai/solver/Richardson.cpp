/**
 * @file Richardson.cpp
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
 * @brief Implementation of methods for the Richardson solver.
 * @author David Schissler
 * @date 17.04.2015
 */

// hpp
#include <scai/solver/Richardson.hpp>

// scai internal libraries
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/norm/L2Norm.hpp>
#include <scai/lama/Vector.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, Richardson<ValueType>::logger, "Solver.Richardson" )

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_Solver* Richardson<ValueType>::create()
{
    return new Richardson<ValueType>( "_genByFactory" );
}

template<typename ValueType>
SolverCreateKeyType Richardson<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "Richardson" );
}

/* ========================================================================= */
/*    Constructors / Destructors                                             */
/* ========================================================================= */

template<typename ValueType>
Richardson<ValueType>::Richardson( const std::string& id ) : 

    OmegaSolver<ValueType>( id, ValueType( -1 ) ) 
{
}

template<typename ValueType>
Richardson<ValueType>::Richardson( const std::string& id, const ValueType omega ) : 

    OmegaSolver<ValueType>( id, omega )
{
}

template<typename ValueType>
Richardson<ValueType>::Richardson( const std::string& id, LoggerPtr logger ) : 

    OmegaSolver<ValueType>( id , ValueType( -1 ), logger )
{
}

template<typename ValueType>
Richardson<ValueType>::Richardson( const std::string& id, ValueType( omega ), LoggerPtr logger ) : 

    OmegaSolver<ValueType>( id, omega, logger )
{
}

template<typename ValueType>
Richardson<ValueType>::Richardson( const Richardson& other ) : 

    OmegaSolver<ValueType>( other )
{
}

template<typename ValueType>
Richardson<ValueType>::~Richardson()
{
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void Richardson<ValueType>::initialize( const lama::Matrix<ValueType>& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )

    IterativeSolver<ValueType>::initialize( coefficients );

    if ( OmegaSolver<ValueType>::getOmega() == ValueType( -1 )  )
    {
        ValueType bound = 2 / coefficients.l2Norm();
        bound *= ValueType( 2 ) / ValueType( 3 );
        OmegaSolver<ValueType>::setOmega( bound );
    }

    hmemo::ContextPtr ctx = coefficients.getContextPtr();

    RichardsonRuntime& runtime = getRuntime();

    runtime.mOldSolution.reset( coefficients.newTargetVector() );
    runtime.mX.reset( coefficients.newTargetVector() );
}

/* ========================================================================= */
/*    solve::init( solution, rhs )                                           */
/* ========================================================================= */

template<typename ValueType>
void Richardson<ValueType>::solveInit( lama::Vector<ValueType>& solution, const lama::Vector<ValueType>& rhs )
{
    if ( solution.getVectorKind() != lama::VectorKind::DENSE )
    {
        COMMON_THROWEXCEPTION( "Richardson solver, solution vector must be dense, but is : " << solution )
    }

    IterativeSolver<ValueType>::solveInit( solution, rhs );
}

/* ========================================================================= */
/*    Solver Iteration                                                       */
/* ========================================================================= */

template<typename ValueType>
void Richardson<ValueType>::iterate()
{
    RichardsonRuntime& runtime = getRuntime();

    const lama::Vector<ValueType>& rhs = *runtime.mRhs;
    const lama::Matrix<ValueType>& A   = *runtime.mCoefficients;

    lama::Vector<ValueType>& oldSolution = *runtime.mOldSolution;
    lama::Vector<ValueType>& solution    = runtime.mSolution.getReference(); // dirty

    // swap old solution and solution, so solution is save for update

    solution.swap( oldSolution );

    lama::Vector<ValueType>& x = *runtime.mX;

    x = A * oldSolution;
    solution = rhs - x;

    ValueType omega = OmegaSolver<ValueType>::getOmega();

    if ( omega != ValueType( 1 ) )
    {
        solution *= omega;
    }

    solution += oldSolution;
}

template<typename ValueType>
typename Richardson<ValueType>::RichardsonRuntime& Richardson<ValueType>::getRuntime()
{
    return mRichardsonRuntime;
}

template<typename ValueType>
const typename Richardson<ValueType>::RichardsonRuntime& Richardson<ValueType>::getRuntime() const
{
    return mRichardsonRuntime;
}

template<typename ValueType>
Richardson<ValueType>* Richardson<ValueType>::copy()
{
    return new Richardson<ValueType>( *this );
}

template<typename ValueType>
void Richardson<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "RichardSon<" << common::TypeTraits<ValueType>::id() << "> ( id = " << Solver<ValueType>::getId()
           << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( Richardson, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
