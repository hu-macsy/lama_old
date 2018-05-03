/**
 * @file solver/examples/myJacobi/MyJacobiModule.cpp
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
 * @brief MyJacobi.cpp
 * @author Kai Buschulte
 * @date 10.08.2011
 */

// hpp
#include "MyJacobiModule.hpp"

// local library
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/DenseVector.hpp>

using namespace scai;
using namespace hmemo;

template<typename ValueType>
MyJacobi<ValueType>::MyJacobi( const std::string& id ) : 

    solver::OmegaSolver<ValueType>( id )
{
}

template<typename ValueType>
MyJacobi<ValueType>::MyJacobi( const std::string& id, const ValueType omega ) : 

    solver::OmegaSolver<ValueType>( id, omega )
{
}

template<typename ValueType>
MyJacobi<ValueType>::MyJacobi( const std::string& id, solver::LoggerPtr logger ) : 

    solver::OmegaSolver<ValueType>( id, logger )
{
}

template<typename ValueType>
MyJacobi<ValueType>::MyJacobi( const std::string& id, const ValueType omega, solver::LoggerPtr logger ) : 

    solver::OmegaSolver<ValueType>( id, omega, logger )
{
}

template<typename ValueType>
MyJacobi<ValueType>::MyJacobi( const MyJacobi<ValueType>& other ) : 

    solver::OmegaSolver<ValueType>( other )
{
}

template<typename ValueType>
MyJacobi<ValueType>::MyJacobiRuntime::MyJacobiRuntime() : 

    solver::IterativeSolver<ValueType>::IterativeSolverRuntime(), 
    mDiagonalTimesLU(), 
    mDiagonalTimesRhs(), 
    mOldSolution()
{
}

template<typename ValueType>
MyJacobi<ValueType>::~MyJacobi()
{
}

template<typename ValueType>
MyJacobi<ValueType>::MyJacobiRuntime::~MyJacobiRuntime()
{
}

template<typename ValueType>
void MyJacobi<ValueType>::initialize( const lama::Matrix<ValueType>& coefficients )
{
    MyJacobiRuntime& runtime = getRuntime();

    // runtime data inherits context of matrix to be solved

    hmemo::ContextPtr ctx = coefficients.getContextPtr();   

    runtime.mDiagonalTimesRhs.setContextPtr( ctx );

    runtime.mDiagonalInverted.setContextPtr( ctx );
    coefficients.getDiagonal( runtime.mDiagonalInverted );
    runtime.mDiagonalInverted.unaryOp( runtime.mDiagonalInverted, common::UnaryOp::RECIPROCAL );

    // mDiagonalTimesLU = ( A - D ) * ( D * rhs )

    runtime.mDiagonalTimesLU.reset( coefficients.copy() );

    lama::Matrix<ValueType>& diagonalTimesLU = *runtime.mDiagonalTimesLU;

    diagonalTimesLU.setDiagonal( ValueType( 0 ) );
    diagonalTimesLU.scaleRows( runtime.mDiagonalInverted );

    runtime.mOldSolution.setContextPtr( ctx );

    solver::Solver<ValueType>::initialize( coefficients );
}

template<typename ValueType>
void MyJacobi<ValueType>::solveInit( lama::Vector<ValueType>& solution, const lama::Vector<ValueType>& rhs )
{
    solver::Solver<ValueType>::solveInit( solution, rhs );

    MyJacobiRuntime& runtime = getRuntime();

    // each iteration needs rhs / D, with D = diagonal(A), so computed it once

    runtime.mDiagonalTimesRhs = runtime.mDiagonalInverted * rhs;  // cwiseProduct
}

template<typename ValueType>
void MyJacobi<ValueType>::iterate()
{
    MyJacobiRuntime& runtime = getRuntime();

    // swap old solution and solution, solution is marked as dirty

    lama::Vector<ValueType>& solution    = runtime.mSolution.getReference();
    lama::Vector<ValueType>& oldSolution = runtime.mOldSolution;

    solution.swap( oldSolution );

    solution = runtime.mDiagonalTimesRhs - *runtime.mDiagonalTimesLU * oldSolution;

    const ValueType omega = this->getOmega();

    if ( omega != 1 )
    {
        solution = omega * solution - ( omega - 1 ) * oldSolution;
    }
}

template<typename ValueType>
typename MyJacobi<ValueType>::MyJacobiRuntime& MyJacobi<ValueType>::getRuntime()
{
    return mMyJacobiRuntime;
}

template<typename ValueType>
const typename MyJacobi<ValueType>::MyJacobiRuntime& MyJacobi<ValueType>::getRuntime() const
{
    return mMyJacobiRuntime;
}

template<typename ValueType>
MyJacobi<ValueType>* MyJacobi<ValueType>::copy()
{
    return new MyJacobi( *this );
}

template<typename ValueType>
void MyJacobi<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "MyJacobi<" << common::TypeTraits<ValueType>::id() << "> ( id = " << this->getId()
           << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
solver::_Solver* MyJacobi<ValueType>::create()
{
    return new MyJacobi<ValueType>( "_genByFactory" );
}

template<typename ValueType>
solver::SolverCreateKeyType MyJacobi<ValueType>::createValue()
{
    return solver::SolverCreateKeyType( common::getScalarType<ValueType>(), "MyJacobi" );
}


/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

template class MyJacobi<DefaultReal>;
