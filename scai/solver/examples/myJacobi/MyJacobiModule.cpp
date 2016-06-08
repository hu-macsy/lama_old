/**
 * @file solver/examples/myJacobi/MyJacobiModule.cpp
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

using namespace scai::hmemo;

MyJacobi::MyJacobi( const std::string& id )
    : scai::solver::OmegaSolver( id )
{
}

MyJacobi::MyJacobi( const std::string& id, const scai::lama::Scalar omega )
    : scai::solver::OmegaSolver( id, omega )
{
}

MyJacobi::MyJacobi( const std::string& id, scai::solver::LoggerPtr logger )
    : scai::solver::OmegaSolver( id, logger )
{
}

MyJacobi::MyJacobi( const std::string& id, const scai::lama::Scalar omega, scai::solver::LoggerPtr logger )
    : scai::solver::OmegaSolver( id, omega, logger )
{
}

MyJacobi::MyJacobi( const MyJacobi& other )
    : scai::solver::OmegaSolver( other )
{
}

MyJacobi::MyJacobiRuntime::MyJacobiRuntime()
    : scai::solver::OmegaSolver::OmegaSolverRuntime(), mDiagonalTimesLU(), mDiagonalTimesRhs(), mOldSolution()
{
}

MyJacobi::~MyJacobi()
{
}

MyJacobi::MyJacobiRuntime::~MyJacobiRuntime()
{
}

void MyJacobi::initialize( const scai::lama::Matrix& coefficients )
{
    MyJacobiRuntime& runtime = getRuntime();

    if ( !runtime.mDiagonalTimesRhs.get() )
    {
        runtime.mDiagonalTimesRhs.reset( coefficients.newDenseVector() );
    }

    coefficients.getDiagonal( *runtime.mDiagonalTimesRhs );
    runtime.mDiagonalTimesRhs->invert();
    runtime.mDiagonalTimesLU.reset( coefficients.copy() );
    runtime.mDiagonalTimesLU->setDiagonal( 0.0 );
    runtime.mDiagonalTimesLU->scale( *runtime.mDiagonalTimesRhs );
    runtime.mDiagonalInverted.reset( coefficients.newMatrix() ); // zero matrix with same storage type
    runtime.mDiagonalInverted->setIdentity( coefficients.getRowDistributionPtr() );
    runtime.mDiagonalInverted->inheritAttributes( coefficients );
    runtime.mDiagonalInverted->setDiagonal( *runtime.mDiagonalTimesRhs );
    runtime.mOldSolution.reset( scai::lama::Vector::create( runtime.mDiagonalTimesRhs->getCreateValue() ) );
    runtime.mOldSolution->setContextPtr( runtime.mDiagonalTimesRhs->getContextPtr() );
    scai::solver::OmegaSolver::initialize( coefficients );
}

void MyJacobi::solve( scai::lama::Vector& solution, const scai::lama::Vector& rhs )
{
    if ( getConstRuntime().mSolveInit )
    {
        std::cout << "Previous initialization of solver found! Will be overriden!" << std::endl;
    }

    solveInit( solution, rhs );
    solveImpl();
    solveFinalize();
}

void MyJacobi::solveInit( scai::lama::Vector& solution, const scai::lama::Vector& rhs )
{
    MyJacobiRuntime& runtime = getRuntime();

    //Check if oldSolution already exists, if not create copy of solution
    if ( !runtime.mOldSolution.get() )
    {
        runtime.mOldSolution.reset( scai::lama::Vector::create( solution.getCreateValue() ) );
    }

    runtime.mProxyOldSolution = runtime.mOldSolution.get();

    if ( !runtime.mDiagonalTimesRhs.get() || !runtime.mDiagonalInverted.get() )
    {
        COMMON_THROWEXCEPTION( "No initialization executed before running solve." )
    }

    *runtime.mDiagonalTimesRhs = *runtime.mDiagonalInverted * rhs;
    IterativeSolver::solveInit( solution, rhs );
}

void MyJacobi::solveFinalize()
{
    MyJacobiRuntime& runtime = getRuntime();

    if ( runtime.mIterations % 2 )
    {
        *runtime.mProxyOldSolution = *runtime.mSolution;
    }
}

template<typename ValueType>
void MyJacobi::iterate()
{
    ValueType omega = mOmega.getValue<ValueType>();
    MyJacobiRuntime& runtime = getRuntime();
    //swap old solution and solution pointer begin
    scai::lama::Vector* ptr_OldSolution = &( *runtime.mProxyOldSolution );
    scai::lama::Vector* ptr_solution = &( *runtime.mSolution );
    runtime.mProxyOldSolution = ptr_solution;
    runtime.mSolution = ptr_OldSolution;
    //swap end now m_proxOldSolution holds the solution of the last iteration
    //and m_solution will be the output of the current iteration
    const scai::lama::Vector& oldSolution = runtime.mProxyOldSolution.getConstReference();
    *runtime.mSolution = *runtime.mDiagonalTimesRhs - *runtime.mDiagonalTimesLU * oldSolution;

    if ( omega != 1.0 )
    {
        *runtime.mSolution = omega * ( *runtime.mSolution ) - ( omega - 1.0 ) * oldSolution;
    }
}

void MyJacobi::iterate()
{
    switch ( getRuntime().mDiagonalTimesLU->getValueType() )
    {
        case scai::common::scalar::FLOAT:
            iterate<float>();
            break;

        case scai::common::scalar::DOUBLE:
            iterate<double>();
            break;

        default:
            COMMON_THROWEXCEPTION( "Unsupported ValueType " << getRuntime().mDiagonalTimesLU->getValueType() )
    }
}

MyJacobi::MyJacobiRuntime& MyJacobi::getRuntime()
{
    return mMyJacobiRuntime;
}

const MyJacobi::MyJacobiRuntime& MyJacobi::getConstRuntime() const
{
    return mMyJacobiRuntime;
}

scai::solver::SolverPtr MyJacobi::copy()
{
    return scai::solver::SolverPtr( new MyJacobi( *this ) );
}

void MyJacobi::writeAt( std::ostream& stream ) const
{
    stream << "MyJacobi ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

std::string MyJacobi::createValue()
{
    return "MyJacobi";
}

scai::solver::Solver* MyJacobi::create( const std::string name )
{
    return new MyJacobi( name );
}
