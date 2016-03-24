/**
 * @file MyJacobi.cpp
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
 * @brief MyJacobi.cpp
 * @author Kai Buschulte
 * @date 10.08.2011
 * @since 1.0.0
 */

// hpp
#include "myJacobi.hpp"

#include <lama.hpp>

// local library
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/DenseVector.hpp>

namespace scai 
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( MyJacobi::logger, "Jacobi.MyJacobi" )

using namespace hmemo;

MyJacobi::MyJacobi( const std::string& id )
    : OmegaSolver( id )
{
}

MyJacobi::MyJacobi( const std::string& id, const lama::Scalar omega )
    : OmegaSolver( id, omega )
{
}

MyJacobi::MyJacobi( const std::string& id, LoggerPtr logger )
    : OmegaSolver( id, logger )
{
}

MyJacobi::MyJacobi( const std::string& id, const lama::Scalar omega, LoggerPtr logger )
    : OmegaSolver( id, omega, logger )
{
}

MyJacobi::MyJacobi( const MyJacobi& other )
    : OmegaSolver( other )
{
}

MyJacobi::MyJacobiRuntime::MyJacobiRuntime()
    : OmegaSolverRuntime(), mDiagonalTimesLU(), mDiagonalTimesRhs(), mOldSolution()
{
}

MyJacobi::~MyJacobi()
{
}

MyJacobi::MyJacobiRuntime::~MyJacobiRuntime()
{
}

void MyJacobi::initialize( const lama::Matrix& coefficients )
{
    MyJacobiRuntime& runtime = getRuntime();

    if( !runtime.mDiagonalTimesRhs.get() )
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

    runtime.mOldSolution.reset( lama::Vector::create( runtime.mDiagonalTimesRhs->getCreateValue() ) );
    runtime.mOldSolution->setContextPtr( runtime.mDiagonalTimesRhs->getContextPtr() );

    OmegaSolver::initialize( coefficients );
}

void MyJacobi::solve( lama::Vector& solution, const lama::Vector& rhs )
{
    if( getConstRuntime().mSolveInit )
    {
        SCAI_LOG_WARN( logger, "Previous initialization of solver found! Will be overriden!" )
    }

    solveInit( solution, rhs );
    solveImpl();
    solveFinalize();
}

void MyJacobi::solveInit( lama::Vector& solution, const lama::Vector& rhs )
{
    MyJacobiRuntime& runtime = getRuntime();

    //Check if oldSolution already exists, if not create copy of solution
    if( !runtime.mOldSolution.get() )
    {
        runtime.mOldSolution.reset( lama::Vector::create( solution.getCreateValue() ) );
    }

    runtime.mProxyOldSolution = runtime.mOldSolution.get();

    if( !runtime.mDiagonalTimesRhs.get() || !runtime.mDiagonalInverted.get() )
    {
        COMMON_THROWEXCEPTION( "No initialization executed before running solve." )
    }

    *runtime.mDiagonalTimesRhs = *runtime.mDiagonalInverted * rhs;

    IterativeSolver::solveInit( solution, rhs );
}

void MyJacobi::solveFinalize()
{
    MyJacobiRuntime& runtime = getRuntime();

    if( runtime.mIterations % 2 )
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
    lama::Vector* ptr_OldSolution = &( *runtime.mProxyOldSolution );
    lama::Vector* ptr_solution = &( *runtime.mSolution );

    runtime.mProxyOldSolution = ptr_solution;
    runtime.mSolution = ptr_OldSolution;

    //swap end now m_proxOldSolution holds the solution of the last iteration
    //and m_solution will be the output of the current iteration

    const lama::Vector& oldSolution = runtime.mProxyOldSolution.getConstReference();

    *runtime.mSolution = *runtime.mDiagonalTimesRhs - *runtime.mDiagonalTimesLU * oldSolution;

    if( omega != 1.0 )
    {
        *runtime.mSolution = omega * ( *runtime.mSolution ) - ( omega - 1.0 ) * oldSolution;
    }
}

void MyJacobi::iterate()
{
    switch( getRuntime().mDiagonalTimesLU->getValueType() )
    {
        case common::scalar::FLOAT:
            iterate<float>();
            break;

        case common::scalar::DOUBLE:
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

SolverPtr MyJacobi::copy()
{
    return SolverPtr( new MyJacobi( *this ) );
}

void MyJacobi::writeAt( std::ostream& stream ) const
{
    stream << "MyJacobi ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

std::string MyJacobi::createValue()
{
	return "MyJacobi";
}

Solver* MyJacobi::create( const std::string name )
{
	return new MyJacobi( name );
}

} /* end namespace solver */

} /* end namespace scai */
