/**
 * @file DefaultJacobi.cpp
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
 * @brief DefaultJacobi.cpp
 * @author Kai Buschulte
 * @date 10.08.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/solver/DefaultJacobi.hpp>

// others
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/DenseVector.hpp>

using namespace memory;

namespace lama
{

SCAI_LOG_DEF_LOGGER( DefaultJacobi::logger, "Jacobi.DefaultJacobi" )

DefaultJacobi::DefaultJacobi( const std::string& id )
    : OmegaSolver( id )
{
}

DefaultJacobi::DefaultJacobi( const std::string& id, const Scalar omega )
    : OmegaSolver( id, omega )
{
}

DefaultJacobi::DefaultJacobi( const std::string& id, LoggerPtr logger )
    : OmegaSolver( id, logger )
{
}

DefaultJacobi::DefaultJacobi( const std::string& id, const Scalar omega, LoggerPtr logger )
    : OmegaSolver( id, omega, logger )
{
}

DefaultJacobi::DefaultJacobi( const DefaultJacobi& other )
    : OmegaSolver( other )
{
}

DefaultJacobi::DefaultJacobiRuntime::DefaultJacobiRuntime()
    : OmegaSolverRuntime(), mDiagonalTimesLU(), mDiagonalTimesRhs(), mOldSolution()
{
}

DefaultJacobi::~DefaultJacobi()
{
}

DefaultJacobi::DefaultJacobiRuntime::~DefaultJacobiRuntime()
{
}

void DefaultJacobi::initialize( const Matrix& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )

    DefaultJacobiRuntime& runtime = getRuntime();

    if( !runtime.mDiagonalTimesRhs.get() )
    {
        runtime.mDiagonalTimesRhs.reset(
            Vector::createVector( coefficients.getValueType(), coefficients.getDistributionPtr() ) );

        SCAI_LOG_DEBUG( logger, "Created diagonalTimesRhs vector = " << *runtime.mDiagonalTimesRhs )
    }

    SCAI_LOG_DEBUG( logger, "Diagonal property of coefficients: " << coefficients.hasDiagonalProperty() )

    coefficients.getDiagonal( *runtime.mDiagonalTimesRhs );
    runtime.mDiagonalTimesRhs->setContext( coefficients.getContextPtr() );
    SCAI_LOG_DEBUG( logger, "Got diagonalTimesRhs = " << *runtime.mDiagonalTimesRhs )
    runtime.mDiagonalTimesRhs->invert();
    SCAI_LOG_DEBUG( logger, "Inverted diagonalTimesRhs = " << *runtime.mDiagonalTimesRhs )

    SCAI_LOG_DEBUG( logger, "Copying main system matrix " << coefficients )
    runtime.mDiagonalTimesLU.reset( coefficients.copy() );
    SCAI_LOG_DEBUG( logger, "Copied main system matrix : " << *runtime.mDiagonalTimesLU )

    SCAI_LOG_DEBUG( logger,
                    "diagonal property of mDiagonalTimesLU: " << runtime.mDiagonalTimesLU->hasDiagonalProperty() )

    runtime.mDiagonalTimesLU->setDiagonal( 0.0 );
    runtime.mDiagonalTimesLU->scale( *runtime.mDiagonalTimesRhs );

    SCAI_LOG_DEBUG( logger, "Create diagonal matrix" )
    runtime.mDiagonalInverted.reset( coefficients.clone() ); // zero matrix with same storage type
    runtime.mDiagonalInverted->setIdentity( coefficients.getDistributionPtr() );
    SCAI_LOG_DEBUG( logger, "identity diagonal matrix = " << *runtime.mDiagonalInverted )
    runtime.mDiagonalInverted->inheritAttributes( coefficients );

    SCAI_LOG_DEBUG( logger,
                    "diagonal property of mDiagonalInverted: " << runtime.mDiagonalInverted->hasDiagonalProperty() )

    runtime.mDiagonalInverted->setDiagonal( *runtime.mDiagonalTimesRhs );

    runtime.mOldSolution.reset( runtime.mDiagonalTimesRhs->clone() );
    runtime.mOldSolution->setContext( runtime.mDiagonalTimesRhs->getContext() );

    OmegaSolver::initialize( coefficients );

    SCAI_LOG_DEBUG( logger, "Initialization performed" )
}

void DefaultJacobi::solve( Vector& solution, const Vector& rhs )
{
    if( getConstRuntime().mSolveInit )
    {
        SCAI_LOG_WARN( logger, "Previous initialization of solver found! Will be overriden!" )
    }

    solveInit( solution, rhs );
    solveImpl();
    solveFinalize();
}

void DefaultJacobi::solveInit( Vector& solution, const Vector& rhs )
{
    DefaultJacobiRuntime& runtime = getRuntime();

    //Check if oldSolution already exists, if not create copy of solution
    if( !runtime.mOldSolution.get() )
    {
        runtime.mOldSolution.reset( solution.clone() );
    }

    runtime.mProxyOldSolution = runtime.mOldSolution.get();

    if( !runtime.mDiagonalTimesRhs.get() || !runtime.mDiagonalInverted.get() )
    {
        COMMON_THROWEXCEPTION( "No initialization executed before running solve." )
    }

    SCAI_LOG_DEBUG( logger, " mDiagonalTimesRhs  =  mDiagonalInverted * rhs , rhs = " << rhs )
    *runtime.mDiagonalTimesRhs = *runtime.mDiagonalInverted * rhs;

    IterativeSolver::solveInit( solution, rhs );
}

void DefaultJacobi::solveFinalize()
{
//    MF: ?????
//    if( &( mProxyOldSolution.getConstReference() ) ==
//        &( mSolution.getConstReference() ) )
    DefaultJacobiRuntime& runtime = getRuntime();

    if( runtime.mIterations % 2 )
    {
        SCAI_LOG_DEBUG( logger, "mProxyOldSolution = *mSolution" )
        *runtime.mProxyOldSolution = *runtime.mSolution;
    }

    SCAI_LOG_DEBUG( logger, " end solve " )
}

template<typename ValueType>
void DefaultJacobi::iterate()
{
    ValueType omega = mOmega.getValue<ValueType>();

    DefaultJacobiRuntime& runtime = getRuntime();

    //swap old solution and solution pointer begin
    Vector* ptr_OldSolution = &( *runtime.mProxyOldSolution );
    Vector* ptr_solution = &( *runtime.mSolution );

    runtime.mProxyOldSolution = ptr_solution;
    runtime.mSolution = ptr_OldSolution;

    //swap end now m_proxOldSolution holds the solution of the last iteration
    //and m_solution will be the output of the current iteration

    const Vector& oldSolution = runtime.mProxyOldSolution.getConstReference();

    SCAI_LOG_DEBUG( logger, " mSolution  =  mDiagonalTimesRhs -  mDiagonalTimesLU * oldSolution " )
    *runtime.mSolution = *runtime.mDiagonalTimesRhs - *runtime.mDiagonalTimesLU * oldSolution;

    if( omega != 1.0 )
    {
        SCAI_LOG_DEBUG( logger, " mSolution = omega * mSolution - (omega - 1.0) * oldSolution " )
        *runtime.mSolution = omega * ( *runtime.mSolution ) - ( omega - 1.0 ) * oldSolution;
    }

    if( SCAI_LOG_TRACE_ON( logger ) )
    {
        SCAI_LOG_TRACE( logger, "Solution " << *runtime.mSolution )
        const DenseVector<ValueType>& sol = dynamic_cast<const DenseVector<ValueType>&>( *runtime.mSolution );
        ReadAccess<ValueType> rsol( sol.getLocalValues() );
        std::cout << "Solution: ";

        for( IndexType i = 0; i < rsol.size(); ++i )
        {
            std::cout << " " << rsol[i];
        }

        std::cout << std::endl;
    }
}

void DefaultJacobi::iterate()
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

DefaultJacobi::DefaultJacobiRuntime& DefaultJacobi::getRuntime()
{
    return mDefaultJacobiRuntime;
}

const DefaultJacobi::DefaultJacobiRuntime& DefaultJacobi::getConstRuntime() const
{
    return mDefaultJacobiRuntime;
}

SolverPtr DefaultJacobi::copy()
{
    return SolverPtr( new DefaultJacobi( *this ) );
}

} // namespace lama
