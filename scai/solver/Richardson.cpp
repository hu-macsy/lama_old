/**
 * @file Richardson.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Richardson.cpp
 * @author David Schissler
 * @date 17.04.2015
 * @since
 */

// hpp
#include <scai/solver/Richardson.hpp>

// local library
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/DenseVector.hpp>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( Richardson::logger, "Solver.Richardson" )

Richardson::Richardson( const std::string& id )
    : OmegaSolver( id, ( lama::Scalar ) - 1.0 ) {}

Richardson::Richardson( const std::string& id, const lama::Scalar omega )
    : OmegaSolver( id, omega ) {}

Richardson::Richardson( const std::string& id, LoggerPtr logger )
    : OmegaSolver( id , ( lama::Scalar ) - 1.0, logger ) {}

Richardson::Richardson( const std::string& id, const lama::Scalar omega, LoggerPtr logger )
    : OmegaSolver( id, omega, logger ) {}

Richardson::Richardson( const Richardson& other )
    : OmegaSolver( other ) {}



Richardson::RichardsonRuntime::RichardsonRuntime()
    : OmegaSolverRuntime() {}

Richardson::~Richardson() {}

Richardson::RichardsonRuntime::~RichardsonRuntime() {}



void Richardson::initialize( const lama::Matrix& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )

    IterativeSolver::initialize( coefficients );

    if ( mOmega == -1.0 )
    {
        lama::L2Norm n;
        lama::Scalar bound = 2.0 / n.apply( coefficients );
        mOmega = ( 2.0 / 3.0 * bound );
    }
}

void Richardson::solveInit( lama::Vector& solution, const lama::Vector& rhs )
{
    RichardsonRuntime& runtime = getRuntime();

    //Check if oldSolution already exists, if not create copy of solution
    if ( !runtime.mOldSolution.get() )
    {
        runtime.mOldSolution.reset( lama::Vector::create( solution.getCreateValue() ) );
    }

    runtime.mProxyOldSolution = runtime.mOldSolution.get();

    IterativeSolver::solveInit( solution, rhs );
}

void Richardson::solveFinalize()
{
    RichardsonRuntime& runtime = getRuntime();

    if ( runtime.mIterations % 2 )
    {
        SCAI_LOG_DEBUG( logger, "mProxyOldSolution = *mSolution" )
        *runtime.mProxyOldSolution = *runtime.mSolution;
    }

    SCAI_LOG_DEBUG( logger, " end solve " )
}

void Richardson::iterate()
{
    RichardsonRuntime& runtime = getRuntime();

    const lama::Vector& rhs = *runtime.mRhs;
    const lama::Matrix& A = *runtime.mCoefficients;
    //swap old solution and solution pointer begin
    lama::Vector* ptr_OldSolution = &( *runtime.mProxyOldSolution );
    lama::Vector* ptr_solution = &( *runtime.mSolution );

    runtime.mProxyOldSolution = ptr_solution;
    runtime.mSolution = ptr_OldSolution;

    const lama::Vector& oldSolution = runtime.mProxyOldSolution.getConstReference();

    lama::Vector* x = lama::Vector::getVector( lama::Vector::DENSE, A.getValueType() );
    lama::Vector& xRef = *x;
    xRef = A * oldSolution;

    *runtime.mSolution = rhs - xRef;

    if ( mOmega != 1.0 )
    {
        *runtime.mSolution = mOmega * ( *runtime.mSolution );
    }

    *runtime.mSolution = oldSolution + ( *runtime.mSolution );
}

Richardson::RichardsonRuntime& Richardson::getRuntime()
{
    return mRichardsonRuntime;
}

const Richardson::RichardsonRuntime& Richardson::getConstRuntime() const
{
    return mRichardsonRuntime;
}

SolverPtr Richardson::copy()
{
    return SolverPtr( new Richardson( *this ) );
}

void Richardson::writeAt( std::ostream& stream ) const
{
    stream << "Richardson ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

std::string Richardson::createValue()
{
	return "Richardson";
}

Solver* Richardson::create( const std::string name )
{
	return new Richardson( name );
}

} /* end namespace solver */

} /* end namespace scai */
