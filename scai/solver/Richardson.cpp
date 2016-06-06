/**
 * @file Richardson.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Richardson.cpp
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
#include <scai/lama/DenseVector.hpp>

#include <scai/common/unique_ptr.hpp>

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

    if( !runtime.mX.get() )
    {
        runtime.mX.reset( lama::Vector::create( solution.getCreateValue() ) );
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

    lama::Vector& xRef = *runtime.mX;
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
