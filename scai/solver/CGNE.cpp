/**
 * @file CGNE.cpp
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
 * @brief CGNE.cpp
 * @author David Schissler
 * @date 27.05.2015
 */

// hpp
#include <scai/solver/CGNE.hpp>

// internal scai libraries
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/common/Constants.hpp>
#include <scai/lama/DenseVector.hpp>

// std
#include <limits>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( CGNE::logger, "Solver.CGNE" )

using lama::Matrix;
using lama::_Vector;
using lama::Scalar;

CGNE::CGNE( const std::string& id )
    : IterativeSolver( id )
{
}


CGNE::CGNE( const std::string& id, LoggerPtr logger )
    : IterativeSolver( id , logger )
{
}

CGNE::CGNE( const CGNE& other )
    : IterativeSolver( other )
{
}

CGNE::CGNERuntime::CGNERuntime()
    : IterativeSolverRuntime()
{
}

CGNE::~CGNE()
{
}

CGNE::CGNERuntime::~CGNERuntime() {}


void CGNE::initialize( const Matrix& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )
    IterativeSolver::initialize( coefficients );
    CGNERuntime& runtime = getRuntime();
    runtime.mEps = Scalar::eps1( coefficients.getValueType() ) * 3.0;
    runtime.mTransposedMat.reset( coefficients.newMatrix() );
    runtime.mTransposedMat->assignTranspose( coefficients );
    runtime.mTransposedMat->conj();
    // get runtime vector with same type / row distribution / context as coefficients
    dmemo::DistributionPtr rowDist( coefficients.getRowDistributionPtr() );
    runtime.mVecP.reset( coefficients.newVector( rowDist ) );
    runtime.mVecZ.reset( coefficients.newVector( rowDist ) );
}


void CGNE::solveInit( _Vector& solution, const _Vector& rhs )
{
    CGNERuntime& runtime = getRuntime();
    runtime.mRhs = &rhs;
    runtime.mSolution = &solution;
    SCAI_ASSERT( runtime.mCoefficients, "solver not initialized" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumRows(), rhs.size(), "size mismatch" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumColumns(), solution.size(), "size mismatch" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getRowDistribution(), rhs.getDistribution(), "distribution mismatch" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getColDistribution(), solution.getDistribution(), "distribution mismatch" )
    // Initialize
    this->getResidual();

    // PRECONDITIONING
    if ( mPreconditioner != NULL )
    {
        runtime.mVecZ->setSameValue( mPreconditioner->getCoefficients().getColDistributionPtr(), 0 );
        mPreconditioner->solve( *runtime.mVecZ, *runtime.mResidual );
    }
    else
    {
        *runtime.mVecZ = *runtime.mResidual;
    }

    *runtime.mVecP = ( *runtime.mTransposedMat ) * ( *runtime.mVecZ );
    runtime.mSolveInit = true;
}

void CGNE::iterate()
{
    CGNERuntime& runtime = getRuntime();
    const Matrix& A = *runtime.mCoefficients;
    const Matrix& transposedA = *runtime.mTransposedMat;
    _Vector& vecP = *runtime.mVecP;
    _Vector& residual = *runtime.mResidual;
    _Vector& solution = *runtime.mSolution;
    _Vector& vecZ = *runtime.mVecZ;
    Scalar alpha;
    Scalar beta;
    Scalar eps = runtime.mEps;
    Scalar scalarProductP = vecP.dotProduct( vecP );
    Scalar scalarProductZR = vecZ.dotProduct( residual );

    SCAI_LOG_INFO( logger, "scalarProductP = " << scalarProductP << ", scalarProductZR = " << scalarProductZR )

    if ( scalarProductP < eps )
    {
        alpha = 0.0;    //norm is small

        SCAI_LOG_INFO( logger, "alpha = 0, as scalarProductP ( " << scalarProductP << " ) < eps ( = " << eps << " )" )
    }
    else
    {
        alpha = scalarProductZR / scalarProductP;

        SCAI_LOG_INFO( logger, "alpha ( = " << alpha << " ) = scalarProductZR ( = " << scalarProductZR
                       << " scalarProductP ( = " << scalarProductP << " )" )
    }

    solution = solution + alpha * vecP;
    residual = residual - alpha * A * vecP;

    // PRECONDITIONING

    if ( mPreconditioner != NULL )
    {
        mPreconditioner->solve( vecZ, residual );
    }
    else
    {
        vecZ = residual;
    }

    if ( scalarProductZR < eps )
    {
        beta = 0.0;    //norm is small

        SCAI_LOG_INFO( logger, "beta = 0, as scalarProductZR ( " << scalarProductZR << " ) < eps ( = " << eps << " )" )
    }
    else
    {
        beta = vecZ.dotProduct( residual ) / scalarProductZR;

        SCAI_LOG_INFO( logger, "beta = " << beta )
    }

    vecP = transposedA * vecZ + beta * vecP;

    // CGNE Implementation End

    mCGNERuntime.mSolution.setDirty( false );
}

SolverPtr CGNE::copy()
{
    return SolverPtr( new CGNE( *this ) );
}

CGNE::CGNERuntime& CGNE::getRuntime()
{
    return mCGNERuntime;
}

const CGNE::CGNERuntime& CGNE::getConstRuntime() const
{
    return mCGNERuntime;
}

std::string CGNE::createValue()
{
    return "CGNE";
}

Solver* CGNE::create( const std::string name )
{
    return new CGNE( name );
}

} /* end namespace solver */

} /* end namespace scai */
