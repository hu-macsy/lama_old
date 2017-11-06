/**
 * @file TFQMR.cpp
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
 * @brief Implementation of methods for the TFQMR solver.
 * @author David Schissler
 * @date 13.05.2015
 */

// hpp
#include <scai/solver/TFQMR.hpp>

// internal scai libraries
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/DenseVector.hpp>

// std
#include <limits>

namespace scai
{

using lama::_Matrix;
using lama::_Vector;
using lama::Scalar;

namespace solver
{

SCAI_LOG_DEF_LOGGER( TFQMR::logger, "Solver.TFQMR" )

TFQMR::TFQMR( const std::string& id )
    : IterativeSolver( id ) {}


TFQMR::TFQMR( const std::string& id, LoggerPtr logger )
    : IterativeSolver( id , logger ) {}

TFQMR::TFQMR( const TFQMR& other )
    : IterativeSolver( other ) {}



TFQMR::TFQMRRuntime::TFQMRRuntime()
    : IterativeSolverRuntime() {}

TFQMR::~TFQMR() {}

TFQMR::TFQMRRuntime::~TFQMRRuntime() {}

void TFQMR::initialize( const _Matrix& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )
    IterativeSolver::initialize( coefficients );
    TFQMRRuntime& runtime = getRuntime();
    runtime.mAlpha = 0.0;
    runtime.mBeta = 0.0;
    runtime.mC = 0.0;
    runtime.mEta = 0.0;
    runtime.mTheta = 0.0;
    runtime.mEps = Scalar::eps1( coefficients.getValueType() ) * 3.0;
    // create dense runtime vectors with same row distribution, type, context as coefficients
    dmemo::DistributionPtr dist = coefficients.getRowDistributionPtr();
    runtime.mVecD.reset( coefficients.newVector( dist ) );
    runtime.mInitialR.reset( coefficients.newVector( dist ) );
    runtime.mVecVEven.reset( coefficients.newVector( dist ) );
    runtime.mVecVOdd.reset( coefficients.newVector( dist ) );
    runtime.mVecW.reset( coefficients.newVector( dist ) );
    runtime.mVecZ.reset( coefficients.newVector( dist ) );
    runtime.mVecVT.reset( coefficients.newVector( dist ) );
}

void TFQMR::solveInit( _Vector& solution, const _Vector& rhs )
{
    TFQMRRuntime& runtime = getRuntime();
    runtime.mRhs = &rhs;
    runtime.mSolution = &solution;
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumRows(), rhs.size(), "mismatch: #rows of matrix, rhs" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumColumns(), solution.size(), "mismatch: #cols of matrix, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getColDistribution(), solution.getDistribution(), "mismatch: matrix col dist, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getRowDistribution(), rhs.getDistribution(), "mismatch: matrix row dist, rhs dist" )
    // Initialize
    this->getResidual();
    const _Matrix& A = *runtime.mCoefficients;
    *runtime.mInitialR = *runtime.mResidual;
    *runtime.mVecVEven = *runtime.mResidual;

    // PRECONDITIONING
    if ( mPreconditioner != NULL )
    {
        runtime.mVecW->setSameValue( runtime.mResidual->getDistributionPtr(), 0 );
        mPreconditioner->solve( *runtime.mVecW , *runtime.mResidual );
    }
    else
    {
        *runtime.mVecW = *runtime.mResidual;
    }

    *runtime.mVecZ = A * ( *runtime.mVecW );
    *runtime.mVecW = *runtime.mResidual;
    // *runtime.mVecD = 0
    runtime.mVecD->setSameValue( A.getRowDistributionPtr(), 0 );
    lama::L2Norm norm;
    runtime.mTau = norm.apply( *runtime.mInitialR );
    runtime.mRhoOld = runtime.mTau * runtime.mTau;
    runtime.mSolveInit = true;
}

void TFQMR::iterationEven()
{
    TFQMRRuntime& runtime = getRuntime();
    const _Vector& vecZ = *runtime.mVecZ;
    const _Vector& initialR = *runtime.mInitialR;
    const _Vector& vecVEven = *runtime.mVecVEven;
    _Vector& vecVOdd = *runtime.mVecVOdd;
    const Scalar& rho = runtime.mRhoOld;
    const Scalar& eps = runtime.mEps;
    const Scalar dotProduct = initialR._dotProduct( vecZ );
    Scalar& alpha = runtime.mAlpha;

    if ( abs( dotProduct ) < eps ) // scalar is small
    {
        alpha = 0.0;
    }
    else
    {
        alpha = rho / dotProduct;
    }

    vecVOdd  = vecVEven - alpha * vecZ;
}

void TFQMR::iterationOdd()
{
    TFQMRRuntime& runtime = getRuntime();
    const _Matrix& A = *runtime.mCoefficients;
    const _Vector& initialR = *runtime.mInitialR;
    const _Vector& vecW = *runtime.mVecW;
    const _Vector& vecVOdd = *runtime.mVecVOdd;
    _Vector& vecVEven = *runtime.mVecVEven;
    Scalar& rhoOld = runtime.mRhoOld;
    Scalar& rhoNew = runtime.mRhoNew;
    Scalar& beta = runtime.mBeta;
    _Vector& vecZ = *runtime.mVecZ;
    _Vector& vecVT = *runtime.mVecVT;
    const Scalar& eps = runtime.mEps;
    rhoNew  = initialR._dotProduct( vecW );

    if ( abs( rhoOld ) < eps )          // scalar is small
    {
        beta = 0.0;
    }
    else
    {
        beta = rhoNew / rhoOld;
    }

    vecVEven = vecW + beta * vecVOdd;

//  PRECONDITIONING
    if ( mPreconditioner != NULL )
    {
        vecVT.setSameValue( mPreconditioner->getCoefficients().getColDistributionPtr(), 0 );
        mPreconditioner->solve( vecVT, vecVEven );
    }
    else
    {
        vecVT = vecVEven;
    }

    vecZ *= beta;
    vecZ = beta * A * vecVOdd + beta * vecZ;
    vecZ = A * vecVT + vecZ;
    rhoOld = rhoNew;
}

void TFQMR::iterate()
{
    TFQMRRuntime& runtime = getRuntime();
    const IndexType& iteration = runtime.mIterations;
    lama::L2Norm norm;
    const _Matrix& A = *runtime.mCoefficients;
    _Vector& vecW = *runtime.mVecW;
    _Vector& vecD = *runtime.mVecD;
    _Vector& solution = *runtime.mSolution;
    const Scalar& alpha = runtime.mAlpha;
    Scalar& c = runtime.mC;
    Scalar& eta = runtime.mEta;
    Scalar& theta = runtime.mTheta;
    Scalar& tau = runtime.mTau;
    Scalar& eps = runtime.mEps;
    const _Vector* vecVp;

    // if iteration = even -> need mVecVEven, else need mVecVOdd
    if ( ( iteration % 2 ) == 0 )
    {
        vecVp = &( *runtime.mVecVEven );
    }
    else
    {
        vecVp = &( *runtime.mVecVOdd );
    }

    const _Vector& vecV = *vecVp;

    if ( ( iteration % 2 ) == 0 )
    {
        iterationEven();
    }

    vecW = vecW - alpha * A * vecV;
    Scalar tempScal;

    if ( abs( alpha ) < eps || abs( theta ) < eps || abs( eta ) < eps ) // scalar is small
    {
        tempScal = 0.0;
    }
    else
    {
        tempScal = ( theta * theta / alpha ) * eta;
    }

    vecD = vecV + tempScal * vecD;

    if ( abs( tau ) < eps ) // scalar is small
    {
        theta = 0.0;
    }
    else
    {
        theta = norm.apply( vecW ) / tau;
    }

    c = 1.0 / sqrt( 1.0 + theta * theta );
    tau = tau * theta * c;
    eta = c * c * alpha;
    solution = solution + eta * vecD;

    if ( ( iteration % 2 ) == 1 )
    {
        iterationOdd();
    }
}

SolverPtr TFQMR::copy()
{
    return SolverPtr( new TFQMR( *this ) );
}

TFQMR::TFQMRRuntime& TFQMR::getRuntime()
{
    return mTFQMRRuntime;
}

const TFQMR::TFQMRRuntime& TFQMR::getConstRuntime() const
{
    return mTFQMRRuntime;
}
void TFQMR::writeAt( std::ostream& stream ) const
{
    stream << "TFQMR ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

std::string TFQMR::createValue()
{
    return "TFQMR";
}

Solver* TFQMR::create( const std::string name )
{
    return new TFQMR( name );
}

} /* end namespace solver */

} /* end namespace scai */
