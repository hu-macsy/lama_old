/**
 * @file MINRES.cpp
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
 * @brief Implementation of methods for the MINRES solver.
 * @author David Schissler
 * @date 13.05.2015
 */

// hpp
#include <scai/solver/MINRES.hpp>

// internal scai libraries
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/common/ScalarType.hpp>

// std
#include <limits>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( MINRES::logger, "Solver.MINRES" )

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

MINRES::MINRES( const std::string& id )
    : IterativeSolver( id ) {}


MINRES::MINRES( const std::string& id, LoggerPtr logger )
    : IterativeSolver( id , logger ) {}

MINRES::MINRES( const MINRES& other )
    : IterativeSolver( other ) {}



MINRES::MINRESRuntime::MINRESRuntime()
    : IterativeSolverRuntime() {}

MINRES::~MINRES() {}

MINRES::MINRESRuntime::~MINRESRuntime() {}

void MINRES::initialize( const Matrix& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )
    IterativeSolver::initialize( coefficients );
    MINRESRuntime& runtime = getRuntime();
    runtime.mBetaNew = 0.0;
    runtime.mC = 1.0;
    runtime.mCNew = 1.0;
    runtime.mS = 0.0;
    runtime.mSNew = 0.0;
    runtime.mVecV.reset( coefficients.newDenseVector() );
    runtime.mVecVOld.reset( coefficients.newDenseVector() );
    runtime.mVecVNew.reset( coefficients.newDenseVector() );
    runtime.mVecP.reset( coefficients.newDenseVector() );
    runtime.mVecPOld.reset( coefficients.newDenseVector() );
    runtime.mVecPNew.reset( coefficients.newDenseVector() );
    runtime.mEps = Scalar::eps1( coefficients.getValueType() ) * 3.0;
}

void MINRES::solveInit( Vector& solution, const Vector& rhs )
{
    MINRESRuntime& runtime = getRuntime();
    runtime.mRhs = &rhs;
    runtime.mSolution = &solution;
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumRows(), rhs.size(), "mismatch: #rows of matrix, rhs" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumColumns(), solution.size(), "mismatch: #cols of matrix, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getColDistribution(), solution.getDistribution(), "mismatch: matrix col dist, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getRowDistribution(), rhs.getDistribution(), "mismatch: matrix row dist, rhs dist" )
    // Initialize
    this->getResidual();
    *runtime.mVecVNew = *runtime.mResidual;
    *runtime.mVecV *= Scalar( 0.0 );
    *runtime.mVecP *= Scalar( 0.0 );
    *runtime.mVecPNew *= Scalar( 0.0 );
    lama::L2Norm norm;
    runtime.mZeta = norm.apply( *runtime.mResidual );
    *runtime.mVecVNew /= runtime.mZeta;
    runtime.mSolveInit = true;
}

void MINRES::Lanczos()
{
    MINRESRuntime& runtime = getRuntime();
    const Matrix& A = *runtime.mCoefficients;
    Vector& vecV = *runtime.mVecV;
    Vector& vecVOld = *runtime.mVecVOld;
    Vector& vecVNew = *runtime.mVecVNew;
    Scalar& alpha = runtime.mAlpha;
    Scalar& betaNew = runtime.mBetaNew;
    Scalar& beta = runtime.mBeta;
    Scalar& eps = runtime.mEps;
    lama::L2Norm norm;
    beta = betaNew;
    vecVOld.swap( vecV );
    vecV.swap( vecVNew );
    vecVNew = A * vecV - beta * vecVOld;
    alpha = vecV.dotProduct( vecVNew );
    vecVNew = vecVNew - alpha * vecV;
    betaNew = norm.apply( vecVNew );

    if ( abs( betaNew ) < eps || 1.0 / abs( betaNew ) < eps )
    {
        vecVNew = Scalar( 0.0 );
    }
    else
    {
        vecVNew = vecVNew / betaNew;
    }
}

void MINRES::applyGivensRotation()
{
    MINRESRuntime& runtime = getRuntime();
    Scalar& c = runtime.mC;
    Scalar& cOld = runtime.mCOld;
    Scalar& cNew = runtime.mCNew;
    Scalar& s = runtime.mS;
    Scalar& sOld = runtime.mSOld;
    Scalar& sNew = runtime.mSNew;
    Scalar& alpha = runtime.mAlpha;
    Scalar& beta = runtime.mBeta;
    Scalar& betaNew = runtime.mBetaNew;
    Scalar rho1, rho2, rho3;
    Scalar& eps = runtime.mEps;
    //Old Givens-rotation
    cOld = c;
    c = cNew;
    sOld = s;
    s = sNew;
    rho1 = sOld * beta;
    rho2 = c * cOld * beta + s * alpha;
    rho3 = c * alpha - s * cOld * beta;
    Scalar tau, nu;
    //New Givens-rotation
    tau = abs( rho3 ) + betaNew;

    if ( abs( tau ) < eps || 1.0 / abs( tau ) < eps )
    {
        nu = 1.0;
        rho3 = 1.0;
        cNew = 0.0;
        sNew = 0.0;
    }
    else
    {
        nu = tau * sqrt( ( rho3 / tau ) * ( rho3 / tau ) + ( betaNew / tau ) * ( betaNew / tau ) );
        cNew = rho3 / nu;
        sNew = betaNew / nu;
        rho3 = nu;
    }

    Vector& vecP = *runtime.mVecP;
    Vector& vecPOld = *runtime.mVecPOld;
    Vector& vecPNew = *runtime.mVecPNew;
    const Vector& vecV = *runtime.mVecV;
    //Update P
    vecPOld.swap( vecP );
    vecP.swap( vecPNew );
    vecPNew = vecV - rho1 * vecPOld;
    vecPNew = vecPNew - rho2 * vecP;
    vecPNew = vecPNew / rho3;
}

void MINRES::iterate()
{
    MINRESRuntime& runtime = getRuntime();
    const Vector& vecPNew = *runtime.mVecPNew;
    Vector& solution = *runtime.mSolution;
    Scalar& cNew = runtime.mCNew;
    Scalar& sNew = runtime.mSNew;
    Scalar& zeta = runtime.mZeta;
    Lanczos();
    applyGivensRotation();
    //New approximation
    solution = solution + cNew * zeta * vecPNew;
    zeta = -sNew * zeta;
    //MINRES Implementation End
}

SolverPtr MINRES::copy()
{
    return SolverPtr( new MINRES( *this ) );
}

MINRES::MINRESRuntime& MINRES::getRuntime()
{
    return mMINRESRuntime;
}

const MINRES::MINRESRuntime& MINRES::getConstRuntime() const
{
    return mMINRESRuntime;
}

void MINRES::writeAt( std::ostream& stream ) const
{
    stream << "MINRES ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

std::string MINRES::createValue()
{
    return "MINRES";
}

Solver* MINRES::create( const std::string name )
{
    return new MINRES( name );
}

} /* end namespace solver */

} /* end namespace scai */
