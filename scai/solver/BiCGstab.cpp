/**
 * @file BiCGstab.cpp
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
 * @brief BiCGstab.cpp
 * @author lschubert
 * @date 06.08.2013
 */

// hpp
#include <scai/solver/BiCGstab.hpp>

// internal scai libraries
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/norm/MaxNorm.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/DenseVector.hpp>

// std
#include <limits>

namespace scai
{

using namespace lama;

namespace solver
{

SCAI_LOG_DEF_LOGGER( BiCGstab::logger, "Solver.BiCGstab" )

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

BiCGstab::BiCGstab( const std::string& id )
    : IterativeSolver( id ) {}


BiCGstab::BiCGstab( const std::string& id, LoggerPtr logger )
    : IterativeSolver( id , logger ) {}

BiCGstab::BiCGstab( const BiCGstab& other )
    : IterativeSolver( other ) {}



BiCGstab::BiCGstabRuntime::BiCGstabRuntime()
    : IterativeSolverRuntime() {}

BiCGstab::~BiCGstab() {}

BiCGstab::BiCGstabRuntime::~BiCGstabRuntime() {}

void BiCGstab::initialize( const Matrix& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )
    IterativeSolver::initialize( coefficients );
    BiCGstabRuntime& runtime = getRuntime();
    runtime.mAlpha  = 1.0;
    runtime.mOmega = 1.0;
    runtime.mRhoOld = 1.0;
    runtime.mResNorm = 1.0;
    runtime.mEps = Scalar::eps1( coefficients.getValueType() ) * 3.0;
    // get runtime vectors with same row distribution / context / type as cofficients matrix
    dmemo::DistributionPtr rowDist = coefficients.getRowDistributionPtr();
    runtime.mRes0.reset( coefficients.newVector( rowDist ) );
    runtime.mVecV.reset( coefficients.newVector( rowDist ) );
    runtime.mVecP.reset( coefficients.newVector( rowDist ) );
    runtime.mVecS.reset( coefficients.newVector( rowDist ) );
    runtime.mVecT.reset( coefficients.newVector( rowDist ) );
    runtime.mVecPT.reset( coefficients.newVector( rowDist ) );
    runtime.mVecST.reset( coefficients.newVector( rowDist ) );
    runtime.mVecTT.reset( coefficients.newVector( rowDist ) );
}

void BiCGstab::solveInit( Vector& solution, const Vector& rhs )
{
    BiCGstabRuntime& runtime = getRuntime();
    runtime.mRhs = &rhs;
    runtime.mSolution = &solution;

    if ( runtime.mCoefficients->getNumRows() != runtime.mRhs->size() )
    {
        COMMON_THROWEXCEPTION(
            "Size of rhs vector " << *runtime.mRhs << " does not match column size of matrix " << *runtime.mCoefficients );
    }

    if ( runtime.mCoefficients->getNumColumns() != solution.size() )
    {
        COMMON_THROWEXCEPTION(
            "Size of solution vector " << solution << " does not match row size of matrix " << *runtime.mCoefficients );
    }

    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getColDistribution(), solution.getDistribution(), "distribution mismatch" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getRowDistribution(), runtime.mRhs->getDistribution(), "distribution mismatch" )
    // Initialize
    this->getResidual();

    *runtime.mRes0 = *runtime.mResidual;
    runtime.mVecV->setSameValue( rhs.getDistributionPtr(), 0 );
    runtime.mVecP->setSameValue( rhs.getDistributionPtr(), 0 );
    runtime.mSolveInit = true;
}

void BiCGstab::iterate()
{
    BiCGstabRuntime& runtime    = getRuntime();
    const Matrix& A = *runtime.mCoefficients;
    const Vector& res0 = *runtime.mRes0;
    Vector& res = *runtime.mResidual;
    Vector& vecV = *runtime.mVecV;
    Vector& vecP = *runtime.mVecP;
    Vector& vecS = *runtime.mVecS;
    Vector& vecT = *runtime.mVecT;
    Vector& solution = *runtime.mSolution;
    Vector& vecPT = *runtime.mVecPT;
    Vector& vecST = *runtime.mVecST;
    Vector& vecTT = *runtime.mVecTT;
    Scalar& alpha = runtime.mAlpha;
    Scalar& beta = runtime.mBeta;
    Scalar& omega = runtime.mOmega;
    Scalar& rhoOld = runtime.mRhoOld;
    Scalar& rhoNew = runtime.mRhoNew;
    const Scalar& eps = runtime.mEps;
    Scalar& resNorm = runtime.mResNorm;
    lama::L2Norm norm;
    rhoNew = res0.dotProduct( res );

    if ( resNorm < eps || rhoOld < eps || omega < eps ) // scalars are small
    {
        beta = 0.0;
    }
    else
    {
        beta = rhoNew / rhoOld * ( alpha / omega );
    }

    SCAI_LOG_INFO( logger, "resNorm = " << resNorm << ", eps = " << eps << ", beta = " << beta )
    vecP = vecP - omega * vecV;
    vecP = res + beta * vecP;

    // PRECONDITIONING
    if ( mPreconditioner != NULL )
    {
        vecPT.setSameValue( A.getRowDistributionPtr(), 0 );
        mPreconditioner->solve( vecPT, vecP );
    }
    else
    {
        vecPT = vecP;
    }

    vecV = A * vecPT;
    Scalar innerProd = res0.dotProduct( vecV );

    if ( resNorm < eps || innerProd < eps ) // scalar is small
    {
        alpha = 0.0;
    }
    else
    {
        alpha = rhoNew / innerProd;
    }

    vecS = res - alpha * vecV;

    // PRECONDITIONING
    if ( mPreconditioner != NULL )
    {
        vecST.setSameValue( A.getColDistributionPtr(), 0 );
        mPreconditioner->solve( vecST, vecS );
        vecT = A * vecST;
        vecTT.setSameValue( A.getColDistributionPtr(), 0 );
        mPreconditioner->solve( vecTT, vecT );
    }
    else
    {
        vecST = vecS;
        vecT = A * vecST;
        vecTT = vecT;
    }

    innerProd = vecTT.dotProduct( vecTT );

    if ( resNorm < eps || innerProd < eps ) //scalar is small
    {
        omega = 0.0;
    }
    else
    {
        omega = vecTT.dotProduct( vecST ) / innerProd;
    }

    solution = solution + alpha * vecP;
    solution = solution + omega * vecS;
    res = vecS - omega * vecT;
    rhoOld = rhoNew;
    resNorm = norm.apply( res );
    //BiCGStab implementation end
    mBiCGstabRuntime.mSolution.setDirty( false );
}

SolverPtr BiCGstab::copy()
{
    return SolverPtr( new BiCGstab( *this ) );
}

BiCGstab::BiCGstabRuntime& BiCGstab::getRuntime()
{
    return mBiCGstabRuntime;
}

const BiCGstab::BiCGstabRuntime& BiCGstab::getConstRuntime() const
{
    return mBiCGstabRuntime;
}

std::string BiCGstab::createValue()
{
    return "BiCGstab";
}

Solver* BiCGstab::create( const std::string name )
{
    return new BiCGstab( name );
}

void BiCGstab::writeAt( std::ostream& stream ) const
{
    stream << "BiCGstab ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

} /* end namespace sovler */

} /* end namespace scai */
