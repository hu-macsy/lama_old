/**
 * @file CGS.cpp
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
 * @brief CGS.cpp
 * @author David Schissler
 * @date 18.05.2015
 */

// hpp
#include <scai/solver/CGS.hpp>

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

namespace solver
{

SCAI_LOG_DEF_LOGGER( CGS::logger, "Solver.CGS" )

using lama::_Matrix;
using lama::_Vector;
using lama::Scalar;

CGS::CGS( const std::string& id )
    : IterativeSolver( id ) {}


CGS::CGS( const std::string& id, LoggerPtr logger )
    : IterativeSolver( id , logger ) {}

CGS::CGS( const CGS& other )
    : IterativeSolver( other ) {}



CGS::CGSRuntime::CGSRuntime()
    : IterativeSolverRuntime() {}

CGS::~CGS() {}

CGS::CGSRuntime::~CGSRuntime() {}

void CGS::initialize( const _Matrix& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )
    IterativeSolver::initialize( coefficients );
    CGSRuntime& runtime = getRuntime();
    runtime.mNormRes = 1.0;
    runtime.mEps = Scalar::eps1( coefficients.getValueType() ) * 3.0;
    dmemo::DistributionPtr rowDist = coefficients.getRowDistributionPtr();
    runtime.mRes0.reset( coefficients.newVector( rowDist ) );
    runtime.mVecT.reset( coefficients.newVector( rowDist ) );
    runtime.mVecP.reset( coefficients.newVector( rowDist ) );
    runtime.mVecQ.reset( coefficients.newVector( rowDist ) );
    runtime.mVecU.reset( coefficients.newVector( rowDist ) );
    runtime.mVecPT.reset( coefficients.newVector( rowDist ) );
    runtime.mVecUT.reset( coefficients.newVector( rowDist ) );
    runtime.mVecTemp.reset( coefficients.newVector( rowDist ) );
}


void CGS::solveInit( _Vector& solution, const _Vector& rhs )
{
    CGSRuntime& runtime = getRuntime();
    runtime.mRhs = &rhs;
    runtime.mSolution = &solution;
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumRows(), rhs.size(), "mismatch: #rows of matrix, rhs" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumColumns(), solution.size(), "mismatch: #cols of matrix, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getColDistribution(), solution.getDistribution(), "mismatch: matrix col dist, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getRowDistribution(), rhs.getDistribution(), "mismatch: matrix row dist, rhs dist" )
    // Initialize
    this->getResidual();
    *runtime.mRes0 = *runtime.mResidual;
    *runtime.mVecP = *runtime.mResidual;
    *runtime.mVecU = *runtime.mResidual;

    // PRECONDITIONING
    if ( mPreconditioner != NULL )
    {
        runtime.mVecPT->setSameValue( runtime.mVecP->getDistributionPtr(), 0 );
        mPreconditioner->solve( *runtime.mVecPT, *runtime.mVecP );
    }
    else
    {
        *runtime.mVecPT = *runtime.mVecP;
    }

    //initial <res,res> inner product;
    runtime.mInnerProdRes = ( *runtime.mRes0 )._dotProduct( *runtime.mRes0 );
    SCAI_LOG_INFO( logger, "solveInit, mInnerProdRes = " << runtime.mInnerProdRes )
    runtime.mSolveInit = true;
}

void CGS::iterate()
{
    CGSRuntime& runtime = getRuntime();
    const _Matrix& A = *runtime.mCoefficients;
    const _Vector& res0 = *runtime.mRes0;
    _Vector& res = *runtime.mResidual;
    _Vector& vecP = *runtime.mVecP;
    _Vector& vecQ = *runtime.mVecQ;
    _Vector& vecU = *runtime.mVecU;
    _Vector& vecT = *runtime.mVecT;
    _Vector& vecPT = *runtime.mVecPT;
    _Vector& vecUT = *runtime.mVecUT;
    _Vector& vecTemp = *runtime.mVecTemp;
    _Vector& solution = *runtime.mSolution;
    Scalar& innerProdRes = runtime.mInnerProdRes;
    Scalar alpha;
    Scalar beta;
    const Scalar& eps = runtime.mEps;
    Scalar& normRes = runtime.mNormRes;
    lama::L2Norm norm;
    vecT = A * vecPT;
    Scalar innerProduct = res0._dotProduct( vecT );

    if ( normRes < eps || innerProduct < eps ) //innerProduct is small
    {
        alpha = 0.0;
    }
    else
    {
        alpha = innerProdRes / innerProduct;
    }

    SCAI_LOG_INFO( logger, "vecQ = vecU - " << alpha << " vecT, normRes = " << normRes
                   << ", innerProdRes = " << innerProdRes << ", innerProduct = " << innerProduct )

    vecQ = vecU - alpha * vecT;

    // PRECONDITIONING
    if ( mPreconditioner != NULL )
    {
        // vecUT = 0;  here the more general approach
        vecUT.setSameValue( mPreconditioner->getCoefficients().getColDistributionPtr(), 0 );
        vecTemp = vecU + vecQ;
        mPreconditioner->solve( vecUT, vecTemp );
    }
    else
    {
        vecUT  = vecU + vecQ;
    }

    solution = solution + alpha * vecUT;
    Scalar innerProdResOld = innerProdRes;
    res = res - alpha * A * vecUT;
    innerProdRes = res0._dotProduct( res );
    normRes = norm.apply( res );

    if ( normRes < eps || innerProdResOld < eps )            // innerProdResOld is small
    {
        beta = 0.0;
    }
    else
    {
        beta = innerProdRes / innerProdResOld ;
    }

    SCAI_LOG_INFO( logger, "beta = " << beta )

    vecU = res + beta * vecQ;
    vecP = vecU + beta * beta * vecP;
    vecP = vecP + beta * vecQ;

    // PRECONDITIONING
    if ( mPreconditioner != NULL )
    {
        vecPT.setSameValue( mPreconditioner->getCoefficients().getColDistributionPtr(), 0 );
        mPreconditioner->solve( vecPT, vecP );
    }
    else
    {
        vecPT = vecP ;
    }

    //End Implementation
    mCGSRuntime.mSolution.setDirty( false );
}


SolverPtr CGS::copy()
{
    return SolverPtr( new CGS( *this ) );
}

CGS::CGSRuntime& CGS::getRuntime()
{
    return mCGSRuntime;
}

const CGS::CGSRuntime& CGS::getConstRuntime() const
{
    return mCGSRuntime;
}

std::string CGS::createValue()
{
    return "CGS";
}

Solver* CGS::create( const std::string name )
{
    return new CGS( name );
}

void CGS::writeAt( std::ostream& stream ) const
{
    stream << "CGS ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

} /* end namespace solver */

} /* end namespace scai */
