/**
 * @file BiCGstab.cpp
 *
 * @license
 * Copyright (c) 2013
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
 * @brief BiCGstab.cpp
 * @author lschubert
 * @date 06.08.2013
 * @since 1.1.0
 */

// hpp
#include <scai/lama/solver/BiCGstab.hpp>

// local library
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

namespace lama
{

SCAI_LOG_DEF_LOGGER( BiCGstab::logger, "Solver.BiCGstab" )

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

    Solver::initialize( coefficients );
    BiCGstabRuntime& runtime = getRuntime();

    runtime.mAlpha  = 1.0;
    runtime.mOmega = 1.0;
    runtime.mRhoOld = 1.0;
    runtime.mResNorm = 1.0;
    runtime.mEps = std::numeric_limits<double>::epsilon() * 3;                  //CAREFUL: No abstract type

    common::scalar::ScalarType type = coefficients.getValueType();

    runtime.mRes0.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecV.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecP.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecS.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecT.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );


    runtime.mRes0->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecV->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecP->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecS->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecT->setContextPtr( coefficients.getContextPtr() );
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

    if ( runtime.mCoefficients->getColDistribution() != solution.getDistribution() )
    {
        COMMON_THROWEXCEPTION(
            "Distribution of lhs " << solution << " = " << solution.getDistribution() << " does not match (row) distribution of " << *runtime.mCoefficients << " = " << runtime.mCoefficients->getColDistribution() );
    }

    if ( runtime.mCoefficients->getDistribution() != runtime.mRhs->getDistribution() )
    {
        COMMON_THROWEXCEPTION(
            "Distribution of old Solution " << *runtime.mRhs << " = " << runtime.mRhs->getDistribution() << " does not match (row) distribution of " << *runtime.mCoefficients << " = " << runtime.mCoefficients->getDistribution() );
    }


    // Initialize
    this->getResidual();

    Vector* initialR = ( *runtime.mResidual ).copy();
    runtime.mRes0.reset( initialR );

    *runtime.mVecV = Scalar( 0 );
    *runtime.mVecP = Scalar( 0 );

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

    Scalar& alpha = runtime.mAlpha;
    Scalar& beta = runtime.mBeta;
    Scalar& omega = runtime.mOmega;
    Scalar& rhoOld = runtime.mRhoOld;
    Scalar& rhoNew = runtime.mRhoNew;

    const Scalar& eps = runtime.mEps;
    Scalar& resNorm = runtime.mResNorm;
    L2Norm norm;

    rhoNew = res0.dotProduct( res );

    if ( resNorm < eps || rhoOld < eps || omega < eps ) // residual is small
    {
        beta = 0.0;
    }
    else
    {
        std::cout << "rhoNew = " << rhoNew << ", rhoOld = " << rhoOld << ", alpha = " << alpha << ", omega = " << omega << std::endl;

        beta = rhoNew / rhoOld * ( alpha / omega );
    }

    SCAI_LOG_INFO( logger, "resNorm = " << resNorm << ", eps = " << eps << ", beta = " << beta )

    vecP = vecP - omega * vecV;
    vecP = res + beta * vecP;
    vecV = A * vecP;

    Scalar innerProd = res0.dotProduct( vecV );

    if ( innerProd < eps )
    {
        alpha = 0.0;
    }
    else
    {
        alpha = rhoNew / innerProd;
    }

    vecS = res - alpha * vecV;
    vecT = A * vecS;

    innerProd = vecT.dotProduct( vecT );

    if ( innerProd < eps ) // residual is small
    {
        omega = 0.0;
    }
    else
    {
        omega = vecT.dotProduct( vecS ) / innerProd;
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

} /* end namespace lama */

} /* end namespace scai */
