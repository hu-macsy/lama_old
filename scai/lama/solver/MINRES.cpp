/**
 * @file MINRES.cpp
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
 * @brief MINRES.cpp
 * @author David schissler
 * @date 13.05.2015
 * @since
 */

// hpp
#include <scai/lama/solver/MINRES.hpp>
// others
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/DenseVector.hpp>

#include <limits>

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( MINRES::logger, "Solver.MINRES" )

MINRES::MINRES( const std::string& id )
    : IterativeSolver(id){}


MINRES::MINRES( const std::string& id, LoggerPtr logger )
    : IterativeSolver(id ,logger){}

MINRES::MINRES( const MINRES& other )
    : IterativeSolver( other ){}



MINRES::MINRESRuntime::MINRESRuntime()
    : IterativeSolverRuntime(){}

MINRES::~MINRES(){}

MINRES::MINRESRuntime::~MINRESRuntime(){}

void MINRES::initialize( const Matrix& coefficients ){
    SCAI_LOG_DEBUG(logger, "Initialization started for coefficients = "<< coefficients)

    Solver::initialize( coefficients );
 	MINRESRuntime& runtime = getRuntime();

    runtime.mBetaNew = 0.0;
    runtime.mC = 1.0;
    runtime.mCNew = 1.0;
    runtime.mS = 0.0;
    runtime.mSNew = 0.0;

    ScalarType type = coefficients.getValueType();
    
    runtime.mVecV.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecVOld.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecVNew.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecP.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecPOld.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecPNew.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );

    runtime.mVecV->setContext( coefficients.getContextPtr() );   
    runtime.mVecVOld->setContext( coefficients.getContextPtr() );
    runtime.mVecVNew->setContext( coefficients.getContextPtr() );
    runtime.mVecP->setContext( coefficients.getContextPtr() );
    runtime.mVecPOld->setContext( coefficients.getContextPtr() );
    runtime.mVecPNew->setContext(coefficients.getContextPtr());
}

void MINRES::solveInit( Vector& solution, const Vector& rhs ){
    MINRESRuntime& runtime = getRuntime();

    runtime.mRhs = &rhs;
    runtime.mSolution = &solution;

    if ( runtime.mCoefficients->getNumRows() != runtime.mRhs->size() ){
        COMMON_THROWEXCEPTION(
            "Size of rhs vector " << *runtime.mRhs << " does not match column size of matrix " << *runtime.mCoefficients );
    }

    if ( runtime.mCoefficients->getNumColumns() != solution.size() ){
        COMMON_THROWEXCEPTION(
            "Size of solution vector " << solution << " does not match row size of matrix " << *runtime.mCoefficients );
    }

    if ( runtime.mCoefficients->getColDistribution() != solution.getDistribution() ){
        COMMON_THROWEXCEPTION(
            "Distribution of lhs " << solution << " = " << solution.getDistribution() << " does not match (row) distribution of " << *runtime.mCoefficients << " = " << runtime.mCoefficients->getColDistribution() );
    }

    if ( runtime.mCoefficients->getDistribution() != runtime.mRhs->getDistribution() ){
        COMMON_THROWEXCEPTION(
            "Distribution of old Solution " << *runtime.mRhs << " = " << runtime.mRhs->getDistribution() << " does not match (row) distribution of " << *runtime.mCoefficients << " = " << runtime.mCoefficients->getDistribution() );
    }


    // Initialize


    this->getResidual();   

    Vector* vecVNew = (*runtime.mResidual).copy();
    runtime.mVecVNew.reset(vecVNew);    

    *runtime.mVecV *= (0.0);
    *runtime.mVecP *= (0.0);
    *runtime.mVecPNew *= (0.0);

    L2Norm norm;
    runtime.mZeta = norm.apply(*runtime.mResidual);
    *runtime.mVecVNew /= runtime.mZeta;

  
    runtime.mSolveInit = true;
}

void MINRES::Lanczos(){
    MINRESRuntime& runtime = getRuntime();
    const Matrix& A = *runtime.mCoefficients;
    Vector& vecV = *runtime.mVecV;
    Vector& vecVOld = *runtime.mVecVOld;
    Vector& vecVNew = *runtime.mVecVNew;
    Scalar& alpha = runtime.mAlpha;
    Scalar& betaNew = runtime.mBetaNew;
    Scalar& beta = runtime.mBeta;
    
    L2Norm norm;
    
    beta=betaNew;

    vecVOld.swap(vecV);
    vecV.swap(vecVNew);
    vecVNew = A*vecV - beta*vecVOld;
    alpha = vecVNew.dotProduct(vecV);
    vecVNew = vecVNew - alpha*vecV;
    betaNew = norm.apply(vecVNew);
    vecVNew = vecVNew/betaNew;

}	

void MINRES::applyGivensRotation(){
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
    //Old Givens-rotation
    cOld = c;
    c=cNew;
    sOld = s;
    s=sNew;

    rho1 = sOld*beta;
    rho2 = c*cOld*beta+s*alpha;
    rho3 = c*alpha-s*cOld*beta;

    Scalar tau, nu;
    //New Givens-rotation
    tau = abs(rho3)+betaNew;
    nu = tau * sqrt((rho3/tau)*(rho3/tau) + (betaNew/tau)*(betaNew/tau));
    cNew = rho3 / nu;
    sNew = betaNew / nu;
    rho3 = nu;

    Vector& vecP = *runtime.mVecP;
    Vector& vecPOld = *runtime.mVecPOld;
    Vector& vecPNew = *runtime.mVecPNew;
    const Vector& vecV = *runtime.mVecV;
    //Update P     
    vecPOld.swap(vecP);
    vecP.swap(vecPNew);
    vecPNew = vecV-rho1*vecPOld;
    vecPNew = vecPNew -rho2*vecP;
    vecPNew = vecPNew/rho3;

}

void MINRES::iterate(){
    MINRESRuntime& runtime = getRuntime();
    
    const Vector& vecPNew = *runtime.mVecPNew;
    Vector& solution = *runtime.mSolution;
    Scalar& cNew = runtime.mCNew;
    Scalar& sNew = runtime.mSNew;
    Scalar& zeta = runtime.mZeta;

    Lanczos();
    applyGivensRotation();    
    
    //New approximation
    solution = solution + cNew*zeta*vecPNew;
    zeta = -sNew*zeta;

    //MINRES Implementation End
}

SolverPtr MINRES::copy(){
    return SolverPtr( new MINRES( *this ) );
}

MINRES::MINRESRuntime& MINRES::getRuntime(){
    return mMINRESRuntime;
}

const MINRES::MINRESRuntime& MINRES::getConstRuntime() const{
    return mMINRESRuntime;
}

} /* end namespace lama */

} /* end namespace scai */
