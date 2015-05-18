/**
 * @file CGS.cpp
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
 * @brief CGS.cpp
 * @author David schissler
 * @date 18.05.2015
 * @since
 */

// hpp
#include <lama/solver/CGS.hpp>
// others
#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixExpressions.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

#include <lama/norm/MaxNorm.hpp>

#include <lama/DenseVector.hpp>

#include <lama/solver/criteria/ResidualThreshold.hpp>


#include <limits>

namespace lama{ 

LAMA_LOG_DEF_LOGGER( CGS::logger, "Solver.CGS" )

CGS::CGS( const std::string& id )
    : IterativeSolver(id){}


CGS::CGS( const std::string& id, LoggerPtr logger )
    : IterativeSolver(id ,logger){}

CGS::CGS( const CGS& other )
    : IterativeSolver( other ){}



CGS::CGSRuntime::CGSRuntime()
    : IterativeSolverRuntime(){}

CGS::~CGS(){}

CGS::CGSRuntime::~CGSRuntime(){}

void CGS::initialize( const Matrix& coefficients ){
    LAMA_LOG_DEBUG(logger, "Initialization started for coefficients = "<< coefficients)

    Solver::initialize( coefficients );
 	CGSRuntime& runtime = getRuntime();

    runtime.mAlpha  = 1.0;


    Scalar::ScalarType type = coefficients.getValueType();

    runtime.mRes0.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mResOld.reset(Vector::createVector( type, coefficients.getDistributionPtr() ));
    runtime.mVecP.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecQ.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecU.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );


    runtime.mRes0->setContext( coefficients.getContextPtr() );   
    runtime.mVecP->setContext( coefficients.getContextPtr() );
    runtime.mVecQ->setContext( coefficients.getContextPtr() );
    runtime.mVecU->setContext( coefficients.getContextPtr() );
}


void CGS::solveInit( Vector& solution, const Vector& rhs ){
    CGSRuntime& runtime = getRuntime();

    runtime.mRhs = &rhs;
    runtime.mSolution = &solution;

    if ( runtime.mCoefficients->getNumRows() != runtime.mRhs->size() ){
        LAMA_THROWEXCEPTION(
            "Size of rhs vector " << *runtime.mRhs << " does not match column size of matrix " << *runtime.mCoefficients );
    }

    if ( runtime.mCoefficients->getNumColumns() != solution.size() ){
        LAMA_THROWEXCEPTION(
            "Size of solution vector " << solution << " does not match row size of matrix " << *runtime.mCoefficients );
    }

    if ( runtime.mCoefficients->getColDistribution() != solution.getDistribution() ){
        LAMA_THROWEXCEPTION(
            "Distribution of lhs " << solution << " = " << solution.getDistribution() << " does not match (row) distribution of " << *runtime.mCoefficients << " = " << runtime.mCoefficients->getColDistribution() );
    }

    if ( runtime.mCoefficients->getDistribution() != runtime.mRhs->getDistribution() ){
        LAMA_THROWEXCEPTION(
            "Distribution of old Solution " << *runtime.mRhs << " = " << runtime.mRhs->getDistribution() << " does not match (row) distribution of " << *runtime.mCoefficients << " = " << runtime.mCoefficients->getDistribution() );
    }


    // Initialize
    this->getResidual();   

    Vector* initialR = (*runtime.mResidual).copy();
    runtime.mRes0.reset( initialR );

    Vector* vecP = (*runtime.mResidual).copy();
    runtime.mVecP.reset(vecP);

    Vector* vecU = (*runtime.mResidual).copy();
    runtime.mVecU.reset(vecU);

    runtime.mSolveInit = true;
}

void CGS::iterate(){
    CGSRuntime& runtime	= getRuntime();
    
    const Matrix& A = *runtime.mCoefficients;
    
    const Vector& res0 = *runtime.mRes0;
    Vector& res = *runtime.mResidual;
    Vector& resOld = *runtime.mResOld;
    Vector& vecP = *runtime.mVecP;
    Vector& vecQ = *runtime.mVecQ;
    Vector& vecU = *runtime.mVecU;
    Vector& solution = *runtime.mSolution;
    Scalar& alpha = runtime.mAlpha;
    Scalar beta;
	
    if(this->getIterationCount() >0){
        resOld = res;
        res = res - alpha*A*vecU;
        res = res - alpha*A*vecQ;
        beta = res.dotProduct(res0)/resOld.dotProduct(res0);
        vecU = res + beta*vecQ;
        vecP = vecU + beta*beta*vecP;
        vecP = vecP + beta*vecQ;
    }

    resOld= A *vecP;                  // CAREFUL: We use resOld as a temporary Vector.

    alpha= res.dotProduct(res0)/res0.dotProduct(resOld);
    vecQ= vecU - alpha*A*vecP;
    solution = solution + alpha*vecU;
    solution = solution +alpha*vecQ;

    //End Implementation
    mCGSRuntime.mSolution.setDirty( false );
}

void CGS::setStoppingCriterion( const CriterionPtr criterion ){
    
    Scalar eps = std::numeric_limits<double>::epsilon()*2;                // NOT ABSTRACT.

    NormPtr norm = NormPtr( new MaxNorm() );
    CriterionPtr rt( new ResidualThreshold( norm, eps, ResidualThreshold::Relative ) );
    LAMA_ASSERT_ERROR( criterion, "Criterion defined is NULL." )
    LAMA_LOG_INFO( logger, "Criteria " << *criterion << " defined." )

    mCriterionRootComponent = ( criterion ||  rt );
}

SolverPtr CGS::copy(){
    return SolverPtr( new CGS( *this ) );
}

CGS::CGSRuntime& CGS::getRuntime(){
    return mCGSRuntime;
}

const CGS::CGSRuntime& CGS::getConstRuntime() const{
    return mCGSRuntime;
}

} /* namespace lama */