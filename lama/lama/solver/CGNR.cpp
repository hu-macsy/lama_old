/**
 * @file CGNR.cpp
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
 * @brief CGNR.cpp
 * @author David schissler
 * @date 27.05.2015
 * @since
 */

// hpp
#include <lama/solver/CGNR.hpp>
// others
#include <lama/expression/VectorExpressions.hpp>
#include <lama/expression/MatrixExpressions.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

#include <lama/norm/L2Norm.hpp>

#include <lama/DenseVector.hpp>

#include <limits>

namespace lama{ 

LAMA_LOG_DEF_LOGGER( CGNR::logger, "Solver.CGNR" )

CGNR::CGNR( const std::string& id )
    : IterativeSolver(id){}


CGNR::CGNR( const std::string& id, LoggerPtr logger )
    : IterativeSolver(id ,logger){}

CGNR::CGNR( const CGNR& other )
    : IterativeSolver( other ){}



CGNR::CGNRRuntime::CGNRRuntime()
    : IterativeSolverRuntime(){}

CGNR::~CGNR(){}

CGNR::CGNRRuntime::~CGNRRuntime(){}


void CGNR::initialize( const Matrix& coefficients ){
    LAMA_LOG_DEBUG(logger, "Initialization started for coefficients = "<< coefficients)

    Solver::initialize(coefficients);
 	CGNRRuntime& runtime = getRuntime();

    common::ScalarType type = coefficients.getValueType();
    runtime.mEps = std::numeric_limits<double>::epsilon()*3;   //CAREFUL: No abstract type 


    runtime.mVecD.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecW.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecZ.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );

    runtime.mVecD->setContext( coefficients.getContextPtr() );   
    runtime.mVecW->setContext( coefficients.getContextPtr() );
    runtime.mVecZ->setContext( coefficients.getContextPtr() );
    

}


void CGNR::solveInit( Vector& solution, const Vector& rhs ){
    CGNRRuntime& runtime = getRuntime();

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
   
    const Matrix& A = *runtime.mCoefficients;
    const Vector& residual = *runtime.mResidual;
    Vector& vecZ = *runtime.mVecZ;
    vecZ = residual*A;

    Vector* vecD = (*runtime.mVecZ).copy();
    runtime.mVecD.reset(vecD);
  
    L2Norm norm;
    runtime.mNormVecZ = norm.apply(vecZ);
    runtime.mSolveInit = true;
}

void CGNR::iterate(){
    CGNRRuntime& runtime = getRuntime();
   
    const Matrix& A = *runtime.mCoefficients;
    Vector& vecW = *runtime.mVecW;
    Vector& vecD = *runtime.mVecD;
    Vector& vecZ = *runtime.mVecZ;
    Vector& residual = *runtime.mResidual;
	Vector& solution = *runtime.mSolution;
    Scalar& normVecZ = runtime.mNormVecZ;
    Scalar alpha;
    Scalar beta;

    L2Norm norm;
    const Scalar& eps = runtime.mEps;    

    vecW= A*vecD;
    Scalar normVecW = norm.apply(vecW);

    if(normVecW< eps)               //norm is small 
        alpha=0.0;
    else alpha = (normVecZ*normVecZ)/(normVecW*normVecW);
    
    solution= solution + alpha*vecD;
    residual = residual - alpha*vecW;
    vecZ = residual*A;
    
    Scalar normVecZNew=norm.apply(vecZ);


    if(normVecZNew < eps)           //norm is small
        beta=0.0;
    else beta = (normVecZNew*normVecZNew)/(normVecZ*normVecZ);

    vecD = vecZ + beta*vecD;
    normVecZ = normVecZNew;
    //CGNR Implementation End
    mCGNRRuntime.mSolution.setDirty(false);
}

SolverPtr CGNR::copy(){
    return SolverPtr( new CGNR( *this ) );
}

CGNR::CGNRRuntime& CGNR::getRuntime(){
    return mCGNRRuntime;
}

const CGNR::CGNRRuntime& CGNR::getConstRuntime() const{
    return mCGNRRuntime;
}

} /* namespace lama */

