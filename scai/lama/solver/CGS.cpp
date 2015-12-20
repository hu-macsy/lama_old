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
#include <scai/lama/solver/CGS.hpp>

// local library
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/norm/MaxNorm.hpp>

#include <scai/lama/DenseVector.hpp>

// std
#include <limits>

namespace scai
{

namespace lama{ 

SCAI_LOG_DEF_LOGGER( CGS::logger, "Solver.CGS" )

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
    SCAI_LOG_DEBUG(logger, "Initialization started for coefficients = "<< coefficients)

    Solver::initialize( coefficients );
 	CGSRuntime& runtime = getRuntime();

    runtime.mNormRes = 1.0;
    runtime.mEps = std::numeric_limits<double>::epsilon()*3;            //CAREFUL: No abstract type

    common::scalar::ScalarType type = coefficients.getValueType();

    runtime.mRes0.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecT.reset(Vector::createVector( type, coefficients.getDistributionPtr() ));
    runtime.mVecP.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecQ.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecU.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );


    runtime.mRes0->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecP->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecQ->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecU->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecT->setContextPtr( coefficients.getContextPtr() );
}


void CGS::solveInit( Vector& solution, const Vector& rhs ){
    CGSRuntime& runtime = getRuntime();

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

    Vector* initialR = (*runtime.mResidual).copy();
    runtime.mRes0.reset( initialR );

    Vector* vecP = (*runtime.mResidual).copy();
    runtime.mVecP.reset(vecP);

    Vector* vecU = (*runtime.mResidual).copy();
    runtime.mVecU.reset(vecU);

    //initial <res,res> inner product;
    runtime.mInnerProdRes = (*runtime.mRes0).dotProduct(*runtime.mRes0);

    runtime.mSolveInit = true;
}

void CGS::iterate(){
    CGSRuntime& runtime	= getRuntime();
    
    const Matrix& A = *runtime.mCoefficients;
    
    const Vector& res0 = *runtime.mRes0;
    Vector& res = *runtime.mResidual;
    Vector& vecP = *runtime.mVecP;
    Vector& vecQ = *runtime.mVecQ;
    Vector& vecU = *runtime.mVecU;
    Vector& vecT = *runtime.mVecT;
    Vector& solution = *runtime.mSolution;
    Scalar& innerProdRes = runtime.mInnerProdRes;

    Scalar alpha;
    Scalar beta;

    const Scalar& eps = runtime.mEps;
    Scalar& normRes = runtime.mNormRes;
	MaxNorm norm;



    vecT= A *vecP;         

    if(normRes< eps)    //residual is small
        alpha=0.0;
    else alpha= innerProdRes/res0.dotProduct(vecT);

    vecQ= vecU - alpha*vecT;
    solution = solution + alpha*vecU;
    solution = solution + alpha*vecQ;

    Scalar innerProdResOld = innerProdRes;

    res = res - alpha*A*vecU;
    res = res - alpha*A*vecQ; 
    innerProdRes = res0.dotProduct(res);

    normRes = norm.apply(res);

    if(normRes < eps)               // residual is small
        beta=0.0;
    else beta = innerProdRes/ innerProdResOld ;

    vecU = res + beta*vecQ;
    vecP = vecU + beta*beta*vecP;
    vecP = vecP + beta*vecQ;

    //End Implementation
    mCGSRuntime.mSolution.setDirty( false );
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

} /* end namespace lama */

} /* end namespace scai */
