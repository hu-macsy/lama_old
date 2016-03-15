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
#include <scai/solver/CGS.hpp>

// local library
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

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

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

void CGS::initialize( const Matrix& coefficients )
{
    SCAI_LOG_DEBUG(logger, "Initialization started for coefficients = "<< coefficients)

    IterativeSolver::initialize( coefficients );
 	CGSRuntime& runtime = getRuntime();

    runtime.mNormRes = 1.0;
    // runtime.mEps = std::numeric_limits<double>::epsilon()*3;            //CAREFUL: No abstract type
    switch(coefficients.getValueType()){
        case common::scalar::FLOAT:
            runtime.mEps = std::numeric_limits<float>::epsilon()*3;
            break;
        case common::scalar::DOUBLE:
            runtime.mEps = std::numeric_limits<double>::epsilon()*3;
            break;
        case common::scalar::LONG_DOUBLE:
            runtime.mEps = std::numeric_limits<long double>::epsilon()*3;
            break;
        case common::scalar::COMPLEX:
            runtime.mEps = std::numeric_limits<float>::epsilon()*3;
            break;
        case common::scalar::DOUBLE_COMPLEX:
            runtime.mEps = std::numeric_limits<double>::epsilon()*3;
            break;
        case common::scalar::LONG_DOUBLE_COMPLEX:
            runtime.mEps = std::numeric_limits<long double>::epsilon()*3;
            break;    
        default:
        SCAI_LOG_INFO(logger,"Valuetype not supported");
        break;
    }
    
    runtime.mRes0.reset( coefficients.newDenseVector() );
    runtime.mVecT.reset( coefficients.newDenseVector() );
    runtime.mVecP.reset( coefficients.newDenseVector() );
    runtime.mVecQ.reset( coefficients.newDenseVector() );
    runtime.mVecU.reset( coefficients.newDenseVector() );
    runtime.mVecPT.reset( coefficients.newDenseVector() );
    runtime.mVecUT.reset( coefficients.newDenseVector() );
    runtime.mVecTemp.reset( coefficients.newDenseVector() );
}


void CGS::solveInit( Vector& solution, const Vector& rhs )
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
    if(mPreconditioner != NULL){
        *runtime.mVecPT  = Scalar(0.0);
        mPreconditioner->solve( *runtime.mVecPT, *runtime.mVecP);      
    } 
    else   *runtime.mVecPT = *runtime.mVecP;



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
    Vector& vecPT = *runtime.mVecPT;
    Vector& vecUT = *runtime.mVecUT;
    Vector& vecTemp = *runtime.mVecTemp;
    Vector& solution = *runtime.mSolution;
    Scalar& innerProdRes = runtime.mInnerProdRes;

    Scalar alpha;
    Scalar beta;

    const Scalar& eps = runtime.mEps;
    Scalar& normRes = runtime.mNormRes;
	lama::L2Norm norm;



    vecT= A *vecPT;

    Scalar innerProduct = res0.dotProduct(vecT);
    if(normRes< eps || innerProduct < eps)    //innerProduct is small
        alpha=0.0;
    else alpha= innerProdRes/innerProduct;

    vecQ= vecU - alpha*vecT;

    // PRECONDITIONING
    if(mPreconditioner != NULL){
        vecUT = Scalar(0.0);
        vecTemp = vecU + vecQ;
        mPreconditioner->solve(vecUT,vecTemp);      
    } 
    else   vecUT  = vecU + vecQ;
    solution = solution + alpha*vecUT;

    Scalar innerProdResOld = innerProdRes;

    res = res - alpha*A*vecUT; 
    innerProdRes = res0.dotProduct(res);

    normRes = norm.apply(res);

    if(normRes < eps || innerProdResOld < eps)               // innerProdResOld is small
        beta=0.0;
    else beta = innerProdRes/ innerProdResOld ;

    vecU = res + beta*vecQ;
    vecP = vecU + beta*beta*vecP;
    vecP = vecP + beta*vecQ;

    // PRECONDITIONING
    if(mPreconditioner != NULL){
        vecPT  = Scalar(0.0);
        mPreconditioner->solve( vecPT, vecP);      
    } 
    else   vecPT = vecP ;

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

} /* end namespace solver */

} /* end namespace scai */
