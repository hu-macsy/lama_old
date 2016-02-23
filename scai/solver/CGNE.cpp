/**
 * @file CGNE.cpp
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
 * @brief CGNE.cpp
 * @author David schissler
 * @date 27.05.2015
 * @since
 */

// hpp
#include <scai/solver/CGNE.hpp>

// local library
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/common/Constants.hpp>
#include <scai/lama/DenseVector.hpp>

// std 
#include <limits>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( CGNE::logger, "Solver.CGNE" )

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

CGNE::CGNE( const std::string& id )
    : IterativeSolver(id){}


CGNE::CGNE( const std::string& id, LoggerPtr logger )
    : IterativeSolver(id ,logger){}

CGNE::CGNE( const CGNE& other )
    : IterativeSolver( other ){}

CGNE::CGNERuntime::CGNERuntime()
    : IterativeSolverRuntime(){}

CGNE::~CGNE(){}

CGNE::CGNERuntime::~CGNERuntime(){}


void CGNE::initialize( const Matrix& coefficients ){
    SCAI_LOG_DEBUG(logger, "Initialization started for coefficients = "<< coefficients)

    IterativeSolver::initialize(coefficients);
    CGNERuntime& runtime = getRuntime();

    common::scalar::ScalarType type = coefficients.getValueType();
    runtime.mEps = std::numeric_limits<double>::epsilon()*3;                  //CAREFUL: No abstract type

    runtime.mTransposedMat.reset( coefficients.newMatrix() );
    runtime.mTransposedMat->assignTranspose( coefficients );
    runtime.mTransposedMat->conj();
    
    runtime.mVecP.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecZ.reset( Vector::createVector( type, coefficients.getDistributionPtr() ) );
    
    runtime.mVecP->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecZ->setContextPtr( coefficients.getContextPtr() );
}


void CGNE::solveInit( Vector& solution, const Vector& rhs ){
    CGNERuntime& runtime = getRuntime();

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
    // PRECONDITIONING
    if(mPreconditioner != NULL){
        *runtime.mVecZ = Scalar(0.0);
        mPreconditioner->solve(*runtime.mVecZ,*runtime.mResidual);
    }
    else *runtime.mVecZ = *runtime.mResidual;

    *runtime.mVecP = (*runtime.mTransposedMat) * (*runtime.mVecZ);
    runtime.mSolveInit = true;
}

void CGNE::iterate(){
    CGNERuntime& runtime = getRuntime();
   
    const Matrix& A = *runtime.mCoefficients;
    const Matrix& transposedA = *runtime.mTransposedMat;

    Vector& vecP= *runtime.mVecP;
    Vector& residual = *runtime.mResidual;
    Vector& solution = *runtime.mSolution;
    Vector& vecZ = *runtime.mVecZ;
    Scalar alpha;
    Scalar beta;
    Scalar eps = runtime.mEps;

    Scalar scalarProductP = vecP.dotProduct(vecP);
    Scalar scalarProductZR = vecZ.dotProduct(residual);
    if(scalarProductP < eps)    alpha=0.0;     //norm is small 
    else    alpha = scalarProductZR/scalarProductP;
    
    solution= solution + alpha*vecP;
    residual = residual - alpha*A*vecP;

    // PRECONDITIONING
    if(mPreconditioner != NULL) mPreconditioner->solve(vecZ,residual);
    else vecZ = residual;

    if(scalarProductZR < eps)     beta=0.0;   //norm is small
    else    beta = vecZ.dotProduct(residual)/scalarProductZR;

    vecP = transposedA*vecZ + beta * vecP;
    //CGNE Implementation End
    mCGNERuntime.mSolution.setDirty(false);
}

SolverPtr CGNE::copy(){
    return SolverPtr( new CGNE( *this ) );
}

CGNE::CGNERuntime& CGNE::getRuntime(){
    return mCGNERuntime;
}

const CGNE::CGNERuntime& CGNE::getConstRuntime() const{
    return mCGNERuntime;
}

std::string CGNE::createValue()
{
    return "CGNE";
}

Solver* CGNE::create( const std::string name )
{
    return new CGNE( name );
}

} /* end namespace solver */

} /* end namespace scai */
