/**
 * @file QMR.cpp
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
 * @brief QMR.cpp
 * @author lschubert
 * @date 06.08.2013
 * @since 1.1.0
 */

// hpp
#include <scai/solver/QMR.hpp>

// internal scai libraries
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/matrix/Matrix.hpp>

#include <scai/lama/Vector.hpp>

// std
#include <limits>
#include <cstddef>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_LOGGER( QMR::logger, "Solver.QMR" )

using lama::Matrix;
using lama::Vector;
using lama::Scalar;

QMR::QMR( const std::string& id )
    : IterativeSolver( id )
{
}


QMR::QMR( const std::string& id, LoggerPtr logger )
    : IterativeSolver( id , logger )
{
}

QMR::QMR( const QMR& other )
    : IterativeSolver( other )
{
}

QMR::~QMR()
{
}

QMR::QMRRuntime::QMRRuntime()
    : IterativeSolverRuntime()
{
}

QMR::QMRRuntime::~QMRRuntime()
{
}

void QMR::initialize( const Matrix& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )

    IterativeSolver::initialize( coefficients );
    QMRRuntime& runtime = getRuntime();

    runtime.mEps = std::numeric_limits<double>::epsilon() * 3;                  //CAREFUL: No abstract type
    common::scalar::ScalarType type = coefficients.getValueType();

    runtime.mTransposeA.reset( coefficients.newMatrix() );
    runtime.mTransposeA->assignTranspose( coefficients );
    runtime.mTransposeA->conj();

    runtime.mVecD.reset( Vector::getDenseVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecP.reset( Vector::getDenseVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecQ.reset( Vector::getDenseVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecS.reset( Vector::getDenseVector( type, coefficients.getDistributionPtr() ) );    
    runtime.mVecV.reset( Vector::getDenseVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecW.reset( Vector::getDenseVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecY.reset( Vector::getDenseVector( type, coefficients.getDistributionPtr() ) );    /*preconditioning 1*/ 
    runtime.mVecZ.reset( Vector::getDenseVector( type, coefficients.getDistributionPtr() ) ); 
    runtime.mVecPT.reset( Vector::getDenseVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecVT.reset( Vector::getDenseVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecWT.reset( Vector::getDenseVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecYT.reset( Vector::getDenseVector( type, coefficients.getDistributionPtr() ) );
    runtime.mVecZT.reset( Vector::getDenseVector( type, coefficients.getDistributionPtr() ) );
   

    runtime.mVecD->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecP->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecQ->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecS->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecV->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecW->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecY->setContextPtr( coefficients.getContextPtr() );      /*preconditioning 1*/ 
    runtime.mVecZ->setContextPtr( coefficients.getContextPtr() ); 
    runtime.mVecVT->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecWT->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecYT->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecZT->setContextPtr( coefficients.getContextPtr() );
    runtime.mVecPT->setContextPtr( coefficients.getContextPtr() );
}


void QMR::solveInit( Vector& solution, const Vector& rhs )
{
    QMRRuntime& runtime = getRuntime();

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

    *runtime.mVecVT = *runtime.mResidual;
    *runtime.mVecWT = *runtime.mResidual;

    runtime.mSolveInit = true;
}

void QMR::iterate(){
    QMRRuntime& runtime    = getRuntime();

    const Matrix& A = *runtime.mCoefficients;
    const Matrix& transposedA = *runtime.mTransposeA;
    Vector& solution = *runtime.mSolution;
    Vector& residual = *runtime.mResidual;
    Vector& vecV = *runtime.mVecV;
    Vector& vecW = *runtime.mVecW;
    Vector& vecP = *runtime.mVecP;
    Vector& vecQ = *runtime.mVecQ;
    Vector& vecS = *runtime.mVecS;
    Vector& vecD = *runtime.mVecD;

    Vector& vecY = *runtime.mVecY;      /*preconditioning*/ 
    Vector& vecZ = *runtime.mVecZ;

    Vector& vecVT = *runtime.mVecVT;
    Vector& vecYT = *runtime.mVecYT;
    Vector& vecZT = *runtime.mVecZT;
    Vector& vecWT = *runtime.mVecWT;
    Vector& vecPT = *runtime.mVecPT;

    Scalar& gamma = runtime.mGamma;
    Scalar& theta = runtime.mTheta;
    Scalar& psi = runtime.mPsi;
    Scalar& rho = runtime.mRho;
    Scalar& epsilon = runtime.mEpsilon;
    Scalar& eta= runtime.mEta;

    Scalar gamma1;
    Scalar theta1;
    Scalar rho1;
    const Scalar& eps = runtime.mEps;
    lama::L2Norm norm;

    if(this->getIterationCount() == 0){
    /*PRECONDITIONING*/
        if(mPreconditioner != NULL)
        {
            vecY = Scalar(0.0);
            mPreconditioner->solve( vecY, vecVT );      
        } 
        else    vecY = vecVT;
        vecZ = vecWT;
        rho = norm(vecY);
        psi = norm(vecZ);
        gamma = 1.0;
        eta = -1.0;
    }
    if( abs(rho) < eps || abs(1.0/rho)<eps || abs(psi) < eps || abs(1.0/psi)<eps)
        return;
    vecV = vecVT/rho;
    vecY = vecY/rho;
    vecW = vecWT/psi;
    vecZ = vecZ/psi;
    Scalar delta = vecZ.dotProduct(vecY);

    if(abs(delta) < eps)
        return;
    /*PRECONDITIONING*/
    vecYT = vecY;
    if(mPreconditioner != NULL){
        vecZT = Scalar(0.0);
        mPreconditioner->solve( vecZT, vecZ );      
    } 
    else vecZT = vecZ;

    if(this->getIterationCount() == 0){
        vecP = vecYT;
        vecQ = vecZT;
    }
    else{
         Scalar pde = psi*delta/epsilon;
        if(abs(pde) < eps || abs(1.0/pde)<eps)
            return;
        Scalar rde = rho* conj(delta/epsilon);
        if(abs(rde) < eps || abs(1.0/rde)<eps)
            return;
        vecP = vecYT - pde*vecP;
        vecQ = vecZT - rde*vecQ;
    }
    vecPT = A*vecP;
    epsilon = vecQ.dotProduct(vecPT);
    if( abs(epsilon) < eps || abs(1.0/eps)<eps )
        return;
    Scalar beta = epsilon/delta;
    if(abs(beta)<eps || abs(1.0/beta)<eps)
        return;
    vecVT = vecPT - beta *vecV;

    /*PRECONDITIONING*/
    if(mPreconditioner != NULL)
    {
        vecY = Scalar(0.0);
        mPreconditioner->solve( vecY, vecVT );      
    }
    else    vecY = vecVT; 
 
    rho1= rho;
    rho = norm(vecY);

    vecWT = transposedA * vecQ;
    vecWT = vecWT - conj(beta)*vecW;

    vecZ = vecWT;
    psi = norm(vecZ);
    if(this->getIterationCount() > 0)
        theta1 = theta;

    theta = rho / (gamma * abs(beta));
    gamma1 = gamma;
    gamma = 1.0 / sqrt(1.0+theta*theta);
    if(abs(gamma) < eps)
        return;
    eta = -eta * rho1* gamma*gamma / (beta * gamma1 * gamma1);
    if(abs(1.0/eta) < eps)
        return;
    if(this->getIterationCount() == 0){
        vecD = eta * vecP;
        vecS = eta * vecPT;
    }
    else{
        vecD = eta * vecP + (theta1*gamma)*(theta1*gamma)*vecD;
        vecS = eta * vecPT + (theta1*gamma)*(theta1*gamma)*vecS;
    }

    solution = solution + vecD;
    residual = residual - vecS;
    mQMRRuntime.mSolution.setDirty( false );
}

SolverPtr QMR::copy()
{
    return SolverPtr( new QMR( *this ) );
}

QMR::QMRRuntime& QMR::getRuntime()
{
    return mQMRRuntime;
}

const QMR::QMRRuntime& QMR::getConstRuntime() const
{
    return mQMRRuntime;
}

std::string QMR::createValue()
{
    return "QMR";
}

Solver* QMR::create( const std::string name )
{
    return new QMR( name );
}

void QMR::writeAt( std::ostream& stream ) const
{
    stream << "QMR ( id = " << mId << ", #iter = " << getConstRuntime().mIterations << " )";
}

} /* end namespace solver */

} /* end namespace scai */
