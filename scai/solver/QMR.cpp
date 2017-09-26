/**
 * @file QMR.cpp
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
 * @brief Implementation of methods for the QMR solver.
 * @author Lauretta Schubert
 * @date 06.08.2013
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
    runtime.mEps = Scalar::eps1( coefficients.getValueType() ) * 3.0;
    runtime.mTransposeA.reset( coefficients.newMatrix() );
    runtime.mTransposeA->assignTranspose( coefficients );
    runtime.mTransposeA->conj();
    dmemo::DistributionPtr rowDist = coefficients.getRowDistributionPtr();
    runtime.mVecD.reset( coefficients.newVector( rowDist ) );
    runtime.mVecP.reset( coefficients.newVector( rowDist ) );
    runtime.mVecQ.reset( coefficients.newVector( rowDist ) );
    runtime.mVecS.reset( coefficients.newVector( rowDist ) );
    runtime.mVecV.reset( coefficients.newVector( rowDist ) );
    runtime.mVecW.reset( coefficients.newVector( rowDist ) );
    runtime.mVecY.reset( coefficients.newVector( rowDist ) );  // preconditioning 1
    runtime.mVecZ.reset( coefficients.newVector( rowDist ) );
    runtime.mVecPT.reset( coefficients.newVector( rowDist ) );
    runtime.mVecVT.reset( coefficients.newVector( rowDist ) );
    runtime.mVecWT.reset( coefficients.newVector( rowDist ) );
    runtime.mVecYT.reset( coefficients.newVector( rowDist ) );
    runtime.mVecZT.reset( coefficients.newVector( rowDist ) );
}

void QMR::solveInit( Vector& solution, const Vector& rhs )
{
    QMRRuntime& runtime = getRuntime();
    runtime.mRhs = &rhs;
    runtime.mSolution = &solution;
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumRows(), rhs.size(), "mismatch: #rows of matrix, rhs" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumColumns(), solution.size(), "mismatch: #cols of matrix, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getColDistribution(), solution.getDistribution(), "mismatch: matrix col dist, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getRowDistribution(), rhs.getDistribution(), "mismatch: matrix row dist, rhs dist" )
    // Initialize
    this->getResidual();
    *runtime.mVecVT = *runtime.mResidual;
    *runtime.mVecWT = *runtime.mResidual;
    runtime.mSolveInit = true;
}

void QMR::iterate()
{
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
    Scalar& eta = runtime.mEta;
    Scalar gamma1;
    Scalar theta1;
    Scalar rho1;
    const Scalar& eps = runtime.mEps;
    lama::L2Norm norm;

    if ( this->getIterationCount() == 0 )
    {
        /*PRECONDITIONING*/
        if ( mPreconditioner != NULL )
        {
            vecY.setSameValue( vecVT.getDistributionPtr(), 0 );
            mPreconditioner->solve( vecY, vecVT );
        }
        else
        {
            vecY = vecVT;
        }

        vecZ = vecWT;
        rho = norm( vecY );
        psi = norm( vecZ );
        gamma = 1.0;
        eta = -1.0;
    }

    if ( abs( rho ) < eps || abs( 1.0 / rho ) < eps || abs( psi ) < eps || abs( 1.0 / psi ) < eps )
    {
        return;
    }

    vecV = vecVT / rho;
    vecY = vecY / rho;
    vecW = vecWT / psi;
    vecZ = vecZ / psi;
    Scalar delta = vecZ.dotProduct( vecY );

    if ( abs( delta ) < eps )
    {
        return;
    }

    /*PRECONDITIONING*/
    vecYT = vecY;

    if ( mPreconditioner != NULL )
    {
        // vecZT = 0, here we make it more safe to be sure about the size
        vecZT.setSameValue( vecZ.getDistributionPtr(), 0 );
        mPreconditioner->solve( vecZT, vecZ );
    }
    else
    {
        vecZT = vecZ;
    }

    if ( this->getIterationCount() == 0 )
    {
        vecP = vecYT;
        vecQ = vecZT;
    }
    else
    {
        Scalar pde = psi * delta / epsilon;

        if ( abs( pde ) < eps || abs( 1.0 / pde ) < eps )
        {
            return;
        }

        Scalar rde = rho * conj( delta / epsilon );

        if ( abs( rde ) < eps || abs( 1.0 / rde ) < eps )
        {
            return;
        }

        vecP = vecYT - pde * vecP;
        vecQ = vecZT - rde * vecQ;
    }

    vecPT = A * vecP;
    epsilon = vecQ.dotProduct( vecPT );

    if ( abs( epsilon ) < eps || abs( 1.0 / eps ) < eps )
    {
        return;
    }

    Scalar beta = epsilon / delta;

    if ( abs( beta ) < eps || abs( 1.0 / beta ) < eps )
    {
        return;
    }

    vecVT = vecPT - beta * vecV;

    /*PRECONDITIONING*/
    if ( mPreconditioner != NULL )
    {
        vecY.setSameValue( vecVT.getDistributionPtr(), 0 );
        mPreconditioner->solve( vecY, vecVT );
    }
    else
    {
        vecY = vecVT;
    }

    rho1 = rho;
    rho = norm( vecY );
    vecWT = transposedA * vecQ;
    vecWT = vecWT - conj( beta ) * vecW;
    vecZ = vecWT;
    psi = norm( vecZ );

    if ( this->getIterationCount() > 0 )
    {
        theta1 = theta;
    }

    theta = rho / ( gamma * abs( beta ) );
    gamma1 = gamma;
    gamma = 1.0 / sqrt( 1.0 + theta * theta );

    if ( abs( gamma ) < eps )
    {
        return;
    }

    eta = -eta * rho1 * gamma * gamma / ( beta * gamma1 * gamma1 );

    if ( abs( 1.0 / eta ) < eps )
    {
        return;
    }

    if ( this->getIterationCount() == 0 )
    {
        vecD = eta * vecP;
        vecS = eta * vecPT;
    }
    else
    {
        vecD = eta * vecP + ( theta1 * gamma ) * ( theta1 * gamma ) * vecD;
        vecS = eta * vecPT + ( theta1 * gamma ) * ( theta1 * gamma ) * vecS;
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
