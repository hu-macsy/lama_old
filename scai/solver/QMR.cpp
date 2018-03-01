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

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

// std
#include <limits>
#include <cstddef>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, QMR<ValueType>::logger, "Solver.IterativeSolver.QMR" )

using lama::Matrix;
using lama::Vector;

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_Solver* QMR<ValueType>::create()
{
    return new QMR<ValueType>( "_genByFactory" );
}

template<typename ValueType>
SolverCreateKeyType QMR<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "QMR" );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
QMR<ValueType>::QMR( const std::string& id ) : 

    IterativeSolver<ValueType>( id )
{
}


template<typename ValueType>
QMR<ValueType>::QMR( const std::string& id, LoggerPtr logger ) : 

    IterativeSolver<ValueType>( id , logger )
{
}

template<typename ValueType>
QMR<ValueType>::QMR( const QMR<ValueType>& other ) : 

    IterativeSolver<ValueType>( other )
{
    // does not copy runtime data
}

template<typename ValueType>
QMR<ValueType>::~QMR()
{
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void QMR<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )

    IterativeSolver<ValueType>::initialize( coefficients );

    QMRRuntime& runtime = getRuntime();

    runtime.mEps = common::TypeTraits<ValueType>::eps1() * ValueType( 3 );

    dmemo::DistributionPtr rowDist = coefficients.getRowDistributionPtr();
    hmemo::ContextPtr ctx = coefficients.getContextPtr();

    runtime.mConjTransposeA.reset( coefficients.newMatrix() );
    runtime.mConjTransposeA->assignTranspose( coefficients );
    runtime.mConjTransposeA->conj();

    runtime.mVecD.reset( coefficients.newTargetVector() );
    runtime.mVecP.reset( coefficients.newTargetVector() );
    runtime.mVecQ.reset( coefficients.newTargetVector() );
    runtime.mVecS.reset( coefficients.newTargetVector() );
    runtime.mVecV.reset( coefficients.newTargetVector() );
    runtime.mVecW.reset( coefficients.newTargetVector() );
    runtime.mVecY.reset( coefficients.newTargetVector() );
    runtime.mVecZ.reset( coefficients.newTargetVector() );
    runtime.mVecPT.reset( coefficients.newTargetVector() );
    runtime.mVecVT.reset( coefficients.newTargetVector() );
    runtime.mVecWT.reset( coefficients.newTargetVector() );
    runtime.mVecYT.reset( coefficients.newTargetVector() );
    runtime.mVecZT.reset( coefficients.newTargetVector() );
}

/* ========================================================================= */
/*    solve : init ( solution, rhs )                                         */
/* ========================================================================= */

template<typename ValueType>
void QMR<ValueType>::solveInit( Vector<ValueType>& solution, const Vector<ValueType>& rhs )
{
    IterativeSolver<ValueType>::solveInit( solution, rhs );

    this->getResidual();

    QMRRuntime& runtime = getRuntime();

    *runtime.mVecVT = *runtime.mResidual;
    *runtime.mVecWT = *runtime.mResidual;
}

/* ========================================================================= */
/*    solve : iterate ( computations for one iteration step )                */
/* ========================================================================= */

template<typename ValueType>
void QMR<ValueType>::iterate()
{
    QMRRuntime& runtime    = getRuntime();

    const Matrix<ValueType>& A = *runtime.mCoefficients;
    const Matrix<ValueType>& Act = *runtime.mConjTransposeA;

    Vector<ValueType>& solution = runtime.mSolution.getReference(); // -> dirty
    Vector<ValueType>& residual = *runtime.mResidual;

    Vector<ValueType>& vecV = *runtime.mVecV;
    Vector<ValueType>& vecW = *runtime.mVecW;
    Vector<ValueType>& vecP = *runtime.mVecP;
    Vector<ValueType>& vecQ = *runtime.mVecQ;
    Vector<ValueType>& vecS = *runtime.mVecS;
    Vector<ValueType>& vecD = *runtime.mVecD;
    Vector<ValueType>& vecY = *runtime.mVecY;      /*preconditioning*/
    Vector<ValueType>& vecZ = *runtime.mVecZ;
    Vector<ValueType>& vecVT = *runtime.mVecVT;
    Vector<ValueType>& vecYT = *runtime.mVecYT;
    Vector<ValueType>& vecZT = *runtime.mVecZT;
    Vector<ValueType>& vecWT = *runtime.mVecWT;
    Vector<ValueType>& vecPT = *runtime.mVecPT;

    ValueType& gamma = runtime.mGamma;
    ValueType& theta = runtime.mTheta;
    ValueType& psi = runtime.mPsi;
    ValueType& rho = runtime.mRho;
    ValueType& epsilon = runtime.mEpsilon;
    ValueType& eta = runtime.mEta;
    ValueType gamma1;
    ValueType theta1 = 0;  // not used in first iteration
    ValueType rho1;

    const RealType<ValueType>& eps = runtime.mEps;

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
        rho = l2Norm( vecY );
        psi = l2Norm( vecZ );
        gamma = 1.0;
        eta = -1.0;
    }

    if ( common::Math::abs( rho ) < eps || 
         common::Math::abs( 1.0 / rho ) < eps || 
         common::Math::abs( psi ) < eps || 
         common::Math::abs( 1.0 / psi ) < eps )
    {
        return;
    }

    vecV = vecVT / rho;
    vecY = vecY / rho;
    vecW = vecWT / psi;
    vecZ = vecZ / psi;
    ValueType delta = vecZ.dotProduct( vecY );

    if ( common::Math::abs( delta ) < eps )
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
        ValueType pde = psi * delta / epsilon;

        if ( common::Math::abs( pde ) < eps || common::Math::abs( 1.0 / pde ) < eps )
        {
            return;
        }

        ValueType rde = rho * common::Math::conj( delta / epsilon );

        if ( common::Math::abs( rde ) < eps || common::Math::abs( 1.0 / rde ) < eps )
        {
            return;
        }

        vecP = vecYT - pde * vecP;
        vecQ = vecZT - rde * vecQ;
    }

    vecPT = A * vecP;
    epsilon = vecQ.dotProduct( vecPT );

    if ( common::Math::abs( epsilon ) < eps || common::Math::abs( 1.0 / eps ) < eps )
    {
        return;
    }

    ValueType beta = epsilon / delta;

    if ( common::Math::abs( beta ) < eps || common::Math::abs( 1.0 / beta ) < eps )
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
    rho = l2Norm( vecY );
    vecWT = Act * vecQ;  // conjTranspose( A ) * vecQ
    vecWT = vecWT - common::Math::conj( beta ) * vecW;
    vecZ = vecWT;
    psi = l2Norm( vecZ );

    if ( this->getIterationCount() > 0 )
    {
        theta1 = theta;
    }

    theta = rho / ( gamma * common::Math::abs( beta ) );
    gamma1 = gamma;
    gamma = ValueType( 1 ) / common::Math::sqrt( ValueType( 1 ) + theta * theta );

    if ( common::Math::abs( gamma ) < eps )
    {
        return;
    }

    eta = -eta * rho1 * gamma * gamma / ( beta * gamma1 * gamma1 );

    if ( common::Math::abs( 1.0 / eta ) < eps )
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
        // Note: theta1 is the value from previous iteration
        vecD = eta * vecP + ( theta1 * gamma ) * ( theta1 * gamma ) * vecD;
        vecS = eta * vecPT + ( theta1 * gamma ) * ( theta1 * gamma ) * vecS;
    }

    solution = solution + vecD;
    residual = residual - vecS;

    mQMRRuntime.mSolution.setDirty( false );  // update of residual already done here
}

/* ========================================================================= */
/*       Getter runtime                                                      */
/* ========================================================================= */

template<typename ValueType>
typename QMR<ValueType>::QMRRuntime& QMR<ValueType>::getRuntime()
{
    return mQMRRuntime;
}

template<typename ValueType>
const typename QMR<ValueType>::QMRRuntime& QMR<ValueType>::getRuntime() const
{
    return mQMRRuntime;
}

/* ========================================================================= */
/*       virtual methods                                                     */
/* ========================================================================= */

template<typename ValueType>
QMR<ValueType>* QMR<ValueType>::copy()
{
    return new QMR( *this );
}

template<typename ValueType>
void QMR<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "QMR<" << this->getValueType() << "> ( id = " << this->getId()
           << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( QMR, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
