/**
 * @file TFQMR.cpp
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
 * @brief Implementation of methods for the TFQMR solver.
 * @author David Schissler
 * @date 13.05.2015
 */

// hpp
#include <scai/solver/TFQMR.hpp>

// internal scai libraries
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

// std
#include <limits>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, TFQMR<ValueType>::logger, "Solver.IterativeSolver.TFQMR" )

using lama::Matrix;
using lama::Vector;

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_Solver* TFQMR<ValueType>::create()
{
    return new TFQMR<ValueType>( "_genByFactory" );
}

template<typename ValueType>
SolverCreateKeyType TFQMR<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "TFQMR" );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
TFQMR<ValueType>::TFQMR( const std::string& id ) : 

    IterativeSolver<ValueType>( id )
{
}

template<typename ValueType>
TFQMR<ValueType>::TFQMR( const std::string& id, LoggerPtr logger ) : 

    IterativeSolver<ValueType>( id , logger )
{
}

template<typename ValueType>
TFQMR<ValueType>::TFQMR( const TFQMR& other ) : 

    IterativeSolver<ValueType>( other )
{
}

template<typename ValueType>
TFQMR<ValueType>::~TFQMR()
{
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void TFQMR<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )
    IterativeSolver<ValueType>::initialize( coefficients );

    TFQMRRuntime& runtime = getRuntime();

    runtime.mAlpha = ValueType( 0 );
    runtime.mBeta = ValueType( 0 );
    runtime.mC = ValueType( 0 );
    runtime.mEta = ValueType( 0 );
    runtime.mTheta = ValueType( 0 );

    runtime.mEps = common::TypeTraits<ValueType>::eps1() * ValueType( 3 );

    // create dense runtime vectors with same row distribution, type, context as coefficients

    runtime.mVecD.reset( coefficients.newTargetVector() );
    runtime.mInitialR.reset( coefficients.newTargetVector() );
    runtime.mVecVEven.reset( coefficients.newTargetVector() );
    runtime.mVecVOdd.reset( coefficients.newTargetVector() );
    runtime.mVecW.reset( coefficients.newTargetVector() );
    runtime.mVecZ.reset( coefficients.newTargetVector() );
    runtime.mVecVT.reset( coefficients.newTargetVector() );
}

/* ========================================================================= */
/*    Solve: init ( solution, rhs )                                          */
/* ========================================================================= */

template<typename ValueType>
void TFQMR<ValueType>::solveInit( Vector<ValueType>& solution, const Vector<ValueType>& rhs )
{
    IterativeSolver<ValueType>::solveInit( solution, rhs );

    TFQMRRuntime& runtime = getRuntime();

    // Initialize

    this->getResidual();

    const Matrix<ValueType>& A = *runtime.mCoefficients;
 
    Vector<ValueType>& initialR = *runtime.mInitialR;
    Vector<ValueType>& residual = *runtime.mResidual;
    Vector<ValueType>& vecW     = *runtime.mVecW;
    Vector<ValueType>& vecD     = *runtime.mVecD;
    Vector<ValueType>& vecZ     = *runtime.mVecZ;
    Vector<ValueType>& vecVEven = *runtime.mVecVEven;

    initialR = residual;
    vecVEven = residual;

    // PRECONDITIONING
    if ( mPreconditioner != NULL )
    {
        vecW = 0;
        mPreconditioner->solve( vecW , residual );
    }
    else
    {
        vecW = residual;
    }

    vecZ = A * vecW;
    vecW = residual;
    vecD = 0;

    runtime.mTau = initialR.l2Norm();
    runtime.mRhoOld = runtime.mTau * runtime.mTau;
}

/* ========================================================================= */
/*    Solve: iterate                                                         */
/* ========================================================================= */

template<typename ValueType>
void TFQMR<ValueType>::iterationEven()
{
    TFQMRRuntime& runtime = getRuntime();

    const Vector<ValueType>& vecZ = *runtime.mVecZ;
    const Vector<ValueType>& initialR = *runtime.mInitialR;
    const Vector<ValueType>& vecVEven = *runtime.mVecVEven;
    Vector<ValueType>& vecVOdd = *runtime.mVecVOdd;

    const ValueType& rho = runtime.mRhoOld;
    const RealType<ValueType>& eps = runtime.mEps;
    const ValueType dotProduct = initialR.dotProduct( vecZ );
    ValueType& alpha = runtime.mAlpha;

    if ( common::Math::abs( dotProduct ) < eps ) // scalar is small
    {
        alpha = 0.0;
    }
    else
    {
        alpha = rho / dotProduct;
    }

    vecVOdd  = vecVEven - alpha * vecZ;
}

template<typename ValueType>
void TFQMR<ValueType>::iterationOdd()
{
    TFQMRRuntime& runtime = getRuntime();
    const Matrix<ValueType>& A = *runtime.mCoefficients;
    const Vector<ValueType>& initialR = *runtime.mInitialR;
    const Vector<ValueType>& vecW = *runtime.mVecW;
    const Vector<ValueType>& vecVOdd = *runtime.mVecVOdd;
    Vector<ValueType>& vecVEven = *runtime.mVecVEven;

    ValueType& rhoOld = runtime.mRhoOld;
    ValueType& rhoNew = runtime.mRhoNew;
    ValueType& beta = runtime.mBeta;
    Vector<ValueType>& vecZ = *runtime.mVecZ;
    Vector<ValueType>& vecVT = *runtime.mVecVT;
    const RealType<ValueType>& eps = runtime.mEps;

    rhoNew  = initialR.dotProduct( vecW );

    if ( common::Math::abs( rhoOld ) < eps )  // scalar is small
    {
        beta = ValueType( 0 );
    }
    else
    {
        beta = rhoNew / rhoOld;
    }

    vecVEven = vecW + beta * vecVOdd;

    //  PRECONDITIONING
 
    if ( mPreconditioner != NULL )
    {
        vecVT.setSameValue( mPreconditioner->getCoefficients().getColDistributionPtr(), 0 );
        mPreconditioner->solve( vecVT, vecVEven );
    }
    else
    {
        vecVT = vecVEven;
    }

    vecZ *= beta;
    vecZ = beta * A * vecVOdd + beta * vecZ;
    vecZ = A * vecVT + vecZ;
    rhoOld = rhoNew;
}

template<typename ValueType>
void TFQMR<ValueType>::iterate()
{
    TFQMRRuntime& runtime = getRuntime();
    const IndexType& iteration = runtime.mIterations;
    const Matrix<ValueType>& A = *runtime.mCoefficients;
    Vector<ValueType>& vecW = *runtime.mVecW;
    Vector<ValueType>& vecD = *runtime.mVecD;
    Vector<ValueType>& solution = runtime.mSolution.getReference();  // dirty
    const ValueType& alpha = runtime.mAlpha;
    ValueType& c = runtime.mC;
    ValueType& eta = runtime.mEta;
    ValueType& theta = runtime.mTheta;
    ValueType& tau = runtime.mTau;

    RealType<ValueType>& eps = runtime.mEps;

    bool isEven = ( iteration % 2 ) == 0;

    const Vector<ValueType>& vecV = isEven ? *runtime.mVecVEven : *runtime.mVecVOdd;

    if ( isEven )
    {
        iterationEven();
    }

    vecW = vecW - alpha * A * vecV;
    ValueType tempScal;

    if (    common::Math::abs( alpha ) < eps 
         || common::Math::abs( theta ) < eps 
         || common::Math::abs( eta ) < eps ) // scalar is small
    {
        tempScal = 0.0;
    }
    else
    {
        tempScal = ( theta * theta / alpha ) * eta;
    }

    vecD = vecV + tempScal * vecD;

    if ( common::Math::abs( tau ) < eps ) // scalar is small
    {
        theta = 0.0;
    }
    else
    {
        theta = l2Norm( vecW ) / tau;
    }

    c = ValueType( 1 ) / common::Math::sqrt( ValueType( 1 ) + theta * theta );
    tau = tau * theta * c;
    eta = c * c * alpha;
    solution = solution + eta * vecD;

    if ( !isEven )
    {
        iterationOdd();
    }
}

template<typename ValueType>
TFQMR<ValueType>* TFQMR<ValueType>::copy()
{
    return new TFQMR( *this );
}

template<typename ValueType>
typename TFQMR<ValueType>::TFQMRRuntime& TFQMR<ValueType>::getRuntime()
{
    return mTFQMRRuntime;
}

template<typename ValueType>
const typename TFQMR<ValueType>::TFQMRRuntime& TFQMR<ValueType>::getRuntime() const
{
    return mTFQMRRuntime;
}

template<typename ValueType>
void TFQMR<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "TFQMR<" << common::TypeTraits<ValueType>::id() << "> ( id = " << Solver<ValueType>::getId()
           << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( TFQMR, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
