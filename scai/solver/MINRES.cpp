/**
 * @file MINRES.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Implementation of methods for the MINRES solver.
 * @author David Schissler
 * @date 13.05.2015
 */

// hpp
#include <scai/solver/MINRES.hpp>

// internal scai libraries
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/Vector.hpp>

#include <scai/common/ScalarType.hpp>
#include <scai/common/macros/instantiate.hpp>

// std
#include <limits>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, MINRES<ValueType>::logger, "Solver.IterativeSolver.MINRES" )

using lama::Matrix;
using lama::Vector;

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_Solver* MINRES<ValueType>::create()
{   
    return new MINRES<ValueType>( "_genByFactory" );
}

template<typename ValueType>
SolverCreateKeyType MINRES<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "MINRES" );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
MINRES<ValueType>::MINRES( const std::string& id ) : 

    IterativeSolver<ValueType>( id ) 
{
}

template<typename ValueType>
MINRES<ValueType>::MINRES( const std::string& id, LoggerPtr logger ) : 

    IterativeSolver<ValueType>( id , logger ) 
{
}

template<typename ValueType>
MINRES<ValueType>::MINRES( const MINRES<ValueType>& other ) : 

    IterativeSolver<ValueType>( other )
{
}


template<typename ValueType>
MINRES<ValueType>::~MINRES() 
{
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void MINRES<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )
    IterativeSolver<ValueType>::initialize( coefficients );
    MINRESRuntime& runtime = getRuntime();
    runtime.mBetaNew = 0.0;
    runtime.mC = 1.0;
    runtime.mCNew = 1.0;
    runtime.mS = 0.0;
    runtime.mSNew = 0.0;

    runtime.mVecV.reset( coefficients.newTargetVector() );
    runtime.mVecVOld.reset( coefficients.newTargetVector() );
    runtime.mVecVNew.reset( coefficients.newTargetVector() );
    runtime.mVecP.reset( coefficients.newTargetVector() );
    runtime.mVecPOld.reset( coefficients.newTargetVector() );
    runtime.mVecPNew.reset( coefficients.newTargetVector() );
    runtime.mEps = common::TypeTraits<ValueType>::eps1() * 3;
}

template<typename ValueType>
void MINRES<ValueType>::solveInit( Vector<ValueType>& solution, const Vector<ValueType>& rhs )
{
    MINRESRuntime& runtime = getRuntime();
    runtime.mRhs = &rhs;
    runtime.mSolution = &solution;
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumRows(), rhs.size(), "mismatch: #rows of matrix, rhs" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getNumColumns(), solution.size(), "mismatch: #cols of matrix, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getColDistribution(), solution.getDistribution(), "mismatch: matrix col dist, solution" )
    SCAI_ASSERT_EQUAL( runtime.mCoefficients->getRowDistribution(), rhs.getDistribution(), "mismatch: matrix row dist, rhs dist" )
    // Initialize
    this->getResidual();

    *runtime.mVecVNew = *runtime.mResidual;
    *runtime.mVecV = 0;
    *runtime.mVecP = 0;
    *runtime.mVecPNew = 0;

    runtime.mZeta = runtime.mResidual->l2Norm();
    *runtime.mVecVNew /= runtime.mZeta;  // *runtime.mVecVNew /= runtime.mZeta
    runtime.mSolveInit = true;
}

template<typename ValueType>
void MINRES<ValueType>::Lanczos()
{
    MINRESRuntime& runtime = getRuntime();

    const Matrix<ValueType>& A = *runtime.mCoefficients;

    Vector<ValueType>& vecV = *runtime.mVecV;
    Vector<ValueType>& vecVOld = *runtime.mVecVOld;
    Vector<ValueType>& vecVNew = *runtime.mVecVNew;

    ValueType& alpha = runtime.mAlpha;
    RealType<ValueType>& betaNew = runtime.mBetaNew;
    RealType<ValueType>& beta = runtime.mBeta;
    RealType<ValueType>& eps = runtime.mEps;

    beta = betaNew;
    vecVOld.swap( vecV );
    vecV.swap( vecVNew );
    vecVNew = A * vecV - ValueType( beta ) * vecVOld;
    alpha = vecV.dotProduct( vecVNew );
    vecVNew = vecVNew - alpha * vecV;
    betaNew = vecVNew.l2Norm();

    if ( common::Math::abs( betaNew ) < eps || 1 / common::Math::abs( betaNew ) < eps )
    {
        vecVNew = ValueType( 0 );
    }
    else
    {
        vecVNew = vecVNew / ValueType( betaNew );
    }
}

template<typename ValueType>
void MINRES<ValueType>::applyGivensRotation()
{
    MINRESRuntime& runtime = getRuntime();
    ValueType& c = runtime.mC;
    ValueType& cOld = runtime.mCOld;
    ValueType& cNew = runtime.mCNew;
    ValueType& s = runtime.mS;
    ValueType& sOld = runtime.mSOld;
    ValueType& sNew = runtime.mSNew;
    ValueType& alpha = runtime.mAlpha;
    RealType<ValueType>& beta = runtime.mBeta;
    RealType<ValueType>& betaNew = runtime.mBetaNew;

    RealType<ValueType>& eps = runtime.mEps;

    //Old Givens-rotation

    cOld = c;
    c = cNew;
    sOld = s;
    s = sNew;

    ValueType rho1 = sOld * beta;
    ValueType rho2 = c * cOld * beta + s * alpha;
    ValueType rho3 = c * alpha - s * cOld * beta;

    //New Givens-rotation

    RealType<ValueType> tau = common::Math::abs( rho3 ) + betaNew;

    ValueType nu;

    if ( common::Math::abs( tau ) < eps || 1.0 / common::Math::abs( tau ) < eps )
    {
        nu = 1.0;
        rho3 = 1.0;
        cNew = 0.0;
        sNew = 0.0;
    }
    else
    {
        nu = tau * common::Math::sqrt( ( rho3 / tau ) * ( rho3 / tau ) + ( betaNew / tau ) * ( betaNew / tau ) );
        cNew = rho3 / nu;
        sNew = betaNew / nu;
        rho3 = nu;
    }

    std::swap( runtime.mVecPOld, runtime.mVecP );
    std::swap( runtime.mVecP, runtime.mVecPNew );

    Vector<ValueType>& vecP    = *runtime.mVecP;
    Vector<ValueType>& vecPOld = *runtime.mVecPOld;
    Vector<ValueType>& vecPNew = *runtime.mVecPNew;

    const Vector<ValueType>& vecV = *runtime.mVecV;

    /*
    // update P(new) -> P -> POld

    vecPOld.swap( vecP );
    vecP.swap( vecPNew );
    */

    vecPNew = vecV - rho1 * vecPOld;
    vecPNew = vecPNew - rho2 * vecP;
    vecPNew = vecPNew / rho3;
}

template<typename ValueType>
void MINRES<ValueType>::iterate()
{
    MINRESRuntime& runtime = getRuntime();

    Lanczos();
    applyGivensRotation();

    const Vector<ValueType>& vecPNew = *runtime.mVecPNew;
    Vector<ValueType>& solution = runtime.mSolution.getReference(); // dirty

    ValueType& cNew = runtime.mCNew;
    ValueType& sNew = runtime.mSNew;
    ValueType& zeta = runtime.mZeta;

    // New approximation

    solution = solution + cNew * zeta * vecPNew;
    zeta = -sNew * zeta;
}

template<typename ValueType>
MINRES<ValueType>* MINRES<ValueType>::copy()
{
    return new MINRES( *this );
}

template<typename ValueType>
typename MINRES<ValueType>::MINRESRuntime& MINRES<ValueType>::getRuntime()
{
    return mMINRESRuntime;
}

template<typename ValueType>
const typename MINRES<ValueType>::MINRESRuntime& MINRES<ValueType>::getRuntime() const
{
    return mMINRESRuntime;
}

template<typename ValueType>
void MINRES<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "MINRES<" << common::TypeTraits<ValueType>::id() << "> ( id = " << Solver<ValueType>::getId()
           << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( MINRES, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
