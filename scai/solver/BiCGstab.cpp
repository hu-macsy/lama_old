/**
 * @file BiCGstab.cpp
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
 * @brief BiCGstab.cpp
 * @author lschubert
 * @date 06.08.2013
 */

// hpp
#include <scai/solver/BiCGstab.hpp>

// internal scai libraries
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/Vector.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

// std
#include <limits>

namespace scai
{

using namespace lama;

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, BiCGstab<ValueType>::logger, "Solver.IterativeSolver.BiCGstab" )

using lama::Matrix;
using lama::Vector;

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_Solver* BiCGstab<ValueType>::create()
{
    return new BiCGstab<ValueType>( "_genByFactory" );
}

template<typename ValueType>
SolverCreateKeyType BiCGstab<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "BiCGstab" );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
BiCGstab<ValueType>::BiCGstab( const std::string& id ) :

    IterativeSolver<ValueType>( id )
{
}

template<typename ValueType>
BiCGstab<ValueType>::BiCGstab( const std::string& id, LoggerPtr logger ) :

    IterativeSolver<ValueType>( id, logger )
{
}

template<typename ValueType>
BiCGstab<ValueType>::BiCGstab( const BiCGstab<ValueType>& other ) :

    IterativeSolver<ValueType>( other )
{
}

template<typename ValueType>
BiCGstab<ValueType>::~BiCGstab()
{
    SCAI_LOG_INFO( logger, "~BiCGstab<" << this->getValueType() << ">" )
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void BiCGstab<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )

    IterativeSolver<ValueType>::initialize( coefficients );

    BiCGstabRuntime& runtime = getRuntime();

    runtime.mAlpha  = ValueType( 1 );
    runtime.mOmega = ValueType( 1 );
    runtime.mRhoOld = ValueType( 1 );
    runtime.mResNorm = ValueType( 1 );

    runtime.mEps = common::TypeTraits<ValueType>::eps1() * ValueType( 3 );

    // get runtime vectors with same row distribution / context / type as cofficients matrix

    runtime.mRes0.reset( coefficients.newTargetVector() );
    runtime.mVecV.reset( coefficients.newTargetVector() );
    runtime.mVecP.reset( coefficients.newTargetVector() );
    runtime.mVecS.reset( coefficients.newTargetVector() );
    runtime.mVecT.reset( coefficients.newTargetVector() );
    runtime.mVecPT.reset( coefficients.newTargetVector() );
    runtime.mVecST.reset( coefficients.newTargetVector() );
    runtime.mVecTT.reset( coefficients.newTargetVector() );
}

/* ========================================================================= */
/*    solve: init ( solution, rhs )                                          */
/* ========================================================================= */

template<typename ValueType>
void BiCGstab<ValueType>::solveInit( Vector<ValueType>& solution, const Vector<ValueType>& rhs )
{
    IterativeSolver<ValueType>::solveInit( solution, rhs );

    BiCGstabRuntime& runtime = getRuntime();

    // Initialize

    this->getResidual();

    *runtime.mRes0 = *runtime.mResidual;    // deep copy of the residual

    runtime.mVecV.reset( rhs.newVector() );
    runtime.mVecP.reset( rhs.newVector() );
}

/* ========================================================================= */
/*    solve: iterate ( one iteration step )                                  */
/* ========================================================================= */

template<typename ValueType>
void BiCGstab<ValueType>::iterate()
{
    BiCGstabRuntime& runtime    = getRuntime();

    const Matrix<ValueType>& A = *runtime.mCoefficients;
    const Vector<ValueType>& res0 = *runtime.mRes0;

    Vector<ValueType>& res = *runtime.mResidual;
    Vector<ValueType>& vecV = *runtime.mVecV;
    Vector<ValueType>& vecP = *runtime.mVecP;
    Vector<ValueType>& vecS = *runtime.mVecS;
    Vector<ValueType>& vecT = *runtime.mVecT;
    Vector<ValueType>& solution = runtime.mSolution.getReference(); // -> dirty
    Vector<ValueType>& vecPT = *runtime.mVecPT;
    Vector<ValueType>& vecST = *runtime.mVecST;
    Vector<ValueType>& vecTT = *runtime.mVecTT;

    ValueType& alpha = runtime.mAlpha;
    ValueType& beta = runtime.mBeta;
    ValueType& omega = runtime.mOmega;
    ValueType& rhoOld = runtime.mRhoOld;
    ValueType& rhoNew = runtime.mRhoNew;
    const RealType<ValueType>& eps = runtime.mEps;

    RealType<ValueType>& resNorm = runtime.mResNorm;

    rhoNew = res0.dotProduct( res );

    if ( resNorm < eps || 
         common::Math::abs( rhoOld ) < eps || 
         common::Math::abs( omega ) < eps ) // scalars are small
    {
        beta = 0.0;
    }
    else
    {
        beta = rhoNew / rhoOld * ( alpha / omega );
    }

    SCAI_LOG_INFO( logger, "resNorm = " << resNorm << ", eps = " << eps << ", beta = " << beta )
    vecP = vecP - omega * vecV;
    vecP = res + beta * vecP;

    // PRECONDITIONING
    if ( mPreconditioner != NULL )
    {
        vecPT.setSameValue( A.getRowDistributionPtr(), 0 );
        mPreconditioner->solve( vecPT, vecP );
    }
    else
    {
        vecPT = vecP;
    }

    vecV = A * vecPT;

    ValueType innerProd = res0.dotProduct( vecV );

    if ( resNorm < eps || common::Math::abs( innerProd ) < eps ) // scalar is small
    {
        alpha = 0.0;
    }
    else
    {
        alpha = rhoNew / innerProd;
    }

    vecS = res - alpha * vecV;

    // PRECONDITIONING
    if ( mPreconditioner != NULL )
    {
        vecST.setSameValue( A.getColDistributionPtr(), 0 );
        mPreconditioner->solve( vecST, vecS );
        vecT = A * vecST;
        vecTT.setSameValue( A.getColDistributionPtr(), 0 );
        mPreconditioner->solve( vecTT, vecT );
    }
    else
    {
        vecST = vecS;
        vecT = A * vecST;
        vecTT = vecT;
    }

    RealType<ValueType> normTT2 = vecTT.dotProduct( vecTT );

    if ( resNorm < eps || normTT2 < eps ) //scalar is small
    {
        omega = 0.0;
    }
    else
    {
        omega = vecTT.dotProduct( vecST ) / normTT2;
    }

    solution = solution + alpha * vecP;
    solution = solution + omega * vecS;

    res = vecS - omega * vecT;
    rhoOld = rhoNew;
    resNorm = l2Norm( res );

    // as we have updated residual already here, we can mark solution as not dirty

    mBiCGstabRuntime.mSolution.setDirty( false );  
}

/* ========================================================================= */
/*       Getter runtime                                                      */
/* ========================================================================= */

template<typename ValueType>
typename BiCGstab<ValueType>::BiCGstabRuntime& BiCGstab<ValueType>::getRuntime()
{
    return mBiCGstabRuntime;
}

template<typename ValueType>
const typename BiCGstab<ValueType>::BiCGstabRuntime& BiCGstab<ValueType>::getRuntime() const
{
    return mBiCGstabRuntime;
}

/* ========================================================================= */
/*       virtual methods                                                     */
/* ========================================================================= */

template<typename ValueType>
BiCGstab<ValueType>* BiCGstab<ValueType>::copy()
{
    return new BiCGstab( *this );
}

template<typename ValueType>
void BiCGstab<ValueType>::writeAt( std::ostream& stream ) const
{
    const char* typeId = common::TypeTraits<ValueType>::id();

    stream << "BiCGstab<" << typeId << "> ( id = " << this->getId()
           << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( BiCGstab, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace sovler */

} /* end namespace scai */
