/**
 * @file CGS.cpp
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
 * @brief CGS.cpp
 * @author David Schissler
 * @date 18.05.2015
 */

// hpp
#include <scai/solver/CGS.hpp>

// internal scai libraries
#include <scai/lama/expression/VectorExpressions.hpp>
#include <scai/lama/expression/MatrixExpressions.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/Vector.hpp>

// common
#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

// std
#include <limits>

namespace scai
{

namespace solver
{

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, CGS<ValueType>::logger, "Solver.IterativeSolver.CGS" )

using lama::Matrix;
using lama::Vector;

/* ========================================================================= */
/*    static methods (for factory)                                           */
/* ========================================================================= */

template<typename ValueType>
_Solver* CGS<ValueType>::create()
{
    return new CGS<ValueType>( "_genByFactory" );
}

template<typename ValueType>
SolverCreateKeyType CGS<ValueType>::createValue()
{
    return SolverCreateKeyType( common::getScalarType<ValueType>(), "CGS" );
}

/* ========================================================================= */
/*    Constructor/Destructor                                                 */
/* ========================================================================= */

template<typename ValueType>
CGS<ValueType>::CGS( const std::string& id ) :

    IterativeSolver<ValueType>( id )
{
}

template<typename ValueType>
CGS<ValueType>::CGS( const std::string& id, LoggerPtr logger ) :

    IterativeSolver<ValueType>( id, logger )
{
}

template<typename ValueType>
CGS<ValueType>::CGS( const CGS& other ) :

    IterativeSolver<ValueType>( other )
{
}

template<typename ValueType>
CGS<ValueType>::~CGS()
{
}

/* ========================================================================= */
/*    Initializaition                                                        */
/* ========================================================================= */

template<typename ValueType>
void CGS<ValueType>::initialize( const Matrix<ValueType>& coefficients )
{
    SCAI_LOG_DEBUG( logger, "Initialization started for coefficients = " << coefficients )

    IterativeSolver<ValueType>::initialize( coefficients );

    CGSRuntime& runtime = getRuntime();

    runtime.mNormRes = 1.0;
    runtime.mEps = common::TypeTraits<ValueType>::eps1() * ValueType( 3 );

    // temporary vectors are distributed corresponding the dist of target(rows) space

    dmemo::DistributionPtr dist = coefficients.getRowDistributionPtr();
    hmemo::ContextPtr ctx = coefficients.getContextPtr();

    runtime.mRes0.reset( coefficients.newSourceVector() );
    runtime.mVecT.reset( coefficients.newTargetVector() );
    runtime.mVecP.reset( coefficients.newTargetVector() );
    runtime.mVecQ.reset( coefficients.newTargetVector() );
    runtime.mVecU.reset( coefficients.newTargetVector() );
    runtime.mVecPT.reset( coefficients.newTargetVector() );
    runtime.mVecUT.reset( coefficients.newTargetVector() );
    runtime.mVecTemp.reset( coefficients.newTargetVector() );
}

/* ========================================================================= */
/*    solve: init                                                            */
/* ========================================================================= */

template<typename ValueType>
void CGS<ValueType>::solveInit( Vector<ValueType>& solution, const Vector<ValueType>& rhs )
{
    IterativeSolver<ValueType>::solveInit( solution, rhs );

    CGSRuntime& runtime = getRuntime();

    this->getResidual();

    const Vector<ValueType>& residual = *runtime.mResidual;

    Vector<ValueType>& res0 = *runtime.mRes0;
    Vector<ValueType>& vecP = *runtime.mVecP;
    Vector<ValueType>& vecU = *runtime.mVecU;
    Vector<ValueType>& vecPT = *runtime.mVecPT;

    // (deep) copy of the residual to res0, vecP, vecU

    res0 = residual;
    vecP = residual;
    vecU = residual;

    // PRECONDITIONING

    if ( mPreconditioner != NULL )
    {
        vecPT = 0;
        mPreconditioner->solve( vecPT, vecP );
    }
    else
    {
        vecPT = vecP;
    }

    //initial <res,res> inner product;

    runtime.mInnerProdRes = res0.dotProduct( res0 );
    SCAI_LOG_INFO( logger, "solveInit, mInnerProdRes = " << runtime.mInnerProdRes )
    runtime.mSolveInit = true;
}

/* ========================================================================= */
/*    solve: iterate                                                         */
/* ========================================================================= */

template<typename ValueType>
void CGS<ValueType>::iterate()
{
    CGSRuntime& runtime = getRuntime();

    const Matrix<ValueType>& A = *runtime.mCoefficients;

    const Vector<ValueType>& res0 = *runtime.mRes0;
    Vector<ValueType>& res = *runtime.mResidual;
    Vector<ValueType>& vecP = *runtime.mVecP;
    Vector<ValueType>& vecQ = *runtime.mVecQ;
    Vector<ValueType>& vecU = *runtime.mVecU;
    Vector<ValueType>& vecT = *runtime.mVecT;
    Vector<ValueType>& vecPT = *runtime.mVecPT;
    Vector<ValueType>& vecUT = *runtime.mVecUT;
    Vector<ValueType>& vecTemp = *runtime.mVecTemp;
    Vector<ValueType>& solution = runtime.mSolution.getReference(); // -> dirty

    ValueType& innerProdRes = runtime.mInnerProdRes;
    ValueType alpha;
    ValueType beta;

    const RealType<ValueType>& eps = runtime.mEps;

    RealType<ValueType>& normRes = runtime.mNormRes;

    vecT = A * vecPT;

    ValueType innerProduct = res0.dotProduct( vecT );

    if ( normRes < eps || common::Math::abs( innerProduct ) < eps )  // innerProduct is small
    {
        alpha = 0.0;
    }
    else
    {
        alpha = innerProdRes / innerProduct;
    }

    SCAI_LOG_INFO( logger, "vecQ = vecU - " << alpha << " vecT, normRes = " << normRes
                   << ", innerProdRes = " << innerProdRes << ", innerProduct = " << innerProduct )

    vecQ = vecU - alpha * vecT;

    // PRECONDITIONING
    if ( mPreconditioner != NULL )
    {
        // vecUT = 0;  here the more general approach
        vecUT.setSameValue( mPreconditioner->getCoefficients().getColDistributionPtr(), 0 );
        vecTemp = vecU + vecQ;
        mPreconditioner->solve( vecUT, vecTemp );
    }
    else
    {
        vecUT  = vecU + vecQ;
    }

    solution = solution + alpha * vecUT;
    ValueType innerProdResOld = innerProdRes;
    res = res - alpha * A * vecUT;
    innerProdRes = res0.dotProduct( res );
    normRes = l2Norm( res );

    if ( normRes < eps || common::Math::abs( innerProdResOld ) < eps )            // innerProdResOld is small
    {
        beta = 0.0;
    }
    else
    {
        beta = innerProdRes / innerProdResOld ;
    }

    SCAI_LOG_INFO( logger, "beta = " << beta )

    vecU = res + beta * vecQ;
    vecP = vecU + beta * beta * vecP;
    vecP = vecP + beta * vecQ;

    // PRECONDITIONING
    if ( mPreconditioner != NULL )
    {
        vecPT.setSameValue( mPreconditioner->getCoefficients().getColDistributionPtr(), 0 );
        mPreconditioner->solve( vecPT, vecP );
    }
    else
    {
        vecPT = vecP ;
    }

    // End Implementation

    mCGSRuntime.mSolution.setDirty( false );
}

/* ========================================================================= */
/*       Getter runtime                                                      */
/* ========================================================================= */

template<typename ValueType>
typename CGS<ValueType>::CGSRuntime& CGS<ValueType>::getRuntime()
{
    return mCGSRuntime;
}

template<typename ValueType>
const typename CGS<ValueType>::CGSRuntime& CGS<ValueType>::getRuntime() const
{
    return mCGSRuntime;
}

/* ========================================================================= */
/*       virtual methods                                                     */
/* ========================================================================= */

template<typename ValueType>
CGS<ValueType>* CGS<ValueType>::copy()
{
    return new CGS<ValueType>( *this );
}

template<typename ValueType>
void CGS<ValueType>::writeAt( std::ostream& stream ) const
{
    const char* typeId = common::TypeTraits<ValueType>::id();

    stream << "CGS<" << typeId << "> ( id = " << this->getId() 
           << ", #iter = " << getRuntime().mIterations << " )";
}

/* ========================================================================= */
/*       Template instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( CGS, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace solver */

} /* end namespace scai */
